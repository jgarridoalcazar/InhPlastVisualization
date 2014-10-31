'''
Created on May 12, 2014

@author: Jesus Garrido (jgarridoalcazar at gmail.com)
'''

import nest
import abc
from CerebellarModel import CerebellarModel
import shutil
import numpy
import os
import logging
from SpikingSimulation.Utils.Utils import WriteConfigFile
from mpi4py import MPI

logger = logging.getLogger('Simulation')

class NestCerebellarModel(CerebellarModel):
    '''
    This class defines an inherited class including all the methods
    needed in order to generate a cerebellum-like network for EDLUT
    simulator.
    '''
    __metaclass__ = abc.ABCMeta
    
    # Cell name translation
    cellNameTranslatorDict = {
         'ConductanceLIF' : 'iaf_cond_exp',
         'ConductanceLIFwIP': 'iaf_cond_exp_ip',
         'ConductanceLIFStowIP': 'iaf_cond_exp_sto_ip',
         'CurrentLIF' : 'iaf_neuron'
    }

    # This dictionary maps the layer configuration parameters into the NEST cell model parameters (used in the keys).
    paramTranslatorDict = { 
        'E_ex':'eexc*1.e3', # Excitatory reversal potential in mV
        'E_in':'einh*1.e3', # Inhibitory reversal potential in mV
        'E_L': 'erest*1.e3', # Resting potential in mV
        'V_reset': 'erest*1.e3', # Reset potential in mV (map to the resting potential)
        'tau_syn_ex': 'texc*1.e3', # Decay time of the excitatory synaptic conductance in ms
        'tau_syn_in': 'tinh*1.e3', # Decay time of the inhibitory synaptic conductance in ms
        'g_L': 'grest*1.e9', # Resting conductance in nS
        'V_th': 'eth*1.e3', # Threshold potential in mV
        't_ref': 'tref*1.e3', # Refractory period in ms
        'C_m': 'cm*1.e12', # Membrane capacitance in pF
        'tau_m': 'cm/grest*1.e3', # Membrane time constant in ms
        'tau_minus':'tau_minus*1.e3', # Time constant of the post-pre part in ms
        # Intrinsic plasticity parameters
        'tau_ip':'tau_ip*1.e3', # Time constant of the intrinsic plasticity in ms
        'beta':'beta_ip*1.', # Beta parameter of the IP (unitless)
        'epsilon_rC': 'epsilon_rc_ip*1.', # epsilon_rC parameter of the IP (unitless)
        'epsilon_rR': 'epsilon_rr_ip*1.', # epsilon_rR parameter of the IP (unitless)
        'r_C': '1./(cm*1.e12)', # Inverse of the membrane capacitance
        # Stochastic IP model parameters
        'ip_rate': 'ip_rate', # IP rate
        'target_firing': 'target_freq', # Target firing frequency
        # Symmetric STDP parameters
        'tau_sym': 'tau_sym*1.e3'
    }
    
    # This dictionary maps the state variables as used in the config file with the state variable names in NEST.
    stateTranslatorDict = {
        'Vm': ['V_m',1.e-3], # Membrane potential
        'Gexc': ['g_ex',1.e-9], # Excitatory conductance
        'Ginh': ['g_in',1.e-9], # Inhibitory conductance
        'rC': ['r_C',1.], # Excitatory conductance
        'gL': ['g_L',1.], # Excitatory conductance
        'Vth': ['V_th',1.e-3], # Threshold potential
        'r0': ['r_0',1.], # Gain frequency
        'ualpha': ['u_alpha',1.e-3], # Alpha parameter
        'refractoriness': ['refractoriness',1.], # Refreactoriness
        'gain': ['gain',1.], # Current frequency
        'firing_probability': ['firing_probability',1.] # Firing probability
    }
    
    # Learning rule translation
    ruleNameTranslatorDict = {
        'STDP' : 'stdp_synapse_hom',
        'STDPSym'   :   'stdp_sym_synapse_hom'
    }
    
    # This dictionary maps the learning rule configuration parameters into the NEST learning rule model parameters (used in the keys).
    ruleTranslatorDict = { 
        'tau_plus':'tau_plus*1.e3', # Time constant of the pre-post part in ms
        'lambda': 'learning_step', # pre-post amplitude (normalized?)
        'alpha': 'minus_plus_ratio', # post-pre/pre-post ratio (normalized?)
        'Wmax': 'max_weight*1e9', # Maximum weight in nS
        'mu_plus': '0.0', # Exponent of the weight dependence (0 to get additive STDP, 1 for multiplicative STDP)
        'mu_minus': '0.0', # Exponent of the weight dependence (0 to get additive STDP, 1 for multiplicative STDP)
        'tau_sym':'tau_sym*1.e3', # Time constant of the symmetric STDP in ms
    
    }
    
    # Connectivity pattern translator
    connectivityNameTranslatorDict = {
        'randomn2one': ('fixed_indegree',['indegree','autapses']),
        'random_with_probability': ('pairwise_bernoulli',['p','autapses'])
    }
    
    # This dictionary maps the connectivity configuration parameters into the NEST connectivity parameters (used in the keys).
    connectivityTranslatorDict = { 
        'indegree':'number_of_source_cells', # Target cell fan-in
        'p': 'connection_probability*1.', # Probability of connection
        'autapses': 'allow_auto_connection' # Connections from one node to itself are allowed
    }
    
    # Weight initialization translator
    distributionNameTranslatorDict = {
        'random': ('uniform',['low','high']),
    }
    
    # This dictionary maps the distribution configuration parameters into the NEST distribution parameters (used in the keys).
    distributionTranslatorDict = { 
        'low':'random_min_weight*1.e9', # Minimum weight
        'high': 'random_max_weight*1.e9', # Maximum weight
    }
    
    
    def _map_cm_parameters(self,neuron_layer):
        '''
        Calculate all the NEST cell model parameters from the values defined in the cfg file. 
        '''
        
        cm_parameters_out=dict()
        
        # Get NEST parameters of the cell model
        cell_model_dict = nest.GetDefaults(neuron_layer.nest_model_name)
        
        for key in cell_model_dict:
            if key in self.paramTranslatorDict:
                try:
                    # Evaluate the expression of the parameter   
                    v_out = eval(self.paramTranslatorDict[key],None, neuron_layer.cell_model_parameters)
                    cm_parameters_out[key] = v_out
                except NameError as e:
                    logger.warning('%s cannot be calculated in layer %s. Variable %s is not defined. Using default value %s', key, neuron_layer.__name__, e.message.split("'")[1], cell_model_dict[key])
                                    
        return cm_parameters_out
    
    def _map_lr_parameters(self,synapsis_layer):
        '''
        Calculate all the NEST synapsis parameters from the values defined in the cfg file. 
        '''
        
        lr_parameters_out=dict()
        
        # Get NEST parameters of the cell model
        rule_model_dict = nest.GetDefaults(synapsis_layer.nest_model_name)
        
        for key in rule_model_dict:
            if key in self.ruleTranslatorDict:
                try:
                    # Evaluate the expression of the parameter   
                    v_out = eval(self.ruleTranslatorDict[key],None, synapsis_layer.learning_rule_parameters)
                    lr_parameters_out[key] = v_out
                except NameError as e:
                    logger.warning('%s cannot be calculated in layer %s. Variable %s is not defined. Using default value %s', key, synapsis_layer.__name__, e.message.split("'")[1], rule_model_dict[key])
                                    
        return lr_parameters_out
    
    
        
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates a new cerebellar model.
        @param config_dict: Dictionary including all the network generation parameters.
        '''
        
        super(NestCerebellarModel, self).__init__(**kwargs)
        
    def initialize_simulation(self):
        '''
        Reset the NEST simulator, build the network and set additional parameters (e.g. the simulation time-step, the number of threads, ...).     
        '''
        
        super(NestCerebellarModel, self).initialize_simulation()
        
        try:
            nest.Install('glplasticitymodule')
        except nest.NESTError:
            print('NEST Error caught on loading user module. Retrying...')
            nest.Install('glplasticitymodule')
        
        nest.sr("M_WARNING setverbosity")
        
        nest.ResetKernel()
        
        self.nest_options = self.config_dict['nest']
        
        nest_options_dict = dict()
        
        if 'num_record_processes' in self.nest_options:
            if self.nest_options['num_record_processes']>0:
                if nest.NumProcesses()>self.nest_options['num_record_processes']:
                    logger.debug('Setting number of record processes to %s',self.nest_options['num_record_processes'])
                    nest.SetNumRecProcesses(self.nest_options['num_record_processes'])
                else:
                    logger.warning('Invalid number of record processes: %s. Ignoring',self.nest_options['num_record_processes'])
               
        if 'number_of_virtual_processes' in self.nest_options:
            nest_options_dict['total_num_virtual_procs'] = self.nest_options['number_of_virtual_processes']
        
        if 'resolution' in self.nest_options:
            nest_options_dict['resolution'] = self.nest_options['resolution']*1e3
        
        nest_options_dict['overwrite_files'] = True # set to True to permit overwriting
        nest_options_dict['print_time'] = False # set to True to print the simulation completed
        
        if not 'record_to_file' in self.simulation_options:
            self.simulation_options['record_to_file'] = False
        
        if self.simulation_options['record_to_file']:        
            if 'data_path' in self.simulation_options:
                nest_options_dict['data_path'] = self.simulation_options['data_path']
            else:
                nest_options_dict['data_path'] = './results'
                
            if 'simulation_name' in self.simulation_options:
                nest_options_dict['data_path'] = nest_options_dict['data_path'] + '/' + self.simulation_options['simulation_name']
                
            data_path = nest_options_dict['data_path']
            
            if self.get_my_process_id()==0 and self.simulation_options['record_to_file']:
                if not os.path.exists(nest_options_dict['data_path']):
                    os.makedirs(nest_options_dict['data_path'])
                else:
                    logger.debug('Removing results folder %s', data_path)
                    try:
                        shutil.rmtree(data_path)
                        os.makedirs(data_path)
                    except OSError:
                        pass
                    
                     
        # Synchronize all the processes to avoid start writting before deleting the folder
        comm = MPI.COMM_WORLD
        
        comm.Barrier()   
            
        
        nest.SetKernelStatus(nest_options_dict)
        
        # Set random seeds
        nest_options_dict = dict()
        nest_options_dict['grng_seed'] = self.simulation_options['seed'] + self.get_number_of_virtual_processes()
        nest_options_dict['rng_seeds'] = range(self.simulation_options['seed'] + self.get_number_of_virtual_processes() + 1,
                                               self.simulation_options['seed'] + 2*self.get_number_of_virtual_processes() + 1)
        nest.SetKernelStatus(nest_options_dict)
        
        # Set record_to_file in spike_detectors to the value    
        nest.SetDefaults('spike_detector', 'to_file', self.simulation_options['record_to_file'])
        nest.SetDefaults('multimeter', 'to_file', self.simulation_options['record_to_file'])
        
        # Initialize ac_current and dc_current stimulation
        self.ac_generator = []
        self.dc_current = []
        
        # Build the network
        self._build_network()
                
        self.simulation_time = 0
        
        if self.simulation_options['record_to_file']:
            # Copy configuration file of this simulation
            WriteConfigFile(self.config_dict, data_path+'/'+'SimulationConfig.cfg')
        
        return
    
    
    def _initialize_weight_recording_buffer(self):
        
        # Check if recording_time_step is above 0
        if self.simulation_options['weight_recording_step'] >= float("inf"):
            return
        
        if 'data_path' in self.simulation_options:
            data_path = self.simulation_options['data_path']
        else:
            data_path = './results'
            
        if 'simulation_name' in self.simulation_options:
            data_path = data_path + '/' + self.simulation_options['simulation_name']

        
        for layer in self.synaptic_layers:
            if layer.weight_recording:
                layer.weight_record = dict()
                
                # logger.debug('Source layer: %s. Target layer: %s',layer.source_layer.nest_layer,layer.target_layer.nest_layer)
                # Store the source and target cells to know the order to be stored
                # We assume GetConnections return only those connections whose target node is local to the process
                layer.weight_record['con'] = nest.GetConnections(source=layer.source_layer.nest_layer,target=layer.target_layer.nest_layer)
                global_connections = numpy.array(nest.GetStatus(layer.weight_record['con'],['source','target']))
                min_source = layer.source_layer.MinIndex
                min_target = layer.target_layer.MinIndex
                layer.weight_record['connections'] = numpy.array([(val[0]-min_source,val[1]-min_target) for val in global_connections])
                
                # logger.debug('Connections created in the weight record of layer %s: %s', layer.__name__, layer.weight_record['connections'])
                
                # Record the initial weights
                scaled_weights = [weight*1.e-9 for weight in nest.GetStatus(layer.weight_record['con'],'weight')]
                layer.weight_record['weights'] = numpy.array([scaled_weights])
                layer.weight_record['time'] = numpy.array([0])
                
                # logger.debug('Initial weights in the weight record of layer %s: %s', layer.__name__, layer.weight_record['weights'])
                
                # Open the file when the weights are going to be stored
                if self.simulation_options['record_to_file']:
                    file_name = data_path + '/' + layer.__name__ + '_weights_' + str(self.get_my_process_id()) + '.dat'
                    layer.weight_record['f_handle'] =  open(file_name,'w')
                    
                    # Write the connections in the first line
                    numpy.savetxt(layer.weight_record['f_handle'], layer.weight_record['connections'], fmt='%1.1u', delimiter=' ', newline=' ')
                    layer.weight_record['f_handle'].write('\n')
                    # Write the initial weights and time
                    numpy.savetxt(layer.weight_record['f_handle'], numpy.append(layer.weight_record['time'][-1:],layer.weight_record['weights'][-1:]), fmt='%1.3e', delimiter=' ', newline=' ')
                    layer.weight_record['f_handle'].write('\n')
                    # Store the IO buffer
                    layer.weight_record['f_handle'].flush()
                    
                
                # print 'Process',self.get_my_process_id(),':','Weight record:',layer.weight_record
            else:
                layer.weight_recording = None

        self.next_weight_step = 1
        
        
        return
    
    def _build_network(self):
        '''
        Generate a NEST model based on the inherited cerebellar model.
        '''
        
        super(NestCerebellarModel, self)._build_network()
        
        # Create nodes in the network
        self._create_nodes()
        
        # Create the connections in the network
        self._create_connections()
        
        # print 'Process:', self.get_my_process_id(),'Network created:',nest.GetStatus(nest.GetConnections(), ['source', 'target', 'weight'])
        
        # print 'Process:', self.get_my_process_id(),'Current netword:',nest.PrintNetwork()
        
        # Initialize weight recording buffer
        self._initialize_weight_recording_buffer()
        
        return
    
    def _create_nodes(self):
        '''
        Generate the nodes of every neuron layer in the model
        '''
        
        for layer in self.neuron_layers:
            # Check if the specified cell model is included
            if layer.cell_model in self.cellNameTranslatorDict:
                layer.nest_model_name = self.cellNameTranslatorDict[layer.cell_model]
            else:
                logger.error('The cell model %s has not been mapped to a NEST model', layer.cell_model)
                raise Exception('Non-MappedCellModel')
        
            # Get the dictionary with the mapped values
            nest_param_dict = self._map_cm_parameters(layer)
            
            # Create the layer nodes and store them in the NuronLayer object
            layer.nest_layer = nest.Create(model=layer.nest_model_name, params=nest_param_dict, n=layer.number_of_neurons)
            
            # Check whether we have to record the activity. If that is the case, create the spike detector
            if layer.register_activity:
                layer.nest_spike_detector = nest.Create(model='spike_detector')
                #nest.ConvergentConnect(layer.nest_layer,layer.nest_spike_detector)
                nest.Connect(layer.nest_layer,layer.nest_spike_detector,conn_spec='all_to_all')
                
                # Set layer registered name
                nest.SetStatus(layer.nest_spike_detector,{'label':layer.__name__+'_spike'})
            else:
                layer.nest_spike_detector = None
                
            # Check wether register any state variable (Vm, gexc, ginh, ...)
            recording_vars = []
            recorded_vars = []
            
            # Get NEST model recordable variables
            available_vars = nest.GetDefaults(layer.nest_model_name)['recordables']
            if (layer.record_vars):
                # If only an string has been used, embed it in an array
                if isinstance(layer.record_vars,str):
                    layer.record_vars = [layer.record_vars]
                    
                # Get recording variables in this cell model
                for var in layer.record_vars:
                    if var in self.stateTranslatorDict:
                        # Check if this variable is recordable in this cell model for NEST
                        if self.stateTranslatorDict[var][0] in available_vars:
                            recording_vars.append(self.stateTranslatorDict[var][0])
                            recorded_vars.append(var)
                        else:
                            logger.warning('%s variable is not recordable in the model %s. Ignoring',self.stateTranslatorDict[var][0],layer.nest_model_name)
                    else:
                        logger.warning('%s state variable is included in the recordable variable map. Ignoring',var)
                
            if recording_vars:
                layer.nest_multimeter = nest.Create(model='multimeter', params = {'withtime': True,
                                                                            'withgid': True, 
                                                                            'interval': layer.record_step*1e3,
                                                                            'record_from': recording_vars})
                #nest.DivergentConnect(layer.nest_multimeter,layer.nest_layer)
                nest.Connect(layer.nest_multimeter,layer.nest_layer,conn_spec='all_to_all')
                
                # Set layer registered name
                nest.SetStatus(layer.nest_multimeter,{'label':layer.__name__+'_mult'})
            else:
                layer.nest_multimeter = None
                
            comm = MPI.COMM_WORLD
            
            # Collect the minimum of the layer indexes
            localMin = numpy.array(numpy.min(layer.nest_layer), dtype='l') 
            layer.MinIndex = numpy.array(0, dtype='l') 
            comm.Allreduce([localMin, MPI.LONG], [layer.MinIndex, MPI.LONG], op=MPI.MIN)
            
            self.config_dict[layer.__name__]['minindex'] = layer.MinIndex
            # Update the record_vars variable to those that effectively can be recorded
            self.config_dict[layer.__name__]['record_vars'] = ','.join(recorded_vars) 
            
            # print 'Process:', self.get_my_process_id(),'Layer:', layer.__name__, 'Collected:', layer.MinIndex
            
            pass

            
        return
    
    def _create_connections(self):
        '''
        Generate the connections of every synaptic layer in the model
        '''
        
        for layer in self.synaptic_layers:
            
            # Check if the specified learning rule is included
            if layer.learning_rule_type:
                if layer.learning_rule_type in self.ruleNameTranslatorDict:
                    layer.nest_model_name = self.ruleNameTranslatorDict[layer.learning_rule_type]
                else:
                    logger.error('The synapsis model %s has not been mapped to a NEST model',layer.learning_rule_type)
                    raise Exception('Non-MappedSynapsisModel')
        
                # Get the dictionary with the mapped values
                nest_param_dict = self._map_lr_parameters(layer)
                nest_param_dict['delay'] = layer.synaptic_delay*1.e3
                #nest_param_dict['min_delay'] = layer.synaptic_delay*1.e3
                #nest_param_dict['max_delay'] = layer.synaptic_delay*1.e3
                
                nest.CopyModel(layer.nest_model_name,layer.__name__,nest_param_dict)
            else:
                # Static synaptic layer
                nest.CopyModel("static_synapse", layer.__name__, {'delay':layer.synaptic_delay*1.e3})
            
            # Search the local nodes
            #node_info = nest.GetStatus(layer.target_layer.nest_layer)
            
            # Create the synaptic connections in NEST
            #presynaptic = [layer.source_layer.nest_layer[layer.source_index[index]] for index in xrange(len(layer.source_index)) if node_info[layer.target_index[index]]['local']] # Get absolute indexes in NEST
            #postsynaptic = [layer.target_layer.nest_layer[layer.target_index[index]] for index in xrange(len(layer.target_index)) if node_info[layer.target_index[index]]['local']] # Get absolute indexes in NEST
            #weights = [layer.weights[index]*1.e9 for index in xrange(len(layer.target_index)) if node_info[layer.target_index[index]]['local']] # Get absolute indexes in NEST
            #delay_values = [layer.synaptic_delay*1.e3] * len(presynaptic)
            
            # Create the synaptic connections in NEST >2.4
            con_dict = self._create_connection_pattern_dict(layer)
            syn_dict = self._create_connection_synapsis_dict(layer)
            
            #print 'Process:', self.get_my_process_id(),'Source cells:', presynaptic,'Target cells:',postsynaptic
            nest.Connect(pre=layer.source_layer.nest_layer, post=layer.target_layer.nest_layer, conn_spec=con_dict, syn_spec=syn_dict)
            #nest.OneToOneConnect(pre=presynaptic, post=postsynaptic, params= weights, delay=delay_values, model=layer.__name__)
            logger.debug('Nest Process: %s. Connections created in layer %s: %s', nest.Rank(), layer.__name__,len(nest.GetConnections(source=layer.source_layer.nest_layer,target=layer.target_layer.nest_layer)))


 
        return
    
    def _create_connection_pattern_dict(self, layer):
        # Check if connectivity_type is defined
        if layer.connectivity_type in self.connectivityNameTranslatorDict:
            dictionary = {'rule':self.connectivityNameTranslatorDict[layer.connectivity_type][0]}

            for param in self.connectivityNameTranslatorDict[layer.connectivity_type][1]:
                dictionary[param] = eval(self.connectivityTranslatorDict[param],None,layer.connectivity_parameters)
                
        else:
            logger.error('Non-mapped connectivity pattern %s in layer %s', layer.connectivity_type, layer.__name__)
            raise Exception('Non-MappedConnectivityPattern')
            
        return dictionary
    
    def _create_connection_synapsis_dict(self, layer):
        # Check if weight_initialization_type is defined
        if layer.weight_initialization_type=='fixed':
            dictionary = {'model':layer.__name__,
                          'weight': layer.weight_initialization_parameters['initial_weight']*1.e9}
            return dictionary;
        elif layer.weight_initialization_type in self.distributionNameTranslatorDict:
            dictionary = {'model':layer.__name__}
            weight_dict = {'distribution':self.distributionNameTranslatorDict[layer.weight_initialization_type][0]}

            for param in self.distributionNameTranslatorDict[layer.weight_initialization_type][1]:
                weight_dict[param] = eval(self.distributionTranslatorDict[param],None,layer.weight_initialization_parameters)
                
            dictionary['weight'] = weight_dict
                
        else:
            logger.error('Non-mapped weight initialization distribution %s in layer %s', layer.weight_initialization_type, layer.__name__)
            raise Exception('Non-MappedWeightInitialization')
        
        return dictionary
            
             
    def add_ac_current(self, **kwargs):
        '''
        Set a new ac current that will be conveyed to the cerebellar mossy fibers according to the parameters.
        @param amplitude Amplitude of the wave (in A)
        @param offset Constant amplitude offset (in pA)
        @param frequency Frequency of the ac_generator (in Hz)
        @param phase Phase of the sine current (0-360 deg)
        '''
        iThreshold = (self.mflayer.cell_model_parameters['eth']-self.mflayer.cell_model_parameters['erest'])*self.mflayer.cell_model_parameters['grest']
        
        # Default values of each parameter in NEST units
        ac_dict = {'amplitude' : 1.*iThreshold, # Amplitude in pA
                   'offset' : 0., # Offset in pA
                   'phase' : 0., # Phase in deg
                   'frequency' : 1. # Frequency in Hz
                   }
        
        if 'amplitude' in kwargs:
            ac_dict['amplitude'] = kwargs.pop('amplitude')*iThreshold*1.e12
            
        if 'offset' in kwargs:
            ac_dict['offset'] = kwargs.pop('offset')*1.e12
            
        if 'phase' in kwargs:
            ac_dict['phase'] = float(kwargs.pop('phase'))
            
        if 'frequency' in kwargs:
            ac_dict['frequency'] = float(kwargs.pop('frequency'))
        
        
        new_ac_generator = nest.Create(model='ac_generator', n=1, params=ac_dict)
        
        # Connect the new ac_generator to the mossy fibers
        #nest.DivergentConnect(pre=new_ac_generator, post=self.mflayer.nest_layer, weight=1., delay=1., model='static_synapse')
        nest.Connect(pre=new_ac_generator, post=self.mflayer.nest_layer, conn_spec='all_to_all', syn_spec={'model':'static_synapse','weight':1.,'delay':1.})
        
        self.ac_generator.append(new_ac_generator)
        
        return
    
    def set_dc_current(self, **kwargs):
        '''
        Set the mossy fiber dc currents to the values specifed in the parameter.
        @param amplitude Amplitude of the current in A. A list with the same length of the number of mossy fibers must be used.
        '''
        
        # If dc_current generators have not been 
        if not self.dc_current:            
            self.dc_current = nest.Create(model='dc_generator', n=self.mflayer.number_of_neurons, params={'amplitude':0.0})
            #weights = [1.]*self.mflayer.number_of_neurons
            #delay = [1.]*self.mflayer.number_of_neurons
            #nest.Connect(pre=self.dc_current, post=self.mflayer.nest_layer, params=weights, delay=delay, model='static_synapse')
            nest.Connect(self.dc_current, self.mflayer.nest_layer, conn_spec='one_to_one', syn_spec={"model": "static_synapse", "weight":1., "delay":1.})
        
        # Default values of each parameter in NEST units
        if 'amplitude' in kwargs:
            iThreshold = (self.mflayer.cell_model_parameters['eth']-self.mflayer.cell_model_parameters['erest'])*self.mflayer.cell_model_parameters['grest']
            
            amp = kwargs.pop('amplitude')*iThreshold
            
            if len(amp)==self.mflayer.number_of_neurons:
                amplitude = amp*1.e12
                nest.SetStatus(nodes=self.dc_current, params='amplitude', val=amplitude)
                # We could check and set status only of the local nodes   
            else:
                logger.error('dc_current amplitude has to be a list with the same length of the number of MFs')
                raise Exception('InvalidDCCurrent')
        else:
            logger.error('Error: dc_current amplitude has to be specified')
            raise Exception('Non-SpecifiedDCCurrent')
        
        return
        
    def _save_weights(self):
        
        # Check if recording_time_step is above 0
        for layer in self.synaptic_layers:
            if layer.weight_recording:
                # Record the weights
                scaled_weights = [weight*1.e-9 for weight in nest.GetStatus(layer.weight_record['con'],'weight')]
                layer.weight_record['weights'] = numpy.append(layer.weight_record['weights'],[scaled_weights], axis=0)
                layer.weight_record['time'] = numpy.append(layer.weight_record['time'],self.next_weight_step * self.simulation_options['weight_recording_step'])
                
                # Open the file when the weights are going to be stored
                if self.simulation_options['record_to_file']:
                    # Write the initial weights and time
                    numpy.savetxt(layer.weight_record['f_handle'], numpy.append(layer.weight_record['time'][-1:],layer.weight_record['weights'][-1:]), fmt='%1.3e', delimiter=' ', newline=' ')
                    layer.weight_record['f_handle'].write('\n')
                    # Store the IO buffer
                    layer.weight_record['f_handle'].flush()
                    
                

        self.next_weight_step += 1
    
    def simulate_network(self,time):
        '''
        Simulate the network for the specified time (in seconds).
        @param time Length to be simulated (in seconds)
        '''
        
        end_time = self.simulation_time + time
        next_weight_recording = self.next_weight_step * self.simulation_options['weight_recording_step'] 
        
        # Simulate until recording weight
        while (end_time-self.simulation_time>=self.nest_options['resolution']):
            next_stop = min(next_weight_recording,end_time)
            sim_time = next_stop-self.simulation_time
            
            # Simulate the step (in ms)
            # Round the simulation time to avoid inconsistent results in nest.
            # nest.Simulate(math.ceil(sim_time*1.e3))
            nest.Simulate(sim_time*1.e3)
            
            # If it is time to record the weights
            if next_stop==next_weight_recording:
                self._save_weights()
                next_weight_recording = self.next_weight_step * self.simulation_options['weight_recording_step'] 
            
            # Get the simulation time from the kernel to avoid inconsistencies between both times
            self.simulation_time = nest.GetKernelStatus('time')/1.e3
            
        return
    
            
    def get_spike_activity(self, **kwargs):
        '''
        Get the spikes been fired in the cells, layer and time specified in the parameters.
        @param neuron_layer Layer name as the section name in the config file (.cfg).
        @param neuron_indexes Indexes of the cell to get the activity. Relative indexes in the layer.
        @param init_time Initial time from which retrieve the activity.
        @param end_time Final time until which retrieve the activity.  
        '''
        
        # Collect all the parameters
        if 'neuron_layer' in kwargs:
            neuron_layer_name = kwargs.pop('neuron_layer')
        else:
            logger.error('Non specified neuron layer in get_spike_activity function')
            raise Exception('Non-DefinedNeuronLayer')
        
        
        if neuron_layer_name in self.layer_map:
            neuron_layer = self.layer_map[neuron_layer_name]
        else:
            logger.error('Invalid neuron layer in get_spike_activity function')
            raise Exception('InvalidNeuronLayer')
        
        if 'neuron_indexes' in kwargs:
            local_indexes = kwargs.pop('neuron_indexes')
            
            if (max(local_indexes)>=(neuron_layer.MinIndex+neuron_layer.number_of_neurons)):
                logger.error('Invalid neuron index in get_spike_activity function')
                raise Exception('InvalidNeuronIndex')
            
            neuron_indexes = [neuron_layer.MinIndex+index for index in local_indexes]            
        else:
            neuron_indexes = neuron_layer.nest_layer
            
        if 'init_time' in kwargs:
            init_time = kwargs.pop('init_time')
        else:
            init_time = 0
            
        if 'end_time' in kwargs:
            end_time = kwargs.pop('end_time')
        else:
            end_time = float('inf')
        
        # Comprobar que la neurona (spike recorder) sea local
        if not neuron_layer.nest_spike_detector: 
            logger.error('Invalid neuron layer in get_spike_activity function. The activity in this layer has not been recorded')
            raise Exception('InvalidNeuronLayer')
        
        # If the spike detector is local to this process
        # logger.debug('Checking locality in layer %s from %s: %s',neuron_layer.__name__, neuron_layer.nest_spike_detector, nest.GetStatus(neuron_layer.nest_spike_detector,'local'))
        if nest.GetStatus(neuron_layer.nest_spike_detector,'local')[0]:
            # logger.debug('Getting spikes in layer %s from %s: %s',neuron_layer.__name__, neuron_layer.nest_spike_detector, nest.GetStatus(neuron_layer.nest_spike_detector,'events'))
            spike_events = nest.GetStatus(neuron_layer.nest_spike_detector,'events')[0]
            time = spike_events['times']*1e-3
            neuron_id = spike_events['senders']
            
            # Select only those events in the specified interval
            index = (time>=init_time) & (time<=end_time) 
            time = time[index]
            neuron_id = neuron_id[index]
            
            # Select only those events in the specified cells
            index = numpy.in1d(neuron_id, neuron_indexes)
            time = time[index]
            neuron_id = neuron_id[index]
        else:
            time = numpy.array([],dtype=numpy.float64)
            neuron_id = numpy.array([],dtype=numpy.float64)
        
        # Send relative neuron_id
        neuron_id = neuron_id-neuron_layer.MinIndex
        
        comm = MPI.COMM_WORLD
        
        # Send the number of elements
        lsum = numpy.array([len(time)], dtype=numpy.uint64)
        if self.get_my_process_id()==0:
            num_spikes = numpy.empty(comm.Get_size(), dtype=numpy.uint64)
        else:
            num_spikes = None
        comm.Gather([lsum, MPI.UNSIGNED_LONG], [num_spikes, MPI.UNSIGNED_LONG], root=0) 
        
        # print 'Process',self.get_my_process_id(),':','Sent number:',lsum,'Collected numbers ->',num_spikes
        
        if self.get_my_process_id()==0:
            # Total number of spikes to be gathered
            gsum = num_spikes.sum()
            num_sent = tuple(num_spikes)
            aux_array = numpy.roll(num_spikes.cumsum(),1)
            aux_array[0] = 0
            offset = tuple(aux_array)
            
            gtime = numpy.empty(gsum, dtype=numpy.float64)
            gneuron_id = numpy.empty(gsum, dtype=numpy.uint64)
        else:
            num_sent = None
            offset = None
            gtime = None
            gneuron_id = None
        
        # Gather the time and neuron_id arrays
        comm.Gatherv(time, [gtime, num_sent, offset, MPI.DOUBLE], root=0)
        comm.Gatherv(neuron_id, [gneuron_id, num_sent, offset, MPI.UNSIGNED_LONG], root=0)
        
#         # Calculate the average activity
#         av_activity = float(gneuron_id.shape[0])/(len(neuron_indexes)*(end_time-init_time))
#         logger.debug('Average firing rate in layer %s: %s',neuron_layer_name, av_activity)
#         av_activity_per_cell = [numpy.count_nonzero(gneuron_id==num_cell)/(8.*(end_time-init_time)) for num_cell in (neuron_indexes-neuron_layer.MinIndex)]
#         logger.debug('Average spikes per cycle: %s', av_activity_per_cell)
#         bins = numpy.arange(-0,4.5)
#         av_activity_histogram = numpy.histogram(av_activity_per_cell,bins)[0]
#         logger.debug('Histogram of spikes per cycle: %s', av_activity_histogram)
        
        
        # print 'Process',self.get_my_process_id(),':','Collected time ->',gtime,'Collected Neurons ->',gneuron_id
        return (gtime,gneuron_id)
    
    def get_state_variable(self, **kwargs):
        '''
        Get the state variables in the cells, layer and during the time specified in the parameters.
        @param neuron_layer Layer name as the section name in the config file (.cfg).
        @param neuron_indexes Indexes of the cell to get the activity.
        @param init_time Initial time from which retrieve the activity.
        @param end_time Final time until which retrieve the activity.
        @param state_variable Name of the state variable to retrieve.   
        '''
        # Collect all the parameters
        if 'neuron_layer' in kwargs:
            neuron_layer_name = kwargs.pop('neuron_layer')
        else:
            logger.error('Non specified neuron layer in get_state_variable function')
            raise Exception('Non-DefinedNeuronLayer')
        
        
        if neuron_layer_name in self.layer_map:
            neuron_layer = self.layer_map[neuron_layer_name]
        else:
            logger.error('Invalid neuron layer in get_state_variable function')
            raise Exception('InvalidNeuronLayer')
        
        if 'neuron_indexes' in kwargs:
            local_indexes = kwargs.pop('neuron_indexes')
            
            if (max(local_indexes)>=(neuron_layer.MinIndex+neuron_layer.number_of_neurons)):
                logger.error('Invalid neuron index in get_state_variable function')
                raise Exception('InvalidNeuronIndex')
            
            neuron_indexes = [neuron_layer.MinIndex+index for index in local_indexes]         
        else:
            neuron_indexes = neuron_layer.nest_layer
            
        if 'state_variable' in kwargs:
            variable_name = kwargs.pop('state_variable')
        else:
            logger.error('Non specified state_variable in get_state_variable function')
            raise Exception('Non-DefinedStateVariable')
        
        if variable_name in self.stateTranslatorDict:
            state_variable_name = self.stateTranslatorDict[variable_name][0]
        else:
            logger.error('Invalid state variable in get_state_variable function')
            raise Exception('InvalidStateVariable')
        
        if not variable_name in neuron_layer.record_vars:
            logger.error('Invalid state variable in get_state_variable function. %s has not been recorded', variable_name)
            raise Exception('InvalidStateVariable')
        
        if 'init_time' in kwargs:
            init_time = kwargs.pop('init_time')
        else:
            init_time = 0
            
        if 'end_time' in kwargs:
            end_time = kwargs.pop('end_time')
        else:
            end_time = float('inf')
        
        # If the spike detector is local to this process
        recording_events = nest.GetStatus(neuron_layer.nest_multimeter,'events')[0]
        # print 'Process',self.get_my_process_id(),': recording_events=',recording_events
        time = recording_events['times']*1e-3
        neuron_id = recording_events['senders']
        value = recording_events[state_variable_name]*self.stateTranslatorDict[variable_name][1]
        
        # Select only those events in the specified interval
        index = (time>=init_time) & (time<=end_time) 
        time = time[index]
        neuron_id = neuron_id[index]
        value = value[index]
        
        # Select only those events in the specified cells
        index = numpy.in1d(neuron_id, neuron_indexes)
        time = time[index]
        neuron_id = neuron_id[index]
        value = value[index]
        
        # Send relative neuron_id
        neuron_id = neuron_id-neuron_layer.MinIndex
        
        comm = MPI.COMM_WORLD
        
        # Send the number of elements
        lsum = numpy.array([len(time)], dtype=numpy.uint64)
        if self.get_my_process_id()==0:
            num_events = numpy.empty(comm.Get_size(), dtype=numpy.uint64)
        else:
            num_events = None
        comm.Gather([lsum, MPI.UNSIGNED_LONG], [num_events, MPI.UNSIGNED_LONG], root=0) 
        
        # print 'Process',self.get_my_process_id(),':','Sent number:',lsum,'Collected numbers ->',num_events
        
        if self.get_my_process_id()==0:
            # Total number of spikes to be gathered
            gsum = num_events.sum()
            num_sent = tuple(num_events)
            aux_array = numpy.roll(num_events.cumsum(),1)
            aux_array[0] = 0
            offset = tuple(aux_array)
            
            gtime = numpy.empty(gsum, dtype=numpy.float64)
            gneuron_id = numpy.empty(gsum, dtype=numpy.uint64)
            gvalue = numpy.empty(gsum, dtype=numpy.float64)
        else:
            num_sent = None
            offset = None
            gtime = None
            gneuron_id = None
            gvalue = None
        
        # Gather the time and neuron_id arrays
        comm.Gatherv(time, [gtime, num_sent, offset, MPI.DOUBLE], root=0)
        comm.Gatherv(neuron_id, [gneuron_id, num_sent, offset, MPI.UNSIGNED_LONG], root=0)
        comm.Gatherv(value, [gvalue, num_sent, offset, MPI.DOUBLE], root=0)
        
        # print 'Process',self.get_my_process_id(),':','Collected time ->',gtime,'Collected Neurons ->',gneuron_id,'Collected values ->',gvalue
        if gtime is not None:
            gtime = gtime
            
        return (gtime,gneuron_id,gvalue)
        
    
    def get_synaptic_weights(self, **kwargs):
        '''
        Get the synaptic weights in the synapses, layer and during the time specified in the parameters.
        @param synaptic_layer Layer name as the section name in the config file (.cfg).
        @param source_indexes Indexes of the source cells of the synapses to get the activity.
        @param target_indexes Indexes of the target cells of the synapses to get the activity.
        @param init_time Initial time from which retrieve the activity.
        @param end_time Final time until which retrieve the activity.
        '''
        # Collect all the parameters
        if 'synaptic_layer' in kwargs:
            synaptic_layer_name = kwargs.pop('synaptic_layer')
        else:
            logger.error('Non specified synaptic layer in get_synaptic_weights function')
            raise Exception('Non-DefinedSynapticLayer')
        
        if synaptic_layer_name in self.layer_map:
            synaptic_layer = self.layer_map[synaptic_layer_name]
        else:
            logger.error('Invalid synaptic layer in get_synaptic_weights function')
            raise Exception('InvalidSynapticLayer')
        
        if 'source_indexes' in kwargs:
            source_indexes = kwargs.pop('source_indexes')
            
            if (max(source_indexes)>=(synaptic_layer.source_layer.number_of_neurons)):
                logger.error('Invalid source index in get_synaptic_weights function')
                raise Exception('InvalidSourceIndex')
        else:
            source_indexes = range(synaptic_layer.source_layer.number_of_neurons)
        
        if 'target_indexes' in kwargs:
            target_indexes = kwargs.pop('target_indexes')
            
            if (max(target_indexes)>=(synaptic_layer.target_layer.number_of_neurons)):
                logger.error('Invalid target index in get_synaptic_weights function')
                raise Exception('InvalidTargetIndex')
        else:
            target_indexes = range(synaptic_layer.target_layer.number_of_neurons)
            
        # print 'Process',self.get_my_process_id(),':','Source indexes:',source_indexes,'Target indexes ->',target_indexes
        
        if 'init_time' in kwargs:
            init_time = kwargs.pop('init_time')
        else:
            init_time = 0
            
        if 'end_time' in kwargs:
            end_time = kwargs.pop('end_time')
        else:
            end_time = float('inf')
        
        
        if not synaptic_layer.weight_recording:
            logger.error('Invalid synaptic layer in get_synaptic_weights function. The weights in this layer have not been recorded')
            raise Exception('InvalidSynapticLayer')
        
        # print 'Process',self.get_my_process_id(),':','Weight record:',synaptic_layer.weight_record
        
        # Calculate selected connection indexes
        connections = synaptic_layer.weight_record['connections']
        connection_indexes = [index for index in range(len(connections)) if (connections[index][0] in source_indexes) and (connections[index][1] in target_indexes)]
        selected_connections = connections[connection_indexes]
        
        # print 'Process',self.get_my_process_id(),':','Selected connections:',selected_connections
        
        # Calculate selected time indexes
        time = synaptic_layer.weight_record['time']
        time_indexes = (time>=init_time) & (time<=end_time) 
        selected_time = time[time_indexes]
        
        # Pick selected weights
        weights = synaptic_layer.weight_record['weights']
        selected_weights = numpy.array([record[connection_indexes] for record in weights[time_indexes]]).transpose()
        
        comm = MPI.COMM_WORLD
        
        # Send the number of elements
        time_elements = len(selected_time) 
        
#         # Gather every connection record individually
#         connection_aux = numpy.array(len(connection_indexes), dtype=numpy.uint64) 
#         connection_elements = numpy.array(0, dtype=numpy.uint64) 
#         comm.Reduce(connection_aux, connection_elements, op=MPI.SUM, root=0)
#         
#         logger.debug('Connection number sent: %s. Connection number collected: %s',connection_aux,connection_elements)
#         
        # Send the number of elements
        connection_aux = numpy.array([len(connection_indexes)], dtype=numpy.uint64)
        if self.get_my_process_id()==0:
            connection_elements = numpy.empty(comm.Get_size(), dtype=numpy.uint64)
        else:
            connection_elements = None
        comm.Gather([connection_aux, MPI.UNSIGNED_LONG], [connection_elements, MPI.UNSIGNED_LONG], root=0) 
        
        # logger.debug('Sent number: %s. Collected numbers: %s',connection_aux,connection_elements)
        
        
        if self.get_my_process_id()==0:
            # Total number of spikes to be gathered
            gsum = connection_elements.sum()
            con_num_sent = tuple(connection_elements*2)
            weight_num_sent = tuple(connection_elements*time_elements)
            aux_array = numpy.roll(connection_elements.cumsum(),1)
            aux_array[0] = 0
            con_offset = tuple(aux_array*2)
            weight_offset = tuple(aux_array*time_elements)
            
            gtime = selected_time
            gconnections = numpy.empty((gsum,2), dtype=numpy.uint64)
            gweights = numpy.empty((gsum, time_elements), dtype=numpy.float64)
        else:
            con_num_sent = None
            weight_num_sent = None
            con_offset = None
            weight_offset = None
            gtime = None
            gconnections = None
            gweights = None
        
        # Gather the time and neuron_id arrays
        #comm.Gatherv(time, [gtime, num_sent, offset, MPI.DOUBLE], root=0)
        comm.Gatherv(selected_connections, [gconnections, con_num_sent, con_offset, MPI.UNSIGNED_LONG], root=0)
        comm.Gatherv(selected_weights.ravel(), [gweights, weight_num_sent, weight_offset, MPI.DOUBLE], root=0)
        
        # print 'Process',self.get_my_process_id(),':','Collected time ->',gtime,'Collected Connections ->',gconnections,'Collected weights ->',gweights
        
        return (gtime,gconnections,gweights)
        
        
    def get_number_of_virtual_processes(self):
        '''
        Return the number of virtual processes. It might be used to decide the number of seeds to be generated.  
        '''
        return nest.GetKernelStatus(['total_num_virtual_procs'])[0]
    
    def get_my_process_id(self):
        '''
        Return the id-number of this process.   
        '''
        return nest.Rank()