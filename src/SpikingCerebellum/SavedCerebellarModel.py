'''
Created on May 12, 2014

@author: Jesus Garrido (jgarridoalcazar at gmail.com)
'''

import abc
from CerebellarModel import CerebellarModel,ConfigSectionMap
import InputLayer
import NeuronLayer
import SynapticLayer
import numpy
import math
import os
import glob
import ConfigParser
from mpi4py import MPI

class SavedCerebellarModel(CerebellarModel):
    '''
    This class defines an inherited class including all the methods
    needed in order to generate a cerebellum-like network with the activity
    data loaded from data-files.
    '''
    __metaclass__ = abc.ABCMeta
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates a new cerebellar model.
        @param config_dict: Dictionary including all the network generation parameters.
        @param simulation_folder: Folder where the simulation results have been stored. [Optional] If not specified, it will be calculated as in NestCerebellarModel.
        '''
        
        if 'simulation_folder' in kwargs:
            self.simulation_folder = kwargs.pop('simulation_folder')
        else:
            self.simulation_folder = None
        
        super(SavedCerebellarModel, self).__init__(**kwargs)
        
    def initialize_simulation(self):
        '''
        Reset the NEST simulator, build the network and set additional parameters (e.g. the simulation time-step, the number of threads, ...).     
        '''
        
        super(SavedCerebellarModel, self).initialize_simulation()
        
        if self.simulation_folder is None:
            if 'data_path' in self.simulation_options:
                self.simulation_folder = self.simulation_options['data_path']
            else:
                self.simulation_folder = './results'
            
            if 'simulation_name' in self.simulation_options:
                self.simulation_folder = self.simulation_folder + '/' + self.simulation_options['simulation_name']
            
        if not os.path.exists(self.simulation_folder):
            print 'Error: Simulation folder', self.simulation_folder, 'does not exist.'
            raise Exception('InvalidSimulationFolder')
        
        
        # Initialize ac_current and dc_current stimulation
        self.ac_generator = []
        self.dc_current = []
        
        # Load the data from the files and build the network
        self._build_network()
                
        self.simulation_time = 0
        
        return
    
    def _file_name_iterator(self, layer, type=None):
        '''
        Return an iterator to the file names of recorded activity/state/weights in the specified layer.
        @param layer Layer to get the files.
        @param type  Type of the files to get in 'spike','state','weight'
        '''
        name_pattern = self.simulation_folder + '/' + layer.__name__
        
        if isinstance(layer,InputLayer.InputLayer):
            if type=='spike':
                name_pattern = name_pattern + '_spike*.gdf'
            elif type=='state':
                name_pattern = name_pattern + '_mult*.dat'
            else:
                print 'Warning: Invalid type of register', type, ' used.'
                raise Exception('InvalidType')
        elif isinstance(layer,SynapticLayer.SynapticLayer):
            if type=='weight' or type is None:
                name_pattern = name_pattern + '_weights*.dat'
            else:
                print 'Warning: Invalid type of register', type, ' used.'
                raise Exception('InvalidType')
        else:
            print 'Warning: Invalid layer type used.'
            raise Exception('InvalidLayerType')
        
        file_list = glob.glob(name_pattern)
        
        for name in file_list:
            yield name
                
    
    
    def _initialize_weight_recording_buffer(self):
        '''
        Initialize and load the weight recording buffer from the files.
        '''
        
        # Check if recording_time_step is above 0
        if self.simulation_options['weight_recording_step'] >= float("inf"):
            return
        
        if not self.simulation_options['record_to_file']:
            return
            
        for layer in self.synaptic_layers:
            if layer.weight_recording:
                layer.weight_record = dict()
            
                # print 'Process',self.get_my_process_id(),':','Source layer:',layer.source_layer.nest_layer,'Target layer:', layer.target_layer.nest_layer
                layer.weight_record['connections'] = None
                layer.weight_record['weights'] = None
                layer.weight_record['time'] = None
                
                
                # Get the weight file list of the layer
                file_iterator = self._file_name_iterator(layer=layer, type='weight')
                for file_name in file_iterator:
                    
                    print 'Process',self.get_my_process_id(),':','Reading weight register file', file_name
                    
                    with open(file_name, 'r') as f:
                        # Get the connections array and transform into an array of tuples
                        connections_array = numpy.array(f.readline().split(' ')[:-1], int)
                        connections_matrix = numpy.reshape(connections_array,(-1,2))
                        
                        if layer.weight_record['connections'] is None:
                            layer.weight_record['connections'] = connections_matrix
                        else:
                            layer.weight_record['connections'] = numpy.append(layer.weight_record['connections'],connections_matrix,axis=0)
                        
                        # Check the number of connections in the file
                        if len(connections_matrix):
                            # Get the rest of the weight file
                            data = numpy.genfromtxt(f, delimiter=' ')
                        
                            # Extract the time and the wiehgts
                            if layer.weight_record['time'] is None:
                                layer.weight_record['time'] = data[:,0]
                        
                            if layer.weight_record['weights'] is None:
                                layer.weight_record['weights'] = data[:,1:]
                            else:
                                layer.weight_record['weights'] = numpy.append(layer.weight_record['weights'],data[:,1:],axis=0)
            else:
                layer.weight_record = None

        return
    
    def _initialize_layer_buffer(self):
        '''
        Initialize and load the spike activity and the state of every layer.
        '''
        
        # Check if record_to_file has been enabled
        if not self.simulation_options['record_to_file']:
            return
            
        for layer in self.neuron_layers:
            # Load registered spike activity
            if layer.register_activity:
                layer.activity_record = dict()
            
                # print 'Process',self.get_my_process_id(),':','Source layer:',layer.source_layer.nest_layer,'Target layer:', layer.target_layer.nest_layer
                layer.activity_record['cell'] = None
                layer.activity_record['time'] = None
                
                # Get the layer min index
                layer.MinIndex = self.config_dict[layer.__name__]['minindex']
                
                # Get the weight file list of the layer
                file_iterator = self._file_name_iterator(layer=layer, type='spike')
                for file_name in file_iterator:
                    
                    print 'Process',self.get_my_process_id(),':','Reading activity register file', file_name
                    
                    # Get the activity file
                    data = numpy.loadtxt(file_name)
                    
                    if len(data):    
                        # Get the connections array and transform into an array of tuples
                        cell_array = data[:,0] - layer.MinIndex
                        time_array = data[:,1]/1.e3
                            
                        if layer.activity_record['cell'] is None:
                            layer.activity_record['cell'] = cell_array
                        else:
                            layer.activity_record['cell'] = numpy.append(layer.activity_record['cell'],cell_array)
                        
                        if layer.activity_record['time'] is None:
                            layer.activity_record['time'] = cell_array
                        else:
                            layer.activity_record['time'] = numpy.append(layer.activity_record['time'],time_array)
            else:
                layer.activity_record = None
                
            
            # Load registered state variables
            if layer.record_vars:
                # If only an string has been used, embed it in an array
                if isinstance(layer.record_vars,str):
                    layer.record_vars = [layer.record_vars]
                
                # Attention: This code assumes that every variable in the list have been recorded in the file. Otherwise, it might crash.
                
                layer.state_record = dict()
            
                # print 'Process',self.get_my_process_id(),':','Source layer:',layer.source_layer.nest_layer,'Target layer:', layer.target_layer.nest_layer
                layer.state_record['cell'] = None
                layer.state_record['time'] = None
                
                for var in layer.record_vars:
                    layer.state_record[var] = None
                
                # Get the layer min index
                layer.MinIndex = self.config_dict[layer.__name__]['minindex']
                
                # Get the weight file list of the layer
                file_iterator = self._file_name_iterator(layer=layer, type='state')
                for file_name in file_iterator:
                    
                    print 'Process',self.get_my_process_id(),':','Reading state recorded file', file_name
                    
                    # Get the activity file
                    data = numpy.loadtxt(file_name)
                        
                    if len(data):
                        # Get the connections array and transform into an array of tuples
                        cell_array = data[:,0] - layer.MinIndex
                        time_array = data[:,1]/1.e3
                            
                        if layer.state_record['cell'] is None:
                            layer.state_record['cell'] = cell_array
                        else:
                            layer.state_record['cell'] = numpy.append(layer.state_record['cell'],cell_array)
                        
                        if layer.state_record['time'] is None:
                            layer.state_record['time'] = time_array
                        else:
                            layer.state_record['time'] = numpy.append(layer.state_record['time'],time_array)
                            
                        for index in range(len(layer.record_vars)):
                            var = layer.record_vars[index]
                            array = data[:,index+2]
                            if layer.state_record[var] is None:
                                layer.state_record[var] = array
                            else:
                                layer.state_record[var] = numpy.append(layer.state_record[var],array)
                     
            else:
                layer.state_record = None

        return
    
    def _build_network(self):
        '''
        Generate a NEST model based on the inherited cerebellar model.
        '''
        
        super(SavedCerebellarModel, self)._build_network()
        
        # Initialize weight recording buffer
        self._initialize_weight_recording_buffer()
        
        # Initialize spike activity and state variables
        self._initialize_layer_buffer()
        
        return
    
        
    def add_ac_current(self, **kwargs):
        '''
        Set a new ac current that will be conveyed to the cerebellar mossy fibers according to the parameters.
        @param amplitude Amplitude of the wave (in A)
        @param offset Constant amplitude offset (in pA)
        @param frequency Frequency of the ac_generator (in Hz)
        @param phase Phase of the sine current (0-360 deg)
        '''
        return
    
    def set_dc_current(self, **kwargs):
        '''
        Set the mossy fiber dc currents to the values specifed in the parameter.
        @param amplitude Amplitude of the current in A. A list with the same length of the number of mossy fibers must be used.
        '''
        return
        
    def simulate_network(self,time):
        '''
        Simulate the network for the specified time (in seconds).
        @param time Length to be simulated (in seconds)
        '''
        
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
            print 'Non specified neuron layer in get_spike_activity function.'
            raise Exception('Non-DefinedNeuronLayer')
        
        
        if neuron_layer_name in self.layer_map:
            neuron_layer = self.layer_map[neuron_layer_name]
        else:
            print 'Invalid neuron layer in get_spike_activity function.'
            raise Exception('InvalidNeuronLayer')
        
        if 'neuron_indexes' in kwargs:
            neuron_indexes = kwargs.pop('neuron_indexes')
            
            if (max(neuron_indexes)>=(neuron_layer.number_of_neurons)):
                print 'Invalid neuron index in get_state_variable function.'
                raise Exception('InvalidNeuronIndex')
        else:
            neuron_indexes = range(neuron_layer.number_of_neurons)
            
        if 'init_time' in kwargs:
            init_time = kwargs.pop('init_time')
        else:
            init_time = 0
            
        if 'end_time' in kwargs:
            end_time = kwargs.pop('end_time')
        else:
            end_time = float('inf')
        
        # Check whether the neuron has been recorded
        if not neuron_layer.register_activity: 
            print 'Invalid neuron layer in get_spike_activity function. The activity in this layer has not been recorded.'
            raise Exception('InvalidNeuronLayer')
        
        # If the spike detector is local to this process
        time = neuron_layer.activity_record['time']
        neuron_id = neuron_layer.activity_record['cell']
            
        # Select only those events in the specified interval
        index = (time>=init_time) & (time<=end_time) 
        time = time[index]
        neuron_id = neuron_id[index]
            
        # Select only those events in the specified cells
        index = numpy.in1d(neuron_id, neuron_indexes)
        time = time[index]
        neuron_id = neuron_id[index]
        
        # print 'Process',self.get_my_process_id(),':','Collected time ->',gtime,'Collected Neurons ->',gneuron_id
        return (time,neuron_id)
    
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
            print 'Non specified neuron layer in get_state_variable function.'
            raise Exception('Non-DefinedNeuronLayer')
        
        
        if neuron_layer_name in self.layer_map:
            neuron_layer = self.layer_map[neuron_layer_name]
        else:
            print 'Invalid neuron layer in get_state_variable function.'
            raise Exception('InvalidNeuronLayer')
        
        if 'neuron_indexes' in kwargs:
            neuron_indexes = kwargs.pop('neuron_indexes')
            
            if (max(neuron_indexes)>=(neuron_layer.number_of_neurons)):
                print 'Invalid neuron index in get_state_variable function.'
                raise Exception('InvalidNeuronIndex')
        else:
            neuron_indexes = range(neuron_layer.number_of_neurons)
            
        if 'state_variable' in kwargs:
            variable_name = kwargs.pop('state_variable')
        else:
            print 'Non specified state_variable in get_state_variable function.'
            raise Exception('Non-DefinedStateVariable')
        
        if not variable_name in neuron_layer.record_vars:
            print 'Invalid state variable in get_state_variable function.',variable_name,'has not been recorded.'
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
        # print 'Process',self.get_my_process_id(),': recording_events=',recording_events
        neuron_id = neuron_layer.state_record['cell']
        time = neuron_layer.state_record['time']
        value = neuron_layer.state_record[variable_name]
        
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
        
        return (time,neuron_id,value)
        
    
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
            print 'Non specified synaptic layer in get_synaptic_weights function.'
            raise Exception('Non-DefinedSynapticLayer')
        
        if synaptic_layer_name in self.layer_map:
            synaptic_layer = self.layer_map[synaptic_layer_name]
        else:
            print 'Invalid synaptic layer in get_synaptic_weights function.'
            raise Exception('InvalidSynapticLayer')
        
        if 'source_indexes' in kwargs:
            source_indexes = kwargs.pop('source_indexes')
            
            if (max(source_indexes)>=(synaptic_layer.source_layer.number_of_neurons)):
                print 'Invalid source index in get_synaptic_weights function.'
                raise Exception('InvalidSourceIndex')
        else:
            source_indexes = range(synaptic_layer.source_layer.number_of_neurons)
        
        if 'target_indexes' in kwargs:
            target_indexes = kwargs.pop('target_indexes')
            
            if (max(target_indexes)>=(synaptic_layer.target_layer.number_of_neurons)):
                print 'Invalid target index in get_synaptic_weights function.'
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
            print 'Invalid synaptic layer in get_synaptic_weights function. The weights in this layer has not been recorded.'
            raise Exception('InvalidSynapticLayer')
        
        # print 'Process',self.get_my_process_id(),':','Weight record:',synaptic_layer.weight_record
        
        # Calculate selected connection indexes
        connections = synaptic_layer.weight_record['connections']
        connection_indexes = [index for index in range(len(connections)) if (connections[index][0] in source_indexes) and (connections[index][1] in target_indexes)]
        gconnections = connections[connection_indexes]
        
        # print 'Process',self.get_my_process_id(),':','Selected connections:',selected_connections
        
        # Calculate selected time indexes
        time = synaptic_layer.weight_record['time']
        time_indexes = (time>=init_time) & (time<=end_time) 
        gtime = time[time_indexes]
        
        # Pick selected weights
        weights = synaptic_layer.weight_record['weights']
        gweights = numpy.array([record[connection_indexes] for record in weights[time_indexes]]).transpose()
        
        return (gtime,gconnections,gweights)
        
        
    def get_number_of_virtual_processes(self):
        '''
        Return the number of virtual processes. It might be used to decide the number of seeds to be generated.  
        '''
        return self.config_dict['nest']['number_of_virtual_processes']
    
    def get_my_process_id(self):
        '''
        Return the id-number of this process.   
        '''
        comm = MPI.COMM_WORLD
        return comm.Get_rank()