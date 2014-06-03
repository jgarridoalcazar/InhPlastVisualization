'''
Created on May 12, 2014

@author: Jesus Garrido (jgarridoalcazar at gmail.com)
'''

import abc
import ConfigParser
import NeuronLayer
import InputLayer
import SynapticLayer
import time
import numpy
from Utils.Utils import ConfigSectionMap

class CerebellarModel(object):
    '''
    This class defines an abstract base class with all the methods
    needed in order to generate a cerebellum-like network. Subclasses
    are thought to implement the translation of the network into the
    target-simulator language (e.g., EDLUT, Nest, PyNN, ...)
    '''
    __metaclass__ = abc.ABCMeta
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates a new cerebellar model.
        @param config_file: Name of the file including all the network generation parameters.
        '''
        
        if ('config_file' in kwargs):
            self.config_file = kwargs.pop('config_file')                         
        else:
            print 'Non-specified cerebellum config file.'
            raise Exception('Non-DefinedCerebellumModel')
        
        # Initialize the cell layer map
        self.layer_map = {}
        
        
        super(CerebellarModel, self).__init__()
        
        return 
        
        
    def _create_network_elements(self):
        '''
        It creates all the cell layers of the cerebellar network.
        '''
        self.neuron_layers = []
        
        # Create cerebellar inputs (mossy fibers and inferior olive)
        mf_options = ConfigSectionMap(config_parser = self.config_parser, section = 'mflayer')
        self.mflayer = NeuronLayer.NeuronLayer(**mf_options)
        self.neuron_layers.append(self.mflayer)
        self.layer_map[self.mflayer.__name__] = self.mflayer
        
#         io_options = ConfigSectionMap(config_parser = self.config_parser, section = 'iolayer')
#         self.iolayer = InputLayer.InputLayer(**io_options)
#         self.neuron_layers.append(self.iolayer)
        
        # Create granule cell layer
        grc_options = ConfigSectionMap(config_parser = self.config_parser, section = 'grclayer')
        self.grclayer = NeuronLayer.NeuronLayer(**grc_options)
        self.neuron_layers.append(self.grclayer)
        self.layer_map[self.grclayer.__name__] = self.grclayer
        
        # Create Golgi cell layer
        goc_options = ConfigSectionMap(config_parser = self.config_parser, section = 'goclayer')
        self.goclayer = NeuronLayer.NeuronLayer(**goc_options)
        self.neuron_layers.append(self.goclayer)
        self.layer_map[self.goclayer.__name__] = self.goclayer
                
        # Create Purkinje cell layer
#         pc_options = ConfigSectionMap(config_parser = self.config_parser, section = 'pclayer')
#         self.pclayer = NeuronLayer.NeuronLayer(**pc_options)
#         self.neuron_layers.append(self.pclayer)
        
        # Create deep cerebellar nuclei cell layer
#         dcn_options = ConfigSectionMap(config_parser = self.config_parser, section = 'dcnlayer')
#         self.dcnlayer = NeuronLayer.NeuronLayer(**dcn_options)
#         self.neuron_layers.append(self.dcnlayer)
#         
        return       
        
    def _create_synapses(self):
        '''
        It creates all the synaptic connections of the cerebellar network.
        '''
        
        self.synaptic_layers = []
        # Get the process id to know which seed has to be used
        process_id = self.get_my_process_id()
        
        # Create MF-GrC synaptic layer
        mfgrc_options = ConfigSectionMap(config_parser = self.config_parser, section = 'mfgrcsynapsis')
        mfgrc_options['source_layer'] = self.mflayer
        mfgrc_options['target_layer'] = self.grclayer
        mfgrc_options['random_generator'] = self.simulation_options['pyrngs'][process_id]
        self.mfgrclayer = SynapticLayer.SynapticLayer(**mfgrc_options)
        self.synaptic_layers.append(self.mfgrclayer)
        self.layer_map[self.mfgrclayer.__name__] = self.mfgrclayer
        
        # Create MF-GoC synaptic layer
        mfgoc_options = ConfigSectionMap(config_parser = self.config_parser, section = 'mfgocsynapsis')
        mfgoc_options['source_layer'] = self.mflayer
        mfgoc_options['target_layer'] = self.goclayer
        mfgoc_options['random_generator'] = self.simulation_options['pyrngs'][process_id]
        self.mfgoclayer = SynapticLayer.SynapticLayer(**mfgoc_options)
        self.synaptic_layers.append(self.mfgoclayer)
        self.layer_map[self.mfgoclayer.__name__] = self.mfgoclayer
        
        # Create GrC-GoC synaptic layer
        grcgoc_options = ConfigSectionMap(config_parser = self.config_parser, section = 'grcgocsynapsis')
        grcgoc_options['source_layer'] = self.grclayer
        grcgoc_options['target_layer'] = self.goclayer
        grcgoc_options['random_generator'] = self.simulation_options['pyrngs'][process_id]
        self.grcgoclayer = SynapticLayer.SynapticLayer(**grcgoc_options)
        self.synaptic_layers.append(self.grcgoclayer)
        self.layer_map[self.grcgoclayer.__name__] = self.grcgoclayer
        
        # Create GoC-GoC synaptic layer
        gocgrc_options = ConfigSectionMap(config_parser = self.config_parser, section = 'gocgrcsynapsis')
        gocgrc_options['source_layer'] = self.goclayer
        gocgrc_options['target_layer'] = self.grclayer
        gocgrc_options['random_generator'] = self.simulation_options['pyrngs'][process_id]
        self.gocgrclayer = SynapticLayer.SynapticLayer(**gocgrc_options)
        self.synaptic_layers.append(self.gocgrclayer)
        self.layer_map[self.gocgrclayer.__name__] = self.gocgrclayer
        
        # Create GoC-GoC synaptic layer
        gocgoc_options = ConfigSectionMap(config_parser = self.config_parser, section = 'gocgocsynapsis')
        gocgoc_options['source_layer'] = self.goclayer
        gocgoc_options['target_layer'] = self.goclayer
        gocgoc_options['random_generator'] = self.simulation_options['pyrngs'][process_id]
        self.gocgoclayer = SynapticLayer.SynapticLayer(**gocgoc_options)
        self.synaptic_layers.append(self.gocgoclayer)
        self.layer_map[self.gocgoclayer.__name__] = self.gocgoclayer
        
        # Create GrC-PC synaptic layer
#         grcpc_options = ConfigSectionMap(config_parser = self.config_parser, section = 'grcpcsynapsis')
#         grcpc_options['source_layer'] = self.grclayer
#         grcpc_options['target_layer'] = self.pclayer
#         self.grcpclayer = SynapticLayer.SynapticLayer(**grcpc_options)
#         self.synaptic_layers.append(self.grcpclayer)
        
        # Create IO-PC synaptic layer
#         iopc_options = ConfigSectionMap(config_parser = self.config_parser, section = 'iopcsynapsis')
#         iopc_options['source_layer'] = self.iolayer
#         iopc_options['target_layer'] = self.pclayer
#         self.iopclayer = SynapticLayer.SynapticLayer(**iopc_options)
#         self.synaptic_layers.append(self.iopclayer)
        
        # Create MF-DCN synaptic layer
#         mfdcn_options = ConfigSectionMap(config_parser = self.config_parser, section = 'mfdcnsynapsis')
#         mfdcn_options['source_layer'] = self.mflayer
#         mfdcn_options['target_layer'] = self.dcnlayer
#         self.mfdcnlayer = SynapticLayer.SynapticLayer(**mfdcn_options)
#         self.synaptic_layers.append(self.mfdcnlayer)
        
        # Create PC-DCN synaptic layer
#         pcdcn_options = ConfigSectionMap(config_parser = self.config_parser, section = 'pcdcnsynapsis')
#         pcdcn_options['source_layer'] = self.pclayer
#         pcdcn_options['target_layer'] = self.dcnlayer
#         self.pcdcnlayer = SynapticLayer.SynapticLayer(**pcdcn_options)
#         self.synaptic_layers.append(self.pcdcnlayer)
        
        return 

    def initialize_simulation(self):
        '''
        Reset the simulator, build the network and set additional parameters (e.g. the simulation time-step, the number of threads, ...).     
        '''
        
        print 'Parsing configuration file ', self.config_file
    
        self.config_parser = ConfigParser.ConfigParser()
        self.config_parser.read(self.config_file)
        
        # Check simulation options
        self.simulation_options = ConfigSectionMap(config_parser = self.config_parser, section = 'simulation')
        if not 'seed' in self.simulation_options:
            self.simulation_options['seed']  = time()
        
        # Generate random number generators for each python virtual process
        self.simulation_options['pyrngs'] = [numpy.random.RandomState(s) for s in range(self.simulation_options['seed'], self.simulation_options['seed']+self.get_number_of_virtual_processes())]
        
        # Check weight recording time step
        if not 'weight_recording_step' in self.simulation_options:
            self.simulation_options['weight_recording_step'] = float("inf")
    
        
        return
    
    def build_network(self):
        '''
        Generate a model based on the inherited cerebellar model.
        '''
        
        # Create the network elements
        self._create_network_elements();
        
        # Connect those elements each others
        self._create_synapses();
        
        return
    
    def get_number_of_elements(self, **kwargs):
        '''
        Retrieve the number of elements (neurons or synapses in a layer.
        @param layer The name of the layer to get the number of elements.
        '''
        
        if 'layer' in kwargs:
            layer_name = kwargs.pop('layer')
        else:
            print 'Non-specified layer name.'
            raise Exception('Non-DefinedLayerName')
        
        if layer_name in self.layer_map:
            layer = self.layer_map[layer_name]
            
            if (isinstance(layer, InputLayer.InputLayer)):
                return layer.number_of_neurons
            elif (isinstance(layer, SynapticLayer.SynapticLayer)):
                return layer.number_of_synapses
            else:
                print 'Invalid layer type in',layer.__name__
                raise Exception('InvalidLayerType')
        else:
            print layer_name,'does not exist in the model'
            raise Exception('InvalidLayerName')
    
    @abc.abstractmethod
    def add_ac_current(self, **kwargs):
        '''
        Set a new ac current that will be conveyed to the cerebellar mossy fibers according to the parameters.
        '''
        return
    
    @abc.abstractmethod
    def set_dc_current(self, **kwargs):
        '''
        Set the mossy fiber dc currents to the values specifed in the parameter.
        '''
        return
        
    @abc.abstractmethod
    def simulate_network(self, time):
        '''
        Simulate the network until the specified time (in seconds).
        '''
        return
    
    @abc.abstractmethod
    def get_spike_activity(self, **kwargs):
        '''
        Get the spikes been fired in the cells, layer and time specified in the parameters.  
        '''
        return
    
    @abc.abstractmethod
    def get_state_variable(self, **kwargs):
        '''
        Get the state variables in the cells, layer and during the time specified in the parameters.  
        '''
        return
    
    @abc.abstractmethod
    def get_synaptic_weights(self, **kwargs):
        '''
        Get the synaptic weights in the synapses, layer and during the time specified in the parameters.  
        '''
        
        return
    
    @abc.abstractmethod
    def get_number_of_virtual_processes(self):
        '''
        Return the number of virtual processes. It might be used to decide the number of   
        '''
        return

    @abc.abstractmethod
    def get_my_process_id(self):
        '''
        Return the id-number of this process.   
        '''
        return 

        
        
        
