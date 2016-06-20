'''
Created on May 12, 2014

@author: Jesus Garrido (jgarridoalcazar at gmail.com)
'''

import abc
import InputLayer
import time
import numpy
import logging
 
 
logger = logging.getLogger('Simulation')

class CerebellarModel(object):
    '''
    This class defines an abstract base class with all the methods
    needed in order to generate a cerebellum-like network. Subclasses
    are thought to implement the translation of the network into the
    target-simulator language (e.g., EDLUT, Nest, PyNN, ...)
    '''
#    __metaclass__ = abc.ABCMeta
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates a new cerebellar model.
        @param config_dict: Dictionary including all the network generation parameters.
        '''
        
        if ('config_dict' in kwargs):
            self.config_dict = kwargs.pop('config_dict')                 
        else:
            logger.error('Non-specified cerebellum config file')
            raise Exception('Non-DefinedCerebellumModel')
        
        if ('network' not in self.config_dict):
            self.config_dict['network'] = {}
            
        if ('save_file' not in self.config_dict['network']):
            self.save_file = None
        else:
            self.save_file = self.config_dict['network']['save_file']
            
        if ('load_file' not in self.config_dict['network']):
            self.load_file = None
        else:
            self.load_file = self.config_dict['network']['load_file']
            
        if ('length' not in self.config_dict['network']):
            self.config_dict['network']['length'] = None
            self.network_size = None
        else:
            # Create the dimensions of the cube from the lenght of the size
            length = self.config_dict['network']['length']
            self.network_size = (length, length, length)
            
        
        # Initialize the cell layer map
        self.layer_map = {}
        
        super(CerebellarModel, self).__init__()
        
        return 
        
        
    def _create_neurons(self):
        '''
        It creates all the cell layers of the cerebellar network.
        '''
        self.neuron_layers = []
        
        if (self.config_dict['simulation']['use_mpi']):
            import NeuronLayerMPI as NeuronLayer
        else:
            import NeuronLayerNoMPI as NeuronLayer
        
        # Check if the network has to be loaded from file
        if (self.load_file is not None):
            import h5py
        
            file = h5py.File(self.load_file)
        else:
            file = None
        
        # Create cerebellar inputs (mossy fibers and inferior olive)
        mf_options = self.config_dict['mflayer']
        mf_options['random_generator'] = self.get_local_py_rng()
        if (self.network_size is not None):
            mf_options['size'] = self.network_size
        node = _search_hdf5_group(file, 'mflayer')
        if node is not None:
            mf_options['load_from_file'] = node
        self.mflayer = NeuronLayer.NeuronLayer(**mf_options)
        self.neuron_layers.append(self.mflayer)
        self.layer_map[self.mflayer.__name__] = self.mflayer
        
#         io_options = self.config_dict['iolayer']
#         self.iolayer = InputLayer.InputLayer(**io_options)
#         self.neuron_layers.append(self.iolayer)
        
        # Create granule cell layer
        grc_options = self.config_dict['grclayer']
        grc_options['random_generator'] = self.get_local_py_rng()
        if (self.network_size is not None):
            grc_options['size'] = self.network_size
        node = _search_hdf5_group(file, 'grclayer')
        if node is not None:
            grc_options['load_from_file'] = node
        self.grclayer = NeuronLayer.NeuronLayer(**grc_options)
        self.neuron_layers.append(self.grclayer)
        self.layer_map[self.grclayer.__name__] = self.grclayer
        
        # Create Golgi cell layer
        goc_options = self.config_dict['goclayer']
        goc_options['random_generator'] = self.get_local_py_rng()
        if (self.network_size is not None):
            goc_options['size'] = self.network_size
        node = _search_hdf5_group(file, 'goclayer')
        if node is not None:
            goc_options['load_from_file'] = node
        self.goclayer = NeuronLayer.NeuronLayer(**goc_options)
        self.neuron_layers.append(self.goclayer)
        self.layer_map[self.goclayer.__name__] = self.goclayer
                
        # Create Purkinje cell layer (not compulsory yet)
        if 'pclayer' in self.config_dict:
            pc_options = self.config_dict['pclayer']
            pc_options['random_generator'] = self.get_local_py_rng()
            if (self.network_size is not None):
                pc_options['size'] = self.network_size
            node = _search_hdf5_group(file, 'pclayer')
            if node is not None:
                pc_options['load_from_file'] = node
            self.pclayer = NeuronLayer.NeuronLayer(**pc_options)
            self.neuron_layers.append(self.pclayer)
            self.layer_map[self.pclayer.__name__] = self.pclayer
        
        # Create deep cerebellar nuclei cell layer
#         dcn_options = self.config_dict['dcnlayer']
#         self.dcnlayer = NeuronLayer.NeuronLayer(**dcn_options)
#         self.neuron_layers.append(self.dcnlayer)
#        
        if file is not None:
            file.close()
        
        return       
        
    def _create_synapses(self):
        '''
        It creates all the synaptic connections of the cerebellar network.
        '''
        
        self.synaptic_layers = []
        
        if (self.config_dict['simulation']['use_mpi']):
            import SynapticLayerMPI as SynapticLayer
        else:
            import SynapticLayerNoMPI as SynapticLayer
        
        # Get the process_id to select the correct random number generator    
        process_id = self.get_my_process_id()
        
        # Check if the network has to be loaded from file
        if (self.load_file is not None):
            import h5py
        
            file = h5py.File(self.load_file)
        else:
            file = None
        
        # Create MF-GrC synaptic layer
        mfgrc_options = self.config_dict['mfgrcsynapsis']
        mfgrc_options['source_layer'] = self.mflayer
        mfgrc_options['target_layer'] = self.grclayer
        mfgrc_options['random_generator'] = self.get_local_py_rng()
        node = _search_hdf5_group(file, 'mfgrcsynapsis')
        if node is not None:
            mfgrc_options['load_from_file'] = node
        self.mfgrclayer = SynapticLayer.SynapticLayer(**mfgrc_options)
        self.synaptic_layers.append(self.mfgrclayer)
        self.layer_map[self.mfgrclayer.__name__] = self.mfgrclayer
        
        # Create MF-GoC synaptic layer
        mfgoc_options = self.config_dict['mfgocsynapsis']
        mfgoc_options['source_layer'] = self.mflayer
        mfgoc_options['target_layer'] = self.goclayer
        mfgoc_options['random_generator'] = self.get_local_py_rng()
        node = _search_hdf5_group(file, 'mfgocsynapsis')
        if node is not None:
            mfgoc_options['load_from_file'] = node
        self.mfgoclayer = SynapticLayer.SynapticLayer(**mfgoc_options)
        self.synaptic_layers.append(self.mfgoclayer)
        self.layer_map[self.mfgoclayer.__name__] = self.mfgoclayer
        
        # Create GrC-GoC synaptic layer
        grcgoc_options = self.config_dict['grcgocsynapsis']
        grcgoc_options['source_layer'] = self.grclayer
        grcgoc_options['target_layer'] = self.goclayer
        grcgoc_options['random_generator'] = self.get_local_py_rng()
        node = _search_hdf5_group(file, 'grcgocsynapsis')
        if node is not None:
            grcgoc_options['load_from_file'] = node
        self.grcgoclayer = SynapticLayer.SynapticLayer(**grcgoc_options)
        self.synaptic_layers.append(self.grcgoclayer)
        self.layer_map[self.grcgoclayer.__name__] = self.grcgoclayer
        
        # Create GoC-GrC synaptic layer
        gocgrc_options = self.config_dict['gocgrcsynapsis']
        gocgrc_options['source_layer'] = self.goclayer
        gocgrc_options['target_layer'] = self.grclayer
        gocgrc_options['intermediate_layer'] = self.mflayer
        gocgrc_options['intermediate_to_target_synaptic_layer'] = self.mfgrclayer
        gocgrc_options['random_generator'] = self.get_local_py_rng()
        node = _search_hdf5_group(file, 'gocgrcsynapsis')
        if node is not None:
            gocgrc_options['load_from_file'] = node
        self.gocgrclayer = SynapticLayer.SynapticLayer(**gocgrc_options)
        self.synaptic_layers.append(self.gocgrclayer)
        self.layer_map[self.gocgrclayer.__name__] = self.gocgrclayer
        
        # Create GoC-GoC synaptic layer
        gocgoc_options = self.config_dict['gocgocsynapsis']
        gocgoc_options['source_layer'] = self.goclayer
        gocgoc_options['target_layer'] = self.goclayer
        gocgoc_options['random_generator'] = self.get_local_py_rng()
        node = _search_hdf5_group(file, 'gocgocsynapsis')
        if node is not None:
            gocgoc_options['load_from_file'] = node
        self.gocgoclayer = SynapticLayer.SynapticLayer(**gocgoc_options)
        self.synaptic_layers.append(self.gocgoclayer)
        self.layer_map[self.gocgoclayer.__name__] = self.gocgoclayer
        
        # Create GrC-PC synaptic layer (optional)
        if 'grcpcsynapsis' in self.config_dict:
            grcpc_options = self.config_dict['grcpcsynapsis']
            grcpc_options['source_layer'] = self.grclayer
            grcpc_options['target_layer'] = self.pclayer
            grcpc_options['random_generator'] = self.get_local_py_rng()
            node = _search_hdf5_group(file, 'grcpcsynapsis')
            if node is not None:
                grcpc_options['load_from_file'] = node
            self.grcpclayer = SynapticLayer.SynapticLayer(**grcpc_options)
            self.synaptic_layers.append(self.grcpclayer)
            self.layer_map[self.grcpclayer.__name__] = self.grcpclayer
            
        # Create PC-PC inhibitory connections (optional)
        if 'pcpcsynapsis' in self.config_dict:
            pcpc_options = self.config_dict['pcpcsynapsis']
            pcpc_options['source_layer'] = self.pclayer
            pcpc_options['target_layer'] = self.pclayer
            pcpc_options['random_generator'] = self.get_local_py_rng()
            node = _search_hdf5_group(file, 'pcpcsynapsis')
            if node is not None:
                pcpc_options['load_from_file'] = node
            self.pcpclayer = SynapticLayer.SynapticLayer(**pcpc_options)
            self.synaptic_layers.append(self.pcpclayer)
            self.layer_map[self.pcpclayer.__name__] = self.pcpclayer
        
        # Create IO-PC synaptic layer
#         iopc_options = self.config_dict['iopcsynapsis']
#         iopc_options['source_layer'] = self.iolayer
#         iopc_options['target_layer'] = self.pclayer
#         self.iopclayer = SynapticLayer.SynapticLayer(**iopc_options)
#         self.synaptic_layers.append(self.iopclayer)
        
        # Create MF-DCN synaptic layer
#         mfdcn_options = self.config_dict['mfdcnsynapsis']
#         mfdcn_options['source_layer'] = self.mflayer
#         mfdcn_options['target_layer'] = self.dcnlayer
#         self.mfdcnlayer = SynapticLayer.SynapticLayer(**mfdcn_options)
#         self.synaptic_layers.append(self.mfdcnlayer)
        
        # Create PC-DCN synaptic layer
#         pcdcn_options = self.config_dict['pcdcnsynapsis']
#         pcdcn_options['source_layer'] = self.pclayer
#         pcdcn_options['target_layer'] = self.dcnlayer
#         self.pcdcnlayer = SynapticLayer.SynapticLayer(**pcdcn_options)
#         self.synaptic_layers.append(self.pcdcnlayer)

        if file is not None:
            file.close()        
        
        return 

    def initialize_simulation(self):
        '''
        Reset the simulator, build the network and set additional parameters (e.g. the simulation time-step, the number of threads, ...).     
        '''
        
        # Check simulation options
        self.simulation_options = self.config_dict['simulation']
        if not 'seed' in self.simulation_options:
            self.simulation_options['seed']  = int(time.time())
        
        
        # Check weight recording time step
        if not 'weight_recording_step' in self.simulation_options:
            self.simulation_options['weight_recording_step'] = float("inf")
    
        # Check activity recording time step
        if not 'activity_recording_step' in self.simulation_options:
            self.simulation_options['activity_recording_step'] = float("inf")
    
    
        logger.debug('Cerebellar simulation initialized')
        
        return
    
    def save_network(self):
        '''
        Save the network to the specified hdf5 file.
        @param filename: Name of the file to store the network.
        '''
        
        if self.save_file is None:
            return
        
        filename = self.save_file
        
        logger.info('Saving network to hdf5 file %s', filename)
        
        # Update the neuron states before saving it to file
        self.update_neuron_states()
       
       
        if self.get_my_process_id()==0:
            import h5py
            
            with h5py.File(filename, 'w') as file:
                
                # Saving neuron layers
                for ind, layer in enumerate(self.neuron_layers):
                    # Show info
                    logger.debug('Writing neuron layer %s',layer.__name__)
                    # Create a subgroup for this neuron layer
                    subgroup = file.create_group('neuron_layer_'+str(ind))
                    # Save the neuron layer in it
                    layer.save_layer(subgroup)
                    
                # Saving synaptic layers
                for ind, layer in enumerate(self.synaptic_layers):
                    # Show info
                    logger.debug('Writing neuron layer %s',layer.__name__)
                    # Create a subgroup for this synaptic layer
                    subgroup = file.create_group('synaptic_layer_'+str(ind))
                    # Save the synaptic layer in it
                    layer.save_layer(subgroup)
                
                file.flush()
        else:
            # This is just in case other processes have to communicate something (neurons or connections)
            for ind, layer in enumerate(self.neuron_layers):
                # Save the neuron layer in it
                layer.save_layer(None)
                    
            for ind, layer in enumerate(self.synaptic_layers):
                # Save the synaptic layer in it
                layer.save_layer(None)
        
        logger.debug('File writing ended')
            
    def visualize_network(self):
        '''
        Visualize the network structure by using matplotlib
        '''
        
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot
        
        colours = { 'mflayer': 'b',
                  'grclayer': 'g',
                  'goclayer': 'r',
                  'mfgrcsynapsis': 'c',
                  'mfgocsynapsis': 'm',
                  'grcgocsynapsis': 'y',
                  'gocgrcsynapsis': 'k',
                  'gocgocsynapsis': 'r'}
        
        ignored_layers = [
#                          'mfgrcsynapsis',
#                          'mfgocsynapsis', 
#                          'grcgocsynapsis',
#                          'gocgrcsynapsis', 
#                          'gocgocsynapsis'
                          ]
        
        figure = matplotlib.pyplot.figure(figsize=[23,14],dpi=80)
        axes = figure.add_subplot(111, projection='3d')
        axes.set_xlabel('X (mm)')
        axes.set_ylabel('Y (mm)')
        axes.set_zlabel('Z (mm)')
        
        axes.set_xlim([0.0,self.network_size[0]])
        axes.set_ylim([0.0,self.network_size[1]])
        axes.set_zlim([0.0,self.network_size[2]])
        
        # Visualize the neuron layers
        for neuron_layer in self.neuron_layers:
            if neuron_layer.__name__ not in ignored_layers:
                logger.debug('Plotting neurons from layer %s',neuron_layer.__name__)
                if neuron_layer.soma_size is None:
                    soma_size = 5
                else:
                    soma_size = neuron_layer.soma_size
            
                neuron_positions = neuron_layer.get_absolute_coordinates()
                
                axes.scatter(neuron_positions[:,0], neuron_positions[:,1], neuron_positions[:,2], \
                             c=colours[neuron_layer.__name__], marker = 'o', s = soma_size*1e4, label=neuron_layer.__name__)
        
        # Visualize the synaptic layers
        for synaptic_layer in self.synaptic_layers:
            if synaptic_layer.__name__ not in ignored_layers:
                logger.debug('Plotting synaptic connections from layer %s',synaptic_layer.__name__)
                
                source_layer_positions = synaptic_layer.source_layer.get_absolute_coordinates()
                target_layer_positions = synaptic_layer.target_layer.get_absolute_coordinates()
                
                for spos, tpos in zip(source_layer_positions[synaptic_layer.source_index],target_layer_positions[synaptic_layer.target_index]):
                    axes.plot([spos[0],tpos[0]],[spos[1],tpos[1]],[spos[2],tpos[2]],\
                              marker = '', color = colours[synaptic_layer.__name__], linewidth=0.2, linestyle = 'solid', label=synaptic_layer.__name__)
    
        logger.debug('Cleaning figure legend')
            
        # Remove repeated legend entries
        handles, labels = axes.get_legend_handles_labels()
        newLabels, newHandles = [], []
        for handle, label in zip(handles, labels):
            if label not in newLabels:
                newLabels.append(label)
                newHandles.append(handle)
        matplotlib.pyplot.legend(newHandles, newLabels, loc='center right', bbox_to_anchor=(1.1,0.5))
        
        logger.debug('Showing the figure')
        matplotlib.pyplot.show()
            
            
            
    
    def get_number_of_elements(self, **kwargs):
        '''
        Retrieve the number of elements (neurons or synapses in a layer.
        @param layer The name of the layer to get the number of elements.
        '''
        
        if 'layer' in kwargs:
            layer_name = kwargs.pop('layer')
        else:
            logger.error('Non-specified layer name')
            raise Exception('Non-DefinedLayerName')
        
        if layer_name in self.layer_map:
            layer = self.layer_map[layer_name]
            
            if (isinstance(layer, InputLayer.InputLayer)):
                return layer.number_of_neurons
            elif (isinstance(layer, SynapticLayer.SynapticLayer)):
                return layer.number_of_synapses
            else:
                logger.error('Invalid layer type in %s',layer.__name__)
                raise Exception('InvalidLayerType')
        else:
            logger.error('%s does not exist in the model',layer_name)
            raise Exception('InvalidLayerName')
    
    @abc.abstractmethod
    def update_network_weights(self):
        '''
        Update the weights of the cerebellar model according to the running implementation.
        '''
        return
    
    @abc.abstractmethod
    def update_neuron_states(self):
        '''
        Update the values of the state variable indicated
        '''
        return
    
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
    def get_synaptic_connections(self, **kwargs):
        '''
        Get the source and target neurons in the synaptic layer.
        '''
    
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

    @abc.abstractmethod
    def get_local_py_rng(self):
        '''
        Return the random number generator of this proccess
        '''
        return
    
    @abc.abstractmethod
    def get_global_py_rng(self):
        '''
        Return the random number generator for serial operations
        '''
        return
        
def _search_hdf5_group(file, name):
    '''
    This function searches the first element in a hdf5 file with the specified name. If it cannot find any element it returns None
    @param file The hdf5 file to explore.
    @param name The name to search.
    '''
    
    if file is None:
        return None
    
    import h5py
    
    for _, elem in file.iteritems():
        if 'name' in elem.attrs and elem.attrs['name']==name:
            return elem
    
    return None
