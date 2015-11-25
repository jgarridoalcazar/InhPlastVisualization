'''
Created on April 23, 2014

@author: Jesus Garrido (jgarridoalcazar at gmail.com)
'''
import numpy
import time
import logging
import scipy.spatial.distance

logger = logging.getLogger('Simulation')

class SynapticLayer(object):
    '''
    This class defines a synaptic layer and the data needed to generate it in a simulator.
    If several MPI processes are used, each process store only the loacl connections (those
    reaching a local cell).
    '''
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates a new synaptic layer.
        @param source_layer: Neuron/input source layer.
        @param target_layer: Neuron target layer.
        @param intermediate_layer: Third layer to be used to create the connections (i.e. GoC-GrC connections). (Optional)
        @param intermediate_to_target_synaptic_layer: Synaptic layer connection the intermediate_layer and the target layer. (Optional)
        @param weight_recording: Whether the weights must be recorded.
        @param weight_normalization: Whether the weights must be normalizated.
        @param weight_sum: Total sum of the incoming connection weights for each target cell.
        @param random_generator: The random number generator to use.
        @param synaptic_delay: Transmission delay in the synapses (in seconds).
        @param connectivity_type: Pattern of connectivity in the connection. Only those patterns included in template_connectivity_parameters are allowed.
        @param connectivity_parameters: Connectivity algorithm parameters as described in the template_connectivity_parameters dictionary.
        @param weight_initialization_type: Weight initialization type in the connection. Only those patterns included in template_weight_parameters are allowed.
        @param weight_initialization_parameters: Weight initialization parameters as described in the template_weight_parameters dictionary.
        @param learning_rule_type: Learning rule type associated to the connections. Only those rules included in template_rule_parameters are allowed.
        @param learning_rule_parameters: Learning rule parameters as described in the template_rule_parameters dictionary.
        '''
        
        # Specific parameters for each learning rule
        self.template_rule_parameters = {
                         'STDP': ['tau_plus','learning_step','minus_plus_ratio','max_weight'],
                         'STDPSym': ['tau_sym','learning_step','minus_plus_ratio','max_weight'],
                         'eSTDP': ['tau_plus','learning_step','minus_plus_ratio','max_weight'],
                         'iSTDP': ['tau_plus','learning_step','minus_plus_ratio','max_weight']
                         }
    
        # Specific parameters of the connectivity algorithm (and function implementing each one)
        self.template_connectivity_parameters = {
                                        'fixedn2one': ['number_of_source_cells', self._generate_fixed_n2one_connections_],
                                        'randomn2one': ['number_of_source_cells', self._generate_random_n2one_connections_],
                                        'random_with_probability': ['connection_probability', self._generate_random_connections_],
                                        'randomn2onelocal': ['average_number_of_source_cells', 'average_dendritic_length', self._generate_random_n2one_local_connections_],
                                        'glomerulin2onelocal': ['average_number_of_source_cells_per_glomerulus', 'average_dendritic_length', self._generate_glomeruli_n2one_local_connections_]
                                        }

        # Specific parameters of the weight initialization (and function implementing each one)
        self.template_weight_parameters = {
                                        'random': ['random_min_weight', 'random_max_weight', self._generate_random_weights_],
                                        'fixed': ['initial_weight', self._generate_fixed_weights_]
                                        }
        # Read name
        self.__name__ = kwargs.pop('name')
        
        # Read source layer
        if ('source_layer' in kwargs):
            self.source_layer = kwargs.pop('source_layer')
        else:
            logger.error('Non-specified source layer')
            raise Exception('Non-DefinedProperty')
        
        # Read target layer
        if ('target_layer' in kwargs):
            self.target_layer = kwargs.pop('target_layer')
        else:
            logger.error('Non-specified target layer')
            raise Exception('Non-DefinedProperty')
        
        # Read intermediate layer
        if ('intermediate_layer' in kwargs):
            self.intermediate_layer = kwargs.pop('intermediate_layer')
        else:
            self.intermediate_layer = None
            
        # Read intermediate_to_target_ synaptic layer
        if ('intermediate_to_target_synaptic_layer' in kwargs):
            self.intermediate_to_target_synaptic_layer = kwargs.pop('intermediate_to_target_synaptic_layer')
        else:
            self.intermediate_to_target_synaptic_layer = None
        
        # Read seed
        if ('random_generator' in kwargs):
            self.random_generator = kwargs.pop('random_generator')
        else:
            logger.warning('Non-specified random generator for synaptic layer generation. Using a new one')
            self.random_generator = numpy.random.RandomState(time())
        
        # Read synaptic type
        if ('weight_recording' in kwargs):
            self.weight_recording = kwargs.pop('weight_recording')
        else:
            self.weight_recording = False
        
        # Read synaptic type
        if ('weight_normalization' in kwargs):
            self.weight_normalization = kwargs.pop('weight_normalization')
            
            if self.weight_normalization:
                if ('weight_sum' in kwargs):
                    self.weight_total_sum = kwargs.pop('weight_sum')
                else:
                    logger.warning('Non-specified weight sum for normalization. Using default value 1e-9')
                    self.weight_total_sum = 1.e-9            
        else:
            self.weight_normalization = False
            
        # Read synaptic delay
        if ('synaptic_delay' in kwargs):
            self.synaptic_delay = kwargs.pop('synaptic_delay')
        else:
            logger.error('Non-specified synaptic delay')
            raise Exception('Non-DefinedProperty')
        
        # Read connectivity type and its properties
        if ('connectivity_type' in kwargs):
            self.connectivity_type = kwargs.pop('connectivity_type')
            if (self.connectivity_type in self.template_connectivity_parameters):
                self.connectivity_parameters = {}
                for param in self.template_connectivity_parameters[self.connectivity_type][:-1]:
                    if (param in kwargs):
                        self.connectivity_parameters[param] = kwargs.pop(param)
                    else:
                        logger.warning('Non-specified connectivity parameter: %s in layer %s. Using default value', param, self.__name__)
                
                if 'allow_auto_connection' in kwargs:
                    self.connectivity_parameters['allow_auto_connection'] = kwargs.pop('allow_auto_connection')
                else:
                    self.connectivity_parameters['allow_auto_connection'] = True
                    
                if 'allow_multiple_connections' in kwargs:
                    self.connectivity_parameters['allow_multiple_connections'] = kwargs.pop('allow_multiple_connections')
                else:
                    self.connectivity_parameters['allow_multiple_connections'] = True
                
                # Generate the individual connections
                self.template_connectivity_parameters[self.connectivity_type][-1]()
            else:
                logger.error('Unknown connectivity type: %s', self.connectivity_type)
                raise Exception('UnknownConnectivityType')                         
        else:
            logger.error('Non-specified connectivity type')
            raise Exception('Non-DefinedProperty')
        
        # Read weight initialization and its properties
        if ('weight_initialization_type' in kwargs):
            self.weight_initialization_type = kwargs.pop('weight_initialization_type')
            if (self.weight_initialization_type in self.template_weight_parameters):
                self.weight_initialization_parameters = {}
                for param in self.template_weight_parameters[self.weight_initialization_type][:-1]:
                    if (param in kwargs):
                        self.weight_initialization_parameters[param] = kwargs.pop(param)
                    else:
                        logger.error('Non-specified weight initialization parameter: %s',param)
                        raise Exception('Non-DefinedWeightInitializationParameter')
                # Generate the individual connections
                self.template_weight_parameters[self.weight_initialization_type][-1]()
            else:
                logger.error('Unknown weight initialization type: %s', self.weight_initialization_type)
                raise Exception('UnknownWeightInitializationType')                         
        else:
            logger.error('Non-specified weight initialization type')
            raise Exception('Non-DefinedProperty')
        
        # Read learning rule and its properties
        if ('learning_rule_type' in kwargs):
            self.learning_rule_type = kwargs.pop('learning_rule_type')
            if (self.learning_rule_type in self.template_rule_parameters):
                self.learning_rule_parameters = {}
                for param in self.template_rule_parameters[self.learning_rule_type]:
                    if (param in kwargs):
                        self.learning_rule_parameters[param] = kwargs.pop(param)
                    else:
                        logger.warning('Non-specified learning rule parameter: %s in layer %s. Using default value', param, kwargs['name'])
            else:
                logger.error('Unknown learning rule type: %s', self.learning_rule_type)
                raise Exception('UnknownLearningRuleType')                         
        else:
            self.learning_rule_type = None
            
        # Gather the synapses to every node
        # It is not efficient in terms of memory consumption, but it makes this class independent of NEST simulator
        #self._share_connections()
        
        # Check whether additional parameters have been used.
        for param in kwargs:
            logger.warning('Unrecognized parameter %s in layer %s',param,self.__name__)
            
        return
    
    def _generate_fixed_n2one_connections_(self):
        '''
        Generate connections between source and target layers n2one following succesive order (s0,s1,...,s(n-1) to t0, sn, s(n+1),...,s(2n-1) to t1,...).
        '''
        if (self.source_layer.number_of_neurons != self.target_layer.number_of_neurons * self.connectivity_parameters['number_of_source_cells']):
            logger.error('Invalid connectivity type. Check source and target layer neuron number')
            raise Exception('InvalidPropertyValue') 
         
        # Get the number of parallel proecesses
        size = 1
        rank = 0
         
        # Calculate the number of synapses this process generates dividing the number of target cells between the number of processes
        local_target_cells = range(rank,self.target_layer.number_of_neurons,size)        
         
        self.number_of_synapses = len(local_target_cells) * self.connectivity_parameters['number_of_source_cells']
         
        # Generate source cell indexes of the connections
        self.source_index = [self.connectivity_parameters['number_of_source_cells']*index for index in local_target_cells]
        self.target_index = numpy.repeat(local_target_cells,self.connectivity_parameters['number_of_source_cells']).tolist()
        return
     
    def _generate_random_n2one_connections_(self):
        '''
        Generate connections between source and target layers n2one with random selection of source cells (sr0,sr1,...,sr(n-1) to t0, srn, s(rn+1),...,s(r2n-1) to t1,...).
        '''
        # Get the number of parallel proecesses
        size = 1
        rank = 0
         
        num_source_cells = self.connectivity_parameters['number_of_source_cells']
         
        # Calculate the number of synapses this process generates dividing the number of target cells between the number of processes
        local_target_cells = range(rank,self.target_layer.number_of_neurons,size)
        self.number_of_synapses = len(local_target_cells) * num_source_cells
         
        # For each target cell generate a permutation of the source indexes
        self.source_index = [None]*self.number_of_synapses
        for i in range(len(local_target_cells)):
            self.source_index[i*num_source_cells:(i+1)*num_source_cells] = self.random_generator.permutation(numpy.arange(self.source_layer.number_of_neurons)).tolist()[0:num_source_cells]
        self.target_index = numpy.repeat(local_target_cells,num_source_cells).tolist()
         
        return
     
    def _generate_random_connections_(self):
        '''
        Generate connections between source and target layers with a given probability.
        '''
        # Get the number of parallel proecesses
        size = 1
        rank = 0
         
        probability = self.connectivity_parameters['connection_probability']
 
        # Calculate the number of synapses this process generates dividing the number of target cells between the number of processes
        local_target_cells = numpy.arange(rank,self.target_layer.number_of_neurons,size)
         
        # Generate num_source_cells * num_target_cells rand numbers
        rand_numbers = self.random_generator.rand(self.source_layer.number_of_neurons,len(local_target_cells))
         
        # Autoconnections are not allowed
        if self.source_layer == self.target_layer:
            for index in range(len(local_target_cells)):
                rand_numbers[local_target_cells[index],index] = 1.0
         
        (source_array,target_array) = numpy.where(rand_numbers<probability)
        self.source_index = source_array
        self.target_index = local_target_cells[target_array]  
        self.number_of_synapses = self.target_index.shape[0]      
        return
    
    def _generate_random_n2one_local_connections_(self):
        '''
        Generate connections between source and target layers n2one with random selection of source cells (sr0,sr1,...,sr(n-1) to t0, srn, s(rn+1),...,s(r2n-1) to t1,...).
        In this case the probability of connections is generated according to an exponential distribution with average length the average_dendritic_length parameter.
        '''
        
        target_coord = self.target_layer.get_absolute_coordinates()
        source_coord = self.source_layer.get_absolute_coordinates()
        
        n_source_cells = source_coord.shape[0]
        n_target_cells = target_coord.shape[0]
        
        lambda_val = 1./self.connectivity_parameters['average_dendritic_length']
                
        distance = scipy.spatial.distance.cdist(source_coord, target_coord, 'euclidean')
        
        # Generate the probability of connection for each pair of source/target neurons following an exponential distribution with average_dendritic_length
        probability = lambda_val * numpy.exp(-lambda_val*distance)
        # For each target neuron normalize the probability and multiply by the average number of source cells
        norm_probability = probability / numpy.sum(probability,axis=0) * self.connectivity_parameters['average_number_of_source_cells']
        
        # Generate random numbers according to uniform distribution and compare with the probability of each connection
        rand_nums = self.random_generator.uniform(0,1,probability.shape)
        
        indexes = numpy.where(rand_nums<norm_probability)
        
        self.source_index = indexes[0]
        self.target_index = indexes[1]
        self.number_of_synapses = self.target_index.shape[0]
        return
     
    def _generate_glomeruli_n2one_local_connections_(self):
        '''
        Generate connections between source and target layers n2one with random selection of source cells (sr0,sr1,...,sr(n-1) to t0, srn, s(rn+1),...,s(r2n-1) to t1,...).
        In this case the probability of connections is generated according to an exponential distribution with average length the average_dendritic_length parameter.
        '''
        
        if self.intermediate_layer is None:
            logger.error('Non-specified intermediate layer. It is required in glomeruli-type connectivity creation.')
            raise Exception('Non-DefinedProperty')
        
        if self.intermediate_to_target_synaptic_layer is None:
            logger.error('Non-specified intermediate_to_target synaptic layer. It is required in glomeruli-type connectivity creation.')
            raise Exception('Non-DefinedProperty')
        
        interm_coord = self.intermediate_layer.get_absolute_coordinates()
        source_coord = self.source_layer.get_absolute_coordinates()
        
        n_source_cells = source_coord.shape[0]
        n_interm_cells = interm_coord.shape[0]
        
        lambda_val = 1./self.connectivity_parameters['average_dendritic_length']
                
        distance = scipy.spatial.distance.cdist(source_coord, interm_coord, 'euclidean')
        
        # Generate the probability of connection for each pair of source/target neurons following an exponential distribution with average_dendritic_length
        probability = lambda_val * numpy.exp(-lambda_val*distance)
        
        # Generate random numbers according to uniform distribution and compare with the probability of each connection
        rand_nums = self.random_generator.uniform(0,1,probability.shape)
        
        norm_probability = numpy.zeros((n_source_cells, n_interm_cells))
        
        self.source_index = []
        self.target_index = []
        
        # For each target neuron normalize the probability and multiply by the average number of source cells
        for glomerulus_index in range(n_interm_cells):
            norm_probability[:,glomerulus_index] = probability[:,glomerulus_index] / numpy.sum(probability[:,glomerulus_index]) * self.connectivity_parameters['average_number_of_source_cells_per_glomerulus']
            
            connected_goc_indexes, = numpy.where(rand_nums[:,glomerulus_index]<norm_probability[:,glomerulus_index])
            
            target_grcs = numpy.unique(numpy.concatenate(_get_target_cell_indexes_([glomerulus_index], self.intermediate_to_target_synaptic_layer)))
            neighboor_glomeruli = numpy.unique(numpy.concatenate(_get_source_cell_indexes_(target_grcs, self.intermediate_to_target_synaptic_layer)))
            
            self.source_index.extend(numpy.repeat(connected_goc_indexes,target_grcs.shape[0]))
            self.target_index.extend(numpy.tile(target_grcs,connected_goc_indexes.shape[0]))
                
            # Update the probability of the following glomeruli to avoid repeating connections
            goc_indexes = numpy.tile(connected_goc_indexes,neighboor_glomeruli.shape[0])
            glomeruli_rep_indexes = numpy.repeat(neighboor_glomeruli,connected_goc_indexes.shape[0])
            probability[goc_indexes,glomeruli_rep_indexes] = 0.0
        
        self.source_index = numpy.array(self.source_index)
        self.target_index = numpy.array(self.target_index)
        self.number_of_synapses = self.target_index.shape[0]
        return
    
    def _generate_random_weights_(self):
        '''
        Generate uniformly distributed random weights between the minimum and maximum values.
        '''
         
        min_weight = self.weight_initialization_parameters['random_min_weight']
        max_weight = self.weight_initialization_parameters['random_max_weight']
         
        rand_numbers = self.random_generator.rand(self.number_of_synapses) * (max_weight-min_weight) + min_weight
         
        self.weights = rand_numbers.tolist()
        return
     
    def _generate_fixed_weights_(self):
        '''
        Initialize all the weights to the specified values.
        '''
         
        initial_weight = self.weight_initialization_parameters['initial_weight']
        self.weights = [initial_weight] * self.number_of_synapses
        return

def _calculate_distance_(source_coord, target_coord):
    '''
    Calculate the euclidean distance between the source and target coordinates
    '''
    return numpy.linalg.norm(numpy.array(source_coord)-numpy.array(target_coord))

def _get_target_cell_indexes_(indexes, interm_to_target_syn_layer):
    '''
    Extract the indexes of the target cells of those cells indicated by indexes in the interm_layer.
    '''
    return numpy.array([interm_to_target_syn_layer.target_index[interm_to_target_syn_layer.source_index==index] for index in indexes])
    
def _get_source_cell_indexes_(indexes, interm_to_target_syn_layer):
    '''
    Extract the indexes in the interm_layer of the sources cells of those cells indicated by indexes in the target_layer.
    '''
    return numpy.array([interm_to_target_syn_layer.source_index[interm_to_target_syn_layer.target_index==index] for index in indexes])
    