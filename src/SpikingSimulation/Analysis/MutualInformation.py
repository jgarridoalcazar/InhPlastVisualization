import numpy
import math
import Analysis
import logging
from mpi4py import MPI
        
logger = logging.getLogger('Simulation')

class MutualInformation(Analysis.Analysis):
    '''
    This class defines the calcule of the mutual information between the
    patterns and the cellular activity.
    '''
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates an object.
        @param data_provider The dataProvider that will be used in the axes to get the data. Obligatory parameter.
        @param pattern_generator the pattern generator to get the pattern happening at every time step. Obligatory parameter.
        @param layer: Name of the layer to analyze. Obligatory parameter.
        @param cell_index: Indexes of the cells to analyze. Optional parameter.
        @param pattern_index: Indexes of the patterns to analyze. Optional parameter.
        @param window_length: Length of the window to perform the analysis. Optinal parameter. If window_lenght is not specified, the simulation length will be used.
        @param time_bin: Discretization time bin for the mutual information analysis. Obligatory parameter
        '''
                
        # Get data_provider parameter 
        if ('data_provider' in kwargs):
            self.data_provider = kwargs.pop('data_provider',None)
        else:
            logger.error('Obligatory data_provider parameter not provided')
            raise Exception('NonProvidedParameter','data_provider')
        
        # Get property name parameter 
        if ('pattern_generator' in kwargs):
            self.pattern_generator = kwargs.pop('pattern_generator',None)
        else:
            logger.error('Obligatory pattern generator not provided')
            raise Exception('NonProvidedParameter','pattern_generator')
        
        # Get layer name parameter
        if ('layer' in kwargs):
            self.layer = kwargs.pop('layer',None)
        else:
            logger.error('Obligatory layer parameter not provided')
            raise Exception('NonProvidedParameter','layer')
        
        if self.layer not in self.data_provider.layer_map:
            logger.error('Invalid cell layer %s in mutual information analysis',self.layer)
            raise Exception('InvalidParameter','layer')
        
        # Get index parameter
        if ('cell_index' in kwargs):
            self.cell_index = kwargs.pop('cell_index',None)
        else:
            self.cell_index = range(self.data_provider.get_number_of_elements(layer=self.layer))
            
        # Get index parameter
        if ('pattern_index' in kwargs):
            self.pattern_index = kwargs.pop('pattern_index',None)
        else:
            self.pattern_index = range(1,self.pattern_generator.number_of_patterns+1)
        
        # Get window_length parameter
        if ('window_length' in kwargs):
            self.window_length = kwargs.pop('window_length',None)
        else:
            self.window_length = None
            
        # Get the time bin
        if 'time_bin' in kwargs:
            self.time_bin = kwargs.pop('time_bin', None)
        else:
            logger.error('Obligatory time_bin not provided')
            raise Exception('NonProvidedParameter','time_bin')
            
        
        super(MutualInformation, self).__init__(**kwargs)
                
                
    def initialize(self):
        '''
        Perform all the required operations needed in order to initialize the analysis.
        '''
        
        super(MutualInformation, self).initialize()
        
        # Time of the last update
        self.data_update = 0
        
        # Time bin inverse
        self.inv_time_bin = 1./self.time_bin
        
        # Generate the time bin matrix
        self.total_time = self.pattern_generator.simulation_time
        bin_time_init = numpy.linspace(0.0, self.total_time-self.time_bin, num=self.total_time*self.inv_time_bin)
        bin_time_end = numpy.linspace(self.time_bin,self.total_time,num=self.total_time*self.inv_time_bin)
        
        # Initialize a matrix
        self.num_patterns = len(self.pattern_index)
        self.num_bins = len(bin_time_init)
        self.num_cells = len(self.cell_index)
        
        # Calculate the time of each pattern interval
        time_end_of_pattern = self.pattern_generator.pattern_length_cum
        time_init_of_pattern = numpy.append([0.0],time_end_of_pattern[:-1])
        
        # Calculate the bin of each pattern interval. Check the round of the last bin to avoid out of range
        bin_end_of_pattern = numpy.floor(time_end_of_pattern * self.inv_time_bin).astype(int)
        if (bin_end_of_pattern[-1]>=self.num_bins):
            bin_end_of_pattern[-1]=self.num_bins-1
        bin_init_of_pattern = numpy.append([0],bin_end_of_pattern[:-1])
        if (bin_init_of_pattern[-1]>=self.num_bins):
            bin_init_of_pattern[-1]=self.num_bins-1
        
        # Final matrix indicating which bins are considered of each pattern
        self.bin_is_pattern = numpy.empty((self.num_patterns, self.num_bins),dtype='bool')
        self.bin_is_pattern[:,:] = False
        self.bin_pattern = numpy.zeros(self.num_bins)
        self.bin_pattern[:] = 0
        
        # Final matrix indicating which bins are registered spikes
        self.bin_has_fired = numpy.empty((self.num_cells, self.num_bins), dtype='bool')
        self.bin_has_fired[:,:] = False
        
        for key, value in enumerate(self.pattern_index):
            if (value!=0):
                time_of_pattern_in_bin = numpy.zeros(self.num_bins)
                for index in self.pattern_generator.pattern_id_index[value]:
                    init_bin = bin_init_of_pattern[index]
                    end_bin = bin_end_of_pattern[index]
                
                    list_of_bins = range(init_bin,end_bin+1)
                
                    # Add the time of the initial bin (if exists)
                    if init_bin!=end_bin:
                        time_of_pattern_in_bin[list_of_bins[0]] += (bin_time_end[init_bin] - time_init_of_pattern[index])
                
                    # Add the time of the intermediate bins (if exist)
                    time_of_pattern_in_bin[list_of_bins[1:-1]] += self.time_bin
                
                    # Add the time of the final bin
                    time_of_pattern_in_bin[list_of_bins[-1]] += (time_end_of_pattern[index] - max(time_init_of_pattern[index],bin_time_init[end_bin]))
            
                # Those bins where the time in the pattern is longer than half of the bin are set to part of that pattern    
                self.bin_is_pattern[key,time_of_pattern_in_bin>(self.time_bin/2.)] = True
                self.bin_pattern[self.bin_is_pattern[key,:]] = value
        
        # Create a map of cells to index
        self.cell_map = dict()
        for key, value in enumerate(self.cell_index):
            self.cell_map[value] = key
            
        # Initialize the mutual information
        self.mutual_information = 0.0
                
        return
    
    def runAtTime(self, simulation_time):
        '''
        This function updates the mutual information analysis.
        @param simulation_time: The simulation end time (in seconds).
        '''
        
        if simulation_time<self.window_length:
            init_time = 0
        else:
            init_time = simulation_time - self.window_length
        
        
        load_data_init = max(init_time, self.data_update)
        
        # Load data from the data provider
        gtime,gcell_id = self.data_provider.get_spike_activity(neuron_layer = self.layer, neuron_indexes = self.cell_index,\
                                                               init_time = load_data_init, end_time = simulation_time)
        
        
        comm = MPI.COMM_WORLD
        
        process_id = comm.Get_rank()
        
        if (process_id==0):
            
            spike_bin_index = numpy.floor(gtime*self.inv_time_bin).astype(int)
            
            self.av_firing_rate = float(len(gtime))/(len(self.cell_index)*(simulation_time-load_data_init))
            logger.info('Average firing rate in MI analysis: %sHz',str(self.av_firing_rate))
            
            cell_index = numpy.array([self.cell_map[value] for value in gcell_id])
            
            # Final matrix indicating which bins are registered spikes
            if len(spike_bin_index) and len(cell_index):
                self.bin_has_fired[cell_index,spike_bin_index] = True
            
            # Calculate mutual information in the time window
            init_bin = int(init_time * self.inv_time_bin)
            end_bin = int(simulation_time * self.inv_time_bin)
            
#             # Calculate probability of cell response
#             cell_response_prob = numpy.sum(self.bin_has_fired[:,init_bin:end_bin],axis=1) / float(end_bin-init_bin)
#             
#             # Calculate probability of pattern
#             pattern_prob = numpy.sum(self.bin_is_pattern[:,init_bin:end_bin],axis=1) / float(end_bin-init_bin)
            
            # Calculate hit matrix
            patterns, hit_matrix = calc_Hit_Matrix(self.bin_has_fired[:,init_bin:end_bin], self.bin_pattern[init_bin:end_bin])
            
            logger.info('Pattern list %s', str(patterns))
            logger.info('Hit matrix')
            logger.info(str(hit_matrix))
            
            # Calculate the firing state of the cell population
            cell_state = calc_Firing_State(self.bin_has_fired[:,init_bin:end_bin])
            
            # Calculate the mutual information
            self.mutual_information, self.max_mutual_information = calc_MI(self.bin_pattern[init_bin:end_bin],cell_state)
            
            logger.info('Mutual information')
            logger.info(str(self.mutual_information))
            
            logger.info('Theoretical maximum of MI')
            logger.info(str(self.max_mutual_information))

            
#             # Calculate probability of hits
#             hit_probability = numpy.array([numpy.sum(self.bin_has_fired[:,init_bin:end_bin] & self.bin_is_pattern[pattern, init_bin:end_bin],axis=1)/float(end_bin-init_bin) \
#                                for pattern in range(len(self.pattern_index))])
#             
#             logger.debug('Hit probability matrix (pattern in rows, cell in columns)')
#             logger.debug(str(hit_probability))
#             
#             
#             # Calculate probability of missed
#             miss_probability = numpy.array([numpy.sum(~self.bin_has_fired[:,init_bin:end_bin] & self.bin_is_pattern[pattern, init_bin:end_bin],axis=1)/float(end_bin-init_bin) \
#                                for pattern in range(len(self.pattern_index))])
#             logger.debug('Misses probability matrix (pattern in rows, cell in columns)')
#             logger.debug(str(miss_probability))
#             
#             
#             # Calculate probability of false alarm
#             false_probability = numpy.array([numpy.sum(self.bin_has_fired[:,init_bin:end_bin] & ~self.bin_is_pattern[pattern, init_bin:end_bin],axis=1)/float(end_bin-init_bin) \
#                                for pattern in range(len(self.pattern_index))])
#             logger.debug('False-alarm probability matrix (pattern in rows, cell in columns)')
#             logger.debug(str(false_probability))
#             
#              
#             # Calculate probability of correct rejection
#             rejection_probability = numpy.array([numpy.sum(~self.bin_has_fired[:,init_bin:end_bin] & ~self.bin_is_pattern[pattern, init_bin:end_bin],axis=1)/float(end_bin-init_bin) \
#                                for pattern in range(len(self.pattern_index))])
#             logger.debug('Correct rejection probability matrix (pattern in rows, cell in columns)')
#             logger.debug(str(rejection_probability))
#             
#             logger.debug('Pattern probability')
#             logger.debug(str(pattern_prob))
#             
#             logger.debug('Cell response probability')
#             logger.debug(str(cell_response_prob))
#             
#             # Print the header of the probability 
#             
#             for key_pat in range(len(self.pattern_index)):
#                 pat_prob = pattern_prob[key_pat]
#                 not_pat_prob = 1. - pat_prob
#                     
#                 for key_cell in range(len(self.cell_index)):
#                     
#                     hit = hit_probability[key_pat,key_cell]
#                     miss = miss_probability[key_pat,key_cell]
#                     false = false_probability[key_pat,key_cell]
#                     reject = rejection_probability[key_pat,key_cell]
#                     
#                     cell_prob = cell_response_prob[key_cell]
#                     not_cell_prob = 1. - cell_prob
#                     
#                     # Calculate the mutual information
#                     self.mutual_information[key_pat,key_cell] = 0
#                     if hit:
#                         self.mutual_information[key_pat,key_cell] += hit*math.log(hit/(pat_prob*cell_prob),2)
#                     if miss:
#                         self.mutual_information[key_pat,key_cell] += miss*math.log(miss/(pat_prob*not_cell_prob),2)
#                     if false:
#                         self.mutual_information[key_pat,key_cell] += false*math.log(false/(not_pat_prob*cell_prob),2)
#                     if reject:
#                         self.mutual_information[key_pat,key_cell] += reject*math.log(reject/(not_pat_prob*not_cell_prob),2)
#             
#             # Calculate maximum mutual information
#             self.max_mutual_information = -pattern_prob*numpy.log2(pattern_prob) - (1.-pattern_prob)*numpy.log2(1.-pattern_prob)
#                                                                 
#             logger.info('Mutual information matrix (pattern in rows, cell in columns)')
#             logger.info(str(self.mutual_information))
#             
#             logger.info('Theoretical maximum of MI (pattern in columns)')
#             logger.info(str(self.max_mutual_information))
            
        self.data_update = simulation_time
        
        return 
    
    def writeToFile(self, file_name):
        '''
        This function writes the results into a file.
        @param file_name Name of the file where the data will be stored
        '''
        
        numpy.savetxt(file_name, [self.mutual_information, self.av_firing_rate], delimiter='\t', newline='\n')

def calc_Hit_Matrix(cell_firing, bin_pattern):
    '''
    Calculate the hit matrix with 1 line for each cell and 1 column for each pattern (including noise).
    @param cell_firing Boolean matrix including 1 line for each cell and 1 column for each time bin.
    @param bin_pattern Array with the index of the pattern for each bin.
    '''
    patterns = numpy.unique(bin_pattern)
    
    hit_matrix = numpy.empty((len(cell_firing),len(patterns)))
    
    for index, pat in enumerate(patterns):
        pattern_bin_index = (bin_pattern==pat)
        
        for cell_index in range(len(cell_firing)):
            hit_matrix[cell_index,index] = numpy.count_nonzero(cell_firing[cell_index,pattern_bin_index])/float(numpy.count_nonzero(pattern_bin_index))   
    
    return patterns, hit_matrix
    
    
     

def calc_Firing_State(fired_matrix):
    '''
    Calculate the state of a cell population based of the cells firing in the bin.
    @param fired_matrix Boolean matrix indicating whether a cell (in rows) fired in a time bin (in columns)
    '''
    
    num_bins = fired_matrix.shape[1]
    
    accumulated = numpy.zeros(num_bins)
    for index, row in enumerate(fired_matrix):
        accumulated += numpy.left_shift(row,index)
        
    return accumulated
        
def calc_MI(X,Y):
    '''
    Calculate the MI between the pattern in the input and the outout in some neuron layer.
    @param X Array with the number of pattern (including background noise) in each time bin.
    @param Y Combination of cells firing in each time bin.
    '''
    patterns = numpy.unique(X)
    cell_states = numpy.unique(Y)
    
    patterns_edges = numpy.append(patterns-0.5,[patterns[-1]+0.5])
    cell_states_edges = numpy.append(cell_states-0.5,[cell_states[-1]+0.5])
     
    c_XY = numpy.histogram2d(X,Y,[patterns_edges, cell_states_edges])[0]
    c_X = numpy.histogram(X,patterns_edges)[0]
    c_Y = numpy.histogram(Y,cell_states_edges)[0]

    H_X = shan_entropy(c_X)
    logger.debug('Shannon entropy of the patterns: %s', H_X)
    
    H_Y = shan_entropy(c_Y)
    logger.debug('Shannon entropy of the population response: %s', H_Y)
    
    H_XY = shan_entropy(c_XY)
    logger.debug('Joint shannon entropy: %s', H_XY)
    
    MI = H_X + H_Y - H_XY
    return MI, H_X

def shan_entropy(c):
    '''
    Calculate Shannon Entropy of the histogram.
    @param c Histogram of the occurrences of the events.
    '''
    c_normalized = c/float(numpy.sum(c))
    c_normalized = c_normalized[numpy.nonzero(c_normalized)]
    H = -sum(c_normalized*numpy.log2(c_normalized))  
    return H
