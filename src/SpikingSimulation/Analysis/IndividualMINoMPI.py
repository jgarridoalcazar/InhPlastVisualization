import numpy
import math
import Analysis
import logging
import os
        
logger = logging.getLogger('Simulation')

class IndividualMI(Analysis.Analysis):
    '''
    This class defines the calcule of the mutual information (average of individual neurons) 
    between a single pattern and the cellular activity.
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
            self.pattern_index = range(0,self.pattern_generator.number_of_patterns)
        
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
            
        
        super(IndividualMI, self).__init__(**kwargs)
                
                
    def initialize(self):
        '''
        Perform all the required operations needed in order to initialize the analysis.
        '''
        
        super(IndividualMI, self).initialize()
        
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
        
        # Final matrix indicating which bins are registered spikes
        self.bin_has_fired = numpy.empty((self.num_cells, self.num_bins), dtype='bool')
        self.bin_has_fired[:,:] = False
        
        for key, value in enumerate(self.pattern_index):
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
            self.bin_pattern[self.bin_is_pattern[key,:]] = self.bin_pattern[self.bin_is_pattern[key,:]] + 2**value
    
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
        
        logger.info('Analyzing mutual information from time %s to %s seconds',str(init_time),str(simulation_time))
        
        load_data_init = max(init_time, self.data_update)
        
        # Load data from the data provider
        gtime,gcell_id = self.data_provider.get_spike_activity(neuron_layer = self.layer, neuron_indexes = self.cell_index,\
                                                               init_time = load_data_init, end_time = simulation_time)
        
        
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
        # Calculate hit matrix for each pattern
        #logger.debug('Firing matrix:')
        #logger.debug('%s', self.bin_has_fired[:,init_bin:end_bin].tolist())
        #logger.debug('Pattern matrix:')
        #logger.debug('%s', self.bin_is_pattern[0,init_bin:end_bin].tolist())
        cr_matrix, hit_matrix, miss_matrix, fa_matrix = calc_Ind_Pattern_Hit_Matrix(self.bin_has_fired[:,init_bin:end_bin], self.bin_is_pattern[:,init_bin:end_bin])
        
        logger.info('Individual pattern hit matrix:')
        logger.info('%s', hit_matrix)
        logger.info('Individual pattern correct rejection matrix:')
        logger.info('%s', cr_matrix)
        logger.info('Individual pattern miss matrix:')
        logger.info('%s', miss_matrix)
        logger.info('Individual pattern false alarm matrix:')
        logger.info('%s', fa_matrix)
        

        # Calculate pattern entropy, cell entropy and joint entropy
        pat_entropy = calc_Entropy(self.bin_is_pattern[:,init_bin:end_bin])
        cell_entropy = calc_Entropy(self.bin_has_fired[:,init_bin:end_bin])
        joint_entropy = calc_Joint_Entropy(self.bin_has_fired[:,init_bin:end_bin], self.bin_is_pattern[:,init_bin:end_bin])

        logger.debug('Shannon entropy of the patterns: %s', pat_entropy)
    
        logger.debug('Shannon entropy of the population response: %s', cell_entropy)
    
        logger.debug('Joint shannon entropy: %s', joint_entropy)

        # Calculate average MI
        Av_MI = calc_MI(pat_entropy, cell_entropy, joint_entropy)

        logger.debug('Average MI of individual cells: %s', Av_MI)

        self.mutual_information, self.max_mutual_information = Av_MI,pat_entropy

        logger.info('Mutual information')
        logger.info('%s',self.mutual_information)
        
        logger.info('Theoretical maximum of MI')
        logger.info('%s',self.max_mutual_information)
           
        self.data_update = simulation_time
        
        return 
    
    def writeToFile(self, file_name):
        '''
        This function writes the results into a file.
        @param file_name Name of the file where the data will be stored
        '''
        
        dirname = os.path.dirname(os.path.abspath(file_name));
        if (not os.path.isdir(dirname)):
            logger.debug(str("Creating directory:")+str(dirname))
            os.makedirs(dirname)
            
        numpy.savetxt(file_name, [self.mutual_information, self.av_firing_rate], delimiter='\t', newline='\n')

def calc_Entropy(bin_matrix):
    '''
    Calculate the entropy for each individual pattern.
    @param bin_matrix Boolean matrix including 1 line for each pattern and 1 column for each time bin.
    '''
    sum_bins = numpy.count_nonzero(bin_matrix,axis=1)
    probability = sum_bins/float(bin_matrix.shape[1])
    entropy = numpy.zeros(probability.shape)
    idx = probability>0.0
    entropy[idx] = -probability[idx]*numpy.log2(probability[idx])
    entropy[idx] -= (1-probability[idx])*numpy.log2(1-probability[idx])
    return entropy


def calc_Joint_Entropy(cell_firing, pattern_present):
    '''
    Calculate the entropy of having both cell firing and pattern present.
    @param cell_firing Boolean matrix including 1 line for each cell and 1 column for each time bin.
    @param pattern_present Boolean matrix including 1 line for each pattern and 1 column for each time bin.
    '''
    hit_matrix = numpy.empty((len(pattern_present),len(cell_firing)))
    cr_matrix = numpy.empty((len(pattern_present),len(cell_firing)))
    miss_matrix = numpy.empty((len(pattern_present),len(cell_firing)))
    fa_matrix = numpy.empty((len(pattern_present),len(cell_firing)))
    
    for index_pat, pattern in enumerate(pattern_present):
        if (numpy.count_nonzero(pattern)):
            for index_cell, firing in enumerate(cell_firing):
                hit_matrix[index_pat,index_cell] = numpy.count_nonzero(firing&pattern)/float(pattern.shape[0])
                cr_matrix[index_pat,index_cell] = numpy.count_nonzero(~firing&~pattern)/float(pattern.shape[0])
                miss_matrix[index_pat,index_cell] = numpy.count_nonzero(~firing&pattern)/float(pattern.shape[0])
                fa_matrix[index_pat,index_cell] = numpy.count_nonzero(firing&~pattern)/float(pattern.shape[0])
        else:
            logger.warning('Pattern %s never occurs. Statistics will not be calculated', index_pat)
    
    entropy = numpy.zeros((len(pattern_present),len(cell_firing)),dtype=numpy.float)           
    idx = hit_matrix>0
    entropy[idx] = -hit_matrix[idx]*numpy.log2(hit_matrix[idx])
    idx = cr_matrix>0
    entropy[idx] -= cr_matrix[idx]*numpy.log2(cr_matrix[idx])
    idx = miss_matrix>0
    entropy[idx] -= miss_matrix[idx]*numpy.log2(miss_matrix[idx])
    idx = fa_matrix>0
    entropy[idx] -= fa_matrix[idx]*numpy.log2(fa_matrix[idx])
    return entropy

def calc_Ind_Pattern_Hit_Matrix(cell_firing, pattern_present):
    '''
    Calculate the correct rejection, hit, miss and false alarm matrisses with 1 line for each cell and 1 column for each pattern (including noise).
    @param cell_firing Boolean matrix including 1 line for each cell and 1 column for each time bin.
    @param pattern_present Boolean matrix including 1 line for each pattern and 1 column for each time bin.
    '''
    hit_matrix = numpy.empty((len(pattern_present),len(cell_firing)))
    cr_matrix = numpy.empty((len(pattern_present),len(cell_firing)))
    miss_matrix = numpy.empty((len(pattern_present),len(cell_firing)))
    fa_matrix = numpy.empty((len(pattern_present),len(cell_firing)))
    
    for index_pat, pattern in enumerate(pattern_present):
        if (numpy.count_nonzero(pattern)):
            for index_cell, firing in enumerate(cell_firing):
                hit_matrix[index_pat,index_cell] = numpy.count_nonzero(firing&pattern)/float(numpy.count_nonzero(pattern))
                cr_matrix[index_pat,index_cell] = numpy.count_nonzero(~firing&~pattern)/float(numpy.count_nonzero(~pattern))
                miss_matrix[index_pat,index_cell] = numpy.count_nonzero(~firing&pattern)/float(numpy.count_nonzero(pattern))
                fa_matrix[index_pat,index_cell] = numpy.count_nonzero(firing&~pattern)/float(numpy.count_nonzero(~pattern))
        else:
            logger.warning('Pattern %s never occurs. Statistics will not be calculated', index_pat)
               
        
    return cr_matrix, hit_matrix, miss_matrix, fa_matrix 
    
def calc_MI(pattern_entropy, cell_entropy, joint_entropy):
    '''
    Calculate the MI between the pattern in the input and each neuron.
    @param pattern_entropy Array with the entropy of each individual pattern.
    @param cell_entropy Array with the entropy of each individual cell.
    @param joint_entropy 2D array with the joint entropy of each cell to each pattern
    '''
    
    MI_matrix = numpy.empty((len(pattern_entropy),len(cell_entropy)))

    for index_pat, pattern_H in enumerate(pattern_entropy):
        MI_matrix[index_pat,:] = pattern_H + cell_entropy[:] - joint_entropy[index_pat,:]

    Av_MI = numpy.average(MI_matrix,axis=1)
    
    return Av_MI

