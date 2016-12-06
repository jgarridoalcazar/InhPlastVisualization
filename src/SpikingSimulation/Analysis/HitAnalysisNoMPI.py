import numpy
import math
import Analysis
import logging
import os
        
logger = logging.getLogger('Simulation')

class HitAnalysis(Analysis.Analysis):
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
            
        
        super(HitAnalysis, self).__init__(**kwargs)
                
                
    def initialize(self):
        '''
        Perform all the required operations needed in order to initialize the analysis.
        '''
        
        super(HitAnalysis, self).initialize()
        
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
        self.av_selected = 0.0
        self.av_nonselected = 0.0
        self.hit_index = 0.0
                
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
        cr_matrix, hit_matrix, miss_matrix, fa_matrix, faa_matrix = calc_Ind_Pattern_Hit_Matrix(self.bin_has_fired[:,init_bin:end_bin], self.bin_is_pattern[:,init_bin:end_bin])
        
        logger.info('Individual pattern hit matrix:')
        logger.info('%s', hit_matrix)
        logger.info('Individual pattern correct rejection matrix:')
        logger.info('%s', cr_matrix)
        logger.info('Individual pattern miss matrix:')
        logger.info('%s', miss_matrix)
        logger.info('Individual pattern false alarm matrix:')
        logger.info('%s', fa_matrix)
        logger.info('Individual pattern false alarm any matrix:')
        logger.info('%s', faa_matrix)
        
        extended_hit = numpy.append(hit_matrix, [faa_matrix], axis=0)
        selected_pattern = numpy.argmax(extended_hit, axis=0)

        # Select those cells mainly responding to one of the patterns
        ind = numpy.where(selected_pattern<self.num_patterns)

        diagonal = hit_matrix[selected_pattern[ind],ind]

        # Create a new matrix setting the diagonal values to zero
        zero_matrix = hit_matrix
        zero_matrix[selected_pattern[ind],ind] = 0.0

        # Calculate the average of the diagonal elements
        av_selected = numpy.average(diagonal)
        av_nonselected = numpy.sum(zero_matrix[:])/(self.num_cells*self.num_patterns)
        logger.debug(str('Average hits of selected elements:')+str(av_selected))
        logger.debug(str('Average hits of non-selected elements:')+str(av_nonselected))
        hit_index = av_selected - av_nonselected
        logger.debug(str('Hit index:')+str(hit_index))
        
        self.av_selected = av_selected
        self.av_nonselected = av_nonselected
        self.hit_index = hit_index 
        
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
    faa_matrix = numpy.empty((len(cell_firing)))
    
    for index_pat, pattern in enumerate(pattern_present):
        if (numpy.count_nonzero(pattern)):
            for index_cell, firing in enumerate(cell_firing):
                hit_matrix[index_pat,index_cell] = numpy.count_nonzero(firing&pattern)/float(numpy.count_nonzero(pattern))
                cr_matrix[index_pat,index_cell] = numpy.count_nonzero(~firing&~pattern)/float(numpy.count_nonzero(~pattern))
                miss_matrix[index_pat,index_cell] = numpy.count_nonzero(~firing&pattern)/float(numpy.count_nonzero(pattern))
                fa_matrix[index_pat,index_cell] = numpy.count_nonzero(firing&~pattern)/float(numpy.count_nonzero(~pattern))
        else:
            logger.warning('Pattern %s never occurs. Statistics will not be calculated', index_pat)
            
    any_pattern = numpy.any(pattern_present, axis=0)
    for index_cell, firing in enumerate(cell_firing):
        faa_matrix[index_cell] = numpy.count_nonzero(firing&~any_pattern)/float(numpy.count_nonzero(~any_pattern))
        
    return cr_matrix, hit_matrix, miss_matrix, fa_matrix, faa_matrix
     

