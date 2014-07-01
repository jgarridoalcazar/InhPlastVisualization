'''
Created on May 27, 2014

@author: Jesus Garrido (jgarridoalcazar at gmail.com)
'''

import time
import numpy
import bisect
import logging

logger = logging.getLogger('Simulation')


class FrequencyPatternGenerator(object):
    '''
    This class defines a frequency pattern generator. It will be normalized
    both in the average frequency in rows and in columns.
    '''
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates a new frequency pattern generator.
        @param seed Seed to be used for the random number generator. If not present, time will be usead instead.
        @param mean_length Average length of each pattern according to an exponential distribution (in seconds)
        @param min_amplitude Minimum current amplitude of each fiber (in A).
        @param max_amplitude Maximum current amplitude of each fiber (in A).
        @param number_of_fibers Number of fibers to be stimulated.
        @param rate_of_fibers_in_pattern Rate of fibers to be considered in the pattern.
        @param rate_of_time_in_pattern Rate of time to be considered in the pattern.
        @param number_of_patterns Number of different patterns to include during the rate_of_time_in_pattern.
        @param simulation_time Total stimulation length (in seconds).
        @param number_of_normalizations Number of iterations to perform in the normalization process.
        '''
        
        if 'seed' in kwargs:
            self.seed = kwargs.pop('seed')                         
        else:
            self.seed = time()
            
        if 'mean_length' in kwargs:
            self.mean_length = kwargs.pop('mean_length')
        else:
            self.mean_length = 0.250
            logger.warning('Non-specified activity mean length.Using default value %s', self.mean_length)
            
        if 'min_amplitude' in kwargs:
            self.min_amplitude = kwargs.pop('min_amplitude')
        else:
            self.min_amplitude = 0
            logger.warning('Non-specified stimulation minimum amplitude. Using default value %s', self.min_amplitude)
            
        if 'max_amplitude' in kwargs:
            self.max_amplitude = kwargs.pop('max_amplitude')
        else:
            self.max_amplitude = 1
            logger.warning('Non-specified stimulation maximum amplitude. Using default value %s', self.max_amplitude)
            
        if 'number_of_fibers' in kwargs:
            self.number_of_fibers = kwargs.pop('number_of_fibers')
        else:
            logger.error('Non-specified activity number of fibers')
            raise Exception('Non-DefinedNumberOfFibers')
        
        if 'rate_of_fibers_in_pattern' in kwargs:
            self.rate_of_fibers_in_pattern = kwargs.pop('rate_of_fibers_in_pattern')
        else:
            self.rate_of_fibers_in_pattern = 0.10
            logger.warning('Non-specified stimulation rate of fibers in pattern. Using default value %s', self.rate_of_fibers_in_pattern)
            
        if 'rate_of_time_in_pattern' in kwargs:
            self.rate_of_time_in_pattern = kwargs.pop('rate_of_time_in_pattern')
        else:
            self.rate_of_time_in_pattern = 0.20
            logger.warning('Non-specified stimulation rate of time in pattern. Using default value %s', self.rate_of_time_in_pattern)
            
        if 'number_of_patterns' in kwargs:
            self.number_of_patterns = kwargs.pop('number_of_patterns')
        else:
            self.number_of_patterns = 1
            logger.warning('Non-specified stimulation number of patterns. Using default value %s', self.number_of_patterns)
            
        if 'number_of_normalizations' in kwargs:
            self.number_of_normalizations = kwargs.pop('number_of_normalizations')
        else:
            self.number_of_normalizations = 10
            logger.warning('Non-specified stimulation number of normalizations. Using default value %s', self.number_of_normalizations)
            
        if 'simulation_time' in kwargs:
            self.simulation_time = kwargs.pop('simulation_time')
        else:
            logger.error('Non-specified stimulation simulation time')
            raise Exception('Non-DefinedSimulationTime')
        
        super(FrequencyPatternGenerator, self).__init__()
        
    def _generate_length(self):
        '''
        Generate the length of each stimulation interval.
        '''
        # Generate the length of every frequency (starting from 2*T/mean_time elements)
        number_of_elements = int(self.simulation_time/self.mean_length)
        total_time = 0
        while (total_time<self.simulation_time):
            number_of_elements = 2*number_of_elements
            # Round to ms the pattern length
            self.pattern_length = numpy.round(self.ran_generator.exponential(scale=self.mean_length, size=number_of_elements),decimals=3)
            total_time = self.pattern_length.sum()
        
        # Calculate the number of elements to include in the simulation
        self.pattern_length_cum = numpy.cumsum(a=self.pattern_length)
        first_index = numpy.where(self.pattern_length_cum>=self.simulation_time)[0][0]
        
        # Select those elements
        self.pattern_length = self.pattern_length[0:(first_index+1)]
        self.pattern_length_cum = self.pattern_length_cum[0:(first_index+1)]
        # Adjust the length of the last element to fit exactly the simulation length
        self.pattern_length[-1] = self.pattern_length[-1] - (self.pattern_length_cum[first_index]-self.simulation_time)
        self.pattern_length_cum[-1] = self.simulation_time
        
        return
    
    def _generate_pattern_id(self):
        '''
        Generate the array of patterns associated to each column (0 indicate no pattern).
        '''
        # If number of patterns is 0, no pattern every column
        if not self.number_of_patterns:
            self.pattern_id = numpy.array([0]*len(self.pattern_length))
            return
        
        probability = numpy.empty(self.number_of_patterns+1)
        
        probability[0] = 1-self.rate_of_time_in_pattern
        probability[1:(self.number_of_patterns+1)] = self.rate_of_time_in_pattern/self.number_of_patterns
        
        # Generate the acumulated probability and find where each random number should be inserted.
        cumsum = numpy.cumsum(probability)
        self.pattern_id = numpy.array([bisect.bisect_right(cumsum, self.ran_generator.uniform()) for _ in range(len(self.pattern_length))])
        self.pattern_id_index = numpy.array([numpy.where(self.pattern_id==num_pattern)[0] for num_pattern in range(self.number_of_patterns+1)])
        
        return
    
    def _generate_fibers_in_pattern(self):
        '''
        Generate the number of fibers that will be used for each pattern.
        '''
        
        self.number_of_selected_fibers = int(self.number_of_fibers*self.rate_of_fibers_in_pattern)
        
        self.fibers_in_pattern = numpy.array([self.ran_generator.permutation(self.number_of_fibers)[0:self.number_of_selected_fibers] for _ in range(self.number_of_patterns)])
        self.fibers_in_pattern.sort(axis=1)
        
        return
    
    def _generate_activation_levels(self):
        '''
        Generate the level of activation normalizing by columns and rows. It includes also the patterns.
        '''
        
        # Generate uniform random levels of activation
        self.activation_levels = self.ran_generator.uniform(size=(len(self.pattern_length),self.number_of_fibers))
        
        # Generate the accumulated levels by rows and columns
        total_columns = self.number_of_fibers*0.5
        total_rows = len(self.pattern_length)*0.5
        total_pattern = total_columns*self.rate_of_fibers_in_pattern
        
        # Create a boolean matrix indicating whether a level of activation is in a pattern
        self.is_in_pattern = numpy.tile(False,(len(self.pattern_length),self.number_of_fibers))
        
        # Create a matrix indicating the length of activation for each pattern
        self.activation_length = numpy.tile(self.pattern_length.reshape((len(self.pattern_length),1)),(1,self.number_of_fibers))
        
        self.pattern_activation = numpy.empty((self.number_of_patterns,self.number_of_selected_fibers))
        
        # Normalize and copy the first realization of every pattern (except 0)
        for index in range(self.number_of_patterns):
            pattern_values = self.activation_levels[self.pattern_id_index[index+1][0],self.fibers_in_pattern[index]]
            norm_pattern_values = pattern_values/sum(pattern_values)*total_pattern
            
            # Store a copy of the pattern
            self.pattern_activation[index] = norm_pattern_values
            
            num_fibers = len(self.fibers_in_pattern[index])
            num_time = len(self.pattern_id_index[index+1])
            num_elements = num_time*num_fibers
            
            pattern_indexes_value = numpy.tile(self.pattern_id_index[index+1],(num_fibers,1)).reshape(num_elements,order='F')
            fibers_indexes_value = numpy.tile(self.fibers_in_pattern[index],(num_time,1)).reshape(num_elements)
            self.activation_levels[pattern_indexes_value,fibers_indexes_value] = numpy.tile(norm_pattern_values,(num_time,1)).reshape(num_elements)
            self.is_in_pattern[pattern_indexes_value, fibers_indexes_value] = True
            
        
        # Calculate the sum of activation levels of pattern elements
        total_pattern_per_row = numpy.sum(self.activation_levels*self.is_in_pattern, axis=0)
        total_non_pattern_per_row = total_rows - total_pattern_per_row
        total_pattern_per_column = numpy.sum(self.activation_levels*self.is_in_pattern, axis=1)
        total_non_pattern_per_column = total_columns - total_pattern_per_column
        
        repeated_total_column = numpy.tile(total_non_pattern_per_column.reshape((len(self.pattern_length),1)), (1,self.number_of_fibers))
        repeated_total_row = numpy.tile(total_non_pattern_per_row, (len(self.pattern_length),1))
            
        
        self.not_is_in_pattern = ~self.is_in_pattern
        
        # Loop for the number of normalizations
        for it in range(self.number_of_normalizations):
            # Normalize per rows
            values = self.activation_levels*self.not_is_in_pattern
            repeated_sum = numpy.tile(numpy.sum(values,axis=0),(len(self.pattern_length),1))
            norm_values = values/repeated_sum*repeated_total_row
            self.activation_levels[self.not_is_in_pattern] = norm_values[self.not_is_in_pattern]
            
            # Saturate those values above 1
            self.activation_levels[self.activation_levels>1.0] = 1.0

            # Normalize per columns
            values = self.activation_levels*self.not_is_in_pattern
            repeated_sum = numpy.tile(numpy.sum(values,axis=1).reshape((len(self.pattern_length),1)),(1,self.number_of_fibers))
            norm_values = values/repeated_sum*repeated_total_column
            self.activation_levels[self.not_is_in_pattern] = norm_values[self.not_is_in_pattern]
            
            # Saturate those values above 1
            self.activation_levels[self.activation_levels>1.0] = 1.0
            
        # Final normalization to set values between 0 and 1
        min_value = numpy.min(self.activation_levels)
        max_value = numpy.max(self.activation_levels)
        self.activation_levels = (self.activation_levels-min_value)/(max_value-min_value)
        
        # Print the sum per rows and per columns
        logger.debug('Normalization on iteration %s',it)
        logger.debug('Average sum per columns %s',numpy.average(numpy.sum(self.activation_levels,axis=1)))
        logger.debug('Average sum per rows %s',numpy.average(numpy.sum(self.activation_levels,axis=0)))
        
        
#         pylab.figure()
#         pylab.hist(self.activation_levels[self.is_in_pattern].ravel(),50)
#         pylab.figure()
#         pylab.hist(self.activation_levels.ravel(),50)
#         pylab.show()
            
        return
        
    def initialize(self):
        '''
        Initialize the pattern generator by iteratively normalizing the firing rates until it gets the desired values.
        '''
        
        # Initialize random number generator
        self.ran_generator = numpy.random.RandomState(self.seed)
        
        # Generate the length of each interval
        self._generate_length()
        
        # Generate pattern ids
        self._generate_pattern_id()
        
        # Choose the fibers that will be part of each pattern
        self._generate_fibers_in_pattern()
        
        # Generate the normalized activation levels
        self._generate_activation_levels()
        
        return
    
    def get_patterns(self):
        '''
        Generate the sequence of stimulation currents.
        '''
        for length, currents in zip(self.pattern_length, self.activation_levels):
            range_current = currents*(self.max_amplitude-self.min_amplitude) + self.min_amplitude
            yield length, range_current
    
    def get_all_patterns(self):
        '''
        Retrieve the whole sequence of stimulation currents.
        '''
        return self.pattern_length,self.activation_levels*(self.max_amplitude-self.min_amplitude)+self.min_amplitude        
