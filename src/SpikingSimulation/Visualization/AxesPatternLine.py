import numpy
import bisect
import AxesPlot
import logging
import math
from mpi4py import MPI

logger = logging.getLogger('Simulation')

class AxesPatternLine(AxesPlot.AxesPlot):
    '''
    This class defines the link between an axes object and the dataProvider for plots.
    '''
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates an object linking the axes with the DataProvider.
        @param pattern_provider: Pattern generator where we can retrieve the pattern times and cells. Optional parameter.
        @param pattern: List of patterns to highlight. 0 is noise. Optional parameter. If nothing is specified, all the patterns will be shown.
        @param show_legend: True if the legend will be shown. Optional parameter. Default: True
        @param visible_data_only: True, it loads only the data that matches within the visualization window. Default=True
        @param load_new_data: True, it loads only the new data from the previus update call. Default=True.
        @param x_length: X-axis window limits in case of trial_per_trial figures. Optinal parameter. If x_lenght is not specified, the trial length will be used.
        '''
        super(AxesPatternLine, self).__init__(**kwargs)
                
        # Get pattern_provider parameter 
        if ('pattern_provider' in kwargs):
            self.pattern_provider = kwargs.pop('pattern_provider',None)
        else:
            logger.error('Obligatory pattern_provider parameter not provided')
            raise Exception('NonProvidedParameter','pattern_provider')
        
        # Get number of patterns to highlight 
        if ('pattern' in kwargs):
            self.pattern = kwargs.pop('pattern',None)
        else:
            self.pattern = None
        
        # Get visible_data_only parameter 
        if ('visible_data_only' in kwargs):
            self.visible_data_only = kwargs.pop('visible_data_only',None)
        else:
            self.visible_data_only = True
            
        # Get load_new_data parameter 
        if ('load_new_data' in kwargs):
            self.load_new_data = kwargs.pop('load_new_data',None)
        else:
            self.load_new_data = True
            
        if ('show_legend' in kwargs):
            self.show_legend = kwargs.pop('show_legend', None)
        else:
            self.show_legend = True
            
        # Get axesType parameter
        if ('x_length' in kwargs):
            self.x_length = kwargs.pop('x_length',None)
        else:
            self.x_length = None
            
        
        # Time of the last update
        self.data_update = 0
        
                
                
    def initialize(self):
        '''
        Perform all the required operations needed in order to initialize the axes.
        It sets the title, axes titles, axes limits, legend and creates the lines.
        The DataProvider object must be initialized before calling this function.
        '''
        self.figure_title = 'Pattern selector'
        self.figure_x_label = 'Time (s)'
        self.figure_y_label = ''
         
        # Set axes lines and legends
        if self.pattern_provider:
            if not self.pattern:
                self.pattern = range(self.pattern_provider.number_of_patterns)
        else:
            self.pattern = []
        
        self.time_bin = 1.e-3
        self.inv_time_bin = 1./self.time_bin
        
        # Initialize the pattern generator information
        # Generate the time bin matrix
        self.total_time = self.pattern_provider.simulation_time
        self.bin_time_init = numpy.linspace(0.0, self.total_time-self.time_bin, num=self.total_time*self.inv_time_bin)
        
        # Initialize a matrix
        self.num_patterns = len(self.pattern)
        self.num_bins = len(self.bin_time_init)
        
        # Calculate the time of each pattern interval
        time_end_of_pattern = self.pattern_provider.pattern_length_cum
        
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
        
        for value in self.pattern:
            for index in self.pattern_provider.pattern_id_index[value]:
                init_bin = bin_init_of_pattern[index]
                end_bin = bin_end_of_pattern[index]
                
                list_of_bins = range(init_bin,end_bin)
                
                # Add the time of the initial bin (if exists)
                self.bin_is_pattern[value,list_of_bins] = True
        
        # Initialize the axes
        data_labels = []
        for pat in self.pattern:
            data_labels.append('Pattern '+str(pat+1))
            
        number_of_lines = len(data_labels)
        
        self.axesLines = []
        
        for _ in range(number_of_lines):
            #newLine = self.axes.scatter([],[],marker='.')
            newLine, = self.axes.plot([], [])
            self.axesLines.append(newLine)
        
        if (self.show_legend):
            self.axes.legend(self.axesLines,data_labels,loc='lower left')
            
        self.axes.set_ylim([-0.05,self.pattern_provider.number_of_patterns + 0.05])  
            
        super(AxesPatternLine, self).initialize()
            
        return
    
   
    def getTimeWindow(self, simulation_time):
        if self.x_length:
            minTime = max(simulation_time-self.x_length,0)
        else:
            minTime = 0
            
        maxTime = max(self.x_length,simulation_time)
        return [minTime,maxTime]
        
        
    def drawAtTime(self, simulation_time):
        '''
        This function updates all the elements of the axes according to the existent data
        until time simulationTime.
        @param simulation_time: The simulation end time (in seconds).
        @return A list with the artist to be updated. 
        '''
        
        time_window = self.getTimeWindow(simulation_time = simulation_time)
        
        # Check if load all the data or only the data within the visualization window
        if (self.visible_data_only):
            data_init_time = time_window[0]
        else:
            data_init_time = self.simulation_limits[0]
            
        # Check if load only new data or all the data
        if (self.load_new_data):
            load_data_init = max(data_init_time, self.data_update)
        else:
            load_data_init = data_init_time
        
        # Load data from the data provider
        init_bin = int(math.floor(load_data_init*self.inv_time_bin))
        end_bin = int(math.floor(simulation_time*self.inv_time_bin))
        
        self.data_update = simulation_time
        
        comm = MPI.COMM_WORLD
        
        process_id = comm.Get_rank()
        
        if (process_id==0):
        
            self.axes.set_xlim(time_window) # Set x_axes limits axes
        
            y_limits = range(2)
            initialized = False
            
            sel_time = self.bin_time_init[init_bin:end_bin]
            
            for index,pattern_idx in enumerate(self.pattern):
                sel_data = self.bin_is_pattern[index,init_bin:end_bin] * (pattern_idx + 1.0)
                
                old_time_data = self.axesLines[index].get_xdata()
                first_index = numpy.searchsorted(old_time_data, data_init_time)
                new_time_data = numpy.append(old_time_data[first_index:],sel_time)
            
                old_signal_data = self.axesLines[index].get_ydata()
                new_signal_data = old_signal_data[first_index:]
                new_signal_data = numpy.append(new_signal_data, sel_data)
                
                self.axesLines[index].set_xdata(new_time_data)
                self.axesLines[index].set_ydata(new_signal_data)
                
#                 if (new_signal_data.size):
#                     if (initialized):
#                         y_limits[0] = min(y_limits[0],numpy.min(new_signal_data))
#                         y_limits[1] = max(y_limits[1],numpy.max(new_signal_data))
#                     else:
#                         initialized = True
#                         y_limits[0] = numpy.min(new_signal_data)
#                         y_limits[1] = numpy.max(new_signal_data)
            
        return self.axesLines