import numpy
import bisect
import AxesPlot
from mpi4py import MPI

class AxesRasterPlot(AxesPlot.AxesPlot):
    '''
    This class defines the link between an axes object and the dataProvider for plots.
    '''
    
    stimulation_layer = 'mflayer'
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates an object linking the axes with the DataProvider.
        @param data_provider The dataProvider that will be used in the axes to get the data. Obligatory parameter.
        @param pattern_provider: Pattern generator where we can retrieve the pattern times and cells. Optional parameter.
        @param layer: Name of the layer to plot. Obligatory parameter.
        @param cell_index: Indexes of the cells to plot. Optional parameter.
        @param pattern: List of patterns to highlight. 0 is noise. Optional parameter. If nothing is specified, all the patterns will be shown.
        @param show_legend: True if the legend will be shown. Optional parameter. Default: True
        @param y_window_lim: Window limits in the Y-axis. Optional parameter.
        @param visible_data_only: True, it loads only the data that matches within the visualization window. Default=True
        @param load_new_data: True, it loads only the new data from the previus update call. Default=True.
        @param x_length: X-axis window limits in case of trial_per_trial figures. Optinal parameter. If x_lenght is not specified, the trial length will be used.
        '''
        super(AxesRasterPlot, self).__init__(**kwargs)
                
        # Get data_provider parameter 
        if ('data_provider' in kwargs):
            self.data_provider = kwargs.pop('data_provider',None)
        else:
            print 'Obligatory data_provider parameter not provided'
            raise Exception('NonProvidedParameter','data_provider')
        
        # Get data_provider parameter 
        if ('pattern_provider' in kwargs):
            self.pattern_provider = kwargs.pop('pattern_provider',None)
        else:
            self.pattern_provider = None
        
        # Get layer name parameter
        if ('layer' in kwargs):
            self.layer = kwargs.pop('layer',None)
        else:
            print 'Obligatory layer parameter not provided'
            raise Exception('NonProvidedParameter','layer')
        
        # Get index parameter
        if ('cell_index' in kwargs):
            self.index = kwargs.pop('cell_index',None)
        else:
            self.index = range(self.data_provider.get_number_of_elements(layer=self.layer))
            
        # Get number of patterns to highlight 
        if ('pattern' in kwargs):
            self.pattern = kwargs.pop('pattern',None)
        else:
            self.pattern = None
        
        # Get y_window_lim parameter 
        if ('y_window_lim' in kwargs):
            self.y_window_lim = kwargs.pop('y_window_lim',None)
        else:
            self.y_window_lim = None
            
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
        self.figure_title = 'Raster plot - ' + self.layer
        self.figure_x_label = 'Time (s)'
        self.figure_y_label = 'Cell id'
         
        # Set axes lines and legends
        if self.pattern_provider:
            if not self.pattern:
                self.pattern = range(self.pattern_provider.number_of_patterns+1)
        else:
            self.pattern = [0]
        
        data_labels = []
        for pat in self.pattern:
            if pat:
                data_labels.append('Pattern '+str(pat))
            else:
                data_labels.append('Noise')
        
        number_of_lines = len(data_labels)
        
        self.axesLines = []
        
        for _ in range(number_of_lines):
            #newLine = self.axes.scatter([],[],marker='.')
            newLine, = self.axes.plot([], [], linestyle='', marker='.', markersize=8, markeredgecolor = 'none')
            self.axesLines.append(newLine)
        
        if (self.show_legend):
            self.axes.legend(self.axesLines,data_labels,loc='lower right')
            
        super(AxesRasterPlot, self).initialize()
            
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
        gtime,gcell_id = self.data_provider.get_spike_activity(neuron_layer = self.layer, neuron_indexes = self.index,\
                                                               init_time = load_data_init, end_time = simulation_time)
        
        self.data_update = simulation_time
        
        comm = MPI.COMM_WORLD
        
        process_id = comm.Get_rank()
        
        if (process_id==0):
        
            self.axes.set_xlim(time_window) # Set x_axes limits axes
        
            y_limits = range(2)
            initialized = False
            
            
            for i in range(len(self.pattern)):
                pat = self.pattern[i]
            
                # Filter the activity for each pattern
                if self.pattern_provider:
                
                    # Find the pattern of the time bin for each spike
                    pattern_of_spike = numpy.array([self.pattern_provider.pattern_id[bisect.bisect_right(self.pattern_provider.pattern_length_cum, time)] for time in gtime])
                
                    if self.layer == self.stimulation_layer:
                        if pat:
                            # Select those spikes fired for the cells in the selected fibers of each pattern
                            indexes = (pattern_of_spike==pat) & numpy.in1d(gcell_id,self.pattern_provider.fibers_in_pattern[pat-1])
                        else:
                            # Select those spikes fired in noise interval or those in pattern bin but are not included in the selected fibers
                            indexes = (pattern_of_spike==pat)
                            for pat1 in range(1,self.pattern_provider.number_of_patterns+1):
                                indexes = indexes | ((pattern_of_spike==pat1) & ~numpy.in1d(gcell_id,self.pattern_provider.fibers_in_pattern[pat1-1]))
                    else:
                        indexes = (pattern_of_spike==pat)
                        
                    sel_time = gtime[indexes]
                    sel_cell = gcell_id[indexes]
                else:
                    sel_time = gtime
                    sel_cell = gcell_id
                    
                old_time_data = self.axesLines[i].get_xdata()
                first_index = numpy.searchsorted(old_time_data, data_init_time)
                new_time_data = numpy.append(old_time_data[first_index:],sel_time)
            
                old_signal_data = self.axesLines[i].get_ydata()
                new_signal_data = old_signal_data[first_index:]
                new_signal_data = numpy.append(new_signal_data, sel_cell)
                
                self.axesLines[i].set_xdata(new_time_data)
                self.axesLines[i].set_ydata(new_signal_data)
                
                if (new_signal_data.size):
                    if (initialized):
                        y_limits[0] = min(y_limits[0],numpy.min(new_signal_data))
                        y_limits[1] = max(y_limits[1],numpy.max(new_signal_data))
                    else:
                        initialized = True
                        y_limits[0] = numpy.min(new_signal_data)
                        y_limits[1] = numpy.max(new_signal_data)
                
            if (self.y_window_lim is not None):
                self.axes.set_ylim(self.y_lim)
            elif initialized:
                self.axes.set_ylim(y_limits)
            else:
                self.axes.set_ylim([0,1])  
    
        return self.axesLines