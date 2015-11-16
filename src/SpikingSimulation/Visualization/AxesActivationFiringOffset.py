import numpy
import bisect
import AxesPlot
import logging
from mpi4py import MPI

logger = logging.getLogger('Simulation')

class AxesActivationFiringOffset(AxesPlot.AxesPlot):
    '''
    This class defines the link between an axes object and the dataProvider for plots.
    '''
    
    stimulation_layer = 'mflayer'
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates an object linking the axes with the DataProvider.
        @param data_provider The dataProvider that will be used in the axes to get the data. Obligatory parameter.
        @param oscillation_freq The frequency of the oscillation, determining the init of each cycle. Obligatory parameter.
        @param pattern_provider: Pattern generator where we can retrieve the pattern times and cells. Optional parameter.
        @param layer: Name of the layer to plot. Obligatory parameter.
        @param cell_index: Indexes of the cells to plot. Optional parameter.
        @param y_window_lim: Window limits in the Y-axis. Optional parameter.
        @param visible_data_only: True, it loads only the data that matches within the visualization window. Default=True
        @param load_new_data: True, it loads only the new data from the previus update call. Default=True.
        '''
        super(AxesActivationFiringOffset, self).__init__(**kwargs)
                
        # Get data_provider parameter 
        if ('data_provider' in kwargs):
            self.data_provider = kwargs.pop('data_provider',None)
        else:
            logger.error('Obligatory data_provider parameter not provided')
            raise Exception('NonProvidedParameter','data_provider')
        
        # Get oscillation_freq parameter 
        if ('oscillation_freq' in kwargs):
            self.oscillation_freq = kwargs.pop('oscillation_freq',None)
        else:
            logger.error('Obligatory oscillation_freq parameter not provided')
            raise Exception('NonProvidedParameter','oscillation_freq')
        
        # Get data_provider parameter 
        if ('pattern_provider' in kwargs):
            self.pattern_provider = kwargs.pop('pattern_provider',None)
        else:
            self.pattern_provider = None
        
        # Get layer name parameter
        if ('layer' in kwargs):
            self.layer = kwargs.pop('layer',None)
        else:
            logger.error('Obligatory layer parameter not provided')
            raise Exception('NonProvidedParameter','layer')
        
        # Get index parameter
        if ('cell_index' in kwargs):
            self.index = kwargs.pop('cell_index',None)
        else:
            self.index = range(self.data_provider.get_number_of_elements(layer=self.layer))
            
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
            
        # Time of the last update
        self.data_update = 0
        
                
                
    def initialize(self):
        '''
        Perform all the required operations needed in order to initialize the axes.
        It sets the title, axes titles, axes limits, legend and creates the lines.
        The DataProvider object must be initialized before calling this function.
        '''
        self.figure_title = 'Activation vs Firing Offset plot - ' + self.layer
        self.figure_x_label = 'Activation'
        self.figure_y_label = 'Offset (s)'
        
        self.oscillation_period = 1./self.oscillation_freq
        
        number_of_lines = 1
        
        self.axesLines = []
        animated_artists = []
        
        for _ in range(number_of_lines):
            #newLine = self.axes.scatter([],[],marker='.')
            if (self.figure.blit):
                newLine, = self.axes.plot([], [], linestyle='', marker='.', markersize=5, animated=True)
                animated_artists.append(newLine)
            else:
                newLine, = self.axes.plot([], [], linestyle='', marker='.', markersize=5)
                
            self.axesLines.append(newLine)
        
        if (self.y_window_lim is not None):
            self.axes.set_ylim(self.y_lim)
        else:
            self.axes.set_ylim([0,self.oscillation_period])  
        
        self.axes.set_xlim([0,1])    
        
        self.animated_artists = tuple(animated_artists)
        
        super(AxesActivationFiringOffset, self).initialize()
            
        return
    
   
    def getTimeWindow(self, simulation_time):
        minTime = max(simulation_time-self.oscillation_period,0)
        maxTime = max(self.oscillation_period,simulation_time)
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
        
            initialized = False

            # Calculate each spike's offset
            if gtime.size:
                sel_offset = numpy.mod(gtime,self.oscillation_period)
                bin_index = numpy.array([bisect.bisect_right(self.pattern_provider.pattern_length_cum, time) for time in gtime])
                sel_activation = self.pattern_provider.activation_levels[bin_index,gcell_id]
            else:
                sel_offset = []
                sel_activation = []
        
            if self.visible_data_only:
                new_time_data = sel_activation
                new_signal_data = sel_offset
            else:
                old_time_data = self.axesLines[0].get_xdata()
                new_time_data = numpy.append(old_time_data,sel_activation)
                old_signal_data = self.axesLines[0].get_ydata()
                new_signal_data = numpy.append(old_signal_data, sel_offset)
                
            self.axesLines[0].set_xdata(new_time_data)
            self.axesLines[0].set_ydata(new_signal_data)
                
                
    
        return self.animated_artists