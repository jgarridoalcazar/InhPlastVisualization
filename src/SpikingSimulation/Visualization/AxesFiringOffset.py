import numpy
import bisect
import AxesPlot
import logging
from mpi4py import MPI

logger = logging.getLogger('Simulation')

class AxesFiringOffset(AxesPlot.AxesPlot):
    '''
    This class defines the link between an axes object and the dataProvider for plots.
    '''
    
    stimulation_layer = 'mflayer'
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates an object linking the axes with the DataProvider.
        @param data_provider The dataProvider that will be used in the axes to get the data. Obligatory parameter.
        @param oscillation_freq The frequency of the oscillation, determining the init of each cycle. Obligatory parameter.
        @param layer: Name of the layer to plot. Obligatory parameter.
        @param cell_index: Indexes of the cells to plot. Optional parameter.
        @param y_window_lim: Window limits in the Y-axis. Optional parameter.
        @param visible_data_only: True, it loads only the data that matches within the visualization window. Default=True
        @param load_new_data: True, it loads only the new data from the previus update call. Default=True.
        @param x_length: X-axis window limits in case of trial_per_trial figures. Optinal parameter. If x_lenght is not specified, the trial length will be used.
        '''
        super(AxesFiringOffset, self).__init__(**kwargs)
                
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
        self.figure_title = 'Firing Offset plot - ' + self.layer
        self.figure_x_label = 'Time (s)'
        self.figure_y_label = 'Offset (s)'
        
        self.oscillation_period = 1./self.oscillation_freq
         
        data_labels = []
        for ind in self.index:
            data_labels.append('Cell '+str(ind))
        
        number_of_lines = len(data_labels) 
        
        self.axesLines = []
        
        for _ in range(number_of_lines):
            #newLine = self.axes.scatter([],[],marker='.')
            newLine, = self.axes.plot([], [], linestyle='', marker='.', markersize=5)
            self.axesLines.append(newLine)
        
        if (self.y_window_lim is not None):
            self.axes.set_ylim(self.y_lim)
        else:
            self.axes.set_ylim([0,self.oscillation_period])  
            
        super(AxesFiringOffset, self).initialize()
            
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
        
            initialized = False

            for ind,cell in enumerate(self.index):
                sel_time = gtime[gcell_id==cell]
                sel_offset = numpy.mod(sel_time,self.oscillation_period)
            
                old_time_data = self.axesLines[ind].get_xdata()
                first_index = numpy.searchsorted(old_time_data, data_init_time)
                new_time_data = numpy.append(old_time_data[first_index:],sel_time)
            
                old_signal_data = self.axesLines[ind].get_ydata()
                new_signal_data = old_signal_data[first_index:]
                new_signal_data = numpy.append(new_signal_data, sel_offset)
                
                self.axesLines[ind].set_xdata(new_time_data)
                self.axesLines[ind].set_ydata(new_signal_data)
                
                
    
        return self.axesLines