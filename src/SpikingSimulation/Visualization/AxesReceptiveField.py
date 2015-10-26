import numpy
import bisect
import AxesPlot
import logging
from mpi4py import MPI
from numpy import source

logger = logging.getLogger('Simulation')

class AxesReceptiveField(AxesPlot.AxesPlot):
    '''
    This class implements a figure showing the level of activation in the inputs when some neurons are firing.
    '''
    
    stimulation_layer = 'mflayer'
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates an object linking the axes with the DataProvider.
        @param data_provider The dataProvider that will be used in the axes to get the data. Obligatory parameter.
        @param pattern_provider: Pattern generator where we can retrieve the pattern times and cells. Optional parameter.
        @param synaptic_layer: Name of the layer to plot. Obligatory parameter.
        @param target_indexes Indexes of the target cells of the synapses to get the activity.
        @param y_window_lim: Window limits in the Y-axis. Optional parameter.
        @param load_new_data: True, it loads only the new data from the previus update call. Default=True.
        @param visible_data_only: True, it loads only the data that matches within the visualization window. Default=True
        '''
        super(AxesReceptiveField, self).__init__(**kwargs)
                
        # Get data_provider parameter 
        if ('data_provider' in kwargs):
            self.data_provider = kwargs.pop('data_provider',None)
        else:
            logger.error('Obligatory data_provider parameter not provided')
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
            logger.error('Obligatory synaptic layer parameter not provided')
            raise Exception('NonProvidedParameter','layer')
        
        # Get target cell index parameter
        if ('target_indexes' in kwargs):
            self.target_indexes = kwargs.pop('target_indexes',None)
        else:
            self.target_indexes = None
            
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
        self.figure_title = 'Receptive Field - ' + self.layer
        self.figure_x_label = 'Time'
        self.figure_y_label = 'Activation'
        
        (self.input_fibers, self.output_fibers) = self.data_provider.get_synaptic_connections(synaptic_layer = self.layer, target_indexes = self.target_indexes)
        
        # Set axes lines and legends
        data_labels = ['Cell ' + str(ind) for ind in self.target_indexes]
        number_of_lines = len(self.target_indexes)
    
        self.axesLines = []
        
        for _ in range(number_of_lines):
            #newLine = self.axes.scatter([],[],marker='.')
            newLine, = self.axes.plot([], [], linestyle='', marker='.', markersize=3)
            self.axesLines.append(newLine)
        
        if (self.show_legend):
            self.axes.legend(self.axesLines,data_labels,loc='lower left')
            
        if (self.y_window_lim is not None):
            self.axes.set_ylim(self.y_window_lim)
        else:
            self.axes.set_ylim([-0.01,1.01])  
                
        
        if (self.y_window_lim is not None):
            self.axes.set_ylim(self.y_lim)
        else:
            self.axes.set_ylim([0,1])  
        
        super(AxesReceptiveField, self).initialize()
            
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
        gtime,gcell_id = self.data_provider.get_spike_activity(neuron_layer = self.data_provider.layer_map[self.layer].target_layer.__name__, neuron_indexes = self.target_indexes,\
                                                               init_time = load_data_init, end_time = simulation_time)
        
        self.data_update = simulation_time
        
        comm = MPI.COMM_WORLD
        
        process_id = comm.Get_rank()
        
        if (process_id==0):
        
            self.axes.set_xlim(time_window) # Set x_axes limits axes
        
            for ind,cell in enumerate(self.target_indexes):
                initialized = False
            
                new_time = gtime[(gcell_id == cell)]
    
                # Calculate each spike's offset
                if new_time.size:
                    sel_output_cells = (self.output_fibers==cell)
                    num_cells = numpy.count_nonzero(sel_output_cells)
                    bin_index = [bisect.bisect_right(self.pattern_provider.pattern_length_cum, time) for time in new_time]
                    sel_cells = numpy.repeat(self.input_fibers[sel_output_cells],len(bin_index))
                    sel_time = numpy.tile(new_time,num_cells)
                    sel_bin = numpy.tile(bin_index,num_cells)
                    sel_activation = self.pattern_provider.activation_levels[sel_bin,sel_cells]
                else:
                    sel_time = []
                    sel_activation = []
            
                old_time_data = self.axesLines[ind].get_xdata()
                first_index = numpy.searchsorted(old_time_data, data_init_time)
                new_time_data = numpy.append(old_time_data[first_index:],sel_time)
            
                old_signal_data = self.axesLines[ind].get_ydata()
                new_signal_data = old_signal_data[first_index:]
                new_signal_data = numpy.append(new_signal_data, sel_activation)
                
                
                self.axesLines[ind].set_xdata(new_time_data)
                self.axesLines[ind].set_ydata(new_signal_data)
                    
                
    
        return self.axesLines