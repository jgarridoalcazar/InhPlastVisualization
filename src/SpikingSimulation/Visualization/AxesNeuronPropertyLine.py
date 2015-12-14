import numpy
import AxesPlot
import logging

logger = logging.getLogger('Simulation')

class AxesNeuronPropertyLine(AxesPlot.AxesPlot):
    '''
    This class defines the link between an axes object and the dataProvider for plots.
    '''
    # Figure title (Figure title, X-axis title, Y-axis title)
    figure_property_map = {'Vm': ['Membrane Potential','Vm (V)'],\
                   'Gexc': ['Excitatory Conductance', 'Gexc (S)'],\
                   'Ginh': ['Inhibitory Conductance','Ginh (S)'],\
                   'gL': ['Leak Conductance', 'Gleak (S)'],\
                   'rC': ['Inverse Membrane Capacitance', 'rC (pF-1)'],\
                   'Vth': ['Threshold Potential', 'Vth (V)'],\
                   'r0': ['Gain Constant', 'Gain (Hz)'],\
                   'ualpha': ['Potential Scale', 'u_alpha (V)'],\
                   'refractoriness': ['Refractoriness', 'Refr.'],\
                   'gain': ['Instantaneous Frequency', 'Freq. (Hz)'],\
                   'firing_probability': ['Firing Probability', 'FProp']
                   }
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates an object linking the axes with the DataProvider.
        @param data_provider The dataProvider that will be used in the axes to get the data. Obligatory parameter.
        @param property: Name of the property to plot. Obligatory parameter.
        @param layer: Name of the layer to plot. Obligatory parameter.
        @param cell_index: Indexes of the cells to plot. Optional parameter.
        @param show_legend: True if the legend will be shown. Optional parameter. Default: True
        @param y_window_lim: Window limits in the Y-axis. Optional parameter.
        @param visible_data_only: True, it loads only the data that matches within the visualization window. Default=True
        @param load_new_data: True, it loads only the new data from the previus update call. Default=True.
        @param x_length: X-axis window limits in case of trial_per_trial figures. Optinal parameter. If x_lenght is not specified, the trial length will be used.
        '''
        super(AxesNeuronPropertyLine, self).__init__(**kwargs)
                
        # Get data_provider parameter 
        if ('data_provider' in kwargs):
            self.data_provider = kwargs.pop('data_provider',None)
        else:
            logger.error('Obligatory data_provider parameter not provided')
            raise Exception('NonProvidedParameter','data_provider')
        
        # Get property name parameter 
        if ('property' in kwargs):
            self.property = kwargs.pop('property',None)
        else:
            logger.error('Obligatory property parameter not provided')
            raise Exception('NonProvidedParameter','property')
        
        if (self.property not in self.figure_property_map):
            logger.error('%s is not implemented as a type %s figure', self.property, type(self))
            raise Exception('FigureNotImplemented')
        
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
        self.figure_title = self.figure_property_map[self.property][0] + ' - ' + self.layer
        self.figure_x_label = 'Time (s)'
        self.figure_y_label = self.figure_property_map[self.property][1]
         
        # Set axes lines and legends
        data_labels = ['Cell ' + str(ind) for ind in self.index]
        number_of_lines = len(self.index)
        
        self.axesLines = []
        animated_artists = []
        
        for _ in range(number_of_lines):
            if (self.figure.blit):
                newLine, = self.axes.plot([], [], animated=True)
                animated_artists.append(newLine)
            else:
                newLine, = self.axes.plot([], [])
            self.axesLines.append(newLine)
        
        if (self.show_legend):
            self.axes.legend(self.axesLines,data_labels,loc='lower left')
         
        self.visual_time_window = self.getVisualTimeWindow(0)
        self.time_window = [0,0]
        
        self.axes.set_xlim(self.time_window)
        
        if (self.figure.blit and not self.x_length):
            self.axes.xaxis.set_animated(True)
            animated_artists.append(self.axes.xaxis)
            
        if (self.y_window_lim is not None):
            self.axes.set_ylim(self.y_window_lim)
        else:
            self.axes.set_ylim([0,1])
            if (self.figure.blit):
                self.axes.yaxis.set_animated(True)
                animated_artists.append(self.axes.yaxis)  
        
        self.animated_artists = tuple(animated_artists)
           
        super(AxesNeuronPropertyLine, self).initialize()
            
        return
    
   
    def getVisualTimeWindow(self, simulation_time):
        if self.x_length:
            minTime = -self.x_length
            maxTime = 0
        else:
            minTime = 0
            maxTime = simulation_time
            
        return [minTime,maxTime]
    
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
        
        old_visual_time_window = self.visual_time_window
        self.visual_time_window = self.getVisualTimeWindow(simulation_time = simulation_time)
        old_time_window = self.time_window
        self.time_window = self.getTimeWindow(simulation_time = simulation_time)
        
        # Check if load all the data or only the data within the visualization window
        if (self.visible_data_only):
            data_init_time = self.time_window[0]
        else:
            data_init_time = self.simulation_limits[0]
            
        # Check if load only new data or all the data
        if (self.load_new_data):
            load_data_init = max(data_init_time, self.data_update)
        else:
            load_data_init = data_init_time
        
        # Load data from the data provider
        gtime,gcell_id,gvalue = self.data_provider.get_state_variable(state_variable = self.property, neuron_layer = self.layer,\
                                                              neuron_indexes = self.index, init_time = load_data_init, end_time = simulation_time)
        
        self.data_update = simulation_time
        
        from mpi4py import MPI
        
        comm = MPI.COMM_WORLD
        
        process_id = comm.Get_rank()
        
        if (process_id==0):
        
            if (self.visual_time_window[0]!=old_time_window[0] or 
                self.visual_time_window[1]!=old_time_window[1]):
                self.axes.set_xlim(self.visual_time_window) # Set x_axes limits axes
        
            y_limits = range(2)
            initialized = False
        
            # Select the values for each cell_id
            for ind in range(len(self.axesLines)):
                cell_index = (gcell_id==self.index[ind])
                time = gtime[cell_index]
                value = gvalue[cell_index]
                
                old_time_data = self.axesLines[ind].get_xdata()
                if (self.x_length is None):
                    first_index = numpy.searchsorted(old_time_data, data_init_time)
                    new_time_data = numpy.append(old_time_data[first_index:],time)
                else:
                    first_index = numpy.searchsorted(old_time_data+old_time_window[1], data_init_time)
                    new_time_data = numpy.append(old_time_data[first_index:],time-self.time_window[1])
                
                old_signal_data = self.axesLines[ind].get_ydata()
                new_signal_data = old_signal_data[first_index:]
                new_signal_data = numpy.append(new_signal_data, value)
                
                self.axesLines[ind].set_xdata(new_time_data)
                self.axesLines[ind].set_ydata(new_signal_data)
                
                if (new_signal_data.size and self.y_window_lim is None):
                    if (initialized):
                        y_limits[0] = min(y_limits[0],numpy.min(new_signal_data))
                        y_limits[1] = max(y_limits[1],numpy.max(new_signal_data))
                    else:
                        initialized = True
                        y_limits[0] = numpy.min(new_signal_data)
                        y_limits[1] = numpy.max(new_signal_data)
                
            if (self.y_window_lim is None):
                self.axes.set_ylim(y_limits)  
    
        return self.animated_artists