import numpy
import AxesPlot
import logging
from mpi4py import MPI

logger = logging.getLogger('Simulation')

class AxesWeightEvolutionLine(AxesPlot.AxesPlot):
    '''
    This class defines the link between an axes object and the dataProvider for plots.
    '''       
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates an object linking the axes with the DataProvider.
        @param data_provider The dataProvider that will be used in the axes to get the data. Obligatory parameter.
        @param layer: Name of the layer to plot. Obligatory parameter.
        @param source_indexes Indexes of the source cells of the synapses to get the activity.
        @param target_indexes Indexes of the target cells of the synapses to get the activity.
        @param show_legend: True if the legend will be shown. Optional parameter. Default: True
        @param y_window_lim: Window limits in the Y-axis. Optional parameter.
        @param visible_data_only: True, it loads only the data that matches within the visualization window. Default=True
        @param load_new_data: True, it loads only the new data from the previus update call. Default=True.
        @param x_length: X-axis window limits in case of trial_per_trial figures. Optinal parameter. If x_lenght is not specified, the trial length will be used.
        '''
        
        
        super(AxesWeightEvolutionLine, self).__init__(**kwargs)
                
        # Get data_provider parameter 
        if ('data_provider' in kwargs):
            self.data_provider = kwargs.pop('data_provider',None)
        else:
            logger.error('Obligatory data_provider parameter not provided')
            raise Exception('NonProvidedParameter','data_provider')
        
        # Get layer name parameter
        if ('layer' in kwargs):
            self.layer = kwargs.pop('layer',None)
        else:
            logger.error('Obligatory layer parameter not provided')
            raise Exception('NonProvidedParameter','layer')
        
        # Get souce cell index parameter
        if ('source_indexes' in kwargs):
            self.source_indexes = kwargs.pop('source_indexes',None)
        else:
            self.source_indexes = None
            
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
        self.figure_title = 'Weight Evolution - ' + self.layer
        self.figure_x_label = 'Time (s)'
        self.figure_y_label = 'Weight (S)'
        
        # Set connection parameters
        self.param = dict()
        
        self.param['synaptic_layer'] = self.layer
        if self.target_indexes:
            self.param['target_indexes'] = self.target_indexes
            target_cells = self.target_indexes
        else:
            target_cells = range(self.data_provider.layer_map[self.layer].target_layer.number_of_neurons)
            
        if self.source_indexes:
            self.param['source_indexes'] = self.source_indexes
            source_cells = self.source_indexes
        else:
            source_cells = range(self.data_provider.layer_map[self.layer].source_layer.number_of_neurons)
        
        self.param['init_time'] = 0
        self.param['end_time'] = 0
        
        # Load data from the data provider
        _,self.connections,_ = self.data_provider.get_synaptic_weights(**self.param)
        
        # self.connections = [[source,target] for source in source_cells for target in target_cells]
        
        # Set axes lines and legends
        data_labels = [str(syn[0]) + ' - ' + str(syn[1]) for syn in self.connections]
        number_of_lines = len(data_labels)
        
        self.axesLines = []
        
        for _ in range(number_of_lines):
            newLine, = self.axes.plot([],[])
            self.axesLines.append(newLine)
        
        if (self.show_legend):
            self.axes.legend(self.axesLines,data_labels,loc='lower right')
            
        super(AxesWeightEvolutionLine, self).initialize()
            
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
        
        self.param['init_time'] = load_data_init
        self.param['end_time'] = simulation_time
        
        # Load data from the data provider
        gtime,gconnections,gvalue = self.data_provider.get_synaptic_weights(**self.param)
        
        self.data_update = simulation_time
        
        comm = MPI.COMM_WORLD
        
        process_id = comm.Get_rank()
        
        if (process_id==0):
        
            self.axes.set_xlim(time_window) # Set x_axes limits axes
        
            y_limits = range(2)
            initialized = False
        
            # Select the values for each cell_id
            for ind in range(len(self.axesLines)):
                indexes = [index for index in range(len(self.connections)) if (len(gconnections)>index) and (gconnections[index][0]==self.connections[ind][0]) and (gconnections[index][1]==self.connections[ind][1])]
                time = gtime
                value = gvalue[indexes]
                
                old_time_data = self.axesLines[ind].get_xdata()
                first_index = numpy.searchsorted(old_time_data, data_init_time)
                new_time_data = numpy.append(old_time_data[first_index:],time)
            
                old_signal_data = self.axesLines[ind].get_ydata()
                new_signal_data = old_signal_data[first_index:]
                new_signal_data = numpy.append(new_signal_data, value)
                
                self.axesLines[ind].set_xdata(new_time_data)
                self.axesLines[ind].set_ydata(new_signal_data)
                
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