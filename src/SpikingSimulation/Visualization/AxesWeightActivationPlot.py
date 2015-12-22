import numpy
import bisect
import AxesPlot
import logging

logger = logging.getLogger('Simulation')

class AxesWeightActivationPlot(AxesPlot.AxesPlot):
    '''
    This class defines the link between an axes object and the dataProvider for plots.
    '''
    
    stimulation_layer = 'mflayer'
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates an object linking the axes with the DataProvider.
        @param data_provider The dataProvider that will be used in the axes to get the data. Obligatory parameter.
        @param pattern_provider: Pattern generator where we can retrieve the pattern times and cells. Obligatory parameter.
        @param layer: Name of the synaptic layer to plot. Obligatory parameter.
        @param source_indexes Indexes of the source cells of the synapses to get the activity.
        @param target_indexes Indexes of the target cells of the synapses to get the activity.
        @param pattern: List of patterns to highlight. From 0 to number of patterns-1. Optional parameter. If nothing is specified, all the patterns will be shown.
        @param y_window_lim: Window limits in the Y-axis. Optional parameter.
        @param show_legend: True if the legend will be shown. Optional parameter. Default: True
        '''
        super(AxesWeightActivationPlot, self).__init__(**kwargs)
                
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
            logger.error('Obligatory pattern_provider parameter not provided')
            raise Exception('NonProvidedParameter','pattern_provider')
        
        if ('show_legend' in kwargs):
            self.show_legend = kwargs.pop('show_legend', None)
        else:
            self.show_legend = True
            
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
        
        # Get number of patterns to highlight 
        if ('pattern' in kwargs):
            self.pattern = kwargs.pop('pattern',None)
        else:
            self.pattern = None
        
        if (self.data_provider.layer_map[self.layer].source_layer!=self.data_provider.layer_map[self.stimulation_layer]):
            logger.error('Invalid connection layer. %s is not connected from %s', self.layer, self.stimulation_layer)
            raise Exception('InvalidConnectionLayer','InvalidLayer')
        
                
                
    def initialize(self):
        '''
        Perform all the required operations needed in order to initialize the axes.
        It sets the title, axes titles, axes limits, legend and creates the lines.
        The DataProvider object must be initialized before calling this function.
        '''
        self.figure_title = 'Activation/Weight plot - ' + self.layer
        self.figure_x_label = 'Normalized Activation'
        self.figure_y_label = 'Weight (nS)'
         
        # Set axes lines and legends
        if not self.pattern:
            self.pattern = range(self.pattern_provider.number_of_patterns)
            
        # Set connection parameters
        self.param = dict()
        
        self.param['synaptic_layer'] = self.layer
        self.param['init_time'] = 0
        self.param['end_time'] = 0
        
        if self.target_indexes:
            self.param['target_indexes'] = self.target_indexes
            
        if self.source_indexes:
            self.param['source_indexes'] = self.source_indexes
        
        
        data_labels = ['Pattern '+str(pat+1) for pat in self.pattern]
        
        number_of_lines = len(data_labels)
        
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
        
        if (self.show_legend):
            self.axes.legend(self.axesLines,data_labels,loc='lower left')
        
        synaptic_layer = self.data_provider.layer_map[self.layer]
        
        if synaptic_layer.learning_rule_type:
            max_weight = synaptic_layer.learning_rule_parameters['max_weight']
        else:
            max_weight = 1.
                
        self.axes.set_xlim([0,1])
        
        if (self.y_window_lim is not None):
            self.axes.set_ylim(self.y_window_lim)
        else:
            self.axes.set_ylim([0,max_weight*1.10])
            
        self.animated_artists = tuple(animated_artists)
            
        super(AxesWeightActivationPlot, self).initialize()
        
        return
    
   
    def drawAtTime(self, simulation_time):
        '''
        This function updates all the elements of the axes according to the existent data
        until time simulationTime.
        @param simulation_time: The simulation end time (in seconds).
        @return A list with the artist to be updated. 
        '''
        
        self.param['end_time'] = simulation_time
        
        # Load data from the data provider
        _,gconnections,gvalue = self.data_provider.get_synaptic_weights(**self.param)
        
        from mpi4py import MPI

        comm = MPI.COMM_WORLD
        
        process_id = comm.Get_rank()
        
        if (process_id==0):
            
            # Select the last weights
            gvalue = gvalue[:,-1]
                        
            for i in range(len(self.pattern)):
                pat = self.pattern[i]
            
                # Find the pattern fibers
                fibers = self.pattern_provider.fibers_in_pattern[pat]
                
                # Search fibers in connections
                connection_indexes = numpy.array([index for index in range(len(gconnections)) if (gconnections[index,0]==fibers).any()])
                
                # Select the weight values
                weight_values = gvalue[connection_indexes]
                
                # Select the number of the fibers
                source_fibers = gconnections[connection_indexes,0]
                
                # Select the activation values
                activation_values = self.pattern_provider.activation_levels[self.pattern_provider.pattern_id_index[pat][0],source_fibers]
                self.axesLines[i].set_xdata(activation_values)
                self.axesLines[i].set_ydata(weight_values)  
    
        return self.animated_artists