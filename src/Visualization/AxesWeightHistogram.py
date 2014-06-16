import numpy
import AxesPlot
from mpi4py import MPI

class AxesWeightHistogram(AxesPlot.AxesPlot):
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
        @param max_weight Maximum weight to consider in the histogram.
        @param num_bins Number of bins to include in the histogram
        '''
        
        
        super(AxesWeightHistogram, self).__init__(**kwargs)
                
        # Get data_provider parameter 
        if ('data_provider' in kwargs):
            self.data_provider = kwargs.pop('data_provider',None)
        else:
            print 'Obligatory data_provider parameter not provided'
            raise Exception('NonProvidedParameter','data_provider')
        
        # Get layer name parameter
        if ('layer' in kwargs):
            self.layer = kwargs.pop('layer',None)
        else:
            print 'Obligatory layer parameter not provided'
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
        
            
        # Get max_weight parameter 
        if ('max_weight' in kwargs):
            self.max_weight = kwargs.pop('max_weight',None)
        else:
            self.max_weight = None
            
        # Get num_bins parameter 
        if ('num_bins' in kwargs):
            self.num_bins = kwargs.pop('num_bins',None)
        else:
            self.num_bins = 10
            
        
                
                
    def initialize(self):
        '''
        Perform all the required operations needed in order to initialize the axes.
        It sets the title, axes titles, axes limits, legend and creates the lines.
        The DataProvider object must be initialized before calling this function.
        '''
        self.figure_title = 'Weight Histogram - ' + self.layer
        self.figure_x_label = 'Weight (S)'
        self.figure_y_label = 'Count'
        
        # Set connection parameters
        self.param = dict()
        
        self.param['synaptic_layer'] = self.layer
        if self.target_indexes:
            self.param['target_indexes'] = self.target_indexes
            
        if self.source_indexes:
            self.param['source_indexes'] = self.source_indexes
        
        synaptic_layer = self.data_provider.layer_map[self.layer]
        
        if self.max_weight is None:
            if synaptic_layer.learning_rule_type:
                self.max_weight = synaptic_layer.learning_rule_parameters['max_weight']
            else:
                self.max_weight = 1.
        
        width = self.max_weight/self.num_bins
        self.positions = numpy.arange(0,self.max_weight,width)
        
        self.param['end_time'] = 0.0
        _,gcon,_ = self.data_provider.get_synaptic_weights(**self.param)
        
        self.axesRect = self.axes.bar(self.positions, [0]*len(self.positions), self.max_weight/self.num_bins)
        
        comm = MPI.COMM_WORLD
        
        process_id = comm.Get_rank()
        
        if (process_id==0):
            self.axes.set_xlim([0,self.max_weight])
            self.axes.set_ylim([0,len(gcon)])
        
        super(AxesWeightHistogram, self).initialize()
            
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
        gtime,_,gvalue = self.data_provider.get_synaptic_weights(**self.param)
        
        comm = MPI.COMM_WORLD
        
        process_id = comm.Get_rank()
        
        if (process_id==0):
            
            weights = gvalue[:,-1]
            frequencies,_ = numpy.histogram(weights, bins=self.positions)
            
            for rect, f in zip(self.axesRect, frequencies):
                rect.set_height(f)
            
        return self.axesRect