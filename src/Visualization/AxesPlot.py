import abc

class AxesPlot(object):
    '''
    This class defines an abstract base class with all the methods
    needed in order to link an axes an its data provider.
    '''
    __metaclass__ = abc.ABCMeta
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates an object linking the axes with the DataProvider.
        @param axes: An axes object with the axes which has just been created to plot the data.
        '''
        
        # Get axesType parameter
        if ('axes' in kwargs):
            self.axes = kwargs.pop('axes',None)
        else:
            print 'Obligatory axes parameter not provided'
            raise Exception('NonProvidedParameter','axes')
        
        return 
               
    
    def initialize(self):
        '''
        Perform all the required operations needed in order to initialize the axes.
        It sets the title, axes titles, axes limits, legend and creates the lines.
        The DataProvider object must be initialized before calling this function.
        '''
        
        # Set axes title
        self.axes.set_title(self.figure_title)
        self.axes.set_ylabel(self.figure_y_label)
        self.axes.set_xlabel(self.figure_x_label)
        
        return
    
    @abc.abstractmethod
    def drawAtTime(self, **kwargs):
        '''
        This function updates all the elements of the axes according to the existent data
        until time simulationTime.
        @return A list with the artist to be updated. 
        '''
        return