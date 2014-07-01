import abc
import logging

logger = logging.getLogger('Simulation')

class Analysis(object):
    '''
    This class defines an abstract base class with all the methods
    needed in order to make some analysis over the simulation data.
    '''
    __metaclass__ = abc.ABCMeta
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates a new object..
        '''
        
        return 
               
    
    def initialize(self):
        '''
        Perform all the required operations needed in order to initialize the analysis.
        '''
        
        return
    
    @abc.abstractmethod
    def runAtTime(self, **kwargs):
        '''
        This function updates the analysis according to the existent data
        until time simulationTime.
        '''
        return
    
    @abc.abstractmethod
    def writeToFile(self, **kwargs):
        '''
        This function writes the results into a file.
        @param file_name Name of the file where the data will be stored
        '''
        return