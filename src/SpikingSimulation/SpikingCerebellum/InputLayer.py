'''
Created on March 22, 2014

@author: Jesus Garrido (jgarridoalcazar at gmail.com)
'''
import logging
import numpy
import time

logger = logging.getLogger('Simulation')

class InputLayer(object):
    '''
    This class defines an input layer to a neural network and the data needed to generate it in a simulator.
    '''
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates a new input layer (composed of fibers). Either the number of neurons or the density of neurons and
        the dimensions (length, width, height)-tuple has to be specified 
        @param number_of_neurons: Number of cells in the layer (optional).
        @param register_activity: Boolean parameter indicating whether the activity in this layer will be registered.
        @param minindex: Index of the first cell in the network (optional).
        @param size: (length, width, height)-tuple with the size of the volume to fill (in mm) (optional).
        @param density_of_neurons: Density of neurons in the volume (optional).
        @param random_generator: The random number generator to use (optional).
        '''
        
        # Read name
        self.__name__ = kwargs.pop('name')
        
        self.number_of_neurons = None
        self.density_of_neurons = None
        self.size = None
        
        # Read number of neurons
        if ('number_of_neurons' in kwargs):
            self.number_of_neurons = kwargs.pop('number_of_neurons')
        else:
            logger.warning('Non-specified number of neurons in layer %s',self.__name__)
            
        # Read volume size
        if ('size' in kwargs):
            self.size = numpy.array(kwargs.pop('size'))
            # Calculate the total volume of the dice
            self.volume = numpy.prod(self.size)
        else:
            logger.warning('Non-specified volume size in layer %s',self.__name__)
            
        # Read volume size
        if ('density_of_neurons' in kwargs):
            self.density_of_neurons = kwargs.pop('density_of_neurons')
        else:
            logger.warning('Non-specified density of neurons in layer %s',self.__name__)
            
        if self.number_of_neurons is not None:
            if self.size is not None:
                if self.density_of_neurons is not None:
                    logger.warning('Both number of neurons and density of neurons specified in layer %s. Density of neurons will be discarded',self.__name__)
                self.density_of_neurons = self.number_of_neurons / self.volume
        elif self.size is not None and self.density_of_neurons is not None: 
            self.number_of_neurons = round(self.volume * self.density_of_neurons)
        else:
            logger.error('Non-specified neither number of neurons nor (density of neurons and size) in layer %s',self.__name__)
            raise Exception('Non-DefinedProperty')
                
        
        # Read register activity parameter
        if ('register_activity' in kwargs):
            self.register_activity = kwargs.pop('register_activity') 
        else:
            logger.warning('Non-specified register_actiity parameter in layer %s. Using default value false',self.__name__)
            self.register_activity = False
            
        # Read minindex parameter
        if ('minindex' in kwargs):
            self.MinIndex = kwargs.pop('minindex') 
        else:
            self.MinIndex = None
            
        # Read seed parameter
        if ('random_generator' in kwargs):
            self.random_generator = kwargs.pop('random_generator') 
        else:
            self.random_generator = numpy.random.RandomState(int(time.time()))
            
        # Check whether additional parameters have been used.
        for param in kwargs:
            logger.warning('Unrecognized parameter %s in layer %s',param,self.__name__)
            
        # Locate the elements if required
        if self.size is not None:
            self._locate_elements_()
            
        return
            
            
    def _locate_elements_(self):
        '''
        This function locate every neuron in its specific place. The size needs to be specified before calling this function
        '''
        
        if self.size is None:
            logger.error('Non-specified size in layer %s. Neurons cannot be spatially located',self.__name__)
            raise Exception('Non-DefinedProperty')
                
        n_dimensions = len(self.size)
        
        self.relative_positions = self.random_generator.uniform(0,1,(self.number_of_neurons, n_dimensions)) 
        
        return
    
    def get_relative_coordinates(self):
        return self.relative_positions
    
    def get_absolute_coordinates(self):
        return self.relative_positions * self.size
                
