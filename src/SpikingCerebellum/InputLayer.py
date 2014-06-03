'''
Created on March 22, 2014

@author: Jesus Garrido (jgarridoalcazar at gmail.com)
'''

class InputLayer(object):
    '''
    This class defines an input layer to a neural network and the data needed to generate it in a simulator.
    '''
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates a new input layer (composed of fibers).
        @param number_of_neurons: Number of cells in the layer.
        @param register_activity: Boolean parameter indicating whether the activity in this layer will be registered.
        '''
        
        # Read name
        self.__name__ = kwargs.pop('name')
        
        # Read number of neurons
        if ('number_of_neurons' in kwargs):
            self.number_of_neurons = kwargs.pop('number_of_neurons')
        else:
            print 'Non-specified number of neurons.'
            raise Exception('Non-DefinedProperty')
        
        # Read register activity parameter
        if ('register_activity' in kwargs):
            self.register_activity = kwargs.pop('register_activity') 
        else:
            print 'Non-specified register_actiity parameter. Using default value false.'
            self.register_activity = False
            
        # Check whether additional parameters have been used.
        for param in kwargs:
            print 'Unrecognized parameter ',param,'in layer',self.__name__
