'''
Created on Oct 9, 2013

@author: jgarrido
'''

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import logging

logger = logging.getLogger('Simulation')

class SimulFigure(object):
    '''
    This figure extends the class matplotlib.pyplot.Figure to provide the auto-location of the new axes
    that are added to the Figure.
    '''


    def __init__(self, **kwargs):
        '''
        Constructor. It creates an empty figure in which you can add subplots that automatically distribute
        on the figure.
        
        Keyword arguments
        @param simulation: The Frequency simulation object to run the simulation as required. 
        @param numcolumns: Number of columns to include in the figure.
        @param numrows: Number of rows to include in the figure.
        '''
        # Get simulation parameter
        if ('simulation' in kwargs):
            self.simulation = kwargs.pop('simulation',None)
        else:
            logger.error('Obligatory simulation parameter not provided')
            raise Exception('NonProvidedParameter','simulation')
        
        # Get numColumns parameter
        if ('numColumns' in kwargs):
            self.numColumns = kwargs.pop('numColumns',None)
        else:
            self.numColumns = 4
        
        # Get numRows parameter
        if ('numRows' in kwargs):
            self.numRows = kwargs.pop('numRows',None)
        else:
            self.numRows = 8
                
        self.axesList = []
        
        self.figure = plt.figure(**kwargs)
        
        self.figure.subplots_adjust(hspace=self.numRows*0.1)
        
        if (issubclass(type(self),animation.TimedAnimation)):
            #super(SimulFigure, self).__init__(self.figure,interval=1000./self.frame_rate,repeat=False,blit=False)
            super(SimulFigure, self).__init__(self.figure,interval=0,repeat=False,blit=self.blit)
        else:
            self.blit = False
            super(SimulFigure, self).__init__()
        
    def add_subplot(self, **kwargs):
        '''
        It adds a new subplot to the figure that is placed in one of the free places.
        
        keyword parameters:
        @param axes_type The name of a class that inherits from AxesType (type of axes). Obligatory parameter.
        @param fig_position The position of the axes in the figure (from 1). Obligatory parameter.
        @param axes_parameters Dictionary with the parameters of the axes_type to be used. Optional parameter.
        '''
        
        # Get axes_type parameter
        if ('axes_type' in kwargs):
            self.axes_type = kwargs.pop('axes_type',None)
        else:
            logger.error('Obligatory axes_type parameter not provided')
            raise Exception('NonProvidedParameter','axes_type')
            
        # Get fig_position position
        if ('fig_position' in kwargs):
            self.position = kwargs.pop('fig_position',None)
        else:
            logger.error('Obligatory fig_position parameter not provided')
            raise Exception('NonProvidedParameter','fig_position')
        
        # Get axes_parameters
        if ('axes_parameters' in kwargs):
            axes_type_arguments = kwargs.pop('axes_parameters',None)
        else:
            axes_type_arguments = dict()
            
        axes_type_arguments['figure'] = self
            
        # Insert position parameters and create the subplot
        args = (self.numRows, self.numColumns, self.position)
        
        ax1 = self.figure.add_subplot(*args,**kwargs)
        
        # Create the link between axes and data
        self.axesList.append(self.axes_type(axes=ax1,**axes_type_arguments))
        
        self.axesList[len(self.axesList)-1].initialize()
        
        
    def update(self, **kwargs):
        '''
        It updates the data of all the axes in the figure. The user will have 
        to call figure.canvas.draw() or pass all these elements to the property
        _drawn_artists (if an animation) to get the figure refreshed.
        
        @param simulation_time: Time until which the simulation will be shown. Optional argument. If none is specified, total simulation time will be used.
        @return All the elements (mainly lines) that have been updated.
        '''
        
        objectsToUpdate = []
        
        for axes in self.axesList:
            objectsToUpdate.extend(axes.drawAtTime(**kwargs))
            
        if ('simulation_time' in kwargs):
            self.figure.canvas.set_window_title(str(kwargs['simulation_time']) + ' seconds')
       
        return objectsToUpdate
    
    def plot_at_time(self, **kwargs):
        """
        Generate the snapshot of the simulation at a certain time.
        
        Keyword arguments:
        @param simulation_time: Time until which the simulation will be shown. Optional argument. If none is specified, total simulation time will be used.
        """
        
        if 'simulation_time' in kwargs:
            simulation_time = kwargs['simulation_time']
        else:
            simulation_time = self.simulation.simulation_time
            kwargs['simulation_time'] = simulation_time
        
        # Run the simulation until that time
        self.simulation.run_simulation(end_time=simulation_time)
        
        self.update(**kwargs)
        
        return
            
            
        
        
