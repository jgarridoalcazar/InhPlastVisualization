#! /usr/bin/env python

import matplotlib.animation as animation
import numpy
import SimulFigure
import time

class SimulAnimation (SimulFigure.SimulFigure,animation.TimedAnimation):
    '''
    This class extends the SimulFigure including the animation behavior.
    '''
       
#     def _blit_draw(self, artists, bg_cache):
#         # Handles blitted drawing, which renders only the artists given instead
#         # of the entire figure.
#         updated_ax = []
#         for a in artists:
#             # If we haven't cached the background for this axes object, do
#             # so now. This might not always be reliable, but it's an attempt
#             # to automate the process.
#             if a.axes not in bg_cache:
#                 # bg_cache[a.axes] = a.figure.canvas.copy_from_bbox(a.axes.bbox)
#                 # change here
#                 bg_cache[a.axes] = a.figure.canvas.copy_from_bbox(a.axes.figure.bbox)
#             a.axes.draw_artist(a)
#             updated_ax.append(a.axes)
#     
#         # After rendering all the needed artists, blit each axes individually.
#         for ax in set(updated_ax):
#             # and here
#             # ax.figure.canvas.blit(ax.bbox)
#             ax.figure.canvas.blit(ax.figure.bbox)
    
    
    def __init__(self, **kwargs):
        """
        Create the figure/animation
        
        Keyword arguments (see CustomFigure and TimedAnimation documentation)
        @param init_time: The initial time of the simulation (default 0)
        @param end_time: The final time of the simulation. Obligatory parameter.
        @param frame_rate: Frame rate (in Hz). 
        """
        
        # Get init_time parameter
        if ('init_time' in kwargs):
            self.init_time = kwargs.pop('init_time',None)
        else:
            self.init_time = 0
            
        # Get end_time parameter
        if ('end_time' in kwargs):
            self.end_time = kwargs.pop('end_time',None)
        else:
            print 'Obligatory end_time parameter not provided'
            raise Exception('NonProvidedParameter','end_time')
        
        # Get frame_rate parameter
        if ('frame_rate' in kwargs):
            self.frame_rate = kwargs.pop('frame_rate',None)
        else:
            self.frame_rate = 10
            
        super(SimulAnimation, self).__init__(**kwargs)
        
        self.last_update = time.time()
        
        return
    
    

    def _draw_frame(self, simulation_time):
        """
        Generate and show the movie steps
        
        Keyword arguments:
        simulation_time: The current relative simulation time.
        """
        new_update = time.time()
        
        print 'Refresh rate:',1./(new_update-self.last_update),'FPS'
        self.last_update = new_update
        
        # Run the simulation until that time
        self.simulation.run_simulation(end_time=simulation_time)
        
        self._drawn_artists = self.update(simulation_time = simulation_time)
                 
        return
    
    def new_frame_seq(self):
        return iter(numpy.arange(self.init_time,self.end_time,1.0/self.frame_rate))
    
    def _init_draw(self):
        # Run the simulation until that time
        self.simulation.run_simulation(end_time=self.init_time)
        
        self.update(init_time=self.init_time,end_time=self.init_time)
    

