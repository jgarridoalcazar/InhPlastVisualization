#! /usr/bin/env python

import sys
sys.path.append('./src')
    
import SpikingSimulation.FrequencySimulation as FrequencySimulation
from SpikingSimulation.Utils.Utils import ReadConfigParameters

arguments = ReadConfigParameters(sys.argv)

if 'config_file' not in arguments:
    print sys.argv
    print 'Error: Configuration file has not been specified. Usage:',sys.argv[0],'-c config_file [section.param=value]'
    sys.exit(1)

system_config_file = arguments.pop('config_file')

simulation = FrequencySimulation.FrequencySimulation(config_file = system_config_file)

for section in arguments.keys():
    for param in arguments[section].keys():
        simulation.config_options[section][param] = arguments[section][param]
            
simulation.initialize()
simulation.run_simulation()

# Plot the sin function
import numpy
import matplotlib.pylab
import matplotlib.pyplot as plt
from cmath import sin,pi

xaxes = numpy.arange(0,1./simulation.config_options['oscillations']['frequency'],1.e-3)
ithreshold = (simulation.config_options['mflayer']['eth']-simulation.config_options['mflayer']['erest'])*simulation.config_options['mflayer']['grest']*simulation.config_options['oscillations']['amplitude']
yaxes = ithreshold*numpy.sin(xaxes*2*pi*simulation.config_options['oscillations']['frequency']+simulation.config_options['oscillations']['phase'])
figureSin = plt.plot(xaxes,yaxes)
matplotlib.pylab.savefig('sinimage.svg', format='svg', dpi=1200)
matplotlib.pylab.show()

def create_figure_at_time(time):
    # Visualize the results
    import SpikingSimulation.Visualization.SimulFigure as SimulFigure
    import SpikingSimulation.Visualization.AxesRasterPlot as AxesRasterPlot
    import SpikingSimulation.Visualization.AxesNeuronPropertyLine as AxesNeuronPropertyLine
    import SpikingSimulation.Visualization.AxesPatternLine as AxesPatternLine
    
    
                
    figureA = SimulFigure.SimulFigure(simulation = simulation, numRows=2,numColumns=2,figsize=[23,14],dpi=80)
    figureA.add_subplot(fig_position=1,axes_type=AxesRasterPlot.AxesRasterPlot,
                         axes_parameters= {'data_provider':simulation.cerebellum,
                                           'pattern_provider':simulation.pattern_generator,
                                           'layer':'mflayer',
                                           'cell_index': range(50),
                                           'visible_data_only':True,
                                           'show_legend':False,
                                           'x_length':1.})
    figureA.add_subplot(fig_position=2,axes_type=AxesRasterPlot.AxesRasterPlot,
                    axes_parameters= {'data_provider':simulation.cerebellum,
                                      'layer':'goclayer',
                                      'visible_data_only':True,
                                      'show_legend':False,
                                      'x_length':1.})
    figureA.add_subplot(fig_position=3,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
                        axes_parameters= {'data_provider':simulation.cerebellum,
                                          'property':'Vm',
                                          'layer':'goclayer',
                                          'visible_data_only':True,
                                          'show_legend':False,
                                          'x_length': 1.,
                                          'y_window_lim': [-80e-3,-30e-3]})
    figureA.add_subplot(fig_position=4,axes_type=AxesPatternLine.AxesPatternLine,
                            axes_parameters= {'pattern_provider':simulation.pattern_generator,
                                              'visible_data_only':True,
                                              'show_legend':False,
                                              'x_length':1.})
     
    figureA.plot_at_time(simulation_time=time)

create_figure_at_time(2.)
matplotlib.pylab.savefig('figure2A.svg', format='svg', dpi=1200) 
matplotlib.pylab.show()

create_figure_at_time(190.)
matplotlib.pylab.savefig('figure2B.svg', format='svg', dpi=1200) 
matplotlib.pylab.show()
        
create_figure_at_time(1503.)
matplotlib.pylab.savefig('figure2C.svg', format='svg', dpi=1200) 
matplotlib.pylab.show()

def create_figure_weights_at_time(time):
    # Visualize the results
    import SpikingSimulation.Visualization.SimulFigure as SimulFigure
    import SpikingSimulation.Visualization.AxesFiringOffset as AxesFiringOffset
    import SpikingSimulation.Visualization.AxesActivationFiringOffset as AxesActivationFiringOffset
    import SpikingSimulation.Visualization.AxesWeightActivationPlot as AxesWeightActivationPlot
    
    
                
    figureA = SimulFigure.SimulFigure(simulation = simulation, numRows=2,numColumns=2,figsize=[23,14],dpi=80)
    figureA.add_subplot(fig_position=1,axes_type=AxesFiringOffset.AxesFiringOffset,
                        axes_parameters= {'data_provider':simulation.cerebellum,
                                          'oscillation_freq':simulation.config_options['oscillations']['frequency'],
                                          'layer':'goclayer',
                                          'visible_data_only':True})
    figureA.add_subplot(fig_position=2,axes_type=AxesActivationFiringOffset.AxesActivationFiringOffset,
                        axes_parameters= {'data_provider':simulation.cerebellum,
                                          'pattern_provider':simulation.pattern_generator,
                                          'oscillation_freq':simulation.config_options['oscillations']['frequency'],
                                          'layer':'mflayer',
                                          'visible_data_only':True})
    figureA.add_subplot(fig_position=3,axes_type=AxesWeightActivationPlot.AxesWeightActivationPlot,
                            axes_parameters= {'data_provider':simulation.cerebellum,
                                              'pattern_provider': simulation.pattern_generator,
                                              'layer':'mfgocsynapsis',
                                              'show_legend':False,
                                              'y_window_lim': [-0.01e-9,1.0e-9]}) 
    figureA.plot_at_time(simulation_time=time)
    
create_figure_weights_at_time(1500.)
matplotlib.pylab.savefig('figure2D.svg', format='svg', dpi=1200) 
matplotlib.pylab.show()