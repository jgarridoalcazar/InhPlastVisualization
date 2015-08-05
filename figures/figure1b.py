#! /usr/bin/env python

import sys
sys.path.append('./src')
    
import SpikingSimulation.FrequencySimulation as FrequencySimulation
from SpikingSimulation.Utils.Utils import ReadConfigParameters
from cmath import sin,pi

if __name__ == "__main__":
    
    
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
    
    # Plot the sin function
    import numpy
    import matplotlib.pylab
    import matplotlib.pyplot as plt
    xaxes = numpy.arange(0,simulation.config_options['simulation']['time'],1.e-3)
    ithreshold = (simulation.config_options['mflayer']['eth']-simulation.config_options['mflayer']['erest'])*simulation.config_options['mflayer']['grest']*simulation.config_options['oscillations']['amplitude']
    yaxes = ithreshold*numpy.sin(xaxes*2*pi*simulation.config_options['oscillations']['frequency']+simulation.config_options['oscillations']['phase'])
    figure8 = plt.plot(xaxes,yaxes)
    matplotlib.pylab.savefig('sinimage.svg', format='svg', dpi=1200)

    simulation.run_simulation()
    
    # Print stimulation currents
    print 'Pattern lengths:'
    print simulation.pattern_generator.pattern_length_cum
    print 'Pattern numbers:'
    print simulation.pattern_generator.pattern_id
    print 'Activation levels:'
    print simulation.pattern_generator.activation_levels
                
    # Visualize the results
    import SpikingSimulation.Visualization.SimulFigure as SimulFigure
    import SpikingSimulation.Visualization.AxesRasterPlot as AxesRasterPlot
    
    
                
    figure7 = SimulFigure.SimulFigure(simulation = simulation, numRows=1,numColumns=1,figsize=[23,14],dpi=80)
    figure7.add_subplot(fig_position=1,axes_type=AxesRasterPlot.AxesRasterPlot,
                         axes_parameters= {'data_provider':simulation.cerebellum,
                                           'pattern_provider':simulation.pattern_generator,
                                           'layer':'mflayer',
                                           'visible_data_only':True,
                                           'show_legend':True,
                                           'x_length':2.})
     
    figure7.plot_at_time()
    #matplotlib.pylab.show()
    matplotlib.pylab.savefig('spikes.svg', format='svg', dpi=1200) 
    
    
        
    