#! /usr/bin/env python

import SpikingSimulation.FrequencySimulation as FrequencySimulation
import sys
from SpikingSimulation.Utils.Utils import ReadConfigParameters

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

    if simulation.config_options['simulation']['visualize_animation']:
        simulation.visualize_animation()
    else:
        simulation.run_simulation()
                
    if simulation.config_options['simulation']['visualize_results']:
        simulation.visualize_results()
    
    simulation.analyze_results()    
