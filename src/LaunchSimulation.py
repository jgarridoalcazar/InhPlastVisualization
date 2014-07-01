#! /usr/bin/env python

import SpikingSimulation.FrequencySimulation as FrequencySimulation
import sys

if __name__ == "__main__":
    
    if len(sys.argv)==1:
        print sys.argv
        print 'Error: Configuration file has not been specified. Usage:',sys.argv[0],'config_file'
        sys.exit(1)
    
    system_config_file = sys.argv[1]
    
    simulation = FrequencySimulation.FrequencySimulation(config_file = system_config_file)
                
    simulation.initialize()

    if simulation.config_options['simulation']['visualize_results']:
        simulation.visualize_results()
    else:
        simulation.run_simulation()
                
    
    simulation.analyze_results()    
