#! /usr/bin/env python

import SpikingSimulation.FrequencySimulation as FrequencySimulation
import sys
import time
import logging
from SpikingSimulation.Utils.Utils import ReadConfigParameters

logger = logging.getLogger('Simulation')

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

    start = time.time()
    
    simulation.initialize()
    
    end_init = time.time()

    logger.info('Simulation init finished: %ss', end_init - start)

    if simulation.new_config_options['simulation']['visualize_animation']:
        simulation.visualize_animation()
    else:
        simulation.run_simulation()
                
    if simulation.new_config_options['simulation']['visualize_results']:
        simulation.visualize_results()

    end_simulation = time.time()

    logger.info('Simulation finished: %ss', end_simulation - end_init)    
    
    # Save network state and activity at the end of the simulation
    simulation.cerebellum.save_activity()
    simulation.cerebellum.save_network_state()

    end_saving = time.time()

    logger.info('Saving network finished: %ss', end_saving - end_simulation)    
    
    simulation.analyze_MI()    
    
    simulation.analyze_Hits()
    
    simulation.analyze_Hits_Top()

    end_analysis = time.time()

    logger.info('Analysis finished: %ss', end_analysis - end_saving)    
    