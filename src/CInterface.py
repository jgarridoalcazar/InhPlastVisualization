#! /usr/bin/env python
import logging
import sys
import os

sys.path.append(os.path.dirname(__file__))

from SpikingSimulation.Utils.Logger import InitializeLogger
from SpikingSimulation.Utils.Utils import ExpandConfigParameters

InitializeLogger('Simulation')

# Get logger with default level to INFO
logger = logging.getLogger('Simulation')

def RunSimulation(kwargs):
    """This function servers as a C interface to run a simulation with
    a given configuration.
    
    :param config_file_name: Name of the file with the basic configuration. Obligatory parameter
    :param Any other configuration parameter in the way dictionary[section][param] = value.
    :returns The mutual information divided by the maximum mutual information.
    """
    
    import SpikingSimulation.FrequencySimulation as FrequencySimulation
    
    config_file = None
    if 'config_file_name' in kwargs:
        config_file = kwargs.pop('config_file_name')
        
    config_options = ExpandConfigParameters(**kwargs)
    
    try:
        if config_file is not None:
            
            simulation = FrequencySimulation.FrequencySimulation(config_file = config_file)
            
            for section in config_options.keys():
                for param in config_options[section].keys():
#                     logger.info('Adding configuration option %s.%s=%s',section,param,config_options[section][param])
                    simulation.config_options[section][param] = config_options[section][param]
        else:
#             logger.info('Initializing simulation with options %s',kwargs)
            simulation = FrequencySimulation.FrequencySimulation(config_options = config_options)
        
        simulation.initialize()
    
        if simulation.config_options['simulation']['visualize_animation']:
            simulation.visualize_animation()
        else:
            simulation.run_simulation()
        
        if simulation.config_options['simulation']['visualize_results']:
            simulation.visualize_results()
    
        [mutual_information] = simulation.analyze_results()
    except KeyboardInterrupt:
        logger.warning('Received SIGNINT signal. Ending simulation')
        import sys
        sys.exit(0)
    except Exception as err:
        mutual_information = 0.0
        logger.debug('Exception caught on individual simulation %s: %s', simulation.config_options, err)
        logger.info('Using default mutual information value %s', mutual_information)
         
    return mutual_information
    