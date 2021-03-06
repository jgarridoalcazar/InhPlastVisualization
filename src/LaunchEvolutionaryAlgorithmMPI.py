#! /usr/bin/env python

import SpikingSimulation.EvolutionaryAlgorithmMPI as EvolutionaryAlgorithm
import sys

if len(sys.argv)==1:
    print sys.argv
    print 'Error: Configuration file has not been specified. Usage:',sys.argv[0],'config_file'
    sys.exit(1)

system_config_file = sys.argv[1]

simulation = EvolutionaryAlgorithm.EvolutionaryAlgorithm(config_file = system_config_file)

simulation.initialize_searcher()

simulation.execute_search()

    
    