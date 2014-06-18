#! /usr/bin/env python

import SpikingSimulation.FrequencySimulation as FrequencySimulation
import time
import sys
from mpi4py import MPI 

if __name__ == "__main__":
    
    if len(sys.argv)==1:
        print sys.argv
        print 'Error: Configuration file has not been specified. Usage:',sys.argv[0],'config_file'
        sys.exit(1)
    
    system_config_file = sys.argv[1]
    
    simulation = FrequencySimulation.FrequencySimulation(config_file = system_config_file)
    
    simulation.initialize()
    
    comm = MPI.COMM_WORLD

    init_time = time.time()
    
    #if comm.Get_size()==1:
    simulation.visualize_results()
    #else:
    #    simulation.run_simulation()
    
    end_time = time.time()
    
    print 'Elapsed time:',end_time-init_time
    
    
    pass
    
