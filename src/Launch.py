#! /usr/bin/env python

import SpikingSimulation.FrequencySimulation as FrequencySimulation
import time
from mpi4py import MPI 

if __name__ == "__main__":
    
    system_config_file = './SimulationConfig.cfg'
    
    simulation = FrequencySimulation.FrequencySimulation(config_file = './SimulationConfig.cfg')
    
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
    
