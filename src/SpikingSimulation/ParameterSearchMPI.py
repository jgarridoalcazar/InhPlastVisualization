
'''
Created on January 8, 2015

@author: Jesus Garrido (jgarridoalcazar at gmail.com)
'''
import time
import logging
import copy
import math
import numpy
import pickle
import itertools
import os.path

from mpi4py import MPI
from Utils.Utils import ReadConfigFile
from Utils.Logger import InitializeLogger, Logger2File
import threading
import Queue

InitializeLogger('ParameterSearch')

# Get logger with default level to INFO
logger = logging.getLogger('ParameterSearch')

SIM_DATA = 0
SIM_EXIT = 1


class ParameterSearch(object):
    '''
    This class implements exhaustive search into the parameter space taken from the
    configuration file passed as a parameter.
    '''
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates a new object.
        @param config_file Name of the file with the options of the model.
        '''
        logger = logging.getLogger('ParameterSearch')
        
        if ('config_file' in kwargs):
            self.config_file = kwargs.pop('config_file')                         
        else:
            logger.error('Non-specified simulation config file')
            raise Exception('Non-DefinedSimulationConfig')
        
        super(ParameterSearch, self).__init__(**kwargs)
        
        return
    
    def initialize_searcher(self):
        '''
        Initialize all the objects needed for running the simulation.
        '''
        
        self.comm = MPI.COMM_WORLD
        self.mpi_rank = self.comm.Get_rank()
        self.mpi_size = self.comm.Get_size()
        
        if (self.mpi_rank==0):
            self._initialize_master()
        else:
            self._initialize_worker()
        
        if (self.mpi_size == 1):
            logger.error("MPI Error. Only one MPI process has been created. No workers defined.")
            raise Exception('MPIError')        
        
        return
    
    def _initialize_logger(self):
        if 'log_file' in self.config_options['algorithm']:
            Logger2File(logger, self.config_options['algorithm']['log_file'])
        
        if 'verbosity' not in self.config_options['algorithm']:
            self.config_options['algorithm']['verbosity'] = 'debug'  
            logger.warning('Non-specified simulation verbosity. Using default value %s',self.config_options['algorithm']['verbosity']) 
        
        numeric_level = getattr(logging, self.config_options['algorithm']['verbosity'].upper(), None)
        if not isinstance(numeric_level, int):
            self.config_options['algorithm']['verbosity'] = 'info'
            numeric_level = getattr(logging, self.config_options['algorithm']['verbosity'].upper(), None)
            logger.warning('Invalid simulation verbosity. Using default value %s',self.config_options['algorithm']['verbosity']) 
            raise ValueError('Invalid log level: %s' % self.config_options['algorithm']['verbosity'])
            
        logger.setLevel(numeric_level)
        
    def _initialize_master(self):
        '''
        This function initializes the master process in the evolutionary algorithm.
        The master is in charge of providing individual to the workers. Thus, the master
        reads the algorithm configuration.
        '''
        
        logger.info('Parsing configuration file %s',self.config_file)
        self.config_options = ReadConfigFile(self.config_file)
        
        if 'algorithm' not in self.config_options:
            self.config_options['algorithm'] = dict()
        
        self._initialize_logger()
        
        if 'number_of_repetitions' not in self.config_options['algorithm']:
            self.config_options['algorithm']['number_of_repetitions'] = 1
        
        # Loading from file
        if 'load_from_file' not in self.config_options['algorithm']:
            self.config_options['algorithm']['load_from_file'] = None
        
        # Initialize the simulation seeds if they have not been initialized    
        if 'simulation' not in self.config_options:
            self.config_options['simulation'] = dict()
        
        if 'seed' not in self.config_options['simulation']:
            self.config_options['simulation']['seed'] = time.time()
        
        self._extract_parameters()
        
        # Saving state parameters
        if 'saving_file' not in self.config_options['algorithm']:
            self.config_options['algorithm']['saving_file'] = None
        else:
            # Remove saving file
            if self.config_options['algorithm']['load_from_file'] is None and os.path.isfile(self.config_options['algorithm']['saving_file']):
                logger.info('Removing existing result file %s', self.config_options['algorithm']['saving_file'])
                os.remove(self.config_options['algorithm']['saving_file'])
                param_names = '# '
                for param_dic in self.parameter_dic:
                    param_names = param_names + param_dic['section'] + '.' + param_dic['parameter'] + '\t'
                    
                param_names = param_names + 'Av.\tStd.\n'
        
                if self.config_options['algorithm']['saving_file'] is not None:
                    with open(self.config_options['algorithm']['saving_file'], 'a') as fileid:
                        logger.debug('Writing file header in the new file')
                        fileid.write(param_names)
        
        self.population = self._generate_config_tuples()
        
        self.population = self._extract_finished_tuples(self.population)
        
        # Initialize communication manager
        self.simulationQueue = Queue.Queue()
        self.completeQueue = Queue.Queue()
        self.managerThread = threading.Thread(target=self._manage_communications)
        
        return
    
    def _initialize_worker(self):
        '''
        This function initializes the worker process in the evolutionary algorithm.
        The workers are in charge of running the simulations with the parameters received from
        the master.
        '''
        
        logger.info('Parsing configuration file %s',self.config_file)
        self.config_options = ReadConfigFile(self.config_file)
        
        if 'algorithm' not in self.config_options:
            self.config_options['algorithm'] = dict()
        
        self._initialize_logger()
        
        if 'simulation' not in self.config_options:
            self.config_options['simulation'] = dict()
        
        # Set important undefined options
        if 'visualize_results' not in self.config_options['simulation']:
            self.config_options['simulation']['visualize_results'] = False
        
        if 'seed' not in self.config_options['simulation']:
            self.config_options['simulation']['seed'] = time.time()
            
        # Extract parameters to explore    
        self._extract_parameters()
        
        # Make a copy of the simulation options and extract the algorithm section
        self.simulation_options = copy.deepcopy(self.config_options)
        self.simulation_options.pop('algorithm')
        
        return
    
    def _extract_parameters(self):
        # Extract every parameter to explore
        self.parameter_keys = [key for key in self.config_options.keys() if key.startswith('parameter')]
        self.parameter_dic = []        
        for key in self.parameter_keys:
            self.parameter_dic.append(self.config_options.pop(key))
            
        for key,parameter in zip(self.parameter_keys,self.parameter_dic):
            # Check if the section and parameter exists
            if not 'section' in parameter:
                logger.error('Parameter section has not been specified in %s',key)
                raise Exception('NonSpecifiedSection')
            
            if parameter['section'] not in self.config_options:
                logger.error('Parameter section %s does not exist',parameter['section'])
                raise Exception('InvalidSection') 
        
            if not 'parameter' in parameter:
                logger.error('Parameter name has not been specified in %s',key)
                raise Exception('NonSpecifiedParameter')
            
            if parameter['parameter'] not in self.config_options[parameter['section']]:
                logger.error('Parameter %s does not exist in section %s',parameter['parameter'],parameter['section'])
                raise Exception('InvalidParameter')

            if not 'min_value' in parameter:
                logger.error('Parameter minimum values has not been specified in %s',key)
                raise Exception('NonSpecifiedMinValue')
            
            if not 'max_value' in parameter:
                logger.error('Parameter maximum values has not been specified in %s',key)
                raise Exception('NonSpecifiedMaxValue')
            
            if not 'num_values' in parameter:
                logger.error('Parameter number of values has not been specified in %s', key)
                raise Exception('NonSpecifiedNumValues')
            
            if not 'type' in parameter:
                logger.error('Parameter evolution type has not been specified in %s',key)
                raise Exception('NonSpecifiedType')
            
            if parameter['type'] not in ['geometric','arithmetic']:
                logger.error('Parameter evolution type %s has not been implemented. Only geometric and arithmetic are allowed so far',parameter['type'])
                raise Exception('InvalidType')
            
    def _generate_config_tuples(self):
        
        # Generate the combinations of values
        value_list = []
        for key,parameter in zip(self.parameter_keys,self.parameter_dic):
            # Arithmetic series
            values = list(numpy.linspace(0.0, 1.0, num=parameter['num_values']))
            
            value_list.append(values)
        
        if len(value_list):
            # Generate the combinations of values
            combinations = list(itertools.product(*value_list))
        else:
            combinations = list()
            
        return combinations

    def _extract_finished_tuples(self,population):
        
        if self.config_options['algorithm']['load_from_file'] is not None and os.path.isfile(self.config_options['algorithm']['load_from_file']):
            loaded_values = [tuple(row) for row in numpy.loadtxt(self.config_options['algorithm']['load_from_file'], usecols=tuple(range(len(self.parameter_keys))))]
            
            unnorm_population = [tuple(self._get_unnormalized_values(individual)) for individual in population]
            
            for row in loaded_values:
                RemIndex = -1;
                for index, ind in enumerate(unnorm_population):
                    if all(abs((val1-val2)/val1)<1e-5 for val1,val2 in zip(row,ind)):
                        logger.debug('%s simulation is already loaded from the file. Removing simulation',row)
                        print index
                        RemIndex = index
                        break
                
                if RemIndex!=-1:
                    del unnorm_population[RemIndex]
                    del population[RemIndex]
                    
        return population
            
            
    def _get_unnormalized_values(self, individual):
        unnorm_values = []
        
        for norm_value, param_dic in zip(individual,self.parameter_dic):
            min_value = param_dic['min_value']
            max_value = param_dic['max_value']
            
            if param_dic['type'] == 'arithmetic':
                value = norm_value*(max_value - min_value) + min_value
            elif param_dic['type'] == 'geometric':
                logmin = math.log10(abs(min_value))
                logmax = math.log10(abs(max_value))
                value = 10.0**(norm_value*(logmax - logmin))*min_value
            
            unnorm_values.append(value)
            
        return unnorm_values
        
        
    
    def _eval_fitness_funct(self,individual, seed):
        
        # Make a copy of the simulation config options
        local_config_options = copy.deepcopy(self.simulation_options)
        
        unnorm_values = self._get_unnormalized_values(individual)
        
        for unnorm, param_dic in zip(unnorm_values,self.parameter_dic):
            local_config_options[param_dic['section']][param_dic['parameter']] = unnorm
        
        local_config_options['simulation']['seed'] = seed
            
        logger.info('Running evaluation with seed %s and parameters %s', seed, self._get_unnormalized_values(individual))
        
        mutual_information = helper_simulation(local_config_options)
        
        logger.info('Mutual information with seed %s and parameters %s: %s', seed, self._get_unnormalized_values(individual), mutual_information)

        return mutual_information
            
    def _evaluate_population(self, population):
        
        # Insert the population into the simulation queue and unlock it
        self.simulationQueue.put(population)
        self.simulationQueue.task_done()
        
        self.end_simulation = self.last_generation
        
        logger.info("Evaluating %i individuals",len(population))
        return self.completeQueue.get()
        
    def _manage_communications(self):
        '''
        This function manages the simulation queue, sending the simulations to the other MPI processes.
        It manages the two simulation queues (SimulationQueue -jobs to be done- and CompleteQueue -jobs finished-).
        This function is thought to be executed in a sepparate thread of the master process.
        '''
        
        # Initialize SimulationMap and RunningDict
        simulationMap = dict()
        availableProcs = range(1,self.mpi_size)
        endedProcs = []
        for ind in availableProcs:
            simulationMap[ind] = None
            
        runningDict = dict()
        
        # List with the simulations to be executed in this "batch"
        simulationList = []
        
        status = MPI.Status()
        population_size = 0
        output_population = []
        data = numpy.empty(1, dtype=numpy.float64)
        
        ########################################
        # Create requests with MPI.Irecv(....)
        ########################################
        request = self.comm.Irecv([data, MPI.DOUBLE], source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG)
        
        # This loop end when every worker process has been finished
        while (len(endedProcs)!=(self.mpi_size-1)):
            Job_finished = request.Test(status)
            
            if Job_finished or (not availableProcs):
                if (not Job_finished):
                    logger.debug('Waiting for something finished')
                    request.Wait(status)
                # There is at least one simulation finished
                mpi_process = status.Get_source()
                tuple_ind = simulationMap[mpi_process]
                logger.debug('%s mutual information has been received from %s: %s', tuple_ind, mpi_process, data[0])
                
                if tuple_ind not in runningDict:
                    logger.warning('Error in data received from process %s',mpi_process)
                    logger.warning('%s not exist in runningDict %s',tuple_ind,runningDict)
                else:
                    runningDict[tuple_ind].append(data[0])
                    
                    # If all the simulations with these parameters are done, get the average and std
                    if (len(runningDict[tuple_ind])==self.config_options['algorithm']['number_of_repetitions']):
                        avMI = numpy.average(runningDict[tuple_ind]), numpy.std(runningDict[tuple_ind])
                        logger.debug('Fitness value calculated for individual %s: %s', tuple_ind, avMI)
                        
                        unnorm_val = self._get_unnormalized_values(tuple_ind)
                        unnorm_val.extend(list(avMI))
                        
                        if self.config_options['algorithm']['saving_file'] is not None:
                            with open(self.config_options['algorithm']['saving_file'], 'a') as fileid:
                                logger.debug('Saving fitness value calculated to %s',self.config_options['algorithm']['saving_file'])
                                numpy.savetxt(fileid, numpy.transpose(unnorm_val), fmt="%.15e", delimiter="\t", newline='\t')
                                fileid.write('\n')
                                    
                        output_population.append(tuple_ind)
                        runningDict.pop(tuple_ind)
                        logger.debug('%s extracted from the running dictionary', tuple_ind)
                        # Check the number of individual to finish before unlocking the EA.
                        population_size -= 1
                        if population_size==0:
                            self.completeQueue.put(output_population)
                            output_population = []
                            self.completeQueue.task_done()
                            logger.debug('Simulation batch has been finished')
            
                
                # Set the process as available
                simulationMap[mpi_process] = None
                availableProcs.append(mpi_process)
                
                ########################################
                # Create requests with MPI.Irecv(....)
                ########################################
                request = self.comm.Irecv([data, MPI.DOUBLE], source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG)
            elif availableProcs:
                # There are available processors
                if simulationList:
                    # There are simulations ready to be launched
                    # Extract the first simulation and the first available process
                    parameters, cur_seed = simulationList.pop(0)                    
                    proc_rank = availableProcs.pop(0)
                    param_tuple = tuple(parameters)
                    
                    if param_tuple not in runningDict:
                        runningDict[param_tuple] = []
                        logger.debug('%s inserted in the running dictionary', param_tuple)
                        
                    simulationMap[proc_rank] = param_tuple
                    
                    ########################################
                    # Send through MPI the parameters and the seed
                    ########################################
                    logger.debug('Sending parameters %s and seed %s to process %s', parameters, cur_seed, proc_rank)
                    data_send = numpy.empty(len(parameters)+1, dtype=numpy.float64)
                    for idx, parameter in enumerate(parameters):
                        data_send[idx] = parameter
                    data_send[-1] = cur_seed
                    self.comm.Send([data_send, MPI.DOUBLE], dest=proc_rank, tag=SIM_DATA)
                elif not self.simulationQueue.empty():
                    # There are available batch simulations in the simulation queue.
                    # Add every individual/seed combination to the simulationList
                    population = self.simulationQueue.get()
                    logger.debug('New population to be evaluated received: %s',population)
                    for ind in population:
                        # Skip those individual already under evaluation
                        population_size += 1
                        for seed in range(self.config_options['simulation']['seed'],self.config_options['simulation']['seed']+self.config_options['algorithm']['number_of_repetitions']):
                            simulationList.append((ind,seed))
                elif self.end_simulation:
                    proc_rank = availableProcs.pop(0)
                    ########################################
                    # Send and ending signal to the worker
                    ########################################
                    logger.debug('Sending ending signal to process %s', proc_rank)
                    data_send = numpy.empty(len(self.parameter_dic)+1, dtype=numpy.float64)
                    self.comm.Send([data_send, MPI.DOUBLE], dest=proc_rank, tag=SIM_EXIT)
                    endedProcs.append(proc_rank)
                else:
                    # Nothing to do
                    logger.debug('Sleeping 1')
                    time.sleep(0.1)
            else:
                # Nothing to do
                logger.debug('Sleeping 2')
                time.sleep(0.1)

    def execute_search(self):
        '''
        Initialize all the objects needed for running the simulation.
        '''
        
        if (self.mpi_rank==0):
            self._execute_search_master()
        else:
            self._execute_search_worker()
        
        return
    
    def _execute_search_master(self):
        '''
        The master node executes the genetic algorithm and provides simulation parameters
        to the workers.
        '''
        
        # If we are in the last generation activate the flag
        self.last_generation = True
        self.end_simulation = False
        
        # Start simulation thread
        self.managerThread.start()
            
        
        logger.debug("Start of simulation")
        
        self.population = self._evaluate_population(self.population)
        
        logger.debug("Simulation ended")
            
        return
    
    def _execute_search_worker(self):
        '''
        Worker nodes receive parameter lists and simulate the network.
        '''
        stay_working = True
        my_status = MPI.Status()
        while stay_working:
            # Receive the simulation parameters and seed
            data_recv = numpy.empty(len(self.parameter_keys)+1, dtype=numpy.float64)
            self.comm.Recv([data_recv, MPI.DOUBLE], source=0, tag=MPI.ANY_TAG, status = my_status)
            tag = my_status.Get_tag()
             
            # Check the tag
            if tag == SIM_EXIT:
                stay_working = False
                continue
            
            if tag != SIM_DATA:
                logger.warning('Unknown tag %s received in worker', tag)
            
            cur_seed = int(data_recv[-1])
            parameters = data_recv[:-1].tolist()
            
            logger.debug('Received parameters %s and seed %s', parameters, cur_seed)
            
            # Launch the simulation with the parameters
            mutual_information = self._eval_fitness_funct(parameters, cur_seed)
            
            logger.debug('Sending mutual information value %s to process root process', mutual_information)
            
            send_array = numpy.array([mutual_information], dtype=numpy.float64)
        
            self.comm.Send([send_array, MPI.DOUBLE], dest=0, tag=SIM_DATA)
                    
        return            

# Function creating a subprocess to launch nest simulation (it avoids NEST getting frozen after raising an exception in the previous simulation)
def helper_simulation(local_config_options):
    import SpikingSimulation.FrequencySimulation as FrequencySimulation
    
    # logger.debug('Simulation parameter dictionary: %s', local_config_options)
        
    try:
        simulation = FrequencySimulation.FrequencySimulation(config_options = local_config_options)
    
        simulation.initialize()
    
        if simulation.config_options['simulation']['visualize_animation']:
            simulation.visualize_animation()
        else:
            simulation.run_simulation()
        
        if simulation.config_options['simulation']['visualize_results']:
            simulation.visualize_results()
    
        parameter_keys_MI = [key for key in simulation.config_options.keys() if key.startswith('mutual_information')]
        parameter_keys_hit = [key for key in simulation.config_options.keys() if key.startswith('hit_analysis')]
        parameter_keys_hit_top = [key for key in simulation.config_options.keys() if key.startswith('hit_top_analysis')]
        if parameter_keys_MI:
            [mutual_information] = simulation.analyze_MI()
        elif parameter_keys_hit:
            [mutual_information] = simulation.analyze_Hits()
        elif parameter_keys_hit_top:
            [mutual_information] = simulation.analyze_Hits_Top()
    except KeyboardInterrupt:
        logger.warning('Received SIGNINT signal. Ending simulation')
        import sys
        sys.exit(0)
    except Exception as err:
        mutual_information = 0.0
        logger.debug('Exception caught on individual simulation %s: %s', local_config_options, err)
        logger.info('Using default mutual information value %s', mutual_information)
         
    return mutual_information



    

    
        
    
    
