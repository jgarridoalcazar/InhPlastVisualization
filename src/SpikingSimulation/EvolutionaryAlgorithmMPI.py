
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

from mpi4py import MPI
from deap import base, creator, tools
from Utils.Utils import ReadConfigFile
from Utils.Logger import InitializeLogger, Logger2File
from deap.gp import mutUniform
import threading
import Queue

InitializeLogger('EvolutionaryAlgorithm')

# Get logger with default level to INFO
logger = logging.getLogger('EvolutionaryAlgorithm')

toolbox = base.Toolbox()

SIM_DATA = 0
SIM_EXIT = 1

def mutUniformCustom(individual, indpb, rand_generator):
    """This function applies a random uniform mutation between 0 and 1
     on the input individual. This mutation expects a
    :term:`sequence` individual composed of real valued attributes.
    The *indpb* argument is the probability of each attribute to be mutated.
    
    :param individual: Individual to be mutated.
    :param indpb: Independent probability for each attribute to be mutated.
    :returns: A tuple of one individual.
    
    
    """
    size = len(individual)
    
    for i in xrange(size):
        if rand_generator.rand() < indpb:
            individual[i] = rand_generator.rand()
    
    return individual,

class EvolutionaryAlgorithm(object):
    '''
    This class implements an evolutionary algorithm where the parameters are taken from the
    configuration file passed as a parameter.
    '''
    
    # Cell name translation
    operatorTranslatorDict = {
         'OnePoint' : tools.cxOnePoint,
         'TwoPoint' : tools.cxTwoPoint,
         'Gaussian' : tools.mutGaussian,
         'MutUniform' : mutUniformCustom,
         'Tournament' : tools.selTournament
    }
    
    operatorParamDict = {
         'OnePoint' : [],
         'TwoPoint' : [],
         'Gaussian' : ['gaussian_mu','gaussian_sigma','gaussian_indpb'],
         'Tournament' : ['tournament_size'],
         'MutUniform' : ['uniform_indpb']
    }
    
    paramTranslatorDict = {
         'gaussian_mu' : 'mu',
         'gaussian_sigma' : 'sigma',
         'gaussian_indpb' : 'indpb',
         'uniform_indpb'  : 'indpb',
         'tournament_size' : 'tournsize'    
    }
    
    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates a new object.
        @param config_file Name of the file with the options of the model.
        '''
        logger = logging.getLogger('EvolutionaryAlgorithm')
        
        if ('config_file' in kwargs):
            self.config_file = kwargs.pop('config_file')                         
        else:
            logger.error('Non-specified simulation config file')
            raise Exception('Non-DefinedSimulationConfig')
        
        super(EvolutionaryAlgorithm, self).__init__(**kwargs)
        
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
        
        # Number of generations    
        if 'number_of_generations' not in self.config_options['algorithm']:
            self.config_options['algorithm']['number_of_generations'] = 1
            logger.warning('Non-specified number_of_generations parameter. Using default value %s', self.config_options['algorithm']['number_of_generations'])
            
        # If number of individual is not defined, the number of available cores will be used.
        if 'number_of_individual' not in self.config_options['algorithm']:
            self.config_options['algorithm']['number_of_individual'] = 64
            logger.warning('Non-specified number_of_individual parameter. Using default value %s', self.config_options['algorithm']['number_of_individual'])
            
        if 'fill_idle_nodes' not in self.config_options['algorithm']:
            self.config_options['algorithm']['fill_idle_nodes'] = False
            logger.warning('Non-specified fill_idle_nodes parameter. Using default value %s', self.config_options['algorithm']['fill_idle_nodes'])
            
        # Crossover probability
        if 'crossover_probability' not in self.config_options['algorithm']:
            self.config_options['algorithm']['crossover_probability'] = 1.
            logger.warning('Non-specified crossover_probability parameter. Using default value %s', self.config_options['algorithm']['crossover_probability'])
        
        # Crossover operator
        if 'crossover_operator' not in self.config_options['algorithm']:
            self.config_options['algorithm']['crossover_operator'] = 'OnePoint'
            logger.warning('Non-specified crossover_operator parameter. Using default value %s', self.config_options['algorithm']['crossover_operator'])
        
        # Mutation probability
        if 'mutation_probability' not in self.config_options['algorithm']:
            self.config_options['algorithm']['mutation_probability'] = 1.
            logger.warning('Non-specified mutation_probability parameter. Using default value %s', self.config_options['algorithm']['mutation_probability'])
        
        # Mutation operator
        if 'mutation_operator' not in self.config_options['algorithm']:
            self.config_options['algorithm']['mutation_operator'] = 'Gaussian'
            logger.warning('Non-specified mutation_operator parameter. Using default value %s', self.config_options['algorithm']['mutation_operator'])
        
        # Selection operator
        if 'selection_operator' not in self.config_options['algorithm']:
            self.config_options['algorithm']['selection_operator'] = 'Tournament'
            logger.warning('Non-specified selection_operator parameter. Using default value %s', self.config_options['algorithm']['selection_operator'])
            
        # Hall of fame size
        if 'hall_of_fame_size' not in self.config_options['algorithm']:
            self.config_options['algorithm']['hall_of_fame_size'] = 1
            logger.warning('Non-specified hall_of_fame_size parameter. Using default value %s', self.config_options['algorithm']['hall_of_fame_size'])
        
        # Loading from file
        if 'load_from_file' not in self.config_options['algorithm']:
            self.config_options['algorithm']['load_from_file'] = None
        
        # Saving state parameters
        if 'saving_file' not in self.config_options['algorithm']:
            self.config_options['algorithm']['saving_file'] = None
            
        if 'saving_step' not in self.config_options['algorithm']:
            self.config_options['algorithm']['saving_step'] = 1
            logger.warning('Non-specified saving_step parameter. Using default value %s', self.config_options['algorithm']['saving_step'])
        
        # Initialize the simulation seeds if they have not been initialized    
        if 'simulation' not in self.config_options:
            self.config_options['simulation'] = dict()
        
        if 'seed' not in self.config_options['simulation']:
            self.config_options['simulation']['seed'] = time.time()
            
        # Extract parameters to explore    
        self._extract_parameters()
        
        # Initialize the evolutionary algorithm
        self._initialize_algorithm()
        
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
            
            if not 'type' in parameter:
                logger.error('Parameter evolution type has not been specified in %s',key)
                raise Exception('NonSpecifiedType')
            
            if parameter['type'] not in ['geometric','arithmetic']:
                logger.error('Parameter evolution type %s has not been implemented. Only geometric and arithmetic are allowed so far',parameter['type'])
                raise Exception('InvalidType')
            
    def _get_operator(self, parameter):
        # Check if the specified operator is included
        if parameter in self.operatorTranslatorDict:
            return self.operatorTranslatorDict[parameter]
        else:
            logger.error('The operator %s has not been mapped to an operator', parameter)
            raise Exception('Non-MappedOperator')
        
    def _get_operator_params(self, parameter, dicAlgorithm):
        # Retrieve the parameters of the operator.
        out_params = list()
        param_dic = dict()
        if parameter in self.operatorParamDict:
            for param in self.operatorParamDict[parameter]:
                if param in dicAlgorithm:
                    out_params.append(dicAlgorithm[param])
                    if param in self.paramTranslatorDict:
                        param_dic[self.paramTranslatorDict[param]] = dicAlgorithm[param]
                    else:
                        logger.error('The required operator parameter %s has not a translation', param)
                        raise Exception('Non-DefinedParameter')
                else:
                    logger.error('The required operator parameter %s has not been set', param)
                    raise Exception('Non-DefinedParameter')
        return param_dic
    
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
        
#         parent_conn, child_conn = multiprocessing.Pipe()
#         
#         p = multiprocessing.Process(target=helper_subprocess_simulation, args=(child_conn,local_config_options))
#         p.start()
#         
#                 
# #         # Catch SIGNINT just in case the parent process is killed before.
# #         import signal
# #         import sys
# #        
# #         def signal_term_handler(signal, frame):
# #             logger.info('Got %s. Killing running subprocesses',signal)
# #             if p.is_alive(): # Child still around?
# #                 p.terminate() # kill it
# #                 p.join()
# #             sys.exit(0)
# #        
# #         signal.signal(signal.SIGUSR2, signal_term_handler)
# #         signal.signal(signal.SIGINT, signal_term_handler)
# # #         signal.signal(signal.SIGKILL, signal_term_handler)
# #         signal.signal(signal.SIGTERM, signal_term_handler)
#         
#         mutual_information = parent_conn.recv()
#         p.join()
        mutual_information = helper_simulation(local_config_options)
        
        logger.info('Mutual information with seed %s and parameters %s: %s', seed, self._get_unnormalized_values(individual), mutual_information)

        return mutual_information
            
    def _initialize_algorithm(self):
        '''
        Initialize the evolutionary algorithm based on the provided parameters.
        '''
        
        self.num_generator = numpy.random.RandomState()
        
        # Create multiobjective optimization (maximize average MI and minimize Std)
        creator.create("FitnessMulti", base.Fitness, weights=(1.0,-1.0e-4))
        
        # Each individual inherits from list and add the FitnessMulti fitness function
        creator.create("Individual", list, fitness=creator.FitnessMulti)
        
        # Attribute generator (each attribute will be the normalized value -or the logartihm-)
        toolbox.register("attr_float", self.num_generator.rand)
        
        # Structure initializers
        toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_float, len(self.parameter_keys))

        # Population initializers        
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)
        
        # Operator registering
        #toolbox.register("evaluate", self._evaluate_population)
        
        # Crossover operator 
        operator = self._get_operator(self.config_options['algorithm']['crossover_operator'])
        paramOperator = self._get_operator_params(self.config_options['algorithm']['crossover_operator'], self.config_options['algorithm'])
        toolbox.register("mate", operator, **paramOperator)
        
        # Mutate operator
        operator = self._get_operator(self.config_options['algorithm']['mutation_operator'])
        paramOperator = self._get_operator_params(self.config_options['algorithm']['mutation_operator'], self.config_options['algorithm'])
        toolbox.register("mutate", operator, **paramOperator)
        toolbox.decorate("mutate", checkBounds())
        
        # Selection operator
        operator = self._get_operator(self.config_options['algorithm']['selection_operator'])
        paramOperator = self._get_operator_params(self.config_options['algorithm']['selection_operator'], self.config_options['algorithm'])
        toolbox.register("select", operator, **paramOperator)
        
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
        individualDict = dict()
        
        # List with the simulations to be executed in this "batch"
        simulationList = []
        
        request = None
        status = MPI.Status()
        population_size = 0
        output_population = toolbox.population(n=0)
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
                    runningDict[tuple_ind][0].append(data[0])
                    
                    # If all the simulations with these parameters are done, get the average and std
                    if (len(runningDict[tuple_ind][0])==self.config_options['algorithm']['number_of_repetitions']):
                        individual = individualDict.pop(tuple_ind)
                        individual.fitness.values = numpy.average(runningDict[tuple_ind][0]), numpy.std(runningDict[tuple_ind][0])
                        logger.debug('Fitness value calculated for individual %s: %s', individual, individual.fitness.values)
                        output_population.append(individual)
                        is_from_queue = runningDict[tuple_ind][1]
                        runningDict.pop(tuple_ind)
                        logger.debug('%s extracted from the running dictionary', tuple_ind)
                        # Check the number of individual to finish before unlocking the EA.
                        if is_from_queue:
                            population_size -= 1
                            if population_size==0:
                                self.completeQueue.put(output_population)
                                output_population = toolbox.population(n=0)
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
                    parameters, cur_seed, from_queue = simulationList.pop(0)                    
                    proc_rank = availableProcs.pop(0)
                    param_tuple = tuple(parameters)
                    
                    if param_tuple not in runningDict:
                        runningDict[param_tuple] = ([],from_queue)
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
                        ind_key = tuple(ind)
                        # Skip those individual already under evaluation
                        if ind_key not in individualDict:
                            population_size += 1
                            individualDict[ind_key] = ind
                            for seed in range(self.config_options['simulation']['seed'],self.config_options['simulation']['seed']+self.config_options['algorithm']['number_of_repetitions']):
                                simulationList.append((ind,seed,True))
                elif self.end_simulation:
                    proc_rank = availableProcs.pop(0)
                    ########################################
                    # Send and ending signal to the worker
                    ########################################
                    logger.debug('Sending ending signal to process %s', proc_rank)
                    data_send = numpy.empty(len(self.parameter_dic)+1, dtype=numpy.float64)
                    self.comm.Send([data_send, MPI.DOUBLE], dest=proc_rank, tag=SIM_EXIT)
                    endedProcs.append(proc_rank)
                elif self.config_options['algorithm']['fill_idle_nodes']:
                    # There is nothing to simulate. Generate a new individual
                    logger.debug('There are idle nodes. Creating random individuals.')
                    individual = toolbox.individual()
                    ind_key = tuple(individual)
                    # Skip those individual already under evaluation
                    if ind_key not in individualDict:
                        individualDict[ind_key] = individual
                        for seed in range(self.config_options['simulation']['seed'],self.config_options['simulation']['seed']+self.config_options['algorithm']['number_of_repetitions']):
                            simulationList.append((individual,seed,False))
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
        
        # Sequence of the parameter names
        param_names = []
        for param_dic in self.parameter_dic:
            param_names.append(param_dic['section'] + '.' + param_dic['parameter'])
        
        # Initialize stats
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", numpy.average, axis=0)
        stats.register("max", numpy.max, axis=0)
        stats.register("min", numpy.min, axis=0)
        stats.register("std", numpy.std, axis=0)
        
        # Initialize mapping to distribute the evaluations to the workers
        #toolbox.register("map", futures.map)
        
        # If we are in the last generation activate the flag
        self.last_generation = False
        self.end_simulation = False
        
        # Start simulation thread
        self.managerThread.start()
            
        
        # Load the previous algorithm state
        if self.config_options['algorithm']['load_from_file']:
            with open(self.config_options['algorithm']['load_from_file'], "r") as cp_file:
                cp = pickle.load(cp_file)
            self.population = cp["population"]
            start_gen = cp["generation"]+1
            halloffame = cp["halloffame"]
            logbook = cp["logbook"]
            self.num_generator.set_state(cp["rndstate"])
        else:
            self.num_generator.seed()
            self.population = toolbox.population(n=self.config_options['algorithm']['number_of_individual'])
            start_gen = 0
            halloffame = tools.HallOfFame(maxsize=self.config_options['algorithm']['hall_of_fame_size'])
            logbook = tools.Logbook()
            logbook.header = "gen", "evals", "fitness"
            logbook.chapters["fitness"].header = "avg", "max", "min", "std"

            # Fill idle nodes with new random individual
            #if self.config_options['algorithm']['fill_idle_nodes']:
            #    population = self._fill_idle_nodes(self.population)
        
            logger.debug("Start of evolution")
        
            self.population = self._evaluate_population(self.population)
            
            halloffame.update(self.population)
            
            logger.info('Parameter sequence: %s', param_names)
            logger.info('Hall of Fame:')
            for ind in halloffame:
                logger.info('Individual: %s. Fitness: %s', self._get_unnormalized_values(ind), ind.fitness.values)

        # Begin the evolution
        for gen in range(start_gen, self.config_options['algorithm']['number_of_generations']):
            logger.info("Generation %i", gen)
            
            # If we are in the last generation activate the flag
            self.last_generation = (gen == self.config_options['algorithm']['number_of_generations']-1)
        
            # Select the next generation individuals
            offspring = toolbox.select(self.population, k=self.config_options['algorithm']['number_of_individual'])
            # Clone the selected individuals
            offspring = list(map(toolbox.clone, offspring))
    
            # Apply crossover and mutation on the offspring
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if self.num_generator.rand() < self.config_options['algorithm']['crossover_probability']:
                    toolbox.mate(child1, child2)
                    del child1.fitness.values
                    del child2.fitness.values

            for mutant in offspring:
                if self.num_generator.rand() < self.config_options['algorithm']['mutation_probability']:
                    toolbox.mutate(mutant, rand_generator=self.num_generator)
                    del mutant.fitness.values
    
            # Fill idle nodes with new random individual
            #if self.config_options['algorithm']['fill_idle_nodes']:
            #    offspring = self._fill_idle_nodes(offspring)
            
            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            
            # Evaluate the population
            invalid_ind = self._evaluate_population(invalid_ind)
        
            # The population is entirely replaced by the offspring
            # population[:] = offspring
            # The population is extended with the offspring
            self.population.extend(invalid_ind)
            
            halloffame.update(invalid_ind)
            record = stats.compile(self.population)
            logbook.record(gen=gen, evals=len(invalid_ind), **record)
            
            # Saving evolution state
            if self.config_options['algorithm']['saving_file'] and gen % self.config_options['algorithm']['saving_step'] == 0:
                # Fill the dictionary using the dict(key=value[, ...]) constructor
                cp = dict(population=self.population, generation=gen, halloffame=halloffame,
                      logbook=logbook, rndstate=self.num_generator.get_state())

                with open(self.config_options['algorithm']['saving_file'], "wb") as cp_file:
                    pickle.dump(cp, cp_file)
                    
                logger.info('Evolution state saved in file %s', self.config_options['algorithm']['saving_file'])
        
            logger.info('Statistics in generation %s. %s evaluations', gen, len(invalid_ind))
            for key,value in record.items():
                logger.info('%s: %s', key, value)
            
            logger.info('Parameter sequence: %s', param_names)
            best_ind = tools.selBest(self.population, 1)[0]
            logger.info("Best individual in current population is %s, %s",self._get_unnormalized_values(best_ind), best_ind.fitness.values)
            logger.info('Hall of Fame:')
            for ind in halloffame:
                logger.info('Individual: %s. Fitness: %s', self._get_unnormalized_values(ind), ind.fitness.values)
    
            logger.info("-- End of (successful) evolution --")
            
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

# Function to check the bounds of the attributes in the interval
def checkBounds():
    def decorator(func):
        def wrapper(*args, **kargs):
            offspring = func(*args, **kargs)
            for child in offspring:
                for i in xrange(len(child)):
                    if child[i] > 1.0:
                        child[i] = 1.0
                    elif child[i] < 0.0:
                        child[i] = 0.0
            return offspring
        return wrapper
    return decorator

# Function creating a subprocess to launch nest simulation (it avoids NEST getting frozen after raising an exception in the previous simulation)
def helper_subprocess_simulation(pipe, local_config_options):
    import SpikingSimulation.FrequencySimulation as FrequencySimulation
    
    #logger.debug('Simulation parameter dictionary: %s', local_config_options)
        
    try:
        simulation = FrequencySimulation.FrequencySimulation(config_options = local_config_options)

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
    except:
        mutual_information = 0.0
        logger.debug('Exception caught on individual simulation %s', local_config_options)
        logger.info('Using default mutual information value %s', mutual_information)
        
    pipe.send(mutual_information)
    pipe.close()
    
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
    
        [mutual_information] = simulation.analyze_results()
    except KeyboardInterrupt:
        logger.warning('Received SIGNINT signal. Ending simulation')
        import sys
        sys.exit(0)
    except Exception as err:
        mutual_information = 0.0
        logger.debug('Exception caught on individual simulation %s: %s', local_config_options, err)
        logger.info('Using default mutual information value %s', mutual_information)
         
    return mutual_information



    

    
        
    
    
