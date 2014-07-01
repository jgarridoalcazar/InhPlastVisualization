'''
Created on May 27, 2014

@author: Jesus Garrido (jgarridoalcazar at gmail.com)
'''
import FrequencySimulation
import numpy
import math
import itertools
import copy
import logging
import os
import time
import subprocess

import scipy.interpolate
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

from popen2 import popen2

from Utils.Utils import ReadConfigFile, WriteConfigFile
from Utils.Logger import InitializeLogger

InitializeLogger('ParamSearcher')

# Get logger with default level to INFO
logger = logging.getLogger('ParamSearcher')

class ParameterSearch(object):
    '''
    This class defines launch succesive simulations to explore one or more parameters.
    '''

    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates a new simulation object.
        @param config_file Name of the file with the options of the model.
        '''
        
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
        
        logger.info('Parsing configuration file %s',self.config_file)
        self.config_options = ReadConfigFile(self.config_file)
        
        if 'simulation' in self.config_options:
            if 'debug' in self.config_options['simulation'] and self.config_options['simulation']['debug']:
                logger.setLevel(logging.DEBUG)
        
        # Set important undefined options
        if 'visualize_results' not in self.config_options['simulation']:
            self.config_options['simulation']['visualize_results'] = False
            
        if 'launcher' not in self.config_options:
            self.config_options['launcher'] = dict()
        
        # By default 1 process is used for the simulations
        if 'num_mpi_processes' not in self.config_options['launcher']:
            self.config_options['launcher']['num_mpi_processes'] = 1
        
        # By default qsub is not enabled    
        if 'use_qsub' in self.config_options['launcher'] and self.config_options['launcher']['use_qsub']:
            self.launch_funct = self._launch_qsub_simulation
        elif self.config_options['launcher']['num_mpi_processes']>1:
            self.config_options['launcher']['use_qsub'] = False
            self.launch_funct = self._launch_mpi_simulation
        else:
            self.config_options['launcher']['use_qsub'] = False
            self.launch_funct = self._launch_serial_simulation
                
        # By default 1 thread is used for the simulations
        if 'num_omp_threads' not in self.config_options['launcher']:
            self.config_options['launcher']['num_omp_threads'] = 1
            
        if 'parallel_environment' not in self.config_options['launcher'] and self.config_options['launcher']['use_qsub'] and \
                    (self.config_options['launcher']['num_omp_threads'] > 1 or self.config_options['launcher']['num_mpi_processes'] > 1):
            logger.error('Non-specified parallel environment for qsub job submission')
            raise Exception('Non-DefinedParallelEnvironment')
            
        if 'queue_name' not in self.config_options['launcher'] and self.config_options['launcher']['use_qsub']:
            logger.error('Non-specified queue name for qsub job submission')
            raise Exception('Non-DefinedQueueName')
            
        if 'mpi_host_file' not in self.config_options['launcher']:
            self.config_options['launcher']['mpi_host_file'] = None
        
        if 'mpi_launcher' not in self.config_options['launcher']:
            self.config_options['launcher']['mpi_launcher'] = 'mpirun'
        
        if 'python_exec' not in self.config_options['launcher']:
            self.config_options['launcher']['python_exec'] = 'python'
            
        if 'nest' not in self.config_options:
            self.config_options['nest'] = dict()
            
        if 'number_of_virtual_processes' not in self.config_options['nest']:
            self.config_options['nest']['number_of_virtual_processes'] = self.config_options['launcher']['num_mpi_processes']
        
        return
    
    def _generate_config_dicts(self):
        
        # Extract every parameter to explore
        self.parameter_keys = [key for key in self.config_options.keys() if key.startswith('parameter')]
        self.parameter_dic = []        
        for key in self.parameter_keys:
            self.parameter_dic.append(self.config_options.pop(key))
        
        # Generate the combinations of values
        value_list = []
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
            
            if not 'step' in parameter:
                logger.error('Parameter step has not been specified in %s',key)
                raise Exception('NonSpecifiedStep')
            
            if not 'type' in parameter:
                logger.error('Parameter evolution type has not been specified in %s',key)
                raise Exception('NonSpecifiedType')
            
            if parameter['type']=='geometric':
                # Geometric series
                min_exponent = math.log(parameter['min_value'],parameter['step'])
                max_exponent = math.log(parameter['max_value'],parameter['step'])
                values = list(numpy.logspace(min_exponent, max_exponent, num=int(max_exponent-min_exponent)+1, base=parameter['step']))
            elif parameter['type']=='arithmetic':
                # Arithmetic series
                values = list(numpy.linspace(parameter['min_value'], parameter['max_value'], num=int((parameter['max_value']-parameter['min_value'])/float(parameter['step']))+1))
            else:
                logger.error('Parameter evolution type %s has not been implemented. Only geometric and arithmetic are allowed so far',parameter['type'])
                raise Exception('InvalidType')
            
            value_list.append(values)
        
        if len(value_list):
            # Generate the combinations of values
            combinations = list(itertools.product(*value_list))
            simulation_options = list() 
            
            for tuple_act in combinations:
                # Copy the dictionary and change every single parameter 
                options_copy = copy.deepcopy(self.config_options)
            
                logger.info('Setting parameters to the following values: %s',tuple_act)
                sim_name = ''
                for param_dic, value in zip(self.parameter_dic,tuple_act):
                    options_copy[param_dic['section']][param_dic['parameter']] = value
                    sim_name += '_' + str(value)
                    
                options_copy['simulation']['simulation_name'] += sim_name
                simulation_options.append(options_copy)
        else:
            simulation_options = [options_copy]
            
        return simulation_options

    
    def execute_search(self):
        '''
        Initialize all the objects needed for running the simulation.
        '''
    
        simulation_options = self._generate_config_dicts()
                        
        self.launch_funct(config_options=simulation_options)
            
    def _launch_serial_simulation(self, config_options):
        '''
        Launch a new simulation according to the proposed method.
        '''
        
        for index,config in enumerate(config_options):
            logger.info('Launching simulation %s of %s', index, len(config_options))
            # Create, initialize and launch the simulation
            logger.debug('Creating the simulation object')
            simulation = FrequencySimulation.FrequencySimulation(config_options = config)
        
            logger.info('Initializing the simulation')
            simulation.initialize()
        
            logger.info('Running the simulation')
            if self.config_options['simulation']['visualize_results']:
                simulation.visualize_results()
            else:
                simulation.run_simulation()
            logger.info('Simulation ended')
            logger.info('Analyzing results')
            simulation.analyze_results()
            logger.info('Analysis ended')
        
    def _launch_mpi_simulation(self, config_options):
        '''
        Launch a new simulation according to the proposed method.
        @param index: Index of the simulation.
        @param config_options List of dictionaries with the configuration to be used for the simulation
        '''
        
        for index,config in enumerate(config_options):
            logger.debug('Writing configuration file for job %s',index)
            file_name = self._save_configuration_file(config)
        
            mpi_command = []
            # mpi_command.append(config_options['launcher']['mpi_launcher'])
            # mpi_command.append('-np')
            # mpi_command.append(str(config_options['launcher']['num_mpi_processes']))
            mpi_command.append(str(config['launcher']['python_exec']))
            mpi_command.append('./src/LaunchSimulation.py')
            mpi_command.append(file_name)
        
            # Create, initialize and launch the simulation
            logger.info('Calling MPI process for simulation %s of %s', index, len(config_options))
            logger.debug(mpi_command)
            subprocess.call(mpi_command)
            logger.info('Simulation ended')
    
    def _get_configuration_file_name(self, config_options):
        '''
        Generate the name of the configuration file.
        '''
        
        # Create configuration file
        if 'data_path' in config_options['simulation']:
            data_path = config_options['simulation']['data_path']
        else:
            data_path = './results'
            config_options['simulation']['data_path'] = data_path
        
        if 'simulation_name' in config_options['simulation']:
            data_path += '/' + config_options['simulation']['simulation_name']
                
        if not os.path.exists(data_path):
            logger.info('Creating result folder %s', data_path)
            os.makedirs(data_path)
        
        file_name = data_path + '/SimulationConfig.cfg'
        
        return file_name
        
    def _save_configuration_file(self, config_options):
        '''
        Create the configuration file according to the config_options parameters.
        @param config_options Dictionary with the configuration to be used for the simulation
        '''
        
        file_name = self._get_configuration_file_name(config_options)
        
        if os.path.isfile(file_name):
            logger.warning('A configuration file %s already exists. It will overwrite that simulation file', file_name)
        
        logger.debug('Writing configuration file %s',file_name)    
        WriteConfigFile(config_options, file_name)
        
        return file_name
            
    def _launch_qsub_simulation(self,config_options):
        '''
        Launch a qsub job to run a simulation by using config_options parameters.
        @param config_options List of dictionaries with the configuration to be used for the simulation
        '''
        
        job_table_file = self.config_options['simulation']['data_path'] + '/' + self.config_options['simulation']['simulation_name'] + '.txt'
        f = open(job_table_file,'w')
        
        for index,config in enumerate(config_options):
            if 'simulation_name' not in config['simulation']:
                config['simulation']['simulation_name'] = 'job' + str(index)
            
            logger.debug('Writing configuration file for job %s',index)
            file_name = self._save_configuration_file(config)
            
            # Writhe the file name into the job array table
            f.write(file_name+'\n') # python will convert \n to os.linesep
        
        f.close() # you can omit in most cases as the destructor will call if

        
        # Create the job submission script
        buf = '#!/bin/sh\n'
        buf += '#$ -S /bin/sh\n'
        buf += '#$ -t 1-' + str(len(config_options)) + '\n'
        buf += '#$ -N ' + self.config_options['simulation']['simulation_name'] + '\n'
        buf += '#$ -o ' + self.config_options['simulation']['data_path'] + '/\n'
        buf += '#$ -M jesusgarrido@ugr.es\n'
        buf += '#$ -m abe\n'
        buf += '#$ -j y\n'
        buf += '#$ -cwd\n'
        buf += '#$ -V\n'
        buf += '#$ -v OMP_NUM_THREADS=' + str(self.config_options['launcher']['num_omp_threads']) + '\n'
        buf += '#$ -q ' + self.config_options['launcher']['queue_name'] + '\n'
        # buf += '#$ -l ' + self.config_options['launcher']['queue_name'] + ',h_rt=2:00:00,h_cpu=2:00:00\n' # Set maximum cpu time to 12 hours
        if (self.config_options['launcher']['num_omp_threads']>1 or self.config_options['launcher']['num_mpi_processes']>1):
            buf += '#$ -pe ' + self.config_options['launcher']['parallel_environment'] + ' ' + str(self.config_options['launcher']['num_omp_threads']*self.config_options['launcher']['num_mpi_processes']) + '\n'
        
        buf += '\nPARAM_FILE=' + job_table_file + '\n'
        buf += 'PARAM=$(cat $PARAM_FILE | head -n $SGE_TASK_ID | tail -n 1)\n\n'
        
        buf += 'mpirun -n ' + str(self.config_options['launcher']['num_mpi_processes']) + ' ' + self.config_options['launcher']['python_exec'] + ' ./src/LaunchSimulation.py $PARAM\n'

        logger.debug('Generated qsub script:')
        logger.debug(buf)        
        logger.info('Launching qsub command')
        # Open a pipe to the qsub command.
        output, inputstr = popen2('qsub')
        
        # Send job_string to qsub
        inputstr.write(buf)
        inputstr.close()
 
        # Print your job and the response to the screen
        logger.info(output.read())
        
    def visualize_results(self):
        '''
        Visualize the results only if the number of parameters is 1 or 2. For the moment it only works
        with two parameters, one pattern and one cell.
        '''

        # Generate the config option dictionaries
        simulation_options = self._generate_config_dicts()
        
        if len(self.parameter_dic)!=2:
            logger.error('Number of parameters is different than 2. It can not be represented')
            raise Exception('InvalidNumberOfParameters')
        
        # Generate the labels for each axis    
        labels = []    
        axis = []
        values = []
        for param in self.parameter_dic:
            labels.append(param['section']+'.'+param['parameter'])
            axis.append([])
    
        # Extract every mutual information to explore
        parameter_keys = [key for key in self.config_options.keys() if key.startswith('mutual_information')]
        for key in parameter_keys:        
            # Generate the configuration options
            for index,config in enumerate(simulation_options):
                logger.debug('Loading configuration file for job %s',index)
                file_name = self._get_configuration_file_name(config)
                
                # Load the MI file
                mi_file = os.path.dirname(os.path.realpath(file_name)) + '/' + key
                value = numpy.loadtxt(mi_file)
                values.append(float(value))
                
                for index,param_dic in enumerate(self.parameter_dic):
                    axis[index].append(config[param_dic['section']][param_dic['parameter']])


                    
            # Create the figure    
            fig = plt.figure()
            #ax = fig.gca(projection='3d')
            ax = fig.gca()
            x = numpy.array(axis[0])*1e9
            y = numpy.array(axis[1])*1.
            z = values

            # Set up a regular grid of interpolation points
            xi, yi = numpy.linspace(min(x), max(x), 100), numpy.linspace(min(y), max(y), 100)
            xi, yi = numpy.meshgrid(xi, yi)

            # Interpolate; there's also method='cubic' for 2-D data such as here
            zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method='linear')

            surf = ax.imshow(zi, vmin=min(z), vmax=max(z), origin='lower', extent=[min(x), max(x), min(y), max(y)])
            #surf = ax.plot_surface(xi, yi, zi, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
            ax.scatter(x, y, c=z)
            fig.colorbar(surf, shrink=0.5, aspect=5)
            # Interpolate the data to generate the mesh
            ax.set_title(key)
            ax.set_xlabel(labels[0])
            ax.set_ylabel(labels[1])
            fig.savefig(key+'.png')
            

            
        #plt.show()
            
    
        
        
    