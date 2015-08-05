'''
Created on May 27, 2014

@author: Jesus Garrido (jgarridoalcazar at gmail.com)
'''
            
import bisect
import ntpath
import logging
from Utils.Utils import ReadConfigFile
from Utils.Logger import InitializeLogger, Logger2File

InitializeLogger('Simulation')

# Get logger with default level to INFO
logger = logging.getLogger('Simulation')

class FrequencySimulation(object):
    '''
    This class defines a simulation where the parameters are taking from the
    configuration file passed as a parameter.
    '''

    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates a new simulation object. 1 of the parameters have to be defined
        @param config_options Dictionary with the parameters of the simulations.
        @param config_file Name of the file where the simulation parameters are stored.
        '''
        
        if 'config_options' in kwargs:
            self.config_options = kwargs.pop('config_options')
            # This code forces exception in the GSLSolver
#             self.config_options['goclayer']['epsilon_rr_ip'] = 5.7832127
#             self.config_options['mfgocsynapsis']['max_weight'] = 1.49209e-9
#             self.config_options['mfgocsynapsis']['learning_step'] = 1.0142578e-5
#             self.config_options['goclayer']['beta_ip'] = 1.713656
#             self.config_options['mfgocsynapsis']['minus_plus_ratio'] = 1.77053
#             self.config_options['goclayer']['epsilon_rc_ip'] = 476.5509369
#             self.config_options['goclayer']['tau_ip'] =111.14285
        elif 'config_file' in kwargs:
            self.config_file = kwargs.pop('config_file')
#             logger.info('Parsing configuration file %s',self.config_file)
            self.config_options = ReadConfigFile(self.config_file)
        else:
            logger.error('Non-specified simulation configuration options or configuration file')
            raise Exception('Non-DefinedSimulationConfig')
        
        super(FrequencySimulation, self).__init__(**kwargs)
        
        return
    
    def initialize(self):
        '''
        Initialize all the objects needed for running the simulation.
        '''
        # Read simulation general options
        if 'simulation' not in self.config_options:
            self.config_options['simulation'] = dict()
            
        if 'log_file' in self.config_options['simulation']:
            Logger2File(logger, self.config_options['simulation']['log_file'])
        
        if 'verbosity' not in self.config_options['simulation']:
            self.config_options['simulation']['verbosity'] = 'debug'  
            logger.warning('Non-specified simulation verbosity. Using default value %s',self.config_options['simulation']['verbosity']) 
    
        numeric_level = getattr(logging, self.config_options['simulation']['verbosity'].upper(), None)
        if not isinstance(numeric_level, int):
            self.config_options['simulation']['verbosity'] = 'info'
            numeric_level = getattr(logging, self.config_options['simulation']['verbosity'].upper(), None)
            logger.warning('Invalid simulation verbosity. Using default value %s',self.config_options['verbosity']['verbosity']) 
            raise ValueError('Invalid log level: %s' % self.config_options['simulation']['verbosity'])
        
        logger.setLevel(numeric_level)
        
        if 'use_mpi' not in self.config_options['simulation']:
            self.config_options['simulation']['use_mpi'] = False
        
        if 'time' in self.config_options['simulation']:
            self.simulation_time = self.config_options['simulation']['time']
        else:
            self.simulation_time = 1
                
                
        if 'visualize_animation' not in self.config_options['simulation']:
            self.config_options['simulation']['visualize_animation'] = False
            
        if 'visualize_results' not in self.config_options['simulation']:
            self.config_options['simulation']['visualize_results'] = False
   
            
        logger.debug('Simulation time fixed to %ss',self.simulation_time)
        
        self.new_config_options = self.config_options
            
        # Initialize cerebellar model
        logger.debug('Creating cerebellum generator')
        if 'run_simulation' in self.config_options['simulation'] and self.config_options['simulation']['run_simulation']:
            # Nest has to be imported before mpi4py
            if self.config_options['simulation']['use_mpi']:
                import SpikingCerebellum.NestCerebellarModel as NestGenerator
            else:
                import SpikingCerebellum.NestCerebellarModelNoMPI as NestGenerator
             
            self.cerebellum = NestGenerator.NestCerebellarModel(config_dict=self.config_options)
        else:
            self.config_options['simulation']['run_simulation'] = False
            # Get the path of the config_file
            import SpikingCerebellum.SavedCerebellarModel as SavedGenerator
            path = self.config_options['simulation']['data_path'] + '/' + self.config_options['simulation']['simulation_name']
            # Read the old configuration file being saved with the simulation and containing specific network information
            self.config_options = ReadConfigFile(path+'/'+'SimulationConfig.cfg')
            self.cerebellum = SavedGenerator.SavedCerebellarModel(config_dict=self.config_options, simulation_folder=path)
    
        logger.debug('Initializing cerebellum generator')
        self.cerebellum.initialize_simulation()
    
        # Initialize oscillatory input current
        if 'oscillations' in self.config_options:
            logger.debug('Creating AC Current generator')
            self.cerebellum.add_ac_current(**self.config_options['oscillations'])
            
        # Initialize frequency stimulation input current
        if 'stimulation' in self.config_options:   
            logger.debug('Creating DC Current generator')
            self.config_options['stimulation']['simulation_time'] = self.simulation_time
            self.config_options['stimulation']['number_of_fibers'] = self.cerebellum.mflayer.number_of_neurons
            self.config_options['stimulation']['seed'] = self.config_options['simulation']['seed']
            
            import Stimulation.PatternGenerator as PatternGenerator
            self.pattern_generator = PatternGenerator.PatternGenerator(**self.config_options['stimulation'])
            self.pattern_generator.initialize()
            
            self.pattern_length, self.pattern_activations = self.pattern_generator.get_all_patterns()
            self.pattern_length_cum = self.pattern_generator.pattern_length_cum
            
        self.current_time = 0.
            
        
    def run_simulation(self, **kwargs):
        '''
        Run the simulation according to the configuration file.
        @param end_time Time until when simulation will be run
        '''
        
        if 'end_time' in kwargs:
            end_time = kwargs.pop('end_time')
            
            if end_time>self.simulation_time:
                logger.warning('Simulation time is shorter than end_time. Simulating %ss',self.simulation_time)
                
            end_time = min(end_time,self.simulation_time)
        else:
            end_time = self.simulation_time
            
        logger.info('Running the simulation from %ss until time %ss',self.current_time, end_time)
            
        if 'stimulation' in self.config_options and self.new_config_options['simulation']['run_simulation']:
            init_index = bisect.bisect_left(self.pattern_length_cum, self.current_time)
            end_index = bisect.bisect_left(self.pattern_length_cum,end_time)
        
            for index in range(init_index,end_index+1):
                sim_time = min(self.pattern_length_cum[index]-self.current_time,end_time-self.current_time)
                
                # Substitution with nest step_current_generator is prefered, but it runs slower. 
                self.cerebellum.set_dc_current(amplitude=self.pattern_activations[index])
            
                logger.debug('Running the simulation %ss until %ss', sim_time, self.cerebellum.simulation_time+sim_time) 
                self.cerebellum.simulate_network(sim_time)
                
                self.current_time = self.cerebellum.simulation_time
        else:
            sim_time = self.simulation_time-self.current_time
            logger.debug('Running the simulation %ss until %ss', sim_time, self.cerebellum.simulation_time+sim_time) 
            self.cerebellum.simulate_network(sim_time)
            self.current_time = self.cerebellum.simulation_time
            
    def visualize_results(self):
        '''
        Visualize the results of the simulation
        '''
        import Visualization.SimulFigure as SimulFigure
        import Visualization.AxesNeuronPropertyLine as AxesNeuronPropertyLine
        import Visualization.AxesPatternLine as AxesPatternLine
        import Visualization.AxesRasterPlot as AxesRasterPlot
        import Visualization.AxesWeightEvolutionLine as AxesWeightEvolutionLine
        import Visualization.AxesWeightHistogram as AxesWeightHistogram
        import Visualization.AxesWeightActivationPlot as AxesWeightActivationPlot
        import Visualization.AxesFiringOffset as AxesFiringOffset
        import Visualization.AxesActivationFiringOffset as AxesActivationFiringOffset
        import matplotlib.pylab
        
                
#         figure7 = SimulFigure.SimulFigure(simulation = self, numRows=4,numColumns=1,figsize=[23,14],dpi=80)
#         figure7.add_subplot(fig_position=1,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'property':'Vm',
#                                               'layer':'goclayer',
#                                               'visible_data_only':True,
#                                               'show_legend':False,
#                                               'x_length': 1.})
#         figure7.add_subplot(fig_position=2,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'property':'Gexc',
#                                               'layer':'goclayer',
#                                               'visible_data_only':True,
#                                               'show_legend':False,
#                                               'x_length': 1.})
#         figure7.add_subplot(fig_position=3,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'property':'Ginh',
#                                               'layer':'goclayer',
#                                               'visible_data_only':True,
#                                               'show_legend':False,
#                                               'x_length': 1.})
#         animation.add_subplot(fig_position=2,axes_type=AxesRasterPlot.AxesRasterPlot,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'pattern_provider':self.pattern_generator,
#                                               'layer':'mflayer',
#                                               'cell_index': range(50),
#                                               'visible_data_only':True,
#                                               'show_legend':True,
#                                               'x_length':1.})
#         figure7.add_subplot(fig_position=2,axes_type=AxesWeightEvolutionLine.AxesWeightEvolutionLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'layer':'mfgocsynapsis',
#                                               'source_indexes': range(100),
#                                               'target_indexes': range(1),
#                                               'visible_data_only':True,
#                                               'show_legend':False})
#         figure7.add_subplot(fig_position=3,axes_type=AxesWeightHistogram.AxesWeightHistogram,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'layer':'mfgocsynapsis',
#                                               'num_bins': 60})
#         figure7.add_subplot(fig_position=4,axes_type=AxesWeightActivationPlot.AxesWeightActivationPlot,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'pattern_provider': self.pattern_generator,
#                                               'layer':'mfgocsynapsis'})
#         figure7.plot_at_time()

        figure8 = SimulFigure.SimulFigure(simulation = self, numRows=4,numColumns=4,figsize=[23,14],dpi=80)
        for i in range(16):
            figure8.add_subplot(fig_position=i+1,axes_type=AxesWeightActivationPlot.AxesWeightActivationPlot,
                            axes_parameters= {'data_provider':self.cerebellum,
                                              'pattern_provider': self.pattern_generator,
                                              'layer':'mfgocsynapsis',
                                              'target_indexes': [i],
                                              'show_legend':False})
        figure8.plot_at_time()
        
        matplotlib.pylab.show() 
        
    def visualize_animation(self):
        '''
        Visualize the results of the simulation
        '''
        
        import Visualization.SimulAnimation as SimulAnimation
        import Visualization.AxesNeuronPropertyLine as AxesNeuronPropertyLine
        import Visualization.AxesPatternLine as AxesPatternLine
        import Visualization.AxesRasterPlot as AxesRasterPlot
        import Visualization.AxesWeightEvolutionLine as AxesWeightEvolutionLine
        import Visualization.AxesWeightHistogram as AxesWeightHistogram
        import Visualization.AxesWeightActivationPlot as AxesWeightActivationPlot
        import Visualization.AxesFiringOffset as AxesFiringOffset
        import Visualization.AxesActivationFiringOffset as AxesActivationFiringOffset
        import matplotlib.pylab
        
                
        # Adjust the frame_rate depending on whether the simulation is running at the same time
        if self.config_options['simulation']['run_simulation']:
            frame_rate = 0.1
        else:
            frame_rate = 0.1
        
        animation = SimulAnimation.SimulAnimation(simulation=self,numRows=3,numColumns=4,end_time=self.simulation_time,frame_rate=frame_rate,figsize=[23,14],dpi=80)
#         animation.add_subplot(fig_position=1,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'property':'Vm',
#                                               'layer':'goclayer',
#                                               'visible_data_only':True,
#                                               'show_legend':True,
#                                               'x_length': 1.})
#         animation.add_subplot(fig_position=2,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'property':'Gexc',
#                                               'layer':'goclayer',
#                                               'visible_data_only':True,
#                                               'show_legend':True,
#                                               'x_length': 1.})
#         animation.add_subplot(fig_position=3,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'property':'Vm',
#                                               'layer':'mflayer',
#                                               'cell_index': range(5),
#                                               'visible_data_only':True,
#                                               'show_legend':True,
#                                               'x_length': 1.})
        animation.add_subplot(fig_position=1,axes_type=AxesRasterPlot.AxesRasterPlot,
                            axes_parameters= {'data_provider':self.cerebellum,
                                              'pattern_provider':self.pattern_generator,
                                              'layer':'mflayer',
                                              'cell_index': range(50),
                                              'visible_data_only':True,
                                              'show_legend':False,
                                              'x_length':1.})
        animation.add_subplot(fig_position=2,axes_type=AxesRasterPlot.AxesRasterPlot,
                            axes_parameters= {'data_provider':self.cerebellum,
                                              'layer':'goclayer',
                                              'visible_data_only':True,
                                              'show_legend':False,
                                              'x_length':1.})
#         animation.add_subplot(fig_position=3,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'property':'Gexc',
#                                               'layer':'goclayer',
#                                               'visible_data_only':True,
#                                               'show_legend':False,
#                                               'x_length': 1.})
#         animation.add_subplot(fig_position=4,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'property':'Ginh',
#                                               'layer':'goclayer',
#                                               'visible_data_only':True,
#                                               'show_legend':False,
#                                               'x_length': 1.})
        animation.add_subplot(fig_position=3,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
                            axes_parameters= {'data_provider':self.cerebellum,
                                              'property':'Vm',
                                              'layer':'goclayer',
                                              'visible_data_only':True,
                                              'show_legend':False,
                                              'x_length': 1.,
                                              'y_window_lim': [-80e-3,-30e-3]})
        animation.add_subplot(fig_position=4,axes_type=AxesPatternLine.AxesPatternLine,
                            axes_parameters= {'pattern_provider':self.pattern_generator,
                                              'visible_data_only':True,
                                              'show_legend':False,
                                              'x_length':1.})
        animation.add_subplot(fig_position=5,axes_type=AxesFiringOffset.AxesFiringOffset,
                            axes_parameters= {'data_provider':self.cerebellum,
                                              'oscillation_freq':self.config_options['oscillations']['frequency'],
                                              'layer':'goclayer',
                                              'visible_data_only':True})
#         animation.add_subplot(fig_position=6,axes_type=AxesActivationFiringOffset.AxesActivationFiringOffset,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'pattern_provider':self.pattern_generator,
#                                               'oscillation_freq':self.config_options['oscillations']['frequency'],
#                                               'layer':'mflayer',
#                                               'visible_data_only':True})
#         animation.add_subplot(fig_position=6,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'property':'rC',
#                                               'layer':'goclayer',
#                                               'visible_data_only':True,
#                                               'show_legend':False})
#         animation.add_subplot(fig_position=7,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'property':'gL',
#                                               'layer':'goclayer',
#                                               'visible_data_only':True,
#                                               'show_legend':False})
        animation.add_subplot(fig_position=8,axes_type=AxesWeightEvolutionLine.AxesWeightEvolutionLine,
                            axes_parameters= {'data_provider':self.cerebellum,
                                              'layer':'gocgocsynapsis',
                                              'visible_data_only':True,
                                              'show_legend':False})
        for i in range(4):
            animation.add_subplot(fig_position=i+9,axes_type=AxesWeightActivationPlot.AxesWeightActivationPlot,
                            axes_parameters= {'data_provider':self.cerebellum,
                                              'pattern_provider': self.pattern_generator,
                                              'layer':'mfgocsynapsis',
                                              'target_indexes': [i],
                                              'show_legend':False})
#             animation.add_subplot(fig_position=2*i+8,axes_type=AxesWeightEvolutionLine.AxesWeightEvolutionLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'layer':'mfgocsynapsis',
#                                               'source_indexes': range(100),
#                                               'target_indexes': [i],
#                                               'visible_data_only':True,
#                                               'show_legend':False})
  
        
        
#         animation.add_subplot(fig_position=6,axes_type=AxesWeightHistogram.AxesWeightHistogram,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'layer':'mfgocsynapsis',
#                                               'num_bins': 60})

#         animation.add_subplot(fig_position=9,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'property':'firing_probability',
#                                               'layer':'goclayer',
#                                               'cell_index': [0],
#                                               'visible_data_only':True,
#                                               'show_legend':False,
#                                               'x_length': 1.})
        
#         animation.add_subplot(fig_position=9,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'property':'gain',
#                                               'layer':'goclayer',
#                                               'visible_data_only':True,
#                                               'show_legend':False,
#                                               'x_length': 1.})
#         animation.add_subplot(fig_position=10,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'property':'Vth',
#                                               'layer':'goclayer',
#                                               'visible_data_only':True,
#                                               'show_legend':False})
#         animation.add_subplot(fig_position=11,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'property':'r0',
#                                               'layer':'goclayer',
#                                               'visible_data_only':True,
#                                               'show_legend':False})
#         animation.add_subplot(fig_position=12,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'property':'ualpha',
#                                               'layer':'goclayer',
#                                               'visible_data_only':True,
#                                               'show_legend':False})
#         animation.add_subplot(fig_position=11,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'property':'refractoriness',
#                                               'layer':'goclayer',
#                                               'visible_data_only':True,
#                                               'show_legend':False})
        


        
                
   
        matplotlib.pylab.show() 
            
    def analyze_results(self):
        '''
        Analyze the estimators that have been set in the configuration file
        '''
        
        if self.config_options['simulation']['use_mpi']:
            import Analysis.MutualInformation as MutualInformation
        else:
            import Analysis.MutualInformationNoMPI as MutualInformation
        
        # Extract every mutual information to explore
        parameter_keys = [key for key in self.config_options.keys() if key.startswith('mutual_information')]
        mutual_information = []
        for key in parameter_keys:
            
            if not 'layer' in self.config_options[key]:
                logger.error('Layer name has not been specified in the mutual information section')
                raise Exception('NonSpecifiedLayer')
            
            if not 'window_length' in self.config_options[key]:
                logger.error('Window length has not been specified in the mutual information section')
                raise Exception('NonSpecifiedWindowLenght')
            
            if not 'time_bin' in self.config_options[key]:
                logger.error('time bin has not been specified in the mutual information section')
                raise Exception('NonSpecifiedTimeBin')
            
            if not 'record_to_file' in self.config_options[key]:
                self.config_options[key]['record_to_file'] = False
            
            
            logger.info('Analyzing mutual information in section %s',key)
            MIAnalysis = MutualInformation.MutualInformation(data_provider=self.cerebellum, pattern_generator=self.pattern_generator, layer=self.config_options[key]['layer'],
                                                             window_length=self.config_options[key]['window_length'], time_bin = self.config_options[key]['time_bin'])
            MIAnalysis.initialize()
            MIAnalysis.runAtTime(self.simulation_time)
            mutual_information.append(MIAnalysis.mutual_information/MIAnalysis.max_mutual_information)
            if self.config_options[key]['record_to_file']:
                filename = self.config_options['simulation']['data_path'] + '/' + self.config_options['simulation']['simulation_name'] + '/' + key 
                logger.debug('Writing mutual information from section %s to file %s',key,filename)
                MIAnalysis.writeToFile(file_name=filename)
                
        return mutual_information
    
 