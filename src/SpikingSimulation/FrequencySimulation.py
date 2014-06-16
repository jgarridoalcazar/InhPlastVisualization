'''
Created on May 27, 2014

@author: Jesus Garrido (jgarridoalcazar at gmail.com)
'''
import SpikingCerebellum.NestCerebellarModel as NestGenerator
from mpi4py import MPI
import Stimulation.FrequencyPatternGenerator as FrequencyPatternGenerator
import matplotlib.pylab
import numpy
import ConfigParser
import bisect
from Utils.Utils import ConfigSectionMap

import Visualization.SimulFigure as SimulFigure
import Visualization.SimulAnimation as SimulAnimation
import Visualization.AxesNeuronPropertyLine as AxesNeuronPropertyLine
import Visualization.AxesRasterPlot as AxesRasterPlot
import Visualization.AxesWeightEvolutionLine as AxesWeightEvolutionLine
import Visualization.AxesWeightHistogram as AxesWeightHistogram
import Visualization.AxesWeightActivationPlot as AxesWeightActivationPlot

class FrequencySimulation(object):
    '''
    This class defines a simulation where the parameters are taking from the
    configuration file passed as a parameter.
    '''

    def __init__(self,**kwargs):
        '''
        Constructor of the class. It creates a new simulation object.
        @param config_file Name of the file with the options of the model.
        '''
        
        if ('config_file' in kwargs):
            self.config_file = kwargs.pop('config_file')                         
        else:
            print 'Non-specified simulation config file.'
            raise Exception('Non-DefinedSimulationConfig')
        
        super(FrequencySimulation, self).__init__(**kwargs)
        
        return
    
    def _read_config_file_(self):
        '''
        Read all the sections in the configuration file.
        '''
        print 'Parsing configuration file ', self.config_file
    
        self.config_parser = ConfigParser.ConfigParser()
        self.config_parser.read(self.config_file)
        
        # Read every section in the file
        sections = self.config_parser.sections()
        self.config_options = dict()
        for sect in sections:
            self.config_options[sect] = ConfigSectionMap(config_parser = self.config_parser, section = sect)
        
        return
    
    def initialize(self):
        '''
        Initialize all the objects needed for running the simulation.
        '''
        comm = MPI.COMM_WORLD
        
        self._read_config_file_()
        
        # Read simulation general options
        if 'simulation' in self.config_options:
            if 'time' in self.config_options['simulation']:
                self.simulation_time = self.config_options['simulation']['time']
            else:
                self.simulation_time = 1
        else:
            self.simulation_time = 1
            
        print 'Simulation time fixed to',self.simulation_time,'s'
        
        # Initialize cerebellar model
        print 'Process:', comm.Get_rank(),'. Creating cerebellum generator'
        self.cerebellum = NestGenerator.NestCerebellarModel(config_dict=self.config_options)
    
        print 'Process:', comm.Get_rank(),'. Initializing cerebellum generator'
        self.cerebellum.initialize_simulation()
    
        print 'Process:', comm.Get_rank(),'. Building the network'
        self.cerebellum.build_network()
        
        
        # Initialize oscillatory input current
        if 'oscillations' in self.config_options:
            print 'Process:', comm.Get_rank(),'. Creating AC Current generator'
            self.cerebellum.add_ac_current(**self.config_options['oscillations'])
            
        # Initialize frequency stimulation input current
        if 'stimulation' in self.config_options:      
            print 'Process:', comm.Get_rank(),'. Creating DC Current generator'
            self.config_options['stimulation']['simulation_time'] = self.simulation_time
            self.config_options['stimulation']['number_of_fibers'] = self.cerebellum.mflayer.number_of_neurons
            self.pattern_generator = FrequencyPatternGenerator.FrequencyPatternGenerator(**self.config_options['stimulation'])
            self.pattern_generator.initialize()
            
            self.pattern_length, self.pattern_activations = self.pattern_generator.get_all_patterns()
            self.pattern_length_cum = self.pattern_generator.pattern_length_cum
            
        self.current_time = 0.
            
        
    def run_simulation(self, **kwargs):
        '''
        Run the simulation according to the configuration file.
        @param end_time Time until when simulation will be run
        '''
        
        comm = MPI.COMM_WORLD
        
        if 'end_time' in kwargs:
            end_time = kwargs.pop('end_time')
            
            if end_time>self.simulation_time:
                print 'Warning: simulation time is shorter than end_time. Simulating', self.simulation_time,'seconds'
                
            end_time = min(end_time,self.simulation_time)
        else:
            end_time = self.simulation_time
            
        if 'stimulation' in self.config_options:
            init_index = bisect.bisect_left(self.pattern_length_cum, self.current_time)
            end_index = bisect.bisect_left(self.pattern_length_cum,end_time)
        
            for index in range(init_index,end_index+1):
                sim_time = min(self.pattern_length_cum[index]-self.current_time,end_time-self.current_time)
                
                self.cerebellum.set_dc_current(amplitude=self.pattern_activations[index])
            
                print 'Process:', comm.Get_rank(),'. Running the simulation',sim_time,'seconds until', self.cerebellum.simulation_time+sim_time,' seconds.'
                self.cerebellum.simulate_network(sim_time)
                
                self.current_time = self.cerebellum.simulation_time
        else:
            sim_time = self.simulation_time-self.current_time
            print 'Process:', comm.Get_rank(),'. Running the simulation',sim_time,'seconds until', self.cerebellum.simulation_time+sim_time,' seconds.'
            self.cerebellum.simulate_network(sim_time)
            self.current_time = self.cerebellum.simulation_time
            
    def visualize_results(self):
        '''
        Visualize the results of the simulation
        '''
                
#         figure7 = SimulFigure.SimulFigure(simulation = self, numRows=4,numColumns=2,figsize=[23,14],dpi=80)
#         figure7.add_subplot(fig_position=1,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'property':'Vm',
#                                               'layer':'goclayer',
#                                               'visible_data_only':True,
#                                               'show_legend':True})
#         figure7.add_subplot(fig_position=2,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'property':'Gexc',
#                                               'layer':'goclayer',
#                                               'visible_data_only':True,
#                                               'show_legend':True})
#         figure7.add_subplot(fig_position=3,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'property':'Vm',
#                                               'layer':'mflayer',
#                                               'cell_index': range(5),
#                                               'visible_data_only':True,
#                                               'show_legend':True})
#         figure7.add_subplot(fig_position=4,axes_type=AxesRasterPlot.AxesRasterPlot,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'pattern_provider':self.pattern_generator,
#                                               'layer':'mflayer',
#                                               'cell_index': range(100),
#                                               'visible_data_only':True,
#                                               'show_legend':True,
#                                               'x_length':10})
#         figure7.add_subplot(fig_position=5,axes_type=AxesWeightEvolutionLine.AxesWeightEvolutionLine,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'layer':'mfgocsynapsis',
#                                               'source_indexes': range(5),
#                                               'target_indexes': range(1),
#                                               'visible_data_only':True,
#                                               'show_legend':True})
#         figure7.add_subplot(fig_position=6,axes_type=AxesWeightHistogram.AxesWeightHistogram,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'layer':'mfgocsynapsis',
#                                               'num_bins': 60})
#         figure7.add_subplot(fig_position=7,axes_type=AxesWeightActivationPlot.AxesWeightActivationPlot,
#                             axes_parameters= {'data_provider':self.cerebellum,
#                                               'pattern_provider': self.pattern_generator,
#                                               'layer':'mfgocsynapsis'})
#         figure7.plot_at_time()
        
        animation = SimulAnimation.SimulAnimation(simulation=self,numRows=5,numColumns=1,end_time=self.simulation_time,frame_rate=0.1,figsize=[23,14],dpi=80)
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
        animation.add_subplot(fig_position=1,axes_type=AxesNeuronPropertyLine.AxesNeuronPropertyLine,
                            axes_parameters= {'data_provider':self.cerebellum,
                                              'property':'Vm',
                                              'layer':'goclayer',
                                              'visible_data_only':True,
                                              'show_legend':False,
                                              'x_length': 1.})
        animation.add_subplot(fig_position=2,axes_type=AxesRasterPlot.AxesRasterPlot,
                            axes_parameters= {'data_provider':self.cerebellum,
                                              'pattern_provider':self.pattern_generator,
                                              'layer':'mflayer',
                                              'cell_index': range(50),
                                              'visible_data_only':True,
                                              'show_legend':True,
                                              'x_length':1.})
        animation.add_subplot(fig_position=3,axes_type=AxesWeightEvolutionLine.AxesWeightEvolutionLine,
                            axes_parameters= {'data_provider':self.cerebellum,
                                              'layer':'mfgocsynapsis',
                                              'source_indexes': range(100),
                                              'target_indexes': range(1),
                                              'visible_data_only':True,
                                              'show_legend':False})
        animation.add_subplot(fig_position=4,axes_type=AxesWeightHistogram.AxesWeightHistogram,
                            axes_parameters= {'data_provider':self.cerebellum,
                                              'layer':'mfgocsynapsis',
                                              'num_bins': 60})
        animation.add_subplot(fig_position=5,axes_type=AxesWeightActivationPlot.AxesWeightActivationPlot,
                            axes_parameters= {'data_provider':self.cerebellum,
                                              'pattern_provider': self.pattern_generator,
                                              'layer':'mfgocsynapsis'})
          

        matplotlib.pylab.show() 
            
            
            
    
        
        
    