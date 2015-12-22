'''
Created on May 12, 2014

@author: Jesus Garrido (jgarridoalcazar at gmail.com)
'''

import abc
from CerebellarModel import CerebellarModel
import shutil
import sys
import numpy
import os
import time
import logging
import NestCerebellarModelNoMPI
import nest

logger = logging.getLogger('Simulation')

class NestCerebellarModel(NestCerebellarModelNoMPI.NestCerebellarModel):
    '''
    This class defines an inherited class including all the methods
    needed in order to generate a cerebellum-like network for EDLUT
    simulator.
    '''
    __metaclass__ = abc.ABCMeta
    
        
    def _synchronize_processes_(self):
        from mpi4py import MPI
        
        comm = MPI.COMM_WORLD

        comm.Barrier()
        return
        
    
    def _create_nodes(self):
        '''
        Generate the nodes of every neuron layer in the model
        '''
        
        super(NestCerebellarModel, self)._create_nodes()
        
        from mpi4py import MPI
        
        comm = MPI.COMM_WORLD

        for layer in self.neuron_layers:
            # Collect the minimum of the layer indexes
            localMin = numpy.array(numpy.min(layer.nest_layer), dtype='l') 
            layer.MinIndex = numpy.array(0, dtype='l') 
            comm.Allreduce([localMin, MPI.LONG], [layer.MinIndex, MPI.LONG], op=MPI.MIN)
            
            self.config_dict[layer.__name__]['minindex'] = layer.MinIndex
            
        return
    
    def get_spike_activity(self, **kwargs):
        '''
        Get the spikes been fired in the cells, layer and time specified in the parameters.
        @param neuron_layer Layer name as the section name in the config file (.cfg).
        @param neuron_indexes Indexes of the cell to get the activity. Relative indexes in the layer.
        @param init_time Initial time from which retrieve the activity.
        @param end_time Final time until which retrieve the activity.  
        '''
        
        # Collect all the parameters
        if 'neuron_layer' in kwargs:
            neuron_layer_name = kwargs.pop('neuron_layer')
        else:
            logger.error('Non specified neuron layer in get_spike_activity function')
            raise Exception('Non-DefinedNeuronLayer')
        
        
        if neuron_layer_name in self.layer_map:
            neuron_layer = self.layer_map[neuron_layer_name]
        else:
            logger.error('Invalid neuron layer in get_spike_activity function')
            raise Exception('InvalidNeuronLayer')
        
        if 'neuron_indexes' in kwargs:
            local_indexes = kwargs.pop('neuron_indexes')
            
            if (max(local_indexes)>=(neuron_layer.MinIndex+neuron_layer.number_of_neurons)):
                logger.error('Invalid neuron index in get_spike_activity function')
                raise Exception('InvalidNeuronIndex')
            
            neuron_indexes = [neuron_layer.MinIndex+index for index in local_indexes]            
        else:
            neuron_indexes = neuron_layer.nest_layer
            
        if 'init_time' in kwargs:
            init_time = kwargs.pop('init_time')
        else:
            init_time = 0
            
        if 'end_time' in kwargs:
            end_time = kwargs.pop('end_time')
        else:
            end_time = float('inf')
        
        # Comprobar que la neurona (spike recorder) sea local
        if not neuron_layer.nest_spike_detector: 
            logger.error('Invalid neuron layer in get_spike_activity function. The activity in this layer has not been recorded')
            raise Exception('InvalidNeuronLayer')
        
        # If the spike detector is local to this process
        # logger.debug('Checking locality in layer %s from %s: %s',neuron_layer.__name__, neuron_layer.nest_spike_detector, nest.GetStatus(neuron_layer.nest_spike_detector,'local'))
        # logger.debug('Getting spikes in layer %s from %s: %s',neuron_layer.__name__, neuron_layer.nest_spike_detector, nest.GetStatus(neuron_layer.nest_spike_detector,'events'))
        spike_events = nest.GetStatus(neuron_layer.nest_spike_detector.tolist(),'events')[0]
        time = spike_events['times']*1e-3
        neuron_id = spike_events['senders']
        
        # Select only those events in the specified interval
        index = (time>=init_time) & (time<=end_time) 
        time = time[index]
        neuron_id = neuron_id[index]
        
        # Select only those events in the specified cells
        index = numpy.in1d(neuron_id, neuron_indexes)
        time = time[index]
        neuron_id = neuron_id[index]
        
        # Send relative neuron_id
        neuron_id = neuron_id-neuron_layer.MinIndex
        
        from mpi4py import MPI
        
        comm = MPI.COMM_WORLD
        
        # Send the number of elements
        lsum = numpy.array([len(time)], dtype=numpy.uint32)
        if self.get_my_process_id()==0:
            num_spikes = numpy.empty(comm.Get_size(), dtype=numpy.uint32)
        else:
            num_spikes = None
        comm.Gather([lsum, MPI.UNSIGNED_INT], [num_spikes, MPI.UNSIGNED_INT], root=0) 
        
        # print 'Process',self.get_my_process_id(),':','Sent number:',lsum,'Collected numbers ->',num_spikes
        
        if self.get_my_process_id()==0:
            # Total number of spikes to be gathered
            gsum = num_spikes.sum()
            num_sent = tuple(num_spikes)
            aux_array = numpy.roll(num_spikes.cumsum(),1)
            aux_array[0] = 0
            offset = tuple(aux_array)
            
            gtime = numpy.empty(gsum, dtype=numpy.float64)
            gneuron_id = numpy.empty(gsum, dtype=numpy.uint64)
        else:
            num_sent = None
            offset = None
            gtime = None
            gneuron_id = None
        
        # Gather the time and neuron_id arrays
        comm.Gatherv(time, [gtime, num_sent, offset, MPI.DOUBLE], root=0)
        comm.Gatherv(neuron_id, [gneuron_id, num_sent, offset, MPI.UNSIGNED_LONG], root=0)
        
        # print 'Process',self.get_my_process_id(),':','Collected time ->',gtime,'Collected Neurons ->',gneuron_id
        return (gtime,gneuron_id)
    
    def get_state_variable(self, **kwargs):
        '''
        Get the state variables in the cells, layer and during the time specified in the parameters.
        @param neuron_layer Layer name as the section name in the config file (.cfg).
        @param neuron_indexes Indexes of the cell to get the activity.
        @param init_time Initial time from which retrieve the activity.
        @param end_time Final time until which retrieve the activity.
        @param state_variable Name of the state variable to retrieve.   
        '''
        # Collect all the parameters
        if 'neuron_layer' in kwargs:
            neuron_layer_name = kwargs.pop('neuron_layer')
        else:
            logger.error('Non specified neuron layer in get_state_variable function')
            raise Exception('Non-DefinedNeuronLayer')
        
        
        if neuron_layer_name in self.layer_map:
            neuron_layer = self.layer_map[neuron_layer_name]
        else:
            logger.error('Invalid neuron layer in get_state_variable function')
            raise Exception('InvalidNeuronLayer')
        
        if 'neuron_indexes' in kwargs:
            local_indexes = kwargs.pop('neuron_indexes')
            
            if (max(local_indexes)>=(neuron_layer.MinIndex+neuron_layer.number_of_neurons)):
                logger.error('Invalid neuron index in get_state_variable function')
                raise Exception('InvalidNeuronIndex')
            
            neuron_indexes = [neuron_layer.MinIndex+index for index in local_indexes]         
        else:
            neuron_indexes = neuron_layer.nest_layer
            
        if 'state_variable' in kwargs:
            variable_name = kwargs.pop('state_variable')
        else:
            logger.error('Non specified state_variable in get_state_variable function')
            raise Exception('Non-DefinedStateVariable')
        
        if variable_name in self.stateTranslatorDict:
            state_variable_name = self.stateTranslatorDict[variable_name][0]
        else:
            logger.error('Invalid state variable in get_state_variable function')
            raise Exception('InvalidStateVariable')
        
        if not variable_name in neuron_layer.record_vars:
            logger.error('Invalid state variable in get_state_variable function. %s has not been recorded', variable_name)
            raise Exception('InvalidStateVariable')
        
        if 'init_time' in kwargs:
            init_time = kwargs.pop('init_time')
        else:
            init_time = 0
            
        if 'end_time' in kwargs:
            end_time = kwargs.pop('end_time')
        else:
            end_time = float('inf')
        
        # If the spike detector is local to this process
        recording_events = nest.GetStatus(neuron_layer.nest_multimeter.tolist(),'events')[0]
        # print 'Process',self.get_my_process_id(),': recording_events=',recording_events
        time = recording_events['times']*1e-3
        neuron_id = recording_events['senders']
        value = recording_events[state_variable_name]*self.stateTranslatorDict[variable_name][1]
        
        # Select only those events in the specified interval
        index = (time>=init_time) & (time<=end_time) 
        time = time[index]
        neuron_id = neuron_id[index]
        value = value[index]
        
        # Select only those events in the specified cells
        index = numpy.in1d(neuron_id, neuron_indexes)
        time = time[index]
        neuron_id = neuron_id[index]
        value = value[index]
        
        # Send relative neuron_id
        neuron_id = neuron_id-neuron_layer.MinIndex
        
        from mpi4py import MPI
        
        comm = MPI.COMM_WORLD
        
        # Send the number of elements
        lsum = numpy.array([len(time)], dtype=numpy.uint64)
        if self.get_my_process_id()==0:
            num_events = numpy.empty(comm.Get_size(), dtype=numpy.uint64)
        else:
            num_events = None
        comm.Gather([lsum, MPI.UNSIGNED_LONG], [num_events, MPI.UNSIGNED_LONG], root=0) 
        
        # print 'Process',self.get_my_process_id(),':','Sent number:',lsum,'Collected numbers ->',num_events
        
        if self.get_my_process_id()==0:
            # Total number of spikes to be gathered
            gsum = num_events.sum()
            num_sent = tuple(num_events)
            aux_array = numpy.roll(num_events.cumsum(),1)
            aux_array[0] = 0
            offset = tuple(aux_array)
            
            gtime = numpy.empty(gsum, dtype=numpy.float64)
            gneuron_id = numpy.empty(gsum, dtype=numpy.uint64)
            gvalue = numpy.empty(gsum, dtype=numpy.float64)
        else:
            num_sent = None
            offset = None
            gtime = None
            gneuron_id = None
            gvalue = None
        
        # Gather the time and neuron_id arrays
        comm.Gatherv(time, [gtime, num_sent, offset, MPI.DOUBLE], root=0)
        comm.Gatherv(neuron_id, [gneuron_id, num_sent, offset, MPI.UNSIGNED_LONG], root=0)
        comm.Gatherv(value, [gvalue, num_sent, offset, MPI.DOUBLE], root=0)
        
        # print 'Process',self.get_my_process_id(),':','Collected time ->',gtime,'Collected Neurons ->',gneuron_id,'Collected values ->',gvalue
        if gtime is not None:
            gtime = gtime
            
        return (gtime,gneuron_id,gvalue)
        
    
    def get_synaptic_weights(self, **kwargs):
        '''
        Get the synaptic weights in the synapses, layer and during the time specified in the parameters.
        @param synaptic_layer Layer name as the section name in the config file (.cfg).
        @param source_indexes Indexes of the source cells of the synapses to get the activity.
        @param target_indexes Indexes of the target cells of the synapses to get the activity.
        @param init_time Initial time from which retrieve the activity.
        @param end_time Final time until which retrieve the activity.
        '''
        # Collect all the parameters
        if 'synaptic_layer' in kwargs:
            synaptic_layer_name = kwargs.pop('synaptic_layer')
        else:
            logger.error('Non specified synaptic layer in get_synaptic_weights function')
            raise Exception('Non-DefinedSynapticLayer')
        
        if synaptic_layer_name in self.layer_map:
            synaptic_layer = self.layer_map[synaptic_layer_name]
        else:
            logger.error('Invalid synaptic layer in get_synaptic_weights function')
            raise Exception('InvalidSynapticLayer')
        
        if 'source_indexes' in kwargs:
            source_indexes = kwargs.pop('source_indexes')
            
            if (max(source_indexes)>=(synaptic_layer.source_layer.number_of_neurons)):
                logger.error('Invalid source index in get_synaptic_weights function')
                raise Exception('InvalidSourceIndex')
        else:
            source_indexes = range(synaptic_layer.source_layer.number_of_neurons)
        
        if 'target_indexes' in kwargs:
            target_indexes = kwargs.pop('target_indexes')
            
            if (max(target_indexes)>=(synaptic_layer.target_layer.number_of_neurons)):
                logger.error('Invalid target index in get_synaptic_weights function')
                raise Exception('InvalidTargetIndex')
        else:
            target_indexes = range(synaptic_layer.target_layer.number_of_neurons)
            
        # print 'Process',self.get_my_process_id(),':','Source indexes:',source_indexes,'Target indexes ->',target_indexes
        
        if 'init_time' in kwargs:
            init_time = kwargs.pop('init_time')
        else:
            init_time = 0
            
        if 'end_time' in kwargs:
            end_time = kwargs.pop('end_time')
        else:
            end_time = float('inf')
        
        
        if not synaptic_layer.weight_recording:
            logger.error('Invalid synaptic layer in get_synaptic_weights function. The weights in this layer have not been recorded')
            raise Exception('InvalidSynapticLayer')
        
        # print 'Process',self.get_my_process_id(),':','Weight record:',synaptic_layer.weight_record
        
        # Calculate selected connection indexes
        connections = synaptic_layer.weight_record['connections']
        # Optimize connection indexes selection
        if (source_indexes == range(synaptic_layer.source_layer.number_of_neurons)):
            if (target_indexes == range(synaptic_layer.target_layer.number_of_neurons)):
                # No selection -> All the connections in the layer
                connection_indexes = range(len(connections))
            else:
                # Target indexes selection
                connection_indexes = [index for index in range(len(connections)) if (connections[index][1] in target_indexes)]
        else:
            if (target_indexes == range(synaptic_layer.target_layer.number_of_neurons)):
                # Source indexes selection
                connection_indexes = [index for index in range(len(connections)) if (connections[index][0] in source_indexes)]
            else:
                # Source and target selection
                connection_indexes = [index for index in range(len(connections)) if (connections[index][0] in source_indexes) and (connections[index][1] in target_indexes)]
        selected_connections = connections[connection_indexes]
        
        #logger.debug('Selected connections: %s',selected_connections)
        
        # Calculate selected time indexes
        time = synaptic_layer.weight_record['time']
        time_indexes = (time>=init_time) & (time<=end_time) 
        selected_time = time[time_indexes]
        
        # Pick selected weights
        weights = synaptic_layer.weight_record['weights']
        selected_weights = numpy.array([record[connection_indexes] for record in weights[time_indexes]]).transpose()
        
        # Send the number of elements
        time_elements = len(selected_time) 
        
        from mpi4py import MPI
        
        comm = MPI.COMM_WORLD
        
#         # Gather every connection record individually
#         connection_aux = numpy.array(len(connection_indexes), dtype=numpy.uint64) 
#         connection_elements = numpy.array(0, dtype=numpy.uint64) 
#         comm.Reduce(connection_aux, connection_elements, op=MPI.SUM, root=0)
#         
#         logger.debug('Connection number sent: %s. Connection number collected: %s',connection_aux,connection_elements)
#         
        # Send the number of elements
        connection_aux = numpy.array([len(connection_indexes)], dtype=numpy.uint64)
        if self.get_my_process_id()==0:
            connection_elements = numpy.empty(comm.Get_size(), dtype=numpy.uint64)
        else:
            connection_elements = None
        comm.Gather([connection_aux, MPI.UNSIGNED_LONG], [connection_elements, MPI.UNSIGNED_LONG], root=0) 
        
        #logger.debug('Sent number: %s. Collected numbers: %s',connection_aux,connection_elements)
        
        if self.get_my_process_id()==0:
            # Total number of spikes to be gathered
            gsum = connection_elements.sum()
            con_num_sent = tuple(connection_elements*2)
            weight_num_sent = tuple(connection_elements*time_elements)
            aux_array = numpy.roll(connection_elements.cumsum(),1)
            aux_array[0] = 0
            con_offset = tuple(aux_array*2)
            weight_offset = tuple(aux_array*time_elements)
            
            gtime = selected_time
            gconnections = numpy.empty(gsum*2, dtype=numpy.uint32)
            gweights = numpy.empty((gsum, time_elements), dtype=numpy.float32)
        else:
            con_num_sent = None
            weight_num_sent = None
            con_offset = None
            weight_offset = None
            gtime = None
            gconnections = None
            gweights = None
        
        # Gather the time and neuron_id arrays
        #print 'Process',self.get_my_process_id(),':','Sent Connections ->',selected_connections.ravel(order='C'),'Sent weights ->',selected_weights
        #print 'Process',self.get_my_process_id(),':','Number Sent: ->',con_num_sent,'Sent offset: ->',con_offset
        #comm.Gatherv(time, [gtime, num_sent, offset, MPI.DOUBLE], root=0)
        comm.Gatherv(selected_connections.ravel(order='C'), [gconnections, con_num_sent, con_offset, MPI.UNSIGNED_INT], root=0)
        comm.Gatherv(selected_weights.ravel(order='C'), [gweights, weight_num_sent, weight_offset, MPI.FLOAT], root=0)
        
        #print 'Process',self.get_my_process_id(),':','Collected time ->',gtime,'Collected Connections ->',gconnections,'Collected weights ->',gweights
        
        if self.get_my_process_id()==0:
            gconnections = numpy.reshape(gconnections,(gsum,2),order='C')
            gweights = numpy.reshape(gweights,(gsum, time_elements), order='C')
        
        return (gtime,gconnections,gweights)
        
        
    def get_number_of_virtual_processes(self):
        '''
        Return the number of virtual processes. It might be used to decide the number of seeds to be generated.  
        '''
        import nest
        
        return nest.GetKernelStatus(['total_num_virtual_procs'])[0]
    
    def get_my_process_id(self):
        '''
        Return the id-number of this process.   
        '''
        import nest
        
        return nest.Rank()
    