'''
Created on December 10, 2015

@author: Jesus Garrido (jgarridoalcazar at gmail.com)
'''
import numpy
import time
import logging
import scipy.spatial.distance
import SynapticLayerNoMPI
from numpy import rank

logger = logging.getLogger('Simulation')

class SynapticLayer(SynapticLayerNoMPI.SynapticLayer):
    '''
    This class defines a synaptic layer and the data needed to generate it in a simulator.
    If several MPI processes are used, each process store only the local connections (those
    reaching a local cell). It inherits all the methods from the NoMPI implementation (since
    the target neurons are assumed to be distributed between parallel processes) except the 
    _generate_glomeruli_n2one_local_connections_ method which needs to be specifically reimplemented.
    '''
    
    def _generate_glomeruli_n2one_local_connections_(self):
        '''
        Generate connections between source and target layers n2one with random selection of source cells (sr0,sr1,...,sr(n-1) to t0, srn, s(rn+1),...,s(r2n-1) to t1,...).
        In this case the probability of connections is generated according to an exponential distribution with average length the average_dendritic_length parameter.
        '''
        
        if self.intermediate_layer is None:
            logger.error('Non-specified intermediate layer. It is required in glomeruli-type connectivity creation.')
            raise Exception('Non-DefinedProperty')
        
        if self.intermediate_to_target_synaptic_layer is None:
            logger.error('Non-specified intermediate_to_target synaptic layer. It is required in glomeruli-type connectivity creation.')
            raise Exception('Non-DefinedProperty')
        
        
        # It is assumed that every neuron layer stores the coordinates of all the neurons (not only the local ones).
        interm_coord = self.intermediate_layer.get_absolute_coordinates()
        source_coord = self.source_layer.get_absolute_coordinates()
        
        n_source_cells = source_coord.shape[0]
        n_interm_cells = interm_coord.shape[0]
        
        lambda_val = 1./self.connectivity_parameters['average_dendritic_length']
                
        # Generate the look-up table of target grcs and neighbouring glomeruli
        target_grcs_list = []
        neighbouring_glomeruli_list = []
        
        logger.debug('Creating glomerulus neighbourhood in layer %s. %s glomeruli', self.__name__, n_interm_cells)
        
        # ------------------------------------------------------------------------------------------------------------------
        # Collect all the intermediate_to_target_synaptic_layer (source_index and target_index in all the MPI processes)
        from mpi4py import MPI
        
        comm = MPI.COMM_WORLD
        
        comm.Barrier()
        
        nprocs = comm.Get_size()
        rank = comm.Get_rank()
        
        # Send the number of elements
        lsyn = numpy.array([self.intermediate_to_target_synaptic_layer.number_of_synapses], dtype=numpy.uint64)
        num_synapses = numpy.empty(comm.Get_size(), dtype=numpy.uint64)
        comm.AllGather([lsyn, MPI.UNSIGNED_LONG], [num_synapses, MPI.UNSIGNED_LONG], root=0) 
        
        print 'Process', rank,':','Sent number:',lsyn,'Collected numbers ->',num_synapses
        
        # Total number of synapses to be gathered
        intermediate_to_target_synaptic_count = num_synapses.sum()
        offset = numpy.roll(numpy.cumsum(num_synapses),1)
        
        intermediate_to_target_source_indexes = numpy.empty(intermediate_to_target_synaptic_count, dtype=numpy.uint32)
        intermediate_to_target_target_indexes = numpy.empty(intermediate_to_target_synaptic_count, dtype=numpy.uint32)
        
        # Gather the time and neuron_id arrays
        comm.Allgatherv(self.intermediate_to_target_synaptic_layer.source_indexes, [intermediate_to_target_source_indexes, num_synapses, offset, MPI.UNSIGNED_INT])
        
        print 'Process', rank,':','Sent number:',self.intermediate_to_target_synaptic_layer.source_indexes,'Collected numbers ->',intermediate_to_target_source_indexes
        
        comm.Allgatherv(self.intermediate_to_target_synaptic_layer.target_indexes, [intermediate_to_target_target_indexes, num_synapses, offset, MPI.UNSIGNED_INT])
        
        print 'Process', rank,':','Sent number:',self.intermediate_to_target_synaptic_layer.target_indexes,'Collected numbers ->',intermediate_to_target_target_indexes
        
        # intermediate_to_target_source_indexes -> Array with all the source indexes in the connections Glomeruli -> GrC
        # intermediate_to_target_target_indexes -> Array with all the target indexes in the connections Glomeruli -> GrC
        # intermediate_to_target_synaptic_count -> Number of synapses in the connections Glomeruli -> GrC
        
        # Sort the intermediate_to_target_synaptic_layer (target_indexes) according to the source layer
        target_order = intermediate_to_target_source_indexes.argsort()
        unique_source,split_indices = numpy.unique(intermediate_to_target_source_indexes[target_order], return_index=True)
        target_ordered = numpy.split(intermediate_to_target_target_indexes[target_order], split_indices[1:])
        
        target_all_ordered = []
        source_index = 0
        for i in range(self.intermediate_layer.number_of_neurons):
            if unique_source[source_index]==i:
                target_all_ordered.append(target_ordered[source_index])
                source_index = source_index + 1
            else:
                target_all_ordered.append(numpy.empty([0]))
                
        target_ordered = numpy.array(target_all_ordered)
        # target_ordered[ind] returns the target elements of the source glomeruli ind
                
        
        # Sort the intermediate_to_target_synaptic_layer (source_indexes) according to the target layer
        source_order = intermediate_to_target_target_indexes.argsort()
        unique_target,split_indices = numpy.unique(intermediate_to_target_target_indexes[source_order], return_index=True)
        source_ordered = numpy.split(intermediate_to_target_source_indexes[source_order], split_indices[1:])
        
        source_all_ordered = []
        target_index = 0
        for i in range(self.target_layer.number_of_neurons):
            if unique_target[target_index]==i:
                source_all_ordered.append(source_ordered[target_index])
                target_index = target_index + 1
            else:
                source_all_ordered.append(numpy.empty([0]))
                
        source_ordered = numpy.array(source_all_ordered)
        # source_ordered[ind] returns the source glomeruli of the target neuron ind
        
        # -------------------------------------------------------------------------------------------------------
        # This section has to be executed only in the first process when run with MPI
        # Number of interm cells in each memory block
        if (rank==0):
            glomerulus_connection_index = []
            goc_connection_index = []
            
            interm_layers_in_block = self.memory_block / self.source_layer.number_of_neurons
            
            # Distribute the intermediate cells in the number of blocks
            interm_indexes = range(0,n_interm_cells, interm_layers_in_block)
            
            for i,_ in enumerate(interm_indexes):
                if i==(len(interm_indexes)-1):
                    glomerulus_indexes_block = range(interm_indexes[i], n_interm_cells)
                else:
                    glomerulus_indexes_block = range(interm_indexes[i], interm_indexes[i+1])
                    
                    
                distance = scipy.spatial.distance.cdist(source_coord, interm_coord[glomerulus_indexes_block], 'euclidean')
            
                # Generate the probability of connection for each pair of source/target neurons following an exponential distribution with average_dendritic_length
                probability = lambda_val * numpy.exp(-lambda_val*distance)
            
                # Generate random numbers according to uniform distribution and compare with the probability of each connection
                rand_nums = self.random_generator.uniform(0,1,probability.shape)
                
                logger.debug('Calculating  neighbourhood in layer %s and block %s.', self.__name__, i+1)
                
                # Generate the neighbouring glomeruli
                if glomerulus_connection_index:
                    target_cells = numpy.unique(numpy.hstack(target_ordered[glomerulus_connection_index]))
                else:
                    target_cells = numpy.array([],dtype=numpy.uint32)
                neighbouring_glomeruli = source_ordered[target_cells]
                
                for index_glom, glomeruli_array in enumerate(neighbouring_glomeruli):
                    current_glomeruli = numpy.logical_and(glomeruli_array>=glomerulus_indexes_block[0], glomeruli_array<=glomerulus_indexes_block[-1])
                    
                    adjusted_neighbour_indexes = numpy.array(glomeruli_array[current_glomeruli]) - glomerulus_indexes_block[0]
                    
                    probability[goc_connection_index[index_glom],adjusted_neighbour_indexes] = 0.0
                
                
                logger.debug('Calculating  glomerulus connections in layer %s and block %s. %s glomeruli', self.__name__, i+1, len(glomerulus_indexes_block))
                # For each target neuron normalize the probability and multiply by the average number of source cells
                # When parallelized this part has to be done sequentially in the first process
                for glomerulus_index in glomerulus_indexes_block:
                    local_index = glomerulus_index-glomerulus_indexes_block[0]
                    
                    norm_probability = probability[:,local_index] / numpy.sum(probability[:,local_index]) * self.connectivity_parameters['average_number_of_source_cells_per_glomerulus']
                
                    connected_goc_indexes, = numpy.where(rand_nums[:,local_index]<norm_probability[:])
                
                    # Store the GoC-Glomerulus connections
                    glomerulus_connection_index.extend([glomerulus_index]*connected_goc_indexes.shape[0])
                    goc_connection_index.extend(connected_goc_indexes)
                
                    neighbour_glomeruli = numpy.unique(numpy.hstack(source_ordered[target_ordered[glomerulus_index]]))
                    current_glomeruli = numpy.logical_and(neighbour_glomeruli>=glomerulus_indexes_block[0],neighbour_glomeruli<=glomerulus_indexes_block[-1])
                    
                    adjusted_neighbour_indexes = numpy.array(neighbour_glomeruli[current_glomeruli]) - glomerulus_indexes_block[0]
                    
                    # Update the probability of the following glomeruli to avoid repeating connections
                    goc_indexes = numpy.tile(connected_goc_indexes,adjusted_neighbour_indexes.shape[0])
                    glomeruli_rep_indexes = numpy.repeat(adjusted_neighbour_indexes,connected_goc_indexes.shape[0])
                    probability[goc_indexes,glomeruli_rep_indexes] = 0.0
                    
                logger.debug('Generated glomerulus-like connections in layer %s (%s of %s). Local: %s', self.__name__,i+1,len(interm_indexes),len(goc_connection_index))
            
            
        # ------------------------------------------------------------------------------------------------------------------
        # Send the GoC-Glomerulus connections lists to all the processes (each one will take only those connections reaching their local neurons)
        
        # Send the number of connections
        if rank==0:    
            ncon = numpy.array([len(glomerulus_connection_index)], dtype=numpy.uint64)
        else:
            ncon = numpy.empty(1, dtyp = numpy.uint64)
        
        comm.Bcast([ncon, MPI.UNSIGNED_LONG], root=0) 
            
        print 'Process', rank,':','Broadcasting number:',ncon
        
        if rank==0:
            glomerulus_connection_index = numpy.array(glomerulus_connection_index)
            goc_connection_index = numpy.array(goc_connection_index)
        else:
            glomerulus_connection_index = numpy.empty(ncon, dtype=numpy.uint32)
            goc_connection_index = numpy.array(ncon, dtype=numpy.uint32)
        
        # Broadcast the glomerulus and goc indexes
        comm.Bcast([glomerulus_connection_index, MPI.UNSIGNED_INT], root=0)
            
        print 'Process', rank,':','Broadcasting list of glomerulus indexes:',glomerulus_connection_index
            
        comm.Bcast([goc_connection_index, MPI.UNSIGNED_INT], root=0)
        
        print 'Process', rank,':','Broadcasting list of goc indexes:',goc_connection_index
            
        # Generate the GoC-GrC connections from the GoC-Glomerulus connections
        self.source_index = []
        self.target_index = []
        
        # Generate the target GrCs of the selected glomeruli
        target_grcs = target_ordered[glomerulus_connection_index]
        local_target_cells, = numpy.where(self.target_layer.is_local_node)
        for i, target_group in enumerate(target_grcs):
            # Select the target grcs which are local to the current processor
            selected_local_grcs, = numpy.where(numpy.in1d(target_group, local_target_cells)) 
            self.source_index.extend(numpy.repeat(goc_connection_index[i],len(selected_local_grcs)))
            self.target_index.extend(selected_local_grcs)
        
        logger.debug('Generated glomerulus-like connections in layer %s. Local: %s', self.__name__,len(self.source_index))
        
        self.source_index = numpy.array(self.source_index, dtype=numpy.uint32)
        self.target_index = numpy.array(self.target_index, dtype=numpy.uint32)
        self.number_of_synapses = self.target_index.shape[0]
        return 
