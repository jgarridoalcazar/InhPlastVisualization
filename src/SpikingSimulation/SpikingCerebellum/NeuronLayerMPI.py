'''
Created on April 22, 2014

@author: Jesus Garrido (jgarridoalcazar at gmail.com)
'''

import NeuronLayerNoMPI
import numpy
import logging

logger = logging.getLogger('Simulation')

class NeuronLayer(NeuronLayerNoMPI.NeuronLayer):
    '''
    This class defines a neuron layer and the data needed to generate it in a simulator.
    It implements the communication to share the relative positions of the neurons in the layer.
    '''
    
    def _share_positions_(self):
        '''
        This function sends the positions stored in the first process to all the MPI processes. This function will be reimplemented in
        inheriting classes. 
        '''
        # ------------------------------------------------------------------------------------------------------------------
        # Broadcast all the positions to the rest of processes.
        from mpi4py import MPI
        
        comm = MPI.COMM_WORLD
        
        nprocs = comm.Get_size()
        rank = comm.Get_rank()
        
        if rank==0:
            positions_array = numpy.reshape(self.relative_positions, self.number_of_neurons*len(self.size), order='C')
        else:
            positions_array = numpy.empty(self.number_of_neurons*len(self.size), dtype=numpy.float32)
        
        # Broadcast the position array
        comm.Bcast([positions_array, MPI.FLOAT], root=0)
            
        #print 'Process', rank,':','Broadcasting position array:',positions_array
            
        self.relative_positions = numpy.reshape(positions_array, (self.number_of_neurons, len(self.size)), order='C')
        
        return
    
    def save_layer(self, root):
        '''
        This function stores the relative positions of the neurons and other attributes in the hdf5 group passed as an argument.
        '''
        
        from mpi4py import MPI
        
        comm = MPI.COMM_WORLD
        
        rank = comm.Get_rank()
        
        if rank==0:
            super(NeuronLayer, self).save_layer(root)
        
        return
    