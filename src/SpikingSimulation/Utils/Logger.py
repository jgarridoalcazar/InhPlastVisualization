# This file should be imported after nest. Otherwise nest might not work as expected
import logging
from mpi4py import MPI

def InitializeLogger(name):
    
    mpirank = MPI.COMM_WORLD.Get_rank()

    class ContextFilter(logging.Filter):
        """
        This is a filter which injects contextual information into the log.
        
        Rather than use actual contextual information, we just use random
        data in this demo.
        """
        def filter(self, record):
            record.mpiid = mpirank
            return True


    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)

    # Create formatter
    formatter = logging.Formatter('%(asctime)s - P%(mpiid)s - %(name)s - %(levelname)s: %(message)s')

    # Create stdout handler
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(formatter)


    # Create filter to add the MPI process id
    filter1 = ContextFilter()

    logger.addFilter(filter1)
    logger.addHandler(handler)

