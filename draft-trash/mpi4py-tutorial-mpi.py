from . import macros

import numpy as np
if (macros.mpi_enabled):
    from mpi4py import MPI


    def Allreduce(arr, op=MPI.SUM, mpitype=None, **kwargs):

        res = np.empty_like(arr)
        comm = MPI.COMM_WORLD
        
        comm.Barrier()
        comm.Allreduce(arr, res, op=op)

        return res
    

    def sum(arr, **kwargs):
        return Allreduce(arr, op=MPI.SUM, **kwargs)
    
    def max(arr, **kwargs):
        return Allreduce(arr, op=MPI.MAX, **kwargs)
    
    def min(arr, **kwargs):
        return Allreduce(arr, op=MPI.MIN, **kwargs)