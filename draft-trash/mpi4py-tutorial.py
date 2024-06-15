from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

data_size = 10000
data_per_process = data_size // size
data = np.arange(rank*data_per_process, (rank+1)*data_per_process, dtype=np.float64)

local_sum = np.sum(data)

global_sum = comm.reduce(local_sum, op=MPI.SUM, root=0)

if rank == 0:
    print("Global sum:", global_sum)
