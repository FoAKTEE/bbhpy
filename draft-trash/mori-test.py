import os 
import pathlib 
import glob
import matplotlib.pyplot as plt 
import numpy as np 
import h5py 
import athena_read


def main_multiple_runs(): 



def main_test_run():



class AthdfAnalyzer:

    def __init__(self, data_prim, savedir, number=None, r0=0.6):
        self.save_dirpath = pathlib.Path(savedir) #save directory (useless info)
        self.pdat = data_prim #primitive variables
        self.number = number 

        #set coordinates
        self.r = data_prim['x1v']*r0
        self.theta = data_prim['x2v']
        
        #set meshgrid in const-phi plane
        self.rr, self.tt = np.meshgrid(self.r, self.theta)
        
        #set cylindrical coordinates
        self.R = self.rr * np.sin(self.tt)
        self.z = self.rr * np.cos(self.tt)


    def overplot_B_stream_2d(self):
        _beta = 
        RR, zz, BR, Bz, beta = 

        strm = 

        plt.colorbar(strm.lines, label=r'log($\beta$)')

    def interp_Rz(self, vals, Rlim=[0,30],zlim=[-15,15],n=301):











def compressed_x2(xmin, xmax, xrat_root, nfaces):

    x2rat = np.abs(xrat_root)
    xmid = 0.5*np.pi

    x = np.arange(nfaces)/(nfaces-1)

    def func(x):
        if x<0.5: 
            ratn = x2rat ** (0.5*nfaces)
            


        else:

    return np.vectorize(func)(x)



if __name__ == "__main__":
    main_test_run()
    #main_multiple_runs()


