# %% [markdown]
# ## Quick Look Snapshoot

# %%
# %load_ext autoreload
# %autoreload 2

import os
import glob
import pathlib
import sys
import argparse
import warnings
import h5py
import numpy as np
import yt
import matplotlib
import time
import scipy
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib import rcParams
import multiprocessing
from multiprocessing import Process, Lock, Array
from pathlib import Path

sys.path.append("/Users/hyw/erm/ppscript/vis/python/") # local
sys.path.append("/home/haiyangw/vis/python/") # server: anta

import bin_convert as bc

sys.path.append("/home/emost/codes/athenak/vis/python/") # server: anta
import athena_read

plt.style.use('ncr-paper.mplstyle')
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "sans-serif",
#     "font.sans-serif": ["Computer Modern Sans Serif"]})

def main(**kwargs):

# %%
    # directory_path = Path('/Users/hyw/Desktop/ERm/athenak-erm/build-mdisk/src/mhdtest/bin/')
    # directory_path = Path('/Users/hyw/Desktop/ERm/athenak-erm/build-mdisk/src/mhdtest/bin/')
    # extension = '.bin'

    # for file_path in directory_path.iterdir():
    #     if file_path.is_file() and file_path.suffix == extension:
    #         print(file_path)
    #         binary_fname = str(file_path)
    #         athdf_fname = binary_fname.replace(".bin", ".athdf")
    #         xdmf_fname = athdf_fname + ".xdmf"
    #         filedata = bc.read_binary(binary_fname)

    # binary_fname = kwargs['data_file']
    # athdf_fname = binary_fname.replace(".bin", ".athdf")
    # xdmf_fname = athdf_fname + ".xdmf"
    # filedata = bc.read_binary(binary_fname)


# %%

    tstathdf = '/Users/hyw/Desktop/ERM/athenak-erm/build-mdisk/src/mhdtest/bin/cbd.mhd_w_bcc.02330.athdf'
    tstathdf = '/Users/hyw/CBD_Summit/cbd.mhd_w_bcc.00020.athdf'
    # tstathdf = athdf_fname
    tstathdf = kwargs['data_file']
    # tstathdf = '/Users/hyw/24.6.1-test/cbd.mhd_w_bcc.00008.athdf'
    # tstathdf = '/Users/hyw/athenak/hydro-test/src/testmb/lowreso/cbd.mhd_w_bcc.00215.athdf'
    
    Rfig=20.0
    tst_athdf = athena_read.athdf(filename=tstathdf,x1_min=-Rfig,x1_max=Rfig, x2_min=-Rfig, x2_max=Rfig, x3_min=-Rfig, x3_max=Rfig)
    tst_athdf.keys()
    
    # %%
    tst = tst_athdf
    
    # initialize the data (read-in)
    # time = tst['time']
    xdat = tst['x1v']
    ydat = tst['x2v']
    zdat = tst['x3v']
    dens = tst['dens']
    eint = tst['eint']
    bcc1 = tst['bcc1']
    bcc2 = tst['bcc2']
    bcc3 = tst['bcc3']
    vx = tst['velx']
    vy = tst['vely']
    vz = tst['velz']
    
    z, y, x = np.meshgrid(zdat, ydat, xdat, indexing='ij')
    print('shape of the data is: ', np.shape(x))
    
    nx, ny, nz = np.shape(x)
    
    
    # index of the center of the domain (index truncation)
    xs, xe, ys, ye, zs, ze = nx//4, nx//4*3, ny//4, ny//4*3, nz//4, nz//4*3
    # xs, xe, ys, ye, zs, ze = nx//8*3, nx//8*5, ny//8*3, ny//8*5, nz//8*3, nz//8*5
    # xs, xe, ys, ye, zs, ze = nx//16*7, nx//16*9, ny//16*7, ny//16*9, nz//16*7, nz//16*9
    # xs, xe, ys, ye, zs, ze = nx//32*15, nx//32*17, ny//32*15, ny//32*17, nz//32*15, nz//32*17
    # xs, xe, ys, ye, zs, ze = nx//64*31, nx//64*33, ny//64*31, ny//64*33, nz//64*31, nz//64*33   
    # xs, xe, ys, ye, zs, ze = nx//128*63, nx//128*65, ny//128*63, ny//128*65, nz//128*63, nz//128*65
    
    nx_tst, ny_tst, nz_tst = xe-xs, ye-ys, ze-zs
    
    # truncated coordinates and mesh data
    x_tst = x[xs:xe,ys:ye,zs:ze]
    y_tst = y[xs:xe,ys:ye,zs:ze]
    z_tst = z[xs:xe,ys:ye,zs:ze]
    
    xdat_tst = xdat[xs:xe]
    ydat_tst = ydat[xs:xe]
    zdat_tst = zdat[xs:xe]
    
    # truncate the data to only the inner region
    dens_tst = dens[xs:xe,ys:ye,zs:ze]
    eint_tst = eint[xs:xe,ys:ye,zs:ze]
    bcc1_tst = bcc1[xs:xe,ys:ye,zs:ze]
    bcc2_tst = bcc2[xs:xe,ys:ye,zs:ze]
    bcc3_tst = bcc3[xs:xe,ys:ye,zs:ze]
    vx_tst = vx[xs:xe,ys:ye,zs:ze]
    vy_tst = vy[xs:xe,ys:ye,zs:ze]
    vz_tst = vz[xs:xe,ys:ye,zs:ze]
    
    # reformulate the data into cyl coordinates
    darray_tst = [dens_tst, eint_tst, bcc1_tst, bcc2_tst, bcc3_tst, vx_tst, vy_tst, vz_tst]
    
    print('shape of the truncated data is: ', np.shape(x_tst))
    # print('rad array: ', rad)
    
    
    # %%
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    
    fig = plt.figure(figsize=(12.5, 12),dpi=75)
    gs = GridSpec(nrows=3, ncols=5, wspace=0.1, hspace=0.1, width_ratios=[1,1,1,0.0,0.05], height_ratios=[1,1,1])
    ax_fig = np.array([[fig.add_subplot(gs[j,i]) for i in range(3)] for j in range(3)])
    # ax_fig = ax_fig.flatten()
    ax_cbar = np.array([fig.add_subplot(gs[i,4]) for i in range(3)]).flatten()
    
    zmid = nz_tst//2
    ymid = ny_tst//2
    xmid = nx_tst//2
    
    vmin,vmax = 1e-5,1.0
    quant = dens_tst
    FigID = r'$\rho$'
    
    RRange_list = np.array([2.5,5.0,10.0])
    for iRange,RRange in enumerate(RRange_list):
        ax = ax_fig[iRange,:]
        
        axz = ax[0]
        axy = ax[1]
        axx = ax[2]
        print(RRange)
        Rticks = np.arange(-RRange+RRange/2,RRange,RRange/2)
        Rticklabels = np.array([str(i) for i in Rticks])
        
        imz = axz.pcolormesh(x_tst[zmid,:,:],y_tst[zmid,:,:],quant[zmid,:,:],cmap='RdBu_r',norm=colors.LogNorm(vmin=vmin, vmax=vmax),label=FigID+r' $\mathrm{x-y \ plane}$')
        imy = axy.pcolormesh(x_tst[:,ymid,:],z_tst[:,ymid,:],quant[:,ymid,:],cmap='RdBu_r',norm=colors.LogNorm(vmin=vmin, vmax=vmax),label=FigID+r' $\mathrm{x-z \ plane}$')
        imx = axx.pcolormesh(y_tst[:,:,xmid],z_tst[:,:,xmid],quant[:,:,xmid],cmap='RdBu_r',norm=colors.LogNorm(vmin=vmin, vmax=vmax),label=FigID+r' $\mathrm{y-z \ plane}$')
        
        patchz = mpatches.Patch(color=None,alpha=0.0, label=FigID+r' $\mathrm{x-y \ plane}$')
        patchy = mpatches.Patch(color=None,alpha=0.0, label=FigID+r' $\mathrm{x-z \ plane}$')
        patchx = mpatches.Patch(color=None,alpha=0.0, label=FigID+r' $\mathrm{y-z \ plane}$')
        
        axx.legend(handles=[patchz],loc='upper right',frameon=False,fontsize=10,labelcolor='black')
        axy.legend(handles=[patchy],loc='upper right',frameon=False,fontsize=10,labelcolor='black')
        axz.legend(handles=[patchx],loc='upper right',frameon=False,fontsize=10,labelcolor='black')
        
        for idirec in range(3):
            ax[idirec].set_xlim(-RRange,RRange)
            ax[idirec].set_ylim(-RRange,RRange)
    
            ax[idirec].set_xticks(Rticks)
            # ax[idirec].set_xticklabels([])
            ax[idirec].set_yticks(Rticks)
            ax[idirec].set_yticklabels([])
            ax[idirec].set_aspect("equal")
    
        cax1 = ax_cbar[iRange]
        cax1.tick_params(which='major',direction='in',labelsize=13,length=7.5,labelleft=True,labelright=False)
        cax1.tick_params(which='minor',direction='in')
        cax1.yaxis.set_minor_locator(mticker.MultipleLocator(0.5))
        colorbar1 = plt.colorbar(imz,cax=cax1)
    

    plt.savefig(kwargs['output_file'],dpi=250)

# %%

# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('data_file',
                        help='name of input file, possibly including path')
    parser.add_argument('output_file',
                        default='./test.png',
                        help=('name of output to be (over)written, possibly including '
                              'path; use "show" to show interactive plot instead'))
    parser.add_argument('--vmin',
                        type=float,
                        default=1e-5, #None,
                        help=('data value to correspond to colormap minimum; use '
                              '--vmin=<val> if <val> has negative sign'))
    parser.add_argument('--vmax',
                        type=float,
                        default=1.0, #None,
                        help=('data value to correspond to colormap maximum; use '
                              '--vmax=<val> if <val> has negative sign'))
    # parser.add_argument('quantity',
    #                     help='name of quantity to be plotted')
#     parser.add_argument('-d', '--direction',
#                         type=int,
#                         choices=(1, 2, 3),
#                         default=3,
#                         help=('direction orthogonal to slice for 3D data'))
#     parser.add_argument('--slice_location',
#                         type=float,
#                         default=None,
#                         help=('coordinate value along which slice is to be taken '
#                               '(default: 0)'))
#     parser.add_argument('-a', '--average',
#                         action='store_true',
#                         help=('flag indicating averaging should be done in orthogonal '
#                               'direction for 3D data'))
#     parser.add_argument('-s', '--sum',
#                         action='store_true',
#                         help=('flag indicating summation should be done in orthogonal '
#                               'direction for 3D data'))
#     parser.add_argument('-l',
#                         '--level',
#                         type=int,
#                         default=None,
#                         help=('refinement level to be used in plotting (default: max '
#                               'level in file)'))
#     parser.add_argument('--x_min',
#                         type=float,
#                         default=None,
#                         help='minimum extent of plot in first plotted direction')
#     parser.add_argument('--x_max',
#                         type=float,
#                         default=None,
#                         help='maximum extent of plot in first plotted direction')
#     parser.add_argument('--y_min',
#                         type=float,
#                         default=None,
#                         help='minimum extent of plot in second plotted direction')
#     parser.add_argument('--y_max',
#                         type=float,
#                         default=None,
#                         help='maximum extent of plot in second plotted direction')
#     parser.add_argument('-f', '--fill',
#                         action='store_true',
#                         help='flag indicating image should fill plot area, even if this '
#                              'distorts the aspect ratio')
#     parser.add_argument('-c',
#                         '--colormap',
#                         default=None,
#                         help=('name of Matplotlib colormap to use instead of default'))
#     parser.add_argument('--logc',
#                         action='store_true',
#                         help='flag indicating data should be colormapped logarithmically')
#     parser.add_argument('--stream',
#                         default=None,
#                         help='name of vector quantity to use to make stream plot')
#     parser.add_argument('--stream_average',
#                         action='store_true',
#                         help='flag indicating stream plot should be averaged in '
#                              'orthogonal direction for 3D data')
#     parser.add_argument('--stream_density',
#                         type=float,
#                         default=1.0,
#                         help='density of stream lines')
    args = parser.parse_args()
    main(**vars(args))


# %%
