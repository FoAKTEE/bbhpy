from bbh_pack import *
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
    
    Rfig=kwargs['rmax'] #the limit of anta mem
    tst = athena_read.athdf(filename=tstathdf,x1_min=-Rfig,x1_max=Rfig, x2_min=-Rfig, x2_max=Rfig, x3_min=-Rfig, x3_max=Rfig)
    tst.keys()
    # time = tst['time']

    # construct the meshgrid
    xdat , ydat , zdat = tst['x1v'], tst['x2v'], tst['x3v']
    z, y, x = np.meshgrid(zdat, ydat, xdat, indexing='ij')
    rad = np.sqrt(x**2 + y**2)
    cosphi, sinphi = x/rad, y/rad
    print('shape of the data is: ', np.shape(x))
    nx, ny, nz = np.shape(x)
    zmid , ymid , xmid = nz//2, ny//2, nx//2

    # setting up the vmin, vmax, FigID, and colormap from the quantity
    # dens, eint, bcc1, bcc2, bcc3, velx, vely, velz
    def select_case(input):
        match input:
            case "dens":
                FigID = r'$\rho$'
                vmin, vmax = 1e-5,1.0
                norm = colors.LogNorm(vmin=vmin, vmax=vmax)

                quant_zmid, quant_ymid, quant_xmid  =  tst['dens'][zmid,:,:], tst['dens'][:,ymid,:], tst['dens'][:,:,xmid]
                return FigID, norm, vmin, vmax, quant_zmid, quant_ymid, quant_xmid
            case "eint":
                FigID = r'$e_{\mathrm{int}}$'
                vmin, vmax = 1e-5,1.0
                norm = colors.LogNorm(vmin=vmin, vmax=vmax)

                quant_zmid, quant_ymid, quant_xmid  =  tst['eint'][zmid,:,:], tst['eint'][:,ymid,:], tst['eint'][:,:,xmid]
                return FigID, norm, vmin, vmax, quant_zmid, quant_ymid, quant_xmid
            case "br":
                FigID = r'$B_r$'
                vmin, vmax = -1,1
                norm = colors.Normalize(vmin=vmin, vmax=vmax)

                bx_zmid = tst['bcc1'][zmid,:,:]
                by_zmid = tst['bcc2'][zmid,:,:]

                bx_ymid = tst['bcc1'][:,ymid,:]
                by_ymid = tst['bcc2'][:,ymid,:]

                bx_xmid = tst['bcc1'][:,:,xmid]
                by_xmid = tst['bcc2'][:,:,xmid]

                br_zmid = bx_zmid*cosphi[zmid,:,:] + by_zmid*sinphi[zmid,:,:]   
                br_ymid = bx_ymid*cosphi[:,ymid,:] + by_ymid*sinphi[:,ymid,:]
                br_xmid = bx_xmid*cosphi[:,:,xmid] + by_xmid*sinphi[:,:,xmid]

                quant_zmid, quant_ymid, quant_xmid  =  br_zmid, br_ymid, br_xmid

                return FigID, norm, vmin, vmax, quant_zmid, quant_ymid, quant_xmid
            case "bphi":
                FigID = r'$B_{\phi}$'
                vmin, vmax = -1,1
                norm = colors.Normalize(vmin=vmin, vmax=vmax)

                bx_zmid = tst['bcc1'][zmid,:,:]
                by_zmid = tst['bcc2'][zmid,:,:]
                
                bx_ymid = tst['bcc1'][:,ymid,:]
                by_ymid = tst['bcc2'][:,ymid,:]

                bx_xmid = tst['bcc1'][:,:,xmid]
                by_xmid = tst['bcc2'][:,:,xmid]
  
                bphi_zmid = -bx_zmid*sinphi[zmid,:,:] + by_zmid*cosphi[zmid,:,:]
                bphi_ymid = -bx_ymid*sinphi[:,ymid,:] + by_ymid*cosphi[:,ymid,:]
                bphi_xmid = -bx_xmid*sinphi[:,:,xmid] + by_xmid*cosphi[:,:,xmid]

                quant_zmid, quant_ymid, quant_xmid  =  bphi_zmid, bphi_ymid, bphi_xmid

                return FigID, norm, vmin, vmax, quant_zmid, quant_ymid, quant_xmid
            case "bz":
                FigID = r'$B_z$'
                vmin, vmax = -10,10
                norm = colors.SymLogNorm(vmin=vmin, vmax=vmax,linthresh=1e-3)

                bz_xmid = tst['bcc3'][:,:,xmid]
                bz_ymid = tst['bcc3'][:,ymid,:]
                bz_zmid = tst['bcc3'][zmid,:,:]

                quant_zmid, quant_ymid, quant_xmid  =  bz_zmid, bz_ymid, bz_xmid

                return FigID, norm, vmin, vmax, quant_zmid, quant_ymid, quant_xmid
            case "vel":
                vmin, vmax = 1e-5,1
                norm = colors.Normalize(vmin=vmin, vmax=vmax)
                return FigID, norm, vmin, vmax, quant_zmid, quant_ymid, quant_xmid
            case "sigma":
                vmin, vmax = 1e-4, 1e3
                norm = colors.LogNorm(vmin=vmin, vmax=vmax)
                return FigID, norm, vmin, vmax, quant_zmid, quant_ymid, quant_xmid
            case "beta":
                vmin, vmax = 1e-4, 1e3
                norm = colors.LogNorm(vmin=vmin, vmax=vmax)
                return FigID, norm, vmin, vmax, quant_zmid, quant_ymid, quant_xmid
            case "total":
                dmin, dmax = 1e-4, 1e3
                norm_dens = colors.LogNorm(vmin=dmin, vmax=dmax)
                bmin, bmax = 1e-4, 1e3

                return FigID, norm, vmin, vmax, quant_zmid, quant_ymid, quant_xmid
            case _:
                return "Invalid case"

    # reading in the data: quantity
    # quant = tst[kwargs['quantity']]
    
    FigID, norm, vmin, vmax, quant_zmid, quant_ymid, quant_xmid = select_case(kwargs['quantity'])
    
    x_zmid, y_zmid = x[zmid,:,:], y[zmid,:,:]
    x_ymid, z_ymid = x[:,ymid,:], z[:,ymid,:]
    y_xmid, z_xmid = y[:,:,xmid], z[:,:,xmid]
    
    # def min_max(vmin_in, vmax_in):

    # choosing the colormap
    if kwargs['norm'] == 'logc':
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    elif kwargs['norm'] == 'linear':   
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
    else:
        norm = norm

    fig = plt.figure(figsize=(12.5, 12),dpi=75)
    gs = GridSpec(nrows=3, ncols=5, wspace=0.1, hspace=0.1, width_ratios=[1,1,1,0.0,0.05], height_ratios=[1,1,1])
    ax_fig = np.array([[fig.add_subplot(gs[j,i]) for i in range(3)] for j in range(3)])
    # ax_fig = ax_fig.flatten()
    ax_cbar = np.array([fig.add_subplot(gs[i,4]) for i in range(3)]).flatten()

    RRange_list = np.array([Rfig/4,Rfig/2,Rfig])
    for iRange,RRange in enumerate(RRange_list):
        ax = ax_fig[iRange,:]

        axz = ax[0]
        axy = ax[1]
        axx = ax[2]
        print(RRange)
        Rticks = np.arange(-RRange+RRange/2,RRange,RRange/2)
        Rticklabels = np.array([str(i) for i in Rticks])

        imz = axz.pcolormesh(x_zmid[:,:],y_zmid[:,:],quant_zmid[:,:],cmap='RdBu_r',norm=norm,label=FigID+r' $\mathrm{x-y \ plane}$')
        imy = axy.pcolormesh(x_ymid[:,:],z_ymid[:,:],quant_ymid[:,:],cmap='RdBu_r',norm=norm,label=FigID+r' $\mathrm{x-z \ plane}$')
        imx = axx.pcolormesh(y_xmid[:,:],z_xmid[:,:],quant_xmid[:,:],cmap='RdBu_r',norm=norm,label=FigID+r' $\mathrm{y-z \ plane}$')

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
    parser.add_argument('quantity',
                        help=('name of quantity to be plotted; choose between '
                             ' dens, eint, bcc1, bcc2, bcc3, velx, vely, velz'))
    parser.add_argument('-c',
                        '--colormap',
                        default=None,
                        help=('name of Matplotlib colormap to use instead of default'))
    parser.add_argument('--norm',
                        help=('flag indicating data should be colormapped logarithmically'
                              ' default is chosen according to the quantity, but can be forced'
                              ' to be: linear, logc'))
    parser.add_argument('--rmax',
                        type=float,
                        default=10.0,
                        help=('maximum radius to be plotted'))
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
