import numpy as np, matplotlib.pyplot as plt
import sys, optparse
o = optparse.OptionParser()
o.add_option('--path2data',default='/Users/saulkohn/Documents/PAPER/powerspectra/nBlGroupExperiment/',help='Path to the top of the directory tree containing the power spectra')
o.add_option('-k','--kbin',default=3,help='Single k_parallel bin to plot against increasing N')
o.add_option('--nblg_min',default=2,help='Minimum number of baseline groups used')
o.add_option('--nblg_max',default=23,help='Maximum number of baseline groups used')
o.add_option('--nlst_min',default=1,help='Minimum number of LST groups used')
o.add_option('--nlst_max',default=10,help='Maximum number of LST groups used')

o.add_option('--plotSingleVsNblg', action='store_true', help='Plot error w.r.t. N_bl-groups for a single k_parallel bin')
o.add_option('--plot3d', action='store_true', help='Plot error w.r.t. N_bl-groups for all k_parallel bins')
o.add_option('--plotAvgVsNblg', action='store_true', help='Plot error w.r.t N_bl-groups for the average of the outer k_parallel bins (first and last 5 bins)')
o.add_option('--plotAvgVsNlst', action='store_true', help='Plot error w.r.t N_lst-groups for the average of the outer k_parallel bins (first and last 5 bins)')
o.add_option('--plotConcat', action='store_true', help='Plot error w.r.t. N samples (N_bl and N_lst)-groups')

o.add_option('--plotAll',action='store_true',help='plot all the things')

opts, args = o.parse_args(sys.argv[1:])

# Option Parsing

path2data = opts.path2data
kbin = int(opts.kbin)
nblg_min,nblg_max = int(opts.nblg_min),int(opts.nblg_max)
nlst_min,nlst_max = int(opts.nlst_min),int(opts.nlst_max)
nblg_L = (nblg_max - nblg_min) + 1
nlst_L = (nlst_max - nlst_min) + 1

pIv_bl,pIn_bl = np.zeros((nblg_L,21)),np.zeros((nblg_L,21)) #21 k_para bins
pIv_lst,pIn_lst = np.zeros((nlst_L,21)),np.zeros((nlst_L,21))
nbl_g,nlst_g = np.zeros((nblg_L)),np.zeros((nlst_L))

# Data Parsing

for i,nBL in enumerate(range(nblg_min,nblg_max+1)):
    d = np.load('%s/nblg%i/inject_sep0,1_0.01/pspec_pk_k3pk.npz'%(path2data,nBL))
    pIv_bl[i,:] = d['pIv_err']
    pIn_bl[i,:] = d['pIn_err']
    nbl_g[i] = d['nbls_g']
    if i==0:
        for j,nLST in enumerate(range(nlst_min,nlst_max+1)):
            dl = np.load('%s/nblg%i/inject_sep0,1_0.01/nlst%i/pspec_pk_k3pk.npz'%(path2data,nBL,nLST))
            pIv_lst[j,:] = dl['pIv_err']
            pIn_lst[j,:] = dl['pIn_err']
            nlst_g[j] = dl['nlsts_g']
        
kpl = d['kpl']
pIn_bl_avg = 0.5*(np.mean(pIn_bl[:,:5],axis=1)+np.mean(pIn_bl[:,16:],axis=1))
pIv_bl_avg = 0.5*(np.mean(pIv_bl[:,:5],axis=1)+np.mean(pIv_bl[:,16:],axis=1))
pIn_lst_avg = 0.5*(np.mean(pIn_lst[:,:5],axis=1)+np.mean(pIn_lst[:,16:],axis=1))
pIv_lst_avg = 0.5*(np.mean(pIv_lst[:,:5],axis=1)+np.mean(pIv_lst[:,16:],axis=1))
# Plotting
if opts.plotSingleVsNblg or opts.plotAll:
    plt.loglog(nbl_g,pIn_bl[:,kbin][-1]*nbl_g[-1]/nbl_g,'k--',label='analytic')
    plt.loglog(nbl_g,pIv_bl[:,kbin],'b-',label='pIv_bl_err')
    plt.loglog(nbl_g,pIn_bl[:,kbin],'g-',label='pIn_bl_err')
    
    plt.xlabel('Nbls/group')
    plt.ylabel('P(k)')
    plt.suptitle('k_para=%f'%kpl[kbin])
    plt.legend()
    plt.show()

if opts.plot3d or opts.plotAll:
    plt.imshow(np.log10(pIv_bl), aspect='auto', interpolation='nearest', extent=(kpl[0],kpl[-1],23,2) )
    plt.colorbar()
    plt.xlabel('k parallel')
    plt.ylabel('Ngroups')
    plt.show()

if opts.plotAvgVsNblg or opts.plotAll:
    plt.loglog(nbl_g,pIn_bl_avg[-1]*nbl_g[-1]/nbl_g,'k--',label='analytic')
    plt.loglog(nbl_g,pIv_bl_avg,'b-',label='pIv_bl_err avg')
    plt.loglog(nbl_g,pIn_bl_avg,'g-',label='pIn_bl_err avg')
    
    plt.xlabel('Nbls/group')
    plt.ylabel('P(k)')
    plt.legend()
    plt.show()

if opts.plotAvgVsNlst or opts.plotAll:
    plt.loglog(nlst_g,pIn_lst_avg[-1]*np.sqrt(nlst_g[-1]/nlst_g),'k--',label='analytic')
    plt.loglog(nlst_g,pIv_lst_avg,'b-',label='pIv_lst_err avg')
    plt.loglog(nlst_g,pIn_lst_avg,'g-',label='pIn_lst_err avg')
    
    plt.xlabel('Nlsts/group')
    plt.ylabel('P(k)')
    plt.legend()
    plt.show()

if opts.plotConcat or opts.plotAll:
    #n_g = np.concatenate((nbl_g,nlst_g))
    plt.loglog(nlst_g[0]+nbl_g,pIn_bl_avg[-1]*nbl_g[-1]/nbl_g,'k--',label='1/N_blGrp')
    plt.loglog(nlst_g,pIn_lst_avg[-1]*np.sqrt(nlst_g[-1]/nlst_g),'k:',label='1/sqrt(N_lstGrp)')
    plt.loglog(nlst_g[0]+nbl_g,pIv_bl_avg,'b-',label='pIv_bl_err avg')
    plt.loglog(nlst_g,pIv_lst_avg,'b--',label='pIv_lst_err avg')
    plt.loglog(nlst_g[0]+nbl_g,pIn_bl_avg,'g-',label='pIn_bl_err avg')
    plt.loglog(nlst_g,pIn_lst_avg,'g--',label='pIn_lst_err avg')
    
    plt.xlabel('N_samples')
    plt.ylabel('P(k)')
    plt.legend()
    plt.show()
