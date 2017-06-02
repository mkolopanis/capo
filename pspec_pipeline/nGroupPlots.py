import numpy as np, matplotlib.pyplot as plt
import sys, optparse
o = optparse.OptionParser()
o.add_option('--path2data',default='/Users/saulkohn/Documents/PAPER/powerspectra/nBlGroupExperiment/',help='Path to the top of the directory tree containing the power spectra')
o.add_option('-k','--kbin',default=3,help='Single k_parallel bin to plot against increasing N')
o.add_option('--nblg_min',default=2,help='Minimum number of baseline groups used')
o.add_option('--nblg_max',default=23,help='Maximum number of baseline groups used')
o.add_option('--nlst_min',default=1,help='Minimum number of LST groups used')
o.add_option('--nlst_max',default=10,help='Maximum number of LST groups used')
o.add_option('--polyOrder',default=None,help='Order of polynomial to fit')
o.add_option('--plotSingleVsNblg', action='store_true', help='Plot error w.r.t. N_bl-groups for a single k_parallel bin')
o.add_option('--plot3d', action='store_true', help='Plot error w.r.t. N_bl-groups for all k_parallel bins')
o.add_option('--plotAvgVsNblg', action='store_true', help='Plot error w.r.t N_bl-groups for the average of the outer k_parallel bins (first and last 5 bins)')
o.add_option('--plotAvgVsNlst', action='store_true', help='Plot error w.r.t N_lst-groups for the average of the outer k_parallel bins (first and last 5 bins)')
o.add_option('--plotConcat', action='store_true', help='Plot error w.r.t. N samples (N_bl and N_lst)-groups')
o.add_option('--plotManyk', action='store_true', help='Plot error w.r.t. N samples for each group type')
o.add_option('--plotAll',action='store_true',help='plot all the things')
o.add_option('--k3',action='store_true')

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
pCv_bl,pCv_lst = np.zeros((nblg_L,21)),np.zeros((nlst_L,21))

nbl_g,nlst_g = np.zeros((nblg_L)),np.zeros((nlst_L))

# Data Parsing

for i,nBL in enumerate(range(nblg_min,nblg_max+1)):
    d = np.load('%s/nblg%i/inject_sep0,1_0.01/pspec_pk_k3pk.npz'%(path2data,nBL))
    pIv_bl[i,:] = d['pIv_err']
    pIn_bl[i,:] = d['pIn_err']
    pCv_bl[i,:] = d['pCv_err']
    nbl_g[i] = d['nbls_g']
    if i==0:
        for j,nLST in enumerate(range(nlst_min,nlst_max+1)):
            dl = np.load('%s/nblg%i/inject_sep0,1_0.01/nlst%i/pspec_pk_k3pk.npz'%(path2data,nBL,nLST))
            pIv_lst[j,:] = dl['pIv_err']
            pIn_lst[j,:] = dl['pIn_err']
            pCv_lst[j,:] = dl['pCv_err']
            nlst_g[j] = dl['nlsts_g']
        
kpl = d['kpl']

nice_minus=[0,1,2,3,4]
nice_plus = [16,17,19,20]
nice = [0,1,2,3,4,16,17,19,20]

pIn_bl_avg = 0.5*(np.mean(pIn_bl[:,nice_minus],axis=1)+np.mean(pIn_bl[:,nice_plus],axis=1))
pIv_bl_avg = 0.5*(np.mean(pIv_bl[:,nice_minus],axis=1)+np.mean(pIv_bl[:,nice_plus],axis=1))
pCv_bl_avg = 0.5*(np.mean(pCv_bl[:,nice_minus],axis=1)+np.mean(pCv_bl[:,nice_plus],axis=1))

pIn_lst_avg = 0.5*(np.mean(pIn_lst[:,nice_minus],axis=1)+np.mean(pIn_lst[:,nice_plus],axis=1))
pIv_lst_avg = 0.5*(np.mean(pIv_lst[:,nice_minus],axis=1)+np.mean(pIv_lst[:,nice_plus],axis=1))
pCv_lst_avg = 0.5*(np.mean(pCv_lst[:,nice_minus],axis=1)+np.mean(pCv_lst[:,nice_plus],axis=1))

if not opts.polyOrder is None:
    order = int(opts.polyOrder)
    # Fitting routines
    coeffs_nbl,C_nbl = np.polyfit(np.log10(nbl_g),np.log10(pIv_bl_avg),deg=order,cov=True)
    coeffs_nlst,C_nlst = np.polyfit(np.log10(nlst_g),np.log10(pIv_lst_avg),deg=order,cov=True)
    
    nblT = np.vstack([np.log10(nbl_g)**(order-i) for i in range(order+1)]).T
    nlstT = np.vstack([np.log10(nlst_g)**(order-i) for i in range(order+1)]).T
    C_nbl_fit = np.dot(nblT,np.dot(C_nbl, nblT.T))
    C_nlst_fit = np.dot(nlstT,np.dot(C_nlst, nlstT.T))
    sig_nbl_fit = np.sqrt(np.diag(C_nbl_fit))
    sig_nlst_fit = np.sqrt(np.diag(C_nlst_fit))
    
    err_nbl_fit = np.mean(sig_nbl_fit)#np.sum([sig**2 for sig in sig_nbl_fit])#/len(sig_nbl_fit)
    err_nlst_fit = np.mean(sig_nlst_fit)#np.sum([sig**2 for sig in sig_nlst_fit])#/len(sig_nlst_fit)
    
    print 'Coeffs for log10(nbl_g) vs log10(pIv_bl_avg): ',coeffs_nbl, '+/-', err_nbl_fit
    print 'Coeffs for log10(nlst_g) vs log10(pIv_lst_avg): ',coeffs_nlst, '+/-', err_nlst_fit
    poly_nbl = np.poly1d(coeffs_nbl)
    poly_nlst= np.poly1d(coeffs_nlst)
        
    f,axarr = plt.subplots(1,2,sharey=True)
    axarr[0].plot(np.log10(nbl_g),np.log10(pIn_bl_avg[-1]*nbl_g[-1]/nbl_g),'k--',label='analytic')
    axarr[0].plot(np.log10(nbl_g),np.log10(pIv_bl_avg),'b-',label='pIv_bl_err avg')
    axarr[0].plot(np.log10(nbl_g),np.log10(pIn_bl_avg),'g-',label='pIn_bl_err avg')
    axarr[0].plot(np.log10(nbl_g),poly_nbl(np.log10(nbl_g)),'c:',label='order %i fit'%order)
    axarr[0].plot(np.log10(nbl_g),np.log10(pCv_bl_avg),'r-',label='pCv_bl_err avg')
    axarr[0].set_xlabel('log10(Nbls/group)')
    axarr[0].set_ylabel('log10(P(k))')
    axarr[0].legend()

    axarr[1].plot(np.log10(nlst_g),np.log10(pIn_lst_avg[-1]*np.sqrt(nlst_g[-1]/nlst_g)),'k--',label='analytic')
    axarr[1].plot(np.log10(nlst_g),np.log10(pIv_lst_avg),'b-',label='pIv_lst_err avg')
    axarr[1].plot(np.log10(nlst_g),np.log10(pIn_lst_avg),'g-',label='pIn_lst_err avg')
    axarr[1].plot(np.log10(nlst_g),np.log10(pCv_lst_avg),'r-',label='pCv_lst_err avg')
    axarr[1].plot(np.log10(nlst_g),poly_nlst(np.log10(nlst_g)),'c:',label='order %i fit'%order)
    axarr[1].set_xlabel('log10(Nlsts/group)')
    axarr[1].legend()
    
    plt.show()

# Plotting
if opts.plotManyk or opts.plotAll:
    f,axarr = plt.subplots(1,2,sharey=True)
    for kbin in nice:
        if opts.k3:
            factor = np.abs(kpl[kbin])**3 / (2*np.pi)
        else:
            factor = 1.
        axarr[0].plot(np.log10(nbl_g),np.log10(factor*pIv_bl[:,kbin]),'b-',alpha=0.2)
        axarr[0].plot(np.log10(nbl_g),np.log10(factor*pIn_bl[:,kbin]),'g-',alpha=0.2)
        axarr[0].plot(np.log10(nbl_g),np.log10(factor*pCv_bl[:,kbin]),'r-',alpha=0.2)
        
        axarr[1].plot(np.log10(nlst_g),np.log10(factor*pIv_lst[:,kbin]),'b-',alpha=0.2)
        axarr[1].plot(np.log10(nlst_g),np.log10(factor*pIn_lst[:,kbin]),'g-',alpha=0.2)
        axarr[1].plot(np.log10(nlst_g),np.log10(factor*pCv_lst[:,kbin]),'r-',alpha=0.2)
    axarr[0].set_xlabel('log10(Nbls/group)')
    if not opts.k3:
        axarr[0].set_ylabel('log10(P(k))')
    else:
        axarr[0].set_ylabel('log10($\Delta^2$(k))')
    axarr[1].set_xlabel('log10(Nlsts/group)')
    plt.show()

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
    f,axarr = plt.subplots(1,2,sharex=True)
    im0=axarr[0].imshow(np.log10(pIv_bl), aspect='auto', interpolation='nearest', extent=(kpl[0],kpl[-1],23,2) )
    f.colorbar(im0,ax=axarr[0])
    axarr[0].set_xlabel('k parallel')
    axarr[0].set_label('# bl groups')
    
    im1=axarr[1].imshow(np.log10(pIv_lst), aspect='auto', interpolation='nearest', extent=(kpl[0],kpl[-1],10,1) )
    f.colorbar(im1,ax=axarr[1])
    axarr[1].set_xlabel('k parallel')
    axarr[1].set_label('# lst groups')

    f2,axarr2 = plt.subplots(1,2,sharex=True)
    
    im02=axarr2[0].imshow(np.log10(pIn_bl), aspect='auto', interpolation='nearest', extent=(kpl[0],kpl[-1],23,2) )
    f2.colorbar(im02,ax=axarr2[0])
    axarr2[0].set_xlabel('k parallel')
    axarr2[0].set_label('# bl groups')
    im12=axarr2[1].imshow(np.log10(pIn_lst), aspect='auto', interpolation='nearest', extent=(kpl[0],kpl[-1],10,1) )
    f2.colorbar(im12,ax=axarr2[1])
    axarr2[1].set_xlabel('k parallel')
    axarr2[1].set_label('# lst groups')
    
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
    f,axarr = plt.subplots(1,2,sharey=True)
    
    axarr[0].loglog(nbl_g,pIn_bl_avg[-1]*nbl_g[-1]/nbl_g,'k--',label='analytic')
    axarr[0].loglog(nbl_g,pIv_bl_avg,'b-',label='pIv_bl_err avg')
    axarr[0].loglog(nbl_g,pIn_bl_avg,'g-',label='pIn_bl_err avg')
    axarr[0].set_xlabel('Nbls/group')
    axarr[0].set_ylabel('P(k)')
    axarr[0].legend()

    axarr[1].loglog(nlst_g,pIn_lst_avg[-1]*np.sqrt(nlst_g[-1]/nlst_g),'k--',label='analytic')
    axarr[1].loglog(nlst_g,pIv_lst_avg,'b-',label='pIv_lst_err avg')
    axarr[1].loglog(nlst_g,pIn_lst_avg,'g-',label='pIn_lst_err avg')
    axarr[1].set_xlabel('Nlsts/group')
    axarr[1].legend()
    
    plt.show()
