#! /usr/bin/env python
"""Compute Signal Loss factor via distribution method v3.
    -Reads in pspec_2d_to_1d.npz (contains PS points for all bootstraps)
    -Finds distribution of Pin vs. Pout curve (using absolute values of points and 2D kernel density estimator) = M
    -Finds distribution of PS points = p (using 1D KDE)
    -Multilpies p element-wise for every column of M (for every Pin)
    -Sums up each row of M, yielding one final distribution along Pin
    -Computes power spectrum points and errors from the distribution
"""
import numpy as n
import pylab as p
from matplotlib import gridspec
from scipy.interpolate import interp1d
from capo.pspec import f2z
import scipy
import glob
import optparse
import sys

o = optparse.OptionParser()
o.set_description(__doc__)
o.add_option('--plot',action='store_true',
            help='Show plots.')
o.add_option('--skip_sigloss',action='store_true',
            help='Save values without correcting for signal loss.')
o.add_option('--sep',default='0,2',
            help='Separation to run script for.')
opts,args = o.parse_args(sys.argv[1:])

# read one pspec_pk_k3pk.npz file for k points
file = n.load(glob.glob('inject_sep'+opts.sep+'*')[0]+'/pspec_pk_k3pk.npz')
kpl = file['kpl']; k = file['k']; kpl_fold = file['kpl_fold']

# loop over injects and read pspec_2d_to_1d outputs
for count in range(2):
    if count == 1: print 'NOISE CASE'
    else: print 'DATA CASE'
    if opts.plot:
        fig1 = p.figure(1, figsize=(15, 7))
        fig2 = p.figure(2, figsize=(15, 7))
    pklo, pkhi = 1e-1, 1e12
    Pouts_fold = {}
    Pouts_I_fold = {}
    Pins_fold = {}
    pCs = {}; pCs_fold = {}
    pIs = {}; pIs_fold = {}
    for inject in glob.glob('inject_sep'+opts.sep+'*'):
        print 'Reading', inject
        file_2d = n.load(inject + '/pspec_2d_to_1d.npz')

        if count == 1: # noise case
            color='b.'
            linecolor='cyan'
            Pout_fold = file_2d['pCs-pCn_fold']
            Pout_I_fold = file_2d['pIs-pIn_fold']
            pC = file_2d['pCn']
            pC_fold = file_2d['pCn_fold']
            pI = file_2d['pIn']
            pI_fold = file_2d['pIn_fold']
        else:
            color='k.'
            linecolor='0.5'
            Pout_fold = file_2d['pCr-pCv_fold'] # shape (#boots, #kpl)
            Pout_I_fold = file_2d['pIr-pIv_fold']  # pI case
            pC = file_2d['pCv']
            pC_fold = file_2d['pCv_fold']
            pI = file_2d['pIv']
            pI_fold = file_2d['pIv_fold']

        Pin_fold = file_2d['pIe_fold']

        for ind in range(len(kpl_fold)): # loop through k for Delta^2(k)
            try:
                Pouts_fold[kpl_fold[ind]].append(Pout_fold[:,ind])
                Pouts_I_fold[kpl_fold[ind]].append(Pout_I_fold[:,ind])
                Pins_fold[kpl_fold[ind]].append(Pin_fold[:,ind])
                pCs_fold[kpl_fold[ind]] = [pC_fold[:,ind]] # no appending because it's the same for every inject
                pIs_fold[kpl_fold[ind]] = [pI_fold[:,ind]]
            except:
                Pouts_fold[kpl_fold[ind]] = [Pout_fold[:,ind]]
                Pouts_I_fold[kpl_fold[ind]] = [Pout_I_fold[:,ind]]
                Pins_fold[kpl_fold[ind]] = [Pin_fold[:,ind]]
                pCs_fold[kpl_fold[ind]] = [pC_fold[:,ind]]
                pIs_fold[kpl_fold[ind]] = [pI_fold[:,ind]]
   
        for ind in range(len(kpl)): # loop through k for P(k)
                pCs[kpl[ind]] = [pC[:,ind]] # no appending because it's the same for every inject
                pIs[kpl[ind]] = [pI[:,ind]]
    

    # Values of interest are contained within dictionaries indexed by kpl_fold:
    #     Pins_fold
    #     Pouts_fold
    #     Pouts_I_fold
    #     pCs, pCs_fold, pIs, pIs_fold  (data/noise values)


    # Plot Pin vs. Pout curves 
    for ind,k in enumerate(Pins_fold.keys()):
       
        if opts.plot:
            p.figure(1) # Pin vs. Pout
            p.subplot(3, 4, ind+1)
            p.plot(Pins_fold[k], Pouts_fold[k], color)  # points
            p.plot(Pins_fold[k], poly[k], 'r-') # polyfit
            #p.plot([pklo, pkhi], [pklo, pkhi], 'k-')  # diagonal line
            p.plot([-pkhi, pkhi], [-pkhi, pkhi], 'k-')  # diagonal line
            p.grid(True)
            p.xscale('symlog',linthreshx=1e6)
            p.yscale('symlog',linthreshy=1e6)
            #p.xlim(pklo, pkhi)
            #p.ylim(pklo, pkhi)
            p.tick_params(axis='both', which='both', labelsize=6)
            p.title('k = '+str("%.4f" % kpl_fold[ind]), fontsize=8)

            p.figure(2) # Pin vs. Pout for I case
            p.subplot(3, 4, ind+1)
            p.plot(Pins_fold[k], Pouts_I_fold[k], color)  # points
            p.plot(Pins_fold[k], poly_I[k], 'r-') # polyfit
            #p.plot([pklo, pkhi], [pklo, pkhi], 'k-')  # diagonal line
            p.plot([-pkhi, pkhi], [-pkhi, pkhi], 'k-')  # diagonal line
            p.grid(True)
            #p.xlim(pklo, pkhi)
            #p.ylim(pklo, pkhi)
            p.xscale('symlog',linthreshx=1e6)
            p.yscale('symlog',linthreshy=1e6)
            p.tick_params(axis='both', which='both', labelsize=6)
            p.title('k = '+str("%.4f" % kpl_fold[ind]), fontsize=8)

    if opts.plot:
        # plot Pin vs. Pout
        p.figure(1)
        p.tight_layout()
        fig1.text(0.5, 0.0, r'$P_{\rm in}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',
                 ha='center')
        fig1.text(0.0, 0.5, r'$P_{\rm out}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',
                 va='center', rotation='vertical')

        # plot Pin vs. Pout for I case
        p.figure(2)
        p.tight_layout()
        fig2.text(0.5, 0.0, r'$P_{\rm in}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',
                 ha='center')
        fig2.text(0.0, 0.5, r'$P_{\rm out,I}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',
                 va='center', rotation='vertical')
        p.show()

    # Function to get bin sizes to use for normalization
    def bin_size(bins): # XXX approximation since they're just divided in two
        bins_size = []
        for bb in range(len(bins)):
            if bb == 0: bins_size.append(bins[bb+1]-bins[bb]) # first bin
            elif bb == len(bins)-1: bins_size.append(bins[bb]-bins[bb-1]) # last bin
            else: bins_size.append((bins[bb+1]-bins[bb-1])/2) # middle bins
        return bins_size

    # Function to compute probability matrix ('transfer function')
    
    def make_bins():
        # lin-log grid-spacing
        nbins = 100 # number of bins for grid upon which to estimate kernels
        dmax = n.max(Pins_fold.values())
        dg = 1 # spacing in linear regime
        G = n.logspace(0,n.log10(dmax),nbins)[-1]/n.logspace(0,n.log10(dmax),nbins)[-2] # way to estimate multiplicative factor in log regime (~1.3) 
        grid = [] # grid spacing that's linear at low regime and log at high regime 
        g = 0
        for i in range(nbins):
            g = G*g + dg
            grid.append(g)
        bins = n.log10(grid) # log-space
        bins_size = bin_size(bins)
        bins_concat = n.concatenate((-10**bins[::-1],10**bins)) # full length, real numbers
        return bins, bins_concat

    def transfer_func(data, identity=False):
        bins_size = bin_size(bins)
        M_matrix = {}
        for kk,k in enumerate(kpl_fold): # only positive k's
            # Use all bootstraps 
            #xs = n.array(Pins_fold[k]).flatten()
            #if identity == True: ys = n.array(Pouts_I_fold[k]).flatten()
            #if identity == False: ys = n.array(Pouts_fold[k]).flatten()
            # Average signal loss curves first
            xs = n.mean(n.array(Pins_fold[k]),axis=1)
            if identity == True: ys = n.mean(n.array(Pouts_I_fold[k]),axis=1)
            if identity == False: ys = n.mean(n.array(Pouts_fold[k]),axis=1)
            # Binning by kernel density estimators
            ygrid,xgrid = n.meshgrid(bins,bins) # create grid on which to sample... note that ygrid and xgrid are flipped for some reason
            positions = n.vstack([xgrid.ravel(),ygrid.ravel()])
            xs = n.abs(xs) # abs all since transfer curve is symmetric
            ys = n.abs(ys) 
            kernel = scipy.stats.gaussian_kde((n.log10(xs),n.log10(ys)),bw_method=0.1)#'scott') # kernel determined by P_in and P_out
            M_matrix[k] = n.reshape(kernel(positions).T,xgrid.shape).T
            # ensure columns add to 1
            for col in range(M_matrix[k].shape[1]):
                if n.sum(M_matrix[k][:,col]) > 0: # avoid nan values
                    M_matrix[k][:,col] /= n.sum(M_matrix[k][:,col]*bins_size)
            # Plot 2D distribution (KDE)
            if False and kk == 0: # just one k-value
                p.plot(n.log10(xs), n.log10(ys),'k.')
                p.pcolormesh(bins_in,bins_out,M_matrix[k]);p.colorbar()
                p.xlim(0,n.log10(dmax));p.ylim(0,n.log10(dmax))
                p.title('k='+str(k))
                p.xlabel('$P_{in}$ (log)');p.ylabel('$P_{out}$ (log)');p.show()
        return M_matrix
        
    # Function to apply transfer function
    def sigloss_func(data, M_matrix):
        new_PS = {} # dictionaries that will hold PS distributions
        old_PS = {}
        for k in data.keys():
            kernel = scipy.stats.gaussian_kde(data[k][0])
            data_dist = kernel(bins_concat)
            bins_concat_size = bin_size(bins_concat)
            data_dist /= n.sum(data_dist*bins_concat_size) # normalize
            data_dist_pos = data_dist[len(bins):] # separate
            data_dist_neg = data_dist[:len(bins)][::-1] # reverse this to apply same M
            # Get new distribution
            try: M = M_matrix[k]
            except: M = M_matrix[-k] # M only exists for positive k's
            
            Mpos = M*(n.resize(data_dist_pos,M.shape).T) # multiply data distribution element-wise per column of M 
            Mneg = M*(n.resize(data_dist_neg,M.shape).T) 
            rowsum_pos = n.zeros(M.shape[0])
            rowsum_neg = n.zeros(M.shape[0])
            for row in Mpos: rowsum_pos += row # add up rows of M
            for row in Mneg: rowsum_neg += row 
            rowsum_combine = n.concatenate((rowsum_neg[::-1],rowsum_pos))
            new_PS[k] = rowsum_combine/n.sum(rowsum_combine*bins_concat_size) # normalize
            old_PS[k] = data_dist # save original distribution too
        return new_PS, old_PS
   
    # Call function for probability matrix
    bins, bins_concat = make_bins() # bins are where to sample distributions
    M = transfer_func(pCs)
    M_I = transfer_func(pIs, identity=True)
   
    # Call function for new PS distributions
    new_pCs, old_pCs = sigloss_func(pCs, M)
    new_pCs_fold, old_pCs_fold = sigloss_func(pCs_fold, M)
    new_pIs, old_pIs = sigloss_func(pIs, M_I)
    new_pIs_fold, old_pIs_fold = sigloss_func(pIs_fold, M_I)

    # Compute PS points and errors
    pC = []; pC_fold = []; pC_err = []; pC_fold_err = [] # final arrays
    pI = []; pI_fold = []; pI_err = []; pI_fold_err = []
   
    def compute_stats(bins,data):
        pts, errs = [], []
        for key in n.sort(data.keys()):
            point = bins[n.argmax(data[key])]
            pts.append(point)
            percents = [n.sum((data[key]*bin_size(bins_concat))[:i]) for i in range(len(data[key]))]
            left = n.interp(0.025,percents,bins)
            right = n.interp(0.975,percents,bins)
            errs.append((right-left)/2) # 1-sigma
        return pts,errs

    pC, pC_err = compute_stats(bins_concat, new_pCs)
    pI, pI_err = compute_stats(bins_concat, new_pIs)
    pC_fold, pC_fold_err = compute_stats(bins_concat, new_pCs_fold)
    pI_fold, pI_fold_err = compute_stats(bins_concat, new_pIs_fold)
    
    # Plot old vs. new distributions for P(k)
    if False:
        num = 1 # weighted
        p.figure(figsize=(15,8))
        for nn in range(len(file['kpl'])):
            p.subplot(3,7,num)
            p.plot(bin_cent,old_pCs[file['kpl'][nn]]/n.max(old_pCs[file['kpl'][nn]]),'b-',label='old')
            p.plot(bin_cent,new_pCs[file['kpl'][nn]]/n.max(new_pCs[file['kpl'][nn]]),'g-',label='new')
            p.xscale('symlog')
            p.title('k = ' + str(file['kpl'][nn]),fontsize=8)
            p.tick_params(axis='both', which='major', labelsize=6)
            p.tick_params(axis='both', which='minor', labelsize=6)
            p.grid()
            num += 1
        p.suptitle('Distribution of Weighted PS, Old (blue) vs. New (green)')
        p.show()
        num = 1 # unweighted
        p.figure(figsize=(15,8))
        for nn in range(len(file['kpl'])):
            p.subplot(3,7,num)
            p.plot(bin_cent_I,old_pIs[file['kpl'][nn]]/n.max(old_pIs[file['kpl'][nn]]),'b-',label='old')
            p.plot(bin_cent_I,new_pIs[file['kpl'][nn]]/n.max(new_pIs[file['kpl'][nn]]),'g-',label='new')
            p.xscale('symlog')
            p.title('k = ' + str(file['kpl'][nn]),fontsize=8)
            p.tick_params(axis='both', which='major', labelsize=6)
            p.tick_params(axis='both', which='minor', labelsize=6)
            p.grid()
            num += 1
        p.suptitle('Distribution of Unweighted PS, Old (blue) vs. New (green)')
        p.show()

    # Save values
    if opts.skip_sigloss:
        print 'Skipping applying signal loss!!!'
        if count == 0: # data case
            pCv = file['pCv']
            pCv_fold = file['pCv_fold']
            pIv = file['pIv']
            pIv_fold = file['pIv_fold']
            pCv_err = file['pCv_err']
            pCv_fold_err = file['pCv_fold_err']
            pIv_err = file['pIv_err']
            pIv_fold_err = file['pIv_fold_err']
        else: # noise case
            pCn = file['pCn']
            pCn_fold = file['pCn_fold']
            pIn = file['pIn']
            pIn_fold = file['pIn_fold']
            pCn_err = file['pCn_err']
            pCn_fold_err = file['pCn_fold_err']
            pIn_err = file['pIn_err']
            pIn_fold_err = file['pIn_fold_err']
        
    if count == 0: # data case
        pCv = pCv    
        pCv_fold = pCv_fold
        pIv = pIv
        pIv_fold = pIv_fold
        pCv_err = pCv_err
        pCv_fold_err = pCv_fold_err
        pIv_err = pIv_err
        pIv_fold_err = pIv_fold_err
    if count == 1: # noise case
        pCn = pCn
        pCn_fold = pCn_fold
        pIn = pIn
        pIn_fold = pIn_fold
        pCn_err = pCn_err
        pCn_fold_err = pCn_fold_err
        pIn_err = pIn_err
        pIn_fold_err = pIn_fold_err

# Write out solutions
outname = 'pspec_final_sep'+opts.sep+'.npz'
print '   Saving', outname

# Save with 1-sigma errors (compatible with pspec_combine.py and plot_pspec_final.py)
n.savez(outname, kpl=kpl, k=file['k'], freq=file['freq'],
        pCv=pCv, pCv_err=pCv_err,
        pCv_fold=pCv_fold, pCv_fold_err=pCv_fold_err,
        pIv=pIv, pIv_err=pIv_err,
        pIv_fold=pIv_fold, pIv_fold_err=pIv_fold_err,
        pCn=pCn, pCn_err=pCn_err,
        pCn_fold=pCn_fold, pCn_fold_err=pCn_fold_err,
        pIn=pIn, pIn_err=pIn_err,
        pIn_fold=pIn_fold, pIn_fold_err=pIn_fold_err,
        prob=0.95,
        ngps=file['ngps'], nbls=file['nbls'], nbls_g=file['nbls_g'], nlsts_g=file['nlsts_g'],
        lsts=file['lsts'], afreqs=file['afreqs'], cnt_eff=file['cnt_eff'],
        frf_inttime=file['frf_inttime'], inttime=file['inttime'], 
        cmd=file['cmd'].item() + ' \n '+' '.join(sys.argv))

