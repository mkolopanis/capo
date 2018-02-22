#! /usr/bin/env python
"""Compute Signal Loss factor via distribution method v5.
    -Reads in pspec_2d_to_1d.npz (contains PS points for aln bootstraps)
    -Uses KDE to smooth out 2D transfer curve (P_in vs. P_out)
    -Takes a horizontal cut at P_out = peak of p_x (data) and sums to 1
    -Uses this to compute new PS points + errors
"""
import numpy as n
import pylab as p
from matplotlib import gridspec
from scipy.interpolate import interp1d
from scipy.stats import norm
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
o.add_option('--nosubtract',action='store_true',
            help='Run for no-subtraction case (P_out = pCr).')
opts,args = o.parse_args(sys.argv[1:])

#-------------
### FUNCTIONS

# Get bin sizes to use for normalization
def bin_size(bins): # XXX approximation since they're just divided in two
    bins_size = []
    for bb in range(len(bins)):
        if bb == 0: bins_size.append(bins[bb+1]-bins[bb]) # first bin
        elif bb == len(bins)-1: bins_size.append(bins[bb]-bins[bb-1]) # last bin
        else: bins_size.append((bins[bb+1]-bins[bb-1])/2) # middle bins
    return bins_size

# Make bins (grid), which is used for sampling distributions 
def make_bins():
    # lin-log grid-spacing
    nbins = 101 # XXX hard-coded number of bins for grid upon which to estimate kernels (must be odd for fftshift to do the right thing)
    dmax = max(n.max(Pins_fold.values()),n.max(pIs.values())) # maximum is determined from whichever is largest: EoR injection or "I" value
    dg = n.min(Pins_fold.values()) #1 # XXX artificially set based on starting P_in values
    """
    G = n.logspace(0,n.log10(dmax),nbins)[-1]/n.logspace(0,n.log10(dmax),nbins)[-2] # way to estimate multiplicative factor in log regime (~1.3) 
    grid = [] # grid spacing that's linear at low regime and log at high regime 
    g = 0
    for i in range(nbins): # XXX nothing is constraining this from going beyond dmax
        g = G*g + dg
        grid.append(g)
    """
    grid = n.logspace(n.log10(dg),n.log10(dmax),nbins) # logspace sampling
    bins = n.log10(grid) # log-space
    bins_concat = n.concatenate((-10**bins[::-1],10**bins)) # full length, real numbers
    binsx, binsy = bins, bins
    binsy_full = n.concatenate((-binsy[::-1],binsy)) # for pos and neg
    return binsx, binsy, binsy_full, bins_concat # bins is n.log10, bins_concat is not

# For all data points (full distributions), smooth the transfer curves for both the C and I cases using KDEs
def smooth_dist(fold=True):
    kde_C = {}
    kde_I = {}
    if fold == True: ks = kpl_fold # k values to loop over
    if fold == False: ks = kpl
    for kk,k in enumerate(ks): 

        if kk == 8 and count == 0 and fold == True: # save values out
            # Save values to use for plotting sigloss terms
            if opts.nosubtract: fn = 'pspec_sigloss_terms_nosubtract.npz'
            else: fn = 'pspec_sigloss_terms.npz'
            print "Saving",fn,"which contains data values for k =",k
            n.savez(fn, k=k, pCv=file['pCv'][n.where(kpl==k)[0][0]], pIv=file['pIv'][n.where(kpl==k)[0][0]], Pins=Pins_fold[k], Pouts=Pouts_fold[k], Pouts_I=Pouts_I_fold[k], pCrs_fold=pCrs_fold[k], pCvs_Cr_fold=pCvs_Cr_fold[k], pCes_Cr_fold=pCes_Cr_fold[k], pCves_fold=pCves_fold[k], pIves_fold=pIves_fold[k])
        if fold == True:
            xs = n.array(Pins_fold[k]).flatten()
            ys_C = n.array(Pouts_fold[k]).flatten()
            ys_I = n.array(Pouts_I_fold[k]).flatten()
        if fold == False:
            xs = n.array(Pins[k]).flatten()
            ys_C = n.array(Pouts[k]).flatten()
            ys_I = n.array(Pouts_I[k]).flatten()
        # Kernel Density Estimator (all values together)
        ygrid,xgrid = n.meshgrid(binsy_full,binsy_full)
        positions = n.vstack([xgrid.ravel(),ygrid.ravel()])
        ind = n.where(xs > 0)[0] # only positive P_in's
        xs = xs[ind] 
        ys_C = ys_C[ind]
        ys_I = ys_I[ind]
        factor = 0.05
        factorI = 0.05
        kernel_C = scipy.stats.gaussian_kde((n.log10(xs),n.sign(ys_C)*n.log10(n.abs(ys_C))),bw_method=factor)
        kernel_I = scipy.stats.gaussian_kde((n.log10(xs),n.sign(ys_I)*n.log10(n.abs(ys_I))),bw_method = factorI)
        kdeC = n.reshape(kernel_C(positions).T,(binsy_full.size,binsy_full.size)).T
        kdeI = n.reshape(kernel_I(positions).T,(binsy_full.size,binsy_full.size)).T
        kdeC = kdeC[:,binsx.size:] # remove negative P_in half (nothing is there)
        kdeI = kdeI[:,binsx.size:]
        # Normalize columns
        for col in range(kdeC.shape[1]):
            if n.sum(kdeC[:,col]) > 0: # avoid nan values
                kdeC[:,col] /= n.sum(kdeC[:,col]*bin_size(binsy_full))
            if n.sum(kdeI[:,col]) > 0:
                kdeI[:,col] /= n.sum(kdeI[:,col]*bin_size(binsy_full))
        # Plot KDE and points
        if count == 1: #opts.plot:
            p.figure(figsize=(10,6))
            p.subplot(121)
            p.pcolormesh(binsx,binsy_full,kdeC,cmap='hot_r')
            if fold == True: dataval = file['pCv_fold'][n.where(ks==k)[0][0]]
            if fold == False: dataval = file['pCv'][n.where(ks==k)[0][0]]
            p.axhline(y=n.sign(dataval)*n.log10(n.abs(dataval)),color='0.5',linewidth=2)
            p.plot(n.log10(xs), n.sign(ys_C)*n.log10(n.abs(ys_C)),'k.')
            p.xlim(binsx.min(), binsx.max())
            p.ylim(binsy_full.min(), binsy_full.max()); p.grid()
            p.subplot(122)
            p.pcolormesh(binsx,binsy_full,kdeI,cmap='hot_r')
            if fold == True: dataval = file['pIv_fold'][n.where(ks==k)[0][0]]
            if fold == False: dataval = file['pIv'][n.where(ks==k)[0][0]]
            p.axhline(y=n.sign(dataval)*n.log10(n.abs(dataval)),color='0.5',linewidth=2)
            p.plot(n.log10(xs), n.sign(ys_I)*n.log10(n.abs(ys_I)),'k.')
            p.xlim(binsx.min(), binsx.max())
            p.ylim(binsy_full.min(), binsy_full.max())
            p.tight_layout(); p.suptitle("k = " + str(k)); p.grid(); p.show()
        kde_C[k] = kdeC # populate dictionary
        kde_I[k] = kdeI
    return kde_C, kde_I

# Get data distribution using 1D KDE
def data_dist(data):
    old_PS = {}
    for k in data.keys():
        kernel = scipy.stats.gaussian_kde(data[k][0]) # KDE for data
        data_dist = kernel(bins_concat)
        bins_concat_size = bin_size(bins_concat)
        data_dist /= n.sum(data_dist*bins_concat_size) # normalize
        old_PS[k] = data_dist
    return old_PS

# Apply transfer function to data
def sigloss_func(pt, M_matrix):
    new_PS = {} # dictionaries that will hold PS distributions
    for k in pt.keys():
        point = pt[k] # value of data point
        M = M_matrix[k].copy()
        peak_ind = n.argmin(n.abs(bins_concat-point)) # matching index
        cut = M[peak_ind] # cut of KDE at peak of data
        #if peak_ind < len(bins_concat)/2: cut = M[len(bins_concat)-peak_ind-1] # XXX positive cut instead
        cut /= n.sum(cut*bin_size(10**binsx))
        cut = n.concatenate((n.zeros_like(cut),cut)) # tack on zeros as negative half (always positive points)
        new_PS[k] = cut
    return new_PS # distribution of bins_concat
 
# Compute PS points and errors
def compute_stats(bins, data, pt):
    pts, errs = [], []
    for key in n.sort(data.keys()):
        if opts.skip_sigloss: point = pt[key] # use no-bootstrapping case
        else: 
            #point = bins[n.argmax(data[key])] # XXX peak of dist
            point = n.sum(bins*data[key])/n.sum(data[key]) # weighted avg
        pts.append(point)
        percents = [n.sum((data[key]*bin_size(bins))[:i]) for i in range(len(data[key]))]
        left = n.interp(0.025,percents,bins)
        right = n.interp(0.975,percents,bins)
        errs.append((right-left)/4) # 1-sigma
    return pts,errs

# Avoid signal gain (i.e. avoid reducing PS values artificially because the transfer curve artificially shows 'signal gain')
def no_gain(pts, errs, old_dist, pt):
    new_errs = []
    new_pts = []
    _, old_errs = compute_stats(bins_concat, old_dist, pt)
    for ee in range(len(errs)):
        if errs[ee] < old_errs[ee]: # if variance after signal loss is smaller than before
            new_errs.append(old_errs[ee])
            try: new_pts.append(n.abs(pt[kpl[ee]])) # un-folded case # XXX abs since points are always positive now?
            except: new_pts.append(n.abs(pt[-kpl[ee]])) # folded case (positive k's only)
        else: 
            new_errs.append(errs[ee])
            new_pts.append(pts[ee])
    return new_pts, new_errs

#-----------------
### Main Program

# read pspec_pk_k3pk.npz file for k points and PS points (from no-boostrapping case)
file = n.load(glob.glob('inject_sep'+opts.sep+'*')[0]+'/pspec_pk_k3pk.npz')
kpl = file['kpl']; k = file['k']; kpl_fold = file['kpl_fold']

# loop over injects and read pspec_2d_to_1d outputs
for count in range(2):
    if count == 1: print 'NOISE CASE'
    else: print 'DATA CASE'
    pklo, pkhi = 1e-1, 1e12
    Pouts_fold = {} # folded case
    Pouts_I_fold = {}
    Pins_fold = {}
    Pouts = {} # un-folded
    Pouts_I = {}
    Pins = {}
    pCs = {}; pCs_fold = {}
    pIs = {}; pIs_fold = {}
    pCes_Cr_fold = {}
    pCvs_Cr_fold = {}
    pCves_fold = {}
    pIves_fold = {}
    pCrs_fold = {}
    pIrs_fold = {}

    for inject in glob.glob('inject_sep'+opts.sep+'*'):
        print 'Reading', inject
        file_2d = n.load(inject + '/pspec_2d_to_1d.npz')

        if count == 1: # noise case
            color='b.'
            linecolor='cyan'
            Pout_fold = file_2d['pCs-pCn_fold']
            Pout_I_fold = file_2d['pIs-pIn_fold']
            Pout = file_2d['pCs-pCn']
            Pout_I = file_2d['pIs-pIn']
            pC = file_2d['pCn']
            pC_fold = file_2d['pCn_fold']
            pI = file_2d['pIn']
            pI_fold = file_2d['pIn_fold']
        else:
            color='k.'
            linecolor='0.5'
            Pout_fold = file_2d['pCr-pCv_fold'] # shape (#boots, #kpl)
            Pout_I_fold = file_2d['pIr-pIv_fold']  # pI case
            Pout = file_2d['pCr-pCv']
            Pout_I = file_2d['pIr-pIv']
            if opts.nosubtract:
                Pout_fold = file_2d['pCr_fold']
                Pout_I_fold = file_2d['pIr_fold']
                Pout = file_2d['pCr']
                Pout_I = file_2d['pIr']
            pCe_Cr_fold = file_2d['pCe_Cr_fold']
            pCv_Cr_fold = file_2d['pCv_Cr_fold']
            try: 
                pCve_fold = file_2d['pCve_fold']
                pIve_fold = file_2d['pIve_fold']
            except: # depending on how recent the oqe run was, this key may not exist
                pCve_fold = n.zeros_like(pCv_Cr_fold)
                pIve_fold = n.zeros_like(pCv_Cr_fold)
            pCr_fold = file_2d['pCr_fold']
            pIr_fold = file_2d['pIr_fold']
            pC = file_2d['pCv']
            pC_fold = file_2d['pCv_fold']
            pI = file_2d['pIv']
            pI_fold = file_2d['pIv_fold']

        Pin_fold = file_2d['pIe_fold']
        Pin = file_2d['pIe']
        
        for ind in range(len(kpl_fold)): # loop through k for Delta^2(k)
            try:
                Pouts_fold[kpl_fold[ind]].append(Pout_fold[:,ind])
                Pouts_I_fold[kpl_fold[ind]].append(Pout_I_fold[:,ind])
                Pins_fold[kpl_fold[ind]].append(Pin_fold[:,ind])
                pCs_fold[kpl_fold[ind]] = [pC_fold[:,ind]] # no appending because it's the same for every inject
                pIs_fold[kpl_fold[ind]] = [pI_fold[:,ind]]
                pCes_Cr_fold[kpl_fold[ind]].append(pCe_Cr_fold[:,ind])
                pCvs_Cr_fold[kpl_fold[ind]].append(pCv_Cr_fold[:,ind])
                pCves_fold[kpl_fold[ind]].append(pCve_fold[:,ind])
                pIves_fold[kpl_fold[ind]].append(pIve_fold[:,ind])
                pCrs_fold[kpl_fold[ind]].append(pCr_fold[:,ind])
                pIrs_fold[kpl_fold[ind]].append(pIr_fold[:,ind])
            except:
                Pouts_fold[kpl_fold[ind]] = [Pout_fold[:,ind]]
                Pouts_I_fold[kpl_fold[ind]] = [Pout_I_fold[:,ind]]
                Pins_fold[kpl_fold[ind]] = [Pin_fold[:,ind]]
                pCs_fold[kpl_fold[ind]] = [pC_fold[:,ind]]
                pIs_fold[kpl_fold[ind]] = [pI_fold[:,ind]]
                pCes_Cr_fold[kpl_fold[ind]] = [pCe_Cr_fold[:,ind]]
                pCvs_Cr_fold[kpl_fold[ind]] = [pCv_Cr_fold[:,ind]]
                pCves_fold[kpl_fold[ind]] = [pCve_fold[:,ind]]
                pIves_fold[kpl_fold[ind]] = [pIve_fold[:,ind]]
                pCrs_fold[kpl_fold[ind]] = [pCr_fold[:,ind]]
                pIrs_fold[kpl_fold[ind]] = [pIr_fold[:,ind]]
        for ind in range(len(kpl)): # loop through k for P(k)
            pCs[kpl[ind]] = [pC[:,ind]] # no appending because it's the same for every inject
            pIs[kpl[ind]] = [pI[:,ind]]
            try:
                Pouts[kpl[ind]].append(Pout[:,ind])
                Pouts_I[kpl[ind]].append(Pout_I[:,ind])
                Pins[kpl[ind]].append(Pin[:,ind])
            except:
                Pouts[kpl[ind]] = [Pout[:,ind]]
                Pouts_I[kpl[ind]] = [Pout_I[:,ind]]
                Pins[kpl[ind]] = [Pin[:,ind]]
                
    # Values of interest are contained within dictionaries indexed by kpl_fold or kpl:
    #     Pins_fold, Pins
    #     Pouts_fold, Pouts
    #     Pouts_I_fold, Pouts_I
    #     pCs, pCs_fold, pIs, pIs_fold  (data/noise values)
 
    ### SIGNAL LOSS CODE

    binsx, binsy, binsy_full, bins_concat = make_bins() # bins are where to sample distributions

    # Get distributions of original data (used to find errors if skip_sigloss)
    old_pCs = data_dist(pCs)
    old_pCs_fold = data_dist(pCs_fold)
    old_pIs = data_dist(pIs)    
    old_pIs_fold = data_dist(pIs_fold)

    # Get original PS points (from no-bootstrapping case)
    if count == 0: # data case
        pt_pC = file['pCv']
        pt_pC_fold = file['pCv_fold']
        pt_pI = file['pIv']
        pt_pI_fold = file['pIv_fold']
    else: # noise case
        pt_pC = file['pCn']
        pt_pC_fold = file['pCn_fold']
        pt_pI = file['pIn']
        pt_pI_fold = file['pIn_fold']
    pt_pCs, pt_pCs_fold, pt_pIs, pt_pIs_fold = {},{},{},{}
    for kk,k in enumerate(kpl):
        pt_pCs[k] = pt_pC[kk]
        pt_pIs[k] = pt_pI[kk]
    for kk,k in enumerate(kpl_fold):
        pt_pCs_fold[k] = pt_pC_fold[kk]
        pt_pIs_fold[k] = pt_pI_fold[kk]
    
    # Compute signal loss and get new distributions 
    if opts.skip_sigloss:
        print "Skipping Signal Loss!"
        new_pCs = old_pCs.copy()
        new_pCs_fold = old_pCs_fold.copy()
        new_pIs = old_pIs.copy()
        new_pIs_fold = old_pIs_fold.copy()
    else:
        kde_C, kde_I = smooth_dist(fold=False) # KDE transfer curves
        kde_C_fold, kde_I_fold = smooth_dist(fold=True)
        # new distributions
        new_pCs = sigloss_func(pt_pCs, kde_C) 
        new_pCs_fold = sigloss_func(pt_pCs_fold, kde_C_fold)
        new_pIs = sigloss_func(pt_pIs, kde_I)
        new_pIs_fold = sigloss_func(pt_pIs_fold, kde_I_fold)
   
    pC, pC_err = compute_stats(bins_concat, new_pCs, pt_pCs)
    pI, pI_err = compute_stats(bins_concat, new_pIs, pt_pIs)
    pC_fold, pC_fold_err = compute_stats(bins_concat, new_pCs_fold, pt_pCs_fold)
    pI_fold, pI_fold_err = compute_stats(bins_concat, new_pIs_fold, pt_pIs_fold)
    """
    if opts.skip_sigloss == None: # XXX artificially disallow signal gain
        pC, pC_err = no_gain(pC, pC_err, old_pCs, pt_pCs)
        pI, pI_err = no_gain(pI, pI_err, old_pIs, pt_pIs)
        pC_fold, pC_fold_err = no_gain(pC_fold, pC_fold_err, old_pCs_fold, pt_pCs_fold)
        pI_fold, pI_fold_err = no_gain(pI_fold, pI_fold_err, old_pIs_fold, pt_pIs_fold)
    """
    # Save values to use for plotting sigloss plots
    if count == 0:
        Ind = -3 # one k-value
        k = kpl[ind]
        if opts.nosubtract: fn = 'pspec_sigloss_nosubtract.npz'
        else: fn = 'pspec_sigloss.npz'
        print "Saving",fn,"which contains data values for k =",k
        n.savez(fn, k=k, binsx=binsx, binsy=binsy, bins_concat=bins_concat, pC=pC[ind], pC_err=pC_err[ind], pI=pI[ind], pI_err=pI_err[ind], new_pCs=new_pCs[k], new_pIs=new_pIs[k], old_pCs=old_pCs[k], old_pIs=old_pIs[k], Pins=Pins_fold[k], Pouts=Pouts_fold[k], Pouts_I=Pouts_I_fold[k])    
       
    if count == 0: # data case
        pCv = pC    
        pCv_fold = pC_fold
        pIv = pI
        pIv_fold = pI_fold
        pCv_err = pC_err
        pCv_fold_err = pC_fold_err
        pIv_err = pI_err
        pIv_fold_err = pI_fold_err
    if count == 1: # noise case
        pCn = pC
        pCn_fold = pC_fold
        pIn = pI
        pIn_fold = pI_fold
        pCn_err = pC_err
        pCn_fold_err = pC_fold_err
        pIn_err = pI_err
        pIn_fold_err = pI_fold_err

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
        prob=0.95, kperp=file['kperp'], sep=opts.sep, kpl_fold=file['kpl_fold'],
        ngps=file['ngps'], nbls=file['nbls'], nbls_g=file['nbls_g'], nlsts_g=file['nlsts_g'],
        lsts=file['lsts'], afreqs=file['afreqs'], cnt_eff=file['cnt_eff'],
        frf_inttime=file['frf_inttime'], inttime=file['inttime'], 
        cmd=file['cmd'].item() + ' \n '+' '.join(sys.argv))

