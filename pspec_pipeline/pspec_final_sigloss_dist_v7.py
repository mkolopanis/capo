#! /usr/bin/env python
"""Compute Signal Loss factor via distribution method v6.
    -Reads in pspec_2d_to_1d.npz (contains PS points for aln bootstraps)
    -Assumes dense sampling in P_in!
    -For bins along the P_in axis, a 1D KDE is used to smooth out the points (along the P_out axis)
    -All KDE's are combined into a 2D transfer curve (P_in vs. P_r)
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
o.add_option('--sep',default='0_2',
            help='Separation to run script for.')
o.add_option('--nosubtract',action='store_true',
            help='Run for no-subtraction case (P_out = pCr).')
opts,args = o.parse_args(sys.argv[1:])

#-------------
### FUNCTIONS

# Calculate normalization factor for log-sampling prior
# Factor is linear width
def prior_factor(bins):
    factors = []
    for bb in range(len(bins)):
        if bb != 0 and bb != len(bins)-1: # skip edge cases
            factor = (bins[bb+1]-bins[bb-1]) #/ (n.log10(bins[bb+1]/bins[bb-1]))
            factors.append(factor)
    factors.append(factors[-1]) # repeat last value for right edge case
    factors = n.insert(factors,0,factors[0]) # repeat first value for left edge case
    return factors

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
    nbins = 101 # XXX hard-coded number of bins for grid upon which to estimate kernels
    dmax = max(n.max(Pins_fold.values()),n.max(pIs.values())) # maximum is determined from whichever is largest: EoR injection or "I" value
    dg = n.min(Pins_fold.values()) #1 # XXX artificially set based on starting P_in values
    grid = n.logspace(n.log10(dg),n.log10(dmax),nbins) # logspace sampling
    binsy_log = n.log10(grid) # log-space
    binsy_lin = n.concatenate((-10**binsy_log[::-1],10**binsy_log)) # real numbers (not log-space)
    binsy_log = n.concatenate((-binsy_log[::-1],binsy_log)) # for pos and neg
    return binsy_log,binsy_lin

# For all data points (full distributions), smooth the transfer curves for both the C and I cases using 1D KDE's per injection level
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

        # Get bins and data for transfer curve
        bins_inj = binsy_log[len(binsy_log)/2:] # positive half of binsy_log will be used as bin boundaries for P_in
        binsy_kde = binsy_log
        binsx_kde = []
        for bb in range(len(bins_inj)-1): # find center P_in bins
            binsx_kde.append((bins_inj[bb+1]+bins_inj[bb])/2)
        binsx_kde = n.insert(binsx_kde,0,binsx_kde[0]-(binsx_kde[1]-binsx_kde[0]))
        binsx_kde = n.insert(binsx_kde,len(binsx_kde),binsx_kde[-1] + (binsx_kde[-1]-binsx_kde[-2]))
        if fold == True:
            xs = n.array(Pins_fold[k])
            ys_C = n.array(Pouts_fold[k])
            ys_I = n.array(Pouts_I_fold[k])
        if fold == False:
            xs = n.array(Pins[k])
            ys_C = n.array(Pouts[k])
            ys_I = n.array(Pouts_I[k])
        xs_kde = (n.log10(xs)).flatten() # log-space
        ys_kde = (n.sign(ys_C)*n.log10(n.abs(ys_C))).flatten()
        ys_I_kde = (n.sign(ys_I)*n.log10(n.abs(ys_I))).flatten()
        binning = n.digitize(xs_kde,bins_inj)
        kdeC = n.zeros((len(binsy_kde),len(binsx_kde)))
        kdeI = n.zeros((len(binsy_kde),len(binsx_kde)))
        for sub_bin in range(len(binsx_kde)): # loop over injection bins
            # get x and y values for each sub_bin
            xs_sub_bin = xs_kde[n.where(binning == sub_bin)[0]] 
            ys_sub_bin = ys_kde[n.where(binning == sub_bin)[0]]
            ys_I_sub_bin = ys_I_kde[n.where(binning == sub_bin)[0]]
            # find distribution for weighted case
            if len(ys_sub_bin) == 1: # if only one value, use delta function
                delta_loc = n.argmin(n.abs(binsy_kde - ys_sub_bin))
                data_dist_C = n.zeros_like(binsy_kde)
                data_dist_C[delta_loc] = 1.0
            elif len(ys_sub_bin) == 0: # if no values, leave as zeros
                data_dist_C = n.zeros_like(binsy_kde)
            else: # if multiple values, use KDE
                kernel_C = scipy.stats.gaussian_kde(ys_sub_bin) 
                ratio = n.sqrt(n.cov(n.abs(ys_sub_bin)) / n.cov(ys_sub_bin))
                factor = ratio * kernel_C.factor 
                kernel_C = scipy.stats.gaussian_kde(ys_sub_bin,bw_method=factor)
                data_dist_C = kernel_C(binsy_kde)
            # find distribution for unweighted case
            if len(ys_I_sub_bin) == 1: 
                delta_loc = n.argmin(n.abs(binsy_kde - ys_I_sub_bin))
                data_dist_I = n.zeros_like(binsy_kde)
                data_dist_I[delta_loc] = 1.0    
            elif len(ys_I_sub_bin) == 0:
                data_dist_I = n.zeros_like(binsy_kde)
            else:
                kernel_I = scipy.stats.gaussian_kde(ys_I_sub_bin) 
                ratioI = n.sqrt(n.cov(n.abs(ys_I_sub_bin)) / n.cov(ys_I_sub_bin))
                factorI = ratioI * kernel_I.factor 
                kernel_I = scipy.stats.gaussian_kde(ys_I_sub_bin,bw_method=factorI)
                data_dist_I = kernel_I(binsy_kde)
            # Normalize and save distribution as a column of the transfer matrix
            if n.sum(data_dist_C) != 0: data_dist_C /= n.sum(data_dist_C*bin_size(binsy_kde))
            if n.sum(data_dist_I) != 0: data_dist_I /= n.sum(data_dist_I*bin_size(binsy_kde))
            kdeC[:,sub_bin] = data_dist_C
            kdeI[:,sub_bin] = data_dist_I
        
        # Save values to use for plotting sigloss plots
        if count == 0 and fold == True and kk == 8:
            if opts.nosubtract: fn = 'pspec_sigloss_nosubtract.npz'
            else: fn = 'pspec_sigloss.npz'
            print "Saving",fn,"which contains data values for k =",k
            n.savez(fn, k=k, binsx=binsx_kde, binsy=binsy_kde, kdeC=kdeC, kdeI=kdeI, xs=xs_kde, ys=ys_kde, ysI=ys_I_kde, pC=file['pCv_fold'][n.where(ks==k)[0][0]], pI=file['pIv_fold'][n.where(ks==k)[0][0]])
        
        # Plot KDE and points
        if opts.plot:
            p.figure(figsize=(10,6))
            p.subplot(121)
            p.pcolormesh(binsx_kde,binsy_kde,kdeC,cmap='hot_r')
            if fold == True and count == 0: dataval = file['pCv_fold'][n.where(ks==k)[0][0]]
            if fold == False and count == 0: dataval = file['pCv'][n.where(ks==k)[0][0]]
            if fold == True and count == 1: dataval = file['pCn_fold'][n.where(ks==k)[0][0]]
            if  fold == False and count == 1: dataval = file['pCn'][n.where(ks==k)[0][0]]
            p.axhline(y=n.sign(dataval)*n.log10(n.abs(dataval)),color='0.5',linewidth=2)
            p.plot(xs_kde, ys_kde,'k.')
            p.xlim(binsx_kde.min(), binsx_kde.max())
            p.ylim(binsy_kde.min(), binsy_kde.max()); p.grid()
            p.subplot(122)
            p.pcolormesh(binsx_kde,binsy_kde,kdeI,cmap='hot_r')
            if fold == True and count == 0: dataval = file['pIv_fold'][n.where(ks==k)[0][0]]
            if fold == False and count == 0: dataval = file['pIv'][n.where(ks==k)[0][0]]
            if fold == True and count == 1: dataval = file['pIn_fold'][n.where(ks==k)[0][0]]
            if fold == False and count == 1: dataval = file['pIn'][n.where(ks==k)[0][0]]
            p.axhline(y=n.sign(dataval)*n.log10(n.abs(dataval)),color='0.5',linewidth=2)
            p.plot(xs_kde,ys_I_kde,'k.')
            p.xlim(binsx_kde.min(), binsx_kde.max())
            p.ylim(binsy_kde.min(), binsy_kde.max())
            p.tight_layout(); p.suptitle("k = " + str(k)); p.grid(); p.show()
        kde_C[k] = kdeC # populate dictionary
        kde_I[k] = kdeI
    return kde_C, kde_I, binsx_kde # return KDE's and P_in bins

# Get data distribution using 1D KDE
def data_dist(data):
    old_PS = {}
    for k in data.keys():
        kernel = scipy.stats.gaussian_kde(data[k][0]) # KDE for data
        data_dist = kernel(binsy_lin)
        bins_size = bin_size(binsy_lin)
        data_dist /= n.sum(data_dist*bins_size) # normalize
        old_PS[k] = data_dist
    return old_PS

# Apply transfer function to data
def sigloss_func(pt, M_matrix):
    new_PS = {} # dictionaries that will hold PS distributions
    for k in pt.keys():
        point = pt[k] # value of data point
        M = M_matrix[k].copy()
        peak_ind = n.argmin(n.abs(binsy_lin-point)) # matching index
        cut = M[peak_ind] # cut of KDE at peak of data
        #if peak_ind < len(binsy_lin)/2: cut = M[len(binsy_lin)-peak_ind-1] # XXX positive cut instead
        cut /= n.sum(cut*bin_size(10**binsx_log)) # normal normalization
        cut * prior_factor(10**binsx_log) # normalize for prior
        new_PS[k] = cut
    return new_PS # distribution of bins_concat

# Compute PS points and errors
def compute_stats(bins, data, pt, old=False):
    pts, errs = [], []
    for key in n.sort(data.keys()):
        if opts.skip_sigloss or old==True:
            point = pt[key] # use no-bootstrapping case
        else:
            #point = bins[n.argmax(data[key])] # XXX peak of dist
            point = n.sum(bins*data[key])/n.sum(data[key]) # weighted avg
        pts.append(point)
        percents = [n.sum((data[key]*bin_size(bins))[:i]) for i in range(len(data[key]))]
        if opts.skip_sigloss or old==True: # find error when dist is centered at PS pt
            left = n.interp(0.025,percents,bins)
            right = n.interp(0.975,percents,bins)
            errs.append((right-left)/4) # 1-sigma error bar
        else:
            up = n.interp(0.975,percents,bins) # upper limit (not error bar) of positive-only posterior
            errs.append(up/2) # divide because it'll be multiplied in the plotting script
    return pts,errs

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

    inj_files = glob.glob('inject_sep'+opts.sep+'*')
    inj_levels = [float(inj.split('_')[-1]) for inj in inj_files]
    inj_inds = n.argsort(inj_levels)
    inj_files = n.take(inj_files, inj_inds)

    for inject in inj_files:
        #if 'moredense' in inject: continue
        print 'Reading', inject
        file_2d = n.load(inject + '/pspec_2d_to_1d.npz')

        if count == 1: # noise case
            color='b.'
            linecolor='cyan'
            Pout_fold = file_2d['pCs-pCn_fold']
            Pout_I_fold = file_2d['pIs-pIn_fold']
            Pout = file_2d['pCs-pCn']
            Pout_I = file_2d['pIs-pIn']
            if opts.nosubtract:
                Pout_fold = file_2d['pCs_fold']
                Pout_I_fold = file_2d['pIs_fold']
                Pout = file_2d['pCs']
                Pout_I = file_2d['pIs']
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

        # estimate variance of inject
        varz = []
        for k_val in range(Pin.shape[1]):
            dat = n.log10(Pin[:,k_val])
            var = n.std(dat)**2
            varz.append(var)
            #print 'k:',kpl[k_val],' ; var:',var
        var_inject = n.mean(varz) # mean over k (should be the same for all injects if seed is the same)

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
    binsy_log, binsy_lin = make_bins() # bins are where to sample distributions
    
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
        binsx_log = binsy_log
    else:
        kde_C, kde_I, binsx_log = smooth_dist(fold=False) # KDE transfer curves
        kde_C_fold, kde_I_fold, binsx_log = smooth_dist(fold=True)
        # new distributions
        new_pCs = sigloss_func(pt_pCs, kde_C)
        new_pCs_fold = sigloss_func(pt_pCs_fold, kde_C_fold)
        new_pIs = sigloss_func(pt_pIs, kde_I)
        new_pIs_fold = sigloss_func(pt_pIs_fold, kde_I_fold)
    
    # Plot un-folded case
    if opts.plot:
        for k in kde_C:
            p.figure(figsize=(10,10))
            xmin,xmax = 0,16
            p.subplot(221)
            p.pcolormesh(binsx_log,binsy_log,kde_C[k],cmap='hot_r')
            if count == 0: dataval = file['pCv'][n.where(kpl==k)[0][0]]
            if count == 1: dataval = file['pCn'][n.where(kpl==k)[0][0]]
            p.axhline(y=n.sign(dataval)*n.log10(n.abs(dataval)),color='0.5',linewidth=2)
            p.xlim(xmin,xmax)
            p.ylim(binsy_log.min(), binsy_log.max()); p.grid()
            p.subplot(222)
            p.pcolormesh(binsx_log,binsy_log,kde_I[k],cmap='hot_r')
            if count == 0: dataval = file['pIv'][n.where(kpl==k)[0][0]]
            if count == 1: dataval = file['pIn'][n.where(kpl==k)[0][0]]
            p.axhline(y=n.sign(dataval)*n.log10(n.abs(dataval)),color='0.5',linewidth=2)
            p.xlim(xmin,xmax)
            p.ylim(binsy_log.min(), binsy_log.max()); p.grid()
            p.subplot(223)
            p.plot(binsx_log,new_pCs[k],'k-'); p.grid()
            p.xlim(xmin,xmax)
            p.subplot(224)
            p.plot(binsx_log,new_pIs[k],'k-'); p.grid()
            p.xlim(xmin,xmax)
            p.tight_layout(); p.suptitle("k = " + str(k)); p.show()

    pC, pC_err = compute_stats(10**n.abs(binsx_log)*n.sign(binsx_log), new_pCs, pt_pCs)
    pI, pI_err = compute_stats(10**n.abs(binsx_log)*n.sign(binsx_log), new_pIs, pt_pIs)
    pC_fold, pC_fold_err = compute_stats(10**binsx_log, new_pCs_fold, pt_pCs_fold)
    pI_fold, pI_fold_err = compute_stats(10**binsx_log, new_pIs_fold, pt_pIs_fold)

    pC_old, pC_err_old = compute_stats(binsy_lin, old_pCs, pt_pCs, old=True)
    pI_old, pI_err_old = compute_stats(binsy_lin, old_pIs, pt_pIs, old=True)
    pC_fold_old, pC_fold_err_old = compute_stats(binsy_lin, old_pCs_fold, pt_pCs_fold, old=True)
    pI_fold_old, pI_fold_err_old = compute_stats(binsy_lin, old_pIs_fold, pt_pIs_fold, old=True)

    if count == 0: # data case
        pCv = pC; pCv_old = pC_old
        pCv_fold = pC_fold; pCv_fold_old = pC_fold_old
        pIv = pI; pIv_old = pI_old
        pIv_fold = pI_fold; pIv_fold_old = pI_fold_old
        pCv_err = pC_err; pCv_err_old = pC_err_old
        pCv_fold_err = pC_fold_err; pCv_fold_err_old = pC_fold_err_old
        pIv_err = pI_err; pIv_err_old = pI_err_old
        pIv_fold_err = pI_fold_err; pIv_fold_err_old = pI_fold_err_old
    if count == 1: # noise case
        pCn = pC; pCn_old = pC_old
        pCn_fold = pC_fold; pCn_fold_old = pC_fold_old
        pIn = pI; pIn_old = pI_old
        pIn_fold = pI_fold; pIn_fold_old = pI_fold_old
        pCn_err = pC_err; pCn_err_old = pC_err_old
        pCn_fold_err = pC_fold_err; pCn_fold_err_old = pC_fold_err_old
        pIn_err = pI_err; pIn_err_old = pI_err_old
        pIn_fold_err = pI_fold_err; pIn_fold_err_old = pI_fold_err_old

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
        #theory_noise = file['theory_noise'],
        #theory_noise_delta2 = file['theory_noise_delta2'],
        prob=0.975, kperp=file['kperp'], sep=opts.sep, kpl_fold=file['kpl_fold'],
        ngps=file['ngps'], nbls=file['nbls'], nbls_g=file['nbls_g'], nlsts_g=file['nlsts_g'],
        lsts=file['lsts'], afreqs=file['afreqs'], cnt_eff=file['cnt_eff'],
        frf_inttime=file['frf_inttime'], inttime=file['inttime'],
        cmd=file['cmd'].item() + ' \n '+' '.join(sys.argv))

outname2 = 'pspec_final_sep'+opts.sep+'_full.npz'
print '   Saving', outname2
# Save with 1-sigma errors (compatible with plot_pspec_final_v2.py)
n.savez(outname2, kpl=kpl, k=file['k'], freq=file['freq'],
        pCv=pCv, pCv_err=pCv_err,
        pCv_old=pCv_old, pCv_err_old=pCv_err_old,
        pCv_fold=pCv_fold, pCv_fold_err=pCv_fold_err,
        pCv_fold_old=pCv_fold_old, pCv_fold_err_old=pCv_fold_err_old,
        pIv=pIv, pIv_err=pIv_err,
        pIv_old=pIv_old, pIv_err_old=pIv_err_old,
        pIv_fold=pIv_fold, pIv_fold_err=pIv_fold_err,
        pIv_fold_old=pIv_fold_old, pIv_fold_err_old=pIv_fold_err_old,
        pCn=pCn, pCn_err=pCn_err,
        pCn_old=pCn_old, pCn_err_old=pCn_err_old,
        pCn_fold=pCn_fold, pCn_fold_err=pCn_fold_err,
        pCn_fold_old=pCn_fold_old, pCn_fold_err_old=pCn_fold_err_old,
        pIn=pIn, pIn_err=pIn_err,
        pIn_old=pIn_old, pIn_err_old=pIn_err_old,
        pIn_fold=pIn_fold, pIn_fold_err=pIn_fold_err,
        pIn_fold_old=pIn_fold_old, pIn_fold_err_old=pIn_fold_err_old,
        #theory_noise = file['theory_noise'],
        #theory_noise_delta2 = file['theory_noise_delta2'],
        prob=0.975, kperp=file['kperp'], sep=opts.sep, kpl_fold=file['kpl_fold'],
        ngps=file['ngps'], nbls=file['nbls'], nbls_g=file['nbls_g'], nlsts_g=file['nlsts_g'],
        lsts=file['lsts'], afreqs=file['afreqs'], cnt_eff=file['cnt_eff'],
        frf_inttime=file['frf_inttime'], inttime=file['inttime'],
        cmd=file['cmd'].item() + ' \n '+' '.join(sys.argv))
