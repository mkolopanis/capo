#! /usr/bin/env python
"""Compute Signal Loss factor via distribution method v9.
    -Reads in pspec_2d_to_1d.npz (contains PS points for aln bootstraps)
    -Assumes dense sampling in P_in
    -Plots all bootstraps per injection level at the same <P_in> value
    -For bins along the P_in axis, a 1D Gaussian is used to smooth out the points (fit in linear-space)
    -All gaussians are combined into a 2D transfer curve (P_in vs. P_r)
    -All gaussians are multiplied by a prior
    -Takes a horizontal cut at P_out = peak of p_x (data) and sums to 1
    -Uses this to compute new PS points + errors
"""
import numpy as n
import pylab as p
from matplotlib import gridspec
from matplotlib.colors import LogNorm
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
    nbins = 101 # XXX hard-coded number of bins for P_out axis
    dmax = max(n.max(Pins_fold.values()),n.max(pIs.values())) # maximum is determined from whichever is largest: EoR injection or "I" value
    dg = n.min(Pins_fold.values()) #1 # XXX artificially set based on starting P_in values
    grid = n.logspace(n.log10(dg),n.log10(dmax),nbins) # logspace sampling
    binsy_log = n.log10(grid) # log-space
    binsy_lin = n.concatenate((-10**binsy_log[::-1],10**binsy_log)) # real numbers (not log-space)
    binsy_log = n.concatenate((-binsy_log[::-1],binsy_log)) # for pos and neg
    return binsy_log,binsy_lin

# Fit gaussian
def make_gauss(bins,mean,sigma):
    A = 1/(sigma*n.sqrt(2*n.pi))
    return A*n.exp((-1/2.)*((bins-mean)/sigma)**2) * bin_size(bins)

# Compute Jeffrey prior
def compute_jeffrey(xs,ys):
    x,y_mean,y_std = [],[],[]
    for i in range(xs.shape[0]): # loop over injections
        x.append(xs[i][0])
        y_mean.append(n.mean(ys[i]))
        y_std.append(n.std(ys[i]))
    x = n.array(x)
    y_mean = n.array(y_mean)
    y_std = n.array(y_std)
    fit_mean = n.poly1d(n.polyfit(x,y_mean,2)) # fit polynomial
    fit_std = n.poly1d(n.polyfit(x,y_std,2))
    fit_mean_d = n.polyder(fit_mean, m=1) # derivative
    fit_std_d = n.polyder(fit_std, m=1)
    dsigma_dx = fit_std_d(x) # partial derivatives
    dy_dx = fit_mean_d(x)
    jeffrey = n.sqrt((1/y_std) * (2*(dsigma_dx)**2 + (dy_dx)**2))
    return jeffrey

# For all data points (full distributions), smooth the transfer curves for both the C and I cases using 1D KDE's per injection level
def smooth_dist(fold=True):
    T_C = {}
    T_I = {}
    binsx = {} # will hold the P_in injection bins per k
    if fold == True: ks = kpl_fold # k values to loop over
    if fold == False: ks = kpl
    for kk,k in enumerate(ks):

        if kk == 8 and count == 0 and fold == True: # save values out
            # Save values to use for plotting sigloss terms
            if opts.nosubtract: fn = 'pspec_sigloss_terms_nosubtract.npz'
            else: fn = 'pspec_sigloss_terms.npz'
            print "Saving",fn,"which contains data values for k =",k
            n.savez(fn, k=k, pCv=file['pCv_fold'][n.where(kpl_fold==k)[0][0]], pIv=file['pIv_fold'][n.where(kpl_fold==k)[0][0]], Pins=Pins_fold[k], Pouts=Pouts_fold[k], Pouts_I=Pouts_I_fold[k], pCrs_fold=pCrs_fold[k], pCvs_Cr_fold=pCvs_Cr_fold[k], pCes_Cr_fold=pCes_Cr_fold[k], pCves_fold=pCves_fold[k], pIves_fold=pIves_fold[k])
        
        binsx_log = [] # will hold the P_in injection bins for this k
        if fold == True:
            xs = n.array(Pins_fold[k])
            ys_C = n.array(Pouts_fold[k])
            ys_I = n.array(Pouts_I_fold[k])
        if fold == False:
            xs = n.array(Pins[k])
            ys_C = n.array(Pouts[k])
            ys_I = n.array(Pouts_I[k])
        # Find average P_in per injection level and build up x-axis using that
        for ll,level in enumerate(xs): 
            Pin_mid = n.mean(xs[ll])
            binsx_log.append(n.log10(Pin_mid)) # 1D array of x-values (log)
            xs[ll] = n.repeat(Pin_mid,len(level)) # same shape as ys
        binsx_log = n.array(binsx_log) # make an array
        binsx[k] = binsx_log
        # Make empty 2D transfer curves that will hold distributions
        TC = n.zeros((len(binsy_log),len(binsx_log))) 
        TI = n.zeros((len(binsy_log),len(binsx_log)))
        # Compute Jeffrey prior
        #jeffrey_prior = compute_jeffrey(xs,ys_C)
        #jeffrey_prior_I = compute_jeffrey(xs,ys_I)
        # loop over injections
        for sub_bin in range(len(binsx_log)): 
            xs_sub_bin = xs[sub_bin] # P_in values for the injection bin
            ys_sub_bin = ys_C[sub_bin] # P_out values for the injection bin
            ys_I_sub_bin = ys_I[sub_bin]
            # fit (normalized) gaussian distributions
            data_dist_C = make_gauss(binsy_lin,n.mean(ys_sub_bin),n.std(ys_sub_bin))
            data_dist_I = make_gauss(binsy_lin,n.mean(ys_I_sub_bin),n.std(ys_I_sub_bin))
            # Multiply the P_in column by the prior
            TC[:,sub_bin] = data_dist_C#*prior_factor(10**binsx_log)[sub_bin]#*jeffrey_prior[sub_bin]
            TI[:,sub_bin] = data_dist_I#*prior_factor(10**binsx_log)[sub_bin]#*jeffrey_prior_I[sub_bin]

        # Save values to use for plotting sigloss plots
        if count == 0 and fold == True and kk == 8:
            if opts.nosubtract: fn = 'pspec_sigloss_nosubtract.npz'
            else: fn = 'pspec_sigloss.npz'
            print "Saving",fn,"which contains data values for k =",k
            n.savez(fn, k=k, binsx=binsx_log, binsy=binsy_log, kdeC=TC, kdeI=TI, xs=xs, ys=ys_C, ysI=ys_I, pC=file['pCv_fold'][n.where(ks==k)[0][0]], pI=file['pIv_fold'][n.where(ks==k)[0][0]])

        # Plot KDE and points
        if opts.plot:
            p.figure(figsize=(10,6))
            p.subplot(121)
            p.pcolormesh(binsx_log,binsy_log,TC,cmap='hot_r',vmax=0.05,vmin=0)
            if fold == True and count == 0: dataval = file['pCv_fold'][n.where(ks==k)[0][0]]
            if fold == False and count == 0: dataval = file['pCv'][n.where(ks==k)[0][0]]
            if fold == True and count == 1: dataval = file['pCn_fold'][n.where(ks==k)[0][0]]
            if  fold == False and count == 1: dataval = file['pCn'][n.where(ks==k)[0][0]]
            p.axhline(y=n.sign(dataval)*n.log10(n.abs(dataval)),color='0.5',linewidth=2)
            p.plot(n.sign(xs)*n.log10(n.abs(xs)), n.sign(ys_C)*n.log10(n.abs(ys_C)), 'k.')
            p.xlim(binsx_log.min(), binsx_log.max())
            p.ylim(binsy_log.min(), binsy_log.max()); p.grid()
            p.subplot(122)
            p.pcolormesh(binsx_log,binsy_log,TI,cmap='hot_r',vmax=0.05,vmin=0)
            if fold == True and count == 0: dataval = file['pIv_fold'][n.where(ks==k)[0][0]]
            if fold == False and count == 0: dataval = file['pIv'][n.where(ks==k)[0][0]]
            if fold == True and count == 1: dataval = file['pIn_fold'][n.where(ks==k)[0][0]]
            if fold == False and count == 1: dataval = file['pIn'][n.where(ks==k)[0][0]]
            p.axhline(y=n.sign(dataval)*n.log10(n.abs(dataval)),color='0.5',linewidth=2)
            p.plot(n.sign(xs)*n.log10(n.abs(xs)), n.sign(ys_I)*n.log10(n.abs(ys_I)), 'k.')
            p.xlim(binsx_log.min(), binsx_log.max())
            p.ylim(binsy_log.min(), binsy_log.max())
            p.tight_layout(); p.suptitle("k = " + str(k)); p.grid(); p.show()

        T_C[k] = TC # populate dictionary
        T_I[k] = TI
    return T_C, T_I, binsx # return KDE's and P_in bins

# Apply transfer function to data
def sigloss_func(pt, M_matrix, bins):
    new_PS = {} # dictionaries that will hold PS distributions
    for k in pt.keys():
        point = pt[k] # value of data point
        M = M_matrix[k].copy()
        peak_ind = n.argmin(n.abs(binsy_lin-point)) # matching index
        cut = M[peak_ind] # cut of KDE at peak of data
        cut /= n.sum(cut) # normalization
        new_PS[k] = cut
    return new_PS # distribution of bins_concat

# Compute PS points and errors based on std(bootstraps)
def boot_stats(pt, boots):
    pts, errs = [],[]
    for key in n.sort(pt.keys()):
        pts.append(pt[key])
        errs.append(n.std(boots[key][0])) # 1-sigma
    return pts, errs

# Compute PS points and errors
def compute_stats(binsx, data, pt):
    pts, errs = [], []
    for key in n.sort(data.keys()):
        #point = bins[n.argmax(data[key])] # XXX peak of dist
        bins = n.sign(binsx[key])*10**(n.abs(binsx[key]))
        point = n.sum(bins*data[key])/n.sum(data[key]) # weighted avg
        pts.append(point)
        percents = [n.sum((data[key])[:i]) for i in range(len(data[key]))]
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

    # Save PS points and errors from bootstrapping
    pC_old, pC_err_old = boot_stats(pt_pCs, pCs)
    pI_old, pI_err_old = boot_stats(pt_pIs, pIs)
    pC_fold_old, pC_fold_err_old = boot_stats(pt_pCs_fold, pCs_fold)
    pI_fold_old, pI_fold_err_old = boot_stats(pt_pIs_fold, pIs_fold)
    
    # Compute signl loss and new limits
    if opts.skip_sigloss: # over-ride new points and errors with old ones
        print "Skipping Signal Loss!"
        pC, pC_err = pC_old, pC_err_old
        pI, pI_err = pI_old, pI_err_old
        pC_fold, pC_fold_err = pC_fold_old, pC_fold_err_old
        pI_fold, pI_fold_err = pI_fold_old, pI_fold_err_old
    else:
        T_C, T_I, binsx_log = smooth_dist(fold=False) # transfer curves
        T_C_fold, T_I_fold, binsx_log_fold = smooth_dist(fold=True)
        # new distributions and stats
        new_pCs = sigloss_func(pt_pCs, T_C, binsx_log)
        new_pCs_fold = sigloss_func(pt_pCs_fold, T_C_fold, binsx_log_fold)
        new_pIs = sigloss_func(pt_pIs, T_I, binsx_log)
        new_pIs_fold = sigloss_func(pt_pIs_fold, T_I_fold, binsx_log_fold)
        pC, pC_err = compute_stats(binsx_log, new_pCs, pt_pCs)
        pI, pI_err = compute_stats(binsx_log, new_pIs, pt_pIs)
        pC_fold, pC_fold_err = compute_stats(binsx_log_fold, new_pCs_fold, pt_pCs_fold)
        pI_fold, pI_fold_err = compute_stats(binsx_log_fold, new_pIs_fold, pt_pIs_fold)
    
    # Plot un-folded case
    if opts.plot:
        for k in T_C:
            p.figure(figsize=(10,10))
            xmin,xmax = 0,16
            p.subplot(221)
            p.pcolormesh(binsx_log[k],binsy_log,T_C[k],cmap='hot_r',vmin=0,vmax=0.05)
            if count == 0: dataval = file['pCv'][n.where(kpl==k)[0][0]]
            if count == 1: dataval = file['pCn'][n.where(kpl==k)[0][0]]
            p.axhline(y=n.sign(dataval)*n.log10(n.abs(dataval)),color='0.5',linewidth=2)
            p.xlim(xmin,xmax)
            p.ylim(binsy_log.min(), binsy_log.max()); p.grid()
            p.subplot(222)
            p.pcolormesh(binsx_log[k],binsy_log,T_I[k],cmap='hot_r',vmin=0,vmax=0.5)
            if count == 0: dataval = file['pIv'][n.where(kpl==k)[0][0]]
            if count == 1: dataval = file['pIn'][n.where(kpl==k)[0][0]]
            p.axhline(y=n.sign(dataval)*n.log10(n.abs(dataval)),color='0.5',linewidth=2)
            p.xlim(xmin,xmax)
            p.ylim(binsy_log.min(), binsy_log.max()); p.grid()
            p.subplot(223)
            p.plot(binsx_log[k],new_pCs[k],'k-'); p.grid()
            p.xlim(xmin,xmax)
            p.subplot(224)
            p.plot(binsx_log[k],new_pIs[k],'k-'); p.grid()
            p.xlim(xmin,xmax)
            p.tight_layout(); p.suptitle("k = " + str(k)); p.show()

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


# create dictionary of extra keys for the save calls
extra_out_dict = {
        'prob':0.975,
        'kperp':file['kperp'],
        'sep':opts.sep,
        'kpl_fold':file['kpl_fold'],
        'ngps':file['ngps'],
        'nbls':file['nbls'],
        'nbls_g':file['nbls_g'],
        'lsts':file['lsts'],
        'afreqs':file['afreqs'],
        'cnt_eff':file['cnt_eff'],
        'frf_inttime':file['frf_inttime'],
        'inttime':file['inttime'],
        'cmd':file['cmd'].item() + ' \n '+' '.join(sys.argv)
        }

if ('theory_noise' and 'theory_noise_delta2') in file.keys():
        extra_out_dict['theory_noise']=file['theory_noise']
        extra_out_dict['theory_noise_delta2']=file['theory_noise_delta2']
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
        **extra_out_dict)

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
        **extra_out_dict)
