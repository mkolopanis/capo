#! /usr/bin/env python
"""Compute Signal Loss factor via distribution method."""
import numpy as n
import pylab as p
from matplotlib import gridspec
from scipy.interpolate import interp1d
from capo.pspec import f2z
import glob
import optparse
import sys

o = optparse.OptionParser()
o.set_description(__doc__)
o.add_option('--plot',action='store_true',
            help='Show plots.')
o.add_option('--median',action='store_true',
            help='Correct by median correction factor.')
opts,args = o.parse_args(sys.argv[1:])

# read one pspec_pk_k3pk.npz file for data points
file = n.load(glob.glob('inject_*')[0]+'/pspec_pk_k3pk.npz')
pCv = n.abs(file['pCv']); pCv_fold = n.abs(file['pCv_fold'])
pCv_err = file['pCv_err']; pCv_fold_err = n.abs(file['pCv_fold_err'])
pIv = n.abs(file['pIv']); pIv_fold = n.abs(file['pIv_fold'])
pIv_err = file['pIv_err']; pIv_fold_err = n.abs(file['pIv_fold_err'])
pCn = n.abs(file['pCn']); pCn_fold = n.abs(file['pCn_fold'])
pCn_err = file['pCn_err']; pCn_fold_err = n.abs(file['pCn_fold_err'])
pIn = n.abs(file['pIn']); pIn_fold = n.abs(file['pIn_fold'])
pIn_err = file['pIn_err']; pIn_fold_err = n.abs(file['pIn_fold_err'])
kpl = file['kpl']; k = file['k']
#  absolute values only used for signal loss estimation
#  both positive/negative values are saved in npz file later

# loop over injects and read pspec_2d_to_1d outputs
for count in range(2):
    if count == 1: print 'NOISE CASE'
    else: print 'DATA CASE'   
    if opts.plot: 
        fig1 = p.figure(1, figsize=(15, 7))
        fig2 = p.figure(2, figsize=(15, 7))
        fig3 = p.figure(3, figsize=(15, 7))
        fig4 = p.figure(4, figsize=(15, 7))
    pklo, pkhi = 1e1, 1e12
    Pouts = {}; Pouts_fold = {}
    Pouts_I = {}; Pouts_I_fold = {}
    pCs = {}; pCs_fold = {}
    pIs = {}; pIs_fold = {}
    alphas = {}; alphas_fold = {}
    alphas_I = {}; alphas_I_fold = {}
    for inject in glob.glob('inject_*'):
        print 'Reading', inject
        file_2d = n.load(inject + '/pspec_2d_to_1d.npz')
        if count == 1: # noise case
            Pout = n.abs(file_2d['pCs-pCn']); Pout_fold = n.abs(file_2d['pCs-pCn_fold'])
            Pout_I = n.abs(file_2d['pIs-pIn']); Pout_I_fold = n.abs(file_2d['pIs-pIn_fold'])
            pC = n.abs(file_2d['pCn']); pC_fold = n.abs(file_2d['pCn_fold'])
            pI = n.abs(file_2d['pIn']); pI_fold = n.abs(file_2d['pIn_fold'])
        else:
            Pout = n.abs(file_2d['pCr-pCv']); Pout_fold = n.abs(file_2d['pCr-pCv_fold']) # shape (#boots, #kpl)
            Pout_I = n.abs(file_2d['pIr-pIv']); Pout_I_fold = n.abs(file_2d['pIr-pIv_fold'])  # pI case
            pC = n.abs(file_2d['pCv']); pC_fold = n.abs(file_2d['pCv_fold'])
            pI = n.abs(file_2d['pIv']); pI_fold = n.abs(file_2d['pIv_fold'])
        Pin = n.abs(file_2d['pIe'])
        Pin_fold = n.abs(file_2d['pIe_fold'])

        for ind in range(len(k)): # loop through k for Delta^2(k)
            try:
                Pouts_fold[k[ind]].append(Pout_fold[:,ind])
                Pouts_I_fold[k[ind]].append(Pout_I_fold[:,ind])
                pCs_fold[k[ind]] = [pC_fold[:,ind]] # no appending because it's the same every injection
                pIs_fold[k[ind]] = [pI_fold[:,ind]]
                alphas_fold[k[ind]].append(n.abs(Pin_fold[:,ind]/Pout_fold[:,ind]))
                alphas_I_fold[k[ind]].append(n.abs(Pin_fold[:,ind]/Pout_I_fold[:,ind]))
            except:
                Pouts_fold[k[ind]] = [Pout_fold[:,ind]]
                Pouts_I_fold[k[ind]] = [Pout_I_fold[:,ind]]
                pCs_fold[k[ind]] = [pC_fold[:,ind]]
                pIs_fold[k[ind]] = [pI_fold[:,ind]]
                alphas_fold[k[ind]] = [n.abs(Pin_fold[:,ind]/Pout_fold[:,ind])]
                alphas_I_fold[k[ind]] = [n.abs(Pin_fold[:,ind]/Pout_I_fold[:,ind])]
        
        for ind in range(len(kpl)):  # loop through kpl  for P(k)
            try:
                Pouts[kpl[ind]].append(Pout[:,ind])
                Pouts_I[kpl[ind]].append(Pout_I[:,ind])
                pCs[kpl[ind]] = [pC[:,ind]] # no appending because it's the same every injection
                pIs[kpl[ind]] = [pI[:,ind]]
                alphas[kpl[ind]].append(n.abs(Pin[:,ind]/Pout[:,ind]))
                alphas_I[kpl[ind]].append(n.abs(Pin[:,ind]/Pout_I[:,ind]))
            except:
                Pouts[kpl[ind]] = [Pout[:,ind]]
                Pouts_I[kpl[ind]] = [Pout_I[:,ind]]
                pCs[kpl[ind]] = [pC[:,ind]]
                pIs[kpl[ind]] = [pI[:,ind]]
                alphas[kpl[ind]] = [n.abs(Pin[:,ind]/Pout[:,ind])]
                alphas_I[kpl[ind]] = [n.abs(Pin[:,ind]/Pout_I[:,ind])]
        
            if opts.plot:
        
                p.figure(1) # Pin vs. Pout
                p.subplot(3, 7, ind)
                p.loglog(Pin[:,ind], Pout[:,ind], 'k.')  # points
                # p.loglog(Pin[ind],Pout_noise[ind],'b.') # noise points
                p.loglog([pklo, pkhi], [pklo, pkhi], 'k-')  # diagonal line
                p.grid(True)
                p.xlim(pklo, pkhi)
                p.ylim(pklo, pkhi)
                p.tick_params(axis='both', which='both', labelsize=6)
                p.title('kpl = '+str("%.4f" % kpl[ind]), fontsize=8)
                 
                p.figure(2) # alpha vs. Pout
                p.subplot(3, 7, ind)
                p.loglog(Pin[:,ind]/Pout[:,ind],Pout[:,ind],'k.') # points
                p.grid(True)
                p.xlim(1e-3,1e7)
                p.ylim(pklo, pkhi)
                p.tick_params(axis='both', which='both', labelsize=6)
                p.title('kpl = '+str("%.4f" % kpl[ind]), fontsize=8)
                
                p.figure(3) # Pin vs. Pout for I case
                p.subplot(3, 7, ind)
                p.loglog(Pin[:,ind], Pout_I[:,ind], 'k.')  # points
                p.loglog([pklo, pkhi], [pklo, pkhi], 'k-')  # diagonal line
                p.grid(True)
                p.xlim(pklo, pkhi)
                p.ylim(pklo, pkhi)
                p.tick_params(axis='both', which='both', labelsize=6)
                p.title('kpl = '+str("%.4f" % kpl[ind]), fontsize=8)

                p.figure(4) # alpha vs. Pout for I case
                p.subplot(3, 7, ind)
                p.loglog(Pin[:,ind]/Pout_I[:,ind],Pout_I[:,ind],'k.') # points
                p.grid(True)
                p.xlim(1e-3,1e7)
                p.ylim(pklo, pkhi)
                p.tick_params(axis='both', which='both', labelsize=6)
                p.title('kpl = '+str("%.4f" % kpl[ind]), fontsize=8)

    if opts.plot:

        # plot Pin vs. Pout
        p.figure(1)
        p.tight_layout()
        fig1.text(0.5, 0.0, r'$P_{\rm in}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',
                 ha='center')
        fig1.text(0.0, 0.5, r'$P_{\rm out}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',
                 va='center', rotation='vertical')

        # plot signal loss factor (alpha) vs. Pout
        p.figure(2)
        p.tight_layout()
        fig2.text(0.0, 0.5, r'$P_{\rm out}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',
                 va='center', rotation='vertical')
        fig2.text(0.5, 0.0, r'$Signal Loss Factors$', ha='center')

        # plot Pin vs. Pout for I case
        p.figure(3)
        p.tight_layout()
        fig3.text(0.5, 0.0, r'$P_{\rm in}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',
                 ha='center')
        fig3.text(0.0, 0.5, r'$P_{\rm out,I}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',
                 va='center', rotation='vertical')

        # plot signal loss factor (alpha) vs. Pout
        p.figure(4)
        p.tight_layout()
        fig4.text(0.0, 0.5, r'$P_{\rm out,I}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',
                 va='center', rotation='vertical')
        fig4.text(0.5, 0.0, r'$Signal Loss Factors$', ha='center')

    if opts.plot: p.show()
    
    for key in alphas:
        alphas[key] = n.array(alphas[key]).flatten()
        alphas_I[key] = n.array(alphas_I[key]).flatten()
        Pouts[key] = n.array(Pouts[key]).flatten()
        Pouts_I[key] = n.array(Pouts_I[key]).flatten()
        pCs[key] = n.array(pCs[key]).flatten()
        pIs[key] = n.array(pIs[key]).flatten()

    for key in alphas_fold:
        alphas_fold[key] = n.array(alphas_fold[key]).flatten()
        alphas_I_fold[key] = n.array(alphas_I_fold[key]).flatten()
        Pouts_fold[key] = n.array(Pouts_fold[key]).flatten()
        Pouts_I_fold[key] = n.array(Pouts_I_fold[key]).flatten()
        pCs_fold[key] = n.array(pCs_fold[key]).flatten()
        pIs_fold[key] = n.array(pIs_fold[key]).flatten()

    # bin Pout, Pout_I, pCv, pIv
    #nbins = Pouts[key].shape[0] / Pout.shape[0] # number of injections
    nbins = 100
    nbins_I = 100
    sigloss_final = {}; sigloss_final_fold = {}
    sigloss_final_I = {}; sigloss_final_I_fold = {}
    
    for ind in range(len(kpl)): # loop for P(k)
        low = n.log10(Pouts[kpl[ind]].min()) # XXX haven't taken into account alphas < 1, which are artificial
        high = n.log10(Pouts[kpl[ind]].max())
        low_I = n.log10(Pouts_I[kpl[ind]].min()) # lowest Pout
        high_I = n.log10(Pouts_I[kpl[ind]].max()) # highest Pout
        bins = n.logspace(low,high,nbins) # bins for Pout
        bins_I = n.logspace(low_I,high_I,nbins_I)
        bin_ind = n.digitize(Pouts[kpl[ind]],bins) # bin indices for all Pout values (ranges from 0 to nbins+1)
        bin_ind_I = n.digitize(Pouts_I[kpl[ind]],bins_I)
        bin_ind_pC = n.digitize(pCs[kpl[ind]],bins) # bin indices for all pC values
        bin_ind_pI = n.digitize(pIs[kpl[ind]],bins_I)
        final_alphas = []
        final_alphas_I = []
        p_sum = []
        p_sum_I = []
        for i in range(len(bins)): # match alpha bins with prob bins
            factors = alphas[kpl[ind]][n.where(bin_ind == i)]
            factors_I = alphas_I[kpl[ind]][n.where(bin_ind_I == i)]
            probs = len(pCs[kpl[ind]][n.where(bin_ind_pC == i)])/float(len(pCs[kpl[ind]]))
            probs_I = len(pIs[kpl[ind]][n.where(bin_ind_pI == i)])/float(len(pIs[kpl[ind]]))
            if len(factors) > 0: p_sum.append(probs)
            if len(factors_I) > 0: p_sum_I.append(probs_I)
            final_alphas = n.concatenate((final_alphas,factors*probs/len(factors))) 
            final_alphas_I = n.concatenate((final_alphas_I,factors_I*probs_I/len(factors_I)))
        sigloss_final[kpl[ind]] = n.sum(final_alphas)/n.sum(p_sum)
        sigloss_final_I[kpl[ind]] = n.sum(final_alphas_I)/n.sum(p_sum_I)
        if True: #opts.plot:
            p.plot(kpl[ind],sigloss_final[kpl[ind]],'k.',label='pC' if ind == 0 else "")
            p.plot(kpl[ind],sigloss_final_I[kpl[ind]],'b.',label='pI' if ind == 0 else "")
    if True: p.xlabel('k');p.ylabel('Signal Loss Factors');p.legend();p.show()
 
    for ind in range(len(k)): # loop for Delta^2(k)
        low = n.log10(Pouts_fold[k[ind]].min()) # XXX haven't taken into account alphas < 1, which are artificial
        high = n.log10(Pouts_fold[k[ind]].max())
        low_I = n.log10(Pouts_I_fold[k[ind]].min()) # lowest Pout
        high_I = n.log10(Pouts_I_fold[k[ind]].max()) # highest Pout
        bins = n.logspace(low,high,nbins) # bins for Pout
        bins_I = n.logspace(low_I,high_I,nbins_I)
        bin_ind = n.digitize(Pouts_fold[k[ind]],bins) # bin indices for all Pout values (ranges from 0 to nbins+1)
        bin_ind_I = n.digitize(Pouts_I_fold[k[ind]],bins_I)
        bin_ind_pC = n.digitize(pCs_fold[k[ind]],bins) # bin indices for all pC values
        bin_ind_pI = n.digitize(pIs_fold[k[ind]],bins_I)
        final_alphas = []
        final_alphas_I = []
        p_sum = []
        p_sum_I = []
        for i in range(len(bins)): # match alpha bins with prob bins
            factors = alphas_fold[k[ind]][n.where(bin_ind == i)]
            factors_I = alphas_I_fold[k[ind]][n.where(bin_ind_I == i)]
            probs = len(pCs_fold[k[ind]][n.where(bin_ind_pC == i)])/float(len(pCs_fold[k[ind]]))
            probs_I = len(pIs_fold[k[ind]][n.where(bin_ind_pI == i)])/float(len(pIs_fold[k[ind]]))
            if len(factors) > 0: p_sum.append(probs)
            if len(factors_I) > 0: p_sum_I.append(probs_I)
            final_alphas = n.concatenate((final_alphas,factors*probs/len(factors))) 
            final_alphas_I = n.concatenate((final_alphas_I,factors_I*probs_I/len(factors_I)))
        sigloss_final_fold[k[ind]] = n.sum(final_alphas)/n.sum(p_sum)
        sigloss_final_I_fold[k[ind]] = n.sum(final_alphas_I)/n.sum(p_sum_I)
    
    if count == 1:
        sigfactors_noise = sigloss_final.values(); sigfactors_noise_fold = sigloss_final_fold.values()
        sigfactors_noise_I = sigloss_final_I.values(); sigfactors_noise_I_fold = sigloss_final_I_fold.values()
    else:
        sigfactors = sigloss_final.values(); sigfactors_fold = sigloss_final_fold.values()
        sigfactors_I = sigloss_final_I.values(); sigfactors_I_fold = sigloss_final_I_fold.values()

# save final values

if opts.median: other_factors = 1/n.log(2)  # median correction factor
else: other_factors = 1
fold_factor = file['k']**3/(2*n.pi**2)

pIv = file['pIv']*sigfactors_I*other_factors
pCv = file['pCv']*sigfactors*other_factors
pIn = file['pIn']*sigfactors_noise_I*other_factors
pCn = file['pCn']*sigfactors_noise*other_factors
pIv_err = file['pIv_err']*sigfactors_I*other_factors
pCv_err = file['pCv_err']*sigfactors*other_factors
pIn_err = file['pIn_err']*sigfactors_noise_I*other_factors
pCn_err = file['pCn_err']*sigfactors_noise*other_factors
pIv_fold = file['pIv_fold']*sigfactors_I_fold*fold_factor*other_factors
pCv_fold = file['pCv_fold']*sigfactors_fold*fold_factor*other_factors
pIn_fold = file['pIn_fold']*sigfactors_noise_I_fold*fold_factor*other_factors
pCn_fold = file['pCn_fold']*sigfactors_noise_fold*fold_factor*other_factors
pIv_fold_err = file['pIv_fold_err']*sigfactors_I_fold*fold_factor*other_factors
pCv_fold_err = file['pCv_fold_err']*sigfactors_fold*fold_factor*other_factors
pIn_fold_err = file['pIn_fold_err']*sigfactors_noise_I_fold*fold_factor*other_factors
pCn_fold_err = file['pCn_fold_err']*sigfactors_noise_fold*fold_factor*other_factors

print '   Saving pspec_final.npz'  # XXX 2-sigma probability is hard-coded
n.savez('pspec_final.npz', kpl=kpl, k=file['k'], freq=file['freq'],
        pC=pCv, pC_up=2 * pCv_err,
        pC_fold=pCv_fold, pC_fold_up=2 * pCv_fold_err,
        pI=pIv, pI_up=2 * pIv_err,
        pI_fold=pIv_fold, pI_fold_up=2 * pIv_fold_err,
        pCn=pCn, pCn_up=2 * pCn_err,
        pCn_fold=pCn_fold, pCn_fold_up=2 * pCn_fold_err,
        pIn=pIn, pIn_up=2 * pIn_err,
        pIn_fold=pIn_fold, pIn_fold_up=2 * pIn_fold_err,
        prob=0.9545, 
        alphaCv=sigfactors, alphaIv=sigfactors_I,
        alphaCn=sigfactors_noise, alphaIn=sigfactors_noise_I,
        alphaCv_fold=sigfactors_fold, alphaIv_fold=sigfactors_I_fold,
        alphaCn_fold=sigfactors_noise_fold, alphaIn_fold=sigfactors_noise_I_fold,
        cmd=' '.join(sys.argv))
