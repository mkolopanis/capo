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
opts,args = o.parse_args(sys.argv[1:])

# read one pspec_pk_k3pk.npz file for data points
file = n.load(glob.glob('inject_*')[0]+'/pspec_pk_k3pk.npz')
pCv = n.abs(file['pCv'])
pCv_err = file['pCv_err']
pIv = n.abs(file['pIv'])
pIv_err = file['pIv_err']
pCn = n.abs(file['pCn'])
pCn_err = file['pCn_err']
pIn = n.abs(file['pIn'])
pIn_err = file['pIn_err']
kpl = file['kpl']
#  absolute values only used for signal loss estimation
#  both positive/negative values are saved in npz file later

# loop over injects and read pspec_2d_to_1d outputs
for count in range(2):
    if count == 1: print 'NOISE CASE'
    else: print 'DATA CASE'    
    fig1 = p.figure(1, figsize=(15, 7))
    fig2 = p.figure(2, figsize=(15, 7))
    fig3 = p.figure(3, figsize=(15, 7))
    fig4 = p.figure(4, figsize=(15, 7))
    pklo, pkhi = 1e1, 1e12
    Pouts = {}
    Pouts_I = {}
    pCs = {}
    pIs = {}
    alphas = {}
    alphas_I = {}
    for inject in glob.glob('inject_*'):
        print 'Reading', inject
        file_2d = n.load(inject + '/pspec_2d_to_1d.npz')
        if count == 1: # noise case
            Pout = n.abs(file_2d['pCs-pCn'])
            Pout_I = n.abs(file_2d['pIs-pIn'])
            pC = n.abs(file_2d['pCn'])
            pI = n.abs(file_2d['pIn'])
        else:
            Pout = n.abs(file_2d['pCr-pCv']) # shape (#boots, #kpl)
            Pout_I = n.abs(file_2d['pIr-pIv'])  # pI case
            pC = n.abs(file_2d['pCv'])
            pI = n.abs(file_2d['pIv'])
        Pin = n.abs(file_2d['pIe'])
        for ind in range(len(kpl)):  # plot for each k
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

    # bin Pout, Pout_I, pCv, pIv
    nbins = 100
    sigloss_final = {}
    sigloss_final_I = {}
    for ind in range(len(kpl)):
        low = n.log10(Pouts[kpl[ind]].min()) # XXX haven't taken into account alphas < 1, which are artificial
        high = n.log10(Pouts[kpl[ind]].max())
        low_I = n.log10(Pouts_I[kpl[ind]].min()) # lowest Pout
        high_I = n.log10(Pouts_I[kpl[ind]].max()) # highest Pout
        bins = n.logspace(low,high,nbins) # bins for Pout
        bins_I = n.logspace(low_I,high_I,nbins)
        bin_ind = n.digitize(Pouts[kpl[ind]],bins) # bin indices for all Pout values (ranges from 0 to nbins+1)
        bin_ind_I = n.digitize(Pouts_I[kpl[ind]],bins_I)
        bin_ind_pC = n.digitize(pCs[kpl[ind]],bins) # bin indices for all pC values
        bin_ind_pI = n.digitize(pIs[kpl[ind]],bins_I)
        """
        xs = n.logspace(low,high,1000) # x values along Pout range
        xs_I = n.logspace(low_I,high_I,1000)
        if count == 1:  
            mu = pCn[ind]
            mu_I = pIn[ind]
            sigma = pCn_err[ind]
            sigma_I = pIn_err[ind]
        else:
            mu = pCv[ind] # data median
            mu_I = pIv[ind] 
            sigma = pCv_err[ind] # data std
            sigma_I = pIv_err[ind] 
        rayleigh = (xs/sigma**2)*n.exp(-(xs**2)/(2*sigma**2)) # bellcurve
        rayleigh_I = (xs_I/sigma_I**2)*n.exp(-(xs_I**2)/(2*sigma_I**2))
        rayleigh = rayleigh / n.sum(rayleigh) # normalize
        rayleigh_I = rayleigh_I / n.sum(rayleigh_I)
        #p.plot(xs,rayleigh); p.xscale('log'); p.show()
        #gaussian = n.exp(-(1/2.)*((xs-mu)**2)/sigma**2) # bellcurve
        #gaussian = gaussian / n.sum(gaussian) # normalize
        probs = [] # probs[i] holds the probability for the ith bin
        probs_I = []
        factors = [] # factors[i] holds all alphas for the ith bin
        factors_I = []
        alphas[kpl[ind]][n.where(alphas[kpl[ind]] < 1)[0]] = 1.0 # overwrite factors < 1 with 1.0
        alphas_I[kpl[ind]][n.where(alphas_I[kpl[ind]] < 1)[0]] = 1.0
        for i,bin in enumerate(bins):
            factors.append(alphas[kpl[ind]][n.where(bin_ind == i)]) 
            factors_I.append(alphas_I[kpl[ind]][n.where(bin_ind_I == i)])
            if i == 0: 
                probs.append(n.sum(rayleigh[:n.argmin(n.abs(bins[i]-xs))]))
                probs_I.append(n.sum(rayleigh_I[:n.argmin(n.abs(bins_I[i]-xs_I))]))
            else: 
                probs.append(n.sum(rayleigh[n.argmin(n.abs(bins[i-1]-xs)):n.argmin(n.abs(bins[i]-xs))]))
                probs_I.append(n.sum(rayleigh_I[n.argmin(n.abs(bins_I[i-1]-xs_I)):n.argmin(n.abs(bins_I[i]-xs_I))]))
        """
        final_alphas = []
        final_alphas_I = []
        for i in range(len(bins)): # match alpha bins with prob bins
            factors = alphas[kpl[ind]][n.where(bin_ind == i)]
            factors_I = alphas_I[kpl[ind]][n.where(bin_ind_I == i)]
            probs = len(pCs[kpl[ind]][n.where(bin_ind_pC == i)])/float(len(pCs[kpl[ind]]))
            probs_I = len(pIs[kpl[ind]][n.where(bin_ind_pI == i)])/float(len(pIs[kpl[ind]]))
            final_alphas = n.concatenate((final_alphas,factors*probs/len(factors))) 
            final_alphas_I = n.concatenate((final_alphas_I,factors_I*probs_I/len(factors_I)))
        sigloss_final[kpl[ind]] = n.sum(final_alphas)
        sigloss_final_I[kpl[ind]] = n.sum(final_alphas_I)
        p.plot(kpl[ind],sigloss_final[kpl[ind]],'k.',label='pC' if ind == 0 else "")
        p.plot(kpl[ind],sigloss_final_I[kpl[ind]],'b.',label='pI' if ind == 0 else "")
    if opts.plot: p.xlabel('k');p.ylabel('Signal Loss Factors');p.legend();p.show()
    if count == 1:
        sigfactors_noise = sigloss_final.values()
        sigfactors_noise_I = sigloss_final_I.values()
    else:
        sigfactors = sigloss_final.values()
        sigfactors_I = sigloss_final_I.values()

# save final values

other_factors = 1/n.log(2)  # median correction factor
# XXX need to multiply other_factors with values below
fold_factor = file['k']**3/(2*n.pi**2)

split_index = n.argmin(n.abs(kpl))
sigfactors_pos = sigfactors[split_index:]
sigfactors_neg = sigfactors[split_index::-1]
sigfactors_I_pos = sigfactors_I[split_index:]
sigfactors_I_neg = sigfactors_I[split_index::-1]
sigfactors_noise_pos = sigfactors_noise[split_index:]
sigfactors_noise_neg = sigfactors_noise[split_index::-1]
sigfactors_noise_I_pos = sigfactors_noise_I[split_index:]
sigfactors_noise_I_neg = sigfactors_noise_I[split_index::-1]
sigfactors_fold = (n.array(sigfactors_pos) + n.array(sigfactors_neg)) / 2
sigfactors_I_fold = (n.array(sigfactors_I_pos) + n.array(sigfactors_I_neg)) / 2
sigfactors_noise_fold = (n.array(sigfactors_noise_pos)
                         + n.array(sigfactors_noise_neg)) / 2
sigfactors_noise_I_fold = (n.array(sigfactors_noise_I_pos)
                         + n.array(sigfactors_noise_I_neg)) / 2
pIv = pIv*sigfactors_I
pCv = pCv*sigfactors
pIn = pIn*sigfactors_noise_I
pCn = pCn*sigfactors_noise
pIv_err = pIv_err*sigfactors_I
pCv_err = pCv_err*sigfactors
pIn_err = pIn_err*sigfactors_noise_I
pCn_err = pCn_err*sigfactors_noise
pIv_fold = n.abs(file['pIv_fold'])*sigfactors_I_fold*fold_factor
pCv_fold = n.abs(file['pCv_fold'])*sigfactors_fold*fold_factor
pIn_fold = n.abs(file['pIn_fold'])*sigfactors_noise_I_fold*fold_factor
pCn_fold = n.abs(file['pCn_fold'])*sigfactors_noise_fold*fold_factor
pIv_fold_err = file['pIv_fold_err']*sigfactors_I_fold*fold_factor
pCv_fold_err = file['pCv_fold_err']*sigfactors_fold*fold_factor
pIn_fold_err = file['pIn_fold_err']*sigfactors_noise_I_fold*fold_factor
pCn_fold_err = file['pCn_fold_err']*sigfactors_noise_fold*fold_factor

neg_ind = n.where(file['pCv'] < 0)
neg_ind_fold = n.where(file['pCv_fold'] < 0)
neg_ind_noise = n.where(file['pCn'] < 0)
neg_ind_noise_fold = n.where(file['pCn_fold'] < 0)

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
        prob=0.9545, neg_ind=neg_ind, neg_ind_fold=neg_ind_fold,
        neg_ind_noise=neg_ind_noise,
        neg_ind_noise_fold=neg_ind_noise_fold,
        cmd=' '.join(sys.argv))
