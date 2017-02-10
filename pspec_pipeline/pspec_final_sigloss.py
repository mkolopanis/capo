#! /usr/bin/env python

import numpy as n, pylab as p
from matplotlib import gridspec
from scipy.interpolate import interp1d
from capo.pspec import f2z
import glob
import optparse, sys

# Signal loss curve for each k

fig = p.figure(1, figsize=(15,7))
pklo,pkhi = 1e1,1e12
Pin_points = {} # for interpolator
Pout_points = {}
Pout_noise_points = {}

file = n.load(glob.glob('inject_*')[0]+'/pspec_pk_k3pk.npz') # read first inject file only
pCv = n.abs(file['pCv']); pCv_err = file['pCv_err'] 
pIv = n.abs(file['pIv']); pIv_err = file['pIv_err']
pCn = n.abs(file['pCn']); pCn_err = file['pCn_err']
pIn = n.abs(file['pIn']); pIn_err = file['pIn_err'] # absolute values only used for signal loss estimation - both positive/negative values are saved in npz file lalr
for inject in glob.glob('inject_*'):
    print 'Reading', inject
    file = n.load(inject + '/pspec_pk_k3pk.npz')
    kpl = file['kpl']
    Pout = n.abs(file['pCr-pCv']); Pout_err = file['pCr-pCv_err']
    Pout_noise = n.abs(file['pCs-pCn']); Pout_noise_err = file['pCs-pCn_err']
    Pin = n.abs(file['pIe']); pIe_err = file['pIe_err']
    for ind in range(len(kpl)): # plot for each k
        p.figure(1)
        p.subplot(3,7,ind)
        try:
            Pin_points[kpl[ind]].append(Pin[ind])
            Pout_points[kpl[ind]].append(Pout[ind]) 
            Pout_noise_points[kpl[ind]].append(Pout_noise[ind])
        except:
            Pin_points[kpl[ind]] = [Pin[ind]]
            Pout_points[kpl[ind]] = [Pout[ind]]
            Pout_noise_points[kpl[ind]] = [Pout_noise[ind]]
        p.loglog(Pin[ind],Pout[ind],'k.') # points
        #p.loglog(Pin[ind],Pout_noise[ind],'b.') # noise points
        p.loglog([pklo,pkhi],[pklo,pkhi], 'k-') # diagonal line
        p.axhline(y=pCv[ind]+2*pCv_err[ind], color='grey', linewidth=1) # 2-sigma upper limit pCv horizontal line (plots for every injection, but should be exactly the same for every injection)
        #p.axhline(y=pCn[ind]+2*pCn_err[ind], color='cyan', linewidth=1) # 2-sigma upper limit pCn horizontal line (plots for every injection, but should be exactly the same for every injection)
        p.grid(True)
        p.xlim(pklo,pkhi)
        p.ylim(pklo,pkhi)
        p.tick_params(axis='both', which='both', labelsize=6)
        p.title('kpl = '+str(kpl[ind]),fontsize=8)

p.tight_layout()
fig.text(0.5, 0.0, r'$P_{\rm in}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$', ha='center')
fig.text(0.0, 0.5, r'$P_{\rm out}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$', va='center', rotation='vertical')
print 'Saving sigloss.png'
fig.savefig('sigloss.png', format='png')

# Obtain signal loss factors vs. k

sigfactors = []
sigfactors_noise = []
for ind in range(len(kpl)): # interpolation for signal loss factors for each k
    order = n.argsort(Pout_points[kpl[ind]]) # build interpolator for data
    Pout_points[kpl[ind]] = n.array(Pout_points[kpl[ind]])[order]
    Pin_points[kpl[ind]] = n.array(Pin_points[kpl[ind]])[order]
    sig_factor_interp = interp1d(Pout_points[kpl[ind]], n.array(Pin_points[kpl[ind]])/n.array(Pout_points[kpl[ind]]),kind='linear',bounds_error=False,fill_value=0)
    sigfactors.append(sig_factor_interp(pCv[ind]+2*pCv_err[ind])) # 2-sigma upper limit for pCv
    order = n.argsort(Pout_noise_points[kpl[ind]]) # build interpolator for noise
    Pout_noise_points[kpl[ind]] = n.array(Pout_noise_points[kpl[ind]])[order]
    Pin_points[kpl[ind]] = n.array(Pin_points[kpl[ind]])[order]
    sig_factor_interp_noise = interp1d(Pout_noise_points[kpl[ind]], n.array(Pin_points[kpl[ind]])/n.array(Pout_noise_points[kpl[ind]]),kind='linear',bounds_error=False,fill_value=0)
    sigfactors_noise.append(sig_factor_interp_noise(pCn[ind]+2*pCn_err[ind])) # 2-sigma upper limit for pCn

# Signal loss factor vs. k

p.figure(2)
for ind in range(len(kpl)):
    p.plot(kpl[ind],sigfactors[ind],'k.',label='Data' if ind==0 else "")
    p.plot(kpl[ind],sigfactors_noise[ind],'b.',label='Noise' if ind==0 else "")
    p.xlabel('kpl')
    p.ylabel('Signal Loss Factor')
p.legend()

split_index = n.argmin(n.abs(kpl))
sigfactors_pos = sigfactors[split_index:]
sigfactors_neg = sigfactors[split_index::-1]
sigfactors_noise_pos = sigfactors_noise[split_index:]
sigfactors_noise_neg = sigfactors_noise[split_index::-1]
sigfactors_fold = (n.array(sigfactors_pos)+n.array(sigfactors_neg))/2
sigfactors_noise_fold = (n.array(sigfactors_noise_pos)+n.array(sigfactors_noise_neg))/2
fold_factor = file['k']**3/(2*n.pi**2)

# Vanilla 2-sigma power spectra plots
"""
p.figure(3)
p.subplot(121)
p.plot(kpl,pIv+2*pIv_err,'k^')
p.errorbar(kpl,pCv*sigfactors,yerr=2*pCv_err*sigfactors,fmt='k.',capsize=3,linewidth=1)
#p.plot(kpl,pIn+2*pIn_err,'b^')
#p.errorbar(kpl,pCn*sigfactors_noise,yerr=2*pCn_err*sigfactors_noise,fmt='b.',capsize=3,linewidth=1)
p.yscale('log',nonposy='clip')
p.ylim(1e-1,1e13)
p.title('P(k)')
p.subplot(122)
p.plot(file['k'],file['pIv_fold']*fold_factor+2*file['pIv_fold_err']*fold_factor,'k^',label='pIv 2-sigma upper limit')
p.errorbar(file['k'],file['pCv_fold']*sigfactors_fold*fold_factor,yerr=2*file['pCv_fold_err']*sigfactors_fold*fold_factor,fmt='k.',capsize=3,linewidth=1,label='pCv')
#p.plot(file['k'],file['pIn_fold']*fold_factor+2*file['pIn_fold_err']*fold_factor,'b^',label='pIn 2-sigma upper limit')
#p.errorbar(file['k'],file['pCn_fold']*sigfactors_noise_fold*fold_factor,yerr=2*file['pCn_fold_err']*sigfactors_fold*fold_factor,fmt='b.',capsize=3,linewidth=1,label='pCn')
p.yscale('log',nonposy='clip')
p.ylim(1e-1,1e13)
p.title('Delta^2(k)')
p.legend(prop={'size':8})
"""
p.show()

# Save values

other_factors = 1/n.log(2) # median correction factor
# XXX need to multiply other_factors with values below

pIv = pIv
pCv = pCv*sigfactors
pIn = pIn
pCn = pCn*sigfactors_noise
pIv_err = pIv_err
pCv_err = pCv_err*sigfactors
pIn_err = pIn_err
pCn_err = pCn_err*sigfactors_noise
pIv_fold = n.abs(file['pIv_fold'])*fold_factor
pCv_fold = n.abs(file['pCv_fold'])*sigfactors_fold*fold_factor
pIn_fold = n.abs(file['pIn_fold'])*fold_factor
pCn_fold = n.abs(file['pCn_fold'])*sigfactors_noise_fold*fold_factor
pIv_fold_err = file['pIv_fold_err']*fold_factor
pCv_fold_err = file['pCv_fold_err']*sigfactors_fold*fold_factor
pIn_fold_err = file['pIn_fold_err']*fold_factor
pCn_fold_err = file['pCn_fold_err']*sigfactors_noise_fold*fold_factor

neg_ind = n.where(file['pCv']<0)
neg_ind_fold = n.where(file['pCv_fold']<0)
neg_ind_noise = n.where(file['pCn']<0)
neg_ind_noise_fold = n.where(file['pCn_fold']<0)

print '   Saving pspec_final.npz' # XXX 2-sigma probability is hard-coded
n.savez('pspec_final.npz', kpl=kpl, k=file['k'], freq=file['freq'],
                        pC=pCv, pC_up=2*pCv_err, 
                        pC_fold=pCv_fold, pC_fold_up=2*pCv_fold_err,
                        pI=pIv, pI_up=2*pIv_err,
                        pI_fold=pIv_fold, pI_fold_up=2*pIv_fold_err,
                        pCn=pCn, pCn_up=2*pCn_err,
                        pCn_fold=pCn_fold, pCn_fold_up=2*pCn_fold_err,
                        pIn=pIn, pIn_up=2*pIn_err,
                        pIn_fold=pIn_fold, pIn_fold_up=2*pIn_fold_err,
                        prob=0.9545, neg_ind=neg_ind, neg_ind_fold=neg_ind_fold,
                        neg_ind_noise=neg_ind_noise, 
                        neg_ind_noise_fold=neg_ind_noise_fold,
                        cmd=' '.join(sys.argv))


 

 
