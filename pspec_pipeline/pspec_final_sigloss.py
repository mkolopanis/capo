#! /usr/bin/env python
"""Compute Signal Loss factor via Ali et. al 2016 method."""
import numpy as n
import pylab as p
from matplotlib import gridspec
from scipy.interpolate import interp1d
from capo.pspec import f2z
import glob
import optparse
import sys

# Signal loss curve for each k

fig = p.figure(1, figsize=(15, 7))
pklo, pkhi = 1e1, 1e12
Pin_points = {}
Pout_points = {}
Pout_noise_points = {}

# read first inject file only
file = n.load(glob.glob('inject_*')[0]+'/pspec_pk_k3pk.npz')
pCv = n.abs(file['pCv'])
pCv_err = file['pCv_err']
pIv = n.abs(file['pIv'])
pIv_err = file['pIv_err']
pCn = n.abs(file['pCn'])
pCn_err = file['pCn_err']
pIn = n.abs(file['pIn'])
pIn_err = file['pIn_err']
#  absolute values only used for signal loss estimation
#  both positive/negative values are saved in npz file later

for inject in glob.glob('inject_*'):
    print 'Reading', inject
    file = n.load(inject + '/pspec_pk_k3pk.npz')
    kpl = file['kpl']
    Pout = n.abs(file['pCr-pCv'])
    Pout_err = n.abs(file['pCr-pCv_err'])
    Pout_noise = n.abs(file['pCs-pCn'])
    Pout_noise_err = file['pCs-pCn_err']
    Pin = n.abs(file['pIe'])
    pIe_err = file['pIe_err']

    for ind in range(len(kpl)):  # plot for each k
        p.figure(1)
        p.subplot(3, 7, ind)
        try:
            Pin_points[kpl[ind]].append(Pin[ind])
            Pout_points[kpl[ind]].append(Pout[ind])
            Pout_noise_points[kpl[ind]].append(Pout_noise[ind])
        except:
            Pin_points[kpl[ind]] = [Pin[ind]]
            Pout_points[kpl[ind]] = [Pout[ind]]
            Pout_noise_points[kpl[ind]] = [Pout_noise[ind]]
        p.loglog(Pin[ind], Pout[ind], 'k.')  # points
        # p.loglog(Pin[ind],Pout_noise[ind],'b.') # noise points
        p.loglog([pklo, pkhi], [pklo, pkhi], 'k-')  # diagonal line
        p.axhline(y=pCv[ind] + 2 * pCv_err[ind], color='grey', linewidth=1)
        # 2-sigma upper limit pCv horizontal line (plots for every injection,
        # but should be exactly the same for every injection)
        # p.axhline(y=pCn[ind]+2*pCn_err[ind], color='cyan', linewidth=1)
        # 2-sigma upper limit pCn horizontal line (plots for every injection,
        # but should be exactly the same for every injection)
        p.grid(True)
        p.xlim(pklo, pkhi)
        p.ylim(pklo, pkhi)
        p.tick_params(axis='both', which='both', labelsize=6)
        p.title('kpl = '+str(kpl[ind]), fontsize=8)

p.tight_layout()
fig.text(0.5, 0.0, r'$P_{\rm in}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',
         ha='center')
fig.text(0.0, 0.5, r'$P_{\rm out}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',
         va='center', rotation='vertical')
print 'Saving sigloss.png'
fig.savefig('sigloss.png', format='png')

# Obtain signal loss factors vs. k

sigfactors = []
sigfactors_noise = []
for ind in range(len(kpl)):  # interpolation for signal loss factors for each k
    order = n.argsort(Pout_points[kpl[ind]])  # build interpolator for data
    Pout_points[kpl[ind]] = n.array(Pout_points[kpl[ind]])[order]
    Pin_points[kpl[ind]] = n.array(Pin_points[kpl[ind]])[order]
    sig_factor_interp = interp1d(Pout_points[kpl[ind]],
                                 n.array(Pin_points[kpl[ind]])
                                 / n.array(Pout_points[kpl[ind]]),
                                 kind='linear', bounds_error=False,
                                 fill_value=0)
    # 2-sigma upper limit for pCv
    sigfactors.append(sig_factor_interp(pCv[ind] + 2 * pCv_err[ind]))
    # build interpolator for noise
    order = n.argsort(Pout_noise_points[kpl[ind]])
    Pout_noise_points[kpl[ind]] = n.array(Pout_noise_points[kpl[ind]])[order]
    Pin_points[kpl[ind]] = n.array(Pin_points[kpl[ind]])[order]
    sig_factor_interp_noise = interp1d(Pout_noise_points[kpl[ind]],
                                       n.array(Pin_points[kpl[ind]])
                                       / n.array(Pout_noise_points[kpl[ind]]),
                                       kind='linear', bounds_error=False,
                                       fill_value=0)
    # 2-sigma upper limit for pCn
    sigfactors_noise.append(sig_factor_interp_noise(pCn[ind]+2*pCn_err[ind]))

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
sigfactors_fold = (n.array(sigfactors_pos) + n.array(sigfactors_neg)) / 2
sigfactors_noise_fold = (n.array(sigfactors_noise_pos)
                         + n.array(sigfactors_noise_neg)) / 2
fold_factor = file['k']**3/(2*n.pi**2)

# Vanilla 2-sigma power spectra plots
#
# p.figure(3)
# p.subplot(121)
# p.plot(kpl,pIv+2*pIv_err,'k^')
# p.errorbar(kpl,pCv*sigfactors,yerr=2*pCv_err*sigfactors,fmt='k.',capsize=3,linewidth=1)
# #p.plot(kpl,pIn+2*pIn_err,'b^')
# #p.errorbar(kpl,pCn*sigfactors_noise,yerr=2*pCn_err*sigfactors_noise,fmt='b.',capsize=3,linewidth=1)
# p.yscale('log',nonposy='clip')
# p.ylim(1e-1,1e13)
# p.title('P(k)')
# p.subplot(122)
# p.plot(file['k'],file['pIv_fold'] * fold_factor
#        + 2 * file['pIv_fold_err'] * fold_factor, 'k^'
#        label='pIv 2-sigma upper limit')
# p.errorbar(file['k'], file['pCv_fold'] * sigfactors_fold * fold_factor,
#            yerr=2 * file['pCv_fold_err'] * sigfactors_fold * fold_factor,
#            fmt='k.', capsize=3, linewidth=1 ,label='pCv')
# #p.plot(file['k'], file['pIn_fold'] * fold_factor
#         + 2 * file['pIn_fold_err'] * fold_factor,'b^',
#         label='pIn 2-sigma upper limit')
# #p.errorbar(file['k'], file['pCn_fold'] * sigfactors_noise_fold
#             * fold_factor, yerr=2 * file['pCn_fold_err'] * sigfactors_fold
#             * fold_factor, fmt='b.' ,capsize=3, linewidth=1,label='pCn')
# p.yscale('log',nonposy='clip')
# p.ylim(1e-1,1e13)
# p.title('Delta^2(k)')
# p.legend(prop={'size':8})

p.show()

# Save values

other_factors = 1/n.log(2)  # median correction factor

pIv = pIv*other_factors
pCv = pCv*sigfactors*other_factors
pIn = pIn*other_factors
pCn = pCn*sigfactors_noise*other_factors
pIv_err = pIv_err*other_factors
pCv_err = pCv_err*sigfactors*other_factors
pIn_err = pIn_err*other_factors
pCn_err = pCn_err*sigfactors_noise*other_factors
pIv_fold = n.abs(file['pIv_fold'])*fold_factor*other_factors
pCv_fold = n.abs(file['pCv_fold'])*sigfactors_fold*fold_factor*other_factors
pIn_fold = n.abs(file['pIn_fold'])*fold_factor*other_factors
pCn_fold = n.abs(file['pCn_fold'])*sigfactors_noise_fold*fold_factor*other_factors
pIv_fold_err = file['pIv_fold_err']*fold_factor*other_factors
pCv_fold_err = file['pCv_fold_err']*sigfactors_fold*fold_factor*other_factors
pIn_fold_err = file['pIn_fold_err']*fold_factor*other_factors
pCn_fold_err = file['pCn_fold_err']*sigfactors_noise_fold*fold_factor*other_factors

neg_ind = n.where(file['pCv'] < 0)
neg_ind_fold = n.where(file['pCv_fold'] < 0)
neg_ind_noise = n.where(file['pCn'] < 0)
neg_ind_noise_fold = n.where(file['pCn_fold'] < 0)

# Get all the meta-data from the npz file and pass it through
flat_power_spectra = [p + x for p in ['pC',  'pI']
                      for x in ['e', 'r', 's', 'v', 'n']]
flat_power_spectra.append('pCr-pCv')
flat_power_spectra.append('pCs-pCn')
folded_power_spectra = [x + '_fold' for x in flat_power_spectra]
flat_errors = [x + '_err' for x in flat_power_spectra]
folded_errors = [x + '_err' for x in folded_power_spectra]

meta_data = {}
generator = (x for x in F.keys()
             if x not in np.concatenate([['kpl', 'k', 'freq', 'k'],
                                         flat_power_spectra, flat_errors,
                                         folded_power_spectra, folded_errors]))
for key in generator:
    if key not in meta_data.keys():
        meta_data[key] = [F[key]]
    else:
        meta_data[key].append(f[key])


print '   Saving pspec_final_median.npz'  # XXX 2-sigma probability is hard-coded
n.savez('pspec_final_median.npz', kpl=kpl, k=file['k'], freq=file['freq'],
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
