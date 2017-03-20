#! /usr/bin/env python
"""
Compute confidence levels for a given power spectrum with injects.

Input a list of signal loss power spectra (as output by pspec_2d_to_1d).
compute a loss-calibrated limit
output as a new pspec file.
"""
import matplotlib as mpl
from glob import glob
import argparse
import os
from capo.eor_results import read_bootstraps_dcj, average_bootstraps
from capo.pspec import dk_du
from capo import cosmo_units
import numpy as np
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

parser = argparse.ArgumentParser(
    description='Calculate a limit from a range of injected pspecs')
parser.add_argument('pspec_files', metavar='pspsec.npz',
                    type=str, nargs='+', default='.',
                    help='Directory containing injections')
parser.add_argument('--outfile', type=str, default='',
                    help='Specifically specify out directory.')
args = parser.parse_args()


if args.outfile and not args.outfile.endswith('/'):
    args.outfile += '/'


def G(x, mx, dx):
    """Compute Gaussian at x with mean mx and std dx."""
    return 1 / (dx * np.sqrt(2 * np.pi)) * np.exp(-1 / 2. * ((x - mx) / dx)**2)


def G_mc(x_lim, dist_x, dist_sig, Ntrials=10000):
    """
    Compute Monte-Carlo integral between x_lim and dist_x, dist_sig.

    Calculate the probability of drawing a point larger than
    x_lim from the gaussian distribution dist_x +/- dist_sig
    """
    mc_x = np.random.normal(dist_x, dist_sig, size=Ntrials)
    return np.sum(mc_x > x_lim) / float(Ntrials)


def injpath_to_injlevel(filepath):
    """Find injection levels from filename."""
    injdir = os.path.split(os.path.dirname(filepath))[-1]
    return float(injdir.split('_')[-1])


def fit_prob_to_tanh(probs, pIs):
    """Fit tanh to probability curve."""
    # tanh fit is VERY sensitive to initial conditions
    p0 = np.zeros(4)
    mean_ind = np.argmin(abs(probs - .5))
    p0[0] = pIs[mean_ind]  # set mean and width to entral value of pspecs
    p0[1] = pIs[mean_ind]  # set mean and width to central value of pspecs
    p0[2] = .25 + np.min(probs)  # offest of guess is min value + .25
    p0[3] = (np.max(probs) - np.min(probs)) * 3 / 4.
    # scale is  3/4 differnce in probabilities
    parms, cov = curve_fit(tanh, pIs, probs, p0=p0, maxfev=50000)
    return parms


def get_pk_k3pk(prob, ks, parms):
    """Compute pk k3pk from prob and ks."""
    nk = len(ks)
    upper_limits = []
    upper_limits.append([atanh(prob, *ps) for ps in parms])
    upper_limits = np.array(upper_limits).squeeze()
    k3pk = np.tile(1e-5, nk)
    k3err = upper_limits * k**3 / (2 * np.pi**2)
    Pk = np.tile(1e-5, nk)
    pkerr = upper_limits
    return Pk, pkerr, k3pk, k3err


def tanh(x, m, s, c=1.0, a=1.0):
    """Compute tanh."""
    return c + a * np.tanh((x - m) / (2 * s))


def atanh(p, m, s, c, a):
    """Compute atanh."""
    return 2 * s * np.arctanh((p - c) / a) + m


# read in the bootstrapped pspecs
pspec_channels = ['pCv_fold', 'pCv_fold_err',  # weighted data pspec
                  'pIv_fold', 'pIv_fold_err',  # unweighted data pspec
                  'pCr_fold', 'pCr_fold_err',  # weighted data+inj pspec
                  'pIe_fold', 'pIe_fold_err',  # unweighted inj pspec
                  'pIn_fold', 'pIn_fold_err',  # unweighted noise
                  'pCn_fold', 'pCn_fold_err',  # weighted noise
                  'pCs_fold', 'pCs_fold_err',  # weighted noise + eor
                  'pIs_fold', 'pIs_fold_err',  # unweighte noise + eor
                  'pIn', 'pIn_err',  # unweighted noise
                  'pCv', 'pCv_err',  # weighted data pspec
                  'pIv', 'pIv_err',  # unweighted data pspec
                  'pCr', 'pCr_err',  # weighted data+inj pspec
                  'pIe', 'pIe_err']  # unweighted inj pspec pos and neg kpls
pspecs = {}
flat_power_spectra = [p + x for p in ['pC',  'pI']
                      for x in ['e', 'r', 's', 'v', 'n']]
flat_power_spectra.append('pCr-pCv')
flat_power_spectra.append('pCs-pCn')
folded_power_spectra = [x + '_fold' for x in flat_power_spectra]
flat_errors = [x + '_err' for x in flat_power_spectra]
folded_errors = [x + '_err' for x in folded_power_spectra]

# sort the input files. makes things easier later
injlevels = [injpath_to_injlevel(filename) for filename in args.pspec_files]
fileorder = np.argsort(injlevels)
filenames = [args.pspec_files[i] for i in fileorder]
for filename in filenames:
    print filename
    F = np.load(filename)
    for pspec_channel in pspec_channels:
        F[pspec_channel]
        try:
            pspecs[pspec_channel].append(F[pspec_channel])
        except(KeyError):
            pspecs[pspec_channel] = [F[pspec_channel]]
for pspec_channel in pspec_channels:
    pspecs[pspec_channel] = np.array(pspecs[pspec_channel])
Ninj, Nk = pspecs[pspec_channel].shape
k = F['k']
kpls = F['kpl_fold']
Nk = len(k)
freq = F['freq']
afreqs = F['afreqs']
print "found injects:", Ninj
print "found k bins:", Nk

# this gets all the the other meta-data in the file and passes it through
meta_data = {}
generator = (x for x in F.keys()
             if x not in np.concatenate([['kpl', 'k', 'freq', 'afreqs', 'k'],
                                         flat_power_spectra, flat_errors,
                                         folded_power_spectra, folded_errors]))
for key in generator:
    if key not in meta_data.keys():
        meta_data[key] = [F[key]]
    else:
        meta_data[key].append(f[key])

# limit option #1. the net probability that pC is above pcV
probs_data = np.zeros((Ninj, Nk))
probs_noise = np.zeros((Ninj, Nk))

for inj in xrange(Ninj):
    for k_ind in xrange(Nk):
        lossy_limit = (pspecs['pCv_fold'][0, k_ind] +
                       pspecs['pCv_fold_err'][0, k_ind])
        if k_ind == 5:
            print "lossy_limit: ", lossy_limit
            print "pCr lower limit", (pspecs['pCr_fold'][inj, k_ind]
                                      - pspecs['pCr_fold_err'][inj, k_ind])
        probs_data[inj, k_ind] = G_mc(lossy_limit,  # limit
                                      pspecs['pCr_fold'][inj, k_ind],
                                      pspecs['pCr_fold_err'][inj, k_ind])

        lossy_limit_n = (pspecs['pCn_fold'][0, k_ind]
                         + pspecs['pCn_fold_err'][0, k_ind])
        probs_noise[inj, k_ind] = G_mc(lossy_limit_n,  # limit
                                       pspecs['pCs_fold'][inj, k_ind],
                                       pspecs['pCs_fold_err'][inj, k_ind])

tanh_parms_data = []
tanh_parms_noise = []
for k_ind in xrange(Nk):
    tanh_parms_data.append(fit_prob_to_tanh(probs_data[:, k_ind],
                                            pspecs['pIe_fold'][:, k_ind]))
    tanh_parms_noise.append(fit_prob_to_tanh(probs_noise[:, k_ind],
                                             pspecs['pIe_fold'][:, k_ind]))


figure()
for k_ind in xrange(Nk):
    plt.semilogx(pspecs['pIe_fold'][:, k_ind], probs_data[:, k_ind],
                 '-', label=k[k_ind])
grid()
legend(loc='best')
xlabel('$P_{inj}$')
ylabel('Probability to find $P_{inj}$')
savefig(args.outfile + 'p_inj_prob.png', format='png')

# pI and pIn values are not dependen on probability
pIn, pIn_up = pspecs['pIn_fold'][0, :], pspecs['pIn_fold_err'][0, :]
pIn_fold = pspecs['pIn_fold'][0, :]
pIn_fold_up = pspecs['pIn_fold_err'][0, :]
pIn_fold = k**3 / (2 * np.pi**2) * pIn_fold
pIn_fold_up = k**3 / (2 * np.pi**2) * pIn_fold_up

pI, pI_up = pspecs['pIv_fold'][0, :], pspecs['pIv_fold_err'][0, :]
pI_fold = pspecs['pIv_fold'][0, :]
pI_fold_up = pspecs['pIv_fold_err'][0, :]
pI_fold = k**3 / (2 * np.pi**2) * pI_fold
pI_fold_up = k**3 / (2 * np.pi**2) * pI_fold_up

prob_limits = [.2, .25, .3, .35, .4, .5, .68, .85, .9, .95, .97, .99]

# all of the probabilities will be the same
# generate a list of all the names of probs with and without folds
# assign each name to per in a dict that can be unraveled in the save
names = ['pI', 'pC', 'pCn', 'pIn']

for per in prob_limits:
    pC, pC_up, pC_fold, pC_fold_up = get_pk_k3pk(per, k, tanh_parms_data)
    pCn, pCn_up, pCn_fold, pCn_fold_up = get_pk_k3pk(per, k, tanh_parms_noise)

    print 'Saving: ' + args.outfile,
    print 'pspec_final_confidence_{0:02d}'.format(int(per * 100)) + '.npz'
    np.savez((args.outfile +
             'pspec_final_confidence_{0:02d}'.format(int(per * 100)) + '.npz'),
             freq=freq, k=k, kpl=kpls, afreqs=afreqs, pC=pC, pC_up=pC_up,
             pC_fold=pC_fold, pC_fold_up=pC_fold_up, pI=pI, pI_up=pI_up,
             pI_fold=pI_fold, pI_fold_up=pI_fold_up, pCn=pCn, pCn_up=pCn_up,
             pCn_fold=pCn_fold, pCn_fold_up=pCn_fold_up, pIn=pIn,
             pIn_up=pIn_up, pIn_fold=pIn_fold, pIn_fold_up=pIn_fold_up,
             prob=per, **meta_data)
