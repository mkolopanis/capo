#! /usr/bin/env python
"""Averages pspec data from multiple separations.

Performs inversre variance weighting averaging to input power spectrum files.
"""
import numpy as np
import argparse
import sys
import os

parser = argparse.ArgumentParser(
    description=('Perform inverse variance weighting on'
                 ' input power spectrum files.')
)

parser.add_argument('files', metavar='<FILE>', type=str, nargs='+',
                    help='List of files to combine')
parser.add_argument('--outfile', type=str, default='pspec_final.npz',
                    help='give custom output file name')
args = parser.parse_args()


def perform_sum(dsum=None, weights=None, kpl=None, values=None, errors=None):
    """Perform inverser variance sum.

    Store answer in dictionaries dsum and weights, with kpl as keys.
    """
    if any(arg is None for arg in [kpl, values, errors]):
        return None
    if dsum is None:
        dsum = {}
    if weights is None:
        weights = {}

    for k, d, err in zip(kpl, values, errors):
        dsum[k] = dsum.get(k, 0) + d / err**2
        weights[k] = weights.get(k, 0) + 1. / err**2

    return dsum, weights


print 'Averaging Separations'

flat_power_spectra = [p + x for p in ['pC',  'pI']
                      for x in ['e', 'r', 's', 'v', 'n']]
flat_power_spectra.append('pCr-pCv')
flat_power_spectra.append('pCs-pCn')
flat_power_spectra.append('pIr-pIv')
flat_power_spectra.append('pIs-pIn')
folded_power_spectra = [x + '_fold' for x in flat_power_spectra]
flat_errors = [x + '_err' for x in flat_power_spectra]
folded_errors = [x + '_err' for x in folded_power_spectra]

summed_pspec = {key: {} for key in np.concatenate([flat_power_spectra,
                                                   folded_power_spectra])}
summed_weights = {key: {} for key in np.concatenate([flat_errors,
                                                     folded_errors])}
single_keys = np.concatenate([['kpl', 'k'], flat_power_spectra, flat_errors,
                              folded_power_spectra, folded_errors])
out_dict = {}
kpls, ks = [], []
for filename in args.files:

    f = np.load(filename)
    # This generator object is all the additional keywords we want to
    # carry through the power spectrum but do not want to average
    # like freq, afreqs, the probability:
    # if these aren't the same you should not be
    # averaging these quantities anyway
    # For all these extra key words, we will keep all of the info for now.
    generator = (x for x in f.keys() if x not in single_keys)
    for key in generator:
        if key not in out_dict.keys():
            out_dict[key] = [f[key]]
        else:
            out_dict[key].append(f[key])

    kpls.append(f['kpl'])
    ks.append(f['k'])
    # Sum all the different pspecs and weights and accumulate in dictionaries
    # specified by the k or kpl values
    for key1, key2 in zip(flat_power_spectra, flat_errors):
        try:
            temp_sum, temp_weight = perform_sum(kpl=f['kpl'], values=f[key1],
                                                errors=f[key2])

            for _k in f['kpl']:
                summed_pspec[key1][_k] = (summed_pspec[key1].get(_k, 0)
                                          + temp_sum[_k])
                summed_weights[key2][_k] = (summed_weights[key2].get(_k, 0)
                                            + temp_weight[_k])
        except(KeyError):
            print 'Cannot find at least one of the',
            print ' keys in your pspec data:', key1, key2
            pass

    for key1, key2 in zip(folded_power_spectra, folded_errors):

        try:
            temp_sum, temp_weight = perform_sum(kpl=f['k'], values=f[key1],
                                                errors=f[key2])
            for _k in f['k']:
                summed_pspec[key1][_k] = (summed_pspec[key1].get(_k, 0)
                                          + temp_sum[_k])
                summed_weights[key2][_k] = (summed_weights[key2].get(_k, 0)
                                            + temp_weight[_k])
        except(KeyError):
            print 'Cannot find at least one of the',
            print ' keys in your pspec data:', key1, key2
            pass

# only keep one copy of the unique k (kpl _fold) and freq (afreqs) values
out_dict['freq'] = np.unique(out_dict['freq'])
out_dict['afreqs'] = np.unique(out_dict['afreqs'])
out_dict['kpl_fold'] = np.unique(out_dict['kpl_fold'])
kpl = np.unique(kpls)
k = np.unique(ks)
# these two statements might seem silly but 0. is in kpl and it becomes
# 0.0 durig the for _k in kpl, so this helps makes sure all keys are the same
kpl = [_k for _k in kpl]
k = [_k for _k in k]

out_dict['k'] = np.array(k)
out_dict['kpl'] = np.array(kpl)

if np.size(out_dict['cmd']) == 1:
    out_dict['cmd'] = out_dict['cmd'].item() + ' \n ' + ' '.join(sys.argv)
else:
    full_cmd = np.concatenate([[_c] for _c in out_dict['cmd']])
    out_dict['cmd'] = ' '.join(full_cmd) + ' \n ' + ' '.join(sys.argv)

for key1, key2 in zip(flat_power_spectra, flat_errors):
    out_dict[key1] = np.array([summed_pspec[key1][_k]
                               / summed_weights[key2][_k] for _k in kpl])

    out_dict[key2] = np.array([1. / np.sqrt(summed_weights[key2][_k])
                               for _k in kpl])

for key1, key2 in zip(folded_power_spectra, folded_errors):
    out_dict[key1] = np.array([summed_pspec[key1][_k]
                               / summed_weights[key2][_k] for _k in k])

    out_dict[key2] = np.array([1. / np.sqrt(summed_weights[key2][_k])
                               for _k in k])
# should we average to closest k_perp value?
ks = np.sqrt(np.min(out_dict['kperp'])**2 + out_dict['kpl_fold']**2)
out_dict['k'] = ks  # rest k values to be the bins
digitized = np.digitize(out_dict['k'], ks)

# take smallest to for the bins
ks = np.sqrt(np.min(out_dict['kperp'])**2 + out_dict['kpl_fold']**2)
digitized = np.digitize(out_dict['k'], ks)

for key1, key2 in zip(folded_power_spectra, folded_errors):
    out_dict[key1] = np.array([np.sum(out_dict[key1][digitized == _k]
                               / out_dict[key2][digitized == _k]**2)
                               for _k in xrange(1, len(ks)+1)])

    out_dict[key2] = np.array([np.sum(1.
                               / out_dict[key2][digitized == _k]**2)
                               for _k in xrange(1, len(ks) + 1)])
    # renormalized summed pspec values and uncertainties
    out_dict[key1] = out_dict[key1]/out_dict[key2]
    out_dict[key2] = np.sqrt(1./out_dict[key2])
# set output k values to mean over kperp
out_dict['k'] = np.sqrt(np.mean(out_dict['kperp'])**2
                        + out_dict['kpl_fold']**2)
gen = (x for x in out_dict if x not in np.concatenate([single_keys,
                                                       ['cmd', 'afreqs', 'k',
                                                        'kpl', 'kpl_fold']]))
for key in gen:
    if np.size(out_dict[key]) == 1:
        continue
    if isinstance(np.squeeze(out_dict[key])[0], str):
        continue
    out_dict[key] = np.mean(out_dict[key], axis=0).squeeze()

print 'Saving output to: '+args.outfile
np.savez(args.outfile, **out_dict)
