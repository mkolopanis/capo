#! /usr/bin/env python
"""Averages pspec data from multiple separations.

Performs inversre variance weighting averaging to input power spectrum files.
"""
import numpy as np
import argparse
import sys
import os

parser = argparse.ArgumentParser(
 description='Perform inverse variance weighting on input power spetrum files.'
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
flat_power_spetra.append(['pCr-pCv'])
flat_power_spetra.append(['pCs-pCn'])
folded_power_spectra = [x + '_fold' for x in flat_power_spectra]
flat_errors = [x + '_err' for x in flat_power_spectra]
folded_errors = [x + '_err' for x in folded_power_spectra]

summed_pspec = {key: {} for key in np.concatenate([flat_power_spectra,
                                                   folded_power_spectra])}
summed_weights = {key: {} for key in np.concatenate([flat_errors,
                                                     folded_errors])}
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
    generator = (x for x in f.keys() if x not in np.concatenate([['kpl', 'k'],
                 flat_power_spectra, flat_errors,
                 folded_power_spectra, folded_errors]))
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

# only keep one copy of the unique k (kpl) and freq (afreqs) values
out_dict['freq'] = np.unique(out_dict['freq'])
out_dict['afreqs'] = np.unique(out_dict['afreqs'])
kpl = np.unique(kpls)
k = np.unique(ks)
# these two statements might seem silly but 0. is in kpl and it becomes
# 0.0 durig the for _k in kpl, so this helps makes sure all keys are the same
kpl = [_k for _k in kpl]
k = [_k for _k in k]

out_dict['k'] = np.array(k)
out_dict['kpl'] = np.array(kpl)

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

print 'Saving output to: '+args.outfile
np.savez(args.outfile, **out_dict)
