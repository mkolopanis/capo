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


def perform_sum(dsum=None, weights=None, kpl=None, values=None,
                errors=None, vsum=None):
    """Perform cumulative sum on dictionaries of pspec values indexed by k.

    dsum = power spectrum mean or medians
    vsum = power spectrum variance
    weights = thing by which the sum needs to be divided to get an answer
              (this is _not_ done in this function)
    kpl = k values to use as keys
    values = new values to add to dsum
    errors = new values to add to variances

    """
    if any(arg is None for arg in [kpl, values, errors]):
        return None
    if dsum is None:
        dsum = {}
    if weights is None:
        weights = {}
    if vsum is None:
        vsum = {}

    for k, d, err in zip(kpl, values, errors):
        # for every k value, add in new data.
        dsum[k] = dsum.get(k, 0) + d  # / err**2
        weights[k] = weights.get(k, 0) + 1.  # / err**2
        vsum[k] = vsum.get(k, 0) + err**2

    return dsum, weights, vsum


print 'Averaging Separations'

flat_power_spectra = [p + x for p in ['pC', 'pI']
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
summed_vars = {key: {} for key in np.concatenate([flat_errors,
                                                  folded_errors])}

# a list of all the keys we expect from the input files
single_keys = np.concatenate([['kpl', 'k'], flat_power_spectra, flat_errors,
                              folded_power_spectra, folded_errors])
out_dict = {}
kpls, ks, kpl_folds = [], [], []
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
    kpl_folds.append(f['kpl_fold'])
    # Sum all the different pspecs and weights and accumulate in dictionaries
    # keyed to the k
    for pspec_key, error_key in zip(flat_power_spectra, flat_errors):
        p, w, v = perform_sum(kpl=f['kpl'], values=f[pspec_key],
                              errors=f[error_key],
                              dsum=summed_pspec[pspec_key],
                              vsum=summed_vars[error_key],
                              weights=summed_weights[error_key])

        summed_pspec[pspec_key] = p
        summed_weights[error_key] = w
        summed_vars[error_key] = v
        del(p, w, v)
    # sum power spectra for the folded pspecs, keyed to kpl_fold
    for pspec_key, error_key in zip(folded_power_spectra, folded_errors):
        p, w, v = perform_sum(kpl=f['kpl_fold'], values=f[pspec_key],
                              errors=f[error_key],
                              dsum=summed_pspec[pspec_key],
                              vsum=summed_vars[error_key],
                              weights=summed_weights[error_key])
        summed_pspec[pspec_key] = p
        summed_weights[error_key] = w
        summed_vars[error_key] = v
        del(p, w, v)
# only keep one copy of the unique k (kpl _fold) and freq (afreqs) values
out_dict['freq'] = np.unique(out_dict['freq'])
out_dict['afreqs'] = np.unique(out_dict['afreqs'])
out_dict['kpl_fold'] = np.unique(out_dict['kpl_fold'])
kpl = np.unique(kpls)
k = np.unique(ks)
kpl_folds = np.unique(kpl_folds)
# these two statements might seem silly but 0. is in kpl and it becomes
# 0.0 durig the for _k in kpl, so this helps makes sure all keys are the same
kpl = [_k for _k in kpl]
k = [_k for _k in k]
# I refuse to add to thix sillyness

out_dict['k'] = np.array(k)
out_dict['kpl'] = np.array(kpl)

# add this step to the history.
if np.size(out_dict['cmd']) == 1:
    out_dict['cmd'] = out_dict['cmd'].item() + ' \n ' + ' '.join(sys.argv)
else:
    full_cmd = np.concatenate([[_c] for _c in out_dict['cmd']])
    out_dict['cmd'] = ' '.join(full_cmd) + ' \n ' + ' '.join(sys.argv)

# divide the summed data and variances by the weights
for pspec_key, error_key in zip(flat_power_spectra, flat_errors):
    out_dict[pspec_key] = np.array([summed_pspec[pspec_key][_k]
                                    / summed_weights[error_key][_k]
                                    for _k in kpl])
    out_dict[error_key] = np.array([np.sqrt(summed_vars[error_key][_k]
                                    / summed_weights[error_key][_k])
                                    for _k in kpl])

# divide the summed folded data and variances by the weights
for pspec_key, error_key in zip(folded_power_spectra, folded_errors):
    out_dict[pspec_key] = np.array([summed_pspec[pspec_key][_k]
                                    / summed_weights[error_key][_k]
                                    for _k in kpl_folds])
    out_dict[error_key] = np.array([np.sqrt(summed_vars[error_key][_k]
                                    / summed_weights[error_key][_k])
                                    for _k in kpl_folds])
# use the mean baseline length as the kperp
out_dict['k'] = np.sqrt(out_dict['kpl_fold']**2 + np.mean(out_dict['kperp'])**2)

# pass on values from input files
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
