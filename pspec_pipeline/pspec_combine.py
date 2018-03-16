#! /usr/bin/env python

import numpy as n
import argparse
import sys, os

parser = argparse.ArgumentParser()
parser.add_argument('files', metavar='<FILE>', type=str, nargs='+',
                    help='List of files to combine.')
parser.add_argument('--weight', type=str, default='uniform',
                    choices=['var','inv_var','uniform'],
                    help='Weights for combining PS data.')
args = parser.parse_args()

# Final dictionaries for all PS files
pCv = {}
pCv_err = {}
pIv = {}
pIv_err = {}
pCv_fold = {}
pCv_fold_err = {}
pIv_fold = {}
pIv_fold_err = {}
pCn = {}
pCn_err = {}
pIn = {}
pIn_err = {}
pCn_fold = {}
pCn_fold_err = {}
pIn_fold = {}
pIn_fold_err = {}

pCv_old = {}
pCv_err_old = {}
pIv_old = {}
pIv_err_old = {}
pCv_fold_old = {}
pCv_fold_err_old = {}
pIv_fold_old = {}
pIv_fold_err_old = {}
pCn_old = {}
pCn_err_old = {}
pIn_old = {}
pIn_err_old = {}
pCn_fold_old = {}
pCn_fold_err_old = {}
pIn_fold_old = {}
pIn_fold_err_old = {}


kperps = []
seps = []
theory_noise = {}
theoury_noise_fold = {}
# These dummy values are to set the "vals" = 0 when combining noise variance
dummy_noise = {}
dummy_noise_fold = {}
# Read in PS for different seps
for file in args.files:  # XXX only compatible with pspec_pk_k3pk.npz files
    print 'Reading', file
    f = n.load(file)
    metakeys = [key for key in f.keys() if not (key.startswith('p') or key.startswith('k') or k.startswith('theory'))]  # all keys except PS values and k,kpl,kpl_fold,kperp
    kperps.append(f['kperp'])
    seps.append(f['sep'])
    for kk, k in enumerate(f['kpl']):  # each dictionary has keys of k's and contains a list of values for that k (3 if there are 3 seps)
        pCv.setdefault(k, []).append(f['pCv'][kk])
        pCv_err.setdefault(k, []).append(f['pCv_err'][kk])
        pIv.setdefault(k, []).append(f['pIv'][kk])
        pIv_err.setdefault(k, []).append(f['pIv_err'][kk])
        pCn.setdefault(k, []).append(f['pCn'][kk])
        pCn_err.setdefault(k, []).append(f['pCn_err'][kk])
        pIn.setdefault(k, []).append(f['pIn'][kk])
        pIn_err.setdefault(k, []).append(f['pIn_err'][kk])

        pCv_old.setdefault(k, []).append(f['pCv_old'][kk])
        pCv_err_old.setdefault(k, []).append(f['pCv_err_old'][kk])
        pIv_old.setdefault(k, []).append(f['pIv_old'][kk])
        pIv_err_old.setdefault(k, []).append(f['pIv_err_old'][kk])
        pCn_old.setdefault(k, []).append(f['pCn_old'][kk])
        pCn_err_old.setdefault(k, []).append(f['pCn_err_old'][kk])
        pIn_old.setdefault(k, []).append(f['pIn_old'][kk])
        pIn_err_old.setdefault(k, []).append(f['pIn_err_old'][kk])

        theory_noise.setdefault(k, []).append(f['theory_noise'][kk])
        dummy_noise.setdefault(k, []).append(0)

    for kk, k in enumerate(f['kpl_fold']):
        pCv_fold.setdefault(k, []).append(f['pCv_fold'][kk])
        pCv_fold_err.setdefault(k, []).append(f['pCv_fold_err'][kk])
        pIv_fold.setdefault(k, []).append(f['pIv_fold'][kk])
        pIv_fold_err.setdefault(k, []).append(f['pIv_fold_err'][kk])
        pCn_fold.setdefault(k, []).append(f['pCn_fold'][kk])
        pCn_fold_err.setdefault(k, []).append(f['pCn_fold_err'][kk])
        pIn_fold.setdefault(k, []).append(f['pIn_fold'][kk])
        pIn_fold_err.setdefault(k, []).append(f['pIn_fold_err'][kk])

        pCv_fold_old.setdefault(k, []).append(f['pCv_fold_old'][kk])
        pCv_fold_err_old.setdefault(k, []).append(f['pCv_fold_err_old'][kk])
        pIv_fold_old.setdefault(k, []).append(f['pIv_fold_old'][kk])
        pIv_fold_err_old.setdefault(k, []).append(f['pIv_fold_err_old'][kk])
        pCn_fold_old.setdefault(k, []).append(f['pCn_fold_old'][kk])
        pCn_fold_err_old.setdefault(k, []).append(f['pCn_fold_err_old'][kk])
        pIn_fold_old.setdefault(k, []).append(f['pIn_fold_old'][kk])
        pIn_fold_err_old.setdefault(k, []).append(f['pIn_fold_err_old'][kk])

        theory_noise_delta2.setdefault(k, []).append(f['theory_noise_delta2'][kk])
        dummy_noise_fold.setdefault(k, []).append(0)


# Combine

def combine_PS(ks, vals, errs):
    """Average vals and errs with given weighting scheme."""
    values = []
    err_values = []
    for kk, k in enumerate(ks):
        var = n.array(errs[k])**2  # square for var
        if args.weight == 'inv_var': weights = 1/var
        elif args.weight == 'var': weights = var
        else: weights = n.ones_like(var)
        value = n.average(vals[k], weights=weights)
        err_value = n.average(var, weights=weights)
        values.append(value)
        err_values.append(n.sqrt(err_value))  # back to std
    return values, err_values


print 'Using', args.weight, 'weighting'
pCv_combine, pCv_err_combine = combine_PS(f['kpl'], pCv, pCv_err)
pIv_combine, pIv_err_combine = combine_PS(f['kpl'], pIv, pIv_err)
pCv_fold_combine, pCv_fold_err_combine = combine_PS(f['kpl_fold'], pCv_fold,
                                                    pCv_fold_err)
pIv_fold_combine, pIv_fold_err_combine = combine_PS(f['kpl_fold'], pIv_fold,
                                                    pIv_fold_err)
pCn_combine, pCn_err_combine = combine_PS(f['kpl'], pCn, pCn_err)
pIn_combine, pIn_err_combine = combine_PS(f['kpl'], pIn, pIn_err)
pCn_fold_combine, pCn_fold_err_combine = combine_PS(f['kpl_fold'], pCn_fold,
                                                    pCn_fold_err)
pIn_fold_combine, pIn_fold_err_combine = combine_PS(f['kpl_fold'], pIn_fold,
                                                    pIn_fold_err)
pCv_old_combine, pCv_err_old_combine = combine_PS(f['kpl'], pCv_old,
                                                  pCv_err_old)
pIv_old_combine, pIv_err_old_combine = combine_PS(f['kpl'], pIv_old,
                                                  pIv_err_old)
pCv_fold_old_combine, pCv_fold_err_old_combine = combine_PS(f['kpl_fold'], pCv_fold_old, pCv_fold_err_old)
pIv_fold_old_combine, pIv_fold_err_old_combine = combine_PS(f['kpl_fold'], pIv_fold_old, pIv_fold_err_old)
pCn_old_combine, pCn_err_old_combine = combine_PS(f['kpl'], pCn_old, pCn_err_old)
pIn_old_combine, pIn_err_old_combine = combine_PS(f['kpl'], pIn_old, pIn_err_old)
pCn_fold_old_combine, pCn_fold_err_old_combine = combine_PS(f['kpl_fold'], pCn_fold_old, pCn_fold_err_old)
pIn_fold_old_combine, pIn_fold_err_old_combine = combine_PS(f['kpl_fold'], pIn_fold_old, pIn_fold_err_old)



_, theory_noise_combine = combine_PS(f['kpl'], dummy_noise, theory_noise)
_, theory_noise_delta2_combine = combine_PS(f['kpl'], dummy_noise_fold,
                                            theory_noise_delta2)
# Save
kperp = n.mean(kperps)  # XXX averaged together k_perps
k = n.sqrt(f['kpl_fold']**2 + kperp**2)  # for folded
kpl = f['kpl']
kpl_fold = f['kpl_fold']
final_file = 'pspec_final_combine.npz'

metadata = {}
for metakey in metakeys: metadata[metakey] = f[metakey]
metadata['cmd'] = metadata['cmd'].item() + '\n ' + ' '.join(sys.argv)
metadata['sep'] = seps

print 'Saving', final_file
n.savez(final_file, k=k, kperp=kperp, kpl=kpl, kpl_fold=kpl_fold,
        pCv=pCv_combine, pCv_err=pCv_err_combine,
        pIv=pIv_combine, pIv_err=pIv_err_combine,
        pCn=pCn_combine, pCn_err=pCn_err_combine,
        pIn=pIn_combine, pIn_err=pIn_err_combine,
        pCv_fold=pCv_fold_combine, pCv_fold_err=pCv_fold_err_combine,
        pIv_fold=pIv_fold_combine, pIv_fold_err=pIv_fold_err_combine,
        pCn_fold=pCn_fold_combine, pCn_fold_err=pCn_fold_err_combine,
        pIn_fold=pIn_fold_combine, pIn_fold_err=pIn_fold_err_combine,
        pIv_fold=pIv_fold_combine, pIv_fold_err=pIv_fold_err_combine,
        pCn_fold=pCn_fold_combine, pCn_fold_err=pCn_fold_err_combine,
        pIn_fold=pIn_fold_combine, pIn_fold_err=pIn_fold_err_combine,
        pCv_old=pCv_old_combine, pCv_err_old=pCv_err_old_combine,
        pIv_old=pIv_old_combine, pIv_err_old=pIv_err_old_combine,
        pCn_old=pCn_old_combine, pCn_err_old=pCn_err_old_combine,
        pIn_old=pIn_old_combine, pIn_err_old=pIn_err_old_combine,
        pCv_fold_old=pCv_fold_old_combine,
        pCv_fold_err_old=pCv_fold_err_old_combine,
        pIv_fold_old=pIv_fold_old_combine,
        pIv_fold_err_old=pIv_fold_err_old_combine,
        pCn_fold_old=pCn_fold_old_combine,
        pCn_fold_err_old=pCn_fold_err_old_combine,
        pIn_fold_old=pIn_fold_old_combine,
        pIn_fold_err_old=pIn_fold_err_old_combine,
        theory_noise=theory_noise, theory_noise_delta2=theory_noise_delta2,
        nseps=len(pCv[pCv.keys()[0]]), **metadata)
