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
kperps = []

# Read in PS for different seps
for file in args.files: # XXX only compatible with pspec_pk_k3pk.npz files 
    f = n.load(file)
    metakeys = [key for key in f.keys() if not (key.startswith('p') or key.startswith('k'))] # all keys except PS values and k,kpl,kpl_fold,kperp
    kperps.append(f['kperp'])
    for kk, k in enumerate(f['kpl']):
        pCv.setdefault(k,[]).append(f['pCv'][kk])
        pCv_err.setdefault(k,[]).append(f['pCv_err'][kk])
        pIv.setdefault(k,[]).append(f['pIv'][kk])
        pIv_err.setdefault(k,[]).append(f['pIv_err'][kk])
        pCn.setdefault(k,[]).append(f['pCn'][kk])
        pCn_err.setdefault(k,[]).append(f['pCn_err'][kk])
        pIn.setdefault(k,[]).append(f['pIn'][kk])
        pIn_err.setdefault(k,[]).append(f['pIn_err'][kk])
    for kk, k in enumerate(f['kpl_fold']):
        pCv_fold.setdefault(k,[]).append(f['pCv_fold'][kk])
        pCv_fold_err.setdefault(k,[]).append(f['pCv_fold_err'][kk])
        pIv_fold.setdefault(k,[]).append(f['pIv_fold'][kk])
        pIv_fold_err.setdefault(k,[]).append(f['pIv_fold_err'][kk])
        pCn_fold.setdefault(k,[]).append(f['pCn_fold'][kk])
        pCn_fold_err.setdefault(k,[]).append(f['pCn_fold_err'][kk])
        pIn_fold.setdefault(k,[]).append(f['pIn_fold'][kk])
        pIn_fold_err.setdefault(k,[]).append(f['pIn_fold_err'][kk])

# Combine

def combine_PS(ks, vals, errs):
    values = []
    err_values = []
    for kk, k in enumerate(ks):
        var = n.array(errs[k])**2 # square for var
        if args.weight == 'inv_var': weights=1/var
        elif args.weight == 'var': weights=var
        else: weights = n.ones_like(var)
        value = n.average(vals[k], weights=weights)
        err_value = n.average(var, weights=weights)
        values.append(value)
        err_values.append(n.sqrt(err_value)) # back to std
    return values, err_values

pCv_combine, pCv_err_combine = combine_PS(f['kpl'], pCv, pCv_err)
pIv_combine, pIv_err_combine = combine_PS(f['kpl'], pIv, pIv_err)
pCv_fold_combine, pCv_fold_err_combine = combine_PS(f['kpl_fold'], pCv_fold, pCv_fold_err)
pIv_fold_combine, pIv_fold_err_combine = combine_PS(f['kpl_fold'], pIv_fold, pIv_fold_err)
pCn_combine, pCn_err_combine = combine_PS(f['kpl'], pCn, pCn_err)
pIn_combine, pIn_err_combine = combine_PS(f['kpl'], pIn, pIn_err)
pCn_fold_combine, pCn_fold_err_combine = combine_PS(f['kpl_fold'], pCn_fold, pCn_fold_err)
pIn_fold_combine, pIn_fold_err_combine = combine_PS(f['kpl_fold'], pIn_fold, pIn_fold_err)

# Save
kperp = n.mean(kperps) # XXX averaged together k_perps
k = n.sqrt(f['kpl_fold']**2 + kperp**2) # for folded
kpl=f['kpl']
kpl_fold=f['kpl_fold']
final_file = 'pspec_final_combine.npz'

metadata = {}
for metakey in metakeys: metadata[metakey] = f[metakey]
metadata['cmd'] = metadata['cmd'].item() + '\n ' + ' '.join(sys.argv)

# TO DO: over-write seps if wanted

n.savez(final_file, k=k, kperp=kperp, kpl=kpl, kpl_fold=kpl_fold,
        pCv=pCv_combine, pCv_err=pCv_err_combine,
        pIv=pIv_combine, pIv_err=pIv_err_combine, 
        pCn=pCn_combine, pCn_err=pCn_err_combine, 
        pIn=pIn_combine, pIn_err=pIn_err_combine,
        pCv_fold=pCv_fold_combine, pCv_fold_err=pCv_fold_err_combine,
        pIv_fold=pIv_fold_combine, pIv_fold_err=pIv_fold_err_combine, 
        pCn_fold=pCn_fold_combine, pCn_fold_err=pCn_fold_err_combine, 
        pIn_fold=pIn_fold_combine, pIn_fold_err=pIn_fold_err_combine, 
        **metadata)




