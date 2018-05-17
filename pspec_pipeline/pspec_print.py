#! /usr/bin/env python

"""Print tabular data of pspec_final_combine.npz or pspec_pk_k3pk.npz data."""

import numpy as np
import sys
import os
import argparse
# import capo.cosmo_units as cosmo
from capo import sensitivity

parser = argparse.ArgumentParser()
parser.add_argument('files', metavar='<FILE>', type=str, nargs='+',
                    help='List of pspec Files to print')
parser.add_argument('--kmags', type=float, nargs='+',
                    default=[.1, .2, .3, .4, .5],
                    help='List of |k| values to print pspec data.')
args = parser.parse_args()

def r(number):
    return '%.2f' % np.round(number, 2)


def f212z(fq):
    return 1420e6/fq - 1

freqs = []
print "Finding Frequency"
for filename in args.files:
    f = np.load(filename)
    try:
        print "npz...",
        freqs.append(f['freq'])
        print "suceess"
    except(KeyError):
        print 'Fail'
        try:
            print "Finding frequency from freqs array"
            freq_array = f['afreqs']
            freqs.append(np.mean(freq_array))
            print "success"
        except(KeyError):
            print "No Frequency Found. Skipping"


files = args.files
inds = np.argsort(freqs)
files = np.take(files, inds)
freqs = np.take(freqs, inds)
redshifts = f212z(freqs*1e9)  # Frequencies are in GHz
for index, file in enumerate(files):
    redshift = redshifts[index]
    pspec_dict = np.load(file)
    inttime = pspec_dict['frf_inttime']
    cnt = pspec_dict['cnt_eff']
    nbls_g = pspec_dict['nbls_g']
    nlsts = len(pspec_dict['lsts']) * pspec_dict['inttime']
    nlsts /= pspec_dict['frf_inttime']
    if pspec_dict['frf_inttime'] == pspec_dict['inttime']:
        omega_eff = .74**2/.32
        # for capo analytical; from T1 of Parsons FRF paper
    else:
        omega_eff = .74**2/.24
    #print 'Using these parameters for the Analytic Model:'
    #print 'Redshift:', redshift
    #print '\tT_int:', inttime
    #print '\tNbls:', pspec_dict['nbls']
    #print '\tNgps:', pspec_dict['ngps']
    #print '\tNdays:', cnt
    #print '\tNlstbins:', nlsts_g
    if 'theory_noise' and 'theory_noise_delta2' not in pspec_dict:
        S = sensitivity.Sense()
        S.z = redshift

        #   Tsys
        # S.Tsys = 551e3  #set to match 21cmsense exactly
        # S.Tsys = 505e3 #Ali et al, at 164MHz
        # set to match noise realization
        S.Tsys = (200 + 180.*(freqs[index]/.180)**-2.55)*1e3
        #print "Tsys = ", S.Tsys

        S.t_int = inttime
        S.Ndays = cnt  # effective number of days
        S.Npols = 2
        S.Nseps = 3
        S.Nblgroups = pspec_dict['ngps']
        S.Omega_eff = omega_eff
        # use the FRF weighted beams listed in
        # T1 of Parsons etal beam sculpting paper
        k = pspec_dict['k']
        S.Nbls = pspec_dict['nbls']
        S.Nlsthours = nlsts
        S.calc()

        pk_noise = {k:S.P_N for k in pspec_dict['k']}
        delta2_noise = {k:S.Delta2_N(k) for k in pspec_dict['k']}
    else:
        pk_noise = {k:pspec_dict['theory_noise'][kk] for kk,k in enumerate(pspec_dict['k'])}
        delta2_noise = {k:pspec_dict['theory_noise_delta2'][kk] for kk,k in enumerate(pspec_dict['k'])}
    #print 'Beggining Table Data: '
    #print '\n'
    for k in args.kmags:
        k_ind = np.argmin(np.abs(k - pspec_dict['k']))
        k_val = pspec_dict['k'][k_ind]
        try: # find from signal loss results
            delta2_value = pspec_dict['pI_fold'][k_ind]
            bootstrap_error = pspec_dict['pI_fold_up'][k_ind]
        except: # find from pspec_2d_to_1d results
            fold_factor = pspec_dict['k'][k_ind]**3/(2*np.pi**2)
            delta2_value = pspec_dict['pIv_fold_old'][k_ind]*fold_factor
            bootstrap_error = pspec_dict['pIv_fold_err_old'][k_ind]*2*fold_factor

        # Formatting here for a LaTEX table
        # Need to print:
        # k, redshift, delta^2, theoretical error, theory sigma, (Det/Ulim)
        #                       bootstrap error, bootstrap sigma, (Det/Ulim)
        #from IPython import embed; embed()
        print k, '&',
        print r(redshift), '&',
        print r(delta2_value), '&',
        print '$\pm$ '+  r(delta2_noise[k_val]*2), '&',
        print r(delta2_value/(delta2_noise[k_val])), '&',
        if delta2_value - (delta2_noise[k_val]*2) > 0:
            print 'Det', '&',
        else:
            print 'ULim', '&',
        print '$\pm$ ' + r(bootstrap_error), '&',
        print r(delta2_value/(bootstrap_error/2.)), '&',
        if delta2_value-bootstrap_error > 0:
            print 'Det',
        else:
            print 'Ulim',
        print '\\tabularnewline'
    print '\t\hline \hline'
