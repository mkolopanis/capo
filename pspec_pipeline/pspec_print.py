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

print "Quick look at best upper lims:"
for redshift,filename in zip(redshifts,files):
    pspec_dict = np.load(filename)
    try: # find from signal loss results
        delta2_value = pspec_dict['pC_fold']
        bootstrap_error = pspec_dict['pC_fold_up']
    except: # find from pspec_2d_to_1d results
        fold_factor = pspec_dict['k']**3/(2*np.pi**2)
        delta2_value = pspec_dict['pCv_fold']*fold_factor
        bootstrap_error = pspec_dict['pCv_fold_err']*2*fold_factor

    good_k = pspec_dict['k'] < .1
    upperlims = delta2_value + bootstrap_error
    upperlims = np.ma.masked_array(upperlims)
    upperlims.mask = good_k
    min_k = upperlims.argmin(fill_value=1e20)

    abs_lims = np.abs(delta2_value) + bootstrap_error
    abs_lims = np.ma.masked_array(abs_lims)
    abs_lims.mask = good_k
    min_k_abs = abs_lims.argmin(fill_value=1e20)

    print "\tz={0:5.2f} at k={1:.2f}: ({2:7.3f} mk)^2 ".format(redshift, pspec_dict['k'][min_k], np.sqrt(np.abs(upperlims))[min_k])
    print "\t{0:<7s}Abs k={1:.2f}: ({2:7.3f} mk)^2".format('', pspec_dict['k'][min_k_abs], np.sqrt(abs_lims)[min_k_abs])

for index, filename in enumerate(files):
    redshift = redshifts[index]
    pspec_dict = np.load(file)
    inttime = pspec_dict['frf_inttime']
    cnt = pspec_dict['cnt_eff']
    nbls_g = pspec_dict['nbls_g']
    nlsts = len(pspec_dict['lsts']) * pspec_dict['inttime']
    nlsts /= pspec_dict['frf_inttime']
    nlsts_g = pspec_dict['nlsts_g']
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
    S.Nlstbins = nlsts_g
    S.calc()
    #print 'Beggining Table Data: '
    #print '\n'
    for k in args.kmags:
        k_ind = np.argmin(np.abs(k - f['k']))
        try: # find from signal loss results
            delta2_value = pspec_dict['pC_fold'][k_ind]
            bootstrap_error = pspec_dict['pC_fold_up'][k_ind]
        except: # find from pspec_2d_to_1d results
            fold_factor = pspec_dict['k'][k_ind]**3/(2*np.pi**2)
            delta2_value = pspec_dict['pCv_fold'][k_ind]*fold_factor
            bootstrap_error = pspec_dict['pCv_fold_err'][k_ind]*2*fold_factor

        # Formatting here for a LaTEX table
        # Need to print:
        # k, redshift, delta^2, theoretical error, theory sigma, (Det/Ulim)
        #                       bootstrap error, bootstrap sigma, (Det/Ulim)
        print k, '&',
        print r(redshift), '&',
        print r(delta2_value), '&',
        print '$\pm$ '+  r(S.Delta2_N(k)*2), '&',
        print r(delta2_value/(S.Delta2_N(k))), '&',
        if delta2_value - (S.Delta2_N(k)*2) > 0:
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
