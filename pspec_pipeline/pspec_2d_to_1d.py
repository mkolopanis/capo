#! /usr/bin/env python

"""Transform 2-Dimensional Power Spectrum to 1-D.

Takes an input a directory full of 2-D power spectrum bootstrap files.
Outpust a power spectrum file with folded and unfoled k//.s
"""

import sys
from glob import glob
import argparse
from capo.eor_results import read_bootstraps_dcj, average_bootstraps
from capo.pspec import dk_du
from capo import cosmo_units
import numpy as np
from IPython import embed

parser = argparse.ArgumentParser(
            description=('Calculate power spectra for a run '
                         'from pspec_oqe_2d.py'))
parser.add_argument('files', metavar='<FILE>', type=str, nargs='+',
                    help='List of files to average')
parser.add_argument('--sub_pCv', action='store_true',
                    help='Remove pCv back from pC before averaging')
parser.add_argument('--output', type=str, default='./',
                    help='Specifically specify out directory.')
parser.add_argument('--nboots', type=int, default=100,
                    help='Number of Bootstraps (averages) default=100')
args = parser.parse_args(sys.argv[1:])

pspecs = read_bootstraps_dcj(args.files)
Nlstbins = np.shape(pspecs['pCr'])[-1]
# get the number of lst integrations in the dataset
try: t_eff = pspecs['frf_inttime']/pspecs['inttime']
except: t_eff = 3887./32. #XXX PSA128
Neff_lst = np.ceil(Nlstbins/t_eff)
# compute the effective number of LST bins
# print Neff_lst
# lets round up because this 'N' is only approximate
for key in pspecs.keys():
    if pspecs[key].dtype not in [np.float]:
        continue
    try:
        pspecs[key] = np.ma.masked_invalid(np.array(pspecs[key]))
    except:
        import ipdb
        ipdb.set_trace()
if args.sub_pCv:
    pspecs['pCr'] -= pspecs['pCv']

# compute Pk vs kpl vs bootstraps
pk_pspecs = average_bootstraps(pspecs, Nt_eff=Neff_lst,
                               Nboots=args.nboots, avg_func=np.median)

# Compute |k|
bl_length = np.linalg.norm(pspecs['uvw'])
wavelength = cosmo_units.c/(pspecs['freq']*1e9)
ubl = bl_length/wavelength
kperp = dk_du(pspecs['freq'])*ubl
print "freq = ", pspecs['freq']
print "kperp = ", kperp
pk_pspecs['k'] = np.sqrt(kperp**2 + pk_pspecs['kpl_fold']**2)
for key in pk_pspecs.keys():
    if pk_pspecs[key].dtype not in [np.float]:
        continue
    try:
        pk_pspecs[key].fill_value = 0 #fills invalid values with 0's
        pk_pspecs[key] = pk_pspecs[key].filled() #returns corrected array
    except:
        import ipdb
        ipdb.set_trace()
print args.output + 'pspec_pk_k3pk.npz'
np.savez(args.output + 'pspec_pk_k3pk.npz', **pk_pspecs)
