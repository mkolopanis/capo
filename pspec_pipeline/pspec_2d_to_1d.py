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

parser = argparse.ArgumentParser(
    description=('Calculate power spectra for a run '
                 'from pspec_oqe_2d.py'))
parser.add_argument('files', metavar='<FILE>', type=str, nargs='+',
                    help='List of files to average')
parser.add_argument('--output', type=str, default='./',
                    help='Specifically specify out directory.')
parser.add_argument('--nboots', type=int, default=100,
                    help='Number of Bootstraps (averages) default=100')
parser.add_argument('--Neff_lst', default=None,
                    help='Number of effective LSTs. If none (default), it is calculated using Nlstbins and t_eff.')
args = parser.parse_args()

np.random.seed(0)
pspecs = read_bootstraps_dcj(args.files)
if args.Neff_lst == None:
    Nlstbins = np.shape(pspecs['pCr'])[-1]
    # get the number of lst integrations in the dataset
    t_eff = pspecs['frf_inttime'] / pspecs['inttime']
    Neff_lst = np.ceil(Nlstbins / t_eff)
else:
    Neff_lst = int(float(args.Neff_lst))

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

pspecs['pCr-pCv'] = pspecs['pCr'] - pspecs['pCv']  # subtracted
pspecs['pCs-pCn'] = pspecs['pCs'] - pspecs['pCn']
pspecs['pIr-pIv'] = pspecs['pIr'] - pspecs['pIv']
pspecs['pIs-pIn'] = pspecs['pIs'] - pspecs['pIn']

# compute Pk vs kpl vs bootstraps
pk_pspecs, vals  = average_bootstraps(pspecs, Nt_eff=Neff_lst,
                                     Nboots=args.nboots, avg_func=np.median)
print 'Saving pspec_2d_to_1d.npz'  # save all values used in bootstrapping
np.savez(args.output + 'pspec_2d_to_1d.npz', **vals)

# Compute |k|
bl_length = np.linalg.norm(pspecs['uvw'])
wavelength = cosmo_units.c / (pspecs['freq'] * 1e9)
ubl = bl_length / wavelength
kperp = dk_du(pspecs['freq']) * ubl
print "freq = ", pspecs['freq']
print "kperp = ", kperp
pk_pspecs['k'] = np.sqrt(kperp**2 + pk_pspecs['kpl_fold']**2)
pk_pspecs['kperp'] = np.ma.masked_invalid(kperp)
pk_pspecs['cmd'] = pk_pspecs['cmd'].item() + ' \n ' + ' '.join(sys.argv)
for key in pk_pspecs.keys():
    if key in ['cmd']:
        continue
    if pk_pspecs[key].dtype not in [np.float]:
        continue
    try:
        pk_pspecs[key].fill_value = 0  # fills invalid values with 0's
        pk_pspecs[key] = pk_pspecs[key].filled()  # returns corrected array
    except:
        import ipdb
        ipdb.set_trace()
print 'Saving', args.output + 'pspec_pk_k3pk.npz'
np.savez(args.output + 'pspec_pk_k3pk.npz', **pk_pspecs)
