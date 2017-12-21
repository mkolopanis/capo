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
                    help=('Number of Bootstraps (averages) default=100 '
                          'if using bootstrap versions 1 or 2.'))
parser.add_argument('--Neff_lst', default=None,
                    help=('Number of effective LSTs '
                          '(if using bootstrap versions 1 or 2). '
                          'If none (default), '
                          'it is calculated using Nlstbins and t_eff.'))
parser.add_argument('--avg_func', default='np.mean', type=str,
                    help='Average function to use (default = np.mean).')
parser.add_argument('-v', '--version', default=4, type=int,
                    help=('Version of bootstrapping method. '
                          'Existing versions are 1,2, or 4.'))
parser.add_argument('--NGPS_LST', default=1, type=int,
                    help='Number of total LST groups to average.')
args = parser.parse_args()

np.random.seed(0)

pspecs = read_bootstraps_dcj(args.files)

for key in pspecs:
    if args.version == 4:
        if key[0] == 'p':
            pspecs[key] = pspecs[key][0] # turn 4-dimensional array back to 3-dim (this is okay since there was only one bootsigloss file read in)
    if args.version == 1 or args.version == 2:
        if key[0] == 'p':
            pspecs[key] = np.average(pspecs[key],axis=1) # average along cross-multiplication axis (bl groups)

if args.Neff_lst is None:
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

# decimate data in time
# for key in pspecs:
#    if key[0] == 'p':
#        shape = pspecs[key].shape
#        indices = np.linspace(0, shape[2], Neff_lst,
#                              dtype='int', endpoint=False)
#        pspecs[key] = pspecs[key][:,:,indices] # down-selecting in time

# average LSTs within groups if given NGPS_LST (default is no averaging)
if args.NGPS_LST == 0 and args.version == 4: NGPS_LST = Neff_lst
elif args.NGPS_LST != 0 and args.version == 4: NGPS_LST = args.NGPS_LST
else: NGPS_LST = 1 # 1 LST group if not version 4

if args.version == 4: # average in LST if version 4
    for pspec in pspecs:
        if pspec[0] == 'p':
            temp_pspec = []
            indices = np.linspace(0, pspecs[pspec].shape[2], NGPS_LST + 1,
                                  endpoint=True, dtype='int')
            for i, index in enumerate(range(len(indices) - 1)):
                temp = pspecs[pspec][:, :, indices[i]:indices[i + 1]]
                temp_pspec.append(np.ma.average(temp, axis=2))
            avg_pspec = np.ma.array(temp_pspec).swapaxes(0, 2).swapaxes(0, 1)
            pspecs[pspec] = avg_pspec

# compute Pk vs kpl vs bootstraps
print "   Bootstrapping..."
pk_pspecs, vals = average_bootstraps(pspecs, Nt_eff=Neff_lst,
                                     Nboots=args.nboots,
                                     avg_func=eval(args.avg_func),
                                     version=args.version)

# Over-write PS points from "pspec_noboot.npz" file
pspec_noboot = np.load('/'.join(args.files[0].split('/')[:-1])+'/pspec_noboot.npz')
for key in pk_pspecs:
    if len(key) == 3 and key[0] == 'p': # "pIn, pCv, etc."
        points = np.mean(pspec_noboot[key],axis=0) # average along cross-multiplication axis
        points = np.mean(points,axis=1) # average along time
        pk_pspecs[key] = points # function of k only

# Save 
outname = 'pspec_2d_to_1d.npz'
print '   Saving', outname  # save all values used in bootstrapping
np.savez(args.output + outname, **vals)

# Compute |k|
bl_length = np.linalg.norm(pspecs['uvw'])
wavelength = cosmo_units.c / (pspecs['freq'] * 1e9)
ubl = bl_length / wavelength
kperp = dk_du(pspecs['freq']) * ubl
print "   freq = ", pspecs['freq']
print "   kperp = ", kperp
pk_pspecs['k'] = np.sqrt(kperp**2 + pk_pspecs['kpl_fold']**2)
pk_pspecs['kperp'] = np.ma.masked_invalid(kperp)
pk_pspecs['cmd'] = pk_pspecs['cmd'].item() + ' \n ' + ' '.join(sys.argv)
if NGPS_LST > 0: pk_pspecs['nlsts_g'] = Neff_lst / NGPS_LST  # number of lsts in one group
else: pk_pspecs['nlsts_g'] = Neff_lst
if args.version == 4: pk_pspecs['nPS'] = pspecs['pCv'].shape[0] * pspecs['pCv'].shape[2]
else: pk_pspecs['nPS'] = pspecs['pCv'].shape[0]

# Important numbers
print "   Total number of bls = ", pk_pspecs['nbls']
print "      number of bl groups = ", pk_pspecs['ngps']
print "      nbls in a group = ", pk_pspecs['nbls'] / pk_pspecs['ngps']
print "   Total number of lsts = ", Neff_lst
print "      number of lst groups = ", NGPS_LST
print "      nlsts in a group  = ", pk_pspecs['nlsts_g']

# Scale for error on error
if pk_pspecs['nPS'] != 1:
    scaling = 1. + (1. / np.sqrt(2 * (pk_pspecs['nPS'] - 1)))
else:
    scaling = 1
    print ('   !!! Warning: Scaling blows up '
           'since there is only one independent sample. '
           'Not applying correction factor !!!')
print '   We have', pk_pspecs['nPS'], 'independent samples...'
print '   ... Therefore, we correct for error on error using factor =', scaling

for key in pk_pspecs:
    if key[0] == 'p' and key[-3:] == 'err':
        pk_pspecs[key] = pk_pspecs[key] * scaling

# Save values
for key in pk_pspecs.keys():
    if isinstance(pk_pspecs[key], np.ma.MaskedArray):
        pk_pspecs[key].fill_value = 0  # fills invalid values with 0's
        pk_pspecs[key] = pk_pspecs[key].filled()  # returns corrected array
filename = args.output + 'pspec_pk_k3pk.npz'
print '   Saving', filename
np.savez(filename, **pk_pspecs)
