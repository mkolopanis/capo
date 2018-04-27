#! /usr/bin/env python

"""Transform 2-Dimensional Power Spectrum to 1-D.

Takes an input a directory full of 2-D power spectrum bootstrap files.
Outpust a power spectrum file with folded and unfoled k//.s
"""

import sys
from glob import glob
import argparse
from capo.eor_results import read_bootstraps_dcj, average_bootstraps, split_stack_kpl
from capo.pspec import dk_du
from capo import cosmo_units
import numpy as np
from capo.cosmo_units import f212z
from capo import sensitivity


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
                          'If none (default), '
                          'it is calculated using Nlstbins and t_eff.'))
parser.add_argument('--avg_func', default='np.mean', type=str,
                    help='Average function to use (default = np.mean).')
args = parser.parse_args()

np.random.seed(0)

pspecs = read_bootstraps_dcj(args.files)

for key in pspecs:
    if key[0] == 'p':
        pspecs[key] = np.average(pspecs[key],axis=1) # average along cross-multiplication axis (bl groups)

if args.Neff_lst is None:
    Nlstbins = np.shape(pspecs['pCr'])[-1]
    # get the number of lst integrations in the dataset
    t_eff = pspecs['frf_inttime'] / pspecs['inttime']
    Neff_lst = np.ceil(Nlstbins / t_eff)
else:
    Neff_lst = int(float(args.Neff_lst))

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
print "   Bootstrapping..."
pk_pspecs, vals = average_bootstraps(pspecs, Nt_eff=Neff_lst,
                                     Nboots=args.nboots,
                                     avg_func=eval(args.avg_func),
                                     version=2)

# Over-write PS points from "pspec_noboot.npz" file
pspec_noboot = read_bootstraps_dcj(['/'.join(args.files[0].split('/')[:-1])+'/pspec_noboot.npz'])
pspec_noboot['pCr-pCv'] = pspec_noboot['pCr'] - pspec_noboot['pCv']  # subtracted
pspec_noboot['pCs-pCn'] = pspec_noboot['pCs'] - pspec_noboot['pCn']
pspec_noboot['pIr-pIv'] = pspec_noboot['pIr'] - pspec_noboot['pIv']
pspec_noboot['pIs-pIn'] = pspec_noboot['pIs'] - pspec_noboot['pIn']

for key in pk_pspecs: # loop over all keys that will get saved in pspec_pk_k3pk file
    overwrite = False
    if len(key) == 3 and key[0] == 'p':  # "pIn", "pCv", etc.
        points = pspec_noboot[key][0, :, :, :]  # only one file so take that one
        overwrite = True
    if key[-4:] == "fold" and key != "kpl_fold":  # "pIn_fold", "pCv_fold", etc.
        _, points = split_stack_kpl(pspec_noboot[key[:-5]][0, :, :, :],
                                    pspec_noboot['kpl'])
        overwrite = True
    if overwrite is True:
        points = np.mean(points, axis=0) # average along cross-multiplication axis
        points = np.mean(points, axis=1)  # average along time
        pk_pspecs[key] = points  # function of k only

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
pk_pspecs['nPS'] = pspecs['pCv'].shape[0]

# Important numbers
print "   Total number of bls = ", pk_pspecs['nbls']
print "      number of bl groups = ", pk_pspecs['ngps']
print "      nbls in a group = ", pk_pspecs['nbls_g']
print "   Total number of lsts = ", Neff_lst

""" # Commented out because not sure if this expression is correct
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
"""

print('Computing 21cmSense_Calc Noise')
print('Using these parameters:')
redshift = f212z(pk_pspecs['freq'] * 1e9)
inttime = pk_pspecs['frf_inttime']
cnt = pk_pspecs['cnt_eff']
nbls_g = pk_pspecs['nbls_g']
nlsts = Neff_lst
if pk_pspecs['frf_inttime'] == pk_pspecs['inttime']:
    omega_eff = .74**2/.32  # for capo analytical; from T1 of Parsons FRF paper
else:
    omega_eff = .74**2/.24
print 'Redshift:', redshift
print '\tT_int:', inttime
print '\tNbls:', pk_pspecs['nbls']
print '\tNgps:', pk_pspecs['ngps']
print '\tNdays:', cnt
print '\tNlstbins:', nlsts
S = sensitivity.Sense()
f = pk_pspecs['freq']
S.z = f212z(f*1e9)

#   Tsys
#S.Tsys = 551e3  #set to match 21cmsense exactly
#S.Tsys = 505e3 #Ali et al, at 164MHz
if 'Trcvr' in pk_pspecs.keys():
    S.Tsys = (pk_pspecs['Trcvr'] + 180.*(f/.180)**-2.55)*1e3
else:
    S.Tsys = ( 144. + 180.*(f/.180)**-2.55)*1e3
    # set to match noise realization
    # calcuation made from correcting Ali 2015 Tsys estimate

print "Tsys = ",S.Tsys

S.t_int = inttime
S.Ndays = cnt  #effective number of days
S.Npols = 2
try: S.Nseps = pk_pspecs['nseps']
except: S.Nseps = 1
print "Nseps = ",S.Nseps
S.Nblgroups = pk_pspecs['ngps']
# use the FRF weighted beams listed in T1 of Parsons etal beam sculpting paper
S.Omega_eff = omega_eff
k = pk_pspecs['k']
S.Nbls = pk_pspecs['nbls']
S.Nlstbins = nlsts
S.calc()
print "capo.sensitivity Pk_noise = ", S.P_N
pk_pspecs['theory_noise'] = np.repeat(S.P_N, len(pk_pspecs['kpl']))
pk_pspecs['theory_noise_delta2'] = S.Delta2_N(k)


# Save values
for key in pk_pspecs.keys():
    if isinstance(pk_pspecs[key], np.ma.MaskedArray):
        pk_pspecs[key].fill_value = 0  # fills invalid values with 0's
        pk_pspecs[key] = pk_pspecs[key].filled()  # returns corrected array
filename = args.output + 'pspec_pk_k3pk.npz'
print '   Saving', filename
np.savez(filename, **pk_pspecs)
