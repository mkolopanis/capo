#!/usr/bin/env python
"""Create simple 2 sigma upper limit plot.

Takes outputs from pspec_final_???.py and creates 2 sigma errorbar plots.
"""
import numpy as np
import sys
import os
import argparse
import py21cmsense as py21cm
import capo
from capo.cosmo_units import f212z, DM, E, ckm, Ho
import matplotlib.pyplot as plt
from matplotlib import gridspec
from capo import sensitivity

parser = argparse.ArgumentParser(
    description=('Calculate power spectra for a run '
                 'from pspec_oqe_2d.py'))
parser.add_argument('files', metavar='<FILE>', type=str, nargs='+',
                    help='List of pspec files to plot')
parser.add_argument('--plot', action='store_true',
                    help='output plot before saving')
parser.add_argument('--analytical', action='store_true',
                    help='Plot analytical noise approximation')
parser.add_argument('--noisefiles', type=str, nargs='*',
                    help='supply 21cmSense files to plot sensitivity')
parser.add_argument('--outfile', type=str, default='pspec_k3pk',
                    help='give custom output file name')
parser.add_argument('--Trcvr', type=float, default=144,
                    help='Receiver Temperature in Kelvin (defualt 144)')
parser.add_argument('--models', type=str,
                    help=('a string that globs a list of 21cmfast pspec'
                          ' output files'))
args = parser.parse_args()


zs = []
for filename in args.files:
    F = np.load(filename)
    freq = F['freq']
    zs.append(f212z(freq * 1e9))

# reverse order to make high redshift of left of plot
zs = np.unique(zs)[::-1]
Nzs = len(zs)
# define some color format sets for the plotter
marker = 'o'
figsize = (5 * (1 + Nzs)/2., 6)

# Create figure and prep subplot sizes for Delta^2
fig = plt.figure(figsize=figsize)
gs = gridspec.GridSpec(1, Nzs)
gs.update(hspace=0.0, wspace=0.2)
ax1 = [plt.subplot(gs[:, i]) for i in range(Nzs)]
plt.subplots_adjust(left=0.2)

# Create figure and prep subplots for P(k)
fig2 = plt.figure(figsize=figsize)
gs = gridspec.GridSpec(1, Nzs)
gs.update(hspace=0.0, wspace=0.2)
ax2 = [plt.subplot(gs[:, i]) for i in range(Nzs)]
plt.subplots_adjust(left=0.2)

# Create figure and prep subplot sizes for Delta^2 Noise
fig3 = plt.figure(figsize=figsize)
gs = gridspec.GridSpec(1, Nzs)
gs.update(hspace=0.0, wspace=0.2)
ax3 = [plt.subplot(gs[:, i]) for i in range(Nzs)]
plt.subplots_adjust(left=0.2)

# Create figure and prep subplot sizes for P(k) Noise
fig4 = plt.figure(figsize=figsize)
gs = gridspec.GridSpec(1, Nzs)
gs.update(hspace=0.0, wspace=0.2)
ax4 = [plt.subplot(gs[:, i]) for i in range(Nzs)]
plt.subplots_adjust(left=0.2)

# group in delta2 and pk for easy looping later
delta2_list = [ax1, ax3]
pk_list = [ax2, ax4]


if args.noisefiles:
    noise_freqs, noise_ks, noises = py21cm.load_noise_files(
        args.noisefiles, polyfit_deg=3)

    freqs = capo.pspec.z2f(zs) * 1e3
    noise_ind = np.array([np.argmin(abs(np.array(noise_freqs) - fq))
                          for fq in freqs])

    noise_freqs = np.take(noise_freqs, noise_ind).tolist()
    noise_ks = np.take(noise_ks, noise_ind, axis=0).tolist()
    noises = np.take(noises, noise_ind, axis=0).tolist()

    # Plot the Noise files on the plot
    for gs_ind in xrange(Nzs):

        d2_n = 2*np.array(noises[gs_ind])
        pk_n = 2*np.pi**2/(np.array(noise_ks[gs_ind])**3)*d2_n
        for ax in pk_list:
            ax[gs_ind].plot(noise_ks[gs_ind], pk_n, '-', color='c')
            ax[gs_ind].plot(-np.array(noise_ks[gs_ind]), pk_n, '-', color='c', label='21cmSense')

        for ax in delta2_list:
            ax[gs_ind].plot(noise_ks[gs_ind], d2_n, 'c-', label='21cmSense')

EoR_MODEL_Delta2 = None
if args.models is not None:
    print "Loading 21cmFAST models in {0}".format(args.models.split('/')[-2])
    from twentyonecmfast_tools import load_andre_models, all_and
    parm_array, model_ks, delta2_array, delta2_err_array = load_andre_models(args.models)
    ZS = np.tile(parm_array[:, 0], (model_ks.shape[1], 1)).T
    EoR_MODEL_Delta2 = interp2d(ZS, model_ks, delta2_array, kind='cubic')


k_max = 0
k_par_max = 0
for filename in args.files:
    pspec_dict = np.load(filename)
    if np.max(pspec_dict['k']) > k_max:
        k_max = np.max(pspec_dict['k'])

    if np.max(np.abs(pspec_dict['kpl'])) > k_par_max:
        k_par_max = np.max(np.abs(pspec_dict['kpl']))

    redshift = f212z(pspec_dict['freq'] * 1e9)
    gs_ind = int(np.where(zs == redshift)[0].item())
    # Find the horizon and do not plot inside it
    z = zs[gs_ind]
    cosmo_factor = DM(z) * E(z)
    cosmo_factor /= (1+z) * (ckm/Ho)
    horizon = 90 * cosmo_factor * np.pi/180 * pspec_dict['kperp']
    horizon_delta2 = np.sqrt(horizon**2 + pspec_dict['kperp']**2)
    horizon /= 2 # XXX divide by 2 because of bug (kperp two times too big?)
    horizon_delta2 /= 2
    for _ax in delta2_list:
        _ax[gs_ind].axvline(horizon_delta2, color='black', linestyle='--',
                            ymin=-.1, ymax=1.1)
    for _ax in pk_list:
        _ax[gs_ind].axvline(-horizon, color='black', linestyle='--',
                            ymin=-.1, ymax=1.1)
        _ax[gs_ind].axvline(horizon, linestyle='--', color='black',
                            ymin=-.1, ymax=1.1)
    k_vals = np.ma.masked_array(pspec_dict['k'])
    kpl_vals = np.ma.masked_array(pspec_dict['kpl'])
    kpl_vals.mask = np.abs(kpl_vals) < horizon
    # kpl_vals.mask = np.abs(kpl_vals) <= horizon
    k_vals = np.ma.masked_less(pspec_dict['k'], horizon_delta2)

    if args.analytical:
        if ('theory_noise' and 'theory_noise_delta2') in pspec_dict.keys():
            print 'Using saved theory_noise and theory_noise_delta2 in pspec filez'
            k3pk_noise = 2*pspec_dict['theory_noise_delta2']
            pk_noise = 2*pspec_dict['theory_noise']

            ax1[gs_ind].plot(k_vals, 2*pspec_dict['theory_noise_delta2'],
                             'g-', label='Analytical 2-sigma')
            ax2[gs_ind].plot(kpl_vals, 2*pspec_dict['theory_noise'], 'g-',
                             label='Analytical 2-sigma')
            ax3[gs_ind].plot(k_vals, 2*pspec_dict['theory_noise_delta2'],
                             'g-', label='Analytical 2-sigma')
            ax4[gs_ind].plot(kpl_vals, 2*pspec_dict['theory_noise'], 'g-',
                             label='Analytical 2-sigma')
        else:
            print 'Using these parameters for the Analytic Model:'

            inttime = pspec_dict['frf_inttime']
            cnt = pspec_dict['cnt_eff']
            nbls_g = pspec_dict['nbls_g']
            nlsts = len(pspec_dict['lsts']) * pspec_dict['inttime']
            nlsts /= pspec_dict['frf_inttime']
            nlsts_g = pspec_dict['nlsts']
            if pspec_dict['frf_inttime'] == pspec_dict['inttime']:
                omega_eff = .74**2/.32  # for capo analytical; from T1 of Parsons FRF paper
            else:
                omega_eff = .74**2/.24
            print 'Redshift:', redshift
            print '\tT_int:', inttime
            print '\tNbls:', pspec_dict['nbls']
            print '\tNgps:', pspec_dict['ngps']
            print '\tNdays:', cnt
            print '\tNlstbins:', nlsts
            S = sensitivity.Sense()
            f = pspec_dict['freq']
            S.z = capo.pspec.f2z(f)

            #   Tsys
            # S.Tsys = 551e3  #set to match 21cmsense exactly
            # S.Tsys = 505e3 #Ali et al, at 164MHz
            S.Tsys = (args.Trcvr + 180.*(f/.180)**-2.55)*1e3  # default set to match noise realization
            print "Tsys = ",S.Tsys

            S.t_int = inttime
            S.Ndays = cnt  # effective number of days
            S.Npols = 2
            try: S.Nseps = pspec_dict['nseps']
            except: S.Nseps = 1
            print "Nseps = ",S.Nseps
            S.Nblgroups = pspec_dict['ngps']
            S.Omega_eff = omega_eff  # use the FRF weighted beams listed in T1 of Parsons etal beam sculpting paper
            k = k_vals
            S.Nbls = pspec_dict['nbls']
            S.Nlstbins = nlsts_g
            S.calc()
            print "capo.sensitivity Pk_noise = ",S.P_N
            k3pk_noise = S.Delta2_N(k)*2
            pk_noise = 2 * S.P_N
            ax5[gs_ind].axhline(S.P_N*2, color='g', marker='_', label='Analytical 2-sigma')
            ax6[gs_ind].axhline(S.P_N*2, color='g', marker='_', label='Analytical 2-sigma')
            ax1[gs_ind].plot(k, S.Delta2_N(k)*2, 'g-', label='Analytical 2-sigma')
            ax2[gs_ind].axhline(S.P_N*2, color='g', marker='_', label='Analytical 2-sigma')
            ax3[gs_ind].plot(k, S.Delta2_N(k)*2, 'g-', label='Analytical 2-sigma')
            ax4[gs_ind].axhline(S.P_N*2, color='g', marker='_', label='Analytical 2-sigma')

    # This should only work with pspec_final_sigloss_dist_v8 outputs (maybe v7)
    fold_factor = pspec_dict['k']**3/(2*np.pi**2)

    # Load upper limits derived from Bayes' Theorem inversion of
    # probability curve
    pCv_limit = pspec_dict['pCv_err']*2
    pCv_fold_limit = pspec_dict['pCv_fold_err']*2*fold_factor
    pIv_limit = pspec_dict['pIv_err']*2
    pIv_fold_limit = pspec_dict['pIv_fold_err']*2*fold_factor
    pIn_limit = pspec_dict['pIn_err']*2
    pIn_fold_limit = pspec_dict['pIn_fold_err']*2*fold_factor
    pCn_limit = pspec_dict['pCn_err']*2
    pCn_fold_limit = pspec_dict['pCn_fold_err']*2*fold_factor

    # Load data points from FFT and C_eff pspecs
    pCv = pspec_dict['pCv_old']
    pIv = pspec_dict['pIv_old']
    pCv_fold = pspec_dict['pCv_fold_old']*fold_factor
    pIv_fold = pspec_dict['pIv_fold_old']*fold_factor
    pCn = pspec_dict['pCn_old']
    pIn = pspec_dict['pIn_old']

    # Load errorbars derived from bootstrap variance
    pCn_fold = pspec_dict['pCn_fold_old']*fold_factor
    pIn_fold = pspec_dict['pIn_fold_old']*fold_factor
    pCv_error = pspec_dict['pCv_err_old']*2
    pIv_error = pspec_dict['pIv_err_old']*2
    pCv_fold_error = pspec_dict['pCv_fold_err_old']*2*fold_factor
    pIv_fold_error = pspec_dict['pIv_fold_err_old']*2*fold_factor
    pCn_error = pspec_dict['pCn_err_old']*2
    pIn_error = pspec_dict['pIn_err_old']*2
    pCn_fold_error = pspec_dict['pCn_fold_err_old']*2*fold_factor
    pIn_fold_error = pspec_dict['pIn_fold_err_old']*2*fold_factor
    prob = pspec_dict['prob']*100

    pos_ind = np.where(pIv >= 0)[0]
    pos_ind_noise = np.where(pIn >= 0)[0]
    pos_ind_fold = np.where(pIv_fold >= 0)[0]
    pos_ind_noise_fold = np.where(pIn_fold >= 0)[0]
    neg_ind = np.where(pIv < 0)[0]
    neg_ind_noise = np.where(pIn < 0)[0]
    neg_ind_fold = np.where(pIv_fold < 0)[0]
    neg_ind_noise_fold = np.where(pIn_fold < 0)[0]

    #  Plot Data delta^2 values
    ax1[gs_ind].errorbar(k_vals[pos_ind_fold],
                         pIv_fold[pos_ind_fold],
                         pIv_fold_error[pos_ind_fold],
                         marker=marker, color='black', linestyle='',
                         label='pI {0:02d}%'.format(int(prob)))
    ax1[gs_ind].errorbar(k_vals[neg_ind_fold],
                         -pIv_fold[neg_ind_fold],
                         pIv_fold_error[neg_ind_fold], linestyle='',
                         color='black', marker=marker, alpha=.5)

    ax1[gs_ind].fill_between(k_vals,
                             np.abs(pIv_fold) - k3pk_noise,
                             np.abs(pIv_fold) + k3pk_noise,
                             color='grey', alpha=.5, step='mid')

    ax1[gs_ind].plot(k_vals.data[k_vals > horizon/1.1],
                     pCv_fold_limit[k_vals > horizon/1.1],'-', color='red',
                     label='C$_{eff}$ Limits', drawstyle='steps-mid')
    ax1[gs_ind].plot(k_vals.data[k_vals > horizon/1.1],
                     pIv_fold_limit[k_vals > horizon/1.1], '-', color='blue',
                     label='Uniform Limits', drawstyle='steps-mid')

    # Plot Data P(K) values
    ax2[gs_ind].errorbar(kpl_vals[pos_ind],
                         pIv[pos_ind],
                         pIv_error[pos_ind],
                         marker=marker, color='black', linestyle='',
                         label='pI {0:02d}%'.format(int(prob)))
    ax2[gs_ind].errorbar(kpl_vals[neg_ind],
                         -pIv[neg_ind],
                         pIv_error[neg_ind], linestyle='',
                         color='black', marker=marker, alpha=.5)

    ax2[gs_ind].fill_between(kpl_vals,
                             np.abs(pIv) - pk_noise,
                             np.abs(pIv) + pk_noise,
                             color='grey', alpha=.5, step='mid')

    #  This splitting of the kpl over +/- makes the plots look nicer
    #  The steps extend over the data points
    for inds in [kpl_vals > horizon/1.1, kpl_vals < -horizon/1.1]:
        ax2[gs_ind].plot(kpl_vals.data[inds],
                         pCv_limit[inds], '-', color='red',
                         label='C$_{eff}$ Limits', drawstyle='steps-mid')
        ax2[gs_ind].plot(kpl_vals[inds], pIv_limit[inds], '-', color='blue',
                         label='Uniform Limits', drawstyle='steps-mid')

    #  Plot noise delta^2 values
    ax3[gs_ind].errorbar(k_vals[pos_ind_noise_fold],
                         pIn_fold[pos_ind_noise_fold],
                         pIn_fold_error[pos_ind_noise_fold],
                         marker=marker, color='black', linestyle='',
                         label='pI {0:02d}%'.format(int(prob)))
    ax3[gs_ind].errorbar(k_vals[neg_ind_noise_fold],
                         -pIn_fold[neg_ind_noise_fold],
                         pIn_fold_error[neg_ind_noise_fold], linestyle='',
                         color='black', marker=marker, alpha=.5)

    ax3[gs_ind].fill_between(k_vals,
                             np.abs(pIn_fold) - k3pk_noise,
                             np.abs(pIn_fold) + k3pk_noise,
                             color='grey', alpha=.5, step='mid')

    ax3[gs_ind].plot(k_vals.data[k_vals > horizon/1.1],
                     pCn_fold_limit[k_vals > horizon/1.1], '-', color='red',
                     label='C$_{eff}$ Limits', drawstyle='steps-mid')
    ax3[gs_ind].plot(k_vals.data[k_vals > horizon/1.1],
                     pIn_fold_limit[k_vals > horizon/1.1], '-', color='blue',
                     label='Uniform Limits', drawstyle='steps-mid')

    #  Plot noise P(K) values
    ax4[gs_ind].errorbar(kpl_vals[pos_ind_noise],
                         pIn_error[pos_ind_noise],
                         pIn[pos_ind_noise],
                         marker=marker, color='black', linestyle='',
                         label='pI {0:02d}%'.format(int(prob)))
    ax4[gs_ind].errorbar(kpl_vals[neg_ind_noise],
                         -pIn[neg_ind_noise],
                         pIn_error[neg_ind_noise], linestyle='',
                         color='black', marker=marker, alpha=.5)

    ax4[gs_ind].fill_between(kpl_vals,
                             np.abs(pIn) - pk_noise,
                             np.abs(pIn) + pk_noise,
                             color='grey', alpha=.5, step='mid')

    #  This splitting of the kpl over +/- makes the plots look nicer
    #  The steps extend over the data points
    for inds in [kpl_vals > horizon/1.1, kpl_vals < -horizon/1.1]:
        ax4[gs_ind].plot(kpl_vals.data[inds],
                         pCn_limit[inds], '-', color='red',
                         label='C$_{eff}$ Limits', drawstyle='steps-mid')
        ax4[gs_ind].plot(kpl_vals.data[inds],
                         pIn_limit[inds], '-', color='blue',
                         label='Uniform Limits', drawstyle='steps-mid')

    if EoR_MODEL_Delta2 is not None:
            for _ax in delta2_list:
                _ax[gs_ind].plot(k_vals,
                                 EoR_MODEL_Delta2(redshift, pspec_dict['k']),
                                 'k', label='Fiducial 21cmFAST model')
            for _ax in pk_list:
                _ax[gs_ind].plot(pspec_dict['kpl_fold'],
                                 EoR_MODEL_Delta2(redshift, pspec_dict['k']) *
                                 (2*np.pi**2) /
                                 np.array(pspec_dict['k'].reshape(-1, 1))**3,
                                 'k', label='Fiducial 21cmFAST model')
                _ax[gs_ind].plot(-pspec_dict['kpl_fold'],
                                 EoR_MODEL_Delta2(redshift, pspec_dict['k']) *
                                 (2*np.pi**2) /
                                 np.array(pspec_dict['k'].reshape(-1, 1))**3, 'k')

# set up some parameters to make the figures pretty
# Loop over axes in each pk and delta 2 in case axes are added or deleted later
for gs_ind in xrange(Nzs):
    # only set ylabel for first plot
    if gs_ind == 0:
        for _ax in delta2_list:
            _ax[gs_ind].set_ylabel('$\\frac{k^{3}}{2\pi^{2}}$ $P(k)$ [mK$^{2}$]', fontsize=16)
        for _ax in pk_list:
            _ax[gs_ind].set_ylabel('$P(k)$ [mK$^{2}(h^{-1}$ Mpc)$^{3}$]', fontsize=16)

    for _ax in delta2_list:
        _ax[gs_ind].set_title('z = {0:.2f}'.format(zs[gs_ind]), fontsize=14)
        _ax[gs_ind].set_yscale('log', nonposy='clip')
        _ax[gs_ind].set_xlabel('$k$ [$h$ Mpc$^{-1}$]', fontsize=16)
        _ax[gs_ind].set_xlim(0, k_max * 1.01)
        _ax[gs_ind].get_shared_y_axes().join(ax1[0], ax1[gs_ind])
        _ax[gs_ind].grid(True)
        _ax[gs_ind].tick_params(axis='both', which='major', labelsize=14)

    for _ax in pk_list:
        _ax[gs_ind].set_yscale('log', nonposy='clip')
        _ax[gs_ind].set_title('z = {0:.2f}'.format(zs[gs_ind]), fontsize=14)
        _ax[gs_ind].set_xlabel('$k_{\\parallel}$ [$h$ Mpc$^{-1}$]', fontsize=16)
        _ax[gs_ind].set_xlim(-1.01 * k_par_max, k_par_max * 1.01)
        _ax[gs_ind].get_shared_y_axes().join(ax2[0], ax2[gs_ind])
        _ax[gs_ind].grid(True)
        _ax[gs_ind].tick_params(axis='both', which='major', labelsize=14)

    # if multi redshift, make shared axes invisible
    if gs_ind > 0:
        for _ax in pk_list + delta2_list:
            plt.setp(_ax[gs_ind].get_yticklabels(), visible=False)
# Check for maximum value along all sub-plots
# use this value to set the ylim for all delta squared plots
delta2_list = [ax1, ax3]
max_d2 = -np.Inf
for ax_set in delta2_list:
    for ax in ax_set:
        for lines in ax.get_lines():
            # must check if there is data in the line
            _data_array = lines.get_ydata()
            if np.size(_data_array):
                max_line_val = np.max(np.abs(_data_array))
                if max_line_val > max_d2:
                    max_d2 = max_line_val

# use this value for all P(k) plots
pk_list = [ax2, ax4]
max_pk = -np.Inf
for ax_set in pk_list:
    for ax in ax_set:
        for lines in ax.get_lines():
            _data_array = lines.get_ydata()
            if np.size(_data_array):
                max_line_val = np.max(np.abs(_data_array))
                if max_line_val > max_pk:
                    max_pk = max_line_val

# round up the power of 10 and add 1
# this will set the ymax value to 10 * the highest value rounded up
# across all plots of each type
max_val_d2 = np.ceil(np.log10(max_d2))
ymax_d2 = np.power(10., max_val_d2)

max_val_pk = np.ceil(np.log10(max_pk))
ymax_pk = np.power(10., max_val_pk)

if ymax_d2 > ymax_pk:  # use highest ymax for both delta^2 and p(k) plots, so that they're both the same
    ymax_pk = ymax_d2.copy()
    max_val_pk = max_val_d2.copy()
if ymax_pk > ymax_d2:
    ymax_d2 = ymax_pk.copy()
    max_val_d2 = max_val_pk.copy()

# ymax_d2 = 4e10 # XXX
# max_val_d2 = np.log10(ymax_d2) # XXX

ax1[0].set_ylim([1e-1, ymax_d2])
ax1[0].set_xlim([0.0, 0.6])
ax2[0].set_ylim([1e-1, ymax_pk])
ax2[0].set_xlim([-0.6, 0.6])
ax3[0].set_ylim([1e-1, ymax_d2])
ax3[0].set_xlim([0.0, 0.6])
ax4[0].set_ylim([1e-1, ymax_pk])
ax4[0].set_xlim([-0.6, 0.6])

# Some matplotlib versions will only give every other power of ten
# The next few lines set the log yticks to be every power of ten
ytick_d2 = np.power(10., np.arange(-1, max_val_d2+1))
ytick_pk = np.power(10., np.arange(-1, max_val_pk+1))

for axes in delta2_list:
    for ax in axes:
        ax.set_yticks(ytick_d2)
        # set to force matplotlib to make more x-ticks
        # some matplotlibs only plot every other .1 kvalues
        ax.locator_params(nbins=10, axis='x')

for axes in pk_list:
    for ax in axes:
        ax.set_yticks(ytick_pk)
        ax.locator_params(nbins=10, axis='x')

# handles, labels = ax1[-1].get_legend_handles_labels()
# ax1[-1].legend(handles, labels, loc='lower right', numpoints=1)
# handles, labels = ax2[-1].get_legend_handles_labels()
# ax2[-1].legend(handles, labels, loc='lower right', numpoints=1)
# handles, labels = ax3[-1].get_legend_handles_labels()
# ax3[-1].legend(handles, labels, loc='lower right', numpoints=1)
# handles, labels = ax4[-1].get_legend_handles_labels()
# ax4[-1].legend(handles, labels, loc='lower right', numpoints=1)

fig.savefig(args.outfile+'.png', format='png')
fig2.savefig(args.outfile+'_pk.png', format='png')
fig3.savefig(args.outfile+'_noise.png', format='png')
fig4.savefig(args.outfile+'_pk_noise.png', format='png')

if args.plot:
    plt.show()
