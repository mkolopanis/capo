#!/usr/bin/env python
"""Create simple 2 sigma upper limit plot.

Takes outputs from pspec_fina_???.py and creates 2 sigma errorbar plots.
"""

import numpy as np
import sys
import os
import argparse
import py21cmsense as py21cm
import capo
from capo.cosmo_units import f212z
import matplotlib.pyplot as plt
from matplotlib import gridspec
import itertools

parser = argparse.ArgumentParser(
    description=('Calculate power spectra for a run '
                 'from pspec_oqe_2d.py'))
parser.add_argument('files', metavar='<FILE>', type=str, nargs='+',
                    help='List of pspec files to plot')
parser.add_argument('--plot', action='store_true',
                    help='output plot before saving')
parser.add_argument('--noisefiles', type=str, nargs='*',
                    help='supply 21cmSense files to plot sensitivity')
parser.add_argument('--outfile', type=str, default='pspec_k3pk',
                    help='give custom output file name')
args = parser.parse_args()


zs = []
for filename in args.files:
    F = np.load(filename)
    freq = F['freq']
    zs.append(f212z(freq * 1e9))

zs = np.unique(zs)
Nzs = len(zs)
# define some color format sets for the plotter
# colors = ['r', 'g', 'k', 'b']
markers = ['o', 'd', '^', 's', 'v']
markers = itertools.cycle(markers)
# Create figure and prep subplot sizes for Delta^2
fig = plt.figure()
gs = gridspec.GridSpec(1, Nzs)
gs.update(hspace=0.0, wspace=0.2)
ax1 = [plt.subplot(gs[:, i]) for i in range(Nzs)]

# Create figure and prep subplots for P(k)
fig2 = plt.figure()
gs = gridspec.GridSpec(1, Nzs)
gs.update(hspace=0.0, wspace=0.2)
ax2 = [plt.subplot(gs[:, i]) for i in range(Nzs)]


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
        ax1[gs_ind].plot(noise_ks[gs_ind], noises[gs_ind], 'k-')

        d2_n = noises[gs_ind]
        pk_n = 2*np.pi**2/(np.array(noise_ks[gs_ind])**3)*d2_n
        if len(kpls[gs_ind]) > len(pk_n):
            ax2[gs_ind].plot(noise_ks[gs_ind], pk_n, '-', color='black')
            ax2[gs_ind].plot(-noise_ks[gs_ind], pk_n, '-', color='black')
        else:
            ax2[gs_ind].plot(noise_ks[gs_ind], pk_n, '-', color='black')


k_max = 0
k_par_max = 0
for filename in args.files:
    pspec_dict = np.load(filename)

    if np.max(pspec_dict['k']) > k_max:
        k_max = np.max(pspec_dict['k'])

    if np.max(np.abs(pspec_dict['kpl'])) > k_par_max:
        k_par_max = np.max(np.abs(pspec_dict['kpl']))

    # get special index for gridspec to plot all psecs on same z value
    redshift = f212z(pspec_dict['freq'] * 1e9)
    gs_ind = np.where(zs == redshift)[0].item()
    marker = markers.next()
    ax1[gs_ind].plot(pspec_dict['k'], pspec_dict['pI_fold'] +
                     pspec_dict['pI_fold_up'], '--', label='unweighted')
    ax1[gs_ind].errorbar(pspec_dict['k'], pspec_dict['pC_fold'],
                         pspec_dict['pC_fold_up'], label='weighted {0}'
                         .format(pspec_dict['pC_fold_prob']),
                         marker=marker)
    ax2[gs_ind].plot(pspec_dict['kpl'], pspec_dict['pI'] +
                     pspec_dict['pI_up'], '--', label='unweighted')
    ax2[gs_ind].errorbar(pspec_dict['kpl'], pspec_dict['pC'],
                         pspec_dict['pC_up'], label='weighted {0}'
                         .format(pspec_dict['pC_prob']),
                         marker=marker)


# set up some parameters to make the figures pretty
for gs_ind in xrange(Nzs):
    ax1[gs_ind].set_title('z = {0:.2f}'.format(zs[gs_ind]))
    ax1[gs_ind].set_ylabel('$\\frac{k^{3}}{2\pi^{2}} P(k) [mK]^{2}$')
    ax1[gs_ind].set_yscale('log')
    ax1[gs_ind].set_xlabel('$k$ [$h$ Mpc$^{-1}$]')
    ax1[gs_ind].set_xlim(0, k_max * 1.01)
    ax1[gs_ind].get_shared_y_axes().join(ax1[0], ax1[gs_ind])
    ax1[gs_ind].grid(True)

    ax2[gs_ind].set_yscale('log')
    ax2[gs_ind].set_title('z = {0:.2f}'.format(zs[gs_ind]))
    ax2[gs_ind].set_ylabel('$ P(k) \\frac{[mK]^{2}}{(hMpc^{-1})^{3}}$')
    ax2[gs_ind].set_xlim(-1.01 * k_par_max, k_par_max * 1.01)
    ax2[gs_ind].get_shared_y_axes().join(ax2[0], ax2[gs_ind])
    ax2[gs_ind].grid(True)

ax1[0].set_ylim([1e0, 1e9])


handles, labels = ax1[-1].get_legend_handles_labels()
ax1[-1].legend(handles, labels, loc='lower right', numpoints=1)
handles, labels = ax2[-1].get_legend_handles_labels()
ax2[-1].legend(handles, labels, loc='lower right', numpoints=1)

fig.savefig(args.outfile+'.png', format='png')
fig2.savefig(args.outfile+'_pk.png', format='png')
if args.plot:
    plt.show(block=True)
