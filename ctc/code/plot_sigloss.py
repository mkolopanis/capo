#! /usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
from scipy import stats

# Read file and get data
file = n.load('pspec_sigloss.npz')
k = file['k']
file2 = n.load('inject_sep0,1_0.001/pspec_pk_k3pk.npz')
file_pts = file2['pCv']
file_pts_I = file2['pIv']
k_pts_ind = n.where(file2['kpl'] == k)[0][0]
#pC = file['pC'] # final PS value
#pI = file['pI']
#pC_err = file['pC_err']
#pI_err = file['pI_err']
#new_pCs = file['new_pCs'] # distributions of PS values (after sigloss correction)
#new_pIs = file['new_pIs']
#old_pCs = file['old_pCs'] # distribution of PS values (before sigloss correction)
#old_pIs = file['old_pIs']
Pins = file['Pins'] # P_in
Pouts = file['Pouts'] # P_out
Pouts_I = file['Pouts_I'] # P_out for I case
bins_concat = file['bins_concat']
bins = file['bins']

# Fit polynomial to mean signal loss transfer curve
xs = n.log10(n.abs(Pins.flatten())) # absolute value since symmetric
ys = n.log10(n.abs(Pouts.flatten()))
ys_I = n.log10(n.abs(Pouts_I.flatten()))
xs = n.append(n.repeat(0,100000),xs) # force fit to go through zero
ys = n.append(n.repeat(0,100000),ys)
ys_I = n.append(n.repeat(0,100000),ys_I)
order = n.argsort(xs) # re-order after padding
xs = xs[order]
ys = ys[order]
ys_I = ys_I[order]
coeff = n.polyfit(xs,ys,8) # coefficients from highest to lowest order
coeff_I = n.polyfit(xs,ys_I,8)

#KDE
ygrid,xgrid = n.meshgrid(bins,bins) # create grid on which to sample
positions = n.vstack([xgrid.ravel(),ygrid.ravel()])
kernel_C = stats.gaussian_kde((n.log10(n.abs(Pins.flatten())),n.log10(n.abs(Pouts.flatten()))),bw_method='scott')
kernel_I = stats.gaussian_kde((n.log10(n.abs(Pins.flatten())),n.log10(n.abs(Pouts_I.flatten()))),bw_method=kernel_C.factor)
M = n.reshape(kernel_C(positions).T,xgrid.shape).T
M_I = n.reshape(kernel_I(positions).T,xgrid.shape).T

# Plot signal loss transfer curves
pklo=1e2
pkhi=1e13

plt.figure(figsize=(12,6))
plt.subplot(121)
#plt.plot(n.abs(Pins.flatten()),n.abs(Pouts.flatten()),'k.',label="All bootstraps")
#plt.plot(10**bins,10**n.polyval(coeff,bins),'r-',label="Polynomial Fit")
plt.pcolormesh(10**bins,10**bins,M,cmap='hot_r')
plt.plot([pklo, pkhi], [pklo, pkhi], 'k-')  # diagonal line
plt.hlines(y=n.abs(file_pts[k_pts_ind]),xmin=pklo,xmax=pkhi,color='0.5',linewidth=3) # peak of original distribution
#plt.legend(numpoints=1,prop={'size':10},loc='best')
#plt.grid()
plt.yscale('log')
plt.xscale('log')
#plt.xlabel('$k_{\\parallel}$ [$h$ Mpc$^{-1}$]')
ttl = plt.title("Inverse Covariance Weighting, k = " + str(n.round(k,3)) + " h Mpc$^{-1}$", fontsize=14)
ttl.set_position([.5, 1.03])
plt.xlabel('$P_{in}$ $[mK^{2}(h^{-1} Mpc)^{3}]$', fontsize=14)
plt.ylabel('$P_{out}$ $[mK^{2}(h^{-1} Mpc)^{3}]$', fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.xlim(pklo,pkhi)
plt.ylim(pklo,pkhi)
plt.grid()

plt.subplot(122)
#plt.plot(n.abs(Pins.flatten()),n.abs(Pouts_I.flatten()),'k.',label="All bootstraps")
#plt.plot(10**bins,10**n.polyval(coeff_I,bins),'r-',label="Polynomial Fit")
plt.hlines(y=file_pts_I[k_pts_ind],xmin=pklo,xmax=pkhi,color='0.5',linewidth=3)
plt.pcolormesh(10**bins,10**bins,M_I,cmap='hot_r')
plt.plot([pklo, pkhi], [pklo, pkhi], 'k-')  # diagonal line
#plt.legend(numpoints=1,prop={'size':10},loc='best')
#plt.grid()
plt.yscale('log')
plt.xscale('log')
#plt.yscale('symlog',linthreshy=5e7)
#plt.xscale('symlog',linthreshx=5e7)
ttl = plt.title("Unweighted, k = " + str(n.round(k,3)) + " h Mpc$^{-1}$", fontsize=14)
ttl.set_position([.5, 1.03])
plt.xlabel('$P_{in}$ $[mK^{2}(h^{-1} Mpc)^{3}]$', fontsize=14)
plt.ylabel('$P_{out}$ $[mK^{2}(h^{-1} Mpc)^{3}]$', fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.xlim(pklo,pkhi)
plt.ylim(pklo,pkhi)
plt.tight_layout()
plt.grid()
plt.show()

"""
# Plot old vs. new distributions
plt.figure(figsize=(12,6))
plt.subplot(121)
plt.plot(bins_concat,old_pCs/n.max(old_pCs),'0.5')
plt.title("Inverse Covariance Weighting, k = " + str(n.round(k,3)))
plt.xlabel('$P(k)$ $[mK^{2}(h^{-1} Mpc)^{3}]$')
plt.plot(bins_concat,new_pCs/n.max(new_pCs),'k-')
#plt.xlim(-1e9,1e9)
plt.xscale('symlog')
plt.subplot(122)
plt.plot(bins_concat,old_pIs/n.max(old_pIs),'0.5',label='Pre-signal loss correction')
plt.plot(bins_concat,new_pIs/n.max(new_pIs),'k-',label='Post-signal loss correction')
plt.title("Unweighted, k = " + str(n.round(k,3)))
plt.xlabel('$P(k)$ $[mK^{2}(h^{-1} Mpc)^{3}]$')
plt.xlim(-1e8,1e8)
plt.legend(numpoints=1,prop={'size':12},loc='best')
plt.tight_layout()
plt.show()
"""


