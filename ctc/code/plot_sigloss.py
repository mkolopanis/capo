#! /usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
from scipy import stats

# Read file and get data
file = n.load('pspec_sigloss.npz')
#file = n.load('pspec_sigloss_nosubtract.npz')
k = file['k']
file2 = n.load('inject_sep0,1_0.001/pspec_pk_k3pk.npz')
file_pts = file2['pCv']
file_pts_I = file2['pIv']
k_pts_ind = n.where(file2['kpl'] == k)[0][0]
Pins = n.array(file['Pins']).flatten() # P_in
Pouts = n.array(file['Pouts']).flatten() # P_out
Pouts_I = n.array(file['Pouts_I']).flatten() # P_out for I case
ind = n.where(Pins > 0)[0] # positive P_ins only
Pins = Pins[ind]
Pouts = Pouts[ind]
Pouts_I = Pouts_I[ind]
bins_concat = file['bins_concat']
binsx = file['binsx']
binsy = file['binsy']
#binsx = file['bins']
#binsy = file['bins']

"""
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
"""

#KDE
Pouts_pos = Pouts[n.where(Pouts > 0)[0]]
Pouts_neg = -Pouts[n.where(Pouts < 0)[0]]
Pins_pos = Pins[n.where(Pouts > 0)[0]]
Pins_neg = Pins[n.where(Pouts < 0)[0]]
Pouts_pos_I = Pouts_I[n.where(Pouts_I > 0)[0]]
Pouts_neg_I = -Pouts_I[n.where(Pouts_I < 0)[0]]
Pins_pos_I = Pins[n.where(Pouts_I > 0)[0]]
Pins_neg_I = Pins[n.where(Pouts_I < 0)[0]]

ygrid,xgrid = n.meshgrid(binsy,binsx) # create grid on which to sample
positions = n.vstack([xgrid.ravel(),ygrid.ravel()])

kernel_C_pos = stats.gaussian_kde((n.log10(Pins_pos),n.log10(Pouts_pos)),bw_method='scott')
factor = kernel_C_pos.factor#+0.2
kernel_C_pos = stats.gaussian_kde((n.log10(Pins_pos),n.log10(Pouts_pos)),bw_method=factor)
M_pos = n.reshape(kernel_C_pos(positions).T,xgrid.shape).T

kernel_I_pos = stats.gaussian_kde((n.log10(Pins_pos_I),n.log10(Pouts_pos_I)),bw_method=factor)
M_pos_I = n.reshape(kernel_I_pos(positions).T,xgrid.shape).T
try: 
    kernel_C_neg = stats.gaussian_kde((n.log10(Pins_neg),n.log10(Pouts_neg)),bw_method=factor)
    M_neg = n.reshape(kernel_C_neg(positions).T,xgrid.shape).T
except: # no negative values
    M_neg = n.zeros((binsx.size,binsy.size))
try:
    kernel_I_neg = stats.gaussian_kde((n.log10(Pins_neg_I),n.log10(Pouts_neg_I)),bw_method=factor)
    M_neg_I = n.reshape(kernel_I_neg(positions).T,xgrid.shape).T
except: # no negative values
    M_neg_I = n.zeros((binsx.size,binsy.size))

M = n.concatenate((M_neg[::-1],M_pos))
M_I = n.concatenate((M_neg_I[::-1],M_pos_I))

for col in range(M.shape[1]): # normalize
    M[:,col] /= n.sum(M[:,col])
    M_I[:,col] /= n.sum(M_I[:,col])

# Plot signal loss transfer curves
pklo=1e2
pkhi=1e13
binsy_full = n.concatenate((-10**binsy[::-1],10**binsy))

plt.figure(figsize=(12,6))
plt.subplot(121)
#plt.plot(n.abs(Pins.flatten()),n.abs(Pouts.flatten()),'k.',label="All bootstraps")
#plt.plot(10**binsx,10**n.polyval(coeff,binsx),'r-',label="Polynomial Fit")
plt.pcolormesh(10**binsx,binsy_full,M,cmap='hot_r')
plt.plot([pklo, pkhi], [pklo, pkhi], 'k:')  # diagonal line
plt.plot([pklo, pkhi], [-pklo,-pkhi], 'k:')  # diagonal line
plt.hlines(y=n.abs(file_pts[k_pts_ind]),xmin=pklo,xmax=pkhi,color='0.5',linewidth=3) # peak of original distribution
#plt.legend(numpoints=1,prop={'size':10},loc='best')
#plt.grid()
plt.yscale('symlog',linthreshy=1e2)
plt.xscale('log')
plt.xlabel('$k_{\\parallel}$ [$h$ Mpc$^{-1}$]')
ttl = plt.title("Inverse Covariance Weighting, k = " + str(n.round(k,3)) + " h Mpc$^{-1}$", fontsize=14)
ttl.set_position([.5, 1.03])
#plt.xlabel('$P_{in}$ $[mK^{2}(h^{-1} Mpc)^{3}]$', fontsize=18)
plt.ylabel('$P_{out}$ $[mK^{2}(h^{-1} Mpc)^{3}]$', fontsize=18)
#plt.ylabel('$\widehat{p}_{r}$ $[mK^{2}(h^{-1} Mpc)^{3}]$', fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.xlim(pklo,pkhi)
plt.ylim(-pkhi,pkhi)
plt.grid()

plt.subplot(122)
#plt.plot(n.abs(Pins.flatten()),n.abs(Pouts_I.flatten()),'k.',label="All bootstraps")
#plt.plot(10**binsx,10**n.polyval(coeff_I,binsx),'r-',label="Polynomial Fit")
plt.hlines(y=file_pts_I[k_pts_ind],xmin=pklo,xmax=pkhi,color='0.5',linewidth=3)
plt.pcolormesh(10**binsx,binsy_full,M_I,cmap='hot_r')
plt.plot([pklo, pkhi], [pklo, pkhi], 'k:')  # diagonal line
plt.plot([pklo, pkhi], [-pklo,-pkhi], 'k:')  # diagonal line
#plt.legend(numpoints=1,prop={'size':10},loc='best')
#plt.grid()
plt.yscale('symlog',linthreshy=1e2)
plt.xscale('log')
#plt.yscale('symlog',linthreshy=5e7)
#plt.xscale('symlog',linthreshx=5e7)
ttl = plt.title("Unweighted, k = " + str(n.round(k,3)) + " h Mpc$^{-1}$", fontsize=14)
ttl.set_position([.5, 1.03])
plt.xlabel('$P_{in}$ $[mK^{2}(h^{-1} Mpc)^{3}]$', fontsize=18)
plt.ylabel('$P_{out}$ $[mK^{2}(h^{-1} Mpc)^{3}]$', fontsize=18)
#plt.ylabel('$\widehat{p}_{r}$ $[mK^{2}(h^{-1} Mpc)^{3}]$', fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.xlim(pklo,pkhi)
plt.ylim(-pkhi,pkhi)
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


