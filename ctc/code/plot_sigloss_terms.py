#! /usr/bin/env python

import numpy as n
import matplotlib.pyplot as p
from scipy import stats

# Read file and get data
file = n.load('pspec_sigloss_terms.npz')
k = file['k']
Pins = n.array(file['Pins']).flatten() # P_in
ind = n.where(Pins > 0)[0] # only pos P_ins
ind2 = n.argsort(Pins[ind]) # sort if using both neg and pos
Pins = Pins[ind][ind2]
Pouts = n.array(file['Pouts']).flatten()[ind][ind2] # P_out
Pouts_I = n.array(file['Pouts_I']).flatten()[ind][ind2] # P_out for I case
pCvs_Cr = n.array(file['pCvs_Cr_fold']).flatten()[ind][ind2]
pCes_Cr = n.array(file['pCes_Cr_fold']).flatten()[ind][ind2]
pCves = n.array(file['pCves_fold']).flatten()[ind][ind2]
pIves = n.array(file['pIves_fold']).flatten()[ind][ind2]
pCrs = n.array(file['pCrs_fold']).flatten()[ind][ind2]
pCv = file['pCv'] # final PS points per k
pIv = file['pIv'] 

# Smooth curves
poly_dim = 4
#"""
    # absolute value
poly_data_Cr = n.polyfit(n.log10(n.abs(Pins)),n.log10(n.abs(pCvs_Cr)),poly_dim)
blue = 10**n.polyval(poly_data_Cr,n.log10(Pins))

red_ind_pos = n.where(pCves>=0)[0] # where positive pCves
red_ind_neg = n.where(pCves<0)[0] # where negative pCves
Pins_red_pos = Pins[red_ind_pos] # adjust Pins
Pins_red_neg = Pins[red_ind_neg]
poly_cross_pos = n.polyfit(n.log10(Pins_red_pos),n.log10(pCves[red_ind_pos]),poly_dim)
poly_cross_neg = n.polyfit(n.log10(Pins_red_neg),n.log10(-pCves[red_ind_neg]),poly_dim)
red_neg = -10**n.polyval(poly_cross_neg,n.log10(Pins_red_neg))
pos_1 = n.where(Pins_red_pos < 1e8)[0] # XXX hack to separate two pos components
pos_2 = n.where(Pins_red_pos > 1e8)[0]
#poly_cross_pos_1 = n.polyfit(n.log10(Pins[red_ind_pos][pos_1]),n.log10(pCves[red_ind_pos][pos_1]),poly_dim)
#poly_cross_pos_2 = n.polyfit(n.log10(Pins[red_ind_pos][pos_2]),n.log10(pCves[red_ind_pos][pos_2]),poly_dim)
red_pos_1 = 10**n.polyval(poly_cross_pos,n.log10(Pins_red_pos[pos_1]))
red_pos_2 = 10**n.polyval(poly_cross_pos,n.log10(Pins_red_pos[pos_2]))

red_I_ind_pos = n.where(pIves>=0)[0] # where positive pIves
red_I_ind_neg = n.where(pIves<0)[0] # where negative pIves
Pins_red_pos_I = Pins[red_I_ind_pos] # adjust Pins
Pins_red_neg_I = Pins[red_I_ind_neg]
poly_cross_I_pos = n.polyfit(n.log10(Pins_red_pos_I),n.log10(pIves[red_I_ind_pos]),poly_dim)
poly_cross_I_neg = n.polyfit(n.log10(Pins_red_neg_I),n.log10(-pIves[red_I_ind_neg]),poly_dim)
red_I_pos = 10**n.polyval(poly_cross_I_pos,n.log10(Pins_red_pos_I))
red_I_neg = -10**n.polyval(poly_cross_I_neg,n.log10(Pins_red_neg_I))

poly_eor_Cr = n.polyfit(n.log10(n.abs(Pins)),n.log10(n.abs(pCes_Cr)),poly_dim)
mag = 10**n.polyval(poly_eor_Cr,n.log10(Pins))

poly_Pouts = n.polyfit(n.log10(n.abs(Pins)),n.log10(n.abs(Pouts)),poly_dim)
black = 10**n.polyval(poly_Pouts,n.log10(Pins))

poly_Pouts_I = n.polyfit(n.log10(n.abs(Pins)),n.log10(n.abs(Pouts_I)),poly_dim)
black_I = 10**n.polyval(poly_Pouts_I,n.log10(Pins))

poly_pCr = n.polyfit(n.log10(n.abs(Pins)),n.log10(n.abs(pCrs)),poly_dim)
yellow = 10**n.polyval(poly_pCr,n.log10(Pins))
#"""
"""
    # log of abs * sign
poly_data_Cr = n.polyfit(n.log10(n.abs(Pins))*n.sign(Pins),n.log10(n.abs(pCvs_Cr))*n.sign(pCvs_Cr),poly_dim)
blue = n.sign(pCvs_Cr)*10**(n.sign(pCvs_Cr)*n.polyval(poly_data_Cr,n.log10(Pins)))
poly_cross = n.polyfit(n.log10(n.abs(Pins))*n.sign(Pins),n.log10(n.abs(pCves))*n.sign(pCves),poly_dim)
red = n.sign(pCves)*10**(n.sign(pCves)*n.polyval(poly_cross,n.log10(Pins)))
poly_cross_I = n.polyfit(n.log10(n.abs(Pins))*n.sign(Pins),n.log10(n.abs(pIves))*n.sign(pIves),poly_dim)
red_I = n.sign(pIves)*10**(n.sign(pIves)*n.polyval(poly_cross_I,n.log10(Pins)))
poly_eor_Cr = n.polyfit(n.log10(n.abs(Pins))*n.sign(Pins),n.log10(n.abs(pCes_Cr))*n.sign(pCes_Cr),poly_dim)
mag = n.sign(pCes_Cr)*10**(n.sign(pCes_Cr)*n.polyval(poly_eor_Cr,n.log10(Pins)))
poly_Pouts = n.polyfit(n.log10(n.abs(Pins))*n.sign(Pins),n.log10(n.abs(Pouts))*n.sign(Pouts),poly_dim)
black = n.sign(Pouts)*10**(n.sign(Pouts)*n.polyval(poly_Pouts,n.log10(Pins)))
poly_Pouts_I = n.polyfit(n.log10(n.abs(Pins))*n.sign(Pins),n.log10(n.abs(Pouts_I))*n.sign(Pouts_I),poly_dim)
black_I = n.sign(Pouts_I)*10**(n.sign(Pouts_I)*n.polyval(poly_Pouts_I,n.log10(Pins)))
poly_pCr = n.polyfit(n.log10(n.abs(Pins))*n.sign(Pins),n.log10(n.abs(pCrs))*n.sign(pCrs),poly_dim)
yellow = n.sign(pCrs)*10**(n.sign(pCrs)*n.polyval(poly_pCr,n.log10(Pins)))
"""
"""
    # convolve with boxcar
window_len = 100 # out of 500
window = n.ones(window_len)
window /= window.sum()
blue_con = n.convolve(window,n.log10(n.abs(pCvs_Cr))*n.sign(pCvs_Cr),mode='same')
blue = n.sign(pCvs_Cr)*10**(n.sign(pCvs_Cr)*blue_con)
red_con = n.convolve(window,n.log10(n.abs(pCves))*n.sign(pCves),mode='same')
red = n.sign(pCves)*10**(n.sign(pCves)*red_con)
red_con_I = n.convolve(window,n.log10(n.abs(pIves))*n.sign(pIves),mode='same')
red_I = n.sign(pIves)*10**(n.sign(pIves)*red_con_I)
mag_con = n.convolve(window,n.log10(n.abs(pCes_Cr))*n.sign(pCes_Cr),mode='same')
mag = n.sign(pCes_Cr)*10**(n.sign(pCes_Cr)*mag_con)
black_con = n.convolve(window,n.log10(n.abs(Pouts))*n.sign(Pouts),mode='same')
black = n.sign(Pouts)*10**(n.sign(Pouts)*black_con)
black_con_I = n.convolve(window,n.log10(n.abs(Pouts_I))*n.sign(Pouts_I),mode='same')
black_I = n.sign(Pouts_I)*10**(n.sign(Pouts_I)*black_con)
yellow_con = n.convolve(window,n.log10(n.abs(pCrs))*n.sign(pCrs),mode='same')
yellow = n.sign(pCrs)*10**(n.sign(pCrs)*yellow_con)
"""



# Plot

#P_out = xx   +  xe +  ex +    ee   - xx (C_x)
#      = blue + red + red + magenta - grey

p.figure(figsize=(14,8))
low=1e2;high=1e16

p.subplot(121)

p.plot(Pins,black,'k-',label='$P_{out}$',linewidth=3)
p.plot(Pins,blue,'b-',linewidth=3,label='$\propto xC_{r}^{-1}QC_{r}^{-1}x$')
#p.plot(Pins_red_pos,red_pos,'r-',linewidth=3,label='$\propto xC_{r}^{-1}QC_{r}^{-1}e$')
p.plot(Pins_red_neg,red_neg,'r-',linewidth=3,label='$\propto xC_{r}^{-1}QC_{r}^{-1}e$')
p.plot(Pins_red_pos[pos_1],red_pos_1,'r-',linewidth=3)
p.plot(Pins_red_pos[pos_2],red_pos_2,'r-',linewidth=3)
p.plot(Pins,mag,'g-',linewidth=3,label='$\propto eC_{r}^{-1}QC_{r}^{-1}e$')
p.axhline(pCv,color='0.5',linestyle='-',linewidth=3,label='$\propto xC_{x}^{-1}QC_{x}^{-1}x$')
#p.plot(Pins, blue + 2*red + magenta - pCv, 'k--')
#p.plot(Pins,yellow,'y-',linewidth=3,label='$\propto rC_{r}^{-1}QC_{r}^{-1}r$')
   # Test all points
#p.plot(Pins,Pouts,'k.')
#p.plot(Pins,pCvs_Cr,'b.')
#p.plot(Pins,pCves,'r.')
#p.plot(Pins,pCes_Cr,'m.')
p.xlabel('$P_{in}$ $[mK^{2}(h^{-1} Mpc)^{3}]$', fontsize=18)
p.ylabel('$P_{out} $ $[mK^{2}(h^{-1} Mpc)^{3}]$', fontsize=18)
p.xscale('log');p.yscale('symlog',linthreshy=1e2)
p.xlim(low,high);p.ylim(-high,high)
p.plot([low,high], [low,high], 'k:')  # diagonal line
p.plot([low,high], [-low,-high], 'k:')  # diagonal line
p.tick_params(axis='both', which='major', labelsize=14)
ttl = p.title("Inverse Covariance Weighting, k = " + str(n.round(k,3)) + " h Mpc$^{-1}$", fontsize=14)
ttl.set_position([.5, 1.03])
p.legend(prop={'size':14},loc=3,numpoints=1,ncol=2)
p.grid()

p.subplot(122)
p.plot(Pins,black_I,'k-',label='$P_{out}$',linewidth=3)
p.plot(Pins_red_pos_I,red_I_pos,'r-',linewidth=3,label='$\propto xIQIe$')
p.plot(Pins_red_neg_I,red_I_neg,'r-',linewidth=3)
p.plot(Pins,Pins,'g-',linewidth=3,label='$\propto eIQIe$')
p.axhline(pIv,color='0.5',linestyle='-',linewidth=3,label='$\propto xIQIx$')
p.xlabel('$P_{in}$ $[mK^{2}(h^{-1} Mpc)^{3}]$', fontsize=18)
p.ylabel('$P_{out} $ $[mK^{2}(h^{-1} Mpc)^{3}]$', fontsize=18)
p.xscale('log');p.yscale('symlog',linthreshy=1e2)
p.xlim(low,high);p.ylim(-high,high)
p.plot([low,high], [low,high], 'k:')  # diagonal line
p.plot([low,high], [-low,-high], 'k:')  # diagonal line
p.tick_params(axis='both', which='major', labelsize=14)
p.legend(prop={'size':14},loc=3,numpoints=1,ncol=2)
p.grid()
ttl = p.title("Unweighted, k = " + str(n.round(k,3)) + " h Mpc$^{-1}$", fontsize=14)
ttl.set_position([.5, 1.03])

p.tight_layout()
p.show()

