#! /usr/bin/env python

import numpy as n
import matplotlib.pyplot as p
from scipy.optimize import curve_fit

# Reads in power spectrum results from projecting out 0,1,2... modes
# Plots power spectrum results before and after signal loss correction as a function of modes removed


# Read data
path = '/home/cacheng/capo/ctc/matt_data/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_RANGEOFPROJECTEDMODES'
nmodes=22
sense=5691450.72021 # XXX from plotting one of the "project_#_modes" directories
PS_i_up = []
PS_f_up = []
PS_i = []
PS_f = []
for mode_num in n.arange(0,nmodes,1):
    filename = path + '/project_' + str(mode_num) + '_modes'
    file_i = n.load(filename + '/inject_sep0,1_0.001/pspec_pk_k3pk.npz')
    file_f = n.load(filename + '/pspec_final_sep0,1.npz')
    kpl = file_i['kpl']
    k_ind = -4 # XXX hard-coded k index for k = 0.344
    k = kpl[k_ind] 
    PS_i_up.append(2*file_i['pCv_err'][k_ind]) # 2-sigma err
    PS_f_up.append(file_f['pC_up'][k_ind]) # 2-sigma err
    PS_i.append(n.abs(file_i['pCv'])[k_ind])
    PS_f.append(n.abs(file_f['pC'])[k_ind])

# Theory from Switzer et al. - first term only
fixmode = 3 # start fix at 3rd mode since first few modes are dominated by systematics
xs = n.arange(fixmode,nmodes,1.) # number of modes removed
err_theory_firstterm = 1./(1 - xs/nmodes)
normalization = PS_f_up[fixmode]/err_theory_firstterm[0]
err_theory_firstterm = err_theory_firstterm*normalization

# Fit N_ind (number of independent modes) in second term
def func(mode_num, N_ind):
    fit = 1./((1-mode_num/nmodes)*(1-mode_num/N_ind))*PS_f_up[fixmode]
    normalization = PS_f_up[fixmode]/fit[0]
    return fit*normalization
N_ind,_ = curve_fit(func, xs, PS_f_up[fixmode:], bounds=(0,1000))
err_theory_fit = 1./((1 - xs/nmodes)*(1 - xs/N_ind))
normalization = PS_f_up[fixmode]/err_theory_fit[0]
err_theory_fit = err_theory_fit*normalization
print "Fit for number of independent modes =", N_ind

# Force fit for full equation
if True:
    N_ind = 15
    err_theory_fit = 1./((1 - xs/nmodes)*(1 - xs/N_ind))
    normalization = PS_f_up[fixmode]/err_theory_fit[0]
    err_theory_fit = err_theory_fit*normalization

# Plot
p.plot(n.arange(0,nmodes,1), PS_i_up, color='0.5', label='Pre-signal loss correction, 2$\sigma$ err')
#p.plot(n.arange(0,nmodes,1), PS_i, color='0.5', marker='.', linestyle="None", label='Mean PS value')
p.plot(n.arange(0,nmodes,1), PS_f_up, 'k-', label='Post-signal loss correction, 2$\sigma$ err')
#p.plot(n.arange(0,nmodes,1), PS_f, 'k.', label='Mean PS value')
p.axhline(sense, color='g', label='Analytical')
p.plot(n.arange(fixmode,nmodes,1), err_theory_firstterm, 'b--', label='Theory from Switzer et al., only frequency modes')
p.plot(n.arange(fixmode,nmodes,1), err_theory_fit, 'b-', label='Theory from Switzer et al., both frequency and time modes')
p.xlabel('Number of modes removed')
p.ylabel('$P(k)$ $[mK^{2}(h^{-1} Mpc)^{3}]$')
p.legend(prop={'size':8}, loc='best')
p.title('k = ' +str(round(k,3))+' & N$_{samples}$ = '+str(N_ind))
p.yscale('log')
p.ylim(1e4,1e10)
p.grid()
p.show()

