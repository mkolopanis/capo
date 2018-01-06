#! /usr/bin/env python

import numpy as n
import matplotlib.pyplot as p
from scipy.optimize import curve_fit

# Reads in power spectrum results from projecting out 0,1,2... modes
# Plots power spectrum results before and after signal loss correction as a function of modes removed

if False: # eigenmodes set to 1
    path = '/data4/paper/ctc/PSA64/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_RANGEOFPROJECTEDMODES'
    startmode=0
    nmodes=21
    deltamode=1
    xlabel='Number of modes down-weighted using inverse covariance weighting'
    f1 = '/project_'
    f2 = '_modes'
    flipx = False

if False: # eigenmodes set to 0
    path = '/data4/paper/ctc/PSA64/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_RANGEOFPROJECTEDMODES_v2'
    startmode=1
    nmodes=20
    deltamode=1
    xlabel='Number of modes down-weighted'
    f1 = '/project_'
    f2 = '_modes'
    flipx = False

if True: # added identity
    path = '/data4/paper/ctc/PSA64/PAPER_METHODS/rangeofaddedidentity_trace'
    startmode=0 
    nmodes=0.2 #20000 
    deltamode=0.01 #1000
    xlabel='Strength of identity added (percentage of Tr(C) added to C)'
    f1 = '/add_'
    f2 = '_identity'
    flipx = True

# Read files
sense=4436767.36822*2 # XXX from plotting one of the "project_#_modes" directories
PS_i_up = []
PS_f_up = []
PS_i = []
PS_f = []
for mode_num in n.arange(startmode,nmodes,deltamode):
    filename = path + f1 + str(mode_num) + f2
    print 'Reading', filename
    print mode_num
    #file_i = n.load(filename + '/inject_sep0,1_0.001/pspec_pk_k3pk.npz')
    file_i = n.load(filename + '/pspec_final_sep0,1_nosigloss.npz')
    file_f = n.load(filename + '/pspec_final_sep0,1.npz')
    kpl = file_i['kpl']
    k_ind = -3 # XXX hard-coded k index for k = 0.344
    k = kpl[k_ind] 
    PS_i_up.append(2*file_i['pCv_err'][k_ind]) # 2-sigma err
    PS_f_up.append(2*file_f['pCv_err'][k_ind]) # 2-sigma err
    PS_i.append(n.abs(file_i['pCv'])[k_ind])
    PS_f.append(n.abs(file_f['pCv'])[k_ind])
    #if mode_num == 0:
    #    PS_f_up[0] = PS_i_up[0]

"""
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
"""

# Best PS (Identity Mult)
f = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_IDENTITYMULTWEIGHT_WEIGHTI/pspec_final_sep0,1.npz')
ps_mult = n.abs(f['pCv'][k_ind]) + 2*f['pCv_err'][k_ind]

# Best PS (Identity Add)
f = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/rangeofaddedidentity/add_2000_identity/pspec_final_sep0,1.npz')
ps_add = n.abs(f['pCv'][k_ind]) + 2*f['pCv_err'][k_ind]

# Plot
p.figure(figsize=(12,8))
#p.plot(n.arange(startmode,nmodes,deltamode), n.array(PS_f), 'r.')
#p.plot(n.arange(startmode,nmodes,deltamode), n.array(PS_i), 'm.')
p.plot(n.arange(startmode,nmodes,deltamode), n.array(PS_i) + n.array(PS_i_up), color='0.5', linewidth=3, label='Pre-signal loss estimation')
p.plot(n.arange(startmode,nmodes,deltamode), n.array(PS_f) + n.array(PS_f_up), 'k-', linewidth=3, label='Post-signal loss estimation')
p.axhline(file_f['pIv'][k_ind]+2*file_f['pIv_err'][k_ind],color='b',linestyle='--',linewidth=3,label='Unweighted')
p.axhline(ps_mult,color='r',linestyle='-',linewidth=3,label='$\hat{C} = \hat{C} \circ I$')
#p.axhline(ps_add,color='c',linestyle='-',linewidth=3,label='$\hat{C} = \hat{C} + 2000I$')
p.axhline(sense,color='g',linestyle='-',linewidth=3,label='Analytical $2\sigma$ Error')
#p.plot(n.arange(fixmode,nmodes,1), err_theory_firstterm, 'b--', label='Theory from Switzer et al., only frequency modes')
#p.plot(n.arange(fixmode,nmodes,1), err_theory_fit, 'b-', label='Theory from Switzer et al., both frequency and time modes')
p.xlabel(xlabel,fontsize=14)
p.ylabel('$P(k)$ $[mK^{2}(h^{-1} Mpc)^{3}]$',fontsize=14)
p.legend(prop={'size':14}, loc='best')
#p.title('k = ' +str(round(k,3))+' & N$_{samples}$ = '+str(N_ind))
p.title('k = ' +str(round(k,3)),fontsize=14)
p.tick_params(axis='both', which='major', labelsize=14)
p.yscale('log')
#p.ylim(1e4,1e10)
if flipx == True: 
    p.gca().invert_xaxis()
p.grid()
p.show()

