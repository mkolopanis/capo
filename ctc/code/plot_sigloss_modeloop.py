#! /usr/bin/env python

import numpy as n
import matplotlib.pyplot as p

# Reads in power spectrum results from projecting out 0,1,2... modes
# Plots power spectrum results before and after signal loss correction as a function of modes removed

path = '/home/cacheng/capo/ctc/matt_data/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_RANGEOFPROJECTEDMODES'
nmodes=22
sense=5691450.72021 # XXX from plotting one of the "project_#_modes" directories
PS_i = []
PS_f = []
for mode_num in n.arange(0,nmodes,1):
    filename = path + '/project_' + str(mode_num) + '_modes'
    file_i = n.load(filename + '/inject_sep0,1_0.001/pspec_pk_k3pk.npz')
    file_f = n.load(filename + '/pspec_final_sep0,1.npz')
    kpl = file_i['kpl']
    k_ind = -4 # XXX hard-coded k index for k = 0.344
    k = kpl[k_ind] 
    PS_i.append(n.abs((file_i['pCv']) + 2*file_i['pCv_err'])[k_ind]) # 2-sigma
    PS_f.append(n.abs((file_f['pC']) + file_f['pC_up'])[k_ind]) # 2-sigma
    
p.plot(n.arange(0,nmodes,1), PS_i, color='0.5', label='Pre-signal loss correction')
p.plot(n.arange(0,nmodes,1), PS_f, 'k-', label='Post-signal loss correction')
p.axhline(sense, color='g', label='Analytical')
p.xlabel('Number of modes removed')
p.ylabel('$P(k)$ $[mK^{2}(h^{-1} Mpc)^{3}]$')
p.legend(prop={'size':8}, loc='best')
p.title('k = ' +str(k))
p.yscale('log')
p.ylim(1e3,1e9)
p.grid()
p.show()

