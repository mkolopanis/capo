#! /usr/bin/env python

import numpy as n
import pylab as p

file_eo = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_JACKKNIFE_EVENODD/pspec_final_sep0,1.npz')
file_bl = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_JACKKNIFE_BASELINES/pspec_final_sep0,1.npz')
file_lst = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_JACKKNIFE_LST_FIRSTLAST/pspec_final_sep0,1.npz')
data = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/PSA64_FRF_RA.5_8.6_CHAN95_115_SEP0,1_IDENTITYMULTWEIGHT/pspec_final_sep0,1_final.npz')


p.errorbar(file_eo['kpl'],file_eo['pCv'],file_eo['pCv_err']*2,linestyle='',marker='.',color='b',label='Even/Odd Null Test')
p.errorbar(file_bl['kpl']+0.005,file_bl['pCv'],file_bl['pCv_err']*2,linestyle='',marker='.',color='g',label='Baselines Null Test')
p.errorbar(file_lst['kpl']+0.01,file_lst['pCv'],file_lst['pCv_err']*2,linestyle='',marker='.',color='m',label='LST Null Test')
p.errorbar(data['kpl']-0.005,data['pCv'],data['pCv_err']*2,linestyle='',marker='.',color='k',label='Original Data')

noise = 4436767.36822 # XXX

p.axhline(0,color='k',linestyle='--')
p.axhspan(-noise*2,noise*2,facecolor='0.5',edgecolor="none",alpha=0.5,label='Estimated $2\sigma$ Error')
p.legend(numpoints=1,prop={'size':14},ncol=2)
p.ylabel('$P(k)$ $[mK^{2}(h^{-1} Mpc)^{3}]$',fontsize=14)
p.xlabel('$k_{\\parallel}$ [$h$ Mpc$^{-1}$]',fontsize=14)
p.tick_params(axis='both', which='major', labelsize=14)
p.ylim(-3.5e8,10e8)
p.grid()
p.show()

