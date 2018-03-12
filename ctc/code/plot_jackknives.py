#! /usr/bin/env python

import numpy as n
import pylab as p

#file_eo = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_evenodd_addI/pspec_final_sep0,1.npz')
#file_bl = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_baselines_addI/pspec_final_sep0,1.npz')
#file_lst = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_lst_addI/pspec_final_sep0,1.npz')
file_eo = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_evenodd_addI/pspec_final_sep0,1_nosigloss.npz')
file_bl = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_baselines_addI/pspec_final_sep0,1_nosigloss.npz')
file_lst = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_lst_addI/pspec_final_sep0,1_nosigloss.npz')
data = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/0.04addedidentity_correcttsys/pspec_final_sep0,1.npz')

noise = 4436767.36822 # XXX

# Figure 1: zoom-out
p.figure()
p.errorbar(file_eo['kpl'],file_eo['pCv'],file_eo['pCv_err']*2,linestyle='',marker='.',color='b',label='Even/Odd Null Test')
p.errorbar(file_bl['kpl']+0.005,file_bl['pCv'],file_bl['pCv_err']*2,linestyle='',marker='.',color='g',label='Baselines Null Test')
p.errorbar(file_lst['kpl']+0.01,file_lst['pCv'],file_lst['pCv_err']*2,linestyle='',marker='.',color='m',label='LST Null Test')
p.errorbar(data['kpl']-0.005,data['pCv'],data['pCv_err']*2,linestyle='',marker='.',color='k',label='Original Data')
p.axhline(0,color='k',linestyle='--')
p.axhspan(-noise*2,noise*2,facecolor='0.5',edgecolor="none",alpha=0.5,label='Estimated $2\sigma$ Error')
p.legend(numpoints=1,prop={'size':14},ncol=2)
p.ylabel('$P(k)$ $[mK^{2}(h^{-1} Mpc)^{3}]$',fontsize=14)
p.xlabel('$k_{\\parallel}$ [$h$ Mpc$^{-1}$]',fontsize=14)
p.tick_params(axis='both', which='major', labelsize=14)
p.ylim(-2e8,1.5e9)
p.grid()

# Figure 2: zoom-in
p.figure()
p.errorbar(file_eo['kpl'],file_eo['pCv'],file_eo['pCv_err']*2,linestyle='',marker='.',color='b',label='Even/Odd Null Test')
p.errorbar(file_bl['kpl']+0.005,file_bl['pCv'],file_bl['pCv_err']*2,linestyle='',marker='.',color='g',label='Baselines Null Test')
p.errorbar(file_lst['kpl']+0.01,file_lst['pCv'],file_lst['pCv_err']*2,linestyle='',marker='.',color='m',label='LST Null Test')
#p.errorbar(data['kpl']-0.005,data['pCv'],data['pCv_err']*2,linestyle='',marker='.',color='k',label='Original Data')
p.axhline(0,color='k',linestyle='--')
p.axhspan(-noise*2,noise*2,facecolor='0.5',edgecolor="none",alpha=0.5,label='Estimated $2\sigma$ Error')
p.legend(numpoints=1,prop={'size':14},ncol=2)
p.ylabel('$P(k)$ $[mK^{2}(h^{-1} Mpc)^{3}]$',fontsize=14)
p.xlabel('$k_{\\parallel}$ [$h$ Mpc$^{-1}$]',fontsize=14)
p.tick_params(axis='both', which='major', labelsize=14)
p.ylim(-1.e8,1.2e8)
p.grid()

p.show()

