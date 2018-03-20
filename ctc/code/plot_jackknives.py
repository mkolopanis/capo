#! /usr/bin/env python

import numpy as n
import pylab as p

# Old jackknives ((e+o) x (e-o) for example)
#file_eo = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_evenodd_addI/pspec_final_sep0,1.npz')
#file_bl = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_baselines_addI/pspec_final_sep0,1.npz')
#file_lst = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_lst_addI/pspec_final_sep0,1.npz')
file_eo = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_evenodd_addI/pspec_final_sep0,1_nosigloss.npz')
file_bl = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_baselines_addI/pspec_final_sep0,1_nosigloss.npz')
file_lst = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_lst_addI/pspec_final_sep0,1_nosigloss.npz')
data = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/rangeofaddedidentity_trace_log_DENSE_NOSEED/add_0.0206913808111_identity/pspec_final_sep0,1.npz')


# New jackknives (differencing of 2 PS)
file_lst_diff = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_lst/pspec_final_sep0,1_full.npz')
file_lst_1 = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_lst/lst_1/pspec_final_sep0,1_full.npz')
file_lst_2 = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_lst/lst_2/pspec_final_sep0,1_full.npz')
file_bls_diff = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_bls/pspec_final_sep0,1_full.npz')
file_bls_1 = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_bls/bls_1/pspec_final_sep0,1_full.npz')
file_bls_2 = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_bls/bls_2/pspec_final_sep0,1_full.npz')
file_eo_diff = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_evenodd/pspec_final_sep0,1_full.npz')
file_eo_1 = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_evenodd/even/pspec_final_sep0,1_full.npz')
file_eo_2 = n.load('/data4/paper/ctc/PSA64/PAPER_METHODS/jackknife_evenodd/odd/pspec_final_sep0,1_full.npz')


"""
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
"""

# LST null test
noise_lst = 7979656.32115
p.figure()
p.errorbar(file_lst_diff['kpl'],file_lst_diff['pIv_old'],file_lst_diff['pIv_err_old']*2,linestyle='',marker='.',color='k')#,label='LST Null Test')
#p.errorbar(file_lst_1['kpl']-0.005,file_lst_1['pIv_old'],file_lst_1['pIv_err_old']*2,linestyle='',marker='.',color='k',label='LST1')
#p.errorbar(file_lst_2['kpl']+0.005,file_lst_2['pIv_old'],file_lst_2['pIv_err_old']*2,linestyle='',marker='.',color='0.5',label='LST2')
p.axhline(0,color='k',linestyle='--')
p.axhspan(-noise_lst*2,noise_lst*2,facecolor='0.5',edgecolor="none",alpha=0.5,label='Estimated $2\sigma$ Error')
#p.legend(numpoints=1,prop={'size':16})
p.ylabel('$P(k)$ $[mK^{2}(h^{-1} Mpc)^{3}]$',fontsize=18)
p.xlabel('$k_{\\parallel}$ [$h$ Mpc$^{-1}$]',fontsize=18)
p.tick_params(axis='both', which='major', labelsize=16)
p.ylim(-1.2e9,1.2e9)
p.grid()
p.title('Null Test: LST',fontsize=18)

# Baselines null test
noise_bl = 8873534.73644
p.figure()
p.errorbar(file_bls_diff['kpl'],file_bls_diff['pIv_old'],file_lst_diff['pIv_err_old']*2,linestyle='',marker='.',color='k')#,label='Baselines Null Test')
#p.errorbar(file_bls_1['kpl']-0.005,file_bls_1['pIv_old'],file_bls_1['pIv_err_old']*2,linestyle='',marker='.',color='k',label='BL1')
#p.errorbar(file_bls_2['kpl']+0.005,file_bls_2['pIv_old'],file_bls_2['pIv_err_old']*2,linestyle='',marker='.',color='0.5',label='BL2')
p.axhline(0,color='k',linestyle='--')
p.axhspan(-noise_bl*2,noise_bl*2,facecolor='0.5',edgecolor="none",alpha=0.5,label='Estimated $2\sigma$ Error')
#p.legend(numpoints=1,prop={'size':16})
#p.ylabel('$P(k)$ $[mK^{2}(h^{-1} Mpc)^{3}]$',fontsize=18)
p.xlabel('$k_{\\parallel}$ [$h$ Mpc$^{-1}$]',fontsize=18)
p.tick_params(axis='both', which='major', labelsize=16)
p.ylim(-1.2e9,1.2e9)
p.grid()
p.title('Null Test: Baselines',fontsize=18)

# Even/Odd null test
noise_eo = 4436767.36822*2 # sensitivity of usual even/odd multiplied by 2, since now we're only doing one quadrant of the grid instead of 2
p.figure()
p.errorbar(file_eo_diff['kpl'],file_eo_diff['pIv_old'],file_eo_diff['pIv_err_old']*2,linestyle='',marker='.',color='k')#,label='Even/Odd Null Test')
#p.errorbar(file_eo_1['kpl']-0.005,file_eo_1['pIv_old'],file_eo_1['pIv_err_old']*2,linestyle='',marker='.',color='k',label='Even')
#p.errorbar(file_eo_2['kpl']+0.005,file_eo_2['pIv_old'],file_eo_2['pIv_err_old']*2,linestyle='',marker='.',color='0.5',label='Odd')
p.axhline(0,color='k',linestyle='--')
p.axhspan(-noise_eo*2,noise_eo*2,facecolor='0.5',edgecolor="none",alpha=0.5,label='Estimated $2\sigma$ Error')
#p.legend(numpoints=1,prop={'size':16})
#p.ylabel('$P(k)$ $[mK^{2}(h^{-1} Mpc)^{3}]$',fontsize=18)
p.xlabel('$k_{\\parallel}$ [$h$ Mpc$^{-1}$]',fontsize=18)
p.tick_params(axis='both', which='major', labelsize=16)
p.ylim(-1.2e9,1.2e9)
p.grid()
p.title('Null Test: Even/Odd',fontsize=18)

p.show()

