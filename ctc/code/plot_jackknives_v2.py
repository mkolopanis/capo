#! /usr/bin/env python

import numpy as n
import pylab as p

"""
Error bars are errors of each PS half added in quadrature
"""

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

# Mask for masking out inner horizon points
mask = n.ones_like(file_lst_diff['kpl'])
mask[9:12] = 0
mask = n.ma.masked_where(mask==0, mask)

# Horizon limit
horizon = 0.06272882578029243 # computed using hera_pspec.conversions.Cosmo_Conversions.tau_to_kpara

#------------------------

# LST null test
noise_lst = n.sqrt(2*7979656.32115**2) # 1-sigma noise added in quadrature (sqrt(n^2 + n^2))
p.figure()
p.errorbar(file_lst_diff['kpl']*mask,file_lst_diff['pIv_old']*mask,n.sqrt(file_lst_1['pIv_err_old']**2 + file_lst_2['pIv_err_old']**2)*2*mask,linestyle='',marker='.',color='k')#,label='LST Null Test')
#p.errorbar(file_lst_1['kpl']-0.005,file_lst_1['pIv_old'],file_lst_1['pIv_err_old']*2,linestyle='',marker='.',color='k',label='LST1')
#p.errorbar(file_lst_2['kpl']+0.005,file_lst_2['pIv_old'],file_lst_2['pIv_err_old']*2,linestyle='',marker='.',color='0.5',label='LST2')
p.axvline(x=horizon,color='k',linestyle='--')
p.axvline(x=-horizon,color='k',linestyle='--')
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
noise_bl = n.sqrt(2*8873534.73644**2)
p.figure()
p.errorbar(file_bls_diff['kpl']*mask,file_bls_diff['pIv_old']*mask,n.sqrt(file_bls_1['pIv_err_old']**2 + file_bls_2['pIv_err_old']**2)*2*mask,linestyle='',marker='.',color='k')#,label='Baselines Null Test')
#p.errorbar(file_bls_1['kpl']-0.005,file_bls_1['pIv_old'],file_bls_1['pIv_err_old']*2,linestyle='',marker='.',color='k',label='BL1')
#p.errorbar(file_bls_2['kpl']+0.005,file_bls_2['pIv_old'],file_bls_2['pIv_err_old']*2,linestyle='',marker='.',color='0.5',label='BL2')
p.axvline(x=horizon,color='k',linestyle='--')
p.axvline(x=-horizon,color='k',linestyle='--')
p.axhline(0,color='k',linestyle='--')
p.axhspan(-noise_bl*2,noise_bl*2,facecolor='0.5',edgecolor="none",alpha=0.5,label='Estimated $2\sigma$ Error')
#p.legend(numpoints=1,prop={'size':16})
p.ylabel('$P(k)$ $[mK^{2}(h^{-1} Mpc)^{3}]$',fontsize=18)
p.xlabel('$k_{\\parallel}$ [$h$ Mpc$^{-1}$]',fontsize=18)
p.tick_params(axis='both', which='major', labelsize=16)
p.ylim(-1.2e9,1.2e9)
p.grid()
p.title('Null Test: Baselines',fontsize=18)

# Even/Odd null test
noise_eo = n.sqrt(2*(4436767.36822*2)**2) # sensitivity of usual even/odd multiplied by 2, since now we're only doing one quadrant of the grid instead of 2
p.figure()
p.errorbar(file_eo_diff['kpl']*mask,file_eo_diff['pIv_old']*mask,n.sqrt(file_eo_1['pIv_err_old']**2 + file_eo_2['pIv_err_old']**2)*2*mask,linestyle='',marker='.',color='k')#,label='Even/Odd Null Test')
#p.errorbar(file_eo_1['kpl']-0.005,file_eo_1['pIv_old'],file_eo_1['pIv_err_old']*2,linestyle='',marker='.',color='k',label='Even')
#p.errorbar(file_eo_2['kpl']+0.005,file_eo_2['pIv_old'],file_eo_2['pIv_err_old']*2,linestyle='',marker='.',color='0.5',label='Odd')
p.axvline(x=horizon,color='k',linestyle='--')
p.axvline(x=-horizon,color='k',linestyle='--')
p.axhline(0,color='k',linestyle='--')
p.axhspan(-noise_eo*2,noise_eo*2,facecolor='0.5',edgecolor="none",alpha=0.5,label='Estimated $2\sigma$ Error')
#p.legend(numpoints=1,prop={'size':16})
p.ylabel('$P(k)$ $[mK^{2}(h^{-1} Mpc)^{3}]$',fontsize=18)
p.xlabel('$k_{\\parallel}$ [$h$ Mpc$^{-1}$]',fontsize=18)
p.tick_params(axis='both', which='major', labelsize=16)
p.ylim(-1.2e9,1.2e9)
p.grid()
p.title('Null Test: Even/Odd',fontsize=18)

p.show()

