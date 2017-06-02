import numpy as np, matplotlib.pyplot as plt
import sys, optparse
o = optparse.OptionParser()

o.add_option('--path2data',default='/Users/saulkohn/Documents/PAPER/powerspectra/nBlGroupExperiment/',help='Path to the top of the directory tree containing the power spectra')
o.add_option('--nblg_min',default=2,help='Minimum number of baseline groups used')
o.add_option('--nblg_max',default=23,help='Maximum number of baseline groups used')
o.add_option('--nlst_min',default=1,help='Minimum number of LST groups used')
o.add_option('--nlst_max',default=10,help='Maximum number of LST groups used')
o.add_option('--fold',action='store_true',help='Fold k and form Delta^2 power spectra')

opts, args = o.parse_args(sys.argv[1:])

# Parse options
path2data = opts.path2data
nblg_min,nblg_max = int(opts.nblg_min),int(opts.nblg_max)
nlst_min,nlst_max = int(opts.nlst_min),int(opts.nlst_max)
nblg_L = (nblg_max - nblg_min) + 1
nlst_L = (nlst_max - nlst_min) + 1

if not opts.fold:
    subs = ''
    nice = [0,1,2,3,4,16,17,19,20]
    Nk = 21
else:
    subs = '_fold'
    nice = [5,6,7,8,9,10]
    Nk = 11

# Data storage
nbl_g,nlst_g = np.zeros((nblg_L)),np.zeros((nlst_L))
pIv_bl,pIn_bl = np.zeros((nblg_L,Nk)),np.zeros((nblg_L,Nk)) #21 k_para bins
pIv_lst,pIn_lst = np.zeros((nlst_L,Nk)),np.zeros((nlst_L,Nk))
pCv_bl,pCv_lst = np.zeros((nblg_L,Nk)),np.zeros((nlst_L,Nk))

pIv_bl,pIn_bl = np.zeros((nblg_L,Nk)),np.zeros((nblg_L,Nk)) #21 k_para bins
pIv_lst,pIn_lst = np.zeros((nlst_L,Nk)),np.zeros((nlst_L,Nk))
pCv_bl,pCv_lst = np.zeros((nblg_L,Nk)),np.zeros((nlst_L,Nk))

# Data Parsing

for i,nBL in enumerate(range(nblg_min,nblg_max+1)):
    d = np.load('%s/nblg%i/inject_sep0,1_0.01/pspec_pk_k3pk.npz'%(path2data,nBL))
    pIv_bl[i,:] = d['pIv%s_err'%subs]
    pIn_bl[i,:] = d['pIn%s_err'%subs]
    pCv_bl[i,:] = d['pCv%s_err'%subs]
    nbl_g[i] = d['nbls_g']
    if i==0:
        for j,nLST in enumerate(range(nlst_min,nlst_max+1)):
            dl = np.load('%s/nblg%i/inject_sep0,1_0.01/nlst%i/pspec_pk_k3pk.npz'%(path2data,nBL,nLST))
            pIv_lst[j,:] = dl['pIv%s_err'%subs]
            pIn_lst[j,:] = dl['pIn%s_err'%subs]
            pCv_lst[j,:] = dl['pCv%s_err'%subs]
            nlst_g[j] = dl['nlsts_g']

kpl = d['kpl']

pIn_bl_avg = np.mean(pIn_bl[:,nice],axis=1)
pIv_bl_avg = np.mean(pIv_bl[:,nice],axis=1)
pCv_bl_avg = np.mean(pCv_bl[:,nice],axis=1)
pIn_lst_avg = np.mean(pIn_lst[:,nice],axis=1)
pIv_lst_avg = np.mean(pIv_lst[:,nice],axis=1)
pCv_lst_avg = np.mean(pCv_lst[:,nice],axis=1)


f,axarr = plt.subplots(1,2,sharey=True)

for nk in nice:
    axarr[0].plot(nbl_g,pIn_bl[:,nk],'g-',alpha=0.1)
    axarr[0].plot(nbl_g,pIv_bl[:,nk],'b-',alpha=0.1)
    axarr[0].plot(nbl_g,pCv_bl[:,nk],'r-',alpha=0.1)

    axarr[1].plot(nlst_g,pIn_lst[:,nk],'g-',alpha=0.1)
    axarr[1].plot(nlst_g,pIv_lst[:,nk],'b-',alpha=0.1)
    axarr[1].plot(nlst_g,pCv_lst[:,nk],'r-',alpha=0.1)
    
    print np.log10(pIv_lst[:,nk])-np.log10(pIn_lst[:,nk])
    
axarr[0].plot(nbl_g,pIn_bl_avg,'g-',lw=2)
axarr[0].plot(nbl_g,pIv_bl_avg,'b-',lw=2)
axarr[0].plot(nbl_g,pCv_bl_avg,'r-',lw=2)

axarr[1].plot(nlst_g,pIn_lst_avg,'g-',lw=2,label='noise')
axarr[1].plot(nlst_g,pIv_lst_avg,'b-',lw=2,label='data+I')
axarr[1].plot(nlst_g,pCv_lst_avg,'r-',lw=2,label='data+C')

axarr[0].set_xlabel('nbl/group')
axarr[1].set_xlabel('nLST/group')
if not opts.fold: ylab = 'P(k)'
else: ylab = r'$\Delta^2(k)$'
axarr[0].set_ylabel(ylab)
axarr[0].set_yscale('log')
plt.legend()

plt.show()
