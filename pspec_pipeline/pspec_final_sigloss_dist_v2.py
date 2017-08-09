#! /usr/bin/env python
"""Compute Signal Loss factor via distribution method v2.
   -Begins with distribution of power spectrum values p (all bootstraps)
   -Computes Pin vs. Pout curve
   -Bins in Pin and Pout
   -Computes probabilities in each bin 
   -Matches them up with probabilities of p
   -End up with rescaled p based Pin values allowed for a given Pout
"""
import numpy as n
import pylab as p
from matplotlib import gridspec
from scipy.interpolate import interp1d
from capo.pspec import f2z
import glob
import optparse
import sys

o = optparse.OptionParser()
o.set_description(__doc__)
o.add_option('--plot',action='store_true',
            help='Show plots.')
o.add_option('--skip_sigloss',action='store_true',
            help='Save values without correcting for signal loss.')
o.add_option('--sep',default='0,2',
            help='Separation to run script for.')
opts,args = o.parse_args(sys.argv[1:])

# read one pspec_pk_k3pk.npz file for k points
file = n.load(glob.glob('inject_sep'+opts.sep+'*')[0]+'/pspec_pk_k3pk.npz')
kpl = file['kpl']; k = file['k']; kpl_fold = file['kpl_fold']

# loop over injects and read pspec_2d_to_1d outputs
for count in range(2):
    if count == 1: print 'NOISE CASE'
    else: print 'DATA CASE'
    if opts.plot:
        fig1 = p.figure(1, figsize=(15, 7))
        fig2 = p.figure(2, figsize=(15, 7))
    pklo, pkhi = 1e-1, 1e12
    Pouts_fold = {}
    Pouts_I_fold = {}
    Pins_fold = {}
    pCs = {}; pCs_fold = {}
    pIs = {}; pIs_fold = {}
    for inject in glob.glob('inject_sep'+opts.sep+'*'):
        print 'Reading', inject
        file_2d = n.load(inject + '/pspec_2d_to_1d.npz')
        if count == 1: # noise case
            color='b.'
            linecolor='cyan'
            Pout_fold = file_2d['pCs-pCn_fold']
            Pout_I_fold = file_2d['pIs-pIn_fold']
            pC = file_2d['pCn']
            pC_fold = file_2d['pCn_fold']
            pI = file_2d['pIn']
            pI_fold = file_2d['pIn_fold']
        else:
            color='k.'
            linecolor='0.5'
            Pout_fold = file_2d['pCr-pCv_fold'] # shape (#boots, #kpl)
            Pout_I_fold = file_2d['pIr-pIv_fold']  # pI case
            pC = file_2d['pCv']
            pC_fold = file_2d['pCv_fold']
            pI = file_2d['pIv']
            pI_fold = file_2d['pIv_fold']

        Pin_fold = file_2d['pIe_fold']
        for ind in range(len(kpl_fold)): # loop through k for Delta^2(k)
            try:
                Pouts_fold[kpl_fold[ind]].append(Pout_fold[:,ind])
                Pouts_I_fold[kpl_fold[ind]].append(Pout_I_fold[:,ind])
                Pins_fold[kpl_fold[ind]].append(Pin_fold[:,ind])
                pCs_fold[kpl_fold[ind]] = [pC_fold[:,ind]] # no appending because it's the same for every inject
                pIs_fold[kpl_fold[ind]] = [pI_fold[:,ind]]
            except:
                Pouts_fold[kpl_fold[ind]] = [Pout_fold[:,ind]]
                Pouts_I_fold[kpl_fold[ind]] = [Pout_I_fold[:,ind]]
                Pins_fold[kpl_fold[ind]] = [Pin_fold[:,ind]]
                pCs_fold[kpl_fold[ind]] = [pC_fold[:,ind]]
                pIs_fold[kpl_fold[ind]] = [pI_fold[:,ind]]

            if opts.plot:

                p.figure(1) # Pin vs. Pout
                p.subplot(3, 4, ind+1)
                p.plot(Pin_fold[:,ind], Pout_fold[:,ind], color)  # points
                #p.plot([pklo, pkhi], [pklo, pkhi], 'k-')  # diagonal line
                p.plot([-pkhi, pkhi], [-pkhi, pkhi], 'k-')  # diagonal line
                p.grid(True)
                p.xscale('symlog',linthreshx=1e6)
                p.yscale('symlog',linthreshy=1e6)
                #p.xlim(pklo, pkhi)
                #p.ylim(pklo, pkhi)
                p.tick_params(axis='both', which='both', labelsize=6)
                p.title('k = '+str("%.4f" % kpl_fold[ind]), fontsize=8)

                p.figure(2) # Pin vs. Pout for I case
                p.subplot(3, 4, ind+1)
                p.plot(Pin_fold[:,ind], Pout_I_fold[:,ind], color)  # points
                #p.plot([pklo, pkhi], [pklo, pkhi], 'k-')  # diagonal line
                p.plot([-pkhi, pkhi], [-pkhi, pkhi], 'k-')  # diagonal line
                p.grid(True)
                #p.xlim(pklo, pkhi)
                #p.ylim(pklo, pkhi)
                p.xscale('symlog',linthreshx=1e6)
                p.yscale('symlog',linthreshy=1e6)
                p.tick_params(axis='both', which='both', labelsize=6)
                p.title('k = '+str("%.4f" % kpl_fold[ind]), fontsize=8)
        
        for ind in range(len(kpl)): # loop through k for P(k)
                pCs[kpl[ind]] = [pC[:,ind]] # no appending because it's the same for every inject
                pIs[kpl[ind]] = [pI[:,ind]]

    if opts.plot:

        # plot Pin vs. Pout
        p.figure(1)
        p.tight_layout()
        fig1.text(0.5, 0.0, r'$P_{\rm in}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',
                 ha='center')
        fig1.text(0.0, 0.5, r'$P_{\rm out}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',
                 va='center', rotation='vertical')

        # plot Pin vs. Pout for I case
        p.figure(2)
        p.tight_layout()
        fig2.text(0.5, 0.0, r'$P_{\rm in}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',
                 ha='center')
        fig2.text(0.0, 0.5, r'$P_{\rm out,I}(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',
                 va='center', rotation='vertical')

    if opts.plot: p.show()

    # Values of interest are contained within dictionaries indexed by kpl_fold:
    #     Pins_fold
    #     Pouts_fold
    #     Pouts_I_fold
    #     pCs, pCs_fold, pIs, pIs_fold  (data/noise values)


    # Function to compute probability matrix ('transfer function')
    def transfer_func(data, identity=False):
   
        # Make bins
        ninjects = n.array(Pins_fold[Pins_fold.keys()[0]]).shape[0]*2
        bins_in = []
        bins_in.append(n.logspace(n.log10(pklo), n.log10(pkhi), ninjects)) # positive
        bins_in.append(-(n.logspace(n.log10(pklo), n.log10(pkhi), ninjects))) # negative 
        bins_in = n.array(bins_in).flatten()
        bins_in.sort()
        #bins_out = []
        #max_val = n.abs(n.array(data.values())).max()
        #min_val = n.abs(n.array(data.values())).min()
        #bins_out.append(n.logspace(n.log10(min_val), n.log10(max_val), ninjects))
        #bins_out.append(-(n.logspace(n.log10(min_val), n.log10(max_val), ninjects)))
        #bins_out = n.array(bins_out).flatten()
        #bins_out.sort()
        bins_out = bins_in.copy()        

        # Compute probabilities in each bin for 'transfer function' Pin vs. Pout
        prob_func = {}
        for k in kpl_fold:
            prob_func[k] = {}
            bin_Pins_fold = n.digitize(n.array(Pins_fold[k]).flatten(),bins_in)
            if identity == False:
                bin_Pouts_fold = n.digitize(n.array(Pouts_fold[k]).flatten(),bins_out)
            if identity == True:
                bin_Pouts_fold = n.digitize(n.array(Pouts_I_fold[k]).flatten(),bins_out)
            for bin_y in range(len(bins_out)+1): # loop through Pout
                prob_Pin = n.zeros(bins_in.size+1)
                index_y = n.where(bin_Pouts_fold == bin_y)
                index_x = bin_Pins_fold[index_y] # indices for bins in Pin that have corresponding Pouts in bin_y
                for x in index_x:
                    prob_Pin[x] += 1
                prob_Pin /= len(index_x)
                prob_func[k][bin_y] = prob_Pin # array of Pin probabilities for this Pout bin
        return bins_in, bins_out, prob_func
        
    # Function to apply transfer function
    def sigloss_func(bins_in, bins_out, data, func):
        # Compute probabilities of PS values in bins
        data_prob = {}
        new_PS = {}
        for k in data.keys():
            data_prob[k] = {}
            bin_data = n.digitize(data[k][0], bins_out) # data[k] is an array in an array for some reason
            for bin_y in range(len(bins_out)+1): # loop through Pout
                index = n.where(bin_data == bin_y)[0]
                prob = float(len(index)) / len(bin_data) # probability of data for a Pout bin
                data_prob[k][bin_y] = prob
            # Match probabilities
            weights = n.zeros(bins_in.size+1)
            for b in data_prob[k].keys(): # loop through Pout bins
                try: val = data_prob[k][b] * func[k][b]
                except: val = data_prob[k][b] * func[-k][b] # if data has negative k's
                val[n.isnan(val)] = 0.0
                weights += val
            new_PS[k] = weights
        return new_PS
   
    # Call function for probability matrix
    bins_in, bins_out, prob_func_pCs = transfer_func(pCs)
    bins_in_I, bins_out_I, prob_func_pIs = transfer_func(pIs, identity=True)
    
    # Call function for PS values
    new_pCs = sigloss_func(bins_in, bins_out, pCs, prob_func_pCs)
    new_pCs_fold = sigloss_func(bins_in, bins_out, pCs_fold, prob_func_pCs)
    new_pIs = sigloss_func(bins_in_I, bins_out_I, pIs, prob_func_pIs)
    new_pIs_fold = sigloss_func(bins_in_I, bins_out_I, pIs_fold, prob_func_pIs)
    
    # Plot histograms for each k
    if opts.plot:        
        for ind,k in enumerate(new_pCs_fold.keys()):
            p.figure(3) # Folded Data Histgraom
            p.subplot(3, 4, ind+1)
            p.plot(bins_in,new_pCs_fold[k][:-1],color,linestyle='-')  # bins are upper limits
            p.grid(True)
            p.xscale('symlog',linthreshx=1e5)
            p.xlim(-1e9, 1e9)
            #p.ylim(pklo, pkhi)
            p.ylim(0,1)
            p.tick_params(axis='both', which='both', labelsize=6)
            p.title('k = '+str("%.4f" % kpl_fold[ind]), fontsize=8)

            p.figure(4) # Folded I Data Histgraom
            p.subplot(3, 4, ind+1)
            p.plot(bins_in_I,new_pIs_fold[k][:-1],color,linestyle='-')  
            p.grid(True)
            p.xscale('symlog',linthreshx=1e5)
            p.xlim(-1e9, 1e9)
            #p.ylim(pklo, pkhi)
            p.ylim(0,1)
            p.tick_params(axis='both', which='both', labelsize=6)
            p.title('k = '+str("%.4f" % kpl_fold[ind]), fontsize=8)
       
        p.figure(3); p.tight_layout()
        p.figure(4); p.tight_layout()
        p.show()

    # Compute PS points 
    bin_cent = [(bins_in[i+1]+bins_in[i])/2.0 for i in range(len(bins_in)-1)]
    bin_cent.append(bin_cent[-1] * 10**(n.log10(bin_cent[-1]/bin_cent[-2]))) # right edge bin
    bin_cent.append(bin_cent[0] * 10**(n.log10(bin_cent[0]/bin_cent[1]))) # left edge bin
    bin_cent.sort()
    bin_cent_I = [(bins_in_I[i+1]+bins_in_I[i])/2.0 for i in range(len(bins_in_I)-1)]
    bin_cent_I.append(bin_cent_I[-1] * 10**(n.log10(bin_cent_I[-1]/bin_cent_I[-2]))) # right edge bin
    bin_cent_I.append(bin_cent_I[0] * 10**(n.log10(bin_cent_I[0]/bin_cent_I[1]))) # left edge bin
    bin_cent_I.sort()

    pC = []; pC_fold = []; pC_err = []; pC_fold_err = []
    pI = []; pI_fold = []; pI_err = []; pI_fold_err = []
    for key in n.sort(new_pCs.keys()):
        pC.append(bin_cent[n.argmax(new_pCs[key])])
        pI.append(bin_cent_I[n.argmax(new_pIs[key])])
    for key in n.sort(new_pCs_fold.keys()):
        pC_fold.append(bin_cent[n.argmax(new_pCs_fold[key])])
        pI_fold.append(bin_cent_I[n.argmax(new_pIs_fold[key])])
    
    # Compute errors
    for key in n.sort(new_pCs.keys()):
        percents = [n.sum(new_pCs[key][:i]) for i in range(len(new_pCs[key]))]
        percents_I = [n.sum(new_pIs[key][:i]) for i in range(len(new_pIs[key]))]
        pC_err.append(n.interp(0.95,percents,bin_cent))
        pI_err.append(n.interp(0.95,percents_I,bin_cent_I))
    for key in n.sort(new_pCs_fold.keys()):
        percents = [n.sum(new_pCs_fold[key][:i]) for i in range(len(new_pCs_fold[key]))]
        percents_I = [n.sum(new_pIs_fold[key][:i]) for i in range(len(new_pIs_fold[key]))]
        pC_fold_err.append(n.interp(0.95,percents,bin_cent))
        pI_fold_err.append(n.interp(0.95,percents_I,bin_cent_I))
    
    # Save values
    if opts.skip_sigloss:
        if count == 0: # data case
            pC = file['pCv']
            pC_fold = file['pCv_fold']
            pI = file['pIv']
            pI_fold = file['pIv_fold']
            pC_err = file['pCv_err']*2
            pC_fold_err = file['pCv_fold_err']*2
            pI_err = file['pIv_err']*2
            pI_fold_err = file['pIv_fold_err']*2
        else: # noise case
            pC = file['pCn']
            pC_fold = file['pCn_fold']
            pI = file['pIn']
            pI_fold = file['pIn_fold']
            pC_err = file['pCn_err']*2
            pC_fold_err = file['pCn_fold_err']*2
            pI_err = file['pIn_err']*2
            pI_fold_err = file['pIn_fold_err']*2
        
    fold_factor = file['k']**3/(2*n.pi**2)
    if count == 0: # data case
        pCv = pC    
        pCv_fold = pC_fold*fold_factor
        pIv = pI
        pIv_fold = pI_fold*fold_factor
        pCv_err = pC_err
        pCv_fold_err = pC_fold_err*fold_factor
        pIv_err = pI_err
        pIv_fold_err = pI_fold_err*fold_factor
    if count == 1: # noise case
        pCn = pC
        pCn_fold = pC_fold*fold_factor
        pIn = pI
        pIn_fold = pI_fold*fold_factor
        pCn_err = pC_err
        pCn_fold_err = pC_fold_err*fold_factor
        pIn_err = pI_err
        pIn_fold_err = pI_fold_err*fold_factor

# Write out solutions
outname = 'pspec_final_sep'+opts.sep+'.npz'
print '   Saving', outname # XXX 2-sigma probability is hard-coded
n.savez(outname, kpl=kpl, k=file['k'], freq=file['freq'],
        pC=pCv, pC_up=pCv_err,
        pC_fold=pCv_fold, pC_fold_up=pCv_fold_err,
        pI=pIv, pI_up=pIv_err,
        pI_fold=pIv_fold, pI_fold_up=pIv_fold_err,
        pCn=pCn, pCn_up=pCn_err,
        pCn_fold=pCn_fold, pCn_fold_up=pCn_fold_err,
        pIn=pIn, pIn_up=pIn_err,
        pIn_fold=pIn_fold, pIn_fold_up=pIn_fold_err,
        prob=0.95,
        ngps=file['ngps'], nbls=file['nbls'], nbls_g=file['nbls_g'], nlsts_g=file['nlsts_g'],
        lsts=file['lsts'], afreqs=file['afreqs'], cnt_eff=file['cnt_eff'],
        frf_inttime=file['frf_inttime'], inttime=file['inttime'], 
        cmd=file['cmd'].item() + ' \n '+' '.join(sys.argv))
