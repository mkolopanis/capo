#! /usr/bin/env python
"""Compute Signal Loss factor via distribution method v3.
    -Reads in pspec_2d_to_1d.npz (contains PS points for all bootstraps)
    -Finds distribution of Pin vs. Pout curve (using absolute values of points and 2D kernel density estimator) = M
    -Finds distribution of PS points = p (using 1D KDE)
    -Multilpies p element-wise for every column of M (for every Pin)
    -Sums up each row of M, yielding one final distribution along Pin
    -Computes power spectrum points and errors from the distribution
"""
import numpy as n
import pylab as p
from matplotlib import gridspec
from scipy.interpolate import interp1d
from capo.pspec import f2z
import scipy
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
   
        for ind in range(len(kpl)): # loop through k for P(k)
                pCs[kpl[ind]] = [pC[:,ind]] # no appending because it's the same for every inject
                pIs[kpl[ind]] = [pI[:,ind]]
    

    # Values of interest are contained within dictionaries indexed by kpl_fold:
    #     Pins_fold
    #     Pouts_fold
    #     Pouts_I_fold
    #     pCs, pCs_fold, pIs, pIs_fold  (data/noise values)


    # Plot Pin vs. Pout curves 
    for ind,k in enumerate(Pins_fold.keys()):
       
        if opts.plot:
            p.figure(1) # Pin vs. Pout
            p.subplot(3, 4, ind+1)
            p.plot(Pins_fold[k], Pouts_fold[k], color)  # points
            p.plot(Pins_fold[k], poly[k], 'r-') # polyfit
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
            p.plot(Pins_fold[k], Pouts_I_fold[k], color)  # points
            p.plot(Pins_fold[k], poly_I[k], 'r-') # polyfit
            #p.plot([pklo, pkhi], [pklo, pkhi], 'k-')  # diagonal line
            p.plot([-pkhi, pkhi], [-pkhi, pkhi], 'k-')  # diagonal line
            p.grid(True)
            #p.xlim(pklo, pkhi)
            #p.ylim(pklo, pkhi)
            p.xscale('symlog',linthreshx=1e6)
            p.yscale('symlog',linthreshy=1e6)
            p.tick_params(axis='both', which='both', labelsize=6)
            p.title('k = '+str("%.4f" % kpl_fold[ind]), fontsize=8)

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
        p.show()

    ### FUNCTIONS

    # Get bin sizes to use for normalization
    def bin_size(bins): # XXX approximation since they're just divided in two
        bins_size = []
        for bb in range(len(bins)):
            if bb == 0: bins_size.append(bins[bb+1]-bins[bb]) # first bin
            elif bb == len(bins)-1: bins_size.append(bins[bb]-bins[bb-1]) # last bin
            else: bins_size.append((bins[bb+1]-bins[bb-1])/2) # middle bins
        return bins_size

    # Make bins (grid), which is used for sampling distributions 
    def make_bins():
        # lin-log grid-spacing
        nbins = 101 # XXX hard-coded number of bins for grid upon which to estimate kernels (must be odd for fftshift to do the right thing)
        dmax = n.max(Pins_fold.values())
        dg = 1 # spacing in linear regime
        G = n.logspace(0,n.log10(dmax),nbins)[-1]/n.logspace(0,n.log10(dmax),nbins)[-2] # way to estimate multiplicative factor in log regime (~1.3) 
        grid = [] # grid spacing that's linear at low regime and log at high regime 
        g = 0
        for i in range(nbins):
            g = G*g + dg
            grid.append(g)
        bins = n.log10(grid) # log-space
        bins_size = bin_size(bins)
        bins_concat = n.concatenate((-10**bins[::-1],10**bins)) # full length, real numbers
        return bins, bins_concat # bins is n.log10, bins_concat is not

    # Fit polynomial to signal loss curve and translate it into a P_in vs. P_out matrix where there's one P_out value for every P_in
    def curve_to_matrix(x,y):
        m = n.zeros((len(bins),len(bins))) # matrix that will be populated
        xs = n.log10(n.abs(x)) # absolute value since symmetric
        ys = n.log10(n.abs(y))
        xs = n.append(n.repeat(0,100000),xs) # force fit to go through zero
        ys = n.append(n.repeat(0,100000),ys)
        order = n.argsort(xs) # re-order after padding
        xs = xs[order]
        ys = ys[order]
        # find polyfit
        coeff = n.polyfit(xs,ys,8) # coefficients from highest to lowest order
        for bb,b in enumerate(bins): # walk through P_in
            y_bin = n.interp(b,xs,n.polyval(coeff,xs)) # get P_out bin for given P_in bin
            y_ind = n.argmin(n.abs(bins-y_bin)) # get P_out index
            m[y_ind][bb] = 1.0 # fill with 1
        return m

    def make_gaussian(sigma,offset):
        x = n.fft.fftshift(n.arange(-bins.size/2,bins.size/2))
        g = n.exp(-x**2/(2*sigma**2))
        g /= n.sum(g)
        g = n.fft.ifft(n.fft.fft(g)*n.exp(-2j*n.pi*x*offset))
        return g

    # Function to fit for a sigma and offset (convolves "I" curve with a gaussian of some sigma to "fatten" it up, and shifts it to best match "C" curve)
    def find_sigma((sigma,offset),curveC,curveI):
        g = make_gaussian(sigma,offset)
        fat_I = n.convolve(g,curveI,mode='same')
        score = n.sum(n.abs(curveC-fat_I)**2) 
        return score # minimize difference

    # Using all bootstrapped points, use kernel density estimators to smooth out both the C and I distributions, and then de-convolve out the extra width that C has on I for every P_in 
    def convolve_kernel():
        M_matrix_C = {}
        M_matrix_I = {}
        convolve_curve = {}
        for kk,k in enumerate(kpl_fold): # only positive k's
            xs = n.abs(n.array(Pins_fold[k]).flatten())
            ys_C = n.abs(n.array(Pouts_fold[k]).flatten())
            ys_I = n.abs(n.array(Pouts_I_fold[k]).flatten())
            # Kernel Density Estimator
            ygrid,xgrid = n.meshgrid(bins,bins) # create grid on which to sample
            positions = n.vstack([xgrid.ravel(),ygrid.ravel()])
            kernel_C = scipy.stats.gaussian_kde((n.log10(xs),n.log10(ys_C)),bw_method='scott')
            kernel_I = scipy.stats.gaussian_kde((n.log10(xs),n.log10(ys_I)),bw_method=kernel_C.factor)
            M_matrix_C[k] = n.reshape(kernel_C(positions).T,xgrid.shape).T
            M_matrix_I[k] = n.reshape(kernel_I(positions).T,xgrid.shape).T
            # ensure columns sum to 1
            for col in range(M_matrix_C[k].shape[1]):
                if n.sum(M_matrix_C[k][:,col]) > 0: # avoid nan values
                    M_matrix_C[k][:,col] /= n.sum(M_matrix_C[k][:,col]*bin_size(bins))
                if n.sum(M_matrix_I[k][:,col]) > 0:
                    M_matrix_I[k][:,col] /= n.sum(M_matrix_I[k][:,col]*bin_size(bins))
            # find convolution
            convolve_curve[k] = n.zeros_like(M_matrix_C[k])
            for col in range(M_matrix_C[k].shape[1]): # doing this for every column of P_in seems overkill and makes the code slow
                curveC = M_matrix_C[k][:,col]
                curveI = M_matrix_I[k][:,col]
                result = scipy.optimize.least_squares(find_sigma,(0,0),bounds=([0,-n.inf],[n.inf,n.inf]),args=(curveC,curveI)) # minimize find_sigma function
                convolve_curve[k][:,col] = n.fft.fftshift(make_gaussian(result.x[0],0))
            # Plot 2D distribution (KDE)
            if False and kk == 0: # just one k-value
                p.plot(n.log10(xs), n.log10(ys),'k.')
                p.pcolormesh(bins_in,bins_out,M_matrix_C[k]);p.colorbar()
                p.title('k='+str(k))
                p.xlabel('$P_{in}$ (log)');p.ylabel('$P_{out}$ (log)');p.show()
        return convolve_curve

    # Average all bootstrapped points together to get one signal loss curve
    # Convolve it based on the extra width "C" has on "I"
    # Return final M matrix (transfer function)
    def transfer_func(identity=False):
        M_matrix = {}
        for kk,k in enumerate(kpl_fold): # only positive k's
            # Average signal loss curves together
            xs = n.mean(n.array(Pins_fold[k]),axis=1)
            if identity == True: ys = n.mean(n.array(Pouts_I_fold[k]),axis=1)
            if identity == False: ys = n.mean(n.array(Pouts_fold[k]),axis=1)
            M_matrix[k] = curve_to_matrix(xs,ys)
            if identity == False: # do convolution here
                for col in range(M_matrix[k].shape[1]):
                    if n.sum(convolve_curve[k][:,col]) > 0: 
                        new_col = n.convolve(M_matrix[k][:,col],convolve_curve[k][:,col],mode='same')
                        M_matrix[k][:,col] = new_col/n.sum(new_col*bin_size(bins))
                    else: M_matrix[k][:,col] /= n.sum(M_matrix[k][:,col]*bin_size(bins))
        return M_matrix
    
    # Get data distribution
    def data_dist(data):
        old_PS = {}
        for k in data.keys():
            kernel = scipy.stats.gaussian_kde(data[k][0]) # KDE for data
            data_dist = kernel(bins_concat)
            bins_concat_size = bin_size(bins_concat)
            data_dist /= n.sum(data_dist*bins_concat_size) # normalize
            old_PS[k] = data_dist
        return old_PS

    # Apply transfer function to data
    def sigloss_func(data_dist, M_matrix):
        new_PS = {} # dictionaries that will hold PS distributions
        for k in data_dist.keys():
            data_dist_pos = data_dist[k][len(bins):] # separate
            data_dist_neg = data_dist[k][:len(bins)][::-1] # reverse this to apply same M
            # Get new distribution
            try: M = M_matrix[k]
            except: M = M_matrix[-k] # M only exists for positive k's
            
            Mpos = M*(n.resize(data_dist_pos,M.shape).T) # multiply data distribution element-wise per column of M 
            Mneg = M*(n.resize(data_dist_neg,M.shape).T) 
            rowsum_pos = n.zeros(M.shape[0])
            rowsum_neg = n.zeros(M.shape[0])
            for row in Mpos: rowsum_pos += row # add up rows of M
            for row in Mneg: rowsum_neg += row 
            rowsum_combine = n.concatenate((rowsum_neg[::-1],rowsum_pos))
            new_PS[k] = rowsum_combine/n.sum(rowsum_combine*bin_size(bins_concat)) # normalize
        return new_PS # distribution of bins_concat
   

    ### SIGNAL LOSS CODE

    bins, bins_concat = make_bins() # bins are where to sample distributions
   
    # Get distributions of original data 
    old_pCs = data_dist(pCs)
    old_pCs_fold = data_dist(pCs_fold)
    old_pIs = data_dist(pIs)    
    old_pIs_fold = data_dist(pIs_fold)
    
    if opts.skip_sigloss:
        print "Skipping Signal Loss!"
        new_pCs = old_pCs.copy()
        new_pCs_fold = old_pCs_fold.copy()
        new_pIs = old_pIs.copy()
        new_pIs_fold = old_pIs_fold.copy()
    else:
        convolve_curve = convolve_kernel() # convolution curve to convolve the "C" signal loss curve by, to account for the extra width that "C" can have on "I"
        M = transfer_func(identity=False) # final signal loss transfer function
        M_I = transfer_func(identity=True)
        new_pCs = sigloss_func(old_pCs, M) # distributions of bins_concat, before and after signal loss correction
        new_pCs_fold = sigloss_func(old_pCs_fold, M)
        new_pIs = sigloss_func(old_pIs, M_I)
        new_pIs_fold = sigloss_func(old_pIs_fold, M_I)
    
    # Compute PS points and errors
    def compute_stats(bins,data):
        pts, errs = [], []
        for key in n.sort(data.keys()):
            point = bins[n.argmax(data[key])]
            pts.append(point)
            percents = [n.sum((data[key]*bin_size(bins_concat))[:i]) for i in range(len(data[key]))]
            left = n.interp(0.025,percents,bins)
            right = n.interp(0.975,percents,bins)
            errs.append((right-left)/2) # 1-sigma
        return pts,errs
    
    pC, pC_err = compute_stats(bins_concat, new_pCs)
    pI, pI_err = compute_stats(bins_concat, new_pIs)
    pC_fold, pC_fold_err = compute_stats(bins_concat, new_pCs_fold)
    pI_fold, pI_fold_err = compute_stats(bins_concat, new_pIs_fold)
    
    # Plot old vs. new distributions for P(k)
    if False:
        num = 1 # weighted
        p.figure(figsize=(15,8))
        for nn in range(len(file['kpl'])):
            p.subplot(3,7,num)
            p.plot(bin_cent,old_pCs[file['kpl'][nn]]/n.max(old_pCs[file['kpl'][nn]]),'b-',label='old')
            p.plot(bin_cent,new_pCs[file['kpl'][nn]]/n.max(new_pCs[file['kpl'][nn]]),'g-',label='new')
            p.xscale('symlog')
            p.title('k = ' + str(file['kpl'][nn]),fontsize=8)
            p.tick_params(axis='both', which='major', labelsize=6)
            p.tick_params(axis='both', which='minor', labelsize=6)
            p.grid()
            num += 1
        p.suptitle('Distribution of Weighted PS, Old (blue) vs. New (green)')
        p.show()
        num = 1 # unweighted
        p.figure(figsize=(15,8))
        for nn in range(len(file['kpl'])):
            p.subplot(3,7,num)
            p.plot(bin_cent_I,old_pIs[file['kpl'][nn]]/n.max(old_pIs[file['kpl'][nn]]),'b-',label='old')
            p.plot(bin_cent_I,new_pIs[file['kpl'][nn]]/n.max(new_pIs[file['kpl'][nn]]),'g-',label='new')
            p.xscale('symlog')
            p.title('k = ' + str(file['kpl'][nn]),fontsize=8)
            p.tick_params(axis='both', which='major', labelsize=6)
            p.tick_params(axis='both', which='minor', labelsize=6)
            p.grid()
            num += 1
        p.suptitle('Distribution of Unweighted PS, Old (blue) vs. New (green)')
        p.show()

    # Save values
    #if opts.skip_sigloss:
    #    print 'Skipping applying signal loss!!!'
    #    pC,pC_err = compute_stats(bins_concat,old_pCs)
    #    pI,pI_err = compute_stats(bins_concat,old_pIs)
    #    pC_fold,pC_fold_err = compute_stats(bins_concat,old_pCs_fold)
    #    pI_fold,pI_fold_err = compute_stats(bins_concat,old_pIs_fold)
        
    if count == 0: # data case
        pCv = pC    
        pCv_fold = pC_fold
        pIv = pI
        pIv_fold = pI_fold
        pCv_err = pC_err
        pCv_fold_err = pC_fold_err
        pIv_err = pI_err
        pIv_fold_err = pI_fold_err
    if count == 1: # noise case
        pCn = pC
        pCn_fold = pC_fold
        pIn = pI
        pIn_fold = pI_fold
        pCn_err = pC_err
        pCn_fold_err = pC_fold_err
        pIn_err = pI_err
        pIn_fold_err = pI_fold_err

# Write out solutions
outname = 'pspec_final_sep'+opts.sep+'.npz'
print '   Saving', outname

# Save with 1-sigma errors (compatible with pspec_combine.py and plot_pspec_final.py)
n.savez(outname, kpl=kpl, k=file['k'], freq=file['freq'],
        pCv=pCv, pCv_err=pCv_err,
        pCv_fold=pCv_fold, pCv_fold_err=pCv_fold_err,
        pIv=pIv, pIv_err=pIv_err,
        pIv_fold=pIv_fold, pIv_fold_err=pIv_fold_err,
        pCn=pCn, pCn_err=pCn_err,
        pCn_fold=pCn_fold, pCn_fold_err=pCn_fold_err,
        pIn=pIn, pIn_err=pIn_err,
        pIn_fold=pIn_fold, pIn_fold_err=pIn_fold_err,
        prob=0.95,
        ngps=file['ngps'], nbls=file['nbls'], nbls_g=file['nbls_g'], nlsts_g=file['nlsts_g'],
        lsts=file['lsts'], afreqs=file['afreqs'], cnt_eff=file['cnt_eff'],
        frf_inttime=file['frf_inttime'], inttime=file['inttime'], 
        cmd=file['cmd'].item() + ' \n '+' '.join(sys.argv))

