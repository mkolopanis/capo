#! /usr/bin/env python
"""Compute Signal Loss factor via distribution method v3.
    -Reads in pspec_2d_to_1d.npz (contains PS points for aln bootstraps)

    -De-convolves out the extra width that M_C can have on M_I by fitting Gaussians to each and calculating a new sigma for the convolution kernel
    -Finds distribution of PS points = p (using 1D KDE)
    -Multiplies p element-wise for every column of M (for every Pin)
    -Sums up each row of M, yielding one final distribution along Pin
    -Computes power spectrum points and errors from the distribution
"""
import numpy as n
import pylab as p
from matplotlib import gridspec
from scipy.interpolate import interp1d
from scipy.stats import norm
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
o.add_option('--nosubtract',action='store_true',
            help='Run for no-subtraction case (P_out = pCr).')
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
    pCes_Cr_fold = {}
    pCvs_Cr_fold = {}
    pCves_fold = {}
    pIves_fold = {}
    pCrs_fold = {}
    pIrs_fold = {}

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
            if opts.nosubtract:
                Pout_fold = file_2d['pCr_fold']
                Pout_I_fold = file_2d['pIr_fold']
            pCe_Cr_fold = file_2d['pCe_Cr_fold']
            pCv_Cr_fold = file_2d['pCv_Cr_fold']
            try: 
                pCve_fold = file_2d['pCve_fold']
                pIve_fold = file_2d['pIve_fold']
            except: 
                pCve_fold = n.zeros_like(pCv_Cr_fold)
                pIve_fold = n.zeros_like(pCv_Cr_fold)
            pCr_fold = file_2d['pCr_fold']
            pIr_fold = file_2d['pIr_fold']
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
                pCes_Cr_fold[kpl_fold[ind]].append(pCe_Cr_fold[:,ind])
                pCvs_Cr_fold[kpl_fold[ind]].append(pCv_Cr_fold[:,ind])
                pCves_fold[kpl_fold[ind]].append(pCve_fold[:,ind])
                pIves_fold[kpl_fold[ind]].append(pIve_fold[:,ind])
                pCrs_fold[kpl_fold[ind]].append(pCr_fold[:,ind])
                pIrs_fold[kpl_fold[ind]].append(pIr_fold[:,ind])
            except:
                Pouts_fold[kpl_fold[ind]] = [Pout_fold[:,ind]]
                Pouts_I_fold[kpl_fold[ind]] = [Pout_I_fold[:,ind]]
                Pins_fold[kpl_fold[ind]] = [Pin_fold[:,ind]]
                pCs_fold[kpl_fold[ind]] = [pC_fold[:,ind]]
                pIs_fold[kpl_fold[ind]] = [pI_fold[:,ind]]
                pCes_Cr_fold[kpl_fold[ind]] = [pCe_Cr_fold[:,ind]]
                pCvs_Cr_fold[kpl_fold[ind]] = [pCv_Cr_fold[:,ind]]
                pCves_fold[kpl_fold[ind]] = [pCve_fold[:,ind]]
                pIves_fold[kpl_fold[ind]] = [pIve_fold[:,ind]]
                pCrs_fold[kpl_fold[ind]] = [pCr_fold[:,ind]]
                pIrs_fold[kpl_fold[ind]] = [pIr_fold[:,ind]]
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
            #p.plot(Pins_fold[k], poly[k], 'r-') # polyfit
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
            #p.plot(Pins_fold[k], poly_I[k], 'r-') # polyfit
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
        dmax = max(n.max(Pins_fold.values()),n.max(pIs.values())) # maximum is determined from whichever is largest: EoR injection or "I" value
        dg = 1 # spacing in linear regime
        G = n.logspace(0,n.log10(dmax),nbins)[-1]/n.logspace(0,n.log10(dmax),nbins)[-2] # way to estimate multiplicative factor in log regime (~1.3) 
        grid = [] # grid spacing that's linear at low regime and log at high regime 
        g = 0
        for i in range(nbins):
            g = G*g + dg
            grid.append(g)
        bins = n.log10(grid) # log-space
        bins_concat = n.concatenate((-10**bins[::-1],10**bins)) # full length, real numbers
        binsx, binsy = bins, bins
        binsy_full = n.concatenate((-binsy[::-1],binsy)) # for pos and neg
        return binsx, binsy, binsy_full, bins_concat # bins is n.log10, bins_concat is not
   
    # Fit polynomial to signal loss curve to translate it into a P_in vs. P_out matrix where there's one P_out value for every P_in
    def curve_to_matrix(x,y):
        m = n.zeros((len(binsx),len(binsy))) # matrix that will be populated
        xs = n.log10(n.abs(x)) # absolute value since symmetric
        ys = n.log10(n.abs(y))
        xs = n.append(n.repeat(0,100000),xs) # force fit to go through zero
        ys = n.append(n.repeat(0,100000),ys) 
        #xs = n.append(0,xs) # force fit to go through zero
        #ys = n.append(0,ys)
        order = n.argsort(xs) # re-order after padding
        xs = xs[order]
        ys = ys[order]
        coeff = n.polyfit(xs,ys,8) # coefficients from highest to lowest order
        for bb,b in enumerate(binsx): # walk through P_in
            if b <= n.max(xs): # only fill out curve if there's data points there
                #y_bin = n.interp(b,xs,ys)
                y_bin = n.interp(b,xs,n.polyval(coeff,xs)) # get P_out bin for given P_in bin
                y_ind = n.argmin(n.abs(binsy-y_bin)) # get P_out index
                m[y_ind][bb] = 1.0 # fill with 1
            else: m[y_ind][bb] = 0.0
        return m

    # Make gaussian with some offset and sigma
    def gauss(x,mu,sigma):
        g = n.exp(-(x-mu)**2/(2*sigma**2))
        g /= n.sum(g*bin_size(x))
        return g

    # Parse data into positive and negative components
    def parse_data(xs, ys_C, ys_I):
        pos_xs = n.where(xs > 0)[0] # positive P_ins only
        xs = xs[pos_xs]
        ys_C = ys_C[pos_xs]
        ys_I = ys_I[pos_xs]
        pos_ys_C = n.where(ys_C > 0)[0]  
        neg_ys_C = n.where(ys_C < 0)[0]
        pos_ys_I = n.where(ys_I > 0)[0]
        neg_ys_I = n.where(ys_I < 0)[0]
        ys_C_pos = ys_C[pos_ys_C]
        ys_C_neg = -ys_C[neg_ys_C]
        ys_I_pos = ys_I[pos_ys_I]
        ys_I_neg = -ys_I[neg_ys_I]
        xs_C_pos = xs[pos_ys_C]
        xs_C_neg = xs[neg_ys_C]
        xs_I_pos = xs[pos_ys_I]
        xs_I_neg = xs[neg_ys_I]
        return xs_C_pos, xs_C_neg, ys_C_pos, ys_C_neg, xs_I_pos, xs_I_neg, ys_I_pos, ys_I_neg

    # For all data points (full distributions), smooth the transfer curves for both the C and I cases, treating positive and negative P_outs separately and combining them together at the end
    def smooth_dist():
        kde_C = {}
        kde_I = {}
        for kk,k in enumerate(kpl_fold): # only positive k's

            if kk == 8 and count == 0: # save values out
                # Save values to use for plotting sigloss terms
                if opts.nosubtract: fn = 'pspec_sigloss_terms_nosubtract.npz'
                else: fn = 'pspec_sigloss_terms.npz'
                print "Saving",fn,"which contains data values for k =",k
                n.savez(fn, k=k, pCv=file['pCv'][n.where(kpl==k)[0][0]], pIv=file['pIv'][n.where(kpl==k)[0][0]], Pins=Pins_fold[k], Pouts=Pouts_fold[k], Pouts_I=Pouts_I_fold[k], pCrs_fold=pCrs_fold[k], pCvs_Cr_fold=pCvs_Cr_fold[k], pCes_Cr_fold=pCes_Cr_fold[k], pCves_fold=pCves_fold[k], pIves_fold=pIves_fold[k])
            xs = n.array(Pins_fold[k]).flatten()
            ys_C = n.array(Pouts_fold[k]).flatten()
            ys_I = n.array(Pouts_I_fold[k]).flatten()
            xs_C_pos,xs_C_neg,ys_C_pos,ys_C_neg,xs_I_pos,xs_I_neg,ys_I_pos,ys_I_neg = parse_data(xs, ys_C, ys_I)
            # Kernel Density Estimator
            ygrid,xgrid = n.meshgrid(binsx,binsy) # create grid on which to sample
            positions = n.vstack([xgrid.ravel(),ygrid.ravel()])
            kernel_C_pos = scipy.stats.gaussian_kde((n.log10(xs_C_pos),n.log10(ys_C_pos)),bw_method='scott')
            factor = kernel_C_pos.factor+0.3 # XXX
            kernel_C_pos = scipy.stats.gaussian_kde((n.log10(xs_C_pos),n.log10(ys_C_pos)),bw_method=factor)
            kernel_I_pos = scipy.stats.gaussian_kde((n.log10(xs_I_pos),n.log10(ys_I_pos)),bw_method=factor)
            if len(xs_C_neg) > 0: # if there are negative points
                kernel_C_neg = scipy.stats.gaussian_kde((n.log10(xs_C_neg),n.log10(ys_C_neg)),bw_method=factor)
                neg_kern = n.reshape(kernel_C_neg(positions).T,(binsx.size,binsy.size)).T[::-1]
            else: # no negative values at all
                neg_kern = n.zeros((binsx.size,binsy.size))
            if n.inf in neg_kern: neg_kern = n.zeros((binsx.size,binsy.size)) # sometimes KDE works but returns inf values if there's too few points
            if len(xs_I_neg) > 0: 
                kernel_I_neg = scipy.stats.gaussian_kde((n.log10(xs_I_neg),n.log10(ys_I_neg)),bw_method=factor)
                neg_kern_I = n.reshape(kernel_I_neg(positions).T,(binsx.size,binsy.size)).T[::-1]
            else: # no negative values at all
                neg_kern_I = n.zeros((binsx.size,binsy.size))
            if n.inf in neg_kern_I: neg_kern_I = n.zeors((binsx.size,binsy.size))
            kde_C[k] = n.concatenate((neg_kern,n.reshape(kernel_C_pos(positions).T,(binsx.size,binsy.size)).T))
            kde_I[k] = n.concatenate((neg_kern_I,n.reshape(kernel_I_pos(positions).T,(binsx.size,binsy.size)).T))
            # ensure columns sum to 1
            for col in range(kde_C[k].shape[1]):
                if n.sum(kde_C[k][:,col]) > 0: # avoid nan values
                    kde_C[k][:,col] /= n.sum(kde_C[k][:,col]*bin_size(binsy_full))
                if n.sum(kde_I[k][:,col]) > 0:
                    kde_I[k][:,col] /= n.sum(kde_I[k][:,col]*bin_size(binsy_full))
        return kde_C, kde_I

    # De-convolve out the extra width that C has on I for every P_in (treat pos and neg halves separately)
    def convolve_kernel(kde_C,kde_I):
        convolve_pos, convolve_neg = {},{}
        for k in kde_C.keys():
            convolve_pos[k] = n.zeros(shape=(binsx.size,binsy.size))
            convolve_neg[k] = n.zeros(shape=(binsx.size,binsy.size))
            neg_MC = kde_C[k][:binsy.size]
            pos_MC = kde_C[k][binsy.size:]
            neg_MI = kde_I[k][:binsy.size]
            pos_MI = kde_I[k][binsy.size:]
            for col in range(kde_C[k].shape[1]): 
                # take a column cut and normalize correctly, fit gaussian
                # Positive half
                curveC_i_pos = pos_MC[:,col] / n.sum(pos_MC[:,col]*bin_size(binsy))
                popt_C_pos,_ = scipy.optimize.curve_fit(gauss,binsy,curveC_i_pos,p0=[10,1]) # popt contains (mu,sigma)  
                curveI_i_pos = pos_MI[:,col] / n.sum(pos_MI[:,col]*bin_size(binsy))
                popt_I_pos,_ = scipy.optimize.curve_fit(gauss,binsy,curveI_i_pos,p0=[10,1])
                if n.abs(popt_C_pos[1]) > n.abs(popt_I_pos[1]): # if width of C is fatter than I
                    curveC = gauss(binsy,popt_C_pos[0],popt_C_pos[1])
                    curveI = gauss(binsy,popt_I_pos[0],popt_I_pos[1])
                    offset = n.argmax(curveC) - n.argmax(curveI) # find offset
                    sigma_diff = n.sqrt(popt_C_pos[1]**2-popt_I_pos[1]**2) # sigmas add in quadrature 
                    convolve_pos[k][:,col] = n.exp(-(n.fft.fftshift(n.arange(binsy.size)))**2/(2*sigma_diff**2)) # convolution gaussian centered at the center
                else: # otherwise just put in a delta function
                    convolve_pos[k][:,col] = n.zeros_like(curveC_i_pos)
                    convolve_pos[k][:,col][curveC_i_pos.size/2] = 1.0
                # Negative half
                if n.all(neg_MC == 0) == False and n.all(neg_MI == 0) == False: # if there are negative values 
                    try: 
                        curveC_i_neg = neg_MC[:,col] / n.sum(neg_MC[:,col]*bin_size(binsy))
                        popt_C_neg,_ = scipy.optimize.curve_fit(gauss,binsy,curveC_i_neg,p0=[10,1]) 
                        curveI_i_neg = neg_MI[:,col] / n.sum(neg_MI[:,col]*bin_size(binsy))
                        popt_I_neg,_ = scipy.optimize.curve_fit(gauss,binsy,curveI_i_neg,p0=[10,1])
                        if n.abs(popt_C_neg[1]) > n.abs(popt_I_neg[1]): 
                            curveC = gauss(binsy,popt_C_neg[0],popt_C_neg[1])
                            curveI = gauss(binsy,popt_I_neg[0],popt_I_neg[1])
                            offset = n.argmax(curveC) - n.argmax(curveI) 
                            sigma_diff = n.sqrt(popt_C_neg[1]**2-popt_I_neg[1]**2) 
                            convolve_neg[k][:,col] = n.exp(-(n.fft.fftshift(n.arange(binsy.size)))**2/(2*sigma_diff**2)) 
                    except: # if curve_fit fails to find a solution
                        convolve_neg[k][:,col] = n.zeros_like(binsx.size)
                        convolve_neg[k][:,col][binsx.size/2] = 1.0
                else: convolve_neg[k] = n.zeros((binsx.size,binsy.size))
        return convolve_pos, convolve_neg

    # Fit polynomial to get one smooth signal loss curve
    # Convolve it based on the extra width "C" has on "I"
    # Return final M matrix (transfer function)
    def transfer_func(identity=False):
        M_matrix = {}
        for kk,k in enumerate(kpl_fold): # only positive k's
            xs = n.array(Pins_fold[k]).flatten() 
            if identity == True: ys = n.array(Pouts_I_fold[k]).flatten() 
            if identity == False: ys = n.array(Pouts_fold[k]).flatten()
            xs_pos,xs_neg,ys_pos,ys_neg, _,_,_,_ = parse_data(xs, ys, ys)
            M_matrix_pos = curve_to_matrix(xs_pos,ys_pos)
            try: M_matrix_neg = curve_to_matrix(xs_neg,ys_neg)
            except: M_matrix_neg = n.zeros_like(M_matrix_pos) # if no neg points
            if identity == False: # do convolution here
                for col in range(M_matrix_pos.shape[1]):
                    new_col_pos = n.convolve(M_matrix_pos[:,col],convolve_pos[k][:,col],mode='same')
                    M_matrix_pos[:,col] = new_col_pos
                    new_col_neg = n.convolve(M_matrix_neg[:,col],convolve_neg[k][:,col],mode='same')
                    M_matrix_neg[:,col] = new_col_neg
            M_matrix[k] = n.concatenate((M_matrix_neg[::-1],M_matrix_pos))
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
            data_dist_pos = data_dist[k][binsy.size:] # separate
            data_dist_neg = data_dist[k][:binsy.size]
            # Get new distribution
            try: M = M_matrix[k]
            except: M = M_matrix[-k] # M only exists for positive k's
           
            Mpos = M[binsy.size:] # separate
            #Mneg = Mpos[::-1] # XXX flipped Mpos
            Mneg = M[:binsy.size]
            """
            if n.all(data_dist_neg==0) == False and n.all(Mneg==0) == True: # rare case when data has negative part of the distribution but there's no transfer function there
                print "WARNING! For k =",k,", you have negative data values but no transfer function for negative P_out's... applying positive transfer function instead."
                M = n.concatenate((Mpos[::-1],Mpos)) # flip Mpos
                for col in range(M.shape[1]): # re-normalize
                    if n.sum(M[:,col]) > 0:
                        M[:,col] /= n.sum(M[:,col]*bin_size(binsy_full))
                Mpos = M[binsy.size:] # separate again
                Mneg = M[:binsy.size]
            """
            # XXX if there are data_dist_neg points but no Mneg, use Mpos instead
            data_neg = n.resize(data_dist_neg,Mneg.shape).T
            Mneg_final = n.zeros_like(Mpos)
            for r in range(len(binsy)):
                if n.all(Mneg[r] == 0) == True and n.all(data_neg[r] == 0) == False: # no Mneg
                    Mneg_final[r] = Mpos[len(binsy)-r-1]*data_neg[r]
                else: Mneg_final[r] = Mneg[r]*data_neg[r]
            Mneg = Mneg_final.copy()
            Mpos *= n.resize(data_dist_pos,Mpos.shape).T # multiply data distribution element-wise per column of M 
            #Mneg = Mneg*(n.resize(data_dist_neg,Mneg.shape).T)
            rowsum_pos = n.zeros(Mpos.shape[0])
            rowsum_neg = n.zeros(Mneg.shape[0])
            for row in Mpos: rowsum_pos += row # add up rows of M
            for row in Mneg: rowsum_neg += row
            rowsum_combine = n.concatenate((rowsum_neg[::-1],rowsum_pos))
            new_PS[k] = rowsum_combine/n.sum(rowsum_combine*bin_size(bins_concat)) # normalize
        return new_PS # distribution of bins_concat
   
    # Shift distribution to new peak
    def shift_dist(data, peakdata, ks):
        shifted_data = {}
        for key in data:
            ind = n.where(ks == key)[0] # ind of peakdata for key
            peak = peakdata[ind]
            old_peak_ind = n.argmax(data[key])
            new_peak_ind = n.argmin(n.abs(bins_concat-peak))
            if n.sign(peak) != n.sign(bins_concat[old_peak_ind]):
                old_peak_ind = n.argmax(data[key][::-1]) # old peak after flipping
                
            diff = new_peak_ind - old_peak_ind
            shifted_bins = bins_concat + diff # shift x-axis
            shifted_data[key] = n.interp(bins_concat, shifted_bins, data[key]) # re-sample bins_concat
            shifted_data[key] /= n.sum(shifted_data[key]*bin_size(bins_concat)) # re-normalize 
        return shifted_data

    ### SIGNAL LOSS CODE

    binsx, binsy, binsy_full, bins_concat = make_bins() # bins are where to sample distributions

    # Get distributions of original data 
    old_pCs = data_dist(pCs)
    old_pCs_fold = data_dist(pCs_fold)
    old_pIs = data_dist(pIs)    
    old_pIs_fold = data_dist(pIs_fold)
     
    # Shift old distributions to peak of no-bootstrapping case
    if count == 0: # data case
        old_pCs = shift_dist(old_pCs, file['pCv'], file['kpl'])
        old_pCs_fold = shift_dist(old_pCs_fold, file['pCv_fold'], file['kpl_fold'])
        old_pIs = shift_dist(old_pIs, file['pIv'], file['kpl'])
        old_pIs_fold = shift_dist(old_pIs_fold, file['pIv_fold'], file['kpl_fold']) 
    else:
        old_pCs = shift_dist(old_pCs, file['pCn'], file['kpl'])
        old_pCs_fold = shift_dist(old_pCs_fold, file['pCn_fold'], file['kpl_fold'])
        old_pIs = shift_dist(old_pIs, file['pIn'], file['kpl'])
        old_pIs_fold = shift_dist(old_pIs_fold, file['pIn_fold'], file['kpl_fold']) 
    
    # Compute signal loss and get new distributions 
    if opts.skip_sigloss:
        print "Skipping Signal Loss!"
        new_pCs = old_pCs.copy()
        new_pCs_fold = old_pCs_fold.copy()
        new_pIs = old_pIs.copy()
        new_pIs_fold = old_pIs_fold.copy()
    else:
        kde_C, kde_I = smooth_dist() # KDE transfer curves
        convolve_pos, convolve_neg = convolve_kernel(kde_C, kde_I) # convolution kernels 
        M = transfer_func(identity=False) # final signal loss transfer function
        M_I = transfer_func(identity=True)
        # new distributions
        new_pCs = sigloss_func(old_pCs, M) 
        new_pCs_fold = sigloss_func(old_pCs_fold, M)
        new_pIs = sigloss_func(old_pIs, M_I)
        new_pIs_fold = sigloss_func(old_pIs_fold, M_I)
    
    # Compute PS points and errors
    def compute_stats(bins,data):
        pts, errs = [], []
        for key in n.sort(data.keys()):
            point = bins[n.argmax(data[key])]
            pts.append(point)
            percents = [n.sum((data[key]*bin_size(bins))[:i]) for i in range(len(data[key]))]
            left = n.interp(0.025,percents,bins)
            right = n.interp(0.975,percents,bins)
            errs.append((right-left)/4) # 1-sigma
        return pts,errs
    
    pC, pC_err = compute_stats(bins_concat, new_pCs)
    pI, pI_err = compute_stats(bins_concat, new_pIs)
    pC_fold, pC_fold_err = compute_stats(bins_concat, new_pCs_fold)
    pI_fold, pI_fold_err = compute_stats(bins_concat, new_pIs_fold)
   
    # Avoid signal gain  (which happens when P_out never goes below the data level)
    def no_gain(pts, errs, old_dist):
        new_errs = []
        new_pts = []
        old_pts, old_errs = compute_stats(bins_concat, old_dist)
        for ee in range(len(errs)):
            if errs[ee] < old_errs[ee]: # if variance after signal loss is smaller than before
                new_errs.append(old_errs[ee])
                new_pts.append(old_pts[ee])
            else: 
                new_errs.append(errs[ee])
                new_pts.append(pts[ee])
        return new_pts, new_errs

    pC, pC_err = no_gain(pC, pC_err, old_pCs)
    pI, pI_err = no_gain(pI, pI_err, old_pIs)
    pC_fold, pC_fold_err = no_gain(pC_fold, pC_fold_err, old_pCs_fold)
    pI_fold, pI_fold_err = no_gain(pI_fold, pI_fold_err, old_pIs_fold)
 
    # Save values to use for plotting sigloss plots
    if count == 0:
        ind = -3 # one k-value
        k = kpl[ind]
        if opts.nosubtract: fn = 'pspec_sigloss_nosubtract.npz'
        else: fn = 'pspec_sigloss.npz'
        print "Saving",fn,"which contains data values for k =",k
        n.savez(fn, k=k, binsx=binsx, binsy=binsy, bins_concat=bins_concat, pC=pC[ind], pC_err=pC_err[ind], pI=pI[ind], pI_err=pI_err[ind], new_pCs=new_pCs[k], new_pIs=new_pIs[k], old_pCs=old_pCs[k], old_pIs=old_pIs[k], Pins=Pins_fold[k], Pouts=Pouts_fold[k], Pouts_I=Pouts_I_fold[k])    
       
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
        pC = pC
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
        prob=0.95, kperp=file['kperp'], sep=opts.sep, kpl_fold=file['kpl_fold'],
        ngps=file['ngps'], nbls=file['nbls'], nbls_g=file['nbls_g'], nlsts_g=file['nlsts_g'],
        lsts=file['lsts'], afreqs=file['afreqs'], cnt_eff=file['cnt_eff'],
        frf_inttime=file['frf_inttime'], inttime=file['inttime'], 
        cmd=file['cmd'].item() + ' \n '+' '.join(sys.argv))

