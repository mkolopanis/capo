import numpy as n,os
from capo import pspec
from capo.cosmo_units import f212z, c
import glob
import ipdb
import matplotlib.pyplot as p
import numpy as np
#measurements

def errorbars(data,axis=1,per=95):
    mean = n.percentile(data,50,axis=axis)
    lower = mean - n.percentile(data, 50-per/2.,axis=axis)
    upper = n.percentile(data,50+per/2.,axis=axis) - mean
    return lower, upper



def PAPER_32_all():
    '''
    Results from  PAPER 32  Jacobs et.al 2014

    outputs results[z] = n.array([k, Delta^2, 2-sigma upper, 2-sigma lower])
    '''

    PAPER_RESULTS_FILES = glob.glob(os.path.dirname(__file__)+'/data/psa32_apj/pspec_*.npz')
    PAPER_RESULTS_FILES.sort()
    freqs= []
    for filename in PAPER_RESULTS_FILES:
        try:
            freqs.append(n.load(filename)['freq']*1e3)
        except(KeyError):
            try:
                dchan = int(filename.split('/')[-1].split('.')[0].split('_')[2])-int(filename.split('/')[-1].split('.')[0].split('_')[1])
                chan = int(filename.split('/')[-1].split('.')[0].split('_')[1]) + dchan/2.
                freqs.append(chan/2. + 100) #a pretty good apprximation of chan 2 freq for 500kHz channels
            except: continue
    freqs = n.array(freqs)
    zs = pspec.f2z(freqs*1e-3)
    results = {}
    for files,z in zip(PAPER_RESULTS_FILES,zs):
        f=n.load(files)
        results[z] = n.array([f['k'],f['k3pk'],f['k3pk']+f['k3err'],f['k3pk']-f['k3err']]).T

    return results

def PAPER_64_all():
    '''
    Results from  PAPER 64  Ali et.al 2015

    outputs results[z] = n.array([k, Delta^2, 2-sigma upper, 2-sigma lower])
    '''

    PAPER_RESULTS_FILES = glob.glob(os.path.dirname(__file__)+'/data/psa64_apj/pspec_*.npz')
    PAPER_RESULTS_FILES.sort()
    freqs = []
    for filename in PAPER_RESULTS_FILES:
        try:
            freqs.append(n.load(filename)['freq']*1e3)
        except(KeyError):
            try:
                dchan = int(filename.split('/')[-1].split('.')[0].split('_')[2])-int(filename.split('.')[0].split('_')[1])
                chan = int(filename.split('/')[-1].split('.')[0].split('_')[1]) + dchan/2.
                freqs.append(chan/2. + 100) #a pretty good apprximation of chan 2 freq for 500kHz channels
            except: continue
    freqs = n.array(freqs)
    zs = pspec.f2z(freqs*1e-3)
    #zs = n.array([8.31])
    results = {}
    for files,z in zip(PAPER_RESULTS_FILES,zs):
        f=n.load(files)
        results[z] = n.array([f['k'],f['k3pk'],f['k3pk']+f['k3err'],f['k3pk']-f['k3err']]).T

    return results

def MWA_32T_all():
    """
    From Josh Dillon April 10, 2014
    Delta^2  |  2sig bottom bar  |  2sig top bar  |  k  |  z | "Detection"

    or

    2sig top bar   |  [ignore this entry]  |  2sig top bar  |  k  |  z | "2sigUL"

    in the units you requested.

    (my request The ideal format would be a text file listing three things:
    Delta^2 [mk^2], k [hMpc^-1], z [z])


    TREATMENT:
    I will store the delta^2 +/- values.  If he doesn't give me a delta^2 or lower limit, I will put
    Delta^2 = 0 and bottom = -top

    return format will be dict[z] = n.array([[k,Delta^2,top,bottom]]) all in mK^2

    """


    MWA_RESULTS_FILE=os.path.dirname(__file__)+'/data/MWA_32T_Results_Summary.txt'
    lines = open(MWA_RESULTS_FILE).readlines()
    result = {}
    for line in lines:
        if line.endswith('2sigUL'):#if its just an upper limit
            l = map(float,line.split()[:-1])
            z = l[-1]
            try:
                result[z].append([l[3],0,l[0],-l[0]])
            except(KeyError):
                result[z] = [[l[3],0,l[0],-l[0]]]
        else:#if its a "detection" which just means its not consistent with 0, not really EoR but probably foregrounds
            #or systematics
            l = map(float,line.split()[:-1])
            z = l[-1]
            try:
                result[z].append([l[3],l[0],l[0],-l[0]])
            except(KeyError):
                result[z] = [[l[3],l[0],l[2],l[1]]]
    for z in result.keys():
        result[z] = n.array(result[z])
    return result

def MWA_128_all():
    '''
    MWSA_128 data from dillion 2015
    return format will be dict[z] = n.array([[k,Delta^2,top,bottom]]) all in mK^2
    '''
    MWA_RESULTS_FILE=glob.glob(os.path.dirname(__file__)+'/data/mwa128/*.dat')
    zs= []
    results = {}
    for files in MWA_RESULTS_FILE:
        name = files.split('/')[-1]
        nums = name.split('=')[1].split('.')[:2]
        z= float(nums[0]+'.'+nums[1])

        data = n.genfromtxt(files, delimiter=' ')
        results[z] = data[:,[0,1,5,4]]

    return results
def z_slice(redshift,pspec_data):
    """
    input a power spectrum data dict output of MWA_32T_all() or GMRT_2014_all()
    returns a slice along k for the input redshift
    example
    z,pspec[k,k3pK] = z_slice(MWA_32T_all())
    """
    zs = n.array(pspec_data.keys())
    closest_z = zs[n.abs(zs-redshift).argmin()]
    return closest_z, pspec_data[closest_z]
def k_slice(k,pspec_data):
    """
    input a power spectrum data dict output of MWA_32T_all() or GMRT_2014_all()
    returns a slice along z for the input redshift
    example
    zs,pspec[k,k3pK] = k_slice(MWA_32T_all())
    """

    zs = n.array(pspec_data.keys())
    k_is = [n.abs(pspec_data[redshift][:,0]-k).argmin() for redshift in zs]
    ks = [pspec_data[redshift][k_i,0] for k_i in k_is]
    power = n.vstack([pspec_data[redshift][k_i,:] for k_i in k_is])
    return zs,power


def MWA_32T_at_z(redshift):
    """
    Input a redshift
    Return the Delta^2(k) power spectrum at the redshift nearest the input redshift
    returns: z,array([k,Delta^2,2sig_upper_limit,2sig_lower_limit])
    """
    data = MWA_32T_all()
    zs = n.array(data.keys())
    closest_z = zs[n.abs(zs-redshift).argmin()]
    return closest_z, data[closest_z]
def MWA_32T_at_k(k):
    """
    Input a k
    Returns the Delta^2(k) power spectrum vs redshift at that k. (all pspec in mk^2)
    return format: z,array([k,Delta^2,2sig_upper_limit,2sig_lower_limit])
    """
    data = MWA_32T_all()
    zs = n.array(data.keys())
    k_is = [n.abs(data[redshift][:,0]-k).argmin() for redshift in zs]
    ks = [data[redshift][k_i,0] for k_i in k_is]
    power = n.vstack([data[redshift][k_i,:] for k_i in k_is])
    return zs,power
def GMRT_2014_all():
    return {8.6:n.array(
            [[0.1, 0,2e5,  -2e5],
            [0.13, 0,4e5,  -4e5],
            [0.16, 0,1e5,  -1e5],
            [0.19, 0,1.9e5,-1.9e5],
            [0.275,0,2e5,  -2e5],
            [0.31, 0,4e5,  -4e5],
            [0.4,  0,6e5,  -6e5],
            [0.5,  0,8e4,  -8e4]])}

def get_pk_from_npz(files=None, verbose=False):
    """
    Load output from plot_pk_k3pk.npz and returns P(k)  spectrum.

    Return lists of k, Pk, Pk_err, Delta^2 ordered by decreasing redshift
    Return format: z, k_parallel, Pk, Pk_err
    """
    if files is None:
        print 'No Files gives for loading'
        return 0, '', '', ''

    if len(n.shape(files)) == 0:
        files = [files]

    freqs = []
    if verbose:
        print "parsing npz file frequencies"
    for filename in files:
        if verbose:
            print filename,
        try:
            if verbose:
                print "npz..",
            freqs.append(n.load(filename)['freq']*1e3)  # load freq in MHz
            if verbose:
                print "[Success]"
        except(KeyError):
            if verbose:
                print "[FAIL]"
            try:
                if verbose:
                    print "looking for path like RUNNAME/chan_chan/I/pspec.npz"
                dchan = (int(filename.split('/')[1].split('_')[1]) -
                         int(filename.split('/')[1].split('_')[0]))
                chan = int(filename.split('/')[1].split('_')[0]) + dchan/2
                freqs.append(chan/2. + 100)
                # a pretty good apprximation of chan 2 freq for 500kHz channels
            except(IndexError):
                if verbose:
                    print "[FAIL] no freq found. Skipping..."

    if len(freqs) == 0:  # check if any files were loaded correctly
        print 'No parsable frequencies found'
        print 'Exiting'
        return 0, '', '', ''

    if verbose:
        print "sorting input files by frequency"
    files = n.array(files)
    files = files[n.argsort(freqs)]
    freqs = n.sort(freqs)
    if verbose:
        print "found freqs"
    freqs = n.array(freqs)
    if verbose:
        print freqs

    z = f212z(freqs*1e6)
    if verbose:
        print "processing redshifts:", z

    kpars = []
    Pks = []
    Pkerr = []
    for i, FILE in enumerate(files):
        F = n.load(FILE)
        if verbose:
            print FILE.split('/')[-1], z[i]
        Pks.append(F['pk'])
        kpars.append(F['kpl'])
        Pkerr.append(F['err'])
    return z, kpars, Pks, Pkerr


def get_k3pk_from_npz(files=None, verbose=False):
    """
    Load output from plot_pk_k3pk.npz and returns Delta^2 spectrum.

    Return lists of k, Delta^2, Delta^2_err ordered by  decreasing redshift
    Return format: z, k_magnitude, Delta^2, Delta^2_err
    """
    if files is None:  # check that files are passed
        print 'No Files given for loading'
        return 0, '', '', ''

    if len(n.shape(files)) == 0:
        files = [files]

    freqs = []
    if verbose:
        print "parsing npz file frequencies"
    for filename in files:
        if verbose:
            print filename,
        try:
            if verbose:
                print "npz..",
            freqs.append(n.load(filename)['freq']*1e3)  # load freq in MHz
            if verbose:
                print "[Success]"
        except(KeyError):
            if verbose:
                print "[FAIL]"
            try:
                if verbose:
                    print "looking for path like RUNNAME/chan_chan/I/pspec.npz"
                dchan = (int(filename.split('/')[1].split('_')[1]) -
                         int(filename.split('/')[1].split('_')[0]))
                chan = int(filename.split('/')[1].split('_')[0]) + dchan/2
                freqs.append(chan/2. + 100)
                # a pretty good apprximation of chan 2 freq for 500kHz channels
            except(IndexError):
                if verbose:
                    print "[FAIL] no freq found. Skipping..."

    if len(freqs) == 0:  # check if any files were loaded correctly
        print 'No parsable frequencies found'
        print 'Exiting'
        return 0, '', '', ''

    if verbose:
        print "sorting input files by frequency"
    files = n.array(files)
    files = files[n.argsort(freqs)]
    freqs = n.sort(freqs)

    if verbose:
        print "found freqs"
    freqs = n.array(freqs)
    if verbose:
        print freqs

    z = f212z(freqs*1e6)
    if verbose:
        print "processing redshifts:", z
    k3Pk = []
    k3err = []
    kmags = []
    for i, FILE in enumerate(files):
        F = n.load(FILE)
        if verbose:
            print FILE.split('/')[-1], z[i]
        k3Pk.append(F['k3pk'])
        k3err.append(F['k3err'])
        kmags.append(F['k'])
    return z, kmags, k3Pk, k3err

def posterior(kpl, pk, err, pkfold=None, errfold=None, f0=.151, umag=16.,
            theo_noise=None,verbose=False):
    import scipy.interpolate as interp
    k0 = n.abs(kpl).argmin()
    kpl = kpl[k0:]
    z = pspec.f2z(f0)
    kpr = pspec.dk_du(z) * umag
    k = n.sqrt(kpl**2 + kpr**2)
    if pkfold is None:
        if verbose: print 'Folding for posterior'
        pkfold = pk[k0:].copy()
        errfold = err[k0:].copy()
        pkpos,errpos = pk[k0+1:].copy(), err[k0+1:].copy()
        pkneg,errneg = pk[k0-1:0:-1].copy(), err[k0-1:0:-1].copy()
        pkfold[1:] = (pkpos/errpos**2 + pkneg/errneg**2) / (1./errpos**2 + 1./errneg**2)
        errfold[1:] = n.sqrt(1./(1./errpos**2 + 1./errneg**2))
        #ind = n.logical_and(kpl>.2, kpl<.5)
    ind = n.logical_and(k>.15, k<.5)
    #ind = n.logical_and(kpl>.12, kpl<.5)
    #print kpl,pk.real,err
    k = k[ind]
    pkfold = pkfold[ind]
    errfold = errfold[ind]
    #if not theo_noise is None:
    #    theo_noise=theo_noise[ind]
    pk= pkfold
    err =errfold
    err_omit = err.copy()
    #s = n.logspace(1,3.5,100)
    s = n.linspace(-5000,5000,10000)
    #    print s
    data = []
    data_omit = []
    for _k, _pk, _err in zip(k, pk, err):
        if verbose: print _k, _pk.real, _err
    #    print '%6.3f    %9.5f     9.5f'%(_k, _pk.real, _err)
    for ss in s:
        data.append(n.exp(-.5*n.sum((pk.real - ss)**2 / err**2)))
        data_omit.append(n.exp(-.5*n.sum((pk.real - ss)**2 / err_omit**2)))
    #    print data[-1]
    data = n.array(data)
    data_omit = n.array(data_omit)
    #print data
    #print s
    #data/=n.sum(data)
    data /= n.max(data)
    data_omit /= n.max(data_omit)
    p.figure(5, figsize=(6.5,5.5))
    p.plot(s, data, 'k', linewidth=2)
#    p.plot(s, data_omit, 'k--', linewidth=1)
    #use a spline interpolator to get the 1 and 2 sigma limits.
    #spline = interp.interp1d(data,s)
    #print spline
    #print max(data), min(data)
    #print spline(.68), spline(.95)
    #p.plot(spline(n.linspace(.0,1,100)),'o')
#    p.plot(s, n.exp(-.5)*n.ones_like(s))
    #    p.plot(s, n.exp(-.5*2**2)*n.ones_like(s))
    data_c = n.cumsum(data)
    data_omit_c = n.cumsum(data_omit)
    data_c /= data_c[-1]
    data_omit_c /= data_omit_c[-1]
    mean = s[n.argmax(data)]
    s1lo,s1hi = s[data_c<0.1586][-1], s[data_c>1-0.1586][0]
    s2lo,s2hi = s[data_c<0.0227][-1], s[data_c>1-0.0227][0]
    if verbose: print 'Posterior: Mean, (1siglo,1sighi), (2siglo,2sighi)'
    if verbose: print 'Posterior:', mean, (s1lo,s1hi), (s2lo,s2hi)
    mean_o = s[n.argmax(data_omit)]
    s1lo_o,s1hi_o = s[data_omit_c<0.1586][-1], s[data_omit_c>1-0.1586][0]
    s2lo_o,s2hi_o = s[data_omit_c<0.0227][-1], s[data_omit_c>1-0.0227][0]
    if verbose: print 'Posterior (omit):', mean_o, (s1lo_o,s1hi_o), (s2lo_o,s2hi_o)

    p.vlines(s1lo,0,1,color=(0,107/255.,164/255.), linewidth=2)
    p.vlines(s1hi,0,1,color=(0,107/255.,164/255.), linewidth=2)

    # limits for data_omit
    p.vlines(s2lo,0,1,color=(1,128/255.,14/255.), linewidth=2)
    p.vlines(s2hi,0,1,color=(1,128/255.,14/255.), linewidth=2)

    if not theo_noise is None:
        s2l_theo=n.sqrt(1./n.mean(1./theo_noise**2))
        p.vlines(s2l_theo,0,1,color='black',linewidth=2)
        if verbose: print('Noise level: {0:0>5.3f} mk^2'.format(s2l_theo))
    p.xlabel(r'$k^3/2\pi^2\ P(k)\ [{\rm mK}^2]$', fontsize='large')
    p.ylabel('Posterior Distribution', fontsize='large')
    p.xlim(0,700)
    p.title('z = {0:.2f}'.format(z))
    if (s2lo > 700) or (s2hi > 1000):
        p.xlim(0,1500)
    p.grid(1)
    p.subplots_adjust(left=.15, top=.95, bottom=.15, right=.95)
    p.savefig('posterior_{0:.2f}.png'.format(z))
    f=open('posterior_{0:.2f}.txt'.format(z), 'w')
    f.write('Posterior: Mean,\t(1siglo,1sighi),\t(2siglo,2sighi)\n')
    f.write('Posterior: {0:.4f},\t({1:.4f},{2:.4f}),\t({3:.4f},{4:.4f})\n'.format( mean, s1lo,s1hi, s2lo,s2hi))
    f.write( 'Posterior (omit): {0:.4f}, ({1:.4f},{2:.4f}),\t({3:.4f},{4:.4f})\n'.format( mean_o, s1lo_o,s1hi_o, s2lo_o,s2hi_o))
    f.write( 'Noise level: {0:0>5.3f} mk^2\n'.format(s2l_theo) )
    f.close()
#Danny's power spectrum bits
def read_bootstraps_dcj(filenames,verbose=False):
    #read in a list of bootstrapped power spectra
    #return a single set of power spectra stacked along the bootstrap dimension
    #only keep the real part!
    """
    ['err_vs_t',    #not sure
     'cmd',         #the command used to generate the file
     'pCv',         #the weighted data power spectrum (no injection) (times,kpls)
     'pC',     #the weighted data power spectrum with injection  (times,kpls)
     'times',       #lsts of data bins
     'scalar',      #conversion from mk^2 to mK^2/h^3Mpc^3 (already applied to data)
     'pI',  #unweighted injected signal
     'freq',        #center frequency of bin in GHz
     'kpl',         #list of k parallels matching the kpl axis of the power spectrum
     'var',  #not sure
     'pIv']             #unweighted power spectrum of data (no injection)

    """
    accumulated_power_spectra = {}
    for filename in filenames:
        F = np.load(filename)
        for thing in F.files:
            try:
                accumulated_power_spectra[thing].append(F[thing])
            except(KeyError):
                accumulated_power_spectra[thing] = [F[thing]]
    power_spectrum_channels = ['pC','pI','err','pCv','var','pIv',
                               'pCe', 'pIe', 'pIr', 'pCr']
    #stack up the various power spectrum channels
    for key in accumulated_power_spectra:
        if key in power_spectrum_channels:
            accumulated_power_spectra[key] = np.real(np.array(accumulated_power_spectra[key]))
        else:    #otherwise just keep the first entry,
                 #   assuming they are all the same but
                 #   for their bootstrapping
            accumulated_power_spectra[key] = accumulated_power_spectra[key][0]
    return accumulated_power_spectra

def average_bootstraps(indata,Nt_eff,avg_func=np.median,Nboots=100):
    """
    "Scramble average" the various power spectrum channels across time
    and then get the mean and standard deviation across scrambles
    compute the error as the standard deviation across scramble
    Assumed axes: (bootstrap,k,time)
    input:
        indata:a dictionary of arrays with names as output by read_bootstraps
        Nt_eff:effective number of independent time samples
    output: a matching dictionary of arrays such as read in by power spectrum
    plotting tools.
    """
    pspec_channels = {  'pCv':'pCv',
                        'pIv':'pIv',
                        'pCe':'pCe',
                        'pIe':'pIe',
                        'pCr':'pCr',
                        'pIr':'pIr'}
    outdata = {}
    for inname in indata:
        if inname in pspec_channels.keys():
            outname = pspec_channels[inname]
            #scramble the times and bootstraps. Average over new time dimension
            #   only draw as many times as we have independent lsts (Nt_eff)
            Z = scramble_avg_bootstrap_array(indata[inname],
                        Nt_eff=Nt_eff,func=avg_func,Nboots=Nboots)
            #power spectrum is the mean and standard dev over scramble dimension
            outdata[outname] = avg_func(Z,axis=0)
            outdata[outname+'_err'] = np.std(Z,axis=0)

            #also do the folded version
            outname += '_fold'
            kpl_fold,X = split_stack_kpl(indata[inname],indata['kpl'])
            Z = scramble_avg_bootstrap_array(X,
                        Nt_eff=Nt_eff,func=avg_func,Nboots=Nboots)
            outdata[outname] = avg_func(Z,axis=0)
            outdata[outname+'_err'] = np.std(Z,axis=0)
            outdata['kpl_fold'] = kpl_fold


        else:
            outdata[inname] = indata[inname]
    return outdata

def scramble_avg_bootstrap_array(X,Nt_eff=10,Nboots=100,func=np.median):
    #choose randomly a time (axis=-1) from a random bootstrap (axis=-2)
    #apply func to the result (default is numpy.median)
    #do for NBOOT iterations
    #assumes input array dimensions (nbootstraps,nks,ntimes)
    bboots = []
    for i in xrange(Nboots):
        times_i = np.random.choice(X.shape[-1],Nt_eff,replace=True)
        bls_i = np.random.choice(X.shape[0],Nt_eff,replace=True)
        bboots.append(X[bls_i,:,times_i].squeeze().T)
    bboots = n.ma.masked_invalid(np.array(bboots))
    return func(bboots,axis=-1)

def split_stack_kpl(X,kpl):
    #split the input X array at kpl=0 and stack along the bootstrap dimension
    #  use in concert with scramble_avg_bootstraps to fold kpls together
    #assumes input dimensions (nbootstraps,nks,ntimes)
    assert(X.shape[1]==len(kpl)) #make sure that kpl matches the kpl axis
    #if theres an odd number of kpls, then there better be a zero value
    assert(len(kpl)%2==0 or np.abs(kpl).min()==0)
    if np.abs(kpl).min()==0:
        kpl0_split_index = np.argmin(np.abs(kpl))
        X_kpos = X[:,kpl0_split_index:,:]
        X_kneg = X[:,kpl0_split_index::-1,:] #select and flip simultaneously
    else:
        kpl0_split_index = np.where(np.logical_and(kpl>0,np.abs(kpl)==np.min(np.abs(kpl))))[0][0]
        X_kpos = X[:,kpl0_split_index:,:]
        X_kneg = X[:,kpl0_split_index-1::-1,:] #select and flip simultaneously
    return kpl[kpl0_split_index:], np.concatenate([X_kpos,X_kneg],axis=0)


#Matt's power spectrum bits
def read_bootstraps(files=None, verbose=False):
    """
    Read bootstrap files and accumulate into dict.

    arguments:
        files: glob of files, or list of file names to be read

    keywords:
        verbose: Print optional output to stdout. Defalt False
    """
    if files is None or not files:
        raise TypeError('Files given are {0}; Must supply input '
                        'files'.format(files))
        return files

    if len(n.shape(files)) == 0:
        files = [files]

    # load the first file to make dummy lists into which boots will aggregate
    npz0 = n.load(files[0])
    num_boots = len(files)
    keys = npz0.keys()
    npz0.close()

    out_dict = {key: [] for key in keys}

    for filename in files:
        if verbose:
            print 'Reading Boots'
        npz = n.load(filename)
        for key in keys:
            out_dict[key].append(n.real(npz[key]))
        npz.close()

    return out_dict

def read_injects(inj_dirs=None):
    '''
    read_inject(inj_dirs)

    iterates over inject directories and runs read_bootstraps on each directory

    arguments:
        inj_dirs: glob or list of inject directories to be loaded

    returns:
        dictionary of outputs from read_bootstraps, keys are the inject_directory names
    '''
    if inj_dirs is None or not inj_dirs:
        raise TypeError('Must supply input list of inject directories')

    ##create list of keys by taking only last party of file name
    ##this could be a problem if you try to read two different channel ranges
    ##wit the same inject values but that sounds like a crazy thing to do.
    keys = [inj.split('/')[-1] for inj in inj_dirs]
    out_dict ={key:{} for key in keys}
    for cnt,key in enumerate(keys):
        in_files = glob.glob( inj_dirs[ cnt ] + '/pspec_boo*.npz' )
        out_dict[key] = read_bootstraps(in_files)
    return out_dict

def random_avg_bootstraps(boot_dict = None,boot_axis=None, time_axis=None,
    outfile=None, verbose=False, nboot=400):
    '''
    random_avg_bootstraps(boot_dict=None, boot_axis=None, time_axis=None
            outfile=None, verbose=False, nboot=400)

    arguments:
        boot_dict: dictionary object of lists with at least 2 dimensions.


        boot_axis: dimension of the lists in boot_dict over which different
                    bootstraps are stored. Required argument. defualt = None

        time_axis: dimension of the lists in boot_dict over which different
                    times are stored. Required argument. defualt = None

    keywords:
        outfile: Optional outfile into which the random averaged bootstraps
                  will be saved.

        nboot: Number of times a random power spectrum will be constructed
                by forming a waterfall with each time element randomly selected
                from a random bootstrap

        verbose: Print optional output to stdout. Default = False.
    '''
    #Check if any required intput is not given
    if boot_dict is None or not boot_dict:
        raise TypeError('Must supply input dictionary of data')
    if boot_axis is None:
        raise TypeError('Must supply boot axis argument')
    if time_axis is None:
        raise TypeError('Must supply time axis argument')

    #check if either axis, if given, is not an integer
    if not isinstance(boot_axis,(int,long)):
        raise TypeError('Expected Integer type for boot_axis '
                            'instead got {0}'.format(type(boot_axis).__name__))
    if not isinstance(time_axis,(int,long)):
        raise TypeError('Expected Integer type for time_axis '
                            'instead got {0}'.format(type(time_axis).__name__))


    # import ipdb; ipdb.set_trace()
    keys = boot_dict.keys()
    out_dict = {key:[] for key in keys}
    # import ipdb; ipdb.set_trace()

    num_times = n.shape(boot_dict[keys[0]])[time_axis]
    num_boots = n.shape(boot_dict[keys[0]])[boot_axis]

    for boot in xrange(nboot):
        if verbose:
            if (boot+1) % 10 == 0:
                    print '   ', boot+1,'/',nboot
        dsum_dict = {key:[] for key in keys}
        # import ipdb; ipdb.set_trace()
        ts = n.random.choice(num_times,num_times)
        bs = n.random.choice(num_boots,num_times)
        for key in keys:
            tmp_dict = n.swapaxes(boot_dict[key],time_axis,-1)
            tmp_dict = n.swapaxes(tmp_dict,boot_axis,-2)
            # import ipdb; ipdb.set_trace()
            dsum_dict[key] = n.array(tmp_dict)[...,bs,ts].real
        # import ipdb; ipdb.set_trace()
        for key in keys:
            tmp = n.median(dsum_dict[key],-1)
            out_dict[key].append(tmp)

    if outfile: n.savez(outfile, **out_dict)
    return out_dict
