#! /usr/bin/env python
"""Compute 2-Dimensional power spectrum from uv files.

Computes the (k// vs t) power spectrum from uv files.
Bootstraps computed over subsets of Baselines in UV files.

Differs from pspec_oqe_2d.py in that it only uses uniform weighting (I)
and doesn't compute all PS channels.
"""

import aipy as a
import numpy as n
import pylab as p
import glob
import optparse
import sys
import random
from capo import zsa, oqe, cosmo_units, fringe
import capo
import time

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True, cal=True)
o.add_option('-b', '--nboot', type='int', default=1,
             help='Number of bootstraps.  Default is 1 (no bootstrapping).')
o.add_option('--plot', action='store_true',
             help='Generate plots')
o.add_option('--window', dest='window', default='blackman-harris',
             help=('Windowing function to use in delay transform. '
                   'Default is blackman-harris. '
                   'Options are: ' + ', '.join(a.dsp.WINDOW_FUNC.keys())))
o.add_option('--sep', default='sep0,1', action='store',
             help='Which separation directory to use for signal loss data.')
o.add_option('-i', '--inject', type='float', default=0.,
             help='EOR injection level.')
o.add_option('--frf', action='store_true',
             help=('The data to be analyzed has been fringe-rate-filtered.'
                   ' Consequently, the injected EoR will be FRF twice '
                   'and noise will be FRF once.'))
o.add_option('--output', type='string', default='',
             help='Output directory for pspec_boot files (default "")')
o.add_option('--weight', type='string', default='L^-1',
             help=('Choice for MC normalization '
                   'Options available L^-1 F^-1/2 I F^-1'))
o.add_option('--Trcvr', type='float', default=200,
             help='Receiver Temperature in Kelvin (defualt 200)')
o.add_option('--rmbls', dest='rmbls', type='string',
             help=('List of baselines (ex:1_4,2_33) '
                   'to remove from the power spectrum analysis.'))
o.add_option('--NGPS', type='int', default=5,
             help='Number of Groups used in bootstrapping (default 5)')
o.add_option('--lmode', type='int', default=None,
             help=('Eigenvalue mode of C (in decreasing order)'
                   ' to be the minimum value used in C^-1'))
o.add_option('--changeC', action='store_true',
             help='Change covariance matrix C.')
o.add_option('--mode_num', default=None,
             help=('Number to dial if changing regularization'
                   ' strength in pspec_batch_modeloop.py'))
o.add_option('--nbls', default='all',
            help=('Number of baselines to use in analysis. Default is all.'))
opts, args = o.parse_args(sys.argv[1:])


# Basic parameters
#random.seed(0)  # for oqe.py (eor generator)
n.random.seed(0)  # for noise generator
POL = opts.pol
if POL == 'xx' or POL == 'yy': NPOL = 1
else: NPOL = 2
DELAY = False
NGPS = opts.NGPS  # number of groups to break the random sampled bls into
PLOT = opts.plot
INJECT_SIG = opts.inject
LMODE = opts.lmode

try:
    rmbls = []
    rmbls_list = opts.rmbls.split(',')
    for bl in rmbls_list:
        i, j = bl.split('_')
        rmbls.append(a.miriad.ij2bl(int(i), int(j)))
    print 'Removing baselines:', rmbls
    # rmbls = map(int, opts.rmbls.split(','))
except:
    rmbls = []


# CLASSES & FUNCTIONS #

class DataSet(oqe.DataSet):
    """Extention of oqe.DataSet to include some covariance regularization."""

def complex_noise(size, noiselev):
    """Generate complex noise of given size and noiselevel."""
    # if noiselev <= 0 or n.isinf(noiselev):
    #     return n.zeros(size)
    noise_real = n.random.normal(size=size, scale=noiselev)/n.sqrt(2)
    noise_imag = n.random.normal(size=size, scale=noiselev)/n.sqrt(2)
    noise = noise_real + 1j*noise_imag
    return noise

def make_noise(d, cnt, inttime, df):  # , freqs, jy2T=None):
    """Create noise with T_rms matching data from T_rcvr and uv['cnt']."""
    # if jy2T is None:
    #    jy2T = capo.pspec.jy2T(freqs)
    Tsys = 180. * n.power(afreqs/0.18, -2.55) + opts.Trcvr  # system temp in K
    Tsys *= 1e3  # system temp in mK
    Trms = Tsys/n.sqrt(df * 1e9 * inttime * cnt * NPOL)  # convert sdf to Hz
    #                     Bchan, inttime, counts (times 2 if even&odd), Npol
    Vrms = Trms  # /jy2T  # jy2T is in units of mK/Jy
    # The following transposes are to create noise correlated in time not
    # frequency. Maybe there is a better way to do it?
    # The masking and filling is to be able to parallelize the noise draws
    # Setting the mask back later and filling zeros out where Vrms ~ inf or < 0
    size = Vrms.shape[0]
    # if opts.frf: # triple size
    #    Vrms = n.repeat(Vrms, 3, axis=0)
    Vrms = n.ma.masked_invalid(Vrms)
    Vrms.mask = n.ma.mask_or(Vrms.mask, Vrms.filled() < 0)
    n.ma.set_fill_value(Vrms, 1e-20)
    noise = n.array([complex_noise(v1.shape, v1.filled()) for v1 in Vrms.T]).T
    noise = n.ma.masked_array(noise)
    noise.mask = Vrms.mask
    n.ma.set_fill_value(noise, 0 + 0j)
    noise = noise.filled()
    # wij = n.ones(noise.shape, dtype=bool)
    # XXX flags are all true (times,freqs)
    # if opts.frf: # FRF noise
    #    noise = fringe_rate_filter(aa, noise, wij, ij[0], ij[1],
    #                               POL, bins, fir)
    # noise = noise[int(size):2*int(size),:]
    # noise.shape = d.shape
    return noise

def fringe_rate_filter(aa, dij, wij, i, j, pol, bins, fir):
    """Apply frf."""
    _d, _w, _, _ = fringe.apply_frf(aa, dij, wij, i, j,
                                    pol=pol, bins=bins, firs=fir)
    return _d

def make_eor(shape):  # Create and fringe rate filter noise
    """Generate White Noise and Apply FRF."""
    dij = oqe.noise(size=shape)
    return dij

def make_PS(keys, dsv, dsn,
            grouping=True):
    """Use OQE formalism to generate power spectrum.

    Output weighted and identity weightings.
    """
    if grouping:
        newkeys, dsIv = dsv.group_data(keys, gps, use_cov=False)
        newkeys, dsIn = dsn.group_data(keys, gps, use_cov=False)
    else:  # no groups (slower)
        newkeys = keys
        dsIv = dsv
        dsIn = dsn
    pIvs = []
    pIns = []
    for k, key1 in enumerate(newkeys):
        if k == 1 and len(newkeys) == 2:
            # NGPS = 1 (skip 'odd' with 'even' if we already did 'even' with 'odd')
            continue
        # print '   ',k+1,'/',len(newkeys)
        for key2 in newkeys[k:]:
            if len(newkeys) > 2 and (key1[0] == key2[0] or key1[1] == key2[1]):
                # NGPS > 1
                continue
            if key1[0] == key2[0]:  # don't do 'even' with 'even', for example
                continue
            #if key1[1] == key2[1] or (key1[0] == 'even' or key2[0] == 'even'): continue # for jackknife with only even or odd (comment out above if statements to use this one instead)
            else:
                FIv = dsIv.get_F(key1, key2, use_cov=False, cov_flagging=False)
                qIv = dsIv.q_hat(key1, key2, use_cov=False, cov_flagging=False)
                MIv, WIv = dsIv.get_MW(FIv, mode='I')
                pIv = dsIv.p_hat(MIv, qIv, scalar=scalar)
                pIvs.append(pIv)
 
                FIn = dsIn.get_F(key1, key2, use_cov=False, cov_flagging=False)
                qIn = dsIn.q_hat(key1, key2, use_cov=False, cov_flagging=False)
                MIn, WIn = dsIn.get_MW(FIn, mode='I')
                pIn = dsIv.p_hat(MIn, qIn, scalar=scalar)
                pIns.append(pIn)
    return n.array(pIvs), n.array(pIns)

# --------------------------------------------------------------------


# Read even&odd data
if 'even' in args[0] or 'odd' in args[0]:
    dsets = {'even': [x for x in args if 'even' in x],
             'odd': [x for x in args if 'odd' in x]
             }
else:
    dsets = {'even': args, 'odd': args}
print dsets

# Get uv file info
WINDOW = opts.window
uv = a.miriad.UV(dsets.values()[0][0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
sdf = uv['sdf']
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
inttime = uv['inttime']
# try to take the frf_inttime from the file
# for old filtered files, need to use parameter still
try: frf_inttime = uv['FRF_NEBW']
except: frf_inttime = inttime
print 'inttime:', inttime
print 'frf_inttime:', frf_inttime
afreqs = freqs.take(chans)
nchan = chans.size
fq = n.average(afreqs)
z = capo.pspec.f2z(fq)
aa = a.cal.get_aa(opts.cal, freqs)  # all freqs
bls, conj = capo.red.group_redundant_bls(aa.ant_layout)
sep2ij, blconj, bl2sep = capo.zsa.grid2ij(aa.ant_layout)
jy2T = capo.pspec.jy2T(afreqs)
window = a.dsp.gen_window(nchan, WINDOW)
if not WINDOW == 'none':
    window.shape = (nchan, 1)

# B = sdf * afreqs.size / capo.pfb.NOISE_EQUIV_BW[WINDOW]
# this is wrong if we aren't inverting
# the window post delay transform
# (or at least dividing out by the gain of the window)
# For windowed data, the FFT divides out by the full bandwidth, B, which is
# then squared. Proper normalization is to multiply by
# B**2 / (B / NoiseEqBand) = B * NoiseEqBand
# XXX NEED TO FIGURE OUT BW NORMALIZATION
B = sdf * afreqs.size * capo.pfb.NOISE_EQUIV_BW[WINDOW]  # proper normalization
etas = n.fft.fftshift(capo.pspec.f2eta(afreqs))
# create etas (fourier dual to frequency)
kpl = etas * capo.pspec.dk_deta(z)
if True:
    bm = n.polyval(capo.pspec.PAPER_BEAM_POLY, fq) * 2.35
    if opts.frf: bm *= 1.3  # correction factor for FRF omega_pp = .32/.24 = 1.3
    # correction for beam^2
    scalar = capo.pspec.X2Y(z) * bm * B
else:
    scalar = 1
if not DELAY:
    # XXX this is a hack
    if WINDOW == 'hamming':
        scalar /= 3.67
    elif WINDOW == 'blackman-harris':
        scalar /= 5.72
print 'Freq:', fq
print 'z:', z
print 'B:', B
print 'scalar:', scalar
sys.stdout.flush()

# Prep FRF Stuff
antstr = 'cross'
_, blconj, _ = zsa.grid2ij(aa.ant_layout)
days = dsets.keys()
s, d, f = capo.miriad.read_files([dsets[days[0]][0]],
                                 antstr=antstr, polstr=POL)  # read first file
ij = d.keys()[0]  # use first baseline
if blconj[a.miriad.ij2bl(ij[0], ij[1])]:
    # makes sure FRP will be the same whether bl is a conjugated one or not
    if ij[0] < ij[1]:
        temp = (ij[1], ij[0])
        ij = temp
sep_type = bl2sep[a.miriad.ij2bl(ij[0], ij[1])]
# convert uvw in light-nanoseconds to m, (cosmo_units.c in m/s)
uvw = aa.get_baseline(ij[0], ij[1], src='z') * cosmo_units.c * 1e-9
bins = fringe.gen_frbins(inttime)
mychan = 101  # XXX use this to match frf_filter.py
frp, bins = fringe.aa_to_fr_profile(aa, ij, mychan, bins=bins)
timebins, firs = fringe.frp_to_firs(frp, bins, aa.get_freqs(),
                                    fq0=aa.get_freqs()[mychan])
# chop firs to frequency range of interest
firs = firs[int(opts.chan.split('_')[0]):int(opts.chan.split('_')[1])+1, :]
fir = {(ij[0], ij[1], POL): firs}
fir_conj = {}  # fir for conjugated baselines
for key in fir:
    fir_conj[key] = n.conj(fir[key])
# aa is now subset of freqs, for use in apply_frf later
aa = a.cal.get_aa(opts.cal, afreqs)

# Acquire data
data_dict_v = {}
data_dict_n = {}
flg_dict = {}
conj_dict = {}
stats, lsts, data, flgs = {}, {}, {}, {}
for k in days:
    stats[k], data[k], flgs[k] = capo.miriad.read_files(dsets[k],
                                                        antstr=antstr,
                                                        polstr=POL,
                                                        verbose=True)
    lsts[k] = n.array(stats[k]['lsts'])
    if rmbls:
        print "    Removing baselines:",
        for bl in rmbls:
            data[k].pop(a.miriad.bl2ij(bl), None)
            flgs[k].pop(a.miriad.bl2ij(bl), None)
            print bl,
        print '\n'
    print 'Generating noise for day: ' + str(k)
    if opts.nbls == 'all': nbls = len(data[k].keys())
    else: nbls = int(opts.nbls)
    for bl in data[k].keys()[:nbls]:
    #for bl in data[k].keys()[:nbls][nbls/2:]: # jackknife in bls
        d = n.array(data[k][bl][POL])[:, chans] * jy2T  # extract freq range
        n_ = make_noise(d, stats[k]['cnt'][:, chans], inttime, sdf)
        flg = n.array(flgs[k][bl][POL])[:, chans]  # extract freq range
        key = (k, bl, POL)
        data_dict_v[key] = d
        data_dict_n[key] = n_
        flg_dict[key] = n.logical_not(flg)
        conj_dict[key[1]] = conj[bl]
keys = data_dict_v.keys()
bls_master = []
for key in keys:  # populate list of baselines
    if key[0] == keys[0][0]:
        bls_master.append(key[1])
print 'Baselines:', len(bls_master)

# Align dataset
inds = oqe.lst_align(lsts)
data_dict_v, flg_dict, lsts = oqe.lst_align_data(inds, dsets=data_dict_v,
                                                 wgts=flg_dict, lsts=lsts)
data_dict_n = oqe.lst_align_data(inds, dsets=data_dict_n)[0]
nlst = data_dict_v[keys[0]].shape[0]
# the lsts given is a dictionary with 'even','odd', etc.
lsts = lsts[lsts.keys()[0]] # same lsts for both even and odd
N = len(bls_master)
nbls = N 

# Calculate effective number of counts in the data
cnts = []
for key in stats.keys(): # loop through even and odd
    cnt_full = stats[key]['cnt'][inds[key]]
    cnt_full = cnt_full[:, chans]
    cnts.append(cnt_full)
cnts = n.array(cnts).flatten()
cnts *= n.sqrt(len(stats.keys())) # multiply by sqrt(# of 'datasets')
cnt_eff = 1./n.sqrt(n.ma.masked_invalid(1./cnts**2).mean())

# Fringe-rate filter noise
if opts.frf:
    print 'Fringe-rate-filtering noise'
    for key in data_dict_n:
        size = data_dict_n[key].shape[0]
        # nij = n.tile(data_dict_n[key].flatten(), 3).reshape((3*size,nchan))
        # add padding
        nij = data_dict_n[key].copy()
        wij = n.ones(nij.shape, dtype=bool)
        fir_size = fir.values()[0].shape[1]
        if size < fir_size:
            diff = fir_size - size
            nij = n.pad(nij, (diff, 0), mode='constant', constant_values=0)
            nij = nij[:, diff:]
            wij = n.ones_like(nij)
        else:
            diff = None
        if conj_dict[key[1]] is True:  # apply frf using the conj of data and the conj fir
            nij_frf = fringe_rate_filter(aa, n.conj(nij), wij, ij[0], ij[1],
                                         POL, bins, fir_conj)
        else:
            nij_frf = fringe_rate_filter(aa, nij, wij, ij[0], ij[1],
                                         POL, bins, fir)
        # data_dict_n[key] = nij_frf[size:2*size,:]
        if diff:
            data_dict_n[key] = nij_frf[diff:].copy()
        else:
            data_dict_n[key] = nij_frf.copy()

# Set data
if opts.changeC:
    dsv = DataSet()
    dsn = DataSet()
else:
    dsv = oqe.DataSet(lmode=LMODE)  # just data
    dsn = oqe.DataSet(lmode=LMODE)  # just noise
dsv.set_data(dsets=data_dict_v, conj=conj_dict, wgts=flg_dict)
dsn.set_data(dsets=data_dict_n, conj=conj_dict, wgts=flg_dict)

"""
# NULL TEST: EVEN/ODD
new_dict = {}
for key in data_dict_v:
    if key[0] == 'even':
        new = data_dict_v[key] + data_dict_v[('odd',key[1],key[2])]
    if key[0] == 'odd':
        new = data_dict_v[('even',key[1],key[2])] - data_dict_v[key]
    new_dict[key] = new
dsv.add_data(dsets=new_dict)
"""
"""
# NULL TEST: BASELINES
new_dict = {}
for k,key in enumerate(data_dict_v):
    bl = key[1]
    num = [bb for bb in range(N) if bls_master[bb] == bl][0]
    if num % 2 == 0: # make negative
        new = -data_dict_v[key]
    else: new = data_dict_v[key] # else leave the same
    new_dict[key] = new
dsv.add_data(dsets=new_dict)
"""
"""
# NULL TEST: LST
new_dict = {}
for key in data_dict_v:
    lst1 = n.append(data_dict_v[('even',key[1],key[2])][:nlst/2],data_dict_v[('odd',key[1],key[2])][:nlst/2],axis=0)
    lst2 = n.append(data_dict_v[('even',key[1],key[2])][nlst/2:],data_dict_v[('odd',key[1],key[2])][nlst/2:],axis=0)
    if key[0] == 'even': # first half of LSTs from both even and odd
        new = lst1+lst2
    if key[0] == 'odd': # last half of LSTs from both even and odd
        new = lst1-lst2
    new_dict[key] = new
dsv.add_data(dsets=new_dict)
"""

"""
n_to_save = {}
v_to_save = {}
for kk in data_dict_n:
    n_to_save[str(kk)] = data_dict_n[kk]
    v_to_save[str(kk)] = data_dict_v[kk]
print 'Saving Noise_Dataset.npz and Data_Dataset.npz'
n.savez('Noise_Dataset.npz', **n_to_save)
n.savez('Data_Dataset.npz', **v_to_save)
sys.exit()
"""
if PLOT and False:
    for key in keys:
        p.subplot(311)
        capo.arp.waterfall(dsv.x[key], mode='real')
        p.colorbar()
        p.title('Data x')
        p.subplot(323)
        capo.arp.waterfall(dsv.C(key))
        p.colorbar()
        p.title('C')
        U, S, V = n.linalg.svd(n.conj(dsv.C(key)))
        p.subplot(324)
        p.plot(S)
        p.yscale('log')
        p.grid()
        p.title('Eigenspectrum')
        # p.subplot(324); p.plot(n.einsum('ij,jk',n.diag(S),V).T.real)
        p.subplot(313)
        capo.arp.waterfall(n.dot(dsv.iC(key), dsv.x[key]),
                           mode='real')  # ,drng=6000,mx=3000)
        p.colorbar()
        p.title('C^-1 x')
        # p.suptitle(key)
        p.tight_layout()
        p.show()

# Bootstrap
for boot in xrange(opts.nboot):
    print '\nBootstrap %d / %d' % (boot + 1, opts.nboot)

    # Make groups
    gps = [bls_master[i::NGPS] for i in range(NGPS)]  # no repeated baselines between or within groups
    nbls_g = n.int(n.round(N/NGPS))  # number of baselines per group
    
    # Compute power spectra
    if NGPS > 1:
        pIv, pIn = make_PS(keys, dsv, dsn, grouping=True)
    elif NGPS == 1:
        pIv, pIn =  make_PS(keys, dsv, dsn, grouping=False)
    
    # Bootstrap
    if opts.nboot > 1: # sample cross-multiplications w/replacement
        n.random.seed(int(time.time())) # random seed for bootstrapping
        sample_ind = n.random.choice(pIv.shape[0],pIv.shape[0],replace=True)
        pIv = pIv[sample_ind,:,:] 
        pIn = pIn[sample_ind,:,:]

    print 'pIv =', n.median(pIv.real)
    print 'pIn =', n.median(pIn.real)

    # Save Output
    if len(opts.output) > 0:
        if opts.nboot == 1: outpath = opts.output + '/pspec_noboot.npz'
        else: outpath = opts.output + '/pspec_bootsigloss%04d.npz' % boot
    else:
        if opts.nboot == 1: outpath = '/pspec_noboot.npz'
        else: outpath = '/pspec_bootsigloss%04d.npz' % boot
    print '   Writing ' + outpath
    n.savez(outpath, kpl=kpl, scalar=scalar, lsts=lsts,
            pIv=n.mean(pIv,axis=0), 
            pIn=n.mean(pIn,axis=0), 
            sep=sep_type, uvw=uvw,
            frf_inttime=frf_inttime, inttime=inttime,
            inject_level=INJECT_SIG, freq=fq, afreqs=afreqs,
            cnt_eff=cnt_eff, nbls=nbls, ngps=NGPS, nbls_g=nbls_g,
            cmd=' '.join(sys.argv))
