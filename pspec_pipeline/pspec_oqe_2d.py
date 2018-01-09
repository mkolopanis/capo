#! /usr/bin/env python
"""Compute 2-Dimensional power spectrum from uv files.

Computes the (k// vs t) power spectrum from uv files.
Bootstraps computed over subsets of Baselines in UV files.
"""

import aipy as a
import numpy as n
import pylab as p
import glob
import optparse
import sys
import random
from capo import zsa, oqe, cosmo_units, frf_conv as fringe
import capo

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
opts, args = o.parse_args(sys.argv[1:])

# Basic parameters
random.seed(0)  # for oqe.py (eor generator)
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

    def iC(self, k, t=None, rcond=1e-12):
        """Regularize covariance before inverting."""
        assert(t is None)
        if k not in self._iC.keys(k):
            C = self.C(k)
            # CHANGE C HERE ###
            # OPTION 1: identity multiplication
            C = C * n.identity(len(C))
            # OPTION 2: identity addition
            # C = C + n.identity(len(C))*n.trace(C)*float(opts.mode_num)
            # C = C + n.identity(len(C))*int(opts.mode_num)
            # OPTION 3: multiplication by identity + 2 diagonals

            # C2 = n.zeros_like(C)
            # for i in range(C.shape[0]):
            #     for j in range(C.shape[1]):
            #         if n.abs(j-i) <= 1: C2[i,j] = 1.0
            # C = C * C2

            # OPTION 4: use 'average C',
            # where it is computed based on average data
            # avgx = self.x[k] - self.avgx
            # C = oqe.cov(avgx, self.w[k])
            # OPTION 5: subtract average C off of C
            # C = C - self.avgC
            # OPTION 6: add identity at strength of median(eigenvalues)
            # U,S,V = n.linalg.svd(C.conj())
            # C = C + n.identity(len(C)) * n.median(S)
            # svd
            U, S, V = n.linalg.svd(C.conj())  # conj in advance of next step
            iS = 1./S
            # OR CHANGE iS HERE ###
            # iS[:3] = 0.0
            # iS[int(opts.mode_num):] = 1.0
            self.set_iC({k: n.einsum('ij,j,jk', V.T, iS, U.T)})
        return self._iC[k]

    def set_data(self, dsets, wgts=None, conj=None):
        """Set data inside of object, also computes average over bls."""
        if type(dsets.values()[0]) == dict:
            dsets, wgts = self.flatten_data(dsets), self.flatten_data(wgts)
        self.x, self.w = {}, {}
        avgx = []
        avgC = []
        for k in dsets:
            self.x[k] = dsets[k].T
            try: self.w[k] = wgts[k].T
            except(TypeError): self.w[k] = n.ones_like(self.x[k])
            try:
                if conj[k[1]]: self.x[k] = n.conj(self.x[k])
            except(TypeError, KeyError): pass
            try:
                avgx.append(self.x[k])
                avgC.append(oqe.cov(self.x[k], self.w[k]))
            except(TypeError, KeyError): pass
        self.avgx = n.average(avgx, axis=0)
        self.avgC = n.average(avgC, axis=0)


def complex_noise(size, noiselev):
    """Generate complex noise of given size and noiselevel."""
    # if noiselev <= 0 or n.isinf(noiselev):
    #     return n.zeros(size)
    noise_real = n.random.normal(size=size, scale=noiselev)/n.sqrt(2)
    noise_imag = n.random.normal(size=size, scale=noiselev)/n.sqrt(2)
    noise = noise_real + 1j*noise_imag
    return noise


def make_noise(d, cnt, inttime, df): #, freqs, jy2T=None):
    """Create noise with T_rms matching data from T_rcvr and uv['cnt']."""
    #if jy2T is None:
    #    jy2T = capo.pspec.jy2T(freqs)
    Tsys = 180. * n.power(afreqs/0.18, -2.55) + opts.Trcvr  # system temp in K
    Tsys *= 1e3  # system temp in mK
    Trms = Tsys/n.sqrt(df * 1e9 * inttime * cnt * NPOL)  # convert sdf to Hz
                        # Bchan, inttime, counts (times 2 if even&odd), Npol
    Vrms = Trms#/jy2T  # jy2T is in units of mK/Jy
    # The following transposes are to create noise correlated in time not
    # frequency. Maybe there is a better way to do it?
    # The masking and filling is to be able to parallelize the noise draws
    # Setting the mask back later and filling zeros out where Vrms ~ inf or < 0
    size = Vrms.shape[0]
    #if opts.frf: # triple size
    #    Vrms = n.repeat(Vrms, 3, axis=0)
    Vrms = n.ma.masked_invalid(Vrms)
    Vrms.mask = n.ma.mask_or(Vrms.mask, Vrms.filled() < 0)
    n.ma.set_fill_value(Vrms, 1e-20)
    noise = n.array([complex_noise(v1.shape, v1.filled()) for v1 in Vrms.T]).T
    noise = n.ma.masked_array(noise)
    noise.mask = Vrms.mask
    n.ma.set_fill_value(noise, 0 + 0j)
    noise = noise.filled()
    #wij = n.ones(noise.shape, dtype=bool) # XXX flags are all true (times,freqs)
    #if opts.frf: # FRF noise
    #    noise = fringe_rate_filter(aa, noise, wij, ij[0], ij[1], POL, bins, fir)
    #noise = noise[int(size):2*int(size),:]
    #noise.shape = d.shape
    return noise


def fringe_rate_filter(aa, dij, wij, i, j, pol, bins, fir):
    """ Apply frf."""
    _d, _w, _, _ = fringe.apply_frf(aa, dij, wij, i, j, pol=pol, bins=bins, firs=fir)
    return _d


def make_eor(shape):  # Create and fringe rate filter noise
    """Generate White Noise and Apply FRF."""
    dij = oqe.noise(size=shape)
    return dij


def make_PS(keys, dsv, dsn, dse, dsr, dss, dse_Cr, dsv_Cr,dsve_Cr,grouping=True):
    """Use OQE formalism to generate power spectrum.

    Output weighted and identity weightings.
    """
    if grouping:
        newkeys, dsCv = dsv.group_data(keys, gps)
        newkeys, dsIv = dsv.group_data(keys, gps, use_cov=False)
        newkeys, dsCn = dsn.group_data(keys, gps)
        newkeys, dsIn = dsn.group_data(keys, gps, use_cov=False)
        newkeys, dsCe = dse.group_data(keys, gps)
        newkeys, dsIe = dse.group_data(keys, gps, use_cov=False)
        newkeys, dsCr = dsr.group_data(keys, gps)
        newkeys, dsIr = dsr.group_data(keys, gps, use_cov=False)
        newkeys, dsCs = dss.group_data(keys, gps)
        newkeys, dsIs = dss.group_data(keys, gps, use_cov=False)
        newkeys, dsCe_Cr = dse_Cr.group_data(keys, gps)
        newkeys, dsCv_Cr = dsv_Cr.group_data(keys, gps)
        newkeys, dsCve_Cr = dsve_Cr.group_data(keys, gps, use_cov=False)
    else:  # no groups (slower)
        newkeys = keys
        dsIv, dsCv = dsv, dsv  # identity and covariance case dataset is the same
        dsIn, dsCn = dsn, dsn
        dsIe, dsCe = dse, dse
        dsIr, dsCr = dsr, dsr
        dsIs, dsCs = dss, dss
        dsCe_Cr = dse_Cr
        dsCv_Cr = dsv_Cr
        dsCve_Cr = dsve_Cr
    pCvs = []; pIvs = []
    pCns = []; pIns = []
    pCes = []; pIes = []
    pCrs = []; pIrs = []
    pCss = []; pIss = []
    pCes_Cr = []
    pCvs_Cr = []
    pCves_Cr = []
    for k, key1 in enumerate(newkeys):
        if k == 1 and len(newkeys) == 2: # NGPS = 1 (skip 'odd' with 'even' if we already did 'even' with 'odd')
            continue
        #print '   ',k+1,'/',len(newkeys)
        for key2 in newkeys[k:]:
            if len(newkeys) > 2 and (key1[0] == key2[0] or key1[1] == key2[1]): # NGPS > 1
                continue
            if key1[0] == key2[0]: # don't do 'even' with 'even', for example
                continue
            else:
                FCv = dsCv.get_F(key1, key2, cov_flagging=False)
                FIv = dsIv.get_F(key1, key2, use_cov=False, cov_flagging=False)
                qCv = dsCv.q_hat(key1, key2, cov_flagging=False)
                qIv = dsIv.q_hat(key1, key2, use_cov=False, cov_flagging=False)
                MCv, WCv = dsCv.get_MW(FCv, mode=opts.weight)
                MIv, WIv = dsIv.get_MW(FIv, mode='I')
                pCv = dsCv.p_hat(MCv, qCv, scalar=scalar)
                pIv = dsIv.p_hat(MIv, qIv, scalar=scalar)
                pCvs.append(pCv)
                pIvs.append(pIv)

                FCn = dsCn.get_F(key1, key2, cov_flagging=False)
                FIn = dsIn.get_F(key1, key2, use_cov=False, cov_flagging=False)
                qCn = dsCn.q_hat(key1, key2, cov_flagging=False)
                qIn = dsIn.q_hat(key1, key2, use_cov=False, cov_flagging=False)
                MCn, WCn = dsCn.get_MW(FCn, mode=opts.weight)
                MIn, WIn = dsIn.get_MW(FIn, mode='I')
                pCn = dsCv.p_hat(MCn, qCn, scalar=scalar)
                pIn = dsIv.p_hat(MIn, qIn, scalar=scalar)
                pCns.append(pCn)
                pIns.append(pIn)

                FCe = dsCe.get_F(key1, key2, cov_flagging=False)
                FIe = dsIe.get_F(key1, key2, use_cov=False, cov_flagging=False)
                qCe = dsCe.q_hat(key1, key2, cov_flagging=False)
                qIe = dsIe.q_hat(key1, key2, use_cov=False, cov_flagging=False)
                MCe, WCe = dsCe.get_MW(FCe, mode=opts.weight)
                MIe, WIe = dsIe.get_MW(FIe, mode='I')
                pCe = dsCe.p_hat(MCe, qCe, scalar=scalar)
                pIe = dsIe.p_hat(MIe, qIe, scalar=scalar)
                pCes.append(pCe)
                pIes.append(pIe)

                FCr = dsCr.get_F(key1, key2, cov_flagging=False)
                FIr = dsIr.get_F(key1, key2, use_cov=False, cov_flagging=False)
                qCr = dsCr.q_hat(key1, key2, cov_flagging=False)
                qIr = dsIr.q_hat(key1, key2, use_cov=False, cov_flagging=False)
                MCr, WCr = dsCr.get_MW(FCr, mode=opts.weight)
                MIr, WIr = dsIr.get_MW(FIr, mode='I')
                pCr = dsCr.p_hat(MCr, qCr, scalar=scalar)
                pIr = dsIr.p_hat(MIr, qIr, scalar=scalar)
                pCrs.append(pCr)
                pIrs.append(pIr)

                FCs = dsCs.get_F(key1, key2, cov_flagging=False)
                FIs = dsIs.get_F(key1, key2, use_cov=False, cov_flagging=False)
                qCs = dsCs.q_hat(key1, key2, cov_flagging=False)
                qIs = dsIs.q_hat(key1, key2, use_cov=False, cov_flagging=False)
                MCs, WCs = dsCs.get_MW(FCs, mode=opts.weight)
                MIs, WIs = dsIs.get_MW(FIs, mode='I')
                pCs = dsCs.p_hat(MCs, qCs, scalar=scalar)
                pIs = dsIs.p_hat(MIs, qIs, scalar=scalar)
                pCss.append(pCs)
                pIss.append(pIs)

                FCe_Cr = dsCe_Cr.get_F(key1, key2, cov_flagging=False)
                qCe_Cr = dsCe_Cr.q_hat(key1, key2, cov_flagging=False)
                MCe_Cr, WCe_Cr = dsCe_Cr.get_MW(FCe_Cr, mode=opts.weight)
                pCe_Cr = dsCe_Cr.p_hat(MCe_Cr, qCe_Cr, scalar=scalar)
                pCes_Cr.append(pCe_Cr)

                FCv_Cr = dsCv_Cr.get_F(key1, key2, cov_flagging=False)
                qCv_Cr = dsCv_Cr.q_hat(key1, key2, cov_flagging=False)
                MCv_Cr, WCv_Cr = dsCv_Cr.get_MW(FCv_Cr, mode=opts.weight)
                pCv_Cr = dsCv_Cr.p_hat(MCv_Cr, qCv_Cr, scalar=scalar)
                pCvs_Cr.append(pCv_Cr)

                FCve_Cr = dsCve_Cr.get_F(key1, key2, cov_flagging=False)
                qCve_Cr = dsCve_Cr.q_hat(key1, key2, cov_flagging=False)
                MCve_Cr, WCve_Cr = dsCve_Cr.get_MW(FCve_Cr, mode=opts.weight)
                pCve_Cr = dsCve_Cr.p_hat(MCve_Cr, qCve_Cr, scalar=scalar)
                pCves_Cr.append(pCve_Cr)

    if PLOT:
        p.subplot(121)
        capo.arp.waterfall(FCv, drng=4)
        p.title('FC')
        p.subplot(122)
        capo.arp.waterfall(FIv, drng=4)
        p.title('FI')
        p.show()
    if PLOT:
        p.subplot(411)
        capo.arp.waterfall(qCv, mode='real')
        p.colorbar(shrink=.5)
        p.title('qC')
        p.subplot(412)
        capo.arp.waterfall(pCv, mode='real')
        p.colorbar(shrink=.5)
        p.title('pC')
        p.subplot(413)
        capo.arp.waterfall(qIv, mode='real')
        p.colorbar(shrink=.5)
        p.title('qI')
        p.subplot(414)
        capo.arp.waterfall(pIv, mode='real')
        p.colorbar(shrink=.5)
        p.title('pI')
        p.show()
    if PLOT:
        p.plot(kpl, n.average(pCv.real, axis=1), 'b.-', label='pC')
        p.plot(kpl, n.average(pIv.real, axis=1), 'k.-', label='pI')
        p.legend()
        p.show()
    return n.array(pCvs), n.array(pIvs), n.array(pCns), n.array(pIns), n.array(pCes), n.array(pIes), n.array(pCrs), n.array(pIrs), n.array(pCss), n.array(pIss), n.array(pCes_Cr), n.array(pCvs_Cr), n.array(pCves_Cr)

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
aa = a.cal.get_aa(opts.cal, freqs) # all freqs
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
    if opts.frf: bm *= 1.3 # correction factor for FRF omega_pp = .32/.24 = 1.3
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
s,d,f = capo.miriad.read_files([dsets[days[0]][0]], antstr=antstr, polstr=POL) # read first file
ij = d.keys()[0] # use first baseline
if blconj[a.miriad.ij2bl(ij[0], ij[1])]:
    # makes sure FRP will be the same whether bl is a conjugated one or not
    if ij[0] < ij[1]:
        temp = (ij[1], ij[0])
        ij = temp
sep_type = bl2sep[a.miriad.ij2bl(ij[0], ij[1])]
# convert uvw in light-nanoseconds to m, (cosmo_units.c in m/s)
uvw = aa.get_baseline(ij[0], ij[1], src='z') * cosmo_units.c * 1e-9
bins = fringe.gen_frbins(inttime)
mychan = 101 # XXX use this to match frf_filter.py
frp, bins = fringe.aa_to_fr_profile(aa, ij, mychan, bins=bins)
timebins, firs = fringe.frp_to_firs(frp, bins, aa.get_freqs(),
                                    fq0=aa.get_freqs()[mychan])
firs = firs[int(opts.chan.split('_')[0]):int(opts.chan.split('_')[1])+1,:] # chop firs to frequency range of interest
fir = {(ij[0], ij[1], POL): firs}
fir_conj = {} # fir for conjugated baselines
for key in fir:
    fir_conj[key] = n.conj(fir[key])
aa = a.cal.get_aa(opts.cal, afreqs) # aa is now subset of freqs, for use in apply_frf later

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
    for bl in data[k]:
        d = n.array(data[k][bl][POL])[:, chans] * jy2T # extract freq range
        n_ = make_noise(d, stats[k]['cnt'][:, chans], inttime, sdf)
        flg = n.array(flgs[k][bl][POL])[:, chans] # extract freq range
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

# Save some information
cnt_full = stats[stats.keys()[0]]['cnt'][inds[stats.keys()[0]]]
cnt_full = cnt_full[:, chans]
lsts = lsts[lsts.keys()[0]]
# calculate the effective number of counts used in the data
cnt_eff = 1./n.sqrt(n.ma.masked_invalid(1./cnt_full**2).mean())
# calculate the effective number of baselines given grouping:
N = len(bls_master)
nbls = N

# Fringe-rate filter noise
if opts.frf:
    print 'Fringe-rate-filtering noise'
    for key in data_dict_n:
        size = data_dict_n[key].shape[0]
        #nij = n.tile(data_dict_n[key].flatten(), 3).reshape((3*size,nchan)) # add padding
        nij = data_dict_n[key].copy()
        wij = n.ones(nij.shape, dtype=bool)
        fir_size = fir.values()[0].shape[1]
        if size < fir_size:
            diff = fir_size - size
            nij = n.pad(nij, (diff,0), mode='constant', constant_values=0)
            nij = nij[:,diff:]
            wij = n.ones_like(nij)
        else:
            diff = None
        if conj_dict[key[1]] is True: # apply frf using the conj of data and the conj fir
            nij_frf = fringe_rate_filter(aa, n.conj(nij), wij, ij[0], ij[1], POL, bins, fir_conj)
        else:
            nij_frf = fringe_rate_filter(aa, nij, wij, ij[0], ij[1], POL, bins, fir)
        #data_dict_n[key] = nij_frf[size:2*size,:]
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
        U,S,V = n.linalg.svd(n.conj(dsv.C(key)))
        p.subplot(324); p.plot(S); p.yscale('log'); p.grid()
        p.title('Eigenspectrum')
        # p.subplot(324); p.plot(n.einsum('ij,jk',n.diag(S),V).T.real)
        p.subplot(313)
        capo.arp.waterfall(n.dot(dsv.iC(key), dsv.x[key]),
                           mode='real')  # ,drng=6000,mx=3000)
        p.colorbar()
        p.title('C^-1 x')
        #p.suptitle(key)
        p.tight_layout()
        p.show()

# Bootstrap
for boot in xrange(opts.nboot):
    print '\nBootstrap %d / %d' % (boot + 1, opts.nboot)

    # Create fake eor signal
    print '  INJECTING SIMULATED SIGNAL @ LEVEL', INJECT_SIG
    eij = make_eor((nlst*3, nchan))
    size = nlst
    #eij = n.tile(eij.flatten(), 3).reshape((3*nlst,nchan)) # add padding
    wij = n.ones(eij.shape, dtype=bool)
    eij_frf = fringe_rate_filter(aa, eij, wij, ij[0], ij[1], POL, bins, fir)
    eij_conj_frf = fringe_rate_filter(aa, n.conj(eij), wij, ij[0], ij[1], POL, bins, fir_conj) # conjugated eor with conjugated FIR
    if opts.frf: # double frf eor
        eij_frf = fringe_rate_filter(aa, eij_frf, wij, ij[0], ij[1], POL, bins, fir)
        eij_conj_frf = fringe_rate_filter(aa, eij_conj_frf, wij, ij[0], ij[1], POL, bins, fir_conj)

    eor = eij_frf[size:2*size,:]*n.abs(INJECT_SIG)
    eor_conj = eij_conj_frf[size:2*size,:]*n.abs(INJECT_SIG)
    data_dict_r = {}
    data_dict_e = {}
    data_dict_s = {}
    for key in data_dict_v:
        if conj_dict[key[1]] is True:
            eorinject = eor_conj
        else:
            eorinject = eor
        # track eor in separate dict
        data_dict_e[key] = eorinject
        # add injected signal to data
        data_dict_r[key] = data_dict_v[key].copy() + eorinject
        # add injected signal to noise
        data_dict_s[key] = data_dict_n[key].copy() + eorinject
        if INJECT_SIG < 0.: # negative injects
            if key[0] == 'even': # XXX hard-coded way to subtract e for negative injects
                data_dict_e[key] = -eorinject
                data_dict_r[key] = data_dict_v[key].copy() - eorinject
                data_dict_s[key] = data_dict_n[key].copy() - eorinject

    # Set data
    if opts.changeC:
        dsr = DataSet()
        dse = DataSet()
        dss = DataSet()
        dse_Cr = DataSet()
        dsv_Cr = DataSet()
        dsve_Cr = DataSet()
    else:
        dsr = oqe.DataSet(lmode=LMODE)  # data + eor
        dse = oqe.DataSet(lmode=LMODE)  # just eor
        dss = oqe.DataSet(lmode=LMODE)  # noise + eor
        dse_Cr = oqe.DataSet(lmode=LMODE) # just eor, but with C_r
        dsv_Cr = oqe.DataSet(lmode=LMODE) # just data, but with C_r
        dsve_Cr = oqe.DataSet(lmode=LMODE) # cross-term x with e, with C_r

    dsr.set_data(dsets=data_dict_r, conj=conj_dict, wgts=flg_dict)
    dse.set_data(dsets=data_dict_e, conj=conj_dict, wgts=flg_dict)
    dss.set_data(dsets=data_dict_s, conj=conj_dict, wgts=flg_dict)
    dse_Cr.set_data(dsets=data_dict_e, conj=conj_dict, wgts=flg_dict)
    dsv_Cr.set_data(dsets=data_dict_v, conj=conj_dict, wgts=flg_dict)

    # Edit data for cross-term
    data_cross = {}
    for key in data_dict_v.keys():
        if key[0] == 'even': # keep as data
            data_cross[key] = data_dict_v[key]
        elif key[0] == 'odd': # change to eor data
            data_cross[key] = data_dict_e[key]
    dsve_Cr.set_data(dsets=data_cross, conj=conj_dict, wgts=flg_dict)

    # Over-write C's for some of the datasets
    for key in data_dict_e.keys():
        C_r = dsr.C(key)
        U,S,V = n.linalg.svd(C_r.conj())
        iS = 1./S
        dse_Cr.set_iC({key:n.einsum('ij,j,jk', V.T, iS, U.T)})
        dsv_Cr.set_iC({key:n.einsum('ij,j,jk', V.T, iS, U.T)})
        dsve_Cr.set_iC({key:n.einsum('ij,j,jk', V.T, iS, U.T)})

    # Make groups
    if opts.nboot > 1: # sample baselines w/replacement
        gps = dsv.gen_gps(bls_master, ngps=NGPS)
        nbls_g = n.int(n.round(N/NGPS))  # number of baselines per group
    elif opts.nboot == 1: # no sampling w/replacement (no bootstrapping)
        gps = [bls_master[i::NGPS] for i in range(NGPS)] # no repeated baselines between or within groups
        nbls_g = n.int(n.round(N/NGPS)) # number of baselines per group
        nbls = nbls_g # over-ride nbls for proper sensitivity calculation later
        NGPS = 1 # over-ride NGPS for proper sensitivity calculation later

    # Compute power spectra
    pCv, pIv, pCn, pIn, pCe, pIe, pCr, pIr, pCs, pIs, pCe_Cr, pCv_Cr, pCve_Cr = make_PS(keys, dsv, dsn, dse, dsr, dss, dse_Cr, dsv_Cr, dsve_Cr, grouping=True)

    print '     Data:         pCv =', n.median(pCv.real),
    print 'pIv =', n.median(pIv.real)
    print '     EoR:          pCe =', n.median(pCe.real),
    print 'pIe =', n.median(pIe.real)
    print '     Noise:        pCn =', n.median(pCn.real),
    print 'pIn =', n.median(pIn.real)
    print '     Data + EoR:   pCr =', n.median(pCr.real),
    print 'pIr =', n.median(pIr.real)
    print '     Noise + EoR:  pCs =', n.median(pCs.real),
    print 'pIs =', n.median(pIs.real)

    print '       Signal Loss Data  ~ pIe/(pCr-pCv) =',
    print n.abs(n.median(pIe.real)) / n.abs(n.median(pCr.real) - n.median(pCv))
    print '       Signal Loss Noise ~ pIe/(pCs-pCn) =',
    print n.abs(n.median(pIe.real)) / n.abs(n.median(pCs.real) - n.median(pCn))

    if PLOT:
        p.plot(kpl, n.average(pCr.real, axis=1), 'b.-')
        p.plot(kpl, n.average(pIr.real, axis=1), 'k.-')
        p.title('Data + EoR')
        p.show()

    # Save Output
    if len(opts.output) > 0:
        if opts.nboot == 1: outpath = opts.output + '/pspec_noboot.npz'
        else: outpath = opts.output + '/pspec_bootsigloss%04d.npz' % boot
    else:
        if opts.nboot == 1: outpath = '/pspec_noboot.npz'
        else: outpath = '/pspec_bootsigloss%04d.npz' % boot
    print '   Writing ' + outpath
    n.savez(outpath, kpl=kpl, scalar=scalar, lsts=lsts,
            pCr=pCr, pIr=pIr, pCv=pCv, pIv=pIv, pCe=pCe, pCe_Cr=pCe_Cr,
            pIe=pIe, pCn=pCn, pIn=pIn, pCs=pCs, pIs=pIs, pCv_Cr=pCv_Cr,
            pCve_Cr=pCve_Cr, sep=sep_type, uvw=uvw,
            frf_inttime=frf_inttime, inttime=inttime,
            inject_level=INJECT_SIG, freq=fq, afreqs=afreqs,
            cnt_eff=cnt_eff, nbls=nbls, ngps=NGPS, nbls_g=nbls_g,
            cmd=' '.join(sys.argv))
