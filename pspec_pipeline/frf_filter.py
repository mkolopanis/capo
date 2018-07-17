#! /usr/bin/env python

import aipy as a
import capo.arp as arp
from capo import fringe
import capo.zsa as zsa
import capo as C
import numpy as n
import pylab as p
import sys
import os
import optparse


def noise_equivalent_bandwidth(t, H):
    # print 'NEQB amp:', 1./n.sum(n.abs(H)*n.diff(t)[0])**2
    return 1./n.max(n.abs(H))**2 * n.sum(n.abs(H)**2)*n.diff(t)[0]


o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, ant=True, pol=True)
o.add_option('-c', '--chan', type='int', default=None,
            help='Channel to use in making the fringe rate filter.')
o.add_option('--alietal', action='store_true',
             help='Uses normalization for alietal frf,(default=False)')
o.add_option('--outpath', action='store',
             help='Output path to write to.')

opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
if opts.outpath:
    outfile = opts.outpath + '/' + args[0] + 'L'
else:
    outfile = args[0]+'L'

if os.path.exists(outfile):
    print 'File exists:', outfile
    sys.exit(0)
nants = uv['nants']
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
pol = a.miriad.pol2str[uv['pol']]

(uvw, t1, (i, j)), d = uv.read()
(uvw, t2, (i, j)), d = uv.read()
while t1 == t2:
    (uvw, t2, (i, j)), d = uv.read()
inttime = (t2-t1) * (3600*24)
#  manually calculate the inttime so frf time bins match data
del(uv)

# set channel to calculate the FRF
if opts.chan != None:
    mychan = int(opts.chan)
else:
    mychan = 101

# Get only the antennas of interest
sep2ij, blconj, bl2sep = zsa.grid2ij(aa.ant_layout)
print "Looking for baselines matching", opts.ant
print "Using normalization for old FRF: ", opts.alietal
ants = [b[0] for b in a.scripting.parse_ants(opts.ant, nants)]
seps = [bl2sep[b] for b in ants]
seps = n.unique(seps)
print 'These are the separations that we are going to use:', seps
print 'This is the channel we are using to build the frf: ', mychan
print 'Current inttime use for gen_frbins: ', inttime

# Get the fir filters for the separation used
bins = fringe.gen_frbins(inttime)
firs = {}
for sep in seps:
    c = 0  # baseline indices
    while True:
        ij = map(int, sep2ij[sep].split(',')[c].split('_'))
        bl = a.miriad.ij2bl(*ij)
        if blconj[bl]:
            c += 1
        else:
            break   # find when conjugation isn't necessary
    frp, bins = fringe.aa_to_fr_profile(aa, ij, mychan, pol=pol,
                                        bins=bins, alietal=opts.alietal)
    timebins, firs[sep] = fringe.frp_to_firs(frp, bins, aa.get_freqs(),
                                             fq0=aa.get_freqs()[mychan])

frps_neqb, frp_freqs = fringe.fir_to_frp(firs[sep], tbins=timebins)
# 1.2 multiplication, Calibrated against square filter of
# equivalent size.
frf_inttime = 1.2/noise_equivalent_bandwidth(frp_freqs, frps_neqb[mychan])

print 'The FRF has a noise equivalent bandwidth [s]: ', frf_inttime

baselines = ''.join(sep2ij[sep] for sep in seps)
times, data, flags = C.miriad.read_files(args, baselines, pol, verbose=True)
# new way of get_dict_of_uv_data
# jds = times['times']
# lsts = [ aa.sidereal_time() for k in map(aa.set_jultime(), jds) ]
lsts = times['lsts']
lst_order = n.argsort(lsts)  # data is not always read in LST order!
lsts = n.array(lsts)[lst_order]
times['times'] = times['times'][lst_order]
for bl in data:  # orders data and flags correctly by LST
    for pol in data[bl]:
        data[bl][pol] = data[bl][pol][lst_order]
        flags[bl][pol] = flags[bl][pol][lst_order]
_d = {}
_w = {}
bins = fringe.gen_frbins(inttime)
for bl in data.keys():  # bl format is (i,j) in data keys
    if bl not in _d.keys():
        _d[bl], _w[bl] = {}, {}
    # get filter which is baseline dependent
    m_bl = a.miriad.ij2bl(bl[0], bl[1])  # miriad bl
    sep = bl2sep[m_bl]
    i, j = bl
    # fir = firs[sep]
    # if blconj[m_bl]: fir = n.conj(fir)
    # print map(int, a.miriad.bl2ij(m_bl)), sep, blconj[m_bl]
    print bl
    for pol in data[bl].keys():
        if blconj[m_bl]:
            fir = {(i, j, pol): n.conj(firs[sep])}  # conjugate fir if needed
        else:
            fir = {(i, j, pol): firs[sep]}
        if pol not in _d[bl].keys():
            _d[bl][pol], _w[bl][pol] = {}, {}
        # _d[bl][pol] = n.zeros_like(data[bl][pol])
        # _w[bl][pol] = n.zeros_like(data[bl][pol])
        dij, wij = data[bl][pol], n.logical_not(flags[bl][pol])
        _d[bl][pol], _w[bl][pol], _, _ = fringe.apply_frf(aa, dij, wij, i, j, pol=pol, bins=bins, firs=fir, alietal=opts.alietal)


def mfunc(uv, p, d, f):
    uvw, t, (i, j) = p
    index = n.where(times['times'] == t)
    # bl = a.miriad.ij2bl(i,j)
    bl = (i, j)
    pol = a.miriad.pol2str[uv['pol']]
    try:
        # The flags are the same flags as input.
        d_ = _d[bl][pol][index, :].reshape(203)
    except(KeyError):
        return p, None, None
    else:
        return p, d_, f


for filename in args:
    if opts.outpath:
        # split on last directory if writing to new directory
        new_name = filename.split('/')[-1]
        outfile = opts.outpath + '/' + new_name + 'L'
    else:
        outfile = filename+'L'
    uvi = a.miriad.UV(filename)
    a.scripting.uv_selector(uvi, ants=baselines)
    uvo = a.miriad.UV(outfile, status='new')
    print 'Writing %s' % (outfile)
    uvo.init_from_uv(uvi, override={'inttime': inttime})
    uvo.add_var('FRF_NEBW', 'r')
    uvo.__setitem__('FRF_NEBW', frf_inttime)
    uvo.pipe(uvi, mfunc=mfunc, append2hist=' '.join(sys.argv)+'\n', raw=True)
