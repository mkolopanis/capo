#! /usr/bin/env python
"""
Create waterfall plot of T_RMS from Miriad UV files.

    Here is some info the help with the legend:
    avg: Data is averaged AFTER  multiplying x.conj * x
    blavg: Data averaged over baseling BEFORE x.conj * x
    e-0: Even - Odd subtraction of data
    var: T_RMS constructed from using the var and cnt info from uv files
    GSM: Global Sky Moden
"""

import aipy as a
import numpy as n
import pylab as p
import sys
import optparse
import glob
import ephem
import capo
import numpy as np

o = optparse.OptionParser()
o.set_usage("plot_trms.py [options]")
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True,
                                 cal=True, cmap=True, chan=True)
o.add_option('--plot', dest='plot', default=False, action='store_true',
             help='Outputs plot to X-Window and saves plot')
o.add_option('--output', type='string', default='',
             help='output directory for image files (default "")')
o.add_option('--vline', default=False, action='store_true',
             help='Emphasizes chosen channel range')
o.add_option('--band', default='0_202', action='store',
             help='Channels from which to Calculate full band Covariance')
o.add_option('--inttime', default=None, action='store', type=float,
             help="Specify inttime for data, else use uv['inttime']")
o.add_option('--Tsky', default=180, type=float,
             help=('rough estimate of the sky temperature '
                   'in this data @ 180MHz. [default 180K]'))
o.add_option('--Trcvr', default=100, type=float,
             help='receiver termp [default 100K]')
o.add_option('--rmbls', type=str,
             help=('list of baselines to NOT use, comma delim, '
                   'integers output by a.miriad.ij2bl'))
o.add_option('--varinttime', type=float,
             help='integration time of the variances (s, default=inttime)')
opts, args = o.parse_args(sys.argv[1:])


POL = opts.pol
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
aa = a.cal.get_aa(opts.cal, .1, .1, 1)
freqs = None
dsets = {
    'even': [x for x in args if 'even' in x],
    'odd': [x for x in args if 'odd' in x]
}

print 'Number of even data sets: {0:d}'.format(len(dsets['even']))
print 'Number of odd data sets: {0:d}'.format(len(dsets['odd']))

uv = a.miriad.UV(dsets.values()[0][0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
sdf = uv['sdf']
print "sdf = ", sdf
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
band_chans = a.scripting.parse_chans(opts.band, uv['nchan'])

if opts.inttime is not None:
    inttime = opts.inttime
else:
    inttime = uv['inttime']
    # * 4  # XXX hack for *E files that have inttime set incorrectly
print 'inttime', inttime
del(uv)

if opts.varinttime is None:
    opts.varinttime = inttime

afreqs = freqs.take(chans)
allfreqs = freqs.take(band_chans)
nchan_band = len(band_chans)
nchan = chans.size
fq = n.average(afreqs)
z = capo.pspec.f2z(fq)

aa = a.cal.get_aa(opts.cal, allfreqs)
bls, conj = capo.red.group_redundant_bls(aa.ant_layout)
jy2T = capo.pspec.jy2T(allfreqs)

antstr = 'cross'
stats = {}
lsts, cnt1, var1, data, flgs = {}, {}, {}, {}, {}
days = dsets.keys()
for k in days:
    stats[k], data[k], flgs[k] = capo.miriad.read_files(dsets[k],
                                                        antstr=antstr,
                                                        polstr=POL,
                                                        verbose=True)
    if rmbls:
        print "    Removing baselines:",
        for bl in rmbls:
            data[k].pop(a.miriad.bl2ij(bl), None)
            flgs[k].pop(a.miriad.bl2ij(bl), None)
            print bl,
        print '\n'
    lsts[k] = stats[k]['lsts']
    cnt1[k] = stats[k]['cnt']
    var1[k] = stats[k]['var']
    print data[k].keys()

# if False:
if True:
    # Align data sets in LST
    print [lsts[k][0] for k in days]
    lstmax = max([lsts[k][0] for k in days])
    for k in days:
        # print k
        for i in xrange(len(lsts[k])):
            # allow for small numerical differences (which shouldn't exist!)
            if lsts[k][i] >= lstmax - .001:
                break
        lsts[k] = lsts[k][i:]
        # cnt1[k], var1[k] , cnt1[k][i:], var1[k][i:]
        for bl in data[k]:
            data[k][bl] = data[k][bl][POL][i:]
            flgs[k][bl] = flgs[k][bl][POL][i:]
            cnt1[k], var1[k] = cnt1[k][i:], var1[k][i:]
    # print [len(lsts[k]) for k in days]
    j = min([len(lsts[k]) for k in days])
    for k in days:
        lsts[k] = lsts[k][:j]
        # cnt1[k], var1[k] = , cnt1[k][:j], var1[k][:j]
        for bl in data[k]:
            data[k][bl] = n.array(data[k][bl][:j])
            flgs[k][bl] = n.array(flgs[k][bl][:j])
            cnt1[k], var1[k] = cnt1[k][:j], var1[k][:j]

lsts = lsts.values()[0]

if len(lsts) == 0:
        print('No data to plot.')
        sys.exit(0)
nlst = len(lsts)
# dlst= n.ceil(nlst/13.)
dlst = 1
x = {}
cnt_o, var_o = {}, {},
cnt_e, var_e = {}, {},
trms_df = {}

# collect count and vars from UV files to create TRMS model
if set(['even', 'odd']) == set(days):
    for bl in data['even']:
        d = (n.copy(data['even'][bl][:, band_chans]
             - data['odd'][bl][:, band_chans]))
        c_e = n.copy(n.array(cnt1['even'])[:, band_chans])
        v_e = n.copy(n.array(var1['even'])[:, band_chans])
        c_o = n.copy(n.array(cnt1['odd'])[:, band_chans])
        v_o = n.copy(n.array(var1['odd'])[:, band_chans])
        if conj[bl]:
            d = n.conj(d)
        x[bl] = n.transpose(d, [1, 0])[band_chans, ::dlst]
        v_e = n.ma.masked_where(c_e < 1, v_e)
        v_o = n.ma.masked_where(c_o < 1, v_o)
        c_e = n.ma.masked_where(c_e < 1, c_e)
        c_o = n.ma.masked_where(c_o < 1, c_o)
        cnt_e[bl] = n.transpose(c_e, [1, 0])[:, ::dlst]
        var_e[bl] = n.transpose(v_e, [1, 0])[:, ::dlst]  # * (inttime*sdf*1e9)
        cnt_o[bl] = n.transpose(c_o, [1, 0])[:, ::dlst]
        var_o[bl] = n.transpose(v_o, [1, 0])[:, ::dlst]  # * (inttime*sdf*1e9)
        D = np.ma.diff(data['even'][bl], axis=1)
        trms_df[bl] = np.sqrt(np.ma.mean(D.conj() * D/2, axis=0)) * jy2T[:-1]
bls_master = x.keys()
nbls = len(bls_master)
bls_master.sort()

trms_data, trms_o_the, trms_e_the = {}, {}, {}
theo_temp = {}
#compute the Trms three ways.  differenced data, variance data, model
for bl in bls_master:
    trms_data[bl] = n.sqrt(n.mean(x[bl].conj() * x[bl], axis=1)/2.).real * jy2T
    trms_e_the[bl] = n.sqrt( n.ma.mean(var_e[bl]/cnt_e[bl]/(inttime/opts.varinttime), axis=1)) * jy2T
    trms_o_the[bl] = n.sqrt( n.ma.mean(var_o[bl]/cnt_o[bl]/(inttime/opts.varinttime), axis=1)) * jy2T
##GSM emission + antenna temp / sqrt(df*dt*cnt)
    theo_temp[bl] = (opts.Tsky*1e3*(allfreqs/.18)**(-2.8)+opts.Trcvr*1e3)/(n.sqrt(inttime * sdf * 1e9)) * n.mean(1./n.sqrt(cnt_e[bl]), axis=1)
# bl_avgeraged counts and vars
var_o_blavg = n.mean([var_o[bl] for bl in bls_master], axis=0)
cnt_o_blavg = n.mean([cnt_o[bl] for bl in bls_master], axis=0)

var_e_blavg = n.mean([var_e[bl] for bl in bls_master], axis=0)
cnt_e_blavg = n.mean([cnt_e[bl] for bl in bls_master], axis=0)
print "nbls = ", nbls
# Creating models from bl averaged data
trms_e_the_blavg = n.sqrt(n.mean(var_e_blavg/(cnt_e_blavg*nbls), axis=1))*jy2T
trms_o_the_blavg = n.sqrt(n.mean(var_o_blavg/(cnt_o_blavg*nbls), axis=1))*jy2T

# the temperature of the bls averaged together
diff_blavg = n.mean([x[bl] for bl in bls_master], axis=0)
trms_blavg = n.sqrt(n.mean(diff_blavg.conj() * diff_blavg, axis=1)/2.) * jy2T

avg_trms = n.mean([trms_data[bl] for bl in bls_master], axis=0)
gsm_data = {}

# for bl in bls_master:
#    gsm_data[bl] = n.ma.masked_invalid(theo_temp[bl]/avg_trms)
# print 'baseline:', bl, 'mean, median, std ratio', n.mean(gsm_data[bl]), n.median(gsm_data[bl]), n.std(gsm_data[bl])

# blavg_ratio=n.ma.masked_invalid(trms_blavg/trms_o_the_blavg)
# print 'BL Averaged Trms/Tvar:','meam', n.mean(blavg_ratio),'meadian', n.median(blavg_ratio)

T_vs_bl = n.zeros((len(bls_master), len(band_chans)))
for i, bl in enumerate(bls_master):
    T_vs_bl[i] = trms_data[bl]
p.subplot(121)
p.plot(bls_master, n.mean(T_vs_bl[:, chans], axis=1))
p.subplot(122)
p.imshow(n.log(T_vs_bl), aspect='auto', interpolation='nearest')
T_vs_bl_mean = n.mean(T_vs_bl[:, chans], axis=1)
hot_bl = bls_master[n.argwhere(T_vs_bl_mean == T_vs_bl_mean.max()).squeeze()]
print "hottest baseline = ", hot_bl, a.miriad.ij2bl(*hot_bl)
p.show(block=True)
p.clf()

for i, bl in enumerate(bls_master):
    print 'Ploting Trms for %d_%d' % bl
    p.subplot(121)
    p.imshow(np.ma.abs(x[bl]).T, aspect='auto',
             vmax=2*np.ma.median(trms_data[bl]))
    p.colorbar()
    p.subplot(122)
    p.plot(band_chans, trms_data[bl], label='$T_{e-o}$')
    p.plot(band_chans, avg_trms, 'm--', label='$T_{e-o,avg}$')
    p.plot(band_chans, trms_o_the[bl], label='$T_{odd,var}$')
    p.plot(band_chans, trms_e_the[bl], label='$T_{even,var}$')
    p.plot(band_chans, trms_blavg, label='$T_{blavg,e-o}$')
    # p.plot(band_chans, trms_o_the_blavg, label='$T_{blavg,odd,var}$')
    # p.plot(band_chans, trms_e_the_blavg, label='$T_{blavg,even,var}$')
    p.plot(band_chans, theo_temp[bl], 'k--', label='GSM')
    p.plot(band_chans, theo_temp[bl]/n.sqrt(nbls), 'k--')
    p.plot(trms_df[bl], label='$Trms_{df}$')
    p.yscale('log')
    p.xlim([band_chans[0]-1, band_chans[-1]+1])
    p.xlabel('chan')
    p.ylabel('$T_{rms}\ [mK]$')
    p.legend(loc='lower center', ncol=4, bbox_to_anchor=(-.15, -.35))
    p.subplots_adjust(bottom=.25)
    p.title('%d\_%d' % bl)
    if opts.plot:
        p.show(block=True)

    outfile = 'Trms_%d_%d.png' % bl
    if not opts.output == '':
        outfile = opts.output+'/'+outfile
        p.savefig(outfile, format='png', bbox_inches='tight')
    p.clf()
