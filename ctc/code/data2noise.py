#! /usr/bin/env python

import os
import aipy
import numpy
import capo
import sys
import optparse

"""
Multiplies data by alternating -1 and 1's along the time axis (random in frequency).
Consequently, if data is FRF, foregrounds will be filtered out since they will have very high fringe rates, and we'll be left with a thermal noise approximation.
"""


### Options ###
o = optparse.OptionParser()
o.set_usage('data2noise.py *uvGA')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])


def mfunc(uv,p,d,f):
    d *= factors[p[2]][p[1]]
    return p,d,f

for file in args:
    # read information from file
    uvi = aipy.miriad.UV(file)
    pol = file.split('.')[-2]
    times,data,_ = capo.miriad.read_files([file],antstr='all',polstr=pol)
    jds = times['times']
    bls = data.keys()
    # build mask
    factors = {}
    for bl in bls:
        factors[bl] = {}
        mask = numpy.random.choice([-1,1],size=data[bl][pol].shape[1]) # random 1 or -1 along frequency
        for jd in jds:
            factors[bl][jd] = mask
            mask = mask*-1
    # create new file
    newfile = file+'n'
    if os.path.exists(newfile): 
        print '   %s exists. Skipping...' % newfile
        continue
    print file, '->', newfile
    uvo = aipy.miriad.UV(newfile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, raw=True, mfunc=mfunc, append2hist='DATA2NOISE: ' + ' '.join(sys.argv) + '\n')
