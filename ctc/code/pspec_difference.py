#! /usr/bin/env python

import numpy as n
import optparse
import sys
import os

o = optparse.OptionParser()
opts, args = o.parse_args(sys.argv[1:])

# Check to make sure 2 files were supplied
if len(args) != 2: 
    print "Must supply two files to difference."
    sys.exit()

# Load files
print "Loading", args[0]
f1 = n.load(args[0])
print"Loading", args[1]
f2 = n.load(args[1])
os.system('mkdir '+args[0].split('/')[-2]) # make inject directory
fname = args[0].split('/')[-2]+'/'+args[0].split('/')[-1] # name of outputted file
fnew = {}

# Difference keys
for key in f1:
    if key[0] == 'p': # PS channel
        fnew[key] = f1[key] - f2[key]        
    else:
        fnew[key] = f1[key] # copy other channels
# Save new file
n.savez(fname, **fnew)
print "   Saving", fname
