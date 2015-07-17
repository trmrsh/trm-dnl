#!/usr/bin/env python

from __future__ import print_function

usage = \
"""
ulog2lc generates a light curve from an ULTRACAM log file applying
time corrections and (optionally) normally by a comparison or the sum
of comparisons. If no comparisons are supplied, just the raw data are
returned. For ULTRASPEC set ccd=1.
"""

import argparse, os
import numpy as np
from astropy.io import fits
from trm import subs
from trm import dnl

parser = argparse.ArgumentParser(description=usage)

# positional
parser.add_argument('ulog', help='name of input ULTRACAM log file')
parser.add_argument('output', help='name of output FITS file (.fits will be added if not present)')
parser.add_argument('name', help='target name for header')
parser.add_argument('tel', help='telescope name')
parser.add_argument('ccd', type=int, help='CCD number (1=red, 2=green, 3=blue)')
parser.add_argument('targ', type=int, help='target aperture number')
parser.add_argument('ra', help='RA in hh:mm:ss.sss form')
parser.add_argument('dec', help='Dec in hh:mm:ss.sss form')
parser.add_argument('comp', type=int, nargs='*',
                    help='comparison aperture number(s) (>1 will be added)')

# optional
parser.add_argument('-a', dest='append', action='store_true',
                    help='Will append to any existing output file rather than clobber it')

parser.add_argument('-c', dest='clobber', action='store_true',
                    help='Clobber existing output files')

parser.add_argument('-f', dest='filter', help='Filter name for header')

# OK, done with arguments.
args = parser.parse_args()

if args.targ < 0:
    print('Target aperture must be > 0')
    exit(1)

for naper in args.comp:
    if naper < 0:
        print('Comparison aperture must be > 0')
        exit(1)

TELNAMES = ['WHT', 'TNT', 'NTT', 'VLT']
if args.tel not in TELNAMES:
    print('Telescope name must be one of: ' + str(TELNAMES))
    exit(1)

if args.ccd < 1 or args.ccd > 3:
    print('CCD number must range from 1 to 3')
    exit(1)

# read log file
ulog = dnl.ulog.Ulog(args.ulog)

# extract target Dset
targ = ulog.tseries(args.ccd, args.targ)

# extract and divide by comparison(s) if specified
if len(args.comp) > 0:
    comp = ulog.tseries(args.ccd, args.comp[0])
    for naper in args.comp[1:]:
        comp += ulog.tseries(args.ccd, naper)
    targ /= comp

# add run name
targ.head['ULOG'] = (args.ulog, 'Name of ULTRACAM/ULTRASPEC log file')

title = args.name + ', ' + args.ulog

# add filter (if specified)
if args.filter is not None:
    targ.head['FILTER'] = (args.filter, 'Filter used')
    title += ', ' + args.filter

# add title for plots
targ.head['TITLE'] = (title, 'Title for plots')

# add positional data
targ.setpos(args.name, args.ra, args.dec)

# add observatory data
telescope,observatory,longitude,latitude,height = subs.observatory(args.tel)
targ.settel(telescope, observatory, longitude, latitude, height)

# convert times
targ.utc2tdb()

# save
head = fits.Header()
head['comment'] = 'File generated from ULTRACAM log file by ulog2lc.py'

fout = args.output if args.output.endswith('.fits') else args.output + '.fits'
dnl.io.wfits(fout, targ, head, args.append, args.clobber)

