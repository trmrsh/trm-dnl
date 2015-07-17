#!/usr/bin/env python

"""
Script to extract particular set of molly spectra given a database file
created by fmolly.py. It stores the name of the original file in the 
molly header.

trm.dnl.molly, trm.simbad, trm.sla and trm.subs
"""

from optparse import OptionParser

usage = """usage: %prog [options] dbase mfile target

Creates a molly file consisting of all the molly spectra that match a given target
or position 

dbase   -- name of FITS file database created by fmolly.py
mfile   -- name of molly file (.mol will be added if need be)
target  -- target name or position. e.g. it could be 'IP Peg' or '01 23 34.5 -34 33 05.2'.
           It will assume that you have specified a target name if any non-numeric character is found. 
           ICRS is assumed.
"""

parser = OptionParser(usage)
parser.add_option("-a", "--append", dest="append", default=False, action="store_true",\
                      help="append spectra to any existing molly file")
parser.add_option("-c", "--clobber", dest="clobber", default=False, action="store_true",\
                      help="clobber files on output")
parser.add_option("-d", "--dist", dest="dist", default=4., type="float",\
                      help="distance in arcminutes for spectrum selection; default = 4")
parser.add_option("-s", "--start", dest="start", default='1900-01-01', \
                      help="earliest date to include spectra YYYY-MM-DD; default = 1900-01-01")
parser.add_option("-e", "--end", dest="end", default='2100-01-01', \
                      help="latest date to include spectra YYYY-MM-DD; default = 2100-01-01")
parser.add_option("-w", "--wave", dest="wave", default=-100., type='float', \
                      help="wavelength that must be covered (angstroms), negative to ignore; default = -100.")
(options, args) = parser.parse_args()

if len(args) < 3:
    print 'Require at least three arguments. Use -h for help'
    exit(1)

import os

dbase = args.pop(0)
if not os.path.exists(dbase):
    print 'Cannot find ',dbase
    exit(1)

mfile = args.pop(0)
if not mfile.endswith('.mol'):
    mfile += '.mol'

clobber = options.clobber
append  = options.append
if not clobber and not append and os.path.exists(mfile):
    print 'File called',mfile,'already exists and will not be overwritten'
    print 'Use the -c option to clobber on output'
    exit(1)
elif not os.path.exists(mfile):
    append = False

# Now try to interpret target. The ultimate aim is 
# to arrive at a position to search around.

target = ' '.join(args)

import re
import trm.subs as subs

if re.search('[a-zA-Z]', target):
    import trm.simbad as simbad
    results = simbad.Query(target).query()
    if len(results):
        print 'Simbad returned coordinates: ',results[0]['Position']
        (ra,dec,system) = subs.str2radec(results[0]['Position'])
else:
    (ra,dec,system) = subs.str2radec(target)

dist  = options.dist
wave  = options.wave

import trm.sla as sla

dre = re.compile('(\d\d\d\d)-(\d\d)-(\d\d)$')
 
m = re.match(dre, options.start) 
if m:
    year  = int(m.group(1))
    month = int(m.group(2))
    day   = int(m.group(3))
    start_mjd = sla.cldj(year, month, day)
else:
    print 'Could not understand the start date:',options.start
    exit(1)

m = re.match(dre, options.end) 
if m:
    year  = int(m.group(1))
    month = int(m.group(2))
    day   = int(m.group(3))
    end_mjd = sla.cldj(year, month, day)
else:
    print 'Could not understand the end date:',options.end
    exit(1)

if start_mjd >= end_mjd:
    print 'Start date must be before end date.'
    exit(1)
   
# OK now read the database file    
import pyfits

hdul   = pyfits.open(dbase)
table  = hdul[1].data

# Get positions, kick out unknown ones
import math

rar     = table.field('RA')
decr    = table.field('Dec')
mjd     = table.field('MJD')
ok      = (rar > -0.01) & (mjd > start_mjd) & (mjd < end_mjd)
if wave > 0.:
    wmin    = table.field('Min. wave')
    wmax    = table.field('Max. wave')
    ok = ok & (wmin < wave) & (wave < wmax)

rar     = 15.*math.pi*rar[ok]/180.
decr    = math.pi*decr[ok]/180.
targets = table.field('Target')[ok]
mfiles  = table.field('File')[ok]
slots   = table.field('Slot')[ok]
hdul.close()

# Convert positions to vectors
import numpy
cosd   = numpy.cos(decr)
v      = numpy.column_stack((cosd*numpy.cos(rar),cosd*numpy.sin(rar),numpy.sin(decr)))

rac    = math.radians(15.*ra)
decc   = math.radians(dec)
cosd   = numpy.cos(decc)
tv     = (cosd*numpy.cos(rac),cosd*numpy.sin(rac),numpy.sin(decc))

# create array of dot products
dot    = tv[0]*v[:,0]+tv[1]*v[:,1]+tv[2]*v[:,2]

# dot product limit
dlim   = math.cos(math.radians(dist/60.))

select = dot > dlim

mfiles = mfiles[select]
slots  = slots[select]

# make sure we read the spectra as efficiently as possible by
# opening each file once only and reading/skipping files as 
# necessary
unique = set(mfiles)

import trm.dnl.molly as molly

if append:
    mout = open(mfile, 'ab')
else:
    mout = open(mfile, 'wb')

nwrote = 0
for mfin in unique:
    slts = slots[mfiles == mfin]
    smax = slts.max()
    try:
        minp  = open(mfin,'rb')
        nspec = 1
        while nspec <= smax:
            if nspec in slts:
                mspec = molly.rmspec(minp)
                mspec['Filename'] = mfin[-36:-4]
                mspec.wmolly(mout, hasarc=True)
            else:
                if not molly.skip_molly(minp):
                    print 'Unexpected error while trying to skip spectrum',nspec,'of file,',mfin
                    break
            nspec += 1
        minp.close()
        print 'Read',len(slts),'spectra from',mfin
        nwrote += len(slts)
    except:
        print 'Error reading file',mfin
mout.close()

print 'Read/wrote',nwrote,'spectra to',mfile
print 'Failed to read',len(slots)-nwrote,'spectra.'
            
            
    

