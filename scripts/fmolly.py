#!/usr/bin/env python

"""
Script to search out all molly files (*.mol) in a directory tree and build a
database in the form of a FITS file with a FITS binary table from them. The
FITS table is most simply viewed with 'fv' but otherwise 'pyfits' in Python
makes it easy to perform searches on it. All coordinates are precessed to
J2000 if possible, or else set to obvious bad values (-1, -100). The script 
also attempts to group spectra by position on the sky by searching recursively 
for all spectra within a certain distance in arcminutes of each other. This is 
to cope with variations in positions. The group classification (an integer 
starting from 0, with -1 as a bad value) as well as the distance in arcminutes 
of each position from the mean of its group are listed in the FITS table that 
is generated.

This script can take a while to run. It is quite robust though, so try just to 
leave it to do its job. e.g. it will read files until it encounters an error but
will then just move onto the next. There are five stages. First, it looks for all 
molly files starting from the top-level directories specified. Second, it reads every 
spectrum in every file. Third, it defines groups of spectra in which every spectrum 
lies at most 'dist' arcminutes from the other spectra in its group. The mean RA and Dec
and the names of the group are printed to the screen. You may want to re-direct these 
results to a file. Fourth, it orders the spectra by group number and then my MJD within 
each group. finally it writes out the data to a FITS file.

It takes about 10 minutes to run this on all of my files and directories, reading 127000
spectra in 2000 files, so it can be updated nightly.

Usage: fmolly.py name dir1 dir2 ..

name   -- name of output fits database file (must end in '.fit' or '.fits'). 
          The routine will not overwrite existing files. 

dist   -- distance in arcminutes to define spectrum groups

dir1   -- first top-level directory to search. Best to supply an absolute
          path

dir2   -- second top-level directory

Any number of directories can be supplied.

Python modules needed: math, datetime, pyfits, numpy,
trm.dnl.molly, trm.sla and trm.subs
"""

from optparse import OptionParser

usage = """usage: %prog [options] dbase dir1 [dir2 dir3 ..]

Generates a FITS database of molly spectra by searching through directory
trees for all molly files, reading them and storing basic data on them in
a FITS file.

dbase   -- name of FITS file that will be created (including .fit or .fits)
dir1    -- top-level of first directory tree to be searched (required)
dir2    -- top-level of second directory tree to be searched (optional);
           any number of such directories can be specified.
"""

parser = OptionParser(usage)
parser.add_option("-a", "--arc", dest="noarc", default=False, action="store_true",\
                  help="ignore arcn.mol and arco.mol files")
parser.add_option("-s", "--sky", dest="nosky", default=False, action="store_true",\
                  help="ignore skyn.mol and skyo.mol files")
parser.add_option("-n", "--nor", dest="nonor", default=False, action="store_true",\
                  help="ignore nor.mol files")
parser.add_option("-d", "--dist", dest="dist", default=4., type="float",\
                  help="distance in arcminutes to define spectrum groups; default = 4")
parser.add_option("-f", "--full", dest="full", default=True, action='store_true',\
                      help="full directory search including dot files")
(options, args) = parser.parse_args()

if len(args) < 2:
    print 'Require at least two arguments. Use -h for help'
    exit(1)

import os

fname = args[0]
if os.path.exists(fname):
    print 'File =',fname,'already exists.'
    exit(1)

if fname.rfind('.fit')  == -1 or (fname.rfind('.fit') != len(fname)-4 and \
                                      fname.rfind('.fits') != len(fname)-5):
    print 'Filename must end in .fit or .fits'
    exit(1)

dist = options.dist

# check list of directories
dirs = args[1:]
for dir in dirs:
    if not os.path.isdir(dir):
        print dir,'is not a directory.'
        exit(1)

import math
import datetime
import pyfits
import numpy
import trm.subs as subs
import trm.dnl.molly as molly
import trm.sla as sla

print 'Looking for molly files ...'

full = options.full

# find all molly files recursively in all directories.
mflist = []    
for dir in dirs:
    for (root,d,files) in os.walk(dir):
        if full or (not root.startswith('.') and not root.find('/.')):
            mflist += [os.path.join(root, f) for f in files if f.endswith('.mol')]  

# apply filters
if options.noarc:
    mflist = [f for f in mflist if not f.endswith('arco.mol') and not f.endswith('arcn.mol')]

if options.nosky:
    mflist = [f for f in mflist if not f.endswith('skyo.mol') and not f.endswith('skyn.mol')]

if options.nonor:
    mflist = [f for f in mflist if not f.endswith('nor.mol')]

if len(mflist) == 0:
    print 'No molly files found. Database composition aborted.'
    exit(1)

print 'Found',len(mflist),'molly files, now reading them ...'

# OK, ready to go

objects = []
ras     = []
decs    = []
mjds    = []
dwells  = []
nrecs   = []
npixs   = []
wmins   = []
wmaxs   = []
nspecs  = []
mfiles  = []

nnoread   = 0
nnopos    = 0
nnowave   = 0
tspec     = 0
nbadeq    = 0
nnojd     = 0
nnoexp    = 0
nnoobject = 0
ntcalc    = 0
nnorec    = 0

for mfile in mflist:
    try:
        for (nspec,mspec) in molly.grmolly(mfile):

            tspec += 1

            # Set all items first. A value is always set although
            # some will be bad as in -1 for RA

            # object name
            if 'Object' in mspec:
                object = mspec['Object']
            else:
                nnoobject += 1
                object = 'UNKNOWN'

            # position
            if 'RA' in mspec and 'Dec' in mspec and 'Equinox' in mspec:
                ra  = mspec['RA']
                dec = mspec['Dec']
                eq  = mspec['Equinox']
                if ra < 0 or ra >= 24. or dec < -90. or dec > 90.:
                    print 'ERROR: RA and/or Dec out of range =' + str(ra) + ' ' + str(dec) + ' file = ' + mfile + ' nspec = ' + str(nspec)
                    ra  = -1.
                    dec = -100.
                    nnopos += 1
                elif abs(eq - 1950) < 0.001:
                    (ra,dec,d1,d2,d3,d4) = sla.fk425(ra,dec)
                elif abs(eq - 2000) > 0.001:
                    ra  = -1.
                    dec = -100.
                    print 'ERROR: unrecognised equinox ' + str(eq) + ' file = ' + mfile + ' nspec = ' + str(nspec)
                    nnopos  += 1
                    nbadeq  += 1
            else:
                ra   = -1.
                dec  = -100.
                nnopos += 1

            # Time
            if 'RJD' in mspec:
                mjd = mspec['RJD'] -2400000.5
            elif 'Day' in mspec and 'Month' in mspec and 'Year' in mspec and 'UTC' in mspec:
                mjd = sla.cldj(mspec['Year'],mspec['Month'],mspec['Day']) + mspec['UTC']/24.
                ntcalc += 1
            else:
                mjd = 0.
                nnojd += 1

            # Dwell
            if 'Dwell' in mspec:
                dwell = mspec['Dwell']
            else:
                dwell = -1.
                nnoexp += 1

            # Record number
            if 'Record' in mspec:
                nrec = mspec['Record']
            else:
                nrec = -1
                nnorec += 1

            # number of pixels
            npix = len(mspec)

            # Min and max wavelengths
            if mspec.narc != 0:
                wmin = mspec.x.dat.min()
                wmax = mspec.x.dat.max()
            else:
                wmin = -1.
                wmax = -1.
                nnowave += 1

            # Finally store
            objects.append(object)
            ras.append(ra)
            decs.append(dec)
            mjds.append(mjd)
            dwells.append(dwell)
            nrecs.append(nrec)
            npixs.append(npix)
            wmins.append(wmin)
            wmaxs.append(wmax)
            nspecs.append(nspec)
            mfiles.append(mfile)

    except:
        print 'ERROR: Failed to read file = ' + mfile 
        nnoread += 1

print '\nA total of',tspec,'spectra were read.'
print 'There were',nnoread,'files which had read errors.'
print 'There were',nnopos,'spectra without coordinates.'
print 'There were',nnowave,'spectra without a wavelength scale.'
print 'There were',nbadeq,'spectra with bad equinoxes.'
print 'There were',nnojd,'spectra without times.'
print 'There were',nnoexp,'spectra without exposure times.'
print 'There were',nnoobject,'spectra without object names.'
print 'There were',nnoobject,'spectra without record numbers.'
    

objects = numpy.array(objects)
ras     = numpy.array(ras,    dtype='float')
decs    = numpy.array(decs,   dtype='float')
mjds    = numpy.array(mjds,   dtype='float')
dwells  = numpy.array(dwells, dtype='float32')
nrecs   = numpy.array(nrecs,  dtype='int32')
npixs   = numpy.array(npixs,  dtype='int32')
wmins   = numpy.array(wmins,  dtype='float32')
wmaxs   = numpy.array(wmins,  dtype='float32')
nspecs  = numpy.array(nspecs, dtype='int32')
mfiles  = numpy.array(mfiles)

print '\nWill now group spectra by position ...\n'

def dot(v1, v2):
    """
    Returns dot product(s) between 3D vectors allowing
    one or other of the arguments to be an array
    """
    if len(v1.shape) == 1 and len(v2.shape) == 1:
        return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
    elif len(v1.shape) == 1 and len(v2.shape) == 2:
        return v1[0]*v2[:,0]+v1[1]*v2[:,1]+v1[2]*v2[:,2]
    elif len(v1.shape) == 2 and len(v2.shape) == 1:
        return v1[:,0]*v2[0]+v1[:,1]*v2[1]+v1[:,2]*v2[2]
    else:
        print 'dot: undefined operation. arguments:',v1,v2
        exit(1)

def find_group(vecc, vec, indx, gindx, thresh):
    """
    Recursive group searcher used below. Given a 3-vector test position 
    and an array of 3-vector positions it searches for all that are closer
    than the threshold defined by thresh. This search is repeatedly
    applied to all fiound in the group so that the group will spread
    to include others not found in the first pass.

    vecc   - test vector (3 element numpy array)
    vec    - vectors to test against
    indx   - indices corresponding to vec
    gindx  - group of indices identified 
    thresh - threshold for dot product

    gindx is returned with updated value of the indices identifying the group members
    vec and indx are returned with any items identified in gindx removed.
    """

    dt     = dot(vecc, vec)
    sgroup = dt > thresh
    sindx  = indx[sgroup]
    if len(sindx) == 0:
        return (gindx,vec,indx)

    # update indices of the group
    gindx = numpy.append(gindx, sindx)

    vg   = vec[sgroup]
    vec  = vec[~sgroup]
    indx = indx[~sgroup]

    for v in vg:
        (gindx,vec,indx) = find_group(v, vec, indx, gindx, thresh)

    return (gindx,vec,indx)


# convert to radians
rar    = 15.*math.pi*ras/180.
decr   = math.pi*decs/180.

# kick out bad positions, set up group number and distance
# arrays
groups = numpy.zeros((len(rar))) - 1
gdist  = numpy.zeros((len(rar)),dtype='float') - 1.
ok     = rar > -0.01
indxok = numpy.arange(len(ok),dtype='int')[ok]
rar    = rar[ok]
decr   = decr[ok]

# Convert positions to vectors
cosd   = numpy.cos(decr)
vecs   = numpy.column_stack((cosd*numpy.cos(rar),cosd*numpy.sin(rar),numpy.sin(decr)))
nvecs  = len(vecs)
cvec   = vecs[:]

# create a mask array set to False, apart from the first 
# this is used to record which spectra have already been
# classified (groups are mutually exclusive).

mask   = cosd > 100
indx   = numpy.arange(len(mask),dtype='int')

thresh = math.cos(math.radians(dist/60.))
total  = 0
ngroup = 0

radict = {}
ginds = []
for i in xrange(nvecs):
    if not mask[i]:
        v     = vecs[0]
        gindx = numpy.array((i,))
        vecs  = vecs[1:]
        indx  = indx[1:]
        (gindx,vecs,indx) = find_group(v, vecs, indx, gindx, thresh)

        # save the indices
        ginds.append(indxok[gindx])

        # store which have been grouped
        mask[gindx] = True
        total += len(gindx)

        # compute mean position of group
        vmean = cvec[gindx,:].mean(0)

        # calculate distance of all members in group from the mean
        gdist[ginds[-1]] = 60.*180.*numpy.arccos(numpy.minimum(1.,dot(vmean, cvec[gindx,:])))/math.pi

        # store the mean RA, Dec of the group, keyed on RA
        mdec = math.degrees(math.asin(vmean[2]))
        mra  = math.degrees(math.atan2(vmean[1], vmean[0]))/15.
        if mra < 0: mra += 24.
        radict[mra] = (ngroup, mra, mdec)
        ngroup += 1

print 'Ordering spectra into groups by RA and by MJD within groups ...\n'

# Next parts are to re-label the groups according to RA order
rakey = radict.keys()
rakey.sort()

# Now relabel each group
for g,key in enumerate(rakey):
    ninds = ginds[radict[key][0]]
    groups[ninds] = g
    (mra,mdec) = radict[key][1:]
    # The next rather horrific line derives a unique set of names from all the names that the object
    # has, and does so case-insensitively to reduce numbers. The conversion back to normal 'str' strings 
    # is to avoid iStr being printed around each name.
    nms = [str(x) for x in set([subs.iStr(x) for x in objects[ninds]])]
    print 'Group %6i, RA,Dec = %s %s, nspec = %5i, dmax = %6.3f, names = ' % (g,subs.d2hms(mra,dp=2),subs.d2hms(mdec,dp=1,sign=True), \
                                                                                  len(ninds),gdist[ninds].max()) + str(nms)

# Sort into group order
order = groups.argsort()

objects = objects[order]
groups  = groups[order]
gdist   = gdist[order]
ras     = ras[order]
decs    = decs[order]
mjds    = mjds[order]
dwells  = dwells[order]
nrecs   = nrecs[order]
npixs   = npixs[order]
wmins   = wmins[order]
wmaxs   = wmaxs[order]
nspecs  = nspecs[order]
mfiles  = mfiles[order]

# Within groups, sort by MJD
for ng in xrange(-1,ngroup):
    gi = groups == ng
    
    order = mjds[gi].argsort()
    if len(order):
        objects[gi] = objects[gi][order]
        groups[gi]  = groups[gi][order]
        gdist[gi]   = gdist[gi][order]
        ras[gi]     = ras[gi][order]
        decs[gi]    = decs[gi][order]
        mjds[gi]    = mjds[gi][order]
        dwells[gi]  = dwells[gi][order]
        nrecs[gi]   = nrecs[gi][order]
        npixs[gi]   = npixs[gi][order]
        wmins[gi]   = wmins[gi][order]
        wmaxs[gi]   = wmaxs[gi][order]
        nspecs[gi]  = nspecs[gi][order]
        mfiles[gi]  = mfiles[gi][order]

print '\nWriting out fits file ...\n'

# Define the columns
cols = []
cols.append(pyfits.Column(name='Target', format='A'+str(objects.itemsize)))
cols.append(pyfits.Column(name='Group', format='J'))
cols.append(pyfits.Column(name='Distance', unit='arcmin', format='E'))
cols.append(pyfits.Column(name='RA',  unit='hours', format='D'))
cols.append(pyfits.Column(name='Dec', unit='degrees', format='D'))
cols.append(pyfits.Column(name='Record', format='J'))
cols.append(pyfits.Column(name='MJD', format='D'))
cols.append(pyfits.Column(name='Dwell', unit='sec', format='E'))
cols.append(pyfits.Column(name='Pixels', format='J'))
cols.append(pyfits.Column(name='Min. wave', unit='A', format='E'))
cols.append(pyfits.Column(name='Max. wave', unit='A', format='E'))
cols.append(pyfits.Column(name='Slot', format='J'))
cols.append(pyfits.Column(name='File', format='A'+str(mfiles.itemsize)))

header = pyfits.Header()
header.update('TIME', str(datetime.datetime.today()), 'when database was compiled')
header.update('DIST', dist, 'distance in arcmin used for defining groups')
hdus = [pyfits.PrimaryHDU(header=header),]

tbhdu = pyfits.new_table(cols, nrows=len(objects))
tbhdu.data.field('Target')[:]      = objects
tbhdu.data.field('Group')[:]       = groups
tbhdu.data.field('Distance')[:]    = gdist
tbhdu.data.field('RA')[:]          = ras
tbhdu.data.field('Dec')[:]         = decs
tbhdu.data.field('MJD')[:]         = mjds
tbhdu.data.field('Dwell')[:]       = dwells
tbhdu.data.field('Record')[:]      = nrecs
tbhdu.data.field('Pixels')[:]      = npixs
tbhdu.data.field('Min. wave')[:]   = wmins
tbhdu.data.field('Max. wave')[:]   = wmaxs
tbhdu.data.field('Slot')[:]        = nspecs
tbhdu.data.field('File')[:]        = mfiles
hdus.append(tbhdu)
hdulist = pyfits.HDUList(hdus)
hdulist.writeto(fname)

print 'All finished! Use "fv" to look at the database file =',fname
