#!/usr/bin/env python

"""
Script to convert a molly file into a FITS file

Usage: mol2fits.py name

name   -- root name of molly file which will also be used for the output fits file

Python modules needed: math, datetime, pyfits, numpy,
trm.dnl.molly, trm.sla and trm.subs
"""

from optparse import OptionParser

usage = """usage: %prog [options] name

Converts a molly file to an equivalent FITS file. The FITS file created
has one table HDU per molly spectrum with columns of wavelength, flux and error 
on the flux.


"""

parser = OptionParser(usage)
parser.add_option("-f", "--flam", dest="flam", default=False, action="store_true",\
                  help="write out as flambda (as opposed to fnu)")

#parser.add_option("-s", "--sky", dest="nosky", default=False, action="store_true",\
#                  help="ignore skyn.mol and skyo.mol files")
#parser.add_option("-n", "--nor", dest="nonor", default=False, action="store_true",\
#                  help="ignore nor.mol files")
#parser.add_option("-d", "--dist", dest="dist", default=4., type="float",\
#                  help="distance in arcminutes to define spectrum groups; default = 4")
#parser.add_option("-f", "--full", dest="full", default=True, action='store_true',\
#                      help="full directory search including dot files")
(options, args) = parser.parse_args()

if len(args) != 1:
    print 'Require an argument. Use -h for help'
    exit(1)

import os, time, math, pyfits, numpy
import trm.subs as subs
import trm.dnl.molly as molly

root  = args[0]
if root.endswith('.mol'):
    root = root[:-4]

mname = root + '.mol'
if not os.path.exists(mname):
    print 'Cannot fine file =',mname
    exit(1)

fname = root + '.fits'
if os.path.exists(fname):
    print 'File =',fname,'already exists.'
    exit(1)

# Read in spectra
specs = []
for nspec, mspec in molly.grmolly(mname):
    specs.append(mspec)

def update_header(mspec, thead):
    """
    Tries to convert as much of a standard molly header as it can into more standard
    FITS (note that molly headers are stored in FITS headers in trm.dnl but mostly
    by using HIERARCH
    """

    mhead = mspec.head

    # Modify the molly header for a few special cases to allow the more
    # generic changes to be nicely ordered,
    if 'Day' in mhead and 'Month' in mhead and 'Year' in mhead:
        mhead['DATE-OBS'] = '{0:4d}-{1:02d}-{2:02d}'.format(mhead['Year'],mhead['Month'],mhead['Day'])
    if 'Equinox' in mhead:
        if mhead['Equinox'] == 1950.:
            mhead['Equinox'] = 'B1950'
        elif mhead['Equinox'] == 2000.:
            mhead['Equinox'] = 'J2000'
    mhead['NARC'] = mspec.narc

    # list of translations (old,new,comment)
    trans = [ \
        ('Object','OBJECT','target name'),
        ('DATE-OBS', 'DATE-OBS', 'Date, YYYY-MM-DD, at centre of exposure'),
        ('UTC','UTC','UTC at centre of exposure (decimal hours)'),
        ('RJD','JD','JD (UTC) at centre of exposure'),
        ('Dwell','EXPOSURE','exposure time (seconds)'),
        ('HJD','HJD','Heliocentric JD (UTC) at centre of exposure'),
        ('Delay','DELAY','JD-HJD (UTC), days'),
        ('Sidereal time', 'SIDEREAL', 'local sidereal time (hours)'),
        ('Hour angle', 'HOURANGL', 'hour angle, hours past meridian'),
        ('Airmass', 'AIRMASS', 'airmass at centre of exposure'),
        ('Vearth','VEARTH','apparent target speed due to Earth (km/s)'),
        ('RA','RA','Right Ascension, hours'),
        ('Dec','DEC','Declination, degrees'),
        ('Equinox','EQUINOX','Equinox of RA, Dec'),
        ('Gal longitude', 'GLONGITU', 'galactic longitude, degrees'),
        ('Gal latitude', 'GLATITUD', 'galactic latitude, degrees'),
        ('Record','RECORD', 'sequential record number'),
        ('Night','NIGHT', 'night number of run'),
        ('Arc(s) used', 'ARCS', 'record number(s) of calibration arc(s)'),
        ('NARC','NARC','number of coefficients for wavelength poly'),
        ('Telescope','TELESCOP','telescope'),
        ('Site','SITE','observing site'),
        ('Longitude','LONGITUD','site longitude, degrees west'),
        ('Latitude','LATITUDE','site latitude, degrees'),
        ('Height','HEIGHT','site height above sea level, metres'),
        ('Extract position','EXTRACT','pixel position along slit of extraction')
        ]

    for old, new, comment in trans:
        if old in mhead: thead[new] = (mhead[old], comment)

    # now copy over anything missed
    skip = [old for old,new,comment in trans] + ['Day','Month','Year']
    for key in mhead:
        if key not in skip and not key.lower().startswith('comment'):
            if len(key) > 8 or key.find(' ') > -1:
                thead['HIERARCH ' + key] = mhead[key]
            else:
                thead[key] = mhead[key]
        elif key == 'COMMENT':
            for val in mhead['COMMENT']:
                thead.add_comment(val)


# Write them out

phdu = pyfits.PrimaryHDU()
head = phdu.header
fpath = os.path.abspath(mname)
start = 0 if len(fpath) < 79 else len(fpath)-78 
head['SOURCE']  = fpath[start:]
head['CREATED'] = (time.asctime(), 'creation date')
head.add_comment(' ')
head.add_comment('This file contains spectra translated from molly by mol2fits.py')
head.add_comment(' ')
head.add_comment('There is one HDU per spectrum each with its own header. The')
head.add_comment('spectra are stored as binary tables with either wavelength,')
head.add_comment('flux, error-in-the-flux and flux/count columns or wavelength,')
head.add_comment('counts and error-in-counts if there is no flux calibration.')
head.add_comment('The absence of a 4th column shows that the latter case is in')
head.add_comment('effect along with the units which are set to ADU. The')
head.add_comment('wavelengths are corrected to the heliocentre and follow')
head.add_comment('smooth polynomials with pixel number. The header item NARC')
head.add_comment('shows the number of coefficients used. If NARC=-2, then the')
head.add_comment('scale was logarithmic. NARC=-2 or 2 indicates that the spectra')
head.add_comment('have been rebinned at least once. Anything else, e.g. NARC=4,')
head.add_comment('means that the data are on their raw wavelength scale.')

hdul = [phdu,]
for spec in specs:
    c1 = pyfits.Column(name='Wavelength', format='D', unit='A', array=spec.x.data)
    counts = spec.cfrat.min() == 1 and spec.cfrat.max() == 1
    if counts:
        c2 = pyfits.Column(name='Flux', format='E', unit='ADU', array=spec.y.data)
        c3 = pyfits.Column(name='Ferr', format='E', unit='ADU', array=spec.y.errors)
    elif options.flam:
        conv = 1.e-15*subs.C/spec.x.dat**2
        c2 = pyfits.Column(name='Flux', format='E', unit='ergs/cm**2/s/A', array=conv*spec.y.data)
        c3 = pyfits.Column(name='Ferr', format='E', unit='ergs/cm**2/s/A', array=conv*spec.y.errors)
        c4 = pyfits.Column(name='Fcratio', format='E', unit='ergs/cm**2/s/A/count', array=conv*spec.cfrat)
    else:
        c2 = pyfits.Column(name='Flux', format='E', unit='mJy', array=spec.y.data)
        c3 = pyfits.Column(name='Ferr', format='E', unit='mJy', array=spec.y.errors)
        c4 = pyfits.Column(name='Fcratio', format='E', unit='mJy/count', array=spec.cfrat)

    if counts:
        coldefs = pyfits.ColDefs([c1,c2,c3])
    else:
        coldefs = pyfits.ColDefs([c1,c2,c3,c4])
    thdu = pyfits.new_table(coldefs)
    update_header(spec, thdu.header)
    hdul.append(thdu)

hdulist = pyfits.HDUList(hdul)
hdulist.writeto(fname)
