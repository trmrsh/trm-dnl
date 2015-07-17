"""
Various miscellaneous pieces of I/O from and to Dsets. Most important
are the rfits and wfits methods for accessing multiple Dset FITS data.
FITS is the "native format" for storage of Dsets.
"""
from __future__ import division, absolute_import, print_function

# built-ins
import os
import copy
import struct

# third party
import numpy as np
from astropy.io import fits

# mine
from trm import subs
from trm.subs import cpp
from trm import sla

# dnl-related stuff
from trm.dnl.core import Axis, Dset, DnlError

def rponto(fobj, endian=''):
    """
    Read from a ponto disk file

    fobj : file pointer
           positioned at start of a ponto data set. For ponto files
           this means that the first 4 bytes have been read and checked.

    endian : string
           endianness. Use cpp.open_ponto to get both fobj and endian
           if unsure. '' for PCs, '>' for the reverse.
    """

    head  = cpp.read_header(fobj, endian)
    (ptype,) = struct.unpack(endian + 'i', fobj.read(4))
    if ptype == 0 or 3:
        dptype = Dset.POINTS
    elif ptype == 1:
        dptype = Dset.LINE
    elif dptype == 2:
        dptype = Dset.BAR

    # read mask
    (npix,) = struct.unpack(endian + 'i', fobj.read(4))
    good    = ~np.fromfile(file=fobj, dtype='b1', count=npix)

    # read axes
    xlabel,xunits,x,xe = cpp.read_ponto_axis(fobj, endian)
    ylabel,yunits,y,ye = cpp.read_ponto_axis(fobj, endian)

    if xe is not None:
        xgood = xe > 0.
    else:
        xgood = None
    if ye is not None:
        ygood = ye > 0.
    else:
        ygood = None

    return Dset(Axis(xlabel,xunits,x,xe,xgood),
                Axis(ylabel,yunits,y,ye,ygood),head,good,ptype=dptype)

def rsdss(fname, air=True, conv=False):
    """sdss = rsdss(fname, air) -- reads an SDSS fits file

    fname : string
          name of FITS file

    air : bool
          if True, mark the wavelengths as air, else vacuum.

    conv : bool
          if True, apply a conversion to get the wavelengths correct. i.e.
          if air=True, and conv=True, the input wavelengths will be assumed
          to be vacuum and converted.  If conv=False, they will be assumed
          to be air from the start. If air=False and conv=True then the
          input wavelengths will be assumed to be air wavelengths and
          converted to vacuum.

    Returns a Dset

    """

    # Read data in one of three ways
    hdulist   = fits.open(fname)
    data      = hdulist[0].data
    if data is not None:
        fhead  = hdulist[0].header
        flux   = data[0,:]
        ferr   = data[2,:]
        lw0    = fhead['COEFF0']
        disp   = fhead['COEFF1']
        wave   = pow(10., lw0 + disp*np.arange(len(flux)))
    else:
        data  = hdulist[1].data
        fhead = hdulist[1].header
        wave  = data.field(0)
        flux  = data.field(1)
        ferr  = data.field(2)
    hdulist.close()

    yaxis  = Axis('f\d\gl\u', '10\u-17\d ergs/s/cm\u2\d/\A', flux, ferr)

    if air:
        if conv:
            wave = subs.vac2air(wave)
        xaxis = Axis('Wavelength', '\A', wave)
    else:
        if conv:
            wave = subs.air2vac(wave)
        xaxis = Axis('Vacuum wavelength', '\A', wave)

    head  = subs.Odict()
    head['ReadFrom'] = fname
    cards = fhead.ascardlist()
    for card in cards:
        head[card.key] = card.value

    return Dset(xaxis, yaxis, head)

def rascii(fname, x=1, xe=0, y=2, ye=3, xlabel='', xunits='',
           ylabel='', yunits=''):
    """
    Reads ASCII column data

    fname : string
          name of ASCII file

    x : int
          number of column representing X

    xe : int
          number of column representing errors on X, 0 to ignore

    y : int
          number of column representing Y

    ye : int
          number of column representing errors on Y, 0 to ignore

    xlabel : string
          label for X axis

    xunits : string
          units for X axis

    ylabel : string
          label for Y axis

    yunits : string
          units for X axis
    """

    data = np.loadtxt(fname)
    if x < 1 or x > data.shape[1]:
        raise DnlError('X column out of range 1 to ' + str(data.shape[1]))
    if y < 1 or y > data.shape[1]:
        raise DnlError('Y column out of range 1 to ' + str(data.shape[1]))
    if xe < 0 or xe > data.shape[1]:
        raise DnlError('X errors column out of range 0 to ' +
                       str(data.shape[1]))
    if ye < 0 or ye > data.shape[1]:
        raise DnlError('Y errors column out of range 0 to ' +
                       str(data.shape[1]))

    if xe == 0:
        xaxis = Axis(xlabel, xunits, data[:,x-1])
    else:
        xaxis = Axis(xlabel, xunits, data[:,x-1], data[:,xe-1])
    if ye == 0:
        yaxis = Axis(ylabel, yunits, data[:,y-1])
    else:
        yaxis = Axis(ylabel, yunits, data[:,y-1], data[:,ye-1])

    head = subs.Odict()
    head['ReadFrom'] = fname
    return Dset(xaxis, yaxis, head)

def wfits(fname, dsets, head=None, append=False, clobber=False):
    """
    writes or appends a set of Dsets to a multi-extension FITS file.
    Each Dset goes into one binary table extension.

    fname : string
          file name

    dsets : iterable
          returns the Dset objects which will be written to the file

    head : astropy.io.fits.Header
          any FITS header to write into the primary HDU.

    append: bool
          if True append if possible, else overwrite.

    clobber: bool
          clobber output files if not appending
    """

    # make iterable
    if isinstance(dsets, Dset): dsets = [dsets,]

    if os.path.exists(fname):
        if not os.path.isfile(fname):
            raise DnlError('io.wfits: ' + fname + ' is not a file')
    else:
        append = False

    if not append:
        if head is None:
            thead  = fits.Header()
        else:
            thead = head.copy()
        thead['comment'] = 'This file contains Dset objects written by trm.dnl.io.wfits'
        thead['comment'] = 'Each one contains X and Y data, and, optionally, associated'
        thead['comment'] = 'errors and logical mask arrays. The mask arrays associated'
        thead['comment'] = 'with X and Y are designed to mark invalid data that should'
        thead['comment'] = 'never be used, but sometimes is useful to keep in order to'
        thead['comment'] = 'synchronise arrays. There is also an optional overall mask'
        thead['comment'] = 'designed to flag, legitimate, but in someway bad, data.'
        thead['comment'] = 'The mask arrays are called "good" to iindicate that True'
        thead['comment'] = 'values mark good data.'

        prihdu = fits.PrimaryHDU(header=thead)
        prihdu.writeto(fname,clobber=clobber)

    fobj = open(fname, 'ab+')
    for dset in dsets:
        tbhdu = dset.toHDU()
        tbhdu.header['EXTNAME'] = 'Dset'
        fits.append(fobj, tbhdu.data, tbhdu.header)
    fobj.close()

def nfits(fname):
    """
    Returns with the number of Dsets in a FITS file written
    by wfits. (Simply returns number of HDUs-1, with no checks.)
    """
    hdul = fits.open(fname)
    ndset = len(hdul)-1
    hdul.close()
    return ndset

def rfits(fname, which=None):
    """
    Reads Dsets in from a multi-extension FITS file.

      fname : string
            file name to read

      which : list / None
            list of extension numbers to load, e.g. [1,5] will load
            the 1st and 5th Dsets in the file. An error will be raised
            if a request is made for a Dset that does not exist. If None,
            all Dsets will be returned
    """

    hdul  = fits.open(fname)
    ndset = len(hdul)
    dsets = []
    if which is None:
        for hdu in hdul[1:]:
            dsets.append(Dset.fromHDU(hdu))
    else:
        for idset in which:
            dsets.append(Dset.fromHDU(hdul[idset]))
    hdul.close()
    return dsets

def grfits(fname, first=1):
    """
    Reads Dsets in from a multi-extension FITS file, generator
    which can be used in expression of form::

      for dset in grfits('test.fits'):

    This is potentially useful when reading very large numbers of
    Dsets in when one would prefer not to store each one but simply
    process one by one.

    Arguments::

      fname : string
            file name to read

      first : int
            first Dset to read, starting at 1
    """

    hdul  = fits.open(fname)
    ndset = len(hdul)
    dsets = []
    for hdu in hdul[first:]:
        yield Dset.fromHDU(hdu)
    hdul.close()

