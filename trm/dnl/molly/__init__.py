"""
molly read/write routines
"""

import math, copy, struct
from scipy import polyfit
import astropy.io.fits as fits
import numpy as np
import trm.subs as subs
import trm.dnl as dnl

class Molly(dnl.Dset):
    """
    Class for molly data. Adds a flux/count ratio array to the usual data
    contained in a Dset
    """

    def __init__(self, xa, ya, head=None, good=None, cfrat=None, narc=0, arc=None):
        """
        Create a new Molly spectrum. Note that header items of REAL*4 type can
        be made by enforcing that they are of numpy.float32 type

        xa     -- X axis
        ya     -- Y axis
        head   -- header
        good   -- Good data mask
        cfrat  -- counts/mJy array
        narc   -- number of arc coefficients
        arc    -- the arc coefficients
        """
        super(Molly, self).__init__(xa,ya,head,good,dnl.Dset.BAR)

        if cfrat == None:
            self.cfrat = np.empty(len(xa),np.float32)
            self.cfrat[:] = 1.0
        else:
            self.cfrat = copy.copy(cfrat)
        self.narc = narc
        self.arc  = arc

    def wmolly(self, mf, namax=6, dpmax=0.01, hasarc=False, output=False, full_output=False):
        """
        Writes out a Molly to an opened molly file

        This assumes that the x-scale is in Angstroms and the
        fluxes are milliJanskys.

        mf : opened file object (binary format)

        namax : maximum number of arc coefficients to try, 1 to force 
                velocity scale attempt only

        dpmax : maximum acceptable deviation in pixels.

        hasarc : indicates that there are arc coefficients that can be
                 used directly.
        """
        if namax < 1:
            raise ValueError('dnl.Dset.wmolly: namax < 1')
        if dpmax <= 0.:
            raise ValueError('dnl.Dset.wmolly: dpmax <= 0')

        fbytes = 44
        if self.y.has_errors:
            fcode = 3
        else:
            fcode = 4
        units  = '%-16s' % 'MILLIJANSKYS'
        npix   = len(self.x)
        xpix   = np.linspace(1,npix,npix)/float(npix)

        if hasarc:
            narc = self.narc
            arc  = self.arc
        else:
            # Attempt to fit arc scale by starting with a linear fit
            # and then successively increasing the number of terms to
            # a preset maximum. Correct to telluric scale according to
            # value of 'Vearth'
            narc   = 1
            dmin = 1.e30
            dmax   = 2.*dpmax
            mdisp  = abs((self.x.data.max()-self.x.data.min()))/(npix-1)
            while dmax > dpmax and narc < namax:
                arc = polyfit(xpix, self.x.data, narc)
                p   = np.poly1d(arc)
                dmax = abs(self.x.data-p(xpix)).max()/mdisp
                dmin = min(dmin, dmax)
                narc += 1
                if full_output:
                    print 'npix, dmax, narc, arc = ',npix,dmax,narc,arc

            if 'Vearth' in self.head:
                arc /= (1.-self.head['Vearth']/(subs.C/1000.))

            if dmax > dpmax:
                print 'trying log scale'
                # One go at a log scale
                lw = np.log(self.x.data)
                mdisp  = abs((lw.max()-lw.min())/(len(lw)-1))
                arc = polyfit(xpix, lw, 2)
                p   = np.poly1d(arc)
                dmax = abs(self.x.data-p(xpix)).max()/mdisp
                dmin = min(dmin, dmax)
                if dmax > dpmax:
                    raise MollyError(
                        'wmolly: could not fit an accurate poly to the X array for output to molly.\n' +
                        'Maximum # coeffs tried = ' + str(namax) + '\n' +
                        'Best achieved = ' + str(dmin) + ' pixels.')
                narc = -2
                if 'Vearth' in self.head:
                    arc -= math.log((1.-self.head['Vearth']/(subs.C/1000.)))
            if output:
                print 'Maximum deviation = ',dmax,'pixels, narc =',narc

            # ensure double precision
            arc = np.cast['float'](arc[::-1])

            if full_output:
                print 'final arc = ',arc

        # split up headers
        chead = subs.Odict()
        nchar = 0
        for k,v in self.head.iteritems():
            if isinstance(v, str):
                chead['%-16s' % k[:min(len(k),16)]] = '%-32s' % v[:min(len(v),32)]
                nchar += 1

        dhead  = subs.Odict()
        ndoub  = 0
        for k,v in self.head.iteritems():
            if isinstance(v, float):
                dhead['%-16s' % k[:min(len(k),16)]] = v
                ndoub += 1

        ihead  = subs.Odict()
        nint   = 0
        for k,v in self.head.iteritems():
            if isinstance(v, int):
                ihead['%-16s' % k[:min(len(k),16)]] = v
                nint += 1

        fhead  = subs.Odict()
        nfloat = 0
        for k,v in self.head.iteritems():
            if isinstance(v, np.float32):
                fhead['%-16s' % k[:min(len(k),16)]] = v
                nfloat += 1

        # first record
        mf.write(struct.pack('2i16s7i',fbytes,fcode,units,npix,narc,nchar,ndoub,nint,nfloat,fbytes))

        # second record
        fbytes = 16*(nchar+ndoub+nint+nfloat)
        mf.write(struct.pack('i',fbytes))
        if chead:
            for key in chead.keys():
                mf.write(struct.pack('16s',key))
        if dhead:
            for key in dhead.keys():
                mf.write(struct.pack('16s',key))
        if ihead:
            for key in ihead.keys():
                mf.write(struct.pack('16s',key))
        if fhead:
            for key in fhead.keys():
                mf.write(struct.pack('16s',key))
        mf.write(struct.pack('i',fbytes))

        # third record
        fbytes = 32*nchar + 8*ndoub + 4*nint + 4*nfloat
        mf.write(struct.pack('i',fbytes))
        if chead:
            for key in chead.keys():
                mf.write(struct.pack('32s',chead[key]))
        if dhead:
            for key in dhead.keys():
                mf.write(struct.pack('d',dhead[key]))
        if ihead:
            for key in ihead.keys():
                mf.write(struct.pack('i',ihead[key]))
        if fhead:
            for key in fhead.keys():
                mf.write(struct.pack('f',fhead[key]))
        mf.write(struct.pack('i',fbytes))

        # fourth record
        fbytes = 8*abs(narc)
        mf.write(struct.pack('i',fbytes))
        arc.tofile(mf)
        mf.write(struct.pack('i',fbytes))

        # fifth record
        if fcode == 3:
            fbytes = 12*npix
            counts = self.y.data * self.cfrat
            errors = self.y.errors * self.cfrat
            flux   = self.y.data.copy()

            flux[counts == 0.] = self.cfrat[counts == 0.]
            mf.write(struct.pack('i',fbytes))
            np.cast['float32'](counts).tofile(mf)
            np.cast['float32'](errors).tofile(mf)
            np.cast['float32'](flux).tofile(mf)
            mf.write(struct.pack('i',fbytes))
        elif fcode == 4:
            fbytes = 4*npix
            mf.write(struct.pack('i',fbytes))
            np.cast['float32'](self.y.data).tofile(mf)
            mf.write(struct.pack('i',fbytes))

        else:
            raise MollyError('fcode = ' + str(fcode) + ' not implemented')

def grmolly(fname):
    """Generator for reading a molly file"""
    mf = open(fname, 'rb')
    nspec = 1
    try:
        mspec = rmspec(mf)
        while mspec != None:
            yield (nspec, mspec)
            mspec = rmspec(mf)
            nspec += 1
        mf.close()
    except MollyError, err:
        mf.close()
        raise MollyError('Error while reading molly file = ' + fname + ': ' + str(err))

def rmolly(fname):
    """Reads entire molly file, a returning list of Molly objects, one for
    each spectrum in the file"""
    mf = open(fname, 'rb')
    mlist = []
    for nspec,mspec in grmolly(fname):
        mlist.append(mspec)
    return mlist

def rmspec(mf):
    """Reads next molly spectrum from an opened file object mf, returns None
    if none found"""
    tup = _read_molly_head(mf)
    if tup != None:
        fcode, head, x, narc, arc, border = tup
        y,cfrat = _read_molly_data(mf,fcode, len(x), border)
        return Molly(x, y, head, cfrat=cfrat, narc=narc, arc=arc)
    else:
        return None

def _read_molly_head(mf):
    """
    Reads headers and arc (if present) of a molly spectrum and sets the X
    array Returns (fcode, head, xaxis, narc, arc, border) where fcode is the
    molly format code needed for reading the data, head is an Odict of the
    molly header, xaxis represents the xaxis, narc is the number of arc
    coefficients, arc are the arc coefficients, and border defines the byte
    order (either '>' or '<'). Returns None if no spectrum is found.
    """

    # If 'fbytes' in the next line comes up blank, we have reached the end of
    # the file
    fbytes = mf.read(4)
    if fbytes == '': return None

    # If it does not start with 44 in either big or little endian form,
    # something is wrong
    (nbyte,)  = struct.unpack('<i', fbytes)
    if nbyte != 44:
        (nbyte,)  = struct.unpack('>i', fbytes)
        if nbyte != 44:
            raise MollyError('_read_molly_header: not a molly spectrum: first 4 bytes = ' + str(nbyte) + ' not 44')
        border = '>'
    else:
        border = '<'

    # Read first line with various format items
    try:
        fcode,units,npix,narc,nchar,ndoub,nint,nfloat = \
            struct.unpack(border + 'i16s6i',mf.read(44))
    except:
        raise MollyError("Failed to read first line of molly spectrum")

    # skip bytes at end of first record and at start of second
    mf.seek(8,1)

    # read names of string header items
    cnames = []
    for i in range(nchar):
        name = mf.read(16).strip()
        cnames.append(name)

    # read names of double header items
    dnames = []
    for i in range(ndoub):
        name = mf.read(16).strip()
        dnames.append(name)

    # read names of integer header items
    inames = []
    for i in range(nint):
        name = mf.read(16).strip()
        inames.append(name)

    # read names of float header items
    fnames = []
    for i in range(nfloat):
        name = mf.read(16).strip()
        fnames.append(name)

    # skip bytes at end of second record and at start of third
    mf.seek(8,1)

    # create header
    head = fits.Header()

    for i in range(nchar):
        value = mf.read(32).strip()
        head['hierarch ' + cnames[i]] = value

    dvals = struct.unpack(border + str(ndoub) + 'd', mf.read(8*ndoub))
    for i in range(ndoub):
        head['hierarch ' + dnames[i]] = dvals[i]

    ivals = struct.unpack(border + str(nint) + 'i', mf.read(4*nint))
    for i in range(nint):
        head['hierarch ' + inames[i]] = ivals[i]

    fvals = struct.unpack(border + str(nfloat) + 'f', mf.read(4*nfloat))
    for i in range(nfloat):
        head['hierarch ' + fnames[i]] = np.float32(fvals[i])

    # skip bytes at end of third record and at start of fourth
    mf.seek(8,1)

    # set X array
    if narc != 0:
        arc = np.fromfile(file=mf, dtype=border + 'f8', count=abs(narc))
        x   = np.polyval(arc[::-1], np.arange(1.,npix+1,1.)/npix)
        if narc < 0:
            x = np.exp(x)
        # correct to heliocentric scale
        if 'Vearth' in head:
            x *= (1.-head['Vearth']/(subs.C/1000.))
            head['comment'] = 'Wavelength scale is heliocentric'
        else:
            head['comment'] = 'Wavelength scale is possibly telluric'
    else:
        x = np.arange(1.,npix+1,1.)
        arc = None

    # skip 4 bytes at end of headers
    mf.seek(4,1)

    return (fcode, head, dnl.Axis('Wavelength', '\A', x), narc, arc, border)

def _read_molly_data(mf, fcode, npix, border):
    """
    (yaxis,fratio) = _read_molly_data(mf, fcode, npix, border)

    mf    -- file object
    fcode -- molly format code
    npix  -- number of pixels
    border-- string defining the byte order

    Reads data of a molly spectrum, assuming the header and arc have been read
    Return y data as an Axis and an array of flux/count ratios (which can be None)
    """
    # skip 4 bytes at start
    mf.seek(4,1)

    cfrat = None

    if fcode == 1:
        y = np.fromfile(file=mf, dtype=border + 'f4', count=npix)
        e = None
        ylabel = 'Counts'
        yunits = ''

    elif fcode == 2:
        y = np.fromfile(file=mf, dtype=border + 'f4', count=npix)
        e = np.fromfile(file=mf, dtype=border + 'f4', count=npix)
        ylabel = 'Counts'
        yunits = ''

    elif fcode == 3:
        counts  = np.fromfile(file=mf, dtype=border + 'f4', count=npix)
        errors  = np.fromfile(file=mf, dtype=border + 'f4', count=npix)
        flux    = np.fromfile(file=mf, dtype=border + 'f4', count=npix)

        cfrat      = np.empty(npix, dtype=border + 'f4')
        mod        = counts == 0.
        cfrat[mod] = flux[mod]
        mod        = counts != 0.
        cfrat[mod] = counts[mod] / flux[mod]

        e       = np.empty_like(errors)
        ok      = cfrat > 0.
        e[ok]   = errors[ok] / cfrat[ok]
        e[~ok]  = -1.
        y       = flux
        y[counts == 0.] = 0.

        ylabel = 'f\d\gn\u'
        yunits = 'mJy'

    elif fcode == 4:
        y = np.fromfile(file=mf, dtype=border + 'f4', count=npix)
        e = None
        ylabel = 'f\d\gn\u'
        yunits = 'mJy'

    elif fcode == 5:
        y = np.fromfile(file=mf, dtype=border + 'f4', count=npix)
        e = np.fromfile(file=mf, dtype=border + 'f4', count=npix)
        ylabel = 'f\d\gn\u'
        yunits = 'mJy'

    else:
        raise MollyError('_read_molly_data: invalid FCODE in molly spectrum = ' + str(fcode))
    
    # skip 4 bytes at end
    mf.seek(4,1)
    
    return (dnl.Axis(ylabel, yunits, y, e), cfrat)

def skip_molly(mf):
    """"
    Skips a molly spectrum, assuming that we are positioned at its start
    """

    # If 'fbytes' in the next line comes up blank, we have reached the end of the file
    fbytes = mf.read(4)
    if fbytes == '': return False
    
    # If it does not start with 44 in either big or little endian form, something is wrong
    (nbyte,)  = struct.unpack('<i', fbytes)
    if nbyte != 44:
        (nbyte,)  = struct.unpack('>i', fbytes)
        if nbyte != 44:
            raise MollyError('skip_molly: not a molly spectrum: first 4 bytes = ' + str(nbyte) + ' not 44')
        border = '>'
    else:
        border = '<'
    
    # Read first line with various format items
    try:
        (fcode,units,npix,narc,nchar,ndoub,nint,nfloat) = struct.unpack(border + 'i16s6i',mf.read(44))
    except:
        raise MollyError("skip_molly: failed to read first line of molly spectrum")

    # compute number of bytes to skip.
    nskip = 36 + 16*(nchar+ndoub+nint+nfloat)+32*nchar+8*ndoub+4*nint+4*nfloat+8*abs(narc)

    if fcode == 1:
        nskip += 4*npix
    elif fcode == 2:
        nskip += 8*npix
    elif fcode == 3:
        nskip += 12*npix
    elif fcode == 5:
        nskip += 8*npix
    else:
        raise MollyError('skip_molly: invalid FCODE in molly spectrum = ' + str(fcode))

    # skip them
    mf.seek(nskip,1)
    return True

# Exception class
class MollyError(dnl.DnlError):
    """For throwing exceptions from the dnl.molly module"""
    pass

