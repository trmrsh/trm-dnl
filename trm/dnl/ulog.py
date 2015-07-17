"""
defines routines useful for handling ULTRACAM log files

Class
=====

Ulog : represents all ULTRACAM data of a series of logfiles

UlogError : exception class

Functions
=========

rulog : reads one star's ultracam data directly from an ASCII log file
        to return a Dset
"""

from __future__ import division, absolute_import, print_function

import numpy as np
import astropy.io.fits as fits
from trm import subs

from .core import Axis, Dset, DnlError

def rulog(fname, nccd, aperture, type='c', emax=5, form='new'):
    """
    Reads ultracam data directly into a Dset.

    NB it is more efficient for multiple Dsets from one log file
    to read the log file first into an Ulog and then to extract
    Dsets from this.

    fname : string
            log file name

    nccd : int
           ccd number, starting from 1.

    aperture : int
               aperture number, starting from 1.

    type : string
           'c' for counts (no other type recognised at present)

    emax : int
           maximum error code

    form : string
           'new' or 'old' format. New style has fewer columns
    """

    if form == 'new':
        off = 0
    elif form == 'old':
        off = 1
    else:
        raise UlogError("rulog: did not recognize value of 'form' " +
                            "which must be either 'new' or 'old'")

    if type == 'c':
        if aperture < 1:
            raise UlogError('rulog: aperture < 1')
        offset = 14*(aperture-1)
        yc     = offset + 15 + off
        yec    = offset + 16 + off
        fc     = offset + 21 + off
        dset   = _rulog(fname, nccd, 'Count rate', 'cts/sec', yc,
                        yec, off, fc, emax)
        dset.y.data /= (2.*86400.)*dset.x.errors
        if dset.y.has_errors:
            dset.y.errors /= (2.*86400.)*dset.x.errors
        dset['NCCD']     = nccd
        dset['Aperture'] = aperture
    else:
        raise UlogError('rulog: type = ' + type + ' unrecognised.')

    dset['Type']     = 'time series'
    dset['Time']     = 'MJD (UTC)'
    return dset

def _rulog(fname, nccd, yaxis, yunits, yc, yec, off, fc, emax):
    """
    Reads ultracam data in, returning a Dset.

    fname : string
            file name

    nccd : int
           ccd number

    yaxis : string
            name for y axis data

    yunits : string
            units for y axis data

    yc : int
         column for y data

    yec : int
          column for y errors. 0 to ignore.

    off : int
          offset, 0 or 1, for new or old log file formats

    fc : int
         column for error flag. 0 to ignore.

    emax : int
           maximum error code
    """

    if yc < 1:
        raise UlogError('_rulog: yc < 1')
    if nccd < 1:
        raise UlogError('_rulog: nccd < 1')

    fin  = open(fname)
    xc   = 1
    xec  = 3+off
    yc  -= 1
    yec -= 1
    fc  -= 1
    x    = []
    xe   = []
    y    = []
    ye   = []
    flag = []
    cmax = max(4+off, yc, yec)
    for line in fin:
        if line[0:1] != '#' and not line.isspace():
            svar = line.split()
            if len(svar) <= cmax:
                raise UlogError('_rulog: too few columns of data in ' +
                                    fname)

            if int(svar[4+off]) == nccd:
                x.append(svar[xc])
                xe.append(svar[xec])
                y.append(svar[yc])
                if yec > -1: ye.append(svar[yec])
                if fc > -1: flag.append(svar[fc])
    fin.close()
    if len(x) == 0:
        raise UlogError('_rulog: no data loaded from ' + fname)

    # set data types, convert to arrays
    x    = np.asarray(x, np.float64)
    xe   = np.asarray(xe, np.float32)/(2.*86400.)
    y    = np.asarray(y, np.float32)
    ye   = np.asarray(ye, np.float32) if yec > -1 else None
    flag = np.asarray(flag, np.int) if fc > -1 else None

    # identify good data
    if yec > -1:
        ygood = ye > 0.
    else:
        ygood = None

    if fc > -1:
        good = flag <= emax
    else:
        good = None

    return Dset(Axis('MJD (UTC)', 'days', x, xe),
                Axis(yaxis, yunits, y, ye, ygood), good=good)

class Ulog(object):

    """
    Class for ULTRACAM data analysis. Contains all data of an ULTRACAM log
    file or files.

    Every attribute is a dictionary keyed by CCD number (e.g. 1, 2 or 3 for
    ULTRACAM). As it is a dictionary, it is possible to have an entry for
    CCD 2 but not 1 for instance.

    The following attributes each return an array, with data type given at
    the end of the description lines::

      utc : dict
           UTC times at mid exposure, in MJD [float64]

      tflag : dict
           reliability of times [bool]

      expose : dict
           exposure time, seconds [float32]

      fwhm : dict
           FWHM in pixels [float32]

      beta : dict
           Moffat beta parameter [float32]

    The next attributes each returns a list of arrays, with one entry per
    aperture used for the CCD in question. Thus x[2][1] gives the x positions
    for the second aperture (usual C-like start at 0 convention for lists) of
    CCD 3.

      x : dict
           x positions used [float32]

      y : dict
          y positions used [float32]

      xm : dict
          x positions measured (0 for invalid or linked apertures) [float32]

      ym : dict
          y positions measured (0 for invalid or linked apertures) [float32]

      exm : dict
          1-sigma uncertainties in measured x positions (-1 if invalid /
          linked) [float32]

      eym : dict
          1-sigma uncertainties in measured y positions (-1 if invalid /
          linked) [float32]

      counts : dict
          counts in the aperture [float32]

      sigma : dict
          1-sigma uncertainties on counts in the aperture [float32]

      sky : dict
          sky background, counts per pixel [float32]

      nsky : dict
          number of sky pixels [int]

      nrej : dict
          number of pixels rejected from the sky [int]

      worst : dict
          worst bad pixel in aperture (exact meaning of this is down to the
          user) [int]

      eflag : dict
          error flag (see ultracam reduce file for meanings) [int]

    Construct using::

      ult = Ulog(flist)

    where flist is a list of log files. These can be either ASCII '.log'
    files or FITS, but must be all of a type.
    """

    def __init__(self, files, form='new'):
        """
        Constructs an Ultracam from a file or files

        files : string / list
                file name or list of file names which will be read in order
                and should match in format, aperture numbers etc.
        """

        if form == 'new':
            off = 0
        else:
            off = 1

        nline = 0
        nccd  = set()
        naper = {}
        self.utc    = {}
        self.tflag  = {}
        self.expose = {}
        self.fwhm   = {}
        self.beta   = {}
        self.x      = {}
        self.y      = {}
        self.xm     = {}
        self.ym     = {}
        self.exm    = {}
        self.eym    = {}
        self.counts = {}
        self.sigma  = {}
        self.sky    = {}
        self.nsky   = {}
        self.nrej   = {}
        self.worst  = {}
        self.eflag  = {}
        found_all_ccds = False

        if isinstance(files, str):
            files = [files]

        if files[0].endswith('.log'):
            ftype = 'log'
        elif files[0].endswith('.fits') or files[0].endswith('.fit') or \
             files[0].endswith('.fits.gz') or files[0].endswith('.fit.gz'):
            ftype = 'fits'
        else:
            raise UlogError('Ulog(files): did not recognize file type of' +
                            ' first one = ' + files[0])

        for fname in files[1:]:
            if files[0].endswith('.log') and ftype == 'fits':
                raise UlogError('Ulog(files): clashing file type. Expected' +
                                ' an ASCII log file but got = ' + fname)
            elif (files[0].endswith('.fits') or files[0].endswith('.fit') or
                  files[0].endswith('.fits.gz') or
                  files[0].endswith('.fit.gz')) and ftype == 'log' :
                raise UlogError('Ulog(files): clashing file type. Expected' +
                                ' a FITS file but got = ' + fname)

        if ftype == 'log':

            for fname in files:

                fin = open(fname)
                for line in fin:
                    nline += 1
                    if line[0:1] != '#' and not line.isspace():
                        svar = line.split()

                        # we accumulate apertures numbers for each new CCD
                        # encounter, but if we re-find a CCD, we check that
                        # aperture numbers match. Also check that extra CCDs
                        # are not found after all were thought to have been
                        # found
                        if (len(svar) - 7 - off ) % 14 > 0:
                            raise UlogError('Ulog.__init__: incorrect number' +
                                            ' of entries in line ' +
                                            str(nline) + ' of ' + fname)
                        nc  = int(svar[4+off])
                        nap = (len(svar) - 7 - off ) // 14
                        if nc in nccd:
                            if nap != naper[nc]:
                                raise UlogError('Ulog.__init__: incorrect' +
                                                ' number of apertures in' +
                                                ' line ' + str(nline) +
                                                ' of ' + fname)
                            found_all_ccds = True
                        elif found_all_ccds:
                            raise UlogError('Ulog.__init__: new CCD was ' +
                                            'found even though all were ' +
                                            'thought to be found in line ' +
                                            str(nline) + ' of ' + fname)
                        else:
                            nccd.add(nc)
                            naper[nc] = nap

                            # initialise the lists for this CCD
                            self.utc[nc]    = []
                            self.tflag[nc]  = []
                            self.expose[nc] = []
                            self.fwhm[nc]   = []
                            self.beta[nc]   = []
                            self.x[nc]      = [[] for i in range(nap)]
                            self.y[nc]      = [[] for i in range(nap)]
                            self.xm[nc]     = [[] for i in range(nap)]
                            self.ym[nc]     = [[] for i in range(nap)]
                            self.exm[nc]    = [[] for i in range(nap)]
                            self.eym[nc]    = [[] for i in range(nap)]
                            self.counts[nc] = [[] for i in range(nap)]
                            self.sigma[nc]  = [[] for i in range(nap)]
                            self.sky[nc]    = [[] for i in range(nap)]
                            self.nsky[nc]   = [[] for i in range(nap)]
                            self.nrej[nc]   = [[] for i in range(nap)]
                            self.worst[nc]  = [[] for i in range(nap)]
                            self.eflag[nc]  = [[] for i in range(nap)]

                        # squirrel the data away
                        self.utc[nc].append(svar[1])
                        self.tflag[nc].append(svar[2])
                        self.expose[nc].append(svar[3+off])
                        self.fwhm[nc].append(svar[5+off])
                        self.beta[nc].append(svar[6+off])
                        for i in range(nap):
                            offset = 14*i + 6 + off
                            self.x[nc][i].append(svar[offset+2])
                            self.y[nc][i].append(svar[offset+3])
                            self.xm[nc][i].append(svar[offset+4])
                            self.ym[nc][i].append(svar[offset+5])
                            self.exm[nc][i].append(svar[offset+6])
                            self.eym[nc][i].append(svar[offset+7])
                            self.counts[nc][i].append(svar[offset+8])
                            self.sigma[nc][i].append(svar[offset+9])
                            self.sky[nc][i].append(svar[offset+10])
                            self.nsky[nc][i].append(svar[offset+11])
                            self.nrej[nc][i].append(svar[offset+12])
                            self.worst[nc][i].append(svar[offset+13])
                            self.eflag[nc][i].append(svar[offset+14])

                fin.close()

            # Transform to numpy arrays of correct type
            for nc in nccd:

                self.utc[nc]    = np.asarray(self.utc[nc], np.float64)
                self.tflag[nc]  = np.asarray(self.tflag[nc], np.bool)
                self.expose[nc] = np.asarray(self.expose[nc], np.float32)
                self.fwhm[nc]   = np.asarray(self.fwhm[nc], np.float32)
                self.beta[nc]   = np.asarray(self.beta[nc], np.float32)

                for nap in range(naper[nc]):
                    self.x[nc][nap]      = np.asarray(self.x[nc][nap],
                                                      np.float32)
                    self.y[nc][nap]      = np.asarray(self.y[nc][nap],
                                                      np.float32)
                    self.xm[nc][nap]     = np.asarray(self.xm[nc][nap],
                                                      np.float32)
                    self.ym[nc][nap]     = np.asarray(self.ym[nc][nap],
                                                      np.float32)
                    self.exm[nc][nap]    = np.asarray(self.exm[nc][nap],
                                                      np.float32)
                    self.eym[nc][nap]    = np.asarray(self.eym[nc][nap],
                                                      np.float32)
                    self.counts[nc][nap] = np.asarray(self.counts[nc][nap],
                                                      np.float32)
                    self.sigma[nc][nap]  = np.asarray(self.sigma[nc][nap],
                                                      np.float32)
                    self.sky[nc][nap]    = np.asarray(self.sky[nc][nap],
                                                      np.float32)
                    self.nsky[nc][nap]   = np.asarray(self.nsky[nc][nap],
                                                      np.int)
                    self.nrej[nc][nap]   = np.asarray(self.nrej[nc][nap],
                                                      np.int)
                    self.worst[nc][nap]  = np.asarray(self.worst[nc][nap],
                                                      np.int)
                    self.eflag[nc][nap]  = np.asarray(self.eflag[nc][nap],
                                                      np.int)

        elif ftype == 'fits':

            for fname in files:
                hdulist = fits.open(fname)

                for n in range(1,len(hdulist)):

                    thead = hdulist[n].header
                    nc  = thead['NCCD']
                    nap = (thead['TFIELDS'] - 5 ) // 13

                    if nc in nccd:
                        # append to lists
                        if nap != naper[nc]:
                            raise UlogError('Ulog.__init__: incorrect' +
                                            ' number of apertures in ' + fname)
                        found_all_ccds = True

                        self.utc[nc] = np.concatenate((self.utc[nc],
                                                       tdata.field('MJD')))
                        self.tflag[nc] = np.concatenate((self.tflag[nc],
                                                         tdata.field('Flag')))
                        self.expose[nc] = np.concatenate((self.expose[nc],
                                                          tdata.field('Expose')))
                        self.fwhm[nc] = np.concatenate((self.fwhm[nc],
                                                        tdata.field('FWHM')))
                        self.beta[nc] = np.concatenate((self.beta[nc],
                                                        tdata.field('beta')))
                        for i in range(nap):
                            lbl = '_' + str(i+1)
                            self.x[nc][i] = np.concatenate(
                                (self.x[nc][i],tdata.field('X' + lbl)))
                            self.y[nc][i] = np.concatenate(
                                (self.y[nc][i], tdata.field('Y' + lbl)))
                            self.xm[nc][i] = np.concatenate(
                                (self.xm[nc][i], tdata.field('XM' + lbl)))
                            self.ym[nc][i] = np.concatenate(
                                (self.ym[nc][i], tdata.field('YM' + lbl)))
                            self.exm[nc][i] = np.concatenate(
                                (self.exm[nc][i],tdata.field('EXM' + lbl)))
                            self.eym[nc][i] = np.concatenate(
                                (self.eym[nc][i],tdata.field('EYM' + lbl)))
                            self.counts[nc][i] = np.concatenate(
                                (self.counts[nc][i],
                                 tdata.field('Counts' + lbl)))
                            self.sigma[nc][i] = np.concatenate(
                                (self.sigma[nc][i],tdata.field('Sigma' + lbl)))
                            self.sky[nc][i]    = np.concatenate(
                                (self.sky[nc][i],tdata.field('Sky' + lbl)))
                            self.nsky[nc][i]   = np.concatenate(
                                (self.nsky[nc][i],tdata.field('Nsky' + lbl)))
                            self.nrej[nc][i]   = np.concatenate(
                                (self.nrej[nc][i],tdata.field('Nsky' + lbl)))
                            self.worst[nc][i]  = np.concatenate(
                                (self.worst[nc][i],tdata.field('Worst' + lbl)))
                            self.eflag[nc][i]  = np.concatenate(
                                (self.eflag[nc][i],tdata.field('Eflag' + lbl)))

                    elif found_all_ccds:
                        raise UlogError('Ulog.__init__: new CCD was found ' +
                                        'even though all were thought to be' +
                                        'found in ' + fname)
                    else:

                        # initialise the lists
                        nccd.add(nc)
                        naper[nc] = nap
                        tdata = hdulist[n].data

                        self.utc[nc]    = tdata.field('MJD')
                        self.tflag[nc]  = tdata.field('Flag')
                        self.expose[nc] = tdata.field('Expose')
                        self.fwhm[nc]   = tdata.field('FWHM')
                        self.beta[nc]   = tdata.field('beta')
                        self.x[nc]      = [tdata.field('X_' + str(i+1))
                                           for i in range(nap)]
                        self.y[nc]      = [tdata.field('Y_' + str(i+1))
                                           for i in range(nap)]
                        self.xm[nc]     = [tdata.field('XM_' + str(i+1))
                                           for i in range(nap)]
                        self.ym[nc]     = [tdata.field('YM_' + str(i+1))
                                           for i in range(nap)]
                        self.exm[nc]    = [tdata.field('EXM_' + str(i+1))
                                           for i in range(nap)]
                        self.eym[nc]    = [tdata.field('EYM_' + str(i+1))
                                           for i in range(nap)]
                        self.counts[nc] = [tdata.field('Counts_' + str(i+1))
                                           for i in range(nap)]
                        self.sigma[nc]  = [tdata.field('Sigma_' + str(i+1))
                                           for i in range(nap)]
                        self.sky[nc]    = [tdata.field('Sky_' + str(i+1))
                                           for i in range(nap)]
                        self.nsky[nc]   = [tdata.field('Nsky_' + str(i+1))
                                           for i in range(nap)]
                        self.nrej[nc]   = [tdata.field('Nsky_' + str(i+1))
                                           for i in range(nap)]
                        self.worst[nc]  = [tdata.field('Worst_' + str(i+1))
                                           for i in range(nap)]
                        self.eflag[nc]  = [tdata.field('Eflag_' + str(i+1))
                                           for i in range(nap)]

                hdulist.close()

    def __iadd__(self, other):
        """
        += in-place addition to add one Ulog onto the end of another

        ccd numbers and aperture numbers for each CCD must match. This is useful
        when reading in a sequence of log files.
        """

        if not isinstance(other, Ulog):
            raise UlogError('Ulog.__iadd__: can' +
                            ' only add another Ulog to an Ulog.')

        nccd = set(self.utc.keys())
        if set(other.utc.keys()) != nccd:
            raise UlogError('Ulog.__iadd__: CCD numbers of inputs do not match')

        for nc in nccd:
            if len(self.x[nc]) != len(other.x[nc]):
                raise UlogError('Ulog.__iadd__: incompatible' +
                                ' aperture numbers for CCD ' + nc)

        # OK, tests passed, add new arrays onto the end of the old ones
        for nc in nccd:

            self.utc[nc] = np.concatenate((self.utc[nc], other.utc[nc]))
            self.tflag[nc] = np.concatenate((self.tflag[nc], other.tflag[nc]))
            self.expose[nc] = np.concatenate((self.expose[nc],other.expose[nc]))
            self.fwhm[nc] = np.concatenate((self.fwhm[nc],other.fwhm[nc]))
            self.beta[nc] = np.concatenate((self.beta[nc],other.beta[nc]))

            for nap in range(len(self.x[nc])):
                self.x[nc][nap] = np.concatenate(
                    (self.x[nc][nap],other.x[nc][nap]))
                self.y[nc][nap] = np.concatenate(
                    (self.y[nc][nap],other.y[nc][nap]))
                self.xm[nc][nap] = np.concatenate(
                    (self.xm[nc][nap],other.xm[nc][nap]))
                self.ym[nc][nap] = np.concatenate(
                    (self.ym[nc][nap],other.ym[nc][nap]))
                self.exm[nc][nap] = np.concatenate(
                    (self.exm[nc][nap],other.exm[nc][nap]))
                self.eym[nc][nap] = np.concatenate(
                    (self.eym[nc][nap],other.eym[nc][nap]))
                self.counts[nc][nap] = np.concatenate(
                    (self.counts[nc][nap],other.counts[nc][nap]))
                self.sigma[nc][nap] = np.concatenate(
                    (self.sigma[nc][nap],other.sigma[nc][nap]))
                self.sky[nc][nap] = np.concatenate(
                    (self.sky[nc][nap],other.sky[nc][nap]))
                self.nsky[nc][nap] = np.concatenate(
                    (self.nsky[nc][nap],other.nsky[nc][nap]))
                self.nrej[nc][nap] = np.concatenate(
                    (self.nrej[nc][nap],other.nrej[nc][nap]))
                self.worst[nc][nap] = np.concatenate(
                    (self.worst[nc][nap],other.worst[nc][nap]))
                self.eflag[nc][nap] = np.concatenate(
                    (self.eflag[nc][nap],other.eflag[nc][nap]))
        return self

    def tseries(self, nccd, naper, ttype='c', wmax=50, emax=5):
        """
        Return a time series as a Dset from an Ulog. The exposure
        times are stored in the "errors" assigned to the X axis
        (after diviasion by 2 so that +/- the error = exposure time).

          nccd : int
                CCD number

          naper : int
                aperture number (starting from 1)

          ttype : string
                'c' counts, 'xm' measured x, 'ym' measured y, 'f' fwhm

          wmax : int
                maximum bad pixel in aperture above which data will be
                flagged as bad (only applies to aperture data).

          emax : int
                maximum error flag, above which data will be flagged as
                bad (only applies to aperture data).
        """
        if nccd not in self.utc:
            raise UlogError('Ulog.tseries: nccd = ' + str(nccd) +
                            ' not found.')
        if naper < 1 or naper > len(self.x[nccd]):
            raise UlogError('Ulog.tseries: naper = ' + str(naper) +
                            ' not found in CCD = ' + str(nccd))
        if ttype != 'c' and ttype != 'xm' and ttype != 'ym' and ttype != 'f':
            raise UlogError('Ulog.tseries: ttype = ' + ttype +
                            ' is not a recognised times series type.')

        nap = naper - 1

        # construct the Dset. First the x-axis, common to all.
        xaxis = Axis('MJD (UTC)', 'days', self.utc[nccd],
                     self.expose[nccd]/86400./2., self.tflag[nccd])

        # common header items
        head = fits.Header()
        head['TYPE']     = ('time series', 'Type of data')
        head['INSTRUME'] = ('ULTRACAM', 'Instrument')
        head['NCCD']     = (nccd, 'CCD number')

        # next the yaxis and data mask arrays
        if ttype == 'c':
            ye    = self.sigma[nccd][nap]
            ygood = ye > 0.
            y     = self.counts[nccd][nap]
            y[ygood]  /= self.expose[nccd][ygood]
            ye[ygood] /= self.expose[nccd][ygood]
            yaxis = Axis('Count rate', 'cts/sec', y, ye, ygood)

            good   = (self.worst[nccd][nap] <= wmax) & \
                     (self.eflag[nccd][nap] <= emax)

            head['NAPERTUR'] = (naper, 'Aperture number')

        elif ttype == 'xm':
            yaxis = Axis('X position', 'pixels', self.xm[nccd][nap],
                         self.xme[nccd][nap])
            good   = (self.worst[nccd][nap] <= wmax) & \
                (self.eflag[nccd][nap] <= emax)

            head['NAPERTUR'] = (naper, 'Aperture number')

        elif ttype == 'ym':
            yaxis = Axis('Y position', 'pixels', self.ym[nccd][nap],
                         self.yme[nccd][nap])
            good   = (self.worst[nccd][nap] <= wmax) & \
                (self.eflag[nccd][nap] <= emax)

            head['NAPERTUR'] = (naper, 'Aperture number')

        elif ttype == 'f':
            yaxis = Axis('FWHM seeing', 'pixels', self.fwhm[nccd])
            good   = (self.worst[nccd][nap] <= wmax) & \
                (self.eflag[nccd][nap] <= emax)

        # create Dset
        dset = Dset(xaxis, yaxis, head, good)

        return dset

# Exception class
class UlogError(DnlError):
    """For throwing exceptions from the dnl.ulog module"""
    pass

