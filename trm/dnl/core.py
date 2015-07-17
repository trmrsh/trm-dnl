"""
This file contains core components of the dnl package such as Axis, Dset
and DnlError.

"""

from __future__ import division, absolute_import, print_function

import copy
import numpy as np
from astropy.io import fits

from trm import subs
from trm import sla

from trm.dnl.base import *
from trm.dnl.plot import *

class Axis(object):

    """Represents an axis with data, errors, array of boolean good/bad values, a
    label and units.  Note that it can often be useful to keep data which is
    known to be bad in order for example to ensure a smooth progression of
    wavelength with pixel number for instance. Axis objects support some
    standard arithmetic operations, propagation errors via assumed
    independence and linear approximation. It is up to the user to ensure that
    these are reasonable assumptions, and one should exercise great care in
    complex cases where it will often be better to handcraft errors rather
    than assume that they propagate safely.

    Attributes::

     data : 1D numpy.ndarray
            Main data array.

     errors : 1D numpy.ndarray or None
              Errors array, if not None must have same length as data

     good : 1D boolean numpy.ndarray or None
            Boolean mask array specifying which pixels are good. Bad data
            are ignored in most operations with the exception of returning
            a length and storing to disk.

     label : string
             Description of quantity the Axis contains, e.g. "Temperature"

     units : string
             Specifies the units of the axis, e.g. "K"

    """

    def __init__(self, *args):
        """
        Initialise an Axis.

        4 versions according to the number of arguments.

        1) 0 arguments -- null constructor

           Creates a dummy Axis with data, errors,
           good all set to None, label and units set blank.
           See later for when this might be useful. Usually it
           should be avoided. Note the object created will fail
           to pass the check method.

        2) 3 to 5 arguments starting with a string::

         label : string
                 Axis label

         units : string
                 Axis units

         data : numpy.ndarray
                1D array of data

         errors : numpy.ndarray or None
                  errors

         good    : OK data (optional, None if not defined)

        3) 2 to 4 arguments starting with an Axis::

         axis    : another axis to provide the label and units
         data    : the data
         errors  : errors (optional, None if not defined)
         good    : OK data (optional, None if not defined)

        4) 1 to 3 arguments starting with an array::

          data    : the data
          errors  : errors (optional, None if not defined)
          good    : OK data (optional, None if not defined)

        NB The data is "data", "errors" and "good" are *copied*, not just
        passed by reference (see the setter routines). This is in order to
        ensure that different Dsets are independent of each other. Obviously
        a price is paid in terms of speed and memory. If this is ever an issue
        then direct access through the "hidden" versions _data, _errors, _good
        might help along with the null constructor from above, but this is
        experts-only territory.
        """

        if len(args) == 0:
            self.label   = ''
            self.units   = ''
            self._data   = None
            self._errors = None
            self._good   = None

        if (len(args)  == 3 or len(args) == 4 or len(args) == 5) and \
                isinstance(args[0], basestring):

            # starts with a label and units
            self.label = args[0]
            self.units = args[1]
            self.data = args[2]

            if len(args) == 4:
                self.errors = args[3]
                self.good = None
            elif len(args) == 5:
                self.errors = args[3]
                self.good  = copy.copy(args[4])
            else:
                self.errors = None
                self.good   = None

        elif (len(args)  == 2 or len(args) == 3) and isinstance(args[0], Axis):
            # starts from another Axis
            self.label = args[0].label
            self.units = args[0].units
            self.data = args[1]

            if len(args) == 3:
                self.errors = args[2]
                self.good = None
            elif len(args) == 4:
                self.errors = args[2]
                self.good = args[3]
            else:
                self.errors = None
                self.good   = None

        elif (len(args)  == 1 or len(args) == 2) and \
                isinstance(args[0], np.ndarray):
            # minimal - just data
            self.label = ''
            self.units = ''
            self.data  = args[0]

            if len(args) == 2:
                self.errors = args[1]
                self.good = None
            elif len(args) == 3:
                self.errors = args[1]
                self.good = args[2]
            else:
                self.errors = None
                self.good   = None
        else:
            raise DnlError('argument list has wrong number or type of arguments')

    def check(self):
        """Diagnostic routine to check that the Axis is in a consistent state
        as far as can be determined. Raises exceptions if not. Examples of
        possible problems are data that is not an array"""

        if not isinstance(self._data, np.ndarray) or self._data.ndim != 1:
            raise DnlError('_data is not a 1D numpy array')

        if self._errors is not None or not isinstance(self._errors, np.ndarray) \
                or self._errors.ndim != 1:
            raise DnlError('_errors is not None or a 1D numpy array')

        if self._errors is not None and len(self._errors) != len(self._data):
            raise DnlError('_data and _errors have conflicting lengths: ' + \
                               str(len(self._data)) + ' vs ' + str(len(self._errors)))

        if self._good is not None or not isinstance(self._good, np.ndarray) \
                or self._good.ndim != 1:
            raise DnlError('_good is not None or a 1D numpy array')

        if self._good is not None and len(self._good) != len(self._data):
            raise DnlError('_data and _good have conflicting lengths: ' + \
                               str(len(self._data)) + ' vs ' + str(len(self._good)))

        if not isinstance(self.label, str):
            raise DnlError('label is not a string')

        if not isinstance(self.units, str):
            raise DnlError('units is not a string')

    def __repr__(self):
        """Returns string representation of an Axis"""
        # Note we show the direct value of good here to avoid conversion to an
        # array
        rep = 'Axis(label=' + repr(self.label) + ', units=' + \
            repr(self.units) + ', data=' + repr(self.data) + \
            ', errors=' + repr(self.errors) + ', good=' + \
            repr(self._good) + ')'
        return rep

    @property
    def data(self):
        """Returns the data of the Axis"""
        return self._data

    @data.setter
    def data(self, data):
        """Sets the data of the Axis"""
        if isinstance(data, np.ndarray) and data.ndim == 1:
            self._data = copy.copy(data)
        elif isinstance(data, (list, tuple)):
            self._data = np.asarray(data)
            if self._data.ndim != 1:
                raise DnlError('data converts to a numpy array but is not 1D')
            dtype = self._data.dtype
            if not issubclass(dtype, np.integer) and not issubclass(dtype, np.float):
                raise DnlError('data converts to a numpy array but has an invalid data type')
        else:
            raise DnlError('data is not a 1D np.ndarray')

    @property
    def good_data(self):
        """Returns the OK data of the Axis"""
        if self.has_mask:
            return self._data[self._good]
        else:
            return self._data

    @property
    def errors(self):
        """Returns the errors of the Axis"""
        return self._errors

    @errors.setter
    def errors(self, errors):
        """Sets the errors of the Axis"""
        if errors is None:
            self._errors = None
        elif isinstance(errors, np.ndarray) and errors.ndim == 1:
            if len(errors) == len(self):
                self._errors = copy.copy(errors)
            else:
                raise DnlError('length of errors array [' + str(len(errors)) +
                               '] differs from data [' + str(len(self)) + ']')
        elif isinstance(errors, (list, tuple)):
            self._errors = np.asarray(errors)
            if self._errors.ndim != 1:
                raise DnlError('errors converts to a numpy array but is not 1D')
            dtype = self._errors.dtype
            if not issubclass(dtype, np.integer) and not issubclass(dtype, np.float):
                raise DnlError('errors converts to a numpy array but has an invalid data type')
        else:
            raise DnlError('errors is not a 1D numpy.ndarray')

    @property
    def good_errors(self):
        """Returns the OK errors of the Axis"""
        if self.has_mask:
            return self._errors[self._good]
        else:
            return self._errors

    @property
    def good(self):
        """Returns the full good array, even if internally it is set to None"""
        if self.has_mask:
            return self._good
        else:
            return self._data == self._data

    @good.setter
    def good(self, good):
        """Sets the good of the Axis"""
        if good is None:
            self._good = None
        elif isinstance(good, np.ndarray) and good.ndim == 1:
            if len(good) == len(self):
                self._good = copy.copy(good)
            else:
                raise DnlError('length of good array [' + str(len(good)) +
                               '] differs from data [' + str(len(self)) + ']')
        elif isinstance(good, (list, tuple)):
            self._good = np.asarray(good, dtype=np.bool)
            if self._good.ndim != 1:
                raise DnlError('good converts to a numpy array but is not 1D')
        else:
            raise DnlError('good is not a 1D numpy.ndarray')

    @property
    def bad(self):
        """Returns inverse of full good array, even if internally it
        is set to None"""
        if self.has_mask:
            return ~self._good
        else:
            return self._data != self._data

    def __len__(self):
        """Returns the length of an Axis, including any bad data"""
        return len(self.data)

    @property
    def has_errors(self):
        """Does the Axis have errors or not"""
        return self._errors is not None

    @property
    def has_mask(self):
        """Does the Axis have a good array or not"""
        return self._good is not None

    # basic arithmetical operations

    # negation
    def __neg__(self):
        ok = self.good
        temp  = copy.deepcopy(self)
        temp.data[ok] = -temp.data[ok]
        return temp

    # absolute value
    def __abs__(self):
        ok = self.good
        temp  = copy.deepcopy(self)
        temp.data[ok] = np.abs(temp.data[ok])
        return temp

    # in place operations
    def __iadd__(self, other):
        """+= in-place addition operator. 'other' can be another Axis, an
        np.ndarray or a constant. Self defines the output label and units.
        Errors are assumed independent"""

        if isinstance(other, Axis):
            ok = self.good & other.good
            if self.has_mask or other.has_mask:
                self.good = ok
            if self.has_errors and other.has_errors:
                self.errors[ok] = np.sqrt(self.errors[ok]**2 +
                                          other.errors[ok]**2)
            elif other.has_errors:
                self.errors[ok] = other.errors[ok]
            self.data[ok] += other.data[ok]
        else:
            ok = self.good
            if isinstance(other, np.ndarray):
                self.data[ok] += other[ok]
            else:
                self.data[ok] += other

        return self

    def __isub__(self, other):
        """-= in-place subtraction operator. 'other' can be another Axis, an
        np.ndarray or a constant. Self defines the output label and units.
        Errors are assumed independent"""

        if isinstance(other, Axis):
            ok = self.good & other.good
            if self.has_mask or other.has_mask:
                self.good = ok
            if self.has_errors and other.has_errors:
                self.errors[ok] = np.sqrt(self.errors[ok]**2 +
                                          other.errors[ok]**2)
            elif other.has_errors:
                self.errors[ok] = other.errors[ok]
            self.data[ok] -= other.data
        else:
            ok = self.good
            if isinstance(other, np.ndarray):
                self.data[ok] -= other[ok]
            else:
                self.data[ok] -= other

        return self

    def __imul__(self, other):
        """*= in-place addition operator. 'other' can be another Axis, an
        np.ndarray or a constant. Self defines the output label and units.
        Errors are assumed independent and combined through linear approx."""

        if isinstance(other, Axis):
            ok = self.good & other.good
            if self.has_mask or other.has_mask:
                self.good = ok
            if self.has_errors and other.has_errors:
                self.errors[ok] = np.sqrt((other.data[ok]*self.errors[ok])**2
                                          + (self.data[ok]*other.errors[ok])**2)
            elif other.has_errors:
                self.errors[ok] = np.abs(self.data[ok])*other.errors[ok]
            else:
                self.errors[ok] *= np.abs(other.data[ok])
            self.data[ok] *= other.data[ok]
        else:
            ok = self.good
            if isinstance(other, np.ndarray):
                self.data[ok] *= other[ok]
                if self.has_errors:
                    self.errors[ok] *= np.abs(other.data[ok])
            else:
                self.data[ok] *= other
                if self.has_errors:
                    self.errors[ok] *= np.abs(other)
        return self

    def __itruediv__(self, other):
        """/= in-place division operator. 'other' can be another Axis, an
        np.ndarray or a constant. Self defines the output label and units.
        Errors are assumed independent"""

        if isinstance(other, Axis):
            ok = self.good & other.good
            if self.has_mask or other.has_mask:
                self.good = ok
            if self.has_errors and other.has_errors:
                self.errors[ok] = np.sqrt(
                    (other.data[ok]*self.errors[ok])**2
                    + (self.data[ok]*other.errors[ok])**2)/other.data[ok]**2
            elif other.has_errors:
                self.errors[ok] = np.abs(self.data[ok])*\
                                  other.errors[ok]/other.data[ok]**2
            else:
                self.errors[ok] /= np.abs(other.data[ok])
            self.data[ok] /= other.data[ok]

        else:
            ok = self.good
            if isinstance(other, np.ndarray):
                self.data[ok] /= other[ok]
                if self.has_errors:
                    self.errors[ok] /= np.abs(other.data[ok])
            else:
                self.data[ok] /= other
                if self.has_errors:
                    self.errors[ok] /= np.abs(other)

        return self

    # temporary for python2.x
    __idiv__ = __itruediv__

    # (self op other) operations
    def __add__(self, other):
        """+ addition operator. 'other' can be another Axis, an np.ndarray
        or a constant. Errors are assumed independent."""
        temp  = copy.deepcopy(self)
        temp += other
        return temp

    def __sub__(self, other):
        """- subtraction operator. 'other' can be another Axis, an np.ndarray
        or a constant. Errors are assumed independent."""
        temp  = copy.deepcopy(self)
        temp -= other
        return temp

    def __mul__(self, other):
        """* multiplication operator. 'other' can be another Axis, an
        np.ndarray or a constant. Errors are assumed independent.
        """
        temp  = copy.deepcopy(self)
        temp *= other
        return temp

    def __truediv__(self, other):
        """/ division operator. 'other' can be another Axis, an np.ndarray
        or a constant. Errors are assumed independent."""
        temp  = copy.deepcopy(self)
        temp /= other
        return temp

    __div__ = __truediv__

    # (other op self) operations
    __radd__ = __add__

    def __rsub__(self, other):
        """- subtraction operator. 'other' can be an np.ndarray
        or a constant. Errors are assumed independent."""
        temp  = -copy.deepcopy(self)
        temp += other
        return temp

    __rmul__ = __mul__

    def __rdiv__(self, other):
        """/ division operator of form 'other / self'. 'other' can
        be an np.ndarray or a constant."""
        temp  = copy.deepcopy(self)
        ok    = self.good
        if isinstance(other, np.ndarray):
            if temp.has_errors:
                temp.errors[ok] = np.abs(other[ok])*\
                                  temp.errors[ok]/temp.data[ok]**2
            temp.data[ok] = other[ok] / temp.data[ok]
        else:
            if temp.has_errors:
                temp.errors[ok] = np.abs(other)*\
                                  temp.errors[ok]/temp.data[ok]**2
            temp.data[ok] = other / temp.data[ok]
        return temp

    # some slightly less common operators
    def log10(self):
        """Takes log10 of the Axis; negative points are masked"""
        ok = self.good & self.data > 0.
        if self.has_errors:
            self.errors[ok]  = self.errors[ok]/self.data[ok]/m.log(10.)
        self.data[ok] = np.log10(self.data[ok])
        if len(self.data[ok]) < len(self.data):
            self.good = ok

    def alog10(self):
        """Takes anti-log10 of the Axis"""
        ok = self.good
        self.data[ok] = 10.**self.data[ok]
        if self.has_errors:
            self.errors[ok] = self.data[ok]*(10.**self.errors[ok]-1.)
        else:
            self.data = 10.**self.data
            if self.has_errors:
                self.errors = self.data*(10.**self.errors-1.)

    # other routines
    def plimits(self, good=None, off=0.):
        """
        Computes plot limits for an Axis

        good -- boolean mask; only these points will be included, along with
                any flagged good internally. Must have full length of data
        off  -- offset to subtract.

        Returns (plo,phi) lower and upper limits adjusted to include errors.
        These will be None if there were no points
        """
        if good is None:
            ok = self.good
        else:
            ok = self.good & good

        if len(self.data[ok]) > 0:
            if self.errors is None:
                plo = self.data[ok].min() - off
                phi = self.data[ok].max() - off
            else:
                plo = (self.data[ok] - self.errors[ok]).min() - off
                phi = (self.data[ok] + self.errors[ok]).max() - off
            extra = (phi - plo)/40.
            plo  -= extra
            phi  += extra
            return (plo,phi)
        else:
            return (None,None)

    def min(self):
        """
        Returns the minimum, ignoring masked points.
        """
        if self._good is None:
            return self.data.min()
        else:
            return self.data[self.good].min()

    def max(self):
        """
        Returns the maximum, ignoring masked points.
        """
        if self._good is None:
            return self.data.max()
        else:
            return self.data[self.good].max()

    def append(self, axis):
        """
        Appends an Axis to this Axis. Each must
        either have or not have errors for this
        to work.
        """
        # both inputs must have errors or neither should
        # have them for this to work
        if self.has_errors and axis.has_errors:
            self._errors = np.concatenate([self._errors,axis._errors])
        elif (self.has_errors and not axis.has_errors) or \
             (not self.has_errors and axis.has_errors):
            raise DnlError(
                'Axis.append: inputs must both have, or both not have, errors')

        self._data = np.concatenate([self._data, axis._data])

        if self.has_mask and axis.has_mask:
            self._good = np.concatenate([self._good, axis._good])
        elif self.has_mask and not axis.has_mask:
            self._good = np.concatenate([self._good,
                                         np.ones(len(axis._data),dtype=bool)])
        elif not self.has_mask and axis.has_mask:
            self._good = np.concatenate([np.ones(len(self._data),dtype=bool),
                                         axis._good])


# Class to represent 1D data

class Dset(object):
    """Represents data with a header, x and y axes and a logical mask.

    This is the key object to represent a set of 1D data of x and y values
    with potentially errors in both sets of values, an array to say which data
    are good. The attributes are::

       x : Axis
           The X Axis

       y : Axis
           The Y Axis

       head : astropy.io.fits.Header
           The header.

       good : 1D numpy.ndarray or None
           mask array saying which data are good, in addition to the internal
           good arrays of the x and y Axis objects. This is designed to be a
           more temporary mask than the x and y mask arrays. The latter flag up
           data that should be ignored in any operation. Thus when the docs
           talk about masked or unmasked points, they mean which points are
           masked or unmasked once those already masked in the x & y arrays are
           removed.

       ptype : int
          default plot type (POINTS, LINE, BAR)

     A number of header keywords are set and used at various points and
     should either not be used, or used with care. Their names and meanings
     are as follows:

     Times (which should all refer to mid-exposure when possible)::

       DATE-OBS : string
            Time. Usual form YYYY-MM-DDTHH:MM:SS.SSS..
            Should refer to mid-exposure where relevant.

       MJD : float / int
            Modified Julian Day (UTC)
            Should refer to mid-exposure where relevant.

       UTC : float
            Decimal hours.
            Should refer to mid-exposure where relevant.

     Target information::

       OBJECT : string
            Astronomical object name

       RA : string / float
            Right Ascension. If string HH:MM:SS.SS, if float decimal degrees.

       DEC : string / float
            Declination. If string [+-]DD:SS.SS, if float decimal degrees.

       SYSTEM : string
            Coordinate system e.g. 'B1950', 'J2000', 'ICRS'

       EPOCH :  float
            Julian epoch of coordinates

       PMRA : float
            Proper motion in RA, arcsec / year

       PARALLAX : float
            Parallax, arcseconds.

       RV : float
            Radial velocity, km/s.

     Observatory information::

        TELESCOP : string
             name of telescope

        OBSERVAT : string
             name of observatory

        LONGITUD : float
             geodetic longitude, degrees

        LATITUDE : float
             geodetic latitude, degrees

        HEIGHT : float
             height in metres

      Instrument::

        INSTRUME : string
             instrument name

        NCCD : int
             CCD id number

        NAPERTUR : int
             Aperture id number

        FILTER : string
             Filter used during observation

      Others::

        TYPE : string
             Type of data, e.g. 'time series', 'spectrum'
    """

    # possible plot types
    POINTS = 0
    LINE   = 1
    BAR    = 2

    def __init__(self,  x, y, head=None, good=None, ptype=POINTS):
        """
        Creates a new Dset.

        Arguments::

        x : Axis or ndarray
            The x-axis data and errors

        y : Axis or ndarray
            The y-axis data and errors

        head : astropy.fits.Header
               Ordered dictionary of header information

        good : 1D boolean ndarray or None
               Mask array to say which data are OK.

        ptype : [POINTS, LINE, BAR]
                default type to use when plotting.
        """

        self.x = x
        self.y = y
        self.head = head
        self.good = good
        self.ptype = ptype

    def check(self):
        """Diagnostic routine to check that a Dset is in a consistent state."""

        self._x.check()
        self._y.check()

        if len(self._x) != len(self._y):
            raise DnlError('X and Y axes have differing lengths: ' + \
                           str(len(self._x)) + ' vs ' + str(len(self._y)))

        if self._good is not None and not isinstance(self._good, np.ndarray) \
           and self._good.ndim != 1:
            raise DnlError('_good is not None and is not a 1D numpy array')

        if self._good is not None and len(self._good) != len(self._x):
            raise DnlError('_good array does not match length of Axes: ' + \
                           str(len(self._good)) + ' vs ' + str(len(self._x)))

    def __len__(self):
        """Returns the length of a Dset"""
        return len(self.x)

    def __repr__(self):
        """Returns string representation of a Dset"""
        rep = 'Dset(x=' + repr(self.x) + ', y=' + repr(self.y) + ', head=' + \
            repr(self.head) + ', good=' + repr(self.good) + \
            ', ptype=' + repr(self.ptype) + ')'
        return rep

    __str__ = __repr__

    @property
    def x(self):
        """Gets the X axis"""
        return self._x

    @x.setter
    def x(self, x):
        """Sets the X axis"""
        if isinstance(x, Axis):
            self._x = copy.deepcopy(x)
        elif isinstance(x, np.ndarray):
            self._x = Axis('X','',x)
        else:
            raise DnlError('dnl.Dset: x is not a dnl.Axis or a numpy.ndarray')

    @property
    def y(self):
        """Gets the Y axis"""
        return self._y

    @y.setter
    def y(self, y):
        """Sets the Y axis"""
        if isinstance(y, Axis):
            if len(y) == len(self):
                self._y = copy.deepcopy(y)
            else:
                raise DnlError('dnl.Dset.y: y length [' + str(len(y)) + \
                               ' does not match x length [' + \
                               str(len(self)) + ']')

        elif isinstance(y, np.ndarray):
            if len(y) == len(self):
                self._y = Axis('Y','',y)
            else:
                raise DnlError('dnl.Dset.y: y length [' + str(len(y)) + \
                               ' does not match x length [' + \
                               str(len(self)) + ']')
        else:
            raise DnlError('dnl.Dset.y: y is not a dnl.Axis' +
                            ' or a numpy.ndarray')

    @property
    def head(self):
        """Gets the header"""
        return self._head

    @head.setter
    def head(self, head):
        """Sets the header"""
        if isinstance(head, fits.Header):
            self._head = copy.copy(head)
        elif head is None:
            self._head = fits.Header()
        else:
            raise DnlError('head must be None or an astropy.io.fits.Header')

    @property
    def good(self):
        """Returns mask of good pixels, with x and y masks applied"""
        if self.has_mask:
            return self._good & self.x.good & self.y.good
        else:
            return self.x.good & self.y.good

    @good.setter
    def good(self, good):
        """Sets the good data"""

        if good is None:
            self._good = None
        elif isinstance(good, np.ndarray) and good.ndim == 1:
            if len(good) == len(self):
                self._good = copy.copy(good)
            else:
                raise DnlError('length of good array [' +
                               str(len(good)) + '] differs from data [' +
                               str(len(self)) + ']')
        elif isinstance(good, (list, tuple)):
            self._good = np.asarray(good, dtype=np.bool)
            if self._good.ndim != 1:
                raise DnlError('good converts to a numpy array but is not 1D')
        else:
            raise DnlError('good is not a 1D numpy.ndarray')

    @property
    def ptype(self):
        """Returns plot type"""
        return self._ptype

    @ptype.setter
    def ptype(self, ptype):
        """Sets the point type"""
        if ptype != Dset.POINTS and ptype != Dset.LINE and ptype != Dset.BAR:
            raise DnlError('ptype invalid: must be one of POINTS, LINE or BAR')
        self._ptype  = ptype

    @property
    def title(self):
        return self.head['TITLE'] if 'TITLE' in self.head else ''

    @property
    def has_mask(self):
        return self._good is not None

    @property
    def bad(self):
        """Returns mask of bad pixels, with x and y masks applied as normal"""

        if self.has_mask:
            return ~self._good | self.x.bad | self.y.bad
        else:
            return self.x.bad | self.y.bad

    @property
    def ngood(self):
        """Returns number of unmasked points (<= value returned by "len")"""
        return len(self.x.data[self.good])

    @property
    def nbad(self):
        """Returns number of masked points (<= value returned by "len")"""
        return len(self.x.data[self.bad])

    def rmask(self):
        """Inverts the mask"""
        if self._good is not None:
            self._good = ~self._good

    @property
    def mean(self):
        """Returns mean y data value of unmasked pixels"""
        return np.mean(self.y.data[self.good])

    # in place operators
    def __iadd__(self, other):
        """+= in-place addition operator. 'other' can be another Dset, an
        Axis, an array or a constant. Addition here is assumed
        to apply to the y array, the x array is untouched."""

        if isinstance(other, Dset):
            if self.has_mask or other.has_mask:
                self.good = self.good & other.good
            self.y += other.y
        elif isinstance(other, Axis):
            if self.has_mask or other.has_mask:
                self.good = self.good & other.good
            self.y += other
        else:
            self.y += other

        return self

    def __isub__(self, other):
        """-= in-place subtraction operator. 'other' can be another Dset, an
        Axis, an array or a constant. Addition here is assumed
        to apply to the y array, the x array is untouched."""

        if isinstance(other, Dset):
            if self.has_mask or other.has_mask:
                self.good = self.good & other.good
            self.y -= other.y
        elif isinstance(other, Axis):
            if self.has_mask or other.has_mask:
                self.good = self.good & other.good
            self.y -= other
        else:
            self.y -= other

        return self

    def __imul__(self, other):
        """*= in-place multiplication operator. 'other' can be another Dset,
        an Axis, an array or a constant"""

        if isinstance(other, Dset):
            if self.has_mask or other.has_mask:
                self.good = self.good & other.good
            self.y *= other.y
        elif isinstance(other, Axis):
            if self.has_mask or other.has_mask:
                self.good = self.good & other.good
            self.y *= other
        else:
            self.y *= other

        return self

    def __itruediv__(self, other):
        """/= in-place division operator. 'other' can be another Dset, an Axis,
        an array or a constant"""
        if isinstance(other, Dset):
            if self.has_mask or other.has_mask:
                self.good = self.good & other.good
            self.y /= other.y
        elif isinstance(other, Axis):
            if self.has_mask or other.has_mask:
                self.good = self.good & other.good
            self.y /= other
        else:
            self.y /= other
        return self

    __idiv__ = __itruediv__

    # standard direction
    def __add__(self, other):
        """+ addition operator. 'other' can be another Dset,
        an Axis, an array or a constant"""
        temp  = copy.deepcopy(self)
        temp += other
        return temp

    def __sub__(self, other):
        """- subtraction operator. 'other' can be another Dset,
        an Axis, an array or a constant"""
        temp  = copy.deepcopy(self)
        temp -= other
        return temp

    def __mul__(self, other):
        """* multiplication operator. 'other' can be another Dset,
        an Axis, an array or a constant"""
        temp  = copy.deepcopy(self)
        temp *= other
        return temp

    def __truediv__(self, other):
        """/ division operator. 'other' can be another Dset, an Axis,
        an array or a constant"""
        temp  = copy.deepcopy(self)
        temp /= other
        return temp

    __div__ = __truediv__

    # reversed
    __radd__ = __add__

    def __rsub__(self, other):
        """- subtraction operator. 'other' can be an Axis, an array
        or a constant"""
        temp   = copy.deepcopy(self)
        temp.y = other - temp.y
        return temp

    __rmul__ = __mul__

    def __rdiv__(self, other):
        """/ division operator. 'other' can be an Axis, an array
        or a constant"""
        temp  = copy.deepcopy(self)
        temp.y = other / temp.y
        return temp

    def bin(self, nx, x1=None, x2=None, weight='u', yerrors='i',
            yemin=0.01, xerrors='i'):
        """
        Bins a Dset into regularly spaced bins.

        Empty bins are not returned; masked input is ignored.

        Arguments::

          nx : int
               number of new pixels.

          x1 : float
               left-hand edge of first new pixel, 'None' for minimum of data

          x2 : float
               right-hand edge of last new pixel, 'None' for maximum of data

          weight : string
                   'u'=uniform, 'v'=inverse y-variance

          yerrors : string
                   'i'=calculate output y errors from weighting the input
                    errors (if any), 'v'=calculate output y errors from the
                    error-on-the-mean given the variance within each bin
                    (requires at least 2 points per bin, and obviously better
                    with more).

          yemin : float
                  lower limit to output errors as multiple of mean error
                  (only applies in y) to prevent overweighting

        Exceptions are raised if there are no output bins.
        """

        ok = self.good

        # make array of weights
        if weight == 'u' or self.y.has_errors:
            weights = np.ones(self.x.data[ok].shape)
        elif weight == 'v':
            weights = 1./np.square(self.y.errors[ok])
        else:
            raise DnlError('weight = ' + weight + ' not recognised.')

        # define minimum number of points required per bin
        if yerrors == 'i' and self.y.has_errors:
            nmin = 1
        elif yerrors == 'v':
            nmin = 2
        else:
            raise DnlError('yerrors = ' + yerrors + ' not recognised.')

        # set start and end values
        if x1 is None:
            x1 = np.min(self.x.data[ok])
        if x2 is None:
            x2 = np.max(self.x.data[ok])

        # create index array
        index = (nx*(self.x.data[ok] - x1)/(x2-x1)).astype(int)

        # count number of points/bin
        nbin = np.bincount(index)

        # only carry on if there are enough bins with enough points
        enough = nbin >= nmin
        if len(nbin[enough]) == 0:
            raise DnlError('dnl.Dset.bin: no points in output binned array')

        # bin the data
        wbin = np.bincount(index, weights)
        xbin = np.bincount(index, weights*self.x.data[ok])
        ybin = np.bincount(index, weights*self.y.data[ok])

        xbin[enough] /= wbin[enough]
        ybin[enough] /= wbin[enough]

        # compute the X errors
        if xerrors == 'i' and self.x.has_errors:
            xebin      = np.bincount(
                index,np.square(weights*self.x.errors[ok]))
            xebin[enough] /= np.square(wbin[enough])
            xebin[enough]  = np.sqrt(xebin[enough])
        elif xerrors == 'v':
            xebin      = np.bincount(
                index,weights*np.square(self.x.data[ok]-xbin[index]))
            xebin[enough] /= wbin[enough]
            xebin[enough]  = np.sqrt(xebin[enough]/nbin[enough])
        else:
            xebin = None

        # compute the Y errors
        if yerrors == 'i' and self.y.has_errors:
            yebin      = np.bincount(
                index,np.square(weights*self.y.errors[ok]))
            yebin[enough] /= np.square(wbin[enough])
            yebin[enough]  = np.sqrt(yebin[enough])
        elif yerrors == 'v':
            yebin      = np.bincount(
                index,weights*np.square(self.y.data[ok]-ybin[index]))
            yebin[enough] /= wbin[enough]
            yebin[enough]  = np.sqrt(yebin[enough]/nbin[enough])
        else:
            yebin = None

        if yebin != None:
            mine = yemin*yebin[enough].mean()
            maxe = 1.1*yebin[enough].max()
            yebin[enough].clip(min=mine,max=maxe)

        newx = Axis(self.x, xbin[enough], xebin[enough]
                    if xebin is not None else None)
        newy = Axis(self.y, ybin[enough], yebin[enough])

        # return the result
        return Dset(newx, newy, self.head)

    def oneLine(self):
        """Prints a one-line summary of a Dset"""
        strg = ''
        if 'OBJECT' in self.head:
            strg += self.head['OBJECT']
        if 'DATE-OBS' in self.head:
            date = self.head['DATE-OBS'].replace('T',' ')
            if strng != '': strng += ', '
            strg += 'time: ' + date

        xmid = (self.x.min()+self.x.max())/2.
        xrng = self.x.max()-self.x.min()
        if strg != '': strg += ', '
        strg += 'x (mid,rng): {0:f} {1:f} [{2:s}]'.format(
            xmid,xrng,self.x.units)

        ymid = (self.y.min()+self.y.max())/2.
        yrng = self.y.max()-self.y.min()
        strg += ', y (mid,rng): {0:f} {1:f} [{2:s}]'.format(
            ymid,yrng,self.y.units)

        strg += ', npts: {0:d}'.format(len(self))
        return strg

    def toHDU(self):
        """Returns a Dset as an astropy.io.fits.BinTableHDU object suitable for
        writing to a FITS file."""

        TFORM = {'int32' : 'J', 'int64' : 'K', 'float32' : 'E',
                 'float64' : 'D', 'bool' : 'L'}

        cols = []
        cols.append(fits.Column(name='X data',
                                format=TFORM[self.x.data.dtype.name],
                                unit=self.x.units, array=self.x.data))

        if self.x.has_errors:
            cols.append(fits.Column(name='X errors',
                                    format=TFORM[self.x.data.dtype.name],
                                    unit=self.x.units, array=self.x.errors))

        if self.x.has_mask:
            cols.append(fits.Column(name='X good',
                                    format=TFORM[self.x.good.dtype.name],
                                    array=self.x.good))

        cols.append(fits.Column(name='Y data',
                                format=TFORM[self.y.data.dtype.name],
                                unit=self.y.units, array=self.y.data))

        if self.y.has_errors:
            cols.append(fits.Column(name='Y errors',
                                    format=TFORM[self.y.data.dtype.name],
                                    unit=self.y.units, array=self.y.errors))
        if self.y.has_mask:
            cols.append(fits.Column(name='Y good',
                                    format=TFORM[self.y.good.dtype.name],
                                    array=self.y.good))

        if self.has_mask:
            cols.append(fits.Column(name='Good',
                                    format=TFORM[self.good.dtype.name],
                                    array=self.good))

        coldefs = fits.ColDefs(cols)
        data = fits.FITS_rec.from_columns(coldefs)

        htemp = self.head.copy()
        htemp['XLABEL'] = self.x.label
        htemp['YLABEL'] = self.y.label

        return fits.BinTableHDU(data,htemp)

    @classmethod
    def fromHDU(cls, hdu):
        """Returns a new Dset given a correctly formatted HDU (as
        written by toHDU)

          hdu : astropy.io.fits.BinTableHDU
             contains X data, Y data
        """
        head    = hdu.header
        table   = hdu.data
        cnames  = table.dtype.names

        xlabel  = head['XLABEL']
        xdata   = table['X data']
        xerrors = table['X errors'] if 'X errors' in cnames else None
        xgood   = table['X good'] if 'X good' in cnames else None
        xdind   = cnames.index('X data') + 1
        xunits  = head['TUNIT' + str(xdind)]
        xaxis   = Axis(xlabel, xunits, xdata, xerrors, xgood)

        ylabel  = head['YLABEL']
        ydata   = table['Y data']
        yerrors = table['Y errors'] if 'Y errors' in cnames else None
        ygood   = table['Y good'] if 'Y good' in cnames else None
        ydind   = cnames.index('Y data') + 1
        yunits  = head['TUNIT' + str(ydind)]
        yaxis   = Axis(ylabel, yunits, ydata, yerrors, ygood)

        good    = table['Good'] if 'Good' in table else None
        return Dset(xaxis, yaxis, head, good, Dset.POINTS)

    def plimits(self, xoff=0., yoff=0.):
        """
        Computes plot limits.

          xoff : float
                 offset to subtract from the X axis
          yoff    -- offset to subtract from the Y axis

        Bad data are ignored.
        """
        ok = self.good
        xlo, xhi = self.x.plimits(ok, xoff)
        ylo, yhi = self.y.plimits(ok, yoff)
        return (xlo,xhi,ylo,yhi)

    def ppoints(self, psymb=17, xoff=0., yoff=0., pbad=False):
        """
        Plots data only as points.

        Arguments::

          psymb : plot symbol. In general this can be a string which
                  will be fed to plstring. This can include escape sequences
                  as detailed in plplot docs. Anything but a string will be
                  sent to plpoin and should therefore be an integer in the
                  range -1 to 127 inclusive to avoid errors
          xoff  : offset in X to subtract
          yoff  : offset in Y to subtract
          pbad  : plot the bad (as opposed to good) data
        """
        nplot = self.nbad if pbad else self.ngood

        if nplot:
            # only plot anything if there are points to plot
            ok = self.bad if pbad else self.good
            if isinstance(psymb, str):
                plstring(self.x.data[ok]-xoff,
                         self.y.data[ok]-yoff, psymb)
            else:
                plpoin(self.x.data[ok]-xoff,
                       self.y.data[ok]-yoff, psymb)

    def pline(self, xoff=0., yoff=0., pbad=False):
        """
        Plots data as line connecting points. There will be gaps
        where points are masked, including those masked within the
        x and y arrays.

        Arguments::

          xoff  : float
                offset in X to subtract
          yoff  : float
                offset in Y to subtract
          pbad  : bool
                plot the bad (as opposed to good) data
        """

        ok = self.bad if pbad else self.good

        for start, end in chunks(ok):
            plline(self.x.data[start:end]-xoff, self.y.data[start:end]-yoff)

    def pbin(self, xoff=0., yoff=0., opt=plg.PL_BIN_CENTRED, pbad=False):
        """
        Plots data as bar chart between points.

          xoff : float
               offset in X to subtract
          yoff : float
               offset in Y to subtract
          opt : int
                 controls positioning of bins
          pbad : bool
               controls positioning of bins
        """

        ok = self.bad if pbad else self.good

        for start, end in chunks(ok):
            plbin(self.x.data[start:end]-xoff, self.y.data[start:end]-yoff,opt)


    def pyebin(self, xoff=0., yoff=0., opt=plg.PL_BIN_CENTRED, pbad=False):
        """
        Plots errors in Y as a bar chart.

          xoff : float
               offset in X to subtract
          yoff : float
               offset in Y to subtract (i.e. not something one
               would usually use)
          opt : int
                 controls positioning of bins
          pbad : bool
               controls positioning of bins
        """

        if self.y.has_errors:
            ok = self.bad if pbad else self.good

            for start, end in chunks(ok):
                plbin(self.x.data[start:end]-xoff,
                      self.y.errors[start:end]-yoff,opt)

    def pxerrors(self, xoff=0., yoff=0., pbad=False, terminal=0.):
        """
        Plots errors bars in X in usual horizontal line style

          xoff : float
               offset in X to subtract
          yoff : float
               offset in Y to subtract (i.e. not something one
               would usually use)
          pbad : bool
               controls positioning of bins
          terminal : float
               length of terminals on error bars
        """
        if self.x.has_errors and ((pbad and self.nbad > 0) or \
                                  (not pbad and self.ngood > 0)):
            ok = self.bad if pbad else self.good
            pgerrx(self.x.data[ok]-xoff-self.x.errors[ok],
                   self.x.data[ok]-xoff+self.x.errors[ok],
                   self.y.data[ok]-yoff, terminal)

    def pyerrors(self, xoff=0., yoff=0., pbad=False, terminal=0.):
        """
        Plots errors bars in X in usual vertical line style

          xoff : float
               offset in X to subtract
          yoff : float
               offset in Y to subtract (i.e. not something one would
               usually use)
          pbad : bool
               controls positioning of bins
          terminal : float
               length of terminals on error bars
        """
        if self.y.has_errors and ((pbad and self.nbad > 0) or \
                                  (not pbad and self.ngood > 0)):
            ok = self.bad if pbad else self.good
            pgerry(self.x.data[ok]-xoff, self.y.data[ok]-yoff-self.y.errors[ok],
                   self.y.data[ok]-yoff+self.y.errors[ok], terminal)

    def plabels(self, xoff=0., yoff=0.):
        """Returns default plot labels (xlabel,ylabel,tlabel)"""

        if xoff == 0.:
            xlabel = self.x.label
        elif xoff < 0.:
            xlabel = self.x.label + ' + ' + str(-xoff)
        else:
            xlabel = self.x.label + ' - ' + str(xoff)
        if self.x.units != '':
            xlabel += ' (' + self.x.units + ')'

        if yoff == 0.:
            ylabel = self.y.label
        elif yoff < 0.:
            ylabel = self.y.label + ' + ' + str(-yoff)
        else:
            ylabel = self.y.label + ' - ' + str(yoff)
        if self.y.units != '':
            ylabel += ' (' + self.y.units + ')'

        return (xlabel,ylabel,self.title)

    def plot(self, xoff=0., yoff=0., psymb=17, dcol=1, ecol=2,
             mdcol=5, mecol=5, perr=True, pbad=True):
        """
        Wraps up a bundle of methods to plot a Dset. You need to have
        opened the plot and defined the scales before using this method.

        xoff  -- offset to subtract from X data, e.g. to make readable
        yoff  -- offset to subtract from Y data
        psymb -- plot symbol for points
        dcol  -- colour for data, colour index 1=white, 2=red, etc
        ecol  -- colour for errors, "     "
        mdcol -- colour for masked data, "    "
        mecol -- colour for masked errors, "    "
        perr  -- plot errors or not
        pbad  -- plot masked data or not
        """

        if self.ptype == Dset.POINTS:
            if pbad:
                if perr:
                    plcol0(mecol)
                    self.pxerrors(xoff,yoff,True)
                    self.pyerrors(xoff,yoff,True)
                plcol0(mdcol)
                self.ppoints(psymb,xoff,yoff,True)

            if perr:
                plcol0(ecol)
                self.pxerrors(xoff,yoff)
                self.pyerrors(xoff,yoff)

            plcol0(dcol)
            self.ppoints(psymb,xoff,yoff)

        elif self.ptype == Dset.LINE:
            if pbad:
                if perr:
                    plcol0(mecol)
                    self.pxerrors(xoff,yoff,True)
                    self.pyerrors(xoff,yoff,True)

                plcol0(mdcol)
                self.pline(xoff,yoff,True)

            if perr:
                plcol0(ecol)
                self.pyerrors(xoff,yoff,False)

            plcol0(dcol)
            self.pline(xoff,yoff,False)

        elif self.ptype == Dset.BAR:
            if masked:
                if errors:
                    plcol0(mecol)
                    self.pyebin(xoff,yoff,True)

                plcol0(mdcol)
                self.pbin(xoff,yoff,True)

            if errors:
                plcol0(ecol)
                self.pyebin(xoff,yoff,False)

            plcol0(dcol)
            self.pbin(xoff,yoff,False)

        else:
            raise DnlError('Dset.plot: unrecognised plot ptype.' +
                           ' Should not have happened')


    def setpos(self, object, ra, dec, system='ICRS', pmra=None, pmdec=None,
               epoch=None, parallax=None, rv=None):
        """
        Set astronomical target data in Dset headers

        Arguments::

          object : string
              target name

          ra : string / float
              Right Ascension. If a string it should be in hh:mm:ss.ss..
              format, if a float it should be in decimal degrees.

          dec : string / float
              Declination. If a string it should be in [-+]dd:mm:ss.ss..
              format, if a float it should be in decimal degrees.

          system : string
              Coordinate system, 'ICRS', 'J2000' or 'B1950'

          pmra : float
              proper motion in RA, arcseconds per year (not seconds of RA)

          pmdec : float
              proper motion in Dec, arcseconds per year

          epoch : float
              Julian epoch of coordinates, e.g. 2000.

          parallax : float
              parallax in arcsec

          rv : float
              radial velocity km/s

        Those which are default 'None' will not be set in the headers unless
        explicitly given.
        """

        # run checks before doing anything
        SYSTEMS = ('ICRS', 'J2000', 'B1950')
        if system not in SYSTEMS:
            raise DnlError('system = ' + str(system) +
                           ' not recognised amongst possibles: ' +
                           str(SYSTEMS))

        same = (pmra is None and pmdec is None and epoch is None) or \
            (pmra is not None and pmdec is not None and epoch is not None)
        if not same:
            raise DnlError('Dset.setpos: pmra, pmdec, and epoch must all be set or not at all')

        self.head['OBJECT'] = (object,'Target name')
        if isinstance(ra, str):
            self.head['RA'] = (ra, 'Right Ascension, hh:mm:ss.ss..')
        else:
            self.head['RA'] = (ra, 'Right Ascension, decimal degrees')
        if isinstance(dec, str):
            self.head['DEC'] = (dec, 'Declination, [-+]dd:mm:ss.ss..')
        else:
            self.head['DEC'] = (dec, 'Declinatio, decimal degrees')
        self.head['SYSTEM'] = (system, 'Coordinate system')

        if epoch is not None:
            self.head['EPOCH'] = (epoch, 'Julian epoch of coordinates')
            self.head['PMRA']  = (pmra, 'Proper motion in RA, arcsec/year')
            self.head['PMDEC'] = (pmdec, 'Proper motion in Dec, arcsec/year')

        if parallax is not None:
            self.head['PARALLAX'] = (parallax,'Parallax to object in arcsec')

        if rv is not None:
            self.head['RV'] = (rv, 'Radial velocity, km/s')

    def getpos(self):
        """
        Get astronomical target data from the headers, if set. Raises
        an Exception if nothing is set, or if any of the keywords 'RA',
        'DEC' or 'SYSTEM' is not set. Other items are optional.

        Returns: (object, ra, dec, system, pmra, pmdec, epoch, parallax, rv)
        where::

          object : string
              target name. Defaults to None if not set.

          ra : float
              Right Ascension, decimal hours. Raises an error if not set.

          dec : float
              Declination, decimal degrees. Raises an error if not set.

          system : string
              Coordinate system, 'ICRS', 'J2000' or 'B1950'. Raises an
              error if not set.

          pmra : float
              proper motion in RA, arcseconds per year (not seconds of RA)
              defaults to 0 if not set.

          pmdec : float
              proper motion in Dec, arcseconds per year
              defaults to 0 if not set.

          epoch : float
              Julian epoch of coordinates
              defaults to 2000. if not set.

          parallax : float
              parallax in arcsec
              defaults to 0 if not set.

          rv : float
              radial velocity km/s
              defaults to 0 if not set.
        """

        if 'RA' not in self.head or 'DEC' not in self.head or \
                'SYSTEM' not in self.head:
            raise DnlError('Dset.getpos: one or more of RA, DEC, SYSTEM not set')

        if 'OBJECT' in self.head:
            object = self.head['OBJECT']
        else:
            object = None

        ra = self.head['RA']
        if isinstance(ra, str): ra = subs.hms2d(ra)
        dec = self.head['DEC']
        if isinstance(dec, str): dec = subs.hms2d(dec)
        system = self.head['SYSTEM']

        epoch = self.head['EPOCH'] if 'EPOCH' in self.head else 2000.
        pmra  = self.head['PMRA'] if 'PMRA' in self.head else 0.
        pmdec  = self.head['PMDEC'] if 'PMDEC' in self.head else 0.
        parallax = self.head['PARALLAX'] if 'PARALLAX' in self.head else 0.
        rv = self.head['RV'] if 'RV' in self.head else 0.

        return (object, ra, dec, system, pmra, pmdec, epoch, parallax, rv)


    def settel(self, telescope, observatory, longitude, latitude, height):
        """
        Set telescope info inside Dset headers

        Arguments::

          telescope : string
               name of telescope

          observatory : string
               name of observatory

          longitude : float
               geodetic longitude, East positive, degrees

          latitude : float
               geodetic latitude, degrees

          height : float
               height in metres
        """

        self.head['TELESCOP'] = (telescope, 'Telescope name')
        self.head['OBSERVAT'] = (observatory, 'Observatory name')
        self.head['LONGITUD'] = (longitude, 'Geodetic longitude, degrees')
        self.head['LATITUDE'] = (latitude, 'Geodetic latitude, degrees')
        self.head['HEIGHT']   = (height, 'Height, metres')

    def gettel(self):
        """
        Gets telescope info from Dset headers, if set. Raises an exception
        if nothing is set, or if either 'LATITUDE' and 'LONGITUD' is not
        set. Other items are optional.

        Returns (telescope, observatory, longitude, latitude, height) where::

          telescope : string
               name of telescope, defaults to None

          observatory : string
               name of observatory, defaults to None

          longitude : float
               geodetic longitude, degrees. Must be set.

          latitude : float
               geodetic latitude, degrees. Must be set.

          height : float
               height in metres, defaults to 0.
        """

        if 'LONGITUD' not in self.head or 'LATITUDE' not in self.head:
            raise DnlError('Dset.getpos: one or both of LONGITUD, LATUTUDE not set')

        if 'TELESCOP' in self.head:
            telescope = self.head['TELESCOP']
        else:
            telescope = None

        if 'OBSERVAT' in self.head:
            observatory = self.head['OBSERVAT']
        else:
            observatory = None

        longitude = self.head['LONGITUD']
        latitude  = self.head['LATITUDE']

        if 'HEIGHT' in self.head:
            height = self.head['HEIGHT']
        else:
            height = 0.

        return (telescope, observatory, longitude, latitude, height)

    def utc2tdb(self):
        """
        Convert X-axis times from MJD (UTC) to BMJD (TDB).
        Target and telescope data must have been set. No checks
        are made for valid MJDs.
        """

        # get positional data
        object, ra, dec, system, pmra, pmdec, epoch, parallax, rv = \
            self.getpos()

        # only ICRS at the moment
        if system != 'ICRS':
            raise DnlError('Dset.utc2tdb: sorry, only ICRS' +
                           ' coordinates supported at present')

        # get observatory information
        telescope, observatory, longitude, latitude, height = \
            self.gettel()

        times = self.x.good_data
        tts,tdbs,btdbs,hutcs,htdbs,vhels,vbars = \
            sla.utc2tdb(times, longitude, latitude, height,
                        ra, dec, pmra, pmdec, parallax, rv)

        # update
        self.x.data[self.x.good] = btdbs
        self.x.label = 'BMJD (TDB)'
        self.x.units = 'days'


    def amass(self):
        """Convert X-axis times from MJD (UTC) to airmass

        Target and observatory data must have been set otherwise a DnlError
        will be thrown.

        """

        # Run checks
        targ = self.checkpos()
        tel  = self.checktel()
        self.checktime()

        # only ICRS at the moment
        if targ['System'] != 'ICRS':
            raise DnlError('Dset.amass: sorry, only' +
                           ' ICRS coordinates supported at present')

        # all tests passed
        for i in xrange(len(self.x.dat)):
            airmass,altmaz,ha,pa,dz = \
                sla.amass(self.x.dat[i],
                          tel['Longitude'], tel['Latitude'], tel['Height'],
                          targ['RA'], targ['Dec'], 0.55, targ['PmRA'],
                          targ['PmDec'], targ['Epoch'], targ['Parallax'],
                          targ['RV'])
            self.x.data[i] = airmass
            self.x.errors  = None

        self.x.label = 'Airmass'
        self.x.units = ''
        del self['Time']

    def xpeak(self, xstart, fwhm, emission=True):
        """
        measures the position of a feature using cross-correlation
        with a gaussian.

        Arguments::

         xstart   : starting position in terms of X
         fwhm     : FWHM of gaussian in terms of X
         emission : emission or absorption

        Returns (xcen,xerr) where::

          xcen : centroid x value
          xerr : 1-sigma uncertainty
        """

        xarr = self.x.data
        index = xarr.searchsorted(xstart)
        if index == 0 or index == len(xarr):
            raise DnlError('Dset.xpeak: xstart out of range of X axis')

        xpos = index -1 + (xstart-xarr[index-1])/(xarr[index]-xarr[index-1])
        yarr = self.y.data
        yerr = self.y.errors
        xcen, xerr = subs.centroid(xpos, fwhm, yarr, emission, yerr)
        index = int(xcen)
        dx    = xarr[index+1]-xarr[index]
        xpos  = xarr[index] + dx*(xcen-index)
        xerr *= abs(dx)
        return (xpos,xerr)

    def fmean(self, filter):
        """
        Computes mean y value folded through a filter. The filter is contained
        in the nx2 2D array filter such that filter[:,0] are the x values,
        while filter[:,1] are the filter throughputs. This is in the form that
        one gets from np.loadtxt given a two column x,y filter table. The
        filter throughputs are interpolated onto the X array of the Dset and
        then the results integrated. The end values will be extended to
        whatever they are at the extremes of xf so you should make sure that
        they go to zero.

        This routine is useful for computing stellar fluxes and magnitudes but
        note that in this case you should have the correct flux for the correct
        x axis (e.g. flambda for wavelength) and you may want to pre-multiply
        your filter response by the wavelengths to get a photon-weighted flux.
        Note that in order to calculate the X-axis decrements, centred
        differences are used and so the first and last pixels are lost.

        No account of units is taken. The integral returned is simply a sum
        over the filter-modified fluxes times the wavelength decrements.
        """

        self.x.data = np.cast[np.float64](self.x.data)
        fint  = subs.linterp(filter[:,0],filter[:,1], self.x.data)
        fflux = (fint*self.y.data)[1:-1]
        dwave = (self.x.data[2:]-self.x.data[:-2])/2.
        return (fflux*dwave).sum()

    def clean(self):
        """Cleans out all bad pixels, defined as those which are bad
        in x and/or y. Leaves in pixels masked at the top-level only.
        """
        ok = self.x.good & self.y.good

        self.x.data = self.x.data[ok]
        if self.x.has_errors: self.x.errors = self.x.errors[ok]
        if self.x.has_mask: self.x.good = self.x.good[ok]

        self.y.data = self.y.data[ok]
        if self.y.has_errors: self.y.errors = self.y.errors[ok]
        if self.y.has_mask: self.y.good = self.y.good[ok]

        if self.has_mask: self._good = self._good[ok]

    def wlcurve(self, fname):
        """
        writes a Dset to an lcurve file. Only data not masked from
        x or y are written. Data masked at the top level are written with
        zero weight.

          fname : string
            name of file to write the data out to.
        """

        if not self.x.has_errors or not self.y.has_errors:
            raise DnlError('require errors in both x & y')

        # only write out good data, set weight on any data
        # masked in the Dset top level to weight=0
        ok   = self.x.good & self.y.good
        good = self.good

        ts    = self.x.data[ok]
        texps = 2.*self.x.errors[ok]
        fs    = self.y.data[ok]
        fes   = self.y.errors[ok]
        goods = good[ok]

        with open(fname,'w') as fp:
            for t,texp,f,fe,good in zip(ts,texps,fs,fes,goods):
                if good:
                    fp.write('{0:15.9f} {1:9.4e} {2:10.4e} {3:9.3e} 1 1\n'.
                             format(t,texp,f,fe))
                else:
                    fp.write('{0:15.9f} {1:9.4e} {2:10.4e} {3:9.3e} 0 1\n'.
                             format(t,texp,f,fe))

    @classmethod
    def rlcurve(cls, fname, data=False):
        """
        reads an lcurve file into a Dset.

          fname : string
            name of file to read the data from

          model : bool
            if True, it is assumed to be a fit and errors are ignored.
        """

        din = np.loadtxt(fname)

        x    = din[:,0]
        xe   = din[:,1]/2. if data else None
        y    = din[:,2]
        ye   = din[:,3] if data else None
        good = din[:,4] > 0. if data else None

        head = fits.Header()
        head['comment'] = 'Data read in from an lcurve file. These files have no headers'
        head['comment'] = 'but the weight and sub-division factor columns are lost.'

        if data:
            xaxis = Axis('Time','days', x, xe)
            yaxis = Axis('Flux','', y, ye)
            ptype = Dset.POINTS
        else:
            xaxis = Axis('Time','days', x)
            yaxis = Axis('Flux','', y)
            ptype = Dset.LINE

        return cls(xaxis, yaxis, head, good, ptype)

    def append(self, dset):
        """
        Appends a Dset to this Dset. Requires the x and
        y Axis components to be errors compatible.
        """
        self.x.append(dset.x)
        self.y.append(dset.y)

        if self.has_mask and dset.has_mask:
            self._good = np.concatenate([self._good, dset._good])
        elif self.has_mask and not dset.has_mask:
            self._good = np.concatenate([self._good,
                                         np.ones(len(dset.x.data),dtype=bool)])
        elif not self.has_mask and dset.has_mask:
            self._good = np.concatenate([np.ones(len(self.x.data),dtype=bool),
                                         dset._good])

def plimits(dsets, x1=None, x2=None, y1=None, y2=None, rmin=0.005):
    """This computes plot limits appropriate for one or more Dsets. It also
    computes offsets to subtract from the data prior to plotting in order to
    avoid overlapping and unreadable axis labels. It produces non-zero offsets
    if the plot limits are much closer to together than the value at the middle
    of an axis.

    The limits are not corrected for the offsets, so it is up to the user to
    do so if the offsets are to be used. The computed limits are determined to
    be the minimum size that encompasses the data.

    Arguments::

       dsets : iterable
           should sequentially produce Dsets

       x1 : float / None
           left-hand X limit

       x2 : float / None
           right-hand X limit

       y1 : float / None
           lower Y limit

       y2 : float / None
           upper Y limit

       rmin : float
           ratio of plot width/mean below which offsets are applied. Should
           be between 0 and 1.

    Returns (xoff,yoff,x1,x2,y1,y2) where xoff and yoff are the suggested
    offsets needed to make plot axes reasonably readable.
    """

    # make iterable
    if isinstance(dsets, Dset): dsets = [dsets,]

    xt1, xt2, yt1, yt2 = None, None, None, None
    for dset in dsets:
        xi1,xi2,yi1,yi2 = dset.plimits()
        if xi1 is not None:
            xt1 = xi1 if xt1 is None else min(xi1, xt1)
            xt2 = xi2 if xt2 is None else max(xi2, xt2)
            yt1 = yi1 if yt1 is None else min(yi1, yt1)
            yt2 = yi2 if yt2 is None else max(yi2, yt2)

    # update x1 etc
    x1 = xt1 if x1 is None else x1
    x2 = xt2 if x2 is None else x2
    y1 = yt1 if y1 is None else y1
    y2 = yt2 if y2 is None else y2

    # calculate offsets
    if abs(x2 - x1) < rmin*abs(x1+x2):
        xoff = x1
    else:
        xoff = 0.

    if abs(y2 - y1) < rmin*abs(y1+y2):
        yoff = y1
    else:
        yoff = 0.

    return (xoff,yoff,x1,x2,y1,y2)

def paxes(device, xoff, yoff, x1, x2, y1, y2,
          xlabel='', ylabel='', tlabel=''):
        """
        Bundles together standard opening of a plot.

         device : string
                  plot device identifier

         xoff : float
                offset to subtract from X values.

         yoff : float
                offset to subtract from Y values

         x1 : float
             left-hand X limit (prior to the application of xoff)

         x2 : float
             right-hand X limit (prior to the application of xoff)

         y1 : float
             lower Y limit (prior to the application of yoff)

         y2 : float
             upper Y limit (prior to the application of yoff)

         xlabel : string
             X-axis label

         ylabel : string
             Y-axis label

         tlabel : string
             Label for top of plot
        """

        # start plot
        pgopen(device)
        plschr(0,1.5)
        plfont(2)
        plwidth(2)
        plcol0(4)
        plenv(x1-xoff,x2-xoff,y1-yoff,y2-yoff,0,0)
        plcol0(2)
        pllab(xlabel, ylabel, tlabel)


