"""
Handles definition and application of data masks
"""

import trm.dnl as dnl

class Mask(object):
    """
    Base class of mask objects. Should never be
    directly implmented. Stores a flag indicating
    whether the object is to mask or unmask.
    """

    def __init__(self, mask):
        self.mask = mask

    def app_mask(self, dset, temp=True):
        """Applies the mask to a given Dset.

        dset  -- the data set to apply the mask to.
        temp  -- if True then the temporary data mask will be altered, else
                 the bad data mask will be set.
        """
        raise MaskError('Mask.app_mask not implemented; must be overloaded in derived classes')

class Xmask(Mask):
    """
    Masks in terms of X ranges

    Example: xmask = Xmask(x1, x2, True) will mask data satisfying x1 < x < x2 
    """

    def __init__(self, x1, x2, mask):
        super(Xmask, self).__init__(mask)
        if x2 <= x1:
            raise MaskError('Xmask.__init__: x2 <= x1')
        self.x1 = x1
        self.x2 = x2

    def app_mask(self, dset, temp=True):
        """Masks a range in X of the dset

        dset  -- the data set to apply the mask to.
        temp  -- if True then the temporary data mask will be altered, else
                 the bad data mask will be set.
        """

#        print 'Setting mask from',self.x1,'to',self.x2,'to',not self.mask
        ok = (dset.x.dat > self.x1) & (dset.x.dat < self.x2)
        if temp:
            dset.mask[ok] = not self.mask
        else:
            dset.good[ok] = not self.mask

class Imask(Mask):
    """
    Masks in terms of pixel indices. Starts from C-style index 0.

    Examples:

    imask = Imask(i1, i2, True) will mask data with indices from satisfying i1 <= i <= i2 
    imask = Imask(i, False) will unmask data with index i
    """

    def __init__(self, *args):
        """
        2 cases:
        i, mask      --- mask pixel i
        i1, i2, mask --- mask pixels i1 to i2
        """
        if len(args) == 2:
            (i1,mask) = args
            i2 = i1
        elif len(args) == 3:
            (i1,i2,mask) = args
            if i2 < i1:
                raise MaskError('Imask.__init__: i2 < i1')
        super(Imask, self).__init__(mask)
        self.i1 = i1
        self.i2 = i2

    def app_mask(self, dset, temp=True):
        """Masks a range of pixels in a dset

        dset  -- the data set to apply the mask to.
        temp  -- if True then the temporary data mask will be altered, else
                 the bad data mask will be set.
        """

        if temp:
            dset.mask[self.i1:self.i2+1] = not self.mask
        else:
            dset.good[self.i1:self.i2+1] = not self.mask

class Gmask(list):
    """
    Group of Masks. This is the class which defines a complete mask
    which can be applied to Dset objects. The idea is to apply the individual
    masks in the order defined by this list so that can one can build up 
    a complete mask as a series of masks and un-masks.
    """
    
    def __init__(self, arg=None):
        if arg == None:
            super(Gmask, self).__init__()
        else:
            super(Gmask, self).__init__(arg)

    def app_mask(self, dset, temp=True):
        """Masks a dset 

        dset  -- the data set to apply the mask to.
        temp  -- if True then the temporary data mask will be altered, else
                 the bad data mask will be set.
        """

        for mask in self:
            mask.app_mask(dset, temp)


# Exception class
class MaskError(dnl.DnlError):
    """For throwing exceptions from the dnl.mask module"""
    pass
