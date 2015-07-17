"""
Basic utility routines and objects.
"""

# This is imported into core, plot etc and should
# be thought of as the lowest level of the import
# hierarchy.

# Exception class
class DnlError(Exception):
    """For throwing exceptions from the dnl module"""
    pass

def chunks(mask):
    """
    Given a 1D boolean array 'mask' this routine generate contiguous
    chunks of 'True' values and returns them one by one as 2-element
    tuples identifying the first and the last+1 indices of each chunk.

    mask : 1D boolean array
           contains True or False values e.g. FFTTTTFFFFTTTTTF

    It should be used in a loop as in::

      for start,end in chunks(mask):

        ... do something with range given ..
    """

    if len(mask.shape) != 1 or mask.dtype != bool or len(mask) < 1:
        raise DnlError('chunks: mask must be a 1D boolean array with at least one element')

    cks = []

    istart = -1
    for i, flag in enumerate(mask):
        if flag and istart < 0:
            istart = i
        elif not flag and istart >= 0:
            yield (istart, i)
            istart = -1

    if istart >= 0:
        yield (istart,len(mask))

if __name__ == '__main__':
    import numpy as np

    a = np.random.uniform(size=1000)
    mask = a < 0.003
    nc = 0
    for start, end in chunks(mask):
        nc += 1
    print 'Found',nc,'chunks'



