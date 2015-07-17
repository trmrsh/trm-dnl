"""
a module to cope with x-y data such as spectra and light-curves.

One often needs to associate two vectors of data representing X-axis
and Y-axis data, sometimes with extra 'header' information, and sometimes
with uncertainties in one or both axes. This module defines a class 'Axis'
to represent each axis and 'Dset' combining both.

Sub-packages
============

dint      -- interactive data analysis (very rough)
mask      -- module for defining masks
molly     -- handles molly data
ulog      -- handles ultracam log files

Classes
=======

Axis     -- class to contain the data of an axis, possibly with errors
DnlError -- exception class for the module
Dset     -- the class that contains complete x, y + header data

Functions
=========

ffit     -- multi-parameter fits to Dsets (not fully developed)
rascii   -- reads in ASCII data to make a Dset
read     -- read a Dset written in binary format
rponto   -- read Dsets from a ponto disk file
rsdss    -- reads an SDSS spectrum fits file into a Dset

"""

from __future__ import division, absolute_import, print_function

# import any C code
from ._dnl import *

# and finally the main stuff of the package
from .core import *

from . import ulog
from . import io


