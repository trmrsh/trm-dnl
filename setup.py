from distutils.core import setup, Extension
import os, numpy

""" Setup script for dnl"""

library_dirs = []
include_dirs = []

# need to direct to where includes and  libraries are
if os.environ.has_key('TRM_SOFTWARE'):
    library_dirs.append(os.path.join(os.environ['TRM_SOFTWARE'], 'lib'))
    include_dirs.append(os.path.join(os.environ['TRM_SOFTWARE'], 'include'))
else:
    print >>sys.stderr, "Environment variable TRM_SOFTWARE pointing to location of shareable libraries and includes not defined!"

include_dirs.append(numpy.get_include())

dnl = Extension('trm.dnl._dnl',
                define_macros        = [('MAJOR_VERSION', '0'),
                                        ('MINOR_VERSION', '1')],
                undef_macros         = ['USE_NUMARRAY'],
                include_dirs         = include_dirs,
                library_dirs         = library_dirs,
                runtime_library_dirs = library_dirs,
                libraries            = ['subs'],
                sources              = [os.path.join('trm', 'dnl', 'dnl.cc')])

setup(name='trm.dnl',
      version='0.1',
      packages = ['trm', 'trm.dnl', 'trm.dnl.dint',
                  'trm.dnl.molly'],

      ext_modules=[dnl],

      scripts=['scripts/danal.py', 'scripts/fmolly.py',
               'scripts/gmolly.py', 'scripts/mol2fits.py',
               'scripts/ulog2lc.py'],

      author='Tom Marsh',
      description='Python module for 1D data analysis',
      author_email='t.r.marsh@warwick.ac.uk',
      url='http://www.astro.warwick.ac.uk/',
      long_description="""
dnl provides an object and a set of routines to help with the analysis of spectra 
and light curves. It handles x-y array data optionally with errors in both x and y.
""",


      )

