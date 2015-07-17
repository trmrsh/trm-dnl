"""
wrapper routines to choose between PGPLOT and plplot equivalents.
Names and arguments mostly match the plplot versions except where
there is no exact equivalent for a PGPLOT routine in which case I
use the pgplot name. These are to ease a possible future transition
to plplot only. See documentation on plplot for details of what the
routines do and their arguments.
"""

PLPLOT = False

try:
    if PLPLOT:
        import plplot as plg
    else:
        import ppgplot as plg
        plg.PL_BIN_CENTRED = 1
except:
    if PLPLOT:
        import ppgplot as plg
        plg.PL_BIN_CENTRED = 1
    else:
        import plplot as plg

# NB. This file is imported into "core" and is thus meant to be more basic
from base import *

def pgopen(device):
    if PLPLOT:
        plg.plsdev(device)
        plg.plinit()
    else:
        plg.pgopen(device)

def plbin(x, y, opt):
    if PLPLOT:
        plg.plbin(x, y, opt)
    else:
        # note rather limited interpretation of opt
        # that can have more than 2 values
        plg.pgbin(x, y, opt == pl.PL_BIN_CENTRED)

def plcol0(col):
    if PLPLOT:
        plg.plcol0(col)
    else:
        plg.pgsci(col)

def plend():
    if PLPLOT:
        plg.plend()
    else:
        plg.pgend()

def plenv(x1, x2, y1, y2, just, axis):
    if PLPLOT:
        plg.plenv(x1,x2,y1,y2,just,axis)
    else:
        plg.pgenv(x1,x2,y1,y2,just,axis)

def pgerrx(x1, x2, y, terminal):
    if PLPLOT:
        plg.plsmin(0,terminal)
        plg.plerrx(x1,x2,y)
    else:
        plg.pgerrx(x1,x2,y,terminal)

def pgerry(x, y1, y2, terminal):
    if PLPLOT:
        plg.plsmin(0,terminal)
        plg.plerrx(x,y1,y2)
    else:
        plg.pgerry(x,y1,y2,terminal)

def plfont(ifont):
    if PLPLOT:
        plg.plfont(ifont)
    else:
        plg.pgscf(ifont)

def pllab(xlabel, ylabel, tlabel):
    if PLPLOT:
        plg.pllab(xlabel, ylabel, tlabel)
    else:
        plg.pglab(xlabel, ylabel, tlabel)

def plline(x,y):
    if PLPLOT:
        plg.plline(x,y)
    else:
        plg.pgline(x,y)

def plpoin(x,y,code):
    if PLPLOT:
        plg.plpoin(x,y,code)
    else:
        plg.pgpt(x,y,code)

def plstring(x,y,string):
    if PLPLOT:
        plg.plstring(x,y,string)
    else:
        raise DnlError("PGPLOT has no equivalent of plplot's plstring")

def plwidth(width):
    if PLPLOT:
        plg.plwidth(col)
    else:
        plg.pgslw(width)

def plschr(default, scale):
    if PLPLOT:
        plg.plschr(default, scale)
    else:
        plg.pgsch(scale)

def psetupb():
    """Before-device open setup"""
    if PLPLOT:
        plg.plscolbg(255,255,255)

def psetupa():
    """After-device open setup"""
    """Pre-device open setup"""
    if PLPLOT:
        plg.plwidth(2)
    else:
        plg.pgscf(2)
        plg.pgslw(2)
