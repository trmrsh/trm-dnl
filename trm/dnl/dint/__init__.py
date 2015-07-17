#!/usr/bin/env python

"""dint is an interactive interface to the dnl routines. It enables a convenient
interface for handling multiple datasets (called 'dsets' in the documentation)
while still allowing access to Python's full command set. The dsets are stored
in slots keyed on integers, held in a sub-class of a dictionary called 'data',
so that e.g. data[10] is the Dset for the 10th slot.

There are a series of routines designed to be called in an interactive context
by passing strings, e.g. 'plot 1', plots the first slot while 'plot 2-5' plots
slots 2 to 5 inclusive. The documentation describes usage in this interactive
context alone. In the in a script context. Arguments listed in [] are hidden
by default but can be revealed by adding 'prompt' when running the command. If
listed as for example [fslot=1], it means that if you do not specify it, it
will return to the value listed, in this case 1.

"""
from __future__ import division, absolute_import, print_function

# core packages
import sys
import cPickle as pickle
import re
import readline

# third party
import numpy as npy
import ppgplot as pg

# mine
import trm.dnl          as dnl
import trm.subs         as subs
import trm.subs.input   as inp

# Module globals
DINT_ENV = None
DINT_DEF = None

# Commands and help strings
COMMANDS = { \
    'airmass'     : 'convert time-series x-axis to airmass', \
    'add'         : 'add slots to each other', \
    'appmask'     : 'apply data mask to slots', \
    'coordinates' : 'set astronomical coordinates', \
    'divide'      : 'divide slots by each other', \
    'help'        : 'list available commands', \
    'lmolly'      : 'load a molly file', \
    'load'        : 'load dsets', \
    'lponto'      : 'load a ponto file', \
    'lultracam'   : 'load an Ultracam file', \
    'mute'        : 'toggle the amount of terminal IO', \
    'verbose'     : 'toggle extra information output in addition to average', \
    'plot'        : 'plot dsets', \
    'setmask'     : 'create a data mask', \
    'show'        : 'show list of dsets', \
    'write'       : 'write dsets to disk', \
    'telescope'   : 'set telescope parameters',
    }

# Command classes as listed by 'help classes'
CLASSES = { \
    #
    'Arithematic' : { \
    'about' : 'arithemetical operations on slots', \
    'list'  : ['divide', 'add']}, \
    #
    'Information' : { \
    'about' : 'commands which give information on slots', \
    'list'  : ['help', 'show']}, \
    #
    'IO'          : { \
    'about' : 'input/output commands', \
    'list'  : ['load', 'lultracam', 'write', 'lponto']}, \
    #
    'Masking'     : { \
    'about' : 'command to mask data', \
    'list'  : ['setmask', 'appmask']}, \
    #
    'Flags'       : { \
    'about' : 'controls look and feel', \
    'list'  : ['mute', 'verbose']}, \
    #
    'Time series' : { \
    'about' : 'time series, light curves', \
    'list'  : ['airmass', 'lultracam']}, \
    #
    'Astronomy' : { \
    'about' : 'astronomical routines', \
    'list'  : ['airmass', 'coordinates', 'lultracam', 'telescope']}, \
    #
    'ULTRACAM' : { \
    'about' : 'all ULTRACAM-associated routines', \
    'list'  : ['lultracam']}, \
    }

def sexec(comms, data=None, group=None, cdata=None, dint_def='.local_defaults'):
    """
    Executes a series of commands in the form of a multi-line string,
    returning the standard data structures (data,group,cdata) in a tuple. This
    is useful at the start of scripts prior to more detailed Python commands
    and thus the 's' at the start of 'sexec'. It can also be used later in
    scripts but you must pass then specify data, group and cdata in the
    function call.
    """

    global DINT_DEF
    if data  is None: data = Data()
    if group is None: group = Group()
    if cdata is None: cdata = Cdata()

    DINT_DEF = dint_def
    if not execute(data,group,cdata,comms):
        print('Aborting execution.')
        exit(1)
    return (data,group,cdata)

def init():
    """
    Initialize standard objects
    """

    # Initialise data, group and cdata
    data   = Data()
    group  = Group()
    cdata  = Cdata()
    return (data,group,cdata)

def execute(data, group, cdata, command):
    """
    Executes a dint command. 'command' is a command string.
    It can be a multiline string in which case it will be split
    on \n and each line will be executed one by one.
    Returns ok meaning whether a command has executed OK
    """

    comlist = command.split('\n')
    comlist = [c for c in comlist if c != '' and not c.isspace()]

    for cmd in comlist:
        valid,nmatch,comm,error = recog(cmd)
        if cdata.verbose:
            print('\n' + cmd)
        if valid:
            try:
                # recognised special commands
                if comm.startswith('airmass'):
                    airmass(comm, data, group, cdata)
                elif comm.startswith('add'):
                    add(comm, data, group, cdata)
                elif comm.startswith('appmask'):
                    appmask(comm, data, group, cdata)
                elif comm.startswith('coordinates'):
                    coordinates(comm, data, group, cdata)
                elif comm.startswith('divide'):
                    divide(comm, data, group, cdata)
                elif comm.startswith('help'):
                    help(comm, cdata)
                elif comm.startswith('load'):
                    load(comm, data, cdata)
                elif comm.startswith('lponto'):
                    lponto(comm, data, cdata)
                elif comm.startswith('lultracam'):
                    lultracam(comm, data, cdata)
                elif comm.startswith('mute'):
                    cdata.mute = not cdata.mute
                    if cdata.mute:
                        print('Terminal output muted')
                    else:
                        print('Full terminal output')
                elif comm.startswith('plot'):
                    plot(comm, data, group, cdata)
                elif comm.startswith('setmask'):
                    setmask(comm, data, cdata)
                elif comm.startswith('show'):
                    data.list()
                elif comm.startswith('telescope'):
                    telescope(comm, data, group, cdata)
                elif comm.startswith('verbose'):
                    cdata.verbose = not cdata.verbose
                    if cdata.verbose:
                        print('When not muted there will be more' +
                              ' terminal output than normal')
                    else:
                        print('Standard terminal output set')
                elif comm.startswith('write'):
                    write(comm, data, group, cdata)

                elif comm == 'exit' or comm == 'quit':
                    print('If you want to quit the program, hit ctrl-D')
                else:
                    raise DintError('Command = "' + comm.split()[0] + '" not yet implemented.')
                ok = True
            except (DintError, dnl.DnlError, subs.SubsError, inp.InputError), err:
                print(err)
                ok = False
                break
        else:

            # can't identify command; have a go with the standard Python interpreter
            try:
                exec(cmd)
                ok = True
            except:
                #            print 'Command "' + command + '" first failed with the error:\n' + error
                #            print 'Trying as a normal Python command then failed with the following error: '
                print(sys.exc_info())
                ok = False
                break

    return ok


class Data(subs.Odict):
    """
    Class to hold the slots. Restricts the type of item stored
    which must be a Dset and the keys which must be positive
    integers
    """

    def __init__(self):
        subs.Odict.__init__(self)

    def __setitem__(self, key, item):
        if not isinstance(key, int):
            raise DintError('dint.Data.__setitem__: key = ' + str(key) +
                            ' is not an integer')
        if key < 1:
            raise DintError('dint.Data.__setitem__: key = ' + str(key) +
                            ' < 1')
        if not isinstance(item, dnl.Dset):
            raise DintError('dint.Data.__setitem__: item is not a Dset')

        subs.Odict.__setitem__(self, key, item)

    def list(self):
        """Lists the slots with some basic information"""
        for key, dset in self.iteritems():
            print('Slot: {0:5d}, {1:s}'.format(key, dset.oneLine()))

class Group(dict):
    """
    Class to hold slot group, which are lists of slot number keyed on the
    group name.
    """

    def __init__(self):
        dict.__init__(self)

    def __setitem__(self, key, item):
        if not isinstance(item, list):
            raise DintError('dint.Group.__setitem__: item is not a list')
        for i in item:
            if not isinstance(i, int) or i < 1:
                raise DintError('dint.Group.__setitem__: element = ' +
                                str(i) + ' is not a positive integer')
        dict.__setitem__(self, key, item)

class Cdata(object):
    """
    Class to hold control data to pass to programs which determines such
    things as the amount of terminal output.

    mute     -- True/False to mute terminal output or not.
    verbose  -- True/False to add extra terminal output or not

    """

    def __init__(self):
        self.mute    = False
        self.verbose = False

    def __str__(self):
        return 'Mute    = ' + str(self.mute)
        return 'Verbose = ' + str(self.verbose)

def recog(command):
    """
    valid,nmatch,comm,message = recog(command) -- command checker.

    Checks whether command is recognised and valid, based on a list of commands

    Arguments:

    command  -- complete command input string

    Returns (valid,nmatch,comm,message)::

      valid : bool
         True/False whether command is valid. This means it is unique
         and recognised so can safely be selected with a simple if/elif
         construction.

      nmatch : int
         number of matches (must = 1 for a valid command).

      comm : string
         returns 'command' but with the command name modified to its
         full length to help operations downstream.

      message : string
         information message, useful in case of error.
    """

    global COMMANDS
    valid   = False
    nmatch  = 0
    cfind   = ''

    m = re.compile('([a-zA-Z0-9]+)').match(command)
    if not m:
        return (valid,nmatch,cfind,
                'No valid command name at start of command = "' +
                command + '"')

    start  = m.group(0)
    for comm in COMMANDS.keys():
        if comm.startswith(start):
            nmatch += 1
            if nmatch > 1:
                break
            cfind = comm

    if nmatch > 1:
        message = 'More than one command matches start of command = "' + \
                  command + '"'
    elif nmatch == 0:
        message = 'No command found matching start of command = "' + \
                  command + '"'
    else:
        message = 'Unique match found for command = "' + command + '"'
        ce = re.compile('[a-zA-Z0-9]+\s*=')
        if ce.match(command):
            message += ' but the equals after the command invalidates it.'
        else:
            valid = True
            # fix up the command string to have the full command at the start.
            command = command.replace(start, cfind, 1)

    return (valid,nmatch,command,message)

def interp_slots(slots, dexist, data, group=None, nfind=None, strict=False):
    """
    slist = interp_slots(slots, group, dexist) -- interprets slot
    ranges/group names

    Arguments::

      slots : string
        slot e.g. '11', range e.g. '21-23' or the name of a group if group
        is not None

      dexist : bool
        True if the slots defined must have data already in which case the
        list returned will only include such slots.

    data : Data
        object containing the slots

    group : Group
        a Group of groups, each entry giving a list of slots. e.g.
        group['ippeg'] might equal [1, 3, 4, 5] for example.

    nfind : int
        set this if you require a specific number of slots. If slot is a
        single number then this will automatically add the right extra number.

    strict : bool
        True if you want an error raised in the case of dexist = True and an
        empty slot being encountered

    Returns::

      slist : list
        the slot numbers

    Throws a DintError if there are problems such as empty slots when dexist
    and strict are True or if no slots are found.
    """

    m = re.compile('([0-9]+)$').match(slots)
    if m:
        if nfind:
            slist = range(int(slots),int(slots)+nfind)
        else:
            slist = [int(slots)]
    else:
        m = re.compile('([0-9]+)-([0-9]+)$').match(slots)
        if m:
            first = int(m.group(1))
            last  = int(m.group(2))
            if first > last:
                raise DintError('interp_slots: invalid slot range or name = ' + slots)
            slist = range(first, last+1)
        elif group is not None:
            if slots.find('-') > -1:
                raise DintError('interp_slots: invalid slot range or name = ' + slots)
            else:
                if slots not in group:
                    raise DintError('interp_slots: could not find group = ' + slots)
                slist = group[slots]
        else:
            raise DintError('interp_slots: no group defined')

    # check for existence of data
    if dexist:
        if strict:
            for slot in slist:
                if slot not in data:
                    raise DintError('interp_slots: slot ' + str(slot) + \
                                    ' does not exist')

        # cut down the list to ones with data
        slist = [slot for slot in slist if slot in data]

    if nfind is not None and len(slist) != nfind:
        raise DintError('interp_slots: found ' + str(len(slist)) + \
                        ' slots when ' + str(nfind) + ' were required.')

    if len(slist) == 0:
        raise DintError('interp_slots: no slots identified')

    return slist

class DintError(Exception):
        """Exception class for interactive data analysis package"""
        pass

# now come a series of possible commands

def airmass(command, data, group, cdata):
    """
    Converts the X axis of astronomical time series slots to airmass, where
    airmass is a measure of the amount of air being looked through. The object
    coordinates must be set with 'coordinates' and the telescope position must
    be set with 'telescope'. If it encounters an error it will abort the loop
    rather than try to continue.

    Interactive usage:

    airmass slots

    Arguments:

    slots     -- slot, slot range or group which will be processed.

    Script usage:

    command  -- contains a string as used interactively
    data     -- a Data containing Dsets.
    group    -- a Group containing slot list groups

    Returns True/False according to whether it ran OK
    """

    # generate arguments
    inpt = inp.Input(DINT_ENV, DINT_DEF, inp.clist(command))

    # register parameters
    inpt.register('slots',  inp.Input.LOCAL,  inp.Input.PROMPT)

    # get inputs
    slots  = inpt.get_value('slots', 'slots to mask', '1-1000')
    slist  = interp_slots(slots, True, data, group)

    for slot in slist:
        data[slot].amass()
        if not cdata.mute:
            print('Converted X-axis of slot ' + str(slot) + ' to airmass.')

def add(command, data, group, cdata):
    """
    Add one set of slots to another

    Interactive usage:

    add inslot1 inslot1 outslot

    Arguments::

      inslot1 : string
        range of slots as in '1-10', or just a single slot '9' to plot, or
        a group name such as 'ippeg'.

      inslot2 : string
        range of slots as in '1-10', or just a single slot '9' to plot, or
        a group name such as 'ippeg', must match number of inslot1 slots

    outslot : int
        first slot for output (sequential from there)

    The headers of the output slots will be taken from the inslot1 slots.
    """

    # generate arguments
    inpt = inp.Input(DINT_ENV, DINT_DEF, inp.clist(command))

    # register parameters
    inpt.register('inslot1',  inp.Input.LOCAL, inp.Input.PROMPT)
    inpt.register('inslot2',  inp.Input.LOCAL, inp.Input.PROMPT)
    inpt.register('outslot',  inp.Input.LOCAL, inp.Input.PROMPT)

    # get inputs
    slots1  = inpt.get_value('inslot1', 'first set of slots', '1')
    slist1  = interp_slots(slots1, True, data, group)
    slots2  = inpt.get_value('inslot2', 'slots to add to the first set', '2')
    slist2  = interp_slots(slots2, True, data, group, nfind=len(slist1))
    slots3  = inpt.get_value('outslot', 'first output slot', '3')
    slist3  = interp_slots(slots3, False, data, group, nfind=len(slist1))

    for i in range(len(slist1)):
        data[slist3[i]] = data[slist1[i]] + data[slist2[i]]
        if not cdata.mute:
            print('Slot ' + str(slist3[i]) + ' = slot ' +
                  str(slist1[i]) + ' +  slot ' + str(slist2[i]))

    return True

def appmask(command, data, group, cdata):
    """
    Applies a data mask from a file created e.g. by 'setmask'

    Interactive usage::

      appmask slots mfile [type=temp]

    Arguments::

      slots : string
         slot, slot range or group which will be masked.

      mfile : string
         mask file

      type : string
         either the temporary data mask will be set ('temp') or the
         bad pixel mask ('bad') will be set, excluding the pixels from
         all plots and processing. See trm.dnl.Dset to understand the
         difference.

    Script usage::

       command : string
          contains a string as used interactively

       data : Data
          contains Dsets.

       group : Group
          contains slot list groups

    Returns True/False according to whether it ran OK
    """

    # generate arguments
    inpt = inp.Input(DINT_ENV, DINT_DEF, inp.clist(command))

    # register parameters
    inpt.register('slots',  inp.Input.LOCAL,  inp.Input.PROMPT)
    inpt.register('mfile',  inp.Input.GLOBAL, inp.Input.PROMPT)
    inpt.register('type',   inp.Input.LOCAL,  inp.Input.HIDE)

    # get inputs
    slots  = inpt.get_value('slots', 'slots to mask', '1-1000')
    slist  = interp_slots(slots, True, data, group)
    mfile  = inpt.get_value('mfile', 'mask file', subs.Fname('mask','.msk'))
    inpt.set_default('type', 'temp')
    mtype  = inpt.get_value('type',  'temp(oraray) or bad mask?', 'temp',
                            lvals=['temp', 'bad'])

    # read and apply mask
    mptr  = open(mfile,'r')
    gmask = pickle.load(mptr)
    mptr.close()

    if mtype == 'temp':
        tmask = True
    else:
        tmask = False

    for i in slist:
        gmask.app_mask(data[i], tmask)
        if not cdata.mute:
            print('Applied mask from ' + str(mfile) +
                  ' to slot ' + str(i))

def coordinates(command, data, group, cdata):
    """
    Adds astronomical position information to slots.

    Interactive usage:

    coordinates slots target position [rapm=0 decpm=0 epoch=2000.0 parallax=0 rv=0]

    Arguments:

    slots     -- slot, slot range or group which will be processed.
    target    -- target name
    position  -- RA, dec as in '12 34 45.56 -00 12 34.2' or decimal degrees '234.5645 -1.2'. Optionally
                 you can add any of 'ICRS', 'J2000' or 'B1950' at the end, but note that only 'ICRS' (the
                 default) is supported as of Aug 2008. See trm.subs.str2radec for full format info.
    epoch     -- Julian epoch of coordinates (only matters if there is proper motion and RV)
    pmra      -- proper motion in RA (arcsec/year)
    pmdec     -- proper motion in Dec (arcsec/year)
    parallax  -- parallax, (arcsec)
    rv        -- radial velocity, (km/s)

    Script usage:

    command  -- contains a string as used interactively
    data     -- a Data containing Dsets.
    group    -- a Group containing slot list groups

    Returns True/False according to whether it ran OK
    """

    # generate arguments
    inpt = inp.Input(DINT_ENV, DINT_DEF, inp.clist(command))

    # register parameters
    inpt.register('slots',     inp.Input.LOCAL,  inp.Input.PROMPT)
    inpt.register('target',    inp.Input.LOCAL,  inp.Input.PROMPT)
    inpt.register('position',  inp.Input.LOCAL,  inp.Input.PROMPT)
    inpt.register('epoch',     inp.Input.LOCAL,  inp.Input.HIDE)
    inpt.register('pmra',      inp.Input.LOCAL,  inp.Input.HIDE)
    inpt.register('pmdec',     inp.Input.LOCAL,  inp.Input.HIDE)
    inpt.register('parallax',  inp.Input.LOCAL,  inp.Input.HIDE)
    inpt.register('rv',        inp.Input.LOCAL,  inp.Input.HIDE)

    # get inputs
    slots    = inpt.get_value('slots', 'slots to mask', '1-1000')
    slist    = interp_slots(slots, True, data, group)
    target   = inpt.get_value('target',  'name of target', 'M31')
    position = inpt.get_value('position',  'celestial coordinates of target', '12 34 45.67 -00 11 23.4')
    ra,dec,system = subs.str2radec(position)

    inpt.set_default('epoch', 2000.0)
    epoch = inpt.get_value('epoch', 'Julian epoch of coordinates', 2000.0)
    inpt.set_default('pmra', 0.)
    pmra = inpt.get_value('pmra', 'proper motion in RA (arcsec/yr)', 0.)
    inpt.set_default('pmdec', 0.)
    pmdec = inpt.get_value('pmdec', 'proper motion in Dec (arcsec/yr)', 0.)
    inpt.set_default('parallax', 0.)
    parallax = inpt.get_value('parallax', 'parallax (arcsec)', 0.)
    inpt.set_default('rv', 0.)
    rv = inpt.get_value('rv', 'radial velocity (km/s)', 0.)

    for slot in slist:
        data[slot].setpos(target, ra, dec, system, pmra, pmdec, epoch, parallax, rv)
        if not cdata.mute:
            print('Added astronomical target data to slot ' + str(slot))

def divide(command, data, group, cdata):
    """
    Divide one set of slots by another

    Interactive usage:

    divide inslot1 inslot1 outslot

    Arguments:

    inslot1 -- range of slots as in '1-10', or just a single slot '9' to plot, or a group
               name such as 'ippeg'.
    inslot2 -- range of slots as in '1-10', or just a single slot '9' to plot, or a group
               name such as 'ippeg', must match number on inslot1 slots
    outslot -- first slot for output (sequential from there)

    The headers of the output slots will be taken from the inslot1 slots.
    """

    # generate arguments
    inpt = inp.Input(DINT_ENV, DINT_DEF, inp.clist(command))

    # register parameters
    inpt.register('inslot1',  inp.Input.LOCAL, inp.Input.PROMPT)
    inpt.register('inslot2',  inp.Input.LOCAL, inp.Input.PROMPT)
    inpt.register('outslot',  inp.Input.LOCAL, inp.Input.PROMPT)

    # get inputs
    slots1  = inpt.get_value('inslot1', 'slots to divide', '1')
    slist1  = interp_slots(slots1, True, data, group)
    slots2  = inpt.get_value('inslot2', 'slots to divide by', '2')
    slist2  = interp_slots(slots2, True, data, group, nfind=len(slist1))
    slots3  = inpt.get_value('outslot', 'first output slots', '3')
    slist3  = interp_slots(slots3, False, data, group, nfind=len(slist1))

    for i in range(len(slist1)):
        data[slist3[i]] = data[slist1[i]] / data[slist2[i]]
        if not cdata.mute:
            print('Slot ' + str(slist3[i]) + ' = slot ' +
                  str(slist1[i]) + ' /  slot ' + str(slist2[i]))

    return True

def help(command, cdata):

    """
    Provides help. Possible ways to call it are:

    'help'          -- full command list
    'help classes'  -- list of command classes
    'help load'     -- help on command 'load'
    'help lo.*'     -- commands matching regular expression 'lo.*'
    'help class IO' -- classes matching regular expression 'lo.*'

    Note that 'help load' would also list help on 'upload' because
    of the regular expression matching.
    """

    global COMMANDS, CLASSES

    try:

        # generate arguments
        inpt = inp.Input(DINT_ENV, DINT_DEF, inp.clist(command))

        # register parameters
        inpt.register('what',   inp.Input.LOCAL, inp.Input.HIDE)
        inpt.register('class',  inp.Input.LOCAL, inp.Input.HIDE)

        # get inputs
        inpt.set_default('what', '')
        what = inpt.get_value('what',  'what do you want help on?', 'classes')

        if what == 'classes':
            keys = CLASSES.keys()
            keys.sort()
            for key in keys:
                print('{0:<10s} -- {1:s}'.format(key, CLASSES[key]['about']))
            print("\nType 'help class classname' for information" +
                  " on a specific class")

        elif what == 'class':
            clss = inpt.get_value('class',  'which class to list (regular expressions supported)', 'IO')
            m = re.compile('clss')
            keys = [key for key in CLASSES.keys() if m.search(key)]
            if len(keys):
                keys.sort()
                for key in keys:
                    found = True
                    comms = CLASSES[key]['list']
                    comms.sort()
                    print('\nCommands of the ' + key +
                          ' class are as follows:\n')
                    for comm in comms:
                        print('{0:<10s} -- {1:s}'.format(key, COMMANDS[key]))
            else:
                print('Sorry; no class matching ' + clss + ' was found.')
        else:
            if what == '':
                comms = COMMANDS.keys()
            else:
                m = re.compile(what)
                comms = [comm for comm in COMMANDS.keys() if m.search(comm)]
            if len(comms):
                comms.sort()
                for comm in comms:
                    print('{0:<11s} -- {1:s}'.format(comm, COMMANDS[comm]))
            else:
                print('Sorry; no command matching ' + what + ' was found.')

        return True

    except inp.InputError, err:
        print('InputError exception raised: ' + str(err))

    return False

def lmolly(command, data, cdata):
    """
    Loads data from a molly file

    Interactive usage:

    lmolly mol slots fslot

    Arguments:

    mol     -- name of molly file
    slots   -- range of slots as in '1-10', or just a single slot '9' to load data into.
    fslot   -- first slot to read in file, starting from 1
    """

    import trm.dnl.molly as molly

    # generate arguments
    inpt = inp.Input(DINT_ENV, DINT_DEF, inp.clist(command))

    # register parameters
    inpt.register('mol',    inp.Input.LOCAL, inp.Input.PROMPT)
    inpt.register('slots',  inp.Input.LOCAL, inp.Input.PROMPT)
    inpt.register('fslot',  inp.Input.LOCAL, inp.Input.PROMPT)

    # get inputs
    mol  = inpt.get_value('mol',  'molly file name', subs.Fname('example', '.mol'))

    slots  = inpt.get_value('slots', 'slots to load data into', '1-1000')
    slist  = interp_slots(slots, False, data)
    fslot  = inpt.get_value('fslot', 'first slot to read from file', 1, 1, 1000000)

    # load molly data
    try:
        mf = open(mol, 'rb')
        nload = 0

        # junk the first fslot-1
        for i in range(fslot-1):
            print('skipping a molly spectrum')
            molly.skip_molly(mf)

        # save the rest
        for slot in slist:
            mspec = molly.rmspec(mf)
            if mspec:
                data[slot] = mspec
                nload += 1
                if not cdata.mute:
                    print('Loaded ponto data into slot',slot)
            else:
                raise EOFError
        mf.close()
    except EOFError:
        print('End-of-file reached.')
        mf.close()

    # loading no slots is flagged as an error
    if nload == 0:
        raise DintError('load: no slots were loaded')
    elif not cdata.mute:
        print(nload,'slots were loaded.')

def load(command, data, cdata):
    """
    Loads dsets. Interactive usage::

      load dset slots fslot

    Arguments::

      dset : string
        name of FITS file containing dsets. '.fits' will be added if
        not already present

      slots : string
        range of slots as in '1-10', or just a single slot '9' to load
        data into. It is OK to specify more slots than the file contains.
        Unused slots will be unchanged.

      fslot : int
        first slot to read from file, starting from 1
    """

    # generate arguments
    inpt = inp.Input(DINT_ENV, DINT_DEF, inp.clist(command))

    # register parameters
    inpt.register('dset',   inp.Input.LOCAL, inp.Input.PROMPT)
    inpt.register('slots',  inp.Input.LOCAL, inp.Input.PROMPT)
    inpt.register('fslot',  inp.Input.LOCAL, inp.Input.PROMPT)

    # Get inputs
    ndset  = inpt.get_value('dset',  'file containing one or more dsets',
                            subs.Fname('example', '.fits'))
    slots  = inpt.get_value('slots', 'slots to load dsets into', '1-1000')
    slist  = interp_slots(slots, False, data)
    fslot  = inpt.get_value('fslot', 'first slot to read from file',
                            1, 1, 1000000)

    # load dsets
    nload = 0
    for i, dset in enumerate(dnl.io.grfits(ndset, fslot)):
        data[slist[i]] = dset
        nload += 1
        if not cdata.mute:
            print('Slot {0:5d}, {1:s}'.format(slist[i],dset.oneLine()))

    # loading no slots is flagged as an error
    if nload == 0:
        raise DintError('load: no slots were loaded')

def lponto(command, data, cdata):
    """
    Loads data from a ponto file

    Interactive usage:

    lponto pnt slots fslot

    Arguments:

    pnt     -- name of ponto file
    slots   -- range of slots as in '1-10', or just a single slot '9' to load data into.
    fslot   -- first slot to read in file, starting from 1
    """

    import trm.subs.cpp as cpp

    # generate arguments
    inpt = inp.Input(DINT_ENV, DINT_DEF, inp.clist(command))

    # register parameters
    inpt.register('pnt',    inp.Input.LOCAL, inp.Input.PROMPT)
    inpt.register('slots',  inp.Input.LOCAL, inp.Input.PROMPT)
    inpt.register('fslot',  inp.Input.LOCAL, inp.Input.PROMPT)

    # get inputs
    pnt  = inpt.get_value('pnt',  'ponto file name', subs.Fname('example', '.pnt'))

    slots  = inpt.get_value('slots', 'slots to load data into', '1-1000')
    slist  = interp_slots(slots, False, data)
    fslot  = inpt.get_value('fslot', 'first slot to read from file', 1, 1, 1000000)

    # load ponto data
    try:
        (fobj,endian) = cpp.open_ponto(str(pnt))
        nload = 0

        # junk the first fslot-1
        for i in range(fslot-1):
            print('skipped a dset')
            dnl.rponto(fobj,endian)

        # save the rest
        for slot in slist:
            data[slot] = dnl.rponto(fobj,endian)
            nload += 1
            if not cdata.mute:
                print('Loaded ponto data into slot',slot)
    except cpp.CppError, err:
        fobj.close()
        raise DintError(str(err))
    except EOFError:
        print('End-of-file reached.')
    fobj.close()

    # loading no slots is flagged as an error
    if nload == 0:
        raise DintError('load: no slots were loaded')
    elif not cdata.mute:
        print(nload,'slots were loaded.')

def lultracam(command, data, cdata):
    """
    Loads dsets from an Ultracam file

    This runs on a file created from concatenating ULTRACAM ASCII log files
    into a single binary file using the script 'cultracam.py'. Such files store
    all data from all apertures and CCDs. This program will try to load all
    apertures from all CCDs, depending on how many slots are being loaded.

    Interactive usage:

    lultracam ult slots what bmax emax

    Arguments:

    ult     -- name of a binary file ending '.ult' created by the script lultracam.py from ultracam log files
    slots   -- range of slots as in '1-10', or just a single slot '9' to load data into.
    what    -- 'c' = count rate light curves, 'f' = fwhm seeing in pixels, 'xm' = measured X positions, 'ym' = measured Y positions
    bmax    -- maximum bad pixel in aperture before data are marked as bad
    emax    -- maximum bad pixel in aperture before data are marked as bad
    """

    # generate arguments
    inpt = inp.Input(DINT_ENV, DINT_DEF, inp.clist(command))

    # register parameters
    inpt.register('ult',    inp.Input.LOCAL, inp.Input.PROMPT)
    inpt.register('slots',  inp.Input.LOCAL, inp.Input.PROMPT)
    inpt.register('what',   inp.Input.LOCAL, inp.Input.PROMPT)
    inpt.register('bmax',   inp.Input.LOCAL, inp.Input.HIDE)
    inpt.register('emax',   inp.Input.LOCAL, inp.Input.HIDE)

    # get parameters
    uname  = inpt.get_value('ult',  'ultracam binary file', subs.Fname('example', '.ult'))
    slots  = inpt.get_value('slots', 'slots to load dsets into', '1-1000')
    slist  = interp_slots(slots, False, data)
    what   = inpt.get_value('what', 'C(ounts), F(whm), XM(easured), YM(easured)', 'c', lvals=['c','C','f','F', 'xm', 'XM', 'ym', 'YM'])
    bmax   = inpt.get_value('bmax', 'maximum bad pixel in aperture', 50)
    emax   = inpt.get_value('emax', 'maximum error flag', 5)

    # read the ultracam binary file
    iptr = open(uname,'rb')
    ult  = pickle.load(iptr)
    iptr.close()

    # extract all Dsets, first work out maximum number of apertures for any one CCD
    # because it is more convenient to extract all CCDs of a given aperture rather than
    # all apertures of a given CCD
    maxaps = npy.array([len(ult.x[nccd]) for nccd in ult.x.keys()]).max()

    nload = 0
    for nap in range(maxaps):
        for nccd in ult.x.keys():
            if nap < len(ult.x[nccd]):
                data[slist[nload]] = ult.tseries(nccd,nap+1,what.lower(),bmax,emax)
                if not cdata.mute:
                    print('Loaded dset for CCD ' + str(nccd) +
                          ' aperture ' + str(nap+1) +
                          ' into slot',slist[nload])
                nload += 1

    # loading no slots is flagged as an error
    if nload == 0:
        raise DintError('load: no slots were loaded')

def plot(command, data, group, cdata):
    """
    Plots dsets from a dictionary called data.

    Interactive usage:

    plot slots [device x1 x2 y1 y2 xoff ysep=0 pmask]

    Arguments:

    slots   -- range of slots as in '1-10', or just a single slot '9' to plot, or a group
               name such as 'ippeg'.
    device  -- plot device (e.g. '/xs', '3/xs', 'hardcopy.ps/cps')
    x1      -- left-hand plot limit
    x2      -- right-hand plot limit. Set = x1 for automatic determination
    y1      -- lower plot limit
    y2      -- upper plot limit.  Set = y1 for automatic determination
    xoff    -- offset to start X axis from, 0 for automatic determination.
    ysep    -- separation in y
    pmask   -- whether to plot masked data
    """

    # generate arguments
    inpt = inp.Input(DINT_ENV, DINT_DEF, inp.clist(command))

    # register parameters
    inpt.register('slots',  inp.Input.LOCAL, inp.Input.PROMPT)
    inpt.register('device', inp.Input.LOCAL, inp.Input.HIDE)
    inpt.register('x1',     inp.Input.LOCAL, inp.Input.HIDE)
    inpt.register('x2',     inp.Input.LOCAL, inp.Input.HIDE)
    inpt.register('y1',     inp.Input.LOCAL, inp.Input.HIDE)
    inpt.register('y2',     inp.Input.LOCAL, inp.Input.HIDE)
    inpt.register('xoff',   inp.Input.LOCAL, inp.Input.HIDE)
    inpt.register('ysep',   inp.Input.LOCAL, inp.Input.HIDE)
    inpt.register('pmask',  inp.Input.LOCAL, inp.Input.HIDE)

    # get inputs
    slots  = inpt.get_value('slots', 'slots to plot', '1')
    slist  = interp_slots(slots, True, data, group)

    device = inpt.get_value('device', 'plot device', '/xs')
    x1   = inpt.get_value('x1', 'left-hand plot limit', 0.0)
    x2   = inpt.get_value('x2', 'right-hand plot limit', 0.0)
    y1   = inpt.get_value('y1', 'lower plot limit', 0.0)
    y2   = inpt.get_value('y2', 'upper plot limit', 0.0)
    xoff = inpt.get_value('xoff', 'X offset', 0.0)
    inpt.set_default('ysep', 0.0)
    ysep = inpt.get_value('ysep', 'vertical separation between successive dsets', 0.0)
    pmask = inpt.get_value('pmask', 'do you want to plot the masked data too?', True)

    # Determine limits automatically if required
    if xoff !=0. or x1 == x2 or y1 == y2:
        xa1 = None
        xa2 = None
        ya1 = None
        ya2 = None
        yadd = 0.
        for i in slist:
            (xi1,xi2,yi1,yi2) = data[i].plimits()
            if xa1 is None:
                xa1 = xi1
            else:
                xa1 = min(xa1, xi1)

            if xa2 is None:
                xa2 = xi2
            else:
                xa2 = max(xa2, xi2)

            if ya1 is None and yi1 is not None:
                ya1 = yi1 + yadd
            elif yi1 is not None:
                ya1 = min(ya1, yi1 + yadd)

            if ya2 is None and yi2 is not None:
                ya2 = yi2 + yadd
            elif yi2 is not None:
                ya2 = max(ya2, yi2 + yadd)

            yadd += ysep

        if xa1 is None or xa2 is None or ya1 is None or ya2 is None:
            raise DintError('plot: no automatic limits could be evaluated; possibly no good data to plot?')

        if xoff == 0.0 and (xa2 - xa1) < (xa1+xa2)/2./100.:
            xoff = xa1
        xa1 -= xoff
        xa2 -= xoff

        if x1 == x2:
            x1 = xa1
            x2 = xa2

        if y1 == y2:
            y1 = ya1
            y2 = ya2

    try:

        pg.pgopen(device)
        pg.pgsci(4)
        pg.pgenv(x1, x2, y1, y2, 0, 0)
        pg.pgsci(2)
        first   = data[slist[0]]
        xlabel  = first.x.label if xoff == 0 else first.x.label + '-' + str(xoff)
        xlabel += ' (' + first.x.units + ')'
        ylabel  = first.y.label + ' (' + first.y.units + ')'
        pg.pglab(xlabel, ylabel, first.title)

        yadd = 0.
        for slot in slist:
            data[slot].plot(xoff,yoff=-yadd,masked=pmask)
            if not cdata.mute:
                print('Plotted slot ' + str(slot))
            yadd += ysep

        pg.pgclos()

    except pg.ioerror, err:
        raise DintError(str(err))


def setmask(command, data, cdata):
    """
    Sets the mask on a dset and dumps a file containing the mask which can be applied
    to other dsets using 'appmask'. This is an interactive routine which will request
    input from the user and is better not used in batch processing. Repeated calls of
    this routine can be used to build complex masks. The masks are always applied in
    the original order so that you can mask then partially unmask for example. Note that
    whatever slot you choose to define the mask will always end up masked; if you don't
    want this you may want to make a copy.

    Interactive usage:

    setmask slot mfile append [device reset x1 x2 y1 y2] mask type

    Arguments:

    slot      -- an example slot to plot.
    mfile     -- mask file
    append    -- append to an old mask file if possible
    device    -- plot device, e.g. '/xs'
    reset     -- rest plot limits or not
    x1        -- left X plot limit
    x2        -- right X plot limit
    y1        -- bottom Y plot limit
    y2        -- top Y plot limit
    mask      -- mask 'M', or unmask 'U' or quit 'Q'.
    type      -- type of mask: 'X' masks using ranges in X

    Mask types:

    X  -- mask a range in X
    Y  -- mask a range in Y
    I  -- mask a range of pixel indices
    P  -- mask in 'phase', i.e. a range that repeats periodically.
    """

    import trm.dnl.mask as mask

    # generate arguments
    inpt = inp.Input(DINT_ENV, DINT_DEF, inp.clist(command))

    # register parameters
    inpt.register('slot',   inp.Input.LOCAL, inp.Input.PROMPT)
    inpt.register('mfile',  inp.Input.GLOBAL,  inp.Input.PROMPT)
    inpt.register('append', inp.Input.LOCAL,  inp.Input.PROMPT)
    inpt.register('device', inp.Input.LOCAL,  inp.Input.HIDE)
    inpt.register('reset',  inp.Input.LOCAL,  inp.Input.HIDE)
    inpt.register('x1',     inp.Input.LOCAL,  inp.Input.HIDE)
    inpt.register('x2',     inp.Input.LOCAL,  inp.Input.HIDE)
    inpt.register('y1',     inp.Input.LOCAL,  inp.Input.HIDE)
    inpt.register('y2',     inp.Input.LOCAL,  inp.Input.HIDE)
    inpt.register('mask',   inp.Input.LOCAL,  inp.Input.PROMPT)
    inpt.register('type',   inp.Input.LOCAL,  inp.Input.PROMPT)

    # get inputs
    slots  = inpt.get_value('slot', 'slot to plot for mask definition', '1')
    slist  = interp_slots(slots, True, data, nfind=1)
    dset   = data[slist[0]]

    device = inpt.get_value('device', 'plot device', '/xs')

    # mask file
    mfile  = inpt.get_value('mfile',  'mask file to save results to', subs.Fname('mask','.msk', subs.Fname.NEW))
    append = inpt.get_value('append', 'add to an old mask file if possible', True)

    if append and mfile.exists():
        mptr  = open(mfile,'rb')
        gmask = pickle.load(mptr)
        gmask.app_mask(dset)
        mptr.close()
    else:
        gmask = mask.Gmask()

    # other parameters
    reset = inpt.get_value('reset',   'reset plot limits automatically?', True)

    # compute default limits
    (x1,x2,y1,y2) = dset.plimits()
    if (x2 - x1) < (x1+x2)/2./100.:
        xoff = x1
        x1   = 0.
        x2  -= xoff
    else:
        xoff = 0.
    yoff = 0.

    if reset:
        inpt.set_default('x1', x1)
        inpt.set_default('x2', x2)
        inpt.set_default('y1', y1)
        inpt.set_default('y2', y2)

    x1 = inpt.get_value('x1', 'left-hand limit of plot', x1)
    x2 = inpt.get_value('x2', 'right-hand limit of plot', x2)
    y1 = inpt.get_value('y1', 'bottom limit of plot', y1)
    y2 = inpt.get_value('y2', 'top limit of plot', y2)

    m_or_u    = inpt.get_value('mask', 'M(ask), U(nmask) or Q(uit)?', 'm', lvals=['m', 'M', 'u', 'U', 'q', 'Q'])
    if m_or_u.upper() == 'M':
        mtext = 'mask'
    else:
        mtext = 'unmask'

    mask_type = inpt.get_value('type', 'X, Y, P(hase), I(ndex) or Q(uit)?', 'x',
                               lvals=['x', 'X', 'y', 'Y', 'p', 'P', 'i', 'I', 'q', 'Q'])

    # initialise plot
    try:
        pg.pgopen(device)
        pg.pgsch(1.5)
        pg.pgscf(2)
        pg.pgslw(2)
        pg.pgsci(4)
        pg.pgenv(x1,x2,y1,y2,0,0)
        (xlabel,ylabel) = dset.plabel(xoff,yoff)
        pg.pgsci(2)
        pg.pglab(xlabel, ylabel, dset.title)

        # plot the dset
        dset.plot(xoff,yoff)

        x = (x1+x2)/2.
        y = (y1+y2)/2.

        # now define masks
        ch = 'X'
        while ch.upper() != 'Q':

            # go through mask options
            if mask_type.upper() == 'X':

                print('Set cursor at the one end of the X range, Q to quit')
                (xm1,y,ch) = pg.pgband(7,0,x,y)
                if ch.upper() != 'Q':
                    print('Set cursor at the other end of ' +
                          'the X range, Q to quit')
                    xm2,y,ch = pg.pgband(7,0,xm1,y)
                    if ch.upper() != 'Q':
                        if xm1 > xm2: xm1,xm2 = xm2,xm1
                        umask = mask.Xmask(xoff+xm1, xoff+xm2, m_or_u.upper() == 'M')

            elif mask_type.upper() == 'I':

                print('Place cursor near a point and click to ' +
                      mtext + ' it, Q to quit')
                x,y,ch = pg.pgband(7,0,x,y)
                if ch.upper() != 'Q':
                    xmm1,xmm2,ymm1,ymm2 = pg.pgqvp(2)
                    xscale  = (xmm2-xmm1)/(x2-x1)
                    yscale  = (ymm2-ymm1)/(y2-y1)

                    # only consider good data of opposite 'polarity' to the
                    # change we are making.
                    ok  = (dset.good == True) & \
                          (dset.mask == (m_or_u.upper() == 'M'))

                    if len(dset.x.dat[ok == True]):
                        # compute physical squared distance of cursor from
                        # points
                        sqdist  = npy.power(
                            xscale*(dset.x.dat[ok]-(xoff+x)),2) + \
                            npy.power(yscale*(dset.y.dat[ok]-(yoff+y)),2)

                        # select the index giving the minimum distance
                        indices = npy.arange(len(dset))[ok]
                        index   = indices[sqdist.min() == sqdist][0]
                        umask   = mask.Imask(index, m_or_u.upper() == 'M')
                    else:
                        print('There seem to be no data to ' + mtext +
                              '; data already ' + mtext + 'ed are ignored.')
                        umask = None

            if ch.upper() != 'Q' and umask is not None:
                gmask.append(umask)
                umask.app_mask(dset)

                print('overplotting data')
                # over-plot the dset
                dset.plot(xoff,yoff)

        pg.pgclos()
    except pg.ioerror, err:
        raise DintError(str(err))

    # Save final mask to disk
    if gmask != None:
        fout = open(mfile,'wb')
        pickle.dump(gmask, fout, -1)
        fout.close()
        print('Mask saved to ' + str(mfile))

def telescope(command, data, group, cdata):
    """
    Adds telescope information to slots.

    Interactive usage::

      telescope slots telescope observatory longitude latitude height

    Arguments::

      slots : string
           slot, slot range or group which will be processed.

      telescope : string
           telescope name. Enter ?? for by-number interactive
           selection from pre-stored data

    observatory : string
           observatory name

    longitude : float
           geodetic longitude, degrees.

    latitude : float
           geodetic latitude, degrees

    height : float
           height above sea level, metres
    """

    # generate arguments
    inpt = inp.Input(DINT_ENV, DINT_DEF, inp.clist(command))

    # register parameters
    inpt.register('slots',       inp.Input.LOCAL,  inp.Input.PROMPT)
    inpt.register('telescope',   inp.Input.LOCAL,  inp.Input.PROMPT)
    inpt.register('observatory', inp.Input.LOCAL,  inp.Input.PROMPT)
    inpt.register('longitude',   inp.Input.LOCAL,  inp.Input.PROMPT)
    inpt.register('latitude',    inp.Input.LOCAL,  inp.Input.PROMPT)
    inpt.register('height',      inp.Input.LOCAL,  inp.Input.PROMPT)

    # get inputs
    slots       = inpt.get_value('slots', 'slots to mask', '1-1000')
    slist       = interp_slots(slots, True, data, group)
    telescope   = inpt.get_value('telescope',
                                 'name of telescope, ?? for a list', '??')
    if telescope == '??':
        try:
            telescope,observatory,longitude,latitude,height = subs.observatory()
        except subs.SubsError, err:
            raise DintError('no observatory data set, ' + str(err))
        inpt.set_default('telescope',   telescope)
        inpt.set_default('observatory', observatory)
        inpt.set_default('longitude',   longitude)
        inpt.set_default('latitude',    latitude)
        inpt.set_default('height',      height)
    else:
        observatory = inpt.get_value('observatory',  'name of observatory',
                                     'La Palma')
        longitude   = inpt.get_value('longitude',
                                     'longitude (degrees, east positive)', 17.0)
        latitude    = inpt.get_value('latitude', 'latitude (degrees)', 28.0)
        height      = inpt.get_value('height', 'height (metres)', 2400.)

    for slot in slist:
        data[slot].settel(telescope, observatory, longitude, latitude, height)
        if not cdata.mute:
            print('Added telescope data to slot ' + str(slot))

def write(command, data, group, cdata):
    """
    Writes dsets to disk

    The data are written using Python's 'cPickle' module to write out
    in binary. If you need long term storage of your data, be aware that
    future changes to the Dset class could render these files unreadable.
    You would be better off writing out in FITS in this case.

    Interactive usage:

    write dfile slots append

    Arguments:

    dfile   -- name of disk file to write to.
    slots   -- range of slots as in '1-10', or just a single slot '9', or the name of a slot group
    append  -- do you want to append the slots to dfile? (careful: this will simply append to *any* file)
    """

    # generate arguments
    inpt = inp.Input(DINT_ENV, DINT_DEF, inp.clist(command))

    # register parameters
    inpt.register('dfile',  inp.Input.LOCAL, inp.Input.PROMPT)
    inpt.register('slots',  inp.Input.LOCAL, inp.Input.PROMPT)
    inpt.register('append', inp.Input.LOCAL, inp.Input.PROMPT)

    # get inputs
    dfile  = inpt.get_value('dfile',  'file to write the data to', subs.Fname('example', '.dnl', subs.Fname.NEW))
    slots  = inpt.get_value('slots', 'slots to load dsets into', '1-1000')
    slist  = interp_slots(slots, True, data)
    append = inpt.get_value('append', 'append to an old file?', False)
    if append and not dfile.exists():
        raise DintError('File = ' + str(dfile) + ' does not exist.')

    if append:
        optr = open(dfile,'ab')
    else:
        optr = open(dfile,'wb')

    # write out data
    for slot in slist:
        data[slot].write(optr)
        if not cdata.mute:
            print('Wrote slot',slot,'into file ' + str(dfile))
    optr.close()
