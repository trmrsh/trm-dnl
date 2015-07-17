#!/usr/bin/env python

"""
danal is an interactive interface to the dnl routines. It enables a simple
command-line interface for handling multiple datasets (called 'Dsets' in the
documentation) while still allowing access to Python's full command set. For
help on Dsets read up on pydoc trm.dnl.Dset; for help on the interactive
commands see pydoc trm.dnl.dint.

Example usage:

load data.dnl 1-10   [loads up to 10 dsets from the file data.dnl
plot 2-3             [plot numbers 2 and 3]
quit                 [quit the program]

'readline' is used to provide a command history and also to enable tab
completion of file names (in the directory the program is running in).
Values of command parameters are also stored in disk files to give the
program memory
"""
from __future__ import division, absolute_import, print_function

#from IPython.Shell import IPShellEmbed
#ipshell = IPShellEmbed(args,
#                       banner = 'Dropping into IPython',
#                       exit_msg = 'Leaving Interpreter, back to program.')

import readline
import os
import sys
from trm.dnl import dint

# Name of environment variable for default storage
dint.DINT_ENV = 'PYTHON_DINT_ENV'

# Name of default directory if DINT_ENV not set
dint.DINT_DEF = '.pydint'

# Initialise command history stuff from previously
# stored file if possible.
if dint.DINT_ENV in os.environ:
    ddir = os.environ[dint.DINT_ENV]
else:
    home = os.environ['HOME']
    ddir = os.path.join(home, dint.DINT_DEF)

if not os.path.lexists(ddir):
    os.mkdir(ddir, 0755)

if not os.path.isdir(ddir):
    raise dint.DintError('Default directory = ' + ddir + ' is not a directory')

chist = os.path.join(ddir, '.history')
if os.path.isfile(chist):
    readline.read_history_file(chist)

# Initialise data, group and cdata
data,group,cdata = dint.init()

print("""
Welcome to danal for interactive data analysis.
Type "help" for a list of commands; ctrl-D to exit.
""")

# Command loop
while 1:

    try:

        command = raw_input('danal> ')
        if command != '':
            dint.execute(data, group, cdata, command)

    except EOFError:

        print('\n')
        reply = raw_input('Do you want to quit? ')
        if reply.lower() == 'y':
            # save command history
            readline.write_history_file(chist)
            exit(0)

    except KeyboardInterrupt:
        print('\nctrl-C intercepted; if you want to quit, hit ctrl-D')
