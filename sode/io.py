# (c) Oscar Benjamin, 2011, under terms of the new BSD license
#
# Author: Oscar Benjamin
# Date: 08 Mar 2011
#
# This module defines the SODE class for numerical solution of Stochastic
# Ordinary Differential Equations with additive noise.

from __future__ import division

import datetime
import re
import itertools
import sys

import numpy as np


def load_csv(fin, parameters=False):
    """Load solution t, Xt from input-file fin

    Usage:
    >>> # Create a system, generate a solution and save it
    >>> from sode.weiners import Weiner
    >>> from sode.algos import solve
    >>> from sode.io import save_csv, load_csv
    >>> system = Weiner()
    >>> x1 = system.get_x0()
    >>> t = np.arange(0, 1, 1e-2)
    >>> Xt = solve(system, x1, t, dtmax=1e-3)
    >>> save_csv(system, t, Xt, open('sys1_1.csv', 'w'))
    >>> # Create a new system
    >>> system = Weiner()
    >>> t, Xt = load_csv(open('sys1_1.csv', 'r'))
    >>> # Also read parameter values
    >>> t, Xt, ret = load_csv(open('sys1_1.csv', 'r'), parameters=True)
    >>> print ret
    ({}, ['mu=0.0', 'sigma=1.0', 'x0=0.0'])

    If parameters is True the file will be checked for the appropriate
    header information. It will be checked for consistency against this
    SODE instance. In particular it must have been saved by an instance of
    the same SODE subclass. Parameters will also be loaded from the file
    into this instance so that it has the same parameter values as used by
    the instance that saved the solution.

    Otherwise, the file is treated as a csv file where lines beginning
    with '#' are comments. The first non-comment line contains the
    (comma-separated) column titles, starting with 'time' and then each of
    the variables in order.  Subsequent lines contain comma-seperated
    values, e.g.:

    # COMMENTS GO HERE
    # COMMENTS HERE
    time, x, y
    0.0, 1.0, 0.0
    0.1, 0.9976547, 0.0132455
    """
    # Doesn't accept string input unless split into lines
    if isinstance(fin, str):
        fin = open(fin, 'r')

    # Use generator to remove newlines if present
    def iter_f():
        for line in fin:
            if line.endswith('\n'):
                yield line[:-1]
            else:
                yield line

    # Create iterator and draw first line
    itf = iter_f()
    line = itf.next()

    # Ignore comment lines
    comment_lines = []
    while line.startswith('#'):
        comment_lines.append(line)
        line = itf.next()

    # Parse variable names
    names = line.split(', ')[1:]
    nvars = len(names)

    # Read parameter values from comment lines
    if parameters:
        ret = load_parameters(comment_lines)

    # Read the rest of the file
    lines = list(itf)

    # Iterate through parsing into arrays
    t = np.zeros((len(lines)))
    Xt = np.zeros((len(lines), nvars))
    for n, l in enumerate(lines):
        nums = [float(s) for s in l.split(', ')]
        t[n], Xt[n, :] = nums[0], nums[1:]

    # And return the result
    if parameters:
        return t, Xt, ret
    else:
        return t, Xt


def save_csv(sodeinst, t, Xt, fout, header=True, titles=True):
    """Save solution t, Xt to output-file fout

    Usage:
    >>> from sode.weiners import Weiner
    >>> from sode.algos import solve
    >>> from sode.io import save_csv
    >>> system = Weiner()
    >>> x1 = system.get_x0()
    >>> t = np.arange(0, 1, 1e-2)
    >>> Xt = solve(system, x1, t, dtmax=1e-3)
    >>> save_csv(system, t, Xt, open('sys1_1.csv', 'w'))

    Writes the solution to fout in csv format where lines beginning with
    '#' are comments. Lines beginning with '##' are used to write
    meta-data used by load_csv. Example:

    ### SODE VERSION 0
    ###.{}
    ###.[]
    ## 1-D linear SODE with equation:
    ##
    ##         dx = - beta x dt + alpha dW
    ##
    ## par beta = 1.0
    ## par alpha = 1.0
    ## var x = 1.0
    # Created: 2011-03-12 18:42:32.040680
    time, x
    0.0, 1.0
    0.1, 1.23475373685
    """
    # First write in comment headers
    if header:
        # Magic line
        fout.write('### SODE VERSION 0\n')

        # _make_sode in Script adds these attributes
        if not hasattr(sodeinst, '_sys_opts'):
            sys_opts, sys_args = {}, []
        else:
            sys_opts, sys_args = sodeinst._sys_opts
        fout.write('##.{0}\n'.format(repr(sys_opts)))
        fout.write('##.{0}\n'.format(repr(sys_args)))

        # Write description
        for line in sodeinst.get_description().split('\n'):
            fout.write('## {0}\n'.format(line))
        dt = datetime.datetime.now()
        fout.write('# Created: {0} by {1}\n'.format(dt, sys.argv[0]))

    if titles:
        # Write column header
        cols = ['time'] + sodeinst.get_variables()
        fout.write(', '.join(cols) + '\n')

    # write raw numbers (csv format)
    for ti, xi in itertools.izip(t, Xt):
        nums = [str(ti)] + [str(x) for x in xi]
        fout.write(', '.join(nums) + '\n')

# Regexes for parsing solution file
_MAGIC_LINE = '^### SODE .*$'
_PV_LINE = '^## ([pv]ar) (\w+) = ([0-9.-]+)'

def load_parameters(fin):
    """Read parameter header information from saved file.

    The information should be formatted as by save_csv"""
    # Create line iterator
    if isinstance(fin, str):
        fin = open(fin, 'r')
    itf = iter(fin)
    line = itf.next()

    # The data we want
    sys_opts = {}
    sys_args = []

    # Check first line
    if not re.match(_MAGIC_LINE, line):
        raise ValueError("Invalid solution file for parameters")
    line = itf.next()

    # See if we have sys_opts
    if line.startswith('##.'):
        sys_opts = eval(line[3:])
        line = itf.next()

    # Attempt to read header
    while line.startswith('##'):
        m = re.match(_PV_LINE, line)
        if m:
            pv, name, val = m.groups()
            val = float(val)
            if pv == 'par':
                a = "{0}={1}".format(name, val)
            elif pv == 'var':
                a = "{0}0={1}".format(name, val)
            sys_args.append(a)
        line = itf.next()

    # Return the parse data
    return sys_opts, sys_args


# Quick test suite to demonstrate usage of sode.py and script.py
if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import os
    os.remove('sys1_1.csv')
