#!/usr/bin/env python


from __future__ import division


import sys
from collections import defaultdict
import re

import opster


GROUP_NAMES = ['controls', 'patients', 'relatives']


def group_name(key):
    """Find group name from key"""
    for g in GROUP_NAMES:
        if key.startswith(g):
            return g
    else:
        raise ValueError("Unrocgnised group in {0}".format(key))

def abs_key(key, bands=False):
    """Get canonical key from raw key"""
    if not bands:
        return key
    else:
        return re.sub(r'_\d\.edf', '', key)

def read_results(fin, bands=False):
    """Compile results from stdin"""
    # Data structure for results
    empty = lambda : (0, 0)
    groups = dict((g, defaultdict(empty)) for g in GROUP_NAMES)

    for line in fin:
        key, T_new, n_new = line.split(', ')
        T_new = float(T_new)
        n_new = float(n_new)
        g = group_name(key)
        k = abs_key(key, bands)
        T_old, n_old = groups[g][k]
        groups[g][k] = (T_old + T_new, n_old + n_new)

    return groups

def write_results(groups, fout):
    """Write data in r format to fout"""
    # Print out info for each key
    for g, d in groups.iteritems():
        for key in sorted(d.iterkeys()):
            T, n = d[key]
            print >> fout, '{0}, {1}, {2}'.format(key, T, n)

def write_stats(groups, fout, lower=False):
    """Write data in r format to fout"""
    # Print data in r format
    for k, d in groups.iteritems():
        rates = [3600 * (n / T) for T, n in d.itervalues()]
        print >> fout, '{0} <- c({1})'.format(k, ' ,'.join(str(r) for r in rates))
        print >> fout, '#', k, (sum(rates) / len(rates)) if rates else 0

    # Hypothesis
    if not lower:
        alt = 'greater'
    else:
        alt = 'less'

    # Add Wilcox tests
    testfmt = 'wilcox.test({0}, {1}, alternative="{2}", exact=FALSE)'
    for k1 in 'patients', 'relatives':
        for k2 in 'relatives', 'controls':
            if k1 != k2:
                print >> fout, testfmt.format(k1, k2, alt)
    # Add fisher exact tests
    fmt = 'contingency <- matrix(c({0}), nr=3, dimnames=list(c({1}), c({2})))'
    mat = []
    NOs, SEs = [], []
    for g in 'patients', 'relatives', 'controls':
        d = groups[g]
        NOs.append(sum(1 for T, n in d.itervalues() if not n))
        SEs.append(sum(1 for T, n in d.itervalues() if n))
    print >> fout, fmt.format(', '.join(str(n) for n in NOs + SEs),
                     '"patients", "relatives", "controls"',
                     '"NO", "SE"')
    print >> fout, "contingency"
    print >> fout, "fisher.test(contingency)"

    # Add means
    for g in 'patients', 'relatives', 'controls':
        print >> fout, 'mean({0})'.format(g)

@opster.command(usage='%name [INFILE [OUTFILE]]')
def main(infile=None,
         outfile=None,
         stats=('s', False, 'Show r format output'),
         bands=('b', False, 'Collect different frequency bands'),
         lower=('', False, 'Show statistically significant higher results')):

    # Open input/output files
    fin = open(infile, 'r') if infile else sys.stdin
    fout = open(outfile, 'w') if outfile else sys.stdout

    # Read input data
    groups = read_results(fin, bands=bands)

    # Write output data
    if stats:
        write_stats(groups, fout, lower)
    else:
        write_results(groups, fout)


if __name__ == "__main__":
    main(argv=sys.argv[1:])
