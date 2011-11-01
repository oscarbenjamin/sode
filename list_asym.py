#!/usr/bin/env python


from __future__ import division


import sys
import os.path

import opster
import numpy as np

from sode.graphviz import GraphParser
from mygraph import Digraph


CHANS = [
    'Fp1',#0
    'Fp2',#1
    'F3' ,#2
    'F4' ,#3
    'C3' ,#4
    'C4' ,#5
    'P3' ,#6
    'P4' ,#7
    'O1' ,#8
    'O2' ,#9
    'F7' ,#10
    'F8' ,#11
    'T3' ,#12
    'T4' ,#13
    'T5' ,#14
    'T6' ,#15
    'Fz' ,#16
    'Cz' ,#17
    'Pz' ,#18
]
LRPAIRS = [
    ('Fp1', 'Fp2'),
    ('F3' , 'F4' ),
    ('C3' , 'C4' ),
    ('P3' , 'P4' ),
    ('F3' , 'F4' ),
    ('O1' , 'O2' ),
    ('F7' , 'F8' ),
    ('T3' , 'T4' ),
    ('T5' , 'T6' ),
]
LRPAIRS_SUB = {
    'all':LRPAIRS,
    'front':[
        ('Fp1', 'Fp2'),
        ('F3' , 'F4' ),
        ('F7' , 'F8' ),
     ],
    'back':[
        ('P3' , 'P4' ),
        ('O1' , 'O2' ),
        ('T5' , 'T6' ),
    ],
    'occ':[
        ('O1' , 'O2' ),
    ],
    'parietal':[
        ('P3' , 'P4' ),
    ],
    'temporal':[
        ('T3' , 'T4' ),
        ('T5' , 'T6' ),
    ],
}


def find_in_directory(a):
    """Find all dot files within directory"""
    return [os.path.join(a, f) for f in os.listdir(a) if f.endswith('.dot')]

def load_graph(p):
    """Obtain adjacency matrix from dot-file name"""
    return Digraph.from_dotfile(p)

# This is the function that measures input/output imbalance
def imbalance(g):
    um = g.get_indegree() - g.get_outdegree()
    return sum(um ** 2)

# This is the function that measures input/output imbalance
def imbalance_mod(g):
    um = g.get_indegree() - g.get_outdegree()
    return sum(abs(um))

# Determine x-value for plot
def graph_factor(g):
    g1st = g.first_transitive_component()
    imfactor = imbalance(g1st)
    Ne = len(g1st.edges())
    if Ne:
        imfactor /= Ne
    return len(g1st.nodes()) - imfactor

# This is the function that measures input/output imbalance
def graph_factor_lr(g, pairs=LRPAIRS):
    nodes = g.nodes()
    total = 0
    for l, r in pairs:
        nl = nodes[CHANS.index(l)]
        nr = nodes[CHANS.index(r)]
        if g.has_edge((nl, nr)) != g.has_edge((nr, nl)):
            total += 1
    return total

# Count the number of edges between pairs of left right electrodes
def graph_factor_lre(g, pairs=LRPAIRS):
    nodes = g.nodes()
    total = 0
    for l, r in pairs:
        nl = nodes[CHANS.index(l)]
        nr = nodes[CHANS.index(r)]
        total += g.has_edge((nl, nr))
        total += g.has_edge((nr, nl))
    return total

# This is the function that measures input/output imbalance
def graph_factor_ws(g):
    return int(g.is_weakly_connected())

opts = [
    ('', 'lr', False, 'Do a left-right comparison'),
    ('', 'lre', False, 'Do a left-right edges comparison'),
    ('', 'ws', False, 'Do a weakly-strongly connected comparison'),
    ('', 'imb', False, 'Compare imbalance'),
    ('', 'mod', False, 'Compare mod imbalance'),
    ('', 'pairs', 'all', '"all", "front", "back"'),
]
@opster.command(usage='%name [OPTS] EEGFILE [EPSFILE]', options=opts)
def main(*args, **opts):
    """Compute and display spectrograms"""
    # Stats counter
    groups = {'controls':[], 'patients':[], 'relatives':[]}

    pairs = LRPAIRS_SUB[opts['pairs']]

    # Append files and edf files from directories
    fnames = []
    for a in args:
        if os.path.isdir(a):
            fnames.extend(find_in_directory(a))
        else:
            fnames.append(a)

    # loop through files
    for fname in fnames:
        # Accumulate stats
        bname = os.path.basename(fname)
        g = load_graph(fname)
        if opts['lr']:
            factor = graph_factor_lr(g, pairs)
        elif opts['lre']:
            factor = graph_factor_lre(g, pairs)
        elif opts['ws']:
            factor = graph_factor_ws(g)
        elif opts['mod']:
            factor = imbalance_mod(g)
        elif opts['imb']:
            factor = imbalance(g)
        else:
            factor = graph_factor(g)
        print '{0}, {1}, {2}'.format(bname, 3600, factor)

    return 0


if __name__ == "__main__":
    main(argv=sys.argv[1:])


