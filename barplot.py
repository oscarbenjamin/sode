#!/usr/bin/env python

from __future__ import division

import sys

import opster
import numpy as np

COLOURS = ['#348ABD',
           '#7A68A6',
           '#A60628',
           '#467821',
           '#CF4457',
           '#188487',
           '#E24A33']
COLOURS = [
    '#A60628',
    '#7A68A6',
    '#348ABD',
]
GROUPS = ('patients', 'relatives', 'controls')

def read_data(infile):
    """Read the data from the input file"""
    data = {}

    for line in infile:
        # Remove newline
        line = line[:-1]

        # Time to start a new data entry
        if line.startswith('##'):
            group_key = line[3:]
            data[group_key] = {}
            dtmp = {}
        elif line.startswith('#'):
            key = line[2:]
            data[group_key][key] = dtmp
            dtmp = {}
        else:
            # Split into name, expr pairs
            name, expr = line.split(':')
            if name == 'degree':
                dtmp['degree'] = int(expr)
            elif name == 'rates':
                dtmp['rates'] = list(eval(expr))
            elif name == 'pvals':
                dtmp['pvals'] = list(eval(expr))

    return data

# Add stars to barplot
def add_stars(ax, rects, pvals):
    for p, r in zip(pvals, rects):
        if p > 0.05:
            continue
        elif 0.01 < p <= 0.05:
            text = r'$\ast$'
        elif 0.001 < p <= 0.01:
            text = r'$\ast\!\ast$'
        elif 0.0001 < p <= 0.001:
            text = r'$\ast\!\ast\!\ast$'
        elif 0.00001 < p <= 0.0001:
            text = r'$\ast\!\ast\!\ast\!\ast$'
        y = r.get_height() + .15
        x = r.get_x() + r.get_width() / 2
        ax.text(x, y, text, ha='center', va='center', size='x-small')

def plot_bars(ax, title, d, legend=False, ylabel=False, xlabels=False,
        groups=GROUPS):
    # patients, relatives, controls
    N = len(groups)
    width = 1 / (N + 1)
    xoffset = width / 2

    # Left sides of patients
    ind = np.arange(len(d))
    order = sorted(d)

    # Keep track of rects for legend
    rects = {}

    # Iterate over patients, relatives, controls
    for n, g in enumerate(groups):
        rates = [d[k]['rates'][GROUPS.index(g)] for k in order]
        lefts = ind + n * width + xoffset
        rs = ax.bar(lefts, rates, width, color=COLOURS[n])
        rects[g] = rs

    # Add significance stars
    pvals = {}
    pvals['patients'] = [d[k]['pvals'][1] for k in order]
    pvals['relatives'] = [d[k]['pvals'][2] for k in order]
    for g in groups:
        if g in pvals:
            add_stars(ax, rects[g], pvals[g])

    # Formatting options
    if legend:
        example_rects = [rects[g][0] for g in groups]
        ax.legend(example_rects, groups)
    if ylabel:
        ax.set_ylabel(r'Frequency ($\mathrm{h}^{-1}$)')
    if xlabels:
        ax.set_xticks(ind + (N/2) * width + xoffset)
        ax.set_xticklabels(order, ha='center')
    else:
        ax.set_xticks([])

    # Disable tick lines
    for t in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    #ax.set_ylim([0,3])
    #ax.set_yticks(range(4))
    ax.set_title(title.capitalize())

opts = [
    ('g', 'groups', '', 'Groups to include'),
    ('l', 'list-groups', False, 'List possible groups and exit'),
]
@opster.command(usage='%name [FILE]', options=opts)
def main(infile=None, epsfile=None, **opts):
    """Read data from infile and plot as bar chart"""
    # Process cl options and arguments
    if opts['list_groups']:
        print ' '.join(GROUPS)
        return 0
    if opts['groups']:
        groups = opts['groups'].split(',')
    else:
        groups = GROUPS
    # Default to stdin
    infile = open(infile) if infile else sys.stdin

    # Read input data
    data = read_data(infile)

    # Create figure
    from pylab import figure, show
    fig = figure(figsize=(4,6))

    # Position 1st (larger) axes
    left, width = 0.15, 0.70
    bottom1 = 0.75
    height1 = 0.2
    ax1 = fig.add_axes([left, bottom1, width, height1])

    # Position other axes
    N = 5
    gap = bottom1 / (N + .5)
    bottom, height = bottom1 - 0.02, 0.075
    axs = []
    for n in range(N):
        bottom -= gap
        ax = fig.add_axes([left, bottom, width, height])
        axs.append(ax)

    # Add barplots and format
    plot_bars(ax1, 'average', data['average'], legend=True, ylabel=True,
            xlabels=True, groups=groups)
    plot_bars(axs[0], r'$\delta$', data['delta'], groups=groups)
    plot_bars(axs[1], r'$\theta$', data['theta'], groups=groups)
    plot_bars(axs[2], r'$\alpha$', data['alpha'], groups=groups)
    plot_bars(axs[3], r'$\beta$' , data['beta'] , groups=groups)
    plot_bars(axs[4], r'$\gamma$', data['gamma'], groups=groups, xlabels=True)

    # Present to the screen
    if not epsfile:
        show()
    else:
        fig.savefig(epsfile)


if __name__ == "__main__":
    main(argv=sys.argv[1:])
