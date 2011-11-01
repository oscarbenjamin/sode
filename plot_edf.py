#!/usr/bin/env python

from __future__ import division


import sys
import os.path
import gc

import numpy
import pylab
import opster

from electrolib import EEGFile

from get_seizures import find_seizures


def plot_eeg(ax, eeg, scale=2000, every=1, mark_times=[]):
    """Plot the EEG file eeg, maring times if appropriate"""

    # Retrieve appropiate data
    sigs = eeg.get_signals(eeg.tstart, eeg.tend)
    t = sigs[0].t[::every]
    xs = [s.x[::every] for s in sigs]
    labels = [s.label for s in sigs]

    # Adjust for units
    mark_times = numpy.array(mark_times, float)
    T = t[-1] - t[0]
    if T / 3600 > 3:
        t /= 3600
        mark_times /= 3600
        tunits = 'hours'
    elif T / 60 > 3:
        t /= 60
        mark_times /= 60
        tunits = 'minutes'
    else:
        tunits = 'seconds'

    # Plot each line one by one
    for n, x in enumerate(xs):
        ax.plot(t, n + x / scale, 'k-')

    # Mark horizontal lines
    ymin, ymax = -1, len(labels)
    for tmark in mark_times:
        ax.plot([tmark, tmark], [ymin, ymax], 'r-', linewidth=2)

    # Sort out formatting
    ax.set_xlim([t[0], t[-1]])
    ax.set_ylim([ymin, ymax])
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_xlabel('Time ({0})'.format(tunits))


opts = [
    ('s', 'scale', 4, 'Approximate scale of signals'),
    ('m', 'mark', True, 'Mark seizures'),
    ('', 'threshold', 1, 'Threshold scale for seizure'),
    ('', 'nmin', 0, 'Minimum number of chans for seizure'),
    ('', 'tmin', 30, 'Minimum time between seizures'),
    ('e', 'every', 1, 'Use 1 in every N time samples'),
    ('O', 'output', False, 'Plot to file using original filename and extension'),
    ('', 'over-write', False, 'Overwrite existing plot files'),
]
@opster.command(usage='%name [OPTS] EEGFILE [EPSFILE]', options=opts)
def main(*args, **opts):
    """Compute and display spectrograms"""
    # Re-use 1 figure instance
    fig = pylab.figure()

    fnames = []
    for a in args:
        if os.path.isdir(a):
            edffiles = [f for f in os.listdir(a) if f.endswith('.edf')]
            fnames.extend(os.path.join(a, f) for f in edffiles)
        else:
            fnames.append(a)

    # loop through files
    for fname in fnames:

        # Clean up
        fig.clear()
        gc.collect()

        # Check if plot file exists already
        if opts['output']:
            epsname = os.path.splitext(fname)[0] + '.eps'
            if not opts['over_write'] and os.path.exists(epsname):
                print "Skipping:", fname
                continue

        # Output for info
        print "Plotting:", fname

        # Open file and find seizures
        eeg = EEGFile(fname)
        if opts['mark']:
            sigs = eeg.get_signals(eeg.tstart, eeg.tend)
            t = sigs[0].t
            xs = [s.x for s in sigs]
            times = find_seizures(t, xs, thresh=opts['threshold'],
                                  nmin=opts['nmin'], tmin=opts['tmin'])
        else:
            times = []

        # Create axes, open EEG file and plot
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        plot_eeg(ax, eeg, mark_times=times, scale=opts['scale'],
                 every=opts['every'])

        # Save to file, or plot to screen
        if opts['output']:
            print "Saving to:", epsname
            fig.savefig(epsname)
        else:
            pylab.show()

    return 0


if __name__ == "__main__":
    main(argv=sys.argv[1:])


