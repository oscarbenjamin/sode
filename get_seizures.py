#!/usr/bin/env python


from __future__ import division


import sys
import os.path
import gc
from cPickle import Pickler, Unpickler
import operator

import numpy
import opster


THRESHOLD = 1
MINTIME = 30


def at_least_n(times, nmin):
    """Rrtuen sorted list of times that appear nmin times"""
    # Special case the empty list
    if not len(times):
        return []

    # Sort so counts can be performed sequentially
    times = sorted(times)

    # Iterate through counting occurrences
    t_nmin = []
    count = 1
    tlast = times[0]
    for tcurrent in times[1:]:
        if tcurrent == tlast:
            count += 1
        else:
            if count >= nmin:
                t_nmin.append(tlast)
            tlast = tcurrent
            count = 1
    if count >= nmin:
        t_nmin.append(tcurrent)

    # Return list in which each entry appears once
    return t_nmin

def find_seizures(t, xs, threshold=THRESHOLD, mintime=MINTIME, nmin=None, **opts):
    """Identify times corresponding to start of seizures"""
    # Default is to use half the number of signals
    if not nmin:
        nmin = len(xs) // 2
    assert nmin

    # Identify all times where a signal exceeds threshold
    t_thresh = numpy.concatenate([t[x > threshold] for x in xs])

    # Only those where nmin signals exceed simultaneously
    t_nmin = at_least_n(t_thresh, nmin)

    # Keep only those times, where at least mintime has passed
    ts = []
    tlast = t[0] - 2 * mintime
    for t in t_nmin:
        if t > tlast + mintime:
            ts.append(t)
            tlast = t

    # Return onset times
    return ts

def find_seizures_file(fname, threshold, mintime, nmin, **opts):
    # Open file and load signal data
    from electrolib import EEGFile
    eeg = EEGFile(fname)
    sigs = eeg.get_signals(eeg.tstart, eeg.tend)
    t = sigs[0].t
    T = t[-1] - t[0]
    xs = [s.x for s in sigs]

    # Find seizures
    times = find_seizures(t, xs, threshold, mintime, nmin)

    # Return seizure times and total Time
    return times, T

def _cache_name_key(fname, **opts):
    # Function shared by load and save cache functions below
    directory, key = os.path.split(fname)
    cache = 'thresh_{threshold}_tmin_{mintime}_nmin_{nmin}.dat'.format(**opts)
    cache_path = os.path.abspath(os.path.join(directory, cache))
    return cache_path, key

__CACHES = {}
__CACHE_COUNT = 0
def _load_cache(fname, **opts):
    # Fail silently
    cache_path, key = _cache_name_key(fname, **opts)
    if cache_path not in __CACHES:
        try:
            data = Unpickler(open(cache_path, 'rb')).load()
            data = dict(data)
        except:
            data = {}
        __CACHES[cache_path] = data
    data = __CACHES[cache_path]
    return data, key

def _save_caches(**opts):
    # Fail silently
    if not opts['cache']:
        return
    for cache_path, data in __CACHES.iteritems():
        try:
            data = sorted(data.items(), key=operator.itemgetter(0))
            Pickler(open(cache_path, 'wb')).dump(data)
        except:
            pass

def _paths_dir(d, **opts):
    # Find keys related to dir
    path = os.path.join(d, '_dummy_')
    cache_path, _ = _cache_name_key(path, **opts)
    _load_cache(path, **opts)
    keys_cached = __CACHES[cache_path].keys()
    keys_files = [p for p in os.listdir(d) if p.endswith('.edf')]
    keys = sorted(set(keys_cached + keys_files))
    paths = [os.path.join(d, k) for k in keys]
    return paths

def get_seizure_count(fname, **opts):
    """Get the seizure count (from cache if available)"""
    # Attempt to use cache
    data, key = _load_cache(fname, **opts)
    if key not in data or opts['force']:
        # Actually compute values
        data[key] = find_seizures_file(fname, **opts)
        # Save if necessary
        global __CACHE_COUNT
        __CACHE_COUNT += 1
        # Save every 10th new entry
        if not __CACHE_COUNT % 10:
            _save_caches(**opts)

    # Just return the value
    return data[key]

def plot_eeg(ax, eeg, scale=1, every=1, mark_times=None, **opts):
    # Retrieve appropiate data
    sigs = eeg.get_signals(eeg.tstart, eeg.tend)

    t = sigs[0].t[::every]
    T = t[-1] - t[0]

    # Adjust for units
    mark_times = numpy.array(mark_times, float)
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

    # Print seizure times
    for tmark in mark_times:
        print "Marking {0} {1}".format(tmark, tunits)

    # Plot each line one by one
    ticklabels = []
    for s in sigs:
        if s.label.endswith('y') or s.label=='lam':
            continue
        offset = len(ticklabels)
        x = offset + s.x[::every] / float(scale)
        ax.plot(t, x, 'k-')
        ticklabels.append(s.label)

    # Mark horizontal lines
    ymin = -1
    ymax = len(ticklabels)
    if mark_times is not None:
        for tmark in mark_times:
            ax.plot([tmark, tmark], [ymin, ymax], 'r-')

    # Sort out formatting
    ax.set_ylim([ymin, ymax])
    ax.set_yticks(range(len(ticklabels)))
    ax.set_yticklabels(ticklabels)
    ax.set_xlabel('Time ({0})'.format(tunits))


opts = [
    ('', 'threshold', 1, 'Threshold scale for seizure'),
    ('', 'nmin', 0, 'Minimum number of chans for seizure'),
    ('', 'mintime', 30, 'Minimum time between seizures'),
    ('s', 'scale', 1, '(with --plot) Approximate scale of signals'),
    ('e', 'every', 10, '(with ---plot) Use 1 in every N time samples'),
    ('p', 'plot', False, 'Plot to screen'),
    ('f', 'force', False, 'Force recomputation (rather than reading cache)'),
    ('c', 'cache', False, 'Save caches'),
]
@opster.command(usage='%name [OPTS] EEGFILE [EPSFILE]', options=opts)
def main(*args, **opts):
    """Compute and display spectrograms"""
    # Stats counter
    groups = {'controls':[], 'patients':[], 'relatives':[]}

    # Append files and edf files from directories
    fnames = []
    for a in args:
        if os.path.isdir(a):
            fnames.extend(_paths_dir(a, **opts))
        else:
            fnames.append(a)

    # loop through files
    for fname in fnames:
        # Find seizures
        try:
            times, T = get_seizure_count(fname, **opts)
        except EOFError:
            continue

        # Accumulate stats
        bname = os.path.basename(fname)
        print '{0}, {1}, {2}'.format(bname, T, len(times))

        # Do Plot
        if opts['plot']:
            # Output for info
            print "Plotting:", fname
            import pylab
            fig = pylab.figure()
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
            plot_eeg(ax, eeg, mark_times=times, **opts)
            pylab.show()

    # Save all data
    _save_caches(**opts)

    return 0


if __name__ == "__main__":
    main(argv=sys.argv[1:])


def test_find_seizures_nothing():
    """Test the find seizures function with zeros"""
    import numpy
    t = numpy.arange(0, 1000, 0.1)
    nc = 10
    xs = [numpy.zeros_like(t) for n in range(nc)]
    assert find_seizures(t, xs) == []

def test_find_seizures_nonoise():
    """Test the find seizures function with zeros"""
    import numpy
    t = numpy.arange(0, 1000, 0.1)
    nc = 10
    xs = [numpy.zeros_like(t) for n in range(nc)]
    ns = len(t)
    samples = [ns // 4, ns // 2, 3 * ns // 4]
    Ts = t[samples]
    dT = 10
    for x in xs:
        for T in Ts:
            x[(t >= T) & (t < T + dT)] = 2
    assert all(find_seizures(t, xs) == Ts)
    assert find_seizures(t, xs, threshold=3) == []

def test_find_seizures_smallnoise():
    """Test the find seizures function with zeros"""
    import numpy
    from numpy.random import randn
    t = numpy.arange(0, 1000, 0.1)
    nc = 10
    xs = [0.1 * randn(*t.shape) for n in range(nc)]
    ns = len(t)
    samples = [ns // 4, ns // 2, 3 * ns // 4]
    Ts = t[samples]
    dT = 10
    for x in xs:
        for T in Ts:
            x[(t >= T) & (t < T + dT)] += 2
    assert all(find_seizures(t, xs) == Ts)
    assert find_seizures(t, xs, threshold=3) == []

def test_find_seizures_roguechan():
    """Test the find seizures function with zeros"""
    import numpy
    from numpy.random import randn
    t = numpy.arange(0, 1000, 0.1)
    nc = 10
    xs = [numpy.zeros_like(t) for n in range(nc)]
    ns = len(t)
    samples = [ns // 4, ns // 2, 3 * ns // 4]
    Ts = t[samples]
    dT = 10
    x = xs[0]
    for T in Ts:
        x[(t >= T) & (t < T + dT)] += 2
    assert find_seizures(t, xs) == []

