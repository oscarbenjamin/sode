#!/usr/bin/env python


from __future__ import division

import datetime
import warnings
import sys
import os.path
import gc

import numpy as np
import opster


# Windows fix
def _fix_mode(fd):
    try:
        import msvcrt, os
        msvcrt.setmode(1, os.O_BINARY)
    except ImportError:
        pass


def csv_to_edf(fin=None, all_chan=False):
    # Open files ready for action
    if fin is None:
        fin = sys.stdin
        fout = sys.stdout
        _fix_mode(fout)
    else:
        name, ext = os.path.splitext(fin)
        fin = open(fin, 'r')
        fout = name + '.edf'
        print fout
        fout = open(fout, 'wb')

    # Begin reading csv file
    itf = iter(fin)
    line = itf.next()

    # Skip comments
    while line.startswith('#'):
        line = itf.next()

    # Read titles
    cols = [s.strip() for s in line.split(',')]
    Ncols = len(cols)

    # Establish channel names
    names = cols[1:]

    # Choose the channels we want to keep
    if not all_chan:
        included = []
        indices = []
        for n, name in enumerate(names):
            if name.endswith('.x'):
                included.append(name[:-2])
                indices.append(n)
    else:
        included = names
        indices = list(range(len(names)))

    Nc = len(included)

    # Read a few lines into buffer
    lines_buffer = []
    for n in range(10):
        lines_buffer.append(itf.next())

    # Read first number of each line
    ts = [float(l.split(',')[0].strip()) for l in lines_buffer]
    dts = np.diff(ts)
    dt = dts.mean()
    err = max(abs(dts - dt)/dt)
    if err > 1e-5:
        warnings.warn("Expected fixed width data: '{0}'".format(err))

    # Attempt to use integer seconds in block
    Hz = 1 / dt
    if ((Hz % 1) / Hz) < 1e-4:
        nsblock = 10 * int(Hz)
        blockdur = 10
    else:
        nsblock = int(10 * Hz)
        blockdur = nsblock * dt
    nblocks = -1

    # These are good for the graph sims
    physmin = -2
    physmax = 2
    digmin = -32768
    digmax = 32767

    # Use now as the date/time
    tstart = datetime.datetime.now()
    nbytes = 256*(1 + Nc)

    # Construct header as string
    header = ''
    header +=  "0".ljust(8)
    header += "simulation".ljust(80)
    header += "graph_file".ljust(80)
    header += tstart.strftime("%d.%m.%y")
    header += tstart.strftime("%H.%M.%S")
    header += str(nbytes).ljust(8)
    header += "".ljust(44)
    header += str(nblocks).ljust(8)
    header += str(blockdur).ljust(8)[:8]
    header += str(Nc).ljust(4)
    assert len(header) == 256

    # Append channel headers
    for name in included:
        header += name.ljust(16)
    header += Nc * ''.ljust(80)
    header += Nc * ''.ljust(8)
    header += Nc * str(physmin).ljust(8)
    header += Nc * str(physmax).ljust(8)
    header += Nc * str(digmin).ljust(8)
    header += Nc * str(digmax).ljust(8)
    header += Nc * ''.ljust(80)
    header += Nc * str(nsblock).ljust(8)
    header += Nc * ''.ljust(32)
    assert len(header) == 256 * (1 + Nc) == nbytes

    # Now write this lot to the file
    fout.write(header)

    # Now write blocks of data
    scale = (digmax - digmin) / (physmax - physmin)
    block = np.zeros((nsblock, Nc))
    iblock = np.zeros((nsblock * Nc), '<i2')

    # Function to write out one block from the line buffer
    def lines_to_block(lines_block):
        assert len(lines_block) == nsblock
        for n, l in enumerate(lines_block):
            numstrs = l.split(',')[1:]
            block[n, :] = [float(numstrs[i]) for i in indices]
        fblock = block.T.flatten() # Reshape
        iblock[:] = (digmin + scale * (fblock - physmin))
        return iblock

    # Iterate over lines
    for line in itf:
        lines_buffer.append(line)
        # Write complete blocks, when possible
        while len(lines_buffer) > nsblock:
            lines_block, lines_buffer = lines_buffer[:nsblock], lines_buffer[nsblock:]
            iblock = lines_to_block(lines_block)
            iblock.tofile(fout)
            fout.flush()
            del lines_block
            gc.collect()
    # Append last partial block with zeros
    if lines_buffer:
        # Add dummy lines_buffer to buffer and write out rest
        blank_line = ','.join(Ncols * ['0'])
        lines_buffer.extend((nsblock - len(lines_buffer)) * [blank_line])
        iblock = lines_to_block(lines_buffer)
        iblock.tofile(fout)
        fout.flush()

    # Finished
    fout.close()


# Script entry point
opts = [
    ('a', 'all-chan', False, 'Include all channels'),
]
@opster.command(usage='%name [OPTS] [FILE1 ...]', options=opts)
def main(*edfnames, **opts):
    """Convert csv files to edf format.

    1) %name [OPTS] FILE1.csv ... : convert FILE1.csv to FILE1.edf and so on
    2) %name [OPTS]           : read cdv from stdin and write edf to stdout
    """
    if edfnames:
        # Iterate over arguments
        for e in edfnames:
            csv_to_edf(e, opts['all_chan'])
    else:
        # stdin to stdout
        csv_to_edf(None, opts['all_chan'])


if __name__ == "__main__":
    main(argv=sys.argv[1:])


