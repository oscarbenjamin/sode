#
# Copyright Oscar Benjamin 2011
# sode/__init__py
#

def get_include():
    import os.path
    import numpy
    sode_base_dir, _ = os.path.split(__file__)
    cdir = os.path.join(sode_base_dir, 'cfiles')
    sode_base_dir, _ = os.path.split(sode_base_dir)
    cydir = sode_base_dir
    return [cdir, cydir, numpy.get_include()]
