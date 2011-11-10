#!/usr/bin/env python
#
# Copyright Oscar Benjamin 2011 under new BSD license

import os.path
from distutils.core import setup
from distutils.extension import Extension

from Cython.Distutils import build_ext
import numpy

import sode

def sode_extension(modname, pyxname):
    return Extension(
        modname,
        [pyxname, os.path.join('sode', 'cfiles', 'randnorm.c')],
        include_dirs=[numpy.get_include(), '.']
    )

ext_modules = [
    sode_extension('sode.cysode',
                   os.path.join('sode', 'cysode.pyx')),
    sode_extension('sode.examples.cyfiles.examples',
                   os.path.join('sode', 'examples', 'cyfiles', 'examples.pyx')),
]

setup(
    name = 'sode',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
)
