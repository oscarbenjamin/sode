#!/usr/bin/env python
#
# Copyright Oscar Benjamin 2011 under new BSD license

import os, os.path, sys
from distutils.core import setup
from distutils.extension import Extension

from Cython.Distutils import build_ext
import numpy

import sode


if not 'win' in sys.platform:
    libraries_cysode = ['m', 'rt']
    libraries_examples = ['m']
else:
    libraries_cysode = libraries_examples = []


ext_modules = [
    Extension(
        'sode.cysode',
        [os.path.join('sode', 'cysode.pyx'),
         os.path.join('sode', 'cfiles', 'randnorm.c')],
        include_dirs=[numpy.get_include(), '.'],
        libraries=libraries_cysode
    ),
    Extension(
        'sode.examples.cyfiles.examples',
        [os.path.join('sode', 'examples', 'cyfiles', 'examples.pyx')],
        include_dirs=[numpy.get_include()],
        libraries = libraries_examples
    )
]

setup(
    name = 'sode',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
)
