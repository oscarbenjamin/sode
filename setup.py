#!/usr/bin/env python
#
# Copyright Oscar Benjamin 2011 under new BSD license

import os.path
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy


def example_extension(name):
    return Extension(
        'sode.examples.cyfiles.{0}'.format(name),
        [os.path.join('sode', 'examples', 'cyfiles', '{0}.pyx'.format(name)),
         os.path.join('sode', 'cfiles', 'randnorm.c')],
        include_dirs=[numpy.get_include(), '.']
    )

ext_modules = [
    Extension(
        'sode.cysode',
        [os.path.join('sode', 'cysode.pyx'),
         os.path.join('sode', 'cfiles', 'randnorm.c')],
        include_dirs=[numpy.get_include()]
    ),
    example_extension('weiner'),
]

setup(
    name = 'sode',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
)
