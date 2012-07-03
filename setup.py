#!/usr/bin/env python
#
# Copyright Oscar Benjamin 2011 under GPLv3+

import os, os.path, sys
from distutils.core import setup
from distutils.extension import Extension
from distutils.command.build_ext import build_ext

import numpy


# Use Cython to generate the c-files if available
#
# The source distribution includes the generated c-files so someone can build
# it without needing cython (although they still need a c-compiler).
#
# Ensuring that the c-files are distributed in the source distribution means
# adding the names of the c-files to MANIFEST.in and ensuring that the cython
# generated c-files are up to date by remembering to do:
#
#    $ python setup.py build_ext
#    $ python setup.py sdist
#
try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True


# Need to link with libm and librt on posix but not on Windows.
if not 'win' in sys.platform:
    libraries_cysode = ['m', 'rt']
    libraries_examples = ['m']
else:
    libraries_cysode = libraries_examples = []


cmdclass = {}
ext_modules = []

# Need to add the pyxname of compiling fully with cython or the cname when
# compiling from sdist without cython.
def add_cython_ext_module(modname, pyxname, cnames=[], **kwargs):
    # Use cython to copmile the pyx file
    if use_cython:
        ext = Extension(modname, [pyxname] + cnames, **kwargs)
    # Use distutils to compile the c file
    else:
        cname = os.path.splitaxt(pyxname)[0] + '.c'
        ext = Extension(modname, [cname] + cnames, **kwargs)
    ext_modules.append(ext)


# This extension module is a core part of sode
add_cython_ext_module('sode.cysode',
    os.path.join('sode', 'cysode.pyx'),
    [os.path.join('sode', 'cfiles', 'randnorm.c')],
    include_dirs = [numpy.get_include(), '.'],
    libraries = libraries_cysode,
)


# This extension module is for the examples
add_cython_ext_module('sode.examples.cyfiles.examples',
    os.path.join('sode', 'examples', 'cyfiles', 'examples.pyx'),
    include_dirs = [numpy.get_include(), '.'],
    libraries = libraries_cysode,
)


cprog_cfiles = ['main.c', 'examples.c', 'randnorm.c', 'solvers.c']
cprog_cfiles = [os.path.join('sode', 'cfiles', p) for p in cprog_cfiles]
cprog_name = os.path.join('sode', 'examples', 'cexamples')


# We also want to build the standalone c program for the examples
class MonkeyPatch_build_ext(build_ext):
    def run(self):
        build_ext.run(self)
        self.compiler.link_executable(cprog_cfiles, cprog_name)


cmdclass['build_ext'] = MonkeyPatch_build_ext


# Use the README.rst file as the front-page on PyPI
#with open('README.rst') as README:
#    LONG_DESCRIPTION = README.read()


# Standard distutils setup
setup(
    # First the metadata
    name = 'sode',
    version = '0.0.1',
    author = 'Oscar Benjamin',
    author_email = 'oscar.j.benjamin@gmail.com',
    url = 'https://github.com/oscarbenjamin/sode',
    description = ('Python/Cython lib for solving ' +
                   'Stochastic Ordinary Differential Equations'),
    # long_description = open('README.rst').read()
    platforms = ['linux2', 'win32', 'darwin'], # tested on these
    license = 'GPLv3+',

    classifiers = [
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

        'Intended Audience :: Science/Research',

        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        'Operating System :: POSIX :: Linux',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: MacOS :: MacOS X',

        'Programming Language :: Cython',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],

    # Package management concepts
    provides = ['sode'],
    requires = ['numpy (>=1.4)', 'cython (>=0.15)'],

    # Now the content
    cmdclass = cmdclass,
    ext_modules = ext_modules,
    packages = ['sode', 'sode.examples'],
    scripts = ['scripts/sode', 'scripts/sode.bat']
)
