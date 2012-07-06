#!/usr/bin/env python
#
# Copyright Oscar Benjamin 2011 under GPLv3+

import os, os.path, sys
from distutils.core import setup
from distutils.extension import Extension
from distutils.command.build_ext import build_ext
from distutils.cygwinccompiler import Mingw32CCompiler # ming32 fix

import numpy


# Use Cython to generate the c-files if available
#
# The source distribution includes the generated c-files so someone can build
# it without needing cython (although they still need a c-compiler).
#
# Ensuring that the c-files are distributed in the source distribution means
# adding the c-files to MANIFEST.in and ensuring that the cython generated
# c-files are up to date by remembering to do:
#
#    $ python setup.py build_ext --inplace
#    $ python setup.py sdist
#
# See here:
# http://stackoverflow.com/questions/4505747/how-should-i-structure-a-python-package-that-contains-cython-code
try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True


# These data structures are conditional on
# 1) Whether cython is installed (ext_modules, cmdclass)
# 2) Whether running on Windows (scripts/*.bat)
cmdclass = {}
ext_modules = []
scripts = []


# :::::::::::::::::::::::::::::
#   Extension modules
# :::::::::::::::::::::::::::::

# Need to add the .pyx name if compiling fully with cython or the .c name when
# compiling from sdist without cython.
def add_cython_ext_module(modname, cname, cnames=[], **kwargs):
    # Use cython to cythonise the .pyx files if possible
    if use_cython:
        cname = os.path.splitext(cname)[0] + '.pyx'
    ext_modules.append(Extension(modname, [cname] + cnames, **kwargs))


# Need to link with libm and librt on posix but not on Windows.
if not 'win' in sys.platform:
    libs_cysode = ['m', 'rt']
    libs_cyexamples = ['m']
    libs_cexamples = ['m', 'rt']
else:
    libs_cysode = []
    libs_cyexamples = []
    libs_cexamples = []


# This extension module is a core part of sode
add_cython_ext_module('sode.cysode',
    os.path.join('sode', 'cysode.c'),
    [os.path.join('sode', 'cfiles', 'randnorm.c')],
    include_dirs = [numpy.get_include(), '.'],
    libraries = libs_cysode,
)


# This extension module is for the examples
add_cython_ext_module('sode.examples.cyexamples.examples',
    os.path.join('sode', 'examples', 'cyexamples', 'examples.c'),
    include_dirs = [numpy.get_include(), '.'],
    libraries = libs_cyexamples,
)


# ::::::::::::::::::::::::::::::::
#   Scripts
# ::::::::::::::::::::::::::::::::

# Executable entry points in scripts
scripts = ['sode', 'sode-pyexamples', 'sode-cyexamples']
scripts_dir = 'scripts'


# These .bat files are needed so that Windows users can run the scripts
# without the .py extension. Giving the scripts .py extensions can cause
# problems with stdin/stdout redirection on Windows
# http://support.microsoft.com/default.aspx?kbid=321788
if 'win' in sys.platform:
    scripts.extend([sname + '.bat' for sname in scripts])

scripts.append('sode-time')

scripts = [os.path.join(scripts_dir, sname) for sname in scripts]


# ::::::::::::::::::::::::::::::::
#   Standalone c program
# ::::::::::::::::::::::::::::::::

# monkey patch build_ext to also build a standalone c program
# http://mail.python.org/pipermail/distutils-sig/2009-September/013216.html
class MonkeyPatch_build_ext(build_ext):
    def run(self):
        global scripts
        build_ext.run(self)
        exe_name = build_examples_cprog(self.compiler)
        scripts.append(exe_name)

cmdclass['build_ext'] = MonkeyPatch_build_ext


# We also need to build the examples c-program
def build_examples_cprog(compiler):
    if isinstance(compiler, Mingw32CCompiler):
        mingw32_compiler_fix(compiler)

    cfiles = ['main.c', 'examples.c', 'randnorm.c', 'solvers.c']
    cfiles = [os.path.join('sode', 'cfiles', p) for p in cfiles]
    exe_name = os.path.join(scripts_dir, 'sode-cexamples')

    compiler.link_executable(cfiles, exe_name, libraries=libs_cexamples)

    # compiler.exe_extension is None on posix, string on Windows
    if compiler.exe_extension is not None:
        exe_name += compiler.exe_extension

    return exe_name


# msvcr90 needs manifest etc. See:
# http://developer.berlios.de/devlog/akruis/2012/06/10/msvcr90dll-and-mingw/
# Since this is a standalone program we can just link against the old version
# of msvcr.
def mingw32_compiler_fix(compiler):
    for n, dllname in enumerate(compiler.dll_libraries):
        if dllname.startswith('msvcr'):
            compiler.dll_libraries[n] = min(dllname, 'msvcr71')
            return


# :::::::::::::::::::::::::::::::
#   Bring it all together
# :::::::::::::::::::::::::::::::

# Use the README.rst file as the front-page on PyPI
#with open('README.rst') as README:
#    LONG_DESCRIPTION = README.read()


# Standard distutils setup
setup(
    # First the metadata
    name = 'sode',
    version = '0.0.3',
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

    # Dependency information
    provides = ['sode'],
    requires = ['numpy (>=1.4)', 'cython (>=0.15)'],

    # Now the content
    cmdclass = cmdclass,
    ext_modules = ext_modules,
    packages = ['sode', 'sode.examples', 'sode.examples.pyexamples',
                                         'sode.examples.cyexamples'],
    scripts = scripts
)
