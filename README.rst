SODE
====

Introduction
------------

SODE is a python/cython library for generating numerical solutions to
Stochastic Ordinary Differential Equations (SODEs)

Setting up the sode-module
--------------------------

1. Compile C and Cython code:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To compile the C and Cython code, you need to type in the terminal::

    $ cd /path/to/sode
    $ python setup.py build_ext --inplace --force

2. Run the sode module
~~~~~~~~~~~~~~~~~~~~~~

Run the following commands (every time), type in the terminal::

    $ cd /path/to/sode
    $ . ./addpath.sh   # note: this line begins with a "." followed by a space " "

3. Run the examples
~~~~~~~~~~~~~~~~~~~

Type in the terminal::

    $ sode solve -p

4. Summary (regular use of the module)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Type in the terminal::

    $ cd sode
    $ . ./addpath.sh
    $ sode-time
    $ sode solve -p

5. Help and specific simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Type in the terminal, for a specific simulation of the tanmult-system::

    $ sode -h
    $ sode-pyexamples solve -s tanmult -p

Troubleshooting
---------------

Check which version of Python is running:

Method 1::

    $ which python

Method 2::

    $ python
    >>> import sys
    >>> sys.version_info
    >>> sys.executable

In case these methods provide refer to different versions, ensure that a
correct version of Python is running by typing in the terminal (before you run
the set-up)::

    $ export PATH="/Library/Frameworks/Python.framework/Versions/7.2/Resources/Python.app/Contents/MacOS:${PATH}"

Check which version of Cython is running::

    $ cython --version
    $ cython.py –version

MAC OSX

1. Opster
~~~~~~~~~

Download Opster for MAC (http://pypi.python.org/pypi/opster), and type in
terminal::

    $ cd /path/to/downloads
    $ tar -xzf opster-3.7.tar.gz
    $ cd opster-3.7
    $ python setup.py install

2. Install C-compiler:
~~~~~~~~~~~~~~~~~~~~~~

http://woss.name/2012/01/24/how-to-install-a-working-set-of-compilers-on-mac-os-x-10-7-lion/

https://github.com/kennethreitz/osx-gcc-installer/downloads

Appendix: background information
--------------------------------
Installing the sode module
~~~~~~~~~~~~~~~~~~~~~~~~~~

Python imports modules that are on its import path. You can see the import
path from a python terminal by printing the value of sys.path::

    >>> import sys
    >>> sys.path

When you download a python package there is usually a file called setup.py in
the root directory of the download. This file is used to install the python
package so that it is on your python path. If you run::

    $ python setup.py install

Then the package will be copied into a directory on your python path so that
its modules can be imported from within python. If the package also has c code
or other code that needs to be compiled, then the setup.py will also compile
those for you. For example, sode has some cython modules that need compiling.
To do this you would run::

    $ python setup.py build

This compiles everything and places the compiled files in a subdirectory
called build/lib along with copies of all the python files. Then when you
run::

    $ python setup.py install

it will copy the compiled files and the python files onto your python path so
that they can be imported from within python.

Using sode without installation If you don't install the files you can
temporarily add them to your python path before running python by setting the
PYTHONPATH environment variable::

     $ export PYTHONPATH=/home/username/work/sode
     $ python
     >>> import sode

If you just want to set PYTHONPATH for one command in the terminal, you can
do::

     $ PYTHONPATH=/home/username/work/sode python
     >>> import sode

This is done automatically by using "addpath.sh", which adds the current
directory, ".", to PYTHONPATH and adds the "./scripts" directory to PATH so
that the shell finds the scripts in there.

Sourcing
~~~~~~~~

Normally when you run a script any changes it makes to environment variables
will only affect programs that are run from within that script. This means
that to run this script you need to source it::

     $ source addpath.sh

or (note the "." at the start)::

     $ . ./addpath.sh

After this any python scripts that use sode should be able to import it. If
you want to use the file "in-place" rather than installing them you will first
need to build the cython files in the current directory (rather than
build/lib) using the command::

     $ python setup.py build_ext --inplace

This places the compiled files in the same place as the other files.
