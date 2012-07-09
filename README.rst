SODE
====

Introduction
------------
SODE is a python/cython library for generating numerical solutions to Stochastic Ordinary
Differential Equations (SODEs)

Setting up the sode-module
--------------------------
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~~~~~~~~~~~~~~~~~~~~

~~~~~~~~~~~~~~~~~~~

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




of Python is running by typing in the terminal (before you run the set-up)::



1. Opster
~~~~~~~~~

    $ cd /path/to/downloads
--------------------------------
~~~~~~~~~~~~~~~~~~~~~~~~~~
python terminal by printing the value of sys.path:
directory of the download. This file is used to install the python package so that it 
is on your python path. If you run::
can be imported from within python. If the package also has c code or other code that needs 
to be compiled, then the setup.py will also compile those for you. For example, sode has 
some cython modules that need compiling. To do this you would run::
along with copies of all the python files. Then when you run::
can be imported from within python. 
running python by setting the PYTHONPATH environment variable::
to PYTHONPATH and adds the "./scripts" directory to PATH so that the shell finds the scripts 
in there.
~~~~~~~~
programs that are run from within that script. This means that to run this script you need to 
source it::
the file "in-place" rather than installing them you will first need to build the cython files 
in the current directory (rather than build/lib) using the command::