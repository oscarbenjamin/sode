#!/usr/bin/env python
# Copyright Oscar Benjamin 2011
#
# This script defines a trivial SODE representing just a Weiner process. This
# script is called when sode is run at the cmd line.


import sys

from sode.cysode import SODE_test, CYSODE_test
from sode.script import Script, MultiScript

main = MultiScript({None: 'sode', 'sode':SODE_test, 'cysode':CYSODE_test}).main
main(sys.argv[1:])
