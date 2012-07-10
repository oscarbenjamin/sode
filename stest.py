#!/usr/bin/env python
# Copyright Oscar Benjamin 2011
#
# This script defines a trivial SODE representing just a Weiner process. This
# script is called when sode is run at the cmd line.


import sys

from sode.cysode import SODE_test, CYSODE_test
from sode.script import Script, MultiScript

class SODE_test_py(SODE_test):
    def _cy_drift(self, a, x, t):
        a[0] = x[0]
        a[1] = x[1]

systems = {
    None: 'sode',
    'sode':SODE_test,
    'cysode':CYSODE_test,
    'py':SODE_test_py,
}
main = MultiScript(systems).main
main(sys.argv[1:])
