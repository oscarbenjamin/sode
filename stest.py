#!/usr/bin/env python
# Copyright Oscar Benjamin 2011
#
# This script defines a trivial SODE representing just a Weiner process. This
# script is called when sode is run at the cmd line.


import sys

from sode.cysode import SODE_test, SODE
from sode.script import Script, MultiScript

class SODE_test_py(SODE):

    def __init__(self):
        self._set_variables(['x', 'y'], [1, 1])

    def _drift(self, a, x, t):
        a[0] = x[1]
        a[1] = - x[0]

    def _diffusion(self, b, x, t):
        b[0] = b[1] = 0.01

systems = {
    None: 'test',
    'test':SODE_test,
    'py':SODE_test_py,
}
main = MultiScript(systems).main
main(sys.argv[1:])
