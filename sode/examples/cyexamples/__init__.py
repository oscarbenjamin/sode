#!/usr/bin/env python
#
# Copyright, Oscar Benjamin 2011
#
# This script brings together the example pure python SODE systems from the 
# pyfiles package for comparison with the faster cython and pure-c
# implementations.

from sode.script import MultiScript

from sode.examples.cyexamples.examples import (
    Weiner, Weiner2, LinearAdditive2D, LinearMultiplicative,
    SinusoidalMultiplicative, TangentMultiplicative)

sysdict = {
    None:'weiner',
    'weiner':Weiner,
    'weiner2':Weiner2,
    'lin2d':LinearAdditive2D,
    'linmult':LinearMultiplicative,
    'sinmult':SinusoidalMultiplicative,
    'tanmult':TangentMultiplicative,
}
examples_script = MultiScript(sysdict)


def main(argv):
    examples_script.main(argv)


if __name__ == "__main__":
    import sys
    main(argv=sys.argv[1:])
