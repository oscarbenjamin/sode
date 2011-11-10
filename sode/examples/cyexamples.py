#!/usr/bin/env python
#
# Copyright, Oscar Benjamin 2011
#
# This script brings together the example pure python SODE systems from the 
# pyfiles package for comparison with the faster cython and pure-c
# implementations.

from sode.script import MultiScript

from cyfiles.weiner import Weiner
#from pyfiles.weiner2 import Weiner2
#from pyfiles.lin2d import LinearAdditive2D
#from pyfiles.linmult import LinearMultiplicative
#from pyfiles.sinmult import SinusoidalMultiplicative
#from pyfiles.tanmult import TangentMultiplicative


sysdict = {
    None:'weiner',
    'weiner':Weiner,
#    'weiner2':Weiner2,
#    'lin2d':LinearAdditive2D,
#    'linmult':LinearMultiplicative,
#    'sinmult':SinusoidalMultiplicative,
#    'tanmult':TangentMultiplicative,
}
examples_script = MultiScript(sysdict)


def main(argv):
    examples_script.main(argv)


if __name__ == "__main__":
    import sys
    main(argv=sys.argv[1:])
