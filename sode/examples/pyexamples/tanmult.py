#!/usr/bin/env python
#
# Author: Oscar Benjamin
# Date: 11 Mar 2011
#
# Example of a script to investigate an SODE with an exactly known solution.


import numpy as np

from sode.pysode  import SODE


class TangentMultiplicative(SODE):
    """Nonlinear 1-D SODE with solutions in tan(Wt)

        dx(t) = alpha^2 x(t) (1 + x(t)^2) dt + alpha (1 + x(t)^2) dW(t)

    This SODE has the exact solution:

        x(t) = tan( alpha W(t) + arctan(x(0)) )
    """
    variables = (('x', 0),)
    parameters = (('alpha', 0.01),)
    def drift(self, alpha, x, t):
        return self.alpha**2 * x * (1 + x**2)
    def diffusion(self, b, x, t):
        return self.alpha * (1 + x**2)
    def exact(self, x0, t, Wt):
        return np.tan(self.alpha*Wt + np.arctan(x0))


if __name__ == "__main__":
    import sys
    from sode.script import Script
    script = Script(TangentMultiplicative)
    script.main(argv=sys.argv[1:])
