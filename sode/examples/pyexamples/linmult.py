#!/usr/bin/env python
#
# Author: Oscar Benjamin
# Date: 11 Mar 2011
#
# Example of a script to investigate an SODE with an exactly known solution.


import numpy as np

from sode.pysode  import SODE


class LinearMultiplicative(SODE):
    """Linear 1-D homogenous SODE with constant coefficients and
    multiplicative noise.

        dx(t) = alpha x(t) dt + beta x(t) dW(t)

    This SODE has the exact solution:

        x(t) = x(0) exp( (alpha - (1/2) beta^2) t + beta W(t) )
    """
    variables = (('x', 1),)
    parameters = (('alpha', -0.0001), ('beta', 1))
    def drift(self, a, x, t):
        return self.alpha * x
    def diffusion(self, b, x, t):
        return self.beta * x
    def exact(self, x0, t, Wt):
        factor = self.alpha
        return x0 * np.exp(factor * t + self.beta * Wt)


if __name__ == "__main__":
    import sys
    from sode.script import Script
    script = Script(LinearMultiplicative)
    script.main(argv=sys.argv[1:])
