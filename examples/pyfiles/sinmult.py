#!/usr/bin/env python
#
# Author: Oscar Benjamin
# Date: 11 Mar 2011
#
# Example of a script to investigate an SODE with an exactly known solution.


import numpy as np

from sode.pysode  import SODE


class SinusoidalMultiplicative(SODE):
    """Nonlinear 1-D SODE with solutions sinusoidal in Wt

        dx(t) = -(1/2) alpha^2 x(t) dt + alpha sqrt(1 - x(t)^2) dW(t)

    This equation has the exact solution

        x(t) = sin( alpha W(t) + arcsin(x(0)) )
    """
    variables = (('x', 0),)
    parameters = (('alpha', 0.01),)
    def drift(self, alpha, x, t):
        return - (self.alpha**2 * x) / 2
    def diffusion(self, b, x, t):
        return self.alpha * np.sqrt(1 - x**2)
    def exact(self, x0, t, Wt):
        return np.sin(self.alpha * Wt + np.arcsin(x0))



if __name__ == "__main__":
    import sys
    from sode.script import Script
    script = Script(SinusoidalMultiplicative)
    script.main(argv=sys.argv[1:])
