#!/usr/bin/env python
#
# Author: Oscar Benjamin
# Date: 11 Mar 2011
#
# Example of a script to investigate an SODE with an exactly known solution.


import numpy as np

from sode import SODE, Script


class SinusoidalMultiplicative(SODE):
    """Nonlinear 1-D SODE with solutions sinusoidal in Wt

        dx(t) = -(1/2) a^2 x(t) dt + a sqrt(1 - x(t)^2) dW(t)

    This equation has the exact solution

        x(t) = sin( a W(t) + arcsin(x(0)) )
    """
    variables = (('x', 0),)
    parameters = (('a', 0.01),)
    def drift(self, a, x, t):
        return - (self.a**2 * x) / 2
    def diffusion(self, b, x, t):
        return self.a * np.sqrt(1 - x**2)
    def exact(self, x0, t, Wt):
        return np.sin(self.a * Wt + np.arcsin(x0))



if __name__ == "__main__":
    import sys
    script = Script(SinusoidalMultiplicative)
    script.main(argv=sys.argv[1:])
