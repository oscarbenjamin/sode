#!/usr/bin/env python
#
# Author: Oscar Benjamin
# Date: 11 Mar 2011
#
# This module uses the framework defined in sode.py to implement a number of
# simple Stochastic Ordinary Differential Equations (SODEs) with known exact
# solutions, for illustration and unit testing.
#
# The examples are mainly taken from "Numerical solution of stochastic
# differential equations" by Peter E. Kloeden and Eckhard Platen


from sode import SODE, Script


__all__ = ['WeinerSODE']


# 1-D SODE subclasses can be written in the simple way.
# arguments a and b can be ignored in this case.
class WeinerSODE(SODE):
    """Weiner process with drift coeff. mu and diffusion coeff. sigma

    dx(t) = mu dt + sigma dW(t)
    """
    variables = (('x', 0),)
    parameters = (('mu', 0), ('sigma', 1))
    def drift(self, a, x, t):
        return self.mu
    def diffusion(self, b, x, t):
        return self.sigma
    def exact(self, x0, t, Wt):
        return (x0 + self.mu * t) + self.sigma * Wt

if __name__ == "__main__":

    import sys

    # Actually run as a script
    script = Script(WeinerSODE)
    script.main(argv=sys.argv[1:])
