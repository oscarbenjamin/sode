#!/usr/bin/env python
#
# Author: Oscar Benjamin
# Date: 11 Mar 2011
#
# Example of a script to investigate a 1-D SODE. This script is a
# reimplementation of weiner.py using a more complex form of the diffusion and
# drift coefficient functions. The more complex form is required for
# higher-dimensional SODEs.


from sode.pysode  import SODE


class Weiner2(SODE):
    """Weiner process with drift coeff. mu and diffusion coeff. sigma

        dx(t) = mu dt + sigma dW(t)

    This trivial SODE has the exact solution:

        x(t) = x(0) + mu t + sigma W(t)
    """
    # SODE subclasses must define the variables and parameters of the system
    # and also the drift and diffusion coefficients as functions of the state
    # x and time t. When the instance is created it is given attributes
    # corresponding to the variable names. These attributes do not represent
    # the value of the variable at any time. They are indices to be used when
    # referencing the arrays a, b and x passed in to the drift and diffusion
    # functions.

    # For each variable, the default initial condition must be provided.
    variables = (('x', 0),)

    # For each parameter, the default value must be provided.
    parameters = (('mu', 0), ('sigma', 1))

    # Drift coefficient. The drift coefficient should initialise the values of
    # the vector a and then return it. 'self.x' is the index of the variable
    # with name 'x'.
    def drift(self, a, x, t):
        a[self.x] = self.mu
        return a

    # Diffusion coefficient. This should also be a vector
    def diffusion(self, b, x, t):
        b[self.x] = self.sigma
        return b

    # The exact solution. If provided, this will be used to check numerical
    # convergence when running the conv command.
    def exact(self, x0, t, Wt):
        return (x0 + self.mu * t) + self.sigma * Wt


if __name__ == "__main__":
    # The Script class imported from sode.pysode allows us to easily turn this file
    # into a script that can be used to investigate the SODE numerically.
    import sys
    from sode.script import Script
    script = Script(Weiner2)
    script.main(argv=sys.argv[1:])
