#!/usr/bin/env python
#
# Author: Oscar Benjamin
# Date: 11 Mar 2011
#
# Example of a script to investigate a 1-D SODE. In this particular case, the
# exact solution is known, as this is a trivial SODE. It serves as an example
# of how to write a 1-D SODE and also to check the numerical performance of
# the sode solver routines.
#
# 1-D SODEs can be written in a particularly simple way. weiner2.py shows how
# to rewrite this same SODE using the syntax needed for higher-dimensional
# SODEs.


from sode.pysode  import SODE


class Weiner(SODE):
    """Weiner process with drift coeff. mu and diffusion coeff. sigma

        dx(t) = mu dt + sigma dW(t)

    This trivial SODE has the exact solution:

        x(t) = x(0) + mu t + sigma W(t)
    """
    # SODE subclasses must define the variables and parameters of the system
    # and also the drift and diffusion coefficients as functions of the state
    # x and time t. When a Weiner instance is created, the parameters
    # values are stored as atributes of self. Note that keyword arguments used
    # when creating an SODE instance can alter parameters from their default
    # values provided in the class definition.

    # For each variable, the default initial condition must be provided.
    variables = (('x', 0),)

    # For each parameter, the default value must be provided.
    parameters = (('mu', 0), ('sigma', 1))

    # The drift coefficient of the SODE. We can access the parameter mu by
    # simply writing 'self.mu'.
    def drift(self, a, x, t):
        return self.mu

    # The diffusion coefficient of the SODE.
    def diffusion(self, b, x, t):
        return self.sigma

    # The exact solution. If provided, this will be used to check numerical
    # convergence when running the conv command.
    def exact(self, x0, t, Wt):
        return (x0 + self.mu * t) + self.sigma * Wt


if __name__ == "__main__":
    # The Script class imported from sode.pysode allows us to easily turn this file
    # into a script that can be used to investigate the SODE numerically.
    import sys
    from sode.script import Script
    script = Script(Weiner)
    script.main(argv=sys.argv[1:])
