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
# This script illustrates how to generate a CYSODE subclass. These are
# implemented in cython and are more efficient than those implemented in pure
# python.


from sode.cysode cimport CYSODE, parameter, variable


# Need to declare the class with cdef
cdef class Weiner(CYSODE):
    """Weiner process with drift coeff. mu and diffusion coeff. sigma

        dx(t) = mu dt + sigma dW(t)

    This trivial SODE has the exact solution:

        x(t) = x(0) + mu t + sigma W(t)
    """
    # CYSODE subclasses must define the variables and parameters of the system
    # and also the drift and diffusion coefficients as functions of the state
    # x and time t. When a Weiner instance is created, the parameters
    # values are stored as atributes of self. Note that keyword arguments used
    # when creating an SODE instance can alter parameters from their default
    # values provided in the class definition.

    # For each variable, the default initial condition must be provided.
    variables = (('x', 0),)

    # For each parameter, the default value must be provided.
    parameters = (('mu', 0), ('sigma', 1))

    # For a CYSODE sublclass we also need to declare the parameters and
    # variables like so:
    cdef public variable x
    cdef public parameter mu, sigma

    def __cinit__(self):
        self.x = 0
        self.mu = 0.0
        self.sigma = 1.0

    # The drift coefficient of the SODE. We can access the parameter mu by
    # simply writing 'self.mu'.
    cdef _drift(self, double* a, double* x, double t):
        a[self.x] = self.mu

    # The diffusion coefficient of the SODE.
    cdef _diffusion(self, double* b, double* x, double t):
        b[self.x] = self.sigma

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
