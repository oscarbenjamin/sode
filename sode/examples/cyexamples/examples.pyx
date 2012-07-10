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


from libc cimport math
import numpy as np
cimport numpy as np
np.import_array()

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
    cdef void _drift(self, double* a, double* x, double t):
        a[self.x] = self.mu

    # The diffusion coefficient of the SODE.
    cdef void _diffusion(self, double* b, double* x, double t):
        b[self.x] = self.sigma

    # The exact solution. If provided, this will be used to check numerical
    # convergence when running the conv command.
    def exact(self, x0, t, Wt):
        return (x0 + self.mu * t) + self.sigma * Wt

# Declare twice for compatibility with pyfiles
Weiner2 = Weiner


cdef class LinearAdditive2D(CYSODE):
    """Linear 2-D homogenous SODE with additive noise

    dx(t) = ( - beta x + alpha y(t) ) dt + dW1(t)
    dy(t) = ( - beta y - alpha x(t) ) dt + dW2(t)
    """
    variables = (('x', 0), ('y', 0))
    parameters = (('alpha', 1), ('beta', .01))
    cdef public variable x, y
    cdef public parameter alpha, beta

    cdef void _drift(self, double* a, double* x, double t):
        a[self.x] = - self.beta * x[self.x] + self.alpha * x[self.y]
        a[self.y] = - self.beta * x[self.y] - self.alpha * x[self.x]

    cdef void _diffusion(self, double* b, double* x, double t):
        b[self.x] = b[self.y] = 1


cdef class LinearMultiplicative(CYSODE):
    """Linear 1-D homogenous SODE with constant coefficients and
    multiplicative noise.

        dx(t) = alpha x(t) dt + beta x(t) dW(t)

    This SODE has the exact solution:

        x(t) = x(0) exp( (alpha - (1/2) beta^2) t + beta W(t) )
    """
    variables = (('x', 1),)
    parameters = (('alpha', -0.0001), ('beta', 1))
    cdef public variable x
    cdef public parameter alpha, beta

    cdef void _drift(self, double* a, double* x, double t):
        a[self.x] = self.alpha * x[self.x]

    cdef void _diffusion(self, double* b, double* x, double t):
        b[self.x] = self.beta * x[self.x]

    def exact(self, x0, t, Wt):
        factor = self.alpha - 0.5 * self.beta ** 2
        return x0 * np.exp(factor * t + self.beta * Wt)


cdef class SinusoidalMultiplicative(CYSODE):
    """Nonlinear 1-D SODE with solutions sinusoidal in Wt

        dx(t) = -(1/2) alpha^2 x(t) dt + alpha sqrt(1 - x(t)^2) dW(t)

    This equation has the exact solution

        x(t) = sin( alpha W(t) + arcsin(x(0)) )
    """
    variables = (('x', 0),)
    parameters = (('alpha', 0.01),)
    cdef public variable x
    cdef public parameter alpha

    cdef void _drift(self, double* a, double* x, double t):
        a[self.x] = - (self.alpha**2 * x[self.x]) / 2

    cdef void _diffusion(self, double* b, double* x, double t):
        b[self.x] = self.alpha * math.sqrt(1 - x[self.x]**2)

    def exact(self, x0, t, Wt):
        return np.sin(self.alpha * Wt + np.arcsin(x0))


cdef class TangentMultiplicative(CYSODE):
    """Nonlinear 1-D SODE with solutions in tan(Wt)

        dx(t) = alpha^2 x(t) (1 + x(t)^2) dt + alpha (1 + x(t)^2) dW(t)

    This SODE has the exact solution:

        x(t) = tan( alpha W(t) + arctan(x(0)) )
    """
    variables = (('x', 0),)
    parameters = (('alpha', 0.01),)
    cdef public variable x
    cdef public parameter alpha

    cdef void _drift(self, double* a, double* x, double t):
        a[self.x] = self.alpha**2 * x[self.x] * (1 + x[self.x]**2)

    cdef void _diffusion(self, double* b, double* x, double t):
        b[self.x] = self.alpha * (1 + x[self.x]**2)

    def exact(self, x0, t, Wt):
        return np.tan(self.alpha*Wt + np.arctan(x0))
