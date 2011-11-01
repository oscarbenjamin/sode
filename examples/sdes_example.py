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


import numpy as np

from sode import SODE


# 1-D SODE subclasses can be written in the simple way (ignore arguments a and
# b)
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

# The more complicated way is required for higher dimensional SODEs
# The same SODE rewritten:
class WeinerSODE2(SODE):
    """Weiner process with drift coeff. mu and diffusion coeff. sigma

    dx(t) = mu dt + sigma dW(t)
    """
    variables = (('x', 0),)
    parameters = (('mu', 0), ('sigma', 1))
    def drift(self, a, x, t):
        a[self.x] = self.mu
        return a
    def diffusion(self, b, x, t):
        b[self.x] = self.sigma
        return b

# 2-D example. self.x and self.y are indices into a, b and x
class LinearAdditive2D(SODE):
    """Linear 2-D homogenous SODE with additive noise

    dx(t) = ( - beta x + alpha y(t) ) dt + dW1(t)
    dy(t) = ( - beta y - alpha x(t) ) dt + dW2(t)
    """
    variables = (('x', 0), ('y', 0))
    parameters = (('alpha', 1), ('beta', .01))
    def drift(self, a, x, t):
        a[self.x] = - self.beta * x[self.x] + self.alpha * x[self.y]
        a[self.y] = - self.beta * x[self.y] - self.alpha * x[self.x]
        return a
    def diffusion(self, b, x, t):
        b[self.x] = b[self.y] = 1
        return b

# What follows are some standard 1-D SODEs with explicit solutions

class LinearAdditive(SODE):
    """Linear 1-D homogenous SODE with additive noise"""
    variables = (('x', 0),)
    parameters = (('alpha', 1), ('beta', 1))
    def drift(self, a, x, t):
        return - self.alpha * x
    def diffusion(self, b, x, t):
        return self.beta

class LinearMultiplicative(SODE):
    """Linear 1-D homogenous SODE with constant coefficients and
    multiplicative noise.

    dXt = alpha Xt dt + beta Xt dWt
    """
    variables = (('x', 1),)
    parameters = (('alpha', -0.0001), ('beta', 1))
    def drift(self, a, x, t):
        return self.alpha * x
    def diffusion(self, b, x, t):
        return self.beta * x
    def exact(self, x0, t, Wt):
        factor = self.alpha - 0.5 * self.beta ** 2
        return x0 * np.exp(factor * t + self.beta * Wt)

# Cannot currently integrate the one below, as it goes outside |x| <= 1

class SinusoidalMultiplicative(SODE):
    """Nonlinear 1-D SODE with solutions sinusoidal in Wt

    dXt = -(1/2) a^2 Xt dt + a sqrt(1 - Xt^2) dWt
    """
    variables = (('x', 0),)
    parameters = (('a', 0.01),)
    def drift(self, a, x, t):
        return - (self.a**2 * x) / 2
    def diffusion(self, b, x, t):
        return self.a * np.sqrt(1 - x**2)
    def exact(self, x0, t, Wt):
        return np.sin(self.a * Wt + np.arcsin(x0))

# Slightly easier one

class TangentMultiplicative(SODE):
    """Nonlinear 1-D SODE with solutions in tan(Wt)

    dXt = a^2 Xt (1 + Xt^2) dt + a (1 + Xt^2) dWt
    """
    variables = (('x', 0),)
    parameters = (('a', 0.01),)
    def drift(self, a, x, t):
        return self.a**2 * x * (1 + x**2)
    def diffusion(self, b, x, t):
        return self.a * (1 + x**2)
    def exact(self, x0, t, Wt):
        return np.tan(self.a*Wt + np.arctan(x0))

SYSTEMS = {
    None:'lin-2d',
    'weiner':WeinerSODE,
    'weiner2':WeinerSODE,
    'lin-add':LinearAdditive,
    'lin-multi':LinearMultiplicative,
    'sin-multi':SinusoidalMultiplicative,
    'tan-multi':TangentMultiplicative,
    'lin-2d':LinearAdditive2D,
}

if __name__ == "__main__":
    # Use Script to turn this into a script
    import sys
    from script import MultiScript

    # Actually run as a script
    script = MultiScript(SYSTEMS)
    script.main(argv=sys.argv[1:])

