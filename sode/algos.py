# (c) Oscar Benjamin, 2011, under terms of the new BSD license
#
# Author: Oscar Benjamin
# Date: 08 Mar 2011
#
# This module defines the SODE class for numerical solution of Stochastic
# Ordinary Differential Equations with additive noise.

from __future__ import division

import numpy as np
from numpy.random import standard_normal


def BrownianIncrements(dt, shape=1):
    """Function to generate an array of Brownian increments

    Usage:
    >>> # Create a single Brownian increment for a 1-d SODE
    >>> dt = 0.001
    >>> dWs = BrownianIncrements(dt)
    >>> # Create nsamples increments for an 1-d SODE covering a period of time
    >>> # given by dt * nsamples.
    >>> nsamples = 1000
    >>> dWs = BrownianIncrements(dt, nsamples)
    >>> # Create nsamples increments for an nvars-dimensional SODE
    >>> nvars = 4
    >>> dWs = BrownianIncrements(dt, (nsamples, nvars))

    Arguments:
    dt       : Sampling interval
    nvars    : Number of variables (width of array)
    'drift', 'diffusion'nsamples : Number of time samples

    Returns:
    dWs      : Array of shape (nvars, nsamples) with values drawn from a
               normal distribution with zero mean and dt variance
    """
    return np.sqrt(dt) * standard_normal(shape)


class BrownianMotion(object):
    """Discretised Brownian motion realisation

    The BrownianMotion class provides an abstraction for defining a particular
    realisation of a series of white noise increments.

    Instances can be used with the solve_bm method of SODE instances to
    generate a numerical solution from a known Broanian path. This can then be
    compared with the exact solution using the same Brownian path, e.g.:

    >>> # Create SODE subclass instance representing equations to solve
    >>> from sode.weiners import Weiner
    >>> from sode.algos import solve_bm
    >>> sodeinst = Weiner()
    >>> # Set up numeric parameters
    >>> x0 = sodeinst.get_x0() # or just: x0 = [0]
    >>> t1 = 0
    >>> t2 = 1
    >>> dt = 0.001
    >>> nsamples = int((t2 - t1)/dt)
    >>> # Create the BrownianMotion realisation and compute the exact solution
    >>> # corresponding to it.
    >>> bm = BrownianMotion(t1, dt, (nsamples, sodeinst.nvars))
    >>> xt_exact = sodeinst.exact(x0, bm.t, bm.Wt)
    >>> # Find numerical solution with same Brownian path
    >>> xt_numerical = solve_bm(sodeinst, x0, bm)
    >>> # Compare
    >>> print max(abs(xt_exact - xt_numerical).flat)
    0.0
    """
    def __init__(self, t1, dt, shape=None, _dWs=None):
        """Realisation of a Brownian process"""
        # Create increments or use provided
        if _dWs is None:
            if shape is None:
                raise ValueError('Need to provide shape')
            _dWs = BrownianIncrements(dt, shape)

        # Could generate an exception
        nsamples, nvars = _dWs.shape

        # Create t and Wt vectors
        Wt = np.zeros((nsamples + 1, nvars))
        Wt[1:, :] = np.cumsum(_dWs, axis=0)
        t = np.zeros_like(Wt)
        t.T[:] = t1 + dt * np.arange(nsamples + 1)

        # Attributes
        self.nvars = nvars
        self.nsamples = nsamples
        self.t = t
        self.dt = dt
        self.t1 = t1
        self.t2 = t1 + dt * nsamples
        self.dWs = _dWs
        self.Wt = Wt

    def coarse_grain(self):
        """Return a BrownianMotion instance describing a path that has been
        coarse grained by a factor of 2.

        >>> # Create a BrownianMotion instance
        >>> t1 = 0
        >>> dt = 0.001
        >>> nsamples = 1000
        >>> nvars = 4
        >>> bm = BrownianMotion(t1, dt, (nsamples, nvars))
        >>> print bm.t1, bm.t2, len(bm.t), bm.dt
        0 1.0 1001 0.001
        >>> # Create the corresponding coarse grain instance
        >>> bm2 = bm.coarse_grain()
        >>> print bm2.t1, bm2.t2, len(bm2.t), bm2.dt
        0 1.0 501 0.002
        """
        # Don't bother to deal with this case
        if self.nsamples % 2:
            raise ValueError("Need multiple of 2 sample to course grain")
        # Coarse grain increments
        dWs = (self.dWs[::2] + self.dWs[1::2])
        dt = 2 * self.dt
        nsamples = self.nsamples // 2
        # Return new instance
        return BrownianMotion(self.t1, dt, _dWs=dWs)


def solveEM(sodeinst, x1, t1, dt, dW):
    """Integrate SODEs numerically using Euler-Maruyama (EM) method.

    Integrate from ti to ti + dt with initial condition xi at ti in a
    single integration step (use solve() to generate a full timeseries)

    Usage:
    >>> from sode.weiners import Weiner
    >>> system = Weiner()
    >>> x1 = system.get_x0()
    >>> t1 = 0
    >>> dt = 0.001
    >>> dW = BrownianIncrements(dt, system.nvars)
    >>> x2, t2 = solveEM(system, x1, t1, dt, dW)

    Arguments:
    xi :    Vector specifying the state of the system at ti.
    ti :    Time corresponding to initial condition xi.
    dt :    Time step to use for EM integration.
    dW :    Change in Brownian motion over the interval ti to ti + dt
            e.diffusion.  dW = W(ti+dt) - W(ti). If W is a Brownian motion
            then dW should be a normally distributed random variable with
            variance dt.

    Returns:
    x2 :    Vector specifying the state of the system at t2
    t2 :    t1 + dt, returned for iterative convenience

    Description:
    The Euler Maruyama method approximates the stochastic integration step
    from Tn to Tn+1

                      / Tn+1              / Tn+1
    X(Tn+1) = X(Tn) + |    a(X(s), s)ds + |    b(X(s), s)dW(s)
                      / Tn                / Tn

    with the approximation

    Xn+1 = Xn + a(Xn, tn) delta_t + b(Xn, tn) delta_Wn

    where delta_Wn is a normally-distributed random variable with expected
    value 0 and variance delta_t. The method is a simple extension of
    Euler's method for deterministic ODEs.
    """
    a = np.zeros(sodeinst.nvars)
    b = np.zeros(sodeinst.nvars)
    a = sodeinst.drift(a, x1, t1)
    b = sodeinst.diffusion(b, x1, t1)
    x2 = x1 + a * dt + b * dW
    t2 = t1 + dt
    return x2, t2

def solveRK4_additive(sodeinst, x1, t1, dt, dW):
    """Integrate SODEs numerically using 4-step Runge-Kutta method.

    Integrate additive SODE from ti to ti + dt with initial condition x1
    at t1 in a using a 4 step Runge-Kutta type method.

    This method is from:
    "Numerical methods for stochastic differential equations" by Joshua
    Wilkie and published in Physical Review E. 2004.

    Usage:
    >>> from sode.weiners import Weiner
    >>> from sode.algos import solveRK4_additive
    >>> system = Weiner()
    >>> x1 = system.get_x0()
    >>> t1 = 0
    >>> dt = 0.001
    >>> dW = BrownianIncrements(dt, system.nvars)
    >>> x2, t2 = solveRK4_additive(system, x1, t1, dt, dW)

    Arguments:
    x1 :    Vector specifying the state of the system at t1.
    t1 :    Time corresponding to initial condition x1.
    dt :    Time step to use for EM integration.
    dW :    Change in Brownian motion over the interval t1 to t1 + dt
            e.diffusion.  dW = W(t1+dt) - W(t1). If W is a Brownian motion
            then dW should be a normally distributed random variable with
            variance dt.

    Returns:
    x2 :    Vector specifying the state of the system at t2
    t2 :    t1 + dt, returned for iterative convenience

    Description:
    We wich to integrate Xn from Tn to Tn+1:

                      / Tn+1              / Tn+1
    X(Tn+1) = X(Tn) + |    a(X(s), s)ds + |    b(X(s), s)dW(s)
                      / Tn                / Tn

    First define the vector function F as:

    Fj(Xn, Tn) = aj(Xn, tn) delta_t + bj(Xn, tn) delta_Wnj

    Then define successive approximations:

    K1 = F(Xn, Tn)
    K2 = F(Xn + 1/2 K1, Tn + 1/2 delta_t)
    K3 = F(Xn + 1/2 K2, Tn + 1/2 delta_t)
    k4 = F(X + K3, Tn + delta_t)

    We can take a weighted average of these increments:

    Xn+1 = Xn + 1/6 (K1 + 2 K2 + 2 K3 + K4)

    The method is a simple extension of the RK4 method for ODEs. However,
    in the case of SDEs convergence is only of order 2.

    The form used here is valid only for additive noise.
    """
    a = np.zeros(sodeinst.nvars)
    b = np.zeros(sodeinst.nvars)
    K1 = _rk4_inc(sodeinst, a, b, x1,         t1        , dt, dW)
    K2 = _rk4_inc(sodeinst, a, b, x1 + K1/2., t1 + dt/2., dt, dW)
    K3 = _rk4_inc(sodeinst, a, b, x1 + K2/2., t1 + dt/2., dt, dW)
    K4 = _rk4_inc(sodeinst, a, b, x1 + K3   , t1 + dt   , dt, dW)
    x2 = x1 + (K1 + 2*K2 + 2*K3 + K4) / 6.
    t2 = t1 + dt
    return x2, t2

def _rk4_inc(sodeinst, a, b, x1, t1, dt, dW):
    return sodeinst.drift(a, x1, t1) * dt + sodeinst.diffusion(b, x1, t1) * dW

def _elementary_method(method):
    method = method or 'EM'
    if method == 'EM':
        return solveEM
    elif method == 'RK4':
        return solveRK4_additive
    else:
        raise ValueError("Unrecognised method '{0}'".format(method))

def largest_dt(T, dtmax):
    """Compute largest dt that gives equal size steps from 0 to T"""
    # Awkward because numpy does not provide round-away-from-zero
    r = T / dtmax
    N = int(np.sign(r) * np.ceil(abs(r)))
    dt = T / N
    return N, dt

def solve(sodeinst, x0, t, dtmax=0.001, method=None):
    """Integrate SODEs numerically to obtain solution at times t.

    Usage:
    >>> from sode.weiners import Weiner
    >>> from sode.algos import solve
    >>> system = Weiner()
    >>> x1 = system.get_x0()
    >>> t = np.arange(0, 1, 1e-2)
    >>> Xt = solve(system, x1, t, dtmax=1e-3)

    Arguments:
    x0    : Vector specifying the state of the system as t[0]. get_x0()
            returns the default initial conditions defined for the
            equations.
    t     : Vector of times at which the solution is desired to be known.
    dtmax : Maximum integration step to use (default 0.001).
    method: 'EM' or 'RK4'. Specifies the elementary integration method.
            'RK4' is valid only for SODEs with additive noise.

    Returns:
    Xt :    2-dimensional array specifying the states of the system at
            times in t. The nth rowth of X (X[n, :]) gives the state of
            the system at time t[n].
    """
    # Choose method
    method = _elementary_method(method)

    # Prepare data for solution
    Nequations = len(x0)
    Ntimes = len(t)
    Xt = np.zeros((Ntimes, Nequations))
    Xt[:] = np.nan

    # Iteratively integrate from x[n], t[n] to x[n+1] t[n+1]
    Xt[0, :] = xi = x0

    if hasattr(sodeinst, 'solveto'):
        for n in range(Ntimes - 1):
            sodeinst.solveto(Xt[n, :], t[n], Xt[n+1,:], t[n+1], dtmax)
        return Xt

    for n in range(Ntimes - 1):

        # Break into substeps and generate BrownianIncrements
        Nsteps, dt = largest_dt(t[n+1] - t[n], dtmax)
        dWs = BrownianIncrements(dt, (Nsteps, Nequations))

        # Iterate over increments with fixed step
        xi, ti = Xt[n, :], t[n]
        for dW in dWs:
            xi, ti = method(sodeinst, xi, ti, dt, dW)


        # ti is incremented from t[n] by dt Nsteps times
        if not np.allclose(ti, t[n+1]):
            raise ValueError("Stepsize underflow")

        # Stop once nans start appearing
        if np.isnan(xi).any():
            warnings.warn('nans in solution. Stopping')
            break

        # Store this state
        Xt[n+1, :] = xi

    # Return array of states
    return Xt

def solve_bm(sodeinst, x0, bm, method=None):
    """Integrate SODEs with fixed step and provided noise realisation.

    Usage:
    >>> # Create an SODE instance
    >>> from sode.weiners import Weiner
    >>> from sode.algos import solve_bm
    >>> system = Weiner()
    >>> t1 = 0
    >>> t2 = 1
    >>> dt = 0.001
    >>> nsamples = int((t2 - t1) / dt)
    >>> bm = BrownianMotion(t1, dt, (nsamples, system.nvars))
    >>> x2 = solve_bm(system, system.get_x0(), bm)

    Arguments:
    x0 :    Vector specifying the state of the system as t[0]. get_x0()
            returns the default initial conditions defined for the
            equations.
    bm :    BrownianMotion instance. Integration is performed with the
            increments from bm. bm is also used to establish the interval
            of integration and the integration step.
    method: 'EM' or 'RK4'. Specifies the elementary integration method.
            'RK4' is valid only for SODEs with additive noise.


    Returns:
    t  :    Vector of time values for the solution
    Xt :    2-dimensional array specifying the state of the system at
            times corresponding to t.

    Notes:
    The integration timestep, dt, is determined by:
    dt = (t2 - t1) / (len(W) - 1)
    """
    # Choose method
    method = _elementary_method(method)

    # Prepare data for solution
    Nequations = len(x0)
    Ntimes = bm.nsamples + 1
    Xt = np.zeros((Ntimes, Nequations))
    Xt[:] = np.nan

    # Iterate over increments with fixed step
    Xt[0, :], ti = x0, bm.t1
    for n, dW in enumerate(bm.dWs):
        Xt[n+1, :], ti = method(sodeinst, Xt[n, :], ti, bm.dt, dW)

    # ti is incremented from t[n] by dt Ntimes times
    if not np.allclose(ti, bm.t2):
        raise ValueError("Stepsize underflow")

    # Return states and times
    return Xt


if __name__ == "__main__":
    import doctest
    doctest.testmod()
