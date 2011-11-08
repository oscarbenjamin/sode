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
    >>> from sode.pysode import Weiner
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
    >>> xt_numerical = sodeinst.solve_bm(x0, bm)
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


if __name__ == "__main__":
    import doctest
    doctest.testmod()
