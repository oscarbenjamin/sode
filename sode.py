#
# Author: Oscar Benjamin
# Date: 08 Mar 2011
#
# This module defines the SODE class for numerical solution of Stochastic
# Ordinary Differential Equations with additive noise.

from __future__ import division

import datetime
import re
import itertools
import warnings
import sys

import numpy as np
from numpy.random import randn

from opster import command


class SODESubclassError(NotImplementedError):
    """Exception raised by subclasses of SODE when the required subclass
    interface has not been defined"""
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


def BrownianIncrements(dt, nsamples=1, nvars=1):
    """Function to generate an array of Brownian increments

    Usage:
    >>> # Create a single Brownian increment for a 1-d SODE
    >>> dWs = BrownianIncrements(dt)
    >>> # Create nsamples increments for an 1-d SODE covering a period of time
    >>> # given by dt * nsamples.
    >>> dWs = BrownianIncrements(dt, nsamples)
    >>> # Create nsamples increments for an nvars-dimensional SODE
    >>> dWs = BrownianIncrements(dt, nsamples, nvars)

    Arguments:
    dt       : Sampling interval
    nvars    : Number of variables (width of array)
    'drift', 'diffusion'nsamples : Number of time samples

    Returns:
    dWs      : Array of shape (nvars, nsamples) with values drawn from a
               normal distribution with zero mean and dt variance
    """
    return np.sqrt(dt) * randn(nsamples, nvars)


class BrownianMotion(object):
    """Discretised Brownian motion realisation

    The BrownianMotion class provides an abstraction for defining a particular
    realisation of a series of white noise increments.

    Instances can be used with the solve_bm method of SODE instances to
    generate a numerical solution from a known Broanian path. This can then be
    compared with the exact solution using the same Brownian path, e.g.:

    >>> # Create Brownian path and determine exact solution
    >>> bm = BrownianMotion(t1, N, dt)
    >>> Xexact = system.exact(x0, bm.t, bm.Wt)
    >>> # Find numerical solution with same Brownian path
    >>> Xnumerical, t = system.solve_bm(x0, bm)
    >>> # Compare
    >>> print max(abs(Xexact - Xnumerical))
    """
    def __init__(self, t1, dt, nsamples=None, nvars=None, dWs=None):
        """Realisation of a Brownian process"""
        # Create increments or use provided
        if dWs is not None:
            if nvars is None:
                nvars = dWs.shape[1]
            if nsamples is None:
                nsamples = dWs.shape[0]
            if dWs.shape != (nsamples, nvars):
                raise ValueError("dWs dimensions != (nsamples, nvars)")
        else:
            dWs = BrownianIncrements(dt, nsamples, nvars)

        # Create t and Wt vectors
        t = t1 + dt * np.arange(nsamples + 1)
        Wt = np.zeros((nsamples + 1, nvars))
        Wt[1:, :] = np.cumsum(dWs, axis=0)

        # Attributes
        self.nvars = nvars
        self.nsamples = nsamples
        self.t = t
        self.dt = dt
        self.t1 = t[0]
        self.t2 = t[-1]
        self.dWs = dWs
        self.Wt = Wt

    def coarse_grain(self):
        """Return a BrownianMotion instance describing a path that has been
        coarse grained by a factor of 2."""
        # Don't bother to deal with this case
        if self.nsamples % 2:
            raise ValueError("Need multiple of 2 sample to course grain")
        # Coarse grain increments
        dWs = (self.dWs[::2] + self.dWs[1::2])
        dt = 2 * self.dt
        nsamples = self.nsamples // 2
        # Return new instance
        return BrownianMotion(self.t1, dt, dWs=dWs)


class SODE(object):
    """The SODE class represents a set of Stochastic Ordinary Differential
    Equations. The SODEs take the form

        dXt = drift(Xt, t)dt + diffusion(Xt, t)dWt

    where X is a vector, drift and diffusion are vector fields and W is a
    vector of Wiener processes, and the equation is interpreted in the Ito
    sense.

    Subclasses should define the functions drift(x, t), diffusion(x, t),
    parameters and variables. drift and diffusion should accept both arguments
    (x, t) where t is scalar and x is a vector. drift and diffusion should
    return vectors of the same size as x.
    """

    # List of regexes that cannot be used for variable names or parameter
    # names. These are reserved because they cannot be used as attributes
    _RES = ['drift', 'diffusion', 'exact']
    _RES_BEG = ['set', 'get', 'has', 'solve', 'make', 'load', 'save']
    _RE_RESERVED = [
        (r'^_.*$',
            "Names must not being with underscore: '{0}'"),
        (r'^({0})$'.format('|'.join(_RES)),
            "'{0}' is a reserved name"),
        (r'^({0})_.*$'.format('|'.join(_RES_BEG)),
            "Names ('{{0}}') cannot begin with {0}".format(','.join(_RES_BEG))),
    ]

    def __init__(self, **kwargs):
        """SODE(**kwargs)

        kwargs: dict specifying initial conditions and parameter values.

        Initial conditions should use variable name followed by '0'.
        Value must be understandable to the float() function.
        Parameters and initial conditions not supplied will have the default
        values from self.parameters and self.variables if not provided in
        kwargs.
        """
        # Check the appropriate variables are defined
        if not (hasattr(self, 'variables') and hasattr(self, 'parameters')):
            msg = "Subclasses should define parameters and variables"
            raise SODESubclassError(msg)

        # Check for reserved or duplicated parameter/variable names, establish
        # number of variables, and create x0
        self._init_pv(self.parameters, self.variables)

        # Read off parameter values specified in kwargs
        self._init_kwargs(**kwargs)

    #
    # The functions below define the SODE system. Since this is an abstract
    # class, we do not define them here.
    #

    # drift coefficient (deterministic derivative)
    def drift(self, a, x, t):
        raise SODESubclassError("Subclasses should override this function")

    # diffusion coefficient (noise amplitude)
    def diffusion(self, b, x, t):
        raise SODESubclassError("Subclasses should override this function")

    # exact solution (used in convergence calculations if provided
    def exact(self, x0, t, Wt):
        raise SODESubclassError("Subclasses should override this function")

    # Code for accessing variables:

    def _init_pv(self, parameters, variables):
        # Check for reserved or duplicated parameter/variable names
        vp_names = set()
        for name, _ in variables + parameters:
            for r, msg in self._RE_RESERVED:
                if re.match(r, name):
                    raise SODESubclassError(msg.format(name))
            if name in vp_names:
                raise ValueError("Duplicate name '{0}'".format(name))
            vp_names.add(name)

        # Prepare scene for variable/parameter access
        # _x0 is the array holding the initial conditions
        # Variable indices are stored as attributes with the same names
        self._variables = [var for var, _ in variables]
        self._parameters = [par for par, _ in parameters]
        self._x0 = np.zeros(len(self._variables))
        for n, var in enumerate(self._variables):
            self._set_var_index(var, n)
        # For convenience access later
        self.nvars = len(self._variables)

        # Initialise default values for parameters and variables
        for var, default in variables:
            self.set_ic(var, default)
        for par, default in parameters:
            self.set_parameter(par, default)

    def _init_kwargs(self, **kwargs):
        # Process parameters and ics provided in kwargs:
        for key, val in kwargs.iteritems():
            # Initial condition or parameter ?
            if key.endswith('0'):
                self.set_ic(key[:-1], val)
            else:
                try:
                    self.set_parameter(key, val)
                except ValueError:
                    pass

    def _shift_indices(self, offset, x0):
        # Copy old initial conditions and update indices
        for var in self._variables:
            old_index = self._get_var_index(var)
            new_index = old_index + offset
            x0[new_index] = self._x0[old_index]
            self._set_var_index(var, new_index)
        # Prevent this system from being used independently
        self._x0 = x0

    def get_variables(self):
        """Return a list of variable names"""
        return list(self._variables)

    def _varchk(self, var):
        if not var in self.get_variables():
            raise ValueError("No such variable {0}".format(var))

    def _get_var_index(self, var):
        return getattr(self, var)

    def _set_var_index(self, var, index):
        self._varchk(var)
        setattr(self, var, index)

    def get_ic(self, var):
        """Returns initial condition. Raises ValueError if unrecognised"""
        self._varchk(var)
        return self._x0[self._get_var_index(var)]

    def set_ic(self, var, val):
        """Set initial condition for var to float(val)"""
        self._varchk(var)
        self._x0[self._get_var_index(var)] = val

    def get_x0(self):
        """Return the initial condition as a vector suitable for drift and diffusion"""
        return self._x0.copy()

    # Code for accessing parameters

    def get_parameters(self):
        """Return a list of parameter names"""
        return list(self._parameters)

    def _parchk(self, par):
        if not par in self.get_parameters():
            raise ValueError("No such parameter {0}".format(par))

    def get_parameter(self, par):
        """Returns parameter value. Raises ValueError if unrecognised"""
        self._parchk(par)
        return getattr(self, par)

    def set_parameter(self, par, val):
        """Set parameter par to float(val)"""
        self._parchk(par)
        setattr(self, par, float(val))

    #
    # Pretty-print system
    #

    def get_description(self):
        """Print parameter values and initial condition to stdout"""
        lines = [self.__doc__]
        for name in self.get_parameters():
            lines.append('par {0} = {1}'.format(name, self.get_parameter(name)))
        for name in self.get_variables():
            lines.append('var {0} = {1}'.format(name, self.get_ic(name)))
        return '\n'.join(lines)

    #
    # Numerical routines
    #

    def largest_dt(self, T, dtmax):
        """Compute largest dt that gives equal size steps from 0 to T"""
        # Awkward because numpy does not provide round-away-from-zero
        r = T / dtmax
        N = int(np.sign(r) * np.ceil(abs(r)))
        dt = T / N
        return N, dt

    def solveEM(self, x1, t1, dt, dW):
        """Integrate SODEs numerically using Euler-Maruyama (EM) method.

        Integrate from ti to ti + dt with initial condition xi at ti in a
        single integration step (use solve() to generate a full timeseries)

        Usage:
        >>> x2, t2 = self.solveEM(x1, t1, dt, dW)

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
        a = np.zeros(self.nvars)
        b = np.zeros(self.nvars)
        a = self.drift(a, x1, t1)
        b = self.diffusion(b, x1, t1)
        x2 = x1 + a * dt + b * dW
        t2 = t1 + dt
        return x2, t2

    def solveRK4_additive(self, x1, t1, dt, dW):
        """Integrate SODEs numerically using 4-step Runge-Kutta method.

        Integrate additive SODE from ti to ti + dt with initial condition x1
        at t1 in a using a 4 step Runge-Kutta type method.

        This method is from:
        "Numerical methods for stochastic differential equations" by Joshua
        Wilkie and published in Physical Review E. 2004.

        Usage:
        >>> x2, t2 = self.solveRK4_additive(x1, t1, dt, dW)

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
        a = np.zeros(self.nvars)
        b = np.zeros(self.nvars)
        K1 = self._rk4_inc(a, b, x1,         t1        , dt, dW)
        K2 = self._rk4_inc(a, b, x1 + K1/2., t1 + dt/2., dt, dW)
        K3 = self._rk4_inc(a, b, x1 + K2/2., t1 + dt/2., dt, dW)
        K4 = self._rk4_inc(a, b, x1 + K3   , t1 + dt   , dt, dW)
        x2 = x1 + (K1 + 2*K2 + 2*K3 + K4) / 6.
        t2 = t1 + dt
        return x2, t2

    def _rk4_inc(self, a, b, x1, t1, dt, dW):
        return self.drift(a, x1, t1) * dt + self.diffusion(b, x1, t1) * dW

    def _elementary_method(self, method):
        method = method or 'EM'
        if method == 'EM':
            return self.solveEM
        elif method == 'RK4':
            return self.solveRK4_additive
        else:
            raise ValueError("Unrecognised method '{0}'".format(method))

    def solve(self, x0, t, dtmax=0.001, method=None):
        """Integrate SODEs numerically to obtain solution at times t.

        Usage:
        >>> Xt = sys.solve(x0, t, dtmax)
        >>> Xt = sys.solve(x0, t)

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
        method = self._elementary_method(method)

        # Prepare data for solution
        Nequations = len(x0)
        Ntimes = len(t)
        Xt = np.zeros((Ntimes, Nequations))
        Xt[:] = np.nan

        # Iteratively integrate from x[n], t[n] to x[n+1] t[n+1]
        Xt[0, :] = xi = x0
        for n in range(Ntimes - 1):

            # Break into substeps and generate BrownianIncrements
            Nsteps, dt = self.largest_dt(t[n+1] - t[n], dtmax)
            dWs = BrownianIncrements(dt, Nsteps, Nequations)

            # Iterate over increments with fixed step
            xi, ti = Xt[n, :], t[n]
            for dW in dWs:
                xi, ti = method(xi, ti, dt, dW)


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

    def solve_bm(self, x0, bm, method=None):
        """Integrate SODEs with fixed step and provided noise realisation.

        Usage:
        >>> x2 = sys.solve(x1, t1, t2, W)

        Arguments:
        x0 :    Vector specifying the state of the system as t[0]. get_x0()
                returns the default initial conditions defined for the
                equations.
        t1 :    Time corresponding to initial condition x0.
        t2 :    Time at end of integration.
        W  :    2-dimensional array specifying the values of the Brownian
                motions at times [t1, t1+dt, d0+2*dt, ...]. Each row of W,
                W[n, :], represents the value of the Brownian motions at time
                t[n].
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
        method = self._elementary_method(method)

        # Prepare data for solution
        Nequations = len(x0)
        Ntimes = bm.nsamples + 1
        Xt = np.zeros((Ntimes, Nequations))
        Xt[:] = np.nan

        # Iterate over increments with fixed step
        Xt[0, :], ti = x0, bm.t1
        for n, dW in enumerate(bm.dWs):
            Xt[n+1, :], ti = method(Xt[n, :], ti, bm.dt, dW)

        # ti is incremented from t[n] by dt Ntimes times
        if not np.allclose(ti, bm.t2):
            raise ValueError("Stepsize underflow")

        # Return states and times
        return Xt

    def make_brownian_path(self, t1, dt, nsamples):
        """Create a discrete Brownian path suitable for use with solve_bm"""
        return BrownianMotion(t1, dt, nsamples, self.nvars)

    @staticmethod
    def load_csv(fin, parameters=False):
        """Load solution t, Xt from input-file fin

        Usage:
        >>> t, Xt = sys1.load_csv(open('sys1_1.txt', 'r'))
        >>> t, Xt = sys1.load_csv(open('sys1_1.txt', 'r'), parameters=True)

        If parameters is True the file will be checked for the appropriate
        header information. It will be checked for consistency against this
        SODE instance. In particular it must have been saved by an instance of
        the same SODE subclass. Parameters will also be loaded from the file
        into this instance so that it has the same parameter values as used by
        the instance that saved the solution.

        Otherwise, the file is treated as a csv file where lines beginning
        with '#' are comments. The first non-comment line contains the
        (comma-separated) column titles, starting with 'time' and then each of
        the variables in order.  Subsequent lines contain comma-seperated
        values, e.g.:

        # COMMENTS GO HERE
        # COMMENTS HERE
        time, x, y
        0.0, 1.0, 0.0
        0.1, 0.9976547, 0.0132455
        """
        # Doesn't accept string input unless split into lines
        if isinstance(fin, str):
            fin = open(fin, 'r')

        # Use generator to remove newlines if present
        def iter_f():
            for line in fin:
                if line.endswith('\n'):
                    yield line[:-1]
                else:
                    yield line

        # Create iterator and draw first line
        itf = iter_f()
        line = itf.next()

        # Ignore comment lines
        comment_lines = []
        while line.startswith('#'):
            comment_lines.append(line)
            line = itf.next()

        # Parse variable names
        names = line.split(', ')[1:]
        nvars = len(names)

        # Read parameter values from comment lines
        if parameters:
            ret = SODE.load_parameters(comment_lines)

        # Read the rest of the file
        lines = list(itf)

        # Iterate through parsing into arrays
        t = np.zeros((len(lines)))
        Xt = np.zeros((len(lines), nvars))
        for n, l in enumerate(lines):
            nums = [float(s) for s in l.split(', ')]
            t[n], Xt[n, :] = nums[0], nums[1:]

        # And return the result
        if parameters:
            return t, Xt, ret
        else:
            return t, Xt

    def save_csv(self, t, Xt, fout, header=True, titles=True):
        """Save solution t, Xt to output-file fout

        Usage:
        >>> sys1.save_csv(t, Xt, open('sys1_1.txt', 'w'))

        Writes the solution to fout in csv format where lines beginning with
        '#' are comments. Lines beginning with '##' are used to write
        meta-data used by load_csv. Example:

        ### SODE VERSION 0
        ###.{}
        ###.[]
        ## 1-D linear SODE with equation:
        ##
        ##         dx = - beta x dt + alpha dW
        ##
        ## par beta = 1.0
        ## par alpha = 1.0
        ## var x = 1.0
        # Created: 2011-03-12 18:42:32.040680
        time, x
        0.0, 1.0
        0.1, 1.23475373685
        """
        # First write in comment headers
        if header:
            # Magic line
            fout.write('### SODE VERSION 0\n')

            # _make_sode in Script adds these attributes
            if not hasattr(self, '_sys_opts'):
                sys_opts, sys_args = {}, []
            else:
                sys_opts, sys_args = self._sys_opts
            fout.write('##.{0}\n'.format(repr(sys_opts)))
            fout.write('##.{0}\n'.format(repr(sys_args)))

            # Write description
            for line in self.get_description().split('\n'):
                fout.write('## {0}\n'.format(line))
            dt = datetime.datetime.now()
            fout.write('# Created: {0} by {1}\n'.format(dt, sys.argv[0]))

        if titles:
            # Write column header
            cols = ['time'] + self.get_variables()
            fout.write(', '.join(cols) + '\n')

        # write raw numbers (csv format)
        for ti, xi in itertools.izip(t, Xt):
            nums = [str(ti)] + [str(x) for x in xi]
            fout.write(', '.join(nums) + '\n')

    # Regexes for parsing solution file
    _MAGIC_LINE = '^### SODE .*$'
    _PV_LINE = '^## ([pv]ar) (\w+) = ([0-9.-]+)'

    @staticmethod
    def load_parameters(fin):
        """Read parameter header information from saved file.

        The information should be formatted as by save_csv"""
        # Create line iterator
        if isinstance(fin, str):
            fin = open(fin, 'r')
        itf = iter(fin)
        line = itf.next()

        # The data we want
        sys_opts = {}
        sys_args = []

        # Check first line
        if not re.match(SODE._MAGIC_LINE, line):
            raise ValueError("Invalid solution file for parameters")
        line = itf.next()

        # See if we have sys_opts
        if line.startswith('##.'):
            sys_opts = eval(line[3:])
            line = itf.next()

        # Attempt to read header
        while line.startswith('##'):
            m = re.match(SODE._PV_LINE, line)
            if m:
                pv, name, val = m.groups()
                val = float(val)
                if pv == 'par':
                    a = "{0}={1}".format(name, val)
                elif pv == 'var':
                    a = "{0}0={1}".format(name, val)
                sys_args.append(a)
            line = itf.next()

        # Return the parse data
        return sys_opts, sys_args


class SODENetwork(SODE):
    """Specialises SODE for a system defined in terms of subsystems"""

    def __init__(self):
        """Initialise before adding subsystems"""
        # Empty subsystem data
        self._subsystems = [] # Remember order
        self._subsysdict = {} # Convenient access (want a sorteddict)

        # Initialise own parameters, variables with no kwargs
        super(SODENetwork, self).__init__()

    def init(self, **kwargs):
        """Flattens subsystems into single x0 and shifts indices"""

        # Determine indices of subsystems within larger state vector
        ntotal = self.nvars
        indices = []
        for name, sys in self._subsystems:
            indices.append((sys, ntotal))
            ntotal += sys.nvars

        # Create larger x0, and copy own values to the beginning
        x0 = np.zeros(ntotal)
        x0[:self.nvars] = self._x0
        self._x0 = x0
        self.nvars = ntotal

        # Now make all subsystems use a view into the larger array
        # Also copys subsystems values into the array
        for sys, n1 in indices:
            sys._shift_indices(n1, x0)

        # Compile list to be used in drift_subsys and diffusion_subsys
        # zero-d systems (just parameters) should be ignored for efficiency
        self._ss_eval = []
        for name, sys in self._subsystems:
            if sys.nvars:
                self._ss_eval.append(sys)

        # Now read off kwargs
        self._init_kwargs(**kwargs)

    # Add subsystem is called before __init__ so need to check everything here
    def add_subsystem(self, name, subsys):
        """Add a subsystem"""
        # Subsystem cannot be owned twice
        if hasattr(subsys, '_owned'):
            msg = "Cannot add already owned subsystem '{0}'".format(name)
            raise ValueError(msg)
        subsys._owned = self

        # Parameter access is non-unique if there are duplicate names
        if name in self._subsysdict:
            msg = "Duplicate subsystem name '{0}'".format(name)
            raise ValueError(msg)

        # Store subsystems and names
        self._subsystems.append((name, subsys))
        self._subsysdict[name] = subsys

    def _shift_indices(self, n, x0):
        # Shift child indices as well as own
        super(SODENetwork, self)._shift_indices(n, x0)
        for name, sys in self._subsystems:
            sys._shift_indices(n, x0)

    def get_parameters(self):
        # Retrieve own parameters
        pars = super(SODENetwork, self).get_parameters()
        # Retrieve subsystem parameters, with name prepended
        for name, sys in self._subsystems:
            for par in sys.get_parameters():
                pars.append('{0}.{1}'.format(name, par))
        return pars

    def get_variables(self):
        # Retrieve own parameters
        vars_ = super(SODENetwork, self).get_variables()
        # Retrieve subsystem parameters, with name prepended
        for name, sys in self._subsystems:
            for var in sys.get_variables():
                vars_.append('{0}.{1}'.format(name, var))
        return vars_

    # Use this to override get_/set_ par ic methods

    def _get_subsystem(self, pv_name):
        names = pv_name.split('.')
        snames, pv = names[:-1], names[-1]
        sys = self
        for sn in snames:
            if not hasattr(sys, '_subsysdict'):
                raise ValueError("No such subsystem '{0}'".format(sn))
            sys = sys._subsysdict[sn]
        return sys, pv

    def get_parameter(self, par):
        """Returns parameter value. Raises ValueError if unrecognised"""
        ss, par = self._get_subsystem(par)
        return SODE.get_parameter(ss, par)

    def set_parameter(self, par, val):
        """Set initial condition for var to float(val)"""
        ss, par = self._get_subsystem(par)
        return SODE.set_parameter(ss, par, val)

    def get_ic(self, var):
        """Returns initial condition. Raises ValueError if unrecognised"""
        ss, var = self._get_subsystem(var)
        return SODE.get_ic(ss, var)

    def set_ic(self, var, val):
        """Set initial condition for var to float(val)"""
        ss, var = self._get_subsystem(var)
        return SODE.set_ic(ss, var, val)

    # Subclasses must call these. They use self._ss_eval for efficiency

    def drift_subsys(self, a, x, t):
        """Compute drift for subsystems. Subclasses must call this in drift"""
        for sys in self._ss_eval:
            sys.drift(a, x, t)
        return a

    def diffusion_subsys(self, b, x, t):
        """Compute diffusion for subsystems. Subclasses must call this in
        diffusion"""
        for sys in self._ss_eval:
            sys.diffusion(b, x, t)
        return b

# A Script instance is a callable that can be used to quickly create a
# script to investigate a particular SODE subclass.
#
# There are two ways to use Script. One is to create a Script instance with
# an SODE subclass as argument. The other is to override the make_sode
# method. This method is called with all command line positional and option
# arguments as *args and **kwargs respectively. These two methods are
# illustrated below.
#
# Example script:
# -----------------------------------
# Define a simple system
# class MySODE(SODE):
#   def f(self, x, t):
#       return x
#   def g(self, x, t):
#       return t
#
# Use script main if run as a script
# if __name__ == "__main__":
#    import sys
#    from script import Script
#    script = Script(MySODE)
#    script(argv=sys.argv[1:])
## -----------------------------------
#
#A more complicated example:
## -----------------------------------
## Define several systems
#class MySODE1(SODE):
#    ...
#class MySODE2(SODE):
#    ...
#
## Collect together in a dict
#SYSTEMS = {'sys1':MySODE1, 'sys2':MySODE2, ...}
#
#if __name__ == "__main__":
#    import sys
#    from script import Script
#
#    # Override Script to choose system type
#    class MyScript(Script):
#        usage = '%name [OPTS] SYSNAME [PAR1=VAL1 ...]
#        def make_sode(self, sysname, *args, **opts):
#            if sysname not in SYSTEMS:
#                msg = "System name should be one of {0}"
#                raise ValueError(msg.format(', '.join(SYSTEMS)))
#            return SYSTEMS[sysname](*args)
#
#    # Actually run as a script
#    script = MyScript()
#    script(argv=sys.argv[1:])
## -----------------------------------


class Script(object):
    """Class for quickly makeing scripts out of SODE objects"""
    def __init__(self, SYSTYPE=None):
        """Stores SYSTYPE as the SODE subclass for this instance"""
        # Store SYSTYPE if provided otherwise assume that make_sode has been
        # overridden and use its docstring
        if SYSTYPE is not None:
            self.SYSTYPE = SYSTYPE

        self.cmdtable = {
            'info'  : (self.info  , self.info_opts  , 'FILE'            ),
            'params': (self.params, self.params_opts, '[SYSARGS]'       ),
            'solve' : (self.solve , self.solve_opts , '[OPTS] [SYSARGS]'),
            'cont'  : (self.cont  , self.cont_opts  , '[OPTS] FILE'     ),
            'plot'  : (self.plot  , self.plot_opts  , '[OPTS] FILE'     ),
            'conv'  : (self.conv  , self.conv_opts  , '[OPTS] [SYSARGS]'),
        }

        for cmd in ['params', 'solve', 'conv']:
            f, opts, u = self.cmdtable[cmd]
            self.cmdtable[cmd] = f, opts + self.sys_opts, u

    def main(self, argv):
        """Parse argv and run as command line script"""
        from opster import dispatch
        dispatch(argv, self.cmdtable)

    # Subclasses can use sys_opts to define options passed to make_sode
    sys_opts = []

    plot_opts = [
        ('p', 'plot', False, 'Plot to screen'),
        ('P', 'plot-file', '', 'File to save plot image to'),
    ]
    solve_opts = sys_opts + plot_opts + [
        ('t', 't1', 0.0, 'time at beginning of solution'),
        ('T', 't2', 1.0, 'time at end of solution'),
        ('d', 'dtmax', 0.001, 'integration step'),
        ('D', 'dtout', 0.01,  'interval between solution samples'),
        ('o', 'output-file', '',  'Output file for solution data'),
        ('' , 'params', False,  'print parameters and exit'),
        ('m', 'method', 'EM',  'integration method'),
    ]
    def solve(self, *args, **opts):
        """Solve numerically using random realisation of noise

        If neither of plot or outputfile is specified, the solution will be
        written to stdout.
        """
        # Default to stdout
        if not (opts['output_file'] or opts['plot'] or opts['plot_file']):
            opts['output_file'] = '-'

        # Create system with args
        sysinst = self._make_sode(*args, **opts)
        if not sysinst:
            return 1
        if opts['params']:
            print sysinst.get_description()
            return 0

        # Create solution
        t = np.arange(opts['t1'], opts['t2'] + opts['dtout'], opts['dtout'])
        x0 = sysinst.get_x0()
        Xt = sysinst.solve(x0, t, opts['dtmax'], method=opts['method'])

        # Ouput file is contained in opts (no-op if not provided)
        self.save_solution(sysinst, t, Xt, **opts)

        # Plot (no-op) if plot options not provided
        fig = self.figure(**opts)
        if fig:
            ax = fig.add_subplot(1, 1, 1)
            self.plot_solution(ax, t, Xt, 'k')
            self.show(fig, **opts)

    cont_opts = [
        ('T', 'T', 10.0, 'time to integrate before saving'),
        ('d', 'dtmax', 0.001, 'integration step'),
        ('D', 'dtout', 0.01,  'interval between solution samples'),
        ('' , 'params', False,  'print parameters and exit'),
        ('m', 'method', 'EM',  'integration method'),
        ('s', 'stdout', False,  'Output file for solution data'),
    ]
    def cont(self, file_name, **opts):
        """Continue previous numerical solution from input file"""
        # Load solution object
        t, Xt, (sysopts, sysargs) = SODE.load_csv(file_name, parameters=True)
        sysinst = self._make_sode(*sysargs, **sysopts)
        if not sysinst:
            return -1
        if opts['params']:
            print sysinst.get_description()
            return 0

        if opts['stdout']:
            self.save_solution(sysinst, t, Xt, '-')
            append_file = sys.stdout
        else:
            append_file = open(file_name, 'a')

        while True:
            # Find final state to use as initial conditions here
            t1 = t[-1]
            x0 = Xt[-1, :]
            del t, Xt

            # Extend solution
            t = np.arange(t1, t1 + opts['T'] + opts['dtout'], opts['dtout'])
            Xt = sysinst.solve(x0, t, opts['dtmax'], method=opts['method'])

            # Remove repeated time
            t = t[1:]
            Xt = Xt[1:, :]

            # Save to same file or write to stdout
            sysinst.save_csv(t, Xt, append_file, header=False, titles=False)
            append_file.flush()

    # Create system, print parameters and exit
    params_opts = []
    def params(self, *args, **opts):
        """Display parameters and exit"""
        sysinst = self._make_sode(*args, **opts)
        if not sysinst:
            return 1
        print sysinst.get_description()

    info_opts = []
    def info(self, input_file=None):
        """Display parameters from solution file"""
        # Load system from file
        if input_file is None:
            input_file = sys.stdin
        sys_opts, sys_args = SODE.load_parameters(input_file)
        sysinst = self._make_sode(*sys_args, **sys_opts)
        if not sysinst:
            return 1
        print sysinst.get_description()

    def plot(self, input_file=None, **opts):
        """Plot previously saved solution file"""
        # Load solution object
        if input_file is None:
            input_file = sys.stdin
        t, Xt = SODE.load_csv(input_file, parameters=False)

        # figure is no-op without plotopts: force plot
        if not opts['plot_file']:
            opts['plot'] = True

        # And plot
        fig = self.figure(**opts)
        if not fig:
            raise Error("Unable to open plot window")

        # Actually plot
        ax = fig.add_subplot(1, 1, 1)
        self.plot_solution(ax, t, Xt, 'k')
        self.show(fig, **opts)

    # Options used by convergence routine
    conv_opts = sys_opts + plot_opts + [
        ('t', 't1', 0.0, 'time at start of solution'),
        ('T', 'T', 0.0, 'duration used for solution'),
        ('d', 'dtmin', 0.001, 'smallest dt'),
        ('n', 'levels', 2, 'Number of levels to coarse-grain by'),
        ('m', 'method', 'EM',  'integration method'),
        ('' , 'params', False,  'print parameters and exit'),
    ]
    def conv(self, *args, **opts):
        """Test convergence in dt

        Starting with dt = dtmin, compare results for 2*dt, 4*dt etc. for
        nlevels. If exact solution is available compares with exact solution,
        else compares all results are compared with dtmin.

        If T is not provided then 10 * (2 ** nlevels) is used so that each dt
        has at least 10 increments.

        Relative error in x_estimate is given as
        rerr = |(x_estimate - x_true) / (x_true + 1e-5)|
        """
        # Create system with args
        sysinst = self._make_sode(*args, **opts)
        if not sysinst:
            return 1
        if opts['params']:
            print sysinst.get_description()
            return 0

        # Create best solution (smallest dt)
        nsamples = abs(opts['T'] / opts['dtmin'])
        # nsamples must be a multiple of 2 ** nlevels for coarse-graining
        nmin = 2 ** opts['levels']
        if not nsamples:
            nsamples = 10 * nmin
        elif nsamples % nmin:
            nsamples = (nsamples // nmin + 1) * nmin

        # First do the best solution
        bp = sysinst.make_brownian_path(opts['t1'], opts['dtmin'], nsamples)
        x0 = sysinst.get_x0()

        # Create exact solution if possible
        try:
            t = bp.t.copy()
            t.resize(bp.Wt.shape)
            Xt_exact = sysinst.exact(x0, t, bp.Wt)
            best = {'t':t, 'Xt':Xt_exact, 'dt':bp.dt, 'label':'exact'}
        except NotImplementedError:
            best = None

        # Create corase-grained solutions and compare
        results = []
        for n in range(opts['levels']):
            r = 2 ** n
            lab = '{0}dt'.format(r)
            Xt = sysinst.solve_bm(x0, bp, method=opts['method'])
            results.append({'t':bp.t, 'Xt':Xt, 'dt':bp.dt, 'r':r, 'label':lab})
            bp = bp.coarse_grain()

        # Use smalles dt if necessary
        if best is None:
            best, results = results[0], results[1:]

        # Coarse grain comparison timeseries
        Xtbest = best['Xt']
        for n, d in enumerate(results):
            d['Xt_comp'] = Xtbest[::d['r'], :]

        # Compute summary
        for d in results:
            Xt = d['Xt']
            if 'Xt_comp' not in d:
                Xt_true = best['Xt']
            else:
                Xt_true = d['Xt_comp']
            # Copmute relative and absolute errors
            err = abs(Xt - Xt_true)
            ferr = err / (abs(Xt_true) + 1e-5)
            # For good matches
            d['exact'] = all(Xt.flat == Xt_true.flat)
            d['close'] = np.allclose(Xt, Xt_true)
            # Absolute errors
            d['err'] = err.sum(axis=1)
            d['mean'] = err.mean()
            d['max'] = err.max()
            # Relative errors
            d['ferr'] = ferr.sum(axis=1)
            d['fmean'] = ferr.mean()
            d['fmax'] = ferr.max()

        # Print output
        print "Coarse grain results from dt = {0}".format(opts['dtmin'])
        for d in reversed(results):
            print "{0} vs. {1} ({2}):".format(best['label'], d['label'], d['dt'])
            print "  Exactly equal   :", d['exact']
            print "  Within eps      :", d['close']
            print "  Mean error      :", d['mean']
            print "  Max. error      :", d['max']
            print "  Mean rel. error :", d['fmean']
            print "  Max. rel. error :", d['fmax']

        # Attempt to identify convergence order
        dts   = np.array([d['dt'  ] for d in results])
        if len(results) > 1:
            for vname, txt in (('max', 'max abs error'),
                               ('fmax', 'max rel error'),
                               ('mean', 'mean abs error'),
                               ('fmean', 'mean rel error')):
                vals  = np.array([d[vname] for d in results])
                inc = (vals != 0) & (dts != 0)
                if sum(inc) >= 2:
                    x = np.log(dts[inc])
                    y = np.log(vals[inc])
                    m, C = np.polyfit(x, y, 1)
                    coeff, order = np.exp(C), m
                    print "{0} behaves as: {1:3.2}*dt^{2:3.2}".format(txt, coeff, order)
                    print "Convergence order: {0:3.2}".format(order)

        # Plot results
        fig = self.figure(**opts)
        if fig:
            # Three plots. Time series for comparison, absolute and relative
            # errors relative to smallest dt.
            ax1 = fig.add_subplot(3, 1, 1)
            ax2 = fig.add_subplot(3, 1, 2, sharex=ax1)
            ax3 = fig.add_subplot(3, 1, 3, sharex=ax2)
            ax1.set_title('Convergence from dt = {0}'.format(opts['dtmin']))
            ax2.set_title('Absolute error from dt = {0}'.format(opts['dtmin']))
            ax3.set_title('Relative error from dt = {0}'.format(opts['dtmin']))
            self.plot_solution(ax1, best['t'], best['Xt'], 'k', label=best['label'])
            for d in results:
                lab = d['label']
                self.plot_solution(ax1, d['t'], d['Xt'], '.', label=lab, **opts)
                ax2.plot(d['t'], d['err'], label=lab)
                ax3.plot(d['t'], d['ferr'], label=lab)
            ax1.legend()
            ax2.legend()
            ax3.legend()
            self.show(fig, **opts)

    #
    # Utility functions used in the command above are defined below
    #

    # Subclasses can override sys_opts to acces cl options in make_sode
    def make_sode(self, syskwargs, **sysopts):
        """Called using the command line arguments to the script.

        Subclasses should override this to choose the system from the cl
        args.

        args : the positional command line arguments.
        opts : the option values provided on the command line.
        """
        return self.SYSTYPE(**syskwargs)

    def _make_sode(self, *args, **opts):
        """Wraps make_sode to parse make_sode_opts from opts"""
        # Parse relevant options from opts
        sysopts = {}
        for _, name, _, _ in self.sys_opts:
            name = name.replace('-', '_')
            if name in opts:
                sysopts[name] = opts[name]

        # Parse args into key val pairs
        syskwargs = dict([a.split('=') for a in args])

        # Create system with opts and parameters and variables
        sysinst = self.make_sode(syskwargs, **sysopts)

        # sysinst.save_csv stores this:
        sysinst._sys_opts = sysopts, args
        return sysinst

    def plot_solution(self, ax, t, Xt, linestr, label=None, **opts):
        """Plot solution in style specified by opts"""
        ax.plot(t, Xt, linestr, label=label)

    def save_solution(self, sysinst, t, Xt, output_file, **opts):
        """Save solution t, Xt to output_file ('-' for stdout)"""
        if not output_file:
            return
        if output_file == '-':
            fout = sys.stdout
        else:
            fout = open(output_file, 'w')
        sysinst.save_csv(t, Xt, fout)

    def figure(self, plot=False, plot_file='', **opts):
        if not (plot or plot_file):
            return None
        from matplotlib.pyplot import figure
        return figure()

    def show(self, fig, plot=False, plot_file='', **opts):
        """Plot to file and/or screen"""
        if plot_file:
            ax.savefig(plot_file)
        if plot:
            from matplotlib.pyplot import show
            show()

class MultiScript(Script):
    """Subclass Script to handle multiple SODE classes"""
    def __init__(self, SYSDICT):
        """Adds system choosing options to Script"""
        self.SYSDICT = SYSDICT
        self.sys_opts = [('s', 'system-type', '', 'Name of SODE system')]
        Script.__init__(self)
        self.cmdtable['list'] = (self.list_, [], 'List possible system names')

    def make_sode(self, syskwargs, system_type):
        """Choose system type from --system-type"""
        if not system_type:
            system_type = self.SYSDICT[None]
        if system_type in self.SYSDICT:
            syscls = self.SYSDICT[system_type]
            return syscls(**syskwargs)
        else:
            msg ="Unrecognised system name '{0}' use list".format(system_type)
            print >> sys.stderr, msg
            return None

    def list_(self, *args, **opts):
        print ', '.join([k for k in self.SYSDICT if k is not None])




# Quick test suite to demonstrate usage of sode.py and script.py
if __name__ == "__main__":

    # Define an example SODE subclass

    class LinearSODE(SODE):
        """1-D linear SODE with equation:

        dx = - beta x dt + alpha dW
        """
        variables = (('x', 1.0),)
        parameters = (('beta', 1.0), ('alpha', 1.0))

        # drift coefficient (deterministic derivative)
        def drift(self, a, x, t):
            a[self.x] = - self.beta * x[self.x]
            return a

        # diffusion coefficient (noise amplitude)
        def diffusion(self, b, x, t):
            b[self.x] = self.alpha
            return b

    import sys
    script = Script(LinearSODE)
    script.main(argv=sys.argv[1:])

