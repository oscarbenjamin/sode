#!/usr/bin/env python

from __future__ import division

import sys
import os.path

import numpy as np

from sode import SODE, SODENetwork, Script
from graphviz import GraphParser


class Region(SODENetwork):
    """2-D SODE with subcritical hopf.

    Complex equation:

    dz = f(z) dt + alpha dW(t)

    where
    f(z) = ((lambda - 1) + i omega) z + 2z|z|^2 - z|z|^4
    """

    # Setup code

    def __init__(self, brain, **kwargs):
        super(Region, self).__init__()
        self._inputs = []
        self._brain = brain
        self.init(**kwargs)

    def add_input_from(self, name, rfrom, **kwargs):
        c = Connection(rfrom, self, **kwargs)
        self._inputs.append(c)
        self.add_subsystem(name, c)
    
    # Define the dynamics

    variables = (('x', 1), ('y', 0))
    parameters = (('omega', 20), ('alpha', 0.01))

    def diffusion(self, b, X, t):
        b[self.x] = b[self.y] = self.alpha
        return b

    def get_lambda(self, X):
        return X[self._brain.lam]

    def drift(self, a, X, t):
        # Uncomment to check consistency:
        #b, c = np.zeros_like(X), np.zeros_like(X)
        #assert np.allclose(self.drift_slow(b, X, t), self.drift_fast(c, X, t))
        return self.drift_fast(a, X, t)

    def drift_slow(self, a, X, t):
        x, y = X[self.x], X[self.y]
        lp, om = self.get_lambda(X) - 1, self.omega
        modz2 = x ** 2 + y ** 2
        a[self.x] = lp*x - om*y + x*(2*modz2 - modz2**2)
        a[self.y] = lp*y + om*x + y*(2*modz2 - modz2**2)
        a[self.x] += self.inputs_x_slow(X)
        a[self.y] += self.inputs_y_slow(X)
        return a

    def inputs_x_slow(self, X):
        return sum([c.input_x(X) for c in self._inputs])

    def inputs_y_slow(self, X):
        return sum([c.input_y(X) for c in self._inputs])

    def drift_fast(self, a, X, t):
        x = X[self.x]
        y = X[self.y]
        lp = self.get_lambda(X)
        om = self.omega
        modz2 = x ** 2 + y ** 2
        factor = lp - 1 + 2 * modz2 - modz2 ** 2
        a[self.x] = x * factor - om * y + + self.inputs_x_fast(X)
        a[self.y] = y * factor + om * x + + self.inputs_y_fast(X)
        return a

    def inputs_x_fast(self, X):
        if not hasattr(self, '_xcache'):
            self._xrfroms = [c.rfrom.x for c in self._inputs]
            self._xrtos = [c.rto.x for c in self._inputs]
            self._xbetas = np.array([c.beta for c in self._inputs])
            self._xcache = True
        return np.dot(self._xbetas, (X[self._xrfroms] - X[self._xrtos]))

    def inputs_y_fast(self, X):
        if not hasattr(self, '_ycache'):
            self._yrfroms = [c.rfrom.y for c in self._inputs]
            self._yrtos = [c.rto.y for c in self._inputs]
            self._ybetas = np.array([c.beta for c in self._inputs])
            self._ycache = True
        return np.dot(self._ybetas, X[self._yrfroms] - X[self._yrtos])

class PathologicalRegion(Region):
    """Like Region but uses lambda -> lambda + delta_lambda
    
    Has a parameter, delta_lambda which is added to lambda before the
    calculations as described in Region"""
    parameters = Region.parameters + (('delta_lam', 0.05),)

    # Overload to increment by delta_lam
    def get_lambda(self, X):
        return Region.get_lambda(self, X) + self.delta_lam


class Connection(SODE):
    """Connection between regions
    
    This is diffusive coupling so that:
    
    dz1dt = f(z1) - beta (z1 - z2)
    """
    
    parameters = (('beta', 1),)
    variables = ()
    
    # Store references, and compute input

    def __init__(self, rfrom, rto, **kwargs):
        self.rfrom = rfrom
        self.rto = rto
        super(Connection, self).__init__(**kwargs)

    def input_x(self, X):
        return self.beta * (X[self.rfrom.x] - X[self.rto.x])
    
    def input_y(self, X):
        return self.beta * (X[self.rfrom.y] - X[self.rto.y])


class Brain(SODENetwork):
    """System describing the Brain composed of a network of Regions
    
    Each region has an equation:
    
    dzi/dt = f(lambda, zi) + g(zi, Z)

    where 

    Lambda rather than a parameter is a global slow variable with equation


    """

    # Setup code

    def __init__(self):
        """Adds 2 nodes A and B"""
        super(Brain, self).__init__()
        self.regions = []

    def add_region(self, name, r):
        self.add_subsystem(name, r)
        self.regions.append(r)

    # Dynamics of the brain

    variables = (('lam', 0.5),)
    parameters = (('lamr', -1), ('taur', 10), ('lame', 1.0), ('taue', 100))

    def g(self, X):
        z2 = 0
        for r in self.regions:
            z2 += X[r.x] ** 2 + X[r.y] ** 2
        return z2 / len(self.regions)

    def drift(self, a, X, t):
        gz = self.g(X)
        a[self.lam] = - (gz / self.taur) * (X[self.lam] - self.lamr)\
                      - ((1 - gz) / self.taue) * (X[self.lam] - self.lame)
        return self.drift_subsys(a, X, t)

    def diffusion(self, b, X, t):
        b[self.lam] = 0
        return self.diffusion_subsys(b, X, t)


def Network(graph_file, syskwargs):
    """Adds the appropriate regions and connections to the brain based on the
    data in the graph_parser file and any parameters provided on the command
    line."""
    # Parse graph file
    graph_parser = GraphParser(graph_file)

    # Global region and connection parameters from command line
    cl_reg_params = {}
    cl_con_params = {}
    if 'alpha' in syskwargs:
        cl_reg_params['alpha'] = syskwargs['alpha']
    if 'beta' in syskwargs:
        cl_con_params['beta'] = syskwargs['beta']
    
    brain = Brain()

    # Create regions (a sorteddict would be usefult here)
    # Regions are initialised using node parameters from graphviz file
    # These are overridden using any global command line parameters
    regs = []
    for name, pars in graph_parser.iter_nodes():
        pars.update(cl_reg_params)
        if 'delta_lam' in pars:
            r = PathologicalRegion(brain, **pars)
        else:
            r = Region(brain, **pars)
        regs.append((name, r))
    regs_d = dict(regs)

    # Connect regions use (region name as connection name)
    # Default parameters from graphviz file. These are overridden by any
    # command line parameters.
    for rfrom, rto, pars in graph_parser.iter_edges():
        # Ignore edges with beta=0 in graph file
        if 'beta' in pars and float(pars['beta']) == 0:
            continue
        pars.update(cl_con_params)
        regs_d[rto].add_input_from(rfrom, regs_d[rfrom], **pars)

    # Add regions and connections to the brain
    for name, region in regs:
        brain.add_region(name, region)

    # Apply graphfile and command lines parameters
    brain_pars = graph_parser.graph_params
    brain_pars.update(syskwargs)
    brain.init(**brain_pars)

    #Return the SODE instance
    return brain


TEST_GRAPH = """
digraph {
    A
}
"""


class GraphScript(Script):
    """Subclass Script to handle multiple SODE classes"""

    def __init__(self, *args, **kwargs):
        """Adds graph file choosing option to Script"""
        self.sys_opts = [
            ('g', 'graph-file', '', 'file to load graph from'),
            ('e', 'eigenvalues', False, 'display eigenavlues of L'),
        ]
        Script.__init__(self)


    def make_sode(self, syskwargs, graph_file, eigenvalues=False):
        """Single node if no graph supplied"""

        # Use test graph
        if not graph_file:
            from StringIO import StringIO
            graph_file = StringIO(TEST_GRAPH)
        
        # Check that graph file exists
        elif not os.path.exists(graph_file):
            msg = "file '{0}' does not exist".format(graph_file)
            print >> sys.stderr, msg
            return None

        if eigenvalues:
            g = GraphParser(graph_file)
            names = g.node_order
            N = len(names)
            M = np.zeros((N, N))

        # And return
        return Network(graph_file, syskwargs)


# Create main function
_script = GraphScript(Brain)
main = _script.main


if __name__ == "__main__":
    main(argv=sys.argv[1:])

