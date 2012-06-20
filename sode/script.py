# Copyright, Oscar Benjamin 2011
#
# script.py
#
# This module exports the Script class. The Script class uses opster.dispatch
# to make a flexible script for generating solutions of SODE classes.

import sys

import numpy as np

from sode.opster import dispatch
from sode.pysode import SODE
from sode.algos import solve, solve_bm, BrownianMotion
from sode.io import load_csv, save_csv

class Script(object):
    """A Script instance is a callable that can be used to quickly create a
    script to investigate a particular SODE subclass.

    There are two ways to use Script. One is to create a Script instance with
    an SODE subclass as argument. The other is to override the make_sode
    method. This method is called with all command line positional and option
    arguments as *args and **kwargs respectively. These two methods are
    illustrated below.

    Example script:
    -----------------------------------
    Define a simple system
    class MySODE(SODE):
      def f(self, x, t):
          return x
      def g(self, x, t):
          return t

    Use script main if run as a script
    if __name__ == "__main__":
        import sys
        from script import Script
        script = Script(MySODE)
        script(argv=sys.argv[1:])
    # -----------------------------------

    A more complicated example:
    # -----------------------------------
    # Define several systems
    class MySODE1(SODE):
        ...
    class MySODE2(SODE):
        ...

    # Collect together in a dict
    SYSTEMS = {'sys1':MySODE1, 'sys2':MySODE2, ...}

    if __name__ == "__main__":
        import sys
        from script import Script

        # Override Script to choose system type
        class MyScript(Script):
           usage = '%name [OPTS] SYSNAME [PAR1=VAL1 ...]
           def make_sode(self, sysname, *args, **opts):
               if sysname not in SYSTEMS:
                   msg = "System name should be one of {0}"
                   raise ValueError(msg.format(', '.join(SYSTEMS)))
               return SYSTEMS[sysname](*args)

        # Actually run as a script
        script = MyScript()
        script(argv=sys.argv[1:])
    """
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
        Xt = solve(sysinst, x0, t, opts['dtmax'], method=opts['method'])

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
        t, Xt, (sysopts, sysargs) = load_csv(file_name, parameters=True)
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
            Xt = solve(sysinst, x0, t, opts['dtmax'], method=opts['method'])

            # Remove repeated time
            t = t[1:]
            Xt = Xt[1:, :]

            # Save to same file or write to stdout
            save_csv(sysinst, t, Xt, append_file, header=False, titles=False)
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
        t, Xt = load_csv(input_file, parameters=False)

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
        bp = BrownianMotion(opts['t1'], opts['dtmin'], (nsamples, sysinst.nvars))
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
            Xt = solve_bm(sysinst, x0, bp, method=opts['method'])
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
        save_csv(sysinst, t, Xt, fout)

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
        sysnames = [k for k in self.SYSDICT if k is not None]
        helptext = 'Name of SODE system: {0}'.format(', '.join(sysnames))
        self.sys_opts = [('s', 'system-type', '', helptext)]
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


if __name__ == "__main__":
    import doctest
    doctest.testmod()
