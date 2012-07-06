# Copyright Oscar Benjamin 2011
#
# This script defines a trivial SODE representing just a Weiner process. This
# script is called when sode is run at the cmd line.


import sys

from sode.pysode import SODE
from sode.script import Script


class Weiner(SODE):
    """Weiner process with drift coeff. mu and diffusion coeff. sigma

        dx(t) = mu dt + sigma dW(t)

    This trivial SODE has the exact solution:

        x(t) = x(0) + mu t + sigma W(t)
    """
    variables = (('x', 0),)
    parameters = (('mu', 0), ('sigma', 1))

    def drift(self, a, x, t):
        return self.mu

    def diffusion(self, b, x, t):
        return self.sigma

    def exact(self, x0, t, Wt):
        return (x0 + self.mu * t) + self.sigma * Wt


# Use sode.script to create a command line interface for this system
main = Script(Weiner).main


if __name__ == "__main__":
    main(sys.argv[1:])
