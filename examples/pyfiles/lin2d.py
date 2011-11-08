#!/usr/bin/env python
#
# Author: Oscar Benjamin
# Date: 11 Mar 2011
#
# Example of a script to investigate a 2-D SODE. In the absence of noise this
# system is a damped oscillator. However, noise causes it to oscillate
# irregularly rather than simply converge on the fixed point.


from sode.pysode  import SODE


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


if __name__ == "__main__":
    import sys
    from sode.script import Script
    script = Script(LinearAdditive2D)
    script.main(argv=sys.argv[1:])
