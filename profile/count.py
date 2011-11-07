#!/usr/bin/env python
import sys, math
nums = [float(l) for l in sys.stdin]
print min(nums), max(nums)

from pylab import *
fig = figure()
ax = fig.add_subplot(1, 1, 1)
xlim = (-5, 5)
ax.hist(nums, bins=100, range=xlim, normed=True)
ax.set_xlim(xlim)
show()


