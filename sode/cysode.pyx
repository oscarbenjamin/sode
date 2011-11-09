# Copyright Oscar Benjamin 2011
#
# Defines the CYSODE Extension type

from libc cimport math

import numpy as np
cimport numpy as np

cdef extern from "numpy/arrayobject.h":
    int _import_array()
    int _import_umath()
_import_array()
_import_umath()

cdef extern from "cfiles/randnorm.h":
    enum:
        RANDNORM_NORMAL
    unsigned long randnorm_jsr, randnorm_jz
    long randnorm_hz
    unsigned long randnorm_iz, randnorm_kn[128]
    double randnorm_wn[128]
    double randnorm_nfix()

cdef class CYSODE:

    nvars = 1
    cdef public _sys_opts

    cpdef drift(self, np.ndarray a, np.ndarray x, double t):
        if a is None or x is None:
            raise ValueError("Invalid args")
        a[0] = 0
        return a

    cpdef diffusion(self, np.ndarray b, np.ndarray x, double t):
        if b is None or x is None:
            raise ValueError("Invalid args")
        b[0] = 1
        return b

    cpdef get_x0(self):
        return np.array([0], np.float64)

    def get_description(self):
        return "CYSODE!!!!!"

    def get_variables(self):
        return ['x']

    cpdef solveto(self, np.ndarray x1, double t1, np.ndarray x2, double t2, double dtmax):
        cdef double t = t1
        cdef double dt
        cdef double tnext
        cdef double sqrtdt = math.sqrt(dtmax)
        cdef np.ndarray a = np.zeros(self.nvars)
        cdef np.ndarray b = np.zeros(self.nvars)
        cdef np.ndarray x = np.zeros(self.nvars)
        while t <= t2:
            tnext = t + dtmax
            if tnext > t2:
                tnext = t2
                dt = t2 - t
                sqrtdt = math.sqrt(dtmax)
            tnext = math.min()
            self.drift(a, x, t)
            self.diffusion(b, x, t)
            for i in range(self.nvars):
                x[i] += a[i] * dt + b[i] * sqrtdt * RANDNORM_NORMAL
            t = tnext
        for i in range(self.nvars):
            x2[i] += x[i]
