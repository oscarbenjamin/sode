# Copyright Oscar Benjamin 2011
#
# Defines the CYSODE Extension type

import cython

from libc cimport math

import numpy as np
cimport numpy as np


cdef extern from "numpy/arrayobject.h":
    int _import_array()
    int _import_umath()
_import_array()
_import_umath()

cdef extern from "cfiles/randnorm.h":
    double RANDNORM_NORMAL()
    unsigned long randnorm_jsr, randnorm_jz
    long randnorm_hz
    unsigned long randnorm_iz, randnorm_kn[128]
    double randnorm_wn[128]
    double randnorm_nfix()
    void randnorm_seed(unsigned int)
    enum:
        RANDNORM_SEED_PID_TIME
randnorm_seed(RANDNORM_SEED_PID_TIME)

DTYPE = np.float64

cdef class CYSODE:

    def __cinit__(self):
        self.nvars = 1

    cpdef drift(self, np.ndarray[DTYPE_t, ndim=1] a,
                      np.ndarray[DTYPE_t, ndim=1] x, double t):
        self._drift(<double*>a.data, <double*>x.data, t)

    cpdef diffusion(self, np.ndarray[DTYPE_t, ndim=1] b,
                          np.ndarray[DTYPE_t, ndim=1] x, double t):
        self._diffusion(<double*>b.data, <double*>x.data, t)

    cdef _drift(self, double* a, double* x, double t):
        raise NotImplementedError("Subclasses need to override this.")

    cdef _diffusion(self, double* b, double* x, double t):
        raise NotImplementedError("Subclasses need to override this.")

    cdef _solve_EM(self, double* x1, double t1, double* x2, double dt,
                         double sqrtdt, double *a, double *b):
        cdef unsigned int i
        self._drift(a, x1, t1)
        self._diffusion(b, x1, t1)
        for i in range(self.nvars):
            x2[i] = x1[i] + a[i] * dt + b[i] * sqrtdt * RANDNORM_NORMAL()

    cpdef get_x0(self):
        return np.array([0], np.float64)

    def get_description(self):
        return "CYSODE!!!!!"

    def get_variables(self):
        return ['x']

    cpdef solveto(self, np.ndarray[DTYPE_t, ndim=1] x1, double t1,
                        np.ndarray[DTYPE_t, ndim=1] x2, double t2, double dtmax):
        cdef double t = t1, tnext = t, dt = dtmax
        cdef double sqrtdt = math.sqrt(dt)
        cdef double temp
        cdef np.ndarray[DTYPE_t, ndim=1] a = np.zeros(self.nvars, DTYPE)
        cdef np.ndarray[DTYPE_t, ndim=1] b = np.zeros(self.nvars, DTYPE)
        cdef np.ndarray[DTYPE_t, ndim=1] x = np.zeros(self.nvars, DTYPE)
        for i in range(self.nvars):
            x[i] = x1[i]
        while t < t2:
            tnext = t + dt
            if tnext > t2:
                tnext = t2
                dt = t2 - t
                sqrtdt = math.sqrt(dt)
            self._solve_EM(<double*>x.data, t, <double*>x.data, dt,
                           sqrtdt, <double*>a.data, <double*>b.data)
            t = tnext
        for i in range(self.nvars):
            x2[i] = x[i]
