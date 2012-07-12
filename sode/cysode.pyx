# Copyright Oscar Benjamin 2011
#
# Defines the CYSODE Extension type

import cython

from libc cimport math
from libc.stdlib cimport malloc, free

import numpy as np
cimport numpy as np

from sode import pysode

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


cdef class Vector:

    def __len__(self):
        return self.length

    def __setitem__(self, Int index, Real val):
        self.data[index] = val

    def __getitem__(self, Int index):
        return self.data[index]


cdef class VectorAlloc(Vector):

    def __cinit__(self, Int length):
        self.length = length
        self.data = <Real*>malloc(length * sizeof(Real))
        if self.data is NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self.data is not NULL:
            free(self.data)
            self.data = NULL


cdef class VectorView(Vector):

    def __cinit__(self, np.ndarray npv):
        self.length = <Int>npv.size
        self.data = <Real*>npv.data


cdef class SODE:

    def __cinit__(self):
        self._nvars = -1
        self._x0 = []
        self._variables = []
        self._init = False

    cpdef _set_variables(self, names, values):
        assert len(names) == len(values), 'len(names) != len(values)'
        self._nvars = <Int>len(names)
        self._variables = list(names)
        self._x0 = list([float(v) for v in values])
        self._init = True

    cpdef _get_variables(self):
        return list(self._variables)

    def get_x0(self):
        return list(self._x0)

    def get_description(self):
        return 'No description for SODE instances'

    def get_variables(self):
        return list(self._variables)

    cpdef _drift(self, Vector a, Vector x, Real t):
        raise NotImplementedError('Subclasses should define this')

    cpdef _diffusion(self, Vector b, Vector x, Real t):
        raise NotImplementedError('Subclasses should define this')

    def _eval_drift(self, x, t):
        cdef np.ndarray[Real, ndim=1] xnp = np.array(x, dtype=DTYPE)
        cdef np.ndarray[Real, ndim=1] anp = np.zeros_like(x)
        cdef xv = VectorView(xnp)
        cdef av = VectorView(anp)
        self._drift(av, xv, t)
        return anp

    def _eval_diffusion(self, x, t):
        cdef np.ndarray[Real, ndim=1] xnp = np.array(x, dtype=DTYPE)
        cdef np.ndarray[Real, ndim=1] bnp = np.zeros_like(x)
        cdef xv = VectorView(xnp)
        cdef bv = VectorView(bnp)
        self._diffusion(bv, xv, t)
        return bnp

    def exact(self, x, t, Wt):
        raise NotImplementedError('Subclasses should define this')

    # This is a temporary stopgap
    cpdef solveto(self, np.ndarray[Real, ndim=1] x1, Real t1,
                        np.ndarray[Real, ndim=1] x2, Real t2, Real dtmax):
        cdef t, tnext, dt, sqrtdt
        cdef Vector a, b, x

        if not self._init:
            raise RuntimeError('Not initialised yet')

        # Allocate temporary memory
        a = VectorAlloc(self._nvars)
        b = VectorAlloc(self._nvars)
        x = VectorAlloc(self._nvars)

        # Prepare for loop
        t = tnext = t1
        dt = dtmax
        sqrtdt = math.sqrt(dt)

        # Copy in the initial state
        for i in range(self._nvars):
            x[i] = x1[i]

        # Main loop
        while t < t2:
            tnext = t + dt
            if tnext > t2:
                tnext = t2
                dt = t2 - t
                sqrtdt = math.sqrt(dt)
            self._solve_EM(x, t, x, dt, sqrtdt, a, b)
            t = tnext

        # Copy out the final state
        for i in range(self._nvars):
            x2[i] = x[i]

    cdef void _solve_EM(self, Vector x1, Real t1, Vector x2,
                              Real dt, Real sqrtdt,
                              Vector a, Vector b):
        cdef Int i
        self._drift(a, x1, t1)
        self._diffusion(b, x1, t1)
        for i in range(self._nvars):
            x2[i] = x1[i] + a[i] * dt + b[i] * sqrtdt * RANDNORM_NORMAL()


cdef class CYSODE(SODE):
    pass


cdef class CYSODENetwork(CYSODE):

    def __cinit__(self):
        self._subsystems = []
        self._subsysdict = {}
        self._ss_eval = []

    cdef void _drift_subsys(self, Vector a, Vector x, Real t):
        """Compute drift for subsystems. Subclasses must call this in drift"""
        cdef unsigned int N = len(self._ss_eval)
        for i in range(N):
            (<CYSODE>self._ss_eval[i])._drift(a, x, t)

    cdef void _diffusion_subsys(self, Vector b, Vector x, Real t):
        """Compute diffusion for subsystems. Subclasses must call this in
        diffusion"""
        cdef unsigned int N = len(self._ss_eval)
        for i in range(N):
            (<CYSODE>self._ss_eval[i])._diffusion(b, x, t)

    add_subsystem  = pysode._n_add_subsystem
    _shift_indices = pysode._n__shift_indices
    get_parameters = pysode._n_get_parameters
    get_variables  = pysode._n_get_variables
    _get_subsystem = pysode._n__get_subsystem
    get_parameter  = pysode._n_get_parameter
    set_parameter  = pysode._n_set_parameter
    get_ic         = pysode._n_get_ic
    set_ic         = pysode._n_set_ic
    _init          = pysode._n__init
