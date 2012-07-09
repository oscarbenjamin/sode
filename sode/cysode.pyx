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

DTYPE = np.float64

ctypedef np.float64_t Real
ctypedef unsigned int Int

cdef class Vector:

    cdef Real* data
    cdef Int length

    def __cinit__(self, Int length):
        self.length = length
        self.data = <Real*>malloc(length * sizeof(Real))
        if self.data is NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self.data is not NULL:
            free(self.data)
            self.data = NULL

    def __setitem__(self, Int index, Real val):
        self.data[index] = val

    def __getitem__(self, Int index):
        return self.data[index]

cdef class SODE:

    cdef readonly Int _nvars
    cdef public dict _sys_opts

    def __cinit__(self):
        self._nvars = -1

    cpdef _cy_drift(self, Vector a, Vector x, Real t):
        raise NotImplementedError('Subclasses should define this')

    cpdef _cy_diffusion(self, Vector b, Vector x, Real t):
        raise NotImplementedError('Subclasses should define this')

    def exact(self, x, t, Wt):
        raise NotImplementedError('Subclasses should define this')

    # This is a temporary stopgap
    cpdef solveto(self, np.ndarray[Real, ndim=1] x1, Real t1,
                        np.ndarray[Real, ndim=1] x2, Real t2, Real dtmax):
        cdef t, tnext, dt, sqrtdt
        cdef Vector a, b, x

        # Allocate temporary memory
        a = Vector(self._nvars)
        b = Vector(self._nvars)
        x = Vector(self._nvars)

        # Prepare for loop
        t = tnext = t1
        dt = dtmax
        sqrtdt = math.sqrt(dt)

        # Copy in the initial state
        for i in range(self.nvars):
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
        for i in range(self.nvars):
            x2[i] = x[i]

    cdef void _solve_EM(self, Vector x1, Real t1, Vector x2,
                              Real dt, Real sqrtdt,
                              Vector a, Vector b):
        cdef Int i
        self._cy_drift(a, x1, t1)
        self._cy_diffusion(b, x1, t1)
        for i in range(self.nvars):
            x2[i] = x1[i] + a[i] * dt + b[i] * sqrtdt * RANDNORM_NORMAL()

cdef class SODE_test(SODE):

    def __cinit__(self):
        self._nvars = 2

    cpdef _cy_drift(self, Vector a, Vector x, Real t):
        a[0] = x[1]
        a[1] = - x[0]

    cpdef _cy_diffusion(self, Vector b, Vector x, Real t):
        b[0] = b[1] = 0.01

cdef class CYSODE:

    def __cinit__(self):
        self.nvars = 0
        self._owned = 0

    def __init__(self, **kwargs):
        self._init_pv(self.parameters, self.variables)
        self._init_kw(**kwargs)

    cpdef drift(self, np.ndarray[DTYPE_t, ndim=1] a,
                      np.ndarray[DTYPE_t, ndim=1] x, double t):
        self._drift(<double*>a.data, <double*>x.data, t)
        return a

    cpdef diffusion(self, np.ndarray[DTYPE_t, ndim=1] b,
                          np.ndarray[DTYPE_t, ndim=1] x, double t):
        self._diffusion(<double*>b.data, <double*>x.data, t)
        return b

    cdef void _drift(self, double* a, double* x, double t):
        raise NotImplementedError("Subclasses need to override this.")

    cdef void _diffusion(self, double* b, double* x, double t):
        raise NotImplementedError("Subclasses need to override this.")

    def exact(self, x, t, Wt):
        raise NotImplementedError("Subclasses need to override this.")

    cdef void _solve_EM(self, double* x1, double t1, double* x2, double dt,
                         double sqrtdt, double *a, double *b):
        cdef unsigned int i
        self._drift(a, x1, t1)
        self._diffusion(b, x1, t1)
        for i in range(self.nvars):
            x2[i] = x1[i] + a[i] * dt + b[i] * sqrtdt * RANDNORM_NORMAL()

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

    # These functions are recycled from pysode
    _init_pv        = pysode._init_pv_
    _init_kw        = pysode._init_kw_
    _shift_indices  = pysode._shift_indices_
    _get_var_index  = pysode._get_var_index_
    _set_var_index  = pysode._set_var_index_
    _varchk         = pysode._varchk_
    _parchk         = pysode._parchk_
    get_variables   = pysode._get_variables
    get_ic          = pysode._get_ic
    set_ic          = pysode._set_ic
    get_x0          = pysode._get_x0
    get_parameters  = pysode._get_parameters
    get_parameter   = pysode._get_parameter
    set_parameter   = pysode._set_parameter
    get_description = pysode._get_description


cdef class CYSODE_test(SODE):

    parameters = ()
    variables = (('x', 1), ('y', 0))

    cdef public parameter x
    cdef public parameter y

    cpdef _cy_drift(self, Vector a, Vector x, Real t):
        a[0] = x[1]
        a[1] = - x[0]

    cpdef _cy_diffusion(self, Vector b, Vector x, Real t):
        b[0] = b[1] = 0.01

cdef class CYSODENetwork(CYSODE):

    def __cinit__(self):
        self._subsystems = []
        self._subsysdict = {}
        self._ss_eval = []

    cdef void _drift_subsys(self, double* a, double* x, double t):
        """Compute drift for subsystems. Subclasses must call this in drift"""
        cdef unsigned int N = len(self._ss_eval)
        for i in range(N):
            (<CYSODE>self._ss_eval[i])._drift(a, x, t)

    cdef void _diffusion_subsys(self, double* b, double* x, double t):
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
