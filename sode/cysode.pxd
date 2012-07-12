
import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ctypedef np.float64_t Real
ctypedef unsigned int Int

ctypedef Real parameter
ctypedef unsigned int variable


cdef class Vector:
    cdef Int length
    cdef Real* data

cdef class SODE:

    cdef readonly Int _nvars
    cdef list _variables
    cdef _x0
    cdef readonly int _init
    cdef public _sys_opts

    cpdef _get_variables(self)
    cpdef _set_variables(self, names, values)

    cpdef _drift(self, Vector a, Vector x, Real t)
    cpdef _diffusion(self, Vector b, Vector x, Real t)

    cdef void _solve_EM(self, Vector x1, Real t1, Vector x2,
                              Real dt, Real sqrtdt, Vector a, Vector b)
    cpdef solveto(self, np.ndarray[Real, ndim=1] x1, Real t1,
                        np.ndarray[Real, ndim=1] x2, Real t2, Real dtmax)

cdef class CYSODE(SODE):
    cdef public _parameters
    cdef public bint _owned

cdef class CYSODENetwork(CYSODE):

    cdef public list _subsystems
    cdef public dict _subsysdict
    cdef public list _ss_eval

    cdef void _drift_subsys(self, Vector a, Vector x, Real t)
    cdef void _diffusion_subsys(self, Vector b, Vector x, Real t)
