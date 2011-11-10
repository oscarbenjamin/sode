
import numpy as np
cimport numpy as np

ctypedef np.float64_t DTYPE_t
ctypedef double parameter
ctypedef unsigned int variable

cdef class CYSODE:
    cdef public unsigned int nvars
    cdef public _sys_opts, _variables, _parameters, _x0
    cdef public bint _owned
    cpdef drift(self, np.ndarray[DTYPE_t, ndim=1] a,
                      np.ndarray[DTYPE_t, ndim=1] x, double t)
    cpdef diffusion(self, np.ndarray[DTYPE_t, ndim=1] b,
                          np.ndarray[DTYPE_t, ndim=1] x, double t)
    cdef void _drift(self, double* a, double* x, double t)
    cdef void _diffusion(self, double* b, double* x, double t)
    cdef void _solve_EM(self, double* x1, double t1, double* x2, double dt,
                              double sqrtdt, double *a, double *b)
    cpdef solveto(self, np.ndarray[DTYPE_t, ndim=1] x1, double t1,
                        np.ndarray[DTYPE_t, ndim=1] x2, double t2, double dtmax)

cdef class CYSODENetwork(CYSODE):

    cdef public list _subsystems
    cdef public dict _subsysdict
    cdef public list _ss_eval

    cdef void _drift_subsys(self, double* a, double* x, double t)
    cdef void _diffusion_subsys(self, double* b, double* x, double t)
