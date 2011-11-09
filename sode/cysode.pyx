# Copyright Oscar Benjamin 2011
#
# Defines the CYSODE Extension type

import numpy as np
cimport numpy as np

cdef extern from "numpy/arrayobject.h":
    int _import_array()
    int _import_umath()
_import_array()
_import_umath()

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
