#!python
#cython: language_level=3
# This file is generated automatically by QuTiP.

import numpy as np
cimport numpy as np
import scipy.special as spe
cimport cython
np.import_array()
cdef extern from "numpy/arrayobject.h" nogil:
    void PyDataMem_NEW_ZEROED(size_t size, size_t elsize)
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)
from qutip.cy.spmatfuncs cimport spmvpy
from qutip.cy.inter cimport _spline_complex_t_second, _spline_complex_cte_second
from qutip.cy.inter cimport _spline_float_t_second, _spline_float_cte_second
from qutip.cy.inter cimport _step_float_cte, _step_complex_cte
from qutip.cy.inter cimport _step_float_t, _step_complex_t
from qutip.cy.interpolate cimport (interp, zinterp)
from qutip.cy.cqobjevo_factor cimport StrCoeff
from qutip.cy.cqobjevo cimport CQobjEvo
from qutip.cy.math cimport erf, zerf
from qutip.qobj import Qobj
cdef double pi = 3.14159265358979323

include 'C:/ProgramData/Anaconda3/lib/site-packages/qutip/cy/complex_math.pxi'

cdef class CompiledStrCoeff(StrCoeff):
    cdef double omega
    cdef double A

    def set_args(self, args):
        self.omega=args['omega']
        self.A=args['A']

    cdef void _call_core(self, double t, complex * coeff):
        cdef double omega = self.omega
        cdef double A = self.A

        coeff[0] = 0.5 * (np.sign(t - 10.000000 + 0.250000/A) + 1)*np.cos(omega*t)* 0.5 * (np.sign(10.000000+0.250000/A-t) + 1)+    0.5 * (np.sign(t - 110.000000 + 0.250000/A)+1)*np.cos(omega*t)*0.5 * (np.sign(110.000000 + 0.250000/A-t)+1)
        coeff[1] = 0.5 * (np.sign(t - 60.000000 + 0.500000/A) + 1)*np.cos(omega*t)*0.5 * (np.sign(60.000000+0.500000/A - t) + 1)
