include "complex.pxd"
import numpy as np
cimport numpy as np
from libc.math cimport exp
from genx.models.lib.complex cimport cexp, I

def surface_lattice_sum(double[::1] x, double[::1] y, double[::1] z,
                        double[::1] h, double[::1] k, double[::1] l,
                        double[::1] u, double[::1] oc, complex[::1, :] f,
                        double[:, :, ::1] Pt, double[::1] dinv):

    cdef double pi = 3.14159265358979311599796346854418516159057617187500
    cdef complex tmp
    cdef complex[::1] ret = np.zeros(h.shape[0], dtype=complex)
    cdef int i, j, m

    #   Loop over all data points
    for i in range(h.shape[0]):
        # Loop over all atoms
        for j in range(oc.shape[0]):
            #Loop over symmetry operations
            tmp = 0.0
            for m in range(Pt.shape[0]):
                tmp += cexp(2.0 * pi * I * (h[i] * (Pt[m, 0, 0] * x[j] + Pt[m, 0, 1]* y[j] + Pt[m, 0, 2]) +
                                            k[i] * (Pt[m, 1, 0] * x[j] + Pt[m, 1, 1] * y[j] + Pt[m, 1, 2]) +
                                            l[i] * z[j]))
            ret[i] += oc[j] * f[i, j] * exp(-2.0 * (pi * dinv[i])**2.0 * u[j]) * tmp


    return np.asarray(ret)
