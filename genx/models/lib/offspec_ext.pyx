include "complex.pxd"
import numpy as np
cimport numpy as np
from libc.math cimport fabs, pow, exp
from genx.models.lib.complex cimport cpow, cexp, conj, creal

def dwba_interdiff_sum(double[::1] qx, complex[:, :, ::1] G,
                       complex[:, :, ::1] q, double eta,
                       double h, complex[::1] sigma,
                       complex[::1] sigmaid, complex[::1] sqn,
                       double[::1] z, complex[::1] table,
                       double[::1] q_min, double[::1] fn,
                       double eta_z):


    cdef double xnew = 0
    cdef int lower = 0
    cdef complex s = 0
    cdef complex ret_temp = 0
    cdef double[::1] ret = np.zeros(qx.shape[0], dtype=float)
    cdef int m, i, j, k, l, p

    for m in range(qx.shape[0]):
        ret_temp = 0
        for i in range(G.shape[1]):
            for j in range(G.shape[1]):
                for k in range(4):
                    for l in range(4):
                        for p in range(fn.shape[0]):
                            xnew=fabs(qx[m]*eta/pow(p+1.0,1.0/2.0/h))
                            #xnew=qx[m]*eta/pow(p+1.0,1.0/2.0/h)
                            lower=int((xnew-q_min[0])/(q_min[1]-q_min[0]))
                            s += 2.0*cpow(q[k, i, m]*conj(q[l, j, m])*sigma[i]*sigma[j]*exp(-fabs(z[i]-z[j])/eta_z),p+1.0)\
                                 /fn[p]*eta/pow(p+1.0,1.0/2.0/h)*((table[lower+1]-table[lower])/(q_min[lower+1]-q_min[lower])*(xnew-q_min[lower])+table[lower])


                        ret_temp += (sqn[i]-sqn[i+1])*conj(sqn[j]-sqn[j+1])*G[k, i, m]*cexp(-0.5*cpow(sigmaid[i]*q[k, i, m],2.0))\
                                 *conj(G[l, j, m]*cexp(-0.5*cpow(sigmaid[j]*q[l, j, m],2.0)))\
                                 *cexp(-0.5*(cpow(q[k, i, m]*sigma[i],2.0)+cpow(conj(q[l, j, m])*sigma[j],2.0)))/q[k, i, m]/conj(q[l, j, m])*s
                        #*cexp(-I*q[k, i, m]*z[i]-conj(I*q[l, j, m])*z[j])
                        s = 0

        ret[m] = creal(ret_temp)
    return np.asarray(ret)

def dwba_sum(double[::1] qx, complex[:, :, ::1] G,
             complex[:, :, ::1] q, double eta,
             double h, complex[::1] sigma,
             complex[::1] sqn, double[::1] z,
             complex[::1] table, double[::1] q_min,
             double[::1] fn, double eta_z):

    cdef double xnew = 0
    cdef int lower = 0
    cdef complex s = 0
    cdef complex ret_temp = 0
    cdef double[::1] ret = np.zeros(G.shape[2], dtype=float)
    cdef int m, i, j, k, l, p


    for m in range(G.shape[2]):
        ret_temp = 0
        for i in range(G.shape[1]):
            for j in range(G.shape[1]):
                for k in range(4):
                    for l in range(4):
                        for p in range(fn.shape[0]):
                            xnew = fabs(qx[m]*eta/pow(p+1.0,1.0/2.0/h))
                            lower=int((xnew-q_min[0])/(q_min[1]-q_min[0]))
                            s+=2.0*cpow(q[k, i, m]*conj(q[l, j, m])*sigma[i]*sigma[j]*exp(-fabs(z[i]-z[j])/eta_z),p+1.0)\
                               /fn[p]*eta/pow(p+1.0,1.0/2.0/h)*((table[lower+1]-table[lower])/(q_min[lower+1]-q_min[lower])*(xnew-q_min[lower])+table[lower])

                        ret_temp += (sqn[i]-sqn[i+1])*conj(sqn[j]-sqn[j+1])*G[k, i, m]*conj(G[l, j, m])\
                                 *cexp(-0.5*(cpow(q[k, i, m]*sigma[i],2.0)+cpow(conj(q[l, j, m])*sigma[j],2.0)))\
                                 /q[k, i, m]/conj(q[l, j, m])*s
                        #*cexp(-I*q[k, i, m]*z[i]-conj(I*q[l, j, m])*z[j])
                        s = 0

        ret[m] = creal(ret_temp)
    return np.asarray(ret)

def born_sum(double[::1] qx, double[::1] qz,
             double eta, double h, complex[::1] sigma,
             complex[::1] sqn, double[::1] z,
             complex[::1] table, double[::1] q_min,
             double[::1] fn, double eta_z):

    cdef double xnew=0
    cdef int lower=0
    cdef complex s = 0
    cdef complex ret_temp = 0
    cdef double[::1] ret = np.zeros(qx.shape[0], dtype=float)
    cdef int m, i, j, p


    for m in range(qx.shape[0]):
        ret_temp = 0
        for i in range(sigma.shape[0]):
            for j in range(sigma.shape[0]):
                for p in range(fn.shape[0]):
                    xnew=qx[m]*eta/pow(p+1.0,1.0/2.0/h)
                    lower=int((xnew-q_min[0])/(q_min[1]-q_min[0]))
                    s += 2.0*cpow(qz[m]*qz[m]*sigma[i]*sigma[j]*exp(-fabs(z[i]-z[j])/eta_z),p+1.0)\
                         /fn[p]*eta/pow(p+1.0,1.0/2.0/h)*((table[lower+1]-table[lower])/(q_min[lower+1]-q_min[lower])*(xnew-q_min[lower])+table[lower])


                ret_temp+=(sqn[i]-sqn[i+1])*conj(sqn[j]-sqn[j+1])*cexp(-0.5*(cpow(qz[m]*sigma[i],2.0)+cpow(qz[m]*sigma[j],2.0)))/qz[m]/qz[m]*s
                s = 0

        ret[m] = creal(ret_temp)

    return np.asarray(ret)
