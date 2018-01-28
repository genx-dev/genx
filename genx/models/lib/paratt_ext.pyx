include "complex.pxd"
import numpy as np
cimport numpy as np
from libc.math cimport pow, cos
from genx.models.lib.complex cimport I, cexp, cpow, conj, cabs

def refl(double[::1] theta, double lamda, complex[::1] n,
                        double[::1] d, double[::1] sigma):

    cdef double k = 2*3.141592/lamda
    cdef double torad = 3.141592/180.0
    cdef double costheta2

    cdef complex Qj, Qj1, r, rp, p
    cdef complex n2amb = cpow(n[sigma.shape[0] - 1], 2.0)

    cdef double[::1] ret = np.zeros(theta.shape[0], dtype=float)

    cdef int i, j

    for i in range(theta.shape[0]):
        costheta2 = pow(cos(theta[i]*torad), 2.0)
        Qj = 2.0 * k * csqrt(n[0]*n[0] - n2amb*costheta2)
        Qj1 = 2.0*k*csqrt(n[1]*n[1] - n2amb*costheta2)
        rp=(Qj - Qj1)/(Qj1 + Qj)*cexp(-Qj1*Qj/2.0*sigma[0]*sigma[0])
        r = rp
        for j in range(1, sigma.shape[0] - 1):
            Qj = Qj1
            Qj1 = 2.0*k*csqrt(n[j+1]*n[j+1] - n2amb*costheta2)
            rp = (Qj-Qj1)/(Qj1+Qj)*cexp(-Qj1*Qj/2.0*sigma[j]*sigma[j])
            p = cexp(I*d[j]*Qj)
            r = (rp+r*p)/(1.0 + rp*r*p)

        ret[i] = cabs(r * conj(r))

    return np.asarray(ret)

def reflq(double[::1] Q, double lamda, complex[::1] n,
                        double[::1] d, double[::1] sigma):

    cdef double Q02 = 16*3.141592*3.141592/lamda/lamda
    cdef complex Q2, Qj, Qj1, r, rp, p
    cdef complex n2amb = cpow(n[sigma.shape[0] - 1], 2.0)

    cdef double[::1] ret = np.zeros(Q.shape[0], dtype=float)

    cdef int i, j

    for i in range(Q.shape[0]):
        Q2 = Q[i]*Q[i]
        Qj = csqrt((n[0]*n[0] - n2amb)* Q02 + n2amb*Q2)
        Qj1 = csqrt((n[1]*n[1] - n2amb)*Q02 + n2amb*Q2)
        rp=(Qj - Qj1)/(Qj1 + Qj)*cexp(-Qj1*Qj/2.0*sigma[0]*sigma[0])
        r = rp
        for j in range(1, sigma.shape[0] - 1):
            Qj = Qj1
            Qj1 = csqrt((n[j+1]*n[j+1] - n2amb)*Q02 + n2amb*Q2)
            rp = (Qj-Qj1)/(Qj1+Qj)*cexp(-Qj1*Qj/2.0*sigma[j]*sigma[j])
            p = cexp(I*d[j]*Qj)
            r = (rp+r*p)/(1.0 + rp*r*p)

        ret[i] = cabs(r * conj (r))

    return np.asarray(ret)

def refl_nvary2(double[::1] theta, double[::1] lamda,
        complex[:, ::1] n, double[::1] d,
                double[::1] sigma):


    cdef double costheta2, k
    cdef double pi = 3.141592
    cdef double torad = pi/180.0
    cdef complex Qj, Qj1, r, rp, p, n2amb

    cdef double[::1] ret = np.zeros(theta.shape[0], dtype=float)

    cdef int i, j
    
    for i in range(theta.shape[0]):
        costheta2 = pow(cos(theta[i]*torad), 2.0)
        n2amb = cpow(n[sigma.shape[0] - 1, i], 2.0)
        k = 4*pi/lamda[i]
        Qj = k*csqrt(n[0, i]*n[0, i]  - n2amb*costheta2)
        Qj1 = k*csqrt(n[1, i]*n[1, i] - n2amb*costheta2)
        rp = (Qj - Qj1)/(Qj1 + Qj)*cexp(-Qj1*Qj/2.0*sigma[0]*sigma[0])
        r = rp
        for j in range(1, sigma.shape[0] - 1):
            Qj = Qj1
            Qj1 = k*csqrt(n[j+1, i]*n[j+1, i] - n2amb*costheta2)
            rp = (Qj - Qj1)/(Qj1 + Qj)*cexp(-Qj1*Qj/2.0*sigma[j]*sigma[j])
            p = cexp(I*d[j]*Qj)
            r = (rp+r*p)/(1.0+rp*r*p)

        ret[i] = cabs(r*conj(r))

    return np.asarray(ret)

def refl_nvary2_nosigma(double[::1] theta, double[::1] lamda,
        complex[:, ::1] n, double[::1] d):

    cdef double costheta2, k
    cdef double pi = 3.141592
    cdef double torad = pi/180.0
    cdef complex Qj, Qj1, r, rp, p, n2amb

    cdef double[::1] ret = np.zeros(theta.shape[0], dtype=float)

    cdef int i, j

    for i in range(theta.shape[0]):
        costheta2 = pow(cos(theta[i]*torad), 2.0)
        n2amb = cpow(n[d.shape[0]-1, i], 2.0)
        k = 4*pi/lamda[i]
        Qj = k*csqrt(n[0, i]*n[0, i]  - n2amb*costheta2)
        Qj1 = k*csqrt(n[1, i]*n[1, i] - n2amb*costheta2)
        rp = (Qj - Qj1)/(Qj1 + Qj)
        r = rp
        for j in range(d.shape[0] - 1):
            Qj = Qj1
            Qj1 = k*csqrt(n[j+1, i]*n[j+1, i] - n2amb*costheta2)
            rp = (Qj - Qj1)/(Qj1 + Qj)
            p = cexp(I*d[j]*Qj)
            r = (rp+r*p)/(1.0+rp*r*p)

        ret[i] = cabs(r*conj(r))

    return np.asarray(ret)
