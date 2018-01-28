cdef extern from "<complex.h>" nogil:
    double complex I
    double complex cexp(double complex z)
    double complex csqrt(double complex z)
    double complex conj(double complex z)
    double cabs(double complex z)
    double creal(double complex z)
    double complex cpow(double complex x, double complex y)
