import numpy as np
cimport numpy as np
cimport cython
from libc.stdlib cimport malloc, free
from libc.math cimport sqrt
from libc.stdint cimport int32_t, uint32_t
from cython.parallel import prange


cdef extern from "mtrandom.h":
    cdef cppclass MTrandoms:
        void seedMT()
        double random0i1i() nogil
        double random0i1e() nogil


cdef MTrandoms rmt
rmt.seedMT()


@cython.boundscheck(False)
def density2points(int32_t n, int32_t npm, double dg, float[:] density, double[:] r0, int32_t ng):
    cdef int32_t ix, iy, iz, offset
    cdef double x0, y0, z0

    cdef np.ndarray[dtype=np.float64_t, ndim=2] _r = np.empty((n, 3), dtype=np.float64, order='C')
    cdef double[:] gd = np.arange(-npm, npm, dtype=np.float64) * dg

    x0 = r0[0]
    y0 = r0[1]
    z0 = r0[2]

    n = 0
    for ix in range(-npm, npm):
        offset = ng * (ng * (ix + npm)) + npm
        for iy in range(-npm, npm):
            for iz in range(-npm, npm):
                if density[iz + offset] > 0:
                    _r[n, 0] = gd[ix + npm] + x0
                    _r[n, 1] = gd[iy + npm] + y0
                    _r[n, 2] = gd[iz + npm] + z0
                    n += 1
            offset += ng
    return _r


@cython.cdivision(True)
@cython.boundscheck(False)
cdef inline c_random_distances(av1, av2, double* distances, uint32_t nSamples):

    cdef uint32_t i, i1, i2
    cdef int32_t lp1, lp2

    cdef double[:, :] p1 = av1.points
    cdef double[:, :] p2 = av2.points

    lp1 = p1.shape[0]
    lp2 = p2.shape[0]

    for i in prange(nSamples, nogil=True):
        i1 = <int>(rmt.random0i1e() * lp1)
        i2 = <int>(rmt.random0i1e() * lp2)
        distances[i] = sqrt(
            (p1[i1, 0] - p2[i2, 0]) * (p1[i1, 0] - p2[i2, 0]) + \
            (p1[i1, 1] - p2[i2, 1]) * (p1[i1, 1] - p2[i2, 1]) + \
            (p1[i1, 2] - p2[i2, 2]) * (p1[i1, 2] - p2[i2, 2])
        )

def random_distances(av1, av2, uint32_t nSamples=10000):
    cdef np.ndarray[ndim=1, dtype=np.float64_t] dist = np.empty(nSamples, dtype=np.float64)
    cdef double* dist_buff = <double*> dist.data
    c_random_distances(av1, av2, dist_buff, nSamples)
    return dist


@cython.cdivision(True)
@cython.boundscheck(False)
def RDAMeanE(av1, av2, double R0=52.0, uint32_t nSamples=50000):
    cdef uint32_t i
    cdef double Esum = 0.0
    cdef double* d = <double*>malloc(nSamples * sizeof(double))
    c_random_distances(av1, av2, d, nSamples)
    for i in prange(nSamples, nogil=True):
        Esum += (1./(1.+(d[i]/R0)**6.0))
    Esum /= nSamples
    free(d)
    return (1./Esum - 1.)**(1./6.) * R0


@cython.cdivision(True)
@cython.boundscheck(False)
def RDAMean(av1, av2, uint32_t nSamples=50000):
    cdef uint32_t i
    cdef double RDA = 0.0

    cdef double* d = <double*>malloc(nSamples * sizeof(double))
    c_random_distances(av1, av2, d, nSamples)
    for i in prange(nSamples, nogil=True):
        RDA += d[i]
    free(d)
    return RDA / nSamples



