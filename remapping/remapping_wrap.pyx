import numpy as np

cdef extern:
    void remap_wrap(int scheme, int n0, double *h0, double *u0,
                    int n1, double *h1, double *u1)

def remap(int scheme,
          double[:] h0 not None,
          double[:] u0 not None,
          double[:] h1 not None):
    cdef int n0, n1
    n0 = len(h0)
    n1 = len(h1)

    u1 = np.empty(n1, dtype=np.double)
    cdef double[:] res = u1

    remap_wrap(scheme, n0, &h0[0], &u0[0], n1, &h1[0], &res[0])

    return u1
