import numpy as np

cdef extern:
    void remap_wrap(int scheme, int n0, double *h0, double *u0,
                    int n1, double *h1, double *u1)

def remap(int scheme,
          double[::1,:,:] h0 not None,
          double[::1,:,:] u0 not None,
          double[::1,:,:] h1 not None):
    cdef unsigned int nx, ny, n0, n1
    cdef unsigned int i, j

    nx = h0.shape[0]
    ny = h0.shape[1]
    n0 = h0.shape[2]
    n1 = h1.shape[2]

    u1 = np.empty((nx, ny, n1), dtype=np.double, order='F')
    cdef double[::1,:,:] res = u1

    ta1 = np.empty(n0, dtype=np.double)
    ta2 = np.empty(n0, dtype=np.double)
    ta3 = np.empty(n1, dtype=np.double)
    ta4 = np.empty(n1, dtype=np.double)
    cdef double[:] t1 = ta1, t2 = ta2, t3 = ta3, t4 = ta4

    for j in range(ny):
        for i in range(nx):
            t1[:] = h0[i,j,:]
            t2[:] = u0[i,j,:]
            t3[:] = h1[i,j,:]

            remap_wrap(scheme, n0, &t1[0], &t2[0], n1, &t3[0], &t4[0])

            u1[i,j,:] = t4

    return u1
