import numpy as np

cdef extern:
    void build_grid_adaptive(double *h, double *t, double *s, double *dz,
                            int nx, int ny, int nz)

# expects Fortran-contiguous arrays
def build_grid(double[::1,:,:] h not None,
               double[::1,:,:] t not None,
               double[::1,:,:] s not None):
    cdef int nx, ny, nz
    nx = h.shape[0]
    ny = h.shape[1]
    nz = h.shape[2]

    dz = np.empty((nx, ny, nz+1), dtype=np.double, order='F')
    cdef double[::1,:,:] res = dz

    build_grid_adaptive(&h[0,0,0], &t[0,0,0], &s[0,0,0], &res[0,0,0],
                        nx, ny, nz)

    return dz
