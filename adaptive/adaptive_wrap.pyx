import numpy as np
cimport numpy as np
import remapping
import sys

cdef extern:
    void build_grid_adaptive(double *h, double *t, double *s,
                             double *dz_a, double *dz_d,
                             double *dk_sig, double *di_sig,
                             int nx, int ny, int nz,
                             double alpha, double max_depth)

    void build_grid_adaptive_tend(double *h, double *t, double *s,
                                  double *dz_a, double *dz_p,
                                  double *dz_i, double *dz_p_i,
                                  double *dk_sig,
                                  int nx, int ny, int nz,
                                  double alpha)

    void remap_wrap(int scheme, int n0, double *h0, double *u0,
                    int n1, double *h1, double *u1)

def build_grid(double[:,:,:] h not None,
               double[:,:,:] t not None,
               double[:,:,:] s not None,
               double alpha, double max_depth):
    cdef int nx, ny, nz
    nx = h.shape[0]
    ny = h.shape[1]
    nz = h.shape[2]

    cdef double[::1,:,:] h_f = h.copy_fortran()
    cdef double[::1,:,:] t_f = t.copy_fortran()
    cdef double[::1,:,:] s_f = s.copy_fortran()

    dz_a = np.empty((nx, ny, nz+1), dtype=np.double, order='F')
    dz_d = np.empty((nx, ny, nz+1), dtype=np.double, order='F')
    dz_p = np.empty((nx, ny, nz+1), dtype=np.double, order='F')
    cdef double[::1,:,:] res_a = dz_a
    cdef double[::1,:,:] res_d = dz_d
    cdef double[::1,:,:] res_p = dz_p

    dk_sig = np.empty((nx, ny, nz), dtype=np.double, order='F')
    di_sig = np.empty((nx+1, ny, nz-1), dtype=np.double, order='F')
    di_pre = np.empty((nx+1, ny, nz-1), dtype=np.double, order='F')
    cdef double[::1,:,:] res_dk = dk_sig
    cdef double[::1,:,:] res_di = di_sig
    cdef double[::1,:,:] res_di_p = di_pre

    # build_grid_adaptive(&h_f[0,0,0], &t_f[0,0,0], &s_f[0,0,0],
    #                     &res_a[0,0,0], &res_d[0,0,0],
    #                     &res_dk[0,0,0], &res_di[0,0,0],
    #                     nx, ny, nz, alpha, max_depth)
    # return dz_a, dz_d, dk_sig, di_sig

    build_grid_adaptive_tend(&h_f[0,0,0], &t_f[0,0,0], &s_f[0,0,0],
                                &res_a[0,0,0], &res_p[0,0,0],
                                &res_di[0,0,0], &res_di_p[0,0,0],
                                &res_dk[0,0,0],
                                nx, ny, nz, 0.0)

    return dz_a, dz_p, dk_sig, di_sig, di_pre

def build_grid_i(np.ndarray[np.double_t, ndim=3] h not None,
                 np.ndarray[np.double_t, ndim=3] t not None,
                 np.ndarray[np.double_t, ndim=3] s not None,
                 int n_iter, int scheme):
    cdef int nx, ny, nz
    nx = h.shape[0]
    ny = h.shape[1]
    nz = h.shape[2]

    cdef double[::1,:,:] h_f = np.asfortranarray(h)
    cdef double[::1,:,:] t_f = np.asfortranarray(t)
    cdef double[::1,:,:] s_f = np.asfortranarray(s)
    h_tmp = np.empty_like(h, order='F')
    cdef double[::1,:,:] h_t = h_tmp

    dz_a = np.empty((nx, ny, nz+1), dtype=np.double, order='F')
    dz_p = np.empty((nx, ny, nz+1), dtype=np.double, order='F')
    cdef double[::1,:,:] res_a = dz_a
    cdef double[::1,:,:] res_p = dz_p

    dk_sig = np.empty((nx, ny, nz), dtype=np.double, order='F')
    di_sig = np.empty((nx+1, ny, nz-1), dtype=np.double, order='F')
    di_pre = np.empty((nx+1, ny, nz-1), dtype=np.double, order='F')
    cdef double[::1,:,:] res_dk = dk_sig
    cdef double[::1,:,:] res_di = di_sig
    cdef double[::1,:,:] res_di_p = di_pre

    ta1 = np.empty(nz, dtype=np.double)
    ta2 = np.empty(nz, dtype=np.double)
    ta3 = np.empty(nz, dtype=np.double)
    ta4 = np.empty(nz, dtype=np.double)
    cdef double[:] h0 = ta1, u0 = ta2, h1 = ta3, u1 = ta4

    for i in range(n_iter):
        if i % 10 == 0:
            sys.stdout.write('.')
        if i % 1000 == 0:
            print 'iteration {}'.format(i)

        build_grid_adaptive_tend(&h_f[0,0,0], &t_f[0,0,0], &s_f[0,0,0],
                                 &res_a[0,0,0], &res_p[0,0,0],
                                 &res_di[0,0,0], &res_di_p[0,0,0],
                                 &res_dk[0,0,0],
                                 nx, ny, nz, 0.0)

        h_tmp = h - np.diff(dz_a)
        t = remapping.remap(scheme, h_f, t_f, h_t)
        s = remapping.remap(scheme, h_f, s_f, h_t)
        h[:] = h_tmp

    return h_f, t_f, s_f
