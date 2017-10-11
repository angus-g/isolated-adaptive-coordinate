module adaptive
  use MOM_EOS, only : calculate_density_derivs, EOS_WRIGHT
  use iso_c_binding, only : c_double, c_int

  implicit none
  private

  real, parameter :: H_subroundoff = 1e-20 * 1e-10
  real, parameter :: H_to_Pa = 9.8 * 1035.0

  logical, parameter :: limiter = .false.

contains

  subroutine build_grid_adaptive(h, t, s, dz, nx, ny, nz) bind(c)
    real(c_double), dimension(nx,ny,nz), intent(in) :: h, t, s
    real(c_double), dimension(nx,ny,nz+1), intent(out) :: dz
    integer(c_int), intent(in), value :: nx, ny, nz

    integer :: i, j, k
    integer :: eos
    real, dimension(nx,ny) :: t_int, t_int_kp1, s_int, s_int_kp1
    real, dimension(nx,ny,nz+1) :: z_int
    real, dimension(nx,ny) :: alpha_int, alpha_int_kp1
    real, dimension(nx,ny) :: beta_int, beta_int_kp1
    real, dimension(nx,ny,nz) :: dk_sig
    real, dimension(0:nx,0:ny+1) :: di_sig
    real, dimension(0:nx+1,0:ny) :: dj_sig
    ! calculated curvature (divergence of flux of d*_sig)
    real :: curv

    eos = EOS_WRIGHT

    ! initially no grid movement
    dz(:,:,:) = 0.

    ! set halo values to zero
    di_sig(:,:) = 0.
    dj_sig(:,:) = 0.

    ! set surface and second interface
    z_int(:,:,1) = 0.
    z_int(:,:,2) = h(:,:,1)

    ! populate data ahead of current interface
    do j = 1, ny
      do i = 1, nx
        t_int_kp1(i,j) = ( &
             t(i,j,1) * (h(i,j,2) + H_subroundoff) + &
             t(i,j,2) * (h(i,j,1) + H_subroundoff)) / &
             (h(i,j,1) + h(i,j,2) + 2*H_subroundoff)
        s_int_kp1(i,j) = ( &
             s(i,j,1) * (h(i,j,2) + H_subroundoff) + &
             s(i,j,2) * (h(i,j,1) + H_subroundoff)) / &
             (h(i,j,1) + h(i,j,2) + 2*H_subroundoff)
      enddo

      call calculate_density_derivs(t(:,j,1), s(:,j,1), z_int(:,j,1) * H_to_Pa, &
           alpha_int(:,j), beta_int(:,j), 1, nx, eos)
      call calculate_density_derivs(t_int_kp1(:,j), s_int_kp1(:,j), z_int(:,j,2) * H_to_Pa, &
           alpha_int_kp1(:,j), beta_int_kp1(:,j), 1, nx, eos)

      do i = 1, nx
        dk_sig(i,j,1) = 0.5 * (alpha_int(i,j) + alpha_int_kp1(i,j)) * (t_int_kp1(i,j) - t(i,j,1)) &
             + 0.5 * (beta_int(i,j) + beta_int_kp1(i,j)) * (s_int_kp1(i,j) - s(i,j,1))
        dk_sig(i,j,1) = max(dk_sig(i,j,1), 0.)
      enddo
    enddo

    ! work on the rest of the interior interfaces (top-down)
    do K = 2, nz
      do j = 1, ny
        do i = 1, nx
          ! calculate next interface position
          z_int(i,j,K+1) = z_int(i,j,K) + h(i,j,k)

          ! copy from interface ahead
          t_int(i,j) = t_int_kp1(i,j)
          s_int(i,j) = s_int_kp1(i,j)
          alpha_int(i,j) = alpha_int_kp1(i,j)
          beta_int(i,j) = beta_int_kp1(i,j)

          if (k == nz) then
            ! use constant value for bottom interface
            t_int_kp1(i,j) = t(i,j,nz)
            s_int_kp1(i,j) = s(i,j,nz)
          else
            ! calculate ahead one interface (interior)
            t_int_kp1(i,j) = ( &
                 t(i,j,k) * (h(i,j,k+1) + H_subroundoff) + &
                 t(i,j,k+1) * (h(i,j,k) + H_subroundoff)) / &
                 (h(i,j,k) + h(i,j,k+1) + 2*H_subroundoff)
            s_int_kp1(i,j) = ( &
                 s(i,j,k) * (h(i,j,k+1) + H_subroundoff) + &
                 s(i,j,k+1) * (h(i,j,k) + H_subroundoff)) / &
                 (h(i,j,k) + h(i,j,k+1) + 2*H_subroundoff)
          endif
        enddo

        call calculate_density_derivs(t_int_kp1(:,j), s_int_kp1(:,j), z_int(:,j,K+1) * H_to_Pa, &
             alpha_int_kp1(:,j), beta_int_kp1(:,j), 1, nx, eos)

        do i = 1, nx
          dk_sig(i,j,k) = 0.5 * (alpha_int(i,j) + alpha_int_kp1(i,j)) * (t_int_kp1(i,j) - t_int(i,j)) &
               + 0.5 * (beta_int(i,j) + beta_int_kp1(i,j)) * (s_int_kp1(i,j) - s_int(i,j))
          dk_sig(i,j,k) = max(dk_sig(i,j,k), 0.)
        enddo
      enddo

      ! u-points
      do j = 1, ny
        do I = 1, nx-1
          di_sig(I,j) = 0.5 * (alpha_int(i,j) + alpha_int(i+1,j)) * (t_int(i+1,j) - t_int(i,j)) &
               + 0.5 * (beta_int(i,j) + beta_int(i+1,j)) * (s_int(i+1,j) - s_int(i,j))

          if (limiter) then
            ! if di_sig is positive, the interface to the right is denser
            ! (and the interface to the left is lighter)
            ! this will cause the right interface to seek lighter water by moving
            ! upwards, and vice versa on the left
            if (di_sig(I,j) > 0.) then
              ! limit downward motion to the left (thickness below left interface)
              di_sig(I,j) = min(di_sig(I,j), h(i,j,k) * dk_sig(i,j,k))
              ! limit upward motion to the right (thickness above right interface)
              di_sig(I,j) = min(di_sig(I,j), h(i+1,j,k-1) * dk_sig(i+1,j,k-1))
            else
              ! limit upward motion to the left (thickness above left interface)
              di_sig(I,j) = max(di_sig(I,j), -h(i,j,k-1) * dk_sig(i,j,k-1))
              ! limit downward motion to the right (thickness below right interface)
              di_sig(I,j) = max(di_sig(I,j), -h(i+1,j,k) * dk_sig(i+1,j,k))
            endif
          endif
        enddo
      enddo

      ! v-points
      do J = 1, ny-1
        do i = 1, nx
          dj_sig(i,J) = 0.5 * (alpha_int(i,j) + alpha_int(i,j+1)) * (t_int(i,j+1) - t_int(i,j)) &
               + 0.5 * (beta_int(i,j) + beta_int(i,j+1)) * (s_int(i,j+1) - s_int(i,j))

          if (limiter) then
            ! if dj_sig is positive, the interface to the north is denser
            ! (and the interface to the south is lighter)
            ! this will cause the north interface to seek lighter water by moving
            ! upwards, and vice versa to the south
            if (dj_sig(i,J) > 0.) then
              ! limit downward motion to the south (thickness below south interface)
              dj_sig(i,J) = min(dj_sig(i,J), h(i,j,k) * dk_sig(i,j,k))
              ! limit upward motion to the north (thickness above north interface)
              dj_sig(i,J) = min(dj_sig(i,J), h(i,j+1,k-1) * dk_sig(i,j+1,k-1))
            else
              ! limit upward motion to the south (thickness above south interface)
              dj_sig(i,J) = max(dj_sig(i,J), -h(i,j,k-1) * dk_sig(i,j,k-1))
              ! limit downward motion to the north (thickness below north interface)
              dj_sig(i,J) = max(dj_sig(i,J), -h(i,j+1,k) * dk_sig(i,j+1,k))
            endif
          endif
        enddo
      enddo

      do j = 1, ny
        do i = 1, nx
          ! calculate divergence of flux
          curv = (di_sig(I,j) - di_sig(I-1,j)) + (dj_sig(i,J) - dj_sig(i,J-1))

          ! choose dk_sig based on sign
          if (curv < 0.) then
            if (dk_sig(i,j,k-1) > 0.) dz(i,j,K) = curv * h(i,j,k-1) / &
                 sqrt(dk_sig(i,j,k-1)**2 + curv**2)
          else
            if (dk_sig(i,j,k) > 0.) dz(i,j,K) = curv * h(i,j,k) / &
                 sqrt(dk_sig(i,j,k)**2 + curv**2)
          endif
        enddo
      enddo
    enddo

  end subroutine build_grid_adaptive
end module adaptive
