module adaptive
  use MOM_EOS, only : calculate_density_derivs, EOS_WRIGHT, EOS_LINEAR
  use iso_c_binding, only : c_double, c_int

  implicit none
  private

  real, parameter :: H_subroundoff = 1e-20 * 1e-10
  real, parameter :: H_to_Pa = 9.8 * 1035.0

  logical, parameter :: limiter = .false.

contains
  subroutine build_adapt_column(z, dz, nz, max_depth)
    real, dimension(nz+1), intent(in) :: z
    real, dimension(nz+1), intent(out) :: dz
    integer, intent(in) :: nz
    real, intent(in) :: max_depth

    integer :: k
    real :: b1, d1, b_denom_1
    real, dimension(nz) :: k_grid, c1
    real, dimension(nz+1) :: z_next

    ! no initial movement of interfaces
    dz(:) = 0.0
    ! copy interfaces
    z_next(:) = z(:)

    do k = 1, nz
      ! background component only
      k_grid(k) = nz**2 / max_depth
    end do

    ! tridiagonal solver
    b1 = 1.0
    d1 = 1.0
    do K = 2, nz
      b_denom_1 = 1. + d1 * k_grid(k-1)
      b1 = 1.0 / (b_denom_1 + k_grid(k))
      c1(K) = k_grid(k) * b1
      d1 = b_denom_1 * b1
      z_next(K) = b1 * (z_next(K) + k_grid(k-1) * z_next(K-1))
    end do
    do K = nz, 2, -1
      z_next(K) = z_next(K) + c1(K) * z_next(K+1)
      ! figure out how much this interface was moved
      dz(K) = z_next(K) - z(K)
      ! don't allow tangling or movement of bottom
      !z(K) = min(z(K), z(K+1))
    end do

  end subroutine build_adapt_column

  subroutine build_grid_adaptive(h, t, s, dz_a, dz_d, dk_sig, di_sig, &
       nx, ny, nz, alpha, max_depth) bind(c)
    real(c_double), dimension(nx,ny,nz), intent(in) :: h, t, s
    real(c_double), dimension(nx,ny,nz+1), intent(out) :: dz_a, dz_d
    integer(c_int), intent(in), value :: nx, ny, nz
    real(c_double), intent(in), value :: alpha, max_depth
    real(c_double), dimension(nx,ny,nz), intent(out) :: dk_sig
    real(c_double), dimension(0:nx,ny,2:nz), intent(out) :: di_sig

    integer :: i, j, k
    integer :: eos
    real, dimension(nx,ny) :: t_int, t_int_kp1, s_int, s_int_kp1
    real, dimension(nx,ny,nz+1) :: z_int
    real, dimension(nx,ny) :: alpha_int, alpha_int_kp1
    real, dimension(nx,ny) :: beta_int, beta_int_kp1
    real, dimension(nx,0:ny) :: dj_sig
    ! calculated curvature (divergence of flux of d*_sig)
    real :: curv

    !eos = EOS_WRIGHT
    eos = EOS_LINEAR

    ! initially no grid movement
    dz_a(:,:,:) = 0.
    dz_d(:,:,:) = 0.

    ! set halo values to zero
    di_sig(:,:,:) = 0.
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
           alpha_int(:,j), beta_int(:,j), 1, nx, eos, 1000.0, -0.2, 0.8)
      call calculate_density_derivs(t_int_kp1(:,j), s_int_kp1(:,j), z_int(:,j,2) * H_to_Pa, &
           alpha_int_kp1(:,j), beta_int_kp1(:,j), 1, nx, eos, 1000.0, -0.2, 0.8)

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
             alpha_int_kp1(:,j), beta_int_kp1(:,j), 1, nx, eos, 1000.0, -0.2, 0.8)

        do i = 1, nx
          dk_sig(i,j,k) = 0.5 * (alpha_int(i,j) + alpha_int_kp1(i,j)) * (t_int_kp1(i,j) - t_int(i,j)) &
               + 0.5 * (beta_int(i,j) + beta_int_kp1(i,j)) * (s_int_kp1(i,j) - s_int(i,j))
          dk_sig(i,j,k) = max(dk_sig(i,j,k), 0.)
        enddo
      enddo

      ! u-points
      do j = 1, ny
        do I = 1, nx-1
          di_sig(I,j,K) = 0.5 * (alpha_int(i,j) + alpha_int(i+1,j)) * (t_int(i+1,j) - t_int(i,j)) &
               + 0.5 * (beta_int(i,j) + beta_int(i+1,j)) * (s_int(i+1,j) - s_int(i,j))

          if (limiter) then
            ! if di_sig is positive, the interface to the right is denser
            ! (and the interface to the left is lighter)
            ! this will cause the right interface to seek lighter water by moving
            ! upwards, and vice versa on the left
            if (di_sig(I,j,K) > 0.) then
              ! limit downward motion to the left (thickness below left interface)
              di_sig(I,j,K) = min(di_sig(I,j,K), h(i,j,k) * dk_sig(i,j,k))
              ! limit upward motion to the right (thickness above right interface)
              di_sig(I,j,K) = min(di_sig(I,j,K), h(i+1,j,k-1) * dk_sig(i+1,j,k-1))
            else
              ! limit upward motion to the left (thickness above left interface)
              di_sig(I,j,K) = max(di_sig(I,j,K), -h(i,j,k-1) * dk_sig(i,j,k-1))
              ! limit downward motion to the right (thickness below right interface)
              di_sig(I,j,K) = max(di_sig(I,j,K), -h(i+1,j,k) * dk_sig(i+1,j,k))
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
          curv = (di_sig(I,j,K) - di_sig(I-1,j,K)) + (dj_sig(i,J) - dj_sig(i,J-1))

          ! choose dk_sig based on sign
          if (curv < 0.) then
            ! negative curvature -- interface should move upward (dz negative)
            if (dk_sig(i,j,k-1) > 0.) dz_a(i,j,K) = curv * h(i,j,k-1) / &
                 sqrt(dk_sig(i,j,k-1)**2 + curv**2)
          else
            ! positive curvature -- interface should move downward (dz positive)
            if (dk_sig(i,j,k) > 0.) dz_a(i,j,K) = curv * h(i,j,k) / &
                 sqrt(dk_sig(i,j,k)**2 + curv**2)
          endif

          z_int(i,j,K) = z_int(i,j,K) + alpha * dz_a(i,j,K)
        enddo
      enddo
    enddo

    do j = 1, ny
      do i = 1, nx
        call build_adapt_column(z_int(i,j,:), dz_d(i,j,:), nz, max_depth)
      enddo
    enddo
  end subroutine build_grid_adaptive

  subroutine build_grid_adaptive_tend(h, t, s, dz_a, dz_p, dz_i, dz_p_i, &
       dk_sig, nx, ny, nz, alpha) bind(c)
    real(c_double), dimension(nx,ny,nz), intent(in) :: h, t, s
    real(c_double), dimension(nx,ny,nz+1), intent(out) :: dz_a, dz_p
    integer(c_int), intent(in), value :: nx, ny, nz
    real(c_double), intent(in), value :: alpha

    real(c_double), dimension(0:nx,ny,2:nz), intent(out) :: dz_i, dz_p_i
    real(c_double), dimension(nx,ny,nz), intent(out) :: dk_sig

    integer :: i, j, k
    integer :: eos
    real, dimension(nx,ny) :: t_int, t_int_kp1, s_int, s_int_kp1
    real, dimension(nx,ny,nz+1) :: z_int
    real, dimension(nx,ny) :: alpha_int, alpha_int_kp1
    real, dimension(nx,ny) :: beta_int, beta_int_kp1
    real, dimension(nx,0:ny) :: dz_j, dz_p_j
    real :: di_sig, dj_sig, di_pre, dj_pre
    real :: extent_left, extent_right

    eos = EOS_LINEAR

    ! initially no grid movement
    dz_a(:,:,:) = 0.
    dz_p(:,:,:) = 0.

    ! set halo values to zero
    dz_i(:,:,:) = 0.
    dz_j(:,:) = 0.

    dz_p_i(:,:,:) = 0.
    dz_p_j(:,:) = 0.

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
           alpha_int(:,j), beta_int(:,j), 1, nx, eos, 1000.0, -0.2, 0.8)
      call calculate_density_derivs(t_int_kp1(:,j), s_int_kp1(:,j), z_int(:,j,2) * H_to_Pa, &
           alpha_int_kp1(:,j), beta_int_kp1(:,j), 1, nx, eos, 1000.0, -0.2, 0.8)

      do i = 1, nx
        dk_sig(i,j,1) = 0.5 * (alpha_int(i,j) + alpha_int_kp1(i,j)) * (t_int_kp1(i,j) - t(i,j,1)) &
             + 0.5 * (beta_int(i,j) + beta_int_kp1(i,j)) * (s_int_kp1(i,j) - s(i,j,1))
        dk_sig(i,j,1) = max(dk_sig(i,j,1), 0.)
      enddo
    enddo

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
             alpha_int_kp1(:,j), beta_int_kp1(:,j), 1, nx, eos, 1000.0, -0.2, 0.8)

        do i = 1, nx
          dk_sig(i,j,k) = 0.5 * (alpha_int(i,j) + alpha_int_kp1(i,j)) * (t_int_kp1(i,j) - t_int(i,j)) &
               + 0.5 * (beta_int(i,j) + beta_int_kp1(i,j)) * (s_int_kp1(i,j) - s_int(i,j))
          dk_sig(i,j,k) = max(dk_sig(i,j,k), 0.)
        enddo
      enddo

      ! u-points
      do j = 1, ny
        do I = 1, nx-1
          ! calculate the density difference across the cell boundary
          di_sig = 0.5 * (alpha_int(i,j) + alpha_int(i+1,j)) * (t_int(i+1,j) - t_int(i,j)) &
               + 0.5 * (beta_int(i,j) + beta_int(i+1,j)) * (s_int(i+1,j) - s_int(i,j))

          if (di_sig < 0) then
            ! left is denser than right
            ! left moves up, right moves down
            !extent_left = h(i,j,k-1) / sqrt(dk_sig(i,j,k-1)**2 + di_sig**2)
            !extent_right = h(i+1,j,k) / sqrt(dk_sig(i+1,j,k)**2 + di_sig**2)
            extent_left = dk_sig(i,j,k-1) / (h(i,j,k-1) + H_subroundoff)
            extent_right = dk_sig(i+1,j,k) / (h(i+1,j,k) + H_subroundoff)
          else
            ! right is denser than left
            ! left moves down, right moves up
            !extent_left = h(i,j,k) / sqrt(dk_sig(i,j,k)**2 + di_sig**2)
            !extent_right = h(i+1,j,k-1) / sqrt(dk_sig(i+1,j,k-1)**2 + di_sig**2)
            extent_left = dk_sig(i,j,k) / (h(i,j,k) + H_subroundoff)
            extent_right = dk_sig(i+1,j,k-1) / (h(i+1,j,k-1) + H_subroundoff)
          end if

          ! we solve for the local change in slope (dz) required to achieve
          ! a change in density of (-alpha * di_sig)
          dz_i(I,j,K) = -0.5 * di_sig / (extent_left + extent_right)

          ! limit slope based on adjacent layers
          if (dz_i(I,j,k) > 0) then
            ! dz_i positive
            dz_i(I,j,K) = min(dz_i(I,j,K), 0.5 * h(i,j,k-1))
            dz_i(I,j,K) = min(dz_i(I,j,K), 0.5 * h(i+1,j,k))
          else
            ! dz_i negative
            dz_i(I,j,K) = max(dz_i(I,j,K), -0.5 * h(i,j,k))
            dz_i(I,j,K) = max(dz_i(I,j,K), -0.5 * h(i+1,j,k-1))
          end if

          ! we also calculate the difference in pressure (interface position)
          di_pre = 0.5 * (z_int(i+1,j,K) - z_int(i,j,K))

          ! limiter
          if (di_pre < 0) then
            ! right higher than left -- left up, right down
            dz_p_i(I,j,K) = max(di_pre, -0.5 * h(i,j,k-1))
            dz_p_i(I,j,K) = max(di_pre, -0.5 * h(i+1,j,k))
          else
            ! left down, right up
            dz_p_i(I,j,K) = min(di_pre, 0.5 * h(i,j,k))
            dz_p_i(I,j,K) = min(di_pre, 0.5 * h(i+1,j,k-1))
          end if
        end do
      end do

      do j = 1, ny
        do i = 1, nx
          dz_a(i,j,K) = (dz_i(I,j,K) - dz_i(I-1,j,K)) + (dz_j(i,J) - dz_j(i,J-1))
          dz_p(i,j,K) = (dz_p_i(I,j,K) - dz_p_i(I-1,j,K)) + (dz_p_j(i,J) - dz_p_j(i,J-1))
        end do
      end do
    end do
  end subroutine build_grid_adaptive_tend
end module adaptive
