!> Provides column-wise vertical remapping functions
module MOM_remapping

! This file is part of MOM6. See LICENSE.md for the license.
! Original module written by Laurent White, 2008.06.09

use MOM_error_handler, only : MOM_error, FATAL
use regrid_edge_values, only : edge_values_explicit_h4, edge_values_implicit_h4
use regrid_edge_values, only : edge_values_implicit_h4, edge_values_implicit_h6
use regrid_edge_slopes, only : edge_slopes_implicit_h3, edge_slopes_implicit_h5
use PCM_functions, only : PCM_reconstruction
use PLM_functions, only : PLM_reconstruction, PLM_boundary_extrapolation
use PPM_functions, only : PPM_reconstruction, PPM_boundary_extrapolation
use PQM_functions, only : PQM_reconstruction, PQM_boundary_extrapolation_v1

implicit none ; private

!> Container for remapping parameters
type, public :: remapping_CS
  !> Determines which reconstruction to use
  integer :: remapping_scheme = -911
  !> Degree of polynomial reconstruction
  integer :: degree = 0
  !> If true, extrapolate boundaries
  logical :: boundary_extrapolation = .true.
  !> If true, reconstructions are checked for consistency.
  logical :: check_reconstruction = .false.
  !> If true, the result of remapping are checked for conservation and bounds.
  logical :: check_remapping = .false.
  !> If true, the intermediate values used in remapping are forced to be bounded.
  logical :: force_bounds_in_subcell = .false.
end type

! The following routines are visible to the outside world
public remapping_core_h

! The following are private parameter constants
integer, parameter  :: REMAPPING_PCM        = 0 !< O(h^1) remapping scheme
integer, parameter  :: REMAPPING_PLM        = 1 !< O(h^2) remapping scheme
integer, parameter  :: REMAPPING_PPM_H4     = 2 !< O(h^3) remapping scheme
integer, parameter  :: REMAPPING_PPM_IH4    = 3 !< O(h^3) remapping scheme
integer, parameter  :: REMAPPING_PQM_IH4IH3 = 4 !< O(h^4) remapping scheme
integer, parameter  :: REMAPPING_PQM_IH6IH5 = 5 !< O(h^5) remapping scheme

integer, parameter  :: INTEGRATION_PCM = 0  !< Piecewise Constant Method
integer, parameter  :: INTEGRATION_PLM = 1  !< Piecewise Linear Method
integer, parameter  :: INTEGRATION_PPM = 3  !< Piecewise Parabolic Method
integer, parameter  :: INTEGRATION_PQM = 5  !< Piecewise Quartic Method

! This CPP macro turns on/off bounding of integrations limits so that they are
! always within the cell. Roundoff can lead to the non-dimensional bounds being
! outside of the range 0 to 1.
#define __USE_ROUNDOFF_SAFE_ADJUSTMENTS__

real, parameter :: h_neglect = 1.E-30 !< A dimensional (H units) number that can be
                                      !! added to thicknesses in a denominator without
                                      !! changing the numerical result, except where
                                      !! a division by zero would otherwise occur.

logical, parameter :: old_algorithm = .false. !< Use the old "broken" algorithm.
                                              !! This is a temporary measure to assist
                                              !! debugging until we delete the old algorithm.

contains

!> Remaps column of values u0 on grid h0 to grid h1
!! assuming the top edge is aligned.
subroutine remapping_core_h(CS,  n0, h0, u0, n1, h1, u1)
  type(remapping_CS),  intent(in)  :: CS !< Remapping control structure
  integer,             intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0), intent(in)  :: h0 !< Cell widths on source grid
  real, dimension(n0), intent(in)  :: u0 !< Cell averages on source grid
  integer,             intent(in)  :: n1 !< Number of cells on target grid
  real, dimension(n1), intent(in)  :: h1 !< Cell widths on target grid
  real, dimension(n1), intent(out) :: u1 !< Cell averages on target grid
  ! Local variables
  integer :: iMethod
  real, dimension(n0,2)           :: ppoly_r_E            !Edge value of polynomial
  real, dimension(n0,2)           :: ppoly_r_S            !Edge slope of polynomial
  real, dimension(n0,CS%degree+1) :: ppoly_r_coefficients !Coefficients of polynomial
  integer :: k
  real :: eps, h0tot, h0err, h1tot, h1err, u0tot, u0err, u0min, u0max, u1tot, u1err, u1min, u1max, uh_err

  call build_reconstructions_1d( CS, n0, h0, u0, ppoly_r_coefficients, ppoly_r_E, ppoly_r_S, iMethod )

  if (CS%check_reconstruction) call check_reconstructions_1d(n0, h0, u0, CS%degree, &
                                   CS%boundary_extrapolation, ppoly_r_coefficients, ppoly_r_E, ppoly_r_S)


  call remap_via_sub_cells( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients, n1, h1, iMethod, &
                            CS%force_bounds_in_subcell, u1, uh_err )

  if (CS%check_remapping) then
    ! Check errors and bounds
    call measure_input_bounds( n0, h0, u0, ppoly_r_E, h0tot, h0err, u0tot, u0err, u0min, u0max )
    call measure_output_bounds( n1, h1, u1, h1tot, h1err, u1tot, u1err, u1min, u1max )
    if (iMethod<5) then ! We except PQM until we've debugged it
    if ( (abs(u1tot-u0tot)>(u0err+u1err)+uh_err .and. abs(h1tot-h0tot)<h0err+h1err) &
        .or. (u1min<u0min .or. u1max>u0max) ) then
      write(0,*) 'iMethod = ',iMethod
      write(0,*) 'H: h0tot=',h0tot,'h1tot=',h1tot,'dh=',h1tot-h0tot,'h0err=',h0err,'h1err=',h1err
      if (abs(h1tot-h0tot)>h0err+h1err) write(0,*) 'H non-conservation difference=',h1tot-h0tot,'allowed err=',h0err+h1err,' <-----!'
      write(0,*) 'UH: u0tot=',u0tot,'u1tot=',u1tot,'duh=',u1tot-u0tot,'u0err=',u0err,'u1err=',u1err,'uh_err=',uh_err
      if (abs(u1tot-u0tot)>(u0err+u1err)+uh_err) write(0,*) 'U non-conservation difference=',u1tot-u0tot,'allowed err=',u0err+u1err+uh_err,' <-----!'
      write(0,*) 'U: u0min=',u0min,'u1min=',u1min
      if (u1min<u0min) write(0,*) 'U minimum overshoot=',u1min-u0min,' <-----!'
      write(0,*) 'U: u0max=',u0max,'u1max=',u1max
      if (u1max>u0max) write(0,*) 'U maximum overshoot=',u1max-u0max,' <-----!'
      write(0,'(a3,6a24)') 'k','h0','left edge','u0','right edge','h1','u1'
      do k = 1, max(n0,n1)
        if (k<=min(n0,n1)) then
          write(0,'(i3,1p6e24.16)') k,h0(k),ppoly_r_E(k,1),u0(k),ppoly_r_E(k,2),h1(k),u1(k)
        elseif (k>n0) then
          write(0,'(i3,96x,1p2e24.16)') k,h1(k),u1(k)
        else
          write(0,'(i3,1p4e24.16)') k,h0(k),ppoly_r_E(k,1),u0(k),ppoly_r_E(k,2)
        endif
      enddo
      write(0,'(a3,2a24)') 'k','u0','Polynomial coefficients'
      do k = 1, n0
        write(0,'(i3,1p6e24.16)') k,u0(k),ppoly_r_coefficients(k,:)
      enddo
      call MOM_error( FATAL, 'MOM_remapping, remapping_core_h: '//&
             'Remapping result is inconsistent!' )
    endif
    endif ! method<5
  endif

end subroutine remapping_core_h

!> Creates polynomial reconstructions of u0 on the source grid h0.
subroutine build_reconstructions_1d( CS, n0, h0, u0, ppoly_r_coefficients, ppoly_r_E, ppoly_r_S, iMethod )
  type(remapping_CS),              intent(in)  :: CS
  integer,                         intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0),             intent(in)  :: h0 !< Cell widths on source grid
  real, dimension(n0),             intent(in)  :: u0 !< Cell averages on source grid
  real, dimension(n0,CS%degree+1), intent(out) :: ppoly_r_coefficients !< Coefficients of polynomial
  real, dimension(n0,2),           intent(out) :: ppoly_r_E !< Edge value of polynomial
  real, dimension(n0,2),           intent(out) :: ppoly_r_S !< Edge slope of polynomial
  integer,                         intent(out) :: iMethod !< Integration method
  ! Local variables
  integer :: local_remapping_scheme

  ! Reset polynomial
  ppoly_r_E(:,:) = 0.0
  ppoly_r_S(:,:) = 0.0
  ppoly_r_coefficients(:,:) = 0.0
  iMethod = -999

  local_remapping_scheme = CS%remapping_scheme
  if (n0<=1) then
    local_remapping_scheme = REMAPPING_PCM
  elseif (n0<=3) then
    local_remapping_scheme = min( local_remapping_scheme, REMAPPING_PLM )
  elseif (n0<=4) then
    local_remapping_scheme = min( local_remapping_scheme, REMAPPING_PPM_H4 )
  endif
  select case ( local_remapping_scheme )
    case ( REMAPPING_PCM )
      call PCM_reconstruction( n0, u0, ppoly_r_E, ppoly_r_coefficients)
      iMethod = INTEGRATION_PCM
    case ( REMAPPING_PLM )
      call PLM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients )
      if ( CS%boundary_extrapolation ) then
        call PLM_boundary_extrapolation( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients)
      end if
      iMethod = INTEGRATION_PLM
    case ( REMAPPING_PPM_H4 )
      call edge_values_explicit_h4( n0, h0, u0, ppoly_r_E )
      call PPM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients )
      if ( CS%boundary_extrapolation ) then
        call PPM_boundary_extrapolation( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients )
      end if
      iMethod = INTEGRATION_PPM
    case ( REMAPPING_PPM_IH4 )
      call edge_values_implicit_h4( n0, h0, u0, ppoly_r_E )
      call PPM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients )
      if ( CS%boundary_extrapolation ) then
        call PPM_boundary_extrapolation( n0, h0, u0, ppoly_r_E, ppoly_r_coefficients )
      end if
      iMethod = INTEGRATION_PPM
    case ( REMAPPING_PQM_IH4IH3 )
      call edge_values_implicit_h4( n0, h0, u0, ppoly_r_E )
      call edge_slopes_implicit_h3( n0, h0, u0, ppoly_r_S )
      call PQM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_S, ppoly_r_coefficients )
      if ( CS%boundary_extrapolation ) then
        call PQM_boundary_extrapolation_v1( n0, h0, u0, ppoly_r_E, ppoly_r_S, ppoly_r_coefficients )
      end if
      iMethod = INTEGRATION_PQM
    case ( REMAPPING_PQM_IH6IH5 )
      call edge_values_implicit_h6( n0, h0, u0, ppoly_r_E )
      call edge_slopes_implicit_h5( n0, h0, u0, ppoly_r_S )
      call PQM_reconstruction( n0, h0, u0, ppoly_r_E, ppoly_r_S, ppoly_r_coefficients )
      if ( CS%boundary_extrapolation ) then
        call PQM_boundary_extrapolation_v1( n0, h0, u0, ppoly_r_E, ppoly_r_S, ppoly_r_coefficients )
      end if
      iMethod = INTEGRATION_PQM
    case default
      call MOM_error( FATAL, 'MOM_remapping, build_reconstructions_1d: '//&
           'The selected remapping method is invalid' )
  end select

end subroutine build_reconstructions_1d

!> Checks that edge values and reconstructions satisfy bounds
subroutine check_reconstructions_1d(n0, h0, u0, deg, boundary_extrapolation, &
                                    ppoly_r_coefficients, ppoly_r_E, ppoly_r_S)
  integer,                  intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0),      intent(in)  :: h0 !< Cell widths on source grid
  real, dimension(n0),      intent(in)  :: u0 !< Cell averages on source grid
  integer,                  intent(in)  :: deg !< Degree of polynomial reconstruction
  logical,                  intent(in)  :: boundary_extrapolation !< Extrapolate at boundaries if true
  real, dimension(n0,deg+1),intent(out) :: ppoly_r_coefficients !< Coefficients of polynomial
  real, dimension(n0,2),    intent(out) :: ppoly_r_E !< Edge value of polynomial
  real, dimension(n0,2),    intent(out) :: ppoly_r_S !< Edge slope of polynomial
  ! Local variables
  integer :: i0, n
  real :: u_l, u_c, u_r ! Cell averages
  real :: u_min, u_max
  logical :: problem_detected

  problem_detected = .false.
  do i0 = 1, n0
    u_l = u0(max(1,i0-1))
    u_c = u0(i0)
    u_r = u0(min(n0,i0+1))
    if (i0 > 1 .or. .not. boundary_extrapolation) then
      u_min = min(u_l, u_c)
      u_max = max(u_l, u_c)
      if (ppoly_r_E(i0,1) < u_min) then
        write(0,'(a,i4,5(x,a,1pe24.16))') 'Left edge undershoot at',i0,'u(i0-1)=',u_l,'u(i0)=',u_c, &
                                          'edge=',ppoly_r_E(i0,1),'err=',ppoly_r_E(i0,1)-u_min
        problem_detected = .true.
      endif
      if (ppoly_r_E(i0,1) > u_max) then
        write(0,'(a,i4,5(x,a,1pe24.16))') 'Left edge overshoot at',i0,'u(i0-1)=',u_l,'u(i0)=',u_c, &
                                          'edge=',ppoly_r_E(i0,1),'err=',ppoly_r_E(i0,1)-u_max
        problem_detected = .true.
      endif
    endif
    if (i0 < n0 .or. .not. boundary_extrapolation) then
      u_min = min(u_c, u_r)
      u_max = max(u_c, u_r)
      if (ppoly_r_E(i0,2) < u_min) then
        write(0,'(a,i4,5(x,a,1pe24.16))') 'Right edge undershoot at',i0,'u(i0)=',u_c,'u(i0+1)=',u_r, &
                                          'edge=',ppoly_r_E(i0,2),'err=',ppoly_r_E(i0,2)-u_min
        problem_detected = .true.
      endif
      if (ppoly_r_E(i0,2) > u_max) then
        write(0,'(a,i4,5(x,a,1pe24.16))') 'Right edge overshoot at',i0,'u(i0)=',u_c,'u(i0+1)=',u_r, &
                                          'edge=',ppoly_r_E(i0,2),'err=',ppoly_r_E(i0,2)-u_max
        problem_detected = .true.
      endif
    endif
    if (i0 > 1) then
      if ( (u_c-u_l)*(ppoly_r_E(i0,1)-ppoly_r_E(i0-1,2)) < 0.) then
        write(0,'(a,i4,5(x,a,1pe24.16))') 'Non-monotonic edges at',i0,'u(i0-1)=',u_l,'u(i0)=',u_c, &
                                          'right edge=',ppoly_r_E(i0-1,2),'left edge=',ppoly_r_E(i0,1)
        write(0,'(5(a,1pe24.16,x))') 'u(i0)-u(i0-1)',u_c-u_l,'edge diff=',ppoly_r_E(i0,1)-ppoly_r_E(i0-1,2)
        problem_detected = .true.
      endif
    endif
    if (problem_detected) then
      write(0,'(a,1p9e24.16)') 'Polynomial coeffs:',ppoly_r_coefficients(i0,:)
      write(0,'(3(a,1pe24.16,x))') 'u_l=',u_l,'u_c=',u_c,'u_r=',u_r
      write(0,'(a4,10a24)') 'i0','h0(i0)','u0(i0)','left edge','right edge','Polynomial coefficients'
      do n = 1, n0
        write(0,'(i4,1p10e24.16)') n,h0(n),u0(n),ppoly_r_E(n,1),ppoly_r_E(n,2),ppoly_r_coefficients(n,:)
      enddo
      call MOM_error(FATAL, 'MOM_remapping, check_reconstructions_1d: '// &
                   'Edge values or polynomial coefficients were inconsistent!')
    endif
  enddo

end subroutine check_reconstructions_1d

!> Remaps column of n0 values u0 on grid h0 to grid h1 with n1 cells by calculating
!! the n0+n1+1 sub-integrals of the intersection of h0 and h1, and the summing the
!! appropriate integrals into the h1*u1 values.
subroutine remap_via_sub_cells( n0, h0, u0, ppoly0_E, ppoly0_coefficients, n1, h1, method, &
                                force_bounds_in_subcell, u1, uh_err, ah_sub, aisub_src, aiss, aise )
  integer,           intent(in)    :: n0     !< Number of cells in source grid
  real,              intent(in)    :: h0(n0)  !< Source grid widths (size n0)
  real,              intent(in)    :: u0(n0)  !< Source cell averages (size n0)
  real,              intent(in)    :: ppoly0_E(n0,2)            !< Edge value of polynomial
  real,              intent(in)    :: ppoly0_coefficients(:,:) !< Coefficients of polynomial
  integer,           intent(in)    :: n1     !< Number of cells in target grid
  real,              intent(in)    :: h1(n1)  !< Target grid widths (size n1)
  integer,           intent(in)    :: method !< Remapping scheme to use
  logical,           intent(in)    :: force_bounds_in_subcell !< Force sub-cell values to be bounded
  real,              intent(out)   :: u1(n1)  !< Target cell averages (size n1)
  real,              intent(out)   :: uh_err !< Estimate of bound on error in sum of u*h
  real, optional,    intent(out)   :: ah_sub(n0+n1+1) !< h_sub
  integer, optional, intent(out)   :: aisub_src(n0+n1+1) !< i_sub_src
  integer, optional, intent(out)   :: aiss(n0) !< isrc_start
  integer, optional, intent(out)   :: aise(n0) !< isrc_ens
  ! Local variables
  integer :: i_sub ! Index of sub-cell
  integer :: i0 ! Index into h0(1:n0), source column
  integer :: i1 ! Index into h1(1:n1), target column
  integer :: i_start0 ! Used to record which sub-cells map to source cells
  integer :: i_start1 ! Used to record which sub-cells map to target cells
  integer :: i_max ! Used to record which sub-cell is the largest contribution of a source cell
  real :: dh_max ! Used to record which sub-cell is the largest contribution of a source cell
  real, dimension(n0+n1+1) :: h_sub ! Width of each each sub-cell
  real, dimension(n0+n1+1) :: uh_sub ! Integral of u*h over each sub-cell
  real, dimension(n0+n1+1) :: u_sub ! Average of u over each sub-cell
  integer, dimension(n0+n1+1) :: isub_src ! Index of source cell for each sub-cell
  integer, dimension(n0) :: isrc_start ! Index of first sub-cell within each source cell
  integer, dimension(n0) :: isrc_end ! Index of last sub-cell within each source cell
  integer, dimension(n0) :: isrc_max ! Index of thickest sub-cell within each source cell
  real, dimension(n0) :: h0_eff ! Effective thickness of source cells
  real, dimension(n0) :: u0_min ! Minimum value of reconstructions in source cell
  real, dimension(n0) :: u0_max ! Minimum value of reconstructions in source cell
  integer, dimension(n1) :: itgt_start ! Index of first sub-cell within each target cell
  integer, dimension(n1) :: itgt_end ! Index of last sub-cell within each target cell
  real :: xa, xb ! Non-dimensional position within a source cell (0..1)
  real :: h0_supply, h1_supply ! The amount of width available for constructing sub-cells
  real :: dh ! The width of the sub-cell
  real :: duh ! The total amount of accumulated stuff (u*h)
  real :: dh0_eff ! Running sum of source cell thickness
  ! For error checking/debugging
  logical, parameter :: force_bounds_in_target = .true. ! To fix round-off issues
  logical, parameter :: adjust_thickest_subcell = .true. ! To fix round-off conservation issues
  logical, parameter :: debug_bounds = .false. ! For debugging overshoots etc.
  integer :: k, i0_last_thick_cell
  real :: h0tot, h0err, h1tot, h1err, h2tot, h2err, u02_err
  real :: u0tot, u0err, u0min, u0max, u1tot, u1err, u1min, u1max, u2tot, u2err, u2min, u2max, u_orig
  logical :: src_has_volume !< True if h0 has not been consumed
  logical :: tgt_has_volume !< True if h1 has not been consumed

  if (old_algorithm) isrc_max(:)=1

  i0_last_thick_cell = 0
  do i0 = 1, n0
    u0_min(i0) = min(ppoly0_E(i0,1), ppoly0_E(i0,2))
    u0_max(i0) = max(ppoly0_E(i0,1), ppoly0_E(i0,2))
    if (h0(i0)>0.) i0_last_thick_cell = i0
  enddo

  ! Initialize algorithm
  h0_supply = h0(1)
  h1_supply = h1(1)
  src_has_volume = .true.
  tgt_has_volume = .true.
  i0 = 1 ; i1 = 1
  i_start0 = 1 ; i_start1 = 1
  i_max = 1
  dh_max = 0.
  dh0_eff = 0.

  ! First sub-cell is always vanished
  h_sub(1) = 0.
  isrc_start(1) = 1
  isrc_end(1) = 1
  isrc_max(1) = 1
  isub_src(1) = 1

  ! Loop over each sub-cell to calculate intersections with source and target grids
  do i_sub = 2, n0+n1+1

    ! This is the width of the sub-cell, determined by which ever column has the least
    ! supply available to consume.
    dh = min(h0_supply, h1_supply)

    ! This is the running sum of the source cell thickness. After summing over each
    ! sub-cell, the sum of sub-cell thickness might differ from the original source
    ! cell thickness due to round off.
    dh0_eff = dh0_eff + min(dh, h0_supply)

    ! Record the source index (i0) that this sub-cell integral belongs to. This
    ! is needed to index the reconstruction coefficients for the source cell
    ! used in the integrals of the sub-cell width.
    isub_src(i_sub) = i0
    h_sub(i_sub) = dh

    ! For recording the largest sub-cell within a source cell.
    if (dh >= dh_max) then
      i_max = i_sub
      dh_max = dh
    endif

    ! Which ever column (source or target) has the least width left to consume determined
    ! the width, dh, of sub-cell i_sub in the expression for dh above.
    if (h0_supply <= h1_supply .and. src_has_volume) then
      ! h0_supply is smaller than h1_supply) so we consume h0_supply and increment the
      ! source cell index.
      h1_supply = h1_supply - dh ! Although this is a difference the result will
                                 ! be non-negative because of the conditional.
      ! Record the sub-cell start/end index that span the source cell i0.
      isrc_start(i0) = i_start0
      isrc_end(i0) = i_sub
      i_start0 = i_sub + 1
      ! Record the sub-cell that is the largest fraction of the source cell.
      isrc_max(i0) = i_max
      i_max = i_sub + 1
      dh_max = 0.
      ! Record the source cell thickness found by summing the sub-cell thicknesses.
      h0_eff(i0) = dh0_eff
      ! Move the source index.
      if (old_algorithm) then
        if (i0 < i0_last_thick_cell) then
          i0 = i0 + 1
          h0_supply = h0(i0)
          dh0_eff = 0.
          do while (h0_supply==0. .and. i0<i0_last_thick_cell)
            ! This loop skips over vanished source cells
            i0 = i0 + 1
            h0_supply = h0(i0)
          enddo
        else
          h0_supply = 1.E30
        endif
      else
        if (i0 < n0) then
          i0 = i0 + 1
          h0_supply = h0(i0)
          dh0_eff = 0.
        else
          h0_supply = 0.
          src_has_volume = .false.
        endif
      endif
    elseif (h0_supply >= h1_supply .and. tgt_has_volume) then
      ! h1_supply is smaller than h0_supply) so we consume h1_supply and increment the
      ! target cell index.
      h0_supply = h0_supply - dh ! Although this is a difference the result will
                                 ! be non-negative because of the conditional.
      ! Record the sub-cell start/end index that span the target cell i1.
      itgt_start(i1) = i_start1
      itgt_end(i1) = i_sub
      i_start1 = i_sub + 1
      ! Move the target index.
      if (i1 < n1) then
        i1 = i1 + 1
        h1_supply = h1(i1)
      else
        if (old_algorithm) then
          h1_supply = 1.E30
        else
          h1_supply = 0.
          tgt_has_volume = .false.
        endif
      endif
    elseif (src_has_volume) then
      ! We ran out of target volume but still have source cells to consume
      h_sub(i_sub) = h0_supply
      ! Record the sub-cell start/end index that span the source cell i0.
      isrc_start(i0) = i_start0
      isrc_end(i0) = i_sub
      i_start0 = i_sub + 1
      ! Record the sub-cell that is the largest fraction of the source cell.
      isrc_max(i0) = i_max
      i_max = i_sub + 1
      dh_max = 0.
      ! Record the source cell thickness found by summing the sub-cell thicknesses.
      h0_eff(i0) = dh0_eff
      if (i0 < n0) then
        i0 = i0 + 1
        h0_supply = h0(i0)
        dh0_eff = 0.
      else
        h0_supply = 0.
        src_has_volume = .false.
      endif
    elseif (tgt_has_volume) then
      ! We ran out of source volume but still have target cells to consume
      h_sub(i_sub) = h1_supply
      ! Record the sub-cell start/end index that span the target cell i1.
      itgt_start(i1) = i_start1
      itgt_end(i1) = i_sub
      i_start1 = i_sub + 1
      ! Move the target index.
      if (i1 < n1) then
        i1 = i1 + 1
        h1_supply = h1(i1)
      else
        h1_supply = 0.
        tgt_has_volume = .false.
      endif
    else
      stop 'remap_via_sub_cells: THIS SHOULD NEVER HAPPEN!'
    endif

  enddo

  ! Loop over each sub-cell to calculate average/integral values within each sub-cell.
  xa = 0.
  dh0_eff = 0.
  uh_sub(1) = 0.
  u_sub(1) = ppoly0_E(1,1)
  u02_err = 0.
  do i_sub = 2, n0+n1

    ! Sub-cell thickness from loop above
    dh = h_sub(i_sub)

    ! Source cell
    i0 = isub_src(i_sub)

    ! Evaluate average and integral for sub-cell i_sub.
    ! Integral is over distance dh but expressed in terms of non-dimensional
    ! positions with source cell from xa to xb  (0 <= xa <= xb <= 1).
    dh0_eff = dh0_eff + dh ! Cumulative thickness within the source cell
    if (h0_eff(i0)>0.) then
      xb = dh0_eff / h0_eff(i0) ! This expression yields xa <= xb <= 1.0
      xb = min(1., xb) ! This is only needed when the total target column is wider than the source column
      u_sub(i_sub) = average_value_ppoly( n0, u0, ppoly0_E, ppoly0_coefficients, method, i0, xa, xb)
    else ! Vanished cell
      xb = 1.
      u_sub(i_sub) = u0(i0)
    endif
    if (debug_bounds) then
      if (method<5 .and.(u_sub(i_sub)<u0_min(i0) .or. u_sub(i_sub)>u0_max(i0))) then
        write(0,*) 'Sub cell average is out of bounds',i_sub,'method=',method
        write(0,*) 'xa,xb: ',xa,xb
        write(0,*) 'Edge values: ',ppoly0_E(i0,:),'mean',u0(i0)
        write(0,*) 'a_c: ',(u0(i0)-ppoly0_E(i0,1))+(u0(i0)-ppoly0_E(i0,2))
        write(0,*) 'Polynomial coeffs: ',ppoly0_coefficients(i0,:)
        write(0,*) 'Bounds min=',u0_min(i0),'max=',u0_max(i0)
        write(0,*) 'Average: ',u_sub(i_sub),'rel to min=',u_sub(i_sub)-u0_min(i0),'rel to max=',u_sub(i_sub)-u0_max(i0)
        call MOM_error( FATAL, 'MOM_remapping, remap_via_sub_cells: '//&
             'Sub-cell average is out of bounds!' )
      endif
    endif
    if (force_bounds_in_subcell) then
      ! These next two lines should not be needed but when using PQM we found roundoff
      ! can lead to overshoots. These lines sweep issues under the rug which need to be
      ! properly .. later. -AJA
      u_orig = u_sub(i_sub)
      u_sub(i_sub) = max( u_sub(i_sub), u0_min(i0) )
      u_sub(i_sub) = min( u_sub(i_sub), u0_max(i0) )
      u02_err = u02_err + dh*abs( u_sub(i_sub) - u_orig )
    endif
    uh_sub(i_sub) = dh * u_sub(i_sub)

    if (isub_src(i_sub+1) /= i0) then
      ! If the next sub-cell is in a different source cell, reset the position counters
      dh0_eff = 0.
      xa = 0.
    else
      xa = xb ! Next integral will start at end of last
    endif

  enddo
  u_sub(n0+n1+1) = ppoly0_E(n0,2)                   ! This value is only needed when total target column
  uh_sub(n0+n1+1) = ppoly0_E(n0,2) * h_sub(n0+n1+1) ! is wider than the source column

  if (adjust_thickest_subcell) then
    ! Loop over each source cell substituting the integral/average for the thickest sub-cell (within
    ! the source cell) with the residual of the source cell integral minus the other sub-cell integrals
    ! aka a genius algorithm for accurate conservation when remapping from Robert Hallberg (@Hallberg-NOAA).
    do i0 = 1, i0_last_thick_cell
      i_max = isrc_max(i0)
      dh_max = h_sub(i_max)
      if (dh_max > 0.) then
        ! duh will be the sum of sub-cell integrals within the source cell except for the thickest sub-cell.
        duh = 0.
        do i_sub = isrc_start(i0), isrc_end(i0)
          if (i_sub /= i_max) duh = duh + uh_sub(i_sub)
        enddo
        uh_sub(i_max) = u0(i0)*h0(i0) - duh
        u02_err = u02_err + max( abs(uh_sub(i_max)), abs(u0(i0)*h0(i0)), abs(duh) )
      endif
    enddo
  endif

  ! Loop over each target cell summing the integrals from sub-cells within the target cell.
  uh_err = 0.
  do i1 = 1, n1
    if (h1(i1) > 0.) then
      duh = 0. ; dh = 0.
      i_sub = itgt_start(i1)
      if (force_bounds_in_target) then
        u1min = u_sub(i_sub)
        u1max = u_sub(i_sub)
      endif
      do i_sub = itgt_start(i1), itgt_end(i1)
        if (force_bounds_in_target) then
          u1min = min(u1min, u_sub(i_sub))
          u1max = max(u1max, u_sub(i_sub))
        endif
        dh = dh + h_sub(i_sub)
        duh = duh + uh_sub(i_sub)
        ! This accumulates the contribution to the error bound for the sum of u*h
        uh_err = uh_err + max(abs(duh),abs(uh_sub(i_sub)))*epsilon(duh)
      enddo
      u1(i1) = duh / dh
      ! This is the contribution from the division to the error bound for the sum of u*h
      uh_err = uh_err + abs(duh)*epsilon(duh)
      if (force_bounds_in_target) then
        u_orig = u1(i1)
        u1(i1) = max(u1min, min(u1max, u1(i1)))
        ! Adjusting to be bounded contributes to the error for the sum of u*h
        uh_err = uh_err + dh*abs( u1(i1)-u_orig )
      endif
    else
      u1(i1) = u_sub(itgt_start(i1))
    endif
  enddo

  ! Check errors and bounds
  if (debug_bounds) then
    call measure_input_bounds( n0, h0, u0, ppoly0_E, h0tot, h0err, u0tot, u0err, u0min, u0max )
    call measure_output_bounds( n1, h1, u1, h1tot, h1err, u1tot, u1err, u1min, u1max )
    call measure_output_bounds( n0+n1+1, h_sub, u_sub, h2tot, h2err, u2tot, u2err, u2min, u2max )
    if (method<5) then ! We except PQM until we've debugged it
    if (     (abs(u1tot-u0tot)>(u0err+u1err)+uh_err+u02_err .and. abs(h1tot-h0tot)<h0err+h1err) &
        .or. (abs(u2tot-u0tot)>u0err+u2err+u02_err .and. abs(h2tot-h0tot)<h0err+h2err) &
        .or. (u1min<u0min .or. u1max>u0max) ) then
      write(0,*) 'method = ',method
      write(0,*) 'Source to sub-cells:'
      write(0,*) 'H: h0tot=',h0tot,'h2tot=',h2tot,'dh=',h2tot-h0tot,'h0err=',h0err,'h2err=',h2err
      if (abs(h2tot-h0tot)>h0err+h2err) write(0,*) 'H non-conservation difference=',h2tot-h0tot,'allowed err=',h0err+h2err,' <-----!'
      write(0,*) 'UH: u0tot=',u0tot,'u2tot=',u2tot,'duh=',u2tot-u0tot,'u0err=',u0err,'u2err=',u2err,'adjustment err=',u02_err
      if (abs(u2tot-u0tot)>u0err+u2err) write(0,*) 'U non-conservation difference=',u2tot-u0tot,'allowed err=',u0err+u2err,' <-----!'
      write(0,*) 'Sub-cells to target:'
      write(0,*) 'H: h2tot=',h2tot,'h1tot=',h1tot,'dh=',h1tot-h2tot,'h2err=',h2err,'h1err=',h1err
      if (abs(h1tot-h2tot)>h2err+h1err) write(0,*) 'H non-conservation difference=',h1tot-h2tot,'allowed err=',h2err+h1err,' <-----!'
      write(0,*) 'UH: u2tot=',u2tot,'u1tot=',u1tot,'duh=',u1tot-u2tot,'u2err=',u2err,'u1err=',u1err,'uh_err=',uh_err
      if (abs(u1tot-u2tot)>u2err+u1err) write(0,*) 'U non-conservation difference=',u1tot-u2tot,'allowed err=',u2err+u1err,' <-----!'
      write(0,*) 'Source to target:'
      write(0,*) 'H: h0tot=',h0tot,'h1tot=',h1tot,'dh=',h1tot-h0tot,'h0err=',h0err,'h1err=',h1err
      if (abs(h1tot-h0tot)>h0err+h1err) write(0,*) 'H non-conservation difference=',h1tot-h0tot,'allowed err=',h0err+h1err,' <-----!'
      write(0,*) 'UH: u0tot=',u0tot,'u1tot=',u1tot,'duh=',u1tot-u0tot,'u0err=',u0err,'u1err=',u1err,'uh_err=',uh_err
      if (abs(u1tot-u0tot)>(u0err+u1err)+uh_err) write(0,*) 'U non-conservation difference=',u1tot-u0tot,'allowed err=',u0err+u1err+uh_err,' <-----!'
      write(0,*) 'U: u0min=',u0min,'u1min=',u1min,'u2min=',u2min
      if (u1min<u0min) write(0,*) 'U minimum overshoot=',u1min-u0min,' <-----!'
      if (u2min<u0min) write(0,*) 'U2 minimum overshoot=',u2min-u0min,' <-----!'
      write(0,*) 'U: u0max=',u0max,'u1max=',u1max,'u2max=',u2max
      if (u1max>u0max) write(0,*) 'U maximum overshoot=',u1max-u0max,' <-----!'
      if (u2max>u0max) write(0,*) 'U2 maximum overshoot=',u2max-u0max,' <-----!'
      write(0,'(a3,6a24,2a3)') 'k','h0','left edge','u0','right edge','h1','u1','is','ie'
      do k = 1, max(n0,n1)
        if (k<=min(n0,n1)) then
          write(0,'(i3,1p6e24.16,2i3)') k,h0(k),ppoly0_E(k,1),u0(k),ppoly0_E(k,2),h1(k),u1(k),itgt_start(k),itgt_end(k)
        elseif (k>n0) then
          write(0,'(i3,96x,1p2e24.16,2i3)') k,h1(k),u1(k),itgt_start(k),itgt_end(k)
        else
          write(0,'(i3,1p4e24.16)') k,h0(k),ppoly0_E(k,1),u0(k),ppoly0_E(k,2)
        endif
      enddo
      write(0,'(a3,2a24)') 'k','u0','Polynomial coefficients'
      do k = 1, n0
        write(0,'(i3,1p6e24.16)') k,u0(k),ppoly0_coefficients(k,:)
      enddo
      write(0,'(a3,3a24,a3,2a24)') 'k','Sub-cell h','Sub-cell u','Sub-cell hu','i0','xa','xb'
      xa = 0.
      dh0_eff = 0.
      do k = 1, n0+n1+1
        dh = h_sub(k)
        i0 = isub_src(k)
        dh0_eff = dh0_eff + dh ! Cumulative thickness within the source cell
        xb = dh0_eff / h0_eff(i0) ! This expression yields xa <= xb <= 1.0
        xb = min(1., xb) ! This is only needed when the total target column is wider than the source column
        write(0,'(i3,1p3e24.16,i3,1p2e24.16)') k,h_sub(k),u_sub(k),uh_sub(k),i0,xa,xb
        if (k<=n0+n1) then
          if (isub_src(k+1) /= i0) then
            dh0_eff = 0.; xa = 0.
          else
            xa = xb
          endif
        endif
      enddo
      call MOM_error( FATAL, 'MOM_remapping, remap_via_sub_cells: '//&
             'Remapping result is inconsistent!' )
    endif
    endif ! method<5
  endif ! debug_bounds

  ! Include the error remapping from source to sub-cells in the estimate of total remapping error
  uh_err = uh_err + u02_err

  if (present(ah_sub)) ah_sub(1:n0+n1+1) = h_sub(1:n0+n1+1)
  if (present(aisub_src)) aisub_src(1:n0+n1+1) = isub_src(1:n0+n1+1)
  if (present(aiss)) aiss(1:n0) = isrc_start(1:n0)
  if (present(aise)) aise(1:n0) = isrc_end(1:n0)

end subroutine remap_via_sub_cells

!> Returns the average value of a reconstruction within a single source cell, i0,
!! between the non-dimensional positions xa and xb (xa<=xb) with dimensional
!! separation dh.
real function average_value_ppoly( n0, u0, ppoly0_E, ppoly0_coefficients, method, i0, xa, xb)
  integer,       intent(in)    :: n0     !< Number of cells in source grid
  real,          intent(in)    :: u0(:)  !< Cell means
  real,          intent(in)    :: ppoly0_E(:,:)            !< Edge value of polynomial
  real,          intent(in)    :: ppoly0_coefficients(:,:) !< Coefficients of polynomial
  integer,       intent(in)    :: method !< Remapping scheme to use
  integer,       intent(in)    :: i0     !< Source cell index
  real,          intent(in)    :: xa     !< Non-dimensional start position within source cell
  real,          intent(in)    :: xb     !< Non-dimensional end position within source cell
  ! Local variables
  real :: u_ave, xa_2, xb_2, xa2pxb2, xapxb
  real, parameter :: r_3 = 1.0/3.0 ! Used in evaluation of integrated polynomials

  real :: mx, a_L, a_R, u_c, Ya, Yb, my, xa2b2ab, Ya2b2ab, a_c

  if (xb > xa) then
    select case ( method )
      case ( INTEGRATION_PCM )
        u_ave = u0(i0)
      case ( INTEGRATION_PLM )
        u_ave = (                                           &
            ppoly0_coefficients(i0,1)                       &
          + ppoly0_coefficients(i0,2) * 0.5 * ( xb + xa ) )
      case ( INTEGRATION_PPM )
        mx = 0.5 * ( xa + xb )
        a_L = ppoly0_E(i0, 1)
        a_R = ppoly0_E(i0, 2)
        u_c = u0(i0)
        a_c = 0.5 * ( ( u_c - a_L ) + ( u_c - a_R ) ) ! a_6 / 6
        if (mx<0.5) then
          ! This integration of the PPM reconstruction is expressed in distances from the left edge
          xa2b2ab = (xa*xa+xb*xb)+xa*xb
          u_ave = a_L + ( ( a_R - a_L ) * mx &
                          + a_c * ( 3. * ( xb + xa ) - 2.*xa2b2ab ) )
        else
          ! This integration of the PPM reconstruction is expressed in distances from the right edge
          Ya = 1. - xa
          Yb = 1. - xb
          my = 0.5 * ( Ya + Yb )
          Ya2b2ab = (Ya*Ya+Yb*Yb)+Ya*Yb
          u_ave = a_R  + ( ( a_L - a_R ) * my &
                           + a_c * ( 3. * ( Yb + Ya ) - 2.*Ya2b2ab ) )
        endif
      case ( INTEGRATION_PQM )
        xa_2 = xa*xa
        xb_2 = xb*xb
        xa2pxb2 = xa_2 + xb_2
        xapxb = xa + xb
        u_ave = (                                                                               &
              ppoly0_coefficients(i0,1)                                                         &
          + ( ppoly0_coefficients(i0,2) * 0.5 * ( xapxb )                                       &
          + ( ppoly0_coefficients(i0,3) * r_3 * ( xa2pxb2 + xa*xb )                             &
          + ( ppoly0_coefficients(i0,4) * 0.25* ( xa2pxb2 * xapxb )                             &
          +   ppoly0_coefficients(i0,5) * 0.2 * ( ( xb*xb_2 + xa*xa_2 ) * xapxb + xa_2*xb_2 ) ) ) ) )
      case default
        call MOM_error( FATAL,'The selected integration method is invalid' )
    end select
  else ! dh == 0.
    select case ( method )
      case ( INTEGRATION_PCM )
        u_ave =        ppoly0_coefficients(i0,1)
      case ( INTEGRATION_PLM )
       !u_ave =        ppoly0_coefficients(i0,1)   &
       !      + xa *   ppoly0_coefficients(i0,2)
        a_L = ppoly0_E(i0, 1)
        a_R = ppoly0_E(i0, 2)
        Ya = 1. - xa
        if (xa < 0.5) then
          u_ave = a_L + xa * ( a_R - a_L )
        else
          u_ave = a_R + Ya * ( a_L - a_R )
        endif
      case ( INTEGRATION_PPM )
       !u_ave =        ppoly0_coefficients(i0,1)   &
       !      + xa * ( ppoly0_coefficients(i0,2)   &
       !      + xa *   ppoly0_coefficients(i0,3) )
        a_L = ppoly0_E(i0, 1)
        a_R = ppoly0_E(i0, 2)
        u_c = u0(i0)
        a_c = 3. * ( ( u_c - a_L ) + ( u_c - a_R ) ) ! a_6
        Ya = 1. - xa
        if (xa < 0.5) then
          u_ave = a_L + xa * ( ( a_R - a_L ) + a_c * Ya )
        else
          u_ave = a_R + Ya * ( ( a_L - a_R ) + a_c * xa )
        endif
      case ( INTEGRATION_PQM )
        u_ave =        ppoly0_coefficients(i0,1)   &
              + xa * ( ppoly0_coefficients(i0,2)   &
              + xa * ( ppoly0_coefficients(i0,3)   &
              + xa * ( ppoly0_coefficients(i0,4)   &
              + xa *   ppoly0_coefficients(i0,5) ) ) )
      case default
        call MOM_error( FATAL,'The selected integration method is invalid' )
    end select
  endif
  average_value_ppoly = u_ave

end function average_value_ppoly

!> Measure totals and bounds on source grid
subroutine measure_input_bounds( n0, h0, u0, ppoly_E, h0tot, h0err, u0tot, u0err, u0min, u0max )
  integer,               intent(in)  :: n0 !< Number of cells on source grid
  real, dimension(n0),   intent(in)  :: h0 !< Cell widths on source grid
  real, dimension(n0),   intent(in)  :: u0 !< Cell averages on source grid
  real, dimension(n0,2), intent(in)  :: ppoly_E !< Cell edge values on source grid
  real,                  intent(out) :: h0tot !< Sum of cell widths
  real,                  intent(out) :: h0err !< Magnitude of round-off error in h0tot
  real,                  intent(out) :: u0tot !< Sum of cell widths times values
  real,                  intent(out) :: u0err !< Magnitude of round-off error in u0tot
  real,                  intent(out) :: u0min !< Minimum value in reconstructions of u0
  real,                  intent(out) :: u0max !< Maximum value in reconstructions of u0
  ! Local variables
  integer :: k
  real :: eps

  eps = epsilon(h0(1))
  h0tot = h0(1)
  h0err = 0.
  u0tot = h0(1) * u0(1)
  u0err = 0.
  u0min = min( ppoly_E(1,1), ppoly_E(1,2) )
  u0max = max( ppoly_E(1,1), ppoly_E(1,2) )
  do k = 2, n0
    h0tot = h0tot + h0(k)
    h0err = h0err + eps * max(h0tot, h0(k))
    u0tot = u0tot + h0(k) * u0(k)
    u0err = u0err + eps * max(abs(u0tot), abs(h0(k) * u0(k)))
    u0min = min( u0min, ppoly_E(k,1), ppoly_E(k,2) )
    u0max = max( u0max, ppoly_E(k,1), ppoly_E(k,2) )
  enddo

end subroutine measure_input_bounds

!> Measure totals and bounds on destination grid
subroutine measure_output_bounds( n1, h1, u1, h1tot, h1err, u1tot, u1err, u1min, u1max )
  integer,               intent(in)  :: n1 !< Number of cells on destination grid
  real, dimension(n1),   intent(in)  :: h1 !< Cell widths on destination grid
  real, dimension(n1),   intent(in)  :: u1 !< Cell averages on destination grid
  real,                  intent(out) :: h1tot !< Sum of cell widths
  real,                  intent(out) :: h1err !< Magnitude of round-off error in h1tot
  real,                  intent(out) :: u1tot !< Sum of cell widths times values
  real,                  intent(out) :: u1err !< Magnitude of round-off error in u1tot
  real,                  intent(out) :: u1min !< Minimum value in reconstructions of u1
  real,                  intent(out) :: u1max !< Maximum value in reconstructions of u1
  ! Local variables
  integer :: k
  real :: eps

  eps = epsilon(h1(1))
  h1tot = h1(1)
  h1err = 0.
  u1tot = h1(1) * u1(1)
  u1err = 0.
  u1min = u1(1)
  u1max = u1(1)
  do k = 2, n1
    h1tot = h1tot + h1(k)
    h1err = h1err + eps * max(h1tot, h1(k))
    u1tot = u1tot + h1(k) * u1(k)
    u1err = u1err + eps * max(abs(u1tot), abs(h1(k) * u1(k)))
    u1min = min(u1min, u1(k))
    u1max = max(u1max, u1(k))
  enddo

end subroutine measure_output_bounds
end module MOM_remapping
