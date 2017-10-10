module MOM_EOS_linear

! This file is part of MOM6. See LICENSE.md for the license.

!***********************************************************************
!*  The subroutines in this file implement a simple linear equation of *
!*  state for sea water with constant coefficients set as parameters.  *
!***********************************************************************

implicit none ; private

public calculate_density_derivs_linear

contains
!> This subroutine calculates the partial derivatives of density    *
!! with potential temperature and salinity.
subroutine calculate_density_derivs_linear(T, S, pressure, drho_dT_out, &
                       drho_dS_out, Rho_T0_S0, dRho_dT, dRho_dS, start, npts)
  real,    intent(in),  dimension(:) :: T           !< Potential temperature relative to the surface
                                                    !! in C.
  real,    intent(in),  dimension(:) :: S           !< Salinity in PSU.
  real,    intent(in),  dimension(:) :: pressure    !< Pressure in Pa.
  real,    intent(out), dimension(:) :: drho_dT_out !< The partial derivative of density with
                                                    !! potential temperature, in kg m-3 K-1.
  real,    intent(out), dimension(:) :: drho_dS_out !< The partial derivative of density with
                                                    !! salinity, in kg m-3 psu-1.
  real,    intent(in)                :: Rho_T0_S0   !< The density at T=0, S=0, in kg m-3.
  real,    intent(in)                :: dRho_dT, dRho_dS !< The derivatives of density with
                                                    !! temperature and salinity, in kg m-3 C-1
                                                    !! and kg m-3 psu-1.
  integer, intent(in)                :: start       !< The starting point in the arrays.
  integer, intent(in)                :: npts        !< The number of values to calculate.

! *   This subroutine calculates the partial derivatives of density    *
! * with potential temperature and salinity.                           *
! *                                                                    *
! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in PSU.                                    *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     drho_dT_out - the partial derivative of density with    *
! *                      potential temperature, in kg m-3 K-1.         *
! *  (out)     drho_dS_out - the partial derivative of density with    *
! *                      salinity, in kg m-3 psu-1.                    *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
! *  (in)      Rho_T0_S0 - The density at T=0, S=0, in kg m-3.         *
! *  (in)      dRho_dT - The derivatives of density with temperature   *
! *  (in)      dRho_dS - and salinity, in kg m-3 C-1 and kg m-3 psu-1. *
  integer :: j

  do j=start,start+npts-1
    drho_dT_out(j) = dRho_dT
    drho_dS_out(j) = dRho_dS
  enddo

end subroutine calculate_density_derivs_linear
end module MOM_EOS_linear
