module MOM_EOS_Wright

! This file is part of MOM6. See LICENSE.md for the license.

!***********************************************************************
!*  The subroutines in this file implement the equation of state for   *
!*  sea water using the formulae given by  Wright, 1997, J. Atmos.     *
!*  Ocean. Tech., 14, 735-740.  Coded by R. Hallberg, 7/00.            *
!***********************************************************************

implicit none ; private

public calculate_density_derivs_wright

!real :: a0, a1, a2, b0, b1, b2, b3, b4, b5, c0, c1, c2, c3, c4, c5
!    One of the two following blocks of values should be commented out.
!  Following are the values for the full range formula.
!
!real, parameter :: a0 = 7.133718e-4, a1 = 2.724670e-7, a2 = -1.646582e-7
!real, parameter :: b0 = 5.613770e8,  b1 = 3.600337e6,  b2 = -3.727194e4
!real, parameter :: b3 = 1.660557e2,  b4 = 6.844158e5,  b5 = -8.389457e3
!real, parameter :: c0 = 1.609893e5,  c1 = 8.427815e2,  c2 = -6.931554
!real, parameter :: c3 = 3.869318e-2, c4 = -1.664201e2, c5 = -2.765195


! Following are the values for the reduced range formula.
real, parameter :: a0 = 7.057924e-4, a1 = 3.480336e-7, a2 = -1.112733e-7
real, parameter :: b0 = 5.790749e8,  b1 = 3.516535e6,  b2 = -4.002714e4
real, parameter :: b3 = 2.084372e2,  b4 = 5.944068e5,  b5 = -9.643486e3
real, parameter :: c0 = 1.704853e5,  c1 = 7.904722e2,  c2 = -7.984422
real, parameter :: c3 = 5.140652e-2, c4 = -2.302158e2, c5 = -3.079464

contains

!> For a given thermodynamic state, return the thermal/haline expansion coefficients
subroutine calculate_density_derivs_wright(T, S, pressure, drho_dT, drho_dS, start, npts)
  real,    intent(in),  dimension(:) :: T        !< Potential temperature relative to the surface
                                                 !! in C.
  real,    intent(in),  dimension(:) :: S        !< Salinity in PSU.
  real,    intent(in),  dimension(:) :: pressure !< Pressure in Pa.
  real,    intent(out), dimension(:) :: drho_dT  !< The partial derivative of density with potential
                                                 !! temperature, in kg m-3 K-1.
  real,    intent(out), dimension(:) :: drho_dS  !< The partial derivative of density with salinity,
                                                 !! in kg m-3 psu-1.
  integer, intent(in)                :: start    !< The starting point in the arrays.
  integer, intent(in)                :: npts     !< The number of values to calculate.

! * Arguments: T - potential temperature relative to the surface in C. *
! *  (in)      S - salinity in PSU.                                    *
! *  (in)      pressure - pressure in Pa.                              *
! *  (out)     drho_dT - the partial derivative of density with        *
! *                      potential temperature, in kg m-3 K-1.         *
! *  (out)     drho_dS - the partial derivative of density with        *
! *                      salinity, in kg m-3 psu-1.                    *
! *  (in)      start - the starting point in the arrays.               *
! *  (in)      npts - the number of values to calculate.               *
  real :: al0, p0, lambda, I_denom2
  integer :: j

  do j=start,start+npts-1
    al0 = (a0 + a1*T(j)) + a2*S(j)
    p0 = (b0 + b4*S(j)) + T(j) * (b1 + T(j)*((b2 + b3*T(j))) + b5*S(j))
    lambda = (c0 +c4*S(j)) + T(j) * (c1 + T(j)*((c2 + c3*T(j))) + c5*S(j))

    I_denom2 = 1.0 / (lambda + al0*(pressure(j) + p0))
    I_denom2 = I_denom2 *I_denom2
    drho_dT(j) = I_denom2 * &
      (lambda* (b1 + T(j)*(2.0*b2 + 3.0*b3*T(j)) + b5*S(j)) - &
       (pressure(j)+p0) * ( (pressure(j)+p0)*a1 + &
        (c1 + T(j)*(c2*2.0 + c3*3.0*T(j)) + c5*S(j)) ))
    drho_dS(j) = I_denom2 * (lambda* (b4 + b5*T(j)) - &
      (pressure(j)+p0) * ( (pressure(j)+p0)*a2 + (c4 + c5*T(j)) ))
  enddo

end subroutine calculate_density_derivs_wright
end module MOM_EOS_Wright
