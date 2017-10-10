!> Provides subroutines for quantities specific to the equation of state
module MOM_EOS

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_EOS_linear, only : calculate_density_derivs_linear
use MOM_EOS_Wright, only : calculate_density_derivs_wright

implicit none ; private

public calculate_density_derivs

! The named integers that might be stored in eqn_of_state_type%form_of_EOS.
integer, parameter, public :: EOS_LINEAR = 1
integer, parameter, public :: EOS_UNESCO = 2
integer, parameter, public :: EOS_WRIGHT = 3
integer, parameter, public :: EOS_TEOS10 = 4
integer, parameter, public :: EOS_NEMO   = 5

contains

!> Calls the appropriate subroutine to calculate density derivatives for 1-D array inputs.
subroutine calculate_density_derivs(T, S, pressure, drho_dT, drho_dS, start, npts, EOS, linear_Rho_T0_S0, linear_dRho_dT, linear_dRho_dS)
  real, dimension(:), intent(in)  :: T !< Potential temperature referenced to the surface (degC)
  real, dimension(:), intent(in)  :: S !< Salinity (PSU)
  real, dimension(:), intent(in)  :: pressure !< Pressure (Pa)
  real, dimension(:), intent(out) :: drho_dT !< The partial derivative of density with potential tempetature, in kg m-3 K-1.
  real, dimension(:), intent(out) :: drho_dS !< The partial derivative of density with salinity, in kg m-3 psu-1.
  integer,            intent(in)  :: start !< Starting index within the array
  integer,            intent(in)  :: npts !< The number of values to calculate
  integer, intent(in)     :: EOS !< Equation of state structure
  real, optional, intent(in) :: linear_Rho_T0_S0, linear_dRho_dT, linear_dRho_dS
  !!
  select case (EOS)
    case (EOS_LINEAR)
      call calculate_density_derivs_linear(T, S, pressure, drho_dT, drho_dS, linear_Rho_T0_S0, linear_dRho_dT, linear_dRho_dS, &
                                           start, npts)
    case (EOS_WRIGHT)
      call calculate_density_derivs_wright(T, S, pressure, drho_dT, drho_dS, start, npts)
  end select

end subroutine calculate_density_derivs
end module MOM_EOS
