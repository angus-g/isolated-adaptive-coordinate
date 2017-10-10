module remapping
  use MOM_remapping, only : remapping_core_h, remapping_CS
  use iso_c_binding, only : c_double, c_int

  implicit none
  private

  integer, parameter  :: REMAPPING_PCM        = 0 !< O(h^1) remapping scheme
  integer, parameter  :: REMAPPING_PLM        = 1 !< O(h^2) remapping scheme
  integer, parameter  :: REMAPPING_PPM_H4     = 2 !< O(h^3) remapping scheme
  integer, parameter  :: REMAPPING_PPM_IH4    = 3 !< O(h^3) remapping scheme
  integer, parameter  :: REMAPPING_PQM_IH4IH3 = 4 !< O(h^4) remapping scheme
  integer, parameter  :: REMAPPING_PQM_IH6IH5 = 5 !< O(h^5) remapping scheme

contains

  subroutine remap_wrap(scheme, n0, h0, u0, n1, h1, u1) bind(c)
    integer(c_int), intent(in), value :: scheme, n0, n1
    real(c_double), dimension(n0), intent(in) :: h0, u0
    real(c_double), dimension(n1), intent(in) :: h1
    real(c_double), dimension(n1), intent(out) :: u1

    type(remapping_CS) :: CS
    CS%remapping_scheme = scheme
    select case (scheme)
    case (REMAPPING_PCM)
      CS%degree = 0
    case (REMAPPING_PLM)
      CS%degree = 1
    case (REMAPPING_PPM_H4, REMAPPING_PPM_IH4)
      CS%degree = 2
    case (REMAPPING_PQM_IH4IH3, REMAPPING_PQM_IH6IH5)
      CS%degree = 4
    end select

    call remapping_core_h(CS, n0, h0, u0, n1, h1, u1)
  end subroutine remap_wrap
end module remapping
