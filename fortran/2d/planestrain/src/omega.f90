
!------------------------------------------------------------------------
!> function Omega(x)
!! evaluates the boxcar function
!------------------------------------------------------------------------
REAL*8 FUNCTION omega(x)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: x

  REAL*8, EXTERNAL :: heaviside
  REAL*8, PARAMETER :: half = 0.500000000000000000000000000000000000_8

  omega=heaviside(x+half)-heaviside(x-half)

END FUNCTION omega

