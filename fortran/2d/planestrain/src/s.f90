
!------------------------------------------------------------------------
!> function S(x)
!! evalutes the shifted boxcar function
!------------------------------------------------------------------------
REAL*8 FUNCTION s(x)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: x

  REAL*8, EXTERNAL :: omega

  s=omega(x-0.500000000000000000_8)

END FUNCTION s

