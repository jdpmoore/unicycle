
!------------------------------------------------------------------------
!> function xLogY
!! computes x*log(y) and enforces 0*log(0)=0 to avoid NaN
!! \author James Moore, June 2016
!------------------------------------------------------------------------
REAL*8 FUNCTION xLogy(x,y)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: x,y

  IF (0 .EQ. x) THEN
     xLogy=0.000000000000000000000000000_8
  ELSE
     xLogy=x*LOG(y)
  END IF

END FUNCTION xLogy

