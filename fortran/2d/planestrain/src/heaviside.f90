
!------------------------------------------------------------------------
!> function heaviside
!! computes the Heaviside function
!------------------------------------------------------------------------
REAL*8 FUNCTION heaviside(x)
  REAL*8, INTENT(IN) :: x

   IF (0 .LT. x) THEN
      heaviside=1.00000000000000000000000_8
   ELSE
      heaviside=0.00000000000000000000000_8
   END IF

END FUNCTION heaviside

