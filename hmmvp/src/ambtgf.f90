! AMB's interface to tgf. Code copied and modified from elh3dtria_omp.f.

SUBROUTINE AMBTGFELD( &
  rlam, rmu, space, triangle, trinorm, obs, isself, slip, u, strain)
  ! A triangle's influence on a point. Slip vector input. isself is 't' if the
  ! obs is the centroid of the triangle and 'f' if obs is outside the
  ! triangle.
  IMPLICIT NONE
  
  CHARACTER, INTENT(IN) :: space, isself
  REAL(8), INTENT(IN) :: rlam, rmu, slip(3), triangle(3,3), obs(3), &
       trinorm(3)
  REAL(8), INTENT(OUT) :: u(3), strain(3,3)
  REAL(8) :: u0(3), strain0(3,3)
  INTEGER :: numfunev

  IF (isself .EQ. 't') then
     CALL eltst3triadirectself(rlam, rmu, 1, 1, triangle, slip, trinorm, &
          obs, u, 1, strain)
  ELSE
     CALL eltst3triadirecttarg(rlam, rmu, 1, triangle, slip, trinorm, &
          obs, u, 1, strain)
  END IF
  IF (space .EQ. 'h') THEN
     ! Originally I was not going to take u0. But some of the loops in their
     ! code check values that aren't set if certain quantities aren't required,
     ! as valgrind pointed out to me. Perhaps it's a benign bug (fortran perhaps
     ! always inits all memory to 0), but I don't want to worry about it.
     CALL elth3triaadap(rlam, rmu, 1, triangle, slip, trinorm, obs, 1, u0, &
          1, strain0, numfunev)
     u = u + u0
     strain = strain + strain0
  END IF
END SUBROUTINE ambtgfeld
