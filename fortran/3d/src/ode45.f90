!-----------------------------------------------------------------------
! Copyright 2017 Nanyang Technological University, Singapore
!
! This file is part of UNICYCLE
!
! UNICYCLE is free software for non commercial usage: 
! you can redistribute it and/or modify it under the terms of the 
! Creative Commons CC BY-NC-SA 4.0 licence, please see
! https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
!
! UNICYCLE is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
! All intellectual property rights and licences for commercial usage 
! are reserved by Nanyang Technological University, please contact
! NTUItive if you are interested in a commericial licence for UNICYCLE.
! https://www.ntuitive.sg/
!
! If you use this code, please cite as James D P Moore, Sylvain Barbot, 
! Lujia Feng, Yu Hang, Valere Lambert, Eric Lindsey, Sagar Masuti,
! Takanori Matsuzawa, Jun Muto, Priyamvada Nanjundiah, Rino Salman,
! Sharadha Sathiakumar, and Harpreet Sethi. (2019, September 25). 
! jdpmoore/unicycle: Unicycle (Version 1.0). Zenodo. 
! http://doi.org/10.5281/zenodo.4471162
!
!-----------------------------------------------------------------------

#include "macros.f90"

MODULE ode45

  IMPLICIT NONE

  PUBLIC
  REAL(8), ALLOCATABLE :: AK2(:),AK3(:),AK4(:),AK5(:),AK6(:)
  REAL(8), ALLOCATABLE :: yrkck(:)

  ! numerical accuracy
  REAL(8) :: epsilon = 1.d-6

  ! maximum time step
  REAL*8 :: maximumTimeStep = 1.7976931348623158d308

CONTAINS

  !------------------------------------------------------------------------------
  !       FIFTH-ORDER RUNGE-KUTTA  : STEPPER ROUTINE
  !       SEE NUMERICAL RECIPES 2ND. ED. P.712
  !------------------------------------------------------------------------------
  SUBROUTINE RKQSm(n,t,y,dydt,yscal,ycand,yerr,ytmp,htry,hdid,hnext,odefun)

    IMPLICIT NONE

    INTEGER, INTENT(IN)   :: n
    REAL*8, INTENT(IN)    :: htry,yscal(n)
    REAL*8, INTENT(INOUT) :: ytmp(n)
    REAL*8, INTENT(INOUT) :: t
    REAL*8, INTENT(INOUT) :: dydt(n),y(n)
    REAL*8, INTENT(INOUT) :: ycand(n),yerr(n)
    REAL*8, INTENT(INOUT) :: hdid,hnext
    EXTERNAL              :: odefun

    INCLUDE 'mpif.h'

    REAL*8, PARAMETER    :: SAFETY=0.7d0,PGROW=-0.2d0,PSHRNK=-0.15d0,ERRCON=1.89d-4
    INTEGER              :: ierr
    REAL*8               :: h,errmax,errmax_thrd,tnew

    h=htry

    DO WHILE (.TRUE.)
      
       CALL RKCKm(n,t,y,dydt,h,ycand,yerr,ytmp,odefun)
       errmax_thrd=MAXVAL(ABS(yerr(1:n)/yscal(1:n)))
       errmax = errmax_thrd

       CALL MPI_ALLREDUCE(errmax_thrd,errmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)

       errmax = errmax/epsilon
       IF (errmax .GT. 1.d0) THEN
          h=h*SAFETY*(errmax**PSHRNK)
          IF (h .LT. 0.5d0*h) THEN
             h=0.5d0*h
          ENDIF
          tnew=t+h
          IF (tnew .EQ. t) THEN
             WRITE(STDERR,'("# stepsize underflow in RKQS")')
             STOP 1
          END IF
       ELSE
          IF (errmax .GT. ERRCON) THEN
             hnext=SAFETY*h*(errmax**PGROW)
          ELSE
             hnext=1.d1*h
          END IF
          IF (hnext .GT. maximumTimeStep) THEN
             hnext=maximumTimeStep
          END IF

          hdid=h
          t=t+h
          y(1:n)=ycand(1:n)
          EXIT 
       ENDIF
    END DO

  END SUBROUTINE RKQSm

  !------------------------------------------------------------------------------
  !       FIFTH-ORDER RUNGE-KUTTA  : ALGORITHM ROUTINE
  !       SEE NUMERICAL RECIPES 2ND. ED. P.713
  !
  ! uses the BLAS library
  !------------------------------------------------------------------------------
  SUBROUTINE RKCKm(n,t,yin,dydt,h,yout,yerr,ytmp,odefun)

    IMPLICIT NONE

    INTEGER, INTENT(IN)   :: n
    REAL*8, INTENT(IN)    :: t,h
    REAL*8, INTENT(INOUT) :: yin(n),dydt(n), &
                             yout(n),yerr(n),ytmp(n)
    EXTERNAL :: odefun

    REAL*8, PARAMETER :: A2=0.2d0, A3=0.3d0, A4=0.6d0, A5=1.d0, A6=0.875d0
    REAL*8, PARAMETER :: B21=0.2d0, B31=3.d0/40.d0, B32=9.d0/40.d0
    REAL*8, PARAMETER :: B41=0.3d0, B42=-0.9d0, B43=1.2d0
    REAL*8, PARAMETER :: B51=-11.d0/54.d0, B52=2.5d0, B53=-70.d0/27.d0, B54=35.d0/27.d0  
    REAL*8, PARAMETER :: B61=1631.d0/55296.d0,   B62=175.d0/512.d0, B63=575.d0/13824.d0, &
                         B64=44275.d0/110592.d0, B65=253.d0/4096.d0
    REAL*8, PARAMETER :: C1=37.d0/378.d0, C3=250.d0/621.d0, C4=125.d0/594.d0, C6=512.d0/1771.d0
    REAL*8, PARAMETER :: DC1=C1-2825.d0/27648.d0,  DC3=C3-18575.d0/48384.d0,        &
                         DC4=C4-13525.d0/55296.d0, DC5=-277.d0/14336.d0, DC6=C6-0.25d0

    !ytmp = yin + h*B21*dydt
    CALL dcopy(n,yin,1,ytmp,1)
    CALL daxpy(n,h*B21,dydt,1,ytmp,1)

    CALL odefun(n,t+A2*h,ytmp,AK2)

    !ytmp = yin + h*(B31*dydt + B32*AK2)
    CALL dcopy(n,dydt,1,    yrkck,1)
    CALL dscal(n,B31,       yrkck,1)
    CALL daxpy(n,B32,AK2 ,1,yrkck,1)
    CALL dcopy(n,yin,1,ytmp,1)
    CALL daxpy(n,h,yrkck,1,ytmp,1)

    CALL odefun(n,t+A3*h,ytmp,AK3)

    !ytmp = yin + h*(B41*dydt + B42*AK2 + B43*AK3)
    CALL dcopy(n,dydt,1,    yrkck,1)
    CALL dscal(n,B41,       yrkck,1)
    CALL daxpy(n,B42,AK2 ,1,yrkck,1)
    CALL daxpy(n,B43,AK3 ,1,yrkck,1)
    CALL dcopy(n,yin,1,ytmp,1)
    CALL daxpy(n,h,yrkck,1,ytmp,1)

    CALL odefun(n,t+A4*h,ytmp,AK4)

    !ytmp = yin + h*(B51*dydt + B52*AK2 + B53*AK3 + B54*AK4)
    CALL dcopy(n,dydt,1,    yrkck,1)
    CALL dscal(n,B51,       yrkck,1)
    CALL daxpy(n,B52,AK2 ,1,yrkck,1)
    CALL daxpy(n,B53,AK3 ,1,yrkck,1)
    CALL daxpy(n,B54,AK4 ,1,yrkck,1)
    CALL dcopy(n,yin,1,ytmp,1)
    CALL daxpy(n,h,yrkck,1,ytmp,1)

    CALL odefun(n,t+A5*h,ytmp,AK5)

    !ytmp = yin + h*(B61*dydt + B62*AK2 + B63*AK3 + B64*AK4 + B65*AK5)
    CALL dcopy(n,dydt,1,    yrkck,1)
    CALL dscal(n,B61,       yrkck,1)
    CALL daxpy(n,B62,AK2 ,1,yrkck,1)
    CALL daxpy(n,B63,AK3 ,1,yrkck,1)
    CALL daxpy(n,B64,AK4 ,1,yrkck,1)
    CALL daxpy(n,B65,AK5 ,1,yrkck,1)
    CALL dcopy(n,yin,1,ytmp,1)
    CALL daxpy(n,h,yrkck,1,ytmp,1)

    CALL odefun(n,t+A6*h,ytmp,AK6)

    !yout = yin + h*(C1*dydt + C3*AK3 + C4*AK4 + C6*AK6)
    CALL dcopy(n,dydt,1,   yrkck,1)
    CALL dscal(n,C1,       yrkck,1)
    CALL daxpy(n,C3,AK3 ,1,yrkck,1)
    CALL daxpy(n,C4,AK4 ,1,yrkck,1)
    CALL daxpy(n,C6,AK6 ,1,yrkck,1)
    CALL dcopy(n,yin,1,yout,1)
    CALL daxpy(n,h,yrkck,1,yout,1)

    !yerr = h*(DC1*dydt + DC3*AK3 + DC4*AK4 + DC5*AK5 + DC6*AK6)
    CALL dcopy(n,dydt,1,    yerr,1)
    CALL dscal(n,DC1,       yerr,1)
    CALL daxpy(n,DC3,AK3 ,1,yerr,1)
    CALL daxpy(n,DC4,AK4 ,1,yerr,1)
    CALL daxpy(n,DC5,AK5 ,1,yerr,1)
    CALL daxpy(n,DC6,AK6 ,1,yerr,1)
    CALL dscal(n,h,yerr,1)

  END SUBROUTINE RKCKm

END MODULE ode45

