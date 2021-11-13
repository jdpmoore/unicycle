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
! http://doi.org/10.5281/zenodo.5688288
!
!-----------------------------------------------------------------------

#include "macros.f90"

MODULE stuart97

  IMPLICIT NONE

CONTAINS

  !------------------------------------------------------------------------
  !> subroutine computeReferenceSystemStuart97
  !! computes the center position and local reference system of triangular
  !! dislocations.
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeReferenceSystemStuart97(nVe,v,ns,i1,i2,i3,sv,dv,nv,xc)
    INTEGER, INTENT(IN) :: nVe,ns
    REAL*8, DIMENSION(3,nVe), INTENT(IN) :: v
    INTEGER, DIMENSION(ns), INTENT(IN) :: i1,i2,i3
    REAL*8, DIMENSION(3,ns), INTENT(OUT) :: sv,dv,nv,xc
    
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: A,B,C
    REAL*8 :: norm
    INTEGER :: j,iostatus

    ALLOCATE(A(3,ns),B(3,ns),C(3,ns),STAT=iostatus)
    IF (0 /= iostatus) STOP "could not allocation the vertices"

    ! vertices
    DO j=1,ns
       A(:,j)=v(:,i1(j))
       B(:,j)=v(:,i2(j))
       C(:,j)=v(:,i3(j))
    END DO

    ! normal vector (B-A) x (C-A)
    DO j=1,ns
       nv(1,j)=(B(2,j)-A(2,j))*(C(3,j)-A(3,j))-(B(3,j)-A(3,j))*(C(2,j)-A(2,j))
       nv(2,j)=(B(3,j)-A(3,j))*(C(1,j)-A(1,j))-(B(1,j)-A(1,j))*(C(3,j)-A(3,j))
       nv(3,j)=(B(1,j)-A(1,j))*(C(2,j)-A(2,j))-(B(2,j)-A(2,j))*(C(1,j)-A(1,j))

       norm=SQRT(nv(1,j)**2+nv(2,j)**2+nv(3,j)**2)
       nv(:,j)=nv(:,j)/norm

       ! choose upward-pointing normal vectors
       IF (0 .LT. nv(3,j)) THEN
          nv(:,j)=-nv(:,j)
       END IF

       ! strike-direction vector
       sv(1,j)= COS(ATAN2(nv(1,j),nv(2,j)))
       sv(2,j)=-SIN(ATAN2(nv(1,j),nv(2,j)))
       sv(3,j)=0._8
                 
       ! dip-direction vector dv = nv x sv
       dv(1,j)=nv(2,j)*sv(3,j)-nv(3,j)*sv(2,j)
       dv(2,j)=nv(3,j)*sv(1,j)-nv(1,j)*sv(3,j)
       dv(3,j)=nv(1,j)*sv(2,j)-nv(2,j)*sv(1,j)

    END DO

    ! center of fault triangular patch
    xc=(A+B+C)/3._8

    DEALLOCATE(A,B,C)

  END SUBROUTINE computeReferenceSystemStuart97

  !------------------------------------------------------------------------
  !> subroutine computeTractionStuart97
  !! calculates the stress associated with a triangular dislocation in an 
  !! elastic half-space using the analytic solution of Comninou and Dundurs
  !! (1975) and the numerical implementation of Stuart (1997).
  !!
  !! INPUT:
  !! @param x1,x2,x3 - coordinates of observation points
  !! @param sv,dv,nv - strike, dip, and normal vectors
  !! @param q1,q2,q2 - coordinates of the triangular dislocation vertices
  !! @param s,d,t    - slip vector components (strike slip, dip slip 
  !!                    and tensile opening)
  !! @param G,lambda - elastic moduli (Lame parameters)
  !!
  !! OUTPUT:
  !! s11,s12,s13,s22,s23,s33 - the stress components
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeTractionStuart97(x1,x2,x3,sv,dv,nv, &
                                    q1,q2,q3, &
                                    s,d,t,G,lambda, &
                                    s11,s12,s13,s22,s23,s33)
    REAL*8, INTENT(IN) :: x1,x2,x3
    REAL*8, DIMENSION(3) :: sv,dv,nv
    REAL*8, INTENT(IN), DIMENSION(3) :: q1,q2,q3
    REAL*8, INTENT(IN) :: s,d,t
    REAL*8, INTENT(IN) :: G,lambda
    REAL*8, INTENT(OUT) :: s11,s12,s13,s22,s23,s33

    REAL*8 :: e11,e12,e13,e22,e23,e33

    REAL*8, DIMENSION(3) :: P1,P2,P3

    P1=(/ q1(1), q2(1), q3(1) /)
    P2=(/ q1(2), q2(2), q3(2) /)
    P3=(/ q1(3), q2(3), q3(3) /)

    ! call C function
    CALL StrainInHalfSpace( &
                     s11,s12,s13,s22,s23,s33,e11,e12,e13,e22,e23,e33, &
                     %VAL(x1),%VAL(x2),%VAL(x3), &
                     P1,P2,P3, &
                     %VAL(s),%VAL(d),%VAL(t), &
                     %VAL(G),%VAL(lambda))

    ! post-process as necessary
    !s11=
  END SUBROUTINE computeTractionStuart97

  !------------------------------------------------------------------------
  !> subroutine computeTractionStuart97
  !! calculates the traction associated with a triangular dislocation in an 
  !! elastic half-space using the analytic solution of Comninou and Dundurs
  !! (1975) and the numerical implementation of Stuart (1997).
  !!
  !! INPUT:
  !! @param x1,x2,x3 - coordinates of observation points
  !! @param q1,q2,q2 - coordinates of the triangular dislocation vertices
  !! @param s,d,t    - slip vector components (strike slip, dip slip 
  !!                    and tensile opening)
  !! @param G,lambda - elastic moduli (Lame parameters)
  !!
  !! OUTPUT:
  !! s11,s12,s13,s22,s23,s33 - the stress components
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeStressStuart97(x1,x2,x3, &
                                    P1,P2,P3, &
                                    ss,d,t,G,lambda, &
                                    s11,s12,s13,s22,s23,s33)
    REAL*8, INTENT(IN) :: x1,x2,x3
    REAL*8, DIMENSION(3), INTENT(IN) :: P1,P2,P3
    REAL*8, INTENT(IN) :: ss,d,t
    REAL*8, INTENT(IN) :: G,lambda
    REAL*8, INTENT(OUT) :: s11,s12,s13,s22,s23,s33

    REAL*8, DIMENSION(3) :: xo,u
    REAL*8, DIMENSION(9) :: tridlc
    REAL*8, DIMENSION(3,3) :: tilt
    REAL*8 :: ds,ekk,nu

    nu=lambda/(lambda+G)/2._8

    ! observation point
    xo=(/ x1,x2,x3 /)

    ! trianglar corner coordinates
    tridlc= (/ P1(1), P1(2), P1(3), &
               P2(1), P2(2), P2(3), &
               P3(1), P3(2), P3(3) /)

    ! normal displacement
    ds=-d

    CALL dstuart(nu,xo,tridlc,ss,ds,t,u,tilt)

    ekk=tilt(1,1)+tilt(2,2)+tilt(3,3)

    ! change from east,north,up to north,east,down
    s11=lambda*ekk+2._8*G*tilt(1,1)
    s12=G*(tilt(1,2)+tilt(2,1))
    s13=G*(tilt(1,3)+tilt(3,1))
    s22=lambda*ekk+2._8*G*tilt(2,2)
    s23=G*(tilt(2,3)+tilt(3,2))
    s33=lambda*ekk+2._8*G*tilt(3,3)

  END SUBROUTINE computeStressStuart97

  !------------------------------------------------------------------------
  !> subroutine computeTractionKernelsStuart97
  !! calculates the traction kernels associated with the motion of
  !! triangular dislocations in an elastic half-space using the analytic
  !! solution of Comninou and Dundurs (1975) and the numerical
  !! implementation of Stuart (1997).
  !!
  !! INPUT:
  !! @param x                - coordinates (NED) of observation point
  !! @param sv,dv,nv         - strike, dip and normal vector
  !! @param ns               - number of sources
  !! @param i1,i2,i3         - index of vertices forming triangles
  !! @param nVe              - number of vertices
  !! @param vertices         - coordinates of triangle vertices
  !! @param s,d,t            - strike slip, dip slip and opening
  !! @param G,lambda         - rigidity, Lame parameter
  !!
  !! OUTPUT:
  !! ts,td,tn                - traction in the strike, dip and normal
  !!                           directions.
  !------------------------------------------------------------------------
  SUBROUTINE computeTractionKernelsStuart97(x,sv,dv,nv, &
                                           ns,i1,i2,i3, &
                                           nVe,vertices, &
                                           s,d,t,lambda,G, &
                                           ts,td,tn)
    REAL*8, DIMENSION(3), INTENT(IN) :: x
    REAL*8, DIMENSION(3), INTENT(IN) :: sv,dv,nv
    INTEGER, INTENT(IN) :: ns
    INTEGER, DIMENSION(ns), INTENT(IN) :: i1,i2,i3
    INTEGER, INTENT(IN) :: nVe
    REAL*8, DIMENSION(3,nVe), INTENT(IN) :: vertices
    REAL*8, INTENT(IN) :: s,d,t,lambda,G
    REAL*8, DIMENSION(ns), INTENT(OUT) :: ts,td,tn

    INTEGER :: i
    REAL*8 :: s11,s12,s13,s22,s23,s33

    ! loop over sources
    DO i=1,ns
       CALL computeStressStuart97( &
                     x(1),x(2),x(3), &
                     vertices(:,i1(i)),vertices(:,i2(i)),vertices(:,i3(i)), &
                     s,d,t,G,lambda, &
                     s11,s12,s13,s22,s23,s33)

       ! rotate to receiver system of coordinates
       ts(i)= ( nv(1)*s11+nv(2)*s12+nv(3)*s13 )*sv(1) &
             +( nv(1)*s12+nv(2)*s22+nv(3)*s23 )*sv(2) &
             +( nv(1)*s13+nv(2)*s23+nv(3)*s33 )*sv(3)
       td(i)= ( nv(1)*s11+nv(2)*s12+nv(3)*s13 )*dv(1) &
             +( nv(1)*s12+nv(2)*s22+nv(3)*s23 )*dv(2) &
             +( nv(1)*s13+nv(2)*s23+nv(3)*s33 )*dv(3)
       tn(i)= ( nv(1)*s11+nv(2)*s12+nv(3)*s13 )*nv(1) &
             +( nv(1)*s12+nv(2)*s22+nv(3)*s23 )*nv(2) &
             +( nv(1)*s13+nv(2)*s23+nv(3)*s33 )*nv(3)

    END DO

  END SUBROUTINE computeTractionKernelsStuart97

  !------------------------------------------------------------------------
  !> subroutine computeStressKernelsStuart97
  !! calculates the stress kernels associated with the motion of triangular
  !! dislocations in an elastic half-space using the analytic solution of 
  !! Comninou and Dundurs (1975) and the numerical implementation of Stuart (1997).
  !!
  !! INPUT:
  !! @param x                - coordinates (NED) of observation point
  !! @param sv,dv,nv         - strike, dip and normal vector
  !! @param ns               - number of sources
  !! @param i1,i2,i3         - index of vertices forming triangles
  !! @param nVe              - number of vertices
  !! @param vertices         - coordinates of triangle vertices
  !! @param s,d,t            - strike slip, dip slip and opening
  !! @param G,lambda         - rigidity, Lame parameter
  !!
  !! OUTPUT:
  !! s11,s12,s13,s22,s23,s33 - stress components in the reference system
  !!                           tied to the receiver.
  !------------------------------------------------------------------------
  SUBROUTINE computeStressKernelsStuart97(x,sv,dv,nv, &
                                           ns,i1,i2,i3, &
                                           nVe,vertices, &
                                           s,d,t,lambda,G, &
                                           s11,s12,s13,s22,s23,s33)
    REAL*8, DIMENSION(3), INTENT(IN) :: x
    REAL*8, DIMENSION(3), INTENT(IN) :: sv,dv,nv
    INTEGER, INTENT(IN) :: ns
    INTEGER, DIMENSION(ns), INTENT(IN) :: i1,i2,i3
    INTEGER, INTENT(IN) :: nVe
    REAL*8, DIMENSION(3,nVe), INTENT(IN) :: vertices
    REAL*8, INTENT(IN) :: s,d,t,lambda,G
    REAL*8, DIMENSION(ns), INTENT(OUT) :: s11,s12,s13,s22,s23,s33

    INTEGER :: i
    REAL*8 :: s11p,s12p,s13p,s22p,s23p,s33p

    ! loop over sources
    DO i=1,ns
       CALL computeStressStuart97( &
                     x(1),x(2),x(3), &
                     vertices(:,i1(i)),vertices(:,i2(i)),vertices(:,i3(i)), &
                     s,d,t,G,lambda, &
                     s11(i),s12(i),s13(i),s22(i),s23(i),s33(i))

       ! rotate to receiver system of coordinates
       s11p= ( sv(1)*s11(i)+sv(2)*s12(i)+sv(3)*s13(i) )*sv(1) &
            +( sv(1)*s12(i)+sv(2)*s22(i)+sv(3)*s23(i) )*sv(2) &
            +( sv(1)*s13(i)+sv(2)*s23(i)+sv(3)*s33(i) )*sv(3)
       s12p= ( sv(1)*s11(i)+sv(2)*s12(i)+sv(3)*s13(i) )*nv(1) &
            +( sv(1)*s12(i)+sv(2)*s22(i)+sv(3)*s23(i) )*nv(2) &
            +( sv(1)*s13(i)+sv(2)*s23(i)+sv(3)*s33(i) )*nv(3)
       s13p=-( sv(1)*s11(i)+sv(2)*s12(i)+sv(3)*s13(i) )*dv(1) &
            -( sv(1)*s12(i)+sv(2)*s22(i)+sv(3)*s23(i) )*dv(2) &
            -( sv(1)*s13(i)+sv(2)*s23(i)+sv(3)*s33(i) )*dv(3)
       s22p= ( nv(1)*s11(i)+nv(2)*s12(i)+nv(3)*s13(i) )*nv(1) &
            +( nv(1)*s12(i)+nv(2)*s22(i)+nv(3)*s23(i) )*nv(2) &
            +( nv(1)*s13(i)+nv(2)*s23(i)+nv(3)*s33(i) )*nv(3)
       s23p=-( nv(1)*s11(i)+nv(2)*s12(i)+nv(3)*s13(i) )*dv(1) &
            -( nv(1)*s12(i)+nv(2)*s22(i)+nv(3)*s23(i) )*dv(2) &
            -( nv(1)*s13(i)+nv(2)*s23(i)+nv(3)*s33(i) )*dv(3)
       s33p=-(-dv(1)*s11(i)-dv(2)*s12(i)-dv(3)*s13(i) )*dv(1) &
            -(-dv(1)*s12(i)-dv(2)*s22(i)-dv(3)*s23(i) )*dv(2) &
            -(-dv(1)*s13(i)-dv(2)*s23(i)-dv(3)*s33(i) )*dv(3)

       s11(i)=s11p
       s12(i)=s12p
       s13(i)=s13p
       s22(i)=s22p
       s23(i)=s23p
       s33(i)=s33p

    END DO

  END SUBROUTINE computeStressKernelsStuart97

  !------------------------------------------------------------------------
  !> subroutine computeDisplacementStuart97
  !! calculates the displacements associated with a triangular dislocation
  !! in an elastic half-space using the analytic solution of Comninou and Dundurs
  !! (1975) and the numerical implementation of Stuart (1997).
  !!
  !! INPUT:
  !! @param x1,x2,x3 - coordinates of observation points
  !! @param q1,q2,q2 - coordinates of the triangular dislocation vertices
  !! @param s,d,t    - slip vector components (strike slip, dip slip 
  !!                    and tensile opening)
  !! @param G,lambda - elastic moduli (Lame parameters)
  !!
  !! OUTPUT:
  !! u1,u2,u3        - the displacement components
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeDisplacementStuart97(x1,x2,x3,q1,q2,q3, &
                                  s,d,t,G,lambda, &
                                  u)
    REAL*8, INTENT(IN) :: x1,x2,x3
    REAL*8, INTENT(IN), DIMENSION(3) :: q1,q2,q3
    REAL*8, INTENT(IN) :: s,d,t
    REAL*8, INTENT(IN) :: G,lambda
    REAL*8, DIMENSION(3), INTENT(OUT) :: u


    ! call C function
    !CALL DispInHalfSpace_()

    ! post-process as necessary
    !u1=
  END SUBROUTINE computeDisplacementStuart97

END MODULE stuart97




