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
! and
!
! Nikkhoo, M., and T. R. Walter (2015), Triangular dislocation: 
! an analytical, artefact-free solution, Geophys. J. Int., 201(2), 
! 1119–1141, https://doi.org/10.1093/gji/ggv035.
!
!-----------------------------------------------------------------------

#include "macros.f90"

MODULE nikkhoo15

  IMPLICIT NONE

CONTAINS

  !------------------------------------------------------------------------
  !> subroutine computeReferenceSystemNikkhoo15
  !! computes the center position and local reference system of triangular
  !! dislocations.
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeReferenceSystemNikkhoo15(nVe,v,ns,i1,i2,i3,sv,dv,nv,xc)
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

  END SUBROUTINE computeReferenceSystemNikkhoo15

  !------------------------------------------------------------------------
  !> subroutine computeStressNikkhoo15
  !! calculates the traction associated with a triangular dislocation in an 
  !! elastic half-space using the analytic solution of Nikkhoo et al. (2015).
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
  SUBROUTINE computeStressNikkhoo15(x1,x2,x3, &
                                    P1,P2,P3, &
                                    ss,d,t,G,lambda, &
                                    s11,s12,s13,s22,s23,s33)
    REAL*8, INTENT(IN) :: x1,x2,x3
    REAL*8, DIMENSION(3), INTENT(IN) :: P1,P2,P3
    REAL*8, INTENT(IN) :: ss,d,t
    REAL*8, INTENT(IN) :: G,lambda
    REAL*8, INTENT(OUT) :: s11,s12,s13,s22,s23,s33

    REAL*8, DIMENSION(3) :: Q1,Q2,Q3
    REAL*8, DIMENSION(6) :: s,e
    REAL*8 :: x,y,z,ds

    ! convert observation to east,north,up
    x=x2
    y=x1
    z=-x3

    ! convert dip slip
    ds=-d

    ! convert vertex position to east,north,up
    Q1=(/ P1(2), P1(1), -P1(3) /)
    Q2=(/ P2(2), P2(1), -P2(3) /)
    Q3=(/ P3(2), P3(1), -P3(3) /)

    ! call C function straininhalfspace (in lower case)
    CALL strainInhalfSpace( &
                     s,e, &
                     %VAL(x),%VAL(y),%VAL(z), &
                     Q1,Q2,Q3, &
                     %VAL(ss),%VAL(ds),%VAL(t), &
                     %VAL(G),%VAL(lambda))

    ! change from east,north,up to north,east,down
    s11=s(4)
    s12=s(2)
    s13=-s(5)
    s22=s(1)
    s23=-s(3)
    s33=s(6)

  END SUBROUTINE computeStressNikkhoo15

  !------------------------------------------------------------------------
  !> subroutine computeTractionKernelsNikkhoo15
  !! calculates the traction kernels associated with the motion of triangular
  !! dislocations in an elastic half space using the analytic solution of
  !!
  !!    Nikkhoo, Mehdi, and Thomas R. Walter. "Triangular dislocation: 
  !!    an analytical, artefact-free solution." Geophys. J. Int. 201(2), 
  !!    1117-1139, 2015.
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
  !! ts,td,tn                - traction in the strike, dip and normal directions
  !------------------------------------------------------------------------
  SUBROUTINE computeTractionKernelsNikkhoo15(x,sv,dv,nv, &
                                           ns,i1,i2,i3, &
                                           nVe,vertices, &
                                           s,d,t,G,lambda, &
                                           ts,td,tn)
    REAL*8, DIMENSION(3), INTENT(IN) :: x
    REAL*8, DIMENSION(3), INTENT(IN) :: sv,dv,nv
    INTEGER, INTENT(IN) :: ns
    INTEGER, DIMENSION(ns), INTENT(IN) :: i1,i2,i3
    INTEGER, INTENT(IN) :: nVe
    REAL*8, DIMENSION(3,nVe), INTENT(IN) :: vertices
    REAL*8, INTENT(IN) :: s,d,t,G,lambda
    REAL*8, DIMENSION(ns), INTENT(OUT) :: ts,td,tn

    INTEGER :: i
    REAL*8 :: s11,s12,s13,s22,s23,s33

    ! loop over sources
    DO i=1,ns
       CALL computeStressNikkhoo15( &
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

  END SUBROUTINE computeTractionKernelsNikkhoo15

  !------------------------------------------------------------------------
  !> subroutine computeStressKernelsNikkhoo15
  !! calculates the stress kernels associated with the motion of triangular
  !! dislocations in an elastic half space using the analytic solution of
  !!
  !!    Nikkhoo, Mehdi, and Thomas R. Walter. "Triangular dislocation: 
  !!    an analytical, artefact-free solution." Geophys. J. Int. 201(2), 
  !!    1117-1139, 2015.
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
  SUBROUTINE computeStressKernelsNikkhoo15(x,sv,dv,nv, &
                                           ns,i1,i2,i3, &
                                           nVe,vertices, &
                                           s,d,t,G,lambda, &
                                           s11,s12,s13,s22,s23,s33)
    REAL*8, DIMENSION(3), INTENT(IN) :: x
    REAL*8, DIMENSION(3), INTENT(IN) :: sv,dv,nv
    INTEGER, INTENT(IN) :: ns
    INTEGER, DIMENSION(ns), INTENT(IN) :: i1,i2,i3
    INTEGER, INTENT(IN) :: nVe
    REAL*8, DIMENSION(3,nVe), INTENT(IN) :: vertices
    REAL*8, INTENT(IN) :: s,d,t,G,lambda
    REAL*8, DIMENSION(ns), INTENT(OUT) :: s11,s12,s13,s22,s23,s33

    INTEGER :: i
    REAL*8 :: s11p,s12p,s13p,s22p,s23p,s33p

    ! loop over sources
    DO i=1,ns
       CALL computeStressNikkhoo15( &
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

  END SUBROUTINE computeStressKernelsNikkhoo15

  !------------------------------------------------------------------------
  !> subroutine computeDisplacementNikkhoo15
  !! calculates the displacements associated with a triangular dislocation
  !! in an elastic half-space using the analytic solution of Nikkhoo et al.
  !! (2015).
  !!
  !! INPUT:
  !! @param x1,x2,x3 - coordinates of observation points
  !! @param q1,q2,q2 - coordinates of the triangular dislocation vertices
  !! @param s,d,t    - slip vector components (strike slip, dip slip 
  !!                   and tensile opening)
  !! @param nu       - Poisson's ratio
  !!
  !! OUTPUT:
  !! u1,u2,u3        - the displacement components
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeDisplacementNikkhoo15(x1,x2,x3, &
                                    P1,P2,P3, &
                                    ss,d,t,G,lambda, &
                                    u1,u2,u3)
    REAL*8, INTENT(IN) :: x1,x2,x3
    REAL*8, DIMENSION(3), INTENT(IN) :: P1,P2,P3
    REAL*8, INTENT(IN) :: ss,d,t
    REAL*8, INTENT(IN) :: G,lambda
    REAL*8, INTENT(OUT) :: u1,u2,u3

    REAL*8, DIMENSION(3) :: Q1,Q2,Q3
    REAL*8, DIMENSION(3) :: u
    REAL*8 :: x,y,z,ds,nu

    ! Poisson ratio
    nu=lambda/(lambda+G)/2d0

    ! convert observation to east,north,up
    x=x2
    y=x1
    z=-x3

    ! convert dip slip
    ds=-d

    ! convert vertex position to east,north,up
    Q1=(/ P1(2), P1(1), -P1(3) /)
    Q2=(/ P2(2), P2(1), -P2(3) /)
    Q3=(/ P3(2), P3(1), -P3(3) /)

    ! call C function dispinhalfspace (in lower case)
    CALL dispInhalfSpace( &
                     u, &
                     %VAL(x),%VAL(y),%VAL(z), &
                     Q1,Q2,Q3, &
                     %VAL(ss),%VAL(ds),%VAL(t), &
                     %VAL(nu))

    ! change from east,north,up to north,east,down
    u1=u(2)
    u2=u(1)
    u3=-u(3)

  END SUBROUTINE computeDisplacementNikkhoo15

  !------------------------------------------------------------------------
  !> subroutine computeDisplacementKernelsNikkhoo15
  !! calculates the displacement kernels associated with the motion of triangular
  !! dislocations in an elastic half space using the analytic solution of
  !!
  !!    Nikkhoo, Mehdi, and Thomas R. Walter. "Triangular dislocation: 
  !!    an analytical, artefact-free solution." Geophys. J. Int. 201(2), 
  !!    1117-1139, 2015.
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
  !! u1,u2,u3                - array displacement in the strike, dip and 
  !!                           normal directions
  !------------------------------------------------------------------------
  SUBROUTINE computeDisplacementKernelsNikkhoo15(x, &
                                           ns,i1,i2,i3, &
                                           nVe,vertices, &
                                           s,d,t,G,lambda, &
                                           u1,u2,u3)
    REAL*8, DIMENSION(3), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: ns
    INTEGER, DIMENSION(ns), INTENT(IN) :: i1,i2,i3
    INTEGER, INTENT(IN) :: nVe
    REAL*8, DIMENSION(3,nVe), INTENT(IN) :: vertices
    REAL*8, INTENT(IN) :: s,d,t,G,lambda
    REAL*8, DIMENSION(ns), INTENT(OUT) :: u1,u2,u3

    INTEGER :: i

    ! loop over sources
    DO i=1,ns
       CALL computeDisplacementNikkhoo15( &
                     x(1),x(2),x(3), &
                     vertices(:,i1(i)),vertices(:,i2(i)),vertices(:,i3(i)), &
                     s,d,t,G,lambda, &
                     u1(i),u2(i),u3(i))

    END DO

  END SUBROUTINE computeDisplacementKernelsNikkhoo15

END MODULE nikkhoo15




