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
! and
!
! Okada, Y. (1992), Internal deformation due to shear and tensile 
! faults in a half-space, Bull. Seism. Soc. Am., 82, 1018–1040.
!
!-----------------------------------------------------------------------

#include "macros.f90"

MODULE okada92

  IMPLICIT NONE

CONTAINS

  !------------------------------------------------------------------------
  !> subroutine computeReferenceSystemOkada92
  !! computes the center position and local reference system tied to the patch
  !!
  !! INPUT:
  !! @param ns
  !! @param x           - upper left coordinate of fault patch (north, east, down)
  !! @param strike      - strike angle
  !! @param dip         - dip angle (radian)
  !! @param L,W         - length and width of the rectangular dislocation
  !!
  !! OUTPUT:
  !! @param sv,dv,nv    - strike, dip and normal vectors of the fault patch
  !! @param xc          - coordinates (north, east, down) of the center
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeReferenceSystemOkada92(ns,x,L,W,strike,dip,sv,dv,nv,xc)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: x
    REAL*8, DIMENSION(ns), INTENT(IN) :: strike,dip,L,W
    REAL*8, DIMENSION(3,ns), INTENT(OUT) :: sv,dv,nv,xc

    ! unit vectors in the strike direction
    sv(1,1:ns)=COS(strike(1:ns))
    sv(2,1:ns)=SIN(strike(1:ns))
    sv(3,1:ns)=0._8
            
    ! unit vectors in the dip direction
    dv(1,1:ns)=+SIN(strike(1:ns))*COS(dip(1:ns))
    dv(2,1:ns)=-COS(strike(1:ns))*COS(dip(1:ns))
    dv(3,1:ns)=-SIN(dip(1:ns))
            
    ! unit vectors in the normal direction
    nv(1,1:ns)=-SIN(strike(1:ns))*SIN(dip(1:ns))
    nv(2,1:ns)=+COS(strike(1:ns))*SIN(dip(1:ns))
    nv(3,1:ns)=-COS(dip(1:ns))
            
    ! center of fault patch
    xc(1,1:ns)=x(1,1:ns)+L(1:ns)/2*sv(1,1:ns)-W(1:ns)/2*dv(1,1:ns)
    xc(2,1:ns)=x(2,1:ns)+L(1:ns)/2*sv(2,1:ns)-W(1:ns)/2*dv(2,1:ns)
    xc(3,1:ns)=x(3,1:ns)+L(1:ns)/2*sv(3,1:ns)-W(1:ns)/2*dv(3,1:ns)
                
                
  END SUBROUTINE computeReferenceSystemOkada92

  !------------------------------------------------------------------------
  !> subroutine computeStressOkada92
  !! calculates the stress associated with a rectangular dislocation in an 
  !! elastic half-space using the analytic solution of Okada et al. (1992).
  !!
  !! INPUT:
  !! @param x1,x2,x3    - coordinates of observation points
  !! @param q1,q2,q3    - coordinates of upper left corner of rectangular
  !!                      dislocation
  !! @param L,W         - length and width of the rectangular dislocation
  !! @param strike,dipd - strike and dip of the rectangular dislocation (radian)
  !! @param s,d,t       - slip vector components (strike slip, dip slip 
  !!                      and tensile opening)
  !! @param G,lambda    - elastic moduli (Lame parameters)
  !!
  !! OUTPUT:
  !! s11,s12,s13,s22,s23,s33 - the stress components
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeStressOkada92(x1,x2,x3, &
                                  q1,q2,q3,L,W,strike,dipd, &
                                  s,d,t,G,lambda, &
                                  s11,s12,s13,s22,s23,s33)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x1,x2,x3,q1,q2,q3
    REAL*8, INTENT(IN) :: L,W,strike,dipd
    REAL*8, INTENT(IN) :: s,d,t,G,lambda
    REAL*8, INTENT(OUT) :: s11,s12,s13,s22,s23,s33

    ! for subroutine dc3d
    INTEGER :: iret
    REAL*4 :: alpha,x,y,z,depth,dip,pot3,pot4, & 
              ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz

    ! rotated displacement gradient
    REAL*8 :: uxxr,uyxr,uzxr,uxyr,uyyr,uzyr,uxzr,uyzr
    ! strain components
    REAL*8 :: e11,e12,e13,e22,e23,e33,ekk

    ! for subroutine dc3d
    REAL*4 :: al1,al2,aw1,aw2,disl1,disl2,disl3

    ! local constants
    REAL*8, PARAMETER :: RAD2DEG=57.29577951308232_8

    REAL*8 :: csst,ssst

    ! receiver and source independent variables
    alpha=SNGL((lambda+G)/(lambda+2.d0*G))
    pot3=0.0
    pot4=0.0
    disl3=0.0
    al1=0.0
    aw2=0.0
    z=REAL(-x3)

    csst=DCOS(strike)
    ssst=DSIN(strike)

    ! rotate reference system to strike direction
    x=SNGL((x1-q1)*csst+(x2-q2)*ssst)
    y=SNGL((x1-q1)*ssst-(x2-q2)*csst)
    depth=SNGL(q3)
    dip=SNGL(dipd*RAD2DEG)

    ! finite source
    al2=SNGL(L)
    aw1=-SNGL(W)
    disl1=SNGL(s)
    disl2=SNGL(d)
    disl3=SNGL(t)

    iret=1

    ! The coordinate system of dc3d is as follows
    !
    !   x is along strike,
    !   z is up
    !   y is such that (x, y, z) is a right-handed coordinate system
    !
    CALL dc3d(alpha,x,y,z,depth,dip,al1,al2,aw1,aw2, &
              disl1,disl2,disl3,ux,uy,uz, &
              uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret)

    ! rotate displacement gradients to Cartesian coordinates
    uxxr=( csst*uxx+ssst*uyx)*csst &
        -(-csst*uxy-ssst*uyy)*ssst
    uxyr=( csst*uxx+ssst*uyx)*ssst &
        +(-csst*uxy-ssst*uyy)*csst
    uxzr=-csst*uxz+ssst*uyz
    uyxr=( ssst*uxx-csst*uyx)*csst &
        -(-ssst*uxy+csst*uyy)*ssst 
    uyyr=( ssst*uxx-csst*uyx)*ssst &
        +(-ssst*uxy+csst*uyy)*csst
    uyzr=-ssst*uxz+csst*uyz
    uzxr=-csst*uzx-ssst*uzy
    uzyr=-ssst*uzx+csst*uzy

    ! strain components
    e11=uxxr
    e12=(uxyr+uyxr)/2._8
    e13=(uxzr+uzxr)/2._8
    e22=uyyr
    e23=(uyzr+uzyr)/2._8
    e33=uzz

    ! isotropic strain
    ekk=e11+e22+e33

    ! stress components
    s11=lambda*ekk+2._8*G*e11
    s12=2._8*G*e12
    s13=2._8*G*e13
    s22=lambda*ekk+2._8*G*e22
    s23=2._8*G*e23
    s33=lambda*ekk+2._8*G*e33

  END SUBROUTINE computeStressOkada92

  !------------------------------------------------------------------------
  !> subroutine computeDisplacementOkada92
  !! calculates the displacements associated with a rectangular dislocation
  !! in an elastic half-space using the analytic solution of Okada (1992).
  !!
  !! INPUT:
  !! @param x1,x2,x3    - coordinates of observation points
  !! @param q1,q2,q3    - coordinates of upper left corner of rectangular
  !!                      dislocation
  !! @param L,W         - length and width of the rectangular dislocation
  !! @param strike,dipd - strike and dip of the rectangular dislocation
  !!                      in degrees.
  !! @param s,d,t       - slip vector components (strike slip, dip slip 
  !!                      and tensile opening)
  !! @param G,lambda    - elastic moduli (Lame parameters)
  !!
  !! OUTPUT:
  !! u1,u2,u3           - the displacement components
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeDisplacementOkada92(x1,x2,x3, &
                                        q1,q2,q3,L,W,strike,dipd, &
                                        s,d,t,G,lambda, &
                                        u1,u2,u3)
    REAL*8, INTENT(IN) :: x1,x2,x3,q1,q2,q3
    REAL*8, INTENT(IN) :: L,W,strike,dipd
    REAL*8, INTENT(IN) :: s,d,t,G,lambda
    REAL*8, INTENT(OUT) :: u1,u2,u3

    ! for subroutine dc3d
    INTEGER :: iret
    REAL*4 :: alpha,x,y,z,depth,dip,pot3,pot4, & 
              ux,uy,uz,uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz
    REAL*4 :: al1,al2,aw1,aw2,disl1,disl2,disl3

    ! local constants
    REAL*8, PARAMETER :: RAD2DEG=57.29577951308232_8

    REAL*8 :: csst,ssst

    ! receiver and source independent variables
    alpha=SNGL((lambda+G)/(lambda+2.d0*G))
    pot3=0.0
    pot4=0.0
    disl3=0.0
    al1=0.0
    aw2=0.0
    z=REAL(-x3)

    csst=COS(strike)
    ssst=SIN(strike)

    ! rotate reference system to strike direction
    x=SNGL((x1-q1)*csst+(x2-q2)*ssst)
    y=SNGL((x1-q1)*ssst-(x2-q2)*csst)
    depth=SNGL(q3)
    dip=SNGL(dipd*RAD2DEG)

    ! finite source
    al2=SNGL(L)
    aw1=-SNGL(W)
    disl1=SNGL(s)
    disl2=SNGL(d)
    disl3=SNGL(t)

    ! The coordinate system of dc3d is as follows
    !
    !   x is along strike,
    !   z is up
    !   y is such that (x, y, z) is a right-handed coordinate system
    !
    iret=1
    CALL dc3d(alpha,x,y,z,depth,dip,al1,al2,aw1,aw2, &
              disl1,disl2,disl3,ux,uy,uz, &
              uxx,uyx,uzx,uxy,uyy,uzy,uxz,uyz,uzz,iret)

    IF (0 .EQ. iret) THEN
       ! rotate to Cartesian coordinate system
       u1=+DBLE(UX)*csst+DBLE(UY)*ssst
       u2=+DBLE(UX)*ssst-DBLE(UY)*csst
       u3=-DBLE(UZ)
    ELSE
       ! singular point at fault edge
       u1=0._8
       u2=0._8
       u3=0._8
    ENDIF

  END SUBROUTINE computeDisplacementOkada92

  !------------------------------------------------------------------------
  !> subroutine computeDisplacementKernelsOkada92
  !! calculates the displacement kernels associated with dislocations in
  !! a three-dimensional elastic half-space using the solution of Okada (1992)
  !!
  !! INPUT:
  !! @param x          - coordinates of observation points
  !! @param ns         - number of sources
  !! @param y          - coordinate (NED) of sources
  !! @param L          - length (along strike) of the fault patch
  !! @param W          - width (down dip) of the fault patch
  !! @param strike     - strike angle (radian)
  !! @param dip        - dip angle (radian)
  !!
  !! OUTPUT:
  !! u1,u2,u3          - array, displacement in the north, east and depth 
  !!                     directions.
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeDisplacementKernelsOkada92(x, &
                        ns,y,L,W,strike,dip, &
                        ss,ds,ts, &
                        G,lambda, &
                        u1,u2,u3)
    REAL*8, DIMENSION(3), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: y
    REAL*8, DIMENSION(ns), INTENT(IN) :: L,W,strike,dip
    REAL*8, INTENT(IN) :: ss,ds,ts
    REAL*8, INTENT(IN) :: G,lambda
    REAL*8, DIMENSION(ns), INTENT(OUT) :: u1,u2,u3

    INTEGER :: i

    DO i=1,ns
       CALL computeDisplacementOkada92( &
                      x(1),x(2),x(3), &
                      y(1,i),y(2,i),y(3,i),L(i),W(i), &
                      strike(i),dip(i), &
                      ss,ds,ts, &
                      G,lambda, &
                      u1(i),u2(i),u3(i))
    END DO

  END SUBROUTINE computeDisplacementKernelsOkada92

  !------------------------------------------------------------------------
  !> subroutine computeTractionKernelsOkada92
  !! calculates the traction kernels associated with a triangular dislocation
  !! in an elastic half-space using the analytic solution of Okada (1992).
  !!
  !! INPUT:
  !! @param x          - coordinates of observation points
  !! @param sv,dv,nv   - strike, dip and normal vector of observation points
  !! @param ns         - number of sources
  !! @param y          - coordinate (NED) of sources
  !! @param L,W        - length and width of the rectangular dislocation
  !! @param strike,dip - strike and dip of the rectangular dislocation
  !! @param s,d,t      - slip vector components (strike slip, dip slip 
  !!                     and tensile opening)
  !! @param G,lambda   - elastic moduli (Lame parameters)
  !!
  !! OUTPUT:
  !! ts,td,tn          - traction in the strike, dip and normal directions.
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeTractionKernelsOkada92(x,sv,dv,nv, &
                        ns,y,L,W,strike,dip, &
                        s,d,t,G,lambda, &
                        ts,td,tn)
    REAL*8, DIMENSION(3), INTENT(IN) :: x,sv,dv,nv
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: y
    REAL*8, DIMENSION(ns), INTENT(IN) :: L,W,strike,dip
    REAL*8, INTENT(IN) :: s,d,t,G,lambda
    REAL*8, DIMENSION(ns), INTENT(OUT) :: ts,td,tn

    INTEGER :: i
    REAL*8 :: s11,s12,s13,s22,s23,s33

    DO i=1,ns
       CALL computeStressOkada92( &
                      x(1),x(2),x(3), &
                      y(1,i),y(2,i),y(3,i),L(i),W(i),strike(i),dip(i), &
                      s,d,t,G,lambda, &
                      s11,s12,s13,s22,s23,s33)

       ! rotate to receiver system of coordinates
       ts(i)= ( nv(1)*s11+nv(2)*s12+nv(3)*s13 )*sv(1) &
             +( nv(1)*s12+nv(2)*s22+nv(3)*s23 )*sv(2) &
             +( nv(1)*s13+nv(2)*s23+nv(3)*s33 )*sv(3)
       td(i)= ( nv(1)*s11+nv(2)*s12+nv(3)*s13 )*dv(1) &
             +( nv(1)*s12+nv(2)*s22+nv(3)*s23 )*dv(2) &
             +( nv(1)*s13+nv(2)*s23+nv(3)*s33 )*dv(3)
       tn(i)=+( nv(1)*s11+nv(2)*s12+nv(3)*s13 )*nv(1) &
             +( nv(1)*s12+nv(2)*s22+nv(3)*s23 )*nv(2) &
             +( nv(1)*s13+nv(2)*s23+nv(3)*s33 )*nv(3)

    END DO

  END SUBROUTINE computeTractionKernelsOkada92
  !------------------------------------------------------------------------
  !> subroutine computeStressKernelsOkada92
  !! calculates the stress kernels associated with a triangular dislocation
  !! in an elastic half-space using the analytic solution of Okada (1992).
  !!
  !! INPUT:
  !! @param x          - coordinates of observation points
  !! @param sv,dv,nv   - strike, dip and normal vector of observation points
  !! @param ns         - number of sources
  !! @param y          - coordinate (NED) of sources
  !! @param L,W        - length and width of the rectangular dislocation
  !! @param strike,dip - strike and dip of the rectangular dislocation
  !! @param s,d,t      - slip vector components (strike slip, dip slip 
  !!                     and tensile opening)
  !! @param G,lambda   - elastic moduli (Lame parameters)
  !!
  !! OUTPUT:
  !! s11,s12,s13,s22,s23,s33 - the stress components
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeStressKernelsOkada92(x,sv,dv,nv, &
                        ns,y,L,W,strike,dip, &
                        s,d,t,G,lambda,s11,s12,s13,s22,s23,s33)
    REAL*8, DIMENSION(3), INTENT(IN) :: x,sv,dv,nv
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: y
    REAL*8, DIMENSION(ns), INTENT(IN) :: L,W,strike,dip
    REAL*8, INTENT(IN) :: s,d,t,G,lambda
    REAL*8, DIMENSION(ns), INTENT(OUT) :: s11,s12,s13,s22,s23,s33

    INTEGER :: i
    REAL*8 :: s11p,s12p,s13p,s22p,s23p,s33p

    DO i=1,ns
       CALL computeStressOkada92( &
                      x(1),x(2),x(3), &
                      y(1,i),y(2,i),y(3,i),L(i),W(i),strike(i),dip(i), &
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

  END SUBROUTINE computeStressKernelsOkada92

END MODULE okada92




