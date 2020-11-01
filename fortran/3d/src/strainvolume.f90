!-----------------------------------------------------------------------
! Copyright (c) 2017 Sylvain Barbot
!
! This file is part of UNICYCLE
!
! UNICYCLE is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! UNICYCLE is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with UNICYCLE.  If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------

#include "macros.f90"

MODULE strainVolume

  IMPLICIT NONE

  PUBLIC

  REAL*8, PRIVATE, PARAMETER :: PI = 3.141592653589793115997963468544185161_8
    
CONTAINS

  !------------------------------------------------------------------------
  !> function Omega(x)
  !! evaluates the boxcar function
  !------------------------------------------------------------------------
  REAL*8 FUNCTION omega(x)
    REAL*8, INTENT(IN) :: x

    omega=heaviside(x+0.5_8)-heaviside(x-0.5_8)

  END FUNCTION omega

  !------------------------------------------------------------------------
  !> function S(x)
  !! evalutes the shifted boxcar function
  !------------------------------------------------------------------------
  REAL*8 FUNCTION s(x)
    REAL*8, INTENT(IN) :: x

    s=omega(x-0.5_8)

  END FUNCTION s

  !------------------------------------------------------------------------
  !> function heaviside
  !! computes the Heaviside function
  !------------------------------------------------------------------------
  REAL*8 FUNCTION heaviside(x)
    REAL*8, INTENT(IN) :: x

     IF (0 .LT. x) THEN
        heaviside=1.0_8
     ELSE
        heaviside=0.0_8
     END IF

  END FUNCTION heaviside

  !------------------------------------------------------------------------
  !> function atan3
  !! computes atan2 with the required value at infinity
  !------------------------------------------------------------------------
  REAL*8 FUNCTION atan3(y,x)
    REAL*8, INTENT(IN) :: y,x

    IF (0 .EQ. x) THEN
       atan3=SIGN(PI/2,y)
    ELSE
       atan3=ATAN(y/x)
    END IF

  END FUNCTION atan3

  !------------------------------------------------------------------------
  !> function xLogY
  !! computes x*log(y) and enforces 0*log(0)=0 to avoid NaN
  !------------------------------------------------------------------------
  REAL*8 FUNCTION xLogy(x,y)
    REAL*8, INTENT(IN) :: x,y

    IF (0 .EQ. x) THEN
       xLogy=0._8
    ELSE
       xLogy=x*log(y)
    END IF

  END FUNCTION xLogy

  !------------------------------------------------------------------------
  !> subroutine computeReferenceSystemVerticalStrainVolume
  !! computes the center position and local reference system tied to the volume
  !!
  !! INPUT:
  !! @param ns
  !! @param x           - upper left coordinate of fault patch (north, east, down)
  !! @param strike      - strike angle
  !! @param L,W         - length and width of the rectangular dislocation
  !!
  !! OUTPUT:
  !! @param sv,dv,nv    - strike, dip and normal vectors of the fault patch
  !! @param xc          - coordinates (north, east, down) of the center
  !!
  !! \author Sylvain Barbot (21/02/17) - original fortran form
  !------------------------------------------------------------------------
  SUBROUTINE computeReferenceSystemVerticalStrainVolume(ns,x,L,W,strike,sv,dv,nv,xc)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: x
    REAL*8, DIMENSION(ns), INTENT(IN) :: strike,L,W
    REAL*8, DIMENSION(3,ns), INTENT(OUT) :: sv,dv,nv,xc

    ! unit vectors in the strike direction
    sv(1,1:ns)=COS(strike(1:ns))
    sv(2,1:ns)=SIN(strike(1:ns))
    sv(3,1:ns)=0._8
            
    ! unit vectors in the dip direction
    dv(1,1:ns)=0._8
    dv(2,1:ns)=0._8
    dv(3,1:ns)=-1._8
            
    ! unit vectors in the normal direction
    nv(1,1:ns)=-SIN(strike(1:ns))
    nv(2,1:ns)=+COS(strike(1:ns))
    nv(3,1:ns)=-0._8
            
    ! center of fault patch
    xc(1,1:ns)=x(1,1:ns)+L(1:ns)/2*sv(1,1:ns)-W(1:ns)/2*dv(1,1:ns)
    xc(2,1:ns)=x(2,1:ns)+L(1:ns)/2*sv(2,1:ns)-W(1:ns)/2*dv(2,1:ns)
    xc(3,1:ns)=x(3,1:ns)+L(1:ns)/2*sv(3,1:ns)-W(1:ns)/2*dv(3,1:ns)
                
                
  END SUBROUTINE computeReferenceSystemVerticalStrainVolume

  !------------------------------------------------------------------------
  !> subroutine ComputeDisplacementVerticalStrainVolume computes the stress
  !! field associated with deforming vertical strain volume using the
  !! analytic solution of 
  !!
  !!   Barbot S., J. D. P. Moore and V. Lambert, Displacement and Stress
  !!   Associated with Distributed Anelastic Deformation in a Half Space,
  !!   Bull. Seism. Soc. Am., 107(2), 10.1785/0120160237, 2017.
  !!
  !! considering the following geometry:
  !!
  !!
  !!                      N (x1p)
  !!                     /
  !!                    /| strike (theta)          E (x2p)
  !!        q1,q2,q3 ->@--------------------------+
  !!                   |                        w |     +
  !!                   |                        i |    /
  !!                   |                        d |   / s
  !!                   |                        t |  / s
  !!                   |                        h | / e
  !!                   |                          |/ n
  !!                   +--------------------------+  k
  !!                   :       l e n g t h       /  c
  !!                   |                        /  i
  !!                   :                       /  h
  !!                   |                      /  t
  !!                   :                     /
  !!                   |                    +
  !!                   Z (x3)
  !!
  !!
  !! INPUT:
  !! @param x1, x2, x3         northing, easting, and depth of the observation point
  !!                           in unprimed system of coordinates.
  !! @param q1, q2, q3         north, east and depth coordinates of the strain volume,
  !! @param L, T, W            length, thickness, and width of the strain volume,
  !! @param theta (radians)    strike of the strain volume,
  !! @param epsijp             anelastic strain component 11, 12, 13, 22, 23 and 33 
  !!                           in the strain volume in the system of reference tied to 
  !!                           the strain volume (primed reference system),
  !! @param G, nu              shear modulus and Poisson's ratio in the half space.
  !!
  !! OUTPUT:
  !! ui                        displacement components in the unprimed reference system.
  !!
  !! \author James D. P. Moore (10/06/16) - J-functions in original form
  !! \author Sylvain Barbot (21/02/17) - original fortran form
  !------------------------------------------------------------------------
  SUBROUTINE computeDisplacementVerticalStrainVolume( &
                          x1,x2,x3,q1,q2,q3,L,T,W,theta, &
                          eps11p,eps12p,eps13p,eps22p,eps23p,eps33p,G,nu, &
                          u1,u2,u3)

    REAL*8, INTENT(IN) :: x1,x2,x3,q1,q2,q3,L,T,W,theta
    REAL*8, INTENT(IN) :: eps11p,eps12p,eps13p,eps22p,eps23p,eps33p
    REAL*8, INTENT(IN) :: G,nu
    REAL*8, INTENT(OUT) :: u1,u2,u3
    
    REAL*8 :: lambda
    REAL*8 :: t1
    REAL*8 :: epskk
    REAL*8 :: x1p,x2p

    ! check valid parameters
    IF ((-1._8 .GT. nu) .OR. (0.5_8 .LT. nu)) THEN
       WRITE (0,'("error: -1<=nu<=0.5, nu=",ES9.2E2," given.")') nu
       STOP 1
    END IF

    IF (0 .GT. x3) THEN
       WRITE (0,'("error: observation depth (x3) must be positive")')
       STOP 1
    END IF

    IF (0 .GT. q3) THEN
       WRITE (0,'("error: source depth (q3) must be positive")')
       STOP 1
    END IF

    ! lame elastic parameters
    lambda=G*2._8*nu/(1-2*nu)

    ! isotropic strain
    epskk=eps11p+eps22p+eps33p

    ! rotate observation points to the strain-volume-centric system of coordinates
    x1p= (x1-q1)*DCOS(theta)+(x2-q2)*DSIN(theta)
    x2p=-(x1-q1)*DSIN(theta)+(x2-q2)*DCOS(theta)

    u1= IU1(L,   T/2,q3+W)-IU1(L,   -T/2,q3+W)+IU1(L,   -T/2,q3)-IU1(L,   T/2,q3) &
       -IU1(0._8,T/2,q3+W)+IU1(0._8,-T/2,q3+W)-IU1(0._8,-T/2,q3)+IU1(0._8,T/2,q3)
    u2= IU2(L,   T/2,q3+W)-IU2(L,   -T/2,q3+W)+IU2(L,   -T/2,q3)-IU2(L,   T/2,q3) &
       -IU2(0._8,T/2,q3+W)+IU2(0._8,-T/2,q3+W)-IU2(0._8,-T/2,q3)+IU2(0._8,T/2,q3)
    u3= IU3(L,   T/2,q3+W)-IU3(L,   -T/2,q3+W)+IU3(L,   -T/2,q3)-IU3(L,   T/2,q3) &
       -IU3(0._8,T/2,q3+W)+IU3(0._8,-T/2,q3+W)-IU3(0._8,-T/2,q3)+IU3(0._8,T/2,q3)

    ! rotate displacement field to reference system of coordinates
    t1=u1*DCOS(theta)-u2*DSIN(theta)
    u2=u1*DSIN(theta)+u2*DCOS(theta)
    u1=t1

  CONTAINS

    !------------------------------------------------------------------------
    !> function r1
    !! computes the distance from the source at y1,y2,y3
    !------------------------------------------------------------------------
    REAL*8 FUNCTION r1(y1,y2,y3) 
      REAL*8, INTENT(IN) :: y1,y2,y3
    
      r1=sqrt((x1p-y1)**2+(x2p-y2)**2+(x3-y3)**2)
   
    END FUNCTION r1

    !------------------------------------------------------------------------
    !> function r2
    !! computes the distance from the image at y1,y2,-y3
    !------------------------------------------------------------------------
    REAL*8 FUNCTION r2(y1,y2,y3) 
      REAL*8, INTENT(IN) :: y1,y2,y3

      r2=sqrt((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)

    END FUNCTION r2

    !---------------------------------------------------------------
    !> function IU1
    !! computes the indefinite integral U1 
    !---------------------------------------------------------------
    REAL*8 FUNCTION IU1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      IU1=(lambda*epskk+2*G*eps11p)*J1123(y1,y2,y3) &
                        +2*G*eps12p*(J1223(y1,y2,y3)+J1113(y1,y2,y3)) &
                        +2*G*eps13p*(J1323(y1,y2,y3)+J1112(y1,y2,y3)) &
         +(lambda*epskk+2*G*eps22p)*J1213(y1,y2,y3) &
                        +2*G*eps23p*(J1212(y1,y2,y3)+J1313(y1,y2,y3)) &
         +(lambda*epskk+2*G*eps33p)*J1312(y1,y2,y3)

    END FUNCTION IU1

    !---------------------------------------------------------------
    !> function IU2
    !! computes the indefinite integral U2
    !---------------------------------------------------------------
    REAL*8 FUNCTION IU2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      IU2=(lambda*epskk+2*G*eps11p)*J2123(y1,y2,y3) &
                        +2*G*eps12p*(J2223(y1,y2,y3)+J2113(y1,y2,y3)) &
                        +2*G*eps13p*(J2323(y1,y2,y3)+J2112(y1,y2,y3)) &
         +(lambda*epskk+2*G*eps22p)*J2213(y1,y2,y3) &
                        +2*G*eps23p*(J2212(y1,y2,y3)+J2313(y1,y2,y3)) &
         +(lambda*epskk+2*G*eps33p)*J2312(y1,y2,y3)

    END FUNCTION IU2

    !---------------------------------------------------------------
    !> function IU3
    !! computes the indefinite integral U3
    !---------------------------------------------------------------
    REAL*8 FUNCTION IU3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      IU3=(lambda*epskk+2*G*eps11p)*J3123(y1,y2,y3) &
                        +2*G*eps12p*(J3223(y1,y2,y3)+J3113(y1,y2,y3)) &
                        +2*G*eps13p*(J3323(y1,y2,y3)+J3112(y1,y2,y3)) &
         +(lambda*epskk+2*G*eps22p)*J3213(y1,y2,y3) &
                        +2*G*eps23p*(J3212(y1,y2,y3)+J3313(y1,y2,y3)) &
         +(lambda*epskk+2*G*eps33p)*J3312(y1,y2,y3)

    END FUNCTION IU3

    !---------------------------------------------------------------
    !> function J1112
    !! computes the J integral J1112
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1112(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1112=((1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(2*lr2**(-1)*x3*( &
            x1p-y1)*(x2p-y2)*y3*((x1p-y1)**2+(x3+y3)**2)**( &
            -1)-4*(-1._8+nu)*((-1)+2*nu)*(x3+y3)*atan3((x2p-y2),(x1p-y1)) &
            -x3*atan2(x3,x1p-y1)-3*x3* &
            atan2(3*x3,x1p-y1)+4*nu*x3*atan2(-nu*x3,x1p- &
            y1)+4*(-1._8+nu)*((-1)+2*nu)*(x3+y3)*atan2(lr2*(-x1p+y1),( &
            x2p-y2)*(x3+y3))-4*(-1._8+nu)*(x3-y3)*atan2(lr1*( &
            x3-y3),(x1p-y1)*(x2p-y2))+3*y3*atan2((-3)*y3, &
            x1p-y1)-y3*atan2(y3,x1p-y1)-4*nu*y3*atan2( &
            nu*y3,x1p-y1)-4*(-1._8+nu)*(x3+y3)*atan2(lr2*(x3+y3),( &
            x1p-y1)*(x2p-y2))+xLogy(-((-3)+4*nu)*(x1p- &
            y1),lr1+x2p-y2)+xLogy((5+4*nu*((-3)+2*nu))*(x1p-y1), &
            lr2+x2p-y2)+xLogy((-4)*(-1._8+nu)*(x2p-y2),lr1+x1p- &
            y1)+xLogy((-4)*(-1._8+nu)*(x2p-y2),lr2+x1p-y1))

    END FUNCTION J1112

    !---------------------------------------------------------------
    !> function J1113
    !! computes the J integral J1113
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1113(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1113=((1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(2*lr2**(-1)*(x1p+( &
            -1)*y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(-((-1)+ &
            nu)*((-1)+2*nu)*lr2**2*(x3+y3)+(-1._8+nu)*((-1)+2*nu)*lr2* &
            y3*(2*x3+y3)+x3*((x1p-y1)**2+(x2p-y2)**2+x3*(x3+y3)) &
            )+x2p*atan2(-x2p,x1p-y1)-3*x2p*atan2(3*x2p,x1p- &
            y1)+4*nu*x2p*atan2(-nu*x2p,x1p-y1)-4*(-1._8+nu)*( &
            x2p-y2)*atan2(lr1*(x2p-y2),(x1p-y1)*(x3-y3) &
            )+4*(-1._8+nu)*(x2p-y2)*atan2(lr2*(x2p-y2),(x1p- &
            y1)*(x3+y3))+3*y2*atan2((-3)*y2,x1p-y1)-y2*atan2( &
            y2,x1p-y1)-4*nu*y2*atan2(nu*y2,x1p-y1)+xLogy((-1) &
            *((-3)+4*nu)*(x1p-y1),lr1+x3-y3)+xLogy(-(3 &
            -6*nu+4*nu**2)*(x1p-y1),lr2+x3+y3)+xLogy((-4)*(-1._8+nu)*( &
            x3-y3),lr1+x1p-y1)+xLogy(4*(-1._8+nu)*(x3+y3),lr2+x1p+( &
            -1)*y1))

    END FUNCTION J1113

    !---------------------------------------------------------------
    !> function J1123
    !! computes the J integral J1123
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1123(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1123=((1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*((-2)*lr2**(-1)*(( &
            x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*((x1p+(-1) &
            *y1)**2+(x3+y3)**2)**(-1)*(x3*((x3**2+(x1p-y1)**2)*( &
            x3**2+(x1p-y1)**2+(x2p-y2)**2)+x3*(3*x3**2+2*(x1p+(-1) &
            *y1)**2+(x2p-y2)**2)*y3+3*x3**2*y3**2+x3*y3**3)-(( &
            -1)+nu)*((-1)+2*nu)*lr2**2*(x3+y3)*((x1p-y1)**2+(x3+y3) &
            **2)+(-1._8+nu)*((-1)+2*nu)*lr2*y3*(2*x3+y3)*((x1p-y1) &
            **2+(x3+y3)**2))+2*(-1._8+nu)*((-1)+2*nu)*(x1p-y1)*atan3((x1p-y1),(x2p-y2)) &
            +x1p*atan2(-x1p,x2p-y2) &
            -3*x1p*atan2(3*x1p,x2p-y2)+4*nu*x1p*atan2(-nu*x1p, &
            x2p-y2)+3*y1*atan2((-3)*y1,x2p-y2)-y1*atan2( &
            y1,x2p-y2)-4*nu*y1*atan2(nu*y1,x2p-y2)+2*((-1)+ &
            2*nu)*(x1p-y1)*atan2(lr1*(-x1p+y1),(x2p-y2)*(x3+ &
            (-1)*y3))+2*(1-2*nu)**2*(x1p-y1)*atan2(lr2*(-x1p+ &
            y1),(x2p-y2)*(x3+y3))+xLogy((-2)*x3,lr2-x2p+y2)+xLogy(( &
            -1)*((-3)+4*nu)*(x2p-y2),lr1+x3-y3)+xLogy(-(3+( &
            -6)*nu+4*nu**2)*(x2p-y2),lr2+x3+y3)+xLogy(-((-3)+4* &
            nu)*(x3-y3),lr1+x2p-y2)+xLogy(-(5+4*nu*((-3)+2* &
            nu))*(x3+y3),lr2+x2p-y2))

    END FUNCTION J1123

    !---------------------------------------------------------------
    !> function J2112
    !! computes the J integral J2112
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2112(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2112=((1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(-lr1+(1+8*(( &
            -1)+nu)*nu)*lr2-2*lr2**(-1)*x3*y3+xLogy((-4)*(-1._8+nu)*(( &
            -1)+2*nu)*(x3+y3),lr2+x3+y3))

    END FUNCTION J2112

    !---------------------------------------------------------------
    !> function J2113
    !! computes the J integral J2113
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2113(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2113=((1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(2*lr2**(-1)*((x1p+ &
            (-1)*y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(-((-1)+ &
            nu)*((-1)+2*nu)*lr2**2*(x3+y3)+(-1._8+nu)*((-1)+2*nu)*lr2* &
            y3*(2*x3+y3)+x3*((x1p-y1)**2+(x2p-y2)**2+x3*(x3+y3)) &
            )+xLogy(-((-1)-2*nu+4*nu**2)*(x2p-y2),lr2+x3+y3)+ &
            xLogy(-x2p+y2,lr1+x3-y3))

    END FUNCTION J2113

    !---------------------------------------------------------------
    !> function J2123
    !! computes the J integral J2123
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2123(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2123=((1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(2*lr2**(-1)*(x1p+( &
            -1)*y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(-((-1)+ &
            nu)*((-1)+2*nu)*lr2**2*(x3+y3)+(-1._8+nu)*((-1)+2*nu)*lr2* &
            y3*(2*x3+y3)+x3*((x1p-y1)**2+(x2p-y2)**2+x3*(x3+y3)) &
            )+xLogy(-((-1)-2*nu+4*nu**2)*(x1p-y1),lr2+x3+y3)+ &
            xLogy(-x1p+y1,lr1+x3-y3))

    END FUNCTION J2123

    !---------------------------------------------------------------
    !> function J3112
    !! computes the J integral J3112
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3112(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3112=(-(1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*((-2)*lr2**(-1)* &
            x3*(x2p-y2)*y3*(x3+y3)*((x1p-y1)**2+(x3+y3)**2)**( &
            -1)+4*(-1._8+nu)*((-1)+2*nu)*(x1p-y1)*atan3((x1p-y1),(x2p-y2)) &
            +4*(-1._8+nu)*((-1)+2*nu)*(x1p-y1)* &
            atan2(lr2*(-x1p+y1),(x2p-y2)*(x3+y3))+xLogy((-4)*((-1)+ &
            nu)*((-1)+2*nu)*(x2p-y2),lr2+x3+y3)+xLogy(x3-y3,lr1+ &
            x2p-y2)+xLogy(-x3-7*y3-8*nu**2*(x3+y3)+8*nu*( &
            x3+2*y3),lr2+x2p-y2))

    END FUNCTION J3112

    !---------------------------------------------------------------
    !> function J3113
    !! computes the J integral J3113
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3113(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3113=(-(1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(lr1+((-1)-8*(( &
            -1)+nu)*nu)*lr2-2*lr2**(-1)*x3*y3+2*((-3)+4*nu)*x3* &
            ACOTH(lr2**(-1)*(x3+y3))+xLogy(2*(3*x3+2*y3-6*nu*(x3+y3)+ &
            4*nu**2*(x3+y3)),lr2+x3+y3))

    END FUNCTION J3113

    !---------------------------------------------------------------
    !> function J3123
    !! computes the J integral J3123
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3123(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3123=(-(1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(2*lr2**(-1)*x3*( &
            x1p-y1)*(x2p-y2)*y3*((x1p-y1)**2+(x3+y3)**2)**( &
            -1)+4*(-1._8+nu)*((-1)+2*nu)*(x3+y3)*atan3((x2p-y2),(x1p-y1)) &
            +4*((-1)+2*nu)*(nu*x3+(-1._8+nu)*y3)*atan2( &
            lr2*(x1p-y1),(x2p-y2)*(x3+y3))+xLogy(x1p-y1,lr1+x2p+ &
            (-1)*y2)+xLogy(-(1+8*(-1._8+nu)*nu)*(x1p-y1),lr2+x2p+( &
            -1)*y2))

    END FUNCTION J3123

    !---------------------------------------------------------------
    !> function J1212
    !! computes the J integral J1212
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1212(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1212=((1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(-lr1+(1+8*(( &
            -1)+nu)*nu)*lr2-2*lr2**(-1)*x3*y3+xLogy((-4)*(-1._8+nu)*(( &
            -1)+2*nu)*(x3+y3),lr2+x3+y3))

    END FUNCTION J1212

    !---------------------------------------------------------------
    !> function J1213
    !! computes the J integral J1213
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1213(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1213=((1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(2*lr2**(-1)*((x1p+ &
            (-1)*y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(-((-1)+ &
            nu)*((-1)+2*nu)*lr2**2*(x3+y3)+(-1._8+nu)*((-1)+2*nu)*lr2* &
            y3*(2*x3+y3)+x3*((x1p-y1)**2+(x2p-y2)**2+x3*(x3+y3)) &
            )+xLogy(-((-1)-2*nu+4*nu**2)*(x2p-y2),lr2+x3+y3)+ &
            xLogy(-x2p+y2,lr1+x3-y3))

    END FUNCTION J1213

    !---------------------------------------------------------------
    !> function J1223
    !! computes the J integral J1223
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1223(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1223=((1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(2*lr2**(-1)*(x1p+( &
            -1)*y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(-((-1)+ &
            nu)*((-1)+2*nu)*lr2**2*(x3+y3)+(-1._8+nu)*((-1)+2*nu)*lr2* &
            y3*(2*x3+y3)+x3*((x1p-y1)**2+(x2p-y2)**2+x3*(x3+y3)) &
            )+xLogy(-((-1)-2*nu+4*nu**2)*(x1p-y1),lr2+x3+y3)+ &
            xLogy(-x1p+y1,lr1+x3-y3))

    END FUNCTION J1223

    !---------------------------------------------------------------
    !> function J2212
    !! computes the J integral J2212
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2212(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2212=((1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(2*lr2**(-1)*x3*( &
            x1p-y1)*(x2p-y2)*y3*((x2p-y2)**2+(x3+y3)**2)**( &
            -1)-4*(-1._8+nu)*((-1)+2*nu)*(x3+y3)*atan3((x1p-y1),(x2p-y2)) &
            -x3*atan2(x3,x1p-y1)-3*x3* &
            atan2(3*x3,x1p-y1)+4*nu*x3*atan2(-nu*x3,x1p- &
            y1)+4*(-1._8+nu)*((-1)+2*nu)*(x3+y3)*atan2(lr2*(-x2p+y2),( &
            x1p-y1)*(x3+y3))-4*(-1._8+nu)*(x3-y3)*atan2(lr1*( &
            x3-y3),(x1p-y1)*(x2p-y2))+3*y3*atan2((-3)*y3, &
            x1p-y1)-y3*atan2(y3,x1p-y1)-4*nu*y3*atan2( &
            nu*y3,x1p-y1)-4*(-1._8+nu)*(x3+y3)*atan2(lr2*(x3+y3),( &
            x1p-y1)*(x2p-y2))+xLogy((-4)*(-1._8+nu)*(x1p-y1), &
            lr1+x2p-y2)+xLogy((-4)*(-1._8+nu)*(x1p-y1),lr2+x2p- &
            y2)+xLogy(-((-3)+4*nu)*(x2p-y2),lr1+x1p-y1)+xLogy( &
            (5+4*nu*((-3)+2*nu))*(x2p-y2),lr2+x1p-y1))

    END FUNCTION J2212

    !---------------------------------------------------------------
    !> function J2213
    !! computes the J integral J2213
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2213(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2213=((1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*((-2)*lr2**(-1)*( &
            x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*((x2p+(-1) &
            *y2)**2+(x3+y3)**2)**(-1)*(x3*((x3**2+(x2p-y2)**2)*( &
            x3**2+(x1p-y1)**2+(x2p-y2)**2)+x3*(3*x3**2+(x1p- &
            y1)**2+2*(x2p-y2)**2)*y3+3*x3**2*y3**2+x3*y3**3)-( &
            (-1)+nu)*((-1)+2*nu)*lr2**2*(x3+y3)*((x2p-y2)**2+(x3+y3) &
            **2)+(-1._8+nu)*((-1)+2*nu)*lr2*y3*(2*x3+y3)*((x2p-y2) &
            **2+(x3+y3)**2))+2*(-1._8+nu)*((-1)+2*nu)*(x2p-y2)*atan3((x2p-y2),(x1p-y1)) &
            +x2p*atan2(-x2p,x1p-y1) &
            -3*x2p*atan2(3*x2p,x1p-y1)+4*nu*x2p*atan2(-nu*x2p, &
            x1p-y1)+3*y2*atan2((-3)*y2,x1p-y1)-y2*atan2( &
            y2,x1p-y1)-4*nu*y2*atan2(nu*y2,x1p-y1)+2*((-1)+ &
            2*nu)*(x2p-y2)*atan2(lr1*(-x2p+y2),(x1p-y1)*(x3+ &
            (-1)*y3))+2*(1-2*nu)**2*(x2p-y2)*atan2(lr2*(-x2p+ &
            y2),(x1p-y1)*(x3+y3))+xLogy((-2)*x3,lr2-x1p+y1)+xLogy(( &
            -1)*((-3)+4*nu)*(x1p-y1),lr1+x3-y3)+xLogy(-(3+( &
            -6)*nu+4*nu**2)*(x1p-y1),lr2+x3+y3)+xLogy(-((-3)+4* &
            nu)*(x3-y3),lr1+x1p-y1)+xLogy(-(5+4*nu*((-3)+2* &
            nu))*(x3+y3),lr2+x1p-y1))

    END FUNCTION J2213

    !---------------------------------------------------------------
    !> function J2223
    !! computes the J integral J2223
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2223(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2223=((1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(2*lr2**(-1)*((x1p+ &
            (-1)*y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(-((-1)+ &
            nu)*((-1)+2*nu)*lr2**2*(x3+y3)+(-1._8+nu)*((-1)+2*nu)*lr2* &
            y3*(2*x3+y3)+x3*((x1p-y1)**2+(x2p-y2)**2+x3*(x3+y3)) &
            )+x1p*atan2(-x1p,x2p-y2)-3*x1p*atan2(3*x1p,x2p- &
            y2)+4*nu*x1p*atan2(-nu*x1p,x2p-y2)-4*(-1._8+nu)*( &
            x1p-y1)*atan2(lr1*(x1p-y1),(x2p-y2)*(x3-y3) &
            )+4*(-1._8+nu)*(x1p-y1)*atan2(lr2*(x1p-y1),(x2p- &
            y2)*(x3+y3))+3*y1*atan2((-3)*y1,x2p-y2)-y1*atan2( &
            y1,x2p-y2)-4*nu*y1*atan2(nu*y1,x2p-y2)+xLogy((-1) &
            *((-3)+4*nu)*(x2p-y2),lr1+x3-y3)+xLogy(-(3 &
            -6*nu+4*nu**2)*(x2p-y2),lr2+x3+y3)+xLogy((-4)*(-1._8+nu)*( &
            x3-y3),lr1+x2p-y2)+xLogy(4*(-1._8+nu)*(x3+y3),lr2+x2p+( &
            -1)*y2))

    END FUNCTION J2223

    !---------------------------------------------------------------
    !> function J3212
    !! computes the J integral J3212
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3212(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3212=(-(1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*((-2)*lr2**(-1)* &
            x3*(x1p-y1)*y3*(x3+y3)*((x2p-y2)**2+(x3+y3)**2)**( &
            -1)+4*(-1._8+nu)*((-1)+2*nu)*(x2p-y2)*atan3((x2p-y2),(x1p-y1)) &
            +4*(-1._8+nu)*((-1)+2*nu)*(x2p-y2)* &
            atan2(lr2*(-x2p+y2),(x1p-y1)*(x3+y3))+xLogy((-4)*((-1)+ &
            nu)*((-1)+2*nu)*(x1p-y1),lr2+x3+y3)+xLogy(x3-y3,lr1+ &
            x1p-y1)+xLogy(-x3-7*y3-8*nu**2*(x3+y3)+8*nu*( &
            x3+2*y3),lr2+x1p-y1))

    END FUNCTION J3212

    !---------------------------------------------------------------
    !> function J3213
    !! computes the J integral J3213
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3213(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3213=(-(1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(2*lr2**(-1)*x3*( &
            x1p-y1)*(x2p-y2)*y3*((x2p-y2)**2+(x3+y3)**2)**( &
            -1)+4*(-1._8+nu)*((-1)+2*nu)*(x3+y3)*atan3((x1p-y1),(x2p-y2)) &
            +4*((-1)+2*nu)*(nu*x3+(-1._8+nu)*y3)*atan2( &
            lr2*(x2p-y2),(x1p-y1)*(x3+y3))+xLogy(x2p-y2,lr1+x1p+ &
            (-1)*y1)+xLogy(-(1+8*(-1._8+nu)*nu)*(x2p-y2),lr2+x1p+( &
            -1)*y1))

    END FUNCTION J3213

    !---------------------------------------------------------------
    !> function J3223
    !! computes the J integral J3223
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3223(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3223=(-(1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(lr1+((-1)-8*(( &
            -1)+nu)*nu)*lr2-2*lr2**(-1)*x3*y3+2*((-3)+4*nu)*x3* &
            ACOTH(lr2**(-1)*(x3+y3))+xLogy(2*(3*x3+2*y3-6*nu*(x3+y3)+ &
            4*nu**2*(x3+y3)),lr2+x3+y3))

    END FUNCTION J3223

    !---------------------------------------------------------------
    !> function J1312
    !! computes the J integral J1312
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1312(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1312=(-(1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(2*lr2**(-1)*x3*( &
            x2p-y2)*y3*(x3+y3)*((x1p-y1)**2+(x3+y3)**2)**(-1)+( &
            -4)*(-1._8+nu)*((-1)+2*nu)*(x1p-y1)*atan3((x1p-y1),(x2p-y2)) &
            +4*(-1._8+nu)*((-1)+2*nu)*(x1p-y1)* &
            atan2(lr2*(x1p-y1),(x2p-y2)*(x3+y3))+xLogy(4*(-1._8+nu) &
            *((-1)+2*nu)*(x2p-y2),lr2+x3+y3)+xLogy(x3-y3,lr1+x2p+( &
            -1)*y2)+xLogy((7+8*((-2)+nu)*nu)*x3+y3+8*(-1._8+nu)*nu*y3, &
            lr2+x2p-y2))

    END FUNCTION J1312

    !---------------------------------------------------------------
    !> function J1313
    !! computes the J integral J1313
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1313(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1313=(-(1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(lr1+lr2**(-1)*((7+ &
            8*((-2)+nu)*nu)*lr2**2+2*x3*y3)+ &
            2*((-3)+4*nu)*x3* &
ACOTH(lr2**(-1)*(x3+y3))+xLogy(2*(-3*x3-2*y3+6*nu*(x3+y3)-4*nu**2*(x3+y3)),lr2+x3+y3))

    END FUNCTION J1313

    !---------------------------------------------------------------
    !> function J1323
    !! computes the J integral J1323
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1323(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1323=(-(1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*((-2)*lr2**(-1)* &
            x3*(x1p-y1)*(x2p-y2)*y3*((x1p-y1)**2+(x3+y3) &
            **2)**(-1)-4*(-1._8+nu)*((-1)+2*nu)*(x3+y3)*atan3((x2p-y2),(x1p-y1)) &
            -4*(-1._8+nu)*((-3)*x3-y3+2* &
            nu*(x3+y3))*atan2(lr2*(x1p-y1),(x2p-y2)*(x3+y3))+ &
            xLogy(x1p-y1,lr1+x2p-y2)+xLogy((7+8*((-2)+nu)*nu)*(x1p+ &
            (-1)*y1),lr2+x2p-y2))

    END FUNCTION J1323

    !---------------------------------------------------------------
    !> function J2312
    !! computes the J integral J2312
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2312(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2312=(-(1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(2*lr2**(-1)*x3*( &
            x1p-y1)*y3*(x3+y3)*((x2p-y2)**2+(x3+y3)**2)**(-1)+( &
            -4)*(-1._8+nu)*((-1)+2*nu)*(x2p-y2)*atan3((x2p-y2),(x1p-y1)) &
            +4*(-1._8+nu)*((-1)+2*nu)*(x2p-y2)* &
            atan2(lr2*(x2p-y2),(x1p-y1)*(x3+y3))+xLogy(4*(-1._8+nu) &
            *((-1)+2*nu)*(x1p-y1),lr2+x3+y3)+xLogy(x3-y3,lr1+x1p+( &
            -1)*y1)+xLogy((7+8*((-2)+nu)*nu)*x3+y3+8*(-1._8+nu)*nu*y3, &
            lr2+x1p-y1))

    END FUNCTION J2312

    !---------------------------------------------------------------
    !> function J2313
    !! computes the J integral J2313
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2313(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2313=(-(1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*((-2)*lr2**(-1)* &
            x3*(x1p-y1)*(x2p-y2)*y3*((x2p-y2)**2+(x3+y3) &
            **2)**(-1)-4*(-1._8+nu)*((-1)+2*nu)*(x3+y3)*atan3((x1p-y1),(x2p-y2)) &
            -4*(-1._8+nu)*((-3)*x3-y3+2* &
            nu*(x3+y3))*atan2(lr2*(x2p-y2),(x1p-y1)*(x3+y3))+ &
            xLogy(x2p-y2,lr1+x1p-y1)+xLogy((7+8*((-2)+nu)*nu)*(x2p+ &
            (-1)*y2),lr2+x1p-y1))

    END FUNCTION J2313

    !---------------------------------------------------------------
    !> function J2323
    !! computes the J integral J2323
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2323(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2323=(-(1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(lr1+lr2**(-1)*((7._8+ &
            8._8*((-2._8)+nu)*nu)*lr2**2+2*x3*y3)+ &
            2._8*((-3._8)+4._8*nu)*x3* &
ACOTH(lr2**(-1)*(x3+y3))+xLogy(2*((-3._8)*x3-2._8*y3+6._8*nu*(x3+y3)-4._8*nu**2*(x3+y3)),lr2+x3+y3))

      END FUNCTION J2323

    !---------------------------------------------------------------
    !> function J3312
    !! computes the J integral J3312
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3312(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3312=((1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(2*lr2**(-1)*x3*( &
            x1p-y1)*(x2p-y2)*y3*((x1p-y1)**2+(x3+y3)**2)**( &
            -1)*((x2p-y2)**2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2+2*(x3+y3)**2)-3*x3*atan2(3*x3,x1p-y1) &
            -5*x3*atan2(5*x3,x2p-y2)+12*nu*x3*atan2((-3)*nu*x3,x2p+( &
            -1)*y2)+4*nu*x3*atan2(-nu*x3,x1p-y1)-8*nu**2* &
            x3*atan2(nu**2*x3,x2p-y2)+3*y3*atan2((-3)*y3,x1p- &
            y1)-5*y3*atan2(5*y3,x2p-y2)+12*nu*y3*atan2((-3)* &
            nu*y3,x2p-y2)-4*nu*y3*atan2(nu*y3,x1p-y1)-8* &
            nu**2*y3*atan2(nu**2*y3,x2p-y2)+2*((-1)+2*nu)*(x3+(-1) &
            *y3)*atan2(lr1*(-x3+y3),(x1p-y1)*(x2p-y2))+2*( &
            1-2*nu)**2*(x3+y3)*atan2(lr2*(x3+y3),(x1p-y1)*(x2p+(-1) &
            *y2))+xLogy(-((-3)+4*nu)*(x1p-y1),lr1+x2p-y2)+ &
            xLogy((5+4*nu*((-3)+2*nu))*(x1p-y1),lr2+x2p-y2)+ &
            xLogy(-((-3)+4*nu)*(x2p-y2),lr1+x1p-y1)+xLogy((5+ &
            4*nu*((-3)+2*nu))*(x2p-y2),lr2+x1p-y1))

    END FUNCTION J3312

    !---------------------------------------------------------------
    !> function J3313
    !! computes the J integral J3313
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3313(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3313=((1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(2*lr2**(-1)*x3*( &
            x1p-y1)*y3*(x3+y3)*((x2p-y2)**2+(x3+y3)**2)**(-1)+5* &
            x2p*atan2((-5)*x2p,x1p-y1)-3*x2p*atan2(3*x2p,x1p-y1) &
            +4*nu*x2p*atan2(-nu*x2p,x1p-y1)-12*nu*x2p*atan2( &
            3*nu*x2p,x1p-y1)+8*nu**2*x2p*atan2(-nu**2*x2p,x1p+(-1) &
            *y1)-4*(-1._8+nu)*(x2p-y2)*atan2(lr1*(x2p-y2),(x1p+ &
            (-1)*y1)*(x3-y3))-8*(-1._8+nu)**2*(x2p-y2)* &
            atan2(lr2*(x2p-y2),(x1p-y1)*(x3+y3))+3*y2*atan2((-3) &
            *y2,x1p-y1)-5*y2*atan2(5*y2,x1p-y1)+12*nu*y2* &
            atan2((-3)*nu*y2,x1p-y1)-4*nu*y2*atan2(nu*y2,x1p+(-1) &
            *y1)-8*nu**2*y2*atan2(nu**2*y2,x1p-y1)+xLogy((-4)* &
            x3,lr2-x1p+y1)+xLogy((-4)*(-1._8+nu)*(x1p-y1),lr1+x3+(-1) &
            *y3)+xLogy((-8)*(-1._8+nu)**2*(x1p-y1),lr2+x3+y3)+xLogy((-1) &
            *((-3)+4*nu)*(x3-y3),lr1+x1p-y1)+xLogy((-7)*x3 &
            -5*y3+12*nu*(x3+y3)-8*nu**2*(x3+y3),lr2+x1p-y1))
          
    END FUNCTION J3313

    !---------------------------------------------------------------
    !> function J3323
    !! computes the J integral J3323
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3323(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3323=((1._8/16._8))*(1-nu)**(-1)*PI**(-1)*G**(-1)*(2*lr2**(-1)*x3*( &
            x2p-y2)*y3*(x3+y3)*((x1p-y1)**2+(x3+y3)**2)**(-1)+5* &
            x1p*atan2((-5)*x1p,x2p-y2)-3*x1p*atan2(3*x1p,x2p-y2) &
            +4*nu*x1p*atan2(-nu*x1p,x2p-y2)-12*nu*x1p*atan2( &
            3*nu*x1p,x2p-y2)+8*nu**2*x1p*atan2(-nu**2*x1p,x2p+(-1) &
            *y2)-4*(-1._8+nu)*(x1p-y1)*atan2(lr1*(x1p-y1),(x2p+ &
            (-1)*y2)*(x3-y3))-8*(-1._8+nu)**2*(x1p-y1)* &
            atan2(lr2*(x1p-y1),(x2p-y2)*(x3+y3))+3*y1*atan2((-3) &
            *y1,x2p-y2)-5*y1*atan2(5*y1,x2p-y2)+12*nu*y1* &
            atan2((-3)*nu*y1,x2p-y2)-4*nu*y1*atan2(nu*y1,x2p+(-1) &
            *y2)-8*nu**2*y1*atan2(nu**2*y1,x2p-y2)+xLogy((-4)* &
            x3,lr2-x2p+y2)+xLogy((-4)*(-1._8+nu)*(x2p-y2),lr1+x3+(-1) &
            *y3)+xLogy((-8)*(-1._8+nu)**2*(x2p-y2),lr2+x3+y3)+xLogy((-1) &
            *((-3)+4*nu)*(x3-y3),lr1+x2p-y2)+xLogy((-7)*x3 &
            -5*y3+12*nu*(x3+y3)-8*nu**2*(x3+y3),lr2+x2p-y2))

    END FUNCTION J3323

  END SUBROUTINE computeDisplacementVerticalStrainVolume

  !------------------------------------------------------------------------
  !> subroutine ComputeStressVerticalStrainVolume computes the stress field 
  !! associated with deforming vertical strain volume using the analytic 
  !! solution of
  !!
  !!   Barbot S., J. D. P. Moore and V. Lambert, Displacement and Stress
  !!   Associated with Distributed Anelastic Deformation in a Half Space,
  !!   Bull. Seism. Soc. Am., 107(2), 10.1785/0120160237, 2017.
  !!
  !! considering the following geometry:
  !!
  !!
  !!                      N (x1p)
  !!                     /
  !!                    /| strike (theta)          E (x2p)
  !!        q1,q2,q3 ->@--------------------------+
  !!                   |                        w |     +
  !!                   |                        i |    /
  !!                   |                        d |   / s
  !!                   |                        t |  / s
  !!                   |                        h | / e
  !!                   |                          |/ n
  !!                   +--------------------------+  k
  !!                   :       l e n g t h       /  c
  !!                   |                        /  i
  !!                   :                       /  h
  !!                   |                      /  t
  !!                   :                     /
  !!                   |                    +
  !!                   Z (x3)
  !!
  !!
  !! INPUT:
  !! @param x1, x2, x3         northing, easting, and depth of the observation point
  !!                           in unprimed system of coordinates.
  !! @param q1, q2, q3         north, east and depth coordinates of the strain volume,
  !! @param L, T, W            length, thickness, and width of the strain volume,
  !! @param theta (degree)     strike of the strain volume,
  !! @param epsijp             anelastic strain component 11, 12, 13, 22, 23 and 33 
  !!                           in the strain volume in the system of reference tied to 
  !!                           the strain volume (primed reference system),
  !! @param G, nu              shear modulus and Poisson's ratio in the half space.
  !!
  !! OUTPUT:
  !! sij                stress components in the unprimed reference system.
  !!
  !! \author James D. P. Moore (10/06/16) - derivatives of J-functions in original form
  !! \author Sylvain Barbot (21/02/17) - original fortran form
  !------------------------------------------------------------------------
  SUBROUTINE computeStressVerticalStrainVolume(x1,x2,x3,q1,q2,q3,L,T,W,theta, &
                                eps11p,eps12p,eps13p,eps22p,eps23p,eps33p,G,nu, &
                                s11,s12,s13,s22,s23,s33)

    REAL*8, INTENT(IN) :: x1,x2,x3,q1,q2,q3,L,T,W,theta
    REAL*8, INTENT(IN) :: eps11p,eps12p,eps13p,eps22p,eps23p,eps33p
    REAL*8, INTENT(IN) :: G,nu
    REAL*8, INTENT(OUT) :: s11,s12,s13,s22,s23,s33
    
    REAL*8 :: lambda
    REAL*8 :: e11,e12,e13,e22,e23,e33,epskk
    REAL*8 :: e11p,e12p,e13p,e22p,e23p
    REAL*8 :: eps11,eps12,eps13,eps22,eps23,eps33
    REAL*8 :: u12,u13,u21,u23,u31,u32
    REAL*8 :: x1p,x2p

    ! check valid parameters
    IF ((-1._8 .GT. nu) .OR. (0.5_8 .LT. nu)) THEN
       WRITE (0,'("error: -1<=nu<=0.5, nu=",ES9.2E2," given.")') nu
       STOP 1
    END IF

    IF (0 .GT. x3) THEN
       WRITE (0,'("error: observation depth (x3) must be positive")')
       STOP 1
    END IF

    IF (0 .GT. q3) THEN
       WRITE (0,'("error: source depth (q3) must be positive")')
       STOP 1
    END IF

    ! lame elastic parameters
    lambda=G*2._8*nu/(1-2*nu)

    ! isotropic strain
    epskk=eps11p+eps22p+eps33p

    ! rotate observation points to the strain-volume-centric system of coordinates
    x1p= (x1-q1)*DCOS(theta)+(x2-q2)*DSIN(theta)
    x2p=-(x1-q1)*DSIN(theta)+(x2-q2)*DCOS(theta)

    ! displacement gradient
    e11p= IU1d1(L,   T/2._8,q3+W)-IU1d1(L,   -T/2._8,q3+W)+IU1d1(L,   -T/2._8,q3)-IU1d1(L,   T/2._8,q3) &
         -IU1d1(0._8,T/2._8,q3+W)+IU1d1(0._8,-T/2._8,q3+W)-IU1d1(0._8,-T/2._8,q3)+IU1d1(0._8,T/2._8,q3)
    u12 = IU1d2(L,   T/2._8,q3+W)-IU1d2(L,   -T/2._8,q3+W)+IU1d2(L,   -T/2._8,q3)-IU1d2(L,   T/2._8,q3) &
         -IU1d2(0._8,T/2._8,q3+W)+IU1d2(0._8,-T/2._8,q3+W)-IU1d2(0._8,-T/2._8,q3)+IU1d2(0._8,T/2._8,q3)
    u13 = IU1d3(L,   T/2._8,q3+W)-IU1d3(L,   -T/2._8,q3+W)+IU1d3(L,   -T/2._8,q3)-IU1d3(L,   T/2._8,q3) &
         -IU1d3(0._8,T/2._8,q3+W)+IU1d3(0._8,-T/2._8,q3+W)-IU1d3(0._8,-T/2._8,q3)+IU1d3(0._8,T/2._8,q3)
    u21 = IU2d1(L,   T/2._8,q3+W)-IU2d1(L,   -T/2._8,q3+W)+IU2d1(L,   -T/2._8,q3)-IU2d1(L,   T/2._8,q3) &
         -IU2d1(0._8,T/2._8,q3+W)+IU2d1(0._8,-T/2._8,q3+W)-IU2d1(0._8,-T/2._8,q3)+IU2d1(0._8,T/2._8,q3)
    e22p= IU2d2(L,   T/2._8,q3+W)-IU2d2(L,   -T/2._8,q3+W)+IU2d2(L,   -T/2._8,q3)-IU2d2(L,   T/2._8,q3) &
         -IU2d2(0._8,T/2._8,q3+W)+IU2d2(0._8,-T/2._8,q3+W)-IU2d2(0._8,-T/2._8,q3)+IU2d2(0._8,T/2._8,q3)
    u23 = IU2d3(L,   T/2._8,q3+W)-IU2d3(L,   -T/2._8,q3+W)+IU2d3(L,   -T/2._8,q3)-IU2d3(L,   T/2._8,q3) &
         -IU2d3(0._8,T/2._8,q3+W)+IU2d3(0._8,-T/2._8,q3+W)-IU2d3(0._8,-T/2._8,q3)+IU2d3(0._8,T/2._8,q3)
    u31 = IU3d1(L,   T/2._8,q3+W)-IU3d1(L,   -T/2._8,q3+W)+IU3d1(L,   -T/2._8,q3)-IU3d1(L,   T/2._8,q3) &
         -IU3d1(0._8,T/2._8,q3+W)+IU3d1(0._8,-T/2._8,q3+W)-IU3d1(0._8,-T/2._8,q3)+IU3d1(0._8,T/2._8,q3)
    u32 = IU3d2(L,   T/2._8,q3+W)-IU3d2(L,   -T/2._8,q3+W)+IU3d2(L,   -T/2._8,q3)-IU3d2(L,   T/2._8,q3) &
         -IU3d2(0._8,T/2._8,q3+W)+IU3d2(0._8,-T/2._8,q3+W)-IU3d2(0._8,-T/2._8,q3)+IU3d2(0._8,T/2._8,q3)
    e33 = IU3d3(L,   T/2._8,q3+W)-IU3d3(L,   -T/2._8,q3+W)+IU3d3(L,   -T/2._8,q3)-IU3d3(L,   T/2._8,q3) &
         -IU3d3(0._8,T/2._8,q3+W)+IU3d3(0._8,-T/2._8,q3+W)-IU3d3(0._8,-T/2._8,q3)+IU3d3(0._8,T/2._8,q3)
         
    ! strain
    e12p=(u12+u21)/2
    e13p=(u13+u31)/2
    e23p=(u23+u32)/2

    ! rotate strain field to reference (unprimed) system of coordinates
    e11=(COS(theta)*e11p-SIN(theta)*e12p)*COS(theta)-(COS(theta)*e12p-SIN(theta)*e22p)*SIN(theta)
    e12=(COS(theta)*e11p-SIN(theta)*e12p)*SIN(theta)+(COS(theta)*e12p-SIN(theta)*e22p)*COS(theta)
    e13= COS(theta)*e13p-SIN(theta)*e23p
    e22=(SIN(theta)*e11p+COS(theta)*e12p)*SIN(theta)+(SIN(theta)*e12p+COS(theta)*e22p)*COS(theta)
    e23= SIN(theta)*e13p+COS(theta)*e23p
 
    ! rotate anelastic strain field to reference (unprimed) system of coordinates
    eps11=(COS(theta)*eps11p-SIN(theta)*eps12p)*COS(theta)-(COS(theta)*eps12p-SIN(theta)*eps22p)*SIN(theta)
    eps12=(COS(theta)*eps11p-SIN(theta)*eps12p)*SIN(theta)+(COS(theta)*eps12p-SIN(theta)*eps22p)*COS(theta)
    eps13= COS(theta)*eps13p-SIN(theta)*eps23p
    eps22=(SIN(theta)*eps11p+COS(theta)*eps12p)*SIN(theta)+(SIN(theta)*eps12p+COS(theta)*eps22p)*COS(theta)
    eps23= SIN(theta)*eps13p+COS(theta)*eps23p
    eps33=eps33p

    ! remove anelastic eigenstrain
    e11=e11-eps11*S(x1p/L)*Omega(x2p/T)*S((x3-q3)/W)
    e12=e12-eps12*S(x1p/L)*Omega(x2p/T)*S((x3-q3)/W)
    e13=e13-eps13*S(x1p/L)*Omega(x2p/T)*S((x3-q3)/W)
    e22=e22-eps22*S(x1p/L)*Omega(x2p/T)*S((x3-q3)/W)
    e23=e23-eps23*S(x1p/L)*Omega(x2p/T)*S((x3-q3)/W)
    e33=e33-eps33*S(x1p/L)*Omega(x2p/T)*S((x3-q3)/W)

    ! stress components
    s11=lambda*(e11+e22+e33)+2*G*e11
    s12=2*G*e12
    s13=2*G*e13
    s22=lambda*(e11+e22+e33)+2*G*e22
    s23=2*G*e23
    s33=lambda*(e11+e22+e33)+2*G*e33

  CONTAINS 

    !------------------------------------------------------------------------
    !> function r1
    !! computes the distance from the source at y1,y2,y3
    !------------------------------------------------------------------------
    REAL*8 FUNCTION r1(y1,y2,y3) 
      REAL*8, INTENT(IN) :: y1,y2,y3
    
      r1=sqrt((x1p-y1)**2+(x2p-y2)**2+(x3-y3)**2)
   
    END FUNCTION r1

    !------------------------------------------------------------------------
    !> function r2
    !! computes the distance from the image at y1,y2,-y3
    !------------------------------------------------------------------------
    REAL*8 FUNCTION r2(y1,y2,y3) 
      REAL*8, INTENT(IN) :: y1,y2,y3

      r2=sqrt((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)

    END FUNCTION r2

    !---------------------------------------------------------------
    !> function IU1d1
    !! computes the derivative of the indefinite integral U1 
    !! with respect to x1p
    !! function IU1d1 is the integrand for displacement gradient u1,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION IU1d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      IU1d1=0
        
      IF ((0 .NE. eps11p) .OR. (0 .NE. epskk)) THEN
         IU1d1=IU1d1+(lambda*epskk+2*G*eps11p)*J1123d1(y1,y2,y3)
      END IF
      IF (0 .NE. eps12p) THEN
         IU1d1=IU1d1+2*G*eps12p*(J1223d1(y1,y2,y3)+J1113d1(y1,y2,y3))
      END IF
      IF (0 .NE. eps13p) THEN
         IU1d1=IU1d1+2*G*eps13p*(J1323d1(y1,y2,y3)+J1112d1(y1,y2,y3))
      END IF
      IF ((0 .NE. eps22p) .OR. (0 .NE. epskk)) THEN
         IU1d1=IU1d1+(lambda*epskk+2*G*eps22p)*J1213d1(y1,y2,y3)
      END IF
      IF (0 .NE. eps23p) THEN
         IU1d1=IU1d1+2*G*eps23p*(J1212d1(y1,y2,y3)+J1313d1(y1,y2,y3))
      END IF
      IF ((0 .NE. eps33p) .OR. (0 .NE. epskk)) THEN
         IU1d1=IU1d1+(lambda*epskk+2*G*eps33p)*J1312d1(y1,y2,y3)
      END IF

    END FUNCTION IU1d1

    !---------------------------------------------------------------
    !> function IU1d2
    !! computes the derivative of the indefinite integral U1 
    !! with respect to x2p
    !! function IU1d2 is the integrand for displacement gradient u1,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION IU1d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      IU1d2=0
        
      IF ((0 .NE. eps11p) .OR. (0 .NE. epskk)) THEN
         IU1d2=IU1d2+(lambda*epskk+2*G*eps11p)*J1123d2(y1,y2,y3)
      END IF
      IF (0 .NE. eps12p) THEN
         IU1d2=IU1d2+2*G*eps12p*(J1223d2(y1,y2,y3)+J1113d2(y1,y2,y3))
      END IF
      IF (0 .NE. eps13p) THEN
         IU1d2=IU1d2+2*G*eps13p*(J1323d2(y1,y2,y3)+J1112d2(y1,y2,y3))
      END IF
      IF ((0 .NE. eps22p) .OR. (0 .NE. epskk)) THEN
         IU1d2=IU1d2+(lambda*epskk+2*G*eps22p)*J1213d2(y1,y2,y3)
      END IF
      IF (0 .NE. eps23p) THEN
         IU1d2=IU1d2+2*G*eps23p*(J1212d2(y1,y2,y3)+J1313d2(y1,y2,y3))
      END IF
      IF ((0 .NE. eps33p) .OR. (0 .NE. epskk)) THEN
         IU1d2=IU1d2+(lambda*epskk+2*G*eps33p)*J1312d2(y1,y2,y3)
      END IF

    END FUNCTION IU1d2

    !---------------------------------------------------------------
    !> function IU1d3
    !! computes the derivative of the indefinite integral U1 
    !! with respect to x3
    !! function IU1d3 is the integrand for displacement gradient u1,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION IU1d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      IU1d3=0
        
      IF ((0 .NE. eps11p) .OR. (0 .NE. epskk)) THEN
         IU1d3=IU1d3+(lambda*epskk+2*G*eps11p)*J1123d3(y1,y2,y3)
      END IF
      IF (0 .NE. eps12p) THEN
         IU1d3=IU1d3+2*G*eps12p*(J1223d3(y1,y2,y3)+J1113d3(y1,y2,y3))
      END IF
      IF (0 .NE. eps13p) THEN
         IU1d3=IU1d3+2*G*eps13p*(J1323d3(y1,y2,y3)+J1112d3(y1,y2,y3))
      END IF
      IF ((0 .NE. eps22p) .OR. (0 .NE. epskk)) THEN
         IU1d3=IU1d3+(lambda*epskk+2*G*eps22p)*J1213d3(y1,y2,y3)
      END IF
      IF (0 .NE. eps23p) THEN
         IU1d3=IU1d3+2*G*eps23p*(J1212d3(y1,y2,y3)+J1313d3(y1,y2,y3))
      END IF
      IF ((0 .NE. eps33p) .OR. (0 .NE. epskk)) THEN
         IU1d3=IU1d3+(lambda*epskk+2*G*eps33p)*J1312d3(y1,y2,y3)
      END IF

    END FUNCTION IU1d3

    !---------------------------------------------------------------
    !> function IU2d1
    !! computes the derivative of the indefinite integral U2 
    !! with respect to x1p
    !! function IU2d1 is the integrand for displacement gradient u2,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION IU2d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      IU2d1=0
        
      IF ((0 .NE. eps11p) .OR. (0 .NE. epskk)) THEN
         IU2d1=IU2d1+(lambda*epskk+2*G*eps11p)*J2123d1(y1,y2,y3)
      END IF
      IF (0 .NE. eps12p) THEN
         IU2d1=IU2d1+2*G*eps12p*(J2223d1(y1,y2,y3)+J2113d1(y1,y2,y3))
      END IF
      IF (0 .NE. eps13p) THEN
         IU2d1=IU2d1+2*G*eps13p*(J2323d1(y1,y2,y3)+J2112d1(y1,y2,y3))
      END IF
      IF ((0 .NE. eps22p) .OR. (0 .NE. epskk)) THEN
         IU2d1=IU2d1+(lambda*epskk+2*G*eps22p)*J2213d1(y1,y2,y3)
      END IF
      IF (0 .NE. eps23p) THEN
         IU2d1=IU2d1+2*G*eps23p*(J2212d1(y1,y2,y3)+J2313d1(y1,y2,y3))
      END IF
      IF ((0 .NE. eps33p) .OR. (0 .NE. epskk)) THEN
         IU2d1=IU2d1+(lambda*epskk+2*G*eps33p)*J2312d1(y1,y2,y3)
      END IF

    END

    !---------------------------------------------------------------
    !> function IU2d2
    !! computes the derivative of the indefinite integral U2 
    !! with respect to x2p
    !! function IU2d2 is the integrand for displacement gradient u2,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION IU2d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      IU2d2=0
        
      IF ((0 .NE. eps11p) .OR. (0 .NE. epskk)) THEN
         IU2d2=IU2d2+(lambda*epskk+2*G*eps11p)*J2123d2(y1,y2,y3)
      END IF
      IF (0 .NE. eps12p) THEN
         IU2d2=IU2d2+2*G*eps12p*(J2223d2(y1,y2,y3)+J2113d2(y1,y2,y3))
      END IF
      IF (0 .NE. eps13p) THEN
         IU2d2=IU2d2+2*G*eps13p*(J2323d2(y1,y2,y3)+J2112d2(y1,y2,y3))
      END IF
      IF ((0 .NE. eps22p) .OR. (0 .NE. epskk)) THEN
         IU2d2=IU2d2+(lambda*epskk+2*G*eps22p)*J2213d2(y1,y2,y3)
      END IF
      IF (0 .NE. eps23p) THEN
         IU2d2=IU2d2+2*G*eps23p*(J2212d2(y1,y2,y3)+J2313d2(y1,y2,y3))
      END IF
      IF ((0 .NE. eps33p) .OR. (0 .NE. epskk)) THEN
         IU2d2=IU2d2+(lambda*epskk+2*G*eps33p)*J2312d2(y1,y2,y3)
      END IF

    END FUNCTION IU2d2

    !---------------------------------------------------------------
    !> function IU2d3
    !! computes the derivative of the indefinite integral U2 
    !! with respect to x3
    !! function IU2d3 is the integrand for displacement gradient u2,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION IU2d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      IU2d3=0
        
      IF ((0 .NE. eps11p) .OR. (0 .NE. epskk)) THEN
         IU2d3=IU2d3+(lambda*epskk+2*G*eps11p)*J2123d3(y1,y2,y3)
      END IF
      IF (0 .NE. eps12p) THEN
         IU2d3=IU2d3+2*G*eps12p*(J2223d3(y1,y2,y3)+J2113d3(y1,y2,y3))
      END IF
      IF (0 .NE. eps13p) THEN
         IU2d3=IU2d3+2*G*eps13p*(J2323d3(y1,y2,y3)+J2112d3(y1,y2,y3))
      END IF
      IF ((0 .NE. eps22p) .OR. (0 .NE. epskk)) THEN
         IU2d3=IU2d3+(lambda*epskk+2*G*eps22p)*J2213d3(y1,y2,y3)
      END IF
      IF (0 .NE. eps23p) THEN
         IU2d3=IU2d3+2*G*eps23p*(J2212d3(y1,y2,y3)+J2313d3(y1,y2,y3))
      END IF
      IF ((0 .NE. eps33p) .OR. (0 .NE. epskk)) THEN
         IU2d3=IU2d3+(lambda*epskk+2*G*eps33p)*J2312d3(y1,y2,y3)
      END IF

    END FUNCTION IU2d3

    !---------------------------------------------------------------
    !> function IU3d1
    !! computes the derivative of the indefinite integral U3
    !! with respect to x1p
    !! function IU3d1 is the integrand for displacement gradient u3,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION IU3d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      IU3d1=0
        
      IF ((0 .NE. eps11p) .OR. (0 .NE. epskk)) THEN
         IU3d1=IU3d1+(lambda*epskk+2*G*eps11p)*J3123d1(y1,y2,y3)
      END IF
      IF (0 .NE. eps12p) THEN
         IU3d1=IU3d1+2*G*eps12p*(J3223d1(y1,y2,y3)+J3113d1(y1,y2,y3))
      END IF
      IF (0 .NE. eps13p) THEN
         IU3d1=IU3d1+2*G*eps13p*(J3323d1(y1,y2,y3)+J3112d1(y1,y2,y3))
      END IF
      IF ((0 .NE. eps22p) .OR. (0 .NE. epskk)) THEN
         IU3d1=IU3d1+(lambda*epskk+2*G*eps22p)*J3213d1(y1,y2,y3)
      END IF
      IF (0 .NE. eps23p) THEN
         IU3d1=IU3d1+2*G*eps23p*(J3212d1(y1,y2,y3)+J3313d1(y1,y2,y3))
      END IF
      IF ((0 .NE. eps33p) .OR. (0 .NE. epskk)) THEN
         IU3d1=IU3d1+(lambda*epskk+2*G*eps33p)*J3312d1(y1,y2,y3)
      END IF

    END FUNCTION IU3d1

    !---------------------------------------------------------------
    !> function IU3d2
    !! computes the derivative of the indefinite integral U3
    !! with respect to x2p
    !! function IU3d2 is the integrand for displacement gradient u3,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION IU3d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      IU3d2=0
        
      IF ((0 .NE. eps11p) .OR. (0 .NE. epskk)) THEN
         IU3d2=IU3d2+(lambda*epskk+2*G*eps11p)*J3123d2(y1,y2,y3)
      END IF
      IF (0 .NE. eps12p) THEN
         IU3d2=IU3d2+2*G*eps12p*(J3223d2(y1,y2,y3)+J3113d2(y1,y2,y3))
      END IF
      IF (0 .NE. eps13p) THEN
         IU3d2=IU3d2+2*G*eps13p*(J3323d2(y1,y2,y3)+J3112d2(y1,y2,y3))
      END IF
      IF ((0 .NE. eps22p) .OR. (0 .NE. epskk)) THEN
         IU3d2=IU3d2+(lambda*epskk+2*G*eps22p)*J3213d2(y1,y2,y3)
      END IF
      IF (0 .NE. eps23p) THEN
         IU3d2=IU3d2+2*G*eps23p*(J3212d2(y1,y2,y3)+J3313d2(y1,y2,y3))
      END IF
      IF ((0 .NE. eps33p) .OR. (0 .NE. epskk)) THEN
         IU3d2=IU3d2+(lambda*epskk+2*G*eps33p)*J3312d2(y1,y2,y3)
      END IF

    END FUNCTION IU3d2

    !---------------------------------------------------------------
    !> function IU3d3
    !! computes the derivative of the indefinite integral U3
    !! with respect to x3
    !! function IU3d3 is the integrand for displacement gradient u3,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION IU3d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      IU3d3=0
        
      IF ((0 .NE. eps11p) .OR. (0 .NE. epskk)) THEN
         IU3d3=IU3d3+(lambda*epskk+2*G*eps11p)*J3123d3(y1,y2,y3)
      END IF
      IF (0 .NE. eps12p) THEN
         IU3d3=IU3d3+2*G*eps12p*(J3223d3(y1,y2,y3)+J3113d3(y1,y2,y3))
      END IF
      IF (0 .NE. eps13p) THEN
         IU3d3=IU3d3+2*G*eps13p*(J3323d3(y1,y2,y3)+J3112d3(y1,y2,y3))
      END IF
      IF ((0 .NE. eps22p) .OR. (0 .NE. epskk)) THEN
         IU3d3=IU3d3+(lambda*epskk+2*G*eps22p)*J3213d3(y1,y2,y3)
      END IF
      IF (0 .NE. eps23p) THEN
         IU3d3=IU3d3+2*G*eps23p*(J3212d3(y1,y2,y3)+J3313d3(y1,y2,y3))
      END IF
      IF ((0 .NE. eps33p) .OR. (0 .NE. epskk)) THEN
         IU3d3=IU3d3+(lambda*epskk+2*G*eps33p)*J3312d3(y1,y2,y3)
      END IF

    END FUNCTION IU3d3

    !---------------------------------------------------------------
    !> function J1112d1
    !! computes the derivative of the J integral J1112,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1112d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1112d1=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(x3**2*(x3**2+(x1p+ &
            (-1)*y1)**2)**(-1)+9._8*x3**2._8*(9._8*x3**2+(x1p-y1)**2)**(-1._8)+ &
            4._8*nu**2*x3**2._8*(nu**2*x3**2+(x1p-y1)**2)**(-1)-4._8*((-1._8) &
            +nu)*(lr1+x1p-y1)**(-1)*(x2p-y2)-4._8*((-1._8)+nu)*lr1**( &
            -1)*(x1p-y1)*(lr1+x1p-y1)**(-1)*(x2p-y2)-4._8*(( &
            -1._8)+nu)*(lr2+x1p-y1)**(-1)*(x2p-y2)-4._8*((-1._8)+nu)* &
            lr2**(-1)*(x1p-y1)*(lr2+x1p-y1)**(-1)*(x2p-y2)+( &
            -1)*((-3._8)+4._8*nu)*lr1**(-1)*(x1p-y1)**2._8*(lr1+x2p-y2)**( &
            -1)+(5._8+4._8*nu*((-3._8)+2._8*nu))*lr2**(-1)*(x1p-y1)**2*(lr2+x2p+( &
            -1)*y2)**(-1)+4._8*((-1._8)+nu)*lr1*(x2p-y2)*((x1p-y1)**2+ &
            (x3-y3)**2)**(-1)*((x2p-y2)**2+(x3-y3)**2)**(-1) &
            *(x3-y3)**2-4._8*((-1._8)+nu)*lr1**(-1)*(x1p-y1)**2*( &
            x2p-y2)*((x1p-y1)**2+(x3-y3)**2)**(-1)*((x2p+(-1) &
            *y2)**2+(x3-y3)**2)**(-1)*(x3-y3)**2+4._8*((-1._8)+nu)*( &
            (-1._8)+2._8*nu)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p- &
            y2)*(x3+y3)+y3**2*((x1p-y1)**2+y3**2)**(-1)+9._8*y3**2._8*((x1p+ &
            (-1)*y1)**2+9._8*y3**2)**(-1)+4._8*nu**2*y3**2*((x1p-y1)**2+ &
            nu**2*y3**2)**(-1)-4._8*lr2**(-1)*x3*(x1p-y1)**2*(x2p+(-1) &
            *y2)*y3*((x1p-y1)**2+(x3+y3)**2)**(-2)+2._8*lr2**(-1)*x3*( &
            x2p-y2)*y3*((x1p-y1)**2+(x3+y3)**2)**(-1)-4._8*((-1._8)+ &
            nu)*((-1._8)+2._8*nu)*lr2*((x1p-y1)**2+(x2p-y2)**2)**(-1)* &
            (x2p-y2)*(x3+y3)**2*((x1p-y1)**2+(x3+y3)**2)**(-1)+( &
            -4._8)*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-1)*(x1p-y1)**2*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(x3+y3)**2*(( &
            x1p-y1)**2+(x3+y3)**2)**(-1)+4._8*((-1._8)+nu)*lr2*(x2p-y2) &
            *(x3+y3)**2*((x1p-y1)**2+(x3+y3)**2)**(-1)*((x2p-y2) &
            **2+(x3+y3)**2)**(-1._8)-4._8*((-1._8)+nu)*lr2**(-1)*(x1p-y1) &
            **2*(x2p-y2)*(x3+y3)**2*((x1p-y1)**2+(x3+y3)**2)**( &
            -1)*((x2p-y2)**2+(x3+y3)**2)**(-1)-2*x3*(x1p-y1) &
            **2*(x2p-y2)*y3*((x1p-y1)**2+(x3+y3)**2)**(-1)*((x1p+ &
            (-1)*y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)+xLogy(3._8-4._8* &
            nu,lr1+x2p-y2)+xLogy(5._8+4._8*nu*((-3._8)+2._8*nu),lr2+x2p-y2))

    END FUNCTION J1112d1

    !---------------------------------------------------------------
    !> function J1112d2
    !! computes the derivative of the J integral J1112,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1112d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1112d2=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*((-4._8)*((-1._8)+nu)* &
            lr1**(-1)*(lr1+x1p-y1)**(-1)*(x2p-y2)**2-4._8*((-1._8)+nu) &
            *lr2**(-1)*(lr2+x1p-y1)**(-1)*(x2p-y2)**2-((-3._8)+ &
            4._8*nu)*(x1p-y1)*(lr1+x2p-y2)**(-1)-((-3._8)+4._8*nu)* &
            lr1**(-1)*(x1p-y1)*(x2p-y2)*(lr1+x2p-y2)**(-1)+(5+ &
            4._8*nu*((-3._8)+2._8*nu))*(x1p-y1)*(lr2+x2p-y2)**(-1)+(5._8+4._8* &
            nu*((-3._8)+2._8*nu))*lr2**(-1)*(x1p-y1)*(x2p-y2)*(lr2+x2p+ &
            (-1)*y2)**(-1)+4._8*((-1._8)+nu)*lr1*(x1p-y1)*((x1p-y1) &
            **2+(x3-y3)**2)**(-1)*((x2p-y2)**2+(x3-y3)**2) &
            **(-1)*(x3-y3)**2-4._8*((-1._8)+nu)*lr1**(-1)*(x1p-y1) &
            *(x2p-y2)**2*((x1p-y1)**2+(x3-y3)**2)**(-1)*(( &
            x2p-y2)**2+(x3-y3)**2)**(-1)*(x3-y3)**2-4._8*(( &
            -1._8)+nu)*((-1._8)+2._8*nu)*(x1p-y1)*((x1p-y1)**2+(x2p- &
            y2)**2)**(-1)*(x3+y3)+2*lr2**(-1)*x3*(x1p-y1)*y3*((x1p+( &
            -1)*y1)**2+(x3+y3)**2)**(-1)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2*(x1p+ &
            (-1)*y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x3+y3)**2* &
            ((x1p-y1)**2+(x3+y3)**2)**(-1)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu) &
            *lr2**(-1)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**( &
            -1)*(x2p-y2)**2*(x3+y3)**2*((x1p-y1)**2+(x3+y3)**2) &
            **(-1)+4._8*((-1._8)+nu)*lr2*(x1p-y1)*(x3+y3)**2*((x1p-y1) &
            **2+(x3+y3)**2)**(-1)*((x2p-y2)**2+(x3+y3)**2)**(-1)-4._8*( &
            (-1._8)+nu)*lr2**(-1)*(x1p-y1)*(x2p-y2)**2*(x3+y3)**2*( &
            (x1p-y1)**2+(x3+y3)**2)**(-1)*((x2p-y2)**2+(x3+y3)**2) &
            **(-1)-2._8*x3*(x1p-y1)*(x2p-y2)**2*y3*((x1p- &
            y1)**2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+ &
            y3)**2)**(-3._8/2._8)+xLogy(4._8-4._8*nu,lr1+x1p-y1)+xLogy(4._8-4._8*nu, &
            lr2+x1p-y1))

    END FUNCTION J1112d2

    !---------------------------------------------------------------
    !> function J1112d3
    !! computes the derivative of the J integral J1112,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1112d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1112d3=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(x3*(x3**2+(x1p+( &
            -1)*y1)**2)**(-1)*(-x1p+y1)+9._8*x3*(9._8*x3**2+(x1p-y1) &
            **2)**(-1)*(-x1p+y1)+4._8*nu**2*x3*(nu**2*x3**2+(x1p- &
            y1)**2)**(-1)*(-x1p+y1)-4._8*((-1._8)+nu)*lr1**(-1)*(lr1+x1p+( &
            -1)*y1)**(-1)*(x2p-y2)*(x3-y3)-((-3._8)+4._8*nu)* &
            lr1**(-1)*(x1p-y1)*(lr1+x2p-y2)**(-1)*(x3-y3)+( &
            -4._8)*((-1._8)+nu)*lr1*(x1p-y1)*(x2p-y2)*((x1p-y1) &
            **2+(x3-y3)**2)**(-1)*((x2p-y2)**2+(x3-y3)**2) &
            **(-1)*(x3-y3)-4._8*((-1._8)+nu)*lr1**(-1)*(x1p-y1)*( &
            x2p-y2)*((x1p-y1)**2+(x3-y3)**2)**(-1)*((x2p+(-1) &
            *y2)**2+(x3-y3)**2)**(-1)*(x3-y3)**3-4._8*((-1._8)+nu) &
            *lr2**(-1)*(lr2+x1p-y1)**(-1)*(x2p-y2)*(x3+y3)+(5._8+4._8* &
            nu*((-3._8)+2._8*nu))*lr2**(-1)*(x1p-y1)*(lr2+x2p-y2)**(-1) &
            *(x3+y3)-4._8*lr2**(-1)*x3*(x1p-y1)*(x2p-y2)*y3*( &
            x3+y3)*((x1p-y1)**2+(x3+y3)**2)**(-2)+2._8*lr2**(-1)*(x1p+(-1) &
            *y1)*(x2p-y2)*y3*((x1p-y1)**2+(x3+y3)**2)**(-1)+4._8*( &
            (-1._8)+nu)*((-1._8)+2._8*nu)*lr2*(x1p-y1)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*(x2p-y2)*(x3+y3)*((x1p-y1)**2+(x3+ &
            y3)**2)**(-1)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-1)*(x1p- &
            y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*( &
            x3+y3)**3*((x1p-y1)**2+(x3+y3)**2)**(-1)-4._8*((-1._8)+nu)* &
            lr2*(x1p-y1)*(x2p-y2)*(x3+y3)*((x1p-y1)**2+(x3+ &
            y3)**2)**(-1)*((x2p-y2)**2+(x3+y3)**2)**(-1)-4._8*((-1._8)+nu) &
            *lr2**(-1)*(x1p-y1)*(x2p-y2)*(x3+y3)**3*((x1p- &
            y1)**2+(x3+y3)**2)**(-1)*((x2p-y2)**2+(x3+y3)**2)**(-1) &
            -2._8*x3*(x1p-y1)*(x2p-y2)*y3*(x3+y3)*((x1p-y1) &
            **2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3) &
            **2)**(-3._8/2._8)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*atan((x1p-y1)**( &
            -1)*(x2p-y2))-atan2(x3,x1p-y1)-3._8*atan2(3._8*x3, &
            x1p-y1)+4._8*nu*atan2(-nu*x3,x1p-y1)+4._8*((-1._8)+nu)* &
            ((-1._8)+2._8*nu)*atan2(lr2*(-x1p+y1),(x2p-y2)*(x3+y3))+(4._8+( &
            -4._8)*nu)*atan2(lr1*(x3-y3),(x1p-y1)*(x2p-y2))+( &
            4._8-4._8*nu)*atan2(lr2*(x3+y3),(x1p-y1)*(x2p-y2)))

    END FUNCTION J1112d3

    !---------------------------------------------------------------
    !> function J1113d1
    !! computes the derivative of the J integral J1113,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1113d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1113d1=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(x2p**2*(x2p**2+(x1p+ &
            (-1)*y1)**2)**(-1)+9._8*x2p**2*(9._8*x2p**2+(x1p-y1)**2)**(-1)+ &
            4._8*nu**2*x2p**2*(nu**2*x2p**2+(x1p-y1)**2)**(-1)+y2**2*(( &
            x1p-y1)**2+y2**2)**(-1)+9._8*y2**2*((x1p-y1)**2+9._8*y2**2) &
            **(-1)+4._8*nu**2*y2**2*((x1p-y1)**2+nu**2*y2**2)**(-1) &
            -4._8*((-1._8)+nu)*(lr1+x1p-y1)**(-1)*(x3-y3)-4._8*((-1._8)+nu) &
            *lr1**(-1)*(x1p-y1)*(lr1+x1p-y1)**(-1)*(x3-y3)+ &
            4._8*((-1._8)+nu)*lr1*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+( &
            -1)*y2)**2*((x2p-y2)**2+(x3-y3)**2)**(-1)*(x3- &
            y3)-4._8*((-1._8)+nu)*lr1**(-1)*(x1p-y1)**2*((x1p-y1) &
            **2+(x2p-y2)**2)**(-1)*(x2p-y2)**2*((x2p-y2)**2+( &
            x3-y3)**2)**(-1)*(x3-y3)-((-3._8)+4._8*nu)*lr1**(-1) &
            *(x1p-y1)**2*(lr1+x3-y3)**(-1)+4._8*((-1._8)+nu)*(lr2+x1p+( &
            -1)*y1)**(-1)*(x3+y3)+4._8*((-1._8)+nu)*lr2**(-1)*(x1p-y1)*( &
            lr2+x1p-y1)**(-1)*(x3+y3)-(3-6*nu+4*nu**2)*lr2**( &
            -1)*(x1p-y1)**2*(lr2+x3+y3)**(-1)-4._8*((-1._8)+nu)*((-1._8)+2._8* &
            nu)*(x1p-y1)**2*((x1p-y1)**2+(x2p-y2)**2)**(-2)* &
            y3*(2*x3+y3)+2._8*((-1._8)+nu)*((-1._8)+2._8*nu)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*y3*(2*x3+y3)+2._8*((-1._8)+nu)*((-1._8)+2._8*nu)* &
            lr2**(-2)*(x1p-y1)**2*((x1p-y1)**2+(x2p-y2)**2)**( &
            -1)*y3*(2*x3+y3)-4._8*((-1._8)+nu)*lr2*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*(x2p-y2)**2*(x3+y3)*((x2p-y2)**2+( &
            x3+y3)**2)**(-1)+4._8*((-1._8)+nu)*lr2**(-1)*(x1p-y1)**2*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)**2*(x3+y3)*(( &
            x2p-y2)**2+(x3+y3)**2)**(-1)+2*(x1p-y1)**2*((x1p- &
            y1)**2+(x2p-y2)**2)**(-1)*((x1p-y1)**2+(x2p-y2) &
            **2+(x3+y3)**2)**(-3._8/2._8)*((-3)*nu*(x3*(x3**2+(x1p-y1)**2+( &
            x2p-y2)**2)+(x3*((-2)*lr2+3*x3)+(x1p-y1)**2+(x2p- &
            y2)**2)*y3-(lr2-3*x3)*y3**2+y3**3)+2*nu**2*(x3*( &
            x3**2+(x1p-y1)**2+(x2p-y2)**2)+(x3*((-2)*lr2+3*x3)+( &
            x1p-y1)**2+(x2p-y2)**2)*y3-(lr2-3*x3)*y3**2+ &
            y3**3)+y3*((x1p-y1)**2+(x2p-y2)**2-(lr2-x3+( &
            -1)*y3)*(2*x3+y3)))+4*lr2**(-1)*(x1p-y1)**2*((x1p- &
            y1)**2+(x2p-y2)**2)**(-2)*((-3._8)*nu*(x3+y3)*((x1p-y1) &
            **2+(x2p-y2)**2+(x3+y3)**2)+2._8*nu**2*(x3+y3)*((x1p-y1) &
            **2+(x2p-y2)**2+(x3+y3)**2)+y3*((x1p-y1)**2+(x2p- &
            y2)**2+(x3+y3)*(2._8*x3+y3)))+lr2**(-1)*((x1p-y1)**2+(x2p+(-1) &
            *y2)**2)**(-1)*(6._8*nu*(x3+y3)*(3._8*(x1p-y1)**2+(x2p- &
            y2)**2+(x3+y3)**2)-4._8*nu**2*(x3+y3)*(3._8*(x1p-y1)**2+(x2p+ &
            (-1)*y2)**2+(x3+y3)**2)-2._8*y3*(3._8*(x1p-y1)**2+(x2p- &
            y2)**2+(x3+y3)*(2._8*x3+y3)))+xLogy(3._8-4._8*nu,lr1+x3-y3)+ &
            xLogy((-3._8)+6._8*nu-4._8*nu**2,lr2+x3+y3))

    END FUNCTION J1113d1

    !---------------------------------------------------------------
    !> function J1113d2
    !! computes the derivative of the J integral J1113,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1113d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1113d2=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(x2p*(x2p**2+(x1p+( &
            -1)*y1)**2)**(-1)*(-x1p+y1)+9._8*x2p*(9._8*x2p**2+(x1p-y1) &
            **2)**(-1)*(-x1p+y1)+4._8*nu**2*x2p*(nu**2*x2p**2+(x1p- &
            y1)**2)**(-1)*(-x1p+y1)-4._8*((-1._8)+nu)*lr1**(-1)*(lr1+x1p+( &
            -1)*y1)**(-1)*(x2p-y2)*(x3-y3)-4._8*((-1._8)+nu)*lr1* &
            (x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+(-1) &
            *y2)*((x2p-y2)**2+(x3-y3)**2)**(-1)*(x3-y3)+( &
            -4._8)*((-1._8)+nu)*lr1**(-1)*(x1p-y1)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*(x2p-y2)**3*((x2p-y2)**2+(x3- &
            y3)**2)**(-1)*(x3-y3)-((-3)+4._8*nu)*lr1**(-1)*(x1p+( &
            -1)*y1)*(x2p-y2)*(lr1+x3-y3)**(-1)+4._8*((-1._8)+nu)* &
            lr2**(-1)*(lr2+x1p-y1)**(-1)*(x2p-y2)*(x3+y3)-(3+ &
            (-6._8)*nu+4._8*nu**2)*lr2**(-1)*(x1p-y1)*(x2p-y2)*(lr2+ &
            x3+y3)**(-1)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*(x1p-y1)*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2)**(-2)*(x2p-y2)*y3*(2*x3+y3) &
            +2._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-2)*(x1p-y1)*((x1p- &
            y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*y3*(2*x3+y3) &
            -4*lr2**(-1)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**( &
            -1)*(x2p-y2)*(y3-3*nu*(x3+y3)+2*nu**2*(x3+y3))+4._8*(( &
            -1)+nu)*lr2*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**( &
            -1)*(x2p-y2)*(x3+y3)*((x2p-y2)**2+(x3+y3)**2)**(-1)+ &
            4._8*((-1._8)+nu)*lr2**(-1)*(x1p-y1)*((x1p-y1)**2+(x2p+(-1) &
            *y2)**2)**(-1)*(x2p-y2)**3*(x3+y3)*((x2p-y2)**2+(x3+ &
            y3)**2)**(-1)+2*(x1p-y1)*((x1p-y1)**2+(x2p-y2) &
            **2)**(-1)*(x2p-y2)*((x1p-y1)**2+(x2p-y2)**2+(x3+ &
            y3)**2)**(-3._8/2._8)*((-3)*nu*(x3*(x3**2+(x1p-y1)**2+(x2p+(-1) &
            *y2)**2)+(x3*((-2)*lr2+3*x3)+(x1p-y1)**2+(x2p-y2)**2) &
            *y3-(lr2-3*x3)*y3**2+y3**3)+2*nu**2*(x3*(x3**2+(x1p+( &
            -1._8)*y1)**2+(x2p-y2)**2)+(x3*((-2._8)*lr2+3*x3)+(x1p-y1) &
            **2+(x2p-y2)**2)*y3-(lr2-3*x3)*y3**2+y3**3)+y3*(( &
            x1p-y1)**2+(x2p-y2)**2-(lr2-x3-y3)*( &
            2*x3+y3)))+4*lr2**(-1)*(x1p-y1)*((x1p-y1)**2+(x2p+(-1) &
            *y2)**2)**(-2)*(x2p-y2)*((-3._8)*nu*(x3+y3)*((x1p-y1) &
            **2+(x2p-y2)**2+(x3+y3)**2)+2._8*nu**2*(x3+y3)*((x1p-y1) &
            **2+(x2p-y2)**2+(x3+y3)**2)+y3*((x1p-y1)**2+(x2p- &
            y2)**2+(x3+y3)*(2._8*x3+y3)))+atan2(x1p-y1,(-1)*x2p)-3._8* &
            atan2(3*x2p,x1p-y1)+4._8*nu*atan2(-nu*x2p,x1p-y1)+( &
            4._8-4._8*nu)*atan2(lr1*(x2p-y2),(x1p-y1)*(x3-y3)) &
            +4._8*((-1._8)+nu)*atan2(lr2*(x2p-y2),(x1p-y1)*(x3+y3)))

    END FUNCTION J1113d2

    !---------------------------------------------------------------
    !> function J1113d3
    !! computes the derivative of the J integral J1113,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1113d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1113d3=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(4._8*((-1._8)+nu)*lr1* &
            (x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+(-1) &
            *y2)**2*((x2p-y2)**2+(x3-y3)**2)**(-1)-4._8*((-1._8)+ &
            nu)*lr1**(-1)*(lr1+x1p-y1)**(-1)*(x3-y3)**2-4._8*(( &
            -1)+nu)*lr1**(-1)*(x1p-y1)*((x1p-y1)**2+(x2p-y2) &
            **2)**(-1)*(x2p-y2)**2*((x2p-y2)**2+(x3-y3)**2) &
            **(-1)*(x3-y3)**2-((-3._8)+4._8*nu)*(x1p-y1)*(lr1+ &
            x3-y3)**(-1)-((-3._8)+4._8*nu)*lr1**(-1)*(x1p-y1)*( &
            x3-y3)*(lr1+x3-y3)**(-1)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*( &
            x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*y3+4._8*(( &
            -1._8)+nu)*lr2**(-1)*(lr2+x1p-y1)**(-1)*(x3+y3)**2-(3+( &
            -6._8)*nu+4._8*nu**2)*(x1p-y1)*(lr2+x3+y3)**(-1)-(3-6* &
            nu+4*nu**2)*lr2**(-1)*(x1p-y1)*(x3+y3)*(lr2+x3+y3)**(-1)+ &
            2._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-2)*(x1p-y1)*((x1p- &
            y1)**2+(x2p-y2)**2)**(-1)*y3*(x3+y3)*(2*x3+y3)-2* &
            lr2**(-1)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1) &
            *(3*nu*((-3)+2*nu)*x3**2+nu*((-3)+2*nu)*((x1p-y1)**2+ &
            (x2p-y2)**2)+2*(2._8-9._8*nu+6._8*nu**2)*x3*y3+(3._8-9._8*nu+6._8* &
            nu**2)*y3**2)-4._8*((-1._8)+nu)*lr2*(x1p-y1)*((x1p-y1) &
            **2+(x2p-y2)**2)**(-1)*(x2p-y2)**2*((x2p-y2)**2+( &
            x3+y3)**2)**(-1)+4._8*((-1._8)+nu)*lr2**(-1)*(x1p-y1)*((x1p+(-1) &
            *y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)**2*(x3+y3)**2*(( &
            x2p-y2)**2+(x3+y3)**2)**(-1)+2*(x1p-y1)*((x1p-y1) &
            **2+(x2p-y2)**2)**(-1)*(x3+y3)*((x1p-y1)**2+(x2p- &
            y2)**2+(x3+y3)**2)**(-3._8/2._8)*((-3)*nu*(x3*(x3**2+(x1p-y1) &
            **2+(x2p-y2)**2)+(x3*((-2)*lr2+3*x3)+(x1p-y1)**2+(x2p+( &
            -1)*y2)**2)*y3-(lr2-3._8*x3)*y3**2+y3**3)+2._8*nu**2*(x3* &
            (x3**2+(x1p-y1)**2+(x2p-y2)**2)+(x3*((-2)*lr2+3*x3)+( &
            x1p-y1)**2+(x2p-y2)**2)*y3-(lr2-3*x3)*y3**2+ &
            y3**3)+y3*((x1p-y1)**2+(x2p-y2)**2-(lr2-x3+( &
            -1)*y3)*(2*x3+y3)))+xLogy(4._8-4._8*nu,lr1+x1p-y1)+xLogy(4._8*( &
            (-1._8)+nu),lr2+x1p-y1))

    END FUNCTION J1113d3

    !---------------------------------------------------------------
    !> function J1123d1
    !! computes the derivative of the J integral J1123,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1123d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1123d1=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(2._8*((-1._8)+nu)*(( &
            -1._8)+2._8*nu)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**( &
            -1)*(x2p-y2)+x1p*(x1p**2+(x2p-y2)**2)**(-1)*(-x2p+ &
            y2)+9._8*x1p*(9._8*x1p**2+(x2p-y2)**2)**(-1)*(-x2p+y2)+4._8* &
            nu**2*x1p*(nu**2*x1p**2+(x2p-y2)**2)**(-1)*(-x2p+y2)+ &
            2._8*lr2**(-1)*x3*(-x1p+y1)*(lr2-x2p+y2)**(-1)-((-3._8) &
            +4._8*nu)*lr1**(-1)*(x1p-y1)*(lr1+x2p-y2)**(-1)*(x3+(-1) &
            *y3)-2*((-1._8)+2._8*nu)*lr1*(x1p-y1)*((x1p-y1)**2+( &
            x2p-y2)**2)**(-1)*(x2p-y2)*((x1p-y1)**2+(x3+(-1) &
            *y3)**2)**(-1)*(x3-y3)-2*((-1._8)+2._8*nu)*lr1**(-1)*(x1p+( &
            -1)*y1)**3*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p- &
            y2)*((x1p-y1)**2+(x3-y3)**2)**(-1)*(x3-y3)+(-1) &
            *((-3._8)+4._8*nu)*lr1**(-1)*(x1p-y1)*(x2p-y2)*(lr1+x3+( &
            -1)*y3)**(-1)-(5._8+4._8*nu*((-3._8)+2*nu))*lr2**(-1)*(x1p- &
            y1)*(lr2+x2p-y2)**(-1)*(x3+y3)-(3._8-6._8*nu+4._8*nu**2)* &
            lr2**(-1)*(x1p-y1)*(x2p-y2)*(lr2+x3+y3)**(-1)+4._8*((-1._8)+ &
            nu)*((-1._8)+2._8*nu)*(x1p-y1)*((x1p-y1)**2+(x2p-y2) &
            **2)**(-2)*(x2p-y2)*y3*(2*x3+y3)-2._8*((-1._8)+nu)*((-1._8)+ &
            2._8*nu)*lr2**(-2)*(x1p-y1)*((x1p-y1)**2+(x2p-y2) &
            **2)**(-1)*(x2p-y2)*y3*(2._8*x3+y3)+4._8*((-1._8)+nu)*((-1._8)+2._8* &
            nu)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*( &
            x2p-y2)*y3*(2*x3+y3)*((x1p-y1)**2+(x3+y3)**2)**(-1)+ &
            (-2._8)*((-1._8)+2._8*nu)*(x1p-y1)*((x1p-y1)**2+(x2p-y2) &
            **2)**(-1)*(x2p-y2)*((x1p-y1)**2+(x3+y3)**2)**(-1)*(( &
            (-1._8)+2._8*nu)*lr2*(x3+y3)+2._8*((-1._8)+nu)*y3*(2._8*x3+y3))+2._8*lr2**(-1) &
            *(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+( &
            -1)*y2)*((x1p-y1)**2+(x3+y3)**2)**(-1)*(-x1p**2*x3+ &
            2._8*x1p*x3*y1-x3*y1**2+3*x1p**2*y3+2._8*x2p**2*y3+8*x3**2* &
            y3-6._8*x1p*y1*y3+3*y1**2*y3-4._8*x2p*y2*y3+2*y2**2*y3+ &
            12._8*x3*y3**2+4._8*y3**3+4._8*nu**2*(x3+y3)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2+2*(x3+y3)**2)-2._8*nu*(x3+y3)*(4._8*(x1p-y1) &
            **2+3._8*((x2p-y2)**2+2*(x3+y3)**2)))-4._8*lr2**(-1)*(x1p+(-1) &
            *y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*( &
            (x1p-y1)**2+(x3+y3)**2)**(-2)*(x1p**4*y3-4._8*x1p**3*y1* &
            y3-3*nu*(x3+y3)*((x1p-y1)**2+(x3+y3)**2)*((x1p- &
            y1)**2+(x2p-y2)**2+(x3+y3)**2)+2*nu**2*(x3+y3)*((x1p- &
            y1)**2+(x3+y3)**2)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2) &
            -2*x1p*y1*y3*((x2p-y2)**2+2._8*(y1**2+(x3+y3)*(2*x3+y3) &
            ))+y3*(2*x3**4+7*x3**3*y3+(y1**2+y3**2)*(y1**2+y2**2+y3**2)+ &
            x3*y3*(6*y1**2+3*y2**2+5._8*y3**2)+x3**2*(4._8*y1**2+2*y2**2+9._8* &
            y3**2)+x2p**2*(y1**2+(x3+y3)*(2*x3+y3))-2*x2p*y2*(y1**2+( &
            x3+y3)*(2*x3+y3)))+x1p**2*y3*((x2p-y2)**2+2*(3*y1**2+( &
            x3+y3)*(2*x3+y3))))-4._8*lr2**(-1)*(x1p-y1)*((x1p- &
            y1)**2+(x2p-y2)**2)**(-2)*(x2p-y2)*((x1p-y1)**2+( &
            x3+y3)**2)**(-1)*(x1p**4*y3-4._8*x1p**3*y1*y3-3._8*nu*(x3+y3) &
            *((x1p-y1)**2+(x3+y3)**2)*((x1p-y1)**2+(x2p-y2) &
            **2+(x3+y3)**2)+2._8*nu**2*(x3+y3)*((x1p-y1)**2+(x3+y3)**2)* &
            ((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)-2*x1p*y1*y3*( &
            (x2p-y2)**2+2*(y1**2+(x3+y3)*(2*x3+y3)))+y3*(2*x3**4+7._8* &
            x3**3*y3+(y1**2+y3**2)*(y1**2+y2**2+y3**2)+x3*y3*(6._8*y1**2+3._8* &
            y2**2+5._8*y3**2)+x3**2*(4._8*y1**2+2._8*y2**2+9._8*y3**2)+x2p**2*(y1**2+ &
            (x3+y3)*(2._8*x3+y3))-2._8*x2p*y2*(y1**2+(x3+y3)*(2*x3+y3)))+ &
            x1p**2*y3*((x2p-y2)**2+2._8*(3._8*y1**2+(x3+y3)*(2._8*x3+y3))))+( &
            -2)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*( &
            x2p-y2)*((x1p-y1)**2+(x3+y3)**2)**(-1)*((x1p-y1) &
            **2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)*(x1p**4*(y3-3._8*nu*( &
            x3+y3)+2._8*nu**2*(x3+y3))-4._8*x1p**3*y1*(y3-3._8*nu*(x3+y3)+ &
            2*nu**2*(x3+y3))-3._8*nu*(y1**2+(x3+y3)**2)*(x3**3+3*x3**2* &
            y3+y3*(y1**2+y2**2-lr2*y3+y3**2)+x3*(y1**2+y2**2-2* &
            lr2*y3+3*y3**2))+2._8*nu**2*(y1**2+(x3+y3)**2)*(x3**3+3*x3**2* &
            y3+y3*(y1**2+y2**2-lr2*y3+y3**2)+x3*(y1**2+y2**2-2* &
            lr2*y3+3*y3**2))+y3*(2._8*x3**4+7._8*x3**3*y3+(y1**2+y3**2)*( &
            y1**2+y2**2+y3**2)+x3*y3*(6._8*y1**2+3._8*y2**2+5._8*y3**2)+x3**2*( &
            4._8*y1**2+2._8*y2**2+9._8*y3**2)-lr2*(2*x3+y3)*(y1**2+(x3+y3) &
            **2))+2._8*x2p*y2*(3._8*nu*(x3+y3)*(y1**2+(x3+y3)**2)-2._8*nu**2* &
            (x3+y3)*(y1**2+(x3+y3)**2)-y3*(y1**2+(x3+y3)*(2*x3+y3))) &
            +x2p**2*((-3._8)*nu*(x3+y3)*(y1**2+(x3+y3)**2)+2._8*nu**2*(x3+y3)* &
            (y1**2+(x3+y3)**2)+y3*(y1**2+(x3+y3)*(2._8*x3+y3)))+2._8*x1p*y1*((( &
            -1._8)+nu)*((-1._8)+2._8*nu)*lr2*y3*(2*x3+y3)+x2p**2*(-y3+3._8*nu* &
            (x3+y3)-2._8*nu**2*(x3+y3))+2._8*x2p*y2*(y3-3._8*nu*(x3+y3)+2._8* &
            nu**2*(x3+y3))+3._8*nu*(x3+y3)*(2._8*y1**2+y2**2+2._8*(x3+y3)**2)+( &
            -2._8)*nu**2*(x3+y3)*(2._8*y1**2+y2**2+2._8*(x3+y3)**2)-y3*(2* &
            y1**2+y2**2+2*(x3+y3)*(2._8*x3+y3)))+x1p**2*(-((-1._8)+nu)*(( &
            -1._8)+2._8*nu)*lr2*y3*(2._8*x3+y3)+2._8*x2p*y2*(-y3+3._8*nu*(x3+y3) &
            -2._8*nu**2*(x3+y3))+x2p**2*(y3-3._8*nu*(x3+y3)+2._8*nu**2*(x3+ &
            y3))-3._8*nu*(x3+y3)*(6._8*y1**2+y2**2+2._8*(x3+y3)**2)+2._8*nu**2*( &
            x3+y3)*(6._8*y1**2+y2**2+2._8*(x3+y3)**2)+y3*(6._8*y1**2+y2**2+2._8*(x3+ &
            y3)*(2._8*x3+y3))))+2._8*((-1._8)+nu)*((-1._8)+2._8*nu)*atan((x1p-y1) &
            *(x2p-y2)**(-1))+atan2(x2p-y2,(-1)*x1p)-3._8*atan2(3._8* &
            x1p,x2p-y2)+4._8*nu*atan2(-nu*x1p,x2p-y2)+((-2._8)+4._8* &
            nu)*atan2(lr1*(-x1p+y1),(x2p-y2)*(x3-y3))+2._8*(1._8+( &
            -2._8)*nu)**2*atan2(lr2*(-x1p+y1),(x2p-y2)*(x3+y3)))

    END FUNCTION J1123d1

    !---------------------------------------------------------------
    !> function J1123d2
    !! computes the derivative of the J integral J1123,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1123d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1123d2=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(x1p**2*(x1p**2+(x2p+ &
            (-1)*y2)**2)**(-1)+9._8*x1p**2*(9._8*x1p**2+(x2p-y2)**2)**(-1)+ &
            4._8*nu**2*x1p**2*(nu**2*x1p**2+(x2p-y2)**2)**(-1)-2._8*((-1._8) &
            +nu)*((-1._8)+2._8*nu)*(x1p-y1)**2*((x1p-y1)**2+(x2p- &
            y2)**2)**(-1)+y1**2*(y1**2+(x2p-y2)**2)**(-1)+9._8*y1**2*(9._8* &
            y1**2+(x2p-y2)**2)**(-1)+4*nu**2*y1**2*(nu**2*y1**2+(x2p+( &
            -1)*y2)**2)**(-1)+2*x3*(lr2-x2p+y2)**(-1)+2*lr2**(-1)*x3* &
            (-x2p+y2)*(lr2-x2p+y2)**(-1)-((-3._8)+4._8*nu)*(lr1+x2p+( &
            -1)*y2)**(-1)*(x3-y3)-((-3._8)+4._8*nu)*lr1**(-1)*(x2p+( &
            -1)*y2)*(lr1+x2p-y2)**(-1)*(x3-y3)+2._8*((-1._8)+2._8*nu)* &
            lr1*(x1p-y1)**2*((x1p-y1)**2+(x2p-y2)**2)**(-1)*( &
            (x1p-y1)**2+(x3-y3)**2)**(-1)*(x3-y3)-2*(( &
            -1)+2*nu)*lr1**(-1)*(x1p-y1)**2*((x1p-y1)**2+(x2p+(-1) &
            *y2)**2)**(-1)*(x2p-y2)**2*((x1p-y1)**2+(x3-y3) &
            **2)**(-1)*(x3-y3)-((-3._8)+4._8*nu)*lr1**(-1)*(x2p- &
            y2)**2*(lr1+x3-y3)**(-1)-(5._8+4._8*nu*((-3)+2*nu))*(lr2+ &
            x2p-y2)**(-1)*(x3+y3)-(5._8+4._8*nu*((-3)+2*nu))*lr2**( &
            -1)*(x2p-y2)*(lr2+x2p-y2)**(-1)*(x3+y3)-(3 &
            -6*nu+4*nu**2)*lr2**(-1)*(x2p-y2)**2*(lr2+x3+y3)**(-1)+4*(( &
            -1)+nu)*((-1._8)+2._8*nu)*((x1p-y1)**2+(x2p-y2)**2)**(-2)* &
            (x2p-y2)**2*y3*(2*x3+y3)-2._8*((-1._8)+nu)*((-1._8)+2._8*nu)* &
            lr2**(-2)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2) &
            **2*y3*(2*x3+y3)+2*((-1._8)+2._8*nu)*((x1p-y1)**2+(x2p- &
            y2)**2)**(-1)*((x1p-y1)**2+(x3+y3)**2)**(-1)*(((-1._8)+2._8*nu) &
            *lr2*(x1p-y1)**2*(x3+y3)-((-1._8)+nu)*y3*(2*x3+y3)*( &
            (x1p-y1)**2+(x3+y3)**2))-4._8*lr2**(-1)*((x1p-y1)**2+( &
            x2p-y2)**2)**(-2)*(x2p-y2)**2*((x1p-y1)**2+(x3+ &
            y3)**2)**(-1)*(x1p**4*y3-4._8*x1p**3*y1*y3-3*nu*(x3+y3)*( &
            (x1p-y1)**2+(x3+y3)**2)*((x1p-y1)**2+(x2p-y2)**2+( &
            x3+y3)**2)+2*nu**2*(x3+y3)*((x1p-y1)**2+(x3+y3)**2)*((x1p+ &
            (-1)*y1)**2+(x2p-y2)**2+(x3+y3)**2)-2*x1p*y1*y3*((x2p+( &
            -1)*y2)**2+2*(y1**2+(x3+y3)*(2*x3+y3)))+y3*(2*x3**4+7* &
            x3**3*y3+(y1**2+y3**2)*(y1**2+y2**2+y3**2)+x3*y3*(6*y1**2+3* &
            y2**2+5*y3**2)+x3**2*(4._8*y1**2+2*y2**2+9._8*y3**2)+x2p**2*(y1**2+ &
            (x3+y3)*(2*x3+y3))-2*x2p*y2*(y1**2+(x3+y3)*(2*x3+y3)))+ &
            x1p**2*y3*((x2p-y2)**2+2._8*(3*y1**2+(x3+y3)*(2*x3+y3))))+( &
            -2)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)**2* &
            ((x1p-y1)**2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p- &
            y2)**2+(x3+y3)**2)**(-3._8/2._8)*(x1p**4*(y3-3._8*nu*(x3+y3)+2* &
            nu**2*(x3+y3))-4._8*x1p**3*y1*(y3-3*nu*(x3+y3)+2._8*nu**2*( &
            x3+y3))-3*nu*(y1**2+(x3+y3)**2)*(x3**3+3._8*x3**2*y3+y3*( &
            y1**2+y2**2-lr2*y3+y3**2)+x3*(y1**2+y2**2-2._8*lr2*y3+3* &
            y3**2))+2*nu**2*(y1**2+(x3+y3)**2)*(x3**3+3._8*x3**2*y3+y3*( &
            y1**2+y2**2-lr2*y3+y3**2)+x3*(y1**2+y2**2-2*lr2*y3+3* &
            y3**2))+y3*(2*x3**4+7*x3**3*y3+(y1**2+y3**2)*(y1**2+y2**2+ &
            y3**2)+x3*y3*(6._8*y1**2+3*y2**2+5*y3**2)+x3**2*(4*y1**2+2* &
            y2**2+9._8*y3**2)-lr2*(2*x3+y3)*(y1**2+(x3+y3)**2))+2._8*x2p* &
            y2*(3._8*nu*(x3+y3)*(y1**2+(x3+y3)**2)-2._8*nu**2*(x3+y3)*( &
            y1**2+(x3+y3)**2)-y3*(y1**2+(x3+y3)*(2._8*x3+y3)))+x2p**2*(( &
            -3)*nu*(x3+y3)*(y1**2+(x3+y3)**2)+2*nu**2*(x3+y3)*(y1**2+( &
            x3+y3)**2)+y3*(y1**2+(x3+y3)*(2._8*x3+y3)))+2._8*x1p*y1*(((-1._8)+nu) &
            *((-1._8)+2._8*nu)*lr2*y3*(2*x3+y3)+x2p**2*(-y3+3*nu*(x3+y3) &
            -2*nu**2*(x3+y3))+2*x2p*y2*(y3-3._8*nu*(x3+y3)+2._8*nu**2* &
            (x3+y3))+3*nu*(x3+y3)*(2*y1**2+y2**2+2*(x3+y3)**2)-2* &
            nu**2*(x3+y3)*(2*y1**2+y2**2+2*(x3+y3)**2)-y3*(2* &
            y1**2+y2**2+2*(x3+y3)*(2*x3+y3)))+x1p**2*(-((-1._8)+nu)*(( &
            -1)+2*nu)*lr2*y3*(2*x3+y3)+2*x2p*y2*(-y3+3*nu*(x3+y3) &
            -2*nu**2*(x3+y3))+x2p**2*(y3-3*nu*(x3+y3)+2._8*nu**2*(x3+ &
            y3))-3*nu*(x3+y3)*(6*y1**2+y2**2+2*(x3+y3)**2)+2*nu**2*( &
            x3+y3)*(6._8*y1**2+y2**2+2*(x3+y3)**2)+y3*(6._8*y1**2+y2**2+2*(x3+ &
            y3)*(2*x3+y3))))+2*lr2**(-1)*((x1p-y1)**2+(x2p-y2) &
            **2)**(-1)*((x1p-y1)**2+(x3+y3)**2)**(-1)*(-x3* &
            y1**2*y2**2+2*x3**4*y3+4*x3**2*y1**2*y3+y1**4._8*y3+6._8*x3**2* &
            y2**2*y3+2*y1**2*y2**2*y3+7*x3**3*y3**2+6._8*x3*y1**2*y3**2+ &
            9._8*x3*y2**2*y3**2+9._8*x3**2*y3**3+2*y1**2*y3**3+3*y2**2* &
            y3**3+5*x3*y3**4+y3**5-nu*(x3+y3)*(3*x3**4+6._8*x3**2* &
            y1**2+3*y1**4+9._8*x3**2*y2**2+5._8*y1**2*y2**2+6._8*x3*(2._8*(x3**2+ &
            y1**2)+3*y2**2)*y3+3._8*(6._8*x3**2+2*y1**2+3._8*y2**2)*y3**2+12* &
            x3*y3**3+3*y3**4)+x1p**4*(y3-3._8*nu*(x3+y3)+2._8*nu**2*(x3+y3) &
            )-4._8*x1p**3*y1*(y3-3._8*nu*(x3+y3)+2._8*nu**2*(x3+y3))+2._8* &
            nu**2*(x3+y3)*(x3**4+y1**4+4._8*x3**3*y3+3*y2**2*y3**2+y3**4+ &
            y1**2*(y2**2+2*y3**2)+x3**2*(2*y1**2+3._8*y2**2+6._8*y3**2)+x3*( &
            4._8*y1**2*y3+6._8*y2**2*y3+4._8*y3**3))+x2p**2*(-x3*y1**2+6._8* &
            x3**2*y3+2*y1**2*y3+9._8*x3*y3**2+3._8*y3**3+2._8*nu**2*(x3+y3)*( &
            y1**2+3._8*(x3+y3)**2)-nu*(x3+y3)*(5._8*y1**2+9._8*(x3+y3)**2))+ &
            2*x2p*y2*(x3*y1**2-6._8*x3**2*y3-2._8*y1**2*y3-9._8*x3* &
            y3**2-3*y3**3-2*nu**2*(x3+y3)*(y1**2+3._8*(x3+y3)**2)+nu* &
            (x3+y3)*(5._8*y1**2+9._8*(x3+y3)**2))+2._8*x1p*y1*(x2p**2*x3-2* &
            x2p*x3*y2+x3*y2**2-2*x2p**2*y3-4._8*x3**2*y3-2*y1**2* &
            y3+4*x2p*y2*y3-2*y2**2*y3-6._8*x3*y3**2-2*y3**3 &
            -2*nu**2*(x3+y3)*((x2p-y2)**2+2*(y1**2+(x3+y3)**2))+nu*( &
            x3+y3)*(5._8*(x2p-y2)**2+6._8*(y1**2+(x3+y3)**2)))+x1p**2*((-1) &
            *x2p**2*x3+2*x2p*x3*y2-x3*y2**2+2*x2p**2*y3+4._8*x3**2* &
            y3+6*y1**2*y3-4._8*x2p*y2*y3+2*y2**2*y3+6._8*x3*y3**2+2* &
            y3**3+2._8*nu**2*(x3+y3)*((x2p-y2)**2+2*(3._8*y1**2+(x3+y3) &
            **2))-nu*(x3+y3)*(5._8*(x2p-y2)**2+6._8*(3*y1**2+(x3+y3) &
            **2))))+xLogy(3._8-4*nu,lr1+x3-y3)+xLogy((-3)+6._8*nu-4._8* &
            nu**2,lr2+x3+y3))

    END FUNCTION J1123d2

    !---------------------------------------------------------------
    !> function J1123d3
    !! computes the derivative of the J integral J1123,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1123d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1123d3=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(2._8*((-1._8)+2._8*nu)* &
            lr1*(x1p-y1)**2*((x1p-y1)**2+(x2p-y2)**2)**(-1)*( &
            x2p-y2)*((x1p-y1)**2+(x3-y3)**2)**(-1)-(( &
            -3._8)+4._8*nu)*lr1**(-1)*(lr1+x2p-y2)**(-1)*(x3-y3)**2+( &
            -2._8)*((-1._8)+2._8*nu)*lr1**(-1)*(x1p-y1)**2*((x1p-y1)**2+( &
            x2p-y2)**2)**(-1)*(x2p-y2)*((x1p-y1)**2+(x3+(-1) &
            *y3)**2)**(-1)*(x3-y3)**2-((-3)+4._8*nu)*(x2p- &
            y2)*(lr1+x3-y3)**(-1)-((-3)+4._8*nu)*lr1**(-1)*(x2p+(-1) &
            *y2)*(x3-y3)*(lr1+x3-y3)**(-1)-2*lr2**(-1)*x3*( &
            lr2-x2p+y2)**(-1)*(x3+y3)-(5._8+4._8*nu*((-3)+2._8*nu))* &
            lr2**(-1)*(lr2+x2p-y2)**(-1)*(x3+y3)**2-(3._8-6._8*nu+4._8* &
            nu**2)*(x2p-y2)*(lr2+x3+y3)**(-1)-(3._8-6._8*nu+4._8* &
            nu**2)*lr2**(-1)*(x2p-y2)*(x3+y3)*(lr2+x3+y3)**(-1)-2._8*( &
            (-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-2._8)*((x1p-y1)**2+(x2p-y2) &
            **2)**(-1)*(x2p-y2)*y3*(x3+y3)*(2*x3+y3)+4._8*((-1._8)+nu)*( &
            (-1._8)+2._8*nu)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p- &
            y2)*y3*(x3+y3)*(2*x3+y3)*((x1p-y1)**2+(x3+y3)**2)**(-1)+ &
            2._8*((-1._8)+2._8*nu)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+( &
            -1)*y2)*((x1p-y1)**2+(x3+y3)**2)**(-1)*(((-1._8)+2._8*nu)*lr2* &
            (x1p-y1)**2-2._8*((-1._8)+nu)*y3*((x1p-y1)**2+(x3+y3)*( &
            3._8*x3+2*y3)))+2*lr2**(-1)*((x1p-y1)**2+(x2p-y2)**2)**( &
            -1._8)*(x2p-y2)*((x1p-y1)**2+(x3+y3)**2)**(-1)*(nu*(( &
            -3._8)+2*nu)*x1p**4+4._8*(3._8-2._8*nu)*nu*x1p**3*y1-x3**2* &
            y1**2+4*x2p**2*x3*y3+8._8*x3**3*y3+6*x3*y1**2*y3-8*x2p* &
            x3*y2*y3+4._8*x3*y2**2*y3+3*x2p**2*y3**2+21*x3**2*y3**2+5._8* &
            y1**2*y3**2-6._8*x2p*y2*y3**2+3*y2**2*y3**2+18*x3*y3**3+5._8* &
            y3**4+2._8*x1p*y1*((1+2*(7-4._8*nu)*nu)*x3**2+(3._8-2._8*nu)* &
            nu*(2*y1**2+(x2p-y2)**2)-2*(3._8+2._8*nu*((-7._8)+4._8*nu))* &
            x3*y3+((-5._8)+2*(7-4*nu)*nu)*y3**2)+x1p**2*(((-1._8)+2._8*nu*(( &
            -7._8)+4._8*nu))*x3**2+nu*((-3._8)+2._8*nu)*(6._8*y1**2+(x2p-y2)**2)+ &
            2._8*(3+2*nu*((-7)+4*nu))*x3*y3+(5+2*nu*((-7._8)+4*nu))*y3**2) &
            -nu*(y1**2+3*(x3+y3)**2)*(3._8*y1**2+3*(x2p-y2)**2+ &
            5._8*(x3+y3)**2)+2*nu**2*(5*x3**4+y1**4+y1**2*y2**2+20*x3**3* &
            y3+4._8*y1**2*y3**2+3*y2**2*y3**2+5._8*y3**4+x3**2*(4._8*y1**2+3* &
            y2**2+30._8*y3**2)+x3*(8._8*y1**2*y3+6._8*y2**2*y3+20._8*y3**3)+x2p**2* &
            (y1**2+3._8*(x3+y3)**2)-2._8*x2p*y2*(y1**2+3*(x3+y3)**2)))-4* &
            lr2**(-1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2) &
            *(x3+y3)*((x1p-y1)**2+(x3+y3)**2)**(-2)*(x1p**4*y3-4._8* &
            x1p**3*y1*y3-3._8*nu*(x3+y3)*((x1p-y1)**2+(x3+y3)**2)*(( &
            x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)+2._8*nu**2*(x3+y3)*(( &
            x1p-y1)**2+(x3+y3)**2)*((x1p-y1)**2+(x2p-y2)**2+( &
            x3+y3)**2)-2*x1p*y1*y3*((x2p-y2)**2+2._8*(y1**2+(x3+y3)* &
            (2*x3+y3)))+y3*(2._8*x3**4+7._8*x3**3*y3+(y1**2+y3**2)*(y1**2+ &
            y2**2+y3**2)+x3*y3*(6._8*y1**2+3._8*y2**2+5._8*y3**2)+x3**2*(4._8* &
            y1**2+2*y2**2+9._8*y3**2)+x2p**2*(y1**2+(x3+y3)*(2*x3+y3))-2* &
            x2p*y2*(y1**2+(x3+y3)*(2*x3+y3)))+x1p**2*y3*((x2p-y2)**2+ &
            2*(3*y1**2+(x3+y3)*(2*x3+y3))))-2._8*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*(x2p-y2)*(x3+y3)*((x1p-y1)**2+(x3+ &
            y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**( &
            -3._8/2._8)*(x1p**4*(y3-3._8*nu*(x3+y3)+2._8*nu**2*(x3+y3))-4* &
            x1p**3*y1*(y3-3*nu*(x3+y3)+2*nu**2*(x3+y3))-3._8*nu*( &
            y1**2+(x3+y3)**2)*(x3**3+3*x3**2*y3+y3*(y1**2+y2**2-lr2* &
            y3+y3**2)+x3*(y1**2+y2**2-2*lr2*y3+3*y3**2))+2._8*nu**2*( &
            y1**2+(x3+y3)**2)*(x3**3+3*x3**2*y3+y3*(y1**2+y2**2-lr2* &
            y3+y3**2)+x3*(y1**2+y2**2-2*lr2*y3+3*y3**2))+y3*(2*x3**4+ &
            7._8*x3**3*y3+(y1**2+y3**2)*(y1**2+y2**2+y3**2)+x3*y3*(6*y1**2+ &
            3._8*y2**2+5*y3**2)+x3**2*(4._8*y1**2+2*y2**2+9._8*y3**2)-lr2*( &
            2._8*x3+y3)*(y1**2+(x3+y3)**2))+2*x2p*y2*(3._8*nu*(x3+y3)*(y1**2+ &
            (x3+y3)**2)-2*nu**2*(x3+y3)*(y1**2+(x3+y3)**2)-y3*( &
            y1**2+(x3+y3)*(2*x3+y3)))+x2p**2*((-3)*nu*(x3+y3)*(y1**2+(x3+ &
            y3)**2)+2*nu**2*(x3+y3)*(y1**2+(x3+y3)**2)+y3*(y1**2+(x3+y3)* &
            (2*x3+y3)))+2._8*x1p*y1*(((-1._8)+nu)*((-1._8)+2._8*nu)*lr2*y3*(2*x3+ &
            y3)+x2p**2._8*(-y3+3*nu*(x3+y3)-2._8*nu**2*(x3+y3))+2*x2p* &
            y2*(y3-3._8*nu*(x3+y3)+2*nu**2*(x3+y3))+3*nu*(x3+y3)*(2* &
            y1**2+y2**2+2._8*(x3+y3)**2)-2._8*nu**2*(x3+y3)*(2*y1**2+y2**2+ &
            2*(x3+y3)**2)-y3*(2*y1**2+y2**2+2._8*(x3+y3)*(2*x3+y3)))+ &
            x1p**2*(-((-1._8)+nu)*((-1._8)+2._8*nu)*lr2*y3*(2*x3+y3)+2*x2p* &
            y2*(-y3+3._8*nu*(x3+y3)-2._8*nu**2*(x3+y3))+x2p**2*(y3 &
            -3*nu*(x3+y3)+2*nu**2*(x3+y3))-3._8*nu*(x3+y3)*(6*y1**2+ &
            y2**2+2._8*(x3+y3)**2)+2._8*nu**2*(x3+y3)*(6*y1**2+y2**2+2*(x3+y3) &
            **2)+y3*(6._8*y1**2+y2**2+2._8*(x3+y3)*(2*x3+y3))))+xLogy((-2._8),lr2+( &
            -1)*x2p+y2)+xLogy(3._8-4._8*nu,lr1+x2p-y2)+xLogy((-5._8)+4._8*(3._8 &
            -2._8*nu)*nu,lr2+x2p-y2))

    END FUNCTION J1123d3

    !---------------------------------------------------------------
    !> function J2112d1
    !! computes the derivative of the J integral J2112,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2112d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2112d1=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*((1._8+8._8*((-1._8)+nu)*nu)*lr2**( &
            -1)*(x1p-y1)+lr1**(-1)*(-x1p+y1)-4._8*((-1._8)+nu)*((-1._8)+ &
            2._8*nu)*lr2**(-1)*(x1p-y1)*(x3+y3)*(lr2+x3+y3)**(-1)+2._8*x3* &
            (x1p-y1)*y3*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8))*G**(-1)

    END FUNCTION J2112d1

    !---------------------------------------------------------------
    !> function J2112d2
    !! computes the derivative of the J integral J2112,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2112d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2112d2=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*((1._8+8*((-1._8)+nu)*nu)*lr2**( &
            -1)*(x2p-y2)+lr1**(-1)*(-x2p+y2)-4._8*((-1._8)+nu)*((-1._8)+ &
            2._8*nu)*lr2**(-1)*(x2p-y2)*(x3+y3)*(lr2+x3+y3)**(-1)+2._8*x3* &
            (x2p-y2)*y3*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8))*G**(-1)

    END FUNCTION J2112d2

    !---------------------------------------------------------------
    !> function J2112d3
    !! computes the derivative of the J integral J2112,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2112d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2112d3=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(lr1**(-1)*(- &
            x3+y3)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*(x3+y3)*(lr2+x3+y3)**(-1)+( &
            -4._8)*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-1)*(x3+y3)**2*(lr2+x3+y3)**( &
            -1._8)+lr2**(-1)*(x3-y3-8._8*nu*(x3+y3)+8._8*nu**2*(x3+y3))+ &
            2._8*x3*y3*(x3+y3)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2) &
            **(-3._8/2._8)+xLogy((-4._8)*((-1._8)+nu)*((-1._8)+2*nu),lr2+x3+y3))

    END FUNCTION J2112d3

    !---------------------------------------------------------------
    !> function J2113d1
    !! computes the derivative of the J integral J2113,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2113d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2113d1=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*(-lr1**(-1)*(x1p- &
            y1)*(x2p-y2)*(lr1+x3-y3)**(-1)-((-1._8)-2._8*nu+ &
            4._8*nu**2)*lr2**(-1)*(x1p-y1)*(x2p-y2)*(lr2+x3+y3)**( &
            -1._8)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*(x1p-y1)*((x1p-y1) &
            **2+(x2p-y2)**2)**(-2)*(x2p-y2)*y3*(2*x3+y3)+2*(( &
            -1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-2._8)*(x1p-y1)*((x1p-y1)**2+ &
            (x2p-y2)**2)**(-1)*(x2p-y2)*y3*(2*x3+y3)-4*lr2**( &
            -1)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*( &
            x2p-y2)*(y3-3*nu*(x3+y3)+2*nu**2*(x3+y3))+2._8*(x1p+(-1._8) &
            *y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*( &
            (x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)*((-3)*nu* &
            (x3*(x3**2+(x1p-y1)**2+(x2p-y2)**2)+(x3*((-2)*lr2+3* &
            x3)+(x1p-y1)**2+(x2p-y2)**2)*y3-(lr2-3*x3)* &
            y3**2+y3**3)+2*nu**2*(x3*(x3**2+(x1p-y1)**2+(x2p-y2) &
            **2)+(x3*((-2)*lr2+3*x3)+(x1p-y1)**2+(x2p-y2)**2)*y3+ &
            (-1._8)*(lr2-3*x3)*y3**2+y3**3)+y3*((x1p-y1)**2+(x2p- &
            y2)**2-(lr2-x3-y3)*(2*x3+y3)))+4._8*lr2**(-1)*( &
            x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-2._8)*(x2p+(-1) &
            *y2)*((-3._8)*nu*(x3+y3)*((x1p-y1)**2+(x2p-y2)**2+(x3+ &
            y3)**2)+2._8*nu**2*(x3+y3)*((x1p-y1)**2+(x2p-y2)**2+(x3+ &
            y3)**2)+y3*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)*(2._8*x3+y3) &
            )))*G**(-1)

    END FUNCTION J2113d1

    !---------------------------------------------------------------
    !> function J2113d2
    !! computes the derivative of the J integral J2113,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2113d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2113d2=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(-lr1**(-1)*( &
            x2p-y2)**2*(lr1+x3-y3)**(-1)-((-1._8)-2._8*nu+4._8* &
            nu**2)*lr2**(-1)*(x2p-y2)**2*(lr2+x3+y3)**(-1)+2._8*((-1._8)+nu) &
            *((-1._8)+2._8*nu)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*y3*( &
            2*x3+y3)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-2._8)*(x2p-y2)**2*y3*(2._8*x3+y3)+2._8*((-1._8)+nu) &
            *((-1._8)+2._8*nu)*lr2**(-2)*((x1p-y1)**2+(x2p-y2)**2)**( &
            -1)*(x2p-y2)**2*y3*(2*x3+y3)+2*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*(x2p-y2)**2*((x1p-y1)**2+(x2p- &
            y2)**2+(x3+y3)**2)**(-3._8/2._8)*((-3._8)*nu*(x3*(x3**2+(x1p-y1) &
            **2+(x2p-y2)**2)+(x3*((-2)*lr2+3*x3)+(x1p-y1)**2+(x2p+( &
            -1)*y2)**2)*y3-(lr2-3*x3)*y3**2+y3**3)+2*nu**2*(x3* &
            (x3**2+(x1p-y1)**2+(x2p-y2)**2)+(x3*((-2._8)*lr2+3*x3)+( &
            x1p-y1)**2+(x2p-y2)**2)*y3-(lr2-3*x3)*y3**2+ &
            y3**3)+y3*((x1p-y1)**2+(x2p-y2)**2-(lr2-x3+( &
            -1)*y3)*(2._8*x3+y3)))+4._8*lr2**(-1)*((x1p-y1)**2+(x2p- &
            y2)**2)**(-2)*(x2p-y2)**2*((-3._8)*nu*(x3+y3)*((x1p- &
            y1)**2+(x2p-y2)**2+(x3+y3)**2)+2._8*nu**2*(x3+y3)*((x1p- &
            y1)**2+(x2p-y2)**2+(x3+y3)**2)+y3*((x1p-y1)**2+(x2p+(-1) &
            *y2)**2+(x3+y3)*(2*x3+y3)))+lr2**(-1)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*(6._8*nu*(x3+y3)*((x1p-y1)**2+3*(x2p+(-1._8) &
            *y2)**2+(x3+y3)**2)-4._8*nu**2*(x3+y3)*((x1p-y1)**2+3._8*( &
            x2p-y2)**2+(x3+y3)**2)-2._8*y3*((x1p-y1)**2+3._8*(x2p+( &
            -1)*y2)**2+(x3+y3)*(2._8*x3+y3)))+log(lr2+x3+y3)+xLogy((-1._8),lr1+x3+( &
            -1)*y3)+xLogy(2._8*(1._8-2._8*nu)*nu,lr2+x3+y3))

    END FUNCTION J2113d2

    !---------------------------------------------------------------
    !> function J2113d3
    !! computes the derivative of the J integral J2113,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2113d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2113d3=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*((-x2p+y2)*(lr1+x3+(-1) &
            *y3)**(-1)-lr1**(-1)*(x2p-y2)*(x3-y3)*(lr1+x3+( &
            -1)*y3)**(-1)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*(x2p-y2)*y3-((-1._8)-2._8*nu+4._8* &
            nu**2)*(x2p-y2)*(lr2+x3+y3)**(-1)-((-1._8)-2._8*nu+4._8* &
            nu**2)*lr2**(-1)*(x2p-y2)*(x3+y3)*(lr2+x3+y3)**(-1)+2._8*(( &
            -1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-2)*((x1p-y1)**2+(x2p-y2) &
            **2)**(-1)*(x2p-y2)*y3*(x3+y3)*(2._8*x3+y3)-2._8*lr2**(-1) &
            *((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(3* &
            nu*((-3._8)+2._8*nu)*x3**2+nu*((-3._8)+2._8*nu)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)+2*(2._8-9._8*nu+6._8*nu**2)*x3*y3+(3._8-9._8*nu+6._8* &
            nu**2)*y3**2)+2*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+( &
            -1)*y2)*(x3+y3)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2) &
            **(-3._8/2._8)*((-3._8)*nu*(x3*(x3**2+(x1p-y1)**2+(x2p-y2) &
            **2)+(x3*((-2)*lr2+3*x3)+(x1p-y1)**2+(x2p-y2)**2)*y3+ &
            (-1)*(lr2-3*x3)*y3**2+y3**3)+2._8*nu**2*(x3*(x3**2+(x1p- &
            y1)**2+(x2p-y2)**2)+(x3*((-2._8)*lr2+3*x3)+(x1p-y1)**2+( &
            x2p-y2)**2)*y3-(lr2-3._8*x3)*y3**2+y3**3)+y3*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2-(lr2-x3-y3)*(2._8*x3+ &
            y3))))*G**(-1)

    END FUNCTION J2113d3

    !---------------------------------------------------------------
    !> function J2123d1
    !! computes the derivative of the J integral J2123,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2123d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2123d1=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(-lr1**(-1)*( &
            x1p-y1)**2*(lr1+x3-y3)**(-1)-((-1._8)-2._8*nu+4._8* &
            nu**2)*lr2**(-1)*(x1p-y1)**2*(lr2+x3+y3)**(-1)-4._8*((-1._8)+ &
            nu)*((-1._8)+2._8*nu)*(x1p-y1)**2*((x1p-y1)**2+(x2p- &
            y2)**2)**(-2)*y3*(2._8*x3+y3)+2._8*((-1._8)+nu)*((-1._8)+2._8*nu)*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2)**(-1)*y3*(2*x3+y3)+2._8*((-1._8)+nu)* &
            ((-1._8)+2._8*nu)*lr2**(-2)*(x1p-y1)**2*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*y3*(2*x3+y3)+2._8*(x1p-y1)**2*((x1p+(-1) &
            *y1)**2+(x2p-y2)**2)**(-1)*((x1p-y1)**2+(x2p-y2) &
            **2+(x3+y3)**2)**(-3._8/2._8)*((-3._8)*nu*(x3*(x3**2+(x1p-y1)**2+( &
            x2p-y2)**2)+(x3*((-2._8)*lr2+3*x3)+(x1p-y1)**2+(x2p- &
            y2)**2)*y3-(lr2-3._8*x3)*y3**2+y3**3)+2._8*nu**2*(x3*( &
            x3**2+(x1p-y1)**2+(x2p-y2)**2)+(x3*((-2._8)*lr2+3*x3)+( &
            x1p-y1)**2+(x2p-y2)**2)*y3-(lr2-3*x3)*y3**2+ &
            y3**3)+y3*((x1p-y1)**2+(x2p-y2)**2-(lr2-x3+( &
            -1)*y3)*(2._8*x3+y3)))+4._8*lr2**(-1)*(x1p-y1)**2*((x1p- &
            y1)**2+(x2p-y2)**2)**(-2._8)*((-3._8)*nu*(x3+y3)*((x1p-y1) &
            **2+(x2p-y2)**2+(x3+y3)**2)+2._8*nu**2*(x3+y3)*((x1p-y1) &
            **2+(x2p-y2)**2+(x3+y3)**2)+y3*((x1p-y1)**2+(x2p- &
            y2)**2+(x3+y3)*(2._8*x3+y3)))+lr2**(-1)*((x1p-y1)**2+(x2p+(-1) &
            *y2)**2)**(-1)*(6._8*nu*(x3+y3)*(3._8*(x1p-y1)**2+(x2p- &
            y2)**2+(x3+y3)**2)-4._8*nu**2*(x3+y3)*(3._8*(x1p-y1)**2+(x2p+ &
            (-1)*y2)**2+(x3+y3)**2)-2._8*y3*(3._8*(x1p-y1)**2+(x2p- &
            y2)**2+(x3+y3)*(2._8*x3+y3)))+log(lr2+x3+y3)+xLogy((-1._8),lr1+x3- &
            y3)+xLogy(2._8*(1._8-2._8*nu)*nu,lr2+x3+y3))

    END FUNCTION J2123d1

    !---------------------------------------------------------------
    !> function J2123d2
    !! computes the derivative of the J integral J2123,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2123d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2123d2=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*(-lr1**(-1)*(x1p- &
            y1)*(x2p-y2)*(lr1+x3-y3)**(-1)-((-1._8)-2._8*nu+ &
            4._8*nu**2)*lr2**(-1)*(x1p-y1)*(x2p-y2)*(lr2+x3+y3)**( &
            -1)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*(x1p-y1)*((x1p-y1) &
            **2+(x2p-y2)**2)**(-2)*(x2p-y2)*y3*(2._8*x3+y3)+2*(( &
            -1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-2)*(x1p-y1)*((x1p-y1)**2+ &
            (x2p-y2)**2)**(-1)*(x2p-y2)*y3*(2._8*x3+y3)-4._8*lr2**( &
            -1._8)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*( &
            x2p-y2)*(y3-3._8*nu*(x3+y3)+2._8*nu**2*(x3+y3))+2*(x1p+(-1) &
            *y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*( &
            (x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)*((-3._8)*nu* &
            (x3*(x3**2+(x1p-y1)**2+(x2p-y2)**2)+(x3*((-2._8)*lr2+3._8* &
            x3)+(x1p-y1)**2+(x2p-y2)**2)*y3-(lr2-3._8*x3)* &
            y3**2+y3**3)+2._8*nu**2*(x3*(x3**2+(x1p-y1)**2+(x2p-y2) &
            **2)+(x3*((-2)*lr2+3._8*x3)+(x1p-y1)**2+(x2p-y2)**2)*y3+ &
            (-1)*(lr2-3*x3)*y3**2+y3**3)+y3*((x1p-y1)**2+(x2p- &
            y2)**2-(lr2-x3-y3)*(2._8*x3+y3)))+4._8*lr2**(-1)*( &
            x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-2)*(x2p+(-1) &
            *y2)*((-3._8)*nu*(x3+y3)*((x1p-y1)**2+(x2p-y2)**2+(x3+ &
            y3)**2)+2._8*nu**2*(x3+y3)*((x1p-y1)**2+(x2p-y2)**2+(x3+ &
            y3)**2)+y3*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)*(2*x3+y3) &
            )))*G**(-1)

    END FUNCTION J2123d2

    !---------------------------------------------------------------
    !> function J2123d3
    !! computes the derivative of the J integral J2123,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2123d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2123d3=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*((-x1p+y1)*(lr1+x3+(-1) &
            *y3)**(-1)-lr1**(-1)*(x1p-y1)*(x3-y3)*(lr1+x3+( &
            -1)*y3)**(-1)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*(x1p-y1)*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2)**(-1)*y3-((-1._8)-2._8*nu+4._8* &
            nu**2)*(x1p-y1)*(lr2+x3+y3)**(-1)-((-1._8)-2._8*nu+4._8* &
            nu**2)*lr2**(-1)*(x1p-y1)*(x3+y3)*(lr2+x3+y3)**(-1)+2._8*(( &
            -1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-2)*(x1p-y1)*((x1p-y1)**2+ &
            (x2p-y2)**2)**(-1)*y3*(x3+y3)*(2*x3+y3)-2*lr2**(-1)*( &
            x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(3._8*nu*( &
            (-3._8)+2._8*nu)*x3**2+nu*((-3._8)+2._8*nu)*((x1p-y1)**2+(x2p- &
            y2)**2)+2._8*(2._8-9._8*nu+6._8*nu**2)*x3*y3+(3._8-9._8*nu+6._8*nu**2)* &
            y3**2)+2._8*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1) &
            *(x3+y3)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8) &
            *((-3._8)*nu*(x3*(x3**2+(x1p-y1)**2+(x2p-y2)**2)+(x3*( &
            (-2._8)*lr2+3*x3)+(x1p-y1)**2+(x2p-y2)**2)*y3-(lr2+( &
            -3._8)*x3)*y3**2+y3**3)+2._8*nu**2*(x3*(x3**2+(x1p-y1)**2+(x2p+ &
            (-1)*y2)**2)+(x3*((-2._8)*lr2+3*x3)+(x1p-y1)**2+(x2p-y2) &
            **2)*y3-(lr2-3*x3)*y3**2+y3**3)+y3*((x1p-y1)**2+( &
            x2p-y2)**2-(lr2-x3-y3)*(2._8*x3+y3))))*G**(-1)

    END FUNCTION J2123d3

    !---------------------------------------------------------------
    !> function J3112d1
    !! computes the derivative of the J integral J3112,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3112d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3112d1=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(4._8*((-1._8)+nu)*(( &
            -1)+2._8*nu)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**( &
            -1)*(x2p-y2)+lr1**(-1)*(x1p-y1)*(lr1+x2p-y2)**(-1) &
            *(x3-y3)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-1)*(x1p+(-1) &
            *y1)*(x2p-y2)*(lr2+x3+y3)**(-1)+4._8*lr2**(-1)*x3*(x1p- &
            y1)*(x2p-y2)*y3*(x3+y3)*((x1p-y1)**2+(x3+y3)**2)**( &
            -2)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2*(x1p-y1)*((x1p- &
            y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(x3+y3)*((x1p+(-1) &
            *y1)**2+(x3+y3)**2)**(-1)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-1) &
            *(x1p-y1)**3*((x1p-y1)**2+(x2p-y2)**2)**(-1)*( &
            x2p-y2)*(x3+y3)*((x1p-y1)**2+(x3+y3)**2)**(-1)+2*x3* &
            (x1p-y1)*(x2p-y2)*y3*(x3+y3)*((x1p-y1)**2+(x3+ &
            y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**( &
            -3._8/2._8)-lr2**(-1)*(x1p-y1)*(lr2+x2p-y2)**(-1)*(x3+ &
            7._8*y3+8._8*nu**2*(x3+y3)-8._8*nu*(x3+2._8*y3))+4._8*((-1._8)+nu)*((-1._8)+ &
            2._8*nu)*atan((x1p-y1)*(x2p-y2)**(-1))+4._8*((-1._8)+nu)*(( &
            -1)+2._8*nu)*atan2(lr2*(-x1p+y1),(x2p-y2)*(x3+y3)))

    END FUNCTION J3112d1

    !---------------------------------------------------------------
    !> function J3112d2
    !! computes the derivative of the J integral J3112,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3112d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3112d2=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*((-4._8)*((-1._8)+nu)*( &
            (-1)+2._8*nu)*(x1p-y1)**2*((x1p-y1)**2+(x2p-y2)**2) &
            **(-1)+(lr1+x2p-y2)**(-1)*(x3-y3)+lr1**(-1)*(x2p- &
            y2)*(lr1+x2p-y2)**(-1)*(x3-y3)-4._8*((-1._8)+nu)*((-1._8)+ &
            2._8*nu)*lr2**(-1)*(x2p-y2)**2*(lr2+x3+y3)**(-1)+4._8*((-1._8)+nu) &
            *((-1._8)+2._8*nu)*lr2*(x1p-y1)**2*((x1p-y1)**2+(x2p- &
            y2)**2)**(-1)*(x3+y3)*((x1p-y1)**2+(x3+y3)**2)**(-1)-4._8* &
            ((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-1)*(x1p-y1)**2*((x1p- &
            y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)**2*(x3+y3)*((x1p+( &
            -1)*y1)**2+(x3+y3)**2)**(-1)-2*lr2**(-1)*x3*y3*(x3+y3)*(( &
            x1p-y1)**2+(x3+y3)**2)**(-1)+2*x3*(x2p-y2)**2*y3*( &
            x3+y3)*((x1p-y1)**2+(x3+y3)**2)**(-1)*((x1p-y1)**2+( &
            x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)-lr2**(-1)*(x2p-y2) &
            *(lr2+x2p-y2)**(-1)*(x3+7._8*y3+8._8*nu**2*(x3+y3)-8._8*nu*( &
            x3+2*y3))+(lr2+x2p-y2)**(-1)*(-x3-7._8*y3-8._8* &
            nu**2*(x3+y3)+8._8*nu*(x3+2._8*y3))+xLogy((-4._8)*((-1._8)+nu)*((-1._8)+2._8* &
            nu),lr2+x3+y3))

    END FUNCTION J3112d2

    !---------------------------------------------------------------
    !> function J3112d3
    !! computes the derivative of the J integral J3112,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3112d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3112d3=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(lr1**(-1)*(lr1+x2p+( &
            -1)*y2)**(-1)*(x3-y3)**2-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*( &
            x2p-y2)*(lr2+x3+y3)**(-1)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)* &
            lr2**(-1)*(x2p-y2)*(x3+y3)*(lr2+x3+y3)**(-1)+4._8*lr2**(-1)* &
            x3*(x2p-y2)*y3*(x3+y3)**2*((x1p-y1)**2+(x3+y3)**2) &
            **(-2)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2*(x1p-y1)**2*((x1p+(-1._8) &
            *y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*((x1p-y1) &
            **2+(x3+y3)**2)**(-1)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-1)*( &
            x1p-y1)**2*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+( &
            -1)*y2)*(x3+y3)**2*((x1p-y1)**2+(x3+y3)**2)**(-1)-2._8* &
            lr2**(-1)*(x2p-y2)*y3*(2._8*x3+y3)*((x1p-y1)**2+(x3+y3) &
            **2)**(-1)+2._8*x3*(x2p-y2)*y3*(x3+y3)**2*((x1p-y1) &
            **2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3) &
            **2)**(-3._8/2._8)-lr2**(-1)*(lr2+x2p-y2)**(-1)*(x3+y3)*(x3+ &
            7._8*y3+8._8*nu**2*(x3+y3)-8._8*nu*(x3+2._8*y3))+log(lr1+x2p-y2)+ &
            xLogy((-1._8)-8._8*((-1._8)+nu)*nu,lr2+x2p-y2))

    END FUNCTION J3112d3

    !---------------------------------------------------------------
    !> function J3113d1
    !! computes the derivative of the J integral J3111,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3113d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3113d1=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*(lr1**(-1)*(x1p-y1)+( &
            -1)*(1+8._8*((-1._8)+nu)*nu)*lr2**(-1)*(x1p-y1)+2._8*lr2**(-1)*( &
            x1p-y1)*(lr2+x3+y3)**(-1)*(3._8*x3+2._8*y3-6._8*nu*(x3+y3)+4._8* &
            nu**2*(x3+y3))+2._8*x3*(x1p-y1)*y3*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2+(x3+y3)**2)**(-3._8/2._8)-2._8*((-3._8)+4._8*nu)*x3*(x1p- &
            y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x3+y3)*((x1p+(-1) &
            *y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-1._8/2._8))*G**(-1)

    END FUNCTION J3113d1

    !---------------------------------------------------------------
    !> function J3113d2
    !! computes the derivative of the J integral J3111,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3113d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3113d2=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*(lr1**(-1)*(x2p-y2)+( &
            -1)*(1._8+8._8*((-1._8)+nu)*nu)*lr2**(-1)*(x2p-y2)+2._8*lr2**(-1)*( &
            x2p-y2)*(lr2+x3+y3)**(-1)*(3._8*x3+2._8*y3-6._8*nu*(x3+y3)+4._8* &
            nu**2*(x3+y3))+2._8*x3*(x2p-y2)*y3*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2+(x3+y3)**2)**(-3._8/2._8)-2._8*((-3._8)+4._8*nu)*x3*((x1p+(-1) &
            *y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(x3+y3)*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-1._8/2._8))*G**(-1)

    END FUNCTION J3113d2

    !---------------------------------------------------------------
    !> function J3113d3
    !! computes the derivative of the J integral J3111,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3113d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3113d3=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(lr1**(-1)*(x3+(-1._8) &
            *y3)+lr2**(-1)*(-x3-3._8*y3+8._8*nu*(x3+y3)-8._8*nu**2*( &
            x3+y3))+2._8*(lr2+x3+y3)**(-1)*(3._8*x3+2._8*y3-6._8*nu*(x3+y3)+4._8* &
            nu**2*(x3+y3))+2._8*lr2**(-1)*(x3+y3)*(lr2+x3+y3)**(-1)*(3*x3+2* &
            y3-6._8*nu*(x3+y3)+4._8*nu**2*(x3+y3))+2._8*x3*y3*(x3+y3)*((x1p+( &
            -1._8)*y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)-2._8*((-3._8)+4._8* &
            nu)*x3*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x3+y3)**2*(( &
            x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-1._8/2._8)+2._8*((-3._8)+4._8* &
            nu)*lr2**(-1)*x3*(1-(x3+y3)**2*((x1p-y1)**2+(x2p+(-1._8) &
            *y2)**2+(x3+y3)**2)**(-1))**(-1)+ &
            2._8*((-3._8)+4._8*nu)*ACOTH(lr2**(-1)*(x3+y3))+ &
            xLogy(6._8+4._8*nu*((-3._8)+2._8*nu),lr2+x3+y3))

    END FUNCTION J3113d3

    !---------------------------------------------------------------
    !> function J3123d1
    !! computes the derivative of the J integral J3123,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3123d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3123d1=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(lr1**(-1)*(x1p+(-1) &
            *y1)**2*(lr1+x2p-y2)**(-1)-(1._8+8._8*((-1._8)+nu)*nu)*lr2**( &
            -1)*(x1p-y1)**2*(lr2+x2p-y2)**(-1)-4._8*((-1._8)+nu)*(( &
            -1)+2._8*nu)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p- &
            y2)*(x3+y3)-4._8*lr2**(-1)*x3*(x1p-y1)**2*(x2p-y2)* &
            y3*((x1p-y1)**2+(x3+y3)**2)**(-2)+2*lr2**(-1)*x3*(x2p+(-1._8) &
            *y2)*y3*((x1p-y1)**2+(x3+y3)**2)**(-1)+4._8*((-1._8)+2._8*nu)* &
            lr2*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(x3+ &
            y3)*(nu*x3+((-1._8)+nu)*y3)*((x1p-y1)**2+(x3+y3)**2)**(-1)+ &
            4._8*((-1._8)+2._8*nu)*lr2**(-1)*(x1p-y1)**2*((x1p-y1)**2+( &
            x2p-y2)**2)**(-1)*(x2p-y2)*(x3+y3)*(nu*x3+((-1._8)+nu) &
            *y3)*((x1p-y1)**2+(x3+y3)**2)**(-1)-2._8*x3*(x1p-y1) &
            **2*(x2p-y2)*y3*((x1p-y1)**2+(x3+y3)**2)**(-1)*((x1p+ &
            (-1)*y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)+log(lr1+x2p- &
            y2)+xLogy((-1._8)-8._8*((-1._8)+nu)*nu,lr2+x2p-y2))

    END FUNCTION J3123d1

    !---------------------------------------------------------------
    !> function J3123d2
    !! computes the derivative of the J integral J3123,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3123d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3123d2=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*((x1p-y1)*(lr1+x2p+(-1) &
            *y2)**(-1)+lr1**(-1)*(x1p-y1)*(x2p-y2)*(lr1+x2p- &
            y2)**(-1)-(1+8._8*((-1._8)+nu)*nu)*(x1p-y1)*(lr2+x2p- &
            y2)**(-1)-(1+8._8*((-1._8)+nu)*nu)*lr2**(-1)*(x1p-y1)*( &
            x2p-y2)*(lr2+x2p-y2)**(-1)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*( &
            x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x3+y3)+ &
            2._8*lr2**(-1)*x3*(x1p-y1)*y3*((x1p-y1)**2+(x3+y3)**2) &
            **(-1)-4._8*((-1._8)+2._8*nu)*lr2*(x1p-y1)*((x1p-y1)**2+( &
            x2p-y2)**2)**(-1)*(x3+y3)*(nu*x3+((-1._8)+nu)*y3)*((x1p+(-1) &
            *y1)**2+(x3+y3)**2)**(-1)+4._8*((-1._8)+2._8*nu)*lr2**(-1)*(x1p- &
            y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)**2* &
            (x3+y3)*(nu*x3+((-1._8)+nu)*y3)*((x1p-y1)**2+(x3+y3)**2)**( &
            -1)-2._8*x3*(x1p-y1)*(x2p-y2)**2*y3*((x1p-y1) &
            **2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3) &
            **2)**(-3._8/2._8))*G**(-1)

    END FUNCTION J3123d2

    !---------------------------------------------------------------
    !> function J3123d3
    !! computes the derivative of the J integral J3123,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3123d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3123d3=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(lr1**(-1)*(x1p+(-1) &
            *y1)*(lr1+x2p-y2)**(-1)*(x3-y3)-(1._8+8._8*((-1._8)+nu) &
            *nu)*lr2**(-1)*(x1p-y1)*(lr2+x2p-y2)**(-1)*(x3+y3)+( &
            -4._8)*lr2**(-1)*x3*(x1p-y1)*(x2p-y2)*y3*(x3+y3)*(( &
            x1p-y1)**2+(x3+y3)**2)**(-2)+2._8*lr2**(-1)*(x1p-y1)*(x2p+ &
            (-1)*y2)*y3*((x1p-y1)**2+(x3+y3)**2)**(-1)-4._8*((-1._8)+2._8* &
            nu)*lr2*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1) &
            *(x2p-y2)*(nu*x3+((-1._8)+nu)*y3)*((x1p-y1)**2+(x3+y3) &
            **2)**(-1)+4._8*((-1._8)+2._8*nu)*lr2**(-1)*(x1p-y1)*((x1p- &
            y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(x3+y3)**2*(nu* &
            x3+((-1._8)+nu)*y3)*((x1p-y1)**2+(x3+y3)**2)**(-1)-2._8*x3*( &
            x1p-y1)*(x2p-y2)*y3*(x3+y3)*((x1p-y1)**2+(x3+ &
            y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**( &
            -3._8/2._8)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*atan((x1p-y1)**(-1)*(x2p+( &
            -1)*y2))+4._8*nu*((-1._8)+2._8*nu)*atan2(lr2*(x1p-y1),(x2p- &
            y2)*(x3+y3)))

    END FUNCTION J3123d3

    !---------------------------------------------------------------
    !> function J1212d1
    !! computes the derivative of the J integral J1212,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1212d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1212d1=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*((1._8+8._8*((-1._8)+nu)*nu)*lr2**( &
            -1)*(x1p-y1)+lr1**(-1)*(-x1p+y1)-4._8*((-1._8)+nu)*((-1._8)+ &
            2._8*nu)*lr2**(-1)*(x1p-y1)*(x3+y3)*(lr2+x3+y3)**(-1)+2._8*x3* &
            (x1p-y1)*y3*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2) &
            **(-3._8/2._8))*G**(-1)

    END FUNCTION J1212d1

    !---------------------------------------------------------------
    !> function J1212d2
    !! computes the derivative of the J integral J1212,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1212d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1212d2=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*((1+8._8*((-1._8)+nu)*nu)*lr2**( &
            -1)*(x2p-y2)+lr1**(-1)*(-x2p+y2)-4._8*((-1._8)+nu)*((-1._8)+ &
            2._8*nu)*lr2**(-1)*(x2p-y2)*(x3+y3)*(lr2+x3+y3)**(-1)+2._8*x3* &
            (x2p-y2)*y3*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2) &
            **(-3._8/2._8))*G**(-1)

    END FUNCTION J1212d2

    !---------------------------------------------------------------
    !> function J1212d3
    !! computes the derivative of the J integral J1212,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1212d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1212d3=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(lr1**(-1)*(- &
            x3+y3)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*(x3+y3)*(lr2+x3+y3)**(-1)+( &
            -4._8)*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-1)*(x3+y3)**2*(lr2+x3+y3)**( &
            -1._8)+lr2**(-1)*(x3-y3-8*nu*(x3+y3)+8*nu**2*(x3+y3))+ &
            2._8*x3*y3*(x3+y3)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2) &
            **(-3._8/2._8)+xLogy((-4._8)*((-1._8)+nu)*((-1._8)+2._8*nu),lr2+x3+y3))

    END FUNCTION J1212d3

    !---------------------------------------------------------------
    !> function J1213d1
    !! computes the derivative of the J integral J1213,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1213d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1213d1=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*(-lr1**(-1)*(x1p- &
            y1)*(x2p-y2)*(lr1+x3-y3)**(-1)-((-1._8)-2._8*nu+ &
            4._8*nu**2)*lr2**(-1)*(x1p-y1)*(x2p-y2)*(lr2+x3+y3)**( &
            -1)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*(x1p-y1)*((x1p-y1) &
            **2+(x2p-y2)**2)**(-2)*(x2p-y2)*y3*(2*x3+y3)+2*(( &
            -1)+nu)*((-1._8)+2._8*nu)*lr2**(-2)*(x1p-y1)*((x1p-y1)**2+ &
            (x2p-y2)**2)**(-1)*(x2p-y2)*y3*(2._8*x3+y3)-4._8*lr2**( &
            -1)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*( &
            x2p-y2)*(y3-3._8*nu*(x3+y3)+2._8*nu**2*(x3+y3))+2._8*(x1p+(-1) &
            *y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*( &
            (x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)*((-3._8)*nu* &
            (x3*(x3**2+(x1p-y1)**2+(x2p-y2)**2)+(x3*((-2._8)*lr2+3._8* &
            x3)+(x1p-y1)**2+(x2p-y2)**2)*y3-(lr2-3*x3)* &
            y3**2+y3**3)+2._8*nu**2*(x3*(x3**2+(x1p-y1)**2+(x2p-y2) &
            **2)+(x3*((-2._8)*lr2+3*x3)+(x1p-y1)**2+(x2p-y2)**2)*y3+ &
            (-1)*(lr2-3*x3)*y3**2+y3**3)+y3*((x1p-y1)**2+(x2p- &
            y2)**2-(lr2-x3-y3)*(2._8*x3+y3)))+4._8*lr2**(-1)*( &
            x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-2)*(x2p+(-1) &
            *y2)*((-3._8)*nu*(x3+y3)*((x1p-y1)**2+(x2p-y2)**2+(x3+ &
            y3)**2)+2._8*nu**2*(x3+y3)*((x1p-y1)**2+(x2p-y2)**2+(x3+ &
            y3)**2)+y3*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)*(2*x3+y3) &
            )))*G**(-1)

    END FUNCTION J1213d1

    !---------------------------------------------------------------
    !> function J1213d2
    !! computes the derivative of the J integral J1213,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1213d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1213d2=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(-lr1**(-1)*( &
            x2p-y2)**2*(lr1+x3-y3)**(-1)-((-1._8)-2._8*nu+4._8* &
            nu**2)*lr2**(-1)*(x2p-y2)**2*(lr2+x3+y3)**(-1)+2._8*((-1._8)+nu) &
            *((-1._8)+2._8*nu)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*y3*( &
            2._8*x3+y3)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-2)*(x2p-y2)**2*y3*(2._8*x3+y3)+2._8*((-1._8)+nu) &
            *((-1._8)+2._8*nu)*lr2**(-2)*((x1p-y1)**2+(x2p-y2)**2)**( &
            -1)*(x2p-y2)**2*y3*(2*x3+y3)+2*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*(x2p-y2)**2*((x1p-y1)**2+(x2p- &
            y2)**2+(x3+y3)**2)**(-3._8/2._8)*((-3)*nu*(x3*(x3**2+(x1p-y1) &
            **2+(x2p-y2)**2)+(x3*((-2)*lr2+3*x3)+(x1p-y1)**2+(x2p+( &
            -1)*y2)**2)*y3-(lr2-3*x3)*y3**2+y3**3)+2*nu**2*(x3* &
            (x3**2+(x1p-y1)**2+(x2p-y2)**2)+(x3*((-2)*lr2+3*x3)+( &
            x1p-y1)**2+(x2p-y2)**2)*y3-(lr2-3*x3)*y3**2+ &
            y3**3)+y3*((x1p-y1)**2+(x2p-y2)**2-(lr2-x3+( &
            -1)*y3)*(2*x3+y3)))+4._8*lr2**(-1)*((x1p-y1)**2+(x2p- &
            y2)**2)**(-2)*(x2p-y2)**2*((-3._8)*nu*(x3+y3)*((x1p- &
            y1)**2+(x2p-y2)**2+(x3+y3)**2)+2*nu**2*(x3+y3)*((x1p- &
            y1)**2+(x2p-y2)**2+(x3+y3)**2)+y3*((x1p-y1)**2+(x2p+(-1) &
            *y2)**2+(x3+y3)*(2._8*x3+y3)))+lr2**(-1)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*(6._8*nu*(x3+y3)*((x1p-y1)**2+3._8*(x2p+(-1) &
            *y2)**2+(x3+y3)**2)-4._8*nu**2*(x3+y3)*((x1p-y1)**2+3._8*( &
            x2p-y2)**2+(x3+y3)**2)-2._8*y3*((x1p-y1)**2+3._8*(x2p+( &
            -1)*y2)**2+(x3+y3)*(2._8*x3+y3)))+log(lr2+x3+y3)+xLogy((-1._8),lr1+x3+( &
            -1)*y3)+xLogy(2._8*(1._8-2._8*nu)*nu,lr2+x3+y3))

    END FUNCTION J1213d2

    !---------------------------------------------------------------
    !> function J1213d3
    !! computes the derivative of the J integral J1213,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1213d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1213d3=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*((-x2p+y2)*(lr1+x3+(-1) &
            *y3)**(-1)-lr1**(-1)*(x2p-y2)*(x3-y3)*(lr1+x3+( &
            -1)*y3)**(-1)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*(x2p-y2)*y3-((-1._8)-2._8*nu+4._8* &
            nu**2)*(x2p-y2)*(lr2+x3+y3)**(-1)-((-1._8)-2._8*nu+4._8* &
            nu**2)*lr2**(-1)*(x2p-y2)*(x3+y3)*(lr2+x3+y3)**(-1)+2*(( &
            -1)+nu)*((-1._8)+2._8*nu)*lr2**(-2)*((x1p-y1)**2+(x2p-y2) &
            **2)**(-1)*(x2p-y2)*y3*(x3+y3)*(2._8*x3+y3)-2*lr2**(-1) &
            *((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(3._8* &
            nu*((-3._8)+2._8*nu)*x3**2+nu*((-3._8)+2._8*nu)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)+2._8*(2._8-9._8*nu+6._8*nu**2)*x3*y3+(3._8-9._8*nu+6._8* &
            nu**2)*y3**2)+2._8*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+( &
            -1)*y2)*(x3+y3)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2) &
            **(-3._8/2._8)*((-3._8)*nu*(x3*(x3**2+(x1p-y1)**2+(x2p-y2) &
            **2)+(x3*((-2._8)*lr2+3*x3)+(x1p-y1)**2+(x2p-y2)**2)*y3+ &
            (-1)*(lr2-3*x3)*y3**2+y3**3)+2._8*nu**2*(x3*(x3**2+(x1p- &
            y1)**2+(x2p-y2)**2)+(x3*((-2._8)*lr2+3*x3)+(x1p-y1)**2+( &
            x2p-y2)**2)*y3-(lr2-3*x3)*y3**2+y3**3)+y3*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2-(lr2-x3-y3)*(2*x3+y3))))*G**(-1)

    END FUNCTION J1213d3

    !---------------------------------------------------------------
    !> function J1223d1
    !! computes the derivative of the J integral J1223,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1223d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1223d1=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(-lr1**(-1)*( &
            x1p-y1)**2*(lr1+x3-y3)**(-1)-((-1._8)-2._8*nu+4._8* &
            nu**2)*lr2**(-1)*(x1p-y1)**2*(lr2+x3+y3)**(-1)-4._8*((-1._8)+ &
            nu)*((-1._8)+2._8*nu)*(x1p-y1)**2*((x1p-y1)**2+(x2p- &
            y2)**2)**(-2)*y3*(2*x3+y3)+2._8*((-1._8)+nu)*((-1._8)+2._8*nu)*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2)**(-1)*y3*(2._8*x3+y3)+2._8*((-1._8)+nu)* &
            ((-1._8)+2._8*nu)*lr2**(-2)*(x1p-y1)**2*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*y3*(2._8*x3+y3)+2*(x1p-y1)**2*((x1p+(-1) &
            *y1)**2+(x2p-y2)**2)**(-1)*((x1p-y1)**2+(x2p-y2) &
            **2+(x3+y3)**2)**(-3._8/2._8)*((-3._8)*nu*(x3*(x3**2+(x1p-y1)**2+( &
            x2p-y2)**2)+(x3*((-2)*lr2+3*x3)+(x1p-y1)**2+(x2p- &
            y2)**2)*y3-(lr2-3*x3)*y3**2+y3**3)+2*nu**2*(x3*( &
            x3**2+(x1p-y1)**2+(x2p-y2)**2)+(x3*((-2._8)*lr2+3*x3)+( &
            x1p-y1)**2+(x2p-y2)**2)*y3-(lr2-3*x3)*y3**2+ &
            y3**3)+y3*((x1p-y1)**2+(x2p-y2)**2-(lr2-x3+( &
            -1)*y3)*(2._8*x3+y3)))+4._8*lr2**(-1)*(x1p-y1)**2*((x1p- &
            y1)**2+(x2p-y2)**2)**(-2)*((-3._8)*nu*(x3+y3)*((x1p-y1) &
            **2+(x2p-y2)**2+(x3+y3)**2)+2._8*nu**2*(x3+y3)*((x1p-y1) &
            **2+(x2p-y2)**2+(x3+y3)**2)+y3*((x1p-y1)**2+(x2p- &
            y2)**2+(x3+y3)*(2*x3+y3)))+lr2**(-1)*((x1p-y1)**2+(x2p+(-1) &
            *y2)**2)**(-1)*(6._8*nu*(x3+y3)*(3._8*(x1p-y1)**2+(x2p- &
            y2)**2+(x3+y3)**2)-4._8*nu**2*(x3+y3)*(3._8*(x1p-y1)**2+(x2p+ &
            (-1)*y2)**2+(x3+y3)**2)-2._8*y3*(3._8*(x1p-y1)**2+(x2p- &
            y2)**2+(x3+y3)*(2._8*x3+y3)))+log(lr2+x3+y3)+xLogy((-1._8),lr1+x3- &
            y3)+xLogy(2._8*(1._8-2._8*nu)*nu,lr2+x3+y3))

    END FUNCTION J1223d1

    !---------------------------------------------------------------
    !> function J1223d2
    !! computes the derivative of the J integral J1223,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1223d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1223d2=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*(-lr1**(-1)*(x1p- &
            y1)*(x2p-y2)*(lr1+x3-y3)**(-1)-((-1._8)-2._8*nu+ &
            4._8*nu**2)*lr2**(-1)*(x1p-y1)*(x2p-y2)*(lr2+x3+y3)**( &
            -1)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*(x1p-y1)*((x1p-y1) &
            **2+(x2p-y2)**2)**(-2)*(x2p-y2)*y3*(2*x3+y3)+2*(( &
            -1)+nu)*((-1._8)+2._8*nu)*lr2**(-2)*(x1p-y1)*((x1p-y1)**2+ &
            (x2p-y2)**2)**(-1)*(x2p-y2)*y3*(2*x3+y3)-4._8*lr2**( &
            -1)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*( &
            x2p-y2)*(y3-3._8*nu*(x3+y3)+2._8*nu**2*(x3+y3))+2._8*(x1p+(-1) &
            *y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*( &
            (x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)*((-3._8)*nu* &
            (x3*(x3**2+(x1p-y1)**2+(x2p-y2)**2)+(x3*((-2._8)*lr2+3._8* &
            x3)+(x1p-y1)**2+(x2p-y2)**2)*y3-(lr2-3*x3)* &
            y3**2+y3**3)+2._8*nu**2*(x3*(x3**2+(x1p-y1)**2+(x2p-y2) &
            **2)+(x3*((-2)*lr2+3._8*x3)+(x1p-y1)**2+(x2p-y2)**2)*y3+ &
            (-1)*(lr2-3*x3)*y3**2+y3**3)+y3*((x1p-y1)**2+(x2p- &
            y2)**2-(lr2-x3-y3)*(2._8*x3+y3)))+4._8*lr2**(-1)*( &
            x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-2)*(x2p+(-1) &
            *y2)*((-3._8)*nu*(x3+y3)*((x1p-y1)**2+(x2p-y2)**2+(x3+ &
            y3)**2)+2._8*nu**2*(x3+y3)*((x1p-y1)**2+(x2p-y2)**2+(x3+ &
            y3)**2)+y3*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)*(2*x3+y3) &
            )))*G**(-1)

    END FUNCTION J1223d2

    !---------------------------------------------------------------
    !> function J1223d3
    !! computes the derivative of the J integral J1223,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1223d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1223d3=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*((-x1p+y1)*(lr1+x3+(-1) &
            *y3)**(-1)-lr1**(-1)*(x1p-y1)*(x3-y3)*(lr1+x3+( &
            -1)*y3)**(-1)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*(x1p-y1)*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2)**(-1)*y3-((-1._8)-2._8*nu+4._8* &
            nu**2)*(x1p-y1)*(lr2+x3+y3)**(-1)-((-1._8)-2._8*nu+4._8* &
            nu**2)*lr2**(-1)*(x1p-y1)*(x3+y3)*(lr2+x3+y3)**(-1)+2._8*(( &
            -1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-2)*(x1p-y1)*((x1p-y1)**2+ &
            (x2p-y2)**2)**(-1)*y3*(x3+y3)*(2._8*x3+y3)-2._8*lr2**(-1)*( &
            x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(3._8*nu*( &
            (-3._8)+2._8*nu)*x3**2+nu*((-3._8)+2._8*nu)*((x1p-y1)**2+(x2p- &
            y2)**2)+2._8*(2._8-9._8*nu+6._8*nu**2)*x3*y3+(3._8-9._8*nu+6._8*nu**2)* &
            y3**2)+2._8*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1) &
            *(x3+y3)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8) &
            *((-3._8)*nu*(x3*(x3**2+(x1p-y1)**2+(x2p-y2)**2)+(x3*( &
            (-2._8)*lr2+3*x3)+(x1p-y1)**2+(x2p-y2)**2)*y3-(lr2+( &
            -3._8)*x3)*y3**2+y3**3)+2._8*nu**2*(x3*(x3**2+(x1p-y1)**2+(x2p+ &
            (-1)*y2)**2)+(x3*((-2._8)*lr2+3*x3)+(x1p-y1)**2+(x2p-y2) &
            **2)*y3-(lr2-3*x3)*y3**2+y3**3)+y3*((x1p-y1)**2+( &
            x2p-y2)**2-(lr2-x3-y3)*(2*x3+y3))))*G**(-1)

    END FUNCTION J1223d3

    !---------------------------------------------------------------
    !> function J1212d1
    !! computes the derivative of the J integral J2212,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2212d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2212d1=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(x3**2*(x3**2+(x1p+ &
            (-1)*y1)**2)**(-1)+9._8*x3**2*(9._8*x3**2+(x1p-y1)**2)**(-1)+ &
            4._8*nu**2*x3**2*(nu**2*x3**2+(x1p-y1)**2)**(-1)-((-3) &
            +4._8*nu)*(lr1+x1p-y1)**(-1)*(x2p-y2)-((-3._8)+4._8*nu) &
            *lr1**(-1)*(x1p-y1)*(lr1+x1p-y1)**(-1)*(x2p-y2)+( &
            5._8+4._8*nu*((-3._8)+2._8*nu))*(lr2+x1p-y1)**(-1)*(x2p-y2)+(5+ &
            4._8*nu*((-3._8)+2._8*nu))*lr2**(-1)*(x1p-y1)*(lr2+x1p-y1)**( &
            -1)*(x2p-y2)-4._8*((-1._8)+nu)*lr1**(-1)*(x1p-y1)**2*( &
            lr1+x2p-y2)**(-1)-4._8*((-1._8)+nu)*lr2**(-1)*(x1p-y1) &
            **2*(lr2+x2p-y2)**(-1)+4._8*((-1._8)+nu)*lr1*(x2p-y2)*((x1p+ &
            (-1)*y1)**2+(x3-y3)**2)**(-1)*((x2p-y2)**2+(x3- &
            y3)**2)**(-1)*(x3-y3)**2-4._8*((-1._8)+nu)*lr1**(-1)*(x1p+( &
            -1)*y1)**2*(x2p-y2)*((x1p-y1)**2+(x3-y3)**2)**( &
            -1)*((x2p-y2)**2+(x3-y3)**2)**(-1)*(x3-y3)**2+( &
            -4._8)*((-1._8)+nu)*((-1._8)+2._8*nu)*((x1p-y1)**2+(x2p-y2)**2) &
            **(-1)*(x2p-y2)*(x3+y3)+y3**2*((x1p-y1)**2+y3**2)**( &
            -1)+9._8*y3**2*((x1p-y1)**2+9._8*y3**2)**(-1)+4._8*nu**2*y3**2*( &
            (x1p-y1)**2+nu**2*y3**2)**(-1)+2*lr2**(-1)*x3*(x2p- &
            y2)*y3*((x2p-y2)**2+(x3+y3)**2)**(-1)+4._8*((-1._8)+nu)*((-1._8)+ &
            2._8*nu)*lr2*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p- &
            y2)*(x3+y3)**2*((x2p-y2)**2+(x3+y3)**2)**(-1)-4._8*((-1._8)+ &
            nu)*((-1._8)+2._8*nu)*lr2**(-1)*(x1p-y1)**2*((x1p-y1)**2+( &
            x2p-y2)**2)**(-1)*(x2p-y2)*(x3+y3)**2*((x2p-y2) &
            **2+(x3+y3)**2)**(-1)+4._8*((-1._8)+nu)*lr2*(x2p-y2)*(x3+y3) &
            **2*((x1p-y1)**2+(x3+y3)**2)**(-1)*((x2p-y2)**2+(x3+ &
            y3)**2)**(-1)-4._8*((-1._8)+nu)*lr2**(-1)*(x1p-y1)**2*(x2p+( &
            -1)*y2)*(x3+y3)**2*((x1p-y1)**2+(x3+y3)**2)**(-1)*((x2p+( &
            -1)*y2)**2+(x3+y3)**2)**(-1)-2._8*x3*(x1p-y1)**2*(x2p+(-1) &
            *y2)*y3*((x2p-y2)**2+(x3+y3)**2)**(-1)*((x1p-y1)**2+ &
            (x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)+xLogy(4._8-4._8*nu,lr1+x2p- &
            y2)+xLogy(4._8-4._8*nu,lr2+x2p-y2))

    END FUNCTION J2212d1

    !---------------------------------------------------------------
    !> function J2212d2
    !! computes the derivative of the J integral J2212,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2212d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2212d2=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(-((-3._8)+4._8*nu) &
            *lr1**(-1)*(lr1+x1p-y1)**(-1)*(x2p-y2)**2+(5._8+4._8*nu*(( &
            -3._8)+2._8*nu))*lr2**(-1)*(lr2+x1p-y1)**(-1)*(x2p-y2)**2+( &
            -4._8)*((-1._8)+nu)*(x1p-y1)*(lr1+x2p-y2)**(-1)-4._8*((-1._8)+ &
            nu)*lr1**(-1)*(x1p-y1)*(x2p-y2)*(lr1+x2p-y2)**( &
            -1)-4._8*((-1._8)+nu)*(x1p-y1)*(lr2+x2p-y2)**(-1)-4._8*( &
            (-1._8)+nu)*lr2**(-1)*(x1p-y1)*(x2p-y2)*(lr2+x2p-y2) &
            **(-1)+4._8*((-1._8)+nu)*lr1*(x1p-y1)*((x1p-y1)**2+(x3+(-1) &
            *y3)**2)**(-1)*((x2p-y2)**2+(x3-y3)**2)**(-1)*(x3+( &
            -1)*y3)**2-4._8*((-1._8)+nu)*lr1**(-1)*(x1p-y1)*(x2p- &
            y2)**2*((x1p-y1)**2+(x3-y3)**2)**(-1)*((x2p-y2) &
            **2+(x3-y3)**2)**(-1)*(x3-y3)**2+4._8*((-1._8)+nu)*((-1._8)+ &
            2._8*nu)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)* &
            (x3+y3)-4._8*lr2**(-1)*x3*(x1p-y1)*(x2p-y2)**2*y3*( &
            (x2p-y2)**2+(x3+y3)**2)**(-2)+2._8*lr2**(-1)*x3*(x1p-y1) &
            *y3*((x2p-y2)**2+(x3+y3)**2)**(-1)-4._8*((-1._8)+nu)*((-1._8)+ &
            2._8*nu)*lr2*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**( &
            -1)*(x3+y3)**2*((x2p-y2)**2+(x3+y3)**2)**(-1)-4._8*((-1._8)+ &
            nu)*((-1._8)+2._8*nu)*lr2**(-1)*(x1p-y1)*((x1p-y1)**2+(x2p+ &
            (-1)*y2)**2)**(-1)*(x2p-y2)**2*(x3+y3)**2*((x2p-y2) &
            **2+(x3+y3)**2)**(-1)+4._8*((-1._8)+nu)*lr2*(x1p-y1)*(x3+y3) &
            **2*((x1p-y1)**2+(x3+y3)**2)**(-1)*((x2p-y2)**2+(x3+ &
            y3)**2)**(-1)-4._8*((-1._8)+nu)*lr2**(-1)*(x1p-y1)*(x2p- &
            y2)**2*(x3+y3)**2*((x1p-y1)**2+(x3+y3)**2)**(-1)*((x2p+(-1) &
            *y2)**2+(x3+y3)**2)**(-1)-2._8*x3*(x1p-y1)*(x2p-y2) &
            **2*y3*((x2p-y2)**2+(x3+y3)**2)**(-1)*((x1p-y1)**2+( &
            x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)+xLogy(3._8-4._8*nu,lr1+x1p- &
            y1)+xLogy(5._8+4._8*nu*((-3._8)+2._8*nu),lr2+x1p-y1))

    END FUNCTION J2212d2

    !---------------------------------------------------------------
    !> function J2212d3
    !! computes the derivative of the J integral J2212,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2212d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2212d3=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(x3*(x3**2+(x1p+( &
            -1)*y1)**2)**(-1)*(-x1p+y1)+9._8*x3*(9._8*x3**2+(x1p-y1) &
            **2)**(-1)*(-x1p+y1)+4._8*nu**2*x3*(nu**2*x3**2+(x1p- &
            y1)**2)**(-1)*(-x1p+y1)-((-3._8)+4._8*nu)*lr1**(-1)*(lr1+x1p+ &
            (-1)*y1)**(-1)*(x2p-y2)*(x3-y3)-4._8*((-1._8)+nu)* &
            lr1**(-1)*(x1p-y1)*(lr1+x2p-y2)**(-1)*(x3-y3)+( &
            -4._8)*((-1._8)+nu)*lr1*(x1p-y1)*(x2p-y2)*((x1p-y1) &
            **2+(x3-y3)**2)**(-1)*((x2p-y2)**2+(x3-y3)**2) &
            **(-1)*(x3-y3)-4._8*((-1._8)+nu)*lr1**(-1)*(x1p-y1)*( &
            x2p-y2)*((x1p-y1)**2+(x3-y3)**2)**(-1)*((x2p+(-1) &
            *y2)**2+(x3-y3)**2)**(-1)*(x3-y3)**3+(5._8+4._8*nu*((-3) &
            +2*nu))*lr2**(-1)*(lr2+x1p-y1)**(-1)*(x2p-y2)*(x3+y3) &
            -4._8*((-1._8)+nu)*lr2**(-1)*(x1p-y1)*(lr2+x2p-y2)**(-1) &
            *(x3+y3)-4._8*lr2**(-1)*x3*(x1p-y1)*(x2p-y2)*y3*( &
            x3+y3)*((x2p-y2)**2+(x3+y3)**2)**(-2)+2*lr2**(-1)*(x1p+(-1) &
            *y1)*(x2p-y2)*y3*((x2p-y2)**2+(x3+y3)**2)**(-1)+4._8*( &
            (-1._8)+nu)*((-1._8)+2._8*nu)*lr2*(x1p-y1)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*(x2p-y2)*(x3+y3)*((x2p-y2)**2+(x3+ &
            y3)**2)**(-1)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-1)*(x1p- &
            y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*( &
            x3+y3)**3*((x2p-y2)**2+(x3+y3)**2)**(-1)-4._8*((-1._8)+nu)* &
            lr2*(x1p-y1)*(x2p-y2)*(x3+y3)*((x1p-y1)**2+(x3+ &
            y3)**2)**(-1)*((x2p-y2)**2+(x3+y3)**2)**(-1)-4._8*((-1._8)+nu) &
            *lr2**(-1)*(x1p-y1)*(x2p-y2)*(x3+y3)**3*((x1p- &
            y1)**2+(x3+y3)**2)**(-1)*((x2p-y2)**2+(x3+y3)**2)**(-1) &
            -2._8*x3*(x1p-y1)*(x2p-y2)*y3*(x3+y3)*((x2p-y2) &
            **2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3) &
            **2)**(-3._8/2._8)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*atan((x1p-y1)*( &
            x2p-y2)**(-1))-atan2(x3,x1p-y1)-3._8*atan2(3*x3, &
            x1p-y1)+4._8*nu*atan2(-nu*x3,x1p-y1)+4._8*((-1._8)+nu)* &
            ((-1._8)+2._8*nu)*atan2(lr2*(-x2p+y2),(x1p-y1)*(x3+y3))+(4._8+( &
            -4._8)*nu)*atan2(lr1*(x3-y3),(x1p-y1)*(x2p-y2))+( &
            4._8-4._8*nu)*atan2(lr2*(x3+y3),(x1p-y1)*(x2p-y2)))

    END FUNCTION J2212d3

    !---------------------------------------------------------------
    !> function J2213d1
    !! computes the derivative of the J integral J2213,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2213d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2213d1=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(x2p**2*(x2p**2+(x1p+ &
            (-1)*y1)**2)**(-1)+9._8*x2p**2*(9._8*x2p**2+(x1p-y1)**2)**(-1)+ &
            4._8*nu**2*x2p**2*(nu**2*x2p**2+(x1p-y1)**2)**(-1)+2._8*x3*(lr2+ &
            (-1)*x1p+y1)**(-1)+2._8*lr2**(-1)*x3*(-x1p+y1)*(lr2-x1p+ &
            y1)**(-1)-2._8*((-1._8)+nu)*((-1._8)+2._8*nu)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*(x2p-y2)**2+y2**2*((x1p-y1)**2+ &
            y2**2)**(-1)+9._8*y2**2*((x1p-y1)**2+9._8*y2**2)**(-1)+4._8* &
            nu**2*y2**2*((x1p-y1)**2+nu**2*y2**2)**(-1)-((-3._8)+ &
            4._8*nu)*(lr1+x1p-y1)**(-1)*(x3-y3)-((-3._8)+4._8*nu)* &
            lr1**(-1)*(x1p-y1)*(lr1+x1p-y1)**(-1)*(x3-y3)+2._8* &
            ((-1._8)+2._8*nu)*lr1*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+( &
            -1)*y2)**2*((x2p-y2)**2+(x3-y3)**2)**(-1)*(x3- &
            y3)-2*((-1._8)+2._8*nu)*lr1**(-1)*(x1p-y1)**2*((x1p-y1) &
            **2+(x2p-y2)**2)**(-1)*(x2p-y2)**2*((x2p-y2)**2+( &
            x3-y3)**2)**(-1)*(x3-y3)-((-3._8)+4._8*nu)*lr1**(-1) &
            *(x1p-y1)**2*(lr1+x3-y3)**(-1)-(5._8+4._8*nu*((-3._8)+ &
            2._8*nu))*(lr2+x1p-y1)**(-1)*(x3+y3)-(5._8+4._8*nu*((-3._8)+2._8* &
            nu))*lr2**(-1)*(x1p-y1)*(lr2+x1p-y1)**(-1)*(x3+y3)+( &
            -1)*(3._8-6._8*nu+4._8*nu**2)*lr2**(-1)*(x1p-y1)**2*(lr2+x3+y3) &
            **(-1)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*(x1p-y1)**2*((x1p- &
            y1)**2+(x2p-y2)**2)**(-2._8)*y3*(2._8*x3+y3)-2._8*((-1._8)+nu)*(( &
            -1._8)+2._8*nu)*lr2**(-2)*(x1p-y1)**2*((x1p-y1)**2+(x2p+(-1._8) &
            *y2)**2)**(-1)*y3*(2*x3+y3)+2._8*((-1._8)+2._8*nu)*((x1p-y1) &
            **2+(x2p-y2)**2)**(-1)*((x2p-y2)**2+(x3+y3)**2)**(-1)* &
            (((-1._8)+2._8*nu)*lr2*(x2p-y2)**2*(x3+y3)-((-1._8)+nu)*y3* &
            (2._8*x3+y3)*((x2p-y2)**2+(x3+y3)**2))-4._8*lr2**(-1)*(x1p+( &
            -1)*y1)**2*((x1p-y1)**2+(x2p-y2)**2)**(-2)*((x2p+(-1) &
            *y2)**2+(x3+y3)**2)**(-1)*(x2p**4*y3+4*x2p**2*x3**2*y3+2._8* &
            x3**4*y3+x2p**2*y1**2*y3+2*x3**2*y1**2*y3-4*x2p**3*y2* &
            y3-8._8*x2p*x3**2*y2*y3-2._8*x2p*y1**2*y2*y3+6._8*x2p**2* &
            y2**2*y3+4._8*x3**2*y2**2*y3+y1**2*y2**2*y3-4._8*x2p*y2**3* &
            y3+y2**4*y3+6._8*x2p**2*x3*y3**2+7._8*x3**3*y3**2+3._8*x3*y1**2* &
            y3**2-12._8*x2p*x3*y2*y3**2+6._8*x3*y2**2*y3**2+2._8*x2p**2* &
            y3**3+9._8*x3**2*y3**3+y1**2*y3**3-4._8*x2p*y2*y3**3+2._8*y2**2* &
            y3**3+5._8*x3*y3**4+y3**5-3._8*nu*(x3+y3)*((x2p-y2)**2+(x3+ &
            y3)**2)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)+2._8*nu**2* &
            (x3+y3)*((x2p-y2)**2+(x3+y3)**2)*((x1p-y1)**2+(x2p+(-1._8) &
            *y2)**2+(x3+y3)**2)+x1p**2*y3*((x2p-y2)**2+(x3+y3)*(2._8*x3+ &
            y3))-2._8*x1p*y1*y3*((x2p-y2)**2+(x3+y3)*(2._8*x3+y3)))+( &
            -2._8)*(x1p-y1)**2*((x1p-y1)**2+(x2p-y2)**2)**(-1)* &
            ((x2p-y2)**2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p- &
            y2)**2+(x3+y3)**2)**(-3._8/2._8)*(x2p**4*y3-2._8*lr2*x2p**2*x3*y3+4._8* &
            x2p**2*x3**2*y3-2._8*lr2*x3**3*y3+2*x3**4*y3+x2p**2*y1**2* &
            y3+2._8*x3**2*y1**2*y3-4._8*x2p**3*y2*y3+4*lr2*x2p*x3*y2*y3+( &
            -8)*x2p*x3**2*y2*y3-2*x2p*y1**2*y2*y3+6._8*x2p**2*y2**2* &
            y3-2._8*lr2*x3*y2**2*y3+4._8*x3**2*y2**2*y3+y1**2*y2**2*y3+( &
            -4._8)*x2p*y2**3*y3+y2**4*y3-lr2*x2p**2*y3**2+6._8*x2p**2*x3* &
            y3**2-5._8*lr2*x3**2*y3**2+7._8*x3**3*y3**2+3._8*x3*y1**2*y3**2+ &
            2._8*lr2*x2p*y2*y3**2-12._8*x2p*x3*y2*y3**2-lr2*y2**2* &
            y3**2+6._8*x3*y2**2*y3**2+2._8*x2p**2*y3**3-4._8*lr2*x3*y3**3+9._8* &
            x3**2*y3**3+y1**2*y3**3-4._8*x2p*y2*y3**3+2._8*y2**2*y3**3+(-1) &
            *lr2*y3**4+5._8*x3*y3**4+y3**5-3._8*nu*(x3*(x3**2+y1**2+(x2p+( &
            -1)*y2)**2)+((-2._8)*lr2*x3+3*x3**2+y1**2+(x2p-y2)**2)*y3+( &
            -1)*(lr2-3._8*x3)*y3**2+y3**3)*((x2p-y2)**2+(x3+y3)**2)+ &
            2._8*nu**2*(x3*(x3**2+y1**2+(x2p-y2)**2)+((-2._8)*lr2*x3+3* &
            x3**2+y1**2+(x2p-y2)**2)*y3-(lr2-3*x3)*y3**2+ &
            y3**3)*((x2p-y2)**2+(x3+y3)**2)+2._8*x1p*y1*(3._8*nu*(x3+y3)* &
            ((x2p-y2)**2+(x3+y3)**2)-2._8*nu**2*(x3+y3)*((x2p-y2) &
            **2+(x3+y3)**2)-y3*((x2p-y2)**2+(x3+y3)*(2._8*x3+y3)))+ &
            x1p**2*((-3)*nu*(x3+y3)*((x2p-y2)**2+(x3+y3)**2)+2* &
            nu**2*(x3+y3)*((x2p-y2)**2+(x3+y3)**2)+y3*((x2p-y2) &
            **2+(x3+y3)*(2._8*x3+y3))))+2*lr2**(-1)*((x1p-y1)**2+(x2p+(-1) &
            *y2)**2)**(-1)*((x2p-y2)**2+(x3+y3)**2)**(-1)*(-x3* &
            y1**2*y2**2+2*x3**4*y3+6._8*x3**2*y1**2*y3+4._8*x3**2*y2**2*y3+ &
            2._8*y1**2*y2**2*y3+y2**4*y3+7._8*x3**3*y3**2+9._8*x3*y1**2*y3**2+ &
            6._8*x3*y2**2*y3**2+9._8*x3**2*y3**3+3*y1**2*y3**3+2*y2**2* &
            y3**3+5._8*x3*y3**4+y3**5+x2p**4*(y3-3._8*nu*(x3+y3)+2._8*nu**2*( &
            x3+y3))-4._8*x2p**3*y2*(y3-3*nu*(x3+y3)+2._8*nu**2*(x3+y3))+ &
            2*x2p*y2*(x3*y1**2-4._8*x3**2*y3-2*y1**2*y3-2* &
            y2**2*y3-6._8*x3*y3**2-2*y3**3-2._8*nu**2*(x3+y3)*( &
            y1**2+2*y2**2+2._8*(x3+y3)**2)+nu*(x3+y3)*(5._8*y1**2+6._8*y2**2+6._8*( &
            x3+y3)**2))+x2p**2*(-x3*y1**2+4._8*x3**2*y3+2*y1**2*y3+6._8* &
            y2**2*y3+6._8*x3*y3**2+2._8*y3**3+2._8*nu**2*(x3+y3)*(y1**2+6._8* &
            y2**2+2*(x3+y3)**2)-nu*(x3+y3)*(5._8*y1**2+18._8*y2**2+6._8*( &
            x3+y3)**2))+x1p**2*(-x2p**2*x3+2._8*x2p*x3*y2-x3*y2**2+ &
            2._8*x2p**2*y3+6._8*x3**2*y3-4._8*x2p*y2*y3+2._8*y2**2*y3+9._8*x3* &
            y3**2+3._8*y3**3+2._8*nu**2*(x3+y3)*((x2p-y2)**2+3._8*(x3+y3)**2) &
            -nu*(x3+y3)*(5._8*(x2p-y2)**2+9._8*(x3+y3)**2))+2._8*x1p* &
            y1*(x2p**2*x3-2._8*x2p*x3*y2+x3*y2**2-2._8*x2p**2*y3-6._8* &
            x3**2*y3+4._8*x2p*y2*y3-2._8*y2**2*y3-9._8*x3*y3**2-3._8* &
            y3**3-2*nu**2*(x3+y3)*((x2p-y2)**2+3*(x3+y3)**2)+nu*( &
            x3+y3)*(5*(x2p-y2)**2+9._8*(x3+y3)**2))-nu*(x3+y3)*( &
            3*x3**4+12*x3**3*y3+3*(y2**2+y3**2)**2+3*x3**2*(3*y1**2+2* &
            y2**2+6*y3**2)+y1**2*(5._8*y2**2+9._8*y3**2)+6._8*x3*y3*(3._8*y1**2+ &
            2*(y2**2+y3**2)))+2*nu**2*(x3+y3)*(x3**4+4*x3**3*y3+(y2**2+ &
            y3**2)**2+y1**2*(y2**2+3*y3**2)+x3**2*(3._8*y1**2+2._8*y2**2+6._8* &
            y3**2)+x3*(6._8*y1**2*y3+4._8*y3*(y2**2+y3**2))))+xLogy(3._8-4._8*nu, &
            lr1+x3-y3)+xLogy((-3._8)+6._8*nu-4._8*nu**2,lr2+x3+y3))

    END FUNCTION J2213d1

    !---------------------------------------------------------------
    !> function J2213d2
    !! computes the derivative of the J integral J2213,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2213d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2213d2=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(x2p*(x2p**2+(x1p+( &
            -1)*y1)**2)**(-1)*(-x1p+y1)+9._8*x2p*(9._8*x2p**2+(x1p-y1) &
            **2)**(-1)*(-x1p+y1)+4._8*nu**2*x2p*(nu**2*x2p**2+(x1p- &
            y1)**2)**(-1)*(-x1p+y1)+2._8*((-1._8)+nu)*((-1._8)+2._8*nu)*(x1p+(-1) &
            *y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)+ &
            2._8*lr2**(-1)*x3*(lr2-x1p+y1)**(-1)*(-x2p+y2)-((-3._8) &
            +4._8*nu)*lr1**(-1)*(lr1+x1p-y1)**(-1)*(x2p-y2)*(x3+(-1) &
            *y3)-2._8*((-1._8)+2._8*nu)*lr1*(x1p-y1)*((x1p-y1)**2+( &
            x2p-y2)**2)**(-1)*(x2p-y2)*((x2p-y2)**2+(x3+(-1) &
            *y3)**2)**(-1)*(x3-y3)-2._8*((-1._8)+2._8*nu)*lr1**(-1)*(x1p+( &
            -1)*y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2) &
            **3*((x2p-y2)**2+(x3-y3)**2)**(-1)*(x3-y3)+(-1) &
            *((-3._8)+4._8*nu)*lr1**(-1)*(x1p-y1)*(x2p-y2)*(lr1+x3+( &
            -1)*y3)**(-1)-(5._8+4._8*nu*((-3._8)+2._8*nu))*lr2**(-1)*(lr2+x1p+( &
            -1)*y1)**(-1)*(x2p-y2)*(x3+y3)-(3._8-6._8*nu+4._8*nu**2) &
            *lr2**(-1)*(x1p-y1)*(x2p-y2)*(lr2+x3+y3)**(-1)+4._8*(( &
            -1._8)+nu)*((-1._8)+2._8*nu)*(x1p-y1)*((x1p-y1)**2+(x2p- &
            y2)**2)**(-2)*(x2p-y2)*y3*(2*x3+y3)-2._8*((-1._8)+nu)*(( &
            -1)+2._8*nu)*lr2**(-2)*(x1p-y1)*((x1p-y1)**2+(x2p- &
            y2)**2)**(-1)*(x2p-y2)*y3*(2._8*x3+y3)+4._8*((-1._8)+nu)*((-1._8)+ &
            2._8*nu)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)* &
            (x2p-y2)*y3*(2._8*x3+y3)*((x2p-y2)**2+(x3+y3)**2)**(-1) &
            -2._8*((-1._8)+2._8*nu)*(x1p-y1)*((x1p-y1)**2+(x2p- &
            y2)**2)**(-1)*(x2p-y2)*((x2p-y2)**2+(x3+y3)**2)**(-1) &
            *(((-1._8)+2._8*nu)*lr2*(x3+y3)+2*((-1._8)+nu)*y3*(2._8*x3+y3))+2* &
            lr2**(-1)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1) &
            *(x2p-y2)*((x2p-y2)**2+(x3+y3)**2)**(-1)*(- &
            x2p**2*x3+2._8*x2p*x3*y2-x3*y2**2+2._8*x1p**2*y3+3._8*x2p**2*y3+ &
            8._8*x3**2*y3-4._8*x1p*y1*y3+2._8*y1**2*y3-6._8*x2p*y2*y3+3._8* &
            y2**2*y3+12._8*x3*y3**2+4._8*y3**3+4._8*nu**2*(x3+y3)*((x1p-y1) &
            **2+(x2p-y2)**2+2._8*(x3+y3)**2)-2._8*nu*(x3+y3)*(3._8*(x1p+( &
            -1._8)*y1)**2+4._8*(x2p-y2)**2+6._8*(x3+y3)**2))-4._8*lr2**(-1)*( &
            x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+(-1) &
            *y2)*((x2p-y2)**2+(x3+y3)**2)**(-2._8)*(x2p**4*y3+4._8*x2p**2* &
            x3**2*y3+2._8*x3**4*y3+x2p**2*y1**2*y3+2*x3**2*y1**2*y3-4._8* &
            x2p**3*y2*y3-8._8*x2p*x3**2*y2*y3-2._8*x2p*y1**2*y2*y3+6._8* &
            x2p**2*y2**2*y3+4._8*x3**2*y2**2._8*y3+y1**2*y2**2*y3-4._8*x2p* &
            y2**3*y3+y2**4*y3+6._8*x2p**2*x3*y3**2+7._8*x3**3*y3**2+3._8*x3* &
            y1**2*y3**2-12._8*x2p*x3*y2*y3**2+6._8*x3*y2**2*y3**2+2._8* &
            x2p**2*y3**3+9._8*x3**2*y3**3+y1**2*y3**3-4._8*x2p*y2*y3**3+2._8* &
            y2**2*y3**3+5._8*x3*y3**4+y3**5-3._8*nu*(x3+y3)*((x2p-y2) &
            **2+(x3+y3)**2)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)+ &
            2._8*nu**2*(x3+y3)*((x2p-y2)**2+(x3+y3)**2)*((x1p-y1) &
            **2+(x2p-y2)**2+(x3+y3)**2)+x1p**2*y3*((x2p-y2)**2+(x3+ &
            y3)*(2*x3+y3))-2*x1p*y1*y3*((x2p-y2)**2+(x3+y3)*(2* &
            x3+y3)))-4._8*lr2**(-1)*(x1p-y1)*((x1p-y1)**2+(x2p+(-1) &
            *y2)**2)**(-2)*(x2p-y2)*((x2p-y2)**2+(x3+y3)**2)**( &
            -1)*(x2p**4*y3+4._8*x2p**2*x3**2*y3+2._8*x3**4*y3+x2p**2*y1**2*y3+ &
            2._8*x3**2*y1**2*y3-4._8*x2p**3*y2*y3-8._8*x2p*x3**2*y2*y3+( &
            -2)*x2p*y1**2*y2*y3+6._8*x2p**2*y2**2*y3+4._8*x3**2*y2**2*y3+ &
            y1**2*y2**2*y3-4*x2p*y2**3*y3+y2**4*y3+6*x2p**2*x3* &
            y3**2+7._8*x3**3*y3**2+3*x3*y1**2*y3**2-12._8*x2p*x3*y2* &
            y3**2+6._8*x3*y2**2*y3**2+2._8*x2p**2*y3**3+9._8*x3**2*y3**3+y1**2* &
            y3**3-4._8*x2p*y2*y3**3+2._8*y2**2*y3**3+5._8*x3*y3**4+y3**5 &
            -3._8*nu*(x3+y3)*((x2p-y2)**2+(x3+y3)**2)*((x1p-y1)**2+( &
            x2p-y2)**2+(x3+y3)**2)+2._8*nu**2*(x3+y3)*((x2p-y2)**2+( &
            x3+y3)**2)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)+x1p**2* &
            y3*((x2p-y2)**2+(x3+y3)*(2*x3+y3))-2._8*x1p*y1*y3*((x2p+ &
            (-1)*y2)**2+(x3+y3)*(2*x3+y3)))-2._8*(x1p-y1)*((x1p+(-1) &
            *y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*((x2p-y2) &
            **2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3) &
            **2)**(-3._8/2._8)*(x2p**4*y3-2._8*lr2*x2p**2*x3*y3+4._8*x2p**2*x3**2._8* &
            y3-2._8*lr2*x3**3*y3+2._8*x3**4*y3+x2p**2*y1**2*y3+2._8*x3**2* &
            y1**2*y3-4*x2p**3*y2*y3+4*lr2*x2p*x3*y2*y3-8*x2p* &
            x3**2*y2*y3-2._8*x2p*y1**2*y2*y3+6._8*x2p**2*y2**2*y3-2._8* &
            lr2*x3*y2**2*y3+4._8*x3**2*y2**2*y3+y1**2*y2**2*y3-4._8*x2p* &
            y2**3*y3+y2**4*y3-lr2*x2p**2*y3**2+6._8*x2p**2*x3*y3**2+( &
            -5._8)*lr2*x3**2*y3**2+7._8*x3**3*y3**2+3._8*x3*y1**2*y3**2+2*lr2* &
            x2p*y2*y3**2-12._8*x2p*x3*y2*y3**2-lr2*y2**2*y3**2+6._8* &
            x3*y2**2*y3**2+2._8*x2p**2*y3**3-4._8*lr2*x3*y3**3+9._8*x3**2* &
            y3**3+y1**2*y3**3-4._8*x2p*y2*y3**3+2._8*y2**2*y3**3-lr2* &
            y3**4+5._8*x3*y3**4+y3**5-3._8*nu*(x3*(x3**2+y1**2+(x2p-y2) &
            **2)+((-2)*lr2*x3+3*x3**2+y1**2+(x2p-y2)**2)*y3-(lr2+ &
            (-3._8)*x3)*y3**2+y3**3)*((x2p-y2)**2+(x3+y3)**2)+2._8*nu**2*( &
            x3*(x3**2+y1**2+(x2p-y2)**2)+((-2)*lr2*x3+3._8*x3**2+y1**2+( &
            x2p-y2)**2)*y3-(lr2-3._8*x3)*y3**2+y3**3)*((x2p+(-1) &
            *y2)**2+(x3+y3)**2)+2._8*x1p*y1*(3._8*nu*(x3+y3)*((x2p-y2) &
            **2+(x3+y3)**2)-2._8*nu**2*(x3+y3)*((x2p-y2)**2+(x3+y3) &
            **2)-y3*((x2p-y2)**2+(x3+y3)*(2*x3+y3)))+x1p**2*(( &
            -3)*nu*(x3+y3)*((x2p-y2)**2+(x3+y3)**2)+2*nu**2*(x3+y3) &
            *((x2p-y2)**2+(x3+y3)**2)+y3*((x2p-y2)**2+(x3+y3)*( &
            2*x3+y3))))+2._8*((-1._8)+nu)*((-1._8)+2._8*nu)*atan((x1p-y1)**(-1) &
            *(x2p-y2))+atan2(x1p-y1,(-1)*x2p)-3*atan2(3*x2p,x1p+( &
            -1)*y1)+4._8*nu*atan2(-nu*x2p,x1p-y1)+((-2._8)+4._8*nu)* &
            atan2(lr1*(-x2p+y2),(x1p-y1)*(x3-y3))+2._8*(1._8-2._8* &
            nu)**2*atan2(lr2*(-x2p+y2),(x1p-y1)*(x3+y3)))

    END FUNCTION J2213d2

    !---------------------------------------------------------------
    !> function J2213d3
    !! computes the derivative of the J integral J2213,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2213d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2213d3=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(2._8*((-1._8)+2._8*nu)* &
            lr1*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+ &
            (-1)*y2)**2*((x2p-y2)**2+(x3-y3)**2)**(-1)-(( &
            -3._8)+4._8*nu)*lr1**(-1)*(lr1+x1p-y1)**(-1)*(x3-y3)**2+( &
            -2._8)*((-1._8)+2._8*nu)*lr1**(-1)*(x1p-y1)*((x1p-y1)**2+(x2p+ &
            (-1)*y2)**2)**(-1)*(x2p-y2)**2*((x2p-y2)**2+(x3+(-1) &
            *y3)**2)**(-1)*(x3-y3)**2-((-3._8)+4._8*nu)*(x1p- &
            y1)*(lr1+x3-y3)**(-1)-((-3._8)+4._8*nu)*lr1**(-1)*(x1p+(-1) &
            *y1)*(x3-y3)*(lr1+x3-y3)**(-1)-2._8*lr2**(-1)*x3*( &
            lr2-x1p+y1)**(-1)*(x3+y3)-(5+4*nu*((-3._8)+2._8*nu))* &
            lr2**(-1)*(lr2+x1p-y1)**(-1)*(x3+y3)**2-(3._8-6._8*nu+4._8* &
            nu**2)*(x1p-y1)*(lr2+x3+y3)**(-1)-(3._8-6._8*nu+4._8* &
            nu**2)*lr2**(-1)*(x1p-y1)*(x3+y3)*(lr2+x3+y3)**(-1)-2._8*( &
            (-1._8)+nu)*((-1._8)+2*nu)*lr2**(-2)*(x1p-y1)*((x1p-y1) &
            **2+(x2p-y2)**2)**(-1)*y3*(x3+y3)*(2._8*x3+y3)+4._8*((-1._8)+nu) &
            *((-1._8)+2._8*nu)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2) &
            **(-1)*y3*(x3+y3)*(2._8*x3+y3)*((x2p-y2)**2+(x3+y3)**2)**( &
            -1)-4._8*lr2**(-1)*(x1p-y1)*((x1p-y1)**2+(x2p-y2) &
            **2)**(-1)*(x3+y3)*((x2p-y2)**2+(x3+y3)**2)**(-2)*(x2p**4* &
            y3+4._8*x2p**2*x3**2*y3+2._8*x3**4*y3+x2p**2*y1**2*y3+2._8*x3**2* &
            y1**2*y3-4._8*x2p**3*y2*y3-8._8*x2p*x3**2*y2*y3-2._8*x2p* &
            y1**2*y2*y3+6._8*x2p**2*y2**2*y3+4._8*x3**2*y2**2*y3+y1**2* &
            y2**2*y3-4._8*x2p*y2**3*y3+y2**4*y3+6._8*x2p**2*x3*y3**2+7* &
            x3**3*y3**2+3._8*x3*y1**2*y3**2-12*x2p*x3*y2*y3**2+6._8*x3* &
            y2**2*y3**2+2._8*x2p**2*y3**3+9._8*x3**2*y3**3+y1**2*y3**3-4._8* &
            x2p*y2*y3**3+2._8*y2**2*y3**3+5._8*x3*y3**4+y3**5-3._8*nu*(x3+y3) &
            *((x2p-y2)**2+(x3+y3)**2)*((x1p-y1)**2+(x2p-y2) &
            **2+(x3+y3)**2)+2*nu**2*(x3+y3)*((x2p-y2)**2+(x3+y3)**2)* &
            ((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)+x1p**2*y3*((x2p+( &
            -1)*y2)**2+(x3+y3)*(2*x3+y3))-2._8*x1p*y1*y3*((x2p-y2) &
            **2+(x3+y3)*(2._8*x3+y3)))+2._8*((-1._8)+2._8*nu)*(x1p-y1)*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2)**(-1)*((x2p-y2)**2+(x3+y3)**2) &
            **(-1)*(((-1._8)+2._8*nu)*lr2*(x2p-y2)**2-2._8*((-1._8)+nu)*y3*( &
            (x2p-y2)**2+(x3+y3)*(3._8*x3+2._8*y3)))+2._8*lr2**(-1)*(x1p- &
            y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*((x2p-y2)**2+ &
            (x3+y3)**2)**(-1)*(nu*((-3._8)+2._8*nu)*x2p**4+4._8*(3._8-2._8*nu)*nu* &
            x2p**3*y2-x3**2*y2**2+4._8*x1p**2*x3*y3+8._8*x3**3*y3-8._8* &
            x1p*x3*y1*y3+4._8*x3*y1**2*y3+6._8*x3*y2**2*y3+3._8*x1p**2*y3**2+ &
            21._8*x3**2*y3**2-6._8*x1p*y1*y3**2+3._8*y1**2*y3**2+5._8*y2**2* &
            y3**2+18._8*x3*y3**3+5._8*y3**4+2._8*x2p*y2*((1._8+2._8*(7._8-4._8*nu)*nu)* &
            x3**2+(3._8-2._8*nu)*nu*((x1p-y1)**2+2._8*y2**2)-2*(3+2* &
            nu*((-7._8)+4._8*nu))*x3*y3+((-5._8)+2._8*(7._8-4._8*nu)*nu)*y3**2)+ &
            x2p**2*(((-1._8)+2._8*nu*((-7._8)+4._8*nu))*x3**2+nu*((-3._8)+2._8*nu)*((x1p+( &
            -1)*y1)**2+6._8*y2**2)+2*(3._8+2._8*nu*((-7._8)+4._8*nu))*x3*y3+(5._8+2._8* &
            nu*((-7._8)+4._8*nu))*y3**2)-nu*(y2**2+3._8*(x3+y3)**2)*(3._8*( &
            x1p-y1)**2+3._8*y2**2+5._8*(x3+y3)**2)+2._8*nu**2*(5._8*x3**4+ &
            y1**2*y2**2+y2**4+20._8*x3**3*y3+3._8*y1**2*y3**2+4._8*y2**2*y3**2+ &
            5._8*y3**4+x3**2*(3._8*y1**2+4._8*y2**2+30._8*y3**2)+x3*(6._8*y1**2*y3+ &
            8._8*y2**2*y3+20._8*y3**3)+x1p**2*(y2**2+3._8*(x3+y3)**2)-2._8*x1p* &
            y1*(y2**2+3._8*(x3+y3)**2)))-2._8*(x1p-y1)*((x1p-y1) &
            **2+(x2p-y2)**2)**(-1)*(x3+y3)*((x2p-y2)**2+(x3+y3) &
            **2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8) &
            *(x2p**4*y3-2._8*lr2*x2p**2*x3*y3+4._8*x2p**2*x3**2*y3-2._8* &
            lr2*x3**3*y3+2*x3**4*y3+x2p**2*y1**2*y3+2*x3**2*y1**2*y3+( &
            -4._8)*x2p**3*y2*y3+4._8*lr2*x2p*x3*y2*y3-8._8*x2p*x3**2*y2*y3+( &
            -2._8)*x2p*y1**2*y2*y3+6._8*x2p**2*y2**2*y3-2._8*lr2*x3*y2**2* &
            y3+4._8*x3**2*y2**2*y3+y1**2*y2**2*y3-4._8*x2p*y2**3*y3+ &
            y2**4*y3-lr2*x2p**2*y3**2+6._8*x2p**2*x3*y3**2-5._8*lr2* &
            x3**2*y3**2+7*x3**3*y3**2+3._8*x3*y1**2*y3**2+2._8*lr2*x2p*y2* &
            y3**2-12._8*x2p*x3*y2*y3**2-lr2*y2**2*y3**2+6._8*x3* &
            y2**2*y3**2+2._8*x2p**2*y3**3-4._8*lr2*x3*y3**3+9._8*x3**2*y3**3+ &
            y1**2*y3**3-4._8*x2p*y2*y3**3+2._8*y2**2*y3**3-lr2*y3**4+ &
            5._8*x3*y3**4+y3**5-3._8*nu*(x3*(x3**2+y1**2+(x2p-y2)**2)+( &
            (-2._8)*lr2*x3+3*x3**2+y1**2+(x2p-y2)**2)*y3-(lr2-3._8* &
            x3)*y3**2+y3**3)*((x2p-y2)**2+(x3+y3)**2)+2._8*nu**2*(x3*( &
            x3**2+y1**2+(x2p-y2)**2)+((-2._8)*lr2*x3+3*x3**2+y1**2+(x2p+( &
            -1)*y2)**2)*y3-(lr2-3._8*x3)*y3**2+y3**3)*((x2p-y2) &
            **2+(x3+y3)**2)+2._8*x1p*y1*(3._8*nu*(x3+y3)*((x2p-y2)**2+(x3+ &
            y3)**2)-2*nu**2*(x3+y3)*((x2p-y2)**2+(x3+y3)**2)- &
            y3*((x2p-y2)**2+(x3+y3)*(2*x3+y3)))+x1p**2*((-3._8)*nu*(x3+ &
            y3)*((x2p-y2)**2+(x3+y3)**2)+2._8*nu**2*(x3+y3)*((x2p- &
            y2)**2+(x3+y3)**2)+y3*((x2p-y2)**2+(x3+y3)*(2._8*x3+y3))))+ &
            xLogy((-2._8),lr2-x1p+y1)+xLogy(3._8-4._8*nu,lr1+x1p-y1)+xLogy( &
            (-5._8)+4._8*(3._8-2._8*nu)*nu,lr2+x1p-y1))

    END FUNCTION J2213d3

    !---------------------------------------------------------------
    !> function J2223d1
    !! computes the derivative of the J integral J2223,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2223d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2223d1=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(x1p*(x1p**2+(x2p+( &
            -1)*y2)**2)**(-1)*(-x2p+y2)+9._8*x1p*(9._8*x1p**2+(x2p-y2) &
            **2)**(-1)*(-x2p+y2)+4._8*nu**2*x1p*(nu**2*x1p**2+(x2p- &
            y2)**2)**(-1)*(-x2p+y2)-4._8*((-1._8)+nu)*lr1**(-1)*(x1p- &
            y1)*(lr1+x2p-y2)**(-1)*(x3-y3)-4._8*((-1._8)+nu)*lr1*( &
            x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+(-1) &
            *y2)*((x1p-y1)**2+(x3-y3)**2)**(-1)*(x3-y3)+( &
            -4._8)*((-1._8)+nu)*lr1**(-1)*(x1p-y1)**3*((x1p-y1)**2+(x2p+ &
            (-1)*y2)**2)**(-1)*(x2p-y2)*((x1p-y1)**2+(x3- &
            y3)**2)**(-1)*(x3-y3)-((-3._8)+4._8*nu)*lr1**(-1)*(x1p+( &
            -1)*y1)*(x2p-y2)*(lr1+x3-y3)**(-1)+4._8*((-1._8)+nu)* &
            lr2**(-1)*(x1p-y1)*(lr2+x2p-y2)**(-1)*(x3+y3)-(3._8+ &
            (-6._8)*nu+4._8*nu**2)*lr2**(-1)*(x1p-y1)*(x2p-y2)*(lr2+ &
            x3+y3)**(-1)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*(x1p-y1)*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2)**(-2)*(x2p-y2)*y3*(2*x3+y3) &
            +2*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-2)*(x1p-y1)*((x1p- &
            y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*y3*(2*x3+y3) &
            -4._8*lr2**(-1)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**( &
            -1)*(x2p-y2)*(y3-3._8*nu*(x3+y3)+2*nu**2*(x3+y3))+4._8*(( &
            -1._8)+nu)*lr2*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**( &
            -1)*(x2p-y2)*(x3+y3)*((x1p-y1)**2+(x3+y3)**2)**(-1)+ &
            4._8*((-1._8)+nu)*lr2**(-1)*(x1p-y1)**3*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*(x2p-y2)*(x3+y3)*((x1p-y1)**2+(x3+ &
            y3)**2)**(-1)+2*(x1p-y1)*((x1p-y1)**2+(x2p-y2) &
            **2)**(-1)*(x2p-y2)*((x1p-y1)**2+(x2p-y2)**2+(x3+ &
            y3)**2)**(-3._8/2._8)*((-3._8)*nu*(x3*(x3**2+(x1p-y1)**2+(x2p+(-1) &
            *y2)**2)+(x3*((-2._8)*lr2+3*x3)+(x1p-y1)**2+(x2p-y2)**2) &
            *y3-(lr2-3*x3)*y3**2+y3**3)+2*nu**2*(x3*(x3**2+(x1p+( &
            -1)*y1)**2+(x2p-y2)**2)+(x3*((-2)*lr2+3*x3)+(x1p-y1) &
            **2+(x2p-y2)**2)*y3-(lr2-3*x3)*y3**2+y3**3)+y3*(( &
            x1p-y1)**2+(x2p-y2)**2-(lr2-x3-y3)*( &
            2._8*x3+y3)))+4._8*lr2**(-1)*(x1p-y1)*((x1p-y1)**2+(x2p+(-1) &
            *y2)**2)**(-2)*(x2p-y2)*((-3._8)*nu*(x3+y3)*((x1p-y1) &
            **2+(x2p-y2)**2+(x3+y3)**2)+2*nu**2*(x3+y3)*((x1p-y1) &
            **2+(x2p-y2)**2+(x3+y3)**2)+y3*((x1p-y1)**2+(x2p- &
            y2)**2+(x3+y3)*(2*x3+y3)))+atan2(x2p-y2,(-1)*x1p)-3._8* &
            atan2(3._8*x1p,x2p-y2)+4._8*nu*atan2(-nu*x1p,x2p-y2)+( &
            4._8-4._8*nu)*atan2(lr1*(x1p-y1),(x2p-y2)*(x3-y3)) &
            +4._8*((-1._8)+nu)*atan2(lr2*(x1p-y1),(x2p-y2)*(x3+y3)))

    END FUNCTION J2223d1

    !---------------------------------------------------------------
    !> function J2223d2
    !! computes the derivative of the J integral J2223,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2223d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2223d2=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(x1p**2*(x1p**2+(x2p+ &
            (-1)*y2)**2)**(-1)+9._8*x1p**2*(9._8*x1p**2+(x2p-y2)**2)**(-1)+ &
            4*nu**2*x1p**2*(nu**2*x1p**2+(x2p-y2)**2)**(-1)+y1**2*( &
            y1**2+(x2p-y2)**2)**(-1)+9._8*y1**2*(9._8*y1**2+(x2p-y2) &
            **2)**(-1)+4._8*nu**2*y1**2*(nu**2*y1**2+(x2p-y2)**2)**(-1)+ &
            (-4._8)*((-1._8)+nu)*(lr1+x2p-y2)**(-1)*(x3-y3)-4._8*((-1._8) &
            +nu)*lr1**(-1)*(x2p-y2)*(lr1+x2p-y2)**(-1)*(x3- &
            y3)+4._8*((-1._8)+nu)*lr1*(x1p-y1)**2*((x1p-y1)**2+(x2p+(-1) &
            *y2)**2)**(-1)*((x1p-y1)**2+(x3-y3)**2)**(-1)*(x3+( &
            -1)*y3)-4._8*((-1._8)+nu)*lr1**(-1)*(x1p-y1)**2*((x1p- &
            y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)**2*((x1p-y1) &
            **2+(x3-y3)**2)**(-1)*(x3-y3)-((-3)+4._8*nu)* &
            lr1**(-1)*(x2p-y2)**2*(lr1+x3-y3)**(-1)+4._8*((-1._8)+nu)*( &
            lr2+x2p-y2)**(-1)*(x3+y3)+4._8*((-1._8)+nu)*lr2**(-1)*(x2p- &
            y2)*(lr2+x2p-y2)**(-1)*(x3+y3)-(3._8-6._8*nu+4*nu**2)* &
            lr2**(-1)*(x2p-y2)**2*(lr2+x3+y3)**(-1)+2._8*((-1._8)+nu)*((-1._8)+ &
            2._8*nu)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*y3*(2*x3+y3)+ &
            (-4._8)*((-1._8)+nu)*((-1._8)+2._8*nu)*((x1p-y1)**2+(x2p-y2)**2) &
            **(-2)*(x2p-y2)**2*y3*(2._8*x3+y3)+2*((-1._8)+nu)*((-1._8)+2._8* &
            nu)*lr2**(-2)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+(-1) &
            *y2)**2*y3*(2._8*x3+y3)-4._8*((-1._8)+nu)*lr2*(x1p-y1)**2*(( &
            x1p-y1)**2+(x2p-y2)**2)**(-1)*(x3+y3)*((x1p-y1) &
            **2+(x3+y3)**2)**(-1)+4._8*((-1._8)+nu)*lr2**(-1)*(x1p-y1)**2*(( &
            x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)**2*(x3+y3) &
            *((x1p-y1)**2+(x3+y3)**2)**(-1)+2._8*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*(x2p-y2)**2*((x1p-y1)**2+(x2p- &
            y2)**2+(x3+y3)**2)**(-3._8/2._8)*((-3._8)*nu*(x3*(x3**2+(x1p-y1) &
            **2+(x2p-y2)**2)+(x3*((-2)*lr2+3*x3)+(x1p-y1)**2+(x2p+( &
            -1)*y2)**2)*y3-(lr2-3*x3)*y3**2+y3**3)+2._8*nu**2*(x3* &
            (x3**2+(x1p-y1)**2+(x2p-y2)**2)+(x3*((-2._8)*lr2+3*x3)+( &
            x1p-y1)**2+(x2p-y2)**2)*y3-(lr2-3*x3)*y3**2+ &
            y3**3)+y3*((x1p-y1)**2+(x2p-y2)**2-(lr2-x3+( &
            -1)*y3)*(2*x3+y3)))+4._8*lr2**(-1)*((x1p-y1)**2+(x2p- &
            y2)**2)**(-2._8)*(x2p-y2)**2*((-3._8)*nu*(x3+y3)*((x1p- &
            y1)**2+(x2p-y2)**2+(x3+y3)**2)+2._8*nu**2*(x3+y3)*((x1p- &
            y1)**2+(x2p-y2)**2+(x3+y3)**2)+y3*((x1p-y1)**2+(x2p+(-1) &
            *y2)**2+(x3+y3)*(2._8*x3+y3)))+lr2**(-1)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*(6._8*nu*(x3+y3)*((x1p-y1)**2+3*(x2p+(-1) &
            *y2)**2+(x3+y3)**2)-4._8*nu**2*(x3+y3)*((x1p-y1)**2+3*( &
            x2p-y2)**2+(x3+y3)**2)-2._8*y3*((x1p-y1)**2+3._8*(x2p+( &
            -1)*y2)**2+(x3+y3)*(2._8*x3+y3)))+xLogy(3._8-4._8*nu,lr1+x3-y3) &
            +xLogy((-3._8)+6._8*nu-4._8*nu**2,lr2+x3+y3))

    END FUNCTION J2223d2

    !---------------------------------------------------------------
    !> function J2223d3
    !! computes the derivative of the J integral J2223,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2223d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2223d3=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(4._8*((-1._8)+nu)*lr1* &
            (x1p-y1)**2*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+( &
            -1)*y2)*((x1p-y1)**2+(x3-y3)**2)**(-1)-4._8*((-1._8)+ &
            nu)*lr1**(-1)*(lr1+x2p-y2)**(-1)*(x3-y3)**2-4._8*(( &
            -1._8)+nu)*lr1**(-1)*(x1p-y1)**2*((x1p-y1)**2+(x2p- &
            y2)**2)**(-1)*(x2p-y2)*((x1p-y1)**2+(x3-y3)**2) &
            **(-1)*(x3-y3)**2-((-3._8)+4._8*nu)*(x2p-y2)*(lr1+ &
            x3-y3)**(-1)-((-3._8)+4._8*nu)*lr1**(-1)*(x2p-y2)*( &
            x3-y3)*(lr1+x3-y3)**(-1)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*( &
            (x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*y3+4._8*(( &
            -1._8)+nu)*lr2**(-1)*(lr2+x2p-y2)**(-1)*(x3+y3)**2-(3._8+( &
            -6._8)*nu+4._8*nu**2)*(x2p-y2)*(lr2+x3+y3)**(-1)-(3._8-6._8* &
            nu+4*nu**2)*lr2**(-1)*(x2p-y2)*(x3+y3)*(lr2+x3+y3)**(-1)+ &
            2._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-2)*((x1p-y1)**2+(x2p- &
            y2)**2)**(-1)*(x2p-y2)*y3*(x3+y3)*(2*x3+y3)-2*lr2**( &
            -1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*( &
            3._8*nu*((-3)+2*nu)*x3**2+nu*((-3._8)+2._8*nu)*((x1p-y1)**2+( &
            x2p-y2)**2)+2._8*(2._8-9._8*nu+6._8*nu**2)*x3*y3+(3._8-9._8*nu+6._8* &
            nu**2)*y3**2)-4._8*((-1._8)+nu)*lr2*(x1p-y1)**2*((x1p- &
            y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*((x1p-y1)**2+( &
            x3+y3)**2)**(-1)+4._8*((-1._8)+nu)*lr2**(-1)*(x1p-y1)**2*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(x3+y3)**2*(( &
            x1p-y1)**2+(x3+y3)**2)**(-1)+2._8*((x1p-y1)**2+(x2p- &
            y2)**2)**(-1)*(x2p-y2)*(x3+y3)*((x1p-y1)**2+(x2p+(-1) &
            *y2)**2+(x3+y3)**2)**(-3._8/2._8)*((-3._8)*nu*(x3*(x3**2+(x1p-y1) &
            **2+(x2p-y2)**2)+(x3*((-2._8)*lr2+3*x3)+(x1p-y1)**2+(x2p+( &
            -1)*y2)**2)*y3-(lr2-3._8*x3)*y3**2+y3**3)+2._8*nu**2*(x3* &
            (x3**2+(x1p-y1)**2+(x2p-y2)**2)+(x3*((-2._8)*lr2+3*x3)+( &
            x1p-y1)**2+(x2p-y2)**2)*y3-(lr2-3*x3)*y3**2+ &
            y3**3)+y3*((x1p-y1)**2+(x2p-y2)**2-(lr2-x3+( &
            -1)*y3)*(2._8*x3+y3)))+xLogy(4._8-4._8*nu,lr1+x2p-y2)+xLogy(4._8*( &
            (-1._8)+nu),lr2+x2p-y2))

    END FUNCTION J2223d3

    !---------------------------------------------------------------
    !> function J3212d1
    !! computes the derivative of the J integral J3212,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3212d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3212d1=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*((-4._8)*((-1._8)+nu)*( &
            (-1._8)+2*nu)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p- &
            y2)**2+(lr1+x1p-y1)**(-1)*(x3-y3)+lr1**(-1)*(x1p- &
            y1)*(lr1+x1p-y1)**(-1)*(x3-y3)-4._8*((-1._8)+nu)*((-1._8)+ &
            2._8*nu)*lr2**(-1)*(x1p-y1)**2*(lr2+x3+y3)**(-1)+4._8*((-1._8)+nu) &
            *((-1._8)+2._8*nu)*lr2*((x1p-y1)**2+(x2p-y2)**2)**(-1)*( &
            x2p-y2)**2*(x3+y3)*((x2p-y2)**2+(x3+y3)**2)**(-1) &
            -4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-1)*(x1p-y1)**2*((x1p- &
            y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)**2*(x3+y3)*((x2p+( &
            -1)*y2)**2+(x3+y3)**2)**(-1)-2._8*lr2**(-1)*x3*y3*(x3+y3)*(( &
            x2p-y2)**2+(x3+y3)**2)**(-1)+2._8*x3*(x1p-y1)**2*y3*( &
            x3+y3)*((x2p-y2)**2+(x3+y3)**2)**(-1)*((x1p-y1)**2+( &
            x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)-lr2**(-1)*(x1p-y1) &
            *(lr2+x1p-y1)**(-1)*(x3+7._8*y3+8._8*nu**2*(x3+y3)-8._8*nu*( &
            x3+2._8*y3))+(lr2+x1p-y1)**(-1)*(-x3-7._8*y3-8._8* &
            nu**2*(x3+y3)+8._8*nu*(x3+2._8*y3))+xLogy((-4._8)*((-1._8)+nu)*((-1._8)+2._8* &
            nu),lr2+x3+y3))

    END FUNCTION J3212d1

    !---------------------------------------------------------------
    !> function J3212d2
    !! computes the derivative of the J integral J3212,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3212d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3212d2=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(4._8*((-1._8)+nu)*(( &
            -1._8)+2._8*nu)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**( &
            -1)*(x2p-y2)+lr1**(-1)*(lr1+x1p-y1)**(-1)*(x2p-y2) &
            *(x3-y3)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-1)*(x1p+(-1) &
            *y1)*(x2p-y2)*(lr2+x3+y3)**(-1)+4._8*lr2**(-1)*x3*(x1p- &
            y1)*(x2p-y2)*y3*(x3+y3)*((x2p-y2)**2+(x3+y3)**2)**( &
            -2)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2*(x1p-y1)*((x1p- &
            y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(x3+y3)*((x2p+(-1) &
            *y2)**2+(x3+y3)**2)**(-1)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-1) &
            *(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+( &
            -1)*y2)**3*(x3+y3)*((x2p-y2)**2+(x3+y3)**2)**(-1)+2._8*x3*( &
            x1p-y1)*(x2p-y2)*y3*(x3+y3)*((x2p-y2)**2+(x3+ &
            y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**( &
            -3._8/2._8)-lr2**(-1)*(lr2+x1p-y1)**(-1)*(x2p-y2)*(x3+ &
            7._8*y3+8._8*nu**2*(x3+y3)-8._8*nu*(x3+2._8*y3))+4._8*((-1._8)+nu)*((-1._8)+ &
            2._8*nu)*atan((x1p-y1)**(-1)*(x2p-y2))+4._8*((-1._8)+nu)*(( &
            -1._8)+2._8*nu)*atan2(lr2*(-x2p+y2),(x1p-y1)*(x3+y3)))

    END FUNCTION J3212d2

    !---------------------------------------------------------------
    !> function J3212d3
    !! computes the derivative of the J integral J3212,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3212d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3212d3=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(lr1**(-1)*(lr1+x1p+( &
            -1)*y1)**(-1)*(x3-y3)**2-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*( &
            x1p-y1)*(lr2+x3+y3)**(-1)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)* &
            lr2**(-1)*(x1p-y1)*(x3+y3)*(lr2+x3+y3)**(-1)+4._8*lr2**(-1)* &
            x3*(x1p-y1)*y3*(x3+y3)**2*((x2p-y2)**2+(x3+y3)**2) &
            **(-2)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2*(x1p-y1)*((x1p- &
            y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)**2*((x2p-y2) &
            **2+(x3+y3)**2)**(-1)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-1)*( &
            x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+(-1) &
            *y2)**2*(x3+y3)**2*((x2p-y2)**2+(x3+y3)**2)**(-1)-2._8* &
            lr2**(-1)*(x1p-y1)*y3*(2*x3+y3)*((x2p-y2)**2+(x3+y3) &
            **2)**(-1)+2._8*x3*(x1p-y1)*y3*(x3+y3)**2*((x2p-y2) &
            **2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3) &
            **2)**(-3._8/2._8)-lr2**(-1)*(lr2+x1p-y1)**(-1)*(x3+y3)*(x3+ &
            7._8*y3+8._8*nu**2*(x3+y3)-8._8*nu*(x3+2*y3))+log(lr1+x1p-y1)+ &
            xLogy((-1._8)-8._8*((-1._8)+nu)*nu,lr2+x1p-y1))

    END FUNCTION J3212d3

    !---------------------------------------------------------------
    !> function J3213d1
    !! computes the derivative of the J integral J3213,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3213d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3213d1=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*((lr1+x1p-y1)**(-1)*( &
            x2p-y2)+lr1**(-1)*(x1p-y1)*(lr1+x1p-y1)**(-1)*(x2p+ &
            (-1)*y2)-(1._8+8._8*((-1._8)+nu)*nu)*(lr2+x1p-y1)**(-1)*(x2p+ &
            (-1)*y2)-(1._8+8._8*((-1._8)+nu)*nu)*lr2**(-1)*(x1p-y1)*( &
            lr2+x1p-y1)**(-1)*(x2p-y2)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*( &
            (x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(x3+y3)+ &
            2._8*lr2**(-1)*x3*(x2p-y2)*y3*((x2p-y2)**2+(x3+y3)**2) &
            **(-1)-4._8*((-1._8)+2._8*nu)*lr2*((x1p-y1)**2+(x2p-y2)**2) &
            **(-1)*(x2p-y2)*(x3+y3)*(nu*x3+((-1._8)+nu)*y3)*((x2p+(-1) &
            *y2)**2+(x3+y3)**2)**(-1)+4._8*((-1._8)+2._8*nu)*lr2**(-1)*(x1p- &
            y1)**2*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)* &
            (x3+y3)*(nu*x3+((-1._8)+nu)*y3)*((x2p-y2)**2+(x3+y3)**2)**( &
            -1)-2._8*x3*(x1p-y1)**2*(x2p-y2)*y3*((x2p-y2) &
            **2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3) &
            **2)**(-3._8/2._8))*G**(-1)

    END FUNCTION J3213d1

    !---------------------------------------------------------------
    !> function J3213d2
    !! computes the derivative of the J integral J3213,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3213d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3213d2=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(lr1**(-1)*(lr1+x1p+( &
            -1)*y1)**(-1)*(x2p-y2)**2-(1._8+8._8*((-1._8)+nu)*nu)*lr2**( &
            -1)*(lr2+x1p-y1)**(-1)*(x2p-y2)**2-4._8*((-1._8)+nu)*(( &
            -1)+2._8*nu)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**( &
            -1)*(x3+y3)-4._8*lr2**(-1)*x3*(x1p-y1)*(x2p-y2)**2* &
            y3*((x2p-y2)**2+(x3+y3)**2)**(-2)+2._8*lr2**(-1)*x3*(x1p+(-1) &
            *y1)*y3*((x2p-y2)**2+(x3+y3)**2)**(-1)+4*((-1._8)+2._8*nu)* &
            lr2*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x3+ &
            y3)*(nu*x3+((-1._8)+nu)*y3)*((x2p-y2)**2+(x3+y3)**2)**(-1)+ &
            4._8*((-1._8)+2._8*nu)*lr2**(-1)*(x1p-y1)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*(x2p-y2)**2*(x3+y3)*(nu*x3+((-1._8)+nu)* &
            y3)*((x2p-y2)**2+(x3+y3)**2)**(-1)-2*x3*(x1p-y1)* &
            (x2p-y2)**2*y3*((x2p-y2)**2+(x3+y3)**2)**(-1)*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)+log(lr1+x1p- &
            y1)+xLogy((-1._8)-8._8*((-1._8)+nu)*nu,lr2+x1p-y1))

    END FUNCTION J3213d2

    !---------------------------------------------------------------
    !> function J3213d3
    !! computes the derivative of the J integral J3213,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3213d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3213d3=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(lr1**(-1)*(lr1+x1p+( &
            -1)*y1)**(-1)*(x2p-y2)*(x3-y3)-(1._8+8._8*((-1._8)+nu) &
            *nu)*lr2**(-1)*(lr2+x1p-y1)**(-1)*(x2p-y2)*(x3+y3)+( &
            -4._8)*lr2**(-1)*x3*(x1p-y1)*(x2p-y2)*y3*(x3+y3)*(( &
            x2p-y2)**2+(x3+y3)**2)**(-2)+2*lr2**(-1)*(x1p-y1)*(x2p+ &
            (-1)*y2)*y3*((x2p-y2)**2+(x3+y3)**2)**(-1)-4._8*((-1._8)+2._8* &
            nu)*lr2*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1) &
            *(x2p-y2)*(nu*x3+((-1._8)+nu)*y3)*((x2p-y2)**2+(x3+y3) &
            **2)**(-1)+4._8*((-1._8)+2._8*nu)*lr2**(-1)*(x1p-y1)*((x1p- &
            y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(x3+y3)**2*(nu* &
            x3+((-1._8)+nu)*y3)*((x2p-y2)**2+(x3+y3)**2)**(-1)-2._8*x3*( &
            x1p-y1)*(x2p-y2)*y3*(x3+y3)*((x2p-y2)**2+(x3+ &
            y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**( &
            -3._8/2._8)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*atan((x1p-y1)*(x2p- &
            y2)**(-1))+4._8*nu*((-1._8)+2._8*nu)*atan2(lr2*(x2p-y2),(x1p- &
            y1)*(x3+y3)))

    END FUNCTION J3213d3

    !---------------------------------------------------------------
    !> function J3223d1
    !! computes the derivative of the J integral J3223,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3223d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3223d1=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*(lr1**(-1)*(x1p-y1)+( &
            -1)*(1._8+8._8*((-1._8)+nu)*nu)*lr2**(-1)*(x1p-y1)+2*lr2**(-1)*( &
            x1p-y1)*(lr2+x3+y3)**(-1)*(3._8*x3+2._8*y3-6._8*nu*(x3+y3)+4._8* &
            nu**2*(x3+y3))+2*x3*(x1p-y1)*y3*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2+(x3+y3)**2)**(-3._8/2._8)-2._8*((-3._8)+4._8*nu)*x3*(x1p- &
            y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x3+y3)*((x1p+(-1) &
            *y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-0.5_8))*G**(-1)

    END FUNCTION J3223d1

    !---------------------------------------------------------------
    !> function J3223d2
    !! computes the derivative of the J integral J3223,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3223d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3223d2=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*(lr1**(-1)*(x2p-y2)+( &
            -1)*(1._8+8._8*((-1._8)+nu)*nu)*lr2**(-1)*(x2p-y2)+2*lr2**(-1)*( &
            x2p-y2)*(lr2+x3+y3)**(-1)*(3._8*x3+2*y3-6._8*nu*(x3+y3)+4._8* &
            nu**2*(x3+y3))+2*x3*(x2p-y2)*y3*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2+(x3+y3)**2)**(-3._8/2._8)-2._8*((-3._8)+4._8*nu)*x3*((x1p+(-1) &
            *y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(x3+y3)*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-0.5_8))*G**(-1)

    END FUNCTION J3223d2

    !---------------------------------------------------------------
    !> function J3223d3
    !! computes the derivative of the J integral J3223,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3223d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3223d3=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(lr1**(-1)*(x3+(-1) &
            *y3)+lr2**(-1)*(-x3-3._8*y3+8._8*nu*(x3+y3)-8._8*nu**2*( &
            x3+y3))+2._8*(lr2+x3+y3)**(-1)*(3._8*x3+2*y3-6._8*nu*(x3+y3)+4._8* &
            nu**2*(x3+y3))+2._8*lr2**(-1)*(x3+y3)*(lr2+x3+y3)**(-1)*(3._8*x3+2._8* &
            y3-6._8*nu*(x3+y3)+4._8*nu**2*(x3+y3))+2._8*x3*y3*(x3+y3)*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)-2._8*((-3._8)+4._8* &
            nu)*x3*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x3+y3)**2*(( &
            x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-1._8/2._8)+2._8*((-3._8)+4._8* &
            nu)*lr2**(-1)*x3*(1-(x3+y3)**2*((x1p-y1)**2+(x2p+(-1) &
            *y2)**2+(x3+y3)**2)**(-1))**(-1)+ &
            2._8*((-3._8)+4._8*nu)*ACOTH(lr2**(-1)*(x3+y3))+ &
            xLogy(6._8+4._8*nu*((-3._8)+2._8*nu),lr2+x3+y3))

    END FUNCTION J3223d3

    !---------------------------------------------------------------
    !> function J1312d1
    !! computes the derivative of the J integral J1312,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1312d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1312d1=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(4._8*((-1._8)+nu)*(( &
            -1._8)+2._8*nu)*(-x1p+y1)*((x1p-y1)**2+(x2p-y2)**2)**( &
            -1)*(x2p-y2)+lr1**(-1)*(x1p-y1)*(lr1+x2p-y2)**(-1) &
            *(x3-y3)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-1)*(x1p- &
            y1)*(x2p-y2)*(lr2+x3+y3)**(-1)+lr2**(-1)*(x1p-y1)*(lr2+ &
            x2p-y2)**(-1)*((7._8+8._8*((-2._8)+nu)*nu)*x3+y3+8._8*((-1._8)+nu)* &
            nu*y3)-4._8*lr2**(-1)*x3*(x1p-y1)*(x2p-y2)*y3*(x3+ &
            y3)*((x1p-y1)**2+(x3+y3)**2)**(-2)+4._8*((-1._8)+nu)*((-1._8)+2._8* &
            nu)*lr2*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1) &
            *(x2p-y2)*(x3+y3)*((x1p-y1)**2+(x3+y3)**2)**(-1)+4._8*( &
            (-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-1)*(x1p-y1)**3*((x1p-y1) &
            **2+(x2p-y2)**2)**(-1)*(x2p-y2)*(x3+y3)*((x1p- &
            y1)**2+(x3+y3)**2)**(-1)-2._8*x3*(x1p-y1)*(x2p-y2)* &
            y3*(x3+y3)*((x1p-y1)**2+(x3+y3)**2)**(-1)*((x1p-y1) &
            **2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)-4._8*((-1._8)+nu)*((-1._8)+ &
            2._8*nu)*atan((x1p-y1)*(x2p-y2)**(-1))+4._8*((-1._8)+nu)*(( &
            -1._8)+2._8*nu)*atan2(lr2*(x1p-y1),(x2p-y2)*(x3+y3)))

    END FUNCTION J1312d1

    !---------------------------------------------------------------
    !> function J1312d2
    !! computes the derivative of the J integral J1312,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1312d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1312d2=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(4._8*((-1._8)+nu)*(( &
            -1._8)+2._8*nu)*(x1p-y1)**2*((x1p-y1)**2+(x2p-y2)**2) &
            **(-1)+(lr1+x2p-y2)**(-1)*(x3-y3)+lr1**(-1)*(x2p- &
            y2)*(lr1+x2p-y2)**(-1)*(x3-y3)+4._8*((-1._8)+nu)*((-1._8)+2._8* &
            nu)*lr2**(-1)*(x2p-y2)**2*(lr2+x3+y3)**(-1)+(lr2+x2p-y2) &
            **(-1)*((7._8+8._8*((-2._8)+nu)*nu)*x3+y3+8._8*((-1._8)+nu)*nu*y3)+lr2**( &
            -1)*(x2p-y2)*(lr2+x2p-y2)**(-1)*((7._8+8._8*((-2._8)+nu)*nu) &
            *x3+y3+8._8*((-1._8)+nu)*nu*y3)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2*( &
            x1p-y1)**2*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x3+ &
            y3)*((x1p-y1)**2+(x3+y3)**2)**(-1)+4._8*((-1._8)+nu)*((-1._8)+2._8* &
            nu)*lr2**(-1)*(x1p-y1)**2*((x1p-y1)**2+(x2p-y2) &
            **2)**(-1)*(x2p-y2)**2*(x3+y3)*((x1p-y1)**2+(x3+y3) &
            **2)**(-1)+2*lr2**(-1)*x3*y3*(x3+y3)*((x1p-y1)**2+(x3+y3) &
            **2)**(-1)-2*x3*(x2p-y2)**2*y3*(x3+y3)*((x1p-y1) &
            **2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3) &
            **2)**(-3._8/2._8)+xLogy(4._8*((-1._8)+nu)*((-1._8)+2._8*nu),lr2+x3+y3))

    END FUNCTION J1312d2

    !---------------------------------------------------------------
    !> function J1312d3
    !! computes the derivative of the J integral J1312,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1312d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1312d3=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(lr1**(-1)*(lr1+x2p+( &
            -1)*y2)**(-1)*(x3-y3)**2+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*(x2p+( &
            -1)*y2)*(lr2+x3+y3)**(-1)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-1)*( &
            x2p-y2)*(x3+y3)*(lr2+x3+y3)**(-1)+lr2**(-1)*(lr2+x2p-y2) &
            **(-1)*(x3+y3)*((7._8+8._8*((-2._8)+nu)*nu)*x3+y3+8._8*((-1._8)+nu)*nu* &
            y3)-4._8*lr2**(-1)*x3*(x2p-y2)*y3*(x3+y3)**2*((x1p- &
            y1)**2+(x3+y3)**2)**(-2)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2*(x1p+( &
            -1)*y1)**2*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p- &
            y2)*((x1p-y1)**2+(x3+y3)**2)**(-1)+4._8*((-1._8)+nu)*((-1._8)+2._8* &
            nu)*lr2**(-1)*(x1p-y1)**2*((x1p-y1)**2+(x2p-y2) &
            **2)**(-1)*(x2p-y2)*(x3+y3)**2*((x1p-y1)**2+(x3+y3) &
            **2)**(-1)+2*lr2**(-1)*(x2p-y2)*y3*(2*x3+y3)*((x1p- &
            y1)**2+(x3+y3)**2)**(-1)-2*x3*(x2p-y2)*y3*(x3+y3)**2* &
            ((x1p-y1)**2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p- &
            y2)**2+(x3+y3)**2)**(-3._8/2._8)+log(lr1+x2p-y2)+xLogy(7._8+8._8*((-2._8)+ &
            nu)*nu,lr2+x2p-y2))

    END FUNCTION J1312d3

    !---------------------------------------------------------------
    !> function J1313d1
    !! computes the derivative of the J integral J1313,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1313d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1313d1=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*(lr1**(-1)*(x1p-y1)+2._8* &
            (7._8+8._8*((-2._8)+nu)*nu)*lr2**(-1)*(x1p-y1)-2*lr2**(-1)*(x1p+ &
            (-1)*y1)*(lr2+x3+y3)**(-1)*(3._8*x3+2*y3-6._8*nu*(x3+y3)+4._8* &
            nu**2*(x3+y3))-2._8*((-3._8)+4._8*nu)*x3*(x1p-y1)*((x1p- &
            y1)**2+(x2p-y2)**2)**(-1)*(x3+y3)*SQRT((x1p-y1)**2+(x2p+( &
            -1)*y2)**2+(x3+y3)**2)**(-1)+(x1p-y1)*((x1p-y1)**2+( &
            x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)*(((-7._8)-8._8*((-2._8)+nu)*nu)* &
            x1p**2+2._8*(7._8+8._8*((-2._8)+nu)*nu)*x1p*y1+((-7._8)-8._8*((-2._8)+nu)*nu)* &
            y1**2-7._8*(x3**2+(x2p-y2)**2)-16._8*x3*y3-7._8*y3**2+( &
            -8._8)*((-2._8)+nu)*nu*((x2p-y2)**2+(x3+y3)**2)))*G**(-1)

    END FUNCTION J1313d1

    !---------------------------------------------------------------
    !> function J1313d2
    !! computes the derivative of the J integral J1313,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1313d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1313d2=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*(lr1**(-1)*(x2p-y2)+2* &
            (7._8+8._8*((-2._8)+nu)*nu)*lr2**(-1)*(x2p-y2)-2*lr2**(-1)*(x2p+ &
            (-1)*y2)*(lr2+x3+y3)**(-1)*(3._8*x3+2*y3-6._8*nu*(x3+y3)+4._8* &
            nu**2*(x3+y3))-2._8*((-3._8)+4._8*nu)*x3*((x1p-y1)**2+(x2p+(-1) &
            *y2)**2)**(-1)*(x2p-y2)*(x3+y3)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2+(x3+y3)**2)**(-0.5_8)+(x2p-y2)*((x1p-y1)**2+( &
            x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)*(((-7._8)-8._8*((-2._8)+nu)*nu)* &
            x1p**2+2._8*(7._8+8._8*((-2._8)+nu)*nu)*x1p*y1+((-7._8)-8._8*((-2._8)+nu)*nu)* &
            y1**2-7._8*(x3**2+(x2p-y2)**2)-16._8*x3*y3-7._8*y3**2+( &
            -8._8)*((-2._8)+nu)*nu*((x2p-y2)**2+(x3+y3)**2)))*G**(-1)

    END FUNCTION J1313d2

    !---------------------------------------------------------------
    !> function J1313d3
    !! computes the derivative of the J integral J1313,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1313d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1313d3=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(lr1**(-1)*(x3+(-1) &
            *y3)+2._8*lr2**(-1)*((7._8+8._8*((-2._8)+nu)*nu)*x3+8._8*((-1._8)+nu)**2*y3)+ &
            2._8*(lr2+x3+y3)**(-1)*((-3._8)*x3-2._8*y3+6._8*nu*(x3+y3)-4._8* &
            nu**2*(x3+y3))-2._8*lr2**(-1)*(x3+y3)*(lr2+x3+y3)**(-1)*(3._8*x3+ &
            2._8*y3-6._8*nu*(x3+y3)+4._8*nu**2*(x3+y3))-2._8*((-3._8)+4._8*nu)* &
            x3*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x3+y3)**2*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-0.5_8)+(x3+y3)*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)*(((-7._8)-8._8*(( &
            -2._8)+nu)*nu)*x1p**2+2._8*(7._8+8._8*((-2._8)+nu)*nu)*x1p*y1+((-7._8)-8._8*(( &
            -2._8)+nu)*nu)*y1**2-7._8*(x3**2+(x2p-y2)**2)-16._8*x3*y3+( &
            -7._8)*y3**2-8._8*((-2._8)+nu)*nu*((x2p-y2)**2+(x3+y3)**2))+2._8* &
            ((-3._8)+4._8*nu)*lr2**(-1)*x3*(1-(x3+y3)**2*((x1p-y1) &
            **2+(x2p-y2)**2+(x3+y3)**2)**(-1))**(-1)+2._8*((-3._8)+4._8*nu)* &
            ACOTH(lr2**(-1)*(x3+y3))+ &
            xLogy((-6._8)+4._8*(3._8-2._8*nu)*nu,lr2+x3+y3))

    END FUNCTION J1313d3

    !---------------------------------------------------------------
    !> function J1323d1
    !! computes the derivative of the J integral J1323,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1323d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1323d1=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(lr1**(-1)*(x1p+(-1) &
            *y1)**2*(lr1+x2p-y2)**(-1)+(7._8+8._8*((-2._8)+nu)*nu)*lr2**(-1)*( &
            x1p-y1)**2*(lr2+x2p-y2)**(-1)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu) &
            *((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(x3+ &
            y3)+4._8*lr2**(-1)*x3*(x1p-y1)**2*(x2p-y2)*y3*((x1p+( &
            -1)*y1)**2+(x3+y3)**2)**(-2)+2._8*lr2**(-1)*x3*(-x2p+y2)*y3* &
            ((x1p-y1)**2+(x3+y3)**2)**(-1)-4._8*((-1._8)+nu)*lr2*((x1p+(-1) &
            *y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(x3+y3)*((-3)* &
            x3-y3+2._8*nu*(x3+y3))*((x1p-y1)**2+(x3+y3)**2)**(-1)+( &
            -4._8)*((-1._8)+nu)*lr2**(-1)*(x1p-y1)**2*((x1p-y1)**2+(x2p+ &
            (-1)*y2)**2)**(-1)*(x2p-y2)*(x3+y3)*((-3._8)*x3-y3+ &
            2._8*nu*(x3+y3))*((x1p-y1)**2+(x3+y3)**2)**(-1)+2._8*x3*(x1p+( &
            -1)*y1)**2*(x2p-y2)*y3*((x1p-y1)**2+(x3+y3)**2)**( &
            -1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)+log( &
            lr1+x2p-y2)+xLogy(7._8+8._8*((-2._8)+nu)*nu,lr2+x2p-y2))

    END FUNCTION J1323d1

    !---------------------------------------------------------------
    !> function J1323d2
    !! computes the derivative of the J integral J1323,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1323d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1323d2=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*((x1p-y1)*(lr1+x2p+(-1) &
            *y2)**(-1)+lr1**(-1)*(x1p-y1)*(x2p-y2)*(lr1+x2p- &
            y2)**(-1)+(7._8+8._8*((-2._8)+nu)*nu)*(x1p-y1)*(lr2+x2p-y2)**( &
            -1)+(7._8+8._8*((-2._8)+nu)*nu)*lr2**(-1)*(x1p-y1)*(x2p-y2)* &
            (lr2+x2p-y2)**(-1)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*(x1p- &
            y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x3+y3)+2*lr2**( &
            -1)*x3*(-x1p+y1)*y3*((x1p-y1)**2+(x3+y3)**2)**(-1)+ &
            4._8*((-1._8)+nu)*lr2*(x1p-y1)*((x1p-y1)**2+(x2p-y2) &
            **2)**(-1)*(x3+y3)*((-3._8)*x3-y3+2._8*nu*(x3+y3))*((x1p+(-1) &
            *y1)**2+(x3+y3)**2)**(-1)-4._8*((-1._8)+nu)*lr2**(-1)*(x1p- &
            y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)**2* &
            (x3+y3)*((-3._8)*x3-y3+2._8*nu*(x3+y3))*((x1p-y1)**2+( &
            x3+y3)**2)**(-1)+2._8*x3*(x1p-y1)*(x2p-y2)**2*y3*((x1p+ &
            (-1)*y1)**2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2) &
            **2+(x3+y3)**2)**(-3._8/2._8))*G**(-1)

    END FUNCTION J1323d2

    !---------------------------------------------------------------
    !> function J1323d3
    !! computes the derivative of the J integral J1323,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J1323d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J1323d3=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(lr1**(-1)*(x1p+(-1) &
            *y1)*(lr1+x2p-y2)**(-1)*(x3-y3)+(7._8+8._8*((-2._8)+nu)*nu) &
            *lr2**(-1)*(x1p-y1)*(lr2+x2p-y2)**(-1)*(x3+y3)+4._8* &
            lr2**(-1)*x3*(x1p-y1)*(x2p-y2)*y3*(x3+y3)*((x1p+(-1) &
            *y1)**2+(x3+y3)**2)**(-2)-2._8*lr2**(-1)*(x1p-y1)*(x2p+(-1) &
            *y2)*y3*((x1p-y1)**2+(x3+y3)**2)**(-1)+4._8*((-1._8)+nu)*lr2*( &
            x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+(-1) &
            *y2)*((-3._8)*x3-y3+2._8*nu*(x3+y3))*((x1p-y1)**2+(x3+ &
            y3)**2)**(-1)-4._8*((-1._8)+nu)*lr2**(-1)*(x1p-y1)*((x1p+(-1) &
            *y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(x3+y3)**2*((-3._8) &
            *x3-y3+2._8*nu*(x3+y3))*((x1p-y1)**2+(x3+y3)**2)**(-1) &
            +2._8*x3*(x1p-y1)*(x2p-y2)*y3*(x3+y3)*((x1p-y1) &
            **2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3) &
            **2)**(-3._8/2._8)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*atan((x1p-y1)**( &
            -1)*(x2p-y2))-4._8*((-1._8)+nu)*((-3._8)+2._8*nu)*atan2(lr2*(x1p+( &
            -1)*y1),(x2p-y2)*(x3+y3)))

    END FUNCTION J1323d3

    !---------------------------------------------------------------
    !> function J2312d1
    !! computes the derivative of the J integral J2312,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2312d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2312d1=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(4._8*((-1._8)+nu)*(( &
            -1._8)+2._8*nu)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p- &
            y2)**2+(lr1+x1p-y1)**(-1)*(x3-y3)+lr1**(-1)*(x1p- &
            y1)*(lr1+x1p-y1)**(-1)*(x3-y3)+4._8*((-1._8)+nu)*((-1._8)+2._8* &
            nu)*lr2**(-1)*(x1p-y1)**2*(lr2+x3+y3)**(-1)+(lr2+x1p-y1) &
            **(-1)*((7._8+8._8*((-2._8)+nu)*nu)*x3+y3+8._8*((-1._8)+nu)*nu*y3)+lr2**( &
            -1)*(x1p-y1)*(lr2+x1p-y1)**(-1)*((7._8+8._8*((-2._8)+nu)*nu) &
            *x3+y3+8._8*((-1._8)+nu)*nu*y3)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2*( &
            (x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)**2*(x3+ &
            y3)*((x2p-y2)**2+(x3+y3)**2)**(-1)+4._8*((-1._8)+nu)*((-1._8)+2._8* &
            nu)*lr2**(-1)*(x1p-y1)**2*((x1p-y1)**2+(x2p-y2) &
            **2)**(-1)*(x2p-y2)**2*(x3+y3)*((x2p-y2)**2+(x3+y3) &
            **2)**(-1)+2*lr2**(-1)*x3*y3*(x3+y3)*((x2p-y2)**2+(x3+y3) &
            **2)**(-1)-2*x3*(x1p-y1)**2*y3*(x3+y3)*((x2p-y2) &
            **2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3) &
            **2)**(-3._8/2._8)+xLogy(4._8*((-1._8)+nu)*((-1._8)+2._8*nu),lr2+x3+y3))

    END FUNCTION J2312d1

    !---------------------------------------------------------------
    !> function J2312d2
    !! computes the derivative of the J integral J2312,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2312d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2312d2=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*((-4._8)*((-1._8)+nu)*( &
            (-1._8)+2._8*nu)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**( &
            -1)*(x2p-y2)+lr1**(-1)*(lr1+x1p-y1)**(-1)*(x2p-y2) &
            *(x3-y3)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-1)*(x1p- &
            y1)*(x2p-y2)*(lr2+x3+y3)**(-1)+lr2**(-1)*(lr2+x1p-y1)**( &
            -1)*(x2p-y2)*((7._8+8._8*((-2._8)+nu)*nu)*x3+y3+8._8*((-1._8)+nu)* &
            nu*y3)-4._8*lr2**(-1)*x3*(x1p-y1)*(x2p-y2)*y3*(x3+ &
            y3)*((x2p-y2)**2+(x3+y3)**2)**(-2)+4._8*((-1._8)+nu)*((-1._8)+2._8* &
            nu)*lr2*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1) &
            *(x2p-y2)*(x3+y3)*((x2p-y2)**2+(x3+y3)**2)**(-1)+4._8*( &
            (-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-1)*(x1p-y1)*((x1p-y1) &
            **2+(x2p-y2)**2)**(-1)*(x2p-y2)**3*(x3+y3)*((x2p+(-1) &
            *y2)**2+(x3+y3)**2)**(-1)-2._8*x3*(x1p-y1)*(x2p-y2) &
            *y3*(x3+y3)*((x2p-y2)**2+(x3+y3)**2)**(-1)*((x1p-y1) &
            **2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)-4._8*((-1._8)+nu)*((-1._8)+ &
            2._8*nu)*atan((x1p-y1)**(-1)*(x2p-y2))+4._8*((-1._8)+nu)*(( &
            -1._8)+2._8*nu)*atan2(lr2*(x2p-y2),(x1p-y1)*(x3+y3)))

    END FUNCTION J2312d2

    !---------------------------------------------------------------
    !> function J2312d3
    !! computes the derivative of the J integral J2312,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2312d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2312d3=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(lr1**(-1)*(lr1+x1p+( &
            -1)*y1)**(-1)*(x3-y3)**2+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*(x1p+( &
            -1)*y1)*(lr2+x3+y3)**(-1)+4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2**(-1)*( &
            x1p-y1)*(x3+y3)*(lr2+x3+y3)**(-1)+lr2**(-1)*(lr2+x1p-y1) &
            **(-1)*(x3+y3)*((7._8+8._8*((-2._8)+nu)*nu)*x3+y3+8*((-1._8)+nu)*nu* &
            y3)-4._8*lr2**(-1)*x3*(x1p-y1)*y3*(x3+y3)**2*((x2p- &
            y2)**2+(x3+y3)**2)**(-2)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*lr2*(x1p+( &
            -1)*y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2) &
            **2*((x2p-y2)**2+(x3+y3)**2)**(-1)+4._8*((-1._8)+nu)*((-1._8)+2._8* &
            nu)*lr2**(-1)*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2) &
            **(-1)*(x2p-y2)**2*(x3+y3)**2*((x2p-y2)**2+(x3+y3) &
            **2)**(-1)+2._8*lr2**(-1)*(x1p-y1)*y3*(2._8*x3+y3)*((x2p- &
            y2)**2+(x3+y3)**2)**(-1)-2._8*x3*(x1p-y1)*y3*(x3+y3)**2* &
            ((x2p-y2)**2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p- &
            y2)**2+(x3+y3)**2)**(-3._8/2._8)+log(lr1+x1p-y1)+xLogy(7._8+8._8*((-2._8)+ &
            nu)*nu,lr2+x1p-y1))

    END FUNCTION J2312d3

    !---------------------------------------------------------------
    !> function J2313d1
    !! computes the derivative of the J integral J2313,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2313d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2313d1=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*((lr1+x1p-y1)**(-1)*( &
            x2p-y2)+lr1**(-1)*(x1p-y1)*(lr1+x1p-y1)**(-1)*(x2p+ &
            (-1)*y2)+(7._8+8._8*((-2._8)+nu)*nu)*(lr2+x1p-y1)**(-1)*(x2p- &
            y2)+(7._8+8._8*((-2._8)+nu)*nu)*lr2**(-1)*(x1p-y1)*(lr2+x1p- &
            y1)**(-1)*(x2p-y2)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*((x1p+(-1) &
            *y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(x3+y3)+2*lr2**( &
            -1)*x3*(-x2p+y2)*y3*((x2p-y2)**2+(x3+y3)**2)**(-1)+ &
            4._8*((-1._8)+nu)*lr2*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+( &
            -1)*y2)*(x3+y3)*((-3._8)*x3-y3+2*nu*(x3+y3))*((x2p- &
            y2)**2+(x3+y3)**2)**(-1)-4._8*((-1._8)+nu)*lr2**(-1)*(x1p-y1) &
            **2*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*( &
            x3+y3)*((-3._8)*x3-y3+2._8*nu*(x3+y3))*((x2p-y2)**2+(x3+ &
            y3)**2)**(-1)+2._8*x3*(x1p-y1)**2*(x2p-y2)*y3*((x2p+( &
            -1)*y2)**2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+ &
            (x3+y3)**2)**(-3._8/2._8))*G**(-1)

    END FUNCTION J2313d1

    !---------------------------------------------------------------
    !> function J2313d2
    !! computes the derivative of the J integral J2313,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2313d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2313d2=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(lr1**(-1)*(lr1+x1p+( &
            -1)*y1)**(-1)*(x2p-y2)**2+(7._8+8._8*((-2._8)+nu)*nu)*lr2**(-1)*( &
            lr2+x1p-y1)**(-1)*(x2p-y2)**2+4._8*((-1._8)+nu)*((-1._8)+2._8*nu) &
            *(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x3+ &
            y3)+4._8*lr2**(-1)*x3*(x1p-y1)*(x2p-y2)**2*y3*((x2p+( &
            -1)*y2)**2+(x3+y3)**2)**(-2)+2*lr2**(-1)*x3*(-x1p+y1)*y3* &
            ((x2p-y2)**2+(x3+y3)**2)**(-1)-4._8*((-1._8)+nu)*lr2*(x1p+(-1) &
            *y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x3+y3)*((-3._8)* &
            x3-y3+2._8*nu*(x3+y3))*((x2p-y2)**2+(x3+y3)**2)**(-1)+( &
            -4._8)*((-1._8)+nu)*lr2**(-1)*(x1p-y1)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*(x2p-y2)**2*(x3+y3)*((-3._8)*x3-y3+ &
            2._8*nu*(x3+y3))*((x2p-y2)**2+(x3+y3)**2)**(-1)+2._8*x3*(x1p+( &
            -1)*y1)*(x2p-y2)**2*y3*((x2p-y2)**2+(x3+y3)**2)**( &
            -1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)+log( &
            lr1+x1p-y1)+xLogy(7._8+8._8*((-2._8)+nu)*nu,lr2+x1p-y1))

    END FUNCTION J2313d2

    !---------------------------------------------------------------
    !> function J2313d3
    !! computes the derivative of the J integral J2313,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2313d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2313d3=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(lr1**(-1)*(lr1+x1p+( &
            -1)*y1)**(-1)*(x2p-y2)*(x3-y3)+(7._8+8._8*((-2._8)+nu)*nu) &
            *lr2**(-1)*(lr2+x1p-y1)**(-1)*(x2p-y2)*(x3+y3)+4* &
            lr2**(-1)*x3*(x1p-y1)*(x2p-y2)*y3*(x3+y3)*((x2p+(-1) &
            *y2)**2+(x3+y3)**2)**(-2)-2._8*lr2**(-1)*(x1p-y1)*(x2p+(-1) &
            *y2)*y3*((x2p-y2)**2+(x3+y3)**2)**(-1)+4._8*((-1._8)+nu)*lr2*( &
            x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+(-1) &
            *y2)*((-3)*x3-y3+2._8*nu*(x3+y3))*((x2p-y2)**2+(x3+ &
            y3)**2)**(-1)-4._8*((-1._8)+nu)*lr2**(-1)*(x1p-y1)*((x1p+(-1) &
            *y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*(x3+y3)**2*((-3) &
            *x3-y3+2._8*nu*(x3+y3))*((x2p-y2)**2+(x3+y3)**2)**(-1) &
            +2*x3*(x1p-y1)*(x2p-y2)*y3*(x3+y3)*((x2p-y2) &
            **2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3) &
            **2)**(-3._8/2._8)-4._8*((-1._8)+nu)*((-1._8)+2._8*nu)*atan((x1p-y1)*( &
            x2p-y2)**(-1))-4._8*((-1._8)+nu)*((-3._8)+2*nu)*atan2(lr2*(x2p+( &
            -1)*y2),(x1p-y1)*(x3+y3)))

    END FUNCTION J2313d3

    !---------------------------------------------------------------
    !> function J2323d1
    !! computes the derivative of the J integral J2323,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2323d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2323d1=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*(lr1**(-1)*(x1p-y1)+2._8* &
            (7._8+8._8*((-2._8)+nu)*nu)*lr2**(-1)*(x1p-y1)-2*lr2**(-1)*(x1p+ &
            (-1)*y1)*(lr2+x3+y3)**(-1)*(3._8*x3+2*y3-6._8*nu*(x3+y3)+4._8* &
            nu**2*(x3+y3))-2._8*((-3._8)+4._8*nu)*x3*(x1p-y1)*((x1p- &
            y1)**2+(x2p-y2)**2)**(-1)*(x3+y3)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2+(x3+y3)**2)**(-0.5_8)+(x1p-y1)*((x1p-y1)**2+( &
            x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)*(((-7._8)-8._8*((-2._8)+nu)*nu)* &
            x1p**2+2._8*(7._8+8._8*((-2._8)+nu)*nu)*x1p*y1+((-7._8)-8._8*((-2._8)+nu)*nu)* &
            y1**2-7._8*(x3**2+(x2p-y2)**2)-16._8*x3*y3-7._8*y3**2+( &
            -8._8)*((-2._8)+nu)*nu*((x2p-y2)**2+(x3+y3)**2)))*G**(-1)

    END FUNCTION J2323d1

    !---------------------------------------------------------------
    !> function J2323d2
    !! computes the derivative of the J integral J2323,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2323d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2323d2=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*(lr1**(-1)*(x2p-y2)+2._8* &
            (7._8+8._8*((-2._8)+nu)*nu)*lr2**(-1)*(x2p-y2)-2*lr2**(-1)*(x2p+ &
            (-1)*y2)*(lr2+x3+y3)**(-1)*(3._8*x3+2._8*y3-6._8*nu*(x3+y3)+4._8* &
            nu**2*(x3+y3))-2._8*((-3._8)+4._8*nu)*x3*((x1p-y1)**2+(x2p+(-1) &
            *y2)**2)**(-1)*(x2p-y2)*(x3+y3)*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2+(x3+y3)**2)**(-0.5_8)+(x2p-y2)*((x1p-y1)**2+( &
            x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)*(((-7._8)-8._8*((-2._8)+nu)*nu)* &
            x1p**2+2._8*(7._8+8._8*((-2._8)+nu)*nu)*x1p*y1+((-7._8)-8._8*((-2._8)+nu)*nu)* &
            y1**2-7._8*(x3**2+(x2p-y2)**2)-16._8*x3*y3-7*y3**2+( &
            -8._8)*((-2._8)+nu)*nu*((x2p-y2)**2+(x3+y3)**2)))*G**(-1)

    END FUNCTION J2323d2

    !---------------------------------------------------------------
    !> function J2323d3
    !! computes the derivative of the J integral J2323,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J2323d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J2323d3=(-1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(lr1**(-1)*(x3+(-1) &
            *y3)+2*lr2**(-1)*((7._8+8._8*((-2._8)+nu)*nu)*x3+8._8*((-1._8)+nu)**2*y3)+ &
            2._8*(lr2+x3+y3)**(-1)*((-3._8)*x3-2._8*y3+6._8*nu*(x3+y3)-4._8* &
            nu**2*(x3+y3))-2*lr2**(-1)*(x3+y3)*(lr2+x3+y3)**(-1)*(3._8*x3+ &
            2._8*y3-6._8*nu*(x3+y3)+4._8*nu**2*(x3+y3))-2._8*((-3._8)+4._8*nu)* &
            x3*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x3+y3)**2*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-0.5_8)+(x3+y3)*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)*(((-7._8)-8._8*(( &
            -2._8)+nu)*nu)*x1p**2+2._8*(7._8+8._8*((-2._8)+nu)*nu)*x1p*y1+((-7._8)-8._8*(( &
            -2._8)+nu)*nu)*y1**2-7._8*(x3**2+(x2p-y2)**2)-16._8*x3*y3+( &
            -7._8)*y3**2-8._8*((-2._8)+nu)*nu*((x2p-y2)**2+(x3+y3)**2))+2* &
            ((-3._8)+4._8*nu)*lr2**(-1)*x3*(1-(x3+y3)**2*((x1p-y1) &
            **2+(x2p-y2)**2+(x3+y3)**2)**(-1))**(-1)+2._8*((-3._8)+4._8*nu)* &
            ACOTH(lr2**(-1)*(x3+y3))+ &
            xLogy((-6._8)+4._8*(3._8-2._8*nu)*nu,lr2+x3+y3))

    END FUNCTION J2323d3

    !---------------------------------------------------------------
    !> function J3312d1
    !! computes the derivative of the J integral J3312,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3312d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3312d1=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(9._8*x3**2._8*(9._8* &
            x3**2+(x1p-y1)**2)**(-1)+4._8*nu**2*x3**2*(nu**2*x3**2+(x1p+( &
            -1)*y1)**2)**(-1)-((-3._8)+4._8*nu)*(lr1+x1p-y1)**(-1)*( &
            x2p-y2)-((-3._8)+4._8*nu)*lr1**(-1)*(x1p-y1)*(lr1+x1p+( &
            -1)*y1)**(-1)*(x2p-y2)+(5._8+4._8*nu*((-3._8)+2._8*nu))*(lr2+x1p+(-1) &
            *y1)**(-1)*(x2p-y2)+(5._8+4._8*nu*((-3._8)+2._8*nu))*lr2**(-1)*(x1p+ &
            (-1)*y1)*(lr2+x1p-y1)**(-1)*(x2p-y2)-((-3._8)+4._8* &
            nu)*lr1**(-1)*(x1p-y1)**2*(lr1+x2p-y2)**(-1)+(5._8+4._8*nu* &
            ((-3._8)+2._8*nu))*lr2**(-1)*(x1p-y1)**2*(lr2+x2p-y2)**(-1)+ &
            2._8*((-1._8)+2._8*nu)*lr1*(x2p-y2)*((x1p-y1)**2+(x3- &
            y3)**2)**(-1)*((x2p-y2)**2+(x3-y3)**2)**(-1)*(x3+(-1) &
            *y3)**2-2*((-1._8)+2._8*nu)*lr1**(-1)*(x1p-y1)**2*(x2p+(-1) &
            *y2)*((x1p-y1)**2+(x3-y3)**2)**(-1)*((x2p-y2) &
            **2+(x3-y3)**2)**(-1)*(x3-y3)**2+9._8*y3**2*((x1p+(-1) &
            *y1)**2+9._8*y3**2)**(-1)+4._8*nu**2*y3**2*((x1p-y1)**2+ &
            nu**2*y3**2)**(-1)-2._8*(1._8-2._8*nu)**2*lr2*(x2p-y2)*(x3+ &
            y3)**2*((x1p-y1)**2+(x3+y3)**2)**(-1)*((x2p-y2)**2+( &
            x3+y3)**2)**(-1)-4._8*lr2**(-1)*x3*(x1p-y1)**2*(x2p- &
            y2)*y3*((x1p-y1)**2+(x3+y3)**2)**(-2)*((x2p-y2)**2+( &
            x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+2*(x3+y3) &
            **2)-2._8*x3*(x1p-y1)**2*(x2p-y2)*y3*((x1p-y1) &
            **2+(x3+y3)**2)**(-1)*((x2p-y2)**2+(x3+y3)**2)**(-1)*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)*((x1p-y1) &
            **2+(x2p-y2)**2+2._8*(x3+y3)**2)+2._8*lr2**(-1)*(x2p-y2)*(( &
            x1p-y1)**2+(x3+y3)**2)**(-1)*((x2p-y2)**2+(x3+y3)**2) &
            **(-1)*(x1p**2*x3**2-2._8*x1p*x3**2*y1+x3**2*y1**2+5._8*x1p**2* &
            x3*y3+x2p**2*x3*y3+2._8*x3**3*y3-10._8*x1p*x3*y1*y3+5._8*x3* &
            y1**2*y3-2._8*x2p*x3*y2*y3+x3*y2**2*y3+x1p**2*y3**2+4._8* &
            x3**2*y3**2-2._8*x1p*y1*y3**2+y1**2*y3**2+2._8*x3*y3**3-4._8* &
            nu*(x1p-y1)**2*(x3+y3)**2+4._8*nu**2*(x1p-y1)**2*(x3+ &
            y3)**2)+xLogy(3._8-4._8*nu,lr1+x2p-y2)+xLogy(5._8+4._8*nu*((-3._8)+2._8* &
            nu),lr2+x2p-y2))

    END FUNCTION J3312d1

    !---------------------------------------------------------------
    !> function J3312d2
    !! computes the derivative of the J integral J3312,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3312d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3312d2=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(25._8*x3**2*(25._8* &
            x3**2+(x2p-y2)**2)**(-1)+36._8*nu**2*x3**2*(9._8*nu**2*x3**2+( &
            x2p-y2)**2)**(-1)+8._8*nu**4*x3**2._8*(nu**4*x3**2+(x2p- &
            y2)**2)**(-1)-((-3._8)+4._8*nu)*lr1**(-1)*(lr1+x1p-y1)**(-1) &
            *(x2p-y2)**2+(5._8+4._8*nu*((-3._8)+2._8*nu))*lr2**(-1)*(lr2+x1p+(-1) &
            *y1)**(-1)*(x2p-y2)**2-((-3._8)+4._8*nu)*(x1p-y1)*( &
            lr1+x2p-y2)**(-1)-((-3._8)+4._8*nu)*lr1**(-1)*(x1p-y1) &
            *(x2p-y2)*(lr1+x2p-y2)**(-1)+(5._8+4._8*nu*((-3._8)+2._8*nu))*( &
            x1p-y1)*(lr2+x2p-y2)**(-1)+(5._8+4._8*nu*((-3._8)+2._8*nu))* &
            lr2**(-1)*(x1p-y1)*(x2p-y2)*(lr2+x2p-y2)**(-1)+2._8* &
            ((-1._8)+2._8*nu)*lr1*(x1p-y1)*((x1p-y1)**2+(x3-y3) &
            **2)**(-1)*((x2p-y2)**2+(x3-y3)**2)**(-1)*(x3- &
            y3)**2-2._8*((-1._8)+2._8*nu)*lr1**(-1)*(x1p-y1)*(x2p-y2) &
            **2*((x1p-y1)**2+(x3-y3)**2)**(-1)*((x2p-y2)**2+ &
            (x3-y3)**2)**(-1)*(x3-y3)**2+25._8*y3**2*((x2p- &
            y2)**2+25._8*y3**2)**(-1)+36._8*nu**2*y3**2*((x2p-y2)**2+9._8* &
            nu**2*y3**2)**(-1)+8._8*nu**4*y3**2*((x2p-y2)**2+nu**4* &
            y3**2)**(-1)-2._8*(1._8-2._8*nu)**2*lr2*(x1p-y1)*(x3+y3) &
            **2*((x1p-y1)**2+(x3+y3)**2)**(-1)*((x2p-y2)**2+(x3+ &
            y3)**2)**(-1)-4._8*lr2**(-1)*x3*(x1p-y1)*(x2p-y2) &
            **2*y3*((x1p-y1)**2+(x3+y3)**2)**(-1)*((x2p-y2)**2+( &
            x3+y3)**2)**(-2)*((x1p-y1)**2+(x2p-y2)**2+2._8*(x3+y3) &
            **2)-2._8*x3*(x1p-y1)*(x2p-y2)**2*y3*((x1p-y1) &
            **2+(x3+y3)**2)**(-1)*((x2p-y2)**2+(x3+y3)**2)**(-1)*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)*((x1p-y1) &
            **2+(x2p-y2)**2+2*(x3+y3)**2)+2*lr2**(-1)*(x1p-y1)*(( &
            x1p-y1)**2+(x3+y3)**2)**(-1)*((x2p-y2)**2+(x3+y3)**2) &
            **(-1)*(x2p**2*x3**2-2._8*x2p*x3**2*y2+x3**2*y2**2+x1p**2*x3* &
            y3+5._8*x2p**2*x3*y3+2._8*x3**3*y3-2._8*x1p*x3*y1*y3+x3*y1**2* &
            y3-10._8*x2p*x3*y2*y3+5._8*x3*y2**2*y3+x2p**2*y3**2+4._8*x3**2* &
            y3**2-2._8*x2p*y2*y3**2+y2**2*y3**2+2._8*x3*y3**3-4._8*nu*(x2p+ &
            (-1)*y2)**2*(x3+y3)**2+4._8*nu**2*(x2p-y2)**2*(x3+y3)**2)+ &
            xLogy(3._8-4._8*nu,lr1+x1p-y1)+xLogy(5._8+4._8*nu*((-3._8)+2._8*nu),lr2+ &
            x1p-y1))

    END FUNCTION J3312d2

    !---------------------------------------------------------------
    !> function J3312d3
    !! computes the derivative of the J integral J3312,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3312d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3312d3=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(9._8*x3*(9._8*x3**2+( &
            x1p-y1)**2)**(-1)*(-x1p+y1)+4._8*nu**2*x3*(nu**2*x3**2+ &
            (x1p-y1)**2)**(-1)*(-x1p+y1)+25._8*x3*(25._8*x3**2+(x2p+(-1) &
            *y2)**2)**(-1)*(-x2p+y2)+36._8*nu**2*x3*(9._8*nu**2*x3**2+( &
            x2p-y2)**2)**(-1)*(-x2p+y2)+8._8*nu**4*x3*(nu**4*x3**2+ &
            (x2p-y2)**2)**(-1)*(-x2p+y2)-((-3._8)+4._8*nu)*lr1**( &
            -1)*(lr1+x1p-y1)**(-1)*(x2p-y2)*(x3-y3)-(( &
            -3)+4._8*nu)*lr1**(-1)*(x1p-y1)*(lr1+x2p-y2)**(-1)*(x3+( &
            -1)*y3)-2*((-1._8)+2._8*nu)*lr1*(x1p-y1)*(x2p-y2)*(( &
            x1p-y1)**2+(x3-y3)**2)**(-1)*((x2p-y2)**2+(x3+( &
            -1)*y3)**2)**(-1)*(x3-y3)-2._8*((-1._8)+2._8*nu)*lr1**(-1)*( &
            x1p-y1)*(x2p-y2)*((x1p-y1)**2+(x3-y3)**2) &
            **(-1)*((x2p-y2)**2+(x3-y3)**2)**(-1)*(x3-y3) &
            **3+(5._8+4._8*nu*((-3._8)+2._8*nu))*lr2**(-1)*(lr2+x1p-y1)**(-1)*( &
            x2p-y2)*(x3+y3)+(5._8+4._8*nu*((-3._8)+2._8*nu))*lr2**(-1)*(x1p+(-1) &
            *y1)*(lr2+x2p-y2)**(-1)*(x3+y3)+2._8*(1._8-2._8*nu)**2*lr2*( &
            x1p-y1)*(x2p-y2)*(x3+y3)*((x1p-y1)**2+(x3+y3) &
            **2)**(-1)*((x2p-y2)**2+(x3+y3)**2)**(-1)-4._8*lr2**(-1)* &
            x3*(x1p-y1)*(x2p-y2)*y3*(x3+y3)*((x1p-y1)**2+( &
            x3+y3)**2)**(-1)*((x2p-y2)**2+(x3+y3)**2)**(-2)*((x1p- &
            y1)**2+(x2p-y2)**2+2*(x3+y3)**2)-4._8*lr2**(-1)*x3*(x1p+( &
            -1)*y1)*(x2p-y2)*y3*(x3+y3)*((x1p-y1)**2+(x3+y3) &
            **2)**(-2)*((x2p-y2)**2+(x3+y3)**2)**(-1)*((x1p-y1) &
            **2+(x2p-y2)**2+2._8*(x3+y3)**2)-2._8*x3*(x1p-y1)*(x2p+( &
            -1)*y2)*y3*(x3+y3)*((x1p-y1)**2+(x3+y3)**2)**(-1)*((x2p+( &
            -1)*y2)**2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+ &
            (x3+y3)**2)**(-3._8/2._8)*((x1p-y1)**2+(x2p-y2)**2+2*(x3+y3) &
            **2)+2._8*lr2**(-1)*(x1p-y1)*(x2p-y2)*((x1p-y1)**2+ &
            (x3+y3)**2)**(-1)*((x2p-y2)**2+(x3+y3)**2)**(-1)*(x3**3+( &
            x1p**2+x2p**2)*y3+9._8*x3**2*y3-2._8*x1p*y1*y3+y1**2*y3-2._8* &
            x2p*y2*y3+y2**2*y3+11._8*x3*y3**2+3._8*y3**3-4._8*nu*(x3+y3)**3+ &
            4._8*nu**2*(x3+y3)**3)-3._8*atan2(3._8*x3,x1p-y1)-5._8*atan2( &
            5._8*x3,x2p-y2)+12._8*nu*atan2((-3._8)*nu*x3,x2p-y2)+4._8*nu* &
            atan2(-nu*x3,x1p-y1)-8._8*nu**2*atan2(nu**2*x3,x2p+( &
            -1)*y2)+((-2._8)+4._8*nu)*atan2(lr1*(-x3+y3),(x1p-y1)*(x2p+ &
            (-1)*y2))+2._8*(1._8-2._8*nu)**2*atan2(lr2*(x3+y3),(x1p-y1)*( &
            x2p-y2)))

    END FUNCTION J3312d3

    !---------------------------------------------------------------
    !> function J3313d1
    !! computes the derivative of the J integral J3313,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3313d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3313d1=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(9._8*x2p**2*(9._8* &
            x2p**2+(x1p-y1)**2)**(-1)+25._8*x2p**2*(25._8*x2p**2+(x1p-y1) &
            **2)**(-1)+4._8*nu**2*x2p**2*(nu**2*x2p**2+(x1p-y1)**2)**(-1)+ &
            36._8*nu**2*x2p**2*(9._8*nu**2*x2p**2+(x1p-y1)**2)**(-1)+8._8* &
            nu**4*x2p**2*(nu**4*x2p**2+(x1p-y1)**2)**(-1)+4._8*x3*(lr2+( &
            -1)*x1p+y1)**(-1)+4._8*lr2**(-1)*x3*(-x1p+y1)*(lr2-x1p+y1) &
            **(-1)+9._8*y2**2*((x1p-y1)**2+9._8*y2**2)**(-1)+25._8*y2**2*(( &
            x1p-y1)**2+25._8*y2**2)**(-1)+4._8*nu**2*y2**2*((x1p-y1) &
            **2+nu**2*y2**2)**(-1)+36._8*nu**2*y2**2*((x1p-y1)**2+9._8* &
            nu**2*y2**2)**(-1)+8._8*nu**4*y2**2*((x1p-y1)**2+nu**4* &
            y2**2)**(-1)-((-3._8)+4._8*nu)*(lr1+x1p-y1)**(-1)*(x3+(-1) &
            *y3)-((-3._8)+4._8*nu)*lr1**(-1)*(x1p-y1)*(lr1+x1p- &
            y1)**(-1)*(x3-y3)+4._8*((-1._8)+nu)*lr1*((x1p-y1)**2+(x2p+( &
            -1)*y2)**2)**(-1)*(x2p-y2)**2*((x2p-y2)**2+(x3- &
            y3)**2)**(-1)*(x3-y3)-4._8*((-1._8)+nu)*lr1**(-1)*(x1p- &
            y1)**2*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2) &
            **2*((x2p-y2)**2+(x3-y3)**2)**(-1)*(x3-y3) &
            -4._8*((-1._8)+nu)*lr1**(-1)*(x1p-y1)**2*(lr1+x3-y3)**(-1)+( &
            -8._8)*((-1._8)+nu)**2*lr2**(-1)*(x1p-y1)**2*(lr2+x3+y3)**(-1)+( &
            lr2+x1p-y1)**(-1)*((-7._8)*x3-5._8*y3+12._8*nu*(x3+y3)-8._8* &
            nu**2*(x3+y3))-lr2**(-1)*(x1p-y1)*(lr2+x1p-y1)**( &
            -1)*(7._8*x3+5._8*y3-12._8*nu*(x3+y3)+8._8*nu**2*(x3+y3))+8._8*((-1._8)+ &
            nu)**2*lr2*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p- &
            y2)**2*(x3+y3)*((x2p-y2)**2+(x3+y3)**2)**(-1)-8._8*((-1._8)+ &
            nu)**2*lr2**(-1)*(x1p-y1)**2*((x1p-y1)**2+(x2p- &
            y2)**2)**(-1)*(x2p-y2)**2*(x3+y3)*((x2p-y2)**2+(x3+ &
            y3)**2)**(-1)+2*lr2**(-1)*x3*y3*(x3+y3)*((x2p-y2)**2+(x3+ &
            y3)**2)**(-1)-2*x3*(x1p-y1)**2*y3*(x3+y3)*((x2p- &
            y2)**2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+ &
            y3)**2)**(-3._8/2._8)+xLogy(4._8-4._8*nu,lr1+x3-y3)+xLogy((-8._8)*((-1._8) &
            +nu)**2,lr2+x3+y3))

    END FUNCTION J3313d1

    !---------------------------------------------------------------
    !> function J3313d2
    !! computes the derivative of the J integral J3313,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3313d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3313d2=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(9._8*x2p*(9._8*x2p**2+( &
            x1p-y1)**2)**(-1)*(-x1p+y1)+25._8*x2p*(25*x2p**2+(x1p+(-1) &
            *y1)**2)**(-1)*(-x1p+y1)+4._8*nu**2*x2p*(nu**2*x2p**2+(x1p+( &
            -1)*y1)**2)**(-1)*(-x1p+y1)+36._8*nu**2*x2p*(9._8*nu**2*x2p**2+ &
            (x1p-y1)**2)**(-1)*(-x1p+y1)+8._8*nu**4*x2p*(nu**4* &
            x2p**2+(x1p-y1)**2)**(-1)*(-x1p+y1)+4._8*lr2**(-1)*x3*(lr2+ &
            (-1)*x1p+y1)**(-1)*(-x2p+y2)-((-3._8)+4._8*nu)*lr1**(-1)*( &
            lr1+x1p-y1)**(-1)*(x2p-y2)*(x3-y3)-4._8*((-1._8)+ &
            nu)*lr1*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1) &
            *(x2p-y2)*((x2p-y2)**2+(x3-y3)**2)**(-1)*(x3+( &
            -1)*y3)-4._8*((-1._8)+nu)*lr1**(-1)*(x1p-y1)*((x1p-y1) &
            **2+(x2p-y2)**2)**(-1)*(x2p-y2)**3*((x2p-y2)**2+( &
            x3-y3)**2)**(-1)*(x3-y3)-4._8*((-1._8)+nu)*lr1**(-1)*( &
            x1p-y1)*(x2p-y2)*(lr1+x3-y3)**(-1)-8._8*((-1._8)+ &
            nu)**2*lr2**(-1)*(x1p-y1)*(x2p-y2)*(lr2+x3+y3)**(-1)+( &
            -1)*lr2**(-1)*(lr2+x1p-y1)**(-1)*(x2p-y2)*(7._8*x3+5* &
            y3-12._8*nu*(x3+y3)+8._8*nu**2*(x3+y3))-4._8*lr2**(-1)*x3*(x1p+( &
            -1)*y1)*(x2p-y2)*y3*(x3+y3)*((x2p-y2)**2+(x3+y3) &
            **2)**(-2)-8._8*((-1._8)+nu)**2*lr2*(x1p-y1)*((x1p-y1) &
            **2+(x2p-y2)**2)**(-1)*(x2p-y2)*(x3+y3)*((x2p- &
            y2)**2+(x3+y3)**2)**(-1)-8._8*((-1._8)+nu)**2*lr2**(-1)*(x1p- &
            y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)**3* &
            (x3+y3)*((x2p-y2)**2+(x3+y3)**2)**(-1)-2._8*x3*(x1p- &
            y1)*(x2p-y2)*y3*(x3+y3)*((x2p-y2)**2+(x3+y3)**2)**( &
            -1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)+5._8* &
            atan2((-5._8)*x2p,x1p-y1)-3._8*atan2(3._8*x2p,x1p-y1)+4._8*nu* &
            atan2(-nu*x2p,x1p-y1)-12._8*nu*atan2(3._8*nu*x2p,x1p+(-1) &
            *y1)+8._8*nu**2*atan2(-nu**2*x2p,x1p-y1)+(4._8-4._8*nu)* &
            atan2(lr1*(x2p-y2),(x1p-y1)*(x3-y3))-8._8*((-1._8)+ &
            nu)**2*atan2(lr2*(x2p-y2),(x1p-y1)*(x3+y3)))

    END FUNCTION J3313d2

    !---------------------------------------------------------------
    !> function J3313d3
    !! computes the derivative of the J integral J3313,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3313d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3313d3=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(4._8*((-1._8)+nu)*lr1* &
            (x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+(-1) &
            *y2)**2*((x2p-y2)**2+(x3-y3)**2)**(-1)-((-3._8)+ &
            4._8*nu)*lr1**(-1)*(lr1+x1p-y1)**(-1)*(x3-y3)**2-4._8*( &
            (-1._8)+nu)*lr1**(-1)*(x1p-y1)*((x1p-y1)**2+(x2p-y2) &
            **2)**(-1)*(x2p-y2)**2*((x2p-y2)**2+(x3-y3)**2) &
            **(-1)*(x3-y3)**2-4._8*((-1._8)+nu)*(x1p-y1)*(lr1+x3+( &
            -1)*y3)**(-1)-4._8*((-1._8)+nu)*lr1**(-1)*(x1p-y1)*(x3+(-1) &
            *y3)*(lr1+x3-y3)**(-1)-4._8*lr2**(-1)*x3*(lr2-x1p+y1) &
            **(-1)*(x3+y3)-8._8*((-1._8)+nu)**2*(x1p-y1)*(lr2+x3+y3)**( &
            -1)-8._8*((-1._8)+nu)**2*lr2**(-1)*(x1p-y1)*(x3+y3)*(lr2+x3+ &
            y3)**(-1)-lr2**(-1)*(lr2+x1p-y1)**(-1)*(x3+y3)*(7*x3+ &
            5._8*y3-12._8*nu*(x3+y3)+8._8*nu**2*(x3+y3))-4._8*lr2**(-1)*x3*( &
            x1p-y1)*y3*(x3+y3)**2*((x2p-y2)**2+(x3+y3)**2)**(-2)+ &
            8._8*((-1._8)+nu)**2*lr2*(x1p-y1)*((x1p-y1)**2+(x2p- &
            y2)**2)**(-1)*(x2p-y2)**2*((x2p-y2)**2+(x3+y3)**2)**( &
            -1)-8._8*((-1._8)+nu)**2*lr2**(-1)*(x1p-y1)*((x1p-y1) &
            **2+(x2p-y2)**2)**(-1)*(x2p-y2)**2*(x3+y3)**2*((x2p+( &
            -1)*y2)**2+(x3+y3)**2)**(-1)+2._8*lr2**(-1)*(x1p-y1)*y3*(2._8* &
            x3+y3)*((x2p-y2)**2+(x3+y3)**2)**(-1)-2._8*x3*(x1p- &
            y1)*y3*(x3+y3)**2*((x2p-y2)**2+(x3+y3)**2)**(-1)*((x1p+( &
            -1)*y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)+xLogy((-4._8),lr2+( &
            -1)*x1p+y1)+xLogy(3._8-4._8*nu,lr1+x1p-y1)+xLogy((-7._8)+4._8*(3._8 &
            -2._8*nu)*nu,lr2+x1p-y1))

    END FUNCTION J3313d3

    !---------------------------------------------------------------
    !> function J3323d1
    !! computes the derivative of the J integral J3323,1
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3323d1(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3323d1=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(9._8*x1p*(9._8*x1p**2+( &
            x2p-y2)**2)**(-1)*(-x2p+y2)+25._8*x1p*(25._8*x1p**2+(x2p+(-1) &
            *y2)**2)**(-1)*(-x2p+y2)+4._8*nu**2*x1p*(nu**2*x1p**2+(x2p+( &
            -1)*y2)**2)**(-1)*(-x2p+y2)+36._8*nu**2*x1p*(9._8*nu**2*x1p**2+ &
            (x2p-y2)**2)**(-1)*(-x2p+y2)+8._8*nu**4*x1p*(nu**4* &
            x1p**2+(x2p-y2)**2)**(-1)*(-x2p+y2)+4._8*lr2**(-1)*x3*(( &
            -1)*x1p+y1)*(lr2-x2p+y2)**(-1)-((-3._8)+4._8*nu)*lr1**(-1)* &
            (x1p-y1)*(lr1+x2p-y2)**(-1)*(x3-y3)-4._8*((-1._8)+ &
            nu)*lr1*(x1p-y1)*((x1p-y1)**2+(x2p-y2)**2)**(-1) &
            *(x2p-y2)*((x1p-y1)**2+(x3-y3)**2)**(-1)*(x3+( &
            -1)*y3)-4._8*((-1._8)+nu)*lr1**(-1)*(x1p-y1)**3*((x1p- &
            y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)*((x1p-y1)**2+( &
            x3-y3)**2)**(-1)*(x3-y3)-4._8*((-1._8)+nu)*lr1**(-1)*( &
            x1p-y1)*(x2p-y2)*(lr1+x3-y3)**(-1)-8._8*((-1._8)+ &
            nu)**2*lr2**(-1)*(x1p-y1)*(x2p-y2)*(lr2+x3+y3)**(-1)+( &
            -1)*lr2**(-1)*(x1p-y1)*(lr2+x2p-y2)**(-1)*(7._8*x3+5._8* &
            y3-12._8*nu*(x3+y3)+8._8*nu**2*(x3+y3))-4._8*lr2**(-1)*x3*(x1p+( &
            -1)*y1)*(x2p-y2)*y3*(x3+y3)*((x1p-y1)**2+(x3+y3) &
            **2)**(-2)-8._8*((-1._8)+nu)**2*lr2*(x1p-y1)*((x1p-y1) &
            **2+(x2p-y2)**2)**(-1)*(x2p-y2)*(x3+y3)*((x1p- &
            y1)**2+(x3+y3)**2)**(-1)-8._8*((-1._8)+nu)**2*lr2**(-1)*(x1p- &
            y1)**3*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)* &
            (x3+y3)*((x1p-y1)**2+(x3+y3)**2)**(-1)-2._8*x3*(x1p- &
            y1)*(x2p-y2)*y3*(x3+y3)*((x1p-y1)**2+(x3+y3)**2)**( &
            -1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)+5._8* &
            atan2((-5._8)*x1p,x2p-y2)-3._8*atan2(3._8*x1p,x2p-y2)+4._8*nu* &
            atan2(-nu*x1p,x2p-y2)-12._8*nu*atan2(3._8*nu*x1p,x2p+(-1) &
            *y2)+8._8*nu**2*atan2(-nu**2*x1p,x2p-y2)+(4._8-4._8*nu)* &
            atan2(lr1*(x1p-y1),(x2p-y2)*(x3-y3))-8._8*((-1._8)+ &
            nu)**2*atan2(lr2*(x1p-y1),(x2p-y2)*(x3+y3)))

    END FUNCTION J3323d1

    !---------------------------------------------------------------
    !> function J3323d2
    !! computes the derivative of the J integral J3323,2
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3323d2(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3323d2=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(9._8*x1p**2*(9._8* &
            x1p**2+(x2p-y2)**2)**(-1)+25._8*x1p**2*(25._8*x1p**2+(x2p-y2) &
            **2)**(-1)+4._8*nu**2*x1p**2*(nu**2*x1p**2+(x2p-y2)**2)**(-1)+ &
            36._8*nu**2*x1p**2*(9._8*nu**2*x1p**2+(x2p-y2)**2)**(-1)+8._8* &
            nu**4*x1p**2*(nu**4*x1p**2+(x2p-y2)**2)**(-1)+9._8*y1**2*(9._8* &
            y1**2+(x2p-y2)**2)**(-1)+25._8*y1**2*(25._8*y1**2+(x2p-y2) &
            **2)**(-1)+4._8*nu**2*y1**2*(nu**2*y1**2+(x2p-y2)**2)**(-1)+ &
            36._8*nu**2*y1**2*(9._8*nu**2*y1**2+(x2p-y2)**2)**(-1)+8._8* &
            nu**4*y1**2*(nu**4*y1**2+(x2p-y2)**2)**(-1)+4._8*x3*(lr2+( &
            -1)*x2p+y2)**(-1)+4._8*lr2**(-1)*x3*(-x2p+y2)*(lr2-x2p+y2) &
            **(-1)-((-3._8)+4._8*nu)*(lr1+x2p-y2)**(-1)*(x3-y3)+( &
            -1)*((-3._8)+4._8*nu)*lr1**(-1)*(x2p-y2)*(lr1+x2p-y2)**(-1) &
            *(x3-y3)+4._8*((-1._8)+nu)*lr1*(x1p-y1)**2*((x1p-y1) &
            **2+(x2p-y2)**2)**(-1)*((x1p-y1)**2+(x3-y3)**2) &
            **(-1)*(x3-y3)-4._8*((-1._8)+nu)*lr1**(-1)*(x1p-y1) &
            **2*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p-y2)**2* &
            ((x1p-y1)**2+(x3-y3)**2)**(-1)*(x3-y3)-4._8*(( &
            -1._8)+nu)*lr1**(-1)*(x2p-y2)**2*(lr1+x3-y3)**(-1)-8._8* &
            ((-1._8)+nu)**2*lr2**(-1)*(x2p-y2)**2*(lr2+x3+y3)**(-1)+(lr2+x2p+ &
            (-1)*y2)**(-1)*((-7._8)*x3-5._8*y3+12._8*nu*(x3+y3)-8._8*nu**2*( &
            x3+y3))-lr2**(-1)*(x2p-y2)*(lr2+x2p-y2)**(-1)*( &
            7._8*x3+5._8*y3-12._8*nu*(x3+y3)+8._8*nu**2*(x3+y3))+8._8*((-1._8)+nu) &
            **2*lr2*(x1p-y1)**2*((x1p-y1)**2+(x2p-y2)**2)**( &
            -1)*(x3+y3)*((x1p-y1)**2+(x3+y3)**2)**(-1)-8._8*((-1._8)+nu) &
            **2*lr2**(-1)*(x1p-y1)**2*((x1p-y1)**2+(x2p-y2) &
            **2)**(-1)*(x2p-y2)**2*(x3+y3)*((x1p-y1)**2+(x3+y3) &
            **2)**(-1)+2*lr2**(-1)*x3*y3*(x3+y3)*((x1p-y1)**2+(x3+y3) &
            **2)**(-1)-2*x3*(x2p-y2)**2*y3*(x3+y3)*((x1p-y1) &
            **2+(x3+y3)**2)**(-1)*((x1p-y1)**2+(x2p-y2)**2+(x3+y3) &
            **2)**(-3._8/2._8)+xLogy(4._8-4._8*nu,lr1+x3-y3)+xLogy((-8._8)*((-1._8)+ &
            nu)**2,lr2+x3+y3))

    END FUNCTION J3323d2

    !---------------------------------------------------------------
    !> function J3323d3
    !! computes the derivative of the J integral J3323,3
    !---------------------------------------------------------------
    REAL*8 FUNCTION J3323d3(y1,y2,y3)
      REAL*8, INTENT(IN) :: y1,y2,y3

      REAL*8 :: lr1,lr2

      lr1=r1(y1,y2,y3)
      lr2=r2(y1,y2,y3)

      J3323d3=(1._8/16._8)*(1._8-nu)**(-1)*PI**(-1)*G**(-1)*(4._8*((-1._8)+nu)*lr1* &
            (x1p-y1)**2*((x1p-y1)**2+(x2p-y2)**2)**(-1)*(x2p+( &
            -1)*y2)*((x1p-y1)**2+(x3-y3)**2)**(-1)-((-3._8)+ &
            4._8*nu)*lr1**(-1)*(lr1+x2p-y2)**(-1)*(x3-y3)**2-4._8*( &
            (-1._8)+nu)*lr1**(-1)*(x1p-y1)**2*((x1p-y1)**2+(x2p- &
            y2)**2)**(-1)*(x2p-y2)*((x1p-y1)**2+(x3-y3)**2) &
            **(-1)*(x3-y3)**2-4._8*((-1._8)+nu)*(x2p-y2)*(lr1+x3+( &
            -1)*y3)**(-1)-4._8*((-1._8)+nu)*lr1**(-1)*(x2p-y2)*(x3+(-1) &
            *y3)*(lr1+x3-y3)**(-1)-4._8*lr2**(-1)*x3*(lr2-x2p+y2) &
            **(-1)*(x3+y3)-8._8*((-1._8)+nu)**2*(x2p-y2)*(lr2+x3+y3)**( &
            -1)-8._8*((-1._8)+nu)**2*lr2**(-1)*(x2p-y2)*(x3+y3)*(lr2+x3+ &
            y3)**(-1)-lr2**(-1)*(lr2+x2p-y2)**(-1)*(x3+y3)*(7._8*x3+ &
            5._8*y3-12._8*nu*(x3+y3)+8._8*nu**2*(x3+y3))-4._8*lr2**(-1)*x3*( &
            x2p-y2)*y3*(x3+y3)**2*((x1p-y1)**2+(x3+y3)**2)**(-2)+ &
            8._8*((-1._8)+nu)**2*lr2*(x1p-y1)**2*((x1p-y1)**2+(x2p+(-1) &
            *y2)**2)**(-1)*(x2p-y2)*((x1p-y1)**2+(x3+y3)**2)**( &
            -1)-8._8*((-1._8)+nu)**2*lr2**(-1)*(x1p-y1)**2*((x1p-y1) &
            **2+(x2p-y2)**2)**(-1)*(x2p-y2)*(x3+y3)**2*((x1p+(-1) &
            *y1)**2+(x3+y3)**2)**(-1)+2*lr2**(-1)*(x2p-y2)*y3*(2._8*x3+ &
            y3)*((x1p-y1)**2+(x3+y3)**2)**(-1)-2._8*x3*(x2p-y2)* &
            y3*(x3+y3)**2*((x1p-y1)**2+(x3+y3)**2)**(-1)*((x1p- &
            y1)**2+(x2p-y2)**2+(x3+y3)**2)**(-3._8/2._8)+xLogy((-4._8),lr2- &
            x2p+y2)+xLogy(3._8-4._8*nu,lr1+x2p-y2)+xLogy((-7._8)+4._8*(3._8-2._8*nu) &
            *nu,lr2+x2p-y2))

    END FUNCTION J3323d3

  END SUBROUTINE computeStressVerticalStrainVolume

  !------------------------------------------------------------------------
  !> subroutine computeTractionKernelsVerticalStrainVolume
  !! calculates the traction kernels associated with distributed anelastic
  !! strain in an elastic half-space using the analytic solution of
  !!
  !!   Barbot S., J. D. P. Moore and V. Lambert, Displacement and Stress
  !!   Associated with Distributed Anelastic Deformation in a Half Space,
  !!   Bull. Seism. Soc. Am., 107(2), 10.1785/0120160237, 2017.
  !!
  !! INPUT:
  !! @param x1,x2,x3   - coordinates of observation points
  !! @param sv,dv,nv   - strike, dip and normal vector of observation points
  !! @param ns         - number of sources
  !! @param y          - position (NED) of vertical strain volume
  !! @param L,T,W      - length, thickness, and width of the strain volume
  !! @param strike     - strike of the vertical strain volume
  !! @param e11p,e12p
  !!        e13p,e22p
  !!        e23p,e33p  - strain components in the primed reference system
  !!                     tied to the strain volume
  !! @param G,lambda   - rigidity and Lame parameter
  !!
  !! OUTPUT:
  !! ts,td,tn          - traction in the strike, dip and normal directions
  !!
  !! \author Sylvain Barbot (21/02/17) - original fortran form
  !------------------------------------------------------------------------
  SUBROUTINE computeTractionKernelsVerticalStrainVolume( &
                         x,sv,dv,nv, &
                         ns,y,L,T,W,strike, &
                         e11p,e12p,e13p,e22p,e23p,e33p,G,lambda, &
                         ts,td,tn)
    REAL*8, DIMENSION(3), INTENT(IN) :: x
    REAL*8, DIMENSION(3), INTENT(IN) :: sv,dv,nv
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: y
    REAL*8, DIMENSION(ns), INTENT(IN) :: L,T,W,strike
    REAL*8, INTENT(IN) :: e11p,e12p,e13p,e22p,e23p,e33p
    REAL*8, INTENT(IN) :: G,lambda
    REAL*8, DIMENSION(ns), INTENT(OUT) :: ts,td,tn

    INTEGER :: i
    REAL*8 :: s11,s12,s13,s22,s23,s33
    REAL*8 :: nu

    nu=lambda/(lambda+G)/2._8

    DO i=1,ns
       CALL computeStressVerticalStrainVolume( &
              x(1),x(2),x(3), &
              y(1,i),y(2,i),y(3,i),L(i),T(i),W(i),strike(i), &
              e11p,e12p,e13p,e22p,e23p,e33p,G,nu, &
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

  END SUBROUTINE computeTractionKernelsVerticalStrainVolume

  !------------------------------------------------------------------------
  !> subroutine computeStressKernelsVerticalStrainVolume
  !! calculates the traction kernels associated with strain in finite
  !! volumes in an elastic half-space using the analytic solution of 
  !!
  !!   Barbot S., J. D. P. Moore and V. Lambert, Displacement and Stress
  !!   Associated with Distributed Anelastic Deformation in a Half Space,
  !!   Bull. Seism. Soc. Am., 107(2), 10.1785/0120160237, 2017.
  !!
  !! INPUT:
  !! @param x1,x2,x3   - coordinates of observation points
  !! @param sv,dv,nv   - strike, dip and normal vector of observation points
  !! @param ns         - number of sources
  !! @param y          - position (NED) of vertical strain volume
  !! @param L,T,W      - length, thickness, and width of the strain volume
  !! @param strike     - strike of the vertical strain volume
  !! @param e11p,e12p
  !!        e13p,e22p
  !!        e23p,e33p  - strain components in the primed reference system
  !!                     tied to the strain volume
  !! @param G,lambda   - rigidity and Lame parameter
  !!
  !! OUTPUT:
  !! s11,s12,s13,s22,s23,s33 - the stress components in the reference
  !!                           system tied to the strain volume.
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeStressKernelsVerticalStrainVolume( &
                         x,sv,dv,nv, &
                         ns,y,L,T,W,strike, &
                         e11p,e12p,e13p,e22p,e23p,e33p,G,lambda, &
                         s11,s12,s13,s22,s23,s33)
    REAL*8, DIMENSION(3), INTENT(IN) :: x
    REAL*8, DIMENSION(3), INTENT(IN) :: sv,dv,nv
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: y
    REAL*8, DIMENSION(ns), INTENT(IN) :: L,T,W,strike
    REAL*8, INTENT(IN) :: e11p,e12p,e13p,e22p,e23p,e33p
    REAL*8, INTENT(IN) :: G,lambda
    REAL*8, DIMENSION(ns), INTENT(OUT) :: s11,s12,s13,s22,s23,s33

    INTEGER :: i
    REAL*8 :: s11p,s12p,s13p,s22p,s23p,s33p
    REAL*8 :: nu

    nu=lambda/(lambda+G)/2._8

    DO i=1,ns
       CALL computeStressVerticalStrainVolume( &
              x(1),x(2),x(3), &
              y(1,i),y(2,i),y(3,i),L(i),T(i),W(i),strike(i), &
              e11p,e12p,e13p,e22p,e23p,e33p,G,nu, &
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

  END SUBROUTINE computeStressKernelsVerticalStrainVolume

  !------------------------------------------------------------------------
  !> subroutine computeDisplacementKernelsVerticalStrainVolume
  !! calculates the displacement kernels associated with distributed
  !! anelastic strain in an elastic half-space using the analytic solution of
  !! of
  !!
  !!   Barbot S., J. D. P. Moore and V. Lambert, Displacement and Stress
  !!   Associated with Distributed Anelastic Deformation in a Half Space,
  !!   Bull. Seism. Soc. Am., 107(2), 10.1785/0120160237, 2017.
  !!
  !! INPUT:
  !! @param x1,x2,x3   - coordinates of observation points
  !! @param sv,dv,nv   - strike, dip and normal vector of observation points
  !! @param ns         - number of sources
  !! @param y          - position (NED) of vertical strain volume
  !! @param L,T,W      - length, thickness, and width of the strain volume
  !! @param strike     - strike of the vertical strain volume
  !! @param e11p,e12p
  !!        e13p,e22p
  !!        e23p,e33p  - strain components in the primed reference system
  !!                     tied to the strain volume
  !! @param G,lambda   - rigidity and Lame parameter
  !!
  !! OUTPUT:
  !! u1,u2,u3          - array. displacement in the strike, dip and normal 
  !!                     directions
  !!
  !! \author Sylvain Barbot (21/02/17) - original fortran form
  !------------------------------------------------------------------------
  SUBROUTINE computeDisplacementKernelsVerticalStrainVolume( &
                         x,&
                         ns,y,L,T,W,strike, &
                         e11p,e12p,e13p,e22p,e23p,e33p,G,lambda, &
                         u1,u2,u3)
    REAL*8, DIMENSION(3), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: y
    REAL*8, DIMENSION(ns), INTENT(IN) :: L,T,W,strike
    REAL*8, INTENT(IN) :: e11p,e12p,e13p,e22p,e23p,e33p
    REAL*8, INTENT(IN) :: G,lambda
    REAL*8, DIMENSION(ns), INTENT(OUT) :: u1,u2,u3

    INTEGER :: i
    REAL*8 :: nu

    nu=lambda/(lambda+G)/2._8

    DO i=1,ns
       CALL computeDisplacementVerticalStrainVolume( &
              x(1),x(2),x(3), &
              y(1,i),y(2,i),y(3,i),L(i),T(i),W(i),strike(i), &
              e11p,e12p,e13p,e22p,e23p,e33p,G,nu, &
              u1(i),u2(i),u3(i))

    END DO

  END SUBROUTINE computeDisplacementKernelsVerticalStrainVolume

END MODULE strainVolume




