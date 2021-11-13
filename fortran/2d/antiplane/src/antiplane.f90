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

MODULE antiplane

  IMPLICIT NONE

  REAL*8, PARAMETER, PRIVATE :: pi = 3.141592653589793_8

CONTAINS

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
  !> function Omega(x)
  !! evaluates the boxcar function
  !------------------------------------------------------------------------
  REAL*8 FUNCTION omega(x)
    REAL*8, INTENT(IN) :: x

    omega=heaviside(x+0.5_8)-heaviside(x-0.5_8)

  END FUNCTION omega

  !------------------------------------------------------------------------
  !> subroutine computeReferenceSystemAntiplane
  !! computes the center position and local reference system tied to the patch
  !!
  !! INPUT:
  !! @param ns
  !! @param x           - upper left coordinate of fault patch (north, east, down)
  !! @param dip         - dip angle (radian)
  !! @param W           - width of the dislocation
  !!
  !! OUTPUT:
  !! @param sv,dv,nv    - strike, dip and normal vectors of the fault patch
  !! @param xc          - coordinates (north, east, down) of the center
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeReferenceSystemAntiplane(ns,x,W,dip,sv,dv,nv,xc)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: x
    REAL*8, DIMENSION(ns), INTENT(IN) :: dip,W
    REAL*8, DIMENSION(3,ns), INTENT(OUT) :: sv,dv,nv,xc

    ! unit vectors in the strike direction
    sv(1,1:ns)=1._8
    sv(2,1:ns)=0._8
    sv(3,1:ns)=0._8
            
    ! unit vectors in the dip direction
    dv(1,1:ns)=+0._8
    dv(2,1:ns)=-COS(dip(1:ns))
    dv(3,1:ns)=-SIN(dip(1:ns))
            
    ! unit vectors in the normal direction
    nv(1,1:ns)=-0._8
    nv(2,1:ns)=+SIN(dip(1:ns))
    nv(3,1:ns)=-COS(dip(1:ns))
            
    ! center of fault patch
    xc(1,1:ns)=x(1,1:ns)-W(1:ns)/2*dv(1,1:ns)
    xc(2,1:ns)=x(2,1:ns)-W(1:ns)/2*dv(2,1:ns)
    xc(3,1:ns)=x(3,1:ns)-W(1:ns)/2*dv(3,1:ns)
                
                
  END SUBROUTINE computeReferenceSystemAntiplane

  !------------------------------------------------------------------------
  !> subroutine computeStressAntiplane
  !! calculates the stress associated with a dislocation in an elastic
  !! elastic half-space in antiplane strain condition.
  !!
  !! INPUT:
  !! @param x2,x3    - coordinates of observation points
  !! @param q2,q3    - coordinates of upper left corner of the dislocation
  !! @param W        - width of the line dislocation
  !! @param dip      - dip of the line dislocation (radian)
  !! @param G        - rigidity
  !!
  !! OUTPUT:
  !! s12,s13         - the stress components
  !!
  !! KNOWN BUGS: 
  !! the dip angle is ignored.
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeStressAntiplane(x2,x3, &
                                  q2,q3,W,dip, &
                                  G,s12,s13)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: x2,x3,q2,q3
    REAL*8, INTENT(IN) :: W,dip
    REAL*8, INTENT(IN) :: G
    REAL*8, INTENT(OUT) :: s12,s13

    REAL*8 :: r2

    r2=x2-q2

    s12=G*( &
        -(x3-q3)/(r2**2+(x3-q3)**2) &
        +(x3+q3)/(r2**2+(x3+q3)**2) &
        +(x3-q3-W)/(r2**2+(x3-q3-W)**2) &
        -(x3+q3+W)/(r2**2+(x3+q3+W)**2))/2._8/pi

    s13=G*( &
         r2/(r2**2+(x3-q3)**2) &
        -r2/(r2**2+(x3+q3)**2) &
        -r2/(r2**2+(x3-q3-W)**2) &
        +r2/(r2**2+(x3+q3+W)**2))/2/pi;

  END SUBROUTINE computeStressAntiplane

  !------------------------------------------------------------------------
  !> subroutine computeDisplacementAntiplane
  !! calculates the displacements associated with a dislocation in an
  !! elastic half-space in antiplane condition.
  !!
  !! INPUT:
  !! @param x2,x3    - coordinates of observation points
  !! @param q2,q3    - coordinates of upper left corner of line dislocation
  !! @param W        - width of the line dislocation
  !! @param dip      - dip of the line dislocation (radian)
  !! @param s        - strike slip
  !!
  !! OUTPUT:
  !! u1              - displacement in the strike direction
  !!
  !! KNOWN BUGS:
  !! the dip angle is ignored.
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeDisplacementAntiplane(x2,x3, &
                                        q2,q3,W,dip, &
                                        s,u1)
    REAL*8, INTENT(IN) :: x2,x3,q2,q3
    REAL*8, INTENT(IN) :: W,dip
    REAL*8, INTENT(IN) :: s
    REAL*8, INTENT(OUT) :: u1

    u1=s*(+atan((x3-q3)/(x2-q2)) &
          -atan((x3+q3)/(x2-q2)) &
          -atan((x3-q3-W)/(x2-q2)) &
          +atan((x3+q3+W)/(x2-q2)))/2._8/pi;

  END SUBROUTINE computeDisplacementAntiplane

  !------------------------------------------------------------------------
  !> subroutine computeDisplacementStrainVolumeAntiplane
  !! calculates the displacements associated with a dislocation in an
  !! elastic half-space in antiplane condition.
  !!
  !! INPUT:
  !! @param x2,x3    - coordinates of observation points
  !! @param q2,q3    - coordinates of upper left corner of line dislocation
  !! @param W        - width of the line dislocation
  !! @param dip      - dip of the line dislocation (radian)
  !! @param e12,e13  - anelastic strain in the horizontal and depth directions
  !!
  !! OUTPUT:
  !! u1              - displacement in the strike direction
  !!
  !! KNOWN BUGS:
  !! the dip angle is ignored.
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeDisplacementStrainVolumeAntiplane(x2,x3, &
                                        q2,q3,T,W,dip, &
                                        e12,e13,u1)
    REAL*8, INTENT(IN) :: x2,x3,q2,q3
    REAL*8, INTENT(IN) :: T,W,dip
    REAL*8, INTENT(IN) :: e12,e13
    REAL*8, INTENT(OUT) :: u1

    REAL*8 :: r2

    r2=x2-q2

    ! displacement due to distributed anelastic strain
    u1=e12/(2*pi)*( &
           (x3-q3-W)*log((r2-T/2)**2+(x3-q3-W)**2) &
          -(x3-q3-W)*log((r2+T/2)**2+(x3-q3-W)**2) &
          -(x3-q3)*log((r2-T/2)**2+(x3-q3)**2) &
          +(x3-q3)*log((r2+T/2)**2+(x3-q3)**2) &
          +2*(r2-T/2)*(atan((x3-q3-W)/(r2-T/2))-atan((x3-q3)/(r2-T/2))) &
          +2*(r2+T/2)*(atan((x3-q3)/(r2+T/2))-atan((x3-q3-W)/(r2+T/2))) &
          +(x3+q3+W)*log((r2+T/2)**2+(x3+q3+W)**2) &
          -(x3+q3+W)*log((r2-T/2)**2+(x3+q3+W)**2) &
          -(x3+q3)*log((r2+T/2)**2+(x3+q3)**2) &
          +(x3+q3)*log((r2-T/2)**2+(x3+q3)**2) &
          +2*(r2+T/2)*(atan((x3+q3+W)/(r2+T/2))-atan((x3+q3)/(r2+T/2))) &
          +2*(r2-T/2)*(atan((x3+q3)/(r2-T/2))-atan((x3+q3+W)/(r2-T/2))) &
       ) + &
       e13/(2*pi)*( &
           (r2-T/2)*log((r2-T/2)**2+(x3-q3-W)**2) &
          -(r2+T/2)*log((r2+T/2)**2+(x3-q3-W)**2) &
          -(r2-T/2)*log((r2-T/2)**2+(x3-q3)**2) &
          +(r2+T/2)*log((r2+T/2)**2+(x3-q3)**2) &
          +2*(x3-W-q3)*(atan((r2-T/2)/(x3-q3-W))-atan((r2+T/2)/(x3-q3-W))) &
          +2*(x3-q3)*(atan((r2+T/2)/(x3-q3))-atan((r2-T/2)/(x3-q3))) &
          +(r2-T/2)*log((r2-T/2)**2+(x3+q3+W)**2) &
          -(r2+T/2)*log((r2+T/2)**2+(x3+q3+W)**2) &
          -(r2-T/2)*log((r2-T/2)**2+(x3+q3)**2) &
          +(r2+T/2)*log((r2+T/2)**2+(x3+q3)**2) &
          +2*(x3+W+q3)*(atan((r2-T/2)/(x3+q3+W))-atan((r2+T/2)/(x3+q3+W))) &
          +2*(x3+q3)*(atan((r2+T/2)/(x3+q3))-atan((r2-T/2)/(x3+q3))) &
       )

  END SUBROUTINE computeDisplacementStrainVolumeAntiplane

  !------------------------------------------------------------------------
  !> subroutine computeDisplacementKernelsAntiplane
  !! calculates the displacement kernels associated with dislocations in
  !! an elastic half-space in antiplane condition.
  !!
  !! INPUT:
  !! @param x          - coordinates of observation points
  !! @param ns         - number of sources
  !! @param y          - coordinate (NED) of sources
  !! @param W          - width of the line dislocation
  !! @param dip        - dip of the line dislocation
  !!
  !! OUTPUT:
  !! u1                - array, displacement in the strike direction.
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeDisplacementKernelsAntiplane(x, &
                        ns,y,W,dip, &
                        u1)
    REAL*8, DIMENSION(3), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: y
    REAL*8, DIMENSION(ns), INTENT(IN) :: W,dip
    REAL*8, DIMENSION(ns), INTENT(OUT) :: u1

    INTEGER :: i

    DO i=1,ns
       CALL computeDisplacementAntiplane( &
                      x(2),x(3), &
                      y(2,i),y(3,i),W(i),dip(i), &
                      1._8,u1(i))
    END DO

  END SUBROUTINE computeDisplacementKernelsAntiplane

  !------------------------------------------------------------------------
  !> subroutine computeDisplacementKernelsStrainVolumeAntiplane
  !! calculates the displacement kernels associated with distributed
  !! anelastic deformation an elastic half-space in antiplane condition.
  !!
  !! INPUT:
  !! @param x          - coordinates of observation points
  !! @param ns         - number of sources
  !! @param y          - coordinate (NED) of sources
  !! @param T,W        - thickness and width of the strain volume
  !! @param dip        - dip of the line dislocation
  !!
  !! OUTPUT:
  !! u1                - array, displacement in the strike direction.
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeDisplacementKernelsStrainVolumeAntiplane(x, &
                        ns,y,T,W,dip, &
                        e12p,e13p,u1)
    REAL*8, DIMENSION(3), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: y
    REAL*8, DIMENSION(ns), INTENT(IN) :: T,W,dip
    REAL*8, INTENT(IN) :: e12p,e13p
    REAL*8, DIMENSION(ns), INTENT(OUT) :: u1

    INTEGER :: i

    DO i=1,ns
       CALL computeDisplacementStrainVolumeAntiplane( &
                      x(2),x(3), &
                      y(2,i),y(3,i),T(i),W(i),dip(i), &
                      e12p,e13p,u1(i))
    END DO

  END SUBROUTINE computeDisplacementKernelsStrainVolumeAntiplane

  !------------------------------------------------------------------------
  !> subroutine computeTractionKernelsAntiplane
  !! calculates the traction kernels associated with a dislocations in an
  !! elastic half-space in antiplane condition.
  !!
  !! INPUT:
  !! @param x          - coordinates of observation points
  !! @param sv,dv,nv   - strike, dip and normal vector of observation points
  !! @param ns         - number of sources
  !! @param y          - coordinate (NED) of sources
  !! @param W          - width of the line dislocation
  !! @param dip        - dip of the line dislocation
  !! @param s          - strike slip
  !! @param G          - rigidity
  !!
  !! OUTPUT:
  !! ts                - traction in the strike direction.
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeTractionKernelsAntiplane(x,sv,dv,nv, &
                        ns,y,W,dip, &
                        G, &
                        ts)
    REAL*8, DIMENSION(3), INTENT(IN) :: x,sv,dv,nv
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: y
    REAL*8, DIMENSION(ns), INTENT(IN) :: W,dip
    REAL*8, INTENT(IN) :: G
    REAL*8, DIMENSION(ns), INTENT(OUT) :: ts

    INTEGER :: i
    REAL*8 :: s12,s13

    DO i=1,ns
       CALL computeStressAntiplane( &
                      x(2),x(3), &
                      y(2,i),y(3,i),W(i),dip(i), &
                      G,s12,s13)

       ! rotate to receiver system of coordinates
       ts(i)= (           nv(2)*s12+nv(3)*s13 )*sv(1) &
             +( nv(1)*s12                     )*sv(2) &
             +( nv(1)*s13                     )*sv(3)

    END DO

  END SUBROUTINE computeTractionKernelsAntiplane

  !------------------------------------------------------------------------
  !> subroutine computeStressKernelsAntiplane
  !! calculates the stress kernels associated with a dislocation in an
  !! in an elastic half-space in antiplane condition.
  !!
  !! INPUT:
  !! @param x          - coordinates of observation points
  !! @param sv,dv,nv   - strike, dip and normal vector of observation points
  !! @param ns         - number of sources
  !! @param y          - coordinate (NED) of sources
  !! @param W          - width of the dislocation
  !! @param dip        - dip angle of the dislocation
  !! @param s          - strike slip
  !! @param G          - rigidity
  !!
  !! OUTPUT:
  !! s12,s13           - the stress components
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeStressKernelsAntiplane(x,sv,dv,nv, &
                        ns,y,W,dip, &
                        G,s12,s13)
    REAL*8, DIMENSION(3), INTENT(IN) :: x,sv,dv,nv
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: y
    REAL*8, DIMENSION(ns), INTENT(IN) :: W,dip
    REAL*8, INTENT(IN) :: G
    REAL*8, DIMENSION(ns), INTENT(OUT) :: s12,s13

    INTEGER :: i
    REAL*8 :: s12p,s13p

    DO i=1,ns
       CALL computeStressAntiplane( &
                      x(2),x(3), &
                      y(2,i),y(3,i),W(i),dip(i), &
                      G,s12(i),s13(i))

       ! rotate to receiver system of coordinates
       s12p= (              sv(2)*s12(i)+sv(3)*s13(i) )*nv(1) &
            +( sv(1)*s12(i)                           )*nv(2) &
            +( sv(1)*s13(i)                           )*nv(3)
       s13p=-(              sv(2)*s12(i)+sv(3)*s13(i) )*dv(1) &
            -( sv(1)*s12(i)                           )*dv(2) &
            -( sv(1)*s13(i)                           )*dv(3)

       s12(i)=s12p
       s13(i)=s13p

    END DO

  END SUBROUTINE computeStressKernelsAntiplane

  !------------------------------------------------------------------------
  !> subroutine ComputeStressStrainVolumeAntiplane computes the stress field 
  !! associated with deforming strain volume using the analytic solution of
  !!
  !!   Lambert, V., and S. Barbot. "Contribution of viscoelastic flow in 
  !!   earthquake cycles within the lithosphere‐asthenosphere system." 
  !!   Geophysical Research Letters 43.19 (2016).
  !!
  !! considering the following geometry:
  !!
  !!                              q2,q3
  !!                   +------------@-------------+---> E (x2)
  !!                   |                          |
  !!                   |                        w |
  !!                   |                        i |
  !!                   |                        d |
  !!                   |                        t |
  !!                   |                        h |
  !!                   |                          |
  !!                   +--------------------------+
  !!                   :     t h i c k n e s s 
  !!                   :
  !!                   |
  !!                   Z (x3)
  !!
  !!
  !! INPUT:
  !! @param x2, x3         easting, and depth of the observation point
  !!                       in unprimed system of coordinates.
  !! @param q2, q3         east and depth coordinates of the strain volume,
  !! @param T, W           thickness, and width of the strain volume,
  !! @param epsijp         anelastic strain component 12 and 13
  !!                       in the strain volume in the system of reference tied to 
  !!                       the strain volume (primed reference system),
  !! @param G              rigidity.
  !!
  !! OUTPUT:
  !! s12,s13               stress components in the unprimed reference system.
  !!
  !! KNOWN BUGS:
  !! the dip angle is ignored.
  !!
  !! \author Sylvain Barbot (06/09/17) - original form
  !------------------------------------------------------------------------
  SUBROUTINE computeStressStrainVolumeAntiplane(x2,x3,q2,q3,T,W,dip, &
                                eps12p,eps13p,G,s12,s13)

    REAL*8, INTENT(IN) :: x2,x3,q2,q3,T,W,dip
    REAL*8, INTENT(IN) :: eps12p,eps13p
    REAL*8, INTENT(IN) :: G
    REAL*8, INTENT(OUT) :: s12,s13
    
    REAL*8 :: r2

    r2=x2-q2

    s12= G/pi*eps12p*( &
           atan((x3-q3  )/(r2+T/2))-atan((x3-q3  )/(r2-T/2)) &
          +atan((x3-q3-W)/(r2-T/2))-atan((x3-q3-W)/(r2+T/2)) &
          -atan((x3+q3+W)/(r2-T/2))-atan((x3+q3  )/(r2+T/2)) &
          +atan((x3+q3  )/(r2-T/2))+atan((x3+q3+W)/(r2+T/2))) &
       + G/(2*pi)*eps13p*( &
           log((r2-T/2)**2+(x3-q3-W)**2) - log((r2+T/2)**2+(x3-q3-W)**2) &
          +log((r2-T/2)**2+(x3+q3+W)**2) - log((r2+T/2)**2+(x3+q3+W)**2) &
          -log((r2-T/2)**2+(x3-q3)**2)   + log((r2+T/2)**2+(x3-q3)**2) &
          -log((r2-T/2)**2+(x3+q3)**2)   + log((r2+T/2)**2+(x3+q3)**2)) &
       -2*G*eps12p*omega(r2/T)*omega((x3-(2*q3+W)/2._8)/W)

    s13= G/(2*pi)*eps12p*( &
           log((r2-T/2)**2+(x3-q3-W)**2) - log((r2+T/2)**2+(x3-q3-W)**2) &
          -log((r2-T/2)**2+(x3+q3+W)**2) + log((r2+T/2)**2+(x3+q3+W)**2) &
          -log((r2-T/2)**2+(x3-q3)**2)   + log((r2+T/2)**2+(x3-q3)**2) &
          +log((r2-T/2)**2+(x3+q3)**2)   - log((r2+T/2)**2+(x3+q3)**2)) &
       + G/pi*eps13p*( &
           atan((r2+T/2)/(x3-q3))  -atan((r2-T/2)/(x3-q3)) &
          -atan((r2+T/2)/(x3-q3-W))+atan((r2-T/2)/(x3-q3-W)) &
          +atan((r2+T/2)/(x3+q3))  -atan((r2-T/2)/(x3+q3)) &
          -atan((r2+T/2)/(x3+q3+W))+atan((r2-T/2)/(x3+q3+W))) &
       -2*G*eps13p*omega(r2/T)*omega((x3-(2*q3+W)/2)/W)

  END SUBROUTINE computeStressStrainVolumeAntiplane

  !------------------------------------------------------------------------
  !> subroutine computeTractionKernelsStrainVolumeAntiplane
  !! calculates the traction kernels associated with a strain volumes
  !! in antiplane strain using the analytic solution of
  !!
  !!   Lambert, V., and S. Barbot. "Contribution of viscoelastic flow in 
  !!   earthquake cycles within the lithosphere‐asthenosphere system." 
  !!   Geophysical Research Letters 43.19 (2016).
  !!
  !! INPUT:
  !! @param x2,x3      - coordinates of observation points
  !! @param sv,dv,nv   - strike, dip and normal vector of observation points
  !! @param ns         - number of sources
  !! @param y          - position (NED) of strain volume
  !! @param T,W        - thickness and width of the strain volume
  !! @param dip        - dip angle (radian) of the strain volume
  !! @param e12p
  !!        e13p       - strain components in the primed reference system
  !!                     tied to the strain volume
  !! @param G          - rigidity
  !!
  !! OUTPUT:
  !! ts                - traction in the strike direction.
  !!
  !! \author Sylvain Barbot (06/09/17) - original form
  !------------------------------------------------------------------------
  SUBROUTINE computeTractionKernelsStrainVolumeAntiplane( &
                         x,sv,dv,nv, &
                         ns,y,T,W,dip, &
                         e12p,e13p,G,ts)
    REAL*8, DIMENSION(3), INTENT(IN) :: x
    REAL*8, DIMENSION(3), INTENT(IN) :: sv,dv,nv
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: y
    REAL*8, DIMENSION(ns), INTENT(IN) :: T,W,dip
    REAL*8, INTENT(IN) :: e12p,e13p
    REAL*8, INTENT(IN) :: G
    REAL*8, DIMENSION(ns), INTENT(OUT) :: ts

    INTEGER :: i
    REAL*8 :: s12,s13

    DO i=1,ns
       CALL computeStressStrainVolumeAntiplane( &
              x(2),x(3), &
              y(2,i),y(3,i),T(i),W(i),dip(i), &
              e12p,e13p,G,s12,s13)

       ! rotate to receiver system of coordinates
       ts(i)= (           nv(2)*s12+nv(3)*s13 )*sv(1) &
             +( nv(1)*s12                     )*sv(2) &
             +( nv(1)*s13                     )*sv(3)

    END DO

  END SUBROUTINE computeTractionKernelsStrainVolumeAntiplane

  !------------------------------------------------------------------------
  !> subroutine computeStressKernelsStrainVolumeAntiplane
  !! calculates the traction kernels associated with strain in finite
  !! volumes in antiplane strain using the analytic solution of 
  !!
  !!   Lambert, V., and S. Barbot. "Contribution of viscoelastic flow in 
  !!   earthquake cycles within the lithosphere‐asthenosphere system." 
  !!   Geophysical Research Letters 43.19 (2016).
  !!
  !! INPUT:
  !! @param x2,x3      - coordinates of observation points
  !! @param sv,dv,nv   - strike, dip and normal vector of observation points
  !! @param ns         - number of sources
  !! @param y          - position (NED) of strain volume
  !! @param T,W        - thickness and width of the strain volume
  !! @param dip        - dip angle (radian) of the strain volume
  !! @param e12p, e13p - strain components in the primed reference system
  !!                     tied to the strain volume
  !! @param G          - rigidity.
  !!
  !! OUTPUT:
  !! s12,s13           - the stress components in the reference
  !!                     system tied to the strain volume.
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeStressKernelsStrainVolumeAntiplane( &
                         x,sv,dv,nv, &
                         ns,y,T,W,dip, &
                         e12p,e13p,G,s12,s13)
    REAL*8, DIMENSION(3), INTENT(IN) :: x
    REAL*8, DIMENSION(3), INTENT(IN) :: sv,dv,nv
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: y
    REAL*8, DIMENSION(ns), INTENT(IN) :: T,W,dip
    REAL*8, INTENT(IN) :: e12p,e13p
    REAL*8, INTENT(IN) :: G
    REAL*8, DIMENSION(ns), INTENT(OUT) :: s12,s13

    INTEGER :: i
    REAL*8 :: s12p,s13p

    DO i=1,ns
       CALL computeStressStrainVolumeAntiplane( &
              x(2),x(3), &
              y(2,i),y(3,i),T(i),W(i),dip(i), &
              e12p,e13p,G,s12(i),s13(i))

       ! rotate to receiver system of coordinates
       s12p= (              sv(2)*s12(i)+sv(3)*s13(i) )*nv(1) &
            +( sv(1)*s12(i)                           )*nv(2) &
            +( sv(1)*s13(i)                           )*nv(3)
       s13p=-(              sv(2)*s12(i)+sv(3)*s13(i) )*dv(1) &
            -( sv(1)*s12(i)                           )*dv(2) &
            -( sv(1)*s13(i)                           )*dv(3)

       s12(i)=s12p
       s13(i)=s13p

    END DO

  END SUBROUTINE computeStressKernelsStrainVolumeAntiplane

END MODULE antiplane




