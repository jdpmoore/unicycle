!-----------------------------------------------------------------------
! Copyright 2017 Sylvain Barbot
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

MODULE planestrain

  IMPLICIT NONE

  REAL*8, PARAMETER, PRIVATE :: pi = 3.141592653589793_8

CONTAINS

  !------------------------------------------------------------------------
  !> subroutine computeReferenceSystemPlaneStrain
  !! computes the center position and local reference system tied to the patch
  !!
  !! INPUT:
  !! @param ns
  !! @param x           - upper left coordinate of fault patch (north, east, down)
  !! @param dip         - dip angle (radian)
  !! @param W           - length and width of the dislocation
  !!
  !! OUTPUT:
  !! @param sv,dv,nv    - strike, dip and normal vectors of the fault patch
  !! @param xc          - coordinates (north, east, down) of the center
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeReferenceSystemPlaneStrain(ns,x,W,dip,sv,dv,nv,xc)
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
                
                
  END SUBROUTINE computeReferenceSystemPlaneStrain

  !------------------------------------------------------------------------
  !> subroutine computeStressPlaneStrain
  !! calculates the stress associated with a dislocation in an elastic 
  !! half-space in plane strain condition.
  !!
  !! INPUT:
  !! @param x2,x3    - coordinates of observation points
  !! @param s        - normal slip
  !! @param q2,q3    - coordinates of upper left corner of rectangular
  !!                   dislocation
  !! @param W        - width of the line dislocation
  !! @param dipd     - dip of the line dislocation (radian)
  !! @param G        - rigidity
  !! @param lambda   - Lame parameter
  !!
  !! OUTPUT:
  !! s22,s23,s33     - the stress components
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeStressPlaneStrain(x2,x3,q2,q3,W,dip,s,G,lambda,s22,s23,s33)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: s
    REAL*8, INTENT(IN) :: x2,x3,q2,q3
    REAL*8, INTENT(IN) :: W,dip
    REAL*8, INTENT(IN) :: G,lambda
    REAL*8, INTENT(OUT) :: s22,s23,s33

    REAL*8 :: T,eps23p

    ! infinitesimal thickness
    T=W/1e4
    ! infinite shear strain, positive for thrust
    eps23p=-s/(2.0_8*T)

    CALL computeStressStrainVolumePlaneStrain( &
            x2,x3,q2,q3,T,W,dip, &
            0._8,eps23p,0._8,G,lambda,s22,s23,s33)

  END SUBROUTINE computeStressPlaneStrain

  !------------------------------------------------------------------------
  !> subroutine computeDisplacementPlaneStrain
  !! calculates the displacements associated with a dislocation in an
  !! elastic half-space in plane strain condition.
  !!
  !! INPUT:
  !! @param x2,x3    - coordinates of observation points
  !! @param s        - normal slip
  !! @param q2,q3    - coordinates of upper left corner of line dislocation
  !! @param W        - width of the line dislocation
  !! @param dip      - dip angle of the line dislocation (radian)
  !! @param s        - dip slip
  !! @param G        - rigidity
  !!
  !! OUTPUT:
  !! u2,u3           - the displacement components in the east and depth
  !!                   directions.
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeDisplacementPlaneStrain(x2,x3, &
                          q2,q3,W,dip,s,G,lambda,u2,u3)
    REAL*8, INTENT(IN) :: x2,x3,q2,q3
    REAL*8, INTENT(IN) :: W,dip
    REAL*8, INTENT(IN) :: s,G,lambda
    REAL*8, INTENT(OUT) :: u2,u3

    REAL*8 :: T,eps23p

    ! infinitesimal thickness
    T=W/1e4
    ! infinite shear strain
    eps23p=-s/(2.0_8*T)

    CALL computeDisplacementStrainVolumePlaneStrain( &
            x2,x3,q2,q3,T,W,dip, &
            0._8,eps23p,0._8,G,lambda,u2,u3)

  END SUBROUTINE computeDisplacementPlaneStrain

  !------------------------------------------------------------------------
  !> subroutine computeDisplacementKernelsPlaneStrain
  !! calculates the displacement kernels associated with dislocations in
  !! an elastic half-space in plane strain condition.
  !!
  !! INPUT:
  !! @param x          - coordinates of observation points
  !! @param ns         - number of sources
  !! @param y          - coordinate (NED) of sources
  !! @param W          - width of the line dislocation
  !! @param dip        - dip of the line dislocation
  !!
  !! OUTPUT:
  !! u2,u3             - array, displacement in the east and depth directions.
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeDisplacementKernelsPlaneStrain(x, &
                        ns,y,W,dip,G,lambda, &
                        u2,u3)
    REAL*8, DIMENSION(3), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: y
    REAL*8, DIMENSION(ns), INTENT(IN) :: W,dip
    REAL*8, INTENT(IN) :: G,lambda
    REAL*8, DIMENSION(ns), INTENT(OUT) :: u2,u3

    INTEGER :: i

    DO i=1,ns
       CALL computeDisplacementPlaneStrain( &
                      x(2),x(3), &
                      y(2,i),y(3,i),W(i),dip(i), &
                      1._8,G,lambda,u2(i),u3(i))
    END DO

  END SUBROUTINE computeDisplacementKernelsPlaneStrain

  !------------------------------------------------------------------------
  !> subroutine computeDisplacementKernelsStrainVolumePlaneStrain
  !! calculates the displacement kernels associated with dislocations in
  !! an elastic half-space in plane strain condition.
  !!
  !! INPUT:
  !! @param x          - coordinates of observation points
  !! @param ns         - number of sources
  !! @param y          - coordinate (NED) of sources
  !! @param W          - width of the line dislocation
  !! @param dip        - dip of the line dislocation
  !!
  !! OUTPUT:
  !! u2,u3             - array, displacement in the east and depth directions.
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeDisplacementKernelsStrainVolumePlaneStrain(x, &
                        ns,y,T,W,dip, &
                        e22p,e23p,e33p, &
                        G,lambda, &
                        u2,u3)
    REAL*8, DIMENSION(3), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: y
    REAL*8, DIMENSION(ns), INTENT(IN) :: T,W,dip
    REAL*8, INTENT(IN) :: e22p,e23p,e33p
    REAL*8, INTENT(IN) :: G,lambda
    REAL*8, DIMENSION(ns), INTENT(OUT) :: u2,u3

    INTEGER :: i

    DO i=1,ns
       CALL computeDisplacementStrainVolumePlaneStrain( &
                      x(2),x(3), &
                      y(2,i),y(3,i),T(i),W(i),dip(i), &
                      e22p,e23p,e33p,G,lambda,u2(i),u3(i))
    END DO

  END SUBROUTINE computeDisplacementKernelsStrainVolumePlaneStrain

  !------------------------------------------------------------------------
  !> subroutine computeTractionKernelsPlaneStrain
  !! calculates the traction kernels associated with dislocations in an
  !! elastic half-space in plane strain condition.
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
  !! @param lambda     - Lame parameter
  !!
  !! OUTPUT:
  !! td,tn             - traction in the dip and normal directions.
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeTractionKernelsPlaneStrain(x,sv,dv,nv,ns,y,W,dip,G,lambda,td,tn)
    REAL*8, DIMENSION(3), INTENT(IN) :: x,sv,dv,nv
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: y
    REAL*8, DIMENSION(ns), INTENT(IN) :: W,dip
    REAL*8, INTENT(IN) :: G,lambda
    REAL*8, DIMENSION(ns), INTENT(OUT) :: td,tn

    INTEGER :: i
    REAL*8 :: s22,s23,s33

    DO i=1,ns
       CALL computeStressPlaneStrain( &
                      x(2),x(3), &
                      y(2,i),y(3,i),W(i),dip(i), &
                      1._8,G,lambda,s22,s23,s33)

       ! rotate to receiver system of coordinates
       td(i)= ( nv(2)*s22+nv(3)*s23 )*dv(2) &
             +( nv(2)*s23+nv(3)*s33 )*dv(3)
       tn(i)= ( nv(2)*s22+nv(3)*s23 )*nv(2) &
             +( nv(2)*s23+nv(3)*s33 )*nv(3)

    END DO

  END SUBROUTINE computeTractionKernelsPlaneStrain

  !------------------------------------------------------------------------
  !> subroutine computeStressKernelsPlaneStrain
  !! calculates the stress kernels associated with a dislocation in an
  !! elastic half-space in plane strain condition.
  !!
  !! INPUT:
  !! @param x          - coordinates of observation points
  !! @param sv,dv,nv   - strike, dip and normal vector of observation points
  !! @param ns         - number of sources
  !! @param y          - coordinate (NED) of sources
  !! @param W          - width of the dislocation
  !! @param dip        - dip of the dislocation
  !! @param s          - dip slip
  !! @param G          - rigidity
  !! @param lambda     - Lame parameter
  !!
  !! OUTPUT:
  !! s22,s23,s33       - the stress components
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeStressKernelsPlaneStrain(x,sv,dv,nv,ns,y,W,dip,lambda,G,s22,s23,s33)
    REAL*8, DIMENSION(3), INTENT(IN) :: x,sv,dv,nv
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: y
    REAL*8, DIMENSION(ns), INTENT(IN) :: W,dip
    REAL*8, INTENT(IN) :: G,lambda
    REAL*8, DIMENSION(ns), INTENT(OUT) :: s22,s23,s33

    INTEGER :: i
    REAL*8 :: s11p,s12p,s13p,s22p,s23p,s33p

    DO i=1,ns
       CALL computeStressPlaneStrain( &
                      x(2),x(3), &
                      y(2,i),y(3,i),W(i),dip(i), &
                      1._8,G,lambda,s22(i),s23(i),s33(i))

       ! rotate to receiver system of coordinates
       s22p= (  nv(2)*s22(i)+nv(3)*s23(i) )*nv(2) &
            +(  nv(2)*s23(i)+nv(3)*s33(i) )*nv(3)
       s23p=-(  nv(2)*s22(i)+nv(3)*s23(i) )*dv(2) &
            -(  nv(2)*s23(i)+nv(3)*s33(i) )*dv(3)
       s33p=-( -dv(2)*s22(i)-dv(3)*s23(i) )*dv(2) &
            -( -dv(2)*s23(i)-dv(3)*s33(i) )*dv(3)

       s22(i)=s22p
       s23(i)=s23p
       s33(i)=s33p

    END DO

  END SUBROUTINE computeStressKernelsPlaneStrain

  !------------------------------------------------------------------------
  !> subroutine computeTractionKernelsStrainVolumePlaneStrain
  !! calculates the traction kernels associated with a strain volumes
  !! in plane strain using the analytic solution of
  !!
  !!   Barbot S., J. D. P. Moore and V. Lambert, Displacement and Stress
  !!   Associated with Distributed Anelastic Deformation in a Half Space,
  !!   Bull. Seism. Soc. Am., 107(2), 10.1785/0120160237, 2017.
  !!
  !! INPUT:
  !! @param x2,x3      - coordinates of observation points
  !! @param sv,dv,nv   - strike, dip and normal vector of observation points
  !! @param ns         - number of sources
  !! @param y          - position (NED) of strain volume
  !! @param T,W        - thickness and width of the strain volume
  !! @param dip        - dip angle (radian) of the strain volume
  !! @param e22p
  !!        e23p
  !!        e33p       - strain components in the primed reference system
  !!                     tied to the strain volume
  !! @param G          - rigidity,
  !! @param lambda     - Lame parameter
  !!
  !! OUTPUT:
  !! td,tn             - traction in the dip and normal directions.
  !!
  !! \author Sylvain Barbot (06/09/17) - original form
  !------------------------------------------------------------------------
  SUBROUTINE computeTractionKernelsStrainVolumePlaneStrain( &
                         x,sv,dv,nv, &
                         ns,y,T,W,dip,e22p,e23p,e33p,G,lambda,td,tn)
    REAL*8, DIMENSION(3), INTENT(IN) :: x
    REAL*8, DIMENSION(3), INTENT(IN) :: sv,dv,nv
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: y
    REAL*8, DIMENSION(ns), INTENT(IN) :: T,W,dip
    REAL*8, INTENT(IN) :: e22p,e23p,e33p
    REAL*8, INTENT(IN) :: G,lambda
    REAL*8, DIMENSION(ns), INTENT(OUT) :: td,tn

    INTEGER :: i
    REAL*8 :: s22,s23,s33

    DO i=1,ns
       CALL computeStressStrainVolumePlaneStrain( &
              x(2),x(3), &
              y(2,i),y(3,i),T(i),W(i),dip(i), &
              e22p,e23p,e33p,G,lambda,s22,s23,s33)

       ! rotate to receiver system of coordinates
       td(i)= ( nv(2)*s22+nv(3)*s23 )*dv(2) &
             +( nv(2)*s23+nv(3)*s33 )*dv(3)
       tn(i)= ( nv(2)*s22+nv(3)*s23 )*nv(2) &
             +( nv(2)*s23+nv(3)*s33 )*nv(3)

    END DO

  END SUBROUTINE computeTractionKernelsStrainVolumePlaneStrain

  !------------------------------------------------------------------------
  !> subroutine computeStressKernelsStrainVolumePlaneStrain
  !! calculates the traction kernels associated with strain in finite
  !! volumes in plane strain using the analytic solution of 
  !!
  !!   Barbot S., J. D. P. Moore and V. Lambert, Displacement and Stress
  !!   Associated with Distributed Anelastic Deformation in a Half Space,
  !!   Bull. Seism. Soc. Am., 107(2), 10.1785/0120160237, 2017.
  !!
  !! INPUT:
  !! @param x2,x3      - coordinates of observation points
  !! @param sv,dv,nv   - strike, dip and normal vector of observation points
  !! @param ns         - number of sources
  !! @param y          - position (NED) of strain volume
  !! @param T,W        - thickness and width of the strain volume
  !! @param dip        - dip angle (radian) of the strain volume
  !! @param e22p, e23p 
  !!        e33p       - strain components in the primed reference system
  !!                     tied to the strain volume
  !! @param G          - rigidity
  !! @param lambda     - Lame parameter.
  !!
  !! OUTPUT:
  !! s22,s23,s33       - the stress components in the reference
  !!                     system tied to the strain volume.
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !------------------------------------------------------------------------
  SUBROUTINE computeStressKernelsStrainVolumePlaneStrain( &
                         x,sv,dv,nv, &
                         ns,y,T,W,dip, &
                         e22p,e23p,e33p,G,lambda,s22,s23,s33)
    REAL*8, DIMENSION(3), INTENT(IN) :: x
    REAL*8, DIMENSION(3), INTENT(IN) :: sv,dv,nv
    INTEGER, INTENT(IN) :: ns
    REAL*8, DIMENSION(3,ns), INTENT(IN) :: y
    REAL*8, DIMENSION(ns), INTENT(IN) :: T,W,dip
    REAL*8, INTENT(IN) :: e22p,e23p,e33p
    REAL*8, INTENT(IN) :: G,lambda
    REAL*8, DIMENSION(ns), INTENT(OUT) :: s22,s23,s33

    INTEGER :: i
    REAL*8 :: s22p,s23p,s33p

    DO i=1,ns
       CALL computeStressStrainVolumePlaneStrain( &
              x(2),x(3), &
              y(2,i),y(3,i),T(i),W(i),dip(i), &
              e22p,e23p,e33p,G,lambda,s22(i),s23(i),s33(i))

       ! rotate to receiver system of coordinates
       s22p= (  nv(2)*s22(i)+nv(3)*s23(i) )*nv(2) &
            +(  nv(2)*s23(i)+nv(3)*s33(i) )*nv(3)
       s23p=-(  nv(2)*s22(i)+nv(3)*s23(i) )*dv(2) &
            -(  nv(2)*s23(i)+nv(3)*s33(i) )*dv(3)
       s33p=-( -dv(2)*s22(i)-dv(3)*s23(i) )*dv(2) &
            -( -dv(2)*s23(i)-dv(3)*s33(i) )*dv(3)

       s22(i)=s22p
       s23(i)=s23p
       s33(i)=s33p

    END DO

  END SUBROUTINE computeStressKernelsStrainVolumePlaneStrain

END MODULE planestrain




