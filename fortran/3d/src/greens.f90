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

MODULE greens

  IMPLICIT NONE

  PUBLIC

  !------------------------------------------------------------------------
  !! The Green's function for traction and stress interaction amongst 
  !! triangular and rectangular dislocations, and strain volumes has the 
  !! following layout
  !!
  !!       / KK KT  KL \
  !!       |           |
  !!   G = | TK TT  TL |
  !!       |           |
  !!       \ LK LT  LL /
  !!
  !! where
  !!
  !!   KK is the matrix for traction on rectangular faults due to rectangular fault slip
  !!   KT is the matrix for traction on rectangular faults due to triangular  fault slip
  !!   TK is the matrix for traction on triangular  faults due to rectangular fault slip
  !!   TT is the matrix for traction on triangular  faults due to triangular  fault slip
  !!
  !!   KL is the matrix for traction on rectangular faults due to strain in volumes
  !!   TL is the matrix for traction on triangular  faults due to strain in volumes
  !!
  !!   LK is the matrix for stress in volumes due to rectangular fault slip
  !!   LT is the matrix for stress in volumes due to triangular  fault slip
  !!   LL is the matrix for stress in volumes due to strain in volumes
  !!
  !! The functions provided in this module computes the matrices KK, KT, KL
  !! TK, TT, TL, LK, LT, LL separately so they can be subsequently combined.
  !------------------------------------------------------------------------


CONTAINS

  !-----------------------------------------------------------------------
  !> subroutine buildG
  !! Builds the stress interaction matrix G
  !----------------------------------------------------------------------
  SUBROUTINE buildG(in,layout,G)
    USE nikkhoo15
    USE okada92
    USE strainvolume
    USE types
    IMPLICIT NONE
    TYPE(SIMULATION_STRUCT), INTENT(IN) :: in
    TYPE(LAYOUT_STRUCT), INTENT(IN) :: layout
    REAL*8, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: G

    INCLUDE 'mpif.h'

    INTEGER :: ierr,rank,csize
    INTEGER :: elementType,elementIndex

    INTEGER :: i,k,n,m

    ! traction kernels
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: KR,KT,KL

    ! stress kernels
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: LR,LT,LL
  
    ! identity matrix
    REAL*8, DIMENSION(6,6), PARAMETER :: &
               eye=RESHAPE( (/ 1._8,0._8,0._8,0._8,0._8,0._8, &
                               0._8,1._8,0._8,0._8,0._8,0._8, &
                               0._8,0._8,1._8,0._8,0._8,0._8, &
                               0._8,0._8,0._8,1._8,0._8,0._8, &
                               0._8,0._8,0._8,0._8,1._8,0._8, &
                               0._8,0._8,0._8,0._8,0._8,1._8 /), (/ 6,6 /))
  
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,csize,ierr)

    ! number of columns
    n=in%rectangularPatch%ns*DGF_PATCH+ &
      in%triangularPatch%ns*DGF_PATCH+ &
      in%strainVolume%ns*DGF_VOLUME

    ! number of rows in current thread
    m=layout%listForceN(1+rank)

    ALLOCATE(G(n,m),STAT=ierr)

    ! three types of traction sources
    ALLOCATE(KR(DGF_PATCH,in%rectangularPatch%ns,DGF_VECTOR), &
             KT(DGF_PATCH,in%triangularPatch%ns,DGF_VECTOR), &
             KL(DGF_VOLUME,in%strainVolume%ns,DGF_VECTOR),STAT=ierr)
    IF (0 /= ierr) STOP "could not allocate the traction kernels"
      
    ! three types of stress sources
    ALLOCATE(LR(DGF_PATCH,in%rectangularPatch%ns,DGF_TENSOR), &
             LT(DGF_PATCH,in%triangularPatch%ns,DGF_TENSOR), &
             LL(DGF_VOLUME,in%strainVolume%ns,DGF_TENSOR),STAT=ierr)
    IF (0 /= ierr) STOP "could not allocate the kernels"
      
    ! loop over elements owned by current thread
    k=1
    DO i=1,SIZE(layout%elementIndex)
       elementType= layout%elementType(i)
       elementIndex=layout%elementIndex(i)

       SELECT CASE (elementType)
       CASE (FLAG_RECTANGLE)
          ! G(:,k:k+layout%elementForceDGF(i)-1) = Green's function
          CALL tractionRows(in%rectangularPatch%xc(:,elementIndex), &
                         in%rectangularPatch%sv(:,elementIndex), &
                         in%rectangularPatch%dv(:,elementIndex), &
                         in%rectangularPatch%nv(:,elementIndex), &
                         layout%elementForceDGF(i),n,G(:,k))
       CASE (FLAG_TRIANGLE)
          ! G(:,k:k+layout%elementForceDGF(i)-1)= Green's function
          CALL tractionRows(in%triangularPatch%xc(:,elementIndex), &
                         in%triangularPatch%sv(:,elementIndex), &
                         in%triangularPatch%dv(:,elementIndex), &
                         in%triangularPatch%nv(:,elementIndex), &
                         layout%elementForceDGF(i),n,G(:,k))
       CASE (FLAG_VOLUME)
          ! G(:,k:k+layout%elementForceDGF(i)-1)= Green's function
          CALL stressRows(in%strainVolume%xc(:,elementIndex), &
                         in%strainVolume%sv(:,elementIndex), &
                         in%strainVolume%dv(:,elementIndex), &
                         in%strainVolume%nv(:,elementIndex), &
                         layout%elementForceDGF(i),n,G(:,k))
       END SELECT

       k=k+layout%elementForceDGF(i)
    END DO

    DEALLOCATE(KR,KT,KL)
    DEALLOCATE(LR,LT,LL)
      
  CONTAINS

    !-----------------------------------------------------------------------
    !> subroutine tractionRows
    !! evaluates all the relevant rows for traction change due to slip (strike
    !! slip and dip slip) in rectangular and triangular patches, and strain
    !! (11, 12, 13, 22, 23, 33) in strain volumes.
    !!
    !! \author Sylvain Barbot (03/08/2017)
    !----------------------------------------------------------------------
    SUBROUTINE tractionRows(xc,sv,dv,nv,m,n,rows)
      IMPLICIT NONE
      REAL*8, DIMENSION(3), INTENT(IN) :: xc,sv,dv,nv
      INTEGER, INTENT(IN) :: m,n
      REAL*8, DIMENSION(n*m), INTENT(OUT) :: rows
  
      INTEGER :: j
  
      ! rectangular patches
      DO j=1,DGF_PATCH
         CALL computeTractionKernelsOkada92(xc,sv,dv,nv, &
                 in%rectangularPatch%ns,in%rectangularPatch%x, &
                 in%rectangularPatch%length,in%rectangularPatch%width, &
                 in%rectangularPatch%strike,in%rectangularPatch%dip, &
                 eye(j,1),eye(j,2),eye(j,3),in%mu,in%lambda, &
                 KR(j,:,1),KR(j,:,2),KR(j,:,3))
      END DO
  
      ! triangular patches
      DO j=1,DGF_PATCH
         CALL computeTractionKernelsNikkhoo15(xc,sv,dv,nv, &
                 in%triangularPatch%ns,in%triangularPatch%i1,in%triangularPatch%i2,in%triangularPatch%i3, &
                 in%triangularPatch%nVe,in%triangularPatch%v, &
                 eye(j,1),eye(j,2),eye(j,3),in%mu,in%lambda, &
                 KT(j,:,1),KT(j,:,2),KT(j,:,3))
      END DO
  
      ! strain volumes
      DO j=1,DGF_VOLUME
         CALL computeTractionKernelsVerticalStrainVolume(xc,sv,dv,nv, &
                 in%strainVolume%ns,in%strainVolume%x, &
                 in%strainVolume%length,in%strainVolume%thickness,in%strainVolume%width, &
                 in%strainVolume%strike, &
                 eye(j,1),eye(j,2),eye(j,3),eye(j,4),eye(j,5),eye(j,6),in%mu,in%lambda, &
                 KL(j,:,1),KL(j,:,2),KL(j,:,3))
      END DO
  
      ! concatenate kernels
      rows=(/ KR(:,:,1),KT(:,:,1),KL(:,:,1), &
              KR(:,:,2),KT(:,:,2),KL(:,:,2), &
              KR(:,:,3),KT(:,:,3),KL(:,:,3) /)
  
    END SUBROUTINE tractionRows
  
    !-----------------------------------------------------------------------
    !> subroutine stressRows
    !! evaluates all the relevant rows for stress change due to slip (strike
    !! slip and dip slip) in rectangular and triangular patches, and strain
    !! (11, 12, 13, 22, 23, 33) in strain volumes.
    !!
    !! \author Sylvain Barbot (03/08/2017)
    !----------------------------------------------------------------------
    SUBROUTINE stressRows(xc,sv,dv,nv,m,n,rows)
      IMPLICIT NONE
      REAL*8, DIMENSION(3), INTENT(IN) :: xc,sv,dv,nv
      INTEGER, INTENT(IN) :: m,n
      REAL*8, DIMENSION(n*m), INTENT(OUT) :: rows
  
      INTEGER :: j
  
      ! rectangular patches
      DO j=1,DGF_PATCH
         CALL computeStressKernelsOkada92(xc,sv,dv,nv, &
                 in%rectangularPatch%ns,in%rectangularPatch%x, &
                 in%rectangularPatch%length,in%rectangularPatch%width, &
                 in%rectangularPatch%strike,in%rectangularPatch%dip, &
                 eye(j,1),eye(j,2),eye(j,3),in%mu,in%lambda, &
                 LR(j,:,1),LR(j,:,2),LR(j,:,3),LR(j,:,4),LR(j,:,5),LR(j,:,6))
      END DO
  
      ! triangular patches
      DO j=1,DGF_PATCH
         CALL computeStressKernelsNikkhoo15(xc,sv,dv,nv, &
                 in%triangularPatch%ns,in%triangularPatch%i1,in%triangularPatch%i2,in%triangularPatch%i3, &
                 in%triangularPatch%nVe,in%triangularPatch%v, &
                 eye(j,1),eye(j,2),eye(j,3),in%mu,in%lambda, &
                 LT(j,:,1),LT(j,:,2),LT(j,:,3),LT(j,:,4),LT(j,:,5),LT(j,:,6))
      END DO
  
      ! strain volumes
      DO j=1,DGF_VOLUME
         CALL computeStressKernelsVerticalStrainVolume(xc,sv,dv,nv, &
                 in%strainVolume%ns,in%strainVolume%x, &
                 in%strainVolume%length,in%strainVolume%thickness,in%strainVolume%width, &
                 in%strainVolume%strike, &
                 eye(j,1),eye(j,2),eye(j,3),eye(j,4),eye(j,5),eye(j,6),in%mu,in%lambda, &
                 LL(j,:,1),LL(j,:,2),LL(j,:,3),LL(j,:,4),LL(j,:,5),LL(j,:,6))
      END DO
  
      ! concatenate kernels
      rows=(/ LR(:,:,1),LT(:,:,1),LL(:,:,1), &
              LR(:,:,2),LT(:,:,2),LL(:,:,2), &
              LR(:,:,3),LT(:,:,3),LL(:,:,3), &
              LR(:,:,4),LT(:,:,4),LL(:,:,4), &
              LR(:,:,5),LT(:,:,5),LL(:,:,5), &
              LR(:,:,6),LT(:,:,6),LL(:,:,6) /)
  
    END SUBROUTINE stressRows
  
  END SUBROUTINE buildG
  
  !-----------------------------------------------------------------------
  !> subroutine buildO
  !! Builds the displacement matrix O following the layout below
  !!
  !!
  !!     /         **  \   +----------------------------------------+
  !!     |         **  |   |                                        |
  !!     |         **  |   |   nObservation                         |
  !! O = |         **  |   | (      *       )  * (  varying size  ) |
  !!     |         **  |   |   DISPLACEMENT                         |
  !!     |         **  |   |   _VECTOR_DGF                          |
  !!     |         **  |   |                                        |
  !!     \         **  /   +----------------------------------------+
  !!
  !! where each thread owns a few columns of O, marked schematically ****.
  !! The number of columns depends of the type of element owned by the
  !! current thread.
  !!
  !! The layout was chosen to minimize the MPI communication assuming
  !! that 
  !!
  !!   nObservation*DISPLACEMENT_VECTOR_DGF < (nPatch*dPatch+nVolume*dVolume).
  !!
  !! \author Sylvain Barbot (06/28/2017)
  !----------------------------------------------------------------------
  SUBROUTINE buildO(in,layout,O,Of,Ol)
    USE nikkhoo15
    USE okada92
    USE strainvolume
    USE types

    IMPLICIT NONE

    TYPE(SIMULATION_STRUCT), INTENT(IN) :: in
    TYPE(LAYOUT_STRUCT), INTENT(IN) :: layout
    REAL*8, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: O,Of,Ol

    INCLUDE 'mpif.h'

    INTEGER :: ierr,rank,csize
    INTEGER :: i,j,k,l,n,m,p
    INTEGER :: elementIndex,elementType

    TYPE(PATCH_ELEMENT_STRUCT) :: patch
    TYPE(STRAINVOLUME_ELEMENT_STRUCT) :: volume
    REAL*8 :: nu

    ! identity matrix
    REAL*8, DIMENSION(6,6), PARAMETER :: &
               eye=RESHAPE( (/ 1._8,0._8,0._8,0._8,0._8,0._8, &
                               0._8,1._8,0._8,0._8,0._8,0._8, &
                               0._8,0._8,1._8,0._8,0._8,0._8, &
                               0._8,0._8,0._8,1._8,0._8,0._8, &
                               0._8,0._8,0._8,0._8,1._8,0._8, &
                               0._8,0._8,0._8,0._8,0._8,1._8 /), (/ 6,6 /))
  
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,csize,ierr)

    ! Poisson's ratio
    nu=in%lambda/(in%lambda+in%mu)/2.d0

    ! number of columns in current thread
    n=layout%listVelocityN(rank+1)

    ! number of rows
    m=in%nObservationPoint*DISPLACEMENT_VECTOR_DGF

    ALLOCATE(O(n,m),Of(n,m),Ol(n,m),STAT=ierr)
    IF (0 /= ierr) STOP "could not allocate the displacement kernels"

    ! zero out partial contribution matrices
    Of=0._8
    Ol=0._8

    ! initiate counter for displacement vector
    p=1
    ! loop over observation points
    DO k=1,in%nObservationPoint

       ! initiate counter for source kinematics
       l=1

       ! loop over elements owned by current thread
       DO i=1,SIZE(layout%elementIndex)
          elementType= layout%elementType(i)
          elementIndex=layout%elementIndex(i)

          SELECT CASE (elementType)
          CASE (FLAG_RECTANGLE)

             patch=in%rectangularPatch%s(elementIndex)

             ! fault patches
             DO j=1,DGF_PATCH
                CALL computeDisplacementOkada92( &
                        in%observationPoint(k)%x(1), &
                        in%observationPoint(k)%x(2), &
                        in%observationPoint(k)%x(3), &
                        in%rectangularPatch%x(1,elementIndex), &
                        in%rectangularPatch%x(2,elementIndex), &
                        in%rectangularPatch%x(3,elementIndex), &
                        in%rectangularPatch%length(elementIndex), &
                        in%rectangularPatch%width(elementIndex), &
                        in%rectangularPatch%strike(elementIndex), &
                        in%rectangularPatch%dip(elementIndex), &
                        eye(j,1),eye(j,2),eye(j,3), &
                        in%mu,in%lambda, &
                        O(l+j-1,p),O(l+j-1,p+1),O(l+j-1,p+2))

                ! contribution from faulting alone
                Of(l+j-1,p  )=O(l+j-1,p)
                Of(l+j-1,p+1)=O(l+j-1,p+1)
                Of(l+j-1,p+2)=O(l+j-1,p+2)
             END DO

             l=l+DGF_PATCH

          CASE (FLAG_TRIANGLE)

             patch=in%triangularPatch%s(elementIndex)

             ! fault patches
             DO j=1,DGF_PATCH
                CALL computeDisplacementNikkhoo15( &
                        in%observationPoint(k)%x(1), &
                        in%observationPoint(k)%x(2), &
                        in%observationPoint(k)%x(3), &
                        in%triangularPatch%v(:,in%triangularPatch%i1(elementIndex)), &
                        in%triangularPatch%v(:,in%triangularPatch%i2(elementIndex)), &
                        in%triangularPatch%v(:,in%triangularPatch%i3(elementIndex)), &
                        eye(j,1),eye(j,2),eye(j,3), &
                        in%mu,in%lambda, &
                        O(l+j-1,p),O(l+j-1,p+1),O(l+j-1,p+2))

                ! contribution from faulting alone
                Of(l+j-1,p  )=O(l+j-1,p)
                Of(l+j-1,p+1)=O(l+j-1,p+1)
                Of(l+j-1,p+2)=O(l+j-1,p+2)
             END DO

             l=l+DGF_PATCH

          CASE (FLAG_VOLUME)

             volume=in%strainVolume%s(elementIndex)

             ! strain volumes
             DO j=1,DGF_VOLUME
                CALL computeDisplacementVerticalStrainVolume( &
                        in%observationPoint(k)%x(1), &
                        in%observationPoint(k)%x(2), &
                        in%observationPoint(k)%x(3), &
                        in%strainVolume%x(1,elementIndex), &
                        in%strainVolume%x(2,elementIndex), &
                        in%strainVolume%x(3,elementIndex), &
                        in%strainVolume%length(elementIndex), &
                        in%strainVolume%thickness(elementIndex), &
                        in%strainVolume%width(elementIndex), &
                        in%strainVolume%strike(elementIndex), &
                        eye(j,1),eye(j,2),eye(j,3), &
                        eye(j,4),eye(j,5),eye(j,6), &
                        in%mu,nu, &
                        O(l+j-1,p),O(l+j-1,p+1),O(l+j-1,p+2))

                ! contribution from flow alone
                Ol(l+j-1,p  )=O(l+j-1,p  )
                Ol(l+j-1,p+1)=O(l+j-1,p+1)
                Of(l+j-1,p+2)=O(l+j-1,p+2)
             END DO

             l=l+DGF_VOLUME

          CASE DEFAULT
             WRITE(STDERR,'("wrong case: this is a bug.")')
             WRITE_DEBUG_INFO(-1)
             STOP -1
          END SELECT

       END DO ! elements owned by current thread

       p=p+DISPLACEMENT_VECTOR_DGF

    END DO ! observation points
  
  END SUBROUTINE buildO
  
END MODULE greens
