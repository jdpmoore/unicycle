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
! Eric Lindsey, Jun Muto, & Sagar Masuti. (2019, September 25). 
! jdpmoore/unicycle: Unicycle (Version 1.0). Zenodo. 
! http://doi.org/10.5281/zenodo.5688288
!
!-----------------------------------------------------------------------

#include "macros.f90"

MODULE greens

  USE antiplane

  IMPLICIT NONE

  PUBLIC

  !------------------------------------------------------------------------
  !! The Green's function for traction and stress interaction amongst 
  !! dislocations and strain volumes has the following layout
  !!
  !!       / KK  KL \
  !!   G = |        |
  !!       \ LK  LL /
  !!
  !! where
  !!
  !!   KK is the matrix for traction on faults due to fault slip
  !!   KL is the matrix for traction on faults due to strain in finite volumes
  !!
  !!   LK is the matrix for stress in volumes due to fault slip
  !!   LL is the matrix for stress in volumes due to strain in finite volumes
  !!
  !! The functions provided in this module computes the matrices KK, KL
  !! LK, and LL separately so they can be subsequently combined.
  !!
  !! The Green's function for observation points interaction amongst 
  !! dislocations and strain volumes has the following layout
  !!
  !!       /        \
  !!   O = | OK  OL |
  !!       \        /
  !!
  !! where
  !!
  !!   OK is the matrix for displacements due to fault slip
  !!   OL is the matrix for displacements due to strain in finite volumes
  !!
  !------------------------------------------------------------------------


CONTAINS

  !-----------------------------------------------------------------------
  !> subroutine buildG
  !! Builds the stress interaction matrix G following the layout
  !! illustrated below
  !!
  !!
  !!     /             \   +------------------------------------------+
  !!     |             |   |                                          |
  !!     |*************|   |                      nPatch  * dPatch    |
  !! G = |*************|   | ( varying size ) * (         +         ) |
  !!     |             |   |                      nVolume * dVolume   |
  !!     |             |   |                                          |
  !!     \             /   +------------------------------------------+
  !!
  !! where each thread owns a few rows of G, marked schematically ****.
  !! The number of rows depend of the type of elements owned by the 
  !! current thread.
  !!
  !! \author Sylvain Barbot (03/08/2017)
  !----------------------------------------------------------------------
  SUBROUTINE buildG(in,layout,G)
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
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: KK,KL

    ! stress kernels
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: LK,LL
  
    ! identity matrix
    REAL*8, DIMENSION(2,2), PARAMETER :: &
               eye=RESHAPE( (/ 1._8,0._8, &
                               0._8,1._8 /), (/ 2,2 /))
  
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,csize,ierr)

    ! number of columns
    n=in%patch%ns*DGF_PATCH+ &
      in%strainVolume%ns*DGF_VOLUME

    ! number of rows in current thread
    m=layout%listForceN(1+rank)

    ALLOCATE(G(n,m),STAT=ierr)

    ! two types of traction sources
    ALLOCATE(KK(DGF_PATCH,in%patch%ns,DGF_VECTOR), &
             KL(DGF_VOLUME,in%strainVolume%ns,DGF_VECTOR),STAT=ierr)
    IF (0 /= ierr) STOP "could not allocate the traction kernels"

    ! two types of stress sources
    ALLOCATE(LK(DGF_PATCH,in%patch%ns,DGF_TENSOR), &
             LL(DGF_VOLUME,in%strainVolume%ns,DGF_TENSOR),STAT=ierr)
    IF (0 /= ierr) STOP "could not allocate the kernels"
      
    ! loop over elements owned by current thread
    k=1
    DO i=1,SIZE(layout%elementIndex)
       elementType= layout%elementType(i)
       elementIndex=layout%elementIndex(i)

       SELECT CASE (elementType)
       CASE (FLAG_PATCH)
          ! G(:,k:k+layout%elementForceDGF(i)-1) = Green's function
          CALL tractionRows(in%patch%xc(:,elementIndex), &
                         in%patch%sv(:,elementIndex), &
                         in%patch%dv(:,elementIndex), &
                         in%patch%nv(:,elementIndex), &
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

    DEALLOCATE(KK,KL)
    DEALLOCATE(LK,LL)
      
  CONTAINS

    !-----------------------------------------------------------------------
    !> subroutine tractionRows
    !! evaluates all the relevant rows for traction change due to strike slip
    !! in patches and strain (12 and 13) in strain volumes.
    !!
    !! \author Sylvain Barbot (03/08/2017)
    !----------------------------------------------------------------------
    SUBROUTINE tractionRows(xc,sv,dv,nv,m,n,rows)
      IMPLICIT NONE
      REAL*8, DIMENSION(3), INTENT(IN) :: xc,sv,dv,nv
      INTEGER, INTENT(IN) :: m,n
      REAL*8, DIMENSION(n*m), INTENT(OUT) :: rows
  
      INTEGER :: j
  
      ! patches
      DO j=1,DGF_PATCH
         CALL computeTractionKernelsAntiplane(xc,sv,dv,nv, &
                 in%patch%ns,in%patch%x, &
                 in%patch%width, &
                 in%patch%dip, &
                 in%mu,KK(j,:,1))
      END DO
  
      ! strain volumes
      DO j=1,DGF_VOLUME
         CALL computeTractionKernelsStrainVolumeAntiplane(xc,sv,dv,nv, &
                 in%strainVolume%ns,in%strainVolume%x, &
                 in%strainVolume%thickness,in%strainVolume%width, &
                 in%strainVolume%dip, &
                 eye(j,1),eye(j,2),in%mu,KL(j,:,1))
      END DO
  
      ! concatenate kernels
      rows=(/ KK(:,:,1),KL(:,:,1) /)
  
    END SUBROUTINE tractionRows
  
    !-----------------------------------------------------------------------
    !> subroutine stressRows
    !! evaluates all the relevant rows for stress change due to strike slip
    !! in patches and strain (12 and 13) in finite volumes.
    !!
    !! \author Sylvain Barbot (06/09/2017)
    !----------------------------------------------------------------------
    SUBROUTINE stressRows(xc,sv,dv,nv,m,n,rows)
      IMPLICIT NONE
      REAL*8, DIMENSION(3), INTENT(IN) :: xc,sv,dv,nv
      INTEGER, INTENT(IN) :: m,n
      REAL*8, DIMENSION(n*m), INTENT(OUT) :: rows
  
      INTEGER :: j
  
      ! patches
      DO j=1,DGF_PATCH
         CALL computeStressKernelsAntiplane(xc,sv,dv,nv, &
                 in%patch%ns,in%patch%x, &
                 in%patch%width, &
                 in%patch%dip,in%mu,LK(j,:,1),LK(j,:,2))
      END DO
  
      ! strain volumes
      DO j=1,DGF_VOLUME
         CALL computeStressKernelsStrainVolumeAntiplane(xc,sv,dv,nv, &
                 in%strainVolume%ns,in%strainVolume%x, &
                 in%strainVolume%thickness,in%strainVolume%width, &
                 in%strainVolume%dip, &
                 eye(j,1),eye(j,2),in%mu,LL(j,:,1),LL(j,:,2))
      END DO
  
      ! concatenate kernels
      rows=(/ LK(:,:,1),LL(:,:,1), &
              LK(:,:,2),LL(:,:,2) /)
  
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

    ! identity matrix
    REAL*8, DIMENSION(2,2), PARAMETER :: &
               eye=RESHAPE( (/ 1._8,0._8, &
                               0._8,1._8 /), (/ 2,2 /))
  
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,csize,ierr)

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
          CASE (FLAG_PATCH)

             patch=in%patch%s(elementIndex)

             ! fault patches
             DO j=1,DGF_PATCH
                CALL computeDisplacementAntiplane( &
                        in%observationPoint(k)%x(2), &
                        in%observationPoint(k)%x(3), &
                        in%patch%x(2,elementIndex), &
                        in%patch%x(3,elementIndex), &
                        in%patch%width(elementIndex), &
                        in%patch%dip(elementIndex), &
                        1._8,O(l+j-1,p))

                ! contribution from faulting alone
                Of(l+j-1,p)=O(l+j-1,p)
             END DO

             l=l+DGF_PATCH

          CASE (FLAG_VOLUME)

             volume=in%strainVolume%s(elementIndex)

             ! strain volumes
             DO j=1,DGF_VOLUME
                CALL computeDisplacementStrainVolumeAntiplane( &
                        in%observationPoint(k)%x(2), &
                        in%observationPoint(k)%x(3), &
                        in%strainVolume%x(2,elementIndex), &
                        in%strainVolume%x(3,elementIndex), &
                        in%strainVolume%thickness(elementIndex), &
                        in%strainVolume%width(elementIndex), &
                        in%strainVolume%dip(elementIndex), &
                        eye(j,1),eye(j,2),O(l+j-1,p))

                ! contribution from flow alone
                Ol(l+j-1,p)=O(l+j-1,p)
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


