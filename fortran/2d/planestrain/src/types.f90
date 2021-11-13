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

MODULE types

  IMPLICIT NONE

  TYPE PATCH_ELEMENT_STRUCT
     SEQUENCE
     REAL*8 :: Vpl
     REAL*8 :: tau0,mu0,sig,a,b,L,Vo,damping
  END TYPE PATCH_ELEMENT_STRUCT

  TYPE PATCH_STRUCT
     SEQUENCE
     INTEGER :: ns
     TYPE(PATCH_ELEMENT_STRUCT), DIMENSION(:), ALLOCATABLE :: s
     REAL*8, DIMENSION(:,:), ALLOCATABLE :: sv,dv,nv
     REAL*8, DIMENSION(:,:), ALLOCATABLE :: x,xc
     REAL*8, DIMENSION(:), ALLOCATABLE :: width,dip
  END TYPE PATCH_STRUCT

  TYPE STRESS_STRUCT
     SEQUENCE
     REAL*8 :: s22,s23,s33
  END TYPE STRESS_STRUCT

  TYPE STRAINVOLUME_ELEMENT_STRUCT
     SEQUENCE
     ! background strain rate
     REAL*8 :: e22,e23,e33
     ! transient (Kelvin) linear properties
     REAL*8 :: Gk,gammadot0k,dok,mk,Qk,Rk
     ! steady-state (Maxwell) linear properties
     REAL*8 :: Gm,gammadot0m,dom,mm,Qm,Rm
     ! transient (Kelvin) nonlinear properties
     REAL*8 :: nGk,ngammadot0k,npowerk,nQk,nRk
     ! steady-state (Maxwell) nonlinear properties
     REAL*8 :: nGm,ngammadot0m,npowerm,nQm,nRm
     ! thermal properties
     REAL*8 :: rhoc,To
  END TYPE STRAINVOLUME_ELEMENT_STRUCT

  TYPE OBSERVATION_POINT_STRUCT
     SEQUENCE
     INTEGER :: file
     REAL*8, DIMENSION(3) :: x
     CHARACTER(LEN=10) :: name
  END TYPE OBSERVATION_POINT_STRUCT

  TYPE STRAINVOLUME_STRUCT
     SEQUENCE
     INTEGER :: ns

     ! number of linear Maxwell volumes
     INTEGER :: nMaxwell

     ! number of linear Kelvin volumes
     INTEGER :: nKelvin

     ! number of nonlinear Kelvin volumes
     INTEGER :: nNonlinearKelvin

     ! number of nonlinear Maxwell volumes
     INTEGER :: nNonlinearMaxwell

     ! array of strain volume properties
     TYPE(STRAINVOLUME_ELEMENT_STRUCT), DIMENSION(:), ALLOCATABLE :: s

     ! strike, dip, and normal vectors
     REAL*8, DIMENSION(:,:), ALLOCATABLE :: sv,dv,nv

     ! upper left central position and center position
     REAL*8, DIMENSION(:,:), ALLOCATABLE :: x,xc

     ! dimension and orientation
     REAL*8, DIMENSION(:), ALLOCATABLE :: thickness,width,dip

     ! initial stress
     TYPE(STRESS_STRUCT), DIMENSION(:), ALLOCATABLE :: s0

  END TYPE STRAINVOLUME_STRUCT

  TYPE :: SIMULATION_STRUCT
     ! rigidity
     REAL*8 :: mu

     ! Lame parameter
     REAL*8 :: lambda

     ! simulation time
     REAL*8 :: interval

     ! number of patches (rectangular and triangular)
     INTEGER :: nPatch

     ! number of degrees of freedom for faults
     INTEGER :: dPatch

     ! number of strain volumes
     INTEGER :: nVolume

     ! number of degrees of freedom for strain volumes
     INTEGER :: dVolume

     ! rectangular patches
     TYPE(PATCH_STRUCT) :: patch

     ! cuboid strain volumes 
     TYPE(STRAINVOLUME_STRUCT) :: strainVolume

     ! output directory
     CHARACTER(256) :: wdir

     ! Greens function directory
     CHARACTER(256) :: greensFunctionDirectory

     ! filenames
     CHARACTER(256) :: timeFilename

     ! number of observation states
     INTEGER :: nObservationState

     ! observation patches
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: observationState

     ! number of observation points
     INTEGER :: nObservationPoint

     ! observation points
     TYPE(OBSERVATION_POINT_STRUCT), DIMENSION(:), ALLOCATABLE :: observationPoint

     ! observation points name
     CHARACTER(LEN=4), DIMENSION(:), ALLOCATABLE :: observationPointName

     ! other options
     LOGICAL :: isdryrun=.FALSE.
     LOGICAL :: isexportgreens=.FALSE.
     LOGICAL :: isexportnetcdf=.FALSE.
     LOGICAL :: ishelp=.FALSE.
     LOGICAL :: isversion=.FALSE.

  END TYPE SIMULATION_STRUCT

  TYPE LAYOUT_STRUCT
     ! list of number of elements in threads
     INTEGER, DIMENSION(:), ALLOCATABLE :: listElements,listOffset

     ! type of elements (FLAG_PATCH, FLAG_VOLUME)
     INTEGER, DIMENSION(:), ALLOCATABLE :: elementType

     ! index of elements
     INTEGER, DIMENSION(:), ALLOCATABLE :: elementIndex

     ! index of state indices
     INTEGER, DIMENSION(:), ALLOCATABLE :: elementStateIndex

     ! degrees of freedom for elements in velocity vector
     !   1 (dip slip) for patches
     INTEGER, DIMENSION(:), ALLOCATABLE :: elementVelocityDGF

     ! degrees of freedom for elements in state vector:
     !   5 (dip and normal traction, dip slip, velocity, state) for patches
     INTEGER, DIMENSION(:), ALLOCATABLE :: elementStateDGF

     ! degrees of freedom for elements in Greens function
     !   2 (dip and normal traction) for patches
     INTEGER, DIMENSION(:), ALLOCATABLE :: elementForceDGF

     ! number of traction and stress element for each thread and array position
     INTEGER, DIMENSION(:), ALLOCATABLE :: listForceN

     ! number of velocity and strain rate components for each thread and array position
     INTEGER, DIMENSION(:), ALLOCATABLE :: listVelocityN,listVelocityOffset

     ! number of state parameters for each thread and array position
     INTEGER, DIMENSION(:), ALLOCATABLE :: listStateN,listStateOffset

  END TYPE LAYOUT_STRUCT

END MODULE types

