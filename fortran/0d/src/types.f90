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

MODULE types

  IMPLICIT NONE

  TYPE PATCH_ELEMENT_STRUCT
     SEQUENCE
     REAL*8 :: Vpl, width
     REAL*8 :: mu0,sig,a,b,L,Vo,damping
  END TYPE PATCH_ELEMENT_STRUCT

  TYPE PATCH_STRUCT
     SEQUENCE
  END TYPE PATCH_STRUCT

  TYPE :: SIMULATION_STRUCT
     ! rigidity
     REAL*8 :: mu

     ! simulation time
     REAL*8 :: interval

     ! rectangular patches
     TYPE(PATCH_ELEMENT_STRUCT) :: s

     ! output directory
     CHARACTER(256) :: wdir

     ! filenames
     CHARACTER(256) :: timeFilename

     ! other options
     LOGICAL :: isdryrun=.FALSE.
     LOGICAL :: ishelp=.FALSE.
     LOGICAL :: isversion=.FALSE.

  END TYPE SIMULATION_STRUCT

END MODULE types

