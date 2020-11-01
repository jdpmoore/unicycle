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

