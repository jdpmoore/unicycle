
#ifndef UNICYCLE_MACROS

#define UNICYCLE_MACROS 1

#define WRITE_DEBUG_INFO(e) WRITE (0,'("error ",I3.3," at line ",I5.5,", rank ",I4.4)') e,__LINE__,rank
#define WRITE_DEBUG_INFO_SERIAL(e) WRITE (0,'("error ",I3.3," at line ",I5.5)') e,__LINE__

! for the Runge-Kutta method
#define TINY 1.d-30

!-------------------
! device numbers
!-------------------

! standard output
#define STDOUT 6

! standard error
#define STDERR 0

! time and time step file
#define FPTIME 14

! data parallelism
#define FLAG_PATCH  1
#define FLAG_VOLUME 2

!-----------------------
! degrees of freedom
!-----------------------

! number of slip components on patches
#define DGF_PATCH 1

! number of strain components in volumes
#define DGF_VOLUME 2

! number of traction components on patches
#define DGF_VECTOR 1

! number of stress components in volumes
#define DGF_TENSOR 2

!-----------------------
! index in state vector
!-----------------------

! patches
#define STATE_VECTOR_SLIP_STRIKE     0
#define STATE_VECTOR_TRACTION_STRIKE 1
#define STATE_VECTOR_STATE_1         2
#define STATE_VECTOR_VELOCITY        3
#define STATE_VECTOR_DGF_PATCH       4

! volumes
#define STATE_VECTOR_E12        0
#define STATE_VECTOR_E13        1
#define STATE_VECTOR_S12        2
#define STATE_VECTOR_S13        3
#define STATE_VECTOR_DGF_VOLUME 4

!-------------------------------------
! index in traction and stress vector
!-------------------------------------

! patches
#define TRACTION_VECTOR_STRIKE 0

! volumes
#define TRACTION_VECTOR_S12     0
#define TRACTION_VECTOR_S13     1

!-------------------------------------
! index in displacement vector
!-------------------------------------

#define DISPLACEMENT_VECTOR_STRIKE 0
#define DISPLACEMENT_VECTOR_DGF    1

#endif


