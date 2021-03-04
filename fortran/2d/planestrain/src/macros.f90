
#ifndef UNICYCLE_MACROS

#define UNICYCLE_MACROS 1

#define WRITE_DEBUG_INFO(e) WRITE (0,'("error ",I3.3," at line ",I5.5,", rank ",I4.4)') e,__LINE__,rank
#define WRITE_DEBUG_INFO_SERIAL(e) WRITE (0,'("error ",I3.3," at line ",I5.5)') e,__LINE__

! for the Runge-Kutta method
#define TINY 1.d-30

! constant variables
#define ZERO 0.000000000000000000000000000000000000000000000000d0
#define  ONE 1.000000000000000000000000000000000000000000000000d0

!-------------------
! device numbers
!-------------------

! standard output
#define STDOUT 6

! standard error
#define STDERR 0

! rank list check output
#define FPRANK 52

! time and time step file
#define FPTIME 14

! data parallelism
#define FLAG_PATCH 1
#define FLAG_VOLUME 2

!-----------------------
! degrees of freedom
!-----------------------

! number of slip components on patches
#define DGF_PATCH 1

! number of strain components in volumes
#define DGF_VOLUME 3

! number of traction components on patches
#define DGF_VECTOR 2

! number of stress components in volumes
#define DGF_TENSOR 3

!-----------------------
! index in state vector
!-----------------------

! patches
#define STATE_VECTOR_SLIP_DIP        0
#define STATE_VECTOR_TRACTION_DIP    1
#define STATE_VECTOR_TRACTION_NORMAL 2
#define STATE_VECTOR_STATE_1         3
#define STATE_VECTOR_VELOCITY        4
#define STATE_VECTOR_DGF_PATCH       5

! volumes
#define STATE_VECTOR_E22          0
#define STATE_VECTOR_E23          1
#define STATE_VECTOR_E33          2
#define STATE_VECTOR_S22          3
#define STATE_VECTOR_S23          4
#define STATE_VECTOR_S33          5
#define STATE_VECTOR_DGF_VOLUME   6

!-------------------------------------
! index in traction and stress vector
!-------------------------------------

! patches
#define TRACTION_VECTOR_DIP     0
#define TRACTION_VECTOR_NORMAL  1

! volumes
#define TRACTION_VECTOR_S22     0
#define TRACTION_VECTOR_S23     1
#define TRACTION_VECTOR_S33     2

!-------------------------------------
! index in displacement vector
!-------------------------------------

#define DISPLACEMENT_VECTOR_EAST 0
#define DISPLACEMENT_VECTOR_DOWN 1
#define DISPLACEMENT_VECTOR_DGF  2

#endif

