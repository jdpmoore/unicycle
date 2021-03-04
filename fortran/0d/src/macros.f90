
#ifndef UNICYCLE_MACROS

#define UNICYCLE_MACROS 1

#define WRITE_DEBUG_INFO(e) WRITE (0,'("error ",I3.3," at line ",I5.5)') e,__LINE__

! for the Runge-Kutta method
#define TINY 1.d-30

! trigonometry
#define PI 3.1415926535897932384626433832795028d0

!-------------------
! device numbers
!-------------------

! standard output
#define STDOUT 6

! standard error
#define STDERR 0

! time and time step file
#define FPTIME 14

! time series of state vector
#define FPSTATE 100

!-----------------------
! index in state vector
!-----------------------

! patches
#define STATE_VECTOR_SLIP     1
#define STATE_VECTOR_TRACTION 2
#define STATE_VECTOR_STATE_1  3
#define STATE_VECTOR_VELOCITY 4
#define STATE_VECTOR_DGF      4

#endif


