
#ifndef UNICYCLE_MACROS

#define UNICYCLE_MACROS 1

#define ACOSH(X) REAL(LOG(X+ZSQRT(DCMPLX(X-1._8,0._8)*(X+1._8))),8)
#define ACOTH(X) 0.5_8*REAL(LOG(DCMPLX(X+1._8,0._8))-LOG(DCMPLX(X-1._8,0._8)),8)

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

! Paraview files
#define FPVTP 15

! data parallelism
#define FLAG_RECTANGLE 1
#define FLAG_TRIANGLE  2
#define FLAG_VOLUME    3

!-----------------------
! degrees of freedom
!-----------------------

! number of slip components on patches
#define DGF_PATCH 2

! number of strain components in volumes
#define DGF_VOLUME 6

! number of traction components on patches
#define DGF_VECTOR 3

! number of stress components on volumes
#define DGF_TENSOR 6

!-----------------------
! index in state vector
!-----------------------

! patches
#define STATE_VECTOR_SLIP_STRIKE     0
#define STATE_VECTOR_SLIP_DIP        1
#define STATE_VECTOR_TRACTION_STRIKE 2
#define STATE_VECTOR_TRACTION_DIP    3
#define STATE_VECTOR_TRACTION_NORMAL 4
#define STATE_VECTOR_STATE_1         5
#define STATE_VECTOR_VELOCITY        6
#define STATE_VECTOR_DGF_PATCH       7

! volumes
#define STATE_VECTOR_E11         0
#define STATE_VECTOR_E12         1
#define STATE_VECTOR_E13         2
#define STATE_VECTOR_E22         3
#define STATE_VECTOR_E23         4
#define STATE_VECTOR_E33         5
#define STATE_VECTOR_S11         6
#define STATE_VECTOR_S12         7
#define STATE_VECTOR_S13         8
#define STATE_VECTOR_S22         9
#define STATE_VECTOR_S23        10
#define STATE_VECTOR_S33        11
#define STATE_VECTOR_DGF_VOLUME 12

!-------------------------------------
! index in traction and stress vector
!-------------------------------------

! patches
#define TRACTION_VECTOR_STRIKE 0
#define TRACTION_VECTOR_DIP    1
#define TRACTION_VECTOR_NORMAL 2

! volumes
#define TRACTION_VECTOR_S11     0
#define TRACTION_VECTOR_S12     1
#define TRACTION_VECTOR_S13     2
#define TRACTION_VECTOR_S22     3
#define TRACTION_VECTOR_S23     4
#define TRACTION_VECTOR_S33     5

!-------------------------------------
! index in displacement vector
!-------------------------------------

#define DISPLACEMENT_VECTOR_NORTH 0
#define DISPLACEMENT_VECTOR_EAST  1
#define DISPLACEMENT_VECTOR_DOWN  2
#define DISPLACEMENT_VECTOR_DGF   3

#endif

