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
! along with RELAX.  If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------

#include "macros.f90"

MODULE input

  USE types

  IMPLICIT NONE

  REAL*8, PARAMETER :: DEG2RAD = 0.01745329251994329547437168059786927_8

CONTAINS

  !---------------------------------------------------------------------
  !> subroutine init
  !! reads simulation parameters from the standard input and initialize
  !! model parameters.
  !!
  !! INPUT:
  !! @param unit - the unit number used to read input data
  !!
  !! OUTPUT:
  !! @param in
  !!
  !! \author Sylvain Barbot (sbarbot@ntu.edu.sg)
  !---------------------------------------------------------------------
  SUBROUTINE init(in)
    USE types
    USE getopt_m

    TYPE(SIMULATION_STRUCT), INTENT(OUT) :: in

    INCLUDE 'mpif.h'

    CHARACTER :: ch
    CHARACTER(512) :: dataline
    CHARACTER(256) :: filename
    INTEGER :: iunit,noptions
!$  INTEGER :: omp_get_num_procs,omp_get_max_threads
    TYPE(OPTION_S) :: opts(4)

    INTEGER :: i,j,k,iostatus,rank,size,ierr,position
    INTEGER :: maxVertexIndex=-1
    INTEGER, PARAMETER :: psize=1024
    INTEGER :: dummy
    CHARACTER, DIMENSION(psize) :: packed

    INTEGER :: nObservationPatch,nObservationVolume
    INTEGER, DIMENSION(:), ALLOCATABLE :: observationPatch,observationVolume
  
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)

    ! define long options, such as --dry-run
    ! parse the command line for options
    opts(1)=OPTION_S("version",.FALSE.,CHAR(21))
    opts(2)=OPTION_S("dry-run",.FALSE.,CHAR(22))
    opts(3)=OPTION_S("epsilon",.TRUE.,'e')
    opts(4)=OPTION_S("help",.FALSE.,'h')

    ! default accuracy (variable epsilon sits in the ode45 module)
    epsilon=1.d-6

    noptions=0;
    DO
       ch=getopt("he:",opts)
       SELECT CASE(ch)
       CASE(CHAR(0))
          EXIT
       CASE(CHAR(21))
          ! option version
          in%isversion=.TRUE.
       CASE(CHAR(22))
          ! option dry-run
          in%isdryrun=.TRUE.
       CASE('e')
          ! numerical accuracy
          READ(optarg,*) epsilon
          noptions=noptions+1
       CASE('h')
          ! option help
          in%ishelp=.TRUE.
       CASE('?')
          WRITE_DEBUG_INFO(100)
          in%ishelp=.TRUE.
          EXIT
       CASE DEFAULT
          WRITE (0,'("unhandled option ", a, " (this is a bug")') optopt
          WRITE_DEBUG_INFO(100)
          STOP 3
       END SELECT
       noptions=noptions+1
    END DO

    IF (in%isversion) THEN
       CALL printversion()
       ! abort parameter input
       RETURN
    END IF

    IF (in%ishelp) THEN
       CALL printhelp()
       ! abort parameter input
       RETURN
    END IF

    in%nPatch=0
    ! minimum number of dynamic variables for patches
    in%dPatch=STATE_VECTOR_DGF_PATCH
    in%nVolume=0
    ! minimum number of dynamic variables for strain volumes
    in%dVolume=STATE_VECTOR_DGF_VOLUME
    in%rectangularPatch%ns=0
    in%triangularPatch%ns=0
    in%triangularPatch%nVe=0
    in%strainVolume%ns=0

    IF (0.eq.rank) THEN
       PRINT 2000
       PRINT '("# UNICYCLE: unified cycle of earthquakes")'
       PRINT '("# quasi-dynamic earthquake simulation in viscoelastic media")'
       PRINT '("# numerical accuracy: epsilon=",ES9.2)', epsilon
!$     PRINT '("#     * parallel OpenMP implementation with ",I3.3,"/",I3.3," threads")', &
!$                  omp_get_max_threads(),omp_get_num_procs()
       PRINT 2000

       IF (noptions .LT. COMMAND_ARGUMENT_COUNT()) THEN
          ! read from input file
          iunit=25
          CALL GET_COMMAND_ARGUMENT(noptions+1,filename)
          OPEN (UNIT=iunit,FILE=filename,IOSTAT=iostatus)
       ELSE
          ! get input parameters from standard input
          iunit=5
       END IF

       PRINT '("# output directory")'
       CALL getdata(iunit,dataline)
       READ (dataline,'(a)') in%wdir
       PRINT '(2X,a)', TRIM(in%wdir)

       in%timeFilename=TRIM(in%wdir)//"/time.txt"

       ! test write permissions on output directory
       OPEN (UNIT=FPTIME,FILE=in%timeFilename,POSITION="APPEND",&
               IOSTAT=iostatus,FORM="FORMATTED")
       IF (iostatus>0) THEN
          WRITE_DEBUG_INFO(102)
          WRITE (STDERR,'("error: unable to access ",a)') TRIM(in%timefilename)
          STOP 1
       END IF
       CLOSE(FPTIME)
   
       PRINT '("# elastic moduli")'
       CALL getdata(iunit,dataline)
       READ  (dataline,*) in%lambda,in%mu
       PRINT '(2ES9.2E1)', in%lambda,in%mu

       IF (0 .GT. in%mu) THEN
          WRITE_DEBUG_INFO(-1)
          WRITE (STDERR,'(a)') TRIM(dataline)
          WRITE (STDERR,'("input error: shear modulus must be positive")')
          STOP -1
       END IF

       in%nu=in%lambda/2._8/(in%lambda+in%mu)
       IF (-1._8 .GT. in%nu) THEN
          WRITE (STDERR,'(a)') TRIM(dataline)
          WRITE (STDERR,'("input error: Poisson''s ratio must be greater than -1.")')
          STOP -1
       END IF
       IF (0.5_8 .LT. in%nu) THEN
          WRITE (STDERR,'(a)') TRIM(dataline)
          WRITE (STDERR,'("input error: Poisson''s ratio must be lower than 0.5.")')
          STOP -1
       END IF

       PRINT '("# time interval")'
       CALL getdata(iunit,dataline)
       READ  (dataline,*) in%interval
       PRINT '(ES20.12E2)', in%interval

       IF (in%interval .LE. 0._8) THEN
          WRITE (STDERR,'("**** error **** ")')
          WRITE (STDERR,'(a)') TRIM(dataline)
          WRITE (STDERR,'("simulation time must be positive. exiting.")')
          STOP 1
       END IF

       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       !       R E C T A N G U L A R   P A T C H E S
       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       PRINT '("# number of rectangular patches")'
       CALL getdata(iunit,dataline)
       READ  (dataline,*) in%rectangularPatch%ns
       PRINT '(I5)', in%rectangularPatch%ns
       IF (in%rectangularPatch%ns .GT. 0) THEN
          ALLOCATE(in%rectangularPatch%s(in%rectangularPatch%ns), &
                   in%rectangularPatch%x(3,in%rectangularPatch%ns), &
                   in%rectangularPatch%xc(3,in%rectangularPatch%ns), &
                   in%rectangularPatch%length(in%rectangularPatch%ns), &
                   in%rectangularPatch%width(in%rectangularPatch%ns), &
                   in%rectangularPatch%strike(in%rectangularPatch%ns), &
                   in%rectangularPatch%dip(in%rectangularPatch%ns), &
                   in%rectangularPatch%sv(3,in%rectangularPatch%ns), &
                   in%rectangularPatch%dv(3,in%rectangularPatch%ns), &
                   in%rectangularPatch%nv(3,in%rectangularPatch%ns),STAT=iostatus)
          IF (iostatus>0) STOP "could not allocate the rectangular patch list"
          PRINT 2000
          PRINT '("# n      Vpl       x1       x2       x3  length   width strike   dip   rake")'
          PRINT 2000
          DO k=1,in%rectangularPatch%ns
             CALL getdata(iunit,dataline)
             READ (dataline,*,IOSTAT=iostatus) i, &
                  in%rectangularPatch%s(k)%Vpl, &
                  in%rectangularPatch%x(1,k), &
                  in%rectangularPatch%x(2,k), &
                  in%rectangularPatch%x(3,k), &
                  in%rectangularPatch%length(k), &
                  in%rectangularPatch%width(k), &
                  in%rectangularPatch%strike(k), &
                  in%rectangularPatch%dip(k), &
                  in%rectangularPatch%s(k)%rake
             in%rectangularPatch%s(k)%opening=0
   
             PRINT '(I3.3,4ES9.2E1,2ES8.2E1,f7.1,f6.1,f7.1)',i, &
                  in%rectangularPatch%s(k)%Vpl, &
                  in%rectangularPatch%x(1,k), &
                  in%rectangularPatch%x(2,k), &
                  in%rectangularPatch%x(3,k), &
                  in%rectangularPatch%length(k), &
                  in%rectangularPatch%width(k), &
                  in%rectangularPatch%strike(k), &
                  in%rectangularPatch%dip(k), &
                  in%rectangularPatch%s(k)%rake
                

             ! convert to radians
             in%rectangularPatch%strike(k)=in%rectangularPatch%strike(k)*DEG2RAD     
             in%rectangularPatch%dip(k)=in%rectangularPatch%dip(k)*DEG2RAD     
             in%rectangularPatch%s(k)%rake=in%rectangularPatch%s(k)%rake*DEG2RAD     

             IF (i .ne. k) THEN
                WRITE (STDERR,'("invalid rectangular patch definition")')
                WRITE (STDERR,'(a)') TRIM(dataline)
                WRITE (STDERR,'("error in input file: unexpected index")')
                STOP 1
             END IF
             IF (MAX(in%rectangularPatch%length(k),in%rectangularPatch%width(k)) .LE. 0._8) THEN
                WRITE (STDERR,'(a)') TRIM(dataline)
                WRITE (STDERR,'("error in input file: patch length and width must be positive.")')
                STOP 1
             END IF
                
          END DO
   
          ! - - - - - - - - - - - - - - - - - - - - - - - - - -
          !        F R I C T I O N   P R O P E R T I E S
          ! - - - - - - - - - - - - - - - - - - - - - - - - - -
          PRINT 2000
          PRINT '("# number of frictional rectangular patches")'
          CALL getdata(iunit,dataline)
          READ  (dataline,*) dummy
          PRINT '(I5)', dummy
          IF (dummy .NE. in%rectangularPatch%ns) THEN
             WRITE_DEBUG_INFO(-1)
             WRITE(STDERR,'("input error: all rectangular patches require frictional properties")')
             STOP -1
          END IF
          PRINT '("# n      mu0      sig        a        b        L       Vo   G/(2Vs)")'
          PRINT 2000
          DO k=1,in%rectangularPatch%ns
             CALL getdata(iunit,dataline)
             READ (dataline,*,IOSTAT=iostatus) i, &
                   in%rectangularPatch%s(k)%mu0, &
                   in%rectangularPatch%s(k)%sig, &
                   in%rectangularPatch%s(k)%a, &
                   in%rectangularPatch%s(k)%b, &
                   in%rectangularPatch%s(k)%L, &
                   in%rectangularPatch%s(k)%Vo, &
                   in%rectangularPatch%s(k)%damping
   
             PRINT '(I3.3,7ES9.2E1)',i, &
                  in%rectangularPatch%s(k)%mu0, &
                  in%rectangularPatch%s(k)%sig, &
                  in%rectangularPatch%s(k)%a, &
                  in%rectangularPatch%s(k)%b, &
                  in%rectangularPatch%s(k)%L, &
                  in%rectangularPatch%s(k)%Vo, &
                  in%rectangularPatch%s(k)%damping
                
             IF (i .ne. k) THEN
                WRITE_DEBUG_INFO(200)
                WRITE (STDERR,'("invalid friction property definition for rectangular patch")')
                WRITE (STDERR,'(a)') TRIM(dataline)
                WRITE (STDERR,'("error in input file: unexpected index")')
                STOP 1
             END IF
          END DO
   
       END IF
          
       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       !        T R I A N G U L A R   P A T C H E S
       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       PRINT '("# number of triangular patches")'
       CALL getdata(iunit,dataline)
       READ  (dataline,*) in%triangularPatch%ns
       PRINT '(I5)', in%triangularPatch%ns
       IF (in%triangularPatch%ns .GT. 0) THEN
          ALLOCATE(in%triangularPatch%s(in%triangularPatch%ns), &
                   in%triangularPatch%i1(in%triangularPatch%ns), &
                   in%triangularPatch%i2(in%triangularPatch%ns), &
                   in%triangularPatch%i3(in%triangularPatch%ns), &
                   in%triangularPatch%sv(3,in%triangularPatch%ns), &
                   in%triangularPatch%dv(3,in%triangularPatch%ns), &
                   in%triangularPatch%nv(3,in%triangularPatch%ns), &
                   in%triangularPatch%xc(3,in%triangularPatch%ns),STAT=iostatus)
          IF (iostatus>0) STOP "could not allocate the triangular patch list"
          PRINT 2000
          PRINT '("# n      Vpl        i1        i2        i3   rake")'
          PRINT 2000
          DO k=1,in%triangularPatch%ns
             CALL getdata(iunit,dataline)
             READ (dataline,*,IOSTAT=iostatus) i, &
                  in%triangularPatch%s(k)%Vpl, &
                  in%triangularPatch%i1(k), &
                  in%triangularPatch%i2(k), &
                  in%triangularPatch%i3(k), &
                  in%triangularPatch%s(k)%rake
             in%triangularPatch%s(k)%opening=0
   
             PRINT '(I3.3,ES9.2E1,3I10,f7.1)',i, &
                  in%triangularPatch%s(k)%Vpl, &
                  in%triangularPatch%i1(k), &
                  in%triangularPatch%i2(k), &
                  in%triangularPatch%i3(k), &
                  in%triangularPatch%s(k)%rake
                
             ! convert to radians
             in%triangularPatch%s(k)%rake=in%triangularPatch%s(k)%rake*DEG2RAD     

             ! check range of vertex index
             IF (in%triangularPatch%i1(k) .GT. maxVertexIndex) THEN
                maxVertexIndex=in%triangularPatch%i1(k)
             END IF
             IF (in%triangularPatch%i2(k) .GT. maxVertexIndex) THEN
                maxVertexIndex=in%triangularPatch%i2(k)
             END IF
             IF (in%triangularPatch%i3(k) .GT. maxVertexIndex) THEN
                maxVertexIndex=in%triangularPatch%i3(k)
             END IF
             IF ((0 .GT. in%triangularPatch%i1(k)) .OR. &
                 (0 .GT. in%triangularPatch%i2(k)) .OR. &
                 (0 .GT. in%triangularPatch%i3(k))) THEN
                WRITE (STDERR,'(a)') TRIM(dataline)
                WRITE (STDERR,'("error in input file: negative index")')
                STOP 1
             END IF

             IF (i .ne. k) THEN
                WRITE (STDERR,'("error in input file: unexpected index")')
                WRITE (STDERR,'(a)') TRIM(dataline)
                WRITE (STDERR,'("invalid triangular patch definition ")')
                STOP 1
             END IF
          END DO
                
          PRINT '("# number of vertices")'
          CALL getdata(iunit,dataline)
          READ  (dataline,*) in%triangularPatch%nVe
          PRINT '(I5)', in%triangularPatch%nVe
          IF (maxVertexIndex .GT. in%triangularPatch%nVe) THEN
             WRITE (STDERR,'("error in input file: not enough vertices")')
             STOP 1
          END IF
          IF (in%triangularPatch%nVe .GT. 0) THEN
             ALLOCATE(in%triangularPatch%v(3,in%triangularPatch%nVe),STAT=iostatus)
             IF (iostatus>0) STOP "could not allocate the list of vertices"
             PRINT 2000
             PRINT '("# n      x1       x2       x3")'
             PRINT 2000
             DO k=1,in%triangularPatch%nVe
                CALL getdata(iunit,dataline)
                READ (dataline,*,IOSTAT=iostatus) i, &
                      in%triangularPatch%v(1,k), &
                      in%triangularPatch%v(2,k), &
                      in%triangularPatch%v(3,k)
   
                PRINT '(I3.3,3ES9.2E1)',i, &
                     in%triangularPatch%v(1,k), &
                     in%triangularPatch%v(2,k), &
                     in%triangularPatch%v(3,k)
                
                IF (i .ne. k) THEN
                   WRITE_DEBUG_INFO(200)
                   WRITE (STDERR,'("invalid vertex definition ")')
                   WRITE (STDERR,'(a)') TRIM(dataline)
                   WRITE (STDERR,'("error in input file: unexpected index")')
                   STOP 1
                END IF
             END DO
                
          END IF
   
          ! - - - - - - - - - - - - - - - - - - - - - - - - - -
          !        F R I C T I O N   P R O P E R T I E S
          ! - - - - - - - - - - - - - - - - - - - - - - - - - -
          PRINT 2000
          PRINT '("# number of frictional triangular patches")'
          CALL getdata(iunit,dataline)
          READ  (dataline,*) dummy
          PRINT '(I5)', dummy
          IF (dummy .NE. in%triangularPatch%ns) THEN
             WRITE_DEBUG_INFO(-1)
             WRITE(STDERR,'("input error: all triangular patches require frictional properties")')
             STOP -1
          END IF

          PRINT '("# n      mu0      sig        a        b        L       Vo    2G/Vs")'
          PRINT 2000
          DO k=1,in%triangularPatch%ns
             CALL getdata(iunit,dataline)
             READ (dataline,*,IOSTAT=iostatus) i, &
                   in%triangularPatch%s(k)%mu0, &
                   in%triangularPatch%s(k)%sig, &
                   in%triangularPatch%s(k)%a, &
                   in%triangularPatch%s(k)%b, &
                   in%triangularPatch%s(k)%L, &
                   in%triangularPatch%s(k)%Vo, &
                   in%triangularPatch%s(k)%damping
   
             PRINT '(I3.3,7ES9.2E1)',i, &
                  in%triangularPatch%s(k)%mu0, &
                  in%triangularPatch%s(k)%sig, &
                  in%triangularPatch%s(k)%a, &
                  in%triangularPatch%s(k)%b, &
                  in%triangularPatch%s(k)%L, &
                  in%triangularPatch%s(k)%Vo, &
                  in%triangularPatch%s(k)%damping
                
             IF (i .ne. k) THEN
                WRITE_DEBUG_INFO(200)
                WRITE (STDERR,'("invalid friction property definition for triangular patch")')
                WRITE (STDERR,'(a)') TRIM(dataline)
                WRITE (STDERR,'("error in input file: unexpected index")')
                STOP 1
             END IF
          END DO
       END IF
          
   
       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       !       S T R A I N   V O L U M E S
       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       PRINT '("# number of cuboid strain volumes")'
       CALL getdata(iunit,dataline)
       READ  (dataline,*) in%strainVolume%ns
       PRINT '(I5)', in%strainVolume%ns
       IF (in%strainVolume%ns .GT. 0) THEN
          ALLOCATE(in%strainVolume%s(in%strainVolume%ns), &
                   in%strainVolume%x(3,in%strainVolume%ns), &
                   in%strainVolume%xc(3,in%strainVolume%ns), &
                   in%strainVolume%length(in%strainVolume%ns), &
                   in%strainVolume%width(in%strainVolume%ns), &
                   in%strainVolume%thickness(in%strainVolume%ns), &
                   in%strainVolume%strike(in%strainVolume%ns), &
                   in%strainVolume%dip(in%strainVolume%ns), &
                   in%strainVolume%sv(3,in%strainVolume%ns), &
                   in%strainVolume%dv(3,in%strainVolume%ns), &
                   in%strainVolume%nv(3,in%strainVolume%ns),STAT=iostatus)
          IF (iostatus>0) STOP "could not allocate the strain volume list"
          PRINT 2000
          PRINT '("# n       e11       e12       e13       e22       e23       e33       ", &
                  "x1       x2       x3  length   width thickness strike")'
          PRINT 2000
          DO k=1,in%strainVolume%ns
             CALL getdata(iunit,dataline)
             READ (dataline,*,IOSTAT=iostatus) i, &
                  in%strainVolume%s(k)%e11, &
                  in%strainVolume%s(k)%e12, &
                  in%strainVolume%s(k)%e13, &
                  in%strainVolume%s(k)%e22, &
                  in%strainVolume%s(k)%e23, &
                  in%strainVolume%s(k)%e33, &
                  in%strainVolume%x(1,k), &
                  in%strainVolume%x(2,k), &
                  in%strainVolume%x(3,k), &
                  in%strainVolume%length(k), &
                  in%strainVolume%width(k), &
                  in%strainVolume%thickness(k), &
                  in%strainVolume%strike(k), &
                  in%strainVolume%dip(k)
   
             PRINT '(I3.3,6ES10.2E2,3ES9.2E1,3ES8.2E1,f7.1,f7.1)',i, &
                  in%strainVolume%s(k)%e11, &
                  in%strainVolume%s(k)%e12, &
                  in%strainVolume%s(k)%e13, &
                  in%strainVolume%s(k)%e22, &
                  in%strainVolume%s(k)%e23, &
                  in%strainVolume%s(k)%e33, &
                  in%strainVolume%x(1,k), &
                  in%strainVolume%x(2,k), &
                  in%strainVolume%x(3,k), &
                  in%strainVolume%length(k), &
                  in%strainVolume%width(k), &
                  in%strainVolume%thickness(k), &
                  in%strainVolume%strike(k), &
                  in%strainVolume%dip(k)
                
             IF (90. .NE. in%strainVolume%dip(k)) THEN
                WRITE_DEBUG_INFO(-1)
                WRITE (STDERR,'(a)') TRIM(dataline)
                WRITE (STDERR,'("input error: strain volume must be vertical (dip=90).")')
                STOP -1
             END IF

             ! convert to radians
             in%strainVolume%strike(k)=in%strainVolume%strike(k)*DEG2RAD     
             in%strainVolume%dip(k)=in%strainVolume%dip(k)*DEG2RAD     

             IF (0 .GT. in%strainVolume%x(3,k)) THEN
                WRITE_DEBUG_INFO(-1)
                WRITE (STDERR,'(a)') TRIM(dataline)
                WRITE (STDERR,'("input error: depth must be positive")')
                STOP -1
             END IF

             IF (i .ne. k) THEN
                WRITE (STDERR,'("invalid strain volume definition")')
                WRITE (STDERR,'(a)') TRIM(dataline)
                WRITE (STDERR,'("input error: unexpected index")')
                STOP 1
             END IF
             IF (MAX(in%strainVolume%length(k),in%strainVolume%width(k),in%strainVolume%thickness(k)) .LE. 0._8) THEN
                WRITE (STDERR,'(a)') TRIM(dataline)
                WRITE (STDERR,'("input error: strain volume dimension must be positive.")')
                STOP 1
             END IF
                
          END DO
   
          ! - - - - - - - - - - - - - - - - - - - - - - - - - -
          !   L I N E A R   K E L V I N   P R O P E R T I E S
          ! - - - - - - - - - - - - - - - - - - - - - - - - - -
          PRINT 2000
          PRINT '("# number of Kelvin strain volumes")'
          CALL getdata(iunit,dataline)
          READ  (dataline,*) in%strainVolume%nKelvin
          PRINT '(I5)', in%strainVolume%nKelvin
          IF (0 .NE. in%strainVolume%nKelvin) THEN
             IF (in%strainVolume%ns .NE. in%strainVolume%nKelvin) THEN
                WRITE_DEBUG_INFO(-1)
                WRITE (STDERR,'(a)') TRIM(dataline)
                WRITE(STDERR,'("input error: Kelvin properties ", &
                               "are for none or all strain volumes")')
                STOP -1
             END IF

             ! plan six dynamic variables per Kelvin strain volume
             in%dVolume=in%dVolume+6

             PRINT '("# n     Gk   gammadot0k    do   m    Q    R")'
             PRINT 2000
             DO k=1,in%strainVolume%ns
                CALL getdata(iunit,dataline)
                READ (dataline,*,IOSTAT=iostatus) i, &
                      in%strainVolume%s(k)%Gk, &
                      in%strainVolume%s(k)%gammadot0k, &
                      in%strainVolume%s(k)%dok, &
                      in%strainVolume%s(k)%mk, &
                      in%strainVolume%s(k)%Qk, &
                      in%strainVolume%s(k)%Rk
   
                PRINT '(I3.3,6ES9.2E1)',i, &
                      in%strainVolume%s(k)%Gk, &
                      in%strainVolume%s(k)%gammadot0k, &
                      in%strainVolume%s(k)%dok, &
                      in%strainVolume%s(k)%mk, &
                      in%strainVolume%s(k)%Qk, &
                      in%strainVolume%s(k)%Rk
                
                IF (i .ne. k) THEN
                   WRITE_DEBUG_INFO(200)
                   WRITE (STDERR,'("invalid property definition for strain volume")')
                   WRITE (STDERR,'(a)') TRIM(dataline)
                   WRITE (STDERR,'("error in input file: unexpected index")')
                   STOP 1
                END IF
             END DO

          END IF
   
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          !  L I N E A R   M A X W E L L   P R O P E R T I E S
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          PRINT 2000
          PRINT '("# number of Maxwell strain volumes")'
          CALL getdata(iunit,dataline)
          READ  (dataline,*) in%strainVolume%nMaxwell
          PRINT '(I5)', in%strainVolume%nMaxwell
          IF (0 .NE. in%strainVolume%nMaxwell) THEN
             IF (in%strainVolume%ns .NE. in%strainVolume%nMaxwell) THEN
                WRITE_DEBUG_INFO(-1)
                WRITE (STDERR,'(a)') TRIM(dataline)
                WRITE(STDERR,'("input error: Maxwell properties ", &
                               "are for none or all strain volumes")')
                STOP -1
             END IF

             PRINT '("# n     Gm   gammadot0m    do   m    Q    R")'
             PRINT 2000
             DO k=1,in%strainVolume%ns
                CALL getdata(iunit,dataline)
                READ (dataline,*,IOSTAT=iostatus) i, &
                      in%strainVolume%s(k)%Gm, &
                      in%strainVolume%s(k)%gammadot0m, &
                      in%strainVolume%s(k)%dom, &
                      in%strainVolume%s(k)%mm, &
                      in%strainVolume%s(k)%Qm, &
                      in%strainVolume%s(k)%Rm
   
                PRINT '(I3.3,6ES9.2E1)',i, &
                      in%strainVolume%s(k)%Gm, &
                      in%strainVolume%s(k)%gammadot0m, &
                      in%strainVolume%s(k)%dom, &
                      in%strainVolume%s(k)%mm, &
                      in%strainVolume%s(k)%Qm, &
                      in%strainVolume%s(k)%Rm
                
                IF (i .ne. k) THEN
                   WRITE_DEBUG_INFO(200)
                   WRITE (STDERR,'("invalid property definition for strain volume")')
                   WRITE (STDERR,'(a)') TRIM(dataline)
                   WRITE (STDERR,'("error in input file: unexpected index")')
                   STOP 1
                END IF
             END DO
          END IF
   
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          !   N O N L I N E A R   K E L V I N   P R O P E R T I E S
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          PRINT 2000
          PRINT '("# number of nonlinear Kelvin strain volumes")'
          CALL getdata(iunit,dataline)
          READ  (dataline,*) in%strainVolume%nNonlinearKelvin
          PRINT '(I5)', in%strainVolume%nNonlinearKelvin
          IF (0 .NE. in%strainVolume%nNonlinearKelvin) THEN
             IF (in%strainVolume%ns .NE. in%strainVolume%nNonlinearKelvin) THEN
                WRITE_DEBUG_INFO(-1)
                WRITE (STDERR,'(a)') TRIM(dataline)
                WRITE(STDERR,'("input error: nonlinear Kelvin properties ", &
                               "are for none or all strain volumes")')
                STOP -1
             END IF

             ! plan six dynamic variables per nonlinear Kelvin strain volume
             in%dVolume=in%dVolume+6

             PRINT '("# n     Gk   gammadot0k    n    Q    R")'
             PRINT 2000
             DO k=1,in%strainVolume%ns
                CALL getdata(iunit,dataline)
                READ (dataline,*,IOSTAT=iostatus) i, &
                      in%strainVolume%s(k)%nGk, &
                      in%strainVolume%s(k)%ngammadot0k, &
                      in%strainVolume%s(k)%npowerk, &
                      in%strainVolume%s(k)%nQk, &
                      in%strainVolume%s(k)%nRk
   
                PRINT '(I3.3,5ES9.2E1)',i, &
                      in%strainVolume%s(k)%nGk, &
                      in%strainVolume%s(k)%ngammadot0k, &
                      in%strainVolume%s(k)%npowerk, &
                      in%strainVolume%s(k)%nQk, &
                      in%strainVolume%s(k)%nRk
                
                IF (i .ne. k) THEN
                   WRITE_DEBUG_INFO(200)
                   WRITE (STDERR,'("invalid property definition for strain volume")')
                   WRITE (STDERR,'(a)') TRIM(dataline)
                   WRITE (STDERR,'("error in input file: unexpected index")')
                   STOP 1
                END IF
             END DO

          END IF
   
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          !   N O N L I N E A R   M A X W E L L   P R O P E R T I E S
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          PRINT 2000
          PRINT '("# number of nonlinear Maxwell strain volumes")'
          CALL getdata(iunit,dataline)
          READ  (dataline,*) in%strainVolume%nNonlinearMaxwell
          PRINT '(I5)', in%strainVolume%nNonlinearMaxwell
          IF (0 .NE. in%strainVolume%nNonlinearMaxwell) THEN
             IF (in%strainVolume%ns .NE. in%strainVolume%nNonlinearMaxwell) THEN
                WRITE_DEBUG_INFO(-1)
                WRITE (STDERR,'(a)') TRIM(dataline)
                WRITE(STDERR,'("input error: nonlinear Maxwell properties ", &
                               "are for none or all strain volumes")')
                STOP -1
             END IF

             PRINT '("# n       Gm  gammadot0m        n        Q        R")'
             PRINT 2000
             DO k=1,in%strainVolume%ns
                CALL getdata(iunit,dataline)
                READ (dataline,*,IOSTAT=iostatus) i, &
                      in%strainVolume%s(k)%nGm, &
                      in%strainVolume%s(k)%ngammadot0m, &
                      in%strainVolume%s(k)%npowerm, &
                      in%strainVolume%s(k)%nQm, &
                      in%strainVolume%s(k)%nRm
   
                PRINT '(I3.3,ES9.2E1,ES11.4E1,3ES9.2E1)',i, &
                      in%strainVolume%s(k)%nGm, &
                      in%strainVolume%s(k)%ngammadot0m, &
                      in%strainVolume%s(k)%npowerm, &
                      in%strainVolume%s(k)%nQm, &
                      in%strainVolume%s(k)%nRm
                
                IF (0 .GE. in%strainVolume%s(k)%nGm) THEN
                   WRITE_DEBUG_INFO(200)
                   WRITE (STDERR,'("invalid property definition for strain volume")')
                   WRITE (STDERR,'(a)') TRIM(dataline)
                   WRITE (STDERR,'("error in input file: sig0 must be positive.")')
                   STOP 1
                END IF

                IF (0 .GE. in%strainVolume%s(k)%nRm) THEN
                   WRITE_DEBUG_INFO(200)
                   WRITE (STDERR,'("invalid property definition for strain volume")')
                   WRITE (STDERR,'(a)') TRIM(dataline)
                   WRITE (STDERR,'("error in input file: R must be positive.")')
                   STOP 1
                END IF

                IF (i .ne. k) THEN
                   WRITE_DEBUG_INFO(200)
                   WRITE (STDERR,'("invalid property definition for strain volume")')
                   WRITE (STDERR,'(a)') TRIM(dataline)
                   WRITE (STDERR,'("error in input file: unexpected index")')
                   STOP 1
                END IF
             END DO

          END IF
   
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          !          T H E R M A L   P R O P E R T I E S
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          PRINT 2000
          PRINT '("# number of thermal strain volumes")'
          CALL getdata(iunit,dataline)
          READ  (dataline,*) dummy
          PRINT '(I5)', dummy
          IF (0 .NE. dummy) THEN
             IF (in%strainVolume%ns .NE. dummy) THEN
                WRITE_DEBUG_INFO(-1)
                WRITE (STDERR,'(a)') TRIM(dataline)
                WRITE(STDERR,'("input error: thermal properties ", &
                               "are for none or all strain volumes")')
                STOP -1
             END IF

             PRINT '("# n     rhoc  temperature")'
             PRINT 2000
             DO k=1,in%strainVolume%ns
                CALL getdata(iunit,dataline)
                READ (dataline,*,IOSTAT=ierr) i, &
                      in%strainVolume%s(k)%rhoc, &
                      in%strainVolume%s(k)%To
   
                PRINT '(I3.3,5ES11.4E1)',i, &
                      in%strainVolume%s(k)%rhoc, &
                      in%strainVolume%s(k)%To
                
                IF (0 .GE. in%strainVolume%s(k)%To) THEN
                   WRITE_DEBUG_INFO(200)
                   WRITE (STDERR,'("invalid property definition for strain volume")')
                   WRITE (STDERR,'(a)') TRIM(dataline)
                   WRITE (STDERR,'("error in input file: To must be positive.")')
                   STOP 1
                END IF

                IF (i .NE. k) THEN
                   WRITE_DEBUG_INFO(200)
                   WRITE (STDERR,'("invalid property definition for strain volume")')
                   WRITE (STDERR,'(a)') TRIM(dataline)
                   WRITE (STDERR,'("error in input file: unexpected index")')
                   STOP 1
                END IF
             END DO

          END IF
   
       END IF
          
       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       !       O B S E R V A T I O N   P A T C H E S
       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       PRINT 2000
       PRINT '("# number of observation patches")'
       CALL getdata(iunit,dataline)
       READ  (dataline,*) nObservationPatch
       PRINT '(I5)', nObservationPatch
       IF (0 .LT. nObservationPatch) THEN
          ALLOCATE(observationPatch(nObservationPatch,2),STAT=ierr)
          IF (ierr>0) STOP "could not allocate the observation patches"
          PRINT 2000
          PRINT '("# n      i")'
          PRINT 2000
          DO k=1,in%nObservationPatch
             CALL getdata(iunit,dataline)
             READ (dataline,*,IOSTAT=ierr) i,observationPatch(k)
             PRINT '(I3.3,X,I6)',i,observationPatch(k)
             IF (i .NE. k) THEN
                WRITE_DEBUG_INFO(200)
                WRITE (STDERR,'(a)') TRIM(dataline)
                WRITE (STDERR,'("error in input file: unexpected index")')
                STOP 1
             END IF
          END DO
       END IF

       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       !       O B S E R V A T I O N   V O L U M E S
       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       PRINT 2000
       PRINT '("# number of observation volumes")'
       CALL getdata(iunit,dataline)
       READ  (dataline,*) nObservationVolume
       PRINT '(I5)', nObservationVolume
       IF (0 .LT. nObservationVolume) THEN
          ALLOCATE(observationVolume(nObservationVolume),STAT=ierr)
          IF (ierr>0) STOP "could not allocate the observation volumes0"
          PRINT 2000
          PRINT '("# n      i")'
          PRINT 2000
          DO k=1,nObservationVolume
             CALL getdata(iunit,dataline)
             READ (dataline,*,IOSTAT=ierr) i,observationVolume(k)
             PRINT '(I3.3,X,I6)',i,observationVolume(k)
             IF (i .NE. k) THEN
                WRITE_DEBUG_INFO(200)
                WRITE (STDERR,'(a)') TRIM(dataline)
                WRITE (STDERR,'("error in input file: unexpected index")')
                STOP 1
             END IF
          END DO
       END IF

       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       !        O B S E R V A T I O N   P O I N T S
       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       PRINT 2000
       PRINT '("# number of observation points")'
       CALL getdata(iunit,dataline)
       READ  (dataline,*) in%nObservationPoint
       PRINT '(I5)', in%nObservationPoint
       IF (0 .LT. in%nObservationPoint) THEN
          ALLOCATE(in%observationPoint(in%nObservationPoint),STAT=ierr)
          IF (ierr>0) STOP "could not allocate the observation points"
          PRINT 2000
          PRINT '("# n name       x1       x2       x3")'
          PRINT 2000
          DO k=1,in%nObservationPoint
             CALL getdata(iunit,dataline)

             READ (dataline,*,IOSTAT=ierr) i, &
                     in%observationPoint(k)%name, &
                     in%observationPoint(k)%x(1), &
                     in%observationPoint(k)%x(2), &
                     in%observationPoint(k)%x(3)

             PRINT '(I3.3,X,a4,3ES9.2E1)',i, &
                     in%observationPoint(k)%name, &
                     in%observationPoint(k)%x(1), &
                     in%observationPoint(k)%x(2), &
                     in%observationPoint(k)%x(3)

             IF (i .NE. k) THEN
                WRITE_DEBUG_INFO(200)
                WRITE (STDERR,'(a)') TRIM(dataline)
                WRITE (STDERR,'("error in input file: unexpected index")')
                STOP 1
             END IF
          END DO
       END IF

       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       !            P E R T U R B A T I O N S
       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       PRINT '("# number of perturbations")'
       CALL getdata(iunit,dataline)
       READ (dataline,*) in%ne
       PRINT '(I5)', in%ne
       IF (in%ne .GT. 0) ALLOCATE(in%event(in%ne),STAT=iostatus)
       IF (iostatus>0) STOP "could not allocate the event list"
       
       DO i=1,in%ne
          IF (1 .NE. i) THEN
             PRINT '("# time of next perturbation")'
             CALL getdata(iunit,dataline)
             READ (dataline,*) in%event(i)%time
             in%event(i)%i=i-1
             PRINT '(ES9.2E1)', in%event(i)%time
   
             IF (in%event(i)%time .LE. in%event(i-1)%time) THEN
                WRITE_DEBUG_INFO(200)
                WRITE (STDERR,'(a)') TRIM(dataline)
                WRITE (STDERR,'(a,a)') "input file error. ", &
                     "timing of perturbations must increase, quiting."
                STOP 1
             END IF
          ELSE
             in%event(1)%time=0._8
             in%event(1)%i=0
          END IF
   
       END DO
   
       ! test the presence of dislocations for coseismic calculation
       IF ((in%rectangularPatch%ns .EQ. 0) .AND. &
           (in%triangularPatch%ns .EQ. 0) .AND. &
           (in%strainVolume%ns .EQ. 0) .OR. &
           (in%interval .LE. 0._8)) THEN
   
          WRITE_DEBUG_INFO(300)
          WRITE (STDERR,'("nothing to do. exiting.")')
          STOP 1
       END IF
   
       PRINT 2000
       ! flush standard output
       CALL FLUSH(6)      

       in%nPatch=in%rectangularPatch%ns+in%triangularPatch%ns

       ! combine observation patches and observation volumes into observation states
       in%nObservationState=nObservationPatch+nObservationVolume
       ALLOCATE(in%observationState(in%nObservationState,2))
       j=1
       DO i=1,nObservationPatch
          in%observationState(j,1)=observationPatch(i)
          j=j+1
       END DO

       DO i=1,nObservationVolume
          in%observationState(j,1)=observationVolume(i)+in%nPatch
          j=j+1
       END DO
       IF (ALLOCATED(observationPatch)) DEALLOCATE(observationPatch)
       IF (ALLOCATED(observationVolume)) DEALLOCATE(observationVolume)

       position=0
       CALL MPI_PACK(in%interval,           1,MPI_REAL8,  packed,psize,position,MPI_COMM_WORLD,ierr)
       CALL MPI_PACK(in%lambda,             1,MPI_REAL8,  packed,psize,position,MPI_COMM_WORLD,ierr)
       CALL MPI_PACK(in%mu,                 1,MPI_REAL8,  packed,psize,position,MPI_COMM_WORLD,ierr)
       CALL MPI_PACK(in%nu,                 1,MPI_REAL8,  packed,psize,position,MPI_COMM_WORLD,ierr)
       CALL MPI_PACK(in%rectangularPatch%ns,1,MPI_INTEGER,packed,psize,position,MPI_COMM_WORLD,ierr)
       CALL MPI_PACK(in%triangularPatch%ns, 1,MPI_INTEGER,packed,psize,position,MPI_COMM_WORLD,ierr)
       CALL MPI_PACK(in%triangularPatch%nVe,1,MPI_INTEGER,packed,psize,position,MPI_COMM_WORLD,ierr)
       CALL MPI_PACK(in%strainVolume%ns,    1,MPI_INTEGER,packed,psize,position,MPI_COMM_WORLD,ierr)
       CALL MPI_PACK(in%ne,                 1,MPI_INTEGER,packed,psize,position,MPI_COMM_WORLD,ierr)
       CALL MPI_PACK(in%nObservationState,  1,MPI_INTEGER,packed,psize,position,MPI_COMM_WORLD,ierr)
       CALL MPI_PACK(in%nObservationPoint,  1,MPI_INTEGER,packed,psize,position,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)

       position=0
       CALL MPI_PACK(in%wdir,256,MPI_CHARACTER,packed,psize,position,MPI_COMM_WORLD,ierr)
       CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)

       ! send the rectangular patches (geometry and friction properties) 
       DO i=1,in%rectangularPatch%ns
          position=0
          CALL MPI_PACK(in%rectangularPatch%s(i)%Vpl,    1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%rectangularPatch%x(1,i),      1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%rectangularPatch%x(2,i),      1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%rectangularPatch%x(3,i),      1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%rectangularPatch%length(i),   1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%rectangularPatch%width(i),    1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%rectangularPatch%strike(i),   1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%rectangularPatch%dip(i),      1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%rectangularPatch%s(i)%rake,   1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%rectangularPatch%s(i)%opening,1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%rectangularPatch%s(i)%mu0,    1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%rectangularPatch%s(i)%sig,    1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%rectangularPatch%s(i)%a,      1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%rectangularPatch%s(i)%b,      1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%rectangularPatch%s(i)%L,      1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%rectangularPatch%s(i)%Vo,     1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%rectangularPatch%s(i)%damping,1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
       END DO

       ! send the triangular patches (geometry and friction properties) 
       DO i=1,in%triangularPatch%ns
          position=0
          CALL MPI_PACK(in%triangularPatch%s(i)%Vpl,    1,MPI_REAL8,  packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%triangularPatch%i1(i),     1,MPI_INTEGER,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%triangularPatch%i2(i),     1,MPI_INTEGER,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%triangularPatch%i3(i),     1,MPI_INTEGER,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%triangularPatch%s(i)%rake,   1,MPI_REAL8,  packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%triangularPatch%s(i)%opening,1,MPI_REAL8,  packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%triangularPatch%s(i)%mu0,    1,MPI_REAL8,  packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%triangularPatch%s(i)%sig,    1,MPI_REAL8,  packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%triangularPatch%s(i)%a,      1,MPI_REAL8,  packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%triangularPatch%s(i)%b,      1,MPI_REAL8,  packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%triangularPatch%s(i)%L,      1,MPI_REAL8,  packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%triangularPatch%s(i)%Vo,     1,MPI_REAL8,  packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%triangularPatch%s(i)%damping,1,MPI_REAL8,  packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
       END DO

       ! send the triangle vertices
       IF (0 .NE. in%triangularPatch%ns) THEN
          DO i=1,in%triangularPatch%nVe
             position=0
             CALL MPI_PACK(in%triangularPatch%v(1,i),1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
             CALL MPI_PACK(in%triangularPatch%v(2,i),1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
             CALL MPI_PACK(in%triangularPatch%v(3,i),1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
             CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
          END DO
       END IF

       ! send the strain volumes
       DO i=1,in%strainVolume%ns
          position=0
          CALL MPI_PACK(in%strainVolume%s(i)%e11,        1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%e12,        1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%e13,        1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%e22,        1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%e23,        1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%e33,        1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%x(1,i),          1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%x(2,i),          1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%x(3,i),          1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%length(i),       1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%width(i),        1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%thickness(i),    1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%strike(i),       1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%dip(i),          1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%Gk,         1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%gammadot0k, 1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%dok,        1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%mk,         1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%Qk,         1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%Rk,         1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%Gm,         1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%gammadot0m, 1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%dom,        1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%mm,         1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%Qm,         1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%Rm,         1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%nGk,        1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%ngammadot0k,1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%npowerk,       1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%nQk,        1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%nRk,        1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%nGm,        1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%ngammadot0m,1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%npowerm,       1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%nQm,        1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%nRm,        1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%rhoc,       1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%strainVolume%s(i)%To,         1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
       END DO

       ! send the observation state
       DO i=1,in%nObservationState
          position=0
          CALL MPI_PACK(in%observationState(i,1),1,MPI_INTEGER,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
       END DO

       ! send the observation points
       DO i=1,in%nObservationPoint
          position=0
          CALL MPI_PACK(in%observationPoint(i)%name,10,MPI_CHARACTER,packed,psize,position,MPI_COMM_WORLD,ierr)
          DO k=1,3
             CALL MPI_PACK(in%observationPoint(i)%x(k),1,MPI_REAL8,packed,psize,position,MPI_COMM_WORLD,ierr)
          END DO
          CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
       END DO

       ! send the perturbation events
       DO k=1,in%ne
          CALL MPI_PACK(in%event(k)%time,1,MPI_REAL8,  packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_PACK(in%event(k)%i,   1,MPI_INTEGER,packed,psize,position,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
       END DO

    ELSE ! IF 0 .NE. rank

       !------------------------------------------------------------------
       ! S L A V E S
       !------------------------------------------------------------------

       position=0
       CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
       CALL MPI_UNPACK(packed,psize,position,in%interval,           1,MPI_REAL8,  MPI_COMM_WORLD,ierr)
       CALL MPI_UNPACK(packed,psize,position,in%lambda,             1,MPI_REAL8,  MPI_COMM_WORLD,ierr)
       CALL MPI_UNPACK(packed,psize,position,in%mu,                 1,MPI_REAL8,  MPI_COMM_WORLD,ierr)
       CALL MPI_UNPACK(packed,psize,position,in%nu,                 1,MPI_REAL8,  MPI_COMM_WORLD,ierr)
       CALL MPI_UNPACK(packed,psize,position,in%rectangularPatch%ns,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
       CALL MPI_UNPACK(packed,psize,position,in%triangularPatch%ns, 1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
       CALL MPI_UNPACK(packed,psize,position,in%triangularPatch%nVe,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
       CALL MPI_UNPACK(packed,psize,position,in%strainVolume%ns,    1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
       CALL MPI_UNPACK(packed,psize,position,in%ne,                 1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
       CALL MPI_UNPACK(packed,psize,position,in%nObservationState,  1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
       CALL MPI_UNPACK(packed,psize,position,in%nObservationPoint,  1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

       position=0
       CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
       CALL MPI_UNPACK(packed,psize,position,in%wdir,256,MPI_CHARACTER,MPI_COMM_WORLD,ierr)

       IF (0 .LT. in%rectangularPatch%ns) &
                    ALLOCATE(in%rectangularPatch%s(in%rectangularPatch%ns), &
                             in%rectangularPatch%x(3,in%rectangularPatch%ns), &
                             in%rectangularPatch%xc(3,in%rectangularPatch%ns), &
                             in%rectangularPatch%length(in%rectangularPatch%ns), &
                             in%rectangularPatch%width(in%rectangularPatch%ns), &
                             in%rectangularPatch%strike(in%rectangularPatch%ns), &
                             in%rectangularPatch%dip(in%rectangularPatch%ns), &
                             in%rectangularPatch%sv(3,in%rectangularPatch%ns), &
                             in%rectangularPatch%dv(3,in%rectangularPatch%ns), &
                             in%rectangularPatch%nv(3,in%rectangularPatch%ns), &
                             STAT=iostatus)
       IF (iostatus>0) STOP "slave could not allocate memory"

       IF (0 .LT. in%triangularPatch%ns) &
                    ALLOCATE(in%triangularPatch%s(in%triangularPatch%ns), &
                             in%triangularPatch%xc(3,in%triangularPatch%ns), &
                             in%triangularPatch%i1(in%triangularPatch%ns), &
                             in%triangularPatch%i2(in%triangularPatch%ns), &
                             in%triangularPatch%i3(in%triangularPatch%ns), &
                             in%triangularPatch%sv(3,in%triangularPatch%ns), &
                             in%triangularPatch%dv(3,in%triangularPatch%ns), &
                             in%triangularPatch%nv(3,in%triangularPatch%ns),STAT=iostatus)
       IF (iostatus>0) STOP "slave could not allocate memory"

       IF (0 .LT. in%triangularPatch%nVe) ALLOCATE(in%triangularPatch%v(3,in%triangularPatch%nVe), STAT=iostatus)
       IF (iostatus>0) STOP "slave could not allocate memory"

       IF (0 .LT. in%strainVolume%ns) &
                    ALLOCATE(in%strainVolume%s(in%strainVolume%ns), &
                             in%strainVolume%x(3,in%strainVolume%ns), &
                             in%strainVolume%xc(3,in%strainVolume%ns), &
                             in%strainVolume%length(in%strainVolume%ns), &
                             in%strainVolume%width(in%strainVolume%ns), &
                             in%strainVolume%thickness(in%strainVolume%ns), &
                             in%strainVolume%strike(in%strainVolume%ns), &
                             in%strainVolume%dip(in%strainVolume%ns), &
                             in%strainVolume%sv(3,in%strainVolume%ns), &
                             in%strainVolume%dv(3,in%strainVolume%ns), &
                             in%strainVolume%nv(3,in%strainVolume%ns),STAT=iostatus)
       IF (iostatus>0) STOP "slave could not allocate memory"

       IF (0 .LT. in%ne) ALLOCATE(in%event(in%ne),STAT=iostatus)
       IF (iostatus>0) STOP "slave could not allocate memory"

       DO i=1,in%rectangularPatch%ns
          position=0
          CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%rectangularPatch%s(i)%Vpl,    1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%rectangularPatch%x(1,i),      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%rectangularPatch%x(2,i),      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%rectangularPatch%x(3,i),      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%rectangularPatch%length(i),   1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%rectangularPatch%width(i),    1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%rectangularPatch%strike(i),   1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%rectangularPatch%dip(i),      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%rectangularPatch%s(i)%rake,   1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%rectangularPatch%s(i)%opening,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%rectangularPatch%s(i)%mu0,    1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%rectangularPatch%s(i)%sig,    1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%rectangularPatch%s(i)%a,      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%rectangularPatch%s(i)%b,      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%rectangularPatch%s(i)%L,      1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%rectangularPatch%s(i)%Vo,     1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%rectangularPatch%s(i)%damping,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
       END DO

       DO i=1,in%triangularPatch%ns
          position=0
          CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%triangularPatch%s(i)%Vpl,    1,MPI_REAL8,  MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%triangularPatch%i1(i),     1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%triangularPatch%i2(i),     1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%triangularPatch%i3(i),     1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%triangularPatch%s(i)%rake,   1,MPI_REAL8,  MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%triangularPatch%s(i)%opening,1,MPI_REAL8,  MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%triangularPatch%s(i)%mu0,    1,MPI_REAL8,  MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%triangularPatch%s(i)%sig,    1,MPI_REAL8,  MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%triangularPatch%s(i)%a,      1,MPI_REAL8,  MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%triangularPatch%s(i)%b,      1,MPI_REAL8,  MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%triangularPatch%s(i)%L,      1,MPI_REAL8,  MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%triangularPatch%s(i)%Vo,     1,MPI_REAL8,  MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%triangularPatch%s(i)%damping,1,MPI_REAL8,  MPI_COMM_WORLD,ierr)
       END DO

       IF (0 .NE. in%triangularPatch%ns) THEN
          DO i=1,in%triangularPatch%nVe
             position=0
             CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
             CALL MPI_UNPACK(packed,psize,position,in%triangularPatch%v(1,i),1,MPI_REAL8,MPI_COMM_WORLD,ierr)
             CALL MPI_UNPACK(packed,psize,position,in%triangularPatch%v(2,i),1,MPI_REAL8,MPI_COMM_WORLD,ierr)
             CALL MPI_UNPACK(packed,psize,position,in%triangularPatch%v(3,i),1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          END DO
       END IF

       DO i=1,in%strainVolume%ns
          position=0
          CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%e11,        1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%e12,        1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%e13,        1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%e22,        1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%e23,        1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%e33,        1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%x(1,i),          1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%x(2,i),          1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%x(3,i),          1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%length(i),       1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%width(i),        1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%thickness(i),    1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%strike(i),       1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%dip(i),          1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%Gk,         1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%gammadot0k, 1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%dok,        1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%mk,         1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%Qk,         1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%Rk,         1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%Gm,         1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%gammadot0m, 1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%dom,        1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%mm,         1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%Qm,         1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%Rm,         1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%nGk,        1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%ngammadot0k,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%npowerk,    1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%nQk,        1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%nRk,        1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%nGm,        1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%ngammadot0m,1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%npowerm,    1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%nQm,        1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%nRm,        1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%rhoc,       1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%strainVolume%s(i)%To,         1,MPI_REAL8,MPI_COMM_WORLD,ierr)
       END DO

       IF (0 .LT. in%nObservationState) &
                    ALLOCATE(in%observationState(in%nObservationState,2),STAT=ierr)
       IF (ierr>0) STOP "slave could not allocate memory for observation states"

       DO i=1,in%nObservationState
          position=0
          CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%observationState(i,1),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
       END DO

       IF (0 .LT. in%nObservationPoint) &
                    ALLOCATE(in%observationPoint(in%nObservationPoint), &
                             STAT=ierr)
       IF (ierr>0) STOP "slave could not allocate memory for observation points"

       DO i=1,in%nObservationPoint
          position=0
          CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%observationPoint(i)%name,10,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
          DO k=1,3
             CALL MPI_UNPACK(packed,psize,position,in%observationPoint(i)%x(k),1,MPI_REAL8,MPI_COMM_WORLD,ierr)
          END DO
       END DO

       DO i=1,in%ne
          position=0
          CALL MPI_BCAST(packed,psize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%event(i)%time,1,MPI_REAL8,  MPI_COMM_WORLD,ierr)
          CALL MPI_UNPACK(packed,psize,position,in%event(i)%i,   1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
       END DO

       CALL FLUSH(6)      

       in%nPatch=in%rectangularPatch%ns+in%triangularPatch%ns

    END IF ! master or slaves

    in%nVolume=in%strainVolume%ns

2000 FORMAT ("# ----------------------------------------------------------------------------")
   
  END SUBROUTINE init
   
  !-----------------------------------------------
  !> subroutine printhelp
  !! displays a help message with master thread.
  !-----------------------------------------------
  SUBROUTINE printhelp()

    INTEGER :: rank,size,ierr
    INCLUDE 'mpif.h'

    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)

    IF (0.EQ.rank) THEN
       PRINT '("usage:")'
       PRINT '("")'
       PRINT '("mpirun -n 2 unicycle-3d [-h] [--dry-run] [--help] [--epsilon 1e-6] [filename]")'
       PRINT '("")'
       PRINT '("options:")'
       PRINT '("   -h                      prints this message and aborts calculation")'
       PRINT '("   --dry-run               abort calculation, only output geometry")'
       PRINT '("   --help                  prints this message and aborts calculation")'
       PRINT '("   --version               print version number and exit")'
       PRINT '("   --epsilon               set the numerical accuracy [1E-6]")'
       PRINT '("")'
       PRINT '("description:")'
       PRINT '("   simulates elasto-dynamics on faults in a viscoelastic")'
       PRINT '("   medium following the radiation-damping approximation")'
       PRINT '("   using the integral method.")'
       PRINT '("")'
       PRINT '("see also: ""man unicycle""")'
       PRINT '("")'
       CALL FLUSH(6)
    END IF

  END SUBROUTINE printhelp

  !-----------------------------------------------
  !> subroutine printversion
  !! displays code version with master thread.
  !-----------------------------------------------
  SUBROUTINE printversion()

    INTEGER :: rank,size,ierr
    INCLUDE 'mpif.h'

    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)

    IF (0.EQ.rank) THEN
       PRINT '("unicycle-3d version 1.0.0, compiled on ",a)', __DATE__
       PRINT '("")'
       CALL FLUSH(6)
    END IF

  END SUBROUTINE printversion

END MODULE input

