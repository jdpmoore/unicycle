!-----------------------------------------------------------------------
!> program RateState simulates evolution of fault slip for a point
!! source with the radiation damping approximation.
!!
!! \mainpage
!! 
!! The time evolution is evaluated numerically using the 4th/5th order
!! Runge-Kutta method with adaptive time steps. The state vector is 
!! as follows:
!!
!!   /    t    \
!!   |    s    |
!!   !  theta* |
!!   \   v*    /
!!
!! where t is the traction, s is the slip, v* is the logarithm of the norm
!! of the velocity vector (v*=log10(V)), and theta* is the logarithm of the state 
!! variable (theta*=log10(theta)) in the rate and state friction framework.
!!
!! \author Sylvain Barbot (2017).
!----------------------------------------------------------------------
PROGRAM ratestate

#include "macros.f90"

  USE ode45
  USE types 

  IMPLICIT NONE

  ! error flag
  INTEGER :: ierr

  CHARACTER(256) :: filename

  ! scaling factor
  REAL*8, PARAMETER :: lg10=LOG(1.d1)

  ! state vector
  REAL*8, DIMENSION(STATE_VECTOR_DGF) :: y
  ! rate of change of state vector
  REAL*8, DIMENSION(STATE_VECTOR_DGF) :: dydt,yscal
  ! temporary variables
  REAL*8, DIMENSION(STATE_VECTOR_DGF) :: ytmp,ytmp1,ytmp2,ytmp3

  ! time
  REAL*8 :: time,t0
  ! time step
  REAL*8 :: dt_try,dt_next,dt_done
  ! time steps
  INTEGER :: i,j
  ! current velocity
  REAL*8 :: velocity

  ! maximum number of time steps (default)
  INTEGER :: maximumIterations=1000000

  ! model parameters
  TYPE(SIMULATION_STRUCT) :: in

  ! start time
  time=0.d0

  ! initial tentative time step
  dt_next=1.0d-3

  ! retrieve input parameters from command line
  CALL init(in)
  CALL FLUSH(STDOUT)

  IF (in%isdryrun) THEN
     PRINT '("dry run: abort calculation")'
  END IF
  IF (in%isdryrun .OR. in%isversion .OR. in%ishelp) THEN
     STOP
  END IF

  PRINT 2000

  ! allocate buffer from ode45 module
  ALLOCATE(AK2(STATE_VECTOR_DGF), &
           AK3(STATE_VECTOR_DGF), &
           AK4(STATE_VECTOR_DGF), &
           AK5(STATE_VECTOR_DGF), &
           AK6(STATE_VECTOR_DGF), &
           yrkck(STATE_VECTOR_DGF), STAT=ierr)
  IF (ierr>0) STOP "could not allocate the AK1-6 work space"

  OPEN (UNIT=FPTIME,FILE=in%timeFilename,IOSTAT=ierr,FORM="FORMATTED")
  IF (ierr>0) THEN
     WRITE_DEBUG_INFO(102)
     WRITE (STDERR,'("error: unable to access ",a)') TRIM(in%timefilename)
     STOP 1
  END IF

  ! initialize the y vector
  CALL initStateVector(STATE_VECTOR_DGF,y,in)

  ! initialize output
  WRITE(STDOUT,'("#       n               time                 dt          v")')
  WRITE(STDOUT,'(I9.9,ES19.12E2,ES19.12E2)') 0,time,dt_next
  WRITE(FPTIME,'("#               time                 dt")')
  WRITE(FPTIME,'(ES19.12E2,ES19.12E2)') 0._8,dt_next

  ! initialize time series of state vector
  WRITE (filename,'(a,"/patch.dat")') TRIM(in%wdir)
  OPEN (UNIT=FPSTATE,FILE=filename,IOSTAT=ierr,FORM="FORMATTED")
  IF (ierr>0) THEN
     WRITE_DEBUG_INFO(102)
     WRITE (STDERR,'("error: unable to access ",a)') TRIM(filename)
     STOP 1
  END IF
  WRITE(FPSTATE,'("# time, slip, stress, log10(state), log10(V), slip rate, stress rate, stress rate / stress, a / V")')

  ! main loop
  DO i=1,maximumIterations

     CALL odefun(STATE_VECTOR_DGF,time,y,dydt)

     CALL export()

     dt_try=dt_next
     yscal(:)=abs(y(:))+abs(dt_try*dydt(:))+TINY

     t0=0.d0
     CALL RKQSm(STATE_VECTOR_DGF,t0,y,dydt, &
               yscal,ytmp1,ytmp2,ytmp3,dt_try,dt_done,dt_next,odefun)

     time=time+dt_done

     ! end calculation
     IF (in%interval .LE. time) THEN
        EXIT
     END IF
   
  END DO

  PRINT '(I9.9," time steps.")', i

  CLOSE(FPTIME)

  ! close observation patch file
  CLOSE(FPSTATE)

  DEALLOCATE(AK2,AK3,AK4,AK5,AK6)

2000 FORMAT ("# ----------------------------------------------------------------------------")
     
CONTAINS

  !-----------------------------------------------------------------------
  !> subroutine export
  ! write the state variables
  !----------------------------------------------------------------------
  SUBROUTINE export()

    IMPLICIT NONE

    ! counter
    INTEGER :: k

    ! index in state vector
    INTEGER :: index

    ! format string
    CHARACTER(1024) :: formatString

    formatString="(ES19.12E2"
    DO k=1,STATE_VECTOR_DGF
       formatString=TRIM(formatString)//",X,ES20.12E3,X,ES20.12E3"
    END DO
    formatString=TRIM(formatString)//")"

    WRITE (FPSTATE,TRIM(formatString)) time, &
              y(index+1:index+STATE_VECTOR_DGF), &
           dydt(index+1:index+STATE_VECTOR_DGF)

    WRITE(FPTIME,'(ES19.12E2,ES19.12E2)') time,dt_done
    IF (0 .EQ. MOD(i,50)) THEN
       WRITE(STDOUT,'(I9.9,ES19.12E2,ES19.12E2,ES11.4E2)') i,time,dt_done,velocity
       CALL FLUSH(STDOUT)
    END IF

  END SUBROUTINE export

  !-----------------------------------------------------------------------
  !> subroutine initStateVector
  ! initialize the state vector
  !
  ! INPUT:
  ! @param n - number of state elements own by current thread
  ! @param y - the state vector (segment owned by currect thread)
  !
  !----------------------------------------------------------------------
  SUBROUTINE initStateVector(n,y,in)

    IMPLICIT NONE

    INTEGER, INTENT(IN)   :: n
    REAL*8, INTENT(OUT)    :: y(n)
    TYPE(SIMULATION_STRUCT), INTENT(IN) :: in

    ! loop over elements owned by current thread
    y(STATE_VECTOR_SLIP) = 0._8

    ! traction
    y(STATE_VECTOR_TRACTION) = in%s%sig*(in%s%mu0+(in%s%a-in%s%b)*log(in%s%Vpl/in%s%Vo)) &
                                         -in%s%damping*in%s%Vpl

    ! state variable log10(theta)
    y(STATE_VECTOR_STATE_1) = log(in%s%L/in%s%Vpl)/lg10

    ! slip velocity log10(V)
    y(STATE_VECTOR_VELOCITY) = log(in%s%Vpl*0.98d0)/lg10

  END SUBROUTINE initStateVector

  !-----------------------------------------------------------------------
  !> subroutine odefun
  ! evalutes the derivative of the state vector
  !
  ! @param n - number of state elements own by current thread
  ! @param m - degrees of freedom
  !
  ! DESCRIPTION:
  !   1- extract slip velocity and strain rate from state vector
  !   2- calculate the rate of traction and rate of stress
  !   3- calculate the rate of remaining state variables
  !----------------------------------------------------------------------
  SUBROUTINE odefun(n,time,y,dydt)

    IMPLICIT NONE

    INTEGER, INTENT(IN)   :: n
    REAL*8, INTENT(IN)    :: time
    REAL*8, INTENT(IN)    :: y(n)
    REAL*8, INTENT(INOUT) :: dydt(n)

    ! traction
    REAL*8 :: t

    ! relative velocity
    REAL*8 :: v

    !--------------------------------------------------------------------
    ! step 1/3 - extract slip velocity and strain rate from state vector
    !--------------------------------------------------------------------

    ! slip velocity
    velocity=EXP(y(STATE_VECTOR_VELOCITY)*lg10)

    ! update state vector (rate of slip components)
    dydt(STATE_VECTOR_SLIP)=velocity

    v=velocity-in%s%Vpl

    !-----------------------------------------------------------------
    ! step 2/3 - calculate the rate of traction and rate of stress
    !-----------------------------------------------------------------

    t=-7._8*PI/16._8*in%mu/in%s%width*v

    !-----------------------------------------------------------------
    ! step 3/3 - calculate the rate of remaining state variables
    !-----------------------------------------------------------------

    ! rate of state
    dydt(STATE_VECTOR_STATE_1)=(EXP(-y(STATE_VECTOR_STATE_1)*lg10)-velocity/in%s%L)/lg10

    ! acceleration
    dydt(STATE_VECTOR_VELOCITY)=(t-in%s%b*in%s%sig*dydt(STATE_VECTOR_STATE_1)*lg10) / &
                 (in%s%a*in%s%sig+in%s%damping*velocity) / lg10

    ! return the traction rate
    dydt(STATE_VECTOR_TRACTION)=t-in%s%damping*velocity*dydt(STATE_VECTOR_VELOCITY)*lg10

  END SUBROUTINE odefun

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
  
    IMPLICIT NONE

    TYPE(SIMULATION_STRUCT), INTENT(OUT) :: in
  
    CHARACTER :: ch
    CHARACTER(256) :: dataline,filename
    INTEGER :: iunit,noptions
    TYPE(OPTION_S) :: opts(6)
  
    INTEGER :: i,ierr,k
    INTEGER :: dummy
  
    ! define long options, such as --dry-run
    ! parse the command line for options
    opts(1)=OPTION_S("version",.FALSE.,CHAR(21))
    opts(2)=OPTION_S("dry-run",.FALSE.,CHAR(22))
    opts(3)=OPTION_S("epsilon",.TRUE.,'e')
    opts(4)=OPTION_S("maximum-step",.TRUE.,'m')
    opts(5)=OPTION_S("maximum-iterations",.TRUE.,'i')
    opts(6)=OPTION_S("help",.FALSE.,'h')
  
    noptions=0
    DO
       ch=getopt("he:i:m:",opts)
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
          ! numerical accuracy (variable epsilon sits in the ode45 module)
          READ(optarg,*) epsilon
          noptions=noptions+1
       CASE('i')
          ! maximum number of iterations
          READ(optarg,*) maximumIterations
          noptions=noptions+1
       CASE('m')
          ! maximum time step (variable maximumTimeStep sits in the ode45 module)
          READ(optarg,*) maximumTimeStep
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
       STOP
    END IF
  
    IF (in%ishelp) THEN
       CALL printhelp()
       ! abort parameter input
       STOP
    END IF
  
    PRINT 2000
    PRINT '("# RATESTATE")'
    PRINT '("# quasi-dynamic slip evolution for a point source")'
    PRINT '("# numerical accuracy: ",ES11.4)', epsilon
    PRINT '("# maximum iterations: ",I11)', maximumIterations
    PRINT '("# maximum time step: ",ES12.4)', maximumTimeStep

    PRINT 2000
  
    IF (noptions .LT. COMMAND_ARGUMENT_COUNT()) THEN
       ! read from input file
       iunit=25
       CALL GET_COMMAND_ARGUMENT(noptions+1,filename)
       OPEN (UNIT=iunit,FILE=filename,IOSTAT=ierr)
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
          IOSTAT=ierr,FORM="FORMATTED")
    IF (ierr>0) THEN
       WRITE_DEBUG_INFO(102)
       WRITE (STDERR,'("error: unable to access ",a)') TRIM(in%timefilename)
       STOP 1
    END IF
    CLOSE(FPTIME)
     
    PRINT '("# shear modulus")'
    CALL getdata(iunit,dataline)
    READ  (dataline,*) in%mu
    PRINT '(ES9.2E1)', in%mu
  
    IF (0 .GT. in%mu) THEN
       WRITE_DEBUG_INFO(-1)
       WRITE (STDERR,'(a)') TRIM(dataline)
       WRITE (STDERR,'("input error: shear modulus must be positive")')
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
    !                 G E O M E T R Y
    ! - - - - - - - - - - - - - - - - - - - - - - - - - -
    PRINT '("#    Vpl      radius")'
    PRINT 2000
    CALL getdata(iunit,dataline)
    READ (dataline,*,IOSTAT=ierr) &
          in%s%Vpl, &
          in%s%width
   
    PRINT '(2ES9.2E1)', &
         in%s%Vpl, &
         in%s%width

    PRINT 2000
    ! flush standard output
    CALL FLUSH(6)      

    ! - - - - - - - - - - - - - - - - - - - - - - - - - -
    !        F R I C T I O N   P R O P E R T I E S
    ! - - - - - - - - - - - - - - - - - - - - - - - - - -
    PRINT '("#    mu0      sig        a        b        L       Vo  G/(2Vs)")'
    PRINT 2000
    CALL getdata(iunit,dataline)
    READ (dataline,*,IOSTAT=ierr) &
          in%s%mu0, &
          in%s%sig, &
          in%s%a, &
          in%s%b, &
          in%s%L, &
          in%s%Vo, &
          in%s%damping
   
    PRINT '(7ES9.2E1)', &
         in%s%mu0, &
         in%s%sig, &
         in%s%a, &
         in%s%b, &
         in%s%L, &
         in%s%Vo, &
         in%s%damping
   
    PRINT 2000
    ! flush standard output
    CALL FLUSH(6)      

2000 FORMAT ("# ----------------------------------------------------------------------------")
     
  END SUBROUTINE init

  !-----------------------------------------------
  !> subroutine printhelp
  !! displays a help message with master thread.
  !-----------------------------------------------
  SUBROUTINE printhelp()
  
    INTEGER :: ierr
    
    PRINT '("usage:")'
    PRINT '("")'
    PRINT '("unicycle-0d-ratestate [-h] [--dry-run] [--help] [--epsilon 1e-6] [filename]")'
    PRINT '("")'
    PRINT '("options:")'
    PRINT '("   -h                      prints this message and aborts calculation")'
    PRINT '("   --dry-run               abort calculation, only output geometry")'
    PRINT '("   --help                  prints this message and aborts calculation")'
    PRINT '("   --version               print version number and exit")'
    PRINT '("   --epsilon               set the numerical accuracy [1E-6]")'
    PRINT '("   --maximum-iterations    set the maximum time step [1000000]")'
    PRINT '("   --maximum-step          set the maximum time step [none]")'
    PRINT '("")'
    PRINT '("description:")'
    PRINT '("   simulates elasto-dynamics for a point source")'
    PRINT '("   following the radiation-damping approximation")'
    PRINT '("   using the integral method.")'
    PRINT '("")'
    PRINT '("   if filename is not provided, reads from standard input.")'
    PRINT '("")'
    PRINT '("see also: ""man unicycle""")'
    PRINT '("")'
    CALL FLUSH(6)
    
  END SUBROUTINE printhelp
  
  !-----------------------------------------------
  !> subroutine printversion
  !! displays code version with master thread.
  !-----------------------------------------------
  SUBROUTINE printversion()
    
    INTEGER :: ierr
    
    PRINT '("unicycle-0d-ratestate version 1.0.0, compiled on ",a)', __DATE__
    PRINT '("")'
    CALL FLUSH(6)
  
  END SUBROUTINE printversion
    
END PROGRAM ratestate

