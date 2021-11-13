!> program Unicycle (unified cycle of earthquakes) simulates evolution
!! of fault slip and distributed strain using the integral method under
!! the radiation damping approximation.
!!
!! \mainpage
!! 
!! The Green's function for traction and stress interaction amongst 
!! triangular and rectangular dislocations has the following layout
!!
!!       / RR RT \
!!   G = |       |
!!       \ TR TT /
!!
!! where
!!
!!   RR is the matrix for traction on rectangular faults due to rectangular fault slip
!!   RT is the matrix for traction on rectangular faults due to triangular  fault slip
!!   TR is the matrix for traction on triangular  faults due to rectangular fault slip
!!   TT is the matrix for traction on triangular  faults due to triangular  fault slip
!!
!! All faults (rectangular or triangular) have two slip directions:
!! strike slip and dip slip. The interaction matrices become
!!
!!        / RRss  RRsd \        / RTss  RRsd \
!!   RR = |            |,  RT = |            |,  
!!        \ RRds  RRdd /        \ RRds  RRdd /
!!
!!        / TRss  TRsd \        / TTss  TTsd \
!!   TR = |            |,  TT = |            |
!!        \ TRds  TRdd /        \ TTds  TTdd /
!!
!! The time evolution is evaluated numerically using the 4th/5th order
!! Runge-Kutta method with adaptive time steps. The state vector is 
!! as follows:
!!
!!    / R1 1       \   +-----------------------+  +-----------------+
!!    | .          |   |                       |  |                 |
!!    | R1 dPatch  |   |                       |  |                 |
!!    | .          |   |  nRectangle * dPatch  |  |                 |
!!    | .          |   |                       |  |                 |
!!    | Rn 1       |   |                       |  |                 |
!!    | .          |   |                       |  |                 |
!!    | Rn dPatch  |   +-----------------------+  |                 |
!!    |            |                              | nPatch * dPatch |
!!    | T1 1       |   +-----------------------+  |                 |
!!    | .          |   |                       |  |                 |
!!    | T1 dPatch  |   |                       |  |                 |
!!    | .          |   |  nTriangle * dPatch   |  |                 |
!!    | .          |   |                       |  |                 |
!!    | Tn 1       |   |                       |  |                 |
!!    | .          |   |                       |  |                 |
!!    \ Tn dPatch  /   +-----------------------+  +-----------------+
!!
!! where nRectangle and nTriangle are the number of rectangular and
!! triangular patches, respectively, and dPatch is the degrees of 
!! freedom for patches.
!!
!! For every rectangular or triangular patch, we have the following 
!! items in the state vector
!!
!!   /  ss   \  1
!!   |  sd   |  .
!!   |  ts   |  .
!!   |  td   |  .
!!   |  tn   |  .
!!   ! theta*|  .
!!   \  v*   /  dPatch
!!
!! where ts, td, and td are the local traction in the strike, dip, and
!! normal directions, ss and sd are the total slip in the strike and dip 
!! directions, v*=log10(v) is the logarithm of the norm of the velocity,
!! and theta*=log10(theta) is the logarithm of the state variable in the 
!! rate and state friction framework.
!!
!! References:<br>
!!
!!   James D P Moore, Sylvain Barbot, Lujia Feng, Yu Hang, Valere Lambert, 
!!   Eric Lindsey, Sagar Masuti, Takanori Matsuzawa, Jun Muto, 
!!   Priyamvada Nanjundiah, Rino Salman, Sharadha Sathiakumar, 
!!   and Harpreet Sethi. (2019, September 25). 
!!   jdpmoore/unicycle: Unicycle (Version 1.0). Zenodo.
!!   http://doi.org/10.5281/zenodo.5688288
!!
!!   Barbot S., J. D.-P. Moore, and V. Lambert, "Displacement and Stress
!!   Associated with Distributed Anelastic Deformation in a Half-Space",
!!   Bull. Seism. Soc. Am., 10.1785/0120160237, 2017.
!!
!!   Jun Muto, James D P Moore, Sylvain Barbot, Iinuma T, Ohta Y, 
!!   Horiuchi S, Hikaru I. Coupled afterslip and transient mantle flow 
!!   after the 2011 Tohoku earthquake, Science Advances 2019.
!!   https://doi.org/10.1126/sciadv.aaw1164
!!
!!   Lambert, V., and S. Barbot. "Contribution of viscoelastic flow in 
!!   earthquake cycles within the lithosphere‐asthenosphere system." 
!!   Geophysical Research Letters 43.19 (2016).
!!
!!   Qiu, Q., E. M. Hill, S. Barbot, J. Hubbard, W. Feng, E. O. Lindsey, 
!!   L. Feng, K. Dai, S. V. Samsonov and P. Tapponnier "The mechanism of 
!!   partial rupture of a locked megathrust: The role of fault morphology." 
!!   Geology 44.10 (2016): 875-878.
!!
!! \author Sylvain Barbot (2017).
!----------------------------------------------------------------------
PROGRAM ratestate

#include "macros.f90"

  USE greens
  USE ode45
  USE types 

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  REAL*8, PARAMETER :: DEG2RAD = 0.01745329251994329547437168059786927_8

  ! MPI rank and size
  INTEGER :: rank,csize

  ! error flag
  INTEGER :: ierr

  CHARACTER(512) :: filename

  ! maximum velocity
  REAL*8 :: vMax

  ! scaling factor
  REAL*8, PARAMETER :: lg10=LOG(1.d1)

  ! state vector
  REAL*8, DIMENSION(:), ALLOCATABLE :: y
  ! rate of change of state vector
  REAL*8, DIMENSION(:), ALLOCATABLE :: dydt,yscal
  ! temporary variables
  REAL*8, DIMENSION(:), ALLOCATABLE :: ytmp,ytmp1,ytmp2,ytmp3

  ! slip velocity vector and strain rate tensor array (temporary space)
  REAL*8, DIMENSION(:), ALLOCATABLE :: v

  ! slip velocity vector and strain rate tensor array (temporary space)
  REAL*8, DIMENSION(:), ALLOCATABLE :: vAll

  ! traction vector and stress tensor array (temporary space)
  REAL*8, DIMENSION(:), ALLOCATABLE :: t

  ! displacement and kinematics vector
  REAL*8, DIMENSION(:), ALLOCATABLE :: d,u,dAll

  ! Green's functions
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: G,O,Of,Ol

  ! time
  REAL*8 :: time,t0
  ! time step
  REAL*8 :: dt_try,dt_next,dt_done
  ! time steps
  INTEGER :: i,j

  ! maximum number of time steps (default)
  INTEGER :: maximumIterations=1000000

  ! layout for parallelism
  TYPE(LAYOUT_STRUCT) :: layout

  ! model parameters
  TYPE(SIMULATION_STRUCT) :: in

  ! initialization
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,csize,ierr)

  ! start time
  time=0.d0

  ! initial tentative time step
  dt_next=1.0d-3

  ! retrieve input parameters from command line
  CALL init(in)
  CALL FLUSH(STDOUT)

  ! export geometrical layout information
  CALL exportGeometry(in)

  IF (in%isdryrun .AND. 0 .EQ. rank) THEN
     PRINT '("dry run: abort calculation")'
  END IF
  IF (in%isdryrun .OR. in%isversion .OR. in%ishelp) THEN
     CALL MPI_FINALIZE(ierr)
     STOP
  END IF

  ! calculate basis vectors
  CALL initGeometry(in)

  ! describe data layout for parallelism
  CALL initParallelism(in,layout)

  ! calculate the stress interaction matrix
  IF (0 .EQ. rank) THEN
     PRINT '("# computing Green''s functions.")'
  END IF

  CALL buildG(in,layout,G)
  CALL buildO(in,layout,O,Of,Ol)
  DEALLOCATE(Of,Ol)

  IF (0 .EQ. rank) THEN
     PRINT 2000
  END IF

  ! velocity vector and strain rate tensor array (t=G*vAll)
  ALLOCATE(u(layout%listVelocityN(1+rank)), &
           v(layout%listVelocityN(1+rank)), &
           vAll(SUM(layout%listVelocityN)),STAT=ierr)
  IF (ierr>0) STOP "could not allocate the velocity and strain rate vector"

  ! traction vector and stress tensor array (t=G*v)
  ALLOCATE(t(layout%listForceN(1+rank)),STAT=ierr)
  IF (ierr>0) STOP "could not allocate the traction and stress vector"

  ! displacement vector (d=O*v)
  ALLOCATE(d   (in%nObservationPoint*DISPLACEMENT_VECTOR_DGF), &
           dAll(in%nObservationPoint*DISPLACEMENT_VECTOR_DGF),STAT=ierr)
  IF (ierr>0) STOP "could not allocate the displacement vector"

  ! state vector
  ALLOCATE(y(layout%listStateN(1+rank)),STAT=ierr)
  IF (ierr>0) STOP "could not allocate the state vector"

  ! rate of state vector
  ALLOCATE(dydt(layout%listStateN(1+rank)), &
           yscal(layout%listStateN(1+rank)),STAT=ierr)
  IF (ierr>0) STOP "could not allocate the state vectors"

  ALLOCATE(ytmp (layout%listStateN(1+rank)), &
           ytmp1(layout%listStateN(1+rank)), &
           ytmp2(layout%listStateN(1+rank)), &
           ytmp3(layout%listStateN(1+rank)),STAT=ierr)
  IF (ierr>0) STOP "could not allocate the RKQS work space"

  ! allocate buffer from ode45 module
  ALLOCATE(AK2(layout%listStateN(1+rank)), &
           AK3(layout%listStateN(1+rank)), &
           AK4(layout%listStateN(1+rank)), &
           AK5(layout%listStateN(1+rank)), &
           AK6(layout%listStateN(1+rank)), &
           yrkck(layout%listStateN(1+rank)), STAT=ierr)
  IF (ierr>0) STOP "could not allocate the AK1-6 work space"

  IF (0 .EQ. rank) THEN
     OPEN (UNIT=FPTIME,FILE=in%timeFilename,IOSTAT=ierr,FORM="FORMATTED")
     IF (ierr>0) THEN
        WRITE_DEBUG_INFO(102)
        WRITE (STDERR,'("error: unable to access ",a)') TRIM(in%timefilename)
        STOP 1
     END IF
  END IF

  ! initialize the y vector
  CALL initStateVector(layout%listStateN(1+rank),y,in)

  ! initialize output
  IF (0 .EQ. rank) THEN
     WRITE(STDOUT,'("#       n               time                 dt       vMax")')
     WRITE(STDOUT,'(I9.9,ES19.12E2,ES19.12E2)') 0,time,dt_next
     WRITE(FPTIME,'("#               time                 dt")')
     WRITE(FPTIME,'(ES19.12E2,ES19.12E2)') 0._8,dt_next
  END IF

  ! initialize observation patch
  DO j=1,in%nObservationState
     IF ((in%observationState(j,1) .GE. layout%listOffset(rank+1)) .AND. &
         (in%observationState(j,1) .LT. layout%listOffset(rank+1)+layout%listElements(rank+1))) THEN

        in%observationState(j,2)=100+j
        WRITE (filename,'(a,"/patch-",I8.8,".dat")') TRIM(in%wdir),in%observationState(j,1)
        OPEN (UNIT=in%observationState(j,2), &
              FILE=filename,IOSTAT=ierr,FORM="FORMATTED")
        IF (ierr>0) THEN
           WRITE_DEBUG_INFO(102)
           WRITE (STDERR,'("error: unable to access ",a)') TRIM(filename)
           STOP 1
        END IF
     END IF
  END DO

  ! initialize observation points
  IF (0 .EQ. rank) THEN
     DO j=1,in%nObservationPoint
        in%observationPoint(j)%file=1000+j
        WRITE (filename,'(a,"/opts-",a,".dat")') TRIM(in%wdir),TRIM(in%observationPoint(j)%name)
        OPEN (UNIT=in%observationPoint(j)%file, &
              FILE=filename,IOSTAT=ierr,FORM="FORMATTED")
        IF (ierr>0) THEN
           WRITE_DEBUG_INFO(102)
           WRITE (STDERR,'("error: unable to access ",a)') TRIM(filename)
           STOP 1
        END IF
     END DO
  END IF

  ! main loop
  DO i=1,maximumIterations

     CALL odefun(layout%listStateN(1+rank),time,y,dydt)

     CALL export()
     CALL exportPoints()

     dt_try=dt_next
     yscal(:)=abs(y(:))+abs(dt_try*dydt(:))+TINY

     t0=0.d0
     CALL RKQSm(layout%listStateN(1+rank),t0,y,dydt, &
               yscal,ytmp1,ytmp2,ytmp3,dt_try,dt_done,dt_next,odefun)

     time=time+dt_done

     ! end calculation
     IF (in%interval .LE. time) THEN
        EXIT
     END IF
   
  END DO

  IF (0 .EQ. rank) THEN
     PRINT '(I9.9," time steps.")', i
  END IF

  IF (0 .EQ. rank) THEN
     CLOSE(FPTIME)
  END IF

  ! close observation state files
  DO j=1,in%nObservationState
     IF ((in%observationState(j,1) .GT. layout%listOffset(rank+1)) .AND. &
         (in%observationState(j,1) .LE. layout%listOffset(rank+1)+layout%listElements(rank+1))) THEN
        CLOSE(in%observationState(j,2))
     END IF
  END DO

  ! close the observation points
  IF (0 .EQ. rank) THEN
     DO j=1,in%nObservationPoint
        CLOSE(in%observationPoint(j)%file)
     END DO
  END IF

  DEALLOCATE(y,dydt,yscal)
  DEALLOCATE(ytmp,ytmp1,ytmp2,ytmp3,yrkck)
  DEALLOCATE(AK2,AK3,AK4,AK5,AK6)
  DEALLOCATE(G,v,vAll,t)
  DEALLOCATE(layout%listForceN)
  DEALLOCATE(layout%listVelocityN,layout%listVelocityOffset)
  DEALLOCATE(layout%listStateN,layout%listStateOffset)
  DEALLOCATE(layout%elementStateIndex)
  DEALLOCATE(layout%listElements,layout%listOffset)
  DEALLOCATE(O,d,u,dAll)

  CALL MPI_FINALIZE(ierr)

2000 FORMAT ("# ----------------------------------------------------------------------------")
     
CONTAINS

  !-----------------------------------------------------------------------
  !> subroutine exportGeometry
  ! export geometry of fault patches
  !----------------------------------------------------------------------
  SUBROUTINE exportGeometry(in)
    IMPLICIT NONE
    TYPE(SIMULATION_STRUCT), INTENT(IN) :: in

    CHARACTER(512) :: filename

    IF (0 .EQ. rank) THEN
       WRITE (filename,'(a,"/patch-rectangular.vtp")') TRIM(in%wdir)
       OPEN (UNIT=FPVTP,FILE=filename,IOSTAT=ierr,FORM="FORMATTED")
       IF (ierr>0) THEN
          WRITE_DEBUG_INFO(102)
          WRITE (STDERR,'("error: unable to access ",a)') TRIM(filename)
          STOP 1
       END IF
       CALL exportvtk_rfaults(in%rectangularPatch%ns, &
                              in%rectangularPatch%x, &
                              in%rectangularPatch%length, &
                              in%rectangularPatch%width, &
                              in%rectangularPatch%strike, &
                              in%rectangularPatch%dip, &
                              in%rectangularPatch%s, &
                              FPVTP)
       CLOSE(FPVTP)
    END IF

  END SUBROUTINE exportGeometry

  !-----------------------------------------------------------------------
  !> subroutine exportPoints
  ! export observation points
  !----------------------------------------------------------------------
  SUBROUTINE exportPoints()

    IMPLICIT NONE

    INTEGER :: i,k,l,ierr
    INTEGER :: elementIndex,elementType
    CHARACTER(1024) :: formatString
    TYPE(PATCH_ELEMENT_STRUCT) :: patch

    !-----------------------------------------------------------------
    ! step 1/4 - gather the kinematics from the state vector
    !-----------------------------------------------------------------

    ! element index in d vector
    k=1
    ! element index in state vector
    l=1
    ! loop over elements owned by current thread
    DO i=1,SIZE(layout%elementIndex)
       elementIndex=layout%elementIndex(i)
       elementType= layout%elementType(i)

       SELECT CASE (elementType)
       CASE (FLAG_RECTANGLE)
             patch=in%rectangularPatch%s(elementIndex)
       CASE (FLAG_TRIANGLE)
             patch=in%triangularPatch%s(elementIndex)
       CASE DEFAULT
          WRITE(STDERR,'("wrong case: this is a bug.")')
          WRITE_DEBUG_INFO(-1)
          STOP -1
       END SELECT

       ! strike slip and dip slip
       u(k:k+layout%elementVelocityDGF(i)-1)= (/ &
               y(l+STATE_VECTOR_SLIP_STRIKE), &
               y(l+STATE_VECTOR_SLIP_DIP) /)

       l=l+in%dPatch

       k=k+layout%elementVelocityDGF(i)
    END DO

    !-----------------------------------------------------------------
    ! step 2/4 - calculate the rate of traction and rate of stress
    !-----------------------------------------------------------------

    ! use the BLAS library to compute the matrix vector product
    CALL DGEMV("T",SIZE(O,1),SIZE(O,2), &
                1._8,O,SIZE(O,1),u,1,0.d0,d,1)

    !-----------------------------------------------------------------
    ! step 3/4 - master thread adds the contribution of all elemets
    !-----------------------------------------------------------------

    CALL MPI_REDUCE(d,dAll,in%nObservationPoint*DISPLACEMENT_VECTOR_DGF, &
                    MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    !-----------------------------------------------------------------
    ! step 4/4 - master thread writes to disk
    !-----------------------------------------------------------------

    IF (0 .EQ. rank) THEN
       formatString="(ES19.12E2 "
       DO i=1,DISPLACEMENT_VECTOR_DGF
          formatString=TRIM(formatString)//" ES19.12E2"
       END DO
       formatString=TRIM(formatString)//")"

       ! element index in d vector
       k=1
       DO i=1,in%nObservationPoint
          WRITE (in%observationPoint(i)%file,TRIM(formatString)) time,dAll(k:k+DISPLACEMENT_VECTOR_DGF-1)
          k=k+DISPLACEMENT_VECTOR_DGF
       END DO
    END IF

  END SUBROUTINE exportPoints

  !-----------------------------------------------------------------------
  !> subroutine export
  ! write the state variables of patch elements, and other information.
  !----------------------------------------------------------------------
  SUBROUTINE export()

    ! degrees of freedom
    INTEGER :: dgf

    ! counters
    INTEGER :: j,k

    ! index in state vector
    INTEGER :: index

    ! maximum velocity
    REAL*8 :: vMaxAll

    ! format string
    CHARACTER(1024) :: formatString

    ! export observation state
    DO j=1,in%nObservationState
       IF ((in%observationState(j,1) .GE. layout%listOffset(rank+1)) .AND. &
           (in%observationState(j,1) .LT. layout%listOffset(rank+1)+layout%listElements(rank+1))) THEN

          dgf=in%dPatch

          formatString="(ES19.12E2 "
          DO k=1,dgf
             formatString=TRIM(formatString)//" ES20.12E3 ES20.12E3"
          END DO
          formatString=TRIM(formatString)//")"

          index=layout%elementStateIndex(in%observationState(j,1)-layout%listOffset(rank+1)+1)-dgf
          WRITE (in%observationState(j,2),TRIM(formatString)) time, &
                    y(index+1:index+dgf), &
                 dydt(index+1:index+dgf)
       END IF
    END DO

    IF (0 .EQ. MOD(i,50)) THEN
       CALL MPI_REDUCE(vMax,vMaxAll,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    END IF

    IF (0 .EQ. rank) THEN
       WRITE(FPTIME,'(ES19.12E2,ES19.12E2)') time,dt_done
       IF (0 .EQ. MOD(i,50)) THEN
          WRITE(STDOUT,'(I9.9,ES19.12E2,ES19.12E2,ES11.4E2)') i,time,dt_done,vMaxAll
          CALL FLUSH(STDOUT)
       END IF
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

    INTEGER :: i,l
    INTEGER :: elementType,elementIndex
    TYPE(PATCH_ELEMENT_STRUCT) :: patch

    REAL*8 :: tau

    ! zero out state vector
    y=0._8

    ! element index in state vector
    l=1
    ! loop over elements owned by current thread
    DO i=1,SIZE(layout%elementIndex)
       elementType= layout%elementType(i)
       elementIndex=layout%elementIndex(i)

       SELECT CASE (elementType)
       CASE (FLAG_RECTANGLE)
             patch=in%rectangularPatch%s(elementIndex)
       CASE (FLAG_TRIANGLE)
             patch=in%triangularPatch%s(elementIndex)
       CASE DEFAULT
          WRITE(STDERR,'("wrong case: this is a bug.")')
          WRITE_DEBUG_INFO(-1)
          STOP -1
       END SELECT

       ! strike slip
       y(l+STATE_VECTOR_SLIP_STRIKE) = 0._8

       ! dip slip
       y(l+STATE_VECTOR_SLIP_DIP) = 0._8

       ! traction
       tau = patch%sig*(patch%mu0+(patch%a-patch%b)*log(patch%Vpl/patch%Vo))-patch%damping*patch%Vpl

       ! traction in strike direction
       y(l+STATE_VECTOR_TRACTION_STRIKE) = tau*COS(patch%rake)

       ! traction in dip direction
       y(l+STATE_VECTOR_TRACTION_DIP) = tau*SIN(patch%rake)

       ! traction in normal direction
       y(l+STATE_VECTOR_TRACTION_NORMAL) = 0._8

       ! state variable log10(theta)
       y(l+STATE_VECTOR_STATE_1) = log(patch%L/patch%Vpl)/lg10

       ! slip velocity log10(V)
       y(l+STATE_VECTOR_VELOCITY) = log(patch%Vpl*0.99d0)/lg10

       l=l+in%dPatch

    END DO

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

    INTEGER :: i,j,k,l,ierr
    INTEGER :: elementType,elementIndex
    TYPE(PATCH_ELEMENT_STRUCT) :: patch
    REAL*8 :: correction

    ! traction components in the strike and dip directions
    REAL*8 :: ts, td

    ! scalar rate of shear traction
    REAL*8 :: dtau

    ! velocity scalar
    REAL*8 :: velocity

    ! slip velocity in the strike and dip directions
    REAL*8 :: vs,vd

    ! rake of traction and velocity
    REAL*8 :: rake

    ! maximum velocity
    vMax=0._8

    ! initialize rate of state vector
    dydt=0._8

    !--------------------------------------------------------------------
    ! step 1/3 - extract slip velocity and strain rate from state vector
    !--------------------------------------------------------------------

    ! element index in v vector
    k=1
    ! element index in state vector
    l=1
    ! loop over elements owned by current thread
    DO i=1,SIZE(layout%elementIndex)
       elementType= layout%elementType(i)
       elementIndex=layout%elementIndex(i)

       SELECT CASE (elementType)
       CASE (FLAG_RECTANGLE)
          patch=in%rectangularPatch%s(elementIndex)
       CASE (FLAG_TRIANGLE)
          patch=in%triangularPatch%s(elementIndex)
       CASE DEFAULT
          WRITE(STDERR,'("wrong case: this is a bug.")')
          WRITE_DEBUG_INFO(-1)
          STOP -1
       END SELECT

       ! traction and rake
       ts=y(l+STATE_VECTOR_TRACTION_STRIKE)
       td=y(l+STATE_VECTOR_TRACTION_DIP)
       rake=ATAN2(td,ts)

       ! slip velocity
       velocity=DEXP(y(l+STATE_VECTOR_VELOCITY)*lg10)
       vs=velocity*COS(rake)
       vd=velocity*SIN(rake)

       ! maximum velocity
       vMax=MAX(velocity,vMax)

       ! update state vector (rate of slip components)
       dydt(l+STATE_VECTOR_SLIP_STRIKE)=vs
       dydt(l+STATE_VECTOR_SLIP_DIP   )=vd

       ! v(k:k+layout%elementVelocityDGF(i)-1) = slip velocity
       v(k:k+layout%elementVelocityDGF(i)-1)=(/ &
                     vs-patch%Vpl*COS(patch%rake), &
                     vd-patch%Vpl*SIN(patch%rake) /)

       l=l+in%dPatch

       k=k+layout%elementVelocityDGF(i)
    END DO

    ! all threads gather velocity from all threads
    CALL MPI_ALLGATHERV(v,layout%listVelocityN(1+rank),MPI_REAL8, &
                        vAll,layout%listVelocityN,layout%listVelocityOffset,MPI_REAL8, &
                        MPI_COMM_WORLD,ierr)

    !-----------------------------------------------------------------
    ! step 2/3 - calculate the rate of traction and rate of stress
    !-----------------------------------------------------------------

    ! use the BLAS library to compute the matrix vector product
    CALL DGEMV("T",SIZE(G,1),SIZE(G,2), &
                1._8,G,SIZE(G,1),vAll,1,0.d0,t,1)

    !-----------------------------------------------------------------
    ! step 3/3 - calculate the rate of remaining state variables
    !-----------------------------------------------------------------

    ! element index in t vector
    j=1
    ! element index in state vector
    l=1
    ! loop over elements owned by current thread
    DO i=1,SIZE(layout%elementIndex)
       elementType= layout%elementType(i)
       elementIndex=layout%elementIndex(i)

       SELECT CASE (elementType)
       CASE (FLAG_RECTANGLE)
          patch=in%rectangularPatch%s(elementIndex)
       CASE (FLAG_TRIANGLE)
          patch=in%triangularPatch%s(elementIndex)
       CASE DEFAULT
          WRITE(STDERR,'("wrong case: this is a bug.")')
          WRITE_DEBUG_INFO(-1)
          STOP -1
       END SELECT

       ! traction and rake
       ts=y(l+STATE_VECTOR_TRACTION_STRIKE)
       td=y(l+STATE_VECTOR_TRACTION_DIP)
       rake=ATAN2(td,ts)

       ! slip velocity
       velocity=DEXP(y(l+STATE_VECTOR_VELOCITY)*lg10)

       ! rate of state
       dydt(l+STATE_VECTOR_STATE_1)=(EXP(-y(l+STATE_VECTOR_STATE_1)*lg10)-velocity/patch%L)/lg10

       ! scalar rate of shear traction
       dtau=t(j+TRACTION_VECTOR_STRIKE)*COS(rake) &
           +t(j+TRACTION_VECTOR_DIP)   *SIN(rake)
          
       ! acceleration
       dydt(l+STATE_VECTOR_VELOCITY)=(dtau-patch%b*patch%sig*dydt(l+STATE_VECTOR_STATE_1)*lg10) / &
                                     (patch%a*patch%sig+patch%damping*velocity) / lg10

       ! correction
       correction=patch%damping*velocity*dydt(l+STATE_VECTOR_VELOCITY)*lg10

       ! traction rate
       dydt(l+STATE_VECTOR_TRACTION_STRIKE)=t(j+TRACTION_VECTOR_STRIKE)-correction*COS(rake)
       dydt(l+STATE_VECTOR_TRACTION_DIP   )=t(j+TRACTION_VECTOR_DIP   )-correction*SIN(rake)
       dydt(l+STATE_VECTOR_TRACTION_NORMAL)=t(j+TRACTION_VECTOR_NORMAL)

       l=l+in%dPatch

       j=j+layout%elementForceDGF(i)
    END DO

  END SUBROUTINE odefun

  !-----------------------------------------------------------------------
  !> subroutine initParallelism()
  !! initialize variables describe the data layout
  !! for parallelism.
  !!
  !! OUTPUT:
  !! layout    - list of receiver type and type index
  !-----------------------------------------------------------------------
  SUBROUTINE initParallelism(in,layout)
    IMPLICIT NONE
    TYPE(SIMULATION_STRUCT), INTENT(IN) :: in
    TYPE(LAYOUT_STRUCT), INTENT(OUT) :: layout

    INTEGER :: i,j,k,ierr,remainder
    INTEGER :: cumulativeIndex,cumulativeVelocityIndex
    INTEGER :: rank,csize
    INTEGER :: nElements

    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,csize,ierr)

    ! list of number of elements in thread
    ALLOCATE(layout%listElements(csize), &
             layout%listOffset(csize),STAT=ierr)
    IF (ierr>0) STOP "could not allocate the list"

    ! total number of elements
    nElements=in%rectangularPatch%ns+in%triangularPatch%ns

    remainder=nElements-INT(nElements/csize)*csize
    IF (0 .LT. remainder) THEN
       layout%listElements(1:(csize-remainder))      =INT(nElements/csize)
       layout%listElements((csize-remainder+1):csize)=INT(nElements/csize)+1
    ELSE
       layout%listElements(1:csize)=INT(nElements/csize)
    END IF

    ! element start index in thread
    j=0
    k=0
    DO i=1,csize
       j=k+1
       k=k+layout%listElements(i)
       layout%listOffset(i)=j
    END DO

    ALLOCATE(layout%elementType       (layout%listElements(1+rank)), &
             layout%elementIndex      (layout%listElements(1+rank)), &
             layout%elementStateIndex (layout%listElements(1+rank)), &
             layout%elementVelocityDGF(layout%listElements(1+rank)), &
             layout%elementVelocityIndex(layout%listElements(1+rank)), &
             layout%elementStateDGF   (layout%listElements(1+rank)), &
             layout%elementForceDGF   (layout%listElements(1+rank)),STAT=ierr)
    IF (ierr>0) STOP "could not allocate the layout elements"

    j=1
    cumulativeIndex=0
    cumulativeVelocityIndex=0
    DO i=1,in%rectangularPatch%ns
       IF ((i .GE. layout%listOffset(1+rank)) .AND. &
           (i .LT. (layout%listOffset(1+rank)+layout%listElements(1+rank)))) THEN
          layout%elementType(j)=FLAG_RECTANGLE
          layout%elementIndex(j)=i
          layout%elementStateIndex(j)=cumulativeIndex+STATE_VECTOR_DGF_PATCH
          cumulativeIndex=layout%elementStateIndex(j)
          layout%elementVelocityDGF(j)=DGF_PATCH
          layout%elementVelocityIndex(j)=cumulativeVelocityIndex+DGF_PATCH
          cumulativeVelocityIndex=layout%elementVelocityIndex(j)
          layout%elementStateDGF(j)=in%dPatch
          layout%elementForceDGF(j)=DGF_VECTOR
          j=j+1
       END IF
    END DO
    DO i=1,in%triangularPatch%ns
       IF (((i+in%rectangularPatch%ns) .GE. layout%listOffset(1+rank)) .AND. &
           ((i+in%rectangularPatch%ns) .LT. (layout%listOffset(1+rank)+layout%listElements(1+rank)))) THEN
          layout%elementType(j)=FLAG_TRIANGLE
          layout%elementIndex(j)=i
          layout%elementStateIndex(j)=cumulativeIndex+STATE_VECTOR_DGF_PATCH
          cumulativeIndex=layout%elementStateIndex(j)
          layout%elementVelocityDGF(j)=DGF_PATCH
          layout%elementVelocityIndex(j)=cumulativeVelocityIndex+DGF_PATCH
          cumulativeVelocityIndex=layout%elementVelocityIndex(j)
          layout%elementStateDGF(j)=in%dPatch
          layout%elementForceDGF(j)=DGF_VECTOR
          j=j+1
       END IF
    END DO

    ! for MPI_ALLGATHERV
    ALLOCATE(layout%listVelocityN(csize), &
             layout%listVelocityOffset(csize), &
             layout%listStateN(csize), &
             layout%listStateOffset(csize), &
             layout%listForceN(csize), &
             STAT=ierr)
    IF (ierr>0) STOP "could not allocate the size list"

    ! share number of elements in threads
    CALL MPI_ALLGATHER(SUM(layout%elementVelocityDGF),1,MPI_INTEGER,layout%listVelocityN,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLGATHER(SUM(layout%elementStateDGF),   1,MPI_INTEGER,layout%listStateN,   1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    CALL MPI_ALLGATHER(SUM(layout%elementForceDGF),   1,MPI_INTEGER,layout%listForceN,   1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

    j=0
    k=0
    DO i=1,csize
       j=k+1
       k=k+layout%listVelocityN(i)
       layout%listVelocityOffset(i)=j-1
    END DO

    j=0
    k=0
    DO i=1,csize
       j=k+1
       k=k+layout%listStateN(i)
       layout%listStateOffset(i)=j-1
    END DO

  END SUBROUTINE initParallelism


  !-----------------------------------------------------------------------
  !> subroutine initGeometry
  ! initializes the position and local reference system vectors
  !
  ! INPUT:
  ! @param in      - input parameters data structure
  !-----------------------------------------------------------------------
  SUBROUTINE initGeometry(in)
    USE stuart97
    USE okada92
    USE types
    TYPE(SIMULATION_STRUCT), INTENT(INOUT) :: in

    IF (0 .LT. in%rectangularPatch%ns) THEN
       CALL computeReferenceSystemOkada92( &
                in%rectangularPatch%ns, &
                in%rectangularPatch%x, &
                in%rectangularPatch%length, &
                in%rectangularPatch%width, &
                in%rectangularPatch%strike, &
                in%rectangularPatch%dip, &
                in%rectangularPatch%sv, &
                in%rectangularPatch%dv, &
                in%rectangularPatch%nv, &
                in%rectangularPatch%xc)
    END IF

    IF (0 .LT. in%triangularPatch%ns) THEN
       CALL computeReferenceSystemStuart97( &
                in%triangularPatch%nVe, &
                in%triangularPatch%v, &
                in%triangularPatch%ns, &
                in%triangularPatch%i1, &
                in%triangularPatch%i2, &
                in%triangularPatch%i3, &
                in%triangularPatch%sv, &
                in%triangularPatch%dv, &
                in%triangularPatch%nv, &
                in%triangularPatch%xc)
    END IF

  END SUBROUTINE initGeometry

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
    TYPE(OPTION_S) :: opts(6)

    INTEGER :: i,k,rank,size,ierr,position
    INTEGER :: maxVertexIndex=-1
    INTEGER, PARAMETER :: psize=1024
    INTEGER :: dummy
    CHARACTER, DIMENSION(psize) :: packed

    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)

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

    in%nPatch=0
    ! minimum number of dynamic variables for patches
    in%dPatch=STATE_VECTOR_DGF_PATCH
    in%rectangularPatch%ns=0
    in%triangularPatch%ns=0
    in%triangularPatch%nVe=0
    in%nVolume=0
    in%dVolume=STATE_VECTOR_DGF_VOLUME
    in%strainVolume%ns=0

    IF (0.eq.rank) THEN
       PRINT 2000
       PRINT '("# RATESTATE")'
       PRINT '("# quasi-dynamic earthquake simulation in three-dimensional media")'
       PRINT '("# numerical accuracy: ",ES11.4)', epsilon
       PRINT '("# maximum iterations: ",I11)', maximumIterations
       PRINT '("# maximum time step: ",ES12.4)', maximumTimeStep
!$     PRINT '("#     * parallel OpenMP implementation with ",I3.3,"/",I3.3," threads")', &
!$                  omp_get_max_threads(),omp_get_num_procs()
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
                   in%rectangularPatch%nv(3,in%rectangularPatch%ns),STAT=ierr)
          IF (ierr>0) STOP "could not allocate the rectangular patch list"
          PRINT 2000
          PRINT '("#    n      Vpl       x1       x2       x3  length   width strike   dip   rake")'
          PRINT 2000
          DO k=1,in%rectangularPatch%ns
             CALL getdata(iunit,dataline)
             READ (dataline,*,IOSTAT=ierr) i, &
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
   
             PRINT '(I6,4ES9.2E1,2ES8.2E1,f7.1,f6.1,f7.1)',i, &
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
             WRITE (STDERR,'(a)') TRIM(dataline)
             WRITE(STDERR,'("input error: all rectangular patches require frictional properties")')
             STOP -1
          END IF
          PRINT '("#    n      mu0      sig        a        b        L       Vo  G/(2Vs)")'
          PRINT 2000
          DO k=1,in%rectangularPatch%ns
             CALL getdata(iunit,dataline)
             READ (dataline,*,IOSTAT=ierr) i, &
                   in%rectangularPatch%s(k)%mu0, &
                   in%rectangularPatch%s(k)%sig, &
                   in%rectangularPatch%s(k)%a, &
                   in%rectangularPatch%s(k)%b, &
                   in%rectangularPatch%s(k)%L, &
                   in%rectangularPatch%s(k)%Vo, &
                   in%rectangularPatch%s(k)%damping
   
             PRINT '(I6,7ES9.2E1)',i, &
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
                   in%triangularPatch%xc(3,in%triangularPatch%ns),STAT=ierr)
          IF (ierr>0) STOP "could not allocate the triangular patch list"
          PRINT 2000
          PRINT '("#    n      Vpl        i1        i2        i3   rake")'
          PRINT 2000
          DO k=1,in%triangularPatch%ns
             CALL getdata(iunit,dataline)
             READ (dataline,*,IOSTAT=ierr) i, &
                  in%triangularPatch%s(k)%Vpl, &
                  in%triangularPatch%i1(k), &
                  in%triangularPatch%i2(k), &
                  in%triangularPatch%i3(k), &
                  in%triangularPatch%s(k)%rake
             in%triangularPatch%s(k)%opening=0
   
             PRINT '(I6,ES9.2E1,3I10,f7.1)',i, &
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
             WRITE (STDERR,'(a)') TRIM(dataline)
             WRITE (STDERR,'("error in input file: not enough vertices")')
             STOP 1
          END IF
          IF (in%triangularPatch%nVe .GT. 0) THEN
             ALLOCATE(in%triangularPatch%v(3,in%triangularPatch%nVe),STAT=ierr)
             IF (ierr>0) STOP "could not allocate the list of vertices"
             PRINT 2000
             PRINT '("# n      x1       x2       x3")'
             PRINT 2000
             DO k=1,in%triangularPatch%nVe
                CALL getdata(iunit,dataline)
                READ (dataline,*,IOSTAT=ierr) i, &
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
             WRITE (STDERR,'(a)') TRIM(dataline)
             WRITE(STDERR,'("input error: all triangular patches require frictional properties")')
             STOP -1
          END IF

          PRINT '("#    n      mu0      sig        a        b        L       Vo    2G/Vs")'
          PRINT 2000
          DO k=1,in%triangularPatch%ns
             CALL getdata(iunit,dataline)
             READ (dataline,*,IOSTAT=ierr) i, &
                   in%triangularPatch%s(k)%mu0, &
                   in%triangularPatch%s(k)%sig, &
                   in%triangularPatch%s(k)%a, &
                   in%triangularPatch%s(k)%b, &
                   in%triangularPatch%s(k)%L, &
                   in%triangularPatch%s(k)%Vo, &
                   in%triangularPatch%s(k)%damping
   
             PRINT '(I6,7ES9.2E1)',i, &
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
       !       O B S E R V A T I O N   P A T C H E S
       ! - - - - - - - - - - - - - - - - - - - - - - - - - -
       PRINT 2000
       PRINT '("# number of observation patches")'
       CALL getdata(iunit,dataline)
       READ  (dataline,*) in%nObservationState
       PRINT '(I5)', in%nObservationState
       IF (0 .LT. in%nObservationState) THEN
          ALLOCATE(in%observationState(in%nObservationState,2),STAT=ierr)
          IF (ierr>0) STOP "could not allocate the observation patches"
          PRINT 2000
          PRINT '("# n      i")'
          PRINT 2000
          DO k=1,in%nObservationState
             CALL getdata(iunit,dataline)
             READ (dataline,*,IOSTAT=ierr) i,in%observationState(k,1)
             PRINT '(I3.3,X,I6)',i,in%observationState(k,1)
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
       IF (in%ne .GT. 0) ALLOCATE(in%event(in%ne),STAT=ierr)
       IF (ierr>0) STOP "could not allocate the event list"
       
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
           (in%triangularPatch%ns .EQ. 0) .OR. &
           (in%interval .LE. 0._8)) THEN
   
          WRITE_DEBUG_INFO(300)
          WRITE (STDERR,'("nothing to do. exiting.")')
          STOP 1
       END IF
   
       PRINT 2000
       ! flush standard output
       CALL FLUSH(6)      

       in%nPatch=in%rectangularPatch%ns+in%triangularPatch%ns

       position=0
       CALL MPI_PACK(in%interval,           1,MPI_REAL8,  packed,psize,position,MPI_COMM_WORLD,ierr)
       CALL MPI_PACK(in%lambda,             1,MPI_REAL8,  packed,psize,position,MPI_COMM_WORLD,ierr)
       CALL MPI_PACK(in%mu,                 1,MPI_REAL8,  packed,psize,position,MPI_COMM_WORLD,ierr)
       CALL MPI_PACK(in%nu,                 1,MPI_REAL8,  packed,psize,position,MPI_COMM_WORLD,ierr)
       CALL MPI_PACK(in%rectangularPatch%ns,1,MPI_INTEGER,packed,psize,position,MPI_COMM_WORLD,ierr)
       CALL MPI_PACK(in%triangularPatch%ns, 1,MPI_INTEGER,packed,psize,position,MPI_COMM_WORLD,ierr)
       CALL MPI_PACK(in%triangularPatch%nVe,1,MPI_INTEGER,packed,psize,position,MPI_COMM_WORLD,ierr)
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

    ELSE ! if 0.ne.rank

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
                             STAT=ierr)
       IF (ierr>0) STOP "slave could not allocate memory"

       IF (0 .LT. in%triangularPatch%ns) &
                    ALLOCATE(in%triangularPatch%s(in%triangularPatch%ns), &
                             in%triangularPatch%xc(3,in%triangularPatch%ns), &
                             in%triangularPatch%i1(in%triangularPatch%ns), &
                             in%triangularPatch%i2(in%triangularPatch%ns), &
                             in%triangularPatch%i3(in%triangularPatch%ns), &
                             in%triangularPatch%sv(3,in%triangularPatch%ns), &
                             in%triangularPatch%dv(3,in%triangularPatch%ns), &
                             in%triangularPatch%nv(3,in%triangularPatch%ns),STAT=ierr)
       IF (ierr>0) STOP "slave could not allocate memory"

       IF (0 .LT. in%triangularPatch%nVe) ALLOCATE(in%triangularPatch%v(3,in%triangularPatch%nVe), STAT=ierr)
       IF (ierr>0) STOP "slave could not allocate memory"

       IF (0 .LT. in%ne) ALLOCATE(in%event(in%ne),STAT=ierr)
       IF (ierr>0) STOP "slave could not allocate memory"

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
       PRINT '("mpirun -n 2 unicycle-3d-ratestate [-h] [--dry-run] [--help] [--epsilon 1e-6] [filename]")'
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
       PRINT '("   simulates elasto-dynamics on faults in three dimensions")'
       PRINT '("   following the radiation-damping approximation")'
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
       PRINT '("unicycle-3d-ratestate version 1.0.0, compiled on ",a)', __DATE__
       PRINT '("")'
       CALL FLUSH(6)
    END IF

  END SUBROUTINE printversion

END PROGRAM ratestate




