!-----------------------------------------------------------------------
!> program Unicycle (unified cycle of earthquakes) simulates evolution
!! of fault slip and distributed strain using the integral method under
!! the radiation damping approximation.
!!
!! \mainpage
!! 
!! The Green's function for traction and stress interaction amongst 
!! triangular and rectangular dislocations, and strain volumes has the 
!! following layout
!!
!!       / RR RT  RL \
!!       |           |
!!   G = | TR TT  TL |
!!       |           |
!!       \ LR LT  LL /
!!
!! where
!!
!!   RR is the matrix for traction on rectangular faults due to rectangular fault slip
!!   RT is the matrix for traction on rectangular faults due to triangular  fault slip
!!   TR is the matrix for traction on triangular  faults due to rectangular fault slip
!!   TT is the matrix for traction on triangular  faults due to triangular  fault slip
!!
!!   RL is the matrix for traction on rectangular faults due to volume strain
!!   TL is the matrix for traction on triangular  faults due to volume strain
!!
!!   LR is the matrix for stress in volumes due to rectangular fault slip
!!   LT is the matrix for stress in volumes due to triangular  fault slip
!!   LL is the matrix for stress in volumes due to strain in volumes
!!
!! The traction and stress interaction matrix can be thought as
!!
!!       / KK  KL \
!!   G = |        |
!!       \ LK  LL /
!!
!! with
!!
!!        / RR  RT \         / RL \
!!   KK = |        |,   KL = |    |,   LK = | LR  LT |
!!        \ TR  TT /         \ TL /         
!!
!! where KK is the traction due to fault slip,
!!       KL is the traction due to volume strain
!!       LK is the stress   due to fault slip
!!       LL is the stress   due to volume strain
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
!! The volumes have six strain directions: e11, e12, e13, e22, e23, e33.
!! The interaction matrices become
!!
!!        / LL1111  LL1112  LL1113  LL1122  LL1123  LL1133 \
!!        |                                                |
!!        | LL1211  LL1212  LL1213  LL1222  LL1223  LL1233 |
!!        |                                                |
!!        | LL1311  LL1312  LL1313  LL1322  LL1323  LL1333 |
!!   LL = |                                                |
!!        | LL2211  LL2212  LL2213  LL2222  LL2223  LL2233 |
!!        |                                                |
!!        | LL2311  LL2312  LL2313  LL2322  LL2323  LL2333 |
!!        |                                                |
!!        \ LL3311  LL3312  LL3313  LL3322  LL3323  LL3333 /
!!
!!        / RLs11  RLs12  RLs13  RLs22  RLs23  RLs33 \
!!   RL = |                                          |
!!        \ RLd11  RLd12  RLd13  RLd22  RLd23  RLd33 /
!!
!!        / TLs11  TLs12  TLs13  TLs22  TLs23  TLs33 \
!!   TL = |                                          |
!!        \ TLd11  TLd12  TLd13  TLd22  TLd23  TLd33 /
!!
!!        / LR11s  LR11d \            / LT11s  LT11d \
!!        |              |            |              |
!!        | LR12s  LR12d |            | LT12s  LT12d |
!!        |              |            |              |
!!        | LR13s  LR13d |            | LT13s  LT13d |
!!   LR = |              |,      LT = |              |
!!        | LR22s  LR22d |            | LT22s  LT22d |
!!        |              |            |              |
!!        | LR23s  LR23d |            | LT23s  LT23d |
!!        |              |            |              |
!!        \ LR33s  LR33d /            \ LT33s  LT33d /
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
!!    | Tn dPatch  |   +-----------------------+  +-----------------+
!!    |            |                           
!!    | V1 1       |   +-----------------------+  +-----------------+
!!    | .          |   |                       |  |                 |
!!    | V1 dVolume |   |                       |  |     nVolume     |
!!    | .          |   |   nVolume * dVolume   |  |        *        |
!!    | .          |   |                       |  |     dVolume     |
!!    | Vn 1       |   |                       |  |                 |
!!    | .          |   |                       |  |                 |
!!    \ Vn dVolume /   +-----------------------+  +-----------------+
!!
!! where nRectangle, nTriangle, and nVolume are the number of rectangular
!! and triangular patches and the number of strain volumes, respectively,
!! and dPatch and dVolume are the degrees of freedom for patches and 
!! volumes. 
!!
!! For every rectangular or triangular patch, we have the following 
!! items in the state vector
!!
!!   /  ts   \  1
!!   |  td   |  .
!!   |  tn   |  .
!!   |  ss   |  .
!!   |  sd   |  .
!!   ! theta |  .
!!   \  v*   /  dPatch
!!
!! where ts, td, and td are the local traction in the strike, dip, and
!! normal directions, ss and sd are the total slip in the strike and dip 
!! directions, v* is the logarithm of the norm of the instantaneous slip 
!! velocity vector (v*=ln(V/Vo)), and theta is the state variable in the
!! rate and state friction framework.
!!
!! For every strain volume, we have the following items in the state vector
!!
!!   /  s11   \  1
!!   |  s12   |  .
!!   |  s13   |  . 
!!   |  s22   |  .
!!   |  s23   |  .
!!   |  s33   |  .
!!   |  e11   |  .
!!   |  e12   |  .
!!   |  e13   |  .
!!   |  e22   |  .
!!   |  e23   |  .
!!   |  e33   |  .
!!   \   d*  /   dVolume
!!
!! where s11, s12, s13, s22, s23, and s33 are the six independent components
!! of the local stress tensor, e11, e12, e13, e22, e23, and e33 are the six
!! independent components of the cumulative anelastic strain tensor, and d*
!! is the logarithm of the grain size (d* = ln(d/do)).
!!
!! References:<br>
!!
!!   James D P Moore, Sylvain Barbot, Lujia Feng, Yu Hang, Valere Lambert, 
!!   Eric Lindsey, Sagar Masuti, Takanori Matsuzawa, Jun Muto, 
!!   Priyamvada Nanjundiah, Rino Salman, Sharadha Sathiakumar, 
!!   and Harpreet Sethi. (2019, September 25). 
!!   jdpmoore/unicycle: Unicycle (Version 1.0). Zenodo.
!!   http://doi.org/10.5281/zenodo.4471162
!!
!!   Moore, J. DP, H. Yu, C.-H. Tang, T. Wang, S. Barbot, D. Peng, S. Masuti, 
!!   J. Dauwels, Y.-J. Hsu, V. Lambert, P. Nanjundiah, S. Wei, E. Lindsey, 
!!   L. Feng and B. Shibazaki. "Imaging the distribution of transient 
!!   viscosity following the 2016 Mw 7.1 Kumamoto earthquake", Science (2017).
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
!!   Feng, Lujia, et al. "Footprints of past earthquakes revealed in the 
!!   afterslip of the 2010 Mw 7.8 Mentawai tsunami earthquake." Geophysical 
!!   Research Letters 43.18 (2016): 9518-9526.
!!
!!   Qiu, Q., E. M. Hill, S. Barbot, J. Hubbard, W. Feng, E. O. Lindsey, 
!!   L. Feng, K. Dai, S. V. Samsonov and P. Tapponnier "The mechanism of 
!!   partial rupture of a locked megathrust: The role of fault morphology." 
!!   Geology 44.10 (2016): 875-878.
!!
!! \author James D P Moore (2016)
!! \author Sylvain Barbot (2017).
!----------------------------------------------------------------------
PROGRAM unicycle

#include "macros.f90"

  USE greens
  USE input 
  USE ode45
  USE rheology
  USE strainvolume
  USE types 

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  TYPE(SIMULATION_STRUCT) :: iparam

  REAL*8, PARAMETER :: pi=3.141592653589793d0

  INTEGER :: ierr,rank,csize
  ! output directory, grd files titles and names.
  CHARACTER(80) :: wdir,title,filename

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

  ! Green's function
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: G

  ! time
  REAL*8 :: time,t0
  ! time step
  REAL*8 :: dt_try,dt_next,dt_done
  ! time steps
  INTEGER :: i,j

  ! maximum number of time steps
  INTEGER, PARAMETER :: maxStep=1000000

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

  ! calculate basis vectors
  CALL initGeometry(in)

  ! describe data layout for parallelism
  CALL initParallelism(in,layout)

  ! calculate the stress interaction matrix
  IF (0 .EQ. rank) THEN
     PRINT '("# computing Green''s functions.")'
  END IF
  CALL buildG(in,layout,G)
  IF (0 .EQ. rank) THEN
     PRINT 2000
  END IF

  ! velocity vector and strain rate tensor array (t=G*vAll)
  ALLOCATE(v(layout%listVelocityN(1+rank)), &
           vAll(SUM(layout%listVelocityN)),STAT=ierr)
  IF (ierr>0) STOP "could not allocate the velocity and strain rate vector"

  ! traction vector and stress tensor array (t=G*v)
  ALLOCATE(t(layout%listForceN(1+rank)),STAT=ierr)
  IF (ierr>0) STOP "could not allocate the traction and stress vector"

  ALLOCATE(y(layout%listStateN(1+rank)),STAT=ierr)
  IF (ierr>0) STOP "could not allocate the state vector"

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
     WRITE(STDOUT,'("#       n               time                 dt")')
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

  ! main loop
  DO i=1,maxStep

     CALL odefun(layout%listStateN(1+rank),time,y,dydt)

     CALL export()

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

  DEALLOCATE(y,dydt,yscal)
  DEALLOCATE(ytmp,ytmp1,ytmp2,ytmp3,yrkck)
  DEALLOCATE(AK2,AK3,AK4,AK5,AK6)
  DEALLOCATE(G,v,vAll,t)
  DEALLOCATE(layout%listForceN)
  DEALLOCATE(layout%listVelocityN,layout%listVelocityOffset)
  DEALLOCATE(layout%listStateN,layout%listStateOffset)
  DEALLOCATE(layout%elementStateIndex)
  DEALLOCATE(layout%listElements,layout%listOffset)

  CALL MPI_FINALIZE(ierr)

2000 FORMAT ("# ----------------------------------------------------------------------------")
     
CONTAINS

  !-----------------------------------------------------------------------
  !> subroutine export
  ! write the state variables of elements, either patch or volume, and 
  ! other information.
  !----------------------------------------------------------------------
  SUBROUTINE export()

    INTEGER :: dgf
    INTEGER :: j,stateIndex

    ! export observation state
    DO j=1,in%nObservationState
       IF ((in%observationState(j,1) .GE. layout%listOffset(rank+1)) .AND. &
           (in%observationState(j,1) .LT. layout%listOffset(rank+1)+layout%listElements(rank+1))) THEN

          SELECT CASE(layout%elementType(in%observationState(j,1)-layout%listOffset(rank+1)+1))
          CASE (FLAG_RECTANGLE,FLAG_TRIANGLE)
             dgf=in%dPatch
          CASE (FLAG_VOLUME)
             dgf=in%dVolume
          CASE DEFAULT
             WRITE(STDERR,'("wrong case: this is a bug.")')
             WRITE_DEBUG_INFO(-1)
             STOP -1
          END SELECT

          stateIndex=layout%elementStateIndex(in%observationState(j,1)-layout%listOffset(rank+1)+1)-dgf
          WRITE (in%observationState(j,2),*) time, &
                    y(stateIndex+1:stateIndex+dgf), &
                 dydt(stateIndex+1:stateIndex+dgf)
       END IF
    END DO

    IF (0 .EQ. rank) THEN
       WRITE(FPTIME,'(ES19.12E2,ES19.12E2)') time,dt_done
       IF (0 .EQ. MOD(i,50)) THEN
          WRITE(STDOUT,'(I9.9,ES19.12E2,ES19.12E2)') i,time,dt_done
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

    INTEGER :: i,l,ierr
    INTEGER :: elementType,elementIndex
    TYPE(PATCH_ELEMENT_STRUCT) :: patch
    TYPE(STRAINVOLUME_ELEMENT_STRUCT) :: volume

    REAL*8 :: eps,tau

    ! zero out state vector
    y=0._8

    ! element index in state vector
    l=1
    ! loop over elements owned by current thread
    DO i=1,SIZE(layout%elementIndex)
       elementType= layout%elementType(i)
       elementIndex=layout%elementIndex(i)

       SELECT CASE (elementType)
       CASE (FLAG_RECTANGLE,FLAG_TRIANGLE)

          IF (elementType .EQ. FLAG_RECTANGLE) THEN
             patch=in%rectangularPatch%s(elementIndex)
          ELSE
             patch=in%triangularPatch%s(elementIndex)
          END IF

          ! strike slip
          y(l+STATE_VECTOR_SLIP_STRIKE) = 0._8

          ! dip slip
          y(l+STATE_VECTOR_SLIP_DIP) = 0._8

          ! traction
          tau = patch%sig*(patch%mu0+(patch%a-patch%b)*log(patch%Vpl/patch%Vo))+patch%damping*patch%Vpl

          ! traction in strike direction
          y(l+STATE_VECTOR_TRACTION_STRIKE) = tau*cos(patch%rake)

          ! traction in dip direction
          y(l+STATE_VECTOR_TRACTION_DIP) = tau*sin(patch%rake)

          ! traction in normal direction
          y(l+STATE_VECTOR_TRACTION_NORMAL) = 0._8

          ! state variable log(theta Vo / L)
          y(l+STATE_VECTOR_STATE_1) = log(patch%Vo/patch%Vpl)

          ! slip velocity in strike direction
          y(l+STATE_VECTOR_VELOCITY_STRIKE) = 0.98d0*patch%Vpl*cos(patch%rake)

          ! slip velocity in dip direction
          y(l+STATE_VECTOR_VELOCITY_DIP) = 0.98d0*patch%Vpl*sin(patch%rake)

          l=l+in%dPatch

       CASE (FLAG_VOLUME)

          volume=in%strainVolume%s(elementIndex)

          ! norm of loading rate
          eps=SQRT((volume%e11**2+ &
                  2*volume%e12**2+ &
                  2*volume%e13**2+ &
                    volume%e22**2+ &
                  2*volume%e23**2+ &
                    volume%e33**2)/2._8)

          ! stress from nonlinear viscosity in Maxwell element
          IF (0._8 .LT. eps) THEN
             tau=volume%nGm*(eps/volume%ngammadot0m*exp(volume%nQm/(volume%nRm*volume%To)))**(1._8/volume%npowerm)

             ! stress components
             y(l+STATE_VECTOR_S11) = tau*volume%e11/eps
             y(l+STATE_VECTOR_S12) = tau*volume%e12/eps
             y(l+STATE_VECTOR_S13) = tau*volume%e13/eps
             y(l+STATE_VECTOR_S22) = tau*volume%e22/eps
             y(l+STATE_VECTOR_S23) = tau*volume%e23/eps
             y(l+STATE_VECTOR_S33) = tau*volume%e33/eps
          ELSE
             ! stress components
             y(l+STATE_VECTOR_S11) = 0._8
             y(l+STATE_VECTOR_S12) = 0._8
             y(l+STATE_VECTOR_S13) = 0._8
             y(l+STATE_VECTOR_S22) = 0._8
             y(l+STATE_VECTOR_S23) = 0._8
             y(l+STATE_VECTOR_S33) = 0._8
          END IF

          ! strain components
          y(l+STATE_VECTOR_E11) = 0._8
          y(l+STATE_VECTOR_E12) = 0._8
          y(l+STATE_VECTOR_E13) = 0._8
          y(l+STATE_VECTOR_E22) = 0._8
          y(l+STATE_VECTOR_E23) = 0._8
          y(l+STATE_VECTOR_E33) = 0._8

          l=l+in%dVolume

       CASE DEFAULT
          WRITE(STDERR,'("wrong case: this is a bug.")')
          WRITE_DEBUG_INFO(-1)
          STOP -1
       END SELECT

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
    TYPE(STRAINVOLUME_ELEMENT_STRUCT) :: volume

    ! traction components in the strike and dip directions
    REAL*8 :: ts, td

    ! scalar rate of shear traction
    REAL*8 :: dtau

    ! velocity and acceleration scalars
    REAL*8 :: velocity,acceleration

    ! slip velocity in the strike and dip directions
    REAL*8 :: vs,vd

    ! independent components of the stress tensor
    REAL*8 :: s11,s12,s13,s22,s23,s33

    ! norm of deviatoric stress
    REAL*8 :: sII

    ! pressure
    REAL*8 :: p

    ! independent components of the strain rate tensor
    REAL*8 :: e11,e12,e13,e22,e23,e33

    ! scalar strain rate
    REAL*8 :: eII

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
       CASE (FLAG_RECTANGLE,FLAG_TRIANGLE)
          ! v(k:k+layout%elementVelocityDGF(i)-1) = slip velocity

          IF (elementType .EQ. FLAG_RECTANGLE) THEN
             patch=in%rectangularPatch%s(elementIndex)
          ELSE
             patch=in%triangularPatch%s(elementIndex)
          END IF

          ! slip velocity
          vs=y(l+STATE_VECTOR_VELOCITY_STRIKE)
          vd=y(l+STATE_VECTOR_VELOCITY_DIP)

          ! update state vector (rate of slip components)
          dydt(l+STATE_VECTOR_SLIP_STRIKE)=vs
          dydt(l+STATE_VECTOR_SLIP_DIP   )=vd

          v(k:k+layout%elementVelocityDGF(i)-1)=(/ &
                        vs-patch%Vpl*COS(patch%rake), &
                        vd-patch%Vpl*SIN(patch%rake) /)

          l=l+in%dPatch

       CASE (FLAG_VOLUME)
          ! v(k:k+layout%elementVelocityDGF(i)-1)= strain rate

          volume=in%strainVolume%s(elementIndex)

          ! stress components
          s11=y(l+STATE_VECTOR_S11)
          s12=y(l+STATE_VECTOR_S12)
          s13=y(l+STATE_VECTOR_S13)
          s22=y(l+STATE_VECTOR_S22)
          s23=y(l+STATE_VECTOR_S23)
          s33=y(l+STATE_VECTOR_S33)

          ! pressure
          p=(s11+s22+s33)/3._8

          ! deviatoric stress components
          s11=s11-p
          s22=s22-p
          s33=s33-p

          ! initiate strain rate
          e11=0._8
          e12=0._8
          e13=0._8
          e22=0._8
          e23=0._8
          e33=0._8

          ! strain rate from linear viscosity in Maxwell element
          IF (0 .LT. in%strainVolume%nMaxwell) THEN
             e11=e11+volume%gammadot0m*(s11/in%mu)
             e12=e12+volume%gammadot0m*(s12/in%mu)
             e13=e13+volume%gammadot0m*(s13/in%mu)
             e22=e22+volume%gammadot0m*(s22/in%mu)
             e23=e23+volume%gammadot0m*(s23/in%mu)
             e33=e33+volume%gammadot0m*(s33/in%mu)
          END IF

          ! strain rate from nonlinear viscosity in Maxwell element
          IF (0 .LT. in%strainVolume%nNonlinearMaxwell) THEN
             sII=SQRT((s11**2+2._8*s12**2+2._8*s13**2+ &
                       s22**2+2._8*s23**2+     s33**2)/2._8)

             eII=volume%ngammadot0m*(sII/in%mu)**(volume%npowerm-1)/in%mu

             e11=e11+eII*s11
             e12=e12+eII*s12
             e13=e13+eII*s13
             e22=e22+eII*s22
             e23=e23+eII*s23
             e33=e33+eII*s33
          END IF

          ! strain rate tensor
          dydt(l+STATE_VECTOR_E11)=e11
          dydt(l+STATE_VECTOR_E12)=e12
          dydt(l+STATE_VECTOR_E13)=e13
          dydt(l+STATE_VECTOR_E22)=e22
          dydt(l+STATE_VECTOR_E23)=e23
          dydt(l+STATE_VECTOR_E33)=e33

          v(k:k+layout%elementVelocityDGF(i)-1)=(/ &
                  e11-volume%e11, &
                  e12-volume%e12, &
                  e13-volume%e13, &
                  e22-volume%e22, &
                  e23-volume%e23, &
                  e33-volume%e33 &
                  /)

          l=l+in%dVolume

       CASE DEFAULT
          WRITE(STDERR,'("wrong case: this is a bug.")')
          WRITE_DEBUG_INFO(-1)
          STOP -1
       END SELECT

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
       CASE (FLAG_RECTANGLE,FLAG_TRIANGLE)

          IF (elementType .EQ. FLAG_RECTANGLE) THEN
             patch=in%rectangularPatch%s(elementIndex)
          ELSE
             patch=in%triangularPatch%s(elementIndex)
          END IF

          ! slip velocity
          vs=y(l+STATE_VECTOR_VELOCITY_STRIKE)
          vd=y(l+STATE_VECTOR_VELOCITY_DIP)

          velocity=SQRT(vs**2+vd**2)

          ! rate of state
          dydt(l+STATE_VECTOR_STATE_1)=(patch%Vo*exp(-y(l+STATE_VECTOR_STATE_1))-velocity) / patch%L

          ! scalar rate of shear traction
          dtau=dydt(l+STATE_VECTOR_TRACTION_STRIKE)*COS(patch%rake) &
              +dydt(l+STATE_VECTOR_TRACTION_DIP)   *SIN(patch%rake)
          
          ! acceleration
          acceleration=((dtau-patch%b*patch%sig*dydt(l+STATE_VECTOR_STATE_1)) / &
                        (patch%a*patch%sig+patch%damping*velocity))*velocity

          ! acceleration vector
          dydt(l+STATE_VECTOR_VELOCITY_STRIKE)=acceleration*COS(patch%rake)
          dydt(l+STATE_VECTOR_VELOCITY_DIP   )=acceleration*SIN(patch%rake)

          ! traction rate
          dydt(l+STATE_VECTOR_TRACTION_STRIKE)=t(j+TRACTION_VECTOR_STRIKE)-patch%damping*dydt(l+STATE_VECTOR_VELOCITY_STRIKE)
          dydt(l+STATE_VECTOR_TRACTION_DIP   )=t(j+TRACTION_VECTOR_DIP   )-patch%damping*dydt(l+STATE_VECTOR_VELOCITY_DIP)
          dydt(l+STATE_VECTOR_TRACTION_NORMAL)=t(j+TRACTION_VECTOR_NORMAL)

          l=l+in%dPatch

       CASE (FLAG_VOLUME)

          ! return the stress rate
          dydt(l+STATE_VECTOR_S11)=t(j+TRACTION_VECTOR_S11)
          dydt(l+STATE_VECTOR_S12)=t(j+TRACTION_VECTOR_S12)
          dydt(l+STATE_VECTOR_S13)=t(j+TRACTION_VECTOR_S13)
          dydt(l+STATE_VECTOR_S22)=t(j+TRACTION_VECTOR_S22)
          dydt(l+STATE_VECTOR_S23)=t(j+TRACTION_VECTOR_S23)
          dydt(l+STATE_VECTOR_S33)=t(j+TRACTION_VECTOR_S33)

          l=l+in%dVolume

       CASE DEFAULT
          WRITE(STDERR,'("wrong case: this is a bug.")')
          WRITE_DEBUG_INFO(-1)
          STOP -1
       END SELECT

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

    INTEGER :: i,j,k,ierr,n,remainder,cumulativeIndex
    INTEGER :: rank,csize
    INTEGER :: nElements,nColumns
    INTEGER :: buffer

    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,csize,ierr)

    ! total number of elements
    nElements=in%rectangularPatch%ns+in%triangularPatch%ns+in%strainVolume%ns

    ! list of number of elements in thread
    ALLOCATE(layout%listElements(csize), &
             layout%listOffset(csize),STAT=ierr)
    IF (ierr>0) STOP "could not allocate the list"

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
             layout%elementStateDGF   (layout%listElements(1+rank)), &
             layout%elementForceDGF   (layout%listElements(1+rank)),STAT=ierr)
    IF (ierr>0) STOP "could not allocate the layout elements"

    j=1
    cumulativeIndex=0
    DO i=1,in%rectangularPatch%ns
       IF ((i .GE. layout%listOffset(1+rank)) .AND. &
           (i .LT. (layout%listOffset(1+rank)+layout%listElements(1+rank)))) THEN
          layout%elementType(j)=FLAG_RECTANGLE
          layout%elementIndex(j)=i
          layout%elementStateIndex(j)=cumulativeIndex+STATE_VECTOR_DGF_PATCH
          cumulativeIndex=layout%elementStateIndex(j)
          layout%elementVelocityDGF(j)=DGF_PATCH
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
          layout%elementStateDGF(j)=in%dPatch
          layout%elementForceDGF(j)=DGF_VECTOR
          j=j+1
       END IF
    END DO
    DO i=1,in%strainVolume%ns
       IF (((i+in%rectangularPatch%ns+in%triangularPatch%ns) .GE. layout%listOffset(1+rank)) .AND. &
           ((i+in%rectangularPatch%ns+in%triangularPatch%ns) .LT. (layout%listOffset(1+rank)+layout%listElements(1+rank)))) THEN
          layout%elementType(j)=FLAG_VOLUME
          layout%elementIndex(j)=i
          layout%elementStateIndex(j)=cumulativeIndex+STATE_VECTOR_DGF_VOLUME
          cumulativeIndex=layout%elementStateIndex(j)
          layout%elementVelocityDGF(j)=DGF_VOLUME
          layout%elementStateDGF(j)=in%dVolume
          layout%elementForceDGF(j)=DGF_TENSOR
          j=j+1
       END IF
    END DO

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
    USE strainvolume
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

    IF (0 .LT. in%strainVolume%ns) THEN
       CALL computeReferenceSystemVerticalStrainVolume( &
                in%strainVolume%ns, &
                in%strainVolume%x, &
                in%strainVolume%length, &
                in%strainVolume%width, &
                in%strainVolume%strike, &
                in%strainVolume%sv, &
                in%strainVolume%dv, &
                in%strainVolume%nv, &
                in%strainVolume%xc)
    END IF

  END SUBROUTINE initGeometry

END PROGRAM unicycle




