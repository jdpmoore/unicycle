!
!  sse3d_mpi_blas.f90
!
!  2014-05-30 Version 0.1  by Takanori Matsuzawa
!  2014-07-08 Version 0.2  by Takanori Matsuzawa
!  2017-01-19 Version 0.3  by Takanori Matsuzawa (Adding comment for sharing this code with Sylvain Barbot)
!
module common_param
  integer, parameter :: MaxStep = 200000000  ! Number of maximum steps
  integer, parameter :: Nall  = 156682       ! Number of subfault elements
  integer, parameter :: NPrmDim = 2          ! Number of parameter in RS-law
  integer, parameter :: NbinoutPrmDim = NPrmDim + 1   ! Number of binary output parameters (Displacement is added as extra output)
  real(8), parameter :: eps = 1.d-6          ! For Runge-Kutta
  real(8), parameter :: TINY = 1.d-30        ! For Runge-Kutta
  integer, parameter :: root_rank0 = 0       ! Number of master thread rank in the original MPI
  integer, parameter :: root_rank  = root_rank0 + 1   ! Number of master thread rank (Our def. is shifted by one from the original value in MPI_COMM_RANK)
  integer, parameter :: Nbinout = 200        ! Output interval for binary output (time steps)
  integer, parameter :: Nascout = 20*Nbinout ! Output interval for ascii output (time steps)
  real(8), parameter :: v_asta = 1.d0        ! Reference slip Velocity [m/s] 
  real(8), parameter :: Ysec =  365.2422d0 * 24.d0 * 60.d0 * 60.d0  ! Second of 1-year
  real(8), parameter :: vcut1 = 1.d0         ! Cut off velocity v1 (faster cut off velocity [m/s])
  real(8), parameter :: rigid = 3.0d4        ! Rigidity [MPa]
  real(8), parameter :: beta = 3.5d3         ! S-wave_Velocity [m/s]
  real(8), parameter :: dc0  = 6.d-2         ! Reference length of Dc [m/s]
  real(8), parameter :: v_init = 1.d-11      ! Initial Velocity [m/s]
  real(8), parameter :: dt_init = 1.d0*v_asta/dc0  ! Initial dt (normalized) for Runge-Kutta 
  real(8), parameter :: DoutIntvYr = 0.1d0   ! Maximum output intervals [year]
  real(8), parameter :: fac1 = rigid/2.d0/beta*v_asta
  real(8), parameter :: fac2 = dc0/v_asta
  real(8), parameter :: theta2in = v_asta/vcut1

  real(8), parameter :: IsnapSW = 0    ! IsnapSW = 0 : start without a snapshot file; =1 : start with a snapshot file
  integer, parameter :: Nversion = 30  ! Version number of this code
  integer, parameter :: Nfricck  = 14  ! Number of output columns for frictional parameter check
  integer, parameter :: FPBO   = 21  ! Device number for binary output
  integer, parameter :: FPAO   = 22  ! Device number for ascii output
  integer, parameter :: FPSSO  = 23  ! Device number for snapshot output
  integer, parameter :: FPCKO  = 51  ! Device number for frictional parameter output
  integer, parameter :: FPRNK  = 52  ! Device number for rank list check output
  integer, parameter :: FPTRM  = 53  ! Device number for tremor list
  integer, parameter :: STDOUT = 6   ! Device number for standard output
  integer, parameter :: STDERR = 0   ! Device number for standard error output

  real(8), parameter :: pi = 3.141592653589793d0
end module common_param

!==============================================
module fric_param
  use common_param 
  real(8), parameter :: randfac = 0.1d0    ! Factor for fluctuation of parameter (b in RS-law)
  real(8), parameter :: vpl0 = 6.d-2/Ysec  ! Plate velocity: 6cm/yr -> m/s
  real(8), parameter :: vpl1 = 2.d-2/Ysec  ! Plate velocity: 2cm/yr -> m/s (at Suruga Trough)
  real(8), parameter :: xvpl0 = 5.2d5      ! X (m) at Kii
  real(8), parameter :: xvpl1 = 7.6d5      ! X (m) at Suruga

!! a, b in RS-law ( depth[m], a, b ) 
  integer, parameter :: n_ab = 7
  real(8), parameter :: zab_rs(n_ab)  = (/ 5.d0,  1.3d4,  2.d4,  2.2d4,  2.7d4,  3.1d4,  4.5d4  /)
  real(8), parameter :: a_rs(n_ab)    = (/ 8.d-3, 8.d-3, 8.d-3,  8.d-3,  8.d-3,  1.2d-2, 2.4d-2 /)
  real(8), parameter :: b_rs(n_ab)    = (/ 0.d0, 1.2d-2, 1.2d-2, 1.2d-2, 1.2d-2, 0.d0,   0.d0   /)
!! a, b in L-SSE region
  real(8), parameter :: zab_sse(n_ab) = (/ 5.d0,  1.3d4,  2.d4,  2.2d4,  3.4d4,  3.4000001d4,  4.5d4  /)

!! a, b in L-SSE region, and tremor (S-SSE) region
  integer, parameter :: n_ab_trm = 2
  real(8), parameter :: zab_trm(n_ab_trm)  = (/ 2.8d4,  4.5d4   /)
  real(8), parameter :: a_trm(n_ab_trm)    = (/ 8.d-3,  8.5d-3  /)
  real(8), parameter :: b_trm(n_ab_trm)    = (/ 1.2d-2, 1.05d-2 /)
   
!! dc in RS-law  ( depth[m], dc[m] ) 
  integer, parameter :: n_dc = 6
  real(8),parameter  :: zdc_rs(n_dc)  = (/ 0.d0,    2.2d4,   2.751020408d4, 2.751020409d4, 2.8d4, 5.2d4 /)
  real(8),parameter  :: dc_rs(n_dc)   = (/ 2.d-2,   2.d-2,   2.d-3,         2.d-3,         4.d-4, 4.d-4 /)
!! dc in L-SSE region, and boundary (edge) region
  real(8),parameter  :: zdc_sse(n_dc) = (/ 0.d0,    2.0d4,   2.2d4,         3.2d4,         3.4d4, 5.2d4 /)
  real(8),parameter  :: zdc_bnd(n_dc) = (/ 0.d0,    1.0d4,   1.5d4,         2.0d4,         2.8d4, 5.2d4 /)

!! Effective normal stress ( depth[m], s_eff[MPa] )
  integer, parameter :: n_seff = 6
  real(8), parameter :: zseff(n_seff)     = (/ 1.5d3,   2.2d4,   2.751020408d4, 2.751020409d4, 2.8d4,   4.5d4    /)
  real(8), parameter :: seff(n_seff)      = (/ 32.18d0, 32.18d0, 3.218d0,       3.218d0,      0.6436d0, 0.6436d0 /)
!! Effective normal stress in L-SSE region, and boundary (edge) region
  real(8),parameter  :: zseff_sse(n_seff) = (/ 1.5d0,   2.0d4,   2.2d4,         3.2d4,         3.4d4,   5.2d4    /)
  real(8),parameter  :: zseff_bnd(n_seff) = (/ 1.5d0,   1.0d4,   1.5d4,         2.0d4,         2.8d4,   5.2d4    /)

!! Cutoff Velocity 
  integer, parameter :: n_vc = 2
  real(8), parameter :: zvc_rs(n_vc) = (/ zseff(2), zseff(3) /)
  real(8), parameter :: vc_rs(n_vc)  = (/ 0.d0,     -6.5d0   /)

!! tremor (S-SSE) region
  real(8), parameter  :: ThTrmdist = 2.5d3   ! Threshold of Distance [m] from hbc tremor
  integer, parameter  :: NtrThresh = 2      ! Minimum number of tremors to set S-SSE region

!! Location of Long-term SSE patch
  real(8), parameter  :: xsse1 = 1.0d5
  real(8), parameter  :: xsse2 = 1.2d5
  real(8), parameter  :: xsse3 = 1.6d5
  real(8), parameter  :: xsse4 = 1.8d5

!! Location of Boundary region (Edge region) in calculation
  real(8), parameter  :: xbnd1  = 0.2d5
  real(8), parameter  :: xbnd2  = 0.6d5
  real(8), parameter  :: xbnd3  = 7.5d5
  real(8), parameter  :: xbnd4  = 7.9d5
  !! b-value transition zone
  real(8), parameter  :: bbnd  = 0.d0    ! b value at x < xb0 and z < zb1
  real(8), parameter  :: zbbnd = 3.1d4   ! bottom of transition zone for b
 
end module fric_param

!==============================================
module common_var
  use mpi
  use common_param
  integer :: my_rank0, my_rank  ! my_rank0: original thread number; my_rank = my_rank9 +1
  integer :: n_rank             ! n_rank: total thread number (size)
  integer :: nc                 ! number of elements in the current thread
  integer :: lstrnk(Nall) ,lstrnk0(Nall)
  real(8) :: vn(Nall)           ! vn: Normalized velocity (all)
  real(8) :: xp(3,Nall)
  
  integer, allocatable :: lstnc(:),lstdispl(:),lstdispl2(:),lstdispl0(:)
  integer, allocatable :: lstnc_binout(:),lstdispl_binout(:),lstdispl0_binout(:)
  real(8), allocatable :: gfn(:,:)
  real(8), allocatable :: AK2(:,:),AK3(:,:),AK4(:,:),AK5(:,:),AK6(:,:)   ! For Runge-Kutta

  real(8), allocatable :: as(:),bs(:),tauref(:),g1(:),g3(:)
  real(8), allocatable :: vnc(:),vncr(:),vnr(:),vnc_pre(:),dispnc(:),yout_nc(:,:),xp_c(:,:),vnpl_c(:) ! vnpl: Normalized plate velocity (all)
  real(8), allocatable :: yderivs(:),yrkck(:,:)

end module common_var

!==============================================
! Main Program
!==============================================
program sse3d_mpi
  use mpi
  use common_param
  use common_var
  implicit none

  external :: DERIVSm
  integer :: istep,ierr,ibinout,iascout,ncbinout,i
  integer :: itm0,itm1,itm2,ntmr,maxtm
  real(8) :: t,t0
  real(8) :: dt_try,dt_next,dt_done
  real(8) :: TmNextWrite
  real(8) :: yout_all(NbinoutPrmDim,Nall)
  real(8), allocatable :: y(:,:),dydx(:,:),yc(:,:),yscal(:,:)
  real(8), allocatable :: ytmp(:,:),ytmp1(:,:),ytmp2(:,:),ytmp3(:,:)
  real(4) :: etmres0,etmres1

!! Initialize MPI
  my_rank0 = 0
  my_rank = 1
  n_rank = 1
  nc = Nall
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_rank, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank0,ierr)
  my_rank = my_rank0 + 1
!! Generate rank index in MPI
  call calcRankIndex
  nc = lstnc(my_rank)   ! number of elements in the current thread
  ncbinout = nc*NbinoutPrmDim

  if ( my_rank0.eq.root_rank0 ) then
    call SYSTEM_CLOCK(itm0,ntmr,maxtm)
    call CPU_TIME(etmres0)
  endif

!! Memory allocations
  call commVarAlloc
  allocate(y(nc,NPrmDim),yc(nc,NprmDim),dydx(nc,NPrmDim),yscal(nc,NprmDim))
  allocate(ytmp(nc,NprmDim),ytmp1(nc,NprmDim),ytmp2(nc,NprmDim),ytmp3(nc,NprmDim))

!! Open Outputfiles
  if ( my_rank0.eq.root_rank0 )  then
    open(FPBO,file="out_sse3d.bin",access= "stream",form="unformatted")
    open(FPAO,file="out_sse3d.asc",access= "stream",form="formatted")
 
    ! Write a header part of binary output
    write(FPBO)Nall,Nversion
    write(FPBO)NPrmDim,NbinoutPrmDim,Nfricck
    write(FPBO)v_asta,vcut1,dc0,fac1,fac2
  endif

!! Read Green's Function
  if ( my_rank0.eq.root_rank0 ) then
    call SYSTEM_CLOCK(itm1,ntmr,maxtm)
    call CPU_TIME(etmres1)
    itm2 = itm1-itm0  
    if ( itm2 .lt. 0 ) then 
      itm2 = itm2 + maxtm
    endif
    write(STDOUT,'(A,X,E12.5,X,E12.5)')"Reading Greens Funtion....  / Time: ",etmres1-etmres0,dble(itm2)/dble(maxtm)
  endif
  call readGFn
  if ( my_rank0.eq.root_rank0 ) then
    call SYSTEM_CLOCK(itm1,ntmr,maxtm)
    call CPU_TIME(etmres1)
    itm2 = itm1-itm0  
    if ( itm2 .lt. 0 ) then 
      itm2 = itm2 + maxtm
    endif
    write(STDOUT,'(A,X,E12.5,X,E12.5)')"Reading Greens Funtion... Done! / Time: ",etmres1-etmres0,dble(itm2)/dble(maxtm)
  endif

!! Set Frictional and Other Parameters
  call setFricParam
  if ( my_rank0.eq.root_rank0 ) then
    write(STDOUT,*)"Friction parameters are set!"
  endif

!! Initial Conditions
  call setInitCond(nc,t,dt_next,y,dispnc)

  ibinout = 0
  iascout = 0

  vnc(:) = v_asta*exp(y(:,1))  ! Unit of vnc in the mainloop is [m/s]
  vnc_pre(:) = vnc(:)

  if ( IsnapSW .eq. 0 ) then
    TmNextWrite = DoutIntvYr
  else
    TmNextWrite = DoutIntvYr*(int(t/DoutIntvYr) + 1)
  endif

  if ( my_rank0.eq.root_rank0 ) then
    call SYSTEM_CLOCK(itm1,ntmr,maxtm)
    call CPU_TIME(etmres1)
    itm2 = itm1-itm0  
    if ( itm2 .lt. 0 ) then 
      itm2 = itm2 + maxtm
    endif
    write(STDOUT,'(A,X,E12.5,X,E12.5)')"Start Main Loop. / Time: ",etmres1-etmres0,dble(itm2)/dble(maxtm)
  endif

!! Main loop
  do istep = 1, MaxStep
    ibinout = ibinout + 1
    iascout = iascout + 1

    call DERIVSm(nc,NPrmDim,t,y,dydx)

    dt_try = dt_next
    yscal(:,:) = abs(y(:,:)) + abs(dt_try * dydx(:,:)) + TINY

    t0 = 0.d0
    call RKQSm(nc,NPrmDim,t0,y,dydx,yscal,ytmp1,ytmp2,ytmp3,dt_try,dt_done,dt_next,DERIVSm)
    t = t + dt_done * dc0/(Ysec * v_asta)

    yc(:,:) = y(:,:)
    vnc(1:nc) = v_asta*exp(y(1:nc,1))
    dispnc(1:nc) = dispnc(1:nc) + 0.5d0*(vnc(1:nc)+vnc_pre(1:nc))*dt_done*dc0/v_asta

    vnc_pre = vnc

!! Output Section
    if ( ibinout.ge.Nbinout .or. t.ge.TmNextWrite ) then
      TmNextWrite = DoutIntvYr*(int(t/DoutIntvYr) + 1)

      if (my_rank0.eq.root_rank0)then
        call CPU_TIME(etmres1)
        call SYSTEM_CLOCK(itm1,ntmr,maxtm)
        itm2 = itm1-itm0  
        if ( itm2 .lt. 0 ) then 
          itm2 = itm2 + maxtm
        endif
        write(6,'(I,X,E21.14,X,E12.5,X,E12.5)')istep,t,etmres1-etmres0,dble(itm2)/dble(ntmr)
      endif

      if ( ibinout.ge.Nbinout ) then
        ibinout = 0
      endif

      yout_nc(1,1:nc) = y(1:nc,1)
      yout_nc(2,1:nc) = y(1:nc,2)
      yout_nc(3,1:nc) = dispnc(1:nc)
      !yout_nc(4,1:nc)    = -as(1:nc)*log(vcut1/vnc(1:nc)+1.d0) +  &
      !                   bs(1:nc)*log(yc(1:nc,2)/g1(1:nc)+1.d0) - tauref(1:nc)

      call MPI_Gatherv(yout_nc,ncbinout,MPI_REAL8,yout_all,lstnc_binout,lstdispl0_binout,MPI_REAL8,root_rank0,MPI_COMM_WORLD,ierr)
      if ( my_rank0.eq.root_rank0 )  then
        write(FPBO) t
        write(FPBO) yout_all

        open(FPSSO,file="ForRestart.snap",access="stream",form="unformatted")
        write(FPSSO)t,dt_next
        write(FPSSO)Nall,NbinoutPrmDim,Nversion
        write(FPSSO)yout_all
        close(FPSSO)

        if ( iascout.ge.Nascout ) then
          ! ASCII output
          iascout = 0
          do i=1,Nall
            write(FPAO,'(E12.5,X,I,5(X,E12.5))')t,i,v_asta*exp(yout_all(1,i)),yout_all(2,i),xp(1,i),xp(3,i),yout_all(3,i)
          end do
        endif

      endif
    endif
  end do
  if ( my_rank0.eq.root_rank0 )  then
    close(FPBO)
    close(FPAO)
  endif
!! Finalize MPI
  call MPI_FINALIZE(ierr)
!!
  stop
end program sse3d_mpi

subroutine calcRankIndex
! Prepare index of ranks
  use mpi
  use common_param
  use common_var
  implicit none
  integer                :: i,j,i1,i2

  allocate(lstnc(n_rank),lstdispl(n_rank),lstdispl2(n_rank),lstdispl0(n_rank))
  allocate(lstdispl_binout(n_rank),lstdispl0_binout(n_rank),lstnc_binout(n_rank))

  j = Nall - int(Nall/n_rank)*n_rank
  if ( j.gt.0 ) then
    lstnc(1:(n_rank-j))        = int(Nall/n_rank)
    lstnc((n_rank-j+1):n_rank) = int(Nall/n_rank) + 1
  else
    lstnc(1:n_rank) = int(Nall/n_rank)
  endif
  lstnc_binout = lstnc*NbinoutPrmDim

  i1 = 0
  i2 = 0
  do i=1,n_rank
    i1 = i2 + 1
    i2 = i2 + lstnc(i)
    lstdispl(i) = i1
    lstdispl_binout(i) = (i1-1)*NbinoutPrmDim + 1
    if ( i2.le.Nall ) then
      lstrnk(i1:i2) = i
      lstdispl2(i) = i2
    else
      write(STDERR,*)"sse3Dmpi: MPI_RANK index checker, i2 exceed Nall!"
      lstrnk(i1:n_rank) = i
      lstdispl2(i) = Nall
    endif
  end do
  lstrnk0(1:Nall) = lstrnk(1:Nall) - 1
  lstdispl0(1:n_rank) = lstdispl(1:n_rank) - 1
  lstdispl0_binout(1:n_rank) = lstdispl_binout(1:n_rank) - 1

  if ( my_rank0.eq.root_rank0) then
    open(FPRNK,file="test_ranklist.txt")
    do i=1,n_rank
      write(FPRNK,'(7(I,X))')i,lstnc(i),lstdispl(i),lstdispl2(i),lstdispl0(i),lstdispl_binout(i),lstdispl0_binout(i)
    end do
    do i=1,Nall
      write(FPRNK,'(3(I,X))')i,lstrnk(i),lstrnk0(i)
    end do

    close(FPRNK)
  endif

  return
end subroutine calcRankIndex

!======================================================================
subroutine commVarAlloc
!! Allocate common variables 
!!       Exception: lstnc, lstdispl, lstdispl2 is allocated in subroutine "calcRankIndex"
  use mpi
  use common_param
  use common_var
  implicit none

  allocate(gfn(Nall,nc))       ! Green's function for normalized calculation
  allocate(AK2(nc,NPrmDim),AK3(nc,NPrmDim),AK4(nc,NPrmDim),AK5(nc,NPrmDim),AK6(nc,NPrmDim))
  allocate(xp_c(3,nc))
  allocate(vnc(nc),vncr(nc),vnc_pre(nc),dispnc(nc),yout_nc(NbinoutPrmDim,nc),vnpl_c(nc))
  allocate(as(nc),bs(nc),tauref(nc),g1(nc),g3(nc))
  allocate(vnr(Nall))
  allocate(yderivs(nc),yrkck(nc,NPrmDim))
  return
end subroutine commVarAlloc

!======================================================================
subroutine readGFn
  use common_param
  use common_var
  use mpi
  implicit none

  integer, parameter :: FP = 51
  integer :: i,j,ilocal,n_elem,n_vertex,ierr,istatus(MPI_STATUS_SIZE)
  real(8), allocatable :: gf(:)

  if (my_rank0.eq.root_rank0) then
    open(FP,file="green_bin.dat",access= "stream",form="unformatted")
    read(FP) n_elem, n_vertex
    if ( n_elem.ne.Nall ) then
      write(STDERR,*)"Warning number of subfaults is different. Nelement/Nall",n_elem,Nall
      write(STDERR,*)"These numbers should be the same. Check the program!!"
    endif
    read(FP) xp  ! read (x,y,z)
  endif 
  xp(3,:) = -xp(3,:)
  call MPI_Bcast(xp,3*Nall,MPI_REAL8,root_rank0,MPI_COMM_WORLD,ierr)

  do i=1,nc
    xp_c(1:3,i) = xp(1:3,i+lstdispl0(my_rank))
  end do

  allocate(gf(Nall))
  do i=1,Nall
    if (my_rank0.eq.root_rank0) then
      read(FP)gf
    else
      gf = 0.d0
    endif

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

!! gf is part of the elastostatic kernel normalized by rigidity. [m^-1]
!! gfn(i,j) : stress change at j-th element caused by the slip on i-th element.
    if ( lstrnk0(i).eq.root_rank0 ) then
      if ( my_rank0.eq.root_rank0 ) then
        ilocal = i - lstdispl0(my_rank)
        if ( ilocal.ge.1.and.ilocal.le.nc ) then
          gfn(1:Nall,ilocal) = gf(1:Nall)
        else
          write(STDERR,*)"Error! ilocal is out of array size.(ilocal, nc)",ilocal,nc
        endif
      endif
    else
      if ( my_rank0.eq.root_rank0) then
        call MPI_Send(gf,Nall,MPI_REAL8,lstrnk0(i),0,MPI_COMM_WORLD,ierr) 
      else if ( my_rank0.eq.lstrnk0(i) ) then
        call MPI_Recv(gf,Nall,MPI_REAL8,root_rank0,0,MPI_COMM_WORLD,istatus,ierr)
        ilocal = i - lstdispl0(my_rank)
        if ( ilocal.ge.1.and.ilocal.le.nc ) then
          gfn(1:Nall,ilocal) = gf(1:Nall)
        else
          write(STDERR,*)"Error! ilocal is out of array size.(ilocal, nc)",ilocal,nc
        endif
      endif
    endif
  end do
  close(FP)
  deallocate(gf)

!! Unit of gfn is changed from [m^-1] to [MPa].
  gfn = gfn*dc0*rigid

  return 
end subroutine readGFn

!======================================================================
subroutine setFricParam
  use common_param
  use common_var
  use fric_param
  use mpi
  implicit none
  integer, parameter :: Ntremormax = 100000000
  integer :: i,j,ierr
  real(8) :: zab(n_ab),zdc(n_dc),zse(n_seff),zvc(n_vc),theta
  real(8) :: arand(Nall)
  real(8), allocatable :: a_c(:),b_c(:),dc_c(:),se_c(:),theta_c(:),buf(:),ckbuf(:,:),xt(:),yt(:),zt(:),ntr_c(:)
  integer(4) :: ntr(Nall)
  real(8) :: trmdist,xtmp,ytmp,fac

  allocate(a_c(nc),b_c(nc),dc_c(nc),theta_c(nc),se_c(nc),xt(nc),yt(nc),zt(nc),ntr_c(nc))

  ntr(:) = 0
  open(FPTRM,file="tremor.rot.pos")
  do i=1,Ntremormax
    read(FPTRM,*,end=11)xtmp,ytmp
    do j=1,Nall
      trmdist = sqrt((xtmp - xp(1,j))**2 + (ytmp-xp(2,j))**2)
      if ( trmdist .le. ThTrmdist ) then
        ntr(j) = ntr(j) + 1
      endif
    end do
  end do
11 close(FPTRM)
  call MPI_Bcast(ntr,Nall,MPI_INTEGER4,root_rank0,MPI_COMM_WORLD,ierr)
 
  do i=1,nc
    ntr_c(i) = ntr(i+lstdispl0(my_rank))
  end do

  ! Calculation of random values
  call arand1(Nall,arand)

  xt(1:nc)  = xp_c(1,1:nc)
  yt(1:nc)  = xp_c(2,1:nc)
  zt(1:nc) = xp_c(3,1:nc)

  do i=1,nc

    if ( xt(i).le.xvpl0 ) then
      vnpl_c(i) = vpl0/v_asta
    else if ( xt(i).gt.xvpl0.and.xt(i).le.xvpl1 ) then
      vnpl_c(i) = ((vpl1-vpl0)*(xt(i)-xvpl0)/(xvpl1-xvpl0) + vpl0)/v_asta
    else
      vnpl_c(i) = vpl1/v_asta
    endif

    zvc(:) = zvc_rs(:)

    !! Set depth parameter for L-SSE region
    if      ( xt(i).gt.xsse1.and.xt(i).le.xsse2 ) then
      fac = (xt(i)-xsse1)/(xsse2-xsse1)
    else if ( xt(i).gt.xsse2.and.xt(i).le.xsse3 ) then
      fac = 1.d0
    else if ( xt(i).gt.xsse3.and.xt(i).le.xsse4 ) then
      fac = (xsse4-xt(i))/(xsse4-xsse3)
    else
      fac = 0.d0
    endif  
    zab(:) = (zab_sse(:)  -zab_rs(:))*fac + zab_rs(:)
    zdc(:) = (zdc_sse(:)  -zdc_rs(:))*fac + zdc_rs(:)
    zse(:) = (zseff_sse(:)-zseff(:) )*fac + zseff(:)

    !! Modification for boundary (edge) region
    if ( xt(i).lt.xbnd2 .or. xt(i).gt.xbnd3 ) then
      if ( xt(i).gt.xbnd1.and.xt(i).lt.xbnd2 ) then
        fac = (xbnd2-xt(i))/(xbnd2-xbnd1)
      else if ( xt(i).gt.xbnd3.and.xt(i).lt.xbnd4 ) then
        fac = (xt(i)-xbnd3)/(xbnd4-xbnd3)
      else
        fac = 1.d0
      endif
      zdc(:) = (zdc_bnd(:)  -zdc_rs(:))*fac + zdc_rs(:)
      zse(:) = (zseff_bnd(:)-zseff(:) )*fac + zseff(:)
    endif

!! Set a and b 
    if ( zt(i).lt.zab(1) ) then
      a_c(i) = a_rs(1) 
      b_c(i) = b_rs(1) 
    else if (zt(i).ge.zab(n_ab)) then
      a_c(i) = a_rs(n_ab) 
      b_c(i) = b_rs(n_ab) 
    else
      do j=1,n_ab-1
        if ( zt(i).ge.zab(j).and.zt(i).lt.zab(j+1)) then
          fac = (zt(i)-zab(j))/(zab(j+1)-zab(j))
          a_c(i) = (a_rs(j+1) -a_rs(j))*fac + a_rs(j)
          b_c(i) = (b_rs(j+1) -b_rs(j))*fac + b_rs(j)
          cycle   
        endif
      end do 
    endif

    !! Modification for tremor region (S-SSE region)
    if (ntr_c(i).ge.NtrThresh .and. zt(i).ge.zab_trm(1) .and. zt(i).le.zab_trm(n_ab_trm) ) then
      do j=1,n_ab_trm -1
        if ( zt(i).ge.zab_trm(j).and.zt(i).lt.zab_trm(j+1)) then
          fac = (zt(i)-zab_trm(j))/(zab_trm(j+1)-zab_trm(j))
          a_c(i) = (a_trm(j+1) -a_trm(j))*fac + a_trm(j)
          b_c(i) = (b_trm(j+1) -b_trm(j))*fac + b_trm(j)
          cycle  
        endif
      end do
    endif

    !! Modification for boundary (edge) region
    if ( zt(i).le.zbbnd ) then
      if ( xt(i).lt.xbnd1.or.xt(i).ge.xbnd4 ) then
        fac = 1.d0
      else if (xt(i).ge.xbnd1.and.xt(i).lt.xbnd2) then
        fac = (xbnd2-xt(i))/(xbnd2-xbnd1)
      else if (xt(i).ge.xbnd3.and.xt(i).lt.xbnd4) then
        fac = (xt(i)-xbnd3)/(xbnd4-xbnd3)
      else
        fac = 0.d0 
      endif
      b_c(i) = (bbnd - b_c(i))*fac + b_c(i)
    endif

!! Fluctuation of b
    b_c(i) = b_c(i)*(1.d0 + randfac*arand(lstdispl0(my_rank)+i))

!! Set dc
    if ( zt(i).lt.zdc(1) ) then
      dc_c(i) = dc_rs(1)
    else if (zt(i).ge.zdc(n_dc)) then
      dc_c(i) = dc_rs(n_dc)
    else
      do j=1,n_dc-1
        if ( zt(i).ge.zdc(j).and.zt(i).lt.zdc(j+1)) then
          dc_c(i) = (dc_rs(j+1)-dc_rs(j))*(zt(i)-zdc(j))/(zdc(j+1)-zdc(j)) + dc_rs(j)
          cycle
        endif
      end do
    endif

!! Set seff 
    if ( zt(i).lt.zse(1) ) then
      se_c(i) = seff(1)
    else if (zt(i).ge.zse(n_seff)) then
      se_c(i) = seff(n_seff)
    else
      do j=1,n_seff-1
        if ( zt(i).ge.zse(j).and.zt(i).lt.zse(j+1)) then
          se_c(i) = (seff(j+1) - seff(j))*(zt(i)-zse(j))/(zse(j+1)-zse(j)) + seff(j)
          cycle
        endif
      end do
    endif
    
!! Set cutoff velocity
   if ( zt(i).lt.zvc(1) ) then
      theta = vc_rs(1)
    else if (zt(i).ge.zvc(n_vc)) then
      theta = vc_rs(n_vc)
    else
      do j=1,n_vc-1
        if ( zt(i).ge.zvc(j).and.zt(i).lt.zvc(j+1)) then
          theta = (vc_rs(j+1)-vc_rs(j))*(zt(i)-zvc(j))/(zvc(j+1)-zvc(j)) + vc_rs(j)
          cycle
        endif
      end do
    endif
    theta_c(i) = 10**(theta)
 
  end do

!! Set variables for DERIVS
  as(1:nc) = a_c(1:nc)*se_c(1:nc) 
  bs(1:nc) = b_c(1:nc)*se_c(1:nc) 
  g1(1:nc) = dc_c(1:nc)/theta_c(1:nc)
  g3(1:nc) = dc0/dc_c(1:nc)

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

!! Check output of frictional parameters

  if ( my_rank0.eq.root_rank0 ) then
    allocate( buf(Nall), ckbuf(Nfricck,Nall) )
    write(STDOUT,*)"Start to gather frictinal parameters."
    do i=1,Nall
      ckbuf(1:3,i) = xp(1:3,i)
    end do
  endif

  call MPI_Gatherv(a_c,nc,MPI_REAL8,buf,lstnc,lstdispl0,MPI_REAL8,root_rank0,MPI_COMM_WORLD,ierr) ! Gather a
  if ( my_rank0.eq.root_rank0 ) then
    ckbuf(4,1:Nall) = buf(1:Nall)
  endif

  call MPI_Gatherv(b_c,nc,MPI_REAL8,buf,lstnc,lstdispl0,MPI_REAL8,root_rank0,MPI_COMM_WORLD,ierr) ! Gather b
  if ( my_rank0.eq.root_rank0 ) then
    ckbuf(5,1:Nall) = buf(1:Nall)
  endif

  call MPI_Gatherv(se_c,nc,MPI_REAL8,buf,lstnc,lstdispl0,MPI_REAL8,root_rank0,MPI_COMM_WORLD,ierr) ! Gather s_eff
  if ( my_rank0.eq.root_rank0 ) then
    ckbuf(6,1:Nall) = buf(1:Nall)
  endif

  call MPI_Gatherv(dc_c,nc,MPI_REAL8,buf,lstnc,lstdispl0,MPI_REAL8,root_rank0,MPI_COMM_WORLD,ierr) ! Gather dc
  if ( my_rank0.eq.root_rank0 ) then
    ckbuf(7,1:Nall) = buf(1:Nall)
  endif

  call MPI_Gatherv(theta_c,nc,MPI_REAL8,buf,lstnc,lstdispl0,MPI_REAL8,root_rank0,MPI_COMM_WORLD,ierr) ! Gather theta_c
  if ( my_rank0.eq.root_rank0 ) then
    ckbuf(8,1:Nall) = buf(1:Nall)
  endif

  call MPI_Gatherv(g1,nc,MPI_REAL8,buf,lstnc,lstdispl0,MPI_REAL8,root_rank0,MPI_COMM_WORLD,ierr) ! Gather g1
  if ( my_rank0.eq.root_rank0 ) then
    ckbuf(9,1:Nall) = buf(1:Nall)
  endif

  call MPI_Gatherv(g3,nc,MPI_REAL8,buf,lstnc,lstdispl0,MPI_REAL8,root_rank0,MPI_COMM_WORLD,ierr) ! Gather g3
  if ( my_rank0.eq.root_rank0 ) then
    ckbuf(10,1:Nall) = buf(1:Nall)
  endif

  call MPI_Gatherv(as,nc,MPI_REAL8,buf,lstnc,lstdispl0,MPI_REAL8,root_rank0,MPI_COMM_WORLD,ierr) ! Gather as
  if ( my_rank0.eq.root_rank0 ) then
    ckbuf(11,1:Nall) = buf(1:Nall)
  endif

  call MPI_Gatherv(bs,nc,MPI_REAL8,buf,lstnc,lstdispl0,MPI_REAL8,root_rank0,MPI_COMM_WORLD,ierr) ! Gather bs
  if ( my_rank0.eq.root_rank0 ) then
    ckbuf(12,1:Nall) = buf(1:Nall)
  endif

  call MPI_Gatherv(vnpl_c,nc,MPI_REAL8,buf,lstnc,lstdispl0,MPI_REAL8,root_rank0,MPI_COMM_WORLD,ierr) ! Gather vnpl_c
  if ( my_rank0.eq.root_rank0 ) then
    ckbuf(13,1:Nall) = buf(1:Nall)
  endif

  tauref(:) = -as(:)*log(vcut1/(v_asta*vnpl_c(:))+1.d0) + bs(:)*log(dc0/(v_asta*vnpl_c(:)*g1(:)*g3(:))+1.d0)
  call MPI_Gatherv(tauref,nc,MPI_REAL8,buf,lstnc,lstdispl0,MPI_REAL8,root_rank0,MPI_COMM_WORLD,ierr) ! Gather tauref
  if ( my_rank0.eq.root_rank0 ) then
    ckbuf(14,1:Nall) = buf(1:Nall)
  endif
  
  deallocate(a_c,b_c,dc_c,theta_c)

  if ( my_rank0.eq.root_rank0 ) then
    open(FPCKO,file="fric_ck.dat")
    do i=1,Nall
      write(FPCKO,'(I,14(E16.8,X))')i,(ckbuf(j,i),j=1,Nfricck)
    end do

    write(FPBO)ckbuf
    deallocate(buf,ckbuf)
  endif
  return
end subroutine setFricParam

!======================================================================
subroutine setInitCond(n,t,dt_next,y,disp)
  use common_param
  use common_var
  use mpi
  implicit none
  integer,intent(in)    :: n
  real(8),intent(inout) :: t,dt_next, y(n,NprmDim),disp(n)
  integer, parameter :: FP =  53
  integer :: Nallcomp,Nprmdimcomp,NVercomp
  real(8) :: t1,t2
  real(8), allocatable :: yinit(:,:)
  integer :: i,ierr

  if ( IsnapSW.eq.0) then
     t = 0.d0
     dt_next = dt_init
     y(1:n,1) = log( v_init/v_asta )
     y(1:n,2) = 0.05d0*dc0/(g3(1:n)*v_init)
     dispnc(1:n) = 0.d0
  else
    allocate(yinit(NbinoutPrmDim,Nall))
    if ( my_rank0.eq.root_rank0 ) then
      open(FP,file="ForRestartIni.snap",access="stream",form="unformatted")
      read(FP) t1,t2
      read(FP) Nallcomp,Nprmdimcomp,NVercomp
      read(FP) yinit
      close(FP)
      t = t1
      dt_next = t2
      if ( Nallcomp .ne. Nall ) then
        write(STDERR,*)"Nall is inconsistent (Restartfile,Program)",Nallcomp,Nall
        stop
      endif
      if ( Nprmdimcomp .ne. NbinoutPrmDim ) then
        write(STDERR,*)"NbinoutPrmDim is inconsistent(Restartfile,Program)",Nprmdimcomp,NbinoutPrmDim
        stop
      endif
      if ( Nvercomp .ne. Nversion ) then
        write(STDERR,*)"Version is inconsistent(Restartfile,Program)",Nvercomp,Nversion
        stop
      endif
    endif

    call MPI_Bcast(t,1,MPI_REAL8,root_rank0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(dt_next,1,MPI_REAL8,root_rank0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(yinit,Nall*NbinoutPrmDim,MPI_REAL8,root_rank0,MPI_COMM_WORLD,ierr)

    do i=1,n
      y(i,1)  = yinit(1,i+lstdispl0(my_rank))
      y(i,2)  = yinit(2,i+lstdispl0(my_rank))
      disp(i) = yinit(3,i+lstdispl0(my_rank))
    end do
    
    deallocate(yinit)
  endif

  return
end subroutine setInitCond

!======================================================================
subroutine DERIVSm(n,m,x,y,dydx)
! Evaluate derivatitves
!   y(:,1) = log(velocity)
!   y(:,2) = theta (i.e., state variable)
!   y(:,3) = displacement (evaluated in the main loop)
  use mpi
  use common_param
  use common_var
  implicit none
  integer, intent(in)    :: n,m
  real(8), intent(in)    :: x,y(n,m)
  real(8), intent(inout) :: dydx(n,m)
  integer :: i,ierr

  vnc(:)  = exp(y(:,1))
  vncr(:) = vnpl_c(:) - vnc(:)
  call MPI_Allgatherv(vncr,n,MPI_REAL8,vnr,lstnc,lstdispl0,MPI_REAL8, MPI_COMM_WORLD,ierr) ! Gather Velocity

  dydx(:,2) = fac2 - y(:,2)*vnc(:)*g3(:)

  !!!!!! If BLAS is used comment out following lines and return to comment the below block !!!!!!
  yderivs(:) = dydx(:,2)*bs(:)/(y(:,2)+g1(:))
  call DGEMV("T",Nall,n,1.d0,gfn,Nall,vnr,1,-1.d0,yderivs,1)
  dydx(:,1) = yderivs(:) / (as(:)/(1.d0+vnc(:)*theta2in)+fac1*vnc(:))
  !!!!!!

  !!!!!! If BLAS is not used, comment out following lines and return to comment the above block !!!!!!
  !do i=1,n
  !  dydx(i,1)     = dot_product(gfn(1:Nall,i),vnr(1:Nall))
  !end do
  !dydx(:,1) = dydx(:,1) - dydx(:,2)*bs(:)/(y(:,2)+g1(:))
  !dydx(:,1) = dydx(:,1) / (as(:)/(1.d0+vnc(:)*theta2in)+fac1*vnc(:))
  !!!!!!

  return
end subroutine DERIVSm

subroutine RKQSm(n,m,x,y,dydx,yscal,ycand,yerr,ytmp,htry,hdid,hnext,DERIVSm)
!------------------------------------------------------------------------------
!       FIFTH-ORDER RUNGE-KUTTA  : STEEPER ROUTINE
!       SEE NUMERICAL RECIPES 2ND. ED. P.712
!------------------------------------------------------------------------------
  use mpi
  use common_param
  use common_var
  implicit none

  integer,intent(in)    :: n,m
  real(8),intent(in)    :: htry,yscal(n,m),ytmp(n,m)
  real(8),intent(inout) :: x,dydx(n,m),y(n,m),ycand(n,m),yerr(n,m)
  real(8),intent(inout) :: hdid,hnext
  external              :: DERIVSm
  real(8), parameter    :: SAFETY=0.7d0,PGROW=-0.2d0,PSHRNK=-0.15d0,ERRCON=1.89d-4
  integer               :: i,ierr
  real(8)               :: h,errmax,errmax_thrd,xnew

  h=htry

  do while (.true.)
    call RKCKm(n,m,x,y,dydx,h,ycand,yerr,ytmp,DERIVSm)
    errmax_thrd=maxval(abs(yerr(1:n,1:m)/yscal(1:n,1:m)))

    errmax = errmax_thrd

    call MPI_Allreduce(errmax_thrd,errmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)

    errmax = errmax/eps
    if (errmax .gt. 1.d0) then
      h=h*SAFETY*(errmax**PSHRNK)
      if (h .lt. 0.5d0*h) then
        h=0.5d0*h
      endif
      xnew=x+h
      if (xnew .eq. x) pause 'Stepsize underflow in RKQS'
    else
      if (errmax.gt.ERRCON) then
        hnext=SAFETY*h*(errmax**PGROW)
      else
        hnext=1.d1*h
      endif
      hdid=h
      x=x+h
      y(1:n,1:m)=ycand(1:n,1:m)
      exit 
    endif
  end do
  return
end subroutine RKQSm

subroutine RKCKm(n,m,x,yin,dydx,h,yout,yerr,ytmp,DERIVSm)
!------------------------------------------------------------------------------
!       FIFTH-ORDER RUNGE-KUTTA  : ALGORITHM ROUTINE
!       SEE NUMERICAL RECIPES 2ND. ED. P.713
!------------------------------------------------------------------------------
  use mpi
  use common_param
  use common_var
  implicit none

  integer, intent(in)    :: n,m
  real(8), intent(in)    :: x,h
  real(8), intent(inout) :: yin(n,m),dydx(n,m),yout(n,m),yerr(n,m),ytmp(n,m)
  external :: DERIVSm

  real(8), parameter :: A2=0.2d0, A3=0.3d0, A4=0.6d0, A5=1.d0, A6=0.875d0
  real(8), parameter :: B21=0.2d0, B31=3.d0/40.d0, B32=9.d0/40.d0
  real(8), parameter :: B41=0.3d0, B42=-0.9d0, B43=1.2d0
  real(8), parameter :: B51=-11.d0/54.d0, B52=2.5d0, B53=-70.d0/27.d0, B54=35.d0/27.d0  
  real(8), parameter :: B61=1631.d0/55296.d0,   B62=175.d0/512.d0, B63=575.d0/13824.d0, &
                        B64=44275.d0/110592.d0, B65=253.d0/4096.d0
  real(8), parameter :: C1=37.d0/378.d0, C3=250.d0/621.d0, C4=125.d0/594.d0, C6=512.d0/1771.d0
  real(8), parameter :: DC1=C1-2825.d0/27648.d0,  DC3=C3-18575.d0/48384.d0,        &
                        DC4=C4-13525.d0/55296.d0, DC5=-277.d0/14336.d0, DC6=C6-0.25d0
  integer :: nm 

  nm = n*m

  call dcopy(nm,yin,1,ytmp,1)
  call daxpy(nm,h*B21,dydx,1,ytmp,1)

  !ytmp = yin + h*B21*dydx
  call DERIVSm(n,m,x+A2*h,ytmp,AK2)

  call dcopy(nm,dydx,1,    yrkck,1)
  call dscal(nm,B31,       yrkck,1)
  call daxpy(nm,B32,AK2 ,1,yrkck,1)
  call dcopy(nm,yin,1,ytmp,1)
  call daxpy(nm,h,yrkck,1,ytmp,1)

  !ytmp = yin + h*(B31*dydx + B32*AK2)
  call DERIVSm(n,m,x+A3*h,ytmp,AK3)

  call dcopy(nm,dydx,1,    yrkck,1)
  call dscal(nm,B41,       yrkck,1)
  call daxpy(nm,B42,AK2 ,1,yrkck,1)
  call daxpy(nm,B43,AK3 ,1,yrkck,1)
  call dcopy(nm,yin,1,ytmp,1)
  call daxpy(nm,h,yrkck,1,ytmp,1)

  !ytmp = yin + h*(B41*dydx + B42*AK2 + B43*AK3)
  call DERIVSm(n,m,x+A4*h,ytmp,AK4)

  call dcopy(nm,dydx,1,    yrkck,1)
  call dscal(nm,B51,       yrkck,1)
  call daxpy(nm,B52,AK2 ,1,yrkck,1)
  call daxpy(nm,B53,AK3 ,1,yrkck,1)
  call daxpy(nm,B54,AK4 ,1,yrkck,1)
  call dcopy(nm,yin,1,ytmp,1)
  call daxpy(nm,h,yrkck,1,ytmp,1)

  !ytmp = yin + h*(B51*dydx + B52*AK2 + B53*AK3 + B54*AK4)
  call DERIVSm(n,m,x+A5*h,ytmp,AK5)

  call dcopy(nm,dydx,1,    yrkck,1)
  call dscal(nm,B61,       yrkck,1)
  call daxpy(nm,B62,AK2 ,1,yrkck,1)
  call daxpy(nm,B63,AK3 ,1,yrkck,1)
  call daxpy(nm,B64,AK4 ,1,yrkck,1)
  call daxpy(nm,B65,AK5 ,1,yrkck,1)
  call dcopy(nm,yin,1,ytmp,1)
  call daxpy(nm,h,yrkck,1,ytmp,1)

  !ytmp = yin + h*(B61*dydx + B62*AK2 + B63*AK3 + B64*AK4 + B65*AK5)
  call DERIVSm(n,m,x+A6*h,ytmp,AK6)

  call dcopy(nm,dydx,1,   yrkck,1)
  call dscal(nm,C1,       yrkck,1)
  call daxpy(nm,C3,AK3 ,1,yrkck,1)
  call daxpy(nm,C4,AK4 ,1,yrkck,1)
  call daxpy(nm,C6,AK6 ,1,yrkck,1)
  call dcopy(nm,yin,1,yout,1)
  call daxpy(nm,h,yrkck,1,yout,1)
  !yout = yin + h*(C1*dydx + C3*AK3 + C4*AK4 + C6*AK6)

  call dcopy(nm,dydx,1,    yerr,1)
  call dscal(nm,DC1,       yerr,1)
  call daxpy(nm,DC3,AK3 ,1,yerr,1)
  call daxpy(nm,DC4,AK4 ,1,yerr,1)
  call daxpy(nm,DC5,AK5 ,1,yerr,1)
  call daxpy(nm,DC6,AK6 ,1,yerr,1)
  call dscal(nm,h,yerr,1)
  !yerr = h*(DC1*dydx + DC3*AK3 + DC4*AK4 + DC5*AK5 + DC6*AK6)

  return
end subroutine RKCKm


!=============================================================
subroutine arand1(n,arand)
  use common_param
  implicit none
  integer, intent(in)    :: n
  real(8), intent(inout) :: arand(n)
  
  integer, parameter :: MBIG=1000000000,MSEED=161803398,MZ=0
  real(8), parameter :: fac=1.d0/dble(MBIG)
  integer            :: idum, iff, inext, inextp
  integer            :: i,k,ii,ik,mj,mk
  integer            :: ma(55)

  iff = 0
  do ik=1,n
    idum=ik
    if (idum.lt.0.or.iff.eq.0) then
      iff=1
      mj=MSEED-abs(idum)
      mj=mod(mj,MBIG)
      ma(55)=mj
      mk=1
      do i=1,54
        ii=mod(21*i,55)
        ma(ii)=mk
        mk=mj-mk
        if (mk.lt.mz) mk=mk+MBIG
        mj=ma(ii)
      enddo
      do k=1,4
        do i=1,55
          ma(i)=ma(i)-ma(1+mod(i+30,55))
          if (ma(i).lt.mz) ma(i)=ma(i)+MBIG
        enddo
      enddo
      inext =0
      inextp=31
      idum  =1
    endif
    inext=inext+1
    if (inext.eq.56) inext=1
    inextp=inextp+1
    if (inextp.eq.56) inextp=1
    mj=ma(inext)-ma(inextp)
    if (mj.lt.mz) mj=mj+MBIG
    ma(inext)=mj
    arand(ik)=mj*fac
  enddo

  return
end subroutine arand1

