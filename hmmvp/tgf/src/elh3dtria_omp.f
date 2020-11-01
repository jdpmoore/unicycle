!/*********************************************************************
!*                                                                    *
!*  CMG Multipole Mediator for Earthquake Simulators                  *
!*                                                                    *
!*  Developed by a collaboration between Brown University, the        *
!*  University of California at Riverside, New York University, and   *
!*  Invisible Software, Inc.  Funding is provided by the National     *
!*  Science Foundation, award DMS-0934711.                            *
!*                                                                    *
!*  PUBLIC DOMAIN -- The software is in the public domain.  You are   *
!*  free to use, copy, modify, and redistribute the software.         *
!*                                                                    *
!*  NO WARRANTY -- The software is provided "AS IS" and without       *
!*  warranty of any kind, express or implied, including warranties    *
!*  of merchantability and fitness for a particular purpose.          *
!*                                                                    *
!*  DISCLAIMER -- In no event will the NSF, Brown University, the     *
!*  University of California at Riverside, New York University,       *
!*  Invisible Software, or any other developer or distributor of the  *
!*  software be liable to you for any damages arising out of the use  *
!*  or inability to use the software (including but not limited to    *
!*  lost profits, lost savings, lost data, inaccurate results,        *
!*  repair costs, incidental or consequential damages, and losses     *
!*  sustained by third parties), even if advised of the possibility   *
!*  of such damages.                                                  *
!*                                                                    *
!*********************************************************************/


!//////////////////////////////////////////////////////////////////////
!//
!// elh3dtria.f: OpenMP subroutines for elh3dtria.f.
!//
!//////////////////////////////////////////////////////////////////////


!/*********************************************************************
!
! Author: Michael Barall, Invisible Software Inc.
!
! Contains code Copyright (C) 2009: Leslie Greengard and Zydrunas Gimbutas
!
!*********************************************************************/








        subroutine elh3dtriadirecttarg_omp_1(
     $     j,
     $     RLAM,RMU,TRIANGLE,TRINORM,NSOURCE,SOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ifptfrc,ptfrc,ifstrain,strain)

        use ifdim_def
        implicit none

        integer, intent(in) :: j
        ! Shared caller's arguments
        real(8), intent(in) :: rlam
        real(8), intent(in) :: rmu
        integer, intent(in) :: nsource
        real(8), intent(in) :: triangle(3,3,nsource)
        real(8), intent(in) :: trinorm(3,nsource)
        real(8), intent(in) :: source(3,nsource)
        integer, intent(in) :: ifsingle
        real(8), intent(in) :: sigma_sl(3,ifdim(nsource,ifsingle))
        integer, intent(in) :: ifdouble
        real(8), intent(in) :: sigma_dl(3,ifdim(nsource,ifdouble))
        integer, intent(in) :: ifptfrc
        real(8), intent(inout) :: ptfrc(3,ifdim(nsource,ifptfrc))
        integer, intent(in) :: ifstrain
        real(8), intent(inout) :: strain(3,3,ifdim(nsource,ifstrain))


        integer :: ione


C_PRIVATE(i,j,ptfrc0,strain0,numfunev)

        integer :: i
!!        integer :: j
        real(8) :: ptfrc0(3)
        real(8) :: strain0(3,3)
        integer :: numfunev


        ione=1


        ! Original OpenMP code

        do i=1,nsource
c
        if (ifsingle .eq. 1 ) then
        if( i .eq. j ) then
        call elust3triadirectself
     $     (rlam,rmu,ione,ione,triangle(1,1,i),
     $     sigma_sl(1,i),
     1     source(1,j),ptfrc0,ifstrain,strain0)
           ! alias OK, repeated ione args are intent(in)
        else
        call elust3triadirecttarg
     $     (rlam,rmu,ione,triangle(1,1,i),
     $     sigma_sl(1,i),
     1     source(1,j),ptfrc0,ifstrain,strain0)
        endif
        if (ifptfrc .eq. 1) then
        ptfrc(1,j)=ptfrc(1,j)+ptfrc0(1)
        ptfrc(2,j)=ptfrc(2,j)+ptfrc0(2)
        ptfrc(3,j)=ptfrc(3,j)+ptfrc0(3)
        endif
        if (ifstrain .eq. 1) then
        strain(1,1,j)=strain(1,1,j)+strain0(1,1)
        strain(2,1,j)=strain(2,1,j)+strain0(2,1)
        strain(3,1,j)=strain(3,1,j)+strain0(3,1)
        strain(1,2,j)=strain(1,2,j)+strain0(1,2)
        strain(2,2,j)=strain(2,2,j)+strain0(2,2)
        strain(3,2,j)=strain(3,2,j)+strain0(3,2)
        strain(1,3,j)=strain(1,3,j)+strain0(1,3)
        strain(2,3,j)=strain(2,3,j)+strain0(2,3)
        strain(3,3,j)=strain(3,3,j)+strain0(3,3)
        endif
        endif
        if (ifdouble .eq. 1) then
        if( i .eq. j ) then 
        call eltst3triadirectself
     $     (rlam,rmu,ione,ione,triangle(1,1,i),
     $     sigma_dl(1,i),trinorm(1,i),
     1     source(1,j),ptfrc0,ifstrain,strain0)
           ! alias OK, repeated ione args are intent(in)
        else
        call eltst3triadirecttarg
     $     (rlam,rmu,ione,triangle(1,1,i),
     $     sigma_dl(1,i),trinorm(1,i),
     1     source(1,j),ptfrc0,ifstrain,strain0)
        endif
        if (ifptfrc .eq. 1) then
        ptfrc(1,j)=ptfrc(1,j)+ptfrc0(1)
        ptfrc(2,j)=ptfrc(2,j)+ptfrc0(2)
        ptfrc(3,j)=ptfrc(3,j)+ptfrc0(3)
        endif
        if (ifstrain .eq. 1) then
        strain(1,1,j)=strain(1,1,j)+strain0(1,1)
        strain(2,1,j)=strain(2,1,j)+strain0(2,1)
        strain(3,1,j)=strain(3,1,j)+strain0(3,1)
        strain(1,2,j)=strain(1,2,j)+strain0(1,2)
        strain(2,2,j)=strain(2,2,j)+strain0(2,2)
        strain(3,2,j)=strain(3,2,j)+strain0(3,2)
        strain(1,3,j)=strain(1,3,j)+strain0(1,3)
        strain(2,3,j)=strain(2,3,j)+strain0(2,3)
        strain(3,3,j)=strain(3,3,j)+strain0(3,3)
        endif
        endif
        enddo
c
c       ... image contribution
c
        if (ifsingle .eq. 1 ) then
        call eluh3triaadap        
     $     (rlam,rmu,nsource,triangle,
     $     sigma_sl,trinorm,source(1,j),
     1     ifptfrc,ptfrc0,ifstrain,strain0,numfunev)
        if( ifptfrc .eq. 1 ) then
        ptfrc(1,j)=ptfrc(1,j)+ptfrc0(1)
        ptfrc(2,j)=ptfrc(2,j)+ptfrc0(2)
        ptfrc(3,j)=ptfrc(3,j)+ptfrc0(3)
        endif
        if( ifstrain .eq. 1 ) then
        strain(1,1,j)=strain(1,1,j)+strain0(1,1)
        strain(2,1,j)=strain(2,1,j)+strain0(2,1)
        strain(3,1,j)=strain(3,1,j)+strain0(3,1)
        strain(1,2,j)=strain(1,2,j)+strain0(1,2)
        strain(2,2,j)=strain(2,2,j)+strain0(2,2)
        strain(3,2,j)=strain(3,2,j)+strain0(3,2)
        strain(1,3,j)=strain(1,3,j)+strain0(1,3)
        strain(2,3,j)=strain(2,3,j)+strain0(2,3)
        strain(3,3,j)=strain(3,3,j)+strain0(3,3)
        endif
        endif

        if (ifdouble .eq. 1 ) then
        call elth3triaadap        
     $     (rlam,rmu,nsource,triangle,
     $     sigma_dl,trinorm,source(1,j),
     1     ifptfrc,ptfrc0,ifstrain,strain0,numfunev)
        if( ifptfrc .eq. 1 ) then
        ptfrc(1,j)=ptfrc(1,j)+ptfrc0(1)
        ptfrc(2,j)=ptfrc(2,j)+ptfrc0(2)
        ptfrc(3,j)=ptfrc(3,j)+ptfrc0(3)
        endif
        if( ifstrain .eq. 1 ) then
        strain(1,1,j)=strain(1,1,j)+strain0(1,1)
        strain(2,1,j)=strain(2,1,j)+strain0(2,1)
        strain(3,1,j)=strain(3,1,j)+strain0(3,1)
        strain(1,2,j)=strain(1,2,j)+strain0(1,2)
        strain(2,2,j)=strain(2,2,j)+strain0(2,2)
        strain(3,2,j)=strain(3,2,j)+strain0(3,2)
        strain(1,3,j)=strain(1,3,j)+strain0(1,3)
        strain(2,3,j)=strain(2,3,j)+strain0(2,3)
        strain(3,3,j)=strain(3,3,j)+strain0(3,3)
        endif
        endif

        
        return
        end subroutine elh3dtriadirecttarg_omp_1








        subroutine elh3dtriadirecttarg_omp_2(
     $     j,
     $     RLAM,RMU,TRIANGLE,TRINORM,NSOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     NTARGET,target,
     $     ifptfrctarg,PTFRCtarg,ifstraintarg,STRAINtarg)

        use ifdim_def
        implicit none

        integer, intent(in) :: j
        ! Shared caller's arguments
        real(8), intent(in) :: rlam
        real(8), intent(in) :: rmu
        integer, intent(in) :: nsource
        real(8), intent(in) :: triangle(3,3,nsource)
        real(8), intent(in) :: trinorm(3,nsource)
        integer, intent(in) :: ifsingle
        real(8), intent(in) :: sigma_sl(3,ifdim(nsource,ifsingle))
        integer, intent(in) :: ifdouble
        real(8), intent(in) :: sigma_dl(3,ifdim(nsource,ifdouble))
        integer, intent(in) :: ntarget
        real(8) :: target(3,ntarget)
        integer, intent(in) :: ifptfrctarg
        real(8), intent(inout)
     &         :: ptfrctarg(3,ifdim(ntarget,ifptfrctarg))
        integer, intent(in) :: ifstraintarg
        real(8), intent(inout)
     &         :: straintarg(3,3,ifdim(ntarget,ifstraintarg))


C_PRIVATE(i,j,ptfrc0,strain0,numfunev)

!!        integer :: i
!!        integer :: j
        real(8) :: ptfrc0(3)
        real(8) :: strain0(3,3)
        integer :: numfunev


        ! Original OpenMP code

        if (ifsingle .eq. 1 ) then
        call elust3triadirecttarg
     $     (rlam,rmu,nsource,triangle,
     $     sigma_sl,
     1     target(1,j),ptfrc0,ifstraintarg,strain0)
        if (ifptfrctarg .eq. 1) then
        ptfrctarg(1,j)=ptfrctarg(1,j)+ptfrc0(1)
        ptfrctarg(2,j)=ptfrctarg(2,j)+ptfrc0(2)
        ptfrctarg(3,j)=ptfrctarg(3,j)+ptfrc0(3)
        endif
        if (ifstraintarg .eq. 1) then
        straintarg(1,1,j)=straintarg(1,1,j)+strain0(1,1)
        straintarg(2,1,j)=straintarg(2,1,j)+strain0(2,1)
        straintarg(3,1,j)=straintarg(3,1,j)+strain0(3,1)
        straintarg(1,2,j)=straintarg(1,2,j)+strain0(1,2)
        straintarg(2,2,j)=straintarg(2,2,j)+strain0(2,2)
        straintarg(3,2,j)=straintarg(3,2,j)+strain0(3,2)
        straintarg(1,3,j)=straintarg(1,3,j)+strain0(1,3)
        straintarg(2,3,j)=straintarg(2,3,j)+strain0(2,3)
        straintarg(3,3,j)=straintarg(3,3,j)+strain0(3,3)
        endif
        endif
        if (ifdouble .eq. 1) then
        call eltst3triadirecttarg
     $     (rlam,rmu,nsource,triangle,
     $     sigma_dl,trinorm,
     1     target(1,j),ptfrc0,ifstraintarg,strain0)
        if (ifptfrctarg .eq. 1) then
        ptfrctarg(1,j)=ptfrctarg(1,j)+ptfrc0(1)
        ptfrctarg(2,j)=ptfrctarg(2,j)+ptfrc0(2)
        ptfrctarg(3,j)=ptfrctarg(3,j)+ptfrc0(3)
        endif
        if (ifstraintarg .eq. 1) then
        straintarg(1,1,j)=straintarg(1,1,j)+strain0(1,1)
        straintarg(2,1,j)=straintarg(2,1,j)+strain0(2,1)
        straintarg(3,1,j)=straintarg(3,1,j)+strain0(3,1)
        straintarg(1,2,j)=straintarg(1,2,j)+strain0(1,2)
        straintarg(2,2,j)=straintarg(2,2,j)+strain0(2,2)
        straintarg(3,2,j)=straintarg(3,2,j)+strain0(3,2)
        straintarg(1,3,j)=straintarg(1,3,j)+strain0(1,3)
        straintarg(2,3,j)=straintarg(2,3,j)+strain0(2,3)
        straintarg(3,3,j)=straintarg(3,3,j)+strain0(3,3)
        endif
        endif
c
c       ... image contribution
c
        if (ifsingle .eq. 1 ) then
        call eluh3triaadap        
     $     (rlam,rmu,nsource,triangle,
     $     sigma_sl,trinorm,target(1,j),
     1     ifptfrctarg,ptfrc0,ifstraintarg,strain0,numfunev)
        if (ifptfrctarg .eq. 1) then
        ptfrctarg(1,j)=ptfrctarg(1,j)+ptfrc0(1)
        ptfrctarg(2,j)=ptfrctarg(2,j)+ptfrc0(2)
        ptfrctarg(3,j)=ptfrctarg(3,j)+ptfrc0(3)
        endif
        if (ifstraintarg .eq. 1) then
        straintarg(1,1,j)=straintarg(1,1,j)+strain0(1,1)
        straintarg(2,1,j)=straintarg(2,1,j)+strain0(2,1)
        straintarg(3,1,j)=straintarg(3,1,j)+strain0(3,1)
        straintarg(1,2,j)=straintarg(1,2,j)+strain0(1,2)
        straintarg(2,2,j)=straintarg(2,2,j)+strain0(2,2)
        straintarg(3,2,j)=straintarg(3,2,j)+strain0(3,2)
        straintarg(1,3,j)=straintarg(1,3,j)+strain0(1,3)
        straintarg(2,3,j)=straintarg(2,3,j)+strain0(2,3)
        straintarg(3,3,j)=straintarg(3,3,j)+strain0(3,3)
        endif
        endif

        if (ifdouble .eq. 1 ) then
        call elth3triaadap        
     $     (rlam,rmu,nsource,triangle,
     $     sigma_dl,trinorm,target(1,j),
     1     ifptfrctarg,ptfrc0,ifstraintarg,strain0,numfunev)
        if (ifptfrctarg .eq. 1) then
        ptfrctarg(1,j)=ptfrctarg(1,j)+ptfrc0(1)
        ptfrctarg(2,j)=ptfrctarg(2,j)+ptfrc0(2)
        ptfrctarg(3,j)=ptfrctarg(3,j)+ptfrc0(3)
        endif
        if (ifstraintarg .eq. 1) then
        straintarg(1,1,j)=straintarg(1,1,j)+strain0(1,1)
        straintarg(2,1,j)=straintarg(2,1,j)+strain0(2,1)
        straintarg(3,1,j)=straintarg(3,1,j)+strain0(3,1)
        straintarg(1,2,j)=straintarg(1,2,j)+strain0(1,2)
        straintarg(2,2,j)=straintarg(2,2,j)+strain0(2,2)
        straintarg(3,2,j)=straintarg(3,2,j)+strain0(3,2)
        straintarg(1,3,j)=straintarg(1,3,j)+strain0(1,3)
        straintarg(2,3,j)=straintarg(2,3,j)+strain0(2,3)
        straintarg(3,3,j)=straintarg(3,3,j)+strain0(3,3)
        endif
        endif

        
        return
        end subroutine elh3dtriadirecttarg_omp_2






