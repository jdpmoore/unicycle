cc Copyright (C) 2009-2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c       
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        This file contains the routines for elastostatic layer
c        potentials in free space in R^3. 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       User-callable routines are:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c      el3dtriadirecttarg - evaluates the elastostatic potential ON OR
c         OFF SURFACE due to a collection of flat triangles with
c         constant single and/or double layer densities using the direct
c         O(N^2) algorithm.
c
c
c      elust3triadirectself - evaluates the elastostatic potential at a
c         single surface point (triangle centroid) due to a collection
c         of flat triangles with piecewise constant single layer density
c         BY DIRECT CALCULATION.
c
c      elust3triadirecttarg - evaluates the elastostatic potential at an
c         (OFF SURFACE) target due to a collection of flat triangles
c         with piecewise constant single layer density BY DIRECT
c         CALCULATION.
c
c      eltst3triadirectself - evaluates the elastostatic potential at a
c         single surface point (triangle centroid) due to a collection
c         of flat triangles with piecewise constant double layer density
c         BY DIRECT CALCULATION.
c     
c      eltst3triadirecttarg - evaluates the elastostatic potential at an
c         (OFF SURFACE) target due to a collection of flat triangles
c         with piecewise constant double layer density BY DIRECT
c         CALCULATION.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine el3dtriadirecttarg(
     $     RLAM,RMU,TRIANGLE,TRINORM,NSOURCE,SOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ifptfrc,ptfrc,ifstrain,strain,NTARGET,
     $     target,ifptfrctarg,PTFRCtarg,
     $     ifstraintarg,STRAINtarg)
c
c
c       Elastostatic interactions in R^3: evaluate all pairwise triangle
c       interactions and interactions with targets using the direct
c       O(N^2) algorithm.
c
c       INPUT:
c
c       rlam,rmu - Lame parameters
c       triangle(3,3,nsource) - array of triangles in standard format
c       trinorm(3,nsource) - array of triangle normals
c       nsource - number of sources
c       source(3,nsource) - source locations
c       ifsingle - single layer computation flag  
c       sigma_sl(3,nsource) - vector strength of nth charge (single layer)
c       ifdouble - double layer computation flag  
c       sigma_dl(3,nsource) - vector strength of nth dipole (double layer)
c       ntarget - number of targets
c       target(3,ntarget) - evaluation target points
c       ifptfrc - displacement computation flag
c       ifstrain - strain computation flag
c       ifptfrctarg - target displacement computation flag
c       ifstraintarg - target strain computation flag
c
c       OUTPUT:
c
c       ptfrc(3,nsource) - displacement at source locations
c       strain(3,3,nsource) - strain at source locations
c       ptfrctarg(3,ntarget) - displacement at target locations
c       straintarg(3,3,ntarget) - strain at target locations
c
c
        use ifdim_def
        use falloc_thread
        implicit none

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
        integer, intent(in) :: ntarget
        real(8), intent(in) :: target(3,ntarget)
        integer, intent(in) :: ifptfrctarg
        real(8), intent(inout)
     &         :: ptfrctarg(3,ifdim(ntarget,ifptfrctarg))
        integer, intent(in) :: ifstraintarg
        real(8), intent(inout)
     &         :: straintarg(3,3,ifdim(ntarget,ifstraintarg))


        integer :: i
        integer :: ier
        integer :: j
c
c
c     NOTE: In the present version, dipole vectors SIGMA_DV must be SET EQUAL
c     to the triangle normal for elt3d direct routines
c
        do i=1,nsource
        if( ifptfrc .eq. 1) then
           ptfrc(1,i)=0
           ptfrc(2,i)=0
           ptfrc(3,i)=0
        endif
        if( ifstrain .eq. 1) then
           strain(1,1,i)=0
           strain(2,1,i)=0
           strain(3,1,i)=0
           strain(1,2,i)=0
           strain(2,2,i)=0
           strain(3,2,i)=0
           strain(1,3,i)=0
           strain(2,3,i)=0
           strain(3,3,i)=0
        endif
        enddo
c       
        do i=1,ntarget
        if( ifptfrctarg .eq. 1) then
           ptfrctarg(1,i)=0
           ptfrctarg(2,i)=0
           ptfrctarg(3,i)=0
        endif
        if( ifstraintarg .eq. 1) then
           straintarg(1,1,i)=0
           straintarg(2,1,i)=0
           straintarg(3,1,i)=0
           straintarg(1,2,i)=0
           straintarg(2,2,i)=0
           straintarg(3,2,i)=0
           straintarg(1,3,i)=0
           straintarg(2,3,i)=0
           straintarg(3,3,i)=0
        endif
        enddo
c
c       ... sources
c
c
        if( ifptfrc .eq. 1 .or. ifstrain .eq. 1 ) then
c
        call falloc_enter_parallel()
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(j)
        do j=1,nsource

        call el3dtriadirecttarg_omp_1(
     $     j,
     $     RLAM,RMU,TRIANGLE,TRINORM,NSOURCE,SOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ifptfrc,ptfrc,ifstrain,strain)

        enddo
C$OMP END PARALLEL DO
        call falloc_exit_parallel(ier)
c
        endif
c
c       ... targets
c
        if( ifptfrctarg .eq. 1 .or. ifstraintarg .eq. 1 ) then
c       
        call falloc_enter_parallel()
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(j)
        do j=1,ntarget

        call el3dtriadirecttarg_omp_2(
     $     j,
     $     RLAM,RMU,TRIANGLE,TRINORM,NSOURCE,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     NTARGET,
     $     target,ifptfrctarg,PTFRCtarg,
     $     ifstraintarg,STRAINtarg)

        enddo
C$OMP END PARALLEL DO
        call falloc_exit_parallel(ier)
c
        endif
c
        return
        end
c
c
c
c
c
c***********************************************************************
c
c       Quadrature routines for elastostatic single and double layers
c       Constant densities on flat triangles
c
c***********************************************************************
c
c
        subroutine elust3triadirecttarg_one
     $     (rlam,rmu,triangle,sigma_sl,ifself,target,
     $     ptfrc,ifstrain,strain)
c
c
c     Direct evaluation of displacement and strain due to constant
c     single layer elastostatic kernel on a flat triangle.
c
c     Double layer elastostatic kernel: constant-density on flat triangles
c
c     Computes displacement and strain at arbitrary point TARGET
c     due to piecewise-constant single layer density on a triangle.
c
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     rlam,rmu             Lame parameters
c     sigma_sl(3)          SLP strengths (constant)
c     triangle(3,3)        vertices of the triangle in standard format
c     target(3)            target location
c     ifself               self interaction flag, 
c                            set ifself=1 if the target is on the triangle
c
c     OUTPUT:
c
c     ptfrc(3)            displacement at TARGET
c     strain(3,3)         strain at TARGET
c
c
c
        implicit none

        real(8), intent(in) :: rlam
        real(8), intent(in) :: rmu
        real(8), intent(in) :: triangle(3,3)
        real(8), intent(in) :: sigma_sl(3)
        integer, intent(in) :: ifself
        real(8), intent(in) :: target(3)
        real(8), intent(out) :: ptfrc(3)
        integer, intent(in) :: ifstrain
        real(8), intent(out) :: strain(3,3)

        real(8) :: w(20)
        real(8) :: vert1(3)
        real(8) :: vert2(3)
        real(8) :: vert3(3)
        real(8) :: vectout(3)
        real(8) :: vertout(3)
        real(8) :: rtable(0:2,0:2,0:2)
        real(8) :: btable(0:4,0:4,0:4)

        real(8) :: C1
        real(8) :: C2
        real(8) :: x0
        real(8) :: y0
        real(8) :: z0
        real(8) :: d
        real(8) :: valx
        real(8) :: valy
        real(8) :: valz

        real(8) :: valxx
        real(8) :: valxy
        real(8) :: valxz
        real(8) :: valyx
        real(8) :: valyy
        real(8) :: valyz
        real(8) :: valzx
        real(8) :: valzy
        real(8) :: valzz
        real(8) :: derx
        real(8) :: dery

        real(8) :: derz
        real(8) :: derxx
        real(8) :: deryy
        real(8) :: derzz
        real(8) :: derxy
        real(8) :: derxz
        real(8) :: deryz

        integer :: i
        integer :: j
        integer :: iquad
        integer :: maxb
        integer :: maxr
        integer :: k
c
        C1 = (rlam+rmu)/(rlam+2*rmu)
        C2 = (rmu)/(rlam+2*rmu)
c
        do i=1,3
        ptfrc(i)=0
        enddo
        do i=1,3
        do j=1,3
        strain(i,j)=0
        enddo
        enddo

        call tri_ini(triangle(1,1),triangle(1,2),
     1                triangle(1,3),w,vert1,vert2,vert3)
            ! alias OK because "triangle" arguments are intent(in)
        call tri_for(w,target,vertout)
        x0 = vertout(1)
        y0 = vertout(2)
        z0 = vertout(3)

ccc        call prin2('vertout=*',vertout,3)
       
        call tri_for_vect(w,sigma_sl,vectout)

ccc        call prin2('vectout=*',vectout,3)

        if( ifself .eq. 1 ) iquad=0
c
        if( ifself .ne. 1 ) then
        iquad = 0
        if (z0.gt.0) iquad = +1
        if (z0.lt.0) iquad = -1
        endif


        if( ifstrain .eq. 0 ) maxb=2
        if( ifstrain .eq. 0 ) maxr=0
        if( ifstrain .eq. 1 ) maxb=3
        if( ifstrain .eq. 1 ) maxr=1

        do k=0,maxb
        do i=0,maxb
        do j=0,maxb
        if( i+j+k .le. maxb ) then
        call triabtable(i,j,k,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        btable(i,j,k)=d
        endif
        enddo
        enddo
        enddo
        
        do k=0,maxr
        do i=0,maxr
        do j=0,maxr
        if( i+j+k .le. maxr ) then
        call triartable(i,j,k,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        rtable(i,j,k)=d
        endif
        enddo
        enddo
        enddo


        do k=1,3

        if ( k .eq. 1 ) then
        
        d=btable(2,0,0)
        valx = -d*C1

        d=btable(1,1,0)
        valy = -d*C1

        d=btable(1,0,1)
        valz = -d*C1

        d=rtable(0,0,0)
        valx = valx+d*2

        if( ifstrain .eq. 1 ) then
        d=btable(3,0,0)
        valxx = -d*C1

        d=btable(2,1,0)
        valxy = -d*C1

        d=btable(2,0,1)
        valxz = -d*C1

        d=btable(2,1,0)
        valyx = -d*C1

        d=btable(1,2,0)
        valyy = -d*C1

        d=btable(1,1,1)
        valyz = -d*C1

        d=btable(2,0,1)
        valzx = -d*C1

        d=btable(1,1,1)
        valzy = -d*C1

        d=btable(1,0,2)
        valzz = -d*C1

        d=rtable(1,0,0)
        valxx = valxx+d*2

        d=rtable(0,1,0)
        valxy = valxy+d*2

        d=rtable(0,0,1)
        valxz = valxz+d*2
        endif

        endif

        if ( k .eq. 2 ) then
        
        d=btable(1,1,0)
        valx = -d*C1

        d=btable(0,2,0)
        valy = -d*C1

        d=btable(0,1,1)
        valz = -d*C1

        d=rtable(0,0,0)
        valy = valy+d*2

        if( ifstrain .eq. 1 ) then
        d=btable(2,1,0)
        valxx = -d*C1

        d=btable(1,2,0)
        valxy = -d*C1

        d=btable(1,1,1)
        valxz = -d*C1

        d=btable(1,2,0)
        valyx = -d*C1

        d=btable(0,3,0)
        valyy = -d*C1

        d=btable(0,2,1)
        valyz = -d*C1

        d=btable(1,1,1)
        valzx = -d*C1

        d=btable(0,2,1)
        valzy = -d*C1

        d=btable(0,1,2)
        valzz = -d*C1

        d=rtable(1,0,0)
        valyx = valyx+d*2

        d=rtable(0,1,0)
        valyy = valyy+d*2

        d=rtable(0,0,1)
        valyz = valyz+d*2
        endif

        endif

        if ( k .eq. 3 ) then
        
        d=btable(1,0,1)
        valx = -d*C1

        d=btable(0,1,1)
        valy = -d*C1

        d=btable(0,0,2)
        valz = -d*C1

        d=rtable(0,0,0)
        valz = valz+d*2

        if( ifstrain .eq. 1 ) then
        d=btable(2,0,1)
        valxx = -d*C1

        d=btable(1,1,1)
        valxy = -d*C1

        d=btable(1,0,2)
        valxz = -d*C1

        d=btable(1,1,1)
        valyx = -d*C1

        d=btable(0,2,1)
        valyy = -d*C1

        d=btable(0,1,2)
        valyz = -d*C1

        d=btable(1,0,2)
        valzx = -d*C1

        d=btable(0,1,2)
        valzy = -d*C1

        d=btable(0,0,3)
        valzz = -d*C1

        d=rtable(1,0,0)
        valzx = valzx+d*2

        d=rtable(0,1,0)
        valzy = valzy+d*2

        d=rtable(0,0,1)
        valzz = valzz+d*2
        endif

        endif

        call rotder3d(w,triangle,valx,valy,valz,derx,dery,derz)
        ptfrc(1)=ptfrc(1)+vectout(k)*derx
        ptfrc(2)=ptfrc(2)+vectout(k)*dery
        ptfrc(3)=ptfrc(3)+vectout(k)*derz
c
c
        if( ifstrain .eq. 1 ) then
c
c       ... symmetrize the derivative matrix, prepare to compute strain
c        
        valxy=(valxy+valyx)/2
        valxz=(valxz+valzx)/2
        valyz=(valyz+valzy)/2

        call rothess3d(w,triangle,
     $      valxx,valyy,valzz,valxy,valxz,valyz,
     $      derxx,deryy,derzz,derxy,derxz,deryz)

        strain(1,1)=strain(1,1)+vectout(k)*derxx
        strain(1,2)=strain(1,2)+vectout(k)*derxy
        strain(1,3)=strain(1,3)+vectout(k)*derxz
        strain(2,2)=strain(2,2)+vectout(k)*deryy
        strain(2,3)=strain(2,3)+vectout(k)*deryz
        strain(3,3)=strain(3,3)+vectout(k)*derzz
        
        strain(2,1)=strain(1,2)
        strain(3,1)=strain(1,3)
        strain(3,2)=strain(2,3)

        endif

        enddo
c
        do i=1,3
        ptfrc(i)=ptfrc(i)/(2*rmu)
        enddo
c
c
        if( ifstrain .eq. 1 ) then
c
        do i=1,3
        do j=1,3
        strain(i,j)=-strain(i,j)
        enddo
        enddo
c
        do i=1,3
        do j=1,3
        strain(i,j)=strain(i,j)/(2*rmu)
        enddo
        enddo
c
        endif
c
        return
        end
c
c
c
c
c
        subroutine elust3triadirecttarg_one_slow
     $     (rlam,rmu,triangle,sigma_sl,ifself,target,ptfrc,strain)
c
c
c     Direct evaluation of displacement and strain due to constant
c     single layer elastostatic kernel on a flat triangle.
c
c     Double layer elastostatic kernel: constant-density on flat triangles
c
c     Computes displacement and strain at arbitrary point TARGET
c     due to piecewise-constant single layer density on a triangle.
c
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     rlam,rmu             Lame parameters
c     ntri                 number of triangles
c     sigma_sl(3)          SLP strengths (constant)
c     triangle(3,3)        vertices of the triangle in standard format
c     target(3)            target location
c     ifself               self interaction flag, 
c                            set ifself=1 if the target is on the triangle
c
c     OUTPUT:
c
c     ptfrc(3)            displacement at TARGET
c     strain(3,3)         strain at TARGET
c
c
c
        implicit none

        real(8), intent(in) :: rlam
        real(8), intent(in) :: rmu
        real(8), intent(in) :: triangle(3,3)
        real(8), intent(in) :: sigma_sl(3)
        integer, intent(in) :: ifself
        real(8), intent(in) :: target(3)
        real(8), intent(out) :: ptfrc(3)
        real(8), intent(out) :: strain(3,3)

        real(8) :: w(20)
        real(8) :: vert1(3)
        real(8) :: vert2(3)
        real(8) :: vert3(3)
        real(8) :: vectout(3)
        real(8) :: vertout(3)

        real(8) :: C1
        real(8) :: C2
        real(8) :: x0
        real(8) :: y0
        real(8) :: z0
        real(8) :: d
        real(8) :: valx
        real(8) :: valy
        real(8) :: valz
        real(8) :: valxx
        real(8) :: valxy
        real(8) :: valxz
        real(8) :: valyx
        real(8) :: valyy
        real(8) :: valyz
        real(8) :: valzx
        real(8) :: valzy
        real(8) :: valzz
        real(8) :: derx
        real(8) :: dery
        real(8) :: derz
        real(8) :: derxx
        real(8) :: deryy
        real(8) :: derzz
        real(8) :: derxy
        real(8) :: derxz
        real(8) :: deryz

        integer :: i
        integer :: j
        integer :: iquad
        integer :: k
c
        C1 = (rlam+rmu)/(rlam+2*rmu)
        C2 = (rmu)/(rlam+2*rmu)
c
        do i=1,3
        ptfrc(i)=0
        enddo
        do i=1,3
        do j=1,3
        strain(i,j)=0
        enddo
        enddo

        call tri_ini(triangle(1,1),triangle(1,2),
     1                triangle(1,3),w,vert1,vert2,vert3)
        call tri_for(w,target,vertout)
        x0 = vertout(1)
        y0 = vertout(2)
        z0 = vertout(3)

ccc        call prin2('vertout=*',vertout,3)
       
        call tri_for_vect(w,sigma_sl,vectout)

ccc        call prin2('vectout=*',vectout,3)

        if( ifself .eq. 1 ) iquad=0
c
        if( ifself .ne. 1 ) then
        iquad = 0
        if (z0.gt.0) iquad = +1
        if (z0.lt.0) iquad = -1
        endif

        do k=1,3

        if ( k .eq. 1 ) then
        
        call triabtable(2,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valx = -d*C1

        call triabtable(1,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valy = -d*C1

        call triabtable(1,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valz = -d*C1

        call triartable(0,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valx = valx+d*2


        call triabtable(3,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = -d*C1

        call triabtable(2,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = -d*C1

        call triabtable(2,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = -d*C1

        call triabtable(2,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = -d*C1

        call triabtable(1,2,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = -d*C1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = -d*C1

        call triabtable(2,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = -d*C1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = -d*C1

        call triabtable(1,0,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = -d*C1

        call triartable(1,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = valxx+d*2

        call triartable(0,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = valxy+d*2

        call triartable(0,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = valxz+d*2

        endif

        if ( k .eq. 2 ) then
        
        call triabtable(1,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valx = -d*C1

        call triabtable(0,2,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valy = -d*C1

        call triabtable(0,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valz = -d*C1

        call triartable(0,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valy = valy+d*2


        call triabtable(2,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = -d*C1

        call triabtable(1,2,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = -d*C1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = -d*C1

        call triabtable(1,2,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = -d*C1

        call triabtable(0,3,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = -d*C1

        call triabtable(0,2,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = -d*C1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = -d*C1

        call triabtable(0,2,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = -d*C1

        call triabtable(0,1,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = -d*C1

        call triartable(1,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = valyx+d*2

        call triartable(0,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = valyy+d*2

        call triartable(0,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = valyz+d*2

        endif

        if ( k .eq. 3 ) then
        
        call triabtable(1,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valx = -d*C1

        call triabtable(0,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valy = -d*C1

        call triabtable(0,0,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valz = -d*C1

        call triartable(0,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valz = valz+d*2


        call triabtable(2,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = -d*C1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = -d*C1

        call triabtable(1,0,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = -d*C1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = -d*C1

        call triabtable(0,2,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = -d*C1

        call triabtable(0,1,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = -d*C1

        call triabtable(1,0,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = -d*C1

        call triabtable(0,1,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = -d*C1

        call triabtable(0,0,3,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = -d*C1

        call triartable(1,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = valzx+d*2

        call triartable(0,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = valzy+d*2

        call triartable(0,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = valzz+d*2


        endif

        call rotder3d(w,triangle,valx,valy,valz,derx,dery,derz)
        ptfrc(1)=ptfrc(1)+vectout(k)*derx
        ptfrc(2)=ptfrc(2)+vectout(k)*dery
        ptfrc(3)=ptfrc(3)+vectout(k)*derz
c
c
c       ... symmetrize the derivative matrix, prepare to compute strain
c        
        valxy=(valxy+valyx)/2
        valxz=(valxz+valzx)/2
        valyz=(valyz+valzy)/2

        call rothess3d(w,triangle,
     $      valxx,valyy,valzz,valxy,valxz,valyz,
     $      derxx,deryy,derzz,derxy,derxz,deryz)

        strain(1,1)=strain(1,1)+vectout(k)*derxx
        strain(1,2)=strain(1,2)+vectout(k)*derxy
        strain(1,3)=strain(1,3)+vectout(k)*derxz
        strain(2,2)=strain(2,2)+vectout(k)*deryy
        strain(2,3)=strain(2,3)+vectout(k)*deryz
        strain(3,3)=strain(3,3)+vectout(k)*derzz
        
        strain(2,1)=strain(1,2)
        strain(3,1)=strain(1,3)
        strain(3,2)=strain(2,3)

        enddo
c
c
c
        do i=1,3
        do j=1,3
        strain(i,j)=-strain(i,j)
        enddo
        enddo
c
c
        do i=1,3
        ptfrc(i)=ptfrc(i)/(2*rmu)
        enddo

        do i=1,3
        do j=1,3
        strain(i,j)=strain(i,j)/(2*rmu)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine elust3triadirecttarg
     $     (rlam,rmu,ntri,triangles,sigma_sl,
     1     target,ptfrc,ifstrain,strain)
c
c     Single layer elastostatic kernel: constant-densities on flat triangles
c
c     Computes displacement and strain at arbitrary point TARGET not lying
c     on the surface due to piecewise-constant single layer density on
c     collection of triangles.
c
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     rlam,rmu             Lame parameters
c     ntri                 number of triangles
c     sigma_sl(3,ntri)     array of SLP strengths (constant)
c     triangles(3,3,ntri)  array of triangles in standard format
c     target(3)            target location
c
c     OUTPUT:
c
c     ptfrc(3)            displacement at TARGET
c     strain(3,3)         strain at TARGET
c
c
        implicit none

        real(8), intent(in) :: rlam
        real(8), intent(in) :: rmu
        integer, intent(in) :: ntri
        real(8), intent(in) :: triangles(3,3,ntri)
        real(8), intent(in) :: sigma_sl(3,ntri)
        real(8), intent(in) :: target(3)
        real(8), intent(out) :: ptfrc(3)
        integer, intent(in) :: ifstrain
        real(8), intent(out) :: strain(3,3)

        real(8) :: ptfrc0(3)
        real(8) :: strain0(3,3)

        integer :: i
        integer :: j
        integer :: k
        integer :: ifself


        do i=1,3
        ptfrc(i)=0
        enddo
        do i=1,3
        do j=1,3
        strain(i,j)=0
        enddo
        enddo

        do k=1,ntri

        ifself=0
        call elust3triadirecttarg_one
     $     (rlam,rmu,triangles(1,1,k),sigma_sl(1,k),
     1     ifself,target,ptfrc0,ifstrain,strain0)

        do i=1,3
        ptfrc(i)=ptfrc(i)+ptfrc0(i)
        enddo
        do i=1,3
        do j=1,3
        strain(i,j)=strain(i,j)+strain0(i,j)
        enddo
        enddo

        enddo

        return
        end
c
c
c
c
c
        subroutine eltst3triadirecttarg_one
     $     (rlam,rmu,triangle,sigma_dl,trinorm,
     $     ifself,target,ptfrc,ifstrain,strain)
c
c
c     Direct evaluation of displacement and strain due to constant
c     double layer elastostatic kernel on a flat triangle.
c
c     Double layer elastostatic kernel: constant-density on flat triangles
c
c     Computes displacement and strain at arbitrary point TARGET 
c     due to piecewise-constant double layer density on a triangle.
c
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     rlam,rmu             Lame parameters
c     sigma_dl(3)          DLP strengths (constant)
c     triangle(3,3)        vertices of the triangle in standard format
c     trianorm(3)          triangle normal
c     target(3)            target location
c     ifself               self interaction flag, 
c                            set ifself=1 if the target is on the triangle
c
c     OUTPUT:
c
c     ptfrc(3)            displacement at TARGET
c     strain(3,3)         strain at TARGET
c
c
c
        implicit none

        real(8), intent(in) :: rlam
        real(8), intent(in) :: rmu
        real(8), intent(in) :: triangle(3,3)
        real(8), intent(in) :: sigma_dl(3)
        real(8), intent(in) :: trinorm(3)  ! not used
        integer, intent(in) :: ifself
        real(8), intent(in) :: target(3)
        real(8), intent(out) :: ptfrc(3)
        integer , intent(in):: ifstrain
        real(8), intent(out) :: strain(3,3)

        real(8) :: w(20)
        real(8) :: vert1(3)
        real(8) :: vert2(3)
        real(8) :: vert3(3)
        real(8) :: vectout(3)
        real(8) :: vertout(3)
        real(8) :: tmatr(3,3)
        real(8) :: umatr(3,3,3)
        real(8) :: smatr(3,3,3)
        real(8) :: dmatr(3,3,3,3)
        real(8) :: rtable(0:2,0:2,0:2)
        real(8) :: btable(0:4,0:4,0:4)

        real(8) :: C1
        real(8) :: C2
        real(8) :: x0
        real(8) :: y0
        real(8) :: z0
        real(8) :: d

        real(8) :: valxx
        real(8) :: valxy
        real(8) :: valxz
        real(8) :: valyx
        real(8) :: valyy
        real(8) :: valyz
        real(8) :: valzx
        real(8) :: valzy
        real(8) :: valzz
        real(8) :: valx
        real(8) :: valy
        real(8) :: valz

        real(8) :: derx
        real(8) :: dery
        real(8) :: derz
        real(8) :: derxx
        real(8) :: deryy
        real(8) :: derzz
        real(8) :: derxy
        real(8) :: derxz
        real(8) :: deryz

        integer :: i
        integer :: j
        integer :: iquad
        integer :: k
        integer :: m
        integer :: maxb
        integer :: maxr
        integer :: ix
        integer :: iy
        integer :: iz
c
        C1 = (rlam+rmu)/(rlam+2*rmu)
        C2 = (rmu)/(rlam+2*rmu)
c
        do i=1,3
        ptfrc(i)=0
        enddo
        do i=1,3
        do j=1,3
        strain(i,j)=0
        enddo
        enddo

        call tri_ini(triangle(1,1),triangle(1,2),
     1                triangle(1,3),w,vert1,vert2,vert3)
            ! alias OK because "triangle" arguments are intent(in)
        call tri_for(w,target,vertout)
        x0 = vertout(1)
        y0 = vertout(2)
        z0 = vertout(3)

ccc        call prin2('vertout=*',vertout,3)
       
ccc        call tri_for_vect(w,trinorm,vectout)
ccc        call prin2('inside eltst3triadirecttarg, trinorm=*',vectout,3)

        call tri_for_vect(w,sigma_dl,vectout)

ccc        call prin2('vectout=*',vectout,3)

        if( ifself .eq. 1 ) iquad=0

        if( ifself .ne. 1 ) then
        iquad = 0
        if (z0.gt.0) iquad = +1
        if (z0.lt.0) iquad = -1
        endif

        do k=1,3
        do i=1,3
        do j=1,3
        umatr(i,j,k)=0
        enddo
        enddo
        enddo
        
        do m=1,3
        do k=1,3
        do i=1,3
        do j=1,3
        dmatr(i,j,k,m)=0
        enddo
        enddo
        enddo
        enddo


        if( ifstrain .eq. 0 ) maxb=3
        if( ifstrain .eq. 0 ) maxr=1
        if( ifstrain .eq. 1 ) maxb=4
        if( ifstrain .eq. 1 ) maxr=2

        do k=0,maxb
        do i=0,maxb
        do j=0,maxb
        if( i+j+k .le. maxb ) then
        call triabtable(i,j,k,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        btable(i,j,k)=d
        endif
        enddo
        enddo
        enddo
        
        do k=0,maxr
        do i=0,maxr
        do j=0,maxr
        if( i+j+k .le. maxr ) then
        call triartable(i,j,k,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        rtable(i,j,k)=d
        endif
        enddo
        enddo
        enddo


        
        do k=1,3

        if ( k .eq. 1 ) then
        
        d=btable(3,0,0)
        valxx = -d*C1

        d=btable(2,1,0)
        valxy = -d*C1

        d=btable(2,0,1)
        valxz = -d*C1

        d=btable(2,1,0)
        valyx = -d*C1

        d=btable(1,2,0)
        valyy = -d*C1

        d=btable(1,1,1)
        valyz = -d*C1

        d=btable(2,0,1)
        valzx = -d*C1

        d=btable(1,1,1)
        valzy = -d*C1

        d=btable(1,0,2)
        valzz = -d*C1

        d=rtable(1,0,0)
        valxx = valxx+d*2

        d=rtable(0,1,0)
        valxy = valxy+d*2

        d=rtable(0,0,1)
        valxz = valxz+d*2

        endif

        if ( k .eq. 2 ) then
        
        d=btable(2,1,0)
        valxx = -d*C1

        d=btable(1,2,0)
        valxy = -d*C1

        d=btable(1,1,1)
        valxz = -d*C1

        d=btable(1,2,0)
        valyx = -d*C1

        d=btable(0,3,0)
        valyy = -d*C1

        d=btable(0,2,1)
        valyz = -d*C1

        d=btable(1,1,1)
        valzx = -d*C1

        d=btable(0,2,1)
        valzy = -d*C1

        d=btable(0,1,2)
        valzz = -d*C1

        d=rtable(1,0,0)
        valyx = valyx+d*2

        d=rtable(0,1,0)
        valyy = valyy+d*2

        d=rtable(0,0,1)
        valyz = valyz+d*2

        endif

        if ( k .eq. 3 ) then
        
        d=btable(2,0,1)
        valxx = -d*C1

        d=btable(1,1,1)
        valxy = -d*C1

        d=btable(1,0,2)
        valxz = -d*C1

        d=btable(1,1,1)
        valyx = -d*C1

        d=btable(0,2,1)
        valyy = -d*C1

        d=btable(0,1,2)
        valyz = -d*C1

        d=btable(1,0,2)
        valzx = -d*C1

        d=btable(0,1,2)
        valzy = -d*C1

        d=btable(0,0,3)
        valzz = -d*C1

        d=rtable(1,0,0)
        valzx = valzx+d*2

        d=rtable(0,1,0)
        valzy = valzy+d*2

        d=rtable(0,0,1)
        valzz = valzz+d*2

        endif

        umatr(1,k,1)=valxx
        umatr(1,k,2)=valxy
        umatr(1,k,3)=valxz
        umatr(2,k,1)=valyx
        umatr(2,k,2)=valyy
        umatr(2,k,3)=valyz
        umatr(3,k,1)=valzx
        umatr(3,k,2)=valzy
        umatr(3,k,3)=valzz
c
        enddo
c
c
        if( ifstrain .eq. 1 ) then
c       ... and the derivatives 
c
        do m=1,3
c
        if( m .eq. 1 ) then
        ix=1
        iy=0
        iz=0
        endif
        if( m .eq. 2 ) then
        ix=0
        iy=1
        iz=0
        endif
        if( m .eq. 3 ) then
        ix=0
        iy=0
        iz=1
        endif

        do k=1,3

        if ( k .eq. 1 ) then
        
        d=btable(3+ix,0+iy,0+iz)
        valxx = -d*C1

        d=btable(2+ix,1+iy,0+iz)
        valxy = -d*C1

        d=btable(2+ix,0+iy,1+iz)
        valxz = -d*C1

        d=btable(2+ix,1+iy,0+iz)
        valyx = -d*C1

        d=btable(1+ix,2+iy,0+iz)
        valyy = -d*C1

        d=btable(1+ix,1+iy,1+iz)
        valyz = -d*C1

        d=btable(2+ix,0+iy,1+iz)
        valzx = -d*C1

        d=btable(1+ix,1+iy,1+iz)
        valzy = -d*C1

        d=btable(1+ix,0+iy,2+iz)
        valzz = -d*C1

        d=rtable(1+ix,0+iy,0+iz)
        valxx = valxx+d*2

        d=rtable(0+ix,1+iy,0+iz)
        valxy = valxy+d*2

        d=rtable(0+ix,0+iy,1+iz)
        valxz = valxz+d*2

        endif

        if ( k .eq. 2 ) then
        
        d=btable(2+ix,1+iy,0+iz)
        valxx = -d*C1

        d=btable(1+ix,2+iy,0+iz)
        valxy = -d*C1

        d=btable(1+ix,1+iy,1+iz)
        valxz = -d*C1
        
        d=btable(1+ix,2+iy,0+iz)
        valyx = -d*C1

        d=btable(0+ix,3+iy,0+iz)
        valyy = -d*C1

        d=btable(0+ix,2+iy,1+iz)
        valyz = -d*C1

        d=btable(1+ix,1+iy,1+iz)
        valzx = -d*C1

        d=btable(0+ix,2+iy,1+iz)
        valzy = -d*C1

        d=btable(0+ix,1+iy,2+iz)
        valzz = -d*C1

        d=rtable(1+ix,0+iy,0+iz)
        valyx = valyx+d*2

        d=rtable(0+ix,1+iy,0+iz)
        valyy = valyy+d*2

        d=rtable(0+ix,0+iy,1+iz)
        valyz = valyz+d*2

        endif

        if ( k .eq. 3 ) then
        
        d=btable(2+ix,0+iy,1+iz)
        valxx = -d*C1

        d=btable(1+ix,1+iy,1+iz)
        valxy = -d*C1

        d=btable(1+ix,0+iy,2+iz)
        valxz = -d*C1

        d=btable(1+ix,1+iy,1+iz)
        valyx = -d*C1

        d=btable(0+ix,2+iy,1+iz)
        valyy = -d*C1

        d=btable(0+ix,1+iy,2+iz)
        valyz = -d*C1

        d=btable(1+ix,0+iy,2+iz)
        valzx = -d*C1

        d=btable(0+ix,1+iy,2+iz)
        valzy = -d*C1

        d=btable(0+ix,0+iy,3+iz)
        valzz = -d*C1

        d=rtable(1+ix,0+iy,0+iz)
        valzx = valzx+d*2

        d=rtable(0+ix,1+iy,0+iz)
        valzy = valzy+d*2

        d=rtable(0+ix,0+iy,1+iz)
        valzz = valzz+d*2

        endif

        dmatr(1,k,1,m)=valxx
        dmatr(1,k,2,m)=valxy
        dmatr(1,k,3,m)=valxz
        dmatr(2,k,1,m)=valyx
        dmatr(2,k,2,m)=valyy
        dmatr(2,k,3,m)=valyz
        dmatr(3,k,1,m)=valzx
        dmatr(3,k,2,m)=valzy
        dmatr(3,k,3,m)=valzz
c
        enddo
        enddo
c
ccc        call prin2('umatr=*',umatr,3*3*3)
c
        endif
c
c
        do i=1,3
        do j=1,3
        tmatr(i,j)=0
        enddo
        enddo
c
        do i=1,3
        tmatr(i,3)=rlam*(umatr(1,i,1)+umatr(2,i,2)+umatr(3,i,3)) 
        enddo
        do i=1,3
        do j=1,3
        tmatr(i,j)=tmatr(i,j)+rmu*(umatr(j,i,3)+umatr(3,i,j))
        enddo
        enddo
c
ccc        call prin2('tmatr=*',tmatr,3*3)
c
        if( ifstrain .eq. 1 ) then
c
        do m=1,3
c
        do i=1,3
        do j=1,3
        smatr(i,j,m)=0
        enddo
        enddo
c
        do i=1,3
        smatr(i,3,m)=
     $     rlam*(dmatr(1,i,1,m)+dmatr(2,i,2,m)+dmatr(3,i,3,m)) 
        enddo
        do i=1,3
        do j=1,3
        smatr(i,j,m)=smatr(i,j,m)+rmu*(dmatr(j,i,3,m)+dmatr(3,i,j,m))
        enddo
        enddo
c
        enddo
c
ccc        call prin2('smatr=*',tmatr,3*3*3)
c
        endif
c
        do k=1,3

        valx=tmatr(1,k)
        valy=tmatr(2,k)
        valz=tmatr(3,k)

        call rotder3d(w,triangle,valx,valy,valz,derx,dery,derz)
        ptfrc(1)=ptfrc(1)+vectout(k)*derx
        ptfrc(2)=ptfrc(2)+vectout(k)*dery
        ptfrc(3)=ptfrc(3)+vectout(k)*derz
c
c
        if( ifstrain .eq. 1 ) then
c
c       ... symmetrize the derivative matrix, prepare to compute strain
c        
        valxx=smatr(1,k,1)
        valyx=smatr(2,k,1)
        valzx=smatr(3,k,1)
        valxy=smatr(1,k,2)
        valyy=smatr(2,k,2)
        valzy=smatr(3,k,2)
        valxz=smatr(1,k,3)
        valyz=smatr(2,k,3)
        valzz=smatr(3,k,3)

        valxy=(valxy+valyx)/2
        valxz=(valxz+valzx)/2
        valyz=(valyz+valzy)/2

        call rothess3d(w,triangle,
     $      valxx,valyy,valzz,valxy,valxz,valyz,
     $      derxx,deryy,derzz,derxy,derxz,deryz)

        strain(1,1)=strain(1,1)+vectout(k)*derxx
        strain(1,2)=strain(1,2)+vectout(k)*derxy
        strain(1,3)=strain(1,3)+vectout(k)*derxz
        strain(2,2)=strain(2,2)+vectout(k)*deryy
        strain(2,3)=strain(2,3)+vectout(k)*deryz
        strain(3,3)=strain(3,3)+vectout(k)*derzz
        
        strain(2,1)=strain(1,2) 
        strain(3,1)=strain(1,3)
        strain(3,2)=strain(2,3)
c
        endif
c
        enddo
c
c
c
        do i=1,3
        ptfrc(i)=ptfrc(i)/(2*rmu)
        enddo
c
        if( ifstrain .eq. 1 ) then
c
        do i=1,3
        do j=1,3
        strain(i,j)=-strain(i,j)
        enddo
        enddo
c
        do i=1,3
        do j=1,3
        strain(i,j)=strain(i,j)/(2*rmu)
        enddo
        enddo
c
        endif
c
        return
        end
c
c
c
c
c
        subroutine eltst3triadirecttarg_one_slow
     $     (rlam,rmu,triangle,sigma_dl,trinorm,
     $     ifself,target,ptfrc,strain)
c
c
c     Direct evaluation of displacement and strain due to constant
c     double layer elastostatic kernel on a flat triangle.
c
c     Double layer elastostatic kernel: constant-density on flat triangles
c
c     Computes displacement and strain at arbitrary point TARGET 
c     due to piecewise-constant double layer density on a triangle.
c
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     rlam,rmu             Lame parameters
c     ntri                 number of triangles
c     sigma_dl(3)          DLP strengths (constant)
c     triangle(3,3)        vertices of the triangle in standard format
c     trianorm(3)          triangle normal
c     target(3)            target location
c     ifself               self interaction flag, 
c                            set ifself=1 if the target is on the triangle
c
c     OUTPUT:
c
c     ptfrc(3)            displacement at TARGET
c     strain(3,3)         strain at TARGET
c
c
c
        implicit none

        real(8), intent(in) :: rlam
        real(8), intent(in) :: rmu
        real(8), intent(in) :: triangle(3,3)
        real(8), intent(in) :: sigma_dl(3)
        real(8), intent(in) :: trinorm(3)   ! not used
        integer, intent(in) :: ifself
        real(8), intent(in) :: target(3)
        real(8), intent(out) :: ptfrc(3)
        real(8), intent(out) :: strain(3,3)

        real(8) :: w(20)
        real(8) :: vert1(3)
        real(8) :: vert2(3)
        real(8) :: vert3(3)
        real(8) :: vectout(3)
        real(8) :: vertout(3)
        real(8) :: tmatr(3,3)
        real(8) :: umatr(3,3,3)
        real(8) :: smatr(3,3,3)
        real(8) :: dmatr(3,3,3,3)

        real(8) :: C1
        real(8) :: C2
        real(8) :: x0
        real(8) :: y0
        real(8) :: z0
        real(8) :: d
        real(8) :: valxx
        real(8) :: valxy
        real(8) :: valxz
        real(8) :: valyx
        real(8) :: valyy
        real(8) :: valyz
        real(8) :: valzx
        real(8) :: valzy
        real(8) :: valzz
        real(8) :: valx
        real(8) :: valy
        real(8) :: valz
        real(8) :: derx
        real(8) :: dery
        real(8) :: derz
        real(8) :: derxx
        real(8) :: deryy
        real(8) :: derzz
        real(8) :: derxy
        real(8) :: derxz
        real(8) :: deryz

        integer :: i
        integer :: j
        integer :: iquad
        integer :: k
        integer :: m
        integer :: ix
        integer :: iy
        integer :: iz
c
        C1 = (rlam+rmu)/(rlam+2*rmu)
        C2 = (rmu)/(rlam+2*rmu)
c
        do i=1,3
        ptfrc(i)=0
        enddo
        do i=1,3
        do j=1,3
        strain(i,j)=0
        enddo
        enddo

        call tri_ini(triangle(1,1),triangle(1,2),
     1                triangle(1,3),w,vert1,vert2,vert3)
        call tri_for(w,target,vertout)
        x0 = vertout(1)
        y0 = vertout(2)
        z0 = vertout(3)

ccc        call prin2('vertout=*',vertout,3)
       
ccc        call tri_for_vect(w,trinorm,vectout)
ccc        call prin2('inside eltst3triadirecttarg, trinorm=*',vectout,3)

        call tri_for_vect(w,sigma_dl,vectout)

ccc        call prin2('vectout=*',vectout,3)

        if( ifself .eq. 1 ) iquad=0

        if( ifself .ne. 1 ) then
        iquad = 0
        if (z0.gt.0) iquad = +1
        if (z0.lt.0) iquad = -1
        endif

        do k=1,3
        do i=1,3
        do j=1,3
        umatr(i,j,k)=0
        enddo
        enddo
        enddo
        
        do m=1,3
        do k=1,3
        do i=1,3
        do j=1,3
        dmatr(i,j,k,m)=0
        enddo
        enddo
        enddo
        enddo
        
        do k=1,3

        if ( k .eq. 1 ) then
        
        call triabtable(3,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = -d*C1

        call triabtable(2,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = -d*C1

        call triabtable(2,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = -d*C1

        call triabtable(2,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = -d*C1

        call triabtable(1,2,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = -d*C1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = -d*C1

        call triabtable(2,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = -d*C1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = -d*C1

        call triabtable(1,0,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = -d*C1

        call triartable(1,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = valxx+d*2

        call triartable(0,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = valxy+d*2

        call triartable(0,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = valxz+d*2

        endif

        if ( k .eq. 2 ) then
        
        call triabtable(2,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = -d*C1

        call triabtable(1,2,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = -d*C1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = -d*C1

        call triabtable(1,2,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = -d*C1

        call triabtable(0,3,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = -d*C1

        call triabtable(0,2,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = -d*C1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = -d*C1

        call triabtable(0,2,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = -d*C1

        call triabtable(0,1,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = -d*C1

        call triartable(1,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = valyx+d*2

        call triartable(0,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = valyy+d*2

        call triartable(0,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = valyz+d*2

        endif

        if ( k .eq. 3 ) then
        
        call triabtable(2,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = -d*C1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = -d*C1

        call triabtable(1,0,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = -d*C1

        call triabtable(1,1,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = -d*C1

        call triabtable(0,2,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = -d*C1

        call triabtable(0,1,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = -d*C1

        call triabtable(1,0,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = -d*C1

        call triabtable(0,1,2,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = -d*C1

        call triabtable(0,0,3,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = -d*C1

        call triartable(1,0,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = valzx+d*2

        call triartable(0,1,0,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = valzy+d*2

        call triartable(0,0,1,iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = valzz+d*2

        endif

        umatr(1,k,1)=valxx
        umatr(1,k,2)=valxy
        umatr(1,k,3)=valxz
        umatr(2,k,1)=valyx
        umatr(2,k,2)=valyy
        umatr(2,k,3)=valyz
        umatr(3,k,1)=valzx
        umatr(3,k,2)=valzy
        umatr(3,k,3)=valzz
c
        enddo
c
c
c       ... and the derivatives 
c
        do m=1,3
c
        if( m .eq. 1 ) then
        ix=1
        iy=0
        iz=0
        endif
        if( m .eq. 2 ) then
        ix=0
        iy=1
        iz=0
        endif
        if( m .eq. 3 ) then
        ix=0
        iy=0
        iz=1
        endif

        do k=1,3

        if ( k .eq. 1 ) then
        
        call triabtable(3+ix,0+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = -d*C1

        call triabtable(2+ix,1+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = -d*C1

        call triabtable(2+ix,0+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = -d*C1

        call triabtable(2+ix,1+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = -d*C1

        call triabtable(1+ix,2+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = -d*C1

        call triabtable(1+ix,1+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = -d*C1

        call triabtable(2+ix,0+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = -d*C1

        call triabtable(1+ix,1+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = -d*C1

        call triabtable(1+ix,0+iy,2+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = -d*C1

        call triartable(1+ix,0+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = valxx+d*2

        call triartable(0+ix,1+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = valxy+d*2

        call triartable(0+ix,0+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = valxz+d*2

        endif

        if ( k .eq. 2 ) then
        
        call triabtable(2+ix,1+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = -d*C1

        call triabtable(1+ix,2+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = -d*C1

        call triabtable(1+ix,1+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = -d*C1

        call triabtable(1+ix,2+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = -d*C1

        call triabtable(0+ix,3+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = -d*C1

        call triabtable(0+ix,2+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = -d*C1

        call triabtable(1+ix,1+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = -d*C1

        call triabtable(0+ix,2+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = -d*C1

        call triabtable(0+ix,1+iy,2+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = -d*C1

        call triartable(1+ix,0+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = valyx+d*2

        call triartable(0+ix,1+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = valyy+d*2

        call triartable(0+ix,0+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = valyz+d*2

        endif

        if ( k .eq. 3 ) then
        
        call triabtable(2+ix,0+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxx = -d*C1

        call triabtable(1+ix,1+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxy = -d*C1

        call triabtable(1+ix,0+iy,2+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valxz = -d*C1

        call triabtable(1+ix,1+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyx = -d*C1

        call triabtable(0+ix,2+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyy = -d*C1

        call triabtable(0+ix,1+iy,2+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valyz = -d*C1

        call triabtable(1+ix,0+iy,2+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = -d*C1

        call triabtable(0+ix,1+iy,2+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = -d*C1

        call triabtable(0+ix,0+iy,3+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = -d*C1

        call triartable(1+ix,0+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzx = valzx+d*2

        call triartable(0+ix,1+iy,0+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzy = valzy+d*2

        call triartable(0+ix,0+iy,1+iz,
     $     iquad,vert1,vert2,vert3,x0,y0,z0,d)
        valzz = valzz+d*2

        endif

        dmatr(1,k,1,m)=valxx
        dmatr(1,k,2,m)=valxy
        dmatr(1,k,3,m)=valxz
        dmatr(2,k,1,m)=valyx
        dmatr(2,k,2,m)=valyy
        dmatr(2,k,3,m)=valyz
        dmatr(3,k,1,m)=valzx
        dmatr(3,k,2,m)=valzy
        dmatr(3,k,3,m)=valzz
c
        enddo
        enddo
c
c
ccc        call prin2('umatr=*',umatr,3*3*3)
c
c
        do i=1,3
        do j=1,3
        tmatr(i,j)=0
        enddo
        enddo
c
        do i=1,3
        tmatr(i,3)=rlam*(umatr(1,i,1)+umatr(2,i,2)+umatr(3,i,3)) 
        enddo
        do i=1,3
        do j=1,3
        tmatr(i,j)=tmatr(i,j)+rmu*(umatr(j,i,3)+umatr(3,i,j))
        enddo
        enddo
c
ccc        call prin2('tmatr=*',tmatr,3*3)
c
        do m=1,3
c
        do i=1,3
        do j=1,3
        smatr(i,j,m)=0
        enddo
        enddo
c
        do i=1,3
        smatr(i,3,m)=
     $     rlam*(dmatr(1,i,1,m)+dmatr(2,i,2,m)+dmatr(3,i,3,m)) 
        enddo
        do i=1,3
        do j=1,3
        smatr(i,j,m)=smatr(i,j,m)+rmu*(dmatr(j,i,3,m)+dmatr(3,i,j,m))
        enddo
        enddo
c
        enddo
c
ccc        call prin2('smatr=*',tmatr,3*3*3)
c
c
        do k=1,3

        valx=tmatr(1,k)
        valy=tmatr(2,k)
        valz=tmatr(3,k)

        call rotder3d(w,triangle,valx,valy,valz,derx,dery,derz)
        ptfrc(1)=ptfrc(1)+vectout(k)*derx
        ptfrc(2)=ptfrc(2)+vectout(k)*dery
        ptfrc(3)=ptfrc(3)+vectout(k)*derz
c
c
c       ... symmetrize the derivative matrix, prepare to compute strain
c        
        valxx=smatr(1,k,1)
        valyx=smatr(2,k,1)
        valzx=smatr(3,k,1)
        valxy=smatr(1,k,2)
        valyy=smatr(2,k,2)
        valzy=smatr(3,k,2)
        valxz=smatr(1,k,3)
        valyz=smatr(2,k,3)
        valzz=smatr(3,k,3)

        valxy=(valxy+valyx)/2
        valxz=(valxz+valzx)/2
        valyz=(valyz+valzy)/2

        call rothess3d(w,triangle,
     $      valxx,valyy,valzz,valxy,valxz,valyz,
     $      derxx,deryy,derzz,derxy,derxz,deryz)

        strain(1,1)=strain(1,1)+vectout(k)*derxx
        strain(1,2)=strain(1,2)+vectout(k)*derxy
        strain(1,3)=strain(1,3)+vectout(k)*derxz
        strain(2,2)=strain(2,2)+vectout(k)*deryy
        strain(2,3)=strain(2,3)+vectout(k)*deryz
        strain(3,3)=strain(3,3)+vectout(k)*derzz
        
        strain(2,1)=strain(1,2) 
        strain(3,1)=strain(1,3)
        strain(3,2)=strain(2,3)

        enddo
c
c
c
        do i=1,3
        do j=1,3
        strain(i,j)=-strain(i,j)
        enddo
        enddo
c
        do i=1,3
        ptfrc(i)=ptfrc(i)/(2*rmu)
        enddo

        do i=1,3
        do j=1,3
        strain(i,j)=strain(i,j)/(2*rmu)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine eltst3triadirecttarg
     $     (rlam,rmu,ntri,triangles,sigma_dl,trinorm,
     1     target,ptfrc,ifstrain,strain)
c
c     Double layer elastostatic kernel: constant-densities on flat triangles
c
c     Computes displacement and strain at arbitrary point TARGET not lying
c     on the surface due to piecewise-constant double layer density on
c     collection of triangles.
c
c     Analytic quadratures are used (see triahquad.f).
c
c     INPUT:
c
c     rlam,rmu             Lame parameters
c     ntri                 number of triangles
c     sigma_dl(3,ntri)     array of DLP strengths (constant)
c     triangles(3,3,ntri)  array of triangles in standard format
c     trianorm(3,ntri)     array of triangle normals
c     target(3)            target location
c
c     OUTPUT:
c
c     ptfrc(3)            displacement at TARGET
c     strain(3,3)         strain at TARGET
c
c
        implicit none

        real(8), intent(in) :: rlam
        real(8), intent(in) :: rmu
        integer, intent(in) :: ntri
        real(8), intent(in) :: triangles(3,3,ntri)
        real(8), intent(in) :: sigma_dl(3,ntri)
        real(8), intent(in) :: trinorm(3,ntri)
        real(8), intent(in) :: target(3)
        real(8), intent(out) :: ptfrc(3)
        integer, intent(in) :: ifstrain
        real(8), intent(out) :: strain(3,3)

        real(8) :: ptfrc0(3)
        real(8) :: strain0(3,3)

        integer :: i
        integer :: j
        integer :: k
        integer :: ifself
c
        do i=1,3
        ptfrc(i)=0
        enddo
        do i=1,3
        do j=1,3
        strain(i,j)=0
        enddo
        enddo

        do k=1,ntri

        ifself=0
        call eltst3triadirecttarg_one
     $     (rlam,rmu,triangles(1,1,k),sigma_dl(1,k),trinorm(1,k),
     1     ifself,target,ptfrc0,ifstrain,strain0)

        do i=1,3
        ptfrc(i)=ptfrc(i)+ptfrc0(i)
        enddo
        if( ifstrain .eq. 1 ) then
        do i=1,3
        do j=1,3
        strain(i,j)=strain(i,j)+strain0(i,j)
        enddo
        enddo
        endif
        enddo

        return
        end
c
c
c
c
c
        subroutine elust3triadirectself
     $     (rlam,rmu,ipatch,ntri,triangles,sigma_sl,
     1     zparts,ptfrc,ifstrain,strain)
c
c     Single layer elastostatic kernel: constant-densities on flat triangles
c
c     Computes displacement and strain at centroid zparts(*,ipatch) 
c     on the surface due to piecewise-constant double layer density on
c     collection of triangles, numbered jpatch = 1,...,ntri. 
c
c     If ipatch equals jpatch, they are assumed to be the same triangle
c     and a singular quadrature rule is used. Otherwise, they are assumed 
c     to be distinct. In either case, analytic quadratures are used
c     (see triahquad.f)
c
c     INPUT:
c
c     rlam,rmu             Lame parameters
c     ntri                 number of triangles
c     sigma_sl(3,ntri)     array of SLP strengths (constant)
c     triangles(3,3,ntri)  array of triangles in standard format
c     zparts(3,1)          array of triangle centroids
c
c     OUTPUT:
c
c     ptfrc(3)            displacement at centroid zparts(*,ipatch)
c     strain(3,3)         strain at centroid zparts(*,ipatch)
c
c     Note: This routine is always called with ipatch=1 and ntri=1.
c
        implicit none

        real(8), intent(in) :: rlam
        real(8), intent(in) :: rmu
        integer, intent(in) :: ipatch
        integer, intent(in) :: ntri
        real(8), intent(in) :: triangles(3,3,ntri)
        real(8), intent(in) :: sigma_sl(3,ntri)
        real(8), intent(in) :: zparts(3,1)
        real(8), intent(out) :: ptfrc(3)
        integer, intent(in) :: ifstrain
        real(8), intent(out) :: strain(3,3)

        real(8) :: ptfrc0(3)
        real(8) :: strain0(3,3)

        integer :: i
        integer :: j
        integer :: k
        integer :: ifself
c
        do i=1,3
        ptfrc(i)=0
        enddo
        do i=1,3
        do j=1,3
        strain(i,j)=0
        enddo
        enddo

        do k=1,ntri

        if( k .eq. ipatch ) ifself=1
        if( k .ne. ipatch ) ifself=0

        call elust3triadirecttarg_one
     $     (rlam,rmu,triangles(1,1,k),sigma_sl(1,k),
     1     ifself,zparts(1,ipatch),ptfrc0,ifstrain,strain0)

        do i=1,3
        ptfrc(i)=ptfrc(i)+ptfrc0(i)
        enddo
        if( ifstrain .eq. 1 ) then
        do i=1,3
        do j=1,3
        strain(i,j)=strain(i,j)+strain0(i,j)
        enddo
        enddo
        endif
        enddo

        return
        end
c
c
c
c
c
        subroutine eltst3triadirectself
     $     (rlam,rmu,ipatch,ntri,triangles,sigma_dl,trinorm,
     1     zparts,ptfrc,ifstrain,strain)
c
c     Double layer elastostatic kernel: constant-densities on flat triangles
c
c     Computes displacement and strain at centroid zparts(*,ipatch) 
c     on the surface due to piecewise-constant double layer density on
c     collection of triangles, numbered jpatch = 1,...,ntri. 
c
c     If ipatch equals jpatch, they are assumed to be the same triangle
c     and a singular quadrature rule is used. Otherwise, they are assumed 
c     to be distinct. In either case, analytic quadratures are used
c     (see triahquad.f)
c
c     INPUT:
c
c     rlam,rmu             Lame parameters
c     ntri                 number of triangles
c     sigma_dl(3,ntri)     array of DLP strengths (constant)
c     triangles(3,3,ntri)  array of triangles in standard format
c     trianorm(3,ntri)     array of triangle normals
c     zparts(3,1)          array of triangle centroids
c
c     OUTPUT:
c
c     ptfrc(3)            displacement at centroid zparts(*,ipatch)
c     strain(3,3)         strain at centroid zparts(*,ipatch)
c
c     Note: This routine is always called with ipatch=1 and ntri=1.
c
        implicit none

        real(8), intent(in) :: rlam
        real(8), intent(in) :: rmu
        integer, intent(in) :: ipatch
        integer, intent(in) :: ntri
        real(8), intent(in) :: triangles(3,3,ntri)
        real(8), intent(in) :: sigma_dl(3,ntri)
        real(8), intent(in) :: trinorm(3,ntri)
        real(8), intent(in) :: zparts(3,1)
        real(8), intent(out) :: ptfrc(3)
        integer, intent(in) :: ifstrain
        real(8), intent(out) :: strain(3,3)

        real(8) :: ptfrc0(3)
        real(8) :: strain0(3,3)

        integer :: i
        integer :: j
        integer :: k
        integer :: ifself
c
        do i=1,3
        ptfrc(i)=0
        enddo
        do i=1,3
        do j=1,3
        strain(i,j)=0
        enddo
        enddo
        
        do k=1,ntri

        if( k .eq. ipatch ) ifself=1
        if( k .ne. ipatch ) ifself=0

        call eltst3triadirecttarg_one
     $     (rlam,rmu,triangles(1,1,k),sigma_dl(1,k),trinorm(1,k),
     1     ifself,zparts(1,ipatch),ptfrc0,ifstrain,strain0)

        do i=1,3
        ptfrc(i)=ptfrc(i)+ptfrc0(i)
        enddo
        if( ifstrain .eq. 1 ) then
        do i=1,3
        do j=1,3
        strain(i,j)=strain(i,j)+strain0(i,j)
        enddo
        enddo
        endif
        enddo

        return
        end










c
c
c
c
c
