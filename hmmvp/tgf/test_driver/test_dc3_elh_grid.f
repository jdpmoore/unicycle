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
c  This is a sample driver for the triangle Green's function code.
c  
c  It computes the Green's function for a triangular dislocation in an
c  elastic half-space, and compares the result to direct numerical
c  integration of the point dislocation kernel over the triangle.
c  
c  This driver constructs a "triangles" array which contains a single
c  triangle.  It constructs a "target" array which consists of a set of
c  target locations, on the surface of a rectangular box surrounding the
c  triangle, in a 12 x 10 x 10 grid.  The code computes the displacement
c  and strain tensor at each target point, for a uniform dislocation on
c  the triangle.  The output includes results computed using both the
c  Green's function and direct numerical integration of the point
c  displacement kernel, and the differences between the two.
c  
c  The entire test is repeated for 7 dip angles from -90 to +90 degrees.
c
c
        subroutine test_dc3_elh_grid ()

        implicit none

        real(8) :: source(3)
        real(8) :: target(3)
        real(8) :: du(3)
        real(8) :: rnorm(3)
        real(8) :: fvec(3)
        real(8) :: strain(3,3)
        real(8) :: triangles(3,3,100)
        real(8) :: trinorm(3,100)
        real(8) :: sigma(3,100)
        real(8) :: sigma_sl(3,100)
        real(8) :: sigma_dl(3,100)
        real(8) :: centroids(3,3)
        real(8) :: fvec0(3)
        real(8) :: ptfrc0(3)
        real(8) :: strain0(3,3)

        real(8) :: rlam
        real(8) :: rmu
        real(8) :: rnu
        real(8) :: alpha
        real(8) :: er1max
        real(8) :: er2max
        real(8) :: depth
        real(8) :: dip
        real(8) :: phi
        real(8) :: pot1
        real(8) :: pot2
        real(8) :: pot3
        real(8) :: pot4
        real(8) :: scale
        real(8) :: x
        real(8) :: y
        real(8) :: z
        real(8) :: sx
        real(8) :: sy
        real(8) :: sz
        real(8) :: er1
        real(8) :: em1
        real(8) :: ev1
        real(8) :: er2
        real(8) :: em2
        real(8) :: ev2

        integer :: idip
        integer :: itri
        integer :: ntri
        integer :: i
        integer :: j
        integer :: ifsingle
        integer :: ifdouble
        integer :: ngridx
        integer :: ngridy
        integer :: ngridz
        integer :: ix
        integer :: iy
        integer :: iz
        integer :: ifstrain
        integer :: ifptfrc
        integer :: ifptfrctarg
        integer :: ifstraintarg
        integer :: ntarget
        integer :: numfunev

        real(8), parameter :: pi = 3.1415926535897932384626433832795d0
c
        rlam = 3.45d0
        rmu = 2.3d0
        rnu = rlam/2/(rlam+rmu)

        alpha = (rlam+rmu)/(rlam+2*rmu)

        er1max=0
        er2max=0

        do 1400 idip=-90,90,30

        depth = 1.50d0
        dip = dble(idip)
        phi = dip/180.0d0*pi

        pot1 = -2
        pot2 = 3
        pot3 = 1
        pot4 = 0

c
        source(1)=0
        source(2)=0
        source(3)=-depth

c
c       convert strike-dip-tensile parameters into cartesian coordinates
c
        du(1)=pot1
        du(2)=pot2*cos(phi)-pot3*sin(phi)
        du(3)=pot2*sin(phi)+pot3*cos(phi)
        rnorm(1)=0
        rnorm(2)=-sin(phi)
        rnorm(3)=+cos(phi)

        scale=1d0
        du(1)=du(1)*scale
        du(2)=du(2)*scale
        du(3)=du(3)*scale

        call prin2('du=*',du,3)
        call prin2('rnorm=*',rnorm,3)

        itri=1
        triangles(1,1,itri) = 0 
        triangles(2,1,itri) = 0
        triangles(3,1,itri) = -depth
        triangles(1,2,itri) = 1
        triangles(2,2,itri) = 0
        triangles(3,2,itri) = -depth
        triangles(1,3,itri) = 1
        triangles(2,3,itri) = cos(phi)
        triangles(3,3,itri) = -depth+sin(phi)

        itri=2
        triangles(1,1,itri) = 0
        triangles(2,1,itri) = 0
        triangles(3,1,itri) = -depth
        triangles(1,2,itri) = 1
        triangles(2,2,itri) = cos(phi)
        triangles(3,2,itri) = -depth+sin(phi)
        triangles(1,3,itri) = 0
        triangles(2,3,itri) = cos(phi)
        triangles(3,3,itri) = -depth+sin(phi)

        ntri = 1

        do i = 1,ntri
        sigma(1,i)=du(1)
        sigma(2,i)=du(2)
        sigma(3,i)=du(3)
        enddo

        call triangle_norm(triangles(1,1,1),trinorm(1,1))
        call triangle_norm(triangles(1,1,2),trinorm(1,2))
        call prin2('triangles=*',triangles,3*3*ntri)
        call prin2('trinorm=*',trinorm,3*ntri)


        ifsingle=0
        ifdouble=1
        do i = 1,ntri
        sigma_sl(1,i)=0
        sigma_sl(2,i)=0
        sigma_sl(3,i)=0
        sigma_dl(1,i)=du(1)
        sigma_dl(2,i)=du(2)
        sigma_dl(3,i)=du(3)
        enddo
c
        do i = 1,ntri
        centroids(1,i)=
     $     (triangles(1,1,i)+triangles(1,2,i)+triangles(1,3,i))/3
        centroids(2,i)=
     $     (triangles(2,1,i)+triangles(2,2,i)+triangles(2,3,i))/3
        centroids(3,i)=
     $     (triangles(3,1,i)+triangles(3,2,i)+triangles(3,3,i))/3
        enddo


        write(15,2030) 'dip=', dip
        write(15,*) 'triangle=['
        write(15,2010) triangles(1,1,1),triangles(2,1,1),
     $     triangles(3,1,1)
        write(15,2010) triangles(1,2,1),triangles(2,2,1),
     $     triangles(3,2,1)
        write(15,2010) triangles(1,3,1),triangles(2,3,1),
     $     triangles(3,3,1)
        write(15,*) ']'
        
        write(15,*) 'trinorm=['
        write(15,2010) trinorm(1,1),trinorm(2,1),
     $     trinorm(3,1)
        write(15,*) ']'

        write(15,*) 'sigma=['
        write(15,2010) du(1),du(2),
     $     du(3)
        write(15,*) ']'
        
        write(15,2020) 'ifsingle=', ifsingle
        write(15,2020) 'ifdouble=', ifdouble
        write(15,2020) 'ifhalfspace=', 1

 2010   format(3(1x,1pe13.5))
 2020   format(A,1x,I3)
 2030   format(A,1x,1pe13.5)

        x=-1.0d0 
        y=-2.0d0
        z= 0.0d0

        ngridx=12
        ngridy=10
        ngridz=10

        sx= 3.0d0
        sy= 4.0d0
        sz=-3.0d0

        write(15,*) 'errors=['

        do 1300 ix=0,ngridx
        do 1200 iy=0,ngridy
        do 1100 iz=0,ngridz

        if (      ix /= 0 .and. ix /= ngridx
     &      .and. iy /= 0 .and. iy /= ngridy
     &      .and. iz /= 0 .and. iz /= ngridz ) then
          cycle
        endif

        target(1)=x+dble(ix)/dble(ngridx)*sx
        target(2)=y+dble(iy)/dble(ngridy)*sy
        target(3)=z+dble(iz)/dble(ngridz)*sz

        ifstrain = 1


        call princ(' *')
        call prin2('target=*',target,3)
        call princ('=== Triangles (analytic+adaptive) ===*')

        
        ifptfrc=0
        ifstrain=0
        ifptfrctarg=1
        ifstraintarg=1
        ntarget=1
        call elh3dtriadirecttarg(
     $     RLAM,RMU,TRIANGLES,TRINORM,NTRI,CENTROIDS,
     $     ifsingle,SIGMA_SL,ifdouble,SIGMA_DL,
     $     ifptfrc,ptfrc0,ifstrain,strain0,NTARGET,
     $     target,ifptfrctarg,fvec,
     $     ifstraintarg,STRAIN)

        do i = 1,3
        fvec(i)=fvec(i)/(4*pi)
        enddo
        do i = 1,3
        do j = 1,3
        strain(i,j)=strain(i,j)/(4*pi)
        enddo
        enddo

        call prin2('after elh3dtriadirecttarg, fvec=*',fvec,3)
        call prin2('after elh3dtriadirecttarg, strain=*',strain,3*3)

        do i=1,3
        fvec0(i)=fvec(i)
        do j=1,3
        strain0(i,j)=strain(i,j)
        enddo
        enddo

        call princ(' *')
        call princ('=== Triangles (adaptive integration) ===*')

        call elthb3triaadap
     $     (rlam,rmu,ntri,triangles,sigma,trinorm,
     1     target,fvec,strain,numfunev)
        call prinfs('numfunev=*',numfunev)

        do i = 1,3
        fvec(i)=fvec(i)/(4*pi)
        enddo
        do i = 1,3
        do j = 1,3
        strain(i,j)=strain(i,j)/(4*pi)
        enddo
        enddo

        call prin2('after elthb3triaadap, fvec=*',fvec,3)
        call prin2('after elthb3triaadap, strain=*',strain,3*3)
c
c       
        er1=0
        em1=0
        ev1=0
        do i=1,3
        er1=er1+(fvec0(i)-fvec(i))**2
        em1=em1+(fvec(i))**2
        ev1=ev1+(fvec0(i))**2
        enddo
        er1=sqrt(er1/em1)
        er1max=max(er1max,er1)
        ev1=sqrt(ev1)
        call prin2s('relative error in displacement=*',er1)

        er2=0
        em2=0
        ev2=0
        do i=1,3
        do j=1,3
        er2=er2+(strain0(i,j)-strain(i,j))**2
        em2=em2+(strain(i,j))**2
        ev2=ev2+(strain0(i,j))**2
        enddo
        enddo
        er2=sqrt(er2/em2)
        er2max=max(er2max,er2)
        ev2=sqrt(ev2)
        call prin2s('relative error in strain=*',er2)

        write(15,2000) target(1),target(2),target(3),ev1,ev2,er1,er2
 2000   format(7(1x,1pe13.5))

 1100   continue
 1200   continue
 1300   continue

        write(15,*) ']'
        write(15,*)
        write(15,*)
        write(15,*)
        write(15,*)

 1400   continue

        write(15,2030) 'max_displacement_error=', er1max
        write(15,2030) 'max_strain_error=', er2max
        write(15,*)

        return
        end
c
c
c
