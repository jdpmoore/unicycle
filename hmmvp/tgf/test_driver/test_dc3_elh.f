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
c  It compares our results to Okada's results.  The comparisions are
c  done for both a point dislocation, and a rectangular dislocation, in
c  an elastic half-space.
c
c  Comparisons are performed at a set of "target" points, which are
c  located on the surface of a rectangular box surrounding the source,
c  in a 7 x 7 x 7 grid.  The code computes the displacement and strain
c  tensor for each target point.  The output includes both sets of
c  results, and the magnitude of the differences between the two.
c
c  For the point dislocation, results are computed using both our code
c  and Okada's code for a point source.
c
c  For the rectangular dislocation, results are computed using both our
c  code an Okada's code for a point source.  Our code constructs a
c  "triangles" array which contains two triangles, which together form
c  the rectangle.  The point displacement kernel is numerically
c  integrated over the two triangles to obtain the result.
c  
c  The entire test is repeated for 7 dip angles from -90 to +90 degrees.
c
c
        subroutine test_dc3_elh ()

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
        real(8) :: fvec0(3)
        real(8) :: strain0(3,3)

        real(8) :: rlam
        real(8) :: rmu
        real(8) :: rnu
        real(8) :: alpha
        real(8) :: per1max
        real(8) :: per2max
        real(8) :: er1max
        real(8) :: er2max
        real(8) :: depth
        real(8) :: dip
        real(8) :: phi
        real(8) :: pot1
        real(8) :: pot2
        real(8) :: pot3
        real(8) :: pot4
        real(8) :: bx
        real(8) :: by
        real(8) :: bz
        real(8) :: sx
        real(8) :: sy
        real(8) :: sz
        real(8) :: x
        real(8) :: y
        real(8) :: z
        real(8) :: ux
        real(8) :: uy
        real(8) :: uz
        real(8) :: uxx
        real(8) :: uxy
        real(8) :: uxz
        real(8) :: uyy
        real(8) :: uyx
        real(8) :: uyz
        real(8) :: uzx
        real(8) :: uzy
        real(8) :: uzz
        real(8) :: per1
        real(8) :: pem1
        real(8) :: pev1
        real(8) :: per2
        real(8) :: pem2
        real(8) :: pev2
        real(8) :: er1
        real(8) :: em1
        real(8) :: ev1
        real(8) :: er2
        real(8) :: em2
        real(8) :: ev2
        real(8) :: al1
        real(8) :: al2
        real(8) :: aw1
        real(8) :: aw2
        real(8) :: disl1
        real(8) :: disl2
        real(8) :: disl3

        integer :: idip
        integer :: itri
        integer :: ntri
        integer :: i
        integer :: ngridx
        integer :: ngridy
        integer :: ngridz
        integer :: ix
        integer :: iy
        integer :: iz
        integer :: IRET
        integer :: j
        integer :: ifstrain
        integer :: numfunev

        real(8), parameter :: pi = 3.1415926535897932384626433832795d0
c
        rlam = 3.45d0
        rmu = 2.3d0
        rnu = rlam/2/(rlam+rmu)

        alpha = (rlam+rmu)/(rlam+2*rmu)

        per1max=0
        per2max=0
        er1max=0
        er2max=0

        do 1400 idip=-90,90,30

        depth = 1.50d0
        dip = dble(idip)
        phi = dip/180.0d0*pi

!!        pot1 = 0
!!        pot2 = 0
!!        pot3 = 1
!!        pot4 = 0

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

        ntri = 2

        do i = 1,ntri
        sigma(1,i)=du(1)
        sigma(2,i)=du(2)
        sigma(3,i)=du(3)
        enddo

        call triangle_norm(triangles(1,1,1),trinorm(1,1))
        call triangle_norm(triangles(1,1,2),trinorm(1,2))
        call prin2('triangles=*',triangles,3*3*ntri)
        call prin2('trinorm=*',trinorm,3*ntri)


        write(15,2030) 'dip=', dip
        write(15,*) 'triangle=['
        do i = 1,ntri
        write(15,2010) triangles(1,1,i),triangles(2,1,i),
     $     triangles(3,1,i)
        write(15,2010) triangles(1,2,i),triangles(2,2,i),
     $     triangles(3,2,i)
        write(15,2010) triangles(1,3,i),triangles(2,3,i),
     $     triangles(3,3,i)
        enddo
        write(15,*) ']'
        
        write(15,*) 'trinorm=['
        do i = 1,ntri
        write(15,2010) trinorm(1,i),trinorm(2,i),
     $     trinorm(3,i)
        enddo
        write(15,*) ']'

        write(15,*) 'sigma=['
        do i = 1,ntri
        write(15,2010) sigma(1,i),sigma(2,i),
     $     sigma(3,i)
        enddo
        write(15,*) ']'

 2010   format(3(1x,1pe13.5))
 2020   format(A,1x,I3)
 2030   format(A,1x,1pe13.5)
 
        bx=-1.0d0 
        by=-2.0d0
        bz=-1.0d-14

        ngridx=7
        ngridy=7
        ngridz=7

!!        ngridx=12
!!        ngridy=10
!!        ngridz=10

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

        x=bx+dble(ix)/dble(ngridx)*sx
        y=by+dble(iy)/dble(ngridy)*sy
        z=bz+dble(iz)/dble(ngridz)*sz




        call princ('=== Point source (Okada) ===*')

        call DC3D0(ALPHA,X,Y,Z,DEPTH,DIP,POT1,POT2,POT3,POT4,  
     $     UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET) 

        fvec(1)=ux
        fvec(2)=uy
        fvec(3)=uz
        strain(1,1)=uxx
        strain(1,2)=(uxy+uyx)/2
        strain(1,3)=(uxz+uzx)/2
        strain(2,1)=(uyx+uxy)/2
        strain(2,2)=uyy
        strain(2,3)=(uyz+uzy)/2
        strain(3,1)=(uzx+uxz)/2
        strain(3,2)=(uzy+uyz)/2
        strain(3,3)=uzz
        call prin2('after dc3d0, fvec=*',fvec,3)
        call prin2('after dc3d0, strain=*',strain,3*3)

        do i=1,3
        fvec0(i)=fvec(i)
        do j=1,3
        strain0(i,j)=strain(i,j)
        enddo
        enddo

        call princ(' *')
        call princ('=== Point source (Mindlin pieces) ===*')


        target(1)=x
        target(2)=y
        target(3)=z

        ifstrain = 1

        call green3elth_eval(rlam,rmu,source,du,rnorm,
     $     target,fvec,ifstrain,strain)

        do i =1,3
        fvec(i)=fvec(i)/(4*pi)
        enddo
        do i =1,3
        do j =1,3
        strain(i,j)=strain(i,j)/(4*pi)
        enddo
        enddo

        call prin2('after green3elth_eval, fvec=*',fvec,3)
        call prin2('after green3elth_eval, strain=*',strain,3*3)
c
c       
        per1=0
        pem1=0
        pev1=0
        do i=1,3
        per1=per1+(fvec0(i)-fvec(i))**2
        pem1=pem1+(fvec(i))**2
        pev1=pev1+(fvec0(i))**2
        enddo
        per1=sqrt(per1/pem1)
        per1max=max(per1max,per1)
        pev1=sqrt(pev1)
        call prin2s('relative error in displacement=*',per1)

        per2=0
        pem2=0
        pev2=0
        do i=1,3
        do j=1,3
        per2=per2+(strain0(i,j)-strain(i,j))**2
        pem2=pem2+(strain(i,j))**2
        pev2=pev2+(strain0(i,j))**2
        enddo
        enddo
        per2=sqrt(per2/pem2)
        per2max=max(per2max,per2)
        pev2=sqrt(pev2)
        call prin2s('relative error in strain=*',per2)
        
        


        call princ(' *')
        call princ('=== Rectangle (Okada) ===*')

        al1 = 0
        al2 = 1
        aw1 = 0
        aw2 = 1
        disl1 = pot1
        disl2 = pot2
        disl3 = pot3

        call  DC3D(ALPHA,X,Y,Z,DEPTH,DIP,                        
     $     AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3,           
     $     UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)

        fvec(1)=ux
        fvec(2)=uy
        fvec(3)=uz
        strain(1,1)=uxx
        strain(1,2)=(uxy+uyx)/2
        strain(1,3)=(uxz+uzx)/2
        strain(2,1)=(uyx+uxy)/2
        strain(2,2)=uyy
        strain(2,3)=(uyz+uzy)/2
        strain(3,1)=(uzx+uxz)/2
        strain(3,2)=(uzy+uyz)/2
        strain(3,3)=uzz

!!        do i = 1,3
!!        fvec(i)=fvec(i)*(4*pi)
!!        enddo
!!        do i = 1,3
!!        do j = 1,3
!!        strain(i,j)=strain(i,j)*(4*pi)
!!        enddo
!!        enddo


        call prin2('after dc3d, fvec=*',fvec,3)
        call prin2('after dc3d, strain=*',strain,3*3)

        do i=1,3
        fvec0(i)=fvec(i)
        do j=1,3
        strain0(i,j)=strain(i,j)
        enddo
        enddo

ccc        return
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

        call prin2('after elth3triaadap, fvec=*',fvec,3)
        call prin2('after elth3triaadap, strain=*',strain,3*3)
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

        write(15,2000) target(1),target(2),target(3),
     &        pev1,pev2,per1,per2,ev1,ev2,er1,er2
 2000   format(11(1x,1pe13.5))
 
 1100   continue
 1200   continue
 1300   continue

        write(15,*) ']'
        write(15,*)
        write(15,*)
        write(15,*)
        write(15,*)

 1400   continue

        write(15,2030) 'max_point_displacement_error=', per1max
        write(15,2030) 'max_point_strain_error=', per2max
        write(15,2030) 'max_rect_displacement_error=', er1max
        write(15,2030) 'max_rect_strain_error=', er2max
        write(15,*)

        return
        end
c
c
