cc Copyright (C) 2009: Leslie Greengard and Zydrunas Gimbutas, Vladimir Rokhlin
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
c    $Date: 2010-01-05 12:57:52 -0500 (Tue, 05 Jan 2010) $
c    $Revision: 782 $
c
c
c        Tensor product triangle quadrature rules for smooth functions.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine triagauc(n,vert1,vert2,vert3,rnodes,
     1      weights,ifinit,w)

        implicit none

        integer, intent(in) :: n
        real(8), intent(in) :: vert1(2)
        real(8), intent(in) :: vert2(2)
        real(8), intent(in) :: vert3(2)
        real(8), intent(inout) :: rnodes(2*n*n)
        real(8), intent(inout) :: weights(n*n)
        integer, intent(in) :: ifinit
        real(8), intent(inout) :: w(4*n+5)

        integer :: it
        integer :: lt
        integer :: iwhts
        integer :: lwhts
        integer :: itvert
        integer :: ltvert
        integer :: iwhtsver
        integer :: lwhtsver
c
c       this subroutine constructs a tensor product Gaussian
c       quadrature formula on the triangle in the plane 
c       specified by its vertices. Note that the one-dimensional 
c       Gaussian quadratures are lined up parallel to the side of 
c       the triangle opposite the vertex vert1; the total number 
c       of nodes and weights to be created is n**2. Actually, this 
c       is simply a memory manager for the subroutine triagau0,
c       which does all the work.
c
c                              input parameters:
c
c  n - the number of nodes in each direction. note that the 
c       total number of nodes to be created is n**2, and the 
c       algebraic order of the quadrature is (n-1)*2
c  vert1,vert2,vert3 - the vertices of the triangle on which the 
c       quadrature rule is to be constructed
c
c                              output parameters:
c
c  rnodes - the array of n**2 nodes in the plane (all inside the
c       user-supplied traiangle)
c  weights - the quadrature weights corresponding to the nodes rnodes
c
c                              work arrays:
c
c  w - must be at least 4*n+5 real(8) locations 
c
c
c        . . . allocate memory in the work array
c
        it=1
        lt=n+1
c
        iwhts=it+lt
        lwhts=n+2
c
        itvert=iwhts+lwhts
        ltvert=n+2
c
        iwhtsver=itvert+ltvert
        lwhtsver=n+1
c
c       construct the quadrature formula on the triangle
c
!!        call triagau0(n,vert1,vert2,vert3,rnodes,weights,
!!     1      ifinit,w(it),w(iwhts),w(itvert),w(iwhtsver) )
        call triagau0(n,vert1,vert2,vert3,rnodes,weights,
     1      ifinit,w(it:it+n-1),w(iwhts:iwhts+n-1),
     &      w(itvert:itvert+n-1),w(iwhtsver:iwhtsver+n-1) )
c
        return
        end
c
c
c
c
c
        subroutine triagau0(n,vert1,vert2,vert3,rnodes,weights,
     1      ifinit,t,whts,tvert,whtsvert)

        implicit none

        integer, intent(in) :: n
        real(8), intent(in) :: vert1(2)
        real(8), intent(in) :: vert2(2)
        real(8), intent(in) :: vert3(2)
        real(8), intent(inout) :: rnodes(2,n,n)
        real(8), intent(inout) :: weights(n,n)
        integer, intent(in) :: ifinit
        real(8), intent(inout) :: t(n)
        real(8), intent(inout) :: whts(n)
        real(8), intent(inout) :: tvert(n)
        real(8), intent(inout) :: whtsvert(n)

        real(8) :: w(12)
        real(8) :: z1(2)
        real(8) :: z2(2)
        real(8) :: z3(2)

        real(8) :: b
        real(8) :: u
        real(8) :: v
        real(8) :: d
        real(8) :: a12
        real(8) :: b12
        real(8) :: c12
        real(8) :: a13
        real(8) :: b13
        real(8) :: c13
        real(8) :: y

        real(8) :: x12
        real(8) :: x13
        real(8) :: thor

        integer :: ifwhts
        integer :: i
        integer :: j
c
c       this subroutine constructs a tensor product Gaussian
c       quadrature formula on the triangle in the plane 
c       specified by its vertices. Note that the one-dimensional 
c       Gaussian quadratures are lined up parallel to the side of 
c       the triangle opposite the vertex vert1; the total number 
c       of nodes and weights to be created is n**2.
c
c                              input parameters:
c
c  n - the number of nodes in each direction. note that the 
c       total number of nodes to be created is n**2, and the 
c       algebraic order of the quadrature is (n-1)*2
c  vert1,vert2,vert3 - the vertices of the triangle on which the 
c       quadrature rule is to be constructed
c
c                              output parameters:
c
c  rnodes - the array of n**2 nodes in the plane (all inside the
c       user-supplied traiangle)
c  weights - the quadrature weights corresponding to the nodes rnodes
c
c                              work arrays:
c
c  t,whts,tvert,whtsvert - must be n+1 real(8) locations each
c                          No, only need n real(8) locations each
c
c
c       . . . if required, construct the Gaussian nodes and 
c             weights on the interval [-1,1]
c
        if(ifinit .eq. 1) then
        ifwhts=1
        call legewhts(n,t,whts,ifwhts)
        endif
c
c        construct affine transformation putting one side
c       of the user-specified triangle on the real axis. The side
c       to be put on the real axis is the one opposite the vertex vert1.
c
        call trianini(vert1,vert2,vert3,w)
c
c       construct the discretization of the version of the triangle 
c       moved so that its side is on the real axis
c
        call trianfor(w,vert1,z1)
        call trianfor(w,vert2,z2)
        call trianfor(w,vert3,z3)
c
c       . . . construct the gaussian discretization of the
c             hight of the triangle
c
        b=z1(2)
        u=b/2
        v=b/2
        d=dabs(b/2)
        do 1200 i=1,n
        tvert(i)=u*t(i)+v
        whtsvert(i)=whts(i)*d
 1200 continue
c
c        find the equations of the two non-horizontal sides of
c        the shifted triangle
c
        call trialine(z1(1),z2(1),z1(2),z2(2),a12,b12,c12)
c
        call trialine(z1(1),z3(1),z1(2),z3(2),a13,b13,c13)
c
c        one after another, construct the horizontal discretizations
c
        do 2000 i=1,n
c
c       find the current this horizontal section
c
        y=tvert(i)
        x12=-(b12*y+c12)/a12
        x13=-(b13*y+c13)/a13

cccc        call prin2s('x12=*',x12)
cccc        call prin2s('x13=*',x13)
c       
        d=dabs(x12-x13)/2
c
        u=(x13-x12)/2        
        v=(x13+x12)/2        
c
cccc        call prin2s('second u=*',u)
cccc        call prin2s('second v=*',v)
c
        do 1400 j=1,n
c
        thor=u*t(j)+v
c
        rnodes(1,j,i)=thor
        rnodes(2,j,i)=y
c
        weights(j,i)=d*whts(j)*whtsvert(i)
c
 1400 continue
 2000 continue
c
c       now, move the nodes in the plane back to the original triangle
c
        do 2400 i=1,n
        do 2200 j=1,n
c
        call trianbak(w,rnodes(1,j,i),z3)
        rnodes(1,j,i)=z3(1)        
        rnodes(2,j,i)=z3(2)        
 2200 continue
 2400 continue
c
        return
        end
c
c
c
c
c
        subroutine trialine(x1,x2,y1,y2,a,b,c)
        implicit none

        real(8), intent(in) :: x1
        real(8), intent(in) :: x2
        real(8), intent(in) :: y1
        real(8), intent(in) :: y2
        real(8), intent(out) :: a
        real(8), intent(out) :: b
        real(8), intent(out) :: c

        real(8) :: dx
        real(8) :: dy
        real(8) :: d
c
c        construct the equation of the line passing through the 
c        two points (x1,y1), (x2,y2)
c
        dx=x2-x1
        dy=y2-y1
c
        if(dabs(dx) .lt. dabs(dy)) goto 1200
c
        b=1
        a=-dy/dx
        c=-(a*x1+b*y1)
c
        goto 1400
c
 1200 continue
c
        a=1
        b=-dx/dy
        c=-(a*x1+b*y1)
c
 1400 continue
c
c       . . . normalize the coefficients
c
        d=a**2+b**2
        d=dsqrt(d)
c
        a=a/d
        b=b/d
        c=c/d
c
        return
        end
c
c
c
c
c
        module trianini_mod
        implicit none

        integer, parameter :: ixshift = 1
        integer, parameter :: iyshift = 2
        integer, parameter :: iumat = 3
        integer, parameter :: ivmat = 7

        end module trianini_mod




        subroutine trianini(vert1,vert2,vert3,w)

        use trianini_mod
        implicit none

        real(8), intent(in) :: vert1(2)
        real(8), intent(in) :: vert2(2)
        real(8), intent(in) :: vert3(2)
        real(8), intent(out) :: w(10)

        real(8) :: xshift
        real(8) :: yshift

c
c       this subroutine constructs an affine transformation putting one side
c       of the user-specified triangle on the real axis. The side
c       to be put on the real axis is the one opposite the vertex vert1.
c       this subroutine also constructs the transformation inverse to the 
c       above one. The actual transformations are applied to particular
c       user-specified points by the entries trianfor and trianbak of this
c       subroutine
c
c                              input parameters:
c
c  vert1,vert2,vert3 - the vertices of the triangle one of whose 
c       sides is to be put on the x axis
c
c                              output parameters:
c
c  w - the real(8) array of length 10 containing the transformations;
c       it is to be used by the entries trianfor, trianbak (see below)
c
        call triambld(vert1(1),vert1(2),vert2(1),vert2(2),vert3(1),
     1      vert3(2),xshift,yshift,
     &      w(iumat:iumat+3),w(ivmat:ivmat+3) )
        w(ixshift) = xshift
        w(iyshift) = yshift
c
        return
        end subroutine trianini
c
c
c
c
        subroutine trianfor(w,zin,zout)

        use trianini_mod
        implicit none

        real(8), intent(in) :: w(10)
        real(8), intent(in) :: zin(2)
        real(8), intent(out) :: zout(2)
c
c       this entry applies to the user-specified point zin \in R^2
c       the first of the transformations constructed by the entry 
c       trianini (see above), obtaining the point zout \in R^2.
c       
c
        call triarotf(zin(1),zin(2),w(ixshift),w(iyshift),
     1      w(iumat),zout(1),zout(2) )
c
cccc        call prin2('in trianfor, zout=*',zout,2)

        return
        end subroutine trianfor
c
c
c
c
        subroutine trianbak(w,zin,zout)

        use trianini_mod
        implicit none

        real(8), intent(in) :: w(10)
        real(8), intent(in) :: zin(2)
        real(8), intent(out) :: zout(2)
c
c       this entry applies to the user-specified point zin \in R^2
c       the second of the transformations constructed by the entry 
c       trianini (see above), obtaining the point zout \in R^2.
c       
        call triarotb(zin(1),zin(2),w(ixshift),w(iyshift),
     1      w(ivmat),zout(1),zout(2) )
c
        return
        end subroutine trianbak
c
c
c
c
c
        subroutine triarotf(x,y,xshift,yshift,umat,xout,yout)

        implicit none

        real(8), intent(in) :: x
        real(8), intent(in) :: y
        real(8), intent(in) :: xshift
        real(8), intent(in) :: yshift
        real(8), intent(in) :: umat(2,2)
        real(8), intent(out) :: xout
        real(8), intent(out) :: yout

        real(8) :: xx
        real(8) :: yy
c
c        apply the shift
c
        xx=x+xshift
        yy=y+yshift

cccc         call prin2s('in triarotf, xx=*',xx)
cccc         call prin2s('in triarotf, yy=*',yy)

c
c        apply the rotation
c
        xout=umat(1,1)*xx+umat(1,2)*yy
        yout=umat(2,1)*xx+umat(2,2)*yy
        return
        end
c
c
c
c
c
        subroutine triarotb(x,y,xshift,yshift,vmat,xout,yout)

        implicit none

        real(8), intent(in) :: x
        real(8), intent(in) :: y
        real(8), intent(in) :: xshift
        real(8), intent(in) :: yshift
        real(8), intent(in) :: vmat(2,2)
        real(8), intent(out) :: xout
        real(8), intent(out) :: yout

        real(8) :: xx
        real(8) :: yy

c
c        apply the rotation
c
        xx=vmat(1,1)*x+vmat(1,2)*y
        yy=vmat(2,1)*x+vmat(2,2)*y
c
cccc         call prin2s('in triarotb, xx=*',xx)
cccc         call prin2s('in triarotb, yy=*',yy)
c
c        apply the shift
c
        xout=xx-xshift
        yout=yy-yshift
        return
        end
c
c
c
c
c
        subroutine triambld(x1,y1,x2,y2,x3,y3,
     1      xshift,yshift,umat,vmat)

        implicit none

        real(8), intent(in) :: x1
        real(8), intent(in) :: y1
        real(8), intent(in) :: x2
        real(8), intent(in) :: y2
        real(8), intent(in) :: x3
        real(8), intent(in) :: y3
        real(8), intent(out) :: xshift
        real(8), intent(out) :: yshift
        real(8), intent(out) :: umat(2,2)
        real(8), intent(out) :: vmat(2,2)

        real(8) :: xx1
        real(8) :: xx2
        real(8) :: xx3
        real(8) :: yy1
        real(8) :: yy2
        real(8) :: yy3
        real(8) :: u
        real(8) :: v
        real(8) :: d
c
c        construct a motion in the plane that will put the side of 
c        the triangle opposite the vertex (x1,y1) on th2 x axis
c
c        . . . the translation
c
        xshift=-x2
        yshift=-y2
c
        xx1=x1+xshift
        xx2=x2+xshift
        xx3=x3+xshift
c
        yy1=y1+yshift
        yy2=y2+yshift
        yy3=y3+yshift
c
c        . . . rotation
c
        u=(x3-x2)
        v=y3-y2
        d=dsqrt(u**2+v**2)
        u=u/d
        v=v/d
c
        umat(1,1)=u
        umat(1,2)=v
        umat(2,1)=-v
        umat(2,2)=u
c
        vmat(1,1)=umat(1,1)
        vmat(2,2)=umat(2,2)
        vmat(1,2)=umat(2,1)
        vmat(2,1)=umat(1,2)
c
        return
        end
c
c
c
c
c
        subroutine triaarrm(x,y,n)
        implicit none

        integer, intent(in) :: n
        real(8), intent(in) :: x(n)
        real(8), intent(inout) :: y(n)

        integer :: i
c
        do 1200 i=1,n
        y(i)=x(i)
 1200 continue
        return
        end
c
c
c
c
c
        subroutine ctriaarrm(x,y,n)
        implicit none

        integer, intent(in) :: n
        complex(8), intent(in) :: x(n)
        complex(8), intent(inout) :: y(n)

        integer :: i
c
        do 1200 i=1,n
        y(i)=x(i)
 1200 continue
        return
        end
c
c
c
c
c
        subroutine triaplot(iw,z)
        implicit none

        integer, intent(in) :: iw
        real(8), intent(in) :: z(2,3)
c
 1200 format(2x,e12.5,2x,e12.5)
c
        write(iw,1200) z(1,1),z(2,1)
        write(iw,1200) z(1,2),z(2,2)
        write(iw,1200) z(1,3),z(2,3)
        write(iw,1200) z(1,1),z(2,1)
c
 1400 format(2x)
        write(iw,1400)

        return
        end      
c
c
c
c
c
        subroutine triadiv(z1,z2,z3,zout)
        implicit none

        real(8), intent(in) :: z1(2)
        real(8), intent(in) :: z2(2)
        real(8), intent(in) :: z3(2)
        real(8), intent(out) :: zout(2,3,4)

        real(8) :: w12(2)
        real(8) :: w13(2)
        real(8) :: w23(2)
c
c        this subroutine subdivides the input triangle 
c        with the nodes z1, z2, z3 into four triangles
c
c        . . . construct the middles of the sides
c
        w12(1)=(z1(1)+z2(1))/2
        w12(2)=(z1(2)+z2(2))/2
c
        w13(1)=(z1(1)+z3(1))/2
        w13(2)=(z1(2)+z3(2))/2
c
        w23(1)=(z2(1)+z3(1))/2
        w23(2)=(z2(2)+z3(2))/2
c
c       . . . construct the new triangles
c
        zout(1,1,1)=w12(1)
        zout(2,1,1)=w12(2)
c
        zout(1,2,1)=z1(1)
        zout(2,2,1)=z1(2)
c
        zout(1,3,1)=w13(1)
        zout(2,3,1)=w13(2)
c      
c      
        zout(1,1,2)=w13(1)
        zout(2,1,2)=w13(2)
c
        zout(1,2,2)=z3(1)
        zout(2,2,2)=z3(2)
c
        zout(1,3,2)=w23(1)
        zout(2,3,2)=w23(2)
c      
c
        zout(1,1,3)=w23(1)
        zout(2,1,3)=w23(2)
c
        zout(1,2,3)=z2(1)
        zout(2,2,3)=z2(2)
c
        zout(1,3,3)=w12(1)
        zout(2,3,3)=w12(2)
c      
c
        zout(1,1,4)=w23(1)
        zout(2,1,4)=w23(2)
c
        zout(1,2,4)=w12(1)
        zout(2,2,4)=w12(2)
c
        zout(1,3,4)=w13(1)
        zout(2,3,4)=w13(2)
c      
        return
        end
c
c
c
c
c
