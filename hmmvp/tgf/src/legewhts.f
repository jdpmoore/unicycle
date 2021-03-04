cc Copyright (C) 2009: Vladimir Rokhlin
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
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        This file contains a set of subroutines for the handling 
c        of Legendre expansions. 
c
c   legewhts:  this subroutine constructs the nodes and the
c        weights of the n-point gaussian quadrature on 
c        the interval [-1,1]
c
c   legepol_sum - evaluates a single Legendre polynomial (together
c         with its derivative) at the user-provided point
c
c***********************************************************************
      subroutine legewhts(n,ts,whts,ifwhts)
c***********************************************************************
      use ifdim_def
      implicit none

      integer, intent(in) :: n
      real(8), intent(inout) :: ts(n)
      integer, intent(in) :: ifwhts
      real(8), intent(inout) :: whts(ifdim(n,ifwhts))

      integer :: i
      integer :: ifout
      integer :: k
      real(8) :: eps
      real(8) :: h
      real(8) :: t
      real(8) :: xk
      real(8) :: delta
      real(8) :: deltold
      real(8) :: der
      real(8) :: pol
      real(8) :: sum

      real(8), parameter :: pi = 3.1415926535897932384626433832795d0
c
c     This subroutine constructs the nodes and the
c     weights of the n-point gaussian quadrature on 
c     the interval [-1,1]
c
c     input parameters:
c
c     n - the number of nodes in the quadrature
c
c     output parameters:
c
c     ts - the nodes of the n-point gaussian quadrature
c     w - the weights of the n-point gaussian quadrature
c
c
c
c     construct the array of initial approximations
c             to the roots of the n-th legendre polynomial
c
      if(n .eq. 0) return

      eps=1.0d-14
      h=pi/(2*n) 
      do i=1,n
         t=(2*i-1)*h
         ts(n-i+1)=dcos(t)
      enddo
c
c     use newton to find all roots of the legendre polynomial
c
      ts(n/2+1)=0
      do i=1,n/2
         xk=ts(i)
         ifout=0
         deltold=1
         do k=1,10
            call legepol_sum(xk,n,pol,der,sum)
            delta=-pol/der
            xk=xk+delta
            if(abs(delta) .lt. eps) ifout=ifout+1
            if(ifout .eq. 3) exit
         enddo
         ts(i)=xk
         ts(n-i+1)=-xk
      enddo
c     
c     construct the weights via the orthogonality relation
c
      if(ifwhts .eq. 0) return
c
      do i=1,(n+1)/2
         call legepol_sum(ts(i),n,pol,der,sum)
         whts(i)=1/sum
         whts(n-i+1)=whts(i)
      enddo
c
      return
      end
c
c
c
c
c
      subroutine legepol_sum(x,n,pol,der,sum)
      implicit none

      real(8), intent(in) :: x
      integer, intent(in) :: n
      real(8), intent(out) :: pol
      real(8), intent(out) :: der
      real(8), intent(out) :: sum

      integer :: k
      real(8) :: pkm1
      real(8) :: pk
      real(8) :: pkp1

      real(8), parameter :: done = 1.0d0
c
c     needs documentation
c
      sum=0 
c
      pkm1=1
      pk=x
      sum=sum+pkm1**2 /2
      sum=sum+pk**2 *(1+done/2)
c
      pk=1
      pkp1=x
c
c     if n=0 or n=1 - exit
c
      if (n .le. 1) then
         sum=0 
         pol=1
         der=0
         sum=sum+pol**2 /2
         if(n .eq. 0) return
         pol=x
         der=1
         sum=sum+pol**2*(1+done/2)
         return
      endif
c
c     n is greater than 1. conduct recursion
c
      do k=1,n-1
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
        sum=sum+pkp1**2*(k+1+done/2)
      enddo
c
c        calculate the derivative
c
      pol=pkp1
      der=n*(x*pkp1-pk)/(x**2-1)
      return
      end
c
