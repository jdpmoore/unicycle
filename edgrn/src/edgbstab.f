	subroutine edgbstab(n)
	implicit none
c
c	First implemented in Potsdam, Feb, 1999
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
	integer n
c
	include 'edgglobal.h'
c
c 	table of J_n(x), dJ_n(x)/dx and n*J_n(x)/x
c	all multiplied by sqrt(x)
c
	integer nnbess,nnbess1
	parameter(nnbess=nbess*ndbess,nnbess1=nnbess+ndbess)
	double precision bsdx,bsfct(0:nnbess1,3)
	common /bessels/ bsdx,bsfct
c
	double precision pi,pi2
	parameter(pi=3.14159265358979d0,pi2=6.28318530717959d0)
	integer i,j
	double precision x,xsqrt,a,b
	double precision bessj0,bessj1,bessj
c
	do j=1,3
	  bsfct(0,j)=0.d0
	enddo
	bsdx=pi2/dble(ndbess)
	if(n.eq.0)then
	  do i=1,nnbess1
	    x=bsdx*dble(i)
	    xsqrt=dsqrt(x)
	    bsfct(i,1)=xsqrt*bessj0(x)
	    bsfct(i,2)=-xsqrt*bessj1(x)
	    bsfct(i,3)=0.d0
	  enddo
	else if(n.eq.1)then
	  do i=1,nnbess1
	    x=bsdx*dble(i)
	    xsqrt=dsqrt(x)
	    bsfct(i,1)=xsqrt*bessj1(x)
	    a=xsqrt*bessj0(x)
	    b=xsqrt*bessj(2,x)
	    bsfct(i,2)=0.5d0*(a-b)
	    bsfct(i,3)=0.5d0*(a+b)
	  enddo
	else if(n.eq.2)then
	  do i=1,nnbess1
	    x=bsdx*dble(i)
	    xsqrt=dsqrt(x)
	    bsfct(i,1)=xsqrt*bessj(2,x)
	    a=xsqrt*bessj1(x)
	    b=xsqrt*bessj(3,x)
	    bsfct(i,2)=0.5d0*(a-b)
	    bsfct(i,3)=0.5d0*(a+b)
	  enddo
	else
	  stop ' Error in edgbstab: check 0<= n <= 2?'
	endif
c
	return
	end
