c*******************************************************************************
c*******************************************************************************
	double precision function bessj(m,x)
	implicit none
c
c	First implemented in Potsdam, Feb, 1999
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
	integer m
	double precision x
c
c	J_n(x), n <= 3
c
	double precision eps
	parameter(eps=1.0d-10)
	integer i,m0
	double precision ax,bx,y0,y1,y,sum
c
	double precision bessj0,bessj1
c
	m0=iabs(m)
	ax=dabs(x)
	if(m0.eq.0)then
	  bessj=bessj0(ax)
	else if(m0.eq.1)then
	  bessj=bessj1(ax)
	else
	  if(ax.gt.dble(m0))then
	    bx=2.d0/ax
	    y0=bessj0(ax)
	    y1=bessj1(ax)
	    do i=1,m0-1
	      y=dble(i)*bx*y1-y0
	      y0=y1
	      y1=y
	    enddo
	    bessj=y
	  else if(ax.eq.0.d0)then
	    bessj=0.d0
	  else
	    y=1.d0
	    bx=0.5d0*ax
	    do i=1,m0
	      y=y*bx/dble(i)
	    enddo
	    sum=y
	    bx=bx*bx
	    i=0
100	    i=i+1
	    y=-y*bx/dble(i)/dble(m0+i)
	    sum=sum+y
	    if(dabs(y).gt.eps*dabs(sum).or.i.lt.10)goto 100
	    bessj=sum
	  endif
	endif
	if(dble(m)*x.lt.0.d0.and.mod(m0,2).eq.1)bessj=-bessj
	return
	end
