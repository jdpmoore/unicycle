	subroutine edgmoment(istype,strength,ro,vp,vs)
	implicit none
c
c	First implemented in Potsdam, Feb, 1999
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
	double precision pi,pi2
	parameter(pi=3.14159265358979d0,pi2=6.28318530717959d0)
c
	integer istype
	double precision ro,vp,vs
	double precision strength
c
c       source parameters
c
        integer ls,ms,ics
        integer kpower(6)
        double precision zs,r0
        double precision sfct(6)
        common /source/ zs,r0,sfct,ls,ms,ics,kpower
c
	integer i
c
	do i=1,6
	  sfct(i)=0.d0
	  kpower(i)=0
	enddo
c
	if(istype.eq.0)then
c
c	  explosion source (m11=m22=m33=M0)
c
	  ms=0
	  ics=1
	  sfct(1)=-strength/(2.d0*pi*ro*vp*vp)
	  sfct(4)=-strength*(vs/vp)**2/pi
	  kpower(4)=1
	else if(istype.eq.1)then
c
c	  strike-slip (m12=m21=M0)
c
	  ms=2
	  ics=-1
	  sfct(4)=strength/(2.d0*pi)
	  sfct(6)=-sfct(4)
	  kpower(4)=1
	  kpower(6)=1
	else if(istype.eq.2)then
c
c	  dip-slip (m13=m31=M0)
c
	  ms=1
	  ics=1
	  sfct(3)=-strength/(2.d0*pi*ro*vs*vs)
	  sfct(5)=sfct(3)
	else if(istype.eq.3)then
c
c	  compensated linear vector dipole (CLVD) (m11=m22=-M0/2, M33=M0)
c
	  ms=0
	  ics=1
	  sfct(1)=-strength/(2.d0*pi*ro*vp*vp)
	  sfct(4)=strength*(3.d0-4.d0*(vs/vp)**2)/(4.d0*pi)
	  kpower(4)=1
	else if(istype.eq.4)then
c
c	  vertical-single-force (fz=F0)
c
	  ms=0
	  ics=1
	  sfct(2)=strength/(2.d0*pi)
	else if(istype.eq.5)then
c
c	  horizontal-single-force (fx=F0)
c
	  ms=1
	  ics=1
	  sfct(4)=strength/(2.d0*pi)
	  sfct(6)=sfct(4)
	else if(istype.eq.6)then
c
c	  airgun in small pool (m11=m22=M0, fz=M0/r_source)
c
	  ms=0
	  ics=1
	  sfct(2)=-strength/(2.d0*pi*r0)
	  sfct(4)=-strength/(2.d0*pi)
	  kpower(4)=1
	else
	  stop ' Error in edgmoment: wrong source type switch no!'
	endif
c
	return
	end


	  
