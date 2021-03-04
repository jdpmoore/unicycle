c*******************************************************************************
c*******************************************************************************
        subroutine edgkern(y,k,ps,sh,eps)
        implicit none
c
c	First implemented in Potsdam, Feb, 1999
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
c	calculation of response function in frequency-wavelength domain
c       y(6): solution vector
c       k: wave number
c       eps: relative accuracy
c
        double precision k,eps
        double precision y(6)
	logical ps,sh
c
	include 'edgglobal.h'
c       
c       model parameter:
c       n0: number of model layers
c
        integer n0
        double precision h(lmax),ro(lmax),vp(lmax),vs(lmax)
        common /model/ h,ro,vp,vs,n0
c
	integer i
c
	do i=1,6
	  y(i)=0.d0
	enddo
c
	if(ps)call edgpsv(y,k,eps)
	if(sh)call edgsh(y,k,eps)
	return
	end	  
