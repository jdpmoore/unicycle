	subroutine edgsublay(resolut,itty,ierr)
	implicit none
c
c	First implemented in Potsdam, Feb, 1999
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
	integer itty,ierr
	double precision resolut(3)
c
	include 'edgglobal.h'
c       
c       model parameter:
c       n0: number of homogeneous layers
c
        integer n0
        double precision h(lmax),ro(lmax),vp(lmax),vs(lmax)
        common /model/ h,ro,vp,vs,n0
c
c	original model parameters
c
	integer l0
        double precision z1(lmax),z2(lmax),ro1(lmax),ro2(lmax)
	double precision vp1(lmax),vp2(lmax),vs1(lmax),vs2(lmax)
	common /model0/z1,z2,ro1,ro2,vp1,vp2,vs1,vs2,l0
c
c	work space
c
	double precision eps
	parameter(eps=1.0d-02)
	integer i,i0,l
	double precision dh,dro,dvp,dvs,z,dz
c
	n0=0
	do l=1,l0-1
	  dz=z2(l)-z1(l)
	  dvp=2.d0*dabs(vp2(l)-vp1(l))/(vp2(l)+vp1(l))
	  dvs=2.d0*dabs(vs2(l)-vs1(l))/(vs2(l)+vs1(l))
	  dro=2.d0*dabs(ro2(l)-ro1(l))/(ro2(l)+ro1(l))
	  i0=idnint(dmax1(dro/resolut(1),
     &       dvp/resolut(2),dvs/resolut(3)))
	  i0=max0(1,i0)
	  dro=(ro2(l)-ro1(l))/dz
	  dvp=(vp2(l)-vp1(l))/dz
	  dvs=(vs2(l)-vs1(l))/dz
	  dh=dz/dble(i0)
	  do i=1,i0
	    n0=n0+1
	    if(n0.ge.lmax)then
	      ierr=1
	      return
	    endif
	    h(n0)=dh
	    z=(dble(i)-0.5d0)*dh
	    ro(n0)=ro1(l)+dro*z
	    vp(n0)=vp1(l)+dvp*z
	    vs(n0)=vs1(l)+dvs*z
	  enddo
	enddo
c
c	last layer is halfspace
c
	n0=n0+1
	h(n0)=0.d0
	ro(n0)=ro1(l0)
	vp(n0)=vp1(l0)
	vs(n0)=vs1(l0)
c
	if(itty.eq.1)then
	  write(*,'(7a)')' no ',' thick(m)  ','   vp(m/s) ',
     &    '  vs(m/s)  ',' ro(kg/m^3)'
	  do i=1,n0
	    write(*,1001)i,h(i),vp(i),vs(i),ro(i)
	  enddo
	endif
1001	format(i4,4f11.4)
	ierr=0
	return
	end
