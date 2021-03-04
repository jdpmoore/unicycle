	subroutine edcdisc(ns,nz,z1,z2,dr,dz,nps)
	implicit none
c
c	First implemented in Potsdam, Feb, 1999
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
	integer ns,nz,nps
	double precision z1,z2,dr,dz
c
c	inputs:
c	ns = total number of source rectangles
c	nz,z1,z2 = number of depth samples, start and end depths used
c		in Green's functions
c	dlength, dwidth = grid size for discretisation
c
c	returned outputs:
c	nps = total number of discrete point sources
c	other outputs through common blocks
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	LOCAL CONSTANTS
c	===============
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double precision DEGTORAD
	parameter(DEGTORAD=1.745329252E-02)
c
	include 'edcglobal.h'
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	RECTANGULAR SOURCE PLANES
c	=========================
c
c	(xs,ys,zs) = coordinates of the start point of strike
c	with x = north, y = east, z = downward.
c	all angles in degree.
c	NSMAX = the max. number of source rectangles
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double precision dislocation(NSMAX)
	double precision xs(NSMAX),ys(NSMAX),zs(NSMAX)
	double precision length(NSMAX),width(NSMAX)
	double precision strike(NSMAX),dip(NSMAX),rake(NSMAX)
c
	common/rectangles/dislocation,xs,ys,zs,length,width,
     &                    strike,dip,rake
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	DISTRETE POINT SOURCES
c	======================
c
c	(xs,ys,zs) = coordinates of the discrete point sources
c	with x = north, y = east, z = downward
c	angles in degree.
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double precision pxs(NPSMAX),pys(NPSMAX),pzs(NPSMAX)
	double precision pmoment(5,NPSMAX)
c
	common/pointsources/pxs,pys,pzs,pmoment
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	OBSERVATION POSITIONS AND OBSERVABLES
c	=====================================
c
c	(xrec(i),yrec(i),zrec0)=coordinates of the observation positions
c	(Note that zrec0 is fixed)
c	disp = the 3 displcement vector components: ux,uy,uz
c	strain = the 6 strain tensor components: exx,eyy,ezz,exy,eyz,ezx
c	NRECMAX = the max. number of observation positions
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double precision xrec(NRECMAX),yrec(NRECMAX)
	double precision zrec0
	double precision disp(NRECMAX,3),strain(NRECMAX,6)
	double precision tilt(NRECMAX,2)
c
	common/obsarray/xrec,yrec,zrec0,disp,strain,tilt
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	WARNING STATISTICS
c	==================
c
c	nwarn = total number of warnings
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	integer nwarn
c
	common/warnings/nwarn
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	LOCAL WORK SPACES
c	=================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	integer is,ix,iy,nx,ny
	double precision x,y,dx,dy,st,di,ra,disarea
	double precision dlength,dwidth
	double precision sm(3,3)
c
c       write(*,'(a)')' ... discretise rectangular plane sources:'
	dlength=dr
	dwidth=dmin1(dr,dz)
	nps=0
	do is=1,ns
c
	  st=strike(is)*DEGTORAD
	  di=dip(is)*DEGTORAD
	  ra=rake(is)*DEGTORAD
c
	  sm(1,1)=-dsin(di)*dcos(ra)*dsin(2.d0*st)
     &            -dsin(2.d0*di)*dsin(ra)*(dsin(st))**2
	  sm(2,2)= dsin(di)*dcos(ra)*sin(2.d0*st)
     &            -dsin(2.d0*di)*dsin(ra)*(dcos(st))**2
	  sm(3,3)=-(sm(1,1)+sm(2,2))
	  sm(1,2)= dsin(di)*dcos(ra)*dcos(2.d0*st)
     &            +0.5d0*dsin(2.d0*di)*dsin(ra)*dsin(2.d0*st)
	  sm(2,1)=sm(1,2)
	  sm(1,3)=-dcos(di)*dcos(ra)*dcos(st)
     &            -dcos(2.d0*di)*dsin(ra)*dsin(st)
	  sm(3,1)=sm(1,3)
	  sm(2,3)=-dcos(di)*dcos(ra)*dsin(st)
     &            +dcos(2.d0*di)*dsin(ra)*dcos(st)
	  sm(3,2)=sm(2,3)
c
	  nx=max0(1,idnint(length(is)/dlength))
	  ny=max0(1,idnint(width(is)/dwidth))
	  dx=length(is)/dble(nx)
	  dy=width(is)/dble(ny)
c
c	  if one of length and width = 0, then it is a line source
c	  if both length and width = 0, then it is a point source
c
	  disarea=dislocation(is)
	  if(dx.gt.0.d0)then
	    disarea=disarea*dx
	  endif
	  if(dy.gt.0.d0)then
	    disarea=disarea*dy
	  endif
c
	  do ix=1,nx
	    x=dx*(dble(ix)-0.5d0)
	    do iy=1,ny
	      y=dy*(dble(iy)-0.5d0)
	      nps=nps+1
	      if(nps.gt.NPSMAX)then
	        print *,' Warning: too large number for discrete ',
     &                  'point sources (i.e., NPSMAX too small)!'
	        nwarn=nwarn+1
	        nps=NPSMAX
	        return
	      endif
	      pxs(nps)=xs(is)+x*dcos(st)-y*dcos(di)*dsin(st)
	      pys(nps)=ys(is)+x*dsin(st)+y*dcos(di)*dcos(st)
	      pzs(nps)=zs(is)+y*dsin(di)
	      if(pzs(nps).lt.z1-dz)then
	        print *,' Warning: parts of source rectangles shallower'
	        print *,'          than the Green function grids!'
	        nwarn=nwarn+1
	      endif
	      if(pzs(nps).gt.z2+dz)then
	        print *,' Warning: parts of source rectangles deeper'
	        print *,'          than the Green function grids!'
	        nwarn=nwarn+1
	      endif
c
c	      1 = weight for strike-slip: m12=m21=1;
c	      2 = weight for dip-slip: m13=m31=1
c	      3 = weight for clvd: m33=-m11=-m22=1
c	      4 = weight for 45 deg strike-slip: m11=-m22=1
c	      5 = weight for 45 deg dip-slip: m23=m32=1
c
	      pmoment(1,nps)=sm(1,2)*disarea
	      pmoment(2,nps)=sm(1,3)*disarea
	      pmoment(3,nps)=sm(3,3)*disarea
	      pmoment(4,nps)=0.5d0*(sm(1,1)-sm(2,2))*disarea
	      pmoment(5,nps)=sm(2,3)*disarea
	    enddo
	  enddo
c         write(*,'(a,i2,a,i6,a)')' the ',is,'. rectangle => ',
c    &                            nx*ny,' point sources.'
	enddo
c       write(*,*)'------------------------------------------------'
c       write(*,'(a,i7)')' the total number of point sources: ',nps
c
	return
	end

