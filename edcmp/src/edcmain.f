	program edcmain
	implicit none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	this program calculates convolution integral (summation)       c
c	of discrete sources to model the elastic deformations induced  c
c       by an eqrthquake.                                              c
c                                                                      c
c	The input data will be read from an input file                 c
c                                                                      c
c	First implemented in Potsdam, Feb, 1999                        c
c	Last modified: Potsdam, Nov, 2001, by R. Wang                  c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	BEGIN DECLARATIONS
c	==================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	GLOBAL CONSTANTS
c	================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
c	OBSERVATION POSITIONS AND OBSERVABLES
c	=====================================
c
c	(xrec(i),yrec(i),zrec0)=coordinates of the observation positions
c	(Note that zrec0 is fixed)
c	disp = the 3 displcement vector components: ux,uy,uz
c	strain = the 6 strain tensor components: exx,eyy,ezz,exy,eyz,ezx
c	tilt = the two vertical tilt components: dux/dz, duy/dz
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
c	ELASTIC PARAMETERS AT OBSERVATION DEPTH
c	=======================================
c
c	lambda,mu = the two Lame constants in pascal
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double precision lambda,mu
c
	common/elasticity/lambda,mu
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
c	MEMORIES FOR OUTPUTS
c	====================
c
c	1 = two displacement components: ux,uy,uz
c	2 = 6 strain components: exx,eyy,ezz,exy,eyz,ezx
c	3 = 6 stress components: sxx,syy,szz,sxy,syz,szx
c	4 = two vertical tilt components: dux/dz,duy/dz
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	character*80 outdir
	integer iouts(NFIELDS)
	character*80 outputs(NFIELDS)
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	LOCAL CONSTANTS
c	==============
c	pi, parameter for transforming degree to radian
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double precision PI,PI2
	parameter(PI=3.14159265,PI2=6.28318531)
	double precision DEGTORAD
	parameter(DEGTORAD=1.745329252E-02)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	LOCAL WORK SPACES
c	=================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	integer i,is,irec,ixyrec,ixrec,iyrec,nxrec,nyrec,imodel
	integer nrec,ns
	double precision xrec1,xrec2,yrec1,yrec2,dxrec,dyrec
	double complex cxyrec1,cxyrec2
	double complex cxyrec(NRECMAX)
	character*80 infile,grndir,grnss,grnds,grncl
	character*180 dataline
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c	END DECLARATIONS
c	================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c00000000000000000000000000000000000000000000000000000000000000000000000
c	BEGIN READ IN INPUT PARAMETERS
c	==============================
c00000000000000000000000000000000000000000000000000000000000000000000000
c
	nwarn=0
c
	print *,'------------------------------------------------------'
	print *,'                    edcmp -'
	print *,'  convolve source and numerical green''s function'
	print *,'------------------------------------------------------'
	write(*,'(a,$)')' Please type the file name of input data: '
	read(*,'(a)')infile
	open(10,file=infile,status='old')
c00000000000000000000000000000000000000000000000000000000000000000000000
c	READ IN PARAMETERS FOR OBSERVATION ARRAY
c	========================================
c00000000000000000000000000000000000000000000000000000000000000000000000
	call getdata(10,dataline)
        read(dataline,*)ixyrec
	if(ixyrec.eq.0)then
c
c	  irregular observation positions
c
	  call getdata(10,dataline)
          read(dataline,*)nrec
	  if(nrec.gt.NRECMAX)then
	    stop ' Error in input file: NRECMAX too small!'
	  endif
	  read(10,*)(cxyrec(irec),irec=1,nrec)
	  do irec=1,nrec
	    xrec(irec)=dreal(cxyrec(irec))
	    yrec(irec)=dimag(cxyrec(irec))
	  enddo
	else if(ixyrec.eq.1)then
c
c	  1D observation profile
c
	  call getdata(10,dataline)
          read(dataline,*)nrec
	  call getdata(10,dataline)
          read(dataline,*)cxyrec1,cxyrec2
	  if(nrec.lt.1)then
	    stop ' Error in input file: wrong input for nrec!'
	  else
	    xrec(1)=dreal(cxyrec1)
	    yrec(1)=dimag(cxyrec1)
	    if(nrec.gt.1)then
	      dxrec=dreal(cxyrec2-cxyrec1)/dble(nrec-1)
	      dyrec=dimag(cxyrec2-cxyrec1)/dble(nrec-1)
	    else
	      dxrec=0.d0
	      dyrec=0.d0
	    endif
	    do irec=1,nrec
	      xrec(irec)=dreal(cxyrec1)+dxrec*dble(irec-1)
	      yrec(irec)=dimag(cxyrec1)+dyrec*dble(irec-1)
	    enddo
	  endif
	else if(ixyrec.eq.2)then
c
c	  2D rectanglar observation array
c
	  call getdata(10,dataline)
          read(dataline,*)nxrec,xrec1,xrec2
	  call getdata(10,dataline)
          read(dataline,*)nyrec,yrec1,yrec2
	  nrec=nxrec*nyrec
	  if(nrec.gt.NRECMAX.or.nrec.lt.1)then
	    stop ' Error in input file: wrong input for nrec!'
	  endif
	  irec=0
	  if(nxrec.gt.1)then
	    dxrec=(xrec2-xrec1)/dble(nxrec-1)
	  else
	    dxrec=0.d0
	  endif
	  if(nyrec.gt.1)then
	    dyrec=(yrec2-yrec1)/dble(nyrec-1)
	  else
	    dyrec=0.d0
	  endif
	  do iyrec=1,nyrec
	    do ixrec=1,nxrec
	      irec=irec+1
	      xrec(irec)=xrec1+dxrec*dble(ixrec-1)
	      yrec(irec)=yrec1+dyrec*dble(iyrec-1)
	    enddo
	  enddo
	else
	  stop' Error in input file: wrong input for ixyrec!'
	endif
c00000000000000000000000000000000000000000000000000000000000000000000000
c	READ IN OUTPUT PARAMETERS
c	=========================
c00000000000000000000000000000000000000000000000000000000000000000000000
	call getdata(10,dataline)
        read(dataline,*)outdir
	call getdata(10,dataline)
        read(dataline,*)(iouts(i),i=1,NFIELDS)
	call getdata(10,dataline)
        read(dataline,*)(outputs(i),i=1,NFIELDS)
c00000000000000000000000000000000000000000000000000000000000000000000000
c	READ IN PARAMETERS FOR RECTANGULAR SOURCES
c	==========================================
c00000000000000000000000000000000000000000000000000000000000000000000000
	call getdata(10,dataline)
        read(dataline,*)ns
	if(ns.gt.NSMAX)then
	  stop ' Error in edcmain: to large no of source rectangles!'
	endif
	do is=1,ns
	  call getdata(10,dataline)
          read(dataline,*)i,dislocation(is),xs(is),ys(is),zs(is),
     &                    length(is),width(is),
     &                    strike(is),dip(is),rake(is)
	  if(length(is).lt.0.d0.or.width(is).lt.0.d0)then
	    stop ' Error in input data: source length or width < 0!'
	  endif
	  if(zs(is).lt.0.d0)then
	    stop ' Error in input data: source depth zs < 0!'
	  endif
	  if(length(is).eq.0.d0.and.width(is).eq.0.d0)then
	    write(*,'(a,i2,a)')' the ',is,'. rectangle is a point.'
	  else if(length(is).gt.0.d0.and.width(is).eq.0.d0)then
	    write(*,'(a,i2,a)')' the ',is,
     &                         '. rectangle is a horizontal line.'
	  else if(length(is).eq.0.d0.and.width(is).gt.0.d0)then
	    write(*,'(a,i2,a)')' the ',is,
     &                         '. rectangle is a vertical line.'
	  endif
	enddo
c00000000000000000000000000000000000000000000000000000000000000000000000
c	READ IN PARAMETERS FOR EARTH MODEL CHOICE
c	=========================================
c00000000000000000000000000000000000000000000000000000000000000000000000
	call getdata(10,dataline)
        read(dataline,*)imodel
	if(imodel.eq.0)then
	  call getdata(10,dataline)
          read(dataline,*)zrec0,lambda,mu
	else if(imodel.eq.1)then
	  call getdata(10,dataline)
          read(dataline,*)grndir,grnss,grnds,grncl
	else
	  stop ' Error in input file: wrong choice of earth model!'
	endif
	close(10)
c00000000000000000000000000000000000000000000000000000000000000000000000
c	END READ IN INPUT PARAMETERS
c	============================

	print *,'... input data successful ...'

c	BEGIN PROCESSING
c	================

	if(imodel.eq.1)then
	  print *,'... layered half-space model considered ...'
	  print *,'... use the Green function approach ...'
	  call edcgrn(ns,nrec,grndir,grnss,grnds,grncl)
	else
	  print *,'... homogeneous half-space model considered ...'
	  print *,'... use analytical solutions of Okada ...'
	  call okada(ns,NSMAX,nrec,NRECMAX,lambda,mu,dislocation,
     &               xs,ys,zs,length,width,strike,dip,rake,
     &               xrec,yrec,zrec0,disp,strain,tilt)
	endif
	print *,'... outputs ...'
	call edcoutput(nrec,infile,outdir,iouts,outputs)

c	END OF STANDARD PROCESSING
c	==========================

	if(nwarn.eq.0)then
	  print *,'----------------------------------------------------'
	  print *,'             successful computation        '
	  print *,'----------------------------------------------------'
	else
	  print *,'----------------------------------------------------'
	  print *,'     Sorry, there have been',nwarn,' warnings.      '
	  print *,'             Results may be inaccurate!             '
	  print *,'----------------------------------------------------'
          stop 1
	endif
c
	stop
	end
