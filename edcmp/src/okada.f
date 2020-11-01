	subroutine okada(ns,NSMAX,nrec,NRECMAX,lambda,mu,dislocations,
     &                 xs,ys,zs,lengths,widths,strikes,dips,rakes,
     &                 xrec,yrec,zrec0,disp,strain,tilt)
	implicit none
c
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
        integer ns,NSMAX,nrec,NRECMAX
        double precision lambda,mu
        double precision dislocations(NSMAX)
        double precision xs(NSMAX),ys(NSMAX),zs(NSMAX)
        double precision lengths(NSMAX),widths(NSMAX)
        double precision strikes(NSMAX),dips(NSMAX),rakes(NSMAX)
        double precision xrec(NRECMAX),yrec(NRECMAX)
        double precision zrec0
        double precision disp(NRECMAX,3)
        double precision strain(NRECMAX,6),tilt(NRECMAX,2)
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c       ns = the really used number of rectangular sources
c       NSMAX = the upper limit of ns
c       nrec = the really used number of observation positions
c       NRECMAX = the upper limit of nrec
c
c       lambda, mu = the two Lame constants in Pascal (SI unit)
c
c       (xs,ys,zs) = coordinates of the start point of strike
c       with x = north, y = east, z = downward.
c       all angles in degree.
c       (xrec,yrec,zrec0) = cartesian coordinates of observations
c             Note zrec0 is a fixed constant
c       disp = 3 displacement components: ux,uy,uz
c       strain = 6 strain components: exx,eyy,ezz,exy,eyz,ezx
c       tilt = 2 vertical tilt components: dux/dz, duy/dz
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c	from Okada's subroutine DC3D0:
c
	INTEGER IRET
	REAL*4 ALPHA,X,Y,Z,DEPTH,DIP,POT1,POT2,POT3,POT4,
     &         UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ
c
c	more from Okada's subroutine DC3D:
c
	REAL*4 AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c	LOCAL CONSTANTS
c	===============
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	double precision degtorad,eps
	parameter(degtorad=1.745329252E-02,eps=1.0d-06)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c	LOCAL WORK SPACES
c	=================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	integer j,is,irec
	double precision st,di,ra
	double precision csst,ssst,csra,ssra,csdi,ssdi
	double precision cs2st,ss2st
	double precision disp0(3),tilt0(2),strain0(6)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c	PROCESSING
c	==========
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c	receiver and source independent variables
c
	ALPHA=sngl((lambda+mu)/(lambda+2.d0*mu))
	POT3=0.0
	POT4=0.0
	DISL3=0.0
	AL1=0.0
	AW2=0.0
	Z=-zrec0
c
	do 901 irec=1,nrec
c
c	  initialization
c
	  do j=1,6
            strain(irec,j)=0.d0
	  enddo
	  do j=1,3
            disp(irec,j)=0.d0
	  enddo
	  do j=1,2
            tilt(irec,j)=0.d0
	  enddo
c
	  do 900 is=1,ns
c
	    st=strikes(is)*degtorad
            csst=dcos(st)
            ssst=dsin(st)
            cs2st=dcos(2.d0*st)
            ss2st=dsin(2.d0*st)
c
	    di=dips(is)*degtorad
            csdi=dcos(di)
            ssdi=dsin(di)
c
	    ra=rakes(is)*degtorad
            csra=dcos(ra)
            ssra=dsin(ra)
c
c	    transform from Aki's to Okada's system
c
            X=sngl((xrec(irec)-xs(is))*csst+(yrec(irec)-ys(is))*ssst)
            Y=sngl((xrec(irec)-xs(is))*ssst-(yrec(irec)-ys(is))*csst)
	    DEPTH=sngl(zs(is))
	    DIP=sngl(dips(is))
c
	    if(lengths(is).eq.0.d0.and.widths(is).eq.0.d0)then
c
c	      point source
c
	      POT1=sngl(dislocations(is)*csra)
	      POT2=sngl(dislocations(is)*ssra)
	      IRET=1
	      call DC3D0(ALPHA,X,Y,Z,DEPTH,DIP,POT1,POT2,POT3,POT4,
     *             UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
c	      if(IRET.eq.1)then
c	        stop ' There is a problem in Okada subroutine!'
c	      endif
	    else
c
c	      finite source
c
	      AL2=sngl(lengths(is))
	      AW1=-sngl(widths(is))
	      DISL1=sngl(dislocations(is)*csra)
	      DISL2=sngl(dislocations(is)*ssra)
	      if(lengths(is).eq.0.d0)then
	        AL2=sngl(widths(is)*eps)
	        DISL1=DISL1/AL2
	        DISL2=DISL2/AL2
	      else if(widths(is).eq.0.d0)then
	        AW1=-sngl(lengths(is)*eps)
	        DISL1=DISL1/(-AW1)
	        DISL2=DISL2/(-AW1)
	      endif
	      IRET=1
	      call DC3D(ALPHA,X,Y,Z,DEPTH,DIP,AL1,AL2,AW1,AW2,
     *             DISL1,DISL2,DISL3,UX,UY,UZ,
     *             UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
c	      if(IRET.eq.1)then
c	        stop ' There is a problem in Okada subroutine!'
c	      endif
	    endif
c
c	    transform from Okada's to Aki's system
c
            disp0(1)=dble(UX)*csst+dble(UY)*ssst
            disp0(2)=dble(UX)*ssst-dble(UY)*csst
	    disp0(3)=-dble(UZ)
c
            tilt0(1)=-(dble(UXZ)*csst+dble(UYZ)*ssst)
            tilt0(2)=-(dble(UXZ)*ssst-dble(UYZ)*csst)
c
            strain0(1)=dble(UXX)*csst*csst+dble(UYY)*ssst*ssst
     &                +0.5d0*dble(UXY+UYX)*ss2st
            strain0(2)=dble(UXX)*ssst*ssst+dble(UYY)*csst*csst
     &                -0.5d0*dble(UXY+UYX)*ss2st
            strain0(3)=dble(UZZ)
            strain0(4)=0.5d0*dble(UXX-UYY)*ss2st
     &                -0.5d0*dble(UXY+UYX)*cs2st
            strain0(5)=-0.5d0*dble(UZX+UXZ)*ssst
     &                 +0.5d0*dble(UYZ+UZY)*csst
            strain0(6)=-0.5d0*dble(UZX+UXZ)*csst
     &                 -0.5d0*dble(UYZ+UZY)*ssst
c
            do j=1,3
              disp(irec,j)=   disp(irec,j)   + disp0(j)
            enddo
            do j=1,2
              tilt(irec,j)=   tilt(irec,j)   + tilt0(j)
            enddo
            do j=1,6
              strain(irec,j)= strain(irec,j) + strain0(j)
            enddo
900	  continue
901     continue
c
	return
	end
