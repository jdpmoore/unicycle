	subroutine edcoutput(nrec,infile,outdir,iouts,outputs)
	implicit none
c
c	First implemented in Potsdam, Feb, 1999
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
	include 'edcglobal.h'
c
	integer nrec,iouts(NFIELDS)
	character*80 outdir,infile,outputs(NFIELDS)
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
	integer i,leninf,lendir,irec
	double precision dilatation
	double precision stress(6)
	character*100 title
	character*160 outfile
c
	lendir=index(outdir,' ')-1
c
	if(lendir.lt.1)then
	  stop ' Error in edcmain: wrong for output directory!'
	endif
c
	leninf=index(infile,' ')-1
c       
c	DATA OUTPUT
c	===========
c
	if(iouts(1).eq.1)then
	  print *,'... output 3 displacement components ...'
c
	  outfile=outdir(1:lendir)//outputs(1)
	  open(31,file=outfile,status='unknown')
	  write(31,'(a)')'# Displacements calculated with edcmp'
	  write(31,'(a)')'# The input data file is '//infile(1:leninf)
	  title='#   X_m         Y_m         '
     &        //'Ux_m        Uy_m        Uz_m'
	  write(31,'(a57)')title
	  do irec=1,nrec
	    write(31,1001)xrec(irec),yrec(irec),(disp(irec,i),i=1,3)
	  enddo
	  close(31)
	endif
c
	if(iouts(2).eq.1)then
	  print *,'... output 6 strain components ...'
c
	  outfile=outdir(1:lendir)//outputs(2)
	  open(32,file=outfile,status='unknown')
	  write(32,'(a)')'# Strains calculated with edcmp'
	  write(32,'(a)')'# The input data file is '
     &                 //infile(1:leninf)
	  title='#   X_m         Y_m         Exx         Eyy'
     &         //'         Ezz         Exy         Eyz         Ezx'
	  write(32,'(a91)')title
	  do irec=1,nrec
	    write(32,1001)xrec(irec),yrec(irec),(strain(irec,i),i=1,6)
	  enddo
	  close(32)
	endif
c
	if(iouts(3).eq.1)then
	  print *,'... output 6 stress components ...'
c
	  outfile=outdir(1:lendir)//outputs(3)
	  open(33,file=outfile,status='unknown')
	  write(33,'(a)')'# Stresses calculated with edcmp'
	  write(33,'(a)')'# The input data file is '//infile(1:leninf)
	  title='#   X_m         Y_m         Sxx_Pa      Syy_Pa'
     &	      //'      Szz_Pa      Sxy_Pa      Syz_Pa      Szx_Pa'
	  write(33,'(a94)')title
	  do irec=1,nrec
	    dilatation=strain(irec,1)+strain(irec,2)+strain(irec,3)
	    stress(1)=lambda*dilatation+2.d0*mu*strain(irec,1)
	    stress(2)=lambda*dilatation+2.d0*mu*strain(irec,2)
	    stress(3)=lambda*dilatation+2.d0*mu*strain(irec,3)
	    stress(4)=2.d0*mu*strain(irec,4)
	    stress(5)=2.d0*mu*strain(irec,5)
	    stress(6)=2.d0*mu*strain(irec,6)
	    if(zrec0.eq.0.d0)then
	      stress(3)=0.d0
	      stress(5)=0.d0
	      stress(6)=0.d0
	    endif
	    write(33,1001)xrec(irec),yrec(irec),(stress(i),i=1,6)
	  enddo
	  close(33)
	endif
c
	if(iouts(4).eq.1)then
	  print *,'... output 2 tilt components ...'
c
	  outfile=outdir(1:lendir)//outputs(4)
	  open(34,file=outfile,status='unknown')
	  write(34,'(a)')'# Vertical tilts calculated with edcmp'
	  write(34,'(a)')'# The input data file is '//infile(1:leninf)
	  title='#   X_m         Y_m         dUx/dz      dUy/dz'
	  write(34,'(a46)')title
	  do irec=1,nrec
	    write(34,1001)xrec(irec),yrec(irec),(tilt(irec,i),i=1,2)
	  enddo
	  close(34)
	endif
1001	format(8E12.4)
	return
	end
