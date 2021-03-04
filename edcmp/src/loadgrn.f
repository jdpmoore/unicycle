        subroutine loadgrn(grndir)
        implicit none

        character*80 grndir

        include 'edcglobal.h'

        double precision lambda,mu
        common/elasticity/lambda,mu

        double precision DEGTORAD,ZSEPS
        parameter(DEGTORAD=1.745329252E-02,ZSEPS=1.0d-02)

        integer nr,nz
        double precision r1,r2,z1,z2
        double precision ssdisp(3,NRMAX+2,NZMAX+2)
        double precision ssstrn(6,NRMAX+2,NZMAX+2)
        double precision ssuzr(NRMAX+2,NZMAX+2)
        double precision dsdisp(3,NRMAX+2,NZMAX+2)
        double precision dsstrn(6,NRMAX+2,NZMAX+2)
        double precision dsuzr(NRMAX+2,NZMAX+2)
        double precision cldisp(2,NRMAX+2,NZMAX+2)
        double precision clstrn(4,NRMAX+2,NZMAX+2)
        double precision cluzr(NRMAX+2,NZMAX+2)
        common/grnfcts/nr,nz,r1,r2,z1,z2,ssdisp,ssstrn,ssuzr,
     &                 dsdisp,dsstrn,dsuzr,cldisp,clstrn,cluzr

        integer i,irec,jrec,ips,nps,ir,iz,lendir,nx,ny,idis
        double precision x1,x2,y1,y2,zobs,si,co,si2,co2,dis,ddis,azi
        double precision dr,dz,dzs,w00,w10,w01,w11
        double precision la1,mu1
        character*160 grnss,grnds,grncl
        character*180 dataline

        double precision xrec(NRECMAX),yrec(NRECMAX)
        double precision zrec0
        double precision disp(NRECMAX,3),strain(NRECMAX,6)
        double precision tilt(NRECMAX,2)
        common/obsarray/xrec,yrec,zrec0,disp,strain,tilt


c	READ IN GREEN'S FUNCTIONS
c	=========================

        lendir=index(grndir,' ')-1
        grnss=grndir(1:lendir)//'.ss'
        grnds=grndir(1:lendir)//'.ds'
        grncl=grndir(1:lendir)//'.cl'

c	READ GREEN'S FUNCTIONS
c	======================

c	for point strike-slip source
c	ssdisp(1-3): Uz, Ur, Ut
c	ssstrn(1-6): Ezz,Err,Ett,Ezr=Erz,Ert=Etr,Etz=Ezt

        print *,'... read Green fcts for the strike-slip source ...'
        print *, 'read ', grnss
        open(21,file=grnss,status='old')
        write (42,*) "open 21 successful"
        call getdata(21,dataline)
        read(dataline,*)nr,r1,r2,nz,z1,z2,zrec0,lambda,mu
        if (nr.gt.NRMAX) then
          stop ' Error in loadgrn: too large no of Green fct distances!'
        else if(nz.gt.NZMAX)then
          stop ' Error in loadgrn: too large no of Green fct depths!'
        endif
        read(21,*)(((ssdisp(i,ir,iz),i=1,3),(ssstrn(i,ir,iz),i=1,6),
     &            ssuzr(ir,iz),ir=1,nr),iz=1,nz)
        close(21)
c
c	for point dip-slip source
c
c	dsdisp(1-3): Uz, Ur, Ut
c	dsstrn(1-6): Ezz,Err,Ett,Ezr=Erz,Ert=Etr,Etz=Ezt
c
        print *,'... read Green fcts for the dip-slip source ...'
	open(22,file=grnds,status='old')
	call getdata(22,dataline)
        read(dataline,*)nx,x1,x2,ny,y1,y2,zobs,la1,mu1
	if(nx.ne.nr.or.x1.ne.r1.or.x2.ne.r2.or.
     &     ny.ne.nz.or.y1.ne.z1.or.y2.ne.z2)then
	  stop ' Error in loadgrn: different grids in Green functions!'
	endif
	if(la1.ne.lambda.or.mu1.ne.mu)then
	  stop ' Error in loadgrn: different Lame const. in Green fcts!'
	endif
	if(zobs.ne.zrec0)then
	  stop ' Error in loadgrn: different obs.depth in Green fcts!'
	endif
	read(22,*)(((dsdisp(i,ir,iz),i=1,3),(dsstrn(i,ir,iz),i=1,6),
     &            dsuzr(ir,iz),ir=1,nr),iz=1,nz)
	close(22)
c
c	for point clvd source
c
c
c	cldisp(1-2): Uz, Ur (Ut=0)
c	clstrn(1-4): Ezz,Err,Ett,Ezr=Erz (Ert=Etr=Etz=Ezt=0)
c
        print *,'... read Green fcts for the clvd source ...'
	open(23,file=grncl,status='old')
	call getdata(23,dataline)
        read(dataline,*)nx,x1,x2,ny,y1,y2,zobs,la1,mu1
	if(nx.ne.nr.or.x1.ne.r1.or.x2.ne.r2.or.
     &     ny.ne.nz.or.y1.ne.z1.or.y2.ne.z2)then
	  stop ' Error in loadgrn: different grids in Green functions!'
	endif
	if(la1.ne.lambda.or.mu1.ne.mu)then
	  stop ' Error in loadgrn: different Lame const. in Green fcts!'
	endif
	if(zobs.ne.zrec0)then
	  stop ' Error in loadgrn: different obs.depth in Green fcts!'
	endif
	read(23,*)(((cldisp(i,ir,iz),i=1,2),(clstrn(i,ir,iz),i=1,4),
     &            cluzr(ir,iz),ir=1,nr),iz=1,nz)
	close(23)
c
	do iz=nz+1,nz+2
	  do ir=nr+1,nr+2
	    do i=1,3
	      ssdisp(i,ir,iz)=0.d0
	      dsdisp(i,ir,iz)=0.d0
	    enddo
	    do i=1,2
	      cldisp(i,ir,iz)=0.d0
	    enddo
	    do i=1,6
	      ssstrn(i,ir,iz)=0.d0
	      dsstrn(i,ir,iz)=0.d0
	    enddo
	    do i=1,4
	      clstrn(i,ir,iz)=0.d0
	    enddo
	    ssuzr(ir,iz)=0.d0
	    dsuzr(ir,iz)=0.d0
	    cluzr(ir,iz)=0.d0
	  enddo
	enddo
c
        print *,'... read in Green fcts successful ...'
c
c	DISCRETISATION OF RECTANGULAR PLANE SOURCES
c	===========================================
c
        print *,'... discretise the finite sources ...'
c
	dr=(r2-r1)/dble(nr-1)
	dz=(z2-z1)/dble(nz-1)

	if(zrec0.ge.z1.and.zrec0.le.z2.and.r1.eq.0.d0)then
	  iz=idnint((zrec0-z1)/dz)+1
	  dzs=(zrec0-(z1+dz*dble(iz-1)))/dz
	  if(dabs(dzs).le.ZSEPS)then
	    do i=1,3
	      ssdisp(i,1,iz)=ssdisp(i,2,iz)
	      dsdisp(i,1,iz)=dsdisp(i,2,iz)
	    enddo
	    do i=1,2
	      cldisp(i,1,iz)=cldisp(i,2,iz)
	    enddo
	    do i=1,6
	      ssstrn(i,1,iz)=ssstrn(i,2,iz)
	      dsstrn(i,1,iz)=dsstrn(i,2,iz)
	    enddo
	    do i=1,4
	      clstrn(i,1,iz)=clstrn(i,2,iz)
	    enddo
	    ssuzr(1,iz)=ssuzr(2,iz)
	    dsuzr(1,iz)=dsuzr(2,iz)
	    cluzr(1,iz)=cluzr(2,iz)
	  endif
	endif
c
        end subroutine loadgrn
