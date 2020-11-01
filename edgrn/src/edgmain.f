        program edgmain
        implicit none
c
c	First implemented in Potsdam, Feb, 1999
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
	include 'edgglobal.h'
c
	double precision pi,pi2
	parameter(pi=3.14159265358979d0,pi2=6.28318530717959d0)
	double precision zseps
	parameter(zseps=1.0d-02)
c
        integer nr
        double precision r(nrmax)
        common /rprofile/ r,nr
c
        integer lp,nno(nzmax)
        double precision hp(nzmax)
        common /sublayer/ hp,lp,nno
c
c       zrec: receiver depth
c       lzrec: sublayer no of receiver
c
        integer lzrec
        double precision zrec
        common /receiver/ zrec,lzrec
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
c       source parameters
c
        integer ls,ms,ics
        integer kpower(6)
        double precision zs,r0
        double precision sfct(6)
        common /source/ zs,r0,sfct,ls,ms,ics,kpower
c
c 	table of J_n(x), dJ_n(x)/dx and n*J_n(x)/x
c	all multiplied by sqrt(x)
c
	integer nnbess,nnbess1
	parameter(nnbess=nbess*ndbess,nnbess1=nnbess+ndbess)
	double precision bsdx,bsfct(0:nnbess1,3)
	common /bessels/ bsdx,bsfct
c
c       work space
c
        integer i,j,l,itty,ierr,izs,nzs,istype,len,unit
        double precision r1,r2,dr,zs1,zs2,dzs,lambda,mu
	double precision srate,moment
	double precision resolut(3)
        double precision u(10,nrmax)
        character*80 inputfile,outdir,grnfile0(3),reports(3)
	character*160 grnfile(3)
	character*180 dataline
c
c       read input file file
c
	print *,'######################################################'
	print *,'#                                                    #'
	print *,'#                  Welcome to                        #'
	print *,'#                                                    #'
	print *,'#                                                    #'
	print *,'#     EEEEE   DDDD     GGGG   RRRR    N   N          #'
	print *,'#     E       D   D   G       R   R   NN  N          #'
	print *,'#     EEEE    D   D   G GGG   RRRR    N N N          #'
	print *,'#     E       D   D   G   G   R R     N  NN          #'
	print *,'#     EEEEE   DDDD     GGGG   R  R    N   N          #'
	print *,'#                                                    #'
	print *,'#                                                    #'
	print *,'#                      by                            #'
	print *,'#                 Rongjiang Wang                     #'
	print *,'#              (wang@gfz-potsdam.de)                 #'
	print *,'#                                                    #'
	print *,'#           GeoForschungsZentrum Potsdam             #'
	print *,'#                   Nov. 2001                        #'
	print *,'######################################################'
	print *,'                                                      '
	write(*,'(a,$)')' Please type the file name of input data: '
        read(*,'(a)')inputfile
        open(10,file=inputfile,status='old')
c
c	source-observation grid parameters
c	==================================
c
	call getdata(10,dataline)
        read(dataline,*)zrec
        call getdata(10,dataline)
	read(dataline,*)nr,r1,r2
	if(nr.le.1)then
	  stop ' Error in edgmain: number of distance samples < 2!'
	else if(nr.gt.nrmax)then
	  stop ' Error in edgmain: no of distance samples > nrmax!'
	else
	  dr=(r2-r1)/dble(nr-1)
	  if(dr.le.0.d0)then
	    stop ' Error in edgmain: wrong distance stepping!'
	  endif
          do i=1,nr
	    r(i)=r1+dr*dble(i-1)
	  enddo
	endif
        call getdata(10,dataline)
	read(dataline,*)nzs,zs1,zs2
	if(nzs.le.1)then
	  stop ' Error in edgmain: number of source depths < 2!'
	else if(nzs.gt.nzsmax)then
	  stop ' Error in edgmain: no of source depths > nzsmax!'
	else
	  dzs=(zs2-zs1)/dble(nzs-1)
	  if(dzs.le.0.d0)then
	    stop ' Error in edgmain: wrong source depth stepping!'
	  endif
	endif
	r0=0.5d0*dmax1(dzs,dr)
	do izs=1,nzs
	  zs=zs1+dble(izs-1)*dzs
	  if(dabs(zrec-zs).gt.0.d0
     &      .and.dabs(zrec-zs).le.zseps*dzs)then
	    zrec=zs
	    print *,'warning: an insignificant depth difference ',
     &              'between source and observation ignored!'
	  endif
	enddo
c
c	wavenumber integration parameters
c	=================================
c
	call getdata(10,dataline)
        read(dataline,*)srate
	if(srate.le.0.d0)srate=10.d0
	itty=1
c
c	output files
c	============
c
	call getdata(10,dataline)
        read(dataline,*)outdir,(grnfile0(i),i=1,3)
c
	len=index(outdir,' ')-1
c
	if(len.lt.1)then
	  stop ' Error in edkmain: wrong for output directory!'
	endif
c
	do i=1,3
	  grnfile(i)=outdir(1:len)//grnfile0(i)
	enddo
c
c	global model parameters
c	=======================
c
	do i=1,3
	  resolut(i)=5.d-02
	enddo
c
c	multilayered model parameters
c	=============================
c
	call getdata(10,dataline)
	read(dataline,*)l
	do i=1,l
	  call getdata(10,dataline)
	  read(dataline,*)j,h(i),vp(i),vs(i),ro(i)
	  if(vs(i).le.0.d0)vs(i)=1.0d-02*vp(i)
	enddo
c
c	end of inputs
c	=============
c
	close(10)
c
c	determine upper und lower parameter values of each layer
c
	l0=1
	do i=2,l
	  if(h(i).gt.h(i-1))then
	    z1(l0)=h(i-1)
	    vp1(l0)=vp(i-1)
	    vs1(l0)=vs(i-1)
	    ro1(l0)=ro(i-1)
c
	    z2(l0)=h(i)
	    vp2(l0)=vp(i)
	    vs2(l0)=vs(i)
	    ro2(l0)=ro(i)
	    l0=l0+1
	  else
	    z1(l0)=h(i)
	    vp1(l0)=vp(i)
	    vs1(l0)=vs(i)
	    ro1(l0)=ro(i)
	  endif
	enddo
	if(z1(l0).eq.0.d0)then
	    z1(l0)=h(l)
	    vp1(l0)=vp(l)
	    vs1(l0)=vs(l)
	    ro1(l0)=ro(l)
	endif
c
c       construction of sublayers at the cutoff frequency
c
	write(*,*)'   '
	if(itty.eq.1)then
	  write(*,*)'the multi-layered model used:'
	endif
	call edgsublay(resolut,itty,ierr)
	if(ierr.eq.1)then
	  stop 'the layer index (lmax) too small defined!'
	endif
c
        call edglayer(zs1,ls,zrec,lzrec,h,n0)
	mu=ro(nno(lzrec))*vs(nno(lzrec))*vs(nno(lzrec))
	lambda=ro(nno(lzrec))*vp(nno(lzrec))*vp(nno(lzrec))-2.d0*mu
	reports(1)='Source type: strike-slip'
	reports(2)='Source type: dip-slip'
	reports(3)='Source type: compensated linear vector dipole'
	len=index(inputfile,' ')-1
c
	do istype=1,3
c
c         open output file
c
	  unit=20+istype
	  write(*,'(a45)')reports(istype)
	  open(unit,file=grnfile(istype),status='unknown')
	  write(unit,'(a)')'# Green functions calculated with '
     &                     //'the program edgrn'
	  write(unit,'(a)')'# The input data file is '
     &                    //inputfile(1:len)
	  write(unit,'(a,a45,a)')'# ',reports(istype),
     &                    '   Dislocation: 1 meter'
	  write(unit,'(a)')'# The 1. data line includes the parameters:'
	  write(unit,'(a)')'#'
	  write(unit,'(a,$)')'# nr, r1[m], r2[m]; '
	  write(unit,'(a,$)')'nzs, zs1[m], zs2[m]; obs.depth[m];'
	  write(unit,'(a)')' lambda[Pa], mu[Pa] (at obs. depth)'
	  write(unit,'(a)')'#'
	  write(unit,'(a)')'# The following data lines are'
	  write(unit,'(a)')'#'
	  if(istype.le.2)then
	    write(unit,'(a)')'#   uz[m]         ur[m]         ut[m] '
     &        //'        ezz           err           ett   '
     &        //'        ezr           ert           etz   '
     &        //'        duz/dr'
	  else
	    write(unit,'(a)')'#   uz[cm]        ur[cm]'
     &        //'        ezz           err           ett   '
     &        //'        ezr           duz/dr'
	  endif
	  write(unit,'(a)')'#'
	  write(unit,1001)nr,r1,r2,nzs,zs1,zs2,zrec,lambda,mu
c
c	  calculate Bessel function tables
c
	  call edgbstab(3-istype)
c
	  do izs=1,nzs
	    zs=zs1+dble(izs-1)*dzs
            call edglayer(zs,ls,zrec,lzrec,h,n0)
	    moment=ro(nno(ls))*vs(nno(ls))*vs(nno(ls))
	    call edgmoment(istype,moment,ro(nno(ls)),
     &                     vp(nno(ls)),vs(nno(ls)))
c
            call edgwvint(u,r,nr,srate,lambda,mu,itty)
	    if(istype.le.2)then
	      do j=1,nr
	        write(unit,'(10E14.6)')(u(i,j),i=1,10)
	      enddo
	    else
	      do j=1,nr
	        write(unit,'(10E14.6)')(u(i,j),i=1,2),(u(i,j),i=4,7)
     &                               ,u(10,j)
	      enddo
	    endif
	  enddo
	  close(unit)
	enddo
	print *,'######################################################'
	print *,'#                                                    #'
	print *,'#        Successful computations with EDGRN          #'
	print *,'#                                                    #'
	print *,'######################################################'
1001	format(2(i4,2E13.5),4E13.5)
 500    stop
        end
