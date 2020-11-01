        subroutine edgwvint(u,r,nr,srate,lambda,mu,itty)
        implicit none
c
c	First implemented in Potsdam, Feb, 1999
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
        integer nr,itty
        double precision srate,lambda,mu
c
	include 'edgglobal.h'
c
        double precision r(nrmax)
        double precision u(10,nrmax)
c
c 	table of J_n(x), dJ_n(x)/dx and n*J_n(x)/x
c	all multiplied by sqrt(x)
c
	integer nnbess,nnbess1
	parameter(nnbess=nbess*ndbess,nnbess1=nnbess+ndbess)
	double precision bsdx,bsfct(0:nnbess1,3)
	common /bessels/ bsdx,bsfct
c
        double precision eps,eps0
        parameter(eps=1.0d-08,eps0=1.0d-03)
	double precision pi,pi2
	parameter(pi=3.14159265358979d0,pi2=6.28318530717959d0)
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
c       n0: number of model layers
c
        integer n0
        double precision h(lmax),ro(lmax),vp(lmax),vs(lmax)
        common /model/ h,ro,vp,vs,n0
c
c       source parameters
c
        integer ls,ms,ics
        integer kpower(6)
        double precision zs,r0
        double precision sfct(6)
        common /source/ zs,r0,sfct,ls,ms,ics,kpower
c
c       parameters for hankel integrations
c
        integer i,ir,ir1,ncall,ik,nk,nx
        double precision k,k0,klimit,dk,x,wl,wr
        double precision cs,fps,fsh,fac,yabs,dyabs,ymax
        double precision y(6),y0(6),u0(6),uk0(6),bs(3),bs1(3)
	double precision r00(nrmax)
        logical ps,sh,analytic
c
c	ics = 1  when the azmuth-factor is cos(ms*theta) for poloidal mode
c	         (psv) and sin(ms*theta) for the toroidal mode (sh);
c	ics = -1 otherwise.
c
        cs=dble(ics)
	fps=0.d0
	do i=1,4
	  fps=fps+dabs(sfct(i))
	enddo
	if(fps.gt.0.d0)then
	  ps=.true.
	else
	  ps=.false.
	endif
	fsh=0.d0
	do i=5,6
	  fsh=fsh+dabs(sfct(i))
	enddo
	if(fsh.gt.0.d0)then
	  sh=.true.
	else
	  sh=.false.
	endif
c
c	u: 1=uz, 2=ur, 3=ut, 4=ezz, 5=err, 6=ett, 7=ezr, 8=ert, 9=etz
c	  10=duz/dr
c	NOTE: uz, ur, ezz, err, ett, ezr duz/dr have the same azimuth-factor
c	      as the poloidal mode (p-sv);
c	      ut, ert and etz have the same azimuth-factor as the
c	      toroidal mode (sh);
c
        do ir=1,nr
	  r00(ir)=dmax1(r0,1.0d-02*dabs(zs-zrec),1.0d-02*r(ir))
          do i=1,10
            u(i,ir)=0.d0
          enddo
        enddo
	do i=1,6
	  u0(i)=0.d0
	  uk0(i)=0.d0
	enddo
	do i=1,6
	  y0(i)=0.d0
	enddo
c
c	determine wavenumber limit
c
c	determine limits of y(i), i=1,...,6
c
        ncall=0
	if(zs.eq.zrec)then
	  k0=eps0*pi2/(r0+dabs(zs-zrec)+r(nr))
          ymax=0.d0
10        yabs=0.d0
          call edgkern(y,k0,ps,sh,eps)
          ncall=ncall+1
          do i=1,5,2
            yabs=yabs+y(i)*y(i)
          enddo
	  yabs=k0*dsqrt(k0*yabs)*dexp(-(k0*r0)**2)
	  ymax=dmax1(ymax,yabs)
          if(yabs.gt.eps0*ymax)then
            k0=1.25d0*k0
            goto 10
          endif
c
	  analytic=.true.
c
	  k=eps0*pi2/(r0+dabs(zs-zrec)+r(nr))
20	  call edgkern(y,k,ps,sh,eps)
          ncall=ncall+1
	  do i=2,6,2
	    y(i)=y(i)/(k*mu)
	  enddo
	  yabs=0.d0
	  dyabs=0.d0
	  do i=1,6
	    yabs=yabs+y(i)**2
	    dyabs=dyabs+(y(i)-y0(i))**2
	    y0(i)=y(i)
	  enddo
	  if(dyabs.gt.eps*yabs)then
	    if(k.ge.k0)then
	      analytic=.false.
	      do i=1,6
	        y0(i)=0.d0
	      enddo
	    else
	      k=1.25d0*k
	      goto 20
	    endif
	  endif
	  do i=2,6,2
	    y0(i)=y0(i)*mu
	  enddo
	else
	  analytic=.false.
	endif
	  
c
	klimit=eps*pi2/(r0+dabs(zs-zrec)+r(nr))
        ymax=0.d0
30      yabs=0.d0
        call edgkern(y,klimit,ps,sh,eps)
        ncall=ncall+1
        do i=1,5,2
          yabs=yabs+(y(i)-y0(i))**2
        enddo
	yabs=klimit*dsqrt(klimit*yabs)*dexp(-(klimit*r00(1))**2)
	ymax=dmax1(ymax,yabs)
        if(yabs.gt.eps0*ymax)then
          klimit=1.2d0*klimit
          goto 30
        endif
c
c	determine wavenumber sampling rate
c
	dk=pi2/(srate*(r00(1)+r(nr))+dabs(zs-zrec))
	nk=500+idnint(klimit/dk)
	dk=klimit/dble(nk)
c
c	too small distances will be treated as r = 0!
c
	if(r(1).gt.0.d0)then
	  ir1=1
	else
	  ir1=2
	endif
c
        do ik=1,nk
          k=dble(ik)*dk
          call edgkern(y,k,ps,sh,eps)
	  if(analytic)then
	    do i=1,5,2
	      y(i)=y(i)-y0(i)
	    enddo
	    do i=2,6,2
	      y(i)=y(i)-y0(i)*k
	    enddo
	  else if(ir1.eq.2)then
c
c	  for r=0
c
	    fac=k*dexp(-(k*r00(1))**2)*dk
	    do i=1,6
	      u0(i)=u0(i)+y(i)*fac
	      uk0(i)=uk0(i)+y(i)*k*fac
	    enddo
	  endif
          do ir=ir1,nr
	    fac=dsqrt(k)*dexp(-(k*r00(ir))**2)
	    if(k*fac.gt.eps)then
	      fac=fac*dk/dsqrt(r(ir))
	      x=k*r(ir)
c
c	      bessels functions from pre-calculated tables
c
	      nx=idint(x/bsdx)
	      wr=x/bsdx-dble(nx)
	      wl=1.d0-wr
	      if(nx.gt.nnbess)then
	        nx=nnbess+mod(nx-nnbess,ndbess)
	        do i=1,3
	          bs(i)=fac*(wl*bsfct(nx,i)+wr*bsfct(nx+1,i))
	        enddo
	        bs(3)=bs(3)*(dble(nx)+wr)*bsdx/x
	      else
	        do i=1,3
	          bs(i)=fac*(wl*bsfct(nx,i)+wr*bsfct(nx+1,i))
	        enddo
	      endif
c
c	      u1-3 are displacement components:
c	      u4 = normal stress: szz
c	      u5 = surface strain: err+ett
c	      u6 will be derived later
c	      u7 = shear stress: szr
c	      u8 = strain component: dut/dr - (dur/dt)/r + ut/r
c	      u9 = shear stress: szt
c	      u10 = tilt: duz/dr
c
	      u(1,ir)=u(1,ir)+y(1)*bs(1)
	      u(2,ir)=u(2,ir)+y(3)*bs(2)+cs*y(5)*bs(3)
	      u(3,ir)=u(3,ir)-cs*y(3)*bs(3)-y(5)*bs(2)
	      u(4,ir)=u(4,ir)+y(2)*bs(1)
	      u(5,ir)=u(5,ir)-y(3)*k*bs(1)
	      u(7,ir)=u(7,ir)+y(4)*bs(2)+cs*y(6)*bs(3)
	      u(8,ir)=u(8,ir)+y(5)*k*bs(1)
	      u(9,ir)=u(9,ir)-cs*y(4)*bs(3)-y(6)*bs(2)
	      u(10,ir)=u(10,ir)+y(1)*k*bs(2)
	    endif
	  enddo
        enddo
c
c       end of total integral
c
        if(itty.eq.1)then
          write(*,'(a,i7,a,i7)')'   wavenumber samples: ',ncall+nk,
     &                          ', really used: ',nk
        endif
c
	if(ir1.eq.2.and..not.analytic)then
c
c	  for very small r including r=0
c
	  if(ms.eq.0)then
	    u(1,1)=u0(1)
	    u(4,1)=u0(2)
	    u(5,1)=-0.5d0*uk0(3)
	    u(6,1)=u(5,1)
	  else if(ms.eq.1)then
	    u(2,1)=0.5d0*(u0(3)+cs*u0(5))
	    u(3,1)=-0.5d0*(cs*u0(3)+u0(5))
	    u(7,1)=0.5d0*(u0(4)+cs*u0(6))
	    u(9,1)=-0.5d0*(cs*u0(4)+u0(6))
	    u(10,1)=0.5d0*uk0(1)
	  else if(ms.eq.2)then
	    u(5,1)=0.25d0*(uk0(3)+cs*uk0(5))
	    u(6,1)=-u(5,1)
	    u(8,1)=-0.25d0*(cs*uk0(3)+uk0(5))
	  endif
	endif
	do ir=ir1,nr
	  if(analytic)then
	    if(ms.eq.0)then
	      bs(1)=0.d0
	      bs(2)=-1.d0/r(ir)**2
	      bs(3)=0.d0
	      bs1(1)=-1.d0/r(ir)**3
	      bs1(2)=0.d0
	      bs1(3)=0.d0
	    else if(ms.eq.1)then
	      bs(1)=1.d0/r(ir)**2
	      bs(2)=-1.d0/r(ir)**2
	      bs(3)=1.d0/r(ir)**2
	      bs1(1)=0.d0
	      bs1(2)=-2.d0/r(ir)**3
	      bs1(3)=1.d0/r(ir)**3
	    else if(ms.eq.2)then
	      bs(1)=2.d0/r(ir)**2
	      bs(2)=-1.d0/r(ir)**2
	      bs(3)=2.d0/r(ir)**2
	      bs1(1)=3.d0/r(ir)**3
	      bs1(2)=-4.d0/r(ir)**3
	      bs1(3)=4.d0/r(ir)**3
	    endif
	    u(1,ir)=u(1,ir)+y0(1)*bs(1)
	    u(2,ir)=u(2,ir)+y0(3)*bs(2)+cs*y0(5)*bs(3)
	    u(3,ir)=u(3,ir)-y0(5)*bs(2)-cs*y0(3)*bs(3)
	    u(4,ir)=u(4,ir)+y0(2)*bs1(1)
	    u(5,ir)=u(5,ir)-y0(3)*bs1(1)
	    u(7,ir)=u(7,ir)+y0(4)*bs1(2)+cs*y0(6)*bs1(3)
	    u(8,ir)=u(8,ir)+y0(5)*bs1(1)
	    u(9,ir)=u(9,ir)-y0(6)*bs1(2)-cs*y0(4)*bs1(3)
	    u(10,ir)=u(10,ir)+y0(1)*bs1(2)
	  endif
c
c	  u6 is ett = ur/r + (dut/dt)/r
c
	  u(6,ir)=(u(2,ir)+cs*dble(ms)*u(3,ir))/r(ir)
c
c	  u5 now is err = u5(before) - ett
c
	  u(5,ir)=u(5,ir)-u(6,ir)
c
c	  u8 now is ert = 0.5 * u8(before) + (dur/dt)/r - ut/r
c	                = 0.5 * (dut/dr + (dur/dt)/r - ut/r)
c
	  u(8,ir)=0.5d0*u(8,ir)-(cs*dble(ms)*u(2,ir)+u(3,ir))/r(ir)
	enddo
	do ir=1,nr
	  u(4,ir)=(u(4,ir)-lambda*(u(5,ir)+u(6,ir)))/(lambda+2.d0*mu)
	  u(7,ir)=u(7,ir)/(2.d0*mu)
	  u(9,ir)=u(9,ir)/(2.d0*mu)
	enddo
c
        return
        end
