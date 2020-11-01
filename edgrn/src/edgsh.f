        subroutine edgsh(y,k,eps)
        implicit none
c
c	First implemented in Potsdam, Feb, 1999
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
c	calculation of response to sh source
c       y(6): solution vector
c       k: wave number
c       eps: relative accuracy
c
        double precision k,eps
        double precision y(6)
c
	include 'edgglobal.h'
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
c       work space
c
        double precision eps0
        parameter(eps0=1.0d-03)
c
        integer i,l,n,lup,llw,key
        double precision fac,exponent
        double precision y0(2),c0(2),b(2)
        double precision y1(2),yup(2),ylw(2)
        double precision hk(2,2),coef(2,2)
c
        do i=5,6
          y(i)=0.d0
        enddo
c
c===============================================================================
c
c       matrix propagation from surface to source
c
        do i=1,2
          c0(i)=0.d0
          y0(i)=0.d0
	  y1(i)=0.d0
          yup(i)=0.d0
        enddo
c
c       yup: the starting solution vector
c
        lup=1
c
	exponent=0.d0
        do l=ls-1,1,-1
          n=nno(l)
	  exponent=exponent-k*hp(l)
          if(exponent.le.dlog(eps))then
            lup=l+1
            if(lup.gt.lzrec)return
            goto 100
          endif
        enddo
100     continue
c
c       determination of starting sublayer for half-space
c
	exponent=0.d0
        llw=lp
        do l=ls,lp-1
          n=nno(l)
	  exponent=exponent-k*hp(l)
          if(exponent.le.dlog(eps))then
            llw=l
            goto 200
          endif
        enddo
200     continue
c
        yup(1)=1.d0
        if(lup.gt.1)then
          n=nno(lup-1)
	  yup(2)=ro(n)*vs(n)**2*k
        endif
        if(lup.eq.lzrec)call memcpy(yup,y0,2)
c
        do l=lup+1,ls
          n=nno(l-1)
c
c         determination of propagation matrix
c
	  call edghask(hk,2,k,hp(l-1),n)
	  call axb(hk,yup,2,2,1,y1)
	  fac=dexp(-k*hp(l-1))
	  do i=1,2
	    yup(i)=y1(i)*fac
	  enddo
c
	  if(l.gt.lzrec)then
	    do i=1,2
	      y0(i)=y0(i)*fac
	    enddo
	  else if(l.eq.lzrec)then
	    call memcpy(yup,y0,2)
	  endif
        enddo
c
c===============================================================================
c
c       matrix propagation from half-space to source
c
        do i=1,2
          c0(i)=0.d0
          y1(i)=0.d0
          ylw(i)=0.d0
        enddo
c
c       ylw: the starting solution vector
c
        n=nno(llw)
        ylw(1)=1.d0
        if(vs(n).lt.eps0*vp(n))then
c
c         the lowest layer is fluid
c
	  ylw(2)=0.d0
        else
c
c         the lowest layer is solid
c
          ylw(2)=-ro(n)*vs(n)**2*k
        endif
        if(llw.gt.ls.and.llw.eq.lzrec)call memcpy(ylw,y0,2)
c
        do l=llw-1,ls,-1
          n=nno(l)
c
c         determination of propagation matrix
c
	  call edghask(hk,2,k,-hp(l),n)
	  call axb(hk,ylw,2,2,1,y1)
	  fac=dexp(-k*hp(l))
	  do i=1,2
	    ylw(i)=y1(i)*fac
	  enddo
c
	  if(l.lt.lzrec)then
	    do i=1,2
	      y0(i)=y0(i)*fac
	    enddo
	  else if(l.gt.ls.and.l.eq.lzrec)then
	    call memcpy(ylw,y0,2)
	  endif
        enddo
c
c===============================================================================
c
c       conditions on the source surface
c
c
c       source function
c
        do i=1,2
          if(kpower(i+4).eq.1)then
            fac=k
          else
            fac=1.d0
          endif
          b(i)=sfct(i+4)*fac
          coef(i,1)=yup(i)
          coef(i,2)=-ylw(i)
        enddo
        key=0
        call gemp(coef,b,2,1,1.d-99,key)
        if(key.eq.0)then
          print *,'warning in edgsh: anormal exit from cgemp!'
          return
        endif
        if(lzrec.le.ls)then
          do i=1,2
            y(i+4)=b(1)*y0(i)
          enddo
        else
          do i=1,2
            y(i+4)=b(2)*y0(i)
          enddo
        endif
        return
        end
