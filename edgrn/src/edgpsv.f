        subroutine edgpsv(y,k,eps)
        implicit none
c
c	First implemented in Potsdam, Feb, 1999
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
c	calculation of response to p-sv source
c       y(6): solution vector
c       k: wave number
c       eps: relative accuracy
c
	include 'edgglobal.h'
c
        double precision k,eps
        double precision y(6)
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
        integer i,j,l,n,lup,llw,key
        double precision h0,fac,exponent,wave,dro
        double precision norm(2)
        double precision y0(4,2),c0(4,2),c1(4,2),b(4)
        double precision y1(4,2),yup(4,2),ylw(4,2)
        double precision ma(4,4),mai(4,4),hk(4,4)
        double precision orth(2,2),coef(4,4)
c
        do i=1,4
          y(i)=0.d0
        enddo
c
c===============================================================================
c
c       matrix propagation from surface to source
c
        do j=1,2
          do i=1,4
            c0(i,j)=0.d0
            yup(i,j)=0.d0
            y0(i,j)=0.d0
          enddo
        enddo
c
c       y0(4,2,1): 2 starting solution vectors
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
        llw=lp
	exponent=0.d0
        do l=ls,lp-1
          n=nno(l)
	  exponent=exponent-k*hp(l)
          if(exponent.le.dlog(eps))then
            llw=l
            if(llw.lt.lzrec)return
            goto 200
          endif
        enddo
200     continue
c
        if(lup.eq.1)then
          yup(1,1)=1.d0
	  dro=ro(1)
	  yup(2,1)=dro*g0*yup(1,1)
          yup(3,2)=1.d0
        else
          c0(1,1)=1.d0
          c0(3,2)=1.d0
	  n=nno(lup-1)
          call edgmatrix(ma,4,k,0.d0,n)
          call axb(ma,c0,4,4,2,yup)
	  dro=ro(nno(lup))-ro(nno(lup-1))
	  do j=1,2
	    yup(2,j)=yup(2,j)+dro*g0*yup(1,j)
	  enddo
        endif
        if(lup.eq.lzrec)call memcpy(yup,y0,8)
c
        do l=lup+1,ls
          h0=hp(l-1)
	  n=nno(l-1)
	  if(k*h0.le.eps0)then
	    call edghask(hk,4,k,h0,n)
	    call axb(hk,yup,4,4,2,y1)
	    call memcpy(y1,yup,8)    
	  else
c
c
c           determination of propagation matrix
c
            call edgmatinv(mai,4,k,0.d0,n)
            call edgmatrix(ma,4,k,h0,n)
            call axb(mai,yup,4,4,2,c0)
	    wave=dexp(-k*h0)
c
c           normalization of all modes
c
	    do j=1,2
              norm(j)=0.d0
              do i=1,4
                norm(j)=norm(j)+c0(i,j)*c0(i,j)
              enddo
	    enddo
	    fac=1.d0/dsqrt(norm(1)*norm(2))
c
c           orthogonalization of the p-sv modes
c
            orth(1,1)=c0(3,2)*fac
            orth(1,2)=-c0(1,2)*fac
            orth(2,1)=-c0(3,1)*fac
            orth(2,2)=c0(1,1)*fac
            call axb(c0,orth,4,2,2,c1)
            if(l.gt.lzrec)then
	      do j=1,2
	        do i=1,2
	          orth(i,j)=orth(i,j)*wave
	        enddo
	      enddo
              call axb(y0,orth,4,2,2,y1)
              call memcpy(y1,y0,8)
            endif
c
c           c1(1,1)=c1(1,1)
            c1(2,1)=c1(2,1)*wave*wave
            c1(3,1)=(0.d0,0.d0)
            c1(4,1)=c1(4,1)*wave*wave
c
            c1(1,2)=(0.d0,0.d0)
            c1(2,2)=c1(2,2)*wave*wave
c           c1(3,2)=c1(3,2)
            c1(4,2)=c1(4,2)*wave*wave
c
            call axb(ma,c1,4,4,2,yup)
	  endif
	  dro=ro(nno(l))-ro(nno(l-1))
	  do j=1,2
	    yup(2,j)=yup(2,j)+dro*g0*yup(1,j)
	  enddo
          if(l.eq.lzrec)call memcpy(yup,y0,8)
        enddo
c
c===============================================================================
c
c       matrix propagation from half-space to source
c
        do i=1,4
          do j=1,2
            c0(i,j)=0.d0
            ylw(i,j)=0.d0
          enddo
        enddo
c
c
c       c0(4,2): 2 coefficient vectors in the half-space
c
	n=nno(llw)
        if(vs(n).lt.eps0*vp(n))then
c
c         the lowest layer is fluid
c
	  ylw(1,1)=1.d0
          ylw(3,2)=1.d0
        else
c
c         the lowest layer is solid
c
          c0(2,1)=1.d0
          c0(4,2)=1.d0
          call edgmatrix(ma,4,k,0.d0,n)
          call axb(ma,c0,4,4,2,ylw)
        endif
        if(llw.gt.ls.and.llw.eq.lzrec)call memcpy(ylw,y0,8)
c
        do l=llw-1,ls,-1
          h0=hp(l)
          n=nno(l)
	  dro=ro(nno(l+1))-ro(nno(l))
	  do j=1,2
	    ylw(2,j)=ylw(2,j)-dro*g0*ylw(1,j)
	  enddo
	  if(k*h0.le.eps0)then
	    call edghask(hk,4,k,-h0,n)
	    call axb(hk,ylw,4,4,2,y1)
	    call memcpy(y1,ylw,8)
	  else
c
c
c           determination of propagation matrix
c
            call edgmatinv(mai,4,k,0.d0,n)
            call edgmatrix(ma,4,k,-h0,n)
            call axb(mai,ylw,4,4,2,c0)
	    wave=dexp(-k*h0)
c
c           normalization of all modes
c
	    do j=1,2
              norm(j)=0.d0
              do i=1,4
                norm(j)=norm(j)+c0(i,j)*c0(i,j)
              enddo
	    enddo
	    fac=1.d0/dsqrt(norm(1)*norm(2))
c
c           orthogonalization of the p-sv modes
c
            orth(1,1)=c0(4,2)*fac
            orth(1,2)=-c0(2,2)*fac
            orth(2,1)=-c0(4,1)*fac
            orth(2,2)=c0(2,1)*fac
            call axb(c0,orth,4,2,2,c1)
            if(l.lt.lzrec)then
	      do j=1,2
	        do i=1,2
	          orth(i,j)=orth(i,j)*wave
	        enddo
	      enddo
              call axb(y0,orth,4,2,2,y1)
              call memcpy(y1,y0,8)
            endif
c
	    c1(1,1)=c1(1,1)*wave*wave
c	    c1(2,1)=c1(2,1)
	    c1(3,1)=c1(3,1)*wave*wave
            c1(4,1)=(0.d0,0.d0)
c
            c1(1,2)=c1(1,2)*wave*wave
            c1(2,2)=(0.d0,0.d0)
            c1(3,2)=c1(3,2)*wave*wave
c           c1(4,2)=c1(4,2)
c
            call axb(ma,c1,4,4,2,ylw)
	  endif
          if(l.gt.ls.and.l.eq.lzrec)call memcpy(ylw,y0,8)
        enddo
c
c===============================================================================
c
c       conditions on the source surface
c
c
c       source function
c
        do i=1,4
          if(kpower(i).eq.1)then
            fac=k
          else
            fac=1.d0
          endif
          b(i)=sfct(i)*fac
          do j=1,2
            coef(i,j)=yup(i,j)
            coef(i,j+2)=-ylw(i,j)
          enddo
        enddo
        key=0
	call gemp(coef,b,4,1,1.d-99,key)
        if(key.eq.0)then
          print *,'warning in edgpsv: anormal exit from cgemp!'
          return
        endif
        if(lzrec.le.ls)then
          do i=1,4
            y(i)=0.d0
            do j=1,2
              y(i)=y(i)+b(j)*y0(i,j)
            enddo
          enddo
        else
          do i=1,4
            y(i)=0.d0
            do j=1,2
              y(i)=y(i)+b(j+2)*y0(i,j)
            enddo
          enddo
        endif
c
        return
        end
