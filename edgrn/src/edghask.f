c*******************************************************************************
c*******************************************************************************
        subroutine edghask(hk,m,k,z,n)
        implicit none
c
c	First implemented in Potsdam, Feb, 1999
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
        integer m,n
        double precision k,z
        double precision hk(m,m)
c
	include 'edgglobal.h'
c
c
c       n0: number of model layer
c
        integer n0
        double precision h(lmax),ro(lmax),vp(lmax),vs(lmax)
        common /model/ h,ro,vp,vs,n0
c
	double precision eps0
	parameter(eps0=1.0d-03)
c
	double precision k2,x,x2,la,mu
        double precision ex,ch,sh,eta
c
        k2=k*k
	x=k*z
	x2=x*x
	mu=ro(n)*vs(n)**2
	la=ro(n)*vp(n)**2-2.d0*mu
	ex=dexp(x)
c
c	ch=(e^x+1/e^x)/2
c	sh=(e^x-1/e^x)/2x
c
	ch=0.5d0*(ex+1.d0/ex)
	if(dabs(x).le.eps0)then
	  sh=1.d0+x2/6.d0*(1.d0+x2/20.d0)
	else
	  sh=0.5d0*(ex-1.d0/ex)/x
	endif
c
	if(m.eq.2)then
c
c	  propagator matrix for SH waves
c
	  hk(1,1)=ch
	  hk(1,2)=z*sh/mu
	  hk(2,1)=k*x*sh*mu
	  hk(2,2)=ch
	else if(m.eq.4)then
c
c	  propagatior matrix for P-SV waves.
c
	  eta=mu/(la+mu)
	  hk(1,1)=ch-x2*sh/(1.d0+eta)
	  hk(1,2)=0.5d0*z*(-ch+(1.d0+2.d0*eta)*sh)/(1.d0+eta)/mu
	  hk(1,3)=x*(ch-eta*sh)/(1.d0+eta)
	  hk(1,4)=x2*sh/(1.d0+eta)
	  hk(2,1)=2.d0*mu*k*x*(-ch+sh)/(1.d0+eta)
	  hk(2,2)=hk(1,1)
	  hk(2,3)=2.d0*mu*k*hk(1,4)
	  hk(2,4)=x*(ch+eta*sh)/(1.d0+eta)
	  hk(3,1)=-hk(2,4)
	  hk(3,2)=-0.5d0*x*z*sh/(1.d0+eta)/mu
	  hk(3,3)=ch+x2*sh/(1.d0+eta)
	  hk(3,4)=0.5d0*z*(ch+(1.d0+2.d0*eta)*sh)/(1.d0+eta)/mu
	  hk(4,1)=-hk(2,3)
	  hk(4,2)=-hk(1,3)
	  hk(4,3)=2.d0*mu*k*x*(ch+sh)/(1.d0+eta)
	  hk(4,4)=hk(3,3)
	else
	  stop ' Error in haskell: m schould be 2 or 4!'
	endif
c
	return
	end
