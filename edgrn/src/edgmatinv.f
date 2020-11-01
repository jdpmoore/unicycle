        subroutine edgmatinv(a,m,k,z,n)
        implicit none
c
c	First implemented in Potsdam, Feb, 1999
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
        integer m,n
        double precision k,z
        double precision a(m,m)
c
	include 'edgglobal.h'
c
c       model parameter:
c       n0: number of model layers
c
        integer n0
        double precision h(lmax),ro(lmax),vp(lmax),vs(lmax)
        common /model/ h,ro,vp,vs,n0
c
        integer i,j
        double precision x,la,mu,eta,alfa,psi
c
        if(m.ne.4)then
          stop ' Error in edgmatinv: n != 4'
        endif
        x=k*z
	mu=ro(n)*vs(n)**2
	la=ro(n)*vp(n)**2-2.d0*mu
        eta=mu/(la+mu)
        alfa=2.d0*mu*k
        psi=1.d0+eta
        a(1,1)=x-eta
        a(1,2)=x/alfa
        a(1,3)=-x+psi
        a(1,4)=(psi+eta-x)/alfa
        a(2,1)=x+eta
        a(2,2)=-x/alfa
        a(2,3)=x+psi
        a(2,4)=-(psi+eta+x)/alfa
        a(3,1)=1.d0
        a(3,2)=1.d0/alfa
        a(3,3)=-1.d0
        a(3,4)=-1.d0/alfa
        a(4,1)=1.d0
        a(4,2)=-1.d0/alfa
        a(4,3)=1.d0
        a(4,4)=-1.d0/alfa
        do i=1,4
          do j=1,4
            a(i,j)=0.5d0*a(i,j)/psi
          enddo
        enddo
c
        return
        end
