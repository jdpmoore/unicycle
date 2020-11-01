	subroutine edglayer(zs,ls,zrec,lzrec,h,n0)
        implicit none
c
c	First implemented in Potsdam, Feb, 1999
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
        integer ls,lzrec,n0
        double precision zs,zrec,h(n0)
c
	include 'edgglobal.h'
c
        integer lp,nno(nzmax)
        double precision hp(nzmax)
        common /sublayer/ hp,lp,nno
c
        integer l,n,li,lp0
        double precision z1
        double precision z(nzmax),z0(nzmax)
c
        lp0=1
        z0(lp0)=0.d0
        do n=1,n0-1
          lp0=lp0+1
          z0(lp0)=z0(lp0-1)+h(n)
        enddo
        lp0=lp0+1
        z0(lp0)=zrec
        lp0=lp0+1
        z0(lp0)=zs
c
c       sort the z0-profile
c
        do l=1,lp0-1
          do li=l+1,lp0
            if(z0(li).lt.z0(l))then
              z1=z0(l)
              z0(l)=z0(li)
              z0(li)=z1
            endif
          enddo
        enddo
c
c       delete duplicates
c
        lp=1
        z(lp)=0.d0
        do l=2,lp0
          if(z0(l).gt.z(lp))then
            hp(lp)=z0(l)-z(lp)
            lp=lp+1
            z(lp)=z0(l)
          endif
        enddo
        hp(lp)=0.d0
c
c       determine ls,lzrec
c
        do l=1,lp
          if(z(l).eq.zs)ls=l
          if(z(l).eq.zrec)lzrec=l
        enddo
c
c       determine layer no of each depth
c
        li=1
        z1=h(1)
        nno(1)=1
        do l=2,lp
          if(z(l).ge.z1.and.li.lt.n0)then
            li=li+1
            z1=z1+h(li)
          endif
          nno(l)=li
        enddo
c
        return
        end
