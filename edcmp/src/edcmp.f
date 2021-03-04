!-----------------------------------------------------------------------
!       subroutine edcmp calculates convolution integral (summation)
!       of discrete sources to model the elastic deformations induced
!       by slip on a fault patch.
!
!        (xrec(i),yrec(i),zrec0)=coordinates of the observation positions
!       (Note that zrec0 is fixed)
!       disp = the 3 displcement vector components: ux,uy,uz
!       strain = the 6 strain tensor components: exx,eyy,ezz,exy,eyz,ezx
!       tilt = the two vertical tilt components: dux/dz, duy/dz
!       NRECMAX = the max. number of observation positions
!
!       (xs,ys,zs) = coordinates of the start point of strike
!       with x = north, y = east, z = downward.
!       all angles in degree.
!       NSMAX = the max. number of source rectangles
!
!       First implemented in Potsdam, Feb, 1999
!       Last modified: Potsdam, Nov, 2001, by R. Wang
!----------------------------------------------------------------------------------------------
        SUBROUTINE edcmp(islip,ixs,iys,izs,iL,iW,istrike,idip,irake,
     &                   nrec,ixrec,iyrec,ux,uy,uz,ierr)
        IMPLICIT NONE
        REAL*8, INTENT(IN) :: islip,ixs,iys,izs,iL,iW,
     &                        istrike,idip,irake
        INTEGER, INTENT(IN) :: nrec
        REAL*8, DIMENSION(nrec), INTENT(IN) :: ixrec,iyrec
        REAL*8, DIMENSION(nrec), INTENT(OUT) :: ux,uy,uz
        INTEGER, INTENT(OUT) :: ierr

        INCLUDE 'edcglobal.h'

        DOUBLE PRECISION dislocation(NSMAX)
        DOUBLE PRECISION xs(NSMAX),ys(NSMAX),zs(NSMAX)
        DOUBLE PRECISION length(NSMAX),width(NSMAX)
        DOUBLE PRECISION strike(NSMAX),dip(NSMAX),rake(NSMAX)
        COMMON/rectangles/dislocation,xs,ys,zs,length,width,
     &                    strike,dip,rake
        
        DOUBLE PRECISION xrec(NRECMAX),yrec(NRECMAX)
        DOUBLE PRECISION zrec0
        DOUBLE PRECISION disp(NRECMAX,3),strain(NRECMAX,6)
        DOUBLE PRECISION tilt(NRECMAX,2)
        COMMON/obsarray/xrec,yrec,zrec0,disp,strain,tilt
        
        INTEGER i,ns

        ierr=0
        IF (nrec .GT. NRECMAX) THEN
           ierr=1
           RETURN
        END IF

        ns=1
        dislocation(1)=islip
        xs(1)=ixs
        ys(1)=iys
        zs(1)=izs
        length(1)=iL
        width(1)=iW
        strike(1)=istrike
        dip(1)=idip
        rake(1)=irake

        DO i=1,nrec
           xrec(i)=ixrec(i)
           yrec(i)=iyrec(i)
        END DO

        CALL edcgrn(ns,nrec)

        DO i=1,nrec
           uy(i)= disp(i,1)
           ux(i)= disp(i,2)
           uz(i)=-disp(i,3)
        END DO

        END

