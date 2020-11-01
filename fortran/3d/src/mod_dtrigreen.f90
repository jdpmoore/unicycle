
!---------------------------------------------------------------------
!> subroutine dstuart
!! Calls Stuart's angle subroutine to get triangle dislocation results.
!! Comninou and Dundurs use a 1,2,3 = NED coordinate system.
!! Their simple angular dislocation is in the 1,3 plane
!! with its tip under the origin and its arm pointing toward +x1.
!!
!!     nu         = poissons ratio
!!     xo(3)      = observation point in global NED=123 coordinate system
!!     tridlc(9)  = x2,x1,x3 coordinates for the three triangle corners.
!!     ss,ds,op   = burg(3)    = burgers vector
!!                    (ss = + for LL; ds = + for normal; op = + for expansion)
!!     u(3)       = displacements = u1,u2,u3
!!     t(9)       = tilts = t11,t12,t13,t21,t22,t23,t31,t32,t33
!!                    (where tij = d(ui)/dxj)
!!
!! \author William Stuart (1997)
!---------------------------------------------------------------------
SUBROUTINE dstuart (nu, xo, tridlc, ss,ds,op, u,t)
  IMPLICIT NONE
  INTEGER :: i,k,iflag
  REAL(8) :: nu,ss,ds,op
  REAL(8) :: xo(3), tridlc(9), burg(3), u(3), t(9)
  REAL(8) :: utmp(3), ttmp(9)
  REAL(8) :: tri(3,4)
  INTEGER :: nwarn
  COMMON/warn0/ nwarn


  REAL(8), PARAMETER :: pi=3.14159265358979323846d0
  REAL(8), PARAMETER :: deg2rad=pi/180.d0
  REAL(8), PARAMETER :: rad2deg=180.d0/pi

  ! initialize displacement and tilt arrays.
  u = 0.d0
  t = 0.d0

  ! fill array with triangle corner coordinates.
  ! (Corner 4 is corner 1 repeated.)
  DO i=1,4
    DO  k=1,3
      tri(k,i) = tridlc(k+3*mod(i-1,3))
    END do
  END DO


  ! express burger's vector in global 123=NED coordinates.
  CALL burger2global (ss,ds,op, tri, burg)
  !PRINT '(a,6f8.3)', 'ss,op,ds,burg =(stu) ', ss,op,ds,burg

  ! for each of 3 sides, calculate fields.
  DO k = 1,3
     CALL twoangles (nu,xo,tri(1,k),tri(1,k+1),burg,utmp,ttmp)
     u = u + utmp
     t = t + ttmp
  END DO

  ! apply additional displacement if pt is under the triangle.
  CALL undertriangle (xo, tri, iflag)

  IF (iflag .EQ. 1) THEN
     ! point lies under interior of triangle.
     u = u - burg
  ELSE IF (iflag .ne. 0) then
     ! point lies on triangle or on edge of prism below triangle.
    nwarn = nwarn + 1
    IF (nwarn .LT. 10) then
      !PRINT *, ' *** Warning, possible singular location'
      !PRINT *, '      for point = ', xo
    ELSE IF (nwarn.EQ.11) THEN
      !PRINT *, ' *** Warning, possible singular location etc...'
    ENDIF
    u = u - 0.5*burg
  ENDIF

  RETURN
END SUBROUTINE dstuart

!---------------------------------------------------------------------
!> subroutine burger2global
!! Express burger's vector in global 123=NED coordinates.
!! Note special definition of ss and ds needed for horizontal dislocations
!! and of ds for vertical dislocations.
!!
!!     tri(3,i) = coordinates of ith corner in 123 coordinate system.
!!                   (i=4 returns to i=1)
!!
!---------------------------------------------------------------------
SUBROUTINE burger2global(ss,ds,op, tri, burg)
  IMPLICIT NONE
  REAL(8) :: ss,ds,op
  REAL(8) :: tri(3,3), burg(3)
  REAL(8) :: side12(3), side13(3)
  REAL(8) :: perp(3), hor(3), dwndip(3), b(3)
  REAL(8), PARAMETER :: pi=3.14159265358979323846d0
  REAL(8), PARAMETER :: deg2rad=pi/180.d0
  REAL(8), PARAMETER :: rad2deg=180.d0/pi
  REAL(8) :: strikerad,strike,dip

  ! for a horizontal triangle define strike- and dip-directions by fiat.
  IF (tri(3,1) .EQ. tri(3,2) .AND. tri(3,1) .EQ. tri(3,3)) THEN
    ! let 1 = ss-axis, 2 = ds-axis, 3 = op-axis
    burg(1) = ss
    burg(2) = ds
    burg(3) = -op
    RETURN
  ENDIF

  ! otherwise, construct a coordinate system attached to the triangle.

  ! calculate a normal vector.
  side12(:) = tri(:,2) - tri(:,1)
  side13(:) = tri(:,3) - tri(:,1)
  CALL cross(side12,side13,perp)
  CALL unit(perp)

  ! Vertical triangle.
  ! By definition, for a vertical triangle, side from which corners
  ! appear clockwise is considered top side.

  ! For dipping triangle...
  IF (perp(3) .LT. 0.d0) THEN
     ! Dipping triangle has CCW order to corners viewed from above.
     ! Switch corners 2 and 3 and reverse perp.
     CALL switch3(tri(1,2),tri(1,3))
     perp(:) = -1.d0*perp(:)
     ! Perp should now be a unit normal vector pointing down,
     ! (away from top), unless triangle is vertical.
  ENDIF

  ! Otherwise define a horizontal vector along strike.
  ! (Strike is to left hand as one stands looking down dip
  ! -- perp points opposite down-dip.)
  strikerad = ATAN2(perp(2),perp(1)) + pi*0.5d0
  hor(1) = COS(strikerad)
  hor(2) = SIN(strikerad)
  hor(3) = 0.d0

  ! Next a vector in the plane of the triangle pointing in the
  ! down-dip direction.
  call cross (perp,hor,dwndip)

  ! Express Burgers vector in this coordinate system.
  burg(:) = ss*hor(:) + ds*dwndip(:) - op*perp(:)

  ! Do a check.
  strike = strikerad * rad2deg
  dip = ACOS(-perp(3)) * rad2deg
  CALL global2triburger (burg, strike, dip, b)

  RETURN
END SUBROUTINE burger2global

!---------------------------------------------------------------------
!> subroutine global2triburger
!! Express burger's vector as b="ss,op,ds" (i.e., in local coordinates
!! attached to dlc.)  Order is here set to match Comninou & Dunders
!! conventions.
!---------------------------------------------------------------------
SUBROUTINE global2triburger (burg, strike, dip, b)
  IMPLICIT NONE
  REAL(8) :: burg(3), b(3)
  REAL(8) :: perp(3), hor(3), dwndip(3)
  REAL(8), parameter :: pi=3.14159265358979323846d0
  REAL(8), parameter :: deg2rad=pi/180.d0
  REAL(8), parameter :: rad2deg=180.d0/pi
  REAL(8) :: strike,dip,strikerad,diprad,dipdirrad

  ! Construct a coordinates system attached to the triangle.

  ! For a horizontal triangle define strike- and dip-directions by fiat.
  IF (dip .EQ. 0.d0) THEN
     ! Let 1 = ss-axis, 2 = op-axis, 3 = ds-axis
     strike = 0.d0
  ENDIF

  ! Define a horizontal vector along strike.
  strikerad = strike * deg2rad
  hor(1) = COS(strikerad)
  hor(2) = SIN(strikerad)
  hor(3) = 0.d0

  ! Next a vector in the plane of the triangle pointing in the
  ! down-dip dipdirection.
  diprad = dip * deg2rad
  dipdirrad = strikerad + pi*0.5
  dwndip(1) = COS(diprad) * COS(dipdirrad)
  dwndip(2) = COS(diprad) * SIN(dipdirrad)
  dwndip(3) = SIN(diprad)

  ! Finally a perpendicular vector
  CALL cross(hor,dwndip,perp)

  ! Get local burgers vector.  
  b(1) =   dot_product(burg,hor)
  b(2) = - dot_product(burg,perp)
  b(3) =   dot_product(burg,dwndip)

  RETURN
END SUBROUTINE global2triburger

!---------------------------------------------------------------------
!> subroutine twoangles
!! Calculates the effect of a vertical side formed by two angles
!! located at points xcorna and xcornb.
!!
!!     nu     = Poissons ratio
!!     xo(3)  = observation point in x1,x2,x3 = NED coordinates.
!!     xca(3) = x1,x2,x3 coordinates of corner A.
!!     xcb(3) = x1,x2,x3 coordinates of corner B.
!---------------------------------------------------------------------
SUBROUTINE twoangles (nu,xo,xcorna,xcornb, burg, u,t)
  IMPLICIT NONE
  REAL(8) :: nu, xo(3), xcorna(3), xcornb(3), burg(3), u(3), t(9)
  REAL(8) :: xop(3), xca(3), xcb(3), b(3)
  REAL(8) :: v(3,3), dv(3,3,3)
  REAL(8) :: va(3,3), dva(3,3,3)
  REAL(8) :: vb(3,3), dvb(3,3,3)

  REAL(8), PARAMETER :: pi=3.14159265358979323846d0
  REAL(8), PARAMETER :: deg2rad=pi/180.d0
  REAL(8), PARAMETER :: rad2deg=180.d0/pi
  INTEGER :: iswitch,j,k,iret
  REAL(8) :: theta,strike,dip,rot,deptha,depthb,xlen,betarad,x1,x2,x3

  ! Copy corners.
  xca(:) = xcorna(:)
  xcb(:) = xcornb(:)

  ! See which is shallower, A or B, and switch if necessary.
  IF (xca(3).gt.xcb(3)) THEN
     ! B is shallower than A
     CALL switch3 (xca, xcb)
     iswitch = 1
  ELSE
     iswitch = 0
  ENDIF

  ! Calculate angle from A to B.
  theta = ATAN2(xcb(2)-xca(2), xcb(1)-xca(1)) * rad2deg

  ! Set parameters.
  strike = theta
  deptha = xca(3)
  depthb = xcb(3)
  xlen = SQRT((xcb(1)-xca(1))**2 + (xcb(2)-xca(2))**2)
  betarad =  ATAN2(xlen, xcb(3)-xca(3))

  ! Locate the observation point relative to the shallower of A or B.
  xop(1) =  xo(1) - xca(1)
  xop(2) =  xo(2) - xca(2)
  xop(3) =  xo(3)
  ! If strike .ne. 0, need to rotate obs points and output fields.
  rot = -strike
  CALL vector_rot3(xop,-rot)
  x1 = xop(1)
  x2 = xop(2)
  x3 = xop(3)

  ! Compute the displacement and tilts.
  CALL comdun2 ( nu, x1,x2,x3, deptha, betarad, va, dva, iret )
  CALL comdun2 ( nu, x1-xlen,x2,x3, depthb, betarad, vb, dvb, iret )

  ! Combine the results from both angles.
  v(:,:) = va(:,:) - vb(:,:)
  dv(:,:,:) = dva(:,:,:) - dvb(:,:,:)

  ! Calculate displacements(u) and tilts(t)
  dip = 90.d0
  CALL global2triburger (burg, strike, dip, b)

  IF (iswitch .EQ. 1)  b(:) = -1d0*b(:)

  u = matmul(b,v)
  DO j=1,3
     DO k=1,3
        t(k + 3*(j-1)) = dot_product(b(:),dv(:,j,k))
     END DO
  END DO

  ! Rotate fields back to global NED=123 coordinate system.
  CALL vector_rot3(u,rot)
  CALL tensor_rot3(t,rot)

  RETURN
END SUBROUTINE twoangles

!---------------------------------------------------------------------
!> subroutine comdun2
!! Performs some checks before calling comdun
!---------------------------------------------------------------------
SUBROUTINE comdun2 ( nu, x1, x2, x3, a, beta, v, dv, iret )
  IMPLICIT NONE
  REAL(8) :: nu,x1,x2,x3,a,beta
  INTEGER :: iret,nwarn
  REAL(8) :: v(3,3), dv(3,3,3)
  REAL(8) :: vp(3,3), dvp(3,3,3)
  REAL(8) :: vm(3,3), dvm(3,3,3)
  REAL(8) :: betatol = 1.d-4
  REAL(8) :: tol = 1.d-4
  REAL(8) :: x2p,x2m,tan2angle

  COMMON/nwarn00/ nwarn

  ! Check for angle beta - if close to zero, return zero result.
  IF (beta .LT. betatol) THEN
     v(:,:) = 0.d0
     dv(:,:,:) = 0.d0
     RETURN
  ENDIF

  ! Check for distance from x1-x3 plane.  If too close, average.
  IF (ABS(x2) .LT. tol) THEN
     x2p = abs(x2) + tol
     x2m = -x2p
     CALL comdun(nu, x1, x2p, x3, a, beta, vp, dvp, iret)
     CALL comdun(nu, x1, x2m, x3, a, beta, vm, dvm, iret)
     v(:,:) = 0.5d0*(vp(:,:)+vm(:,:))
     dv(:,:,:) = 0.5d0*(dvp(:,:,:)+dvm(:,:,:))
  ELSE
     ! Normal call to comdun
     CALL comdun ( nu, x1, x2, x3, a, beta, v, dv, iret )
  ENDIF

  ! Check to see if point is inside angle.
  tan2angle = ATAN2(x1,x3-a)
  IF (ABS(x2) .LT.tol .AND. (tan2angle .LE. beta .AND. tan2angle .GE. 0.d0)) THEN
     ! Close to the displacement singularity.
     nwarn = nwarn + 1
     IF (nwarn.lt.10) THEN
        !PRINT *, '  *** Displacement singularity...'
     ELSE IF (nwarn .EQ. 11) THEN
       !PRINT *, '  *** Displacement singularity...etc,etc'
     ENDIF
     iret = 1
  ENDIF

  RETURN
END SUBROUTINE comdun2

!---------------------------------------------------------------------
!> subroutine twoangles
!!
!! Returns a flag = 1 if the point xo lies under the triangle interior.
!! Returns a flag = 2 if the point xo lies under a triangle edge or corner.
!! Returns a flag = 3 if the point xo lies on a triangle edge or corner.
!! Returns a flag = 0 otherwise.
!---------------------------------------------------------------------
SUBROUTINE undertriangle(xo, tri, iflag)
  IMPLICIT NONE
  INTEGER :: iflag
  REAL(8) :: xo(3), tri(3,4)
  REAL(8) :: xc(3), yc(3)
  REAL(8) :: side12(3), side13(3), perp(3)
  REAL(8) :: tol = 1.d-6
  REAL(8) :: xpt,ypt,zpt,zontri,ht
  INTEGER :: inside

  iflag=0
  ! Make arrays of x,y coords of triangle corners.
  xc(1:3) = tri(1,1:3)
  yc(1:3) = tri(2,1:3)

  xpt = xo(1)
  ypt = xo(2)
  zpt = xo(3)
  iflag = inside (xpt,ypt, xc,yc,3)

  ! Point is outside vertical prism....
  IF (iflag .EQ. 0) RETURN

  ! Point is inside vertical prism.... is it below triangle?

  ! Form perpendicular to the triangle.
  side12(:) = tri(:,2) - tri(:,1)
  side13(:) = tri(:,3) - tri(:,1)
  CALL cross(side12,side13,perp)
  CALL unit(perp)

  IF (perp(3) .EQ. 0.0) THEN
     ! Triangle is vertical, return edge flag.
     iflag = 2
     RETURN
  ENDIF

  ! Calculate height (+) below or (-) above triangle surface.
  zontri = (dot_product(tri(:,1),perp(:)) - xpt*perp(1) - ypt*perp(2))/perp(3)
  ht = zpt - zontri
  IF (abs(ht) .LT. tol) THEN
     iflag = 3
  ELSE IF (ht .LT. -tol) THEN
     iflag = 0
  ELSE
     iflag = abs(iflag)
  ENDIF

  RETURN
END SUBROUTINE undertriangle

!---------------------------------------------------------------------
!> function inside
!!
!!     From J.Berger, et al. BSSA v. 74, p. 1849-1862, October 1984.
!!     x0,y0 = point to test
!!     px,py,n = corners of polygon of n sides.
!!     Return value = 0 if point is outside
!!                  = +/-1 if point is inside
!!                  = 2 if point is on an edge or vertex.
!---------------------------------------------------------------------
INTEGER FUNCTION inside(x0,y0, px,py,n)
  INTEGER :: i
  INTEGER :: n
  REAL(8) :: x0,y0,px(n), py(n)

  inside = 0
  DO i=1,n-1
     isicr=ksicr(px(i)-x0,py(i)-y0,px(i+1)-x0,py(i+1)-y0)
     IF (isicr .EQ. 4) THEN
        inside = 2
        RETURN
     ENDIF
     inside = inside + isicr
  END DO

  isicr=ksicr(px(n)-x0,py(n)-y0,px(1)-x0,py(1)-y0)

  IF (isicr .EQ. 4) THEN
     inside = 2
     RETURN
  ENDIF

  inside = (inside + isicr)/2.
  RETURN
END FUNCTION inside


!---------------------------------------------------------------------
!> function ksicr
!---------------------------------------------------------------------
FUNCTION ksicr(x1,y1,x2,y2)
  IMPLICIT NONE
  INTEGER :: ksicr
  REAL(8) :: x1,x2,y1,y2
  IF (y1*y2 .GT. 0.0) GOTO 600

  IF (x1*y2 .NE. x2*y1 .OR. x1*x2 .GT. 0.0) GOTO 100

  ksicr = 4
  RETURN

100 IF (y1*y2 .LT. 0.0) GOTO 300
  IF (y2.eq.0.0) GOTO 200
  IF (x1.gt.0.0) GOTO 600
  IF (y2.gt.0.0) GOTO 700
  GOTO 800

200 IF (y1 .EQ. 0.0 .OR. x2 .GT. 0.0) GOTO 600
  IF (y1 .GT. 0.0) GOTO 800
  GOTO 700

300 IF (y1 .GT. 0.0) GOTO 400
  IF (x1*y2 .ge. y1*x2) GOTO 600
  ksicr = 2
  RETURN

400 IF (y1*x2 .GE. x1*y2) GOTO 600
  ksicr = -2
  RETURN

600 ksicr = 0
  RETURN

700 ksicr = 1
  RETURN

800 ksicr = -1
  RETURN

END FUNCTION ksicr

!---------------------------------------------------------------------
!> subroutine vector_rot3
!! Modified from Laurie Erickson's version of the DIS3D program.
!! Rotates vector u around the x3 axis in a CW sense by angle phi 
!! in the x1-x2 plane.
!---------------------------------------------------------------------
SUBROUTINE vector_rot3(u,phi)
  IMPLICIT NONE
  REAL(8) :: phi
  REAL(8) :: u(3), utmp(3)
  REAL(8), PARAMETER :: pi=3.14159265358979323846d0
  REAL(8), PARAMETER :: deg2rad=pi/180.d0
  REAL(8), PARAMETER :: rad2deg=180.d0/pi

  REAL(8) :: phirad,cp,sp

  IF (phi .EQ. 0.0) RETURN

  phirad = phi*deg2rad
  cp = COS(phirad)
  sp = SIN(phirad)

  utmp(1:2) = u(1:2)

  u(1) =  cp*utmp(1) + sp*utmp(2)
  u(2) = -sp*utmp(1) + cp*utmp(2)

  RETURN
END SUBROUTINE vector_rot3


!---------------------------------------------------------------------
!> subroutine tensor_rot3
!! Modified from Laurie Erickson's version of the DIS3D program.
!!
!! Rotates a general tensor t by the angle phi around the
!! x3 axis (within the x1-x2 plane).
!---------------------------------------------------------------------
SUBROUTINE tensor_rot3(t,phi)
  REAL(8) :: phi
  REAL(8) :: t(9), ttmp(9)
  REAL(8), PARAMETER :: pi=3.14159265358979323846d0
  REAL(8), PARAMETER :: deg2rad=pi/180.d0
  REAL(8), PARAMETER :: rad2deg=180.d0/pi
  REAL(8) :: phirad,cp,sp

  IF (phi .EQ. 0.d0) RETURN

  phirad = phi*deg2rad
  cp = cos(phirad)
  sp = sin(phirad)

  ttmp(1:8) = t(1:8)

  t(1) =  (cp**2)*ttmp(1)+cp*sp*ttmp(2)+cp*sp*ttmp(4)+(sp**2)*ttmp(5)
  t(2) = -cp*sp*ttmp(1)+(cp**2)*ttmp(2)-(sp**2)*ttmp(4)+sp*cp*ttmp(5)
  t(3) =  cp*ttmp(3) + sp*ttmp(6)
  t(4) = -sp*cp*ttmp(1)-(sp**2)*ttmp(2)+(cp**2)*ttmp(4)+cp*sp*ttmp(5)
  t(5) =  (sp**2)*ttmp(1)-sp*cp*ttmp(2)-sp*cp*ttmp(4)+(cp**2)*ttmp(5)
  t(6) = -sp*ttmp(3) + cp*ttmp(6)
  t(7) =  cp*ttmp(7) + sp*ttmp(8)
  t(8) = -sp*ttmp(7) + cp*ttmp(8)

  RETURN
END SUBROUTINE tensor_rot3

!---------------------------------------------------------------------
!> subroutine switch
!! Switches real numbers.
!---------------------------------------------------------------------
SUBROUTINE switch(a, b)
  IMPLICIT NONE
  REAL(8) :: a,b,temp
  temp = a
  a = b
  b = temp

  RETURN
END SUBROUTINE switch

!---------------------------------------------------------------------
!> subroutine switch3
!! Interchanges two vectors
!---------------------------------------------------------------------
SUBROUTINE switch3(xca, xcb)
  IMPLICIT NONE
  REAL(8) :: xca(3), xcb(3), xtmp(3)

  xtmp(:) = xca(:)
  xca(:) = xcb(:)
  xcb(:) = xtmp(:)

  RETURN
END SUBROUTINE switch3


!---------------------------------------------------------------------
!> subroutine cross
!! Calculates the cross product, acrsb, of vectors a and b.
!---------------------------------------------------------------------
SUBROUTINE cross(a,b,acrsb)
  IMPLICIT NONE
  REAL(8) :: a(3), b(3), acrsb(3)

  acrsb(1)= a(2)*b(3) - a(3)*b(2)
  acrsb(2)= a(3)*b(1) - a(1)*b(3)
  acrsb(3)= a(1)*b(2) - a(2)*b(1)

  RETURN
END SUBROUTINE cross

!---------------------------------------------------------------------
!> subroutine unit
!! Makes a unit vector out of x.
!---------------------------------------------------------------------
SUBROUTINE unit(x)
  IMPLICIT NONE
  REAL(8) :: x(3)
  REAL(8) :: d
  d = SQRT(dot_product(x,x))

  IF (d .GT. 0.d0) THEN
     x(:) = x(:)/d
  ELSE
     PRINT *, ' ** Warning... zero unit vector in sub unit.'
     STOP
  ENDIF

  RETURN
END SUBROUTINE unit

