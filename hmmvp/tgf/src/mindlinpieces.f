cc Copyright (C) 2009-2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Direct calculation of Mindlin B and C parts for half space
c       elastostatic Green's functions in R^3.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c***********************************************************************
      subroutine intker_mindlinb(icomp,rlame,source,ifsingle,
     1           sigma_sl,ifdouble,dipstr,dipvec,ztrg,pot,iffld,fld)
c***********************************************************************
c
c     This subroutine computes the displacement due to the 
c     harmonic part of the Mindlin B solution and all of its
c     derivatives, from which stress and strain can be extracted.
c     (See okada_imageanalysis for details.)
c     It doesn't include the 1/rmu scaling.
c
c     INPUT:
c
c     icomp      desired component of displacement vector
c     rlame(2)   Lame coefficients supplied in the form
c                rlame(1) = rlam, rlame(2) = rmu
c     source(3)  image source location.
c                This source must be the reflected image in the 
c                upper half-space, NOT the original source in the 
c                lower half-space.
c     ifsingle   (the single layer kernel flag)
c                0 means ignore the <<force>> vector 
c                1 means a <<force>> vector is supplied 
c     sigma_sl   force vector
c     ifdouble   (the double layer kernel flag)
c                0 means ignore the <<dislocation>> 
c                1 means a <<dislocation vector>> is supplied 
c                  with surface normal <<dipvec>>
c     dipstr     dislocation vector
c     dipvec     surface orientation
c                
c     ztrg       target location
c     pot        component <<icomp>> of displacement vector
c     iffld      fld computation flag: 
c                0 means don't compute
c                1 means compute fld = -gradient
c
c     fld        -gradient(pot) at ztrg
c-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: icomp
      real(8), intent(in) :: rlame(2)
      real(8), intent(in) :: source(3)
      integer, intent(in) :: ifsingle
      real(8), intent(in) :: sigma_sl(3)
      integer, intent(in) :: ifdouble
      real(8), intent(in) :: dipstr(3)
      real(8), intent(in) :: dipvec(3)
      real(8), intent(in) :: ztrg(3)
      complex(8), intent(out) :: pot
      integer, intent(in) :: iffld
      complex(8), intent(inout) :: fld(3)

!!      integer :: nder(3)
      integer :: j
      integer :: k
!!      integer :: l
      integer :: n
      real(8) :: uder(3,3,3)
      real(8) :: uuu
      real(8) :: uuuder(3)
      real(8) :: udertot(0:4,0:4,0:4)
      real(8) :: denom1
      real(8) :: denom2
      real(8) :: rr1
      real(8) :: rr2
      real(8) :: rr3
      real(8) :: cd
      real(8) :: rlam
      real(8) :: rmu
      real(8) :: scale
      real(8) :: sgn
      real(8) :: u11
      real(8) :: u12
      real(8) :: u13
      real(8) :: u22
      real(8) :: u23
      real(8) :: u33
!!      real(8) :: uder1
!!      real(8) :: uder2
!!      real(8) :: uder4
      real(8) :: alpha
      real(8) :: d1
      real(8) :: d2
      real(8) :: dotprod
      real(8) :: dx
      real(8) :: dy
      real(8) :: dz
      real(8) :: cd__2
      real(8) :: cd__3
      real(8) :: cd__5
      real(8) :: rr1__2
      real(8) :: rr2__2
      real(8) :: rr3__2
      real(8) :: rr3__3
      real(8) :: rr3__4
      real(8) :: rr3__5
      real(8) :: rr3__6
      real(8) :: rr3__7
      real(8) :: dx__2
      real(8) :: dx__3
      real(8) :: dx__4
      real(8) :: dx__6
      real(8) :: dy__2
      real(8) :: dy__3
      real(8) :: dy__4
      real(8) :: dy__6
      real(8) :: term1
      real(8) :: term1__2
      real(8) :: term1__3
      real(8) :: term1__4
      real(8) :: term2
      real(8) :: term2__3
      real(8) :: term3
      real(8) :: term3__2
      real(8) :: dipcross(3,3)
c      
c    
c
      rlam = rlame(1)
      rmu = rlame(2)
      alpha = (rlam+rmu)/(rlam+2*rmu)
      scale = (1.0d0-alpha)/alpha
c
      pot = 0.0d0
      if (iffld.eq.1) then
         fld(1) = 0.0d0
         fld(2) = 0.0d0
         fld(3) = 0.0d0
      endif
c
      dx = ztrg(1)-source(1)
      dy = ztrg(2)-source(2)
      dz = ztrg(3)-source(3)
      rr1 = dx
      rr2 = dy
      rr3 = -dz
!!      cd = sqrt(dx*dx+dy*dy+dz*dz)

      cd__2 = dx*dx+dy*dy+dz*dz
      cd = sqrt(cd__2)
      cd__3 = cd * cd__2
      cd__5 = cd__3 * cd__2

      rr1__2 = rr1 * rr1

      rr2__2 = rr2 * rr2

      rr3__2 = rr3 * rr3
      rr3__3 = rr3__2 * rr3
      rr3__4 = rr3__3 * rr3
      rr3__5 = rr3__4 * rr3
      rr3__6 = rr3__5 * rr3
      rr3__7 = rr3__6 * rr3

      dx__2 = dx * dx
      dx__3 = dx__2 * dx
      dx__4 = dx__2 * dx__2
      dx__6 = dx__3 * dx__3

      dy__2 = dy * dy
      dy__3 = dy__2 * dy
      dy__4 = dy__2 * dy__2
      dy__6 = dy__3 * dy__3

      term1 = cd+rr3
      term1__2 = term1 * term1
      term1__3 = term1__2 * term1
      term1__4 = term1__3 * term1

      term2 = dx__2+dy__2+rr3__2
      term2__3 = term2 * term2 * term2

      term3 = cd-dz
      term3__2 = term3 * term3
c
      if ( (ifdouble.eq.1) .or. (iffld.eq.1) ) then
c
c     Compute three derivatives of harmonic function 
c     U = R_3 log(R+R_3) - R
c     with respect to R_1, R_2 and R_3. We will use these
c     derivatives to compute both the gradient of the SLP
c     and the DLP.
c
c
!!         denom1 = (term3*cd)**3
         denom1 = (term3*cd)
         denom1 = denom1 * denom1 * denom1
         denom2 = (term3__2)*cd__3
         uder(1,1,1) = 3*dx*cd*(term3*cd-dx__2) + dz*dx__3
         uder(1,1,1) = uder(1,1,1)/denom1
         uder(1,1,2) = term3*dy*(cd__2-dx__2)-2*cd*dy*dx__2
         uder(1,1,2) = uder(1,1,2)/denom1
         uder(1,1,3) = term3*(cd__2-dx__2) - dx__2*cd
         uder(1,1,3) = uder(1,1,3)/denom2
         uder(1,2,1) = uder(1,1,2)       
         uder(1,2,2) = term3*dx*(cd__2-dy__2)-2*cd*dx*dy__2
         uder(1,2,2) = uder(1,2,2)/denom1
         uder(1,2,3) = -dx*dy*(2*cd-dz)
         uder(1,2,3) = uder(1,2,3)/denom2
         uder(1,3,1) = uder(1,1,3)       
         uder(1,3,2) = uder(1,2,3)       
         uder(1,3,3) = -dx*term3__2
         uder(1,3,3) = uder(1,3,3)/denom2
c
         if (icomp.ne.1) then
            uder(2,1,1) = uder(1,1,2)
            uder(2,1,2) = uder(1,2,2)
            uder(2,1,3) = uder(1,2,3)
            uder(2,2,1) = uder(1,2,2)
            uder(2,2,2) = 3*dy*cd*(term3*cd-dy__2) + dz*dy__3
            uder(2,2,2) = uder(2,2,2)/denom1
            uder(2,2,3) = term3*(cd__2-dy__2) - dy__2*cd
            uder(2,2,3) = uder(2,2,3)/denom2
            uder(2,3,1) = uder(1,2,3)
            uder(2,3,2) = uder(2,2,3)
            uder(2,3,3) = -dy*term3__2
            uder(2,3,3) = uder(2,3,3)/denom2
         endif
c
         if (icomp.eq.3) then
            uder(3,1,1) = uder(1,1,3)
            uder(3,1,2) = uder(1,2,3)
            uder(3,1,3) = uder(1,3,3)
            uder(3,2,1) = uder(1,2,3)
            uder(3,2,2) = uder(2,2,3)
            uder(3,2,3) = uder(2,3,3)
            uder(3,3,1) = uder(1,3,3)
            uder(3,3,2) = uder(2,3,3)
            uder(3,3,3) = dz/(cd__3)
         endif
      endif
c
      if (ifsingle.eq.1) then
c
c
c     Compute two derivatives of harmonic function 
c     U = R_3 log(R+R_3) - R
c     with respect to R_1, R_2 and R_3. These play a role only 
c     in the SLP.
c
         d1 = 1.0d0/term1
         d2 = 1.0d0/(cd*term1)
         u11 = d1*(-1.0d0 + d2*rr1__2)
         u12 = d1*d2*rr1*rr2
         u13 = d2*rr1
         u22 = d1*(-1.0d0 + d2*rr2__2)
         u23 = d2*rr2
         u33 = 1.0d0/cd
c
c     See notes (okada_imageanalysis) for explanation.
c
         if (icomp.eq.1) then
            pot = -sigma_sl(1)*u11 - sigma_sl(2)*u12 + sigma_sl(3)*u13
            pot = pot*scale
         else if (icomp.eq.2) then
            pot = -sigma_sl(1)*u12 - sigma_sl(2)*u22 + sigma_sl(3)*u23
            pot = pot*scale
         else 
            pot = -sigma_sl(1)*u13 - sigma_sl(2)*u23 + sigma_sl(3)*u33
            pot = pot*scale
         endif 
c
c      compute contribution to gradient of SLP for component icomp
c
c      Because Laplace FMM is based on FIELD rather than GRADIENT,
c      we compute FIELD = -GRADIENT here. We just take (-gradient)
c      of preceding formula. Since
c      d/dx_1 = d/dR_1, d/dx_2 = d/dR_2 and d/dx_3 = -d/dR_3
c      fld(1) and fld(2) and fld(3) use different signs below.
c
         if (iffld.eq.1) then
            sgn = scale
            do k = 1,3
               if (k.eq.3) sgn = -sgn
               fld(k) = fld(k) + sgn*(
     1          uder(icomp,1,k)*sigma_sl(1) +
     1          uder(icomp,2,k)*sigma_sl(2) -
     1          uder(icomp,3,k)*sigma_sl(3))
            enddo
         endif
      endif
c
      if (ifdouble.eq.1) then
c
c     See notes (okada_imageanalysis) for explanation of sign
c     flip with j in loop below.
c
         do j = 1,3
         do k = 1,3
            uder(icomp,j,k)=uder(icomp,j,k)*scale
            if (j.eq.3) uder(icomp,j,k) = -uder(icomp,j,k)
         enddo
         enddo
c
         uuu = 0
         do n=1,3
            uuu = uuu+uder(icomp,n,n)
         enddo
c
         do j=1,3
         do k=1,3
            dipcross(j,k) = dipstr(j)*dipvec(k)
         enddo
         enddo

         do j=1,3
         do k=1,3
            pot = pot + 
     1        rmu*(uder(icomp,j,k)+uder(icomp,k,j))*dipcross(j,k)
         enddo
         enddo
         dotprod = dipcross(1,1) + 
     1             dipcross(2,2) + 
     1             dipcross(3,3) 
         pot = pot + rlam*uuu*dotprod
         if (iffld.eq.1) then
            udertot(4,0,0) = scale*(
     #-3*(2*rr3*dx__6+6*dx__4*rr3__3+4*cd*dx__4*dy__2+6*rr3
     #*dx__4*dy__2+4*rr3__2*dx__4*cd+2*rr3*dx__2*dy__4+4*dx__2*
     #dy__2*rr3__3+3*rr3__4*dx__2*cd+2*dx__2*rr3__5+6*rr3__2*dx__2*
     #cd*dy__2+3*cd*dx__2*dy__4-cd*
     #dy__6-2*rr3__6*cd-4*rr3__2*cd*
     #dy__4-5*rr3__4*cd*dy__2-2*rr3*dy__6-6*dy__4*rr3__3
     #-6*dy__2*rr3__5-2*rr3__7)/term1__4/term2__3)
         udertot(3,1,0) = scale*(
     # 3*(2*cd*dx__4-dx__2*cd*dy__2
     #-4*rr3__3*dx__2-4*rr3*dx__2*dy__2-rr3__2*dx__2*cd-4*rr3*dy__4
     #-8*rr3__3*dy__2-7*rr3__2*cd*dy__2-4*rr3__4*cd
     #-4*rr3__5-3*cd*dy__4)*dx*dy/term1__4/term2__3)
         udertot(3,0,1) = scale*(
     #dx*(2*cd*dx__4-9*rr3*dx__2*dy__2-9*rr3__3*dx__2-4*
     #dx__2*cd*dy__2-4*rr3__2*dx__2*cd-6*cd*dy__4-9*rr3__5-
     #9*rr3__4*cd-15*rr3__2
     #*cd*dy__2-9*rr3*dy__4-18*rr3__3*dy__2)/
     #term1__3/term2__3)
         udertot(2,2,0) = scale*(
     #-(2*cd*dx__6+2*rr3*dx__6-9*cd
     #*dx__4*dy__2+3*rr3__2*dx__4*cd+2*dx__4*rr3__3-6*rr3*dx__4*
     #dy__2-12*rr3__2*dx__2*cd*dy__2-8*dx__2*dy__2*rr3__3-rr3__4*
     #dx__2*cd-9*cd*dx__2*dy__4-6*rr3*dx__2
     #*dy__4-2*dx__2*rr3__5-2*dy__2*rr3__5-2*rr3__6*cd-rr3__4*cd
     #*dy__2+2*cd*dy__6+2*dy__4*rr3__3+2*rr3
     #*dy__6+3*rr3__2*cd*dy__4-2*rr3__7)/
     #term1__4/term2__3)
         udertot(2,1,1) = scale*(
     #dy*(6*cd*dx__4+6*rr3*dx__4+3*rr3__3*dx__2+6*rr3__2*
     #dx__2*cd+3*rr3*dx__2*dy__2+4*dx__2*cd*
     #dy__2-3*rr3__4*cd-2*cd*dy__4-5*rr3__2
     #*cd*dy__2-6*rr3__3*dy__2-3*rr3__5-3*rr3*dy__4)/
     #term1__3/term2__3)
         udertot(2,0,2) = scale*(
     #(2*cd*dx__4+4*rr3*dx__4+2*rr3*dx__2*dy__2+3*rr3__2*
     #dx__2*cd+2*rr3__3*dx__2+dx__2*cd*dy__2
     #-2*rr3__4*cd-3*rr3__2*cd*dy__2-2*rr3*
     #dy__4-4*rr3__3*dy__2-2*rr3__5-cd*dy__4)/
     #term2__3/term1__2)
         udertot(1,3,0) = scale*(
     #-3*(4*rr3*dx__4+3*cd*dx__4+dx__2*cd*
     #dy__2+7*rr3__2*dx__2*cd+4*rr3*dx__2*dy__2+8*rr3__3*
     #dx__2+rr3__2*cd*dy__2-2*cd*dy__4+4*
     #rr3__4*cd+4*rr3__3*dy__2+4*rr3__5)*dx*dy/term1__4/
     #term2__3)
         udertot(1,2,1) = scale*(
     # -(2*cd*dx__4+3*rr3*dx__4+5*rr3__2*dx__2*cd+6*rr3__3*dx__2-
     #4*dx__2*cd*dy__2-3*rr3*dx__2*dy__2
     #+3*rr3__5-6*rr3*dy__4+3*rr3__4*cd-6*cd*
     #dy__4-6*rr3__2*cd*dy__2-3*rr3__3*dy__2)*dx/term1__3/
     #term2__3)
         udertot(1,1,2) = scale*(
     #3*(cd*dx__2+2*rr3*dx__2+2*cd*
     #rr3__2+2*rr3*dy__2+cd*dy__2+2*rr3__3)*dx*dy/
     #term2__3/term1__2)
         udertot(1,0,3) = scale*3/cd__5*dx*rr3
         udertot(0,4,0) = scale*(
     #3*(cd*dx__6+2*rr3*dx__6+6*dx__4*rr3__3+4*rr3__2*dx__4
     #*cd-3*cd*dx__4*dy__2-2*rr3*dx__4*
     #dy__2-6*rr3__2*dx__2*cd*dy__2+6*dx__2*rr3__5+5*rr3__4*dx__2*
     #sqrt(term2)-6*rr3*dx__2*dy__4-4*cd*dx__2*dy__4
     #-4*dx__2*dy__2*rr3__3-6*dy__4*rr3__3-2*rr3*dy__6-4*rr3__2*cd
     #*dy__4+2*rr3__6*cd+2*rr3__7-3*rr3__4*cd
     #*dy__2-2*dy__2*rr3__5)/term1__4/
     #term2__3)
         udertot(0,3,1) = scale*(
     # -(9*rr3*dx__4+6*cd*dx__4+18*rr3__3*dx__2+15*rr3__2*
     #dx__2*cd+9*rr3*dx__2*dy__2+4*dx__2*cd
     #*dy__2+9*rr3__5+9*rr3__4*cd-2*cd*
     #dy__4+4*rr3__2*cd*dy__2+9*rr3__3*dy__2)*dy/term1__3/
     #term2__3)
         udertot(0,2,2) =  scale*(
     #-(2*rr3*dx__4+cd*dx__4+4*rr3__3*dx__2+3*rr3__2*dx__2*
     #cd-2*rr3*dx__2*dy__2-dx__2*cd*dy__2-
     #2*cd*dy__4-4*rr3*dy__4+2*rr3__5-3*rr3__2*cd
     #*dy__2+2*rr3__4*cd-2*rr3__3*dy__2)/
     #term2__3/term1__2)
         udertot(0,1,3) =  scale*3/cd__5*dy*rr3
         udertot(0,0,4) = -scale*(dx__2+dy__2-2*rr3__2)/cd__5
c
!!         do l = 1,3
!!            uuuder(l) = 0
!!            do n=1,3
!!               nder(1) = 0
!!               nder(2) = 0
!!               nder(3) = 0
!!               nder(icomp) = nder(icomp)+1
!!               nder(l) = nder(l)+1
!!               nder(n) = nder(n)+2
!!c
!!c              uder4 should be d_{x_l} (d u_{icomp}^n / d \xi_n)
!!c
!!               uder4 = -udertot(nder(1),nder(2),nder(3))
!!               if (n.eq.3) uder4 = -uder4
!!               if (l.eq.3) uder4 = -uder4
!!               uuuder(l) = uuuder(l)+uder4
!!ccc               write(6,*)' icomp,l,n',icomp,l,n
!!ccc               write(6,*)' nder',nder
!!ccc               write(6,*)' uder4 =',uder4
!!            enddo
!!         enddo


!        In the following, where you see subscripts like 0+1+2,
!        the first number comes from icomp, the second from l,
!        and the third from n, where n gets a 2 and the others a 1.

         if (icomp.eq.1) then

         ! l = 1

         uuuder(1) =
     &     - udertot(1+1+2,0+0+0,0+0+0)
     &     - udertot(1+1+0,0+0+2,0+0+0)
     &     + udertot(1+1+0,0+0+0,0+0+2)

         ! l = 2

         uuuder(2) =
     &     - udertot(1+0+2,0+1+0,0+0+0)
     &     - udertot(1+0+0,0+1+2,0+0+0)
     &     + udertot(1+0+0,0+1+0,0+0+2)

         ! l = 3

         uuuder(3) =
     &       udertot(1+0+2,0+0+0,0+1+0)
     &     + udertot(1+0+0,0+0+2,0+1+0)
     &     - udertot(1+0+0,0+0+0,0+1+2)

         else if (icomp.eq.2) then

         ! l = 1

         uuuder(1) =
     &     - udertot(0+1+2,1+0+0,0+0+0)
     &     - udertot(0+1+0,1+0+2,0+0+0)
     &     + udertot(0+1+0,1+0+0,0+0+2)

         ! l = 2

         uuuder(2) =
     &     - udertot(0+0+2,1+1+0,0+0+0)
     &     - udertot(0+0+0,1+1+2,0+0+0)
     &     + udertot(0+0+0,1+1+0,0+0+2)

         ! l = 3

         uuuder(3) =
     &       udertot(0+0+2,1+0+0,0+1+0)
     &     + udertot(0+0+0,1+0+2,0+1+0)
     &     - udertot(0+0+0,1+0+0,0+1+2)

         else if (icomp.eq.3) then

         ! l = 1

         uuuder(1) =
     &     - udertot(0+1+2,0+0+0,1+0+0)
     &     - udertot(0+1+0,0+0+2,1+0+0)
     &     + udertot(0+1+0,0+0+0,1+0+2)

         ! l = 2

         uuuder(2) =
     &     - udertot(0+0+2,0+1+0,1+0+0)
     &     - udertot(0+0+0,0+1+2,1+0+0)
     &     + udertot(0+0+0,0+1+0,1+0+2)

         ! l = 3

         uuuder(3) =
     &       udertot(0+0+2,0+0+0,1+1+0)
     &     + udertot(0+0+0,0+0+2,1+1+0)
     &     - udertot(0+0+0,0+0+0,1+1+2)

         endif


c
!!         do l=1,3
!!         do j=1,3
!!         do k=1,3
!!            nder(1) = 0
!!            nder(2) = 0
!!            nder(3) = 0
!!            nder(icomp) = nder(icomp)+1
!!            nder(l) = nder(l)+1
!!            nder(j) = nder(j)+1
!!            nder(k) = nder(k)+1
!!            uder1 = -udertot(nder(1),nder(2),nder(3))
!!            uder2 = -udertot(nder(1),nder(2),nder(3))
!!            if (k.eq.3) uder2 = -uder2
!!            if (j.eq.3) uder1 = -uder1
!!            if (l.eq.3) uder1 = -uder1
!!            if (l.eq.3) uder2 = -uder2
!!            fld(l) = fld(l) + 
!!     1               rmu*(uder1+uder2)*dipcross(j,k)
!!         enddo
!!         enddo
!!         fld(l)=fld(l)+rlam*uuuder(l)*dotprod
!!         enddo


!        In the following, where you see subscripts like 0+1+2,
!        the first number comes from icomp, the second from l,
!        and the third from j and k.

         if (icomp.eq.1) then

         ! l = 1
         
         fld(1) = fld(1) - rmu*2.0d0*(
     &      udertot(1+1+2,0+0+0,0+0+0)*dipcross(1,1)
     &    + udertot(1+1+1,0+0+1,0+0+0)*(dipcross(1,2) + dipcross(2,1))
     &    + udertot(1+1+0,0+0+2,0+0+0)*dipcross(2,2)
     &    - udertot(1+1+0,0+0+0,0+0+2)*dipcross(3,3) )
     &    + rlam*uuuder(1)*dotprod

         ! l = 2
         
         fld(2) = fld(2) - rmu*2.0d0*(
     &      udertot(1+0+2,0+1+0,0+0+0)*dipcross(1,1)
     &    + udertot(1+0+1,0+1+1,0+0+0)*(dipcross(1,2) + dipcross(2,1))
     &    + udertot(1+0+0,0+1+2,0+0+0)*dipcross(2,2)
     &    - udertot(1+0+0,0+1+0,0+0+2)*dipcross(3,3) )
     &    + rlam*uuuder(2)*dotprod

         ! l = 3
         
         fld(3) = fld(3) + rmu*2.0d0*(
     &      udertot(1+0+2,0+0+0,0+1+0)*dipcross(1,1)
     &    + udertot(1+0+1,0+0+1,0+1+0)*(dipcross(1,2) + dipcross(2,1))
     &    + udertot(1+0+0,0+0+2,0+1+0)*dipcross(2,2)
     &    - udertot(1+0+0,0+0+0,0+1+2)*dipcross(3,3) )
     &    + rlam*uuuder(3)*dotprod

         else if (icomp.eq.2) then

         ! l = 1
         
         fld(1) = fld(1) - rmu*2.0d0*(
     &      udertot(0+1+2,1+0+0,0+0+0)*dipcross(1,1)
     &    + udertot(0+1+1,1+0+1,0+0+0)*(dipcross(1,2) + dipcross(2,1))
     &    + udertot(0+1+0,1+0+2,0+0+0)*dipcross(2,2)
     &    - udertot(0+1+0,1+0+0,0+0+2)*dipcross(3,3) )
     &    + rlam*uuuder(1)*dotprod

         ! l = 2
         
         fld(2) = fld(2) - rmu*2.0d0*(
     &      udertot(0+0+2,1+1+0,0+0+0)*dipcross(1,1)
     &    + udertot(0+0+1,1+1+1,0+0+0)*(dipcross(1,2) + dipcross(2,1))
     &    + udertot(0+0+0,1+1+2,0+0+0)*dipcross(2,2)
     &    - udertot(0+0+0,1+1+0,0+0+2)*dipcross(3,3) )
     &    + rlam*uuuder(2)*dotprod

         ! l = 3
         
         fld(3) = fld(3) + rmu*2.0d0*(
     &      udertot(0+0+2,1+0+0,0+1+0)*dipcross(1,1)
     &    + udertot(0+0+1,1+0+1,0+1+0)*(dipcross(1,2) + dipcross(2,1))
     &    + udertot(0+0+0,1+0+2,0+1+0)*dipcross(2,2)
     &    - udertot(0+0+0,1+0+0,0+1+2)*dipcross(3,3) )
     &    + rlam*uuuder(3)*dotprod

         else if (icomp.eq.3) then

         ! l = 1
         
         fld(1) = fld(1) - rmu*2.0d0*(
     &      udertot(0+1+2,0+0+0,1+0+0)*dipcross(1,1)
     &    + udertot(0+1+1,0+0+1,1+0+0)*(dipcross(1,2) + dipcross(2,1))
     &    + udertot(0+1+0,0+0+2,1+0+0)*dipcross(2,2)
     &    - udertot(0+1+0,0+0+0,1+0+2)*dipcross(3,3) )
     &    + rlam*uuuder(1)*dotprod

         ! l = 2
         
         fld(2) = fld(2) - rmu*2.0d0*(
     &      udertot(0+0+2,0+1+0,1+0+0)*dipcross(1,1)
     &    + udertot(0+0+1,0+1+1,1+0+0)*(dipcross(1,2) + dipcross(2,1))
     &    + udertot(0+0+0,0+1+2,1+0+0)*dipcross(2,2)
     &    - udertot(0+0+0,0+1+0,1+0+2)*dipcross(3,3) )
     &    + rlam*uuuder(2)*dotprod

         ! l = 3
         
         fld(3) = fld(3) + rmu*2.0d0*(
     &      udertot(0+0+2,0+0+0,1+1+0)*dipcross(1,1)
     &    + udertot(0+0+1,0+0+1,1+1+0)*(dipcross(1,2) + dipcross(2,1))
     &    + udertot(0+0+0,0+0+2,1+1+0)*dipcross(2,2)
     &    - udertot(0+0+0,0+0+0,1+1+2)*dipcross(3,3) )
     &    + rlam*uuuder(3)*dotprod

         endif


         endif
      endif
      return
      end
c
c
c
c
c
c***********************************************************************
      subroutine intker_mindlinc(icomp,rlame,source,ifsingle,
     1           sigma_sl,ifdouble,dipstr,dipvec,ztrg,pot,iffld,fld)
c***********************************************************************
c
c     This subroutine computes the displacement due to the 
c     harmonic Mindlin C image solution and all of its
c     derivatives, from which stress and strain can be extracted.
c     (See okada_imageanalysis for details.)
c     It doesn't include the 1/rmu scaling.
c
c     INPUT:
c
c     icomp      desired component of displacement vector
c     rlame(2)   Lame coefficients supplied in the form
c                rlame(1) = rlam, rlame(2) = rmu
c     source(3)  image source location.
c                This source must be the reflected image in the 
c                upper half-space, NOT the original source in the 
c                lower half-space.
c     ifsingle   (the single layer kernel flag)
c                0 means ignore the <<force>> vector 
c                1 means a <<force>> vector is supplied 
c     sigma_sl   force vector
c     ifdouble   (the double layer kernel flag)
c                0 means ignore the <<dislocation>> 
c                1 means a <<dislocation vector>> is supplied 
c                  with surface normal <<dipvec>>
c     dipstr     dislocation vector
c     dipvec     surface orientation
c                
c     ztrg       target location
c     pot        component <<icomp>> of displacement vector
c     fld        -gradient(pot) at ztrg
c-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: icomp
      real(8), intent(in) :: rlame(2)
      real(8), intent(in) :: source(3)
      integer, intent(in) :: ifsingle
      real(8), intent(in) :: sigma_sl(3)
      integer, intent(in) :: ifdouble
      real(8), intent(in) :: dipstr(3)
      real(8), intent(in) :: dipvec(3)
      real(8), intent(in) :: ztrg(3)
      complex(8), intent(out) :: pot
      integer, intent(in) :: iffld
      complex(8), intent(inout) :: fld(3)

!!      integer :: nder(3)
      integer :: i
      integer :: j
      integer :: k
!!      integer :: l
      integer :: n
!!      real(8) :: delta(3,3)
      real(8) :: uder(3,3,3)
      real(8) :: uuuder(3)
      real(8) :: v(3,3,3)
      real(8) :: vdertot(0:4,0:4,0:4)
      real(8) :: alpha
      real(8) :: d3
      real(8) :: d5
      real(8) :: d7
      real(8) :: d9
      real(8) :: dotprod
      real(8) :: dx
      real(8) :: dy
      real(8) :: dz
      real(8) :: rlam
      real(8) :: rmu
      real(8) :: rr1
      real(8) :: rr2
      real(8) :: rr3
      real(8) :: scale2c
      real(8) :: sgn
!!      real(8) :: uder1
!!      real(8) :: uder2
!!      real(8) :: uder4
      real(8) :: uuu
      real(8) :: v1
      real(8) :: v2
      real(8) :: v3
      real(8) :: v11
      real(8) :: v12
      real(8) :: v13
      real(8) :: v22
      real(8) :: v23
      real(8) :: v33
      real(8) :: scale2
      real(8) :: cd
      real(8) :: cd__2
      real(8) :: cd__3
      real(8) :: cd__5
      real(8) :: cd__7
      real(8) :: cd__9
      real(8) :: rr1__2
      real(8) :: rr1__3
      real(8) :: rr2__2
      real(8) :: rr2__3
      real(8) :: rr3__2
      real(8) :: rr3__3
      real(8) :: rr3__4
      real(8) :: dx__2
      real(8) :: dx__4
      real(8) :: dy__2
      real(8) :: dy__4
      real(8) :: dipcross(3,3)
      real(8) :: talsz
c      
      rlam = rlame(1)
      rmu = rlame(2)
      alpha = (rlam+rmu)/(rlam+2*rmu)
      scale2 = -(2.0d0-alpha)
ccc      write(6,*)' rlam =',rlam
ccc      write(6,*)' rmu =',rmu
ccc      write(6,*)' alpha =',alpha
!!        do i=1,3
!!        do j=1,3
!!        delta(i,j)=0
!!        enddo
!!        enddo
!!c
!!        do i=1,3
!!        delta(i,i)=1
!!        enddo
c
      dx = ztrg(1)-source(1)
      dy = ztrg(2)-source(2)
      dz = ztrg(3)-source(3)
      rr1 = dx
      rr2 = dy
      rr3 = -dz
!!      cd = sqrt(dx*dx+dy*dy+dz*dz)

      cd__2 = dx*dx+dy*dy+dz*dz
      cd = sqrt(cd__2)
      cd__3 = cd * cd__2
      cd__5 = cd__3 * cd__2
      cd__7 = cd__5 * cd__2
      cd__9 = cd__7 * cd__2

      rr1__2 = rr1 * rr1
      rr1__3 = rr1__2 * rr1

      rr2__2 = rr2 * rr2
      rr2__3 = rr2__2 * rr2

      rr3__2 = rr3 * rr3
      rr3__3 = rr3__2 * rr3

!!      if (iffld.eq.1) then

      rr3__4 = rr3__3 * rr3

      dx__2 = dx * dx
      dx__4 = dx__2 * dx__2

      dy__2 = dy * dy
      dy__4 = dy__2 * dy__2

!!      endif

c
      pot = 0.0d0
      if (iffld.eq.1) then
         fld(1) = 0.0d0
         fld(2) = 0.0d0
         fld(3) = 0.0d0
      endif
c
c
c     Compute three derivatives of harmonic function 
c     V = 1/R
c     with respect to R_1, R_2 and R_3. These are V1,V11,
c     V(1,1,1) etc. We will use these
c     derivatives to compute both the gradient of the SLP
c     and the DLP.
c
c     In this subroutine, it is convenient to first compute
c     uder(i,j,k) = -gradient with respect to R_k of
c     SLP kernel (assuming xi_3 is fixed). 
c     When used to evaluate gradient of SLP below, 
c     uder(i,j,k) needs sign flip for component k=3, since 
c
c     d/dx_1 = d/dR_1, d/dx_2 = d/dR_2, d/dx_3 = - d/dR_3
c
c     It is, however, correctly structured for use in the DLP
c     kernel (if present).
c
      d9 = 1.0d0/(cd__9)
      d7 = 1.0d0/(cd__7)
      d5 = 1.0d0/(cd__5)
      d3 = 1.0d0/(cd__3)
      v1 = -rr1*d3
      v2 = -rr2*d3
      v3 = -rr3*d3
      v11 = -d3 + 3*d5*rr1__2
      v12 = 3*d5*rr1*rr2
      v13 = 3*d5*rr1*rr3
      v22 = -d3 + 3*d5*rr2__2
      v23 = 3*d5*rr2*rr3
      v33 = -d3 + 3*d5*rr3__2
      if ( (ifdouble.eq.1) .or. (iffld.eq.1) ) then
         v(1,1,1) = 9*rr1*d5 - 15*d7*rr1__3
         v(1,1,2) = 3*rr2*d5 - 15*d7*rr2*rr1__2
         v(1,1,3) = 3*rr3*d5 - 15*d7*rr3*rr1__2
         v(1,2,2) = 3*rr1*d5 - 15*d7*rr1*rr2__2
         v(1,2,3) = - 15*d7*rr1*rr2*rr3
         v(1,3,3) = 3*rr1*d5 - 15*d7*rr1*rr3__2
         uder(1,1,1) = v(1,1,1)
         uder(1,1,2) = v(1,1,2)
         uder(1,2,1) = v(1,1,2)
         v(1,2,1) = v(1,1,2)
         uder(1,1,3) = v(1,1,3)
         uder(1,3,1) = v(1,1,3)
         v(1,3,1) = v(1,1,3)
         uder(1,3,3) = v(1,3,3)
         uder(1,2,3) = v(1,2,3)
         uder(1,3,2) = v(1,2,3)
         v(1,3,2) = v(1,2,3)
         uder(1,2,2) = v(1,2,2)
         if (icomp.ne.1) then
            uder(2,1,1) = v(1,1,2)
            v(2,1,1) = v(1,1,2)
            v(2,2,2)  = 9*rr2*d5 - 15*d7*rr2__3
            v(2,2,3)  = 3*rr3*d5 - 15*d7*rr3*rr2__2
            v(2,3,3)  = 3*rr2*d5 - 15*d7*rr2*rr3__2
            uder(2,2,2) = v(2,2,2)
            uder(2,1,2) = v(1,2,2)
            v(2,1,2) = v(1,2,2)
            uder(2,2,1) = v(1,2,2)
            v(2,2,1) = v(1,2,2)
            uder(2,2,3) = v(2,2,3)
            uder(2,3,2) = v(2,2,3)
            v(2,3,2) = v(2,2,3)
            uder(2,3,3) = v(2,3,3)
            uder(2,1,3) = v(1,2,3)
            v(2,1,3) = v(1,2,3)
            uder(2,3,1) = v(1,2,3)
            v(2,3,1) = v(1,2,3)
         endif
         if (icomp.eq. 3) then
            uder(3,1,1) = v(1,1,3)
            v(3,1,1) = v(1,1,3)
            uder(3,2,2) = v(2,2,3)
            v(3,2,2) = v(2,2,3)
            v(3,3,3) = 9*rr3*d5 - 15*d7*rr3__3
            uder(3,3,3) = v(3,3,3)
            uder(3,1,3) = v(1,3,3)
            v(3,1,3) = v(1,3,3)
            uder(3,3,1) = v(1,3,3)
            v(3,3,1) = v(1,3,3)
            uder(3,2,3) = v(2,3,3)
            v(3,2,3) = v(2,3,3)
            uder(3,3,2) = v(2,3,3)
            v(3,3,2) = v(2,3,3)
            uder(3,1,2) = v(1,2,3)
            v(3,1,2) = v(1,2,3)
            uder(3,2,1) = v(1,2,3)
            v(3,2,1) = v(1,2,3)
         endif
         do j = 1,3
            do k = 1,3
               uder(icomp,j,k) = -uder(icomp,j,k)*alpha*source(3)
               if (icomp.eq.3) uder(icomp,j,k) = -uder(icomp,j,k)
            enddo
         enddo
         if (icomp.eq.1) then
            uder(1,3,1)=uder(1,3,1)-scale2*v11
            uder(1,3,2)=uder(1,3,2)-scale2*v12
            uder(1,3,3)=uder(1,3,3)-scale2*v13
         endif
         if (icomp.eq.2) then
            uder(2,3,1)=uder(2,3,1)-scale2*v12
            uder(2,3,2)=uder(2,3,2)-scale2*v22
            uder(2,3,3)=uder(2,3,3)-scale2*v23
         endif
         if (icomp.eq.3) then
            uder(3,1,1)=uder(3,1,1)-scale2*v11
            uder(3,1,2)=uder(3,1,2)-scale2*v12
            uder(3,1,3)=uder(3,1,3)-scale2*v13
            uder(3,2,1)=uder(3,2,1)-scale2*v12
            uder(3,2,2)=uder(3,2,2)-scale2*v22
            uder(3,2,3)=uder(3,2,3)-scale2*v23
         endif
      endif
c
      if (ifsingle.eq.1) then
c
c     compute contribution to SLP for component icomp
c
         scale2c = -(2.0d0-alpha)*sigma_sl(3)
         if (icomp.eq.1) then
            pot = scale2c*v1 + source(3)* 
     1         alpha*(sigma_sl(1)*v11+sigma_sl(2)*v12+sigma_sl(3)*v13)
         else if (icomp.eq.2) then
            pot = scale2c*v2 + source(3)* 
     1         alpha*(sigma_sl(1)*v12+sigma_sl(2)*v22+sigma_sl(3)*v23)
         else 
            pot = scale2*(v1*sigma_sl(1)+v2*sigma_sl(2)) - source(3)* 
     1         alpha*(sigma_sl(1)*v13 + sigma_sl(2)*v23 + 
     1         sigma_sl(3)*v33)
         endif 
c
c
c     compute contribution to -gradient of SLP for component icomp
c
        if (iffld.eq.1) then
           do i=icomp,icomp
           do k=1,3
           sgn = 1.0d0
              if (k.eq.3) sgn = -1.0d0
              do j=1,3
                 fld(k) = fld(k) + 
     1               sgn*uder(i,j,k)*sigma_sl(j)
              enddo
           enddo
           enddo
         endif
      endif
ccc      write(6,*)' after slp'
ccc      write(6,*)(dreal(fld(iii)),iii=1,3)
c
c
      if (ifdouble.eq.1) then
c
c     Compute the double layer kernel.
c
c     Recall that above, we computed
c     uder(i,j,k) = -gradient with respect to R_k of
c     SLP kernel (assuming xi_3 is fixed). 
c     Since
c
c     d/dxi_1 = -d/dR_1, d/dxi_2 = -d/dR_2, d/dxi_3 = -d/dR_3
c
c     uder(i,j,k) can now be interpreted as DLP kernel.  
c
c     However, from C image formula, we clearly need
c     an additional contribution because xi_3
c     appears in definition of kernel and that was held fixed
c     above when we differentiated w.r.t. x_3.
c
         if (icomp.eq.1) then
            uder(1,1,3) = uder(1,1,3) - alpha*v11
            uder(1,2,3) = uder(1,2,3) - alpha*v12
            uder(1,3,3) = uder(1,3,3) - alpha*v13
         else if (icomp.eq.2) then
            uder(2,1,3) = uder(2,1,3) - alpha*v12
            uder(2,2,3) = uder(2,2,3) - alpha*v22
            uder(2,3,3) = uder(2,3,3) - alpha*v23
         else if (icomp.eq.3) then
            uder(3,1,3) = uder(3,1,3) + alpha*v13
            uder(3,2,3) = uder(3,2,3) + alpha*v23
            uder(3,3,3) = uder(3,3,3) + alpha*v33
         endif
c
c     Compute vdertot to hold fourth derivatives of V.
c
      if (iffld.eq.1) then
      vdertot(4,0,0) = 3*(8*dx__4-24*dx__2*dy__2-24*dx__2*rr3__2+
     #3*dy__4+6*dy__2*rr3__2+3*rr3__4)*d9
      vdertot(3,1,0) = 15*dx*dy*(4*dx__2-3*dy__2-3*rr3__2)*d9
      vdertot(2,2,0) = -3*(-27*dx__2*dy__2+4*dx__4+3*dx__2*rr3__2+
     #4*dy__4+3*dy__2*rr3__2-rr3__4)*d9
      vdertot(2,1,1) = 15*dy*rr3*(6*dx__2-dy__2-rr3__2)*d9
      vdertot(3,0,1) = 15*dx*rr3*(4*dx__2-3*dy__2-3*rr3__2)*d9
      vdertot(2,0,2) = -3*(-27*dx__2*rr3__2+4*dx__4+3*dx__2*dy__2+
     #3*dy__2*rr3__2+4*rr3__4-dy__4)*d9
      vdertot(1,3,0) = -15*dx*dy*(-4*dy__2+3*dx__2+3*rr3__2)*d9
      vdertot(1,2,1) = -15*rr3*dx*(-6*dy__2+dx__2+rr3__2)*d9
      vdertot(1,1,2) = -15*dx*dy*(-6*rr3__2+dx__2+dy__2)*d9
      vdertot(1,0,3) = -15*dx*rr3*(-4*rr3__2+3*dx__2+3*dy__2)*d9
      vdertot(0,4,0) = 3*(8*dy__4-24*dx__2*dy__2-24*dy__2*rr3__2+
     #3*dx__4+6*dx__2*rr3__2+3*rr3__4)*d9
      vdertot(0,3,1) = -15*dy*rr3*(-4*dy__2+3*dx__2+3*rr3__2)*d9
      vdertot(0,2,2) = 3*(27*dy__2*rr3__2-3*dx__2*dy__2-4*dy__4-
     #3*dx__2*rr3__2-4*rr3__4+dx__4)*d9
      vdertot(0,1,3) = -15*dy*rr3*(-4*rr3__2+3*dx__2+3*dy__2)*d9
      vdertot(0,0,4) = 3*(8*rr3__4-24*dx__2*rr3__2-24*dy__2*rr3__2+
     #3*dx__4+6*dx__2*dy__2+3*dy__4)*d9
      endif
c
c    compute uuu and uuuder 
c    uuu = sum_n du_{icomp}^n/d\xi_n
c    uuuder(l) = sum_n -d/dx_l du_{icomp}^n/d\xi_n
c
c
         uuu = 0
         do n=1,3
            uuu = uuu+uder(icomp,n,n)
         enddo
c
         if( iffld .eq. 1 ) then

!!         do l = 1,3
!!            uuuder(l) = 0
!!            do n=1,3
!!               nder(1) = 0
!!               nder(2) = 0
!!               nder(3) = 0
!!               nder(icomp) = nder(icomp)+1
!!               nder(l) = nder(l)+1
!!               nder(n) = nder(n)+2
!!c
!!c     nder counts number of derivs w.r.t R_1, R_2, R_3.
!!c
!!               uder4= -vdertot(nder(1),nder(2),nder(3))*alpha*source(3)
!!               if (icomp.eq.3) uder4 = -uder4
!!               uuuder(l) = uuuder(l)+uder4
!!               if (icomp.eq.3) uuuder(l) = uuuder(l) - scale2*v(n,n,l)
!!            enddo
!!            if (icomp.ne.3) uuuder(l) = uuuder(l) - scale2*v(icomp,3,l)
!!            if (icomp.eq.3) uuuder(l) = uuuder(l) + scale2*v(3,3,l)
!!            if (icomp.eq.3) uuuder(l) = uuuder(l) + alpha*v(3,3,l)
!!            if (icomp.ne.3) uuuder(l) = uuuder(l) - alpha*v(icomp,3,l)
!!            if (l.ne.3) uuuder(l) = -uuuder(l)
!!         enddo


!        In the following, where you see subscripts like 0+1+2,
!        the first number comes from icomp, the second from l,
!        and the third from n, where n gets a 2 and the others a 1.

         if (icomp.eq.1) then

         ! l = 1

         uuuder(1) = - (alpha*source(3)*(
     &       - vdertot(1+1+2,0+0+0,0+0+0)
     &       - vdertot(1+1+0,0+0+2,0+0+0)
     &       - vdertot(1+1+0,0+0+0,0+0+2)
     &     )
     &     - (scale2 + alpha)*v(1,3,1)
     &   )

         ! l = 2

         uuuder(2) = - (alpha*source(3)*(
     &       - vdertot(1+0+2,0+1+0,0+0+0)
     &       - vdertot(1+0+0,0+1+2,0+0+0)
     &       - vdertot(1+0+0,0+1+0,0+0+2)
     &     )
     &     - (scale2 + alpha)*v(1,3,2)
     &   )

         ! l = 3

         uuuder(3) =   (alpha*source(3)*(
     &       - vdertot(1+0+2,0+0+0,0+1+0)
     &       - vdertot(1+0+0,0+0+2,0+1+0)
     &       - vdertot(1+0+0,0+0+0,0+1+2)
     &     )
     &     - (scale2 + alpha)*v(1,3,3)
     &   )

         else if (icomp.eq.2) then

         ! l = 1

         uuuder(1) = - (alpha*source(3)*(
     &       - vdertot(0+1+2,1+0+0,0+0+0)
     &       - vdertot(0+1+0,1+0+2,0+0+0)
     &       - vdertot(0+1+0,1+0+0,0+0+2)
     &     )
     &     - (scale2 + alpha)*v(2,3,1)
     &   )

         ! l = 2

         uuuder(2) = - (alpha*source(3)*(
     &       - vdertot(0+0+2,1+1+0,0+0+0)
     &       - vdertot(0+0+0,1+1+2,0+0+0)
     &       - vdertot(0+0+0,1+1+0,0+0+2)
     &     )
     &     - (scale2 + alpha)*v(2,3,2)
     &   )

         ! l = 3

         uuuder(3) =   (alpha*source(3)*(
     &       - vdertot(0+0+2,1+0+0,0+1+0)
     &       - vdertot(0+0+0,1+0+2,0+1+0)
     &       - vdertot(0+0+0,1+0+0,0+1+2)
     &     )
     &     - (scale2 + alpha)*v(2,3,3)
     &   )

         else if (icomp.eq.3) then

         ! l = 1

         uuuder(1) = - (alpha*source(3)*(
     &         vdertot(0+1+2,0+0+0,1+0+0)
     &       + vdertot(0+1+0,0+0+2,1+0+0)
     &       + vdertot(0+1+0,0+0+0,1+0+2)
     &     )
     &     - scale2*(v(1,1,1) + v(2,2,1)) + alpha*v(3,3,1)
     &   )

         ! l = 2

         uuuder(2) = - (alpha*source(3)*(
     &         vdertot(0+0+2,0+1+0,1+0+0)
     &       + vdertot(0+0+0,0+1+2,1+0+0)
     &       + vdertot(0+0+0,0+1+0,1+0+2)
     &     )
     &     - scale2*(v(1,1,2) + v(2,2,2)) + alpha*v(3,3,2)
     &   )

         ! l = 3

         uuuder(3) =   (alpha*source(3)*(
     &         vdertot(0+0+2,0+0+0,1+1+0)
     &       + vdertot(0+0+0,0+0+2,1+1+0)
     &       + vdertot(0+0+0,0+0+0,1+1+2)
     &     )
     &     - scale2*(v(1,1,3) + v(2,2,3)) + alpha*v(3,3,3)
     &   )

         endif


         endif
c
c     Compute contribution to displacement component.
c
         do j=1,3
         do k=1,3
            dipcross(j,k) = dipstr(j)*dipvec(k)
         enddo
         enddo

         dotprod = dipcross(1,1) + 
     1             dipcross(2,2) + 
     1             dipcross(3,3) 
         do j=1,3
         do k=1,3
            pot = pot + 
     1      rmu*(uder(icomp,j,k)+uder(icomp,k,j))*dipcross(j,k)
         enddo
         enddo
         pot = pot + rlam*uuu*dotprod
        if (iffld.eq.1) then
c
c     Compte contributions to -gradient of displacement compoment.
c
!!         do l=1,3
!!         do j=1,3
!!         do k=1,3
!!            nder(1) = 0
!!            nder(2) = 0
!!            nder(3) = 0
!!            nder(icomp) = nder(icomp)+1
!!            nder(l) = nder(l)+1
!!            nder(j) = nder(j)+1
!!            nder(k) = nder(k)+1
!!c
!!c     nder counts number of derivs w.r.t R_1, R_2, R_3.
!!c     uder1 will hold d/dR_l du_{icomp}^j/d\xi_k
!!c     uder2 will hold d/dR_l du_{icomp}^k/d\xi_j
!!c
!!c     lots of if statements to account for the numerous
!!c     special case contributions to kernel.
!!c
!!            uder1 = -vdertot(nder(1),nder(2),nder(3))*alpha*source(3)
!!            uder2 = -vdertot(nder(1),nder(2),nder(3))*alpha*source(3)
!!            if (icomp.eq.3) uder1 = -uder1
!!            if (icomp.eq.3) uder2 = -uder2
!!            if ((icomp.ne.3).and.(j.eq.3)) uder1 = uder1 -
!!     1                         scale2*v(icomp,k,l)
!!            if ((icomp.ne.3).and.(k.eq.3)) uder2 = uder2 -
!!     1                         scale2*v(icomp,j,l)
!!            if ((icomp.eq.3).and.(j.ne.3)) 
!!     1                     uder1 = uder1 - scale2*v(j,k,l)
!!            if ((icomp.eq.3).and.(k.ne.3)) 
!!     1                     uder2 = uder2 - scale2*v(j,k,l)
!!c
!!            if ((icomp.ne.3).and.(k.eq.3)) uder1 = uder1 -
!!     1                         alpha*v(icomp,j,l)
!!            if ((icomp.ne.3).and.(j.eq.3)) uder2 = uder2 -
!!     1                         alpha*v(icomp,k,l)
!!            if ((icomp.eq.3).and.(k.eq.3)) uder1 = uder1 +
!!     1                         alpha*v(icomp,j,l)
!!            if ((icomp.eq.3).and.(j.eq.3)) uder2 = uder2 +
!!     1                         alpha*v(icomp,k,l)
!!c
!!c     want -d/dx_l, so flip sign for first and second
!!c     component. (d/dR_3 is already - d/dx_3.)
!!c
!!            if (l.ne.3) uder1 = -uder1
!!            if (l.ne.3) uder2 = -uder2
!!            fld(l) = fld(l) + 
!!     1               rmu*(uder1+uder2)*dipcross(j,k)
!!            if (j.eq.k) fld(l)=fld(l)+rlam*uuuder(l)*dipcross(j,k)
!!         enddo
!!         enddo
!!         enddo


         talsz = alpha*source(3)*2.0d0

!!         do l=1,3
!!         do j=1,3
!!         do k=1,3
!!            nder(1) = 0
!!            nder(2) = 0
!!            nder(3) = 0
!!            nder(icomp) = nder(icomp)+1
!!            nder(l) = nder(l)+1
!!            nder(j) = nder(j)+1
!!            nder(k) = nder(k)+1
!!
!!            uder1 = -vdertot(nder(1),nder(2),nder(3))*talsz
!!            if (icomp.eq.3) uder1 = -uder1
!!
!!            if ((icomp.ne.3).and.(j.eq.3)) uder1 =
!!     &        uder1 - (scale2 + alpha)*v(icomp,k,l)
!!            if ((icomp.ne.3).and.(k.eq.3)) uder1 =
!!     &        uder1 - (scale2 + alpha)*v(icomp,j,l)
!!
!!            if ((icomp.eq.3).and.(j.ne.3)) uder1 =
!!     &        uder1 - scale2*v(j,k,l)
!!            if ((icomp.eq.3).and.(j.eq.3)) uder1 =
!!     &        uder1 + alpha*v(icomp,k,l)
!!            if ((icomp.eq.3).and.(k.ne.3)) uder1 =
!!     &        uder1 - scale2*v(j,k,l)
!!            if ((icomp.eq.3).and.(k.eq.3)) uder1 =
!!     &        uder1 + alpha*v(icomp,j,l)
!!            
!!            if (l.ne.3) uder1 = -uder1
!!            fld(l) = fld(l) + rmu*uder1*dipcross(j,k)
!!         enddo
!!         enddo
!!         fld(l) = fld(l) + rlam*uuuder(l)*dotprod
!!         enddo


!        In the following, where you see subscripts like 0+0+1+0,
!        the first number comes from icomp, the second from l,
!        the third from j, and the fourth from k.

         if (icomp.eq.1) then

         ! l = 1

         fld(1) = fld(1) + rlam*uuuder(1)*dotprod - rmu*(
     &     dipcross(1,1)*(-vdertot(1+1+1+1,0+0+0+0,0+0+0+0)*talsz)
     &   + dipcross(1,2)*(-vdertot(1+1+1+0,0+0+0+1,0+0+0+0)*talsz)
     &   + dipcross(1,3)*(-vdertot(1+1+1+0,0+0+0+0,0+0+0+1)*talsz
     &                     - (scale2+alpha)*v(1,1,1))
     &   + dipcross(2,1)*(-vdertot(1+1+0+1,0+0+1+0,0+0+0+0)*talsz)
     &   + dipcross(2,2)*(-vdertot(1+1+0+0,0+0+1+1,0+0+0+0)*talsz)
     &   + dipcross(2,3)*(-vdertot(1+1+0+0,0+0+1+0,0+0+0+1)*talsz
     &                     - (scale2+alpha)*v(1,2,1))
     &   + dipcross(3,1)*(-vdertot(1+1+0+1,0+0+0+0,0+0+1+0)*talsz
     &                     - (scale2+alpha)*v(1,1,1))
     &   + dipcross(3,2)*(-vdertot(1+1+0+0,0+0+0+1,0+0+1+0)*talsz
     &                     - (scale2+alpha)*v(1,2,1))
     &   + dipcross(3,3)*(-vdertot(1+1+0+0,0+0+0+0,0+0+1+1)*talsz
     &                     - 2.0d0*(scale2+alpha)*v(1,3,1))
     &   )

         ! l = 2

         fld(2) = fld(2) + rlam*uuuder(2)*dotprod - rmu*(
     &     dipcross(1,1)*(-vdertot(1+0+1+1,0+1+0+0,0+0+0+0)*talsz)
     &   + dipcross(1,2)*(-vdertot(1+0+1+0,0+1+0+1,0+0+0+0)*talsz)
     &   + dipcross(1,3)*(-vdertot(1+0+1+0,0+1+0+0,0+0+0+1)*talsz
     &                     - (scale2+alpha)*v(1,1,2))
     &   + dipcross(2,1)*(-vdertot(1+0+0+1,0+1+1+0,0+0+0+0)*talsz)
     &   + dipcross(2,2)*(-vdertot(1+0+0+0,0+1+1+1,0+0+0+0)*talsz)
     &   + dipcross(2,3)*(-vdertot(1+0+0+0,0+1+1+0,0+0+0+1)*talsz
     &                     - (scale2+alpha)*v(1,2,2))
     &   + dipcross(3,1)*(-vdertot(1+0+0+1,0+1+0+0,0+0+1+0)*talsz
     &                     - (scale2+alpha)*v(1,1,2))
     &   + dipcross(3,2)*(-vdertot(1+0+0+0,0+1+0+1,0+0+1+0)*talsz
     &                     - (scale2+alpha)*v(1,2,2))
     &   + dipcross(3,3)*(-vdertot(1+0+0+0,0+1+0+0,0+0+1+1)*talsz
     &                     - 2.0d0*(scale2+alpha)*v(1,3,2))
     &   )

         ! l = 3

         fld(3) = fld(3) + rlam*uuuder(3)*dotprod + rmu*(
     &     dipcross(1,1)*(-vdertot(1+0+1+1,0+0+0+0,0+1+0+0)*talsz)
     &   + dipcross(1,2)*(-vdertot(1+0+1+0,0+0+0+1,0+1+0+0)*talsz)
     &   + dipcross(1,3)*(-vdertot(1+0+1+0,0+0+0+0,0+1+0+1)*talsz
     &                     - (scale2+alpha)*v(1,1,3))
     &   + dipcross(2,1)*(-vdertot(1+0+0+1,0+0+1+0,0+1+0+0)*talsz)
     &   + dipcross(2,2)*(-vdertot(1+0+0+0,0+0+1+1,0+1+0+0)*talsz)
     &   + dipcross(2,3)*(-vdertot(1+0+0+0,0+0+1+0,0+1+0+1)*talsz
     &                     - (scale2+alpha)*v(1,2,3))
     &   + dipcross(3,1)*(-vdertot(1+0+0+1,0+0+0+0,0+1+1+0)*talsz
     &                     - (scale2+alpha)*v(1,1,3))
     &   + dipcross(3,2)*(-vdertot(1+0+0+0,0+0+0+1,0+1+1+0)*talsz
     &                     - (scale2+alpha)*v(1,2,3))
     &   + dipcross(3,3)*(-vdertot(1+0+0+0,0+0+0+0,0+1+1+1)*talsz
     &                     - 2.0d0*(scale2+alpha)*v(1,3,3))
     &   )

         else if (icomp.eq.2) then

         ! l = 1

         fld(1) = fld(1) + rlam*uuuder(1)*dotprod - rmu*(
     &     dipcross(1,1)*(-vdertot(0+1+1+1,1+0+0+0,0+0+0+0)*talsz)
     &   + dipcross(1,2)*(-vdertot(0+1+1+0,1+0+0+1,0+0+0+0)*talsz)
     &   + dipcross(1,3)*(-vdertot(0+1+1+0,1+0+0+0,0+0+0+1)*talsz
     &                     - (scale2+alpha)*v(2,1,1))
     &   + dipcross(2,1)*(-vdertot(0+1+0+1,1+0+1+0,0+0+0+0)*talsz)
     &   + dipcross(2,2)*(-vdertot(0+1+0+0,1+0+1+1,0+0+0+0)*talsz)
     &   + dipcross(2,3)*(-vdertot(0+1+0+0,1+0+1+0,0+0+0+1)*talsz
     &                     - (scale2+alpha)*v(2,2,1))
     &   + dipcross(3,1)*(-vdertot(0+1+0+1,1+0+0+0,0+0+1+0)*talsz
     &                     - (scale2+alpha)*v(2,1,1))
     &   + dipcross(3,2)*(-vdertot(0+1+0+0,1+0+0+1,0+0+1+0)*talsz
     &                     - (scale2+alpha)*v(2,2,1))
     &   + dipcross(3,3)*(-vdertot(0+1+0+0,1+0+0+0,0+0+1+1)*talsz
     &                     - 2.0d0*(scale2+alpha)*v(2,3,1))
     &   )

         ! l = 2

         fld(2) = fld(2) + rlam*uuuder(2)*dotprod - rmu*(
     &     dipcross(1,1)*(-vdertot(0+0+1+1,1+1+0+0,0+0+0+0)*talsz)
     &   + dipcross(1,2)*(-vdertot(0+0+1+0,1+1+0+1,0+0+0+0)*talsz)
     &   + dipcross(1,3)*(-vdertot(0+0+1+0,1+1+0+0,0+0+0+1)*talsz
     &                     - (scale2+alpha)*v(2,1,2))
     &   + dipcross(2,1)*(-vdertot(0+0+0+1,1+1+1+0,0+0+0+0)*talsz)
     &   + dipcross(2,2)*(-vdertot(0+0+0+0,1+1+1+1,0+0+0+0)*talsz)
     &   + dipcross(2,3)*(-vdertot(0+0+0+0,1+1+1+0,0+0+0+1)*talsz
     &                     - (scale2+alpha)*v(2,2,2))
     &   + dipcross(3,1)*(-vdertot(0+0+0+1,1+1+0+0,0+0+1+0)*talsz
     &                     - (scale2+alpha)*v(2,1,2))
     &   + dipcross(3,2)*(-vdertot(0+0+0+0,1+1+0+1,0+0+1+0)*talsz
     &                     - (scale2+alpha)*v(2,2,2))
     &   + dipcross(3,3)*(-vdertot(0+0+0+0,1+1+0+0,0+0+1+1)*talsz
     &                     - 2.0d0*(scale2+alpha)*v(2,3,2))
     &   )

         ! l = 3

         fld(3) = fld(3) + rlam*uuuder(3)*dotprod + rmu*(
     &     dipcross(1,1)*(-vdertot(0+0+1+1,1+0+0+0,0+1+0+0)*talsz)
     &   + dipcross(1,2)*(-vdertot(0+0+1+0,1+0+0+1,0+1+0+0)*talsz)
     &   + dipcross(1,3)*(-vdertot(0+0+1+0,1+0+0+0,0+1+0+1)*talsz
     &                     - (scale2+alpha)*v(2,1,3))
     &   + dipcross(2,1)*(-vdertot(0+0+0+1,1+0+1+0,0+1+0+0)*talsz)
     &   + dipcross(2,2)*(-vdertot(0+0+0+0,1+0+1+1,0+1+0+0)*talsz)
     &   + dipcross(2,3)*(-vdertot(0+0+0+0,1+0+1+0,0+1+0+1)*talsz
     &                     - (scale2+alpha)*v(2,2,3))
     &   + dipcross(3,1)*(-vdertot(0+0+0+1,1+0+0+0,0+1+1+0)*talsz
     &                     - (scale2+alpha)*v(2,1,3))
     &   + dipcross(3,2)*(-vdertot(0+0+0+0,1+0+0+1,0+1+1+0)*talsz
     &                     - (scale2+alpha)*v(2,2,3))
     &   + dipcross(3,3)*(-vdertot(0+0+0+0,1+0+0+0,0+1+1+1)*talsz
     &                     - 2.0d0*(scale2+alpha)*v(2,3,3))
     &   )

         else if (icomp.eq.3) then

         ! l = 1

         fld(1) = fld(1) + rlam*uuuder(1)*dotprod - rmu*(
     &     dipcross(1,1)*( vdertot(0+1+1+1,0+0+0+0,1+0+0+0)*talsz
     &                     - 2.0d0*scale2*v(1,1,1))
     &   + dipcross(1,2)*( vdertot(0+1+1+0,0+0+0+1,1+0+0+0)*talsz
     &                     - 2.0d0*scale2*v(1,2,1))
     &   + dipcross(1,3)*( vdertot(0+1+1+0,0+0+0+0,1+0+0+1)*talsz
     &                     - scale2*v(1,3,1) + alpha*v(3,1,1))
     &   + dipcross(2,1)*( vdertot(0+1+0+1,0+0+1+0,1+0+0+0)*talsz
     &                     - 2.0d0*scale2*v(2,1,1))
     &   + dipcross(2,2)*( vdertot(0+1+0+0,0+0+1+1,1+0+0+0)*talsz
     &                     - 2.0d0*scale2*v(2,2,1))
     &   + dipcross(2,3)*( vdertot(0+1+0+0,0+0+1+0,1+0+0+1)*talsz
     &                     - scale2*v(2,3,1) + alpha*v(3,2,1))
     &   + dipcross(3,1)*( vdertot(0+1+0+1,0+0+0+0,1+0+1+0)*talsz
     &                     + (alpha - scale2)*v(3,1,1))
     &   + dipcross(3,2)*( vdertot(0+1+0+0,0+0+0+1,1+0+1+0)*talsz
     &                     + (alpha - scale2)*v(3,2,1))
     &   + dipcross(3,3)*( vdertot(0+1+0+0,0+0+0+0,1+0+1+1)*talsz
     &                     + 2.0d0*alpha*v(3,3,1))
     &   )

         ! l = 2

         fld(2) = fld(2) + rlam*uuuder(2)*dotprod - rmu*(
     &     dipcross(1,1)*( vdertot(0+0+1+1,0+1+0+0,1+0+0+0)*talsz
     &                     - 2.0d0*scale2*v(1,1,2))
     &   + dipcross(1,2)*( vdertot(0+0+1+0,0+1+0+1,1+0+0+0)*talsz
     &                     - 2.0d0*scale2*v(1,2,2))
     &   + dipcross(1,3)*( vdertot(0+0+1+0,0+1+0+0,1+0+0+1)*talsz
     &                     - scale2*v(1,3,2) + alpha*v(3,1,2))
     &   + dipcross(2,1)*( vdertot(0+0+0+1,0+1+1+0,1+0+0+0)*talsz
     &                     - 2.0d0*scale2*v(2,1,2))
     &   + dipcross(2,2)*( vdertot(0+0+0+0,0+1+1+1,1+0+0+0)*talsz
     &                     - 2.0d0*scale2*v(2,2,2))
     &   + dipcross(2,3)*( vdertot(0+0+0+0,0+1+1+0,1+0+0+1)*talsz
     &                     - scale2*v(2,3,2) + alpha*v(3,2,2))
     &   + dipcross(3,1)*( vdertot(0+0+0+1,0+1+0+0,1+0+1+0)*talsz
     &                     + (alpha - scale2)*v(3,1,2))
     &   + dipcross(3,2)*( vdertot(0+0+0+0,0+1+0+1,1+0+1+0)*talsz
     &                     + (alpha - scale2)*v(3,2,2))
     &   + dipcross(3,3)*( vdertot(0+0+0+0,0+1+0+0,1+0+1+1)*talsz
     &                     + 2.0d0*alpha*v(3,3,2))
     &   )

         ! l = 3

         fld(3) = fld(3) + rlam*uuuder(3)*dotprod + rmu*(
     &     dipcross(1,1)*( vdertot(0+0+1+1,0+0+0+0,1+1+0+0)*talsz
     &                     - 2.0d0*scale2*v(1,1,3))
     &   + dipcross(1,2)*( vdertot(0+0+1+0,0+0+0+1,1+1+0+0)*talsz
     &                     - 2.0d0*scale2*v(1,2,3))
     &   + dipcross(1,3)*( vdertot(0+0+1+0,0+0+0+0,1+1+0+1)*talsz
     &                     - scale2*v(1,3,3) + alpha*v(3,1,3))
     &   + dipcross(2,1)*( vdertot(0+0+0+1,0+0+1+0,1+1+0+0)*talsz
     &                     - 2.0d0*scale2*v(2,1,3))
     &   + dipcross(2,2)*( vdertot(0+0+0+0,0+0+1+1,1+1+0+0)*talsz
     &                     - 2.0d0*scale2*v(2,2,3))
     &   + dipcross(2,3)*( vdertot(0+0+0+0,0+0+1+0,1+1+0+1)*talsz
     &                     - scale2*v(2,3,3) + alpha*v(3,2,3))
     &   + dipcross(3,1)*( vdertot(0+0+0+1,0+0+0+0,1+1+1+0)*talsz
     &                     + (alpha - scale2)*v(3,1,3))
     &   + dipcross(3,2)*( vdertot(0+0+0+0,0+0+0+1,1+1+1+0)*talsz
     &                     + (alpha - scale2)*v(3,2,3))
     &   + dipcross(3,3)*( vdertot(0+0+0+0,0+0+0+0,1+1+1+1)*talsz
     &                     + 2.0d0*alpha*v(3,3,3))
     &   )

         endif


        endif
      endif
ccc      write(6,*)' after dlp'
ccc      write(6,*)(dreal(fld(iii)),iii=1,3)
      return
      end

