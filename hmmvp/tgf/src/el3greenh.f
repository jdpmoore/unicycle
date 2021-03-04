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
c       Direct calculation of half-space elastostatic 
c       Green's functions in R^3 (the Mindlin solution),
c       that satisfy zero normal stress.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine green3elth_eval(rlam,rmu,source,du,rnorm,
     $     target,fvec,ifstrain,strain)
c
c       Elastostatic double layer Green's function 
c
c       This subroutine evaluates the displacement vector fvec and the
c       strain tensor strain at the location xyz due to the double force
c       (aka slip, displacement jump) du at the origin with the
c       orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while displacement can be 
c       in any direction.
c       
c       ... elastostatic Green's function for the double force source
c
c       Okada 1992
c
c
c       INPUT:
c
c       rlam, rmu       Lame parameters
c       source(3)       Source location in lower half-space (z<0).
c       du(3)           Strength of double force source
c       rnorm(3)        Orientation vector of the double force source
c       target(3)       Target point in lower half-space (z<0).
c       ifstrain        Strain flag (1 means compute, 0 means do not).
c
c       OUTPUT:
c
c       fvec(3)         Displacement at the target location
c       strain(3,3)     Strain at the target
c
c
      implicit none

      real(8), intent(in) :: rlam
      real(8), intent(in) :: rmu
      real(8), intent(in) :: source(3)
      real(8), intent(in) :: du(3)
      real(8), intent(in) :: rnorm(3)
      real(8), intent(in) :: target(3)
      real(8), intent(out) :: fvec(3)
      integer, intent(in) :: ifstrain
      real(8), intent(inout) :: strain(3,3)

      real(8) :: fvecloc(3)
      real(8) :: strainloc(3,3)
      
      integer :: i
      integer :: j

c
c***********************************************************************
c     PART 1 === direct arrival
c***********************************************************************
c
      call green3elt_eval(rlam,rmu,source,du,rnorm,target,fvec,
     1     ifstrain,strain) 
ccc      if (2.ne.3) return
c      fvec(1)=0
c      fvec(2)=0
c      fvec(3)=0
c      if (ifstrain.eq.1) then
c         do i=1,3
c         do j=1,3
c            strain(i,j)=0
c         enddo
c         enddo
c      endif  
c
c***********************************************************************
c     PARTS 2, 3, and 4 === all image contributions
c     (in Okada 1992).
c***********************************************************************
c
      call green3elth_image_eval(rlam,rmu,source,du,rnorm,
     $     target,fvecloc,ifstrain,strainloc)

      fvec(1)=fvec(1)+fvecloc(1)
      fvec(2)=fvec(2)+fvecloc(2)
      fvec(3)=fvec(3)+fvecloc(3)
      if (ifstrain.eq.1) then
         do i=1,3
         do j=1,3
            strain(i,j)=strain(i,j)+strainloc(i,j)
         enddo
         enddo
      endif  

      return
      end
c
c
c
c
        subroutine green3elth_image_eval(rlam,rmu,source,du,rnorm,
     $     target,fvec,ifstrain,strain)
c
c       Elastostatic double layer Green's function 
c
c       This subroutine evaluates the displacement vector fvec and the
c       strain tensor strain at the location xyz due to the double force
c       (aka slip, displacement jump) du at the origin with the
c       orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while displacement can be 
c       in any direction.
c       
c       ... elastostatic Green's function for the double force source
c
c       Okada 1992
c
c
c       INPUT:
c
c       rlam, rmu       Lame parameters
c       source(3)       Source location in lower half-space (z<0).
c       du(3)           Strength of double force source
c       rnorm(3)        Orientation vector of the double force source
c       target(3)       Target point in lower half-space (z<0).
c       ifstrain        Strain flag (1 means compute, 0 means do not).
c
c       OUTPUT:
c
c       fvec(3)         Displacement at the target location
c       strain(3,3)     Strain at the target
c
c
      implicit none

      real(8), intent(in) :: rlam
      real(8), intent(in) :: rmu
      real(8), intent(in) :: source(3)
      real(8), intent(in) :: du(3)
      real(8), intent(in) :: rnorm(3)
      real(8), intent(in) :: target(3)
      real(8), intent(out) :: fvec(3)
      integer, intent(in) :: ifstrain
      real(8), intent(inout) :: strain(3,3)

      real(8) :: xyz(3)
      real(8) :: dudx(3,3)
      real(8) :: fvecloc(3)
      real(8) :: sourceim(3)
      real(8) :: targim(3)
      real(8) :: rlame(2)
      real(8) :: strainloc(3,3)
      real(8) :: df(3)

      real(8) :: alpha
      real(8) :: alphaim
      real(8) :: rlamim
      real(8) :: rmuim

      integer :: i
      integer :: j
      integer :: ifsingle
      integer :: ifdouble

      complex(8) :: pot
      complex(8) :: fld(3)
c
c
c***********************************************************************
c     PART 2 === elastic <<anti-image>>   [see notes]
c     merges A image with first part of Stresslet-like B image 
c     (in Okada 1992).
c***********************************************************************
c                        
      targim(1) = target(1)
      targim(2) = target(2)
      targim(3) = -target(3)
      alpha = (rlam+rmu)/(rlam+2*rmu)
c
      rlamim = rlam+4*rmu
      rmuim = -rmu
      alphaim = 2-alpha
ccc      if (2.ne.3) return
c
      xyz(1)=target(1)-source(1)
      xyz(2)=target(2)-source(2)
      xyz(3)=-target(3)-source(3)
c
      call green3elt_strain_image(rlam,rmu,alpha,source,du,rnorm,
     1     targim,fvecloc,ifstrain,strainloc) 
c
      fvec(1)=fvecloc(1)
      fvec(2)=fvecloc(2)
      fvec(3)=fvecloc(3)
c
      if (ifstrain.eq.1) then
         do i=1,3
         do j=1,3
            strain(i,j)=strainloc(i,j)
         enddo
         enddo
      endif  
ccc      if (2.ne.3) return
c
c***********************************************************************
c       === PART 3 ===
c       Harmonic contributions from image source.
c***********************************************************************
c
      rlame(1) = rlam
      rlame(2) = rmu
      sourceim(1) = source(1)
      sourceim(2) = source(2)
      sourceim(3) = -source(3)
c
      ifsingle = 0
      ifdouble = 1
      do i=1,3
         call intker_mindlinb(i,rlame,sourceim,ifsingle,
     1           df,ifdouble,du,rnorm,target,pot,ifstrain,fld)
         fvec(i) = fvec(i)+dreal(pot)/rmu
         if (ifstrain.eq.1) then
            dudx(i,1) = -dreal(fld(1))/rmu
            dudx(i,2) = -dreal(fld(2))/rmu
            dudx(i,3) = -dreal(fld(3))/rmu
         endif
      enddo
      if (ifstrain.eq.1) then
         do i = 1,3
         do j = 1,3
            strainloc(i,j) = (dudx(i,j)+dudx(j,i))/2
         enddo
         enddo
c
         do i=1,3
         do j=1,3
            strain(i,j)=strain(i,j)+strainloc(i,j)
         enddo
         enddo
      endif
ccc      if (2.ne.3) return
      do i=1,3
         call intker_mindlinc(i,rlame,sourceim,ifsingle,
     1           df,ifdouble,du,rnorm,target,pot,ifstrain,fld)
         fvec(i) = fvec(i)+target(3)*dreal(pot)/rmu
         if (ifstrain.eq.1) then
            dudx(i,1) = -target(3)*dreal(fld(1))/rmu
            dudx(i,2) = -target(3)*dreal(fld(2))/rmu
            dudx(i,3) = -target(3)*dreal(fld(3))/rmu
     1               +dreal(pot)/rmu
         endif
      enddo
c
      if (ifstrain.eq.1) then
         do i = 1,3
         do j = 1,3
            strainloc(i,j) = (dudx(i,j)+dudx(j,i))/2
         enddo
         enddo
c
         do i=1,3
         do j=1,3
            strain(i,j)=strain(i,j)+strainloc(i,j)
         enddo
         enddo
      endif
      return
      end
c
c
c
c
        subroutine green3elth_image_a_eval(rlam,rmu,source,du,rnorm,
     $     target,fvec,ifstrain,strain)
c
c       Elastostatic double layer Green's function 
c
c       This subroutine evaluates the displacement vector fvec and the
c       strain tensor strain at the location xyz due to the double force
c       (aka slip, displacement jump) du at the origin with the
c       orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while displacement can be 
c       in any direction.
c       
c       ... elastostatic Green's function for the double force source
c
c       Okada 1992
c
c
c       INPUT:
c
c       rlam, rmu       Lame parameters
c       source(3)       Source location in lower half-space (z<0).
c       du(3)           Strength of double force source
c       rnorm(3)        Orientation vector of the double force source
c       target(3)       Target point in lower half-space (z<0).
c       ifstrain        Strain flag (1 means compute, 0 means do not).
c
c       OUTPUT:
c
c       fvec(3)         Displacement at the target location
c       strain(3,3)     Strain at the target
c
c
      implicit none

      real(8), intent(in) :: rlam
      real(8), intent(in) :: rmu
      real(8), intent(in) :: source(3)
      real(8), intent(in) :: du(3)
      real(8), intent(in) :: rnorm(3)
      real(8), intent(in) :: target(3)
      real(8), intent(out) :: fvec(3)
      integer, intent(in) :: ifstrain
      real(8), intent(inout) :: strain(3,3)

      real(8) :: xyz(3)
      real(8) :: dudx(3,3)
      real(8) :: fvecloc(3)
      real(8) :: sourceim(3)
      real(8) :: targim(3)
      real(8) :: rlame(2)
      real(8) :: strainloc(3,3)
      real(8) :: df(3)

      real(8) :: alpha
      real(8) :: alphaim
      real(8) :: rlamim
      real(8) :: rmuim

      integer :: i
      integer :: j
      integer :: ifsingle
      integer :: ifdouble

      complex(8) :: pot
      complex(8) :: fld(3)
c
c
c***********************************************************************
c     PART 2 === elastic <<anti-image>>   [see notes]
c     merges A image with first part of Stresslet-like B image 
c     (in Okada 1992).
c***********************************************************************
c                        
      targim(1) = target(1)
      targim(2) = target(2)
      targim(3) = -target(3)
      alpha = (rlam+rmu)/(rlam+2*rmu)
c
      rlamim = rlam+4*rmu
      rmuim = -rmu
      alphaim = 2-alpha
ccc      if (2.ne.3) return
c
      xyz(1)=target(1)-source(1)
      xyz(2)=target(2)-source(2)
      xyz(3)=-target(3)-source(3)
c
      call green3elt_strain_image(rlam,rmu,alpha,source,du,rnorm,
     1     targim,fvecloc,ifstrain,strainloc) 
c
      fvec(1)=fvecloc(1)
      fvec(2)=fvecloc(2)
      fvec(3)=fvecloc(3)
c
      if (ifstrain.eq.1) then
         do i=1,3
         do j=1,3
            strain(i,j)=strainloc(i,j)
         enddo
         enddo
      endif  
      if (2.ne.3) return
c
c***********************************************************************
c       === PART 3 ===
c       Harmonic contributions from image source.
c***********************************************************************
c
      rlame(1) = rlam
      rlame(2) = rmu
      sourceim(1) = source(1)
      sourceim(2) = source(2)
      sourceim(3) = -source(3)
c
      ifsingle = 0
      ifdouble = 1
      do i=1,3
         call intker_mindlinb(i,rlame,sourceim,ifsingle,
     1           df,ifdouble,du,rnorm,target,pot,ifstrain,fld)
         fvec(i) = fvec(i)+dreal(pot)/rmu
         if (ifstrain.eq.1) then
            dudx(i,1) = -dreal(fld(1))/rmu
            dudx(i,2) = -dreal(fld(2))/rmu
            dudx(i,3) = -dreal(fld(3))/rmu
         endif
      enddo
      if (ifstrain.eq.1) then
         do i = 1,3
         do j = 1,3
            strainloc(i,j) = (dudx(i,j)+dudx(j,i))/2
         enddo
         enddo
c
         do i=1,3
         do j=1,3
            strain(i,j)=strain(i,j)+strainloc(i,j)
         enddo
         enddo
      endif
ccc      if (2.ne.3) return
      do i=1,3
         call intker_mindlinc(i,rlame,sourceim,ifsingle,
     1           df,ifdouble,du,rnorm,target,pot,ifstrain,fld)
         fvec(i) = fvec(i)+target(3)*dreal(pot)/rmu
         if (ifstrain.eq.1) then
            dudx(i,1) = -target(3)*dreal(fld(1))/rmu
            dudx(i,2) = -target(3)*dreal(fld(2))/rmu
            dudx(i,3) = -target(3)*dreal(fld(3))/rmu
     1               +dreal(pot)/rmu
         endif
      enddo
c
      if (ifstrain.eq.1) then
         do i = 1,3
         do j = 1,3
            strainloc(i,j) = (dudx(i,j)+dudx(j,i))/2
         enddo
         enddo
c
         do i=1,3
         do j=1,3
            strain(i,j)=strain(i,j)+strainloc(i,j)
         enddo
         enddo
      endif
      return
      end
c
c
c
c
        subroutine green3elth_mindlin_eval(rlam,rmu,source,du,rnorm,
     $     target,fvec,ifstrain,strain)
c
c       Elastostatic double layer Green's function 
c
c       This subroutine evaluates the displacement vector fvec and the
c       strain tensor strain at the location xyz due to the double force
c       (aka slip, displacement jump) du at the origin with the
c       orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while displacement can be 
c       in any direction.
c       
c       ... elastostatic Green's function for the double force source
c
c       Okada 1992
c
c       Mindlin image part only, parts B and C
c
c       INPUT:
c
c       rlam, rmu       Lame parameters
c       source(3)       Source location in lower half-space (z<0).
c       du(3)           Strength of double force source
c       rnorm(3)        Orientation vector of the double force source
c       target(3)       Target point in lower half-space (z<0).
c       ifstrain        Strain flag (1 means compute, 0 means do not).
c
c       OUTPUT:
c
c       fvec(3)         Displacement at the target location
c       strain(3,3)     Strain at the target
c
c
      implicit none

      real(8), intent(in) :: rlam
      real(8), intent(in) :: rmu
      real(8), intent(in) :: source(3)
      real(8), intent(in) :: du(3)
      real(8), intent(in) :: rnorm(3)
      real(8), intent(in) :: target(3)
      real(8), intent(out) :: fvec(3)
      integer, intent(in) :: ifstrain
      real(8), intent(inout) :: strain(3,3)

      real(8) :: xyz(3)
      real(8) :: dudx(3,3)
      real(8) :: sourceim(3)
      real(8) :: targim(3)
      real(8) :: rlame(2)
      real(8) :: strainloc(3,3)
      real(8) :: df(3)

      real(8) :: alpha
      real(8) :: rlamim
      real(8) :: rmuim
      real(8) :: alphaim

      integer :: i
      integer :: j
      integer :: ifsingle
      integer :: ifdouble

      complex(8) :: pot
      complex(8) :: fld(3)
c
c
c***********************************************************************
c     PART 2 === elastic <<anti-image>>   [see notes]
c     merges A image with first part of Stresslet-like B image 
c     (in Okada 1992).
c***********************************************************************
c                        
      targim(1) = target(1)
      targim(2) = target(2)
      targim(3) = -target(3)
      alpha = (rlam+rmu)/(rlam+2*rmu)
c
      rlamim = rlam+4*rmu
      rmuim = -rmu
      alphaim = 2-alpha
ccc      if (2.ne.3) return
c
      xyz(1)=target(1)-source(1)
      xyz(2)=target(2)-source(2)
      xyz(3)=-target(3)-source(3)
c
ccc      call green3elt_strain_image(rlam,rmu,alpha,source,du,rnorm,
ccc     1     targim,fvecloc,ifstrain,strainloc) 
c
      fvec(1)=0
      fvec(2)=0
      fvec(3)=0
c
      if (ifstrain.eq.1) then
         do i=1,3
         do j=1,3
            strain(i,j)=0
         enddo
         enddo
      endif  
ccc      if (2.ne.3) return
c
c***********************************************************************
c       === PART 3 ===
c       Harmonic contributions from image source.
c***********************************************************************
c
      rlame(1) = rlam
      rlame(2) = rmu
      sourceim(1) = source(1)
      sourceim(2) = source(2)
      sourceim(3) = -source(3)
c
      ifsingle = 0
      ifdouble = 1
      do i=1,3
         call intker_mindlinb(i,rlame,sourceim,ifsingle,
     1           df,ifdouble,du,rnorm,target,pot,ifstrain,fld)
         fvec(i) = fvec(i)+dreal(pot)/rmu
         if (ifstrain.eq.1) then
            dudx(i,1) = -dreal(fld(1))/rmu
            dudx(i,2) = -dreal(fld(2))/rmu
            dudx(i,3) = -dreal(fld(3))/rmu
         endif
      enddo
      if (ifstrain.eq.1) then
         do i = 1,3
         do j = 1,3
            strainloc(i,j) = (dudx(i,j)+dudx(j,i))/2
         enddo
         enddo
c
         do i=1,3
         do j=1,3
            strain(i,j)=strain(i,j)+strainloc(i,j)
         enddo
         enddo
      endif
ccc      if (2.ne.3) return
      do i=1,3
         call intker_mindlinc(i,rlame,sourceim,ifsingle,
     1           df,ifdouble,du,rnorm,target,pot,ifstrain,fld)
         fvec(i) = fvec(i)+target(3)*dreal(pot)/rmu
         if (ifstrain.eq.1) then
            dudx(i,1) = -target(3)*dreal(fld(1))/rmu
            dudx(i,2) = -target(3)*dreal(fld(2))/rmu
            dudx(i,3) = -target(3)*dreal(fld(3))/rmu
     1               +dreal(pot)/rmu
         endif
      enddo
c
      if (ifstrain.eq.1) then
         do i = 1,3
         do j = 1,3
            strainloc(i,j) = (dudx(i,j)+dudx(j,i))/2
         enddo
         enddo
c
         do i=1,3
         do j=1,3
            strain(i,j)=strain(i,j)+strainloc(i,j)
         enddo
         enddo
      endif
      return
      end
c
c
c
c
c***********************************************************************
        subroutine green3eluh_eval(rlam,rmu,source,df,target,
     $     fvec,ifstrain,strain)
c***********************************************************************
c
c       This subroutine evaluates the displacement vector <fvec> and 
c       strain tensor <strain> at the
c       location <target> due to the static force vector <df> at 
c       location <source>: the Mindlin solution
c
c       Mindlin's solution corresponds to a half-space solution
c       (z .leq. 0), which satisfies the condition that the normal 
c       component of stress vanishes at z=0.
c
c
c       INPUT:
c
c       rlam, rmu       Lame parameters
c       source(3)       Source location in lower half-space (z<0).
c       df(3)           Strength of single force 
c       target(3)       Target point in lower half-space (z<0).
c       ifstrain        Strain flag (1 means compute, 0 means do not).
c
c
c       OUTPUT:
c
c       fvec (real(8)) - the displacement at the target
c       strain (real(8)) - the strain at the target
c
c
      implicit none

      real(8), intent(in) :: rlam
      real(8), intent(in) :: rmu
      real(8), intent(in) :: source(3)
      real(8), intent(in) :: df(3)
      real(8), intent(in) :: target(3)
      real(8), intent(out) :: fvec(3)
      integer, intent(in) :: ifstrain
      real(8), intent(inout) :: strain(3,3)

      real(8) :: fvecloc(3)
      real(8) :: strainloc(3,3)

      integer :: i
      integer :: j

c
c***********************************************************************
c     PART 1 === direct arrival
c***********************************************************************
c
      call green3elu_eval(rlam,rmu,source,df,target,fvec,
     1     ifstrain,strain) 
ccc      if (2.ne.3) return
c
c***********************************************************************
c     PARTS 2, 3, and 4 === all image contributions
c     (in Okada 1992).
c***********************************************************************
c
      call green3eluh_image_eval(rlam,rmu,source,df,target,
     $     fvecloc,ifstrain,strainloc)

      fvec(1)=fvec(1)+fvecloc(1)
      fvec(2)=fvec(2)+fvecloc(2)
      fvec(3)=fvec(3)+fvecloc(3)
      if (ifstrain.eq.1) then
         do i=1,3
         do j=1,3
            strain(i,j)=strain(i,j)+strainloc(i,j)
         enddo
         enddo
      endif  

c
      return
      end
c
c
c
c
c
c***********************************************************************
        subroutine green3eluh_image_eval(rlam,rmu,source,df,target,
     $     fvec,ifstrain,strain)
c***********************************************************************
c
c       This subroutine evaluates the displacement vector <fvec> and 
c       strain tensor <strain> at the
c       location <target> due to the static force vector <df> at 
c       location <source>: the Mindlin solution
c
c       Mindlin's solution corresponds to a half-space solution
c       (z .leq. 0), which satisfies the condition that the normal 
c       component of stress vanishes at z=0.
c
c
c       INPUT:
c
c       rlam, rmu       Lame parameters
c       source(3)       Source location in lower half-space (z<0).
c       df(3)           Strength of single force 
c       target(3)       Target point in lower half-space (z<0).
c       ifstrain        Strain flag (1 means compute, 0 means do not).
c
c
c       OUTPUT:
c
c       fvec (real(8)) - the displacement at the target
c       strain (real(8)) - the strain at the target
c
c
      implicit none

      real(8), intent(in) :: rlam
      real(8), intent(in) :: rmu
      real(8), intent(in) :: source(3)
      real(8), intent(in) :: df(3)
      real(8), intent(in) :: target(3)
      real(8), intent(out) :: fvec(3)
      integer, intent(in) :: ifstrain
      real(8), intent(inout) :: strain(3,3)

      real(8) :: fvecloc(3)
      real(8) :: sourceim(3)
      real(8) :: targim(3)
      real(8) :: rlame(2)
      real(8) :: dudx(3,3)
      real(8) :: strainloc(3,3)
      real(8) :: dipstr(3)
      real(8) :: dipvec(3)

      real(8) :: rlamim
      real(8) :: rmuim

      integer :: i
      integer :: j
      integer :: ifsingle
      integer :: ifdouble

      complex(8) :: pot
      complex(8) :: fld(3)
c
c
c***********************************************************************
c     PART 2 === elastic <<anti-image>>   [see notes]
c     merges A image with first part of Stokeslet-like B image 
c     (in Okada 1992).
c***********************************************************************
c                        
      targim(1) = target(1)
      targim(2) = target(2)
      targim(3) = -target(3)
      rlamim = rlam+4*rmu
      rmuim = -rmu
c
      call green3elu_strain_image(rlamim,rmuim,source,df,targim,
     1     fvecloc,ifstrain,strainloc) 
c
c     inside the preceding subroutine,
c     scaling is by 1/(2*rmuim) and it should be by 1/(2*rmu).
c     This is fixed by flipping the sign (subtracting) in contributions
c     to displacement and strain here.
c
      fvec(1)=-fvecloc(1)
      fvec(2)=-fvecloc(2)
      fvec(3)=-fvecloc(3)
      if (ifstrain.eq.1) then
         do i=1,3
         do j=1,3
            strain(i,j)=-strainloc(i,j)
         enddo
         enddo
      endif  
cc      if (2.ne.3) return
c
c***********************************************************************
c       === PART 3 ===
c       Harmonic contributions from image source.
c***********************************************************************
c
c
      rlame(1) = rlam
      rlame(2) = rmu
      sourceim(1) = source(1)
      sourceim(2) = source(2)
      sourceim(3) = -source(3)
c
      ifsingle = 1
      ifdouble = 0
      do i=1,3
         call intker_mindlinb(i,rlame,sourceim,ifsingle,
     1           df,ifdouble,dipstr,dipvec,target,pot,ifstrain,fld)
         fvec(i) = fvec(i)+dreal(pot)/rmu
         if (ifstrain.eq.1) then
            dudx(i,1) = -dreal(fld(1))/rmu
            dudx(i,2) = -dreal(fld(2))/rmu
            dudx(i,3) = -dreal(fld(3))/rmu
         endif
      enddo
      if (ifstrain.eq.1) then
         do i = 1,3
         do j = 1,3
            strainloc(i,j) = (dudx(i,j)+dudx(j,i))/2
         enddo
         enddo
c
         do i=1,3
         do j=1,3
            strain(i,j)=strain(i,j)+strainloc(i,j)
         enddo
         enddo
      endif
      do i=1,3
         call intker_mindlinc(i,rlame,sourceim,ifsingle,
     1           df,ifdouble,dipstr,dipvec,target,pot,ifstrain,fld)
         if (ifstrain.eq.1) then
            dudx(i,1) = dreal(fld(1))/rmu
            dudx(i,2) = dreal(fld(2))/rmu
            dudx(i,3) = dreal(fld(3))/rmu
         endif
         fvec(i) = fvec(i)+target(3)*dreal(pot)/rmu
         if (ifstrain.eq.1) then
            dudx(i,1) = -target(3)*dreal(fld(1))/rmu
            dudx(i,2) = -target(3)*dreal(fld(2))/rmu
            dudx(i,3) = -target(3)*dreal(fld(3))/rmu
     1               +dreal(pot)/rmu
         endif
      enddo
c
      if (ifstrain.eq.1) then
         do i = 1,3
         do j = 1,3
            strainloc(i,j) = (dudx(i,j)+dudx(j,i))/2
         enddo
         enddo
c
         do i=1,3
         do j=1,3
            strain(i,j)=strain(i,j)+strainloc(i,j)
         enddo
         enddo
      endif
      return
      end
c
c
c
c
c
c***********************************************************************
        subroutine green3eluh_image_a_eval(rlam,rmu,source,df,target,
     $     fvec,ifstrain,strain)
c***********************************************************************
c
c       This subroutine evaluates the displacement vector <fvec> and 
c       strain tensor <strain> at the
c       location <target> due to the static force vector <df> at 
c       location <source>: the Mindlin solution
c
c       Mindlin's solution corresponds to a half-space solution
c       (z .leq. 0), which satisfies the condition that the normal 
c       component of stress vanishes at z=0.
c
c
c       INPUT:
c
c       rlam, rmu       Lame parameters
c       source(3)       Source location in lower half-space (z<0).
c       df(3)           Strength of single force 
c       target(3)       Target point in lower half-space (z<0).
c       ifstrain        Strain flag (1 means compute, 0 means do not).
c
c
c       OUTPUT:
c
c       fvec (real(8)) - the displacement at the target
c       strain (real(8)) - the strain at the target
c
c
      implicit none

      real(8), intent(in) :: rlam
      real(8), intent(in) :: rmu
      real(8), intent(in) :: source(3)
      real(8), intent(in) :: df(3)
      real(8), intent(in) :: target(3)
      real(8), intent(out) :: fvec(3)
      integer, intent(in) :: ifstrain
      real(8), intent(inout) :: strain(3,3)

      real(8) :: fvecloc(3)
      real(8) :: sourceim(3)
      real(8) :: targim(3)
      real(8) :: rlame(2)
      real(8) :: dudx(3,3)
      real(8) :: strainloc(3,3)
      real(8) :: dipstr(3)
      real(8) :: dipvec(3)

      real(8) :: rlamim
      real(8) :: rmuim

      integer :: i
      integer :: j
      integer :: ifsingle
      integer :: ifdouble

      complex(8) :: pot
      complex(8) :: fld(3)
c
c
c***********************************************************************
c     PART 2 === elastic <<anti-image>>   [see notes]
c     merges A image with first part of Stokeslet-like B image 
c     (in Okada 1992).
c***********************************************************************
c                        
      targim(1) = target(1)
      targim(2) = target(2)
      targim(3) = -target(3)
      rlamim = rlam+4*rmu
      rmuim = -rmu
c
      call green3elu_strain_image(rlamim,rmuim,source,df,targim,
     1     fvecloc,ifstrain,strainloc) 
c
c     inside the preceding subroutine,
c     scaling is by 1/(2*rmuim) and it should be by 1/(2*rmu).
c     This is fixed by flipping the sign (subtracting) in contributions
c     to displacement and strain here.
c
      fvec(1)=-fvecloc(1)
      fvec(2)=-fvecloc(2)
      fvec(3)=-fvecloc(3)
      if (ifstrain.eq.1) then
         do i=1,3
         do j=1,3
            strain(i,j)=-strainloc(i,j)
         enddo
         enddo
      endif  
      if (2.ne.3) return
c
c***********************************************************************
c       === PART 3 ===
c       Harmonic contributions from image source.
c***********************************************************************
c
c
      rlame(1) = rlam
      rlame(2) = rmu
      sourceim(1) = source(1)
      sourceim(2) = source(2)
      sourceim(3) = -source(3)
c
      ifsingle = 1
      ifdouble = 0
      do i=1,3
         call intker_mindlinb(i,rlame,sourceim,ifsingle,
     1           df,ifdouble,dipstr,dipvec,target,pot,ifstrain,fld)
         fvec(i) = fvec(i)+dreal(pot)/rmu
         if (ifstrain.eq.1) then
            dudx(i,1) = -dreal(fld(1))/rmu
            dudx(i,2) = -dreal(fld(2))/rmu
            dudx(i,3) = -dreal(fld(3))/rmu
         endif
      enddo
      if (ifstrain.eq.1) then
         do i = 1,3
         do j = 1,3
            strainloc(i,j) = (dudx(i,j)+dudx(j,i))/2
         enddo
         enddo
c
         do i=1,3
         do j=1,3
            strain(i,j)=strain(i,j)+strainloc(i,j)
         enddo
         enddo
      endif
      do i=1,3
         call intker_mindlinc(i,rlame,sourceim,ifsingle,
     1           df,ifdouble,dipstr,dipvec,target,pot,ifstrain,fld)
         if (ifstrain.eq.1) then
            dudx(i,1) = dreal(fld(1))/rmu
            dudx(i,2) = dreal(fld(2))/rmu
            dudx(i,3) = dreal(fld(3))/rmu
         endif
         fvec(i) = fvec(i)+target(3)*dreal(pot)/rmu
         if (ifstrain.eq.1) then
            dudx(i,1) = -target(3)*dreal(fld(1))/rmu
            dudx(i,2) = -target(3)*dreal(fld(2))/rmu
            dudx(i,3) = -target(3)*dreal(fld(3))/rmu
     1               +dreal(pot)/rmu
         endif
      enddo
c
      if (ifstrain.eq.1) then
         do i = 1,3
         do j = 1,3
            strainloc(i,j) = (dudx(i,j)+dudx(j,i))/2
         enddo
         enddo
c
         do i=1,3
         do j=1,3
            strain(i,j)=strain(i,j)+strainloc(i,j)
         enddo
         enddo
      endif
      return
      end
c
c
c
c
c
c***********************************************************************
        subroutine green3eluh_mindlin_eval(rlam,rmu,source,df,target,
     $     fvec,ifstrain,strain)
c***********************************************************************
c
c       This subroutine evaluates the displacement vector <fvec> and 
c       strain tensor <strain> at the
c       location <target> due to the static force vector <df> at 
c       location <source>: the Mindlin solution
c
c       Mindlin's solution corresponds to a half-space solution
c       (z .leq. 0), which satisfies the condition that the normal 
c       component of stress vanishes at z=0.
c
c       Mindlin image part only, parts B and C
c
c       INPUT:
c
c       rlam, rmu       Lame parameters
c       source(3)       Source location in lower half-space (z<0).
c       df(3)           Strength of single force 
c       target(3)       Target point in lower half-space (z<0).
c       ifstrain        Strain flag (1 means compute, 0 means do not).
c
c
c       OUTPUT:
c
c       fvec (real(8)) - the displacement at the target
c       strain (real(8)) - the strain at the target
c
c
      implicit none

      real(8), intent(in) :: rlam
      real(8), intent(in) :: rmu
      real(8), intent(in) :: source(3)
      real(8), intent(in) :: df(3)
      real(8), intent(in) :: target(3)
      real(8), intent(out) :: fvec(3)
      integer, intent(in) :: ifstrain
      real(8), intent(inout) :: strain(3,3)

      real(8) :: sourceim(3)
      real(8) :: targim(3)
      real(8) :: rlame(2)
      real(8) :: dudx(3,3)
      real(8) :: strainloc(3,3)
      real(8) :: dipstr(3)
      real(8) :: dipvec(3)

      real(8) :: rlamim
      real(8) :: rmuim

      integer :: i
      integer :: j
      integer :: ifsingle
      integer :: ifdouble

      complex(8) :: pot
      complex(8) :: fld(3)
c
c
c***********************************************************************
c     PART 2 === elastic <<anti-image>>   [see notes]
c     merges A image with first part of Stokeslet-like B image 
c     (in Okada 1992).
c***********************************************************************
c                        
      targim(1) = target(1)
      targim(2) = target(2)
      targim(3) = -target(3)
      rlamim = rlam+4*rmu
      rmuim = -rmu
c
ccc      call green3elu_strain_image(rlamim,rmuim,source,df,targim,
ccc     1     fvecloc,ifstrain,strainloc) 
c
c     inside the preceding subroutine,
c     scaling is by 1/(2*rmuim) and it should be by 1/(2*rmu).
c     This is fixed by flipping the sign (subtracting) in contributions
c     to displacement and strain here.
c
      fvec(1)=0
      fvec(2)=0
      fvec(3)=0
      if (ifstrain.eq.1) then
         do i=1,3
         do j=1,3
            strain(i,j)=0
         enddo
         enddo
      endif  
cc      if (2.ne.3) return
c
c***********************************************************************
c       === PART 3 ===
c       Harmonic contributions from image source.
c***********************************************************************
c
c
      rlame(1) = rlam
      rlame(2) = rmu
      sourceim(1) = source(1)
      sourceim(2) = source(2)
      sourceim(3) = -source(3)
c
      ifsingle = 1
      ifdouble = 0
      do i=1,3
         call intker_mindlinb(i,rlame,sourceim,ifsingle,
     1           df,ifdouble,dipstr,dipvec,target,pot,ifstrain,fld)
         fvec(i) = fvec(i)+dreal(pot)/rmu
         if (ifstrain.eq.1) then
            dudx(i,1) = -dreal(fld(1))/rmu
            dudx(i,2) = -dreal(fld(2))/rmu
            dudx(i,3) = -dreal(fld(3))/rmu
         endif
      enddo
      if (ifstrain.eq.1) then
         do i = 1,3
         do j = 1,3
            strainloc(i,j) = (dudx(i,j)+dudx(j,i))/2
         enddo
         enddo
c
         do i=1,3
         do j=1,3
            strain(i,j)=strain(i,j)+strainloc(i,j)
         enddo
         enddo
      endif
      do i=1,3
         call intker_mindlinc(i,rlame,sourceim,ifsingle,
     1           df,ifdouble,dipstr,dipvec,target,pot,ifstrain,fld)
         if (ifstrain.eq.1) then
            dudx(i,1) = dreal(fld(1))/rmu
            dudx(i,2) = dreal(fld(2))/rmu
            dudx(i,3) = dreal(fld(3))/rmu
         endif
         fvec(i) = fvec(i)+target(3)*dreal(pot)/rmu
         if (ifstrain.eq.1) then
            dudx(i,1) = -target(3)*dreal(fld(1))/rmu
            dudx(i,2) = -target(3)*dreal(fld(2))/rmu
            dudx(i,3) = -target(3)*dreal(fld(3))/rmu
     1               +dreal(pot)/rmu
         endif
      enddo
c
      if (ifstrain.eq.1) then
         do i = 1,3
         do j = 1,3
            strainloc(i,j) = (dudx(i,j)+dudx(j,i))/2
         enddo
         enddo
c
         do i=1,3
         do j=1,3
            strain(i,j)=strain(i,j)+strainloc(i,j)
         enddo
         enddo
      endif
      return
      end
c
c
c
c
c
c**********************************************************************C
        subroutine green3eluh_brute(rlam,rmu,source,target,df,fout)
c**********************************************************************C
c
c
c       Naive implementation of elastic SLP with zero normal stress
c       on lower half-space (Mindlin solution) using Okada formula.
c
c       Half-space boundary condition is assumed, i.e. normal stress
c       component is zero at z=0.
c
c       INPUT:
c
c       rlam, rmu       Lame parameters
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       df(3)           Strength of single force source
c       rnorm(3)        Orientation vector of the double force source
c
c       OUTPUT:
c
c       fout (real(8)) - the displacement at the target
c
        implicit none

        real(8), intent(in) :: rlam
        real(8), intent(in) :: rmu
        real(8), intent(in) :: source(3)
        real(8), intent(in) :: target(3)
        real(8), intent(in) :: df(3)
        real(8), intent(out) :: fout(3)

        real(8) :: fvec(3)
        real(8) :: delta(3,3)
        real(8) :: rr(3)

        real(8) :: alpha
        real(8) :: dx
        real(8) :: dy
        real(8) :: dz
        real(8) :: cd

        integer :: i
        integer :: j

c
        alpha=(rlam+rmu)/(rlam+2*rmu)
c
        do i=1,3
        do j=1,3
        delta(i,j)=0
        enddo
        enddo
c
        do i=1,3
        delta(i,i)=1
        enddo
c
c
c       === PART 1 ===
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)-source(3)
        rr(1)=dx
        rr(2)=dy
        rr(3)=dz
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
        fvec(1)=0
        fvec(2)=0
        fvec(3)=0
c
        do i=1,3
        do j=1,3
        fvec(i)=fvec(i)+((2-alpha)*delta(i,j)/cd) * df(j)
        fvec(i)=fvec(i)+(alpha*rr(i)*rr(j)/cd**3) * df(j)
        enddo
        enddo
c
c       ... scale by 1/(2*rmu)
c
        fvec(1)=fvec(1)/(2*rmu)
        fvec(2)=fvec(2)/(2*rmu)
        fvec(3)=fvec(3)/(2*rmu)
c
        fout(1)=fvec(1)
        fout(2)=fvec(2)
        fout(3)=fvec(3)
c
ccc        if (2.ne.3) return
c
c
c       === PART 2 ===
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)+source(3)
        rr(1)=dx
        rr(2)=dy
        rr(3)=-dz
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
        fvec(1)=0
        fvec(2)=0
        fvec(3)=0
c
        do i=1,3
        do j=1,3
        fvec(i)=fvec(i)+((2-alpha)*delta(i,j)/cd) * df(j)
        fvec(i)=fvec(i)+(alpha*rr(i)*rr(j)/cd**3) * df(j)
        enddo
        enddo
ccc        fvec(1) = 0
ccc        fvec(2) = 0
ccc        fvec(3) = 0
c
c       ... scale by 1/(2*rmu)
c
        fvec(1)=fvec(1)/(2*rmu)
        fvec(2)=fvec(2)/(2*rmu)
        fvec(3)=fvec(3)/(2*rmu)
c
        fvec(1)=-fvec(1)
        fvec(2)=-fvec(2)
        fvec(3)=-fvec(3)
c
        fout(1)=fout(1)+fvec(1)
        fout(2)=fout(2)+fvec(2)
        fout(3)=fout(3)+fvec(3)
c
ccc        if (2.ne.3) return
c
c
c       === PART 3 ===
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)+source(3)
        rr(1)=dx
        rr(2)=dy
        rr(3)=-dz
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
        fvec(1)=0
        fvec(2)=0
        fvec(3)=0
c
        do i=1,3
        do j=1,3
c
        fvec(i)=fvec(i)+(delta(i,j)/cd+rr(i)*rr(j)/cd**3) *df(j)
c
        fvec(i)=fvec(i)+(1-alpha)/alpha*
     $   ( delta(i,j)/(cd+rr(3)) +
     $     (rr(i)*delta(j,3)-rr(j)*delta(i,3)*(1-delta(j,3)))
     $     /cd/(cd+rr(3))
     $     -rr(i)*rr(j)*(1-delta(j,3))*(1-delta(i,3))/cd/(cd+rr(3))**2 )
     $     *df(j) 
c
        enddo
        enddo
c
c       ... scale by 1/(rmu)
c
        fvec(1)=fvec(1)/(rmu)
        fvec(2)=fvec(2)/(rmu)
        fvec(3)=fvec(3)/(rmu)
c
        fout(1)=fout(1)+fvec(1)
        fout(2)=fout(2)+fvec(2)
        fout(3)=fout(3)+fvec(3)
c
c
c       === PART 4 ===
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)+source(3)
        rr(1)=dx
        rr(2)=dy
        rr(3)=-dz
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
        fvec(1)=0
        fvec(2)=0
        fvec(3)=0
c
        do i=1,3
        do j=1,3
c
        fvec(i)=fvec(i)+(1-2*delta(i,3))*
     $     ( (2-alpha)*(rr(i)*delta(j,3)-rr(j)*delta(i,3))/cd**3
     $     + alpha*source(3)*(delta(i,j)/cd**3-3*rr(i)*rr(j)/cd**5))
     $     * df(j)
c
        enddo
        enddo
c
c       ... scale by 1/(rmu)
c
        fvec(1)=fvec(1)/(rmu)
        fvec(2)=fvec(2)/(rmu)
        fvec(3)=fvec(3)/(rmu)
c
        fout(1)=fout(1)+target(3)*fvec(1)
        fout(2)=fout(2)+target(3)*fvec(2)
        fout(3)=fout(3)+target(3)*fvec(3)
c
c
        return
        end
c
c
c
c
c
c 
c**********************************************************************C
        subroutine green3elth_brute(rlam,rmu,source,target,
     1           du,rnorm,fout)
c**********************************************************************C
c
c
c       Naive implementation of elastic DLP with zero normal stress
c       on lower half-space (Mindlin solution) using Okada formula.
c
c       INPUT:
c
c       rlam, rmu       Lame parameters
c       source(3)       Source location in lower half-space (z<0).
c       target(3)       Target point in lower half-space (z<0).
c       du(3)           Strength of double force source
c       rnorm(3)        Orientation vector of the double force source
c
c       OUTPUT:
c
c       fout (real(8)) - the displacement at the target
c
        implicit none

        real(8), intent(in) :: rlam
        real(8), intent(in) :: rmu
        real(8), intent(in) :: source(3)
        real(8), intent(in) :: target(3)
        real(8), intent(in) :: du(3)
        real(8), intent(in) :: rnorm(3)
        real(8), intent(out) :: fout(3)

        real(8) :: xyz(3)
        real(8) :: fvec(3)
        real(8) :: delta(3,3)
        real(8) :: rr(3)
        real(8) :: uder(3,3,3)
        real(8) :: uuu(3)

        real(8) :: alpha
        real(8) :: cd

        integer :: i
        integer :: j
        integer :: k
        integer :: n
c
        alpha=(rlam+rmu)/(rlam+2*rmu)
c
        xyz(1) = target(1)-source(1)
        xyz(2) = target(2)-source(2)
        xyz(3) = target(3)-source(3)
        do i=1,3
        do j=1,3
        delta(i,j)=0
        enddo
        enddo
c
        do i=1,3
        delta(i,i)=1
        enddo
c
        rr(1) = xyz(1)
        rr(2) = xyz(2)
        rr(3) = xyz(3)
        cd=sqrt(rr(1)**2+rr(2)**2+rr(3)**2)
c
        do i = 1,3
        do j = 1,3
        do k = 1,3
           uder(i,j,k)= (2-alpha)*rr(k)*delta(i,j)/cd**3
     1          -alpha*(rr(i)*delta(j,k)+rr(j)*delta(i,k))/cd**3
     2          +3*alpha*rr(i)*rr(j)*rr(k)/cd**5
        enddo
        enddo
        enddo
c
        fvec(1)=0
        fvec(2)=0
        fvec(3)=0
c
        do i = 1,3
           uuu(i) = 0
           do n=1,3
              uuu(i) = uuu(i)+uder(i,n,n)
           enddo
        enddo
c
        do i=1,3
        do j=1,3
        do k=1,3
           fvec(i) = fvec(i) + 
     1               rmu*(uder(i,j,k)+uder(i,k,j))*du(j)*rnorm(k)
           fvec(i) = fvec(i) + rlam*delta(j,k)*uuu(i)*du(j)*rnorm(k)
        enddo
        enddo
        enddo
c
c       ... scale by 1/(2*rmu)
c
        fout(1)=fvec(1)/(2*rmu)
        fout(2)=fvec(2)/(2*rmu)
        fout(3)=fvec(3)/(2*rmu)
c
ccc       call prin2(' after direct arrival fout is *',fout,3)
ccc       if (2.ne.3) return
ccc       if (2.ne.3) goto 111
c
c      subtract A image
c
c
        rr(1) = xyz(1)
        rr(2) = xyz(2)
        rr(3) = -target(3)-source(3)
        cd=sqrt(rr(1)**2+rr(2)**2+rr(3)**2)
c
        do i = 1,3
        do j = 1,3
        do k = 1,3
           uder(i,j,k)= (2-alpha)*rr(k)*delta(i,j)/cd**3
     1          -alpha*(rr(i)*delta(j,k)+rr(j)*delta(i,k))/cd**3
     2          +3*alpha*rr(i)*rr(j)*rr(k)/cd**5
        enddo
        enddo
        enddo
c
        fvec(1)=0
        fvec(2)=0
        fvec(3)=0
c
        do i = 1,3
           uuu(i) = 0
           do n=1,3
              uuu(i) = uuu(i)+uder(i,n,n)
           enddo
        enddo
c
        do i=1,3
        do j=1,3
        do k=1,3
           fvec(i) = fvec(i) + 
     1               rmu*(uder(i,j,k)+uder(i,k,j))*du(j)*rnorm(k)
           fvec(i) = fvec(i) + rlam*delta(j,k)*uuu(i)*du(j)*rnorm(k)
        enddo
        enddo
        enddo
c
c       ... scale by 1/(2*rmu)
c
        fvec(1)=fvec(1)/(2*rmu)
        fvec(2)=fvec(2)/(2*rmu)
        fvec(3)=fvec(3)/(2*rmu)
        fout(1)=fout(1)-fvec(1)
        fout(2)=fout(2)-fvec(2)
        fout(3)=fout(3)-fvec(3)
c
c
ccc       if (2.ne.3) return
c
ccc       call prin2(' A image arrival fvec is *',fvec,3)
c
c      add first part of B image
c
c
        rr(1) = xyz(1)
        rr(2) = xyz(2)
        rr(3) = -target(3)-source(3)
        cd=sqrt(rr(1)**2+rr(2)**2+rr(3)**2)
c
        do i = 1,3
        do j = 1,3
        do k = 1,3
           uder(i,j,k)= rr(k)*delta(i,j)/cd**3
     1          -(rr(i)*delta(j,k)+rr(j)*delta(i,k))/cd**3
     2          +3*rr(i)*rr(j)*rr(k)/cd**5
        enddo
        enddo
        enddo
c
c
        fvec(1)=0
        fvec(2)=0
        fvec(3)=0
c
        do i = 1,3
           uuu(i) = 0
           do n=1,3
              uuu(i) = uuu(i)+uder(i,n,n)
           enddo
        enddo
c
        do i=1,3
        do j=1,3
        do k=1,3
           fvec(i) = fvec(i) + 
     1               rmu*(uder(i,j,k)+uder(i,k,j))*du(j)*rnorm(k)
           fvec(i) = fvec(i) + rlam*delta(j,k)*uuu(i)*du(j)*rnorm(k)
        enddo
        enddo
        enddo
c
c       ... scale by 1/(2*rmu)
c
        fvec(1)=fvec(1)/(rmu)
        fvec(2)=fvec(2)/(rmu)
        fvec(3)=fvec(3)/(rmu)
        fout(1)=fout(1)+fvec(1)
        fout(2)=fout(2)+fvec(2)
        fout(3)=fout(3)+fvec(3)
c
ccc       call prin2(' B image arrival fvec is *',fvec,3)
c
111    continue
ccc       if (2.ne.3) return
c
c      second part of B image
c
c
        rr(1) = xyz(1)
        rr(2) = xyz(2)
        rr(3) = -target(3)-source(3)
        cd=sqrt(rr(1)**2+rr(2)**2+rr(3)**2)
c
        do i = 1,3
        do j = 1,3
        do k = 1,3
           uder(i,j,k)=(delta(3,k)*cd+rr(k))*delta(i,j)/
     1         (cd*(cd+rr(3))**2)
     1    -(delta(i,k)*delta(j,3)-delta(j,k)*delta(i,3)*(1-delta(j,3)))
     1         /(cd*(cd+rr(3)))
     1    +( rr(i)*delta(j,3)-rr(j)*delta(i,3)*(1-delta(j,3)))*
     1    ( delta(3,k)*cd**2 + rr(k)*(2*cd+rr(3)))/
     1               ((cd**3)*(cd+rr(3))**2) 
     1    +( (rr(i)*delta(j,k)+rr(j)*delta(i,k))/(cd*(cd+rr(3))**2)-
     1     rr(i)*rr(j)*
     1     ( 2*delta(3,k)*cd**2+rr(k)*(3*cd+rr(3)))/
     1          ((cd**3)*(cd+rr(3))**3))*(1-delta(i,3))*(1-delta(j,3))
           uder(i,j,k)=uder(i,j,k)*(1-alpha)/alpha
        enddo
        enddo
        enddo
c
c
        fvec(1)=0
        fvec(2)=0
        fvec(3)=0
c
        do i = 1,3
           uuu(i) = 0
           do n=1,3
              uuu(i) = uuu(i)+uder(i,n,n)
           enddo
        enddo
c
        do i=1,3
        do j=1,3
        do k=1,3
           fvec(i) = fvec(i) + 
     1               rmu*(uder(i,j,k)+uder(i,k,j))*du(j)*rnorm(k)
           fvec(i) = fvec(i) + rlam*delta(j,k)*uuu(i)*du(j)*rnorm(k)
        enddo
        enddo
        enddo
c
c       ... scale by 1/(rmu)
c
        fvec(1)=fvec(1)/(rmu)
        fvec(2)=fvec(2)/(rmu)
        fvec(3)=fvec(3)/(rmu)
        fout(1)=fout(1)+fvec(1)
        fout(2)=fout(2)+fvec(2)
        fout(3)=fout(3)+fvec(3)
ccc        fout(1)=fvec(1)
ccc        fout(2)=fvec(2)
ccc        fout(3)=fvec(3)
c
ccc       call prin2(' second part of B image fvec is *',fvec,3)
ccc       if (2.ne.3) return
c
c      C image
c
c
        rr(1) = xyz(1)
        rr(2) = xyz(2)
        rr(3) = -target(3)-source(3)
        cd=sqrt(rr(1)**2+rr(2)**2+rr(3)**2)
c
ccc           write(6,*) 'from Okada'
ccc           write(13,*) 'from Okada'
        do i = 1,3
        do j = 1,3
        do k = 1,3
           uder(i,j,k)= (1-2*delta(i,3))*(
     1       (2-alpha)*( 
     1         (delta(j,k)*delta(i,3)-delta(i,k)*delta(j,3))/(cd**3)+
     1         3*rr(k)*( rr(i)*delta(j,3) - rr(j)*delta(i,3))/(cd**5))
     1    + alpha*(delta(i,j)/(cd**3)-3*rr(i)*rr(j)/(cd**5))*delta(3,k)
     1    + 3*alpha*source(3)*(
     1       (rr(i)*delta(j,k)+rr(j)*delta(i,k)+
     1                 rr(k)*delta(i,j))/(cd**5)
     1       - 5*rr(i)*rr(j)*rr(k)/(cd**7) ))
           uder(i,j,k)= uder(i,j,k)*target(3)
ccc           write(6,*) i,j,k,uder(i,j,k)
ccc           write(13,*) i,j,k,uder(i,j,k)
        enddo
        enddo
        enddo
c
c
c
        fvec(1)=0
        fvec(2)=0
        fvec(3)=0
c
        do i = 1,3
           uuu(i) = 0
           do n=1,3
              uuu(i) = uuu(i)+uder(i,n,n)
           enddo
        enddo
c
        do i=1,3
        do j=1,3
        do k=1,3
           fvec(i) = fvec(i) + 
     1               rmu*(uder(i,j,k)+uder(i,k,j))*du(j)*rnorm(k)
           fvec(i) = fvec(i) + rlam*delta(j,k)*uuu(i)*du(j)*rnorm(k)
        enddo
        enddo
        enddo
c
c       ... scale by 1/(rmu)
c
        fvec(1)=fvec(1)/(rmu)
        fvec(2)=fvec(2)/(rmu)
        fvec(3)=fvec(3)/(rmu)
        fout(1)=fout(1)+fvec(1)
        fout(2)=fout(2)+fvec(2)
        fout(3)=fout(3)+fvec(3)
ccc        fout(1)=fvec(1)
ccc        fout(2)=fvec(2)
ccc        fout(3)=fvec(3)
c
ccc       call prin2(' C image fvec is *',fvec,3)

        return
        end
c
c
