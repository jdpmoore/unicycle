cc Copyright (C) 2009: Leslie Greengard and Zydrunas Gimbutas
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
c      This file contains simple geometry I/O and processing routines.
c
c      atrireadchk:  reads simple (single component) .tri format file 
c                 and stores in fortran arrays
c                 includes some dimension checking and error flags
c
c      atriread:  reads simple (single component) .tri format file 
c                 and stores in fortran arrays
c                 No dimension checking or error flags
c
c      atriwrite: writes geometry description (single component) to 
c                 file in .tri format 
c
c      atriwrite2: writes geometry description (single component) to 
c                 file in .tri format in reordered format. (CHECK) 
c
c      aoffwrite: writes geometry description (single component) to 
c                 file in GEOMVIEW/OOGL .off format 
c
c      buildtri:  created triangles/normals/areas from verts.
c
c      triangle_norm:  computes normals from triangles.
c      triangle_area:  computes areas from triangles.
c
c      tri_ini: constructs a mapping that translates and
c       rotates the triangle T defined by the vertices
c       vert1, vert2, vert3 (in R^3) onto the
c       x-y plane in such a way that vert1 is mapped to the origin and 
c       the segment vert1-vert2 lies on the (positive) x-axis. 
c
c       Calls rotshift3d
c
c      rotder3d: rotates gradient from triangle coords to lab frame.
c                Needs initial call to tri_ini
c      tri_for: applies transl/rot to arbitrary point
c      tri_bak: applies inverse transformation
c
c
c      tri_for_vect: vector forward rotation
c      tri_bak_vect: vector backward rotation
c      rothess3d: hessian rotation 
c
c***********************************************************************
        subroutine atrireadchk(ir,verts,lv,nverts,ifaces,li,nfaces,ier)
c***********************************************************************
        implicit none

        integer, intent(in) :: ir
        real(8), intent(inout) :: verts(3,*)
        integer, intent(in) :: lv
        integer, intent(out) :: nverts
        integer, intent(inout) :: ifaces(3,*)
        integer, intent(in) :: li
        integer, intent(out) :: nfaces
        integer, intent(out) :: ier

        integer :: istat
        integer :: j
c
c       read geometry description file in (.tri) format. 
c
        ier = 0
        open(unit=ir, status="OLD", iostat=istat)
        if (istat .ne. 0) then
          write(*,*)"in atriread, OPEN statement iostat = ", istat
          ier = 1
          return
        endif
c
        read(ir,*) nVerts, nFaces
        if (lv.lt.nVerts) then
          ier = 2
          return
        endif
c
        if (li.lt.nFaces) then
          ier = 3
          return
        endif
c
        read(ir,*) (verts(1,j),verts(2,j),verts(3,j),j=1,nVerts) 
        read(ir,*) (ifaces(1,j),ifaces(2,j),ifaces(3,j),j=1,nFaces) 
c
        return
        end
c
c
c
c
c***********************************************************************
        subroutine atriread(ir,verts,nverts,ifaces,nfaces)
c***********************************************************************
        implicit none

        integer, intent(in) :: ir
        real(8), intent(inout) :: verts(3,*)
        integer, intent(out) :: nverts
        integer, intent(inout) :: ifaces(3,*)
        integer, intent(out) :: nfaces

        integer :: istat
        integer :: j
c
c       read geometry description file in (.tri) format. 
c
        open(unit=ir, status="OLD", iostat=istat)
        if (istat .ne. 0) then
          write(*,*)"in atriread, OPEN statement iostat = ", istat
          stop
        endif
c
        read(ir,*) nVerts, nFaces
c
        read(ir,*) (verts(1,j),verts(2,j),verts(3,j),j=1,nVerts) 
        read(ir,*) (ifaces(1,j),ifaces(2,j),ifaces(3,j),j=1,nFaces) 
c
        return
        end
c
c
c
c
c***********************************************************************
        subroutine atriwrite(iw,verts,nverts,ifaces,nfaces)
c***********************************************************************
        implicit none

        integer :: iw
        integer :: nverts
        real(8) :: verts(3,nverts)
        integer :: nfaces
        integer :: ifaces(3,nfaces)

        integer :: j
c
c       write geometry description file in cart3d (.tri) format. 
c
        write(iw,*) nverts, nfaces
        write(iw,1800) (verts(1,j),verts(2,j),verts(3,j),j=1,nverts) 
        write(iw,1900) (ifaces(1,j),ifaces(2,j),ifaces(3,j),j=1,nfaces) 
 1800   format(3(1x,e23.16))
 1900   format(3(1x,i7))
c
        return
        end


c***********************************************************************
        subroutine atriwrite2(iw,verts,nverts,ifaces,itrinew,nfaces)
c***********************************************************************
        implicit none

        integer, intent(in) :: iw
        integer, intent(in) :: nverts
        real(8), intent(in) :: verts(3,nverts)
        integer, intent(in) :: nfaces
        integer, intent(in) :: ifaces(3,*)
        integer, intent(in) :: itrinew(nfaces)

        integer :: j
c
c       write geometry description file in cart3d (.tri) format. 
c

        write(iw,*) nVerts, nFaces
        write(iw,1200) (verts(1,j),verts(2,j),verts(3,j),j=1,nVerts) 
        write(iw,1400) (ifaces(1,itrinew(j)),ifaces(2,itrinew(j)),
     1      ifaces(3,itrinew(j)),j=1,nFaces) 
 1200   format(3(1x,e23.16))
 1400   format(3(1x,i7))

c
        return
        end
c
c
c
c***********************************************************************
      subroutine aoffwrite(iw,verts,nverts,ifaces,nfaces)
c***********************************************************************
c     
c     write triangles into an OOGL OFF file
c     
      implicit none

      integer, intent(in) :: iw
      integer, intent(in) :: nverts
      real(8), intent(in) :: verts(3,nverts)
      integer, intent(in) :: nfaces
      integer, intent(in) :: ifaces(3,nfaces)

      integer :: nedges

      integer :: j
      
ccc   nb of edges of a triangle is 3
      nedges = 3
c     
      write(iw,*) nVerts, nFaces
      write(iw,1200) (verts(1,j),verts(2,j),verts(3,j),j=1,nVerts) 
      write(iw,1400) (nedges,ifaces(1,j),ifaces(2,j),
     1     ifaces(3,j),j=1,nFaces) 
 1200 format(3(1x,e23.16))
 1400 format(4(1x,i7))
c     
      return
      end


c***********************************************************************
      subroutine buildtri(verts,nverts,ntri,itrivert,
     1     triangle,centroids,trinorm,triarea)
c***********************************************************************
c
c     Build triangles, normals and centroids from the native (TRI) 
c     file format has been factorized in this subroutine 
c
c     INPUT: 
c     verts(3,nverts)       array of vertices
c     ntri                  number of triangles
c     itrivert(3,ntri)      indices of triangle vertices
c     
c     OUTPUT:
c     triangle(3,3,ntri)    array of triangles in standard format
c     centroids(3,ntri)     array of triangle centroids
c     trinorm(3,ntri)       array of triangle normals
c     triarea(ntri)         array of triangle areas
c
      implicit none

      integer, intent(in) :: nverts
      real(8), intent(in) :: verts(3,nverts)
      integer, intent(in) :: ntri
      integer, intent(in) :: itrivert(3,ntri)
      real(8), intent(inout) :: triangle(3,3,ntri)
      real(8), intent(inout) :: centroids(3,ntri)
      real(8), intent(inout) :: trinorm(3,ntri)
      real(8), intent(inout) :: triarea(ntri)

      integer :: itri
      real(8) :: v11
      real(8) :: v21
      real(8) :: v31
      real(8) :: v12
      real(8) :: v22
      real(8) :: v32
      real(8) :: v13
      real(8) :: v23
      real(8) :: v33
c
      do itri = 1,ntri
         v11 = verts(1,itrivert(1,itri))
         v21 = verts(2,itrivert(1,itri))
         v31 = verts(3,itrivert(1,itri))
         v12 = verts(1,itrivert(2,itri))
         v22 = verts(2,itrivert(2,itri))
         v32 = verts(3,itrivert(2,itri))
         v13 = verts(1,itrivert(3,itri))
         v23 = verts(2,itrivert(3,itri))
         v33 = verts(3,itrivert(3,itri))

         triangle(1,1,itri) = v11
         triangle(2,1,itri) = v21
         triangle(3,1,itri) = v31
         triangle(1,2,itri) = v12
         triangle(2,2,itri) = v22
         triangle(3,2,itri) = v32
         triangle(1,3,itri) = v13
         triangle(2,3,itri) = v23
         triangle(3,3,itri) = v33
         centroids(1,itri) = (v11+v12+v13)/3
         centroids(2,itri) = (v21+v22+v23)/3
         centroids(3,itri) = (v31+v32+v33)/3
         call triangle_norm(triangle(1,1,itri),trinorm(1,itri))
         call triangle_area(triangle(1,1,itri),triarea(itri))
      enddo
      return
      end
c 
c 
c
c
c
c
c***********************************************************************
        subroutine triangle_norm(triang,trinorm)
c***********************************************************************
        implicit none

        real(8), intent(in) :: triang(3,3)
        real(8), intent(out) :: trinorm(3)

        real(8) :: x(3)
        real(8) :: y(3)
        real(8) :: scale
c
c       constructs normal vector of positively
c       oriented triangle
c
c       INPUT: 
c   
c       triang(3,3)     triangle vertices in standard format
c
c       OUTPUT: 
c   
c       trinorm(3)     triangle vertices in standard format
c
c
        x(1)=triang(1,2)-triang(1,1)
        x(2)=triang(2,2)-triang(2,1)
        x(3)=triang(3,2)-triang(3,1)
c
        y(1)=triang(1,3)-triang(1,2)
        y(2)=triang(2,3)-triang(2,2)
        y(3)=triang(3,3)-triang(3,2)
c
        trinorm(1)= x(2)*y(3)-x(3)*y(2)
        trinorm(2)= x(3)*y(1)-x(1)*y(3)
        trinorm(3)= x(1)*y(2)-x(2)*y(1)
c
        scale=dsqrt(trinorm(1)**2+trinorm(2)**2+trinorm(3)**2)
        trinorm(1)=trinorm(1)/scale
        trinorm(2)=trinorm(2)/scale
        trinorm(3)=trinorm(3)/scale
c
        return
        end
c
c
c
c
c
c***********************************************************************
        subroutine triangle_area(triang,area)
c***********************************************************************
        implicit none

        real(8), intent(in) :: triang(3,3)
        real(8), intent(out) :: area

        real(8) :: x(3)
        real(8) :: y(3)
        real(8) :: z(3)
c
c       compute area of positively oriented triangle in R^3
c
c       INPUT: 
c   
c       triang(3,3)     triangle vertices in standard format
c
c       OUTPUT: 
c   
c       area            triangle area
c
        x(1)=triang(1,2)-triang(1,1)
        x(2)=triang(2,2)-triang(2,1)
        x(3)=triang(3,2)-triang(3,1)
c
        y(1)=triang(1,3)-triang(1,1)
        y(2)=triang(2,3)-triang(2,1)
        y(3)=triang(3,3)-triang(3,1)
c
        z(1)= x(2)*y(3)-x(3)*y(2)
        z(2)= x(3)*y(1)-x(1)*y(3)
        z(3)= x(1)*y(2)-x(2)*y(1)
c
        area=dsqrt(z(1)**2+z(2)**2+z(3)**2)
        area=area*0.5d0
c
        return
        end
c
c
c
c
c
c***********************************************************************
        subroutine rotder3d(w,triang,valx,valy,valz,derx,dery,derz)
c***********************************************************************
        implicit none

        real(8), intent(in) :: w(12)
        real(8), intent(in) :: triang(3,3)
        real(8), intent(in) :: valx
        real(8), intent(in) :: valy
        real(8), intent(in) :: valz
        real(8), intent(out) :: derx
        real(8), intent(out) :: dery
        real(8), intent(out) :: derz

        real(8) :: derin(3)
        real(8) :: derout(3)
c
c       transform the gradient in the rotated triangle coordinate 
c       system to the gradient in the original Cartesian coordinate 
c       system.
c
c       The array w must be initialized by a preceeding call to
c       tri_ini(triang(1,1),triang(1,2),triang(1,3),w)
c
        derin(1)=valx
        derin(2)=valy
        derin(3)=valz
        call tri_bak(w,derin,derout)
        derx=derout(1)-triang(1,1)
        dery=derout(2)-triang(2,1)
        derz=derout(3)-triang(3,1)
        return
        end
c
c
c
c
c
c
c
c***********************************************************************
        subroutine rotder3db(w,triang,valx,valy,valz,derx,dery,derz)
c***********************************************************************
        implicit none

        real(8), intent(in) :: w(12)
        real(8), intent(in) :: triang(3,3)
        real(8), intent(in) :: valx
        real(8), intent(in) :: valy
        real(8), intent(in) :: valz
        real(8), intent(out) :: derx
        real(8), intent(out) :: dery
        real(8), intent(out) :: derz

        real(8) :: derin(3)
        real(8) :: derout(3)
c
c       transform the gradient in the rotated coordinate system to the
c       gradient in the original coordinate system.
c
c       The array w must be initialized by a preceeding call to
c       tri_ini(triang(1,1),triang(1,2),triang(1,3),w)
c
        derin(1)=valx
        derin(2)=valy
        derin(3)=valz
        call tri_bak(w,derin,derout)
        derx=derout(1)-triang(1,1)
        dery=derout(2)-triang(2,1)
        derz=derout(3)-triang(3,1)
        return
        end
c
c
c
c
c
c
c
c***********************************************************************
        subroutine rotder3df(w,triang,valx,valy,valz,derx,dery,derz)
c***********************************************************************
        implicit none

        real(8), intent(in) :: w(12)
        real(8), intent(in) :: triang(3,3)
        real(8), intent(in) :: valx
        real(8), intent(in) :: valy
        real(8), intent(in) :: valz
        real(8), intent(out) :: derx
        real(8), intent(out) :: dery
        real(8), intent(out) :: derz

        real(8) :: derin(3)
        real(8) :: derout(3)
c
c       transform the gradient in the original coordinate system to the
c       gradient in the rotated coordinate system.
c
c       The array w must be initialized by a preceeding call to
c       tri_ini(triang(1,1),triang(1,2),triang(1,3),w)
c
        derin(1)=valx
        derin(2)=valy
        derin(3)=valz
        call tri_for(w,derin,derout)
        derx=derout(1)-triang(1,1)
        dery=derout(2)-triang(2,1)
        derz=derout(3)-triang(3,1)
        return
        end
c
c
c
c
c
c
c
c***********************************************************************
        subroutine rotder3dr(u,z1,z2,z3,znew1,znew2,znew3)
c***********************************************************************
        implicit none

        real(8), intent(in) :: u(3,3)
        real(8), intent(in) :: z1
        real(8), intent(in) :: z2
        real(8), intent(in) :: z3
        real(8), intent(out) :: znew1
        real(8), intent(out) :: znew2
        real(8), intent(out) :: znew3
c
c       transform the gradient in the rotated triangle coordinate 
c       system to the gradient in the original Cartesian coordinate 
c       system.
c
c       u is rotation matrix - multiply by its transpose
c
        znew1=u(1,1)*z1+u(1,2)*z2+u(1,3)*z3
        znew2=u(2,1)*z1+u(2,2)*z2+u(2,3)*z3
        znew3=u(3,1)*z1+u(3,2)*z2+u(3,3)*z3
        return
        end
c
c       
c***********************************************************************
        subroutine rotder3dz(u,z1,z2,z3,znew1,znew2,znew3)
c***********************************************************************
        implicit none

        real(8), intent(in) :: u(3,3)
        complex(8), intent(in) :: z1
        complex(8), intent(in) :: z2
        complex(8), intent(in) :: z3
        complex(8), intent(out) :: znew1
        complex(8), intent(out) :: znew2
        complex(8), intent(out) :: znew3
c
c       transform the gradient in the rotated triangle coordinate 
c       system to the gradient in the original Cartesian coordinate 
c       system.
c
c       u is rotation matrix - multiply by its transpose
c
        znew1=u(1,1)*z1+u(1,2)*z2+u(1,3)*z3
        znew2=u(2,1)*z1+u(2,2)*z2+u(2,3)*z3
        znew3=u(3,1)*z1+u(3,2)*z2+u(3,3)*z3
        return
        end
c
c       
c***********************************************************************
        subroutine rotshift3d(vert1,vert2,vert3,rotmat,shift,v1,v2,v3)
c***********************************************************************
        implicit none

        real(8), intent(in) :: vert1(3)
        real(8), intent(in) :: vert2(3)
        real(8), intent(in) :: vert3(3)
        real(8), intent(out) :: rotmat(3,3)
        real(8), intent(out) :: shift(3)
        real(8), intent(out) :: v1(2)
        real(8), intent(out) :: v2(2)
        real(8), intent(out) :: v3(2)

        real(8) :: vec1(3)
        real(8) :: vec2(3)
        real(8) :: vec3(3)
        real(8) :: dd
c
c       rotate and translate triangle in space so that 
c       vert1-vert2 lies along the positive x-axis.
c
c       Step 1) Shift triangle to origin.
c
        v1(1) = 0.0d0
        v1(2) = 0.0d0
        vec1(1)=vert2(1)-vert1(1)
        vec1(2)=vert2(2)-vert1(2)
        vec1(3)=vert2(3)-vert1(3)
        dd = dsqrt(vec1(1)**2 + vec1(2)**2 + vec1(3)**2)
        v2(1) = dd 
        v2(2) = 0.0d0 
c
c       Step 2) Normalize the first side
c
        vec1(1)=vec1(1)/dd
        vec1(2)=vec1(2)/dd
        vec1(3)=vec1(3)/dd
        rotmat(1,1)=vec1(1)
        rotmat(2,1)=vec1(2)
        rotmat(3,1)=vec1(3)
c
c       Step 3) Get v1->v3  vector
c
        vec2(1) = vert3(1) - vert1(1)
        vec2(2) = vert3(2) - vert1(2)
        vec2(3) = vert3(3) - vert1(3)
c
c       Compute cross-product and normalize - constitutes third
c       column of rotation matrix
c
        vec3(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
        vec3(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
        vec3(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
c
        dd = dsqrt(vec3(1)**2 + vec3(2)**2 + vec3(3)**2)
        vec3(1)=vec3(1)/dd
        vec3(2)=vec3(2)/dd
        vec3(3)=vec3(3)/dd
        rotmat(1,3)=vec3(1)
        rotmat(2,3)=vec3(2)
        rotmat(3,3)=vec3(3)
c
c       Step 4) Generate orthogonal basis vector for second column
c       of rotation matrix by cross product of 1st and 3rd orthogonal
c       vectors. To preserve positive orientation, do (vec3 X vec1).
c
        rotmat(1,2) = vec3(2)*vec1(3) - vec3(3)*vec1(2)
        rotmat(2,2) = vec3(3)*vec1(1) - vec3(1)*vec1(3)
        rotmat(3,2) = vec3(1)*vec1(2) - vec3(2)*vec1(1)
c
c       Step 5) Construct shift vector (using fact that vert1 maps
c       to origin.
c
        shift(1)=-vert1(1)
        shift(2)=-vert1(2)
        shift(3)=-vert1(3)
c
c       map vert3 to xy plane using vec2 = vert3-vert1.
c
        v3(1)=rotmat(1,1)*vec2(1)+rotmat(2,1)*vec2(2)+
     1        rotmat(3,1)*vec2(3)
        v3(2)=rotmat(1,2)*vec2(1)+rotmat(2,2)*vec2(2)+
     1        rotmat(3,2)*vec2(3)
c
        return
        end
c
c
c
c
c
c***********************************************************************
        subroutine tri_ini(vert1,vert2,vert3,w,v1,v2,v3)     
c***********************************************************************
        implicit none

        real(8), intent(in) :: vert1(3)
        real(8), intent(in) :: vert2(3)
        real(8), intent(in) :: vert3(3)
        real(8), intent(out) :: w(12)
        real(8), intent(out) :: v1(2)
        real(8), intent(out) :: v2(2)
        real(8), intent(out) :: v3(2)

c              
c       this subroutine constructs a mapping that translates and
c       rotates the triangle T defined by the vertices
c       vert1, vert2, vert3 (all vertices are in R^3) onto the
c       x-y plane in such a way that vert1 is mapped to the origin and 
c       the segment vert1-vert2 lies on the (positive) x-axis. 
c
c       Once initialized by this routine, the transformations can be 
c       applied to an arbitrary point in R^3 via the 
c       entries 
c
c       tri_for:  to apply the transformation itself
c       tri_bak:  to apply the inverse.
c       
c       INPUT:
c       
c       vert1,vert2,vert3   the vertices of the triangle in the original
c                           coordinate system.
c       
c       OUTPUT:
c       
c       w                   an array of length 12 containing the 
c                           necessary rotation matrix and translation
c                           vector for subsequent use by 
c                           the entries tri_for, tri_bak.
c       v1,v2,v3   the vertices of the triangle in the plane.
c                  v1 is at origin, v2 is on the x-axis, etc.
c       
c
c       contruct the transformation
c       
        call rotshift3d(vert1,vert2,vert3,w(1:9),w(10:12),v1,v2,v3) 
        return
        end subroutine tri_ini
c       
c
c***********************************************************************
        subroutine tri_for(w,zin,zout)
c***********************************************************************
        implicit none

        real(8), intent(in) :: w(12)
        real(8), intent(in) :: zin(3)
        real(8), intent(out) :: zout(3)
c       
c       this entry applies the mapping to the point zin \in R^3,
c       returning zout.
c       
        call tri_rotf(zin,w(1:9),w(10:12),zout)
        return
        end subroutine tri_for
c       
c       
c       
c***********************************************************************
        subroutine tri_bak(w,zin,zout)
c***********************************************************************
        implicit none

        real(8), intent(in) :: w(12)
        real(8), intent(in) :: zin(3)
        real(8), intent(out) :: zout(3)
c       
c       this entry applies the inverse mapping to the point zin in R^3,
c       returning the pont zout in R^3.
c      
        call tri_rotb(zin,w(1:9),w(10:12),zout)
c       
        return
        end subroutine tri_bak
c       
c       
c       
c
c       
c***********************************************************************
      subroutine tri_rotf(z,u,shift,y)     
c***********************************************************************
      implicit none

      real(8), intent(in) :: z(3)
      real(8), intent(in) :: u(3,3)
      real(8), intent(in) :: shift(3)
      real(8), intent(out) :: y(3)

      real(8) :: temp(3)
c
c     translate by vector shift and rotate using matrix u^T.
c
      temp(1)=z(1)+shift(1)
      temp(2)=z(2)+shift(2)
      temp(3)=z(3)+shift(3)
      y(1)=u(1,1)*temp(1)+u(2,1)*temp(2)+u(3,1)*temp(3)
      y(2)=u(1,2)*temp(1)+u(2,2)*temp(2)+u(3,2)*temp(3)
      y(3)=u(1,3)*temp(1)+u(2,3)*temp(2)+u(3,3)*temp(3)
      return
      end
c       
c       
c       
c       
c       
c***********************************************************************
      subroutine tri_rotb(z,u,shift,y)     
c***********************************************************************
        implicit none

        real(8), intent(in) :: z(3)
        real(8), intent(in) :: u(3,3)
        real(8), intent(in) :: shift(3)
        real(8), intent(out) :: y(3)

        real(8) :: temp(3)
c
c     rotate using matrix u and translate by vector -shift.
c
      temp(1)=u(1,1)*z(1)+u(1,2)*z(2)+u(1,3)*z(3)
      temp(2)=u(2,1)*z(1)+u(2,2)*z(2)+u(2,3)*z(3)
      temp(3)=u(3,1)*z(1)+u(3,2)*z(2)+u(3,3)*z(3)
      y(1)=temp(1)-shift(1)
      y(2)=temp(2)-shift(2)
      y(3)=temp(3)-shift(3)
      return
      end
c       
c

c***********************************************************************
        subroutine tri_crotf(z,u,shift,y)     
c***********************************************************************
        implicit none

        complex(8), intent(in) :: z(3)
        real(8), intent(in) :: u(3,3)
        real(8), intent(in) :: shift(3)
        complex(8), intent(out) :: y(3)

        complex(8) :: temp(3)
c
c     translate by vector shift and rotate using matrix u^T.
c       complex(8) version
c
        temp(1)=z(1)+shift(1)
        temp(2)=z(2)+shift(2)
        temp(3)=z(3)+shift(3)
        y(1)=u(1,1)*temp(1)+u(2,1)*temp(2)+u(3,1)*temp(3)
        y(2)=u(1,2)*temp(1)+u(2,2)*temp(2)+u(3,2)*temp(3)
        y(3)=u(1,3)*temp(1)+u(2,3)*temp(2)+u(3,3)*temp(3)
        return
        end
c       
c       
c       
c       
c       
c***********************************************************************
        subroutine tri_crotb(z,u,shift,y)     
c***********************************************************************
        implicit none

        complex(8), intent(in) :: z(3)
        real(8), intent(in) :: u(3,3)
        real(8), intent(in) :: shift(3)
        complex(8), intent(out) :: y(3)

        complex(8) :: temp(3)
c
c     rotate using matrix u and translate by vector -shift.
c       complex(8) version
c       
        temp(1)=u(1,1)*z(1)+u(1,2)*z(2)+u(1,3)*z(3)
        temp(2)=u(2,1)*z(1)+u(2,2)*z(2)+u(2,3)*z(3)
        temp(3)=u(3,1)*z(1)+u(3,2)*z(2)+u(3,3)*z(3)
        y(1)=temp(1)-shift(1)
        y(2)=temp(2)-shift(2)
        y(3)=temp(3)-shift(3)
        return
        end
c       
c
c
c***********************************************************************
c
c                    Vector rotation routines
c
c***********************************************************************
c
c
c
        subroutine tri_for_vect(w,zin,zout)
        implicit none

        real(8), intent(in) :: w(9)
        real(8), intent(in) :: zin(3)
        real(8), intent(out) :: zout(3)

        real(8) :: shift(3)
c       
c       this entry applies the mapping to the vector zin \in R^3,
c       returning zout.
c       
        shift(1)=0
        shift(2)=0
        shift(3)=0
        call tri_rotf(zin,w(1),shift,zout)
        return
        end
c       
c       
c       
        subroutine tri_bak_vect(w,zin,zout)
        implicit none

        real(8), intent(in) :: w(9)
        real(8), intent(in) :: zin(3)
        real(8), intent(out) :: zout(3)

        real(8) :: shift(3)
c       
c       this entry applies the mapping to the vector zin \in R^3,
c       returning zout.
c       
        shift(1)=0
        shift(2)=0
        shift(3)=0
        call tri_rotb(zin,w(1),shift,zout)
        return
        end
c       
c       
        subroutine tri_cfor_vect(w,zin,zout)
        implicit none

        real(8), intent(in) :: w(9)
        complex(8), intent(in) :: zin(3)
        complex(8), intent(out) :: zout(3)

        real(8) :: shift(3)
c       
c       this entry applies the mapping to the complex vector zin \in R^3,
c       returning zout. complex(8) version
c       
        shift(1)=0
        shift(2)=0
        shift(3)=0
        call tri_crotf(zin,w(1),shift,zout)
        return
        end
c       
c       
c       
        subroutine tri_cbak_vect(w,zin,zout)
        implicit none

        real(8), intent(in) :: w(9)
        complex(8), intent(in) :: zin(3)
        complex(8), intent(out) :: zout(3)

        real(8) :: shift(3)
c       
c       this entry applies the mapping to the complex vector zin \in R^3,
c       returning zout. complex(8) version
c       
        shift(1)=0
        shift(2)=0
        shift(3)=0
        call tri_crotb(zin,w(1),shift,zout)
        return
        end
c       
c       
c       
c***********************************************************************
c
c                    Hessian rotation routines
c
c***********************************************************************
        subroutine rothess3d(w,triangle,
     $      valxx,valyy,valzz,valxy,valxz,valyz,
     $      derxx,deryy,derzz,derxy,derxz,deryz)

        implicit none

        real(8), intent(in) :: w(12)
        real(8), intent(in) :: triangle(3,3)   ! not used
        real(8), intent(in) :: valxx
        real(8), intent(in) :: valyy
        real(8), intent(in) :: valzz
        real(8), intent(in) :: valxy
        real(8), intent(in) :: valxz
        real(8), intent(in) :: valyz
        real(8), intent(out) :: derxx
        real(8), intent(out) :: deryy
        real(8), intent(out) :: derzz
        real(8), intent(out) :: derxy
        real(8), intent(out) :: derxz
        real(8), intent(out) :: deryz

        real(8) :: hmatr(3,3)
        real(8) :: hders(3,3)
c        
        hmatr(1,1)=valxx
        hmatr(2,2)=valyy
        hmatr(3,3)=valzz
        hmatr(1,2)=valxy
        hmatr(1,3)=valxz
        hmatr(2,3)=valyz
        hmatr(2,1)=hmatr(1,2)
        hmatr(3,1)=hmatr(1,3)
        hmatr(3,2)=hmatr(2,3)
c
        call rothess3dmatr(hmatr,w(1),hders)
c
        derxx=hders(1,1)
        deryy=hders(2,2)
        derzz=hders(3,3)
        derxy=hders(1,2)
        derxz=hders(1,3)
        deryz=hders(2,3)
c
        return
        end
c
c
c
        subroutine rothess3dmatr(h,rotmat,d)
        implicit none

        real(8), intent(in) :: h(3,3)
        real(8), intent(in) :: rotmat(3,3)
        real(8), intent(out) :: d(3,3)

        real(8) :: temp(3,3)
c
        call rothess3drot1(h,rotmat,temp)
        call rothess3drot2(temp,rotmat,d)
c
        return
        end
c
c
c
        subroutine rothess3drot1(h,u,d)
        implicit none

        real(8), intent(in) :: h(3,3)
        real(8), intent(in) :: u(3,3)
        real(8), intent(out) :: d(3,3)

        integer :: k
c
        do k=1,3
        d(1,k)=u(1,1)*h(1,k)+u(1,2)*h(2,k)+u(1,3)*h(3,k)
        d(2,k)=u(2,1)*h(1,k)+u(2,2)*h(2,k)+u(2,3)*h(3,k)
        d(3,k)=u(3,1)*h(1,k)+u(3,2)*h(2,k)+u(3,3)*h(3,k)
        enddo
c
        return
        end
c
c
c
        subroutine rothess3drot2(h,u,d)
        implicit none

        real(8), intent(in) :: h(3,3)
        real(8), intent(in) :: u(3,3)
        real(8), intent(out) :: d(3,3)

        integer :: k
c
        do k=1,3
        d(1,k)=u(1,1)*h(k,1)+u(1,2)*h(k,2)+u(1,3)*h(k,3)
        d(2,k)=u(2,1)*h(k,1)+u(2,2)*h(k,2)+u(2,3)*h(k,3)
        d(3,k)=u(3,1)*h(k,1)+u(3,2)*h(k,2)+u(3,3)*h(k,3)
        enddo

c        do k=1,3
c        d(1,k)=u(1,1)*h(1,k)+u(2,1)*h(2,k)+u(3,1)*h(3,k)
c        d(2,k)=u(1,2)*h(1,k)+u(2,2)*h(2,k)+u(3,2)*h(3,k)
c        d(3,k)=u(1,3)*h(1,k)+u(2,3)*h(2,k)+u(3,3)*h(3,k)
c        enddo
c
        return
        end
c
c
c
c
c
