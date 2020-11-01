Triangle Green's Functions for Half-Space and Full-Space

----------------------------------------------------------------------

FILES

Files are supplied in three directories:

  source -- Source files for evaluating layer potentials 
            (integrating the Green's function over triangles), 
            for both half-space and full-space.

  test_driver -- Example code which illustrates how to call the evaluation
                 routines, and performs a series of tests to demonstrate
                 their accuracy.

  docs -- Documentation (this file).

----------------------------------------------------------------------

COMPILING THE CODE

The code is written in Fortran 95, and should compile using any recent
Fortran compiler.  

We recommend that you compile for a 64-bit machine architecture, and that you
set compiler options so that the Fortran INTEGER data type is 64 bits.

Source file mod_tgf.f must be compiled first, because it contains Fortran
modules that are referenced by other files.  The remaining files can be
compiled in any order.

The code supports OpenMP.  If you set compiler options to enable OpenMP, then
the code uses multiple CPU cores when you compute Green's function values
at multiple target locations.

----------------------------------------------------------------------

PRECISION CONTROL

Precision is controlled by two parameters which are defined in mod_tgf.f.

For 6 digits of precision, use the following settings.  This is adequate
for ordinary usage.

      integer, parameter :: elhtriaadap_nq = 6
      real(8), parameter :: elhtriaadap_eps = 1e-6

For 12 digits of precision, use the following settings.  This is useful
for testing the accuracy of the Green's functions.

      integer, parameter :: elhtriaadap_nq = 20
      real(8), parameter :: elhtriaadap_eps = 1e-12

----------------------------------------------------------------------

PROBLEM DESCRIPTION

The problem solved by the Green's function code is the following.  The setting
is an elastic medium, which can be either a half-space or a full-space.  We
use (x,y,z) coordinates.  In the case of a half-space, the elastic medium
occupies the half-space z<=0, and the free surface is z==0.  Coordinates are
given in units of meters.

The elastic medium is characterized by its two Lame moduli:  the volumetric
modulus lamdba, and the shear modulus mu.  The Lame moduli are given in units
of Pascal.

We are given a set of triangles within the elastic medium.  For each triangle,
we are given a slip vector, which is defined to be the displacement on the
"positive" side of the triangle minus the displacement on the "negative" side
of the triangle.  The slip vector is assumed to be uniform on each triangle
(hence the slip is discontinuous at the triangle edges).  The slip vector is
not required to be tangent to the triangle, that is, we allow opening.

Optionally, we are also given a set of observation points within the elastic
medium, which are called targets.

A target is specified by giving its (x,y,z) coordinates.  A triangle is
specified by giving the (x,y,z) coordinates of each of its three vertices.

For the half-space case, all coordinates must satisfy z<=0.  If you wish to
place a vertex or target exactly on the free surface z==0, we recommend that
you use a small negative value such as z=-1.0e-14, rather than setting z
exactly equal to zero.

Given this setup, the code can calculate four things:

* The displacement at the centroid of each triangle.

* The strain tensor at the centroid of each triangle.

* The displacement at each target location.

* The strain tensor at each target location.

When you call the Green's function code, you supply four flags that specify
which of these you wish to compute.

A target may not be located within a triangle (nor within a triangle edge,
nor coincident with a triangle vertex).

----------------------------------------------------------------------

CALLING THE GREEN'S FUNCTION CODE

For the half-space case, the following subroutine evaluates layer potentials
using the Green's function:

        subroutine elh3dtriadirecttarg(
     $     rlam,rmu,triangle,trinorm,nsource,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntarget,target,ifptfrctarg,ptfrctarg,
     $     ifstraintarg,straintarg)

For the full-space case, the following subroutine evaluates layer potentials
using the Green's function:

        subroutine el3dtriadirecttarg(
     $     rlam,rmu,triangle,trinorm,nsource,source,
     $     ifsingle,sigma_sl,ifdouble,sigma_dl,
     $     ifptfrc,ptfrc,ifstrain,strain,
     $     ntarget,target,ifptfrctarg,ptfrctarg,
     $     ifstraintarg,straintarg)

Both subroutines use the same set of arguments, which are declared as:

        real(8), intent(in) :: rlam
        real(8), intent(in) :: rmu
        integer, intent(in) :: nsource
        real(8), intent(in) :: triangle(3,3,nsource)
        real(8), intent(in) :: trinorm(3,nsource)
        real(8), intent(in) :: source(3,nsource)
        integer, intent(in) :: ifsingle
        real(8), intent(in) :: sigma_sl(3,ifdim(nsource,ifsingle))
        integer, intent(in) :: ifdouble
        real(8), intent(in) :: sigma_dl(3,ifdim(nsource,ifdouble))
        integer, intent(in) :: ifptfrc
        real(8), intent(inout) :: ptfrc(3,ifdim(nsource,ifptfrc))
        integer, intent(in) :: ifstrain
        real(8), intent(inout) :: strain(3,3,ifdim(nsource,ifstrain))
        integer, intent(in) :: ntarget
        real(8), intent(in) :: target(3,ntarget)
        integer, intent(in) :: ifptfrctarg
        real(8), intent(inout)
     &         :: ptfrctarg(3,ifdim(ntarget,ifptfrctarg))
        integer, intent(in) :: ifstraintarg
        real(8), intent(inout)
     &         :: straintarg(3,3,ifdim(ntarget,ifstraintarg))

The function ifdim is defined as follows.  It is used to indicate whether
or not an array is used, depending on the setting of a flag.

      pure integer function ifdim (n, ifn)
      implicit none
      integer, intent(in) :: n
      integer, intent(in) :: ifn

      if (ifn /= 0) then
        ifdim = n
      else
        ifdim = 0
      endif

      end function ifdim

We now describe each of the arguments to the subroutines.


rlam

Lame volumetric modulus lambda.


rmu

Lame shear modulus mu.


nsource

The number of triangles.  It must be a positive integer.


triangle(3,3,nsource)

The vertices of each triangle.  It is defined as:

  triangle(1,j,n) == x coordinate of the j-th vertex of the n-th triangle.
  triangle(2,j,n) == y coordinate of the j-th vertex of the n-th triangle.
  triangle(3,j,n) == z coordinate of the j-th vertex of the n-th triangle.


trinorm(3,nsource)

Unit normal vector of each triangle, pointing toward the "positive" side of
the triangle.  It is defined as:

  trinorm(1,n) == x component of the unit normal to the n-th triangle.
  trinorm(2,n) == y component of the unit normal to the n-th triangle.
  trinorm(3,n) == z component of the unit normal to the n-th triangle.

Important:  The "positive" side of the triangle is defined to the be the side
from which the vertices appear in counterclockwise order (assuming a right-
handed coordinate system).  To calculate the unit normal correctly, we
recommend that you use the subroutine triangle_norm, like this:

        do i = 1,nsource
        call triangle_norm(triangle(1,1,i),trinorm(1,i))
        enddo


source(3,nsource)

Centroid of each triangle.  It is defined as:

  source(1,n) == x coordinate of the centroid of the n-th triangle.
  source(2,n) == y coordinate of the centroid of the n-th triangle.
  source(3,n) == z coordinate of the centroid of the n-th triangle.

Important:  The coordinates of the centroid must equal the average of the
coordinates of the three vertices, like this:

        do i = 1,nsource
        source(1,i)=
     $     (triangle(1,1,i)+triangle(1,2,i)+triangle(1,3,i))/3
        source(2,i)=
     $     (triangle(2,1,i)+triangle(2,2,i)+triangle(2,3,i))/3
        source(3,i)=
     $     (triangle(3,1,i)+triangle(3,2,i)+triangle(3,3,i))/3
        enddo


ifsingle

A flag indicating if you are supplying a "single layer potential" for each
triangle (which can be thought of as an external force acting on the
triangle).  The value must be either 0 or 1:

  ifsingle == 0  -->  No single layer potential is supplied.
  ifsingle == 1  -->  A single layer potential is supplied.

Ordinarily, you should set ifsingle=0.


sigma_sl(3,ifdim(nsource,ifsingle))

If ifsingle==1, then sigma_sl contains the "single layer potential" for each
triangle, as follows:

  sigma_sl(1,n) == x component of the single layer potential of the
                   n-th triangle.
  sigma_sl(2,n) == y component of the single layer potential of the
                   n-th triangle.
  sigma_sl(3,n) == z component of the single layer potential of the
                   n-th triangle.

If ifsingle==0, then sigma_sl is not used.


ifdouble

A flag indicating if you are supplying a slip vector for each triangle.  The
slip vector is the displacement on the positive side of the triangle, minus
the displacement on the negative side of the triangle, and is assumed to be
uniform over the triangle.  The value must be either 0 or 1:

  ifdouble == 0  -->  No slip vector is supplied.
  ifdouble == 1  -->  A slip vector is supplied.

Ordinarily, you should set ifdouble=1.


sigma_dl(3,ifdim(nsource,ifdouble))

If ifdouble==1, then sigma_dl contains the slip vector for each triangle, as
follows:

  sigma_dl(1,n) == x component of the slip vector of the n-th triangle.
  sigma_dl(2,n) == y component of the slip vector of the n-th triangle.
  sigma_dl(3,n) == z component of the slip vector of the n-th triangle.

If ifdouble==0, then sigma_dl is not used.



ifptfrc

A flag indicating if you want to compute the displacement at each triangle
centroid.  The value must be either 0 or 1:

  ifptfrc == 0  -->  Do not compute displacments at triangle centroids.
  ifptfrc == 1  -->  Compute displacments at triangle centroids.


ptfrc(3,ifdim(nsource,ifptfrc))

If ifptfrc==1, then ptfrc is used to return the computed displacment at each
triangle centroid, as follows:

  ptfrc(1,n) == x component of displacment at the n-th triangle centroid.
  ptfrc(2,n) == y component of displacment at the n-th triangle centroid.
  ptfrc(3,n) == z component of displacment at the n-th triangle centroid.

Note:  After calling the Green's function code, you must divide all the
elements of ptfrc by 4*pi.

If ifptfrc==0, then ptfrc is not used.



ifstrain

A flag indicating if you want to compute the strain tensor at each triangle
centroid.  The value must be either 0 or 1:

  ifstrain == 0  -->  Do not compute strain tensors at triangle centroids.
  ifstrain == 1  -->  Compute strain tensors at triangle centroids.



strain(3,3,ifdim(nsource,ifstrain))

If ifstrain==1, then strain is used to return the computed strain tensor at
each triangle centroid, as follows:

  strain(1,1,n) == xx component of strain tensor at the n-th triangle centroid.
  strain(2,1,n) == xy component of strain tensor at the n-th triangle centroid.
  strain(3,1,n) == xz component of strain tensor at the n-th triangle centroid.
  strain(1,2,n) == yx component of strain tensor at the n-th triangle centroid.
  strain(2,2,n) == yy component of strain tensor at the n-th triangle centroid.
  strain(3,2,n) == yz component of strain tensor at the n-th triangle centroid.
  strain(1,3,n) == zx component of strain tensor at the n-th triangle centroid.
  strain(2,3,n) == zy component of strain tensor at the n-th triangle centroid.
  strain(3,3,n) == zz component of strain tensor at the n-th triangle centroid.

Note:  After calling the Green's function code, you must divide all the
elements of strain by 4*pi.

If ifstrain==0, then strain is not used.



ntarget

The number of targets.  It must be a non-negative integer.



target(3,ntarget)

Coordinates of each target.  It is defined as:

  target(1,n) == x coordinate of the n-th target.
  target(2,n) == y coordinate of the n-th target.
  target(3,n) == z coordinate of the n-th target.



ifptfrctarg

A flag indicating if you want to compute the displacement at each target.
The value must be either 0 or 1:

  ifptfrctarg == 0  -->  Do not compute displacments at targets.
  ifptfrctarg == 1  -->  Compute displacments at targets.



ptfrctarg(3,ifdim(nsource,ifptfrctarg))

If ifptfrctarg==1, then ptfrctarg is used to return the computed displacment
at each target, as follows:

  ptfrctarg(1,n) == x component of displacment at the n-th target.
  ptfrctarg(2,n) == y component of displacment at the n-th target.
  ptfrctarg(3,n) == z component of displacment at the n-th target.

Note:  After calling the Green's function code, you must divide all the
elements of ptfrctarg by 4*pi.

If ifptfrctarg==0, then ptfrctarg is not used.



ifstraintarg

A flag indicating if you want to compute the strain tensor at each target.
The value must be either 0 or 1:

  ifstraintarg == 0  -->  Do not compute strain tensors at targets.
  ifstraintarg == 1  -->  Compute strain tensors at targets.



straintarg(3,3,ifdim(nsource,ifstraintarg))

If ifstraintarg==1, then straintarg is used to return the computed strain
tensor at each target, as follows:

  straintarg(1,1,n) == xx component of strain tensor at the n-th target.
  straintarg(2,1,n) == xy component of strain tensor at the n-th target.
  straintarg(3,1,n) == xz component of strain tensor at the n-th target.
  straintarg(1,2,n) == yx component of strain tensor at the n-th target.
  straintarg(2,2,n) == yy component of strain tensor at the n-th target.
  straintarg(3,2,n) == yz component of strain tensor at the n-th target.
  straintarg(1,3,n) == zx component of strain tensor at the n-th target.
  straintarg(2,3,n) == zy component of strain tensor at the n-th target.
  straintarg(3,3,n) == zz component of strain tensor at the n-th target.

Note:  After calling the Green's function code, you must divide all the
elements of straintarg by 4*pi.

If ifstraintarg==0, then straintarg is not used.


