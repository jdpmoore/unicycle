!-----------------------------------------------------------------------
! Copyright 2017 Sylvain Barbot
!
! This file is part of UNICYCLE
!
! UNICYCLE is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! UNICYCLE is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with UNICYCLE.  If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------

MODULE exportnetcdf

  USE netcdf

  IMPLICIT NONE

  PUBLIC

  INTEGER, PRIVATE, PARAMETER :: NDIMS = 2

  CHARACTER (LEN = *), PRIVATE, PARAMETER :: ACTUAL_RANGE = "actual_range"

CONTAINS

  !------------------------------------------------------------------------------
  ! subroutine writeNetcdf
  ! writes a netcdf file compatible with the GMT .grd format
  ! use 0 for gridline registration and 1 for pixel node registration
  !------------------------------------------------------------------------------
  SUBROUTINE writeNetcdf(filename,nx,x,ny,y,z,registration)
    CHARACTER(LEN=256), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: nx,ny
    REAL*8, INTENT(IN) :: x(nx),y(ny),z(nx,ny)
    INTEGER, INTENT(IN) :: registration

    CHARACTER (LEN = *), PARAMETER :: X_NAME = "x"
    CHARACTER (LEN = *), PARAMETER :: Y_NAME = "y"
    CHARACTER (LEN = *), PARAMETER :: Z_NAME = "z"

    CHARACTER (LEN = *), PARAMETER :: LONG_NAME =    "long_name"
    CHARACTER (LEN = *), PARAMETER :: NODE_OFFSET =    "node_offset"
    CHARACTER (LEN = *), PARAMETER :: TITLE = "title"
    CHARACTER (LEN = *), PARAMETER :: HISTORY = "Greens function"
    CHARACTER (LEN = *), PARAMETER :: DESCRIPTION = "unicycle"
    CHARACTER (LEN = *), PARAMETER :: CONVENTIONS = "Conventions"
    CHARACTER (LEN = *), PARAMETER :: COARDS = "COARDS/CF-1.0"

    INTEGER :: ncid
    INTEGER :: x_dimid, y_dimid
    INTEGER :: x_varid, y_varid, z_varid
    INTEGER :: dimids(NDIMS)

    ! create the netcdf file
    CALL check(nf90_create(filename,nf90_clobber,ncid))

    ! dimensions
    CALL check(nf90_def_dim(ncid,X_NAME,nx,x_dimid))
    CALL check(nf90_def_dim(ncid,Y_NAME,ny,y_dimid))

    ! define the coordinate variables
    CALL check(nf90_def_var(ncid,X_NAME,NF90_DOUBLE,x_dimid,x_varid))
    CALL check(nf90_def_var(ncid,Y_NAME,NF90_DOUBLE,y_dimid,y_varid))

    dimids = (/ x_dimid, y_dimid /)
  
    ! prepare data variable
    CALL check(nf90_def_var(ncid,Z_NAME,NF90_REAL,dimids,z_varid))

    ! attributes
    CALL check(nf90_put_att(ncid,x_varid,LONG_NAME,X_NAME))
    CALL check(nf90_put_att(ncid,y_varid,LONG_NAME,Y_NAME))
    CALL check(nf90_put_att(ncid,z_varid,LONG_NAME,Z_NAME))

    IF (0 .EQ. registration) THEN
       CALL check(nf90_put_att(ncid,x_varid,ACTUAL_RANGE,(/ x(1), x(nx) /)))
       CALL check(nf90_put_att(ncid,y_varid,ACTUAL_RANGE,(/ y(1), y(ny) /)))
    ELSE
       CALL check(nf90_put_att(ncid,x_varid,ACTUAL_RANGE,(/ x(1)-0.5, x(nx)+0.5 /)))
       CALL check(nf90_put_att(ncid,y_varid,ACTUAL_RANGE,(/ y(1)-0.5, y(ny)+0.5 /)))
       CALL check(nf90_put_att(ncid,NF90_GLOBAL,NODE_OFFSET,1))
    END IF
    CALL check(nf90_put_att(ncid,z_varid,ACTUAL_RANGE,(/ MINVAL(z), MAXVAL(z) /)))

    ! global attributes
    CALL check(nf90_put_att(ncid,NF90_GLOBAL,TITLE,TITLE))
    CALL check(nf90_put_att(ncid,NF90_GLOBAL,HISTORY,DESCRIPTION))
    CALL check(nf90_put_att(ncid,NF90_GLOBAL,CONVENTIONS,COARDS))

    ! end define mode
    CALL check(nf90_enddef(ncid))
 
    ! write coordinate system
    CALL check(nf90_put_var(ncid,x_varid,x))
    CALL check(nf90_put_var(ncid,y_varid,y))

    ! write the data
    CALL check(nf90_put_var(ncid,z_varid,z))

    ! flush
    CALL check(nf90_sync(ncid))

    ! close file
    CALL check(nf90_close(ncid))

  END SUBROUTINE writeNetcdf

  !------------------------------------------------------------------------------
  ! subroutine openNetcdfUnlimited
  ! open a netcdf file and prepare to write time series
  !------------------------------------------------------------------------------
  SUBROUTINE openNetcdfUnlimited(filename,n,x,ncid,y_varid,z_varid)
    CHARACTER(LEN=256), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: n
    REAL*8, INTENT(IN) :: x(n)
    INTEGER, INTENT(OUT) :: ncid,y_varid,z_varid

    CHARACTER (LEN = *), PARAMETER :: X_NAME = "x"
    CHARACTER (LEN = *), PARAMETER :: Y_NAME = "y"
    CHARACTER (LEN = *), PARAMETER :: Z_NAME = "z"

    CHARACTER (LEN = *), PARAMETER :: LONG_NAME =    "long_name"
    CHARACTER (LEN = *), PARAMETER :: TITLE = "title"
    CHARACTER (LEN = *), PARAMETER :: HISTORY = "plane strain simulation"
    CHARACTER (LEN = *), PARAMETER :: DESCRIPTION = "unicycle"

    INTEGER :: x_dimid, y_dimid
    INTEGER :: x_varid
    INTEGER :: dimids(NDIMS)

    ! create the netcdf file
    CALL check(nf90_create(filename,nf90_clobber,ncid))

    ! time dimension has unlimited length
    CALL check(nf90_def_dim(ncid,X_NAME,n,             x_dimid))
    CALL check(nf90_def_dim(ncid,Y_NAME,NF90_UNLIMITED,y_dimid))

    ! define the coordinate variables
    CALL check(nf90_def_var(ncid,X_NAME,NF90_DOUBLE,x_dimid,x_varid))
    CALL check(nf90_def_var(ncid,Y_NAME,NF90_DOUBLE,y_dimid,y_varid))

    dimids = (/ x_dimid, y_dimid /)
  
    ! prepare data variable
    CALL check(nf90_def_var(ncid,Z_NAME,NF90_REAL,dimids,z_varid))

    ! attributes
    CALL check(nf90_put_att(ncid,x_varid,LONG_NAME,X_NAME))
    CALL check(nf90_put_att(ncid,y_varid,LONG_NAME,Y_NAME))
    CALL check(nf90_put_att(ncid,z_varid,LONG_NAME,Z_NAME))

    CALL check(nf90_put_att(ncid,x_varid,ACTUAL_RANGE,(/ x(1), x(n) /)))
    CALL check(nf90_put_att(ncid,y_varid,ACTUAL_RANGE,(/ 1, -1 /)))
    CALL check(nf90_put_att(ncid,z_varid,ACTUAL_RANGE,(/ 0, 0 /)))

    ! global attributes
    CALL check(nf90_put_att(ncid,NF90_GLOBAL,TITLE,TITLE))
    CALL check(nf90_put_att(ncid,NF90_GLOBAL,HISTORY,DESCRIPTION))

    ! end define mode
    CALL check(nf90_enddef(ncid))
 
    ! write coordinate system
    CALL check(nf90_put_var(ncid,x_varid,x))

    ! flush
    CALL check(nf90_sync(ncid))

  END SUBROUTINE openNetcdfUnlimited

  !------------------------------------------------------------------------------
  ! subroutine writeNetcdfUnlimited
  ! write d(n) in the variable varid of netcdf file ncid.
  !------------------------------------------------------------------------------
  SUBROUTINE writeNetcdfUnlimited(ncid,y_varid,z_varid,index,n,z)
    INTEGER, INTENT(IN) :: ncid,y_varid,z_varid
    INTEGER, INTENT(IN) :: n
    REAL*4, DIMENSION(n), INTENT(IN) :: z
    INTEGER, INTENT(IN) :: index

    INTEGER :: i,ierr
    INTEGER :: start(NDIMS), count(NDIMS)

    REAL*8, DIMENSION(:), ALLOCATABLE :: y

    count = (/ n, 1 /)
    start = (/ 1, index /)
    
    CALL check(nf90_put_var(ncid,z_varid,z,START=start,COUNT=count))

    ! flush every so often
    IF (0 .EQ. MOD(index,100)) THEN

       ALLOCATE(y(index),STAT=ierr)
       IF (0/=ierr) STOP "could not allocate the index array"

       DO i=1,index
          y(i)=REAL(i,8)
       END DO

       ! update the time index
       CALL check(nf90_put_var(ncid,y_varid,y))
       CALL check(nf90_put_att(ncid,y_varid,ACTUAL_RANGE,(/ 1, index /)))
       CALL check(nf90_sync(ncid))
    END IF

  END SUBROUTINE writeNetcdfUnlimited

  !------------------------------------------------------------------------------
  ! subroutine closeNetcdfUnlimited
  ! close the current netcdf file
  !------------------------------------------------------------------------------
  SUBROUTINE closeNetcdfUnlimited(ncid,y_varid,z_varid,count)
    INTEGER, INTENT(IN) :: ncid,y_varid,z_varid,count

    INTEGER :: i,ierr

    REAL*8, DIMENSION(:), ALLOCATABLE :: y

    ALLOCATE(y(count),STAT=ierr)
    IF (0/=ierr) STOP "could not allocate the index array"

    DO i=1,count
       y(i)=REAL(i,8)
    END DO

    ! update the time index
    CALL check(nf90_put_var(ncid,y_varid,y))
    CALL check(nf90_put_att(ncid,y_varid,ACTUAL_RANGE,(/ 1, count /)))

    CALL check(nf90_close(ncid))

  END SUBROUTINE closeNetcdfUnlimited

  !------------------------------------------------------------------------------
  ! subroutine check
  ! checks the return of nf90_ functions
  !------------------------------------------------------------------------------
  SUBROUTINE check(status)
    INTEGER, INTENT(IN) :: status
          
    IF (status /= nf90_noerr) THEN 
       PRINT *, TRIM(nf90_strerror(status))
       STOP "netcdf write stopped"
    END IF
  END SUBROUTINE check

END MODULE exportnetcdf

