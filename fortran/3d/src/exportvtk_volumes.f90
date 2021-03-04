
!------------------------------------------------------------------
!> subroutine ExportVTK_Volumes
!! creates a .vtp file (in the VTK PolyData XML format) containing
!! the strain volumes.
!!
!! \author sylvain barbot 06/24/09 - original form
!------------------------------------------------------------------
SUBROUTINE exportvtk_volumes(ns,x,length,thickness,width,strike,dip,id)
  USE types 
  INTEGER, INTENT(IN) :: ns
  REAL*8, DIMENSION(3,ns), INTENT(IN) :: x
  REAL*8, DIMENSION(ns), INTENT(IN) :: length
  REAL*8, DIMENSION(ns), INTENT(IN) :: thickness
  REAL*8, DIMENSION(ns), INTENT(IN) :: width
  REAL*8, DIMENSION(ns), INTENT(IN) :: strike
  REAL*8, DIMENSION(ns), INTENT(IN) :: dip
  INTEGER, INTENT(IN) :: id

  INTEGER :: k
  CHARACTER :: q

  REAL*8 :: x1,x2,x3,cstrike,sstrike,cdip,sdip,L,T,W
         
  REAL*8, DIMENSION(3) :: s,d,n

  ! double-quote character
  q=char(34)

  WRITE (id,'("<?xml version=",a,"1.0",a,"?>")') q,q
  WRITE (id,'("<VTKFile type=",a,"PolyData",a," version=",a,"0.1",a,">")') q,q,q,q
  WRITE (id,'("  <PolyData>")')

  DO k=1,ns

     ! volume dimension
     L=length(k)
     W=width(k)
     T=thickness(k)

     ! fault top-corner position
     x1=x(1,k)
     x2=x(2,k)
     x3=x(3,k)

     cstrike=COS(strike(k))
     sstrike=SIN(strike(k))
     cdip=COS(dip(k))
     sdip=SIN(dip(k))
 
     ! strike-slip unit direction
     s(1)=cstrike
     s(2)=sstrike
     s(3)=0._8

     ! dip-slip unit direction
     d(1)=+sstrike*cdip
     d(2)=-cstrike*cdip
     d(3)=-sdip

     ! surface normal vector components
     n(1)=-sdip*sstrike
     n(2)=+sdip*cstrike
     n(3)=-cdip

     WRITE (id,'("    <Piece NumberOfPoints=",a,"8",a," NumberOfPolys=",a,"1",a,">")') q,q,q,q
     WRITE (id,'("      <Points>")')
     WRITE (id,'("        <DataArray type=",a,"Float32",a, &
                         & " Name=",a,"Fault Patch",a, &
                         & " NumberOfComponents=",a,"3",a, &
                         & " format=",a,"ascii",a,">")') q,q,q,q,q,q,q,q

     ! fault edge coordinates
     WRITE (id,'(24ES11.2)') &
                   x1              -n(1)*T/2.d0,x2              -n(2)*T/2.d0,x3              -n(3)*T/2.d0, &
                   x1       +s(1)*L-n(1)*T/2.d0,x2       +s(2)*L-n(2)*T/2.d0,x3       +s(3)*L-n(3)*T/2.d0, &
                   x1-d(1)*W+s(1)*L-n(1)*T/2.d0,x2-d(2)*W+s(2)*L-n(2)*T/2.d0,x3-d(3)*W+s(3)*L-n(3)*T/2.d0, &
                   x1-d(1)*W       -n(1)*T/2.d0,x2-d(2)*W       -n(2)*T/2.d0,x3-d(3)*W       -n(3)*T/2.d0, &
                   x1-d(1)*W       +n(1)*T/2.d0,x2-d(2)*W       +n(2)*T/2.d0,x3-d(3)*W       +n(3)*T/2.d0, &
                   x1              +n(1)*T/2.d0,x2              +n(2)*T/2.d0,x3              +n(3)*T/2.d0, &
                   x1       +s(1)*L+n(1)*T/2.d0,x2       +s(2)*L+n(2)*T/2.d0,x3       +s(3)*L+n(3)*T/2.d0, &
                   x1-d(1)*W+s(1)*L+n(1)*T/2.d0,x2-d(2)*W+s(2)*L+n(2)*T/2.d0,x3-d(3)*W+s(3)*L+n(3)*T/2.d0

     WRITE (id,'("        </DataArray>")')
     WRITE (id,'("      </Points>")')
     WRITE (id,'("      <Polys>")')
     WRITE (id,'("        <DataArray type=",a,"Int32",a, &
                            & " Name=",a,"connectivity",a, &
                            & " format=",a,"ascii",a, &
                            & " RangeMin=",a,"0",a, &
                            & " RangeMax=",a,"7",a,">")') q,q,q,q,q,q,q,q,q,q
     WRITE (id,'("7 4 5 6 7 4 3 2 7 2 1 6")')
     WRITE (id,'("        </DataArray>")')
     WRITE (id,'("        <DataArray type=",a,"Int32",a, &
                                 & " Name=",a,"offsets",a, &
                                 & " format=",a,"ascii",a, &
                                 & " RangeMin=",a,"12",a, &
                                 & " RangeMax=",a,"12",a,">")') q,q,q,q,q,q,q,q,q,q
     WRITE (id,'("          12")')
     WRITE (id,'("        </DataArray>")')
     WRITE (id,'("      </Polys>")')

     WRITE (id,'("    </Piece>")')

     WRITE (id,'("    <Piece NumberOfPoints=",a,"8",a," NumberOfPolys=",a,"1",a,">")') q,q,q,q
     WRITE (id,'("      <Points>")')
     WRITE (id,'("        <DataArray type=",a,"Float32",a, &
                        & " Name=",a,"Weak Zone",a, &
                        & " NumberOfComponents=",a,"3",a, &
                        & " format=",a,"ascii",a,">")') q,q,q,q,q,q,q,q

     ! edge coordinates
     WRITE (id,'(24ES11.2)') &
                   x1              -n(1)*T/2.d0,x2              -n(2)*T/2.d0,x3              -n(3)*T/2.d0, &
                   x1       +s(1)*L-n(1)*T/2.d0,x2       +s(2)*L-n(2)*T/2.d0,x3       +s(3)*L-n(3)*T/2.d0, &
                   x1-d(1)*W+s(1)*L-n(1)*T/2.d0,x2-d(2)*W+s(2)*L-n(2)*T/2.d0,x3-d(3)*W+s(3)*L-n(3)*T/2.d0, &
                   x1-d(1)*W       -n(1)*T/2.d0,x2-d(2)*W       -n(2)*T/2.d0,x3-d(3)*W       -n(3)*T/2.d0, &
                   x1-d(1)*W       +n(1)*T/2.d0,x2-d(2)*W       +n(2)*T/2.d0,x3-d(3)*W       +n(3)*T/2.d0, &
                   x1              +n(1)*T/2.d0,x2              +n(2)*T/2.d0,x3              +n(3)*T/2.d0, &
                   x1       +s(1)*L+n(1)*T/2.d0,x2       +s(2)*L+n(2)*T/2.d0,x3       +s(3)*L+n(3)*T/2.d0, &
                   x1-d(1)*W+s(1)*L+n(1)*T/2.d0,x2-d(2)*W+s(2)*L+n(2)*T/2.d0,x3-d(3)*W+s(3)*L+n(3)*T/2.d0

     WRITE (id,'("        </DataArray>")')
     WRITE (id,'("      </Points>")')
     WRITE (id,'("      <Polys>")')
     WRITE (id,'("        <DataArray type=",a,"Int32",a, &
                         & " Name=",a,"connectivity",a, &
                         & " format=",a,"ascii",a, &
                         & " RangeMin=",a,"0",a, &
                         & " RangeMax=",a,"7",a,">")') q,q,q,q,q,q,q,q,q,q
     WRITE (id,'("0 1 2 3 0 5 4 3 0 1 6 5")')
     WRITE (id,'("        </DataArray>")')
     WRITE (id,'("        <DataArray type=",a,"Int32",a, &
                              & " Name=",a,"offsets",a, &
                              & " format=",a,"ascii",a, &
                              & " RangeMin=",a,"12",a, &
                              & " RangeMax=",a,"12",a,">")') q,q,q,q,q,q,q,q,q,q
     WRITE (id,'("          12")')
     WRITE (id,'("        </DataArray>")')
     WRITE (id,'("      </Polys>")')
     WRITE (id,'("    </Piece>")')
  END DO

  WRITE (id,'("  </PolyData>")')
  WRITE (id,'("</VTKFile>")')

END SUBROUTINE exportvtk_volumes

