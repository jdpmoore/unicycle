
!------------------------------------------------------------------
!> subroutine ExportVTK_RFaults
!! creates a .vtp file (in the VTK PolyData XML format) containing
!! the rectangular faults. The faults are characterized with a set
!! of subsegments (rectangles) each associated with a slip vector. 
!!
!! \author sylvain barbot 06/24/09 - original form
!------------------------------------------------------------------
SUBROUTINE exportvtk_rfaults(ns,x,length,width,strike,dip,patch,id)
  USE types 
  INTEGER, INTENT(IN) :: ns
  REAL*8, DIMENSION(3,ns), INTENT(IN) :: x
  REAL*8, DIMENSION(ns), INTENT(IN) :: length
  REAL*8, DIMENSION(ns), INTENT(IN) :: width
  REAL*8, DIMENSION(ns), INTENT(IN) :: strike
  REAL*8, DIMENSION(ns), INTENT(IN) :: dip
  TYPE(PATCH_ELEMENT_STRUCT), DIMENSION(ns), INTENT(IN) :: patch
  INTEGER, INTENT(IN) :: id

  INTEGER :: k
  CHARACTER :: q

  REAL*8 :: x1,x2,x3,cstrike,sstrike,cdip,sdip,L,W
         
  REAL*8, DIMENSION(3) :: s,d,n

  ! double-quote character
  q=char(34)

  WRITE (id,'("<?xml version=",a,"1.0",a,"?>")') q,q
  WRITE (id,'("<VTKFile type=",a,"PolyData",a," version=",a,"0.1",a,">")') q,q,q,q
  WRITE (id,'("  <PolyData>")')

  DO k=1,ns

     ! fault dimension
     L=length(k)
     W=width(k)

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

     WRITE (id,'("    <Piece NumberOfPoints=",a,"4",a," NumberOfPolys=",a,"1",a,">")') q,q,q,q
     WRITE (id,'("      <Points>")')
     WRITE (id,'("        <DataArray type=",a,"Float32",a, &
                         & " Name=",a,"Fault Patch",a, &
                         & " NumberOfComponents=",a,"3",a, &
                         & " format=",a,"ascii",a,">")') q,q,q,q,q,q,q,q

     ! fault edge coordinates
     WRITE (id,'(12ES11.2)') &
                   x1              ,x2              ,x3, &
                   x1       +s(1)*L,x2       +s(2)*L,x3       +s(3)*L, &
                   x1-d(1)*W+s(1)*L,x2-d(2)*W+s(2)*L,x3-d(3)*W+s(3)*L, &
                   x1-d(1)*W       ,x2-d(2)*W       ,x3-d(3)*W

     WRITE (id,'("        </DataArray>")')
     WRITE (id,'("      </Points>")')
     WRITE (id,'("      <Polys>")')
     WRITE (id,'("        <DataArray type=",a,"Int32",a, &
                            & " Name=",a,"connectivity",a, &
                            & " format=",a,"ascii",a, &
                            & " RangeMin=",a,"0",a, &
                            & " RangeMax=",a,"3",a,">")') q,q,q,q,q,q,q,q,q,q
     WRITE (id,'("0 1 2 3")')
     WRITE (id,'("        </DataArray>")')
     WRITE (id,'("        <DataArray type=",a,"Int32",a, &
                                 & " Name=",a,"offsets",a, &
                                 & " format=",a,"ascii",a, &
                                 & " RangeMin=",a,"4",a, &
                                 & " RangeMax=",a,"4",a,">")') q,q,q,q,q,q,q,q,q,q
     WRITE (id,'("          4")')
     WRITE (id,'("        </DataArray>")')
     WRITE (id,'("      </Polys>")')

     WRITE (15,'("      <CellData Normals=",a,"Normals",a,">")') q,q

     WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                               & " Name=",a,"a",a, &
                               & " NumberOfComponents=",a,"1",a, &
                               & " format=",a,"ascii",a,">")') q,q,q,q,q,q,q,q
     WRITE (15,'(1ES11.2)') patch(k)%a
     WRITE (15,'("        </DataArray>")')

     WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                               & " Name=",a,"b",a, &
                               & " NumberOfComponents=",a,"1",a, &
                               & " format=",a,"ascii",a,">")') q,q,q,q,q,q,q,q
     WRITE (15,'(1ES11.2)') patch(k)%b
     WRITE (15,'("        </DataArray>")')

     WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                               & " Name=",a,"Vpl",a, &
                               & " NumberOfComponents=",a,"1",a, &
                               & " format=",a,"ascii",a,">")') q,q,q,q,q,q,q,q
     WRITE (15,'(1ES11.2)') patch(k)%Vpl
     WRITE (15,'("        </DataArray>")')

     WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                               & " Name=",a,"rake",a, &
                               & " NumberOfComponents=",a,"1",a, &
                               & " format=",a,"ascii",a,">")') q,q,q,q,q,q,q,q
     WRITE (15,'(1ES11.2)') patch(k)%rake
     WRITE (15,'("        </DataArray>")')

     WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                               & " Name=",a,"mu0",a, &
                               & " NumberOfComponents=",a,"1",a, &
                               & " format=",a,"ascii",a,">")') q,q,q,q,q,q,q,q
     WRITE (15,'(1ES11.2)') patch(k)%mu0
     WRITE (15,'("        </DataArray>")')

     WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                               & " Name=",a,"sig",a, &
                               & " NumberOfComponents=",a,"1",a, &
                               & " format=",a,"ascii",a,">")') q,q,q,q,q,q,q,q
     WRITE (15,'(1ES11.2)') patch(k)%sig
     WRITE (15,'("        </DataArray>")')

     WRITE (15,'("        <DataArray type=",a,"Float32",a, &
                               & " Name=",a,"L",a, &
                               & " NumberOfComponents=",a,"1",a, &
                               & " format=",a,"ascii",a,">")') q,q,q,q,q,q,q,q
     WRITE (15,'(1ES11.2)') patch(k)%L
     WRITE (15,'("        </DataArray>")')

     WRITE (15,'("      </CellData>")')

     WRITE (id,'("    </Piece>")')

  END DO

  WRITE (id,'("  </PolyData>")')
  WRITE (id,'("</VTKFile>")')

END SUBROUTINE exportvtk_rfaults

