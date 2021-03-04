!//////////////////////////////////////////////////////////////////////
!//
!// test_tgf_main.f: Main function for testing triangle Green's functions.
!//
!//////////////////////////////////////////////////////////////////////





      program tgf_01
      
      implicit none

      integer :: ioerr




c     Initialize simple printing routines.
c     Calling PRINI(6,13) causes printing to screen and Fortran unit 13.     

      call prini(6,13)
      



c     Test the half-space triangle Green's function, by comparing its
c     results to direct numerical integration of the point dislocation.

      open(unit=13, status='replace', action='write',
     &     access='sequential', form='formatted',
     &     iostat=ioerr, file='test_tgf_log_h.txt')

      open(unit=15, status='replace', action='write',
     &     access='sequential', form='formatted',
     &     iostat=ioerr, file='test_tgf_out_h.txt')

      call test_dc3_elh_grid()

      close(unit=15)

      close(unit=13)
      

      

c     Test the full-space triangle Green's function, by comparing its
c     results to direct numerical integration of the point dislocation.

      open(unit=13, status='replace', action='write',
     &     access='sequential', form='formatted',
     &     iostat=ioerr, file='test_tgf_log_t.txt')

      open(unit=15, status='replace', action='write',
     &     access='sequential', form='formatted',
     &     iostat=ioerr, file='test_tgf_out_t.txt')

      call test_dc3_elt_grid()

      close(unit=15)

      close(unit=13)
      

      

c     Test the half-space point dislocation formula by comparing its
c     results to the Okada formula for a point dislocation.
c     Also, test the numerical integration of the point dislocation formula
c     by comparing its results to the Okada formula for a rectangular dislocation.

      open(unit=13, status='replace', action='write',
     &     access='sequential', form='formatted',
     &     iostat=ioerr, file='test_tgf_log_o.txt')

      open(unit=15, status='replace', action='write',
     &     access='sequential', form='formatted',
     &     iostat=ioerr, file='test_tgf_out_o.txt')

      call test_dc3_elh()

      close(unit=15)

      close(unit=13)



      
      stop      
      end program tgf_01








 