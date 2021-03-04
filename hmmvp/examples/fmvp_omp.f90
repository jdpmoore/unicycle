! An example of how to use Hmat.hpp in Fortran code by using the SFHmat.h
! interface.
!   Run it as
!     ./fmvp_omp H-matrix-filename number-omp-threads repeat
! where
!     H-matrix-filename
! is the name of an H-matrix file (*.hm) produced by dc3dm or hmmvp,
!     number-omp-threads
! is the number of OpenMP threads to use in the test, and
!     repeat
! is the number of times to repeat the matrix-vector product to get good timing
! results.

module util
contains
  function elapsed_time(c1, c2) result (et)
    ! c1 and c2 are 'count' from system_clock.
    real(8) et
    integer, intent(in) :: c1, c2
    integer :: c3, clock_rate

    call system_clock(c3, clock_rate)
    et = real(c2 - c1) / real(clock_rate)
  end function elapsed_time
end module util

program main
  use util
  implicit none

  real(8), dimension(:), allocatable :: x, y
  character(256) :: arg, filename
  integer :: repeat, c1, c2, i
  ! The integers that are passed to SFHmat must agree in size.
  integer(8) :: m, n, nthreads, ncol

  ! Process command-line arguments.
  if (iargc() /= 3) then
     write(*,*) 'Requires 2 args: H-matrix-filename number-omp-threads repeat'
     call exit(-1)
  end if

  call getarg(1, filename)
  call getarg(2, arg)
  read(arg, '(i10)'), nthreads
  call getarg(3, arg)
  read(arg, '(i10)'), repeat

  ! Initialize the H-matrix.
  ncol = 1
  write(*,*) 'nthreads: ', nthreads
  call sfhmat_init(trim(filename), ncol, nthreads)
  
  ! Make x and y in y = A*x.
  call sfhmat_get_size(m, n)
  write(*,*) 'm: ', m
  write(*,*) 'n: ', n

  allocate (x(n), y(m))
  do i = 1, n
     x(i) = 1.0 / i
  end do

  ! Do the matrix-vector product repeatedly while timing it.
  call system_clock(c1)
  do i = 1, repeat
     call sfhmat_mvp(x, y, ncol)
  end do
  call system_clock(c2)
  write(*,*) 'average time per MVP: ', elapsed_time(c1, c2) / repeat

  ! Deallocate.
  call sfhmat_cleanup()
  deallocate (x, y)

end program main
