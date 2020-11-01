!//////////////////////////////////////////////////////////////////////
!//
!// mod_tgf.f: Fortran modules for triangle Green's functions.
!//
!//////////////////////////////////////////////////////////////////////








c----------------------------------------------------------------------
c
c     Modules for use in multiple source files.








c----------------------------------------------------------------------
c
c     tria3ifun_enum - Enumeration that selects function to be numerically
c     integrated by tria3adam.

      module tria3ifun_enum
      implicit none

      integer, parameter :: i_fun3eluh_eval = 1
      integer, parameter :: i_fun3elth_eval = 2

      integer, parameter :: i_fun3elt_eval = 3
      integer, parameter :: i_fun3elthb_eval = 4

      end module tria3ifun_enum








c----------------------------------------------------------------------
c
c     elhtriaadap_mod - Parameters for adaptive integration.
c     This controls the precision of the adaptive integration routines
c      in elh3dtria.f, which is used for integrating half-space image
c      contributions.
c     The first values are adequate for ordinary usage.  The second
c      higher-precision values are suitable for code testing.

      module elhtriaadap_mod
      implicit none

      ! elhtriaadap_nq - Quadrature order to use on each sub-triangle.

!!      integer, parameter :: elhtriaadap_nq = 6
      integer, parameter :: elhtriaadap_nq = 20

      ! elhtriaadap_eps - Accuracy with which integrals are evaluated.

!!      real(8), parameter :: elhtriaadap_eps = 1e-6
      real(8), parameter :: elhtriaadap_eps = 1e-12

      end module elhtriaadap_mod








c----------------------------------------------------------------------
c
c     triaq_mod - Parameters for triangle quadrature.
c     This is provided for callers of triasymq.

      module triaq_mod
      implicit none

      ! triaq_maxnodes - Maximum number of quadrature nodes.
      ! It is guaranteed that the number of nodes returned by triasymq
      !  is less than or equal to triaq_maxnodes.
      ! The arrays passed in to triasymq must have size at least as
      !  large as triaq_maxnodes.

      integer, parameter :: triaq_maxnodes = 500

      ! triaq_worksize - Workspace array size.
      ! The workspace array passed in to triasymq must have size at
      !  least as large as triaq_worksize.

      integer, parameter :: triaq_worksize = 9600

      ! triaq_maxdegree - Maximum quadrature degree.

      integer, parameter :: triaq_maxdegree = 50

      ! triaq_nodes - Number of quadrature nodes.
      ! For degree n, where n ranges from 1 to triaq_maxdegree, the
      !  number of nodes returned by triasymq is triaq_nodes(n).

      integer, parameter :: triaq_nodes(50) = (/
     &     1,    3,    6,    6,    7,   12,   15,   16,   19,   25,
     &    28,   33,   37,   42,   49,   55,   60,   67,   73,   79,
     &    87,   96,  103,  112,  120,  130,  141,  150,  159,  171,
     &   181,  193,  204,  214,  228,  243,  252,  267,  282,  295,
     &   309,  324,  339,  354,  370,  385,  399,  423,  435,  453  /)

      end module triaq_mod








c----------------------------------------------------------------------
c
c     ifdim_def - Function for dimension of conditionally-present arrays.
c
c     The function ifdim is a pure function, which is used to specify
c     the dimension of an array that is a subroutine dummy argument.
c     The presence of the array is given by the flag ifn, which is 1 if
c     the array is present, or 0 if it is not.  If the array is present,
c     its dimension is n.  If not present, we supply a dimension of 0
c     (although it might make sense to supply 1).

      module ifdim_def
      implicit none

      contains

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

      end module ifdim_def








c----------------------------------------------------------------------
c
c     falloc_thread - OpenMP control functions.

      module falloc_thread
      implicit none

      contains




c     falloc_enter_parallel - Enter a parallel region.
c     This function is called before entering an OpenMP parallel region.
c     This is useful if any code needs to be aware of when it is executing
c      in an OpenMP parallel region.

      subroutine falloc_enter_parallel()

      implicit none

      return
      end subroutine falloc_enter_parallel




c     falloc_exit_parallel - Exit a parallel region.
c     Parameters:
c      reduced_ier = Variable to receive error code.
c     This function is called after exiting an OpenMP parallel region.
c     This is useful if any code needs to be aware of when it is executing
c      in an OpenMP parallel region.

      subroutine falloc_exit_parallel(reduced_ier)

      implicit none

      integer, intent(out) :: reduced_ier

      reduced_ier = 0

      return
      end subroutine falloc_exit_parallel




      end module falloc_thread








c----------------------------------------------------------------------
c
c     elfmm_thunk - Memory allocation.
c
c     These routines are responsible for allocating working memory.
c     In this version, we use the system memory allocator.
c     They could be replaced by routines using a different allocation method.
c
c     Note that these routines must work properly within an OpenMP parallel
c     region.  You can use falloc_enter_parallel and falloc_exit_parallel
c     if you need to be aware of when you are in a parallel region.
c
c     (The word "thunk" means a small piece of code, executed before the
c     start of a subroutine, to establish proper conditions for the execution
c     of the subroutine.)




c     elfmm_thunk - Interfaces to thunks written in C.

      module elfmm_thunk
      implicit none

      contains




      subroutine tria3nrem_thunk(ier,s_stack,vert1,vert2,vert3,i_fun,nm,
     1    par1,par2,m,s_vals,nnmax,eps,
     2    rint,maxdepth,maxrec,nrec,ifrel,
     $    numint,iw,numfunev,
     3    s_value1,s_value2,s_value3,s_value4,s_vals2,
     &    s_irec,s_rnodes,s_weights,s_work)

      implicit none

      integer, intent(out) :: ier
      integer, intent(in) :: s_stack
      real(8), intent(in) :: vert1(3)
      real(8), intent(in) :: vert2(3)
      real(8), intent(in) :: vert3(3)
      integer, intent(in) :: i_fun
      integer, intent(in) :: nm
      real(8), intent(in) :: par1(*)
      real(8), intent(in) :: par2(*)
      integer, intent(in) :: m
      integer, intent(in) :: s_vals
      integer, intent(in) :: nnmax
      real(8), intent(in) :: eps
      real(8), intent(inout) :: rint(nm)
      integer, intent(in) :: maxdepth
      integer, intent(inout) :: maxrec
      integer, intent(in) :: nrec
      integer, intent(in) :: ifrel
      integer, intent(out) :: numint
      integer, intent(in) :: iw
      integer, intent(out) :: numfunev
      integer, intent(in) :: s_value1
      integer, intent(in) :: s_value2
      integer, intent(in) :: s_value3
      integer, intent(in) :: s_value4
      integer, intent(in) :: s_vals2
      integer, intent(in) :: s_irec
      integer, intent(in) :: s_rnodes
      integer, intent(in) :: s_weights
      integer, intent(in) :: s_work

      real(8), allocatable :: stack(:)
      real(8), allocatable :: vals(:)
      real(8), allocatable :: value1(:)
      real(8), allocatable :: value2(:)
      real(8), allocatable :: value3(:)
      real(8), allocatable :: value4(:)
      real(8), allocatable :: vals2(:)
      integer, allocatable :: irec(:)
      real(8), allocatable :: rnodes(:)
      real(8), allocatable :: weights(:)
      real(8), allocatable :: work(:)

      allocate( stack(s_stack) )
      allocate( vals(s_vals) )
      allocate( value1(s_value1) )
      allocate( value2(s_value2) )
      allocate( value3(s_value3) )
      allocate( value4(s_value4) )
      allocate( vals2(s_vals2) )
      allocate( irec(s_irec) )
      allocate( rnodes(s_rnodes) )
      allocate( weights(s_weights) )
      allocate( work(s_work) )

      call tria3nrem(ier,stack,vert1,vert2,vert3,i_fun,nm,
     1    par1,par2,m,vals,nnmax,eps,
     2    rint,maxdepth,maxrec,nrec,ifrel,
     $    numint,iw,numfunev,
     3    value1,value2,value3,value4,vals2,
     &    irec,rnodes,weights,work)

      return
      end subroutine tria3nrem_thunk




      end module elfmm_thunk







 