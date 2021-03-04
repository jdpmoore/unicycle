c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the 
c        start of the actual quadrature routines.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
!!        subroutine tria3adam(ier,vert1,vert2,vert3,i_fun,nm,
!!     1      par1,par2,m,eps,rint,maxrec,numfunev,w)
        subroutine tria3adam(ier,vert1,vert2,vert3,i_fun,nm,
     1      par1,par2,m,eps,rint,maxrec,numfunev)

        implicit none

        integer, intent(out) :: ier
        real(8), intent(in) :: vert1(3)
        real(8), intent(in) :: vert2(3)
        real(8), intent(in) :: vert3(3)
        integer, intent(in) :: i_fun
        integer, intent(in) :: nm
        real(8), intent(in) :: par1(*)
        real(8), intent(in) :: par2(*)
        integer, intent(in) :: m
        real(8), intent(in) :: eps
        real(8), intent(inout) :: rint(nm)
        integer, intent(inout) :: maxrec
        integer, intent(out) :: numfunev
!!        real(8), intent(inout) :: w(*)

        integer :: nrec
        integer :: ifrel
c
        nrec=200
        ifrel=1
c
        call tria3adamar(ier,vert1,vert2,vert3,i_fun,nm,
     1      par1,par2,m,eps,rint,maxrec,numfunev,nrec,ifrel)

        return
        end
c
c
c
!!        subroutine tria3adamar(ier,vert1,vert2,vert3,i_fun,nm,
!!     1      par1,par2,m,eps,rint,maxrec,numfunev,w,nrec,ifrel)
        subroutine tria3adamar(ier,vert1,vert2,vert3,i_fun,nm,
     1      par1,par2,m,eps,rint,maxrec,numfunev,nrec,ifrel)

        use triaq_mod  ! triaq_maxnodes, triaq_worksize
        use elfmm_thunk
        implicit none

        integer, intent(out) :: ier
        real(8), intent(in) :: vert1(3)
        real(8), intent(in) :: vert2(3)
        real(8), intent(in) :: vert3(3)
        integer, intent(in) :: i_fun
        integer, intent(in) :: nm
        real(8), intent(in) :: par1(*)
        real(8), intent(in) :: par2(*)
        integer, intent(in) :: m
        real(8), intent(in) :: eps
        real(8), intent(inout) :: rint(nm)
        integer, intent(inout) :: maxrec
        integer, intent(out) :: numfunev
!!        real(8), intent(inout) :: w(*)
        integer, intent(in) :: nrec
        integer, intent(in) :: ifrel

        integer :: nnmax
        integer :: maxdepth

!!        integer :: istack
!!        integer :: lstack
!!        integer :: iw
!!        integer :: lw
!!        integer :: ivals
!!        integer :: lvals
!!        integer :: iww
!!        integer :: lww
!!        integer :: ivalue1
!!        integer :: lvalue1
!!        integer :: ivalue2
!!        integer :: lvalue2
!!        integer :: ivalue3
!!        integer :: lvalue3
!!        integer :: ivalue4
!!        integer :: lvalue4
!!        integer :: ivals2
!!        integer :: lvals2

        integer :: iiiw
        integer :: numint

        integer :: s_stack
        integer :: s_vals
        integer :: s_value1
        integer :: s_value2
        integer :: s_value3
        integer :: s_value4
        integer :: s_vals2
        integer :: s_irec
        integer :: s_rnodes
        integer :: s_weights
        integer :: s_work
c
c       This subroutine uses the adaptive symmetric quadrature (for
c       n.le.50) or the adaptive tensor product gaussian quadrature (for
c       n.ge.51) to integrate a collection of user-supplied 
c       functions R^2 \to R^1 on a triangle in R^3 (also 
c       user-supplied). In fact, this is simply a memory management 
c       routine; all actual work is done by the subroutine tria3nrem 
c       (see below).
c
c                       input parameters:
c
c  vert1,vert2,vert3 - the vertices in R^3 of the triangle
c       over which the function is to be integrated
c  i_fun - the user-supplied function to be integrated. the calling
c       sequence of fun must be 
c
c        call fun(x,y,z,par1,par2,f).                                (1)
c
c        in (1), (x,y) is a point in R^3 where 
c        the function is to be evaluated, and par1, par2
c        are two parameters to be used by fun; they can be 
c        variables or arrays, real or integer, as desired. 
c        f is assumed to be a real(8) vector of length nm
c  nm - the length of the vectors rint in the calling sequence of the
c        subroutine triaadam (see above), and of the vector f in (1) 
c        above
c  par1, par2 - parameters to be used by the user-supplied 
c       subroutine fun (see above)
c  m - the order of the quadrature to me used on each subinterval
c  eps - the accuracy (relative/absolute) to which the integrals will be 
c       evaluated
c  nrec - the depth of recursion, must not exceed 200
c  ifrel - relative/absolute computation flag
c            ifrel=0 - absolute precision
c            ifrel=1 - relative precision
c
c                       output parameters:
c
c  ier - error return code. 
c          ier=0 means normal conclusion
c          ier=4 means that at some point, the depth of recursion
c                reached nrec. 
c          ier=8 means that at some point, the depth of recursion
c                reached 200. this is a fatal error.
c          ier=16 means that the total number of subtriangles in the
c                adaptive subdivision of [a,b] turned out to be greater 
c                than nnmax*4.  this is a fatal error.
c                
c  rint - the integrals as evaluated (nm of them)
c  maxrec - the maximum depth to which the recursion went at its 
c         deepest point. can not be greater than 200, since at that
c         point ier is set to 8 and the execution of the subroutine
c         terminated.
c  numfunev - the total number of function evaluations (calls to the
c         user-supplied subroutine fun) that has occured in the 
c         calculation
c         
c                         work arrays:
c
c  w - must be at least 3*m**2+4*m+1500 real(8) elements long
c
c        . . . integrate the user-supplied function using the 
c              adaptive gaussian quadratures
c
        nnmax=100000
        maxdepth=200
c
c        allocate memory for the subroutine tria3nrem
c
!!        istack=1
!!        lstack=1207
!!c
!!        iw=istack+lstack
!!        lw=3*m**2+50
!!c
!!        ivals=iw+lw
!!        lvals=207*nm+1000
!!c
!!cccc        call prinfs('lvals=*',lvals)
!!c
!!        iww=ivals+lvals
!!        lww=4*m+10
!!c
!!        ivalue1=iww+lww
!!        lvalue1=nm+4
!!c
!!        ivalue2=ivalue1+lvalue1
!!        lvalue2=nm+4
!!c
!!        ivalue3=ivalue2+lvalue2
!!        lvalue3=nm+4
!!c
!!        ivalue4=ivalue3+lvalue3
!!        lvalue4=nm+4
!!c
!!        ivals2=ivalue4+lvalue4
!!        lvals2=nm+4
c
        iiiw=21

        s_stack = 2*3*maxdepth
        s_vals = nm*maxdepth
        s_value1 = nm
        s_value2 = nm
        s_value3 = nm
        s_value4 = nm
        s_vals2 = nm
        s_irec = maxdepth
        s_rnodes = 2*max(triaq_maxnodes,m*m)
        s_weights = max(triaq_maxnodes,m*m)
        s_work = max(triaq_worksize,4*m+5)
c
c       . . . integrate
c
!!        call tria3nrem(ier,w(istack:istack+lstack-1),
!!     &      vert1,vert2,vert3,i_fun,nm,
!!     1      par1,par2,w(iw:iw+lw-1),m,w(ivals:ivals+lvals-1),nnmax,eps,
!!     2      rint,maxdepth,maxrec,nrec,ifrel,
!!     $      numint,iiiw,numfunev,w(iww:iww+lww-1),
!!     3      w(ivalue1:ivalue1+lvalue1-1),w(ivalue2:ivalue2+lvalue2-1),
!!     &      w(ivalue3:ivalue3+lvalue3-1),w(ivalue4:ivalue4+lvalue4-1),
!!     4      w(ivals2:ivals2+lvals2-1) )

        call tria3nrem_thunk(ier,s_stack,
     &      vert1,vert2,vert3,i_fun,nm,
     1      par1,par2,m,s_vals,nnmax,eps,
     2      rint,maxdepth,maxrec,nrec,ifrel,
     $      numint,iiiw,numfunev,
     3      s_value1,s_value2,
     &      s_value3,s_value4,
     4      s_vals2,s_irec,s_rnodes,s_weights,s_work )
c
        return
        end
c
c
c
c
c
!!        subroutine tria3nrem(ier,stack,vert1,vert2,vert3,i_fun,nm,
!!     1      par1,par2,w,m,vals,nnmax,eps,
!!     2      rint,maxdepth,maxrec,nrec,ifrel,
!!     $      numint,iw,numfunev,ww,
!!     3      value1,value2,value3,value4,vals2)
        subroutine tria3nrem(ier,stack,vert1,vert2,vert3,i_fun,nm,
     1      par1,par2,m,vals,nnmax,eps,
     2      rint,maxdepth,maxrec,nrec,ifrel,
     $      numint,iw,numfunev,
     3      value1,value2,value3,value4,vals2,
     &      irec,rnodes,weights,work)
c
        use triaq_mod  ! triaq_maxnodes, triaq_worksize
        implicit none

!thunked
        integer, intent(out) :: ier
        integer, intent(in) :: maxdepth
        real(8), intent(inout) :: stack(2,3,maxdepth)
        real(8), intent(in) :: vert1(3)
        real(8), intent(in) :: vert2(3)
        real(8), intent(in) :: vert3(3)
        integer, intent(in) :: i_fun
        integer, intent(in) :: nm
        real(8), intent(in) :: par1(*)
        real(8), intent(in) :: par2(*)
!!        real(8), intent(inout) :: w(*)
        integer, intent(in) :: m
        real(8), intent(inout) :: vals(nm,maxdepth)
        integer, intent(in) :: nnmax
        real(8), intent(in) :: eps
        real(8), intent(inout) :: rint(nm)
        integer, intent(inout) :: maxrec
        integer, intent(in) :: nrec
        integer, intent(in) :: ifrel
        integer, intent(out) :: numint
        integer, intent(in) :: iw
        integer, intent(out) :: numfunev
!!        real(8), intent(inout) :: ww(*)
        real(8), intent(inout) :: value1(nm)
        real(8), intent(inout) :: value2(nm)
        real(8), intent(inout) :: value3(nm)
        real(8), intent(inout) :: value4(nm)
        real(8), intent(inout) :: vals2(nm)
        integer, intent(inout) :: irec(maxdepth)
        real(8), intent(inout) :: rnodes(2,max(triaq_maxnodes,m*m))
        real(8), intent(inout) :: weights(max(triaq_maxnodes,m*m))
        real(8), intent(inout) :: work(max(triaq_worksize,4*m+5))

        real(8) :: z1(2)
        real(8) :: z2(2)
        real(8) :: z3(2)
        real(8) :: zout(2,3,10)
        integer :: numnodes
        integer :: j
        integer :: ij
        integer :: i
        integer :: ifsubdivide
        integer :: ifdone
        integer :: irectmp
        real(8) :: dmax
c
c       This subroutine uses the adaptive symmetric quadrature (for
c       n.le.50) or the adaptive tensor product gaussian quadrature (for
c       n.ge.51) to integrate a user-supplied function R^2 \to R^1
c       on a triangle in R^3 (also user-supplied).
c
c                       input parameters:
c
c  vert1,vert2,vert3 - the vertices in R^3 of the triangle
c       over which the function is to be integrated
c  i_fun - the user-supplied function to be integrated. the calling
c       sequence of fun must be 
c
c       sequence of fun must be 
c
c        call fun(x,y,z,par1,par2,f).                                (1)
c
c        in (1), (x,y) is a point in R^3 where 
c        the function is to be evaluated, and par1, par2
c        are two parameters to be used by fun; they can be 
c        variables or arrays, real or integer, as desired. 
c        f is assumed to be a real(8) vector of length nm
c  nm - the length of the vectors rint in the calling sequence of the
c        subroutine triaadam (see above), and of the vector f in (1) 
c        above
c  par1, par2 - parameters to be used by the user-supplied 
c       subroutine fun (see above)
c  m - the order of the quadrature to be used on each subinterval
c  nnmax - maximum permitted number of subdivisions
c  eps - the accuracy (relative/absolute) to which the integrals will be 
c       evaluated
c  maxrec - maximum permitted depth of recursion
c  nrec - the depth of recursion (user-requested), must not exceed 200
c  ifrel - relative/absolute computation flag
c            ifrel=0 - absolute precision
c            ifrel=1 - relative precision
c
c                       output parameters:
c
c  ier - error return code. 
c          ier=0 means normal conclusion
c          ier=4 means that at some point, the depth of recursion
c                reached nrec. 
c          ier=8 means that at some point, the depth of recursion
c                reached 200. this is a fatal error.
c          ier=16 means that the total number of subtriangles in the
c                adaptive subdivision of the user-supplied one turned 
c                out to be greater than nnmax*4.  this is a fatal error.
c  rint - the integrals as evaluated (nm of them)
c  maxrec - the maximum depth to which the recursion went at its 
c         deepest point. can not be greater than 200, since at that
c         point ier is set to 8 and the execution of the subroutine
c         terminated.
c  numint - the total number of triangles in the subdivision divided by 
c         four. can not be greater than nnmax,  since at that
c         point ier is set to 16 and the execution of the subroutine
c         terminated.
c  numfunev - the total number of calls to the subroutine fun that
c         have been performed ("the number of function evaluations")
c
c                         work arrays:
c
c  stack - must be at least 1206 real(8) elements long
c  w - must be at least 3*m**2+50 real(8) elements long
c  vals - must be at least 200 real(8) elements long
c  ww - must be at leats m*4+10 real(8) elements long
c
c       . . . start the recursion
c
        numfunev=0
c
cccc        call prinfs('in tria3nrem, nm=*',nm)
c
        z1(1)=0
        z1(2)=0
        z2(1)=1
        z2(2)=0
        z3(1)=0
        z3(2)=1        
c
        call triaarrm(z1,stack(1,1,1),2)
        call triaarrm(z2,stack(1,2,1),2)
        call triaarrm(z3,stack(1,3,1),2)
c
c        
c        plot the top level triangle
c
        if( iw .le. 0) goto 1000
cccc        call triaplot(iw,stack(1,1,1) )
c
 1000   continue
c
        call tria3inm1(m,vert1,vert2,vert3,z1,z2,z3,i_fun,
     1      par1,par2,vals(1,1),rnodes,weights,work,nm,vals2,numnodes)
        numfunev=numfunev+numnodes
c
c       recursively integrate the thing
c
        j=1
        do 1200 ij=1,nm
        rint(ij)=0
 1200 continue
c
        irec(j)=1
c
cccc        rint(1)=0
        ier=0
        maxrec=0
        do 3000 i=1,nnmax
        numint=i
        if(j .gt. maxrec) maxrec=j
c
c       subdivide the current triangle
c
        ifsubdivide=1
c        
        if( irec(j) .ge. nrec ) then 
        ier=4 
        ifsubdivide=0
        endif
c
        if( ifsubdivide .eq. 1 ) then
c
        call triadiv(stack(1,1,j),stack(1,2,j),stack(1,3,j),zout)
        ! alias OK because arguments are intent(in)
c
        call tria3inm1(m,vert1,vert2,vert3,
     $     zout(1,1,1),zout(1,2,1),zout(1,3,1),i_fun,
     1      par1,par2,value1,rnodes,weights,work,nm,vals2,numnodes)
        ! alias OK because arguments are intent(in)
        numfunev=numfunev+numnodes
c
        call tria3inm1(m,vert1,vert2,vert3,
     $     zout(1,1,2),zout(1,2,2),zout(1,3,2),i_fun,
     1      par1,par2,value2,rnodes,weights,work,nm,vals2,numnodes)
        ! alias OK because arguments are intent(in)
        numfunev=numfunev+numnodes
c
        call tria3inm1(m,vert1,vert2,vert3,
     $     zout(1,1,3),zout(1,2,3),zout(1,3,3),i_fun,
     1      par1,par2,value3,rnodes,weights,work,nm,vals2,numnodes)
        ! alias OK because arguments are intent(in)
        numfunev=numfunev+numnodes
c
        call tria3inm1(m,vert1,vert2,vert3,
     $     zout(1,1,4),zout(1,2,4),zout(1,3,4),i_fun,
     1      par1,par2,value4,rnodes,weights,work,nm,vals2,numnodes)
        ! alias OK because arguments are intent(in)
        numfunev=numfunev+numnodes
c
        ifdone=1
c
c       ... estimate the maximum magnitude of integrals to be evaluated,
c       currently, this is done for each subdivision separately, which
c       may cause some difficulties if all functions are very small
c       inside a subregion. This should not be a big problem for
c       non-oscillatory kernels.
c
        dmax = 0
        do ij=1,nm
        if( dmax .lt. abs(vals(ij,j)) ) dmax = abs(vals(ij,j))
        enddo
c
        do 1600 ij=1,nm
c
c       ... check if absolute precision has been reached
c
        if( ifrel .eq. 0 ) then
        if( dabs(value1(ij)+value2(ij)+value3(ij)
     1      +value4(ij)-vals(ij,j)) .gt. eps) ifdone=0
        endif
c
c       ... check if relative precision has been reached
c
        if( ifrel .eq. 1 ) then
        if( dabs(value1(ij)+value2(ij)+value3(ij)
     1      +value4(ij)-vals(ij,j)) .gt. eps*dmax) ifdone=0
        endif
c
 1600 continue
c
c       if the function on this subinterval has been 
c       integrated with sufficient accuracy - add the 
c       value to that of the global integral and move up
c       in the stack
c
        if(ifdone  .eq. 0) goto 2000
c
        do 1800 ij=1,nm
        rint(ij)=rint(ij)+value1(ij)+value2(ij)+value3(ij)+value4(ij)
 1800 continue
c
cccc        rint(1)=rint(1)+value1(1)+value2(1)+value3(1)+value4(1)
        j=j-1
c
        endif
c
        if( ifsubdivide .eq. 0 ) then
c
c       ... maximum requested depth has been reached, accumulate sums
c
        do 1900 ij=1,nm
        rint(ij)=rint(ij)+vals(ij,j)
 1900 continue
c
        j=j-1
c
        endif
c
c        if the whole thing has been integrated - return
c
        if(j .eq. 0) return
        goto 3000
c
 2000 continue
c     
c       if the depth of the recursion has become excessive - bomb
c
        if((j+3) .gt. maxdepth) then
          ier=8
          return
        endif
c        
c       if the function on this subinterval has not been 
c       integrated with sufficient accuracy - move 
c       down the stack
c
        call triaarrm(zout(1,1,1),stack(1,1,j),6)
        call triaarrm(zout(1,1,2),stack(1,1,j+1),6)
        call triaarrm(zout(1,1,3),stack(1,1,j+2),6)
        call triaarrm(zout(1,1,4),stack(1,1,j+3),6)
c
        do 2200 ij=1,nm
c
        vals(ij,j)=value1(ij)
        vals(ij,j+1)=value2(ij)
        vals(ij,j+2)=value3(ij)
        vals(ij,j+3)=value4(ij)
 2200 continue
c
cccc        vals(j)=value1(1)
cccc        vals(j+1)=value2(1)
cccc        vals(j+2)=value3(1)
cccc        vals(j+3)=value4(1)
c
        irectmp=irec(j)+1
        irec(j)=irectmp
        irec(j+1)=irectmp
        irec(j+2)=irectmp
        irec(j+3)=irectmp
c
        j=j+3
c
c        plot the newly created triangles
c
        if(iw .le. 0) goto 2300  
cccc        call triaplot(iw,zout(1,1,1) )
cccc        call triaplot(iw,zout(1,1,2) )
cccc        call triaplot(iw,zout(1,1,3) )
cccc        call triaplot(iw,zout(1,1,4) )
c
 2300 continue
c     
c       if the depth of the recursion has become excessive - bomb
c
!!        if(j .le. maxdepth) goto 3000
!!        ier=8
!!        return
 3000 continue
        ier=16
        return
        end
c
c
c
c
c
!!        subroutine tria3inm1(n,vert1,vert2,vert3,z1,z2,z3,i_fun,
!!     1      par1,par2,rint,w,work,nm,vals,numnodes)
        subroutine tria3inm1(n,vert1,vert2,vert3,z1,z2,z3,i_fun,
     1      par1,par2,rint,rnodes,weights,work,nm,vals,numnodes)

        use triaq_mod  ! triaq_maxnodes, triaq_worksize
        implicit none

        integer, intent(in) :: n
        real(8), intent(in) :: vert1(3)
        real(8), intent(in) :: vert2(3)
        real(8), intent(in) :: vert3(3)
        real(8), intent(in) :: z1(2)
        real(8), intent(in) :: z2(2)
        real(8), intent(in) :: z3(2)
        integer, intent(in) :: i_fun
        real(8), intent(in) :: par1(*)
        real(8), intent(in) :: par2(*)
        integer, intent(in) :: nm
        real(8), intent(inout) :: rint(nm)
!!        real(8), intent(inout) :: w(*)
!!        real(8), intent(inout) :: work(450)
        real(8), intent(inout) :: rnodes(2,max(triaq_maxnodes,n*n))
        real(8), intent(inout) :: weights(max(triaq_maxnodes,n*n))
        real(8), intent(inout) :: work(max(triaq_worksize,4*n+5))
        real(8), intent(inout) :: vals(nm)
        integer, intent(out) :: numnodes

        integer :: ifinit

!!        integer :: irnodes
!!        integer :: lrnodes
!!        integer :: iweights
!!        integer :: lweights

        integer, parameter :: nold = -10
c
c        construct the quadrature formula on this triangle
c
        ifinit=1
        if(n .eq. nold) ifinit=0
c
!!        irnodes=1
!!        lrnodes=n**2 *2 +10
!!c
!!        iweights=irnodes+lrnodes
!!        lweights=n**2+5
c
        if( n .le. 50 ) then
!!        call triasymq(n,z1,z2,z3,w(irnodes:irnodes+lrnodes-1),
!!     1      w(iweights:iweights+lweights-1),numnodes)
        call triasymq(n,z1,z2,z3,rnodes,
     1      weights,numnodes,work)
        endif
c
        if( n .gt. 50 ) then
!!        call triagauc(n,z1,z2,z3,w(irnodes:irnodes+lrnodes-1),
!!     1      w(iweights:iweights+lweights-1),ifinit,work)
        call triagauc(n,z1,z2,z3,rnodes,
     1      weights,ifinit,work)
        numnodes=n*n
        endif
c
c       integrate the user-specified function 
c
!!        call tria3inm0(vert1,vert2,vert3,
!!     $     numnodes,w(irnodes),w(iweights),i_fun,par1,par2,
!!     1     rint,nm,vals)
!!        ! alias OK because arrays are intent(in)
        call tria3inm0(vert1,vert2,vert3,
     $     numnodes,rnodes,weights,i_fun,par1,par2,
     1     rint,nm,vals)
c
        return
        end
c
c
c
c
c
        subroutine tria3inm0(vert1,vert2,vert3,
     $     numnodes,rnodes,weights,i_fun,par1,par2,
     1     rint,nm,vals)

        use tria3ifun_enum
        implicit none

        real(8), intent(in) :: vert1(3)
        real(8), intent(in) :: vert2(3)
        real(8), intent(in) :: vert3(3)
        integer, intent(in) :: numnodes
        real(8), intent(in) :: rnodes(2,numnodes)
        real(8), intent(in) :: weights(numnodes)
        integer, intent(in) :: i_fun
        real(8), intent(in) :: par1(*)
        real(8), intent(in) :: par2(*)
        integer, intent(in) :: nm
        real(8), intent(inout) :: rint(nm)
        real(8), intent(inout) :: vals(nm)

        integer :: ij
        integer :: i
        real(8) :: area
        real(8) :: u
        real(8) :: v
        real(8) :: x
        real(8) :: y
        real(8) :: z

c       this subroutine integrates the user-specified function 
c       fun over a triangle in R^3; it uses the quadrature 
c       formula with numnodes nodes rnodes and weights weights; 

        do ij=1,nm
          rint(ij)=0
        enddo

        call tria3area(vert1,vert2,vert3,area)

        select case (i_fun)

        case (i_fun3eluh_eval)

          do i=1,numnodes

            u=rnodes(1,i)
            v=rnodes(2,i)
            x=vert1(1)+u*(vert2(1)-vert1(1))+v*(vert3(1)-vert1(1))
            y=vert1(2)+u*(vert2(2)-vert1(2))+v*(vert3(2)-vert1(2))
            z=vert1(3)+u*(vert2(3)-vert1(3))+v*(vert3(3)-vert1(3))

            call fun3eluh_eval(x,y,z,par1,par2,vals)

            do ij=1,nm
              rint(ij)=rint(ij)+vals(ij)*weights(i)*2*area
            enddo

          enddo

        case (i_fun3elth_eval)

          do i=1,numnodes

            u=rnodes(1,i)
            v=rnodes(2,i)
            x=vert1(1)+u*(vert2(1)-vert1(1))+v*(vert3(1)-vert1(1))
            y=vert1(2)+u*(vert2(2)-vert1(2))+v*(vert3(2)-vert1(2))
            z=vert1(3)+u*(vert2(3)-vert1(3))+v*(vert3(3)-vert1(3))

            call fun3elth_eval(x,y,z,par1,par2,vals)

            do ij=1,nm
              rint(ij)=rint(ij)+vals(ij)*weights(i)*2*area
            enddo

          enddo

        case (i_fun3elt_eval)

          do i=1,numnodes

            u=rnodes(1,i)
            v=rnodes(2,i)
            x=vert1(1)+u*(vert2(1)-vert1(1))+v*(vert3(1)-vert1(1))
            y=vert1(2)+u*(vert2(2)-vert1(2))+v*(vert3(2)-vert1(2))
            z=vert1(3)+u*(vert2(3)-vert1(3))+v*(vert3(3)-vert1(3))

            call fun3elt_eval(x,y,z,par1,par2,vals)

            do ij=1,nm
              rint(ij)=rint(ij)+vals(ij)*weights(i)*2*area
            enddo

          enddo

        case (i_fun3elthb_eval)

          do i=1,numnodes

            u=rnodes(1,i)
            v=rnodes(2,i)
            x=vert1(1)+u*(vert2(1)-vert1(1))+v*(vert3(1)-vert1(1))
            y=vert1(2)+u*(vert2(2)-vert1(2))+v*(vert3(2)-vert1(2))
            z=vert1(3)+u*(vert2(3)-vert1(3))+v*(vert3(3)-vert1(3))

            call fun3elthb_eval(x,y,z,par1,par2,vals)

            do ij=1,nm
              rint(ij)=rint(ij)+vals(ij)*weights(i)*2*area
            enddo

          enddo

        end select

        return
        end
c
c
c
c
c
c
        subroutine tria3area(vert1,vert2,vert3,area)

        implicit none

        real(8), intent(in) :: vert1(3)
        real(8), intent(in) :: vert2(3)
        real(8), intent(in) :: vert3(3)
        real(8), intent(out) :: area

        real(8) :: x(3)
        real(8) :: y(3)
        real(8) :: z(3)
c
c       compute area of an oriented triangle in R^3
c
c       INPUT: 
c   
c       vert1(3),vert2(3),vert3(3)   -  triangle vertices
c
c       OUTPUT: 
c   
c       area          -  triangle area
c
        x(1)=vert2(1)-vert1(1)
        x(2)=vert2(2)-vert1(2)
        x(3)=vert2(3)-vert1(3)
c
        y(1)=vert3(1)-vert1(1)
        y(2)=vert3(2)-vert1(2)
        y(3)=vert3(3)-vert1(3)
c
        z(1)=x(2)*y(3)-x(3)*y(2)
        z(2)=x(3)*y(1)-x(1)*y(3)
        z(3)=x(1)*y(2)-x(2)*y(1)
c
        area=dsqrt(z(1)**2+z(2)**2+z(3)**2)
        area=area*0.5d0
c
        return
        end
c
c
