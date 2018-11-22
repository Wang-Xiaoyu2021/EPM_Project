cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c******** SUBROUTINE dqedhdjwl  ******************************
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine dqedhdjwl(x, fj, ldfj, igo, iopt, ropt)

!*****************************************************************************80!
!! DQEDHD evaluates the functions and derivatives for DQED.
!
!  Discussion:
!
!    The user problem has MCON constraint functions,
!    MEQUA least squares equations, and involves NVARS
!    unknown variables.
!
!    When this subprogram is entered, the general (near)
!    linear constraint partial derivatives, the derivatives
!    for the least squares equations, and the associated
!    function values are placed into the array FJ(*,*).
!
!    All partials and functions are evaluated at the point
!    in X(*).  Then the subprogram returns to the calling
!    program unit.  Typically one could do the following
!    steps:
!
!      step 1. Place the partials of the i-th constraint
!      function with respect to variable j in the
!      array FJ(i,j), i = 1,...,MCON, j=1,...,NVARS.
!
!      step 2. Place the values of the i-th constraint
!      equation into FJ(i,NVARS+1).
!
!      step 3. Place the partials of the i-th least squares
!      equation with respect to variable j in the
!      array FJ(MCON+i,j), i = 1,...,MEQUA,
!      j = 1,...,NVARS.
!
!      step 4. Place the value of the i-th least squares
!      equation into FJ(MCON+i,NVARS+1).
!
!      step 5. Return to the calling program unit.
!
      implicit none

      include 'epm.inc'

      integer ldfj
      integer mcon
      integer mequa
      integer nvars
      integer i,j
      integer igo
      real*8 fj(ldfj,*)
      real*8 x(*)
c       real*8, allocatable ::  x(:)
      real*8 ropt(*)
      integer iopt(*)
      external diffor, funjwl
c local
      real*8,allocatable :: fjt(:,:)
      real*8,allocatable :: fxt(:)

      nvars=nparamf
      mcon=0     ! no bounding constraint linear equation
      mequa=1    ! just one nonlinear equation

      allocate(fjt(ldfj,nvars))
      allocate(fxt(mequa+mcon))

c  get function values and derivatives
c
      call diffor(fjt,funjwl,fxt,iopt,ldfj,mcon,mequa,nvars,ropt,x)

      do j=1,mcon+mequa
      do i=1,nvars
      fj(j,i)=fjt(j,i)
      enddo
      enddo

      do i=1,mequa
      fj(mcon+i,nvars+1)=fxt(i)
      enddo

      deallocate(fjt)
      deallocate(fxt)

      return
      end subroutine dqedhdjwl

