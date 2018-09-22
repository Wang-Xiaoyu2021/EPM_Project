      function funklocal (X0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc The EPMFIT code original was written by Lin-Wang Wang.
ccc Gabriel Bester at March 2007 introduced the dqed code package
ccc as the Minimization method which work well for EPMFIT.
ccc Jun-Wei Luo at April 2007 adjusts the code in order to work well
ccc with multiepm.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'epm.inc'

      REAL*8   X0(*),funklocal
      external funx,monitor,funjwl

c DQED variables
      integer, parameter :: liwork =2000
      integer, parameter :: lwork = 20000
      integer mcon
      integer mequa
      integer ldfj
      integer, allocatable :: ind(:)
      real ( kind = 8 ), allocatable :: bl(:)
      real ( kind = 8 ), allocatable :: bu(:)
      real ( kind = 8 ), allocatable :: fj(:,:)
      real ( kind = 8 ), allocatable :: x(:)
      real ( kind = 8 ), allocatable :: fx(:)
      real ( kind = 8 ) work(lwork)
      real ( kind = 8 ) ropt(3)
      real ( kind = 8 ) fnorm
      integer igo
      integer iopt(24)
      integer iwork(liwork)
      external dqedhdjwl

c local variables
      integer i,j,k, iempty, ig,nvars
      integer ipgrid(nparamf)
      real*8 X1(nparamf),dlta
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccc  initialize for DQED
c
      nvars=nparamf
      mcon=0     ! no linear constrain equation on unknowns
      mequa=1    ! just one nonlinear equation
      ldfj=mcon+mequa

      allocate(ind(nvars+mcon))
      allocate(bl(nvars+mcon))
      allocate(bu(nvars+mcon))
      allocate(fj(ldfj,nvars+1))
      allocate(x(nparamf))
      allocate(fx(ldfj))

c
c Define the bounding constraints on the variables X(i) >=0, and
c define the initial values of the variables.
c
      do i=1,nvars
      dlta=amax(i)-amin(i)
      X1(i)=mod(X0(i)-amin(i)+5*dlta,dlta)+amin(i)
      enddo

      iempty=0
      do i=1,nvars
      ig=X1(i)*npgrid(i)*0.99999+1
      ipgrid(i)=ig
      if(ncpgrid(ig,i).eq.0) iempty=1
      ncpgrid(ig,i)=ncpgrid(ig,i)+1
      ind(i)=3
      bl(i)=10**agmin(ig,i)
      bu(i)=10**agmax(ig,i)
      x(i)=0.5*(bl(i)+bu(i))
      enddo
c
c if this grid has been called, return.
c
      if(iempty.eq.0) then
      do 50 j=1,ilocal
      do i=1,nvars
      if(ipgrid(i).ne.gridlocal(j,i)) goto 50
      enddo
      write(6,*) "debug, return"
      funklocal=flocal(j)
      return
50    continue
      endif
   
      ilocal=ilocal+1
      do i=1,nvars
      gridlocal(ilocal,i)=ipgrid(i)
      enddo
c
c  Tell how much storage we gave the solver.
c
      iwork(1) = lwork
      iwork(2) = liwork

c
c Change the number of iterations, default is 75.
c
      iopt(1)=2
      iopt(2)=75
c
c Set the option to change the value of TOLF.
c 
      iopt(3)=4
c
c The next entry points to the place in ROPT where
c the new value of TOLF islocated.
c
      iopt(4)=1
c
c This is the new value of TOLF
c
      ropt(1)=1.d-2
c
c This next option is a signal that there are no more options.
c
      iopt(5)=99

c
c Call the program
c
      call dqed ( dqedhdjwl, mequa, nvars, mcon, ind, bl, bu, x, fj,
     &   ldfj, fnorm, igo, iopt, ropt, iwork, work )

      write (6,'(a,i6)') 'Output flag from DQED, IGO = ', igo
      write (6, '(a)' ) ' '
      write (6, '(a)' ) 'Computed X:'
      write (6, '(a)' ) ' '
      do i=1,nvars
      write (6, '(4(g14.6,2x))' ) x(i), bl(i), bu(i), X1(i)
      enddo
      write (6, '(a)' ) ' '
      write (6, '(a,g14.6)' ) 'L2 norm of the residual, FNORM = ',fnorm
      write (6, '(a)' ) ' '
      write (6, '(a,i6)' ) 'EPMFIT_DQED, NO.', ilocal
      write (6, '(a)' ) 'Normal end of execution.'

c
c print the last EPM parameters and comparison table
c
      numiter1=iprt-1
      call funjwl(fx,iopt,mcon,mequa,nvars,ropt,x)
      funklocal=fx(1)
      flocal(ilocal)=funklocal
c
c Save the best solutions
c
      do i=1,nbest
      if(funklocal.ge.fbest(i)) goto 60
      enddo

60    continue

      i=i-1
      if(i.ge.2) then
      do j=1,i-1
      fbest(j)=fbest(j+1)
      do k=1,nvars
      xbest(k,j)=xbest(k,j+1)
      enddo
      enddo
      fbest(i)=funklocal
      do k=1,nvars
      xbest(k,i)=x(k)
      enddo
      endif   ! if(i.ge.2)

      if(i.eq.1) then
      fbest(1)=funklocal
      do k=1,nvars
      xbest(k,1)=x(k)
      enddo
      endif
c
c Deallocate array
c

      deallocate(ind)
      deallocate(bl)
      deallocate(bu)
      deallocate(fj)
      deallocate(x)
      deallocate(fx)

      return
      end

