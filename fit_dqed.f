      program fit
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc The EPMFIT code original was written by Lin-Wang Wang.
ccc Gabriel Bester at March 2007 introduced the dqed code package
ccc as the Minimization method which work well for EPMFIT.
ccc Jun-Wei Luo at April 2007 adjusts the code in order to work well
ccc with multiepm.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      
      include 'epm.inc'
      
      integer nparam,mparam
      external funx,monitor,funjwl

c DQED variables
      integer, parameter :: liwork =2000
      integer, parameter :: lwork = 30000
      integer mcon
      integer mequa
      integer nvars
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
      integer i,j,k,ips
      real( kind = 8 ), allocatable :: aminl(:)
      real( kind = 8 ), allocatable :: amaxl(:)
      real( kind = 8 ) :: tmp
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call timestamp( )

ccc read input files      
c
      call readdata(nparam)
      nparamf=nparam

      numiter=0
      numiter1=0
c
cccc  initialize for DQED
c
      mcon=0     ! no bounding constraint linear equation
      mequa=1    ! just one nonlinear equation
      nvars=nparam
      ldfj=mcon+mequa

      allocate(ind(nvars+mcon))
      allocate(bl(nvars+mcon))
      allocate(bu(nvars+mcon))
      allocate(fj(ldfj,nvars+1))
      allocate(x(nvars))
      allocate(fx(ldfj))
      allocate(aminl(nvars))
      allocate(amaxl(nvars))
c
c Define the bounding constraints on the variables X(i) >=0, and
c define the initial values of the variables.
c
      open(11,file="dqed.in")
      do i=1,nvars
      read(11,*) aminl(i), amaxl(i)
      enddo
      close(11)

      k=0
      do 20 ips=1,itotps
         if(moveSO(ips).eq.1) then
            k=k+1
            tmp=dabs(psSO0L(ips)/psSO0(ips))
            if(tmp.gt.aminl(k)) aminl(k)=tmp
            tmp=dabs(psSO0U(ips)/psSO0(ips))
            if(tmp.lt.amaxl(k)) amaxl(k)=tmp
         endif
         if(mvpsf(ips).eq.1) then
            k=k+1
            tmp=dabs(psf0L(ips)/psf0(ips))
            if(tmp.gt.aminl(k)) aminl(k)=tmp
            tmp=dabs(psf0U(ips)/psf0(ips))
            if(tmp.lt.amaxl(k)) amaxl(k)=tmp
         endif
         if(mvpsbeta(ips).eq.1) then
            k=k+1
            tmp=dabs(psbeta0L(ips)/psbeta0(ips))
            if(tmp.gt.aminl(k)) aminl(k)=tmp
            tmp=dabs(psbeta0U(ips)/psbeta0(ips))
            if(tmp.lt.amaxl(k)) amaxl(k)=tmp
         endif

         do 30 j=1,ngauss(ips)
            if (mvpsa(j,ips).eq.1) then
               k=k+1
               tmp=dabs(psa0L(j,ips)/psa0(j,ips))
               if(tmp.gt.aminl(k)) aminl(k)=tmp
               tmp=dabs(psa0U(j,ips)/psa0(j,ips))
               if(tmp.lt.amaxl(k)) amaxl(k)=tmp
            end if
            if (mvpsb(j,ips).eq.1) then
               k=k+1
               tmp=dabs(psb0L(j,ips)/psb0(j,ips))
               if(tmp.gt.aminl(k)) aminl(k)=tmp
               tmp=dabs(psb0U(j,ips)/psb0(j,ips))
               if(tmp.lt.amaxl(k)) amaxl(k)=tmp
            end if
            if (mvpsc(j,ips).eq.1) then
               k=k+1
               tmp=dabs(psc0L(j,ips)/psc0(j,ips))
               if(tmp.gt.aminl(k)) aminl(k)=tmp
               tmp=dabs(psc0U(j,ips)/psc0(j,ips))
               if(tmp.lt.amaxl(k)) amaxl(k)=tmp
            end if
 30      continue
 20   continue

      if (ifit_sig==1) then
      k=k+1
      aminl(k)=1.0
      amaxl(k)=2.0
      end if

      if(nvars.ne.k) then
      write(*,*) 'nvars != nparam, stop',nvars,k
      stop
      endif

      do i=1,nvars
      ind(i)=3
      bl(i)=aminl(i)
      bu(i)=amaxl(i)
      x(i)=1.d0
      enddo
c
c  Tell how much storage we gave the solver.
c
      iwork(1) = lwork
      iwork(2) = liwork

c Change the number of iterations, default is 75
      iopt(1)=2
      iopt(2)=975
      iopt(3)=99

c check the initial EPM parameters
      write(6,*) "Initial EPM:"
      numiter1=iprt-1
      call funjwl(fx,iopt,mcon,mequa,nvars,ropt,x)

c Call the program
c 
      call dqed ( dqedhdjwl, mequa, nvars, mcon, ind, bl, bu, x, fj, 
     &   ldfj, fnorm, igo, iopt, ropt, iwork, work )

      write (6,'(a,i6)') 'Output flag from DQED, IGO = ', igo
      write (6, '(a)' ) ' ' 
      write (6, '(a)' ) 'Computed X:'
      write (6, '(a)' ) ' '
      write (6, '(g14.6)' ) x(1:nvars)
      write (6, '(a)' ) ' '
      write (6, '(a,g14.6)' ) 'L2 norm of the residual, FNORM = ',fnorm
      write (6, '(a)' ) ' '
      write (6, '(a)' ) 'EPMFIT_DQED'
      write (6, '(a)' ) 'Normal end of execution.'

c print the last EPM parameters and comparison table
      numiter1=iprt-1
      call funjwl(fx,iopt,mcon,mequa,nvars,ropt,x)
  
      deallocate(ind)
      deallocate(bl)
      deallocate(bu)
      deallocate(fj)
      deallocate(x)
      deallocate(fx)
      deallocate(aminl)
      deallocate(amaxl)

      write(6,*) "Ending at"
      call timestamp()
      end

