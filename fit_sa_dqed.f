      program fit

      implicit none
      
      include 'epm.inc'
      
      integer nparam,mparam,lj,lv,liw,lw
      double precision temper0, temper, tempdec
      integer ndecay, istep, maxcal1

      parameter(mparam=24,lj=mproperty,lv=mparam,liw=1)
      parameter(lw=6*mparam+mparam*mproperty+2*mproperty+
     &     mparam*(mparam-1)/2)
      
      integer niter,nf,iw(liw),ifail,maxcal,iprint
      double precision fsumsq,fvecc(mproperty),fjacc(lj,mparam)
      double precision v(lv,mparam),w(lw),eta,xtol,s(mparam)
      double precision stepmx,xc(mparam)

c dqed variables
      integer, parameter :: mcon=0
      integer, parameter :: mequa=1
      integer iopt(24)
      real*8 ropt(3),fx(5)


      integer ip,jp,ips,i,j,k
      double precision lambda,xp(1:mnpar,1:mnpar+1)
      character*20 fileps
      double precision pamoeba(mparam+1,mparam), yamoeba(mparam+1)
      double precision pbest(mparam), ybest, dx
      double precision funk,funklocal , rnd
      external funjwl,monitor,funk,funklocal
     
      
      call readdata(nparam)

      nparamf=nparam
c
c read input file
c
      open (unit=4,file="sa.input",status='old')
      read (4,*) temper0, tempdec, ndecay, maxcal
      read (4,*) dx,stepmx,xtol
      do i=1,nparamf
      read (4,*) amin(i), amax(i), npgrid(i)
      enddo
      close (unit=4)
c
c Divide the epm parameters into subcell
c
      ncpgrid=0
      gridlocal=0
      ilocal=0
      agmin=0.d0
      agmax=0.d0
      fbest=1.d20
      xbest=0.d0

      do i=1,nparamf
      dx=(AMAX(i)-AMIN(i))/npgrid(i)
      do j=1,npgrid(i)
      agmin(j,i)=AMIN(i)+(j-1)*dx
      agmax(j,i)=AMIN(i)+j*dx
      enddo
      enddo

c
c Initializtion for simulation anneling
c 
       pamoeba =0.0d0
       yamoeba=0.0d0

c       pamoeba(nparam+1,1:nparam)=1.0d0
       pamoeba(nparam+1,1:nparam)=0.0d0
       do i=1,nparam
           pamoeba(i,1:nparam)=pamoeba(1+nparam,1:nparam)
           pamoeba(i,i) = pamoeba(i,i) + dx
       enddo

       numiter = 0
       numiter1 = 49
       do i=nparam+1,1,-1
           yamoeba(i)=funklocal(pamoeba(i,1:nparam))
       enddo
       pbest(1:nparam)=1.0d0
       ybest = 1.0d30
       ifail=0
       iprint=1
!       maxcal=8000*nparam
       eta=0.5d0
c       xtol=0.1d-12
c       stepmx=100.0d0

       call random_seed()

       numiter = 0
       numiter1 = 0
c
c call simulation annealing
c
       do istep = 0, ndecay
       temper = temper0*exp(-istep/tempdec)
       maxcal1 =maxcal
       write(*,'(a,1x,i3,10x,a,1x,f12.5)') "Step= ", istep, "T=", temper
       call amebsa(pamoeba(1:nparam+1,1:nparam),yamoeba(1:nparam+1),
     & nparam+1,nparam,nparam,pbest(1:nparam), 
     & ybest, xtol,funklocal,maxcal1, temper)
       write(6,*) 'The number of iterations performed = ',maxcal-maxcal1
       write(*,*) "MC Step",  istep,"finished with:"
       numiter1=iprt-1
       call funjwl(fx,iopt,mcon,mequa,nparamf,ropt,pbest)
       enddo
c
c print the last EPM parameters and comparison table
c
      write(6,*)
      write(6,*) "Finish global optimiation"
      do i=1,nbest
      numiter1=iprt-1
      call funjwl(fx,iopt,mcon,mequa,nparamf,ropt,xbest(1,i))
      enddo

      call timestamp()
      end
