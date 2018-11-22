      subroutine funjwl(fx,iopt,mcon,mequa,nvars,ropt,xc)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc changed from funk to funjwl by jwluo. 
ccc It can be called from subrouinte difcen of dqed code package 
ccc which is a subroutine to estimate jacobian matrix.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include 'epm.inc'
      real*8 fx(*),ropt(*),xc(*)
      integer iopt(*),mcon,mequa,nvars
      
      double precision funk, fvecc(nproperty)
      double precision s(nparamf),fjacc(1,nparamf)
      integer iflag, iw, liw, w, lw
      integer i
      integer ljc,igrade,niter,nf
      external monitor, funx

      ljc=1 
      iw=0
      liw=0
      w=0
      lw=0
      iflag=0

      if(nvars.ne.nparamf) then
      write(6,*) "funjwl: number of parameters are not same, stop"
      write(6,*) nvars,nparamf
      stop 
      endif

      if(mequa.ne.1) then
      write(6,*) "mequa/=1, stop",mequa
      stop
      endif

      call funx(iflag,nproperty,nparamf,xc,fvecc,iw,liw,w,lw)

      funk=0.0d0
      do i=1,nproperty
      funk = funk + fvecc(i)*fvecc(i)
      enddo

      fx(1)=funk

      numiter = numiter + 1
      numiter1 = numiter1 + 1
      if (numiter1 .eq. iprt) then
          call monitor(nproperty,nparamf, xc, fvecc,fjacc,
     & ljc, s, igrade, numiter, nf, iw, liw, w, lw)
          numiter1 = 0
      endif
      return
      end

