
      function funk(xc)
      include 'epm.inc'
      double precision funk, xc(nparamf),fvecc(nproperty)
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

      call funx(iflag,nproperty,nparamf,xc,fvecc,iw,liw,w,lw)

      funk=0.0d0
      do i=1,nproperty
      funk = funk + fvecc(i)*fvecc(i)
      enddo

      numiter = numiter + 1
      numiter1 = numiter1 + 1
      if (numiter1 .eq. 50) then
          call monitor(nproperty,nparamf, xc, fvecc,fjacc,
     & ljc, s, igrade, numiter, nf, iw, liw, w, lw)
          numiter1 = 0
      endif
      return
      end
