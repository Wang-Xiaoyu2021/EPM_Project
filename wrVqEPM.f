      program wrVqEPM
cccccccccc this program writes out the potential from a
cccccc potential parameter file, to a vq.atom file
cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc vq(kk), vq(kk+1) are two cations for heterostructure
cccc vq(kk+2) is the anion of vq(kk)
cccc vq(kk+3) is the anion of vq(kk+1)
cccc vq is in the unit of Ryd.
cccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision (a-h,o-z)
      parameter (mtype=8,nq=1000)
      real*8 a1(mtype),a2(mtype),a3(mtype),a4(mtype)
      real*8 b1(mtype),b2(mtype),b3(mtype),b4(mtype)
      real*8 c1(mtype),c2(mtype),c3(mtype),c4(mtype)
      real*8 f(mtype),beta(mtype),vol(mtype),y(mtype)
      real*8 fa(mtype)
      real*8 param1(3,20,mtype)
      real*8 coef1(30),coef2(30)
      character*4 atom(mtype)
      integer iepm
      integer iatom(mtype)
      open(10,file='pseudo.input')
      rewind(10)

      read(10,*) num,iepm

      if(num.gt.mtype) then
      write(6,*) "num.gt.mtype, stop", num,mtype
      stop
      endif

      if(num.lt.mtype) then
      do i=num+1,mtype
      iatom(i)=0
      enddo
      endif


      if(iepm .eq. 0) then
         write(6,*) 'EPM: 4 Gaussians'
         do ii=1,num
            read(10,*) atom(ii),iatom(ii)
            read(10,*) vol(ii),f(ii),beta(ii)
            read(10,*) a1(ii), b1(ii), c1(ii)
            read(10,*) a2(ii), b2(ii), c2(ii)
            read(10,*) a3(ii), b3(ii), c3(ii)
            read(10,*) a4(ii), b4(ii), c4(ii)
c            read(10,*) a1(ii), c1(ii), b1(ii)
c            read(10,*) a2(ii), c2(ii), b2(ii)
c            read(10,*) a3(ii), c3(ii), b3(ii)
c            read(10,*) a4(ii), c4(ii), b4(ii)
            fa(ii)=0.d0
         enddo
      endif

c      if(iepm.eq.1) then
c         write(6,*) 'EPM: new Gaussian function'
c         do ii=1,num
c            read(10,*) atom(ii),iatom(ii)
c            read(10,*) vol(ii),fa(ii),beta(ii)
c            read(10,*) a1(ii), c1(ii), b1(ii)
c            read(10,*) a2(ii), c2(ii), b2(ii)
c         enddo
c      endif

      if(iepm.eq.1) then
         write(6,*) 'EPM: new Gaussian function'
         read(10,*) Ecut,Smth,scalkin
         do ii=1,num
            read(10,*)
            read(10,*) atom(ii),iatom(ii)
            read(10,*) vol(ii)
            read(10,*) fSO
            read(10,*) fa(ii)
            read(10,*) beta(ii)
            read(10,*) a1(ii), b1(ii), c1(ii)
            read(10,*) a2(ii), b2(ii), c2(ii)
            read(10,*) a3(ii), b3(ii), c3(ii)
            read(10,*) a4(ii), b4(ii), c4(ii)
            f(ii)=0.d0
         enddo
      endif


      if(iepm.eq.5) then
         write(6,*) 'EPM: type 5, new form from the multifit.x'
         read(10,*) Ecut,Smth,scalkin
         do ii=1,num
            read(10,*)
            read(10,*) atom(ii),iatom(ii)
            read(10,*) fSO
            read(10,*) fa(ii)
            read(10,*) a1(ii)
            read(10,*) c1(ii)
            read(10,*) b1(ii)
            read(10,*) a2(ii)
            vol(ii)=1.d0
            beta(ii)=0.d0
            c2(ii)=2.d0
            b2(ii)=0.d0
         enddo
      endif
      if(iepm.eq.51) then
        write(6,*) 'EPM: type 51, formula of Cohen'
        read(10,*) Ecut,Smth,scalkin
        do ii=1,num
        read(10,*)
        read(10,*) atom(ii),iatom(ii)
         read(10,*) fSO
         read(10,*) fa(ii)
         read(10,*) a1(ii)
         read(10,*) c1(ii)
         read(10,*) b1(ii)
         read(10,*) a2(ii)
            vol(ii)=1.d0
            beta(ii)=0.d0
            c2(ii)=2.d0
            b2(ii)=0.d0
         enddo
       endif
       if(iepm.eq.52) then
        write(6,*) 'EPM: type 52, formula of Schluter'
        read(10,*) Ecut,Smth,scalkin
        do ii=1,num
        read(10,*)
        read(10,*) atom(ii),iatom(ii)
         read(10,*) fSO
         read(10,*) fa(ii)
         read(10,*) a1(ii)
         read(10,*) c1(ii)
         read(10,*) b1(ii)
         read(10,*) a2(ii)
            vol(ii)=1.d0
            beta(ii)=0.d0
            c2(ii)=2.d0
            b2(ii)=0.d0
         read(10,*) b2(ii)
         read(10,*) c2(ii)
         enddo
       endif
c
c pstype = 54, exponential formula of Goano from JAP88,6467(2000).
c
       if(iepm.eq.54) then
        write(6,*) 'EPM: type 54, formula of Goano'
        read(10,*) Ecut,Smth,scalkin
        do ii=1,num
        read(10,*)
        read(10,*) atom(ii),iatom(ii)
         read(10,*) fSO
         read(10,*) fa(ii)
         read(10,*) beta(ii)
         read(10,*) a1(ii)
         read(10,*) b1(ii)
         read(10,*) c1(ii)
         read(10,*) a2(ii)
         read(10,*) b2(ii)
         read(10,*) c2(ii)
            vol(ii)=1.d0
         enddo
       endif

      if(iepm.eq.2) then
        write(6,*) "EPM: SEPM 20 Gaussians, special"

        if(num.ne.2) then
        write(6,*) "num.ne.2, stop"
        stop
        endif
        do i=1,20
        read(10,*) param1(1,i,1),param1(2,i,1),param1(3,i,1)
        enddo
        read(10,*) iatom(1),iatom(2)
        do i=1,20
        read(10,*) param1(1,i,2),param1(2,i,2),param1(3,i,2)
        enddo
        fa(1)=0.d0
        fa(2)=0.d0
      endif


      if(iepm.eq.13) then
        write(6,*) "EPM: style 1 and syle 3"
        if(num.ne.2) then
        write(6,*) "num.ne.2, stop"
        stop
        endif

ccccccc for the first atom, style 1
            read(10,*) atom(1),iatom(1)
            read(10,*) vol(1),fa(1),beta(1)
            read(10,*) a1(1), c1(1), b1(1)
            read(10,*) a2(1), c2(1), b2(1)
cccccc for the second atom, style 3
            read(10,*) atom(2),iatom(2)
            read(10,*) npower,qcut
            read(10,*) (coef1(i),i=1,npower)
            read(10,*) (coef2(i),i=1,npower)
            fa(1)=0.d0
            fa(2)=0.d0
      endif


      close(10)

      Ecut=20.d0/2.d0
c      Ecut=45.d0/2.d0
      qm=2*dsqrt(2*Ecut)


      do i=1,num
      open(10+i,file="vq.atom"//char(i+48))
      rewind(10+i)
      write(6,*) "header for ",
     & "vq.atom"//char(i+48), " : ", nq+1,iatom(i),fa(i)
      enddo

      do 100 k=0,nq

       q=k*qm/nq

      if(iepm .eq. 0) then
         do ii=1,num
            y(ii)=a1(ii)*dexp(-c1(ii)*(q-b1(ii))**2)+
     &           a2(ii)*dexp(-c2(ii)*(q-b2(ii))**2)+
     &           a3(ii)*dexp(-c3(ii)*(q-b3(ii))**2)+
     &           a4(ii)*dexp(-c4(ii)*(q-b4(ii))**2)
            y(ii)=vol(ii)*y(ii)*(1+f(ii)*dexp(-beta(ii)*q**2))
         enddo
      endif

c      if(iepm.eq.1) then
c         do ii=1,num
c            y(ii)=
c     &  vol(ii)*a1(ii)*(q*q-c1(ii))/(b1(ii)*dexp(a2(ii)*q*q)-1.0d0)
c         enddo
c      endif

      if(iepm.eq.1) then
          do ii=1,num
            y(ii)=a1(ii)*dexp(-c1(ii)*(q-b1(ii))**2)+
     &           a2(ii)*dexp(-c2(ii)*(q-b2(ii))**2)+
     &           a3(ii)*dexp(-c3(ii)*(q-b3(ii))**2)+
     &           a4(ii)*dexp(-c4(ii)*(q-b4(ii))**2)
            y(ii)=vol(ii)*y(ii)*(1+f(ii)*dexp(-beta(ii)*q**2))
         enddo
      endif


      if(iepm.eq.2) then
ccccccc special 
         y1=0.d0
         y2=0.d0
         do m=1,20
         y1=y1+param1(3,m,1)*
     &     dexp(-((q-param1(1,m,1))/param1(2,m,1))**2)
         y2=y2+param1(3,m,2)*
     &     dexp(-((q-param1(1,m,2))/param1(2,m,2))**2)
         enddo
ccccccc make it in the unit of Ryd.
         y(1)=(y1+y2)/2*2
         y(2)=(y1-y2)/2*2
         y(3)=0.d0
         y(4)=0.d0
      endif

      if(iepm.eq.13) then
        y(1)=
     &  vol(1)*a1(1)*(q*q-c1(1))/(b1(1)*dexp(a2(1)*q*q)-1.0d0)


        s1=0.d0
        s2=0.d0
        do i=1,npower
        s1=s1+coef1(i)*q**(i-1)
        s2=s2+coef2(i)/(q+1.D-10)**(i-1)
        enddo
        if(q.le.qcut) y(2)=s1
        if(q.ge.qcut) y(2)=s2

        y(1)=y(1)*2     ! to change to Ryd
        y(2)=y(2)*2
        y(3)=0.d0
        y(4)=0.d0

      endif

      if(iepm.eq.5) then
      do i=1,num
      y(i)=a1(i)*(q*q-c1(i))/(b1(i)*exp(a2(i)*q*q)-1.d0)
      enddo
      endif
      
      if(iepm.eq.51) then
      do i=1,num
c      y(i)=a1(i)*(q*q-c1(i))/(b1(i)*exp(a2(i)*q*q)-1.d0)  
      y(i)=a1(i)*(q*q-c1(i))/(exp(b1(i)*(q*q-a2(i)))+1.d0)
      enddo
      endif

      if(iepm.eq.52) then
      do i=1,num
c      y(i)=a1(i)*(q*q-c1(i))/(b1(i)*exp(a2(i)*q*q)-1.d0)
      y(i)=a1(i)*(q*q-c1(i))/(exp(b1(i)*(q*q-a2(i)))+1.d0)*
     &     (0.5d0*tanh((b2(i)-q*q)/c2(i))+0.5d0)
      enddo
      endif

      if(iepm.eq.54) then
      do i=1,num
         q2=q*q
         qmax2=beta(i)*beta(i)
         if(q2.le.qmax2) then
           vsr = 1.d0-dexp(-(b1(i)*qmax2)**c1(i))
           vsr = 1.d0/vsr
           vlr = 1.d0-dexp( -(b1(i)*(qmax2-q2))**c1(i) )
           y(i) = -a1(i)+(a2(i)+a1(i))*(1.d0-vsr*vlr)
         else
           vlr = dexp( -(b2(i)*(q2-qmax2))**c2(i) )
           y(i) = a2(i)*vlr
         endif
      enddo
      endif

      do i=1,num
      y(i)=y(i)/2     ! change the potential to unit of Hartree
      if(dabs(y(i)).le.1.d-99) y(i)=0.d0
      enddo

      do i=1,num
      write(10+i,200) q,y(i)
      enddo

100   continue

      do i=1,num
      close(10+i)
      enddo

200   format(2(E12.6,1x))

300   format(5(f10.5,1x))


      stop
      end
