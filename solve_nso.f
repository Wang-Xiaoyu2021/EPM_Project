      subroutine solve_nso(kpc,vr,work,inG,ev,proj,
     &    iproj,zz_st)

       implicit none

       include 'epm.inc'

       integer lwork,iproj,iflag,j1,j2,j
       parameter (lwork=10000)
       complex*16 workx(lwork)
       real*8 workrx(3*mnpw)

       real*8 proj(mnpw*2,mnpw*2)
       complex*16 zz_st(mnpw*2,mnpw*2)
       complex*16 zz(mnpw,mnpw)

       integer nl1,nh1,nl2,nh2,nl3,nh3,info
       integer i,ipw,jpw,ih,nd1,nd2,nd3,inG,naux
       double precision kpc(3),vr(1:mngrid),gkk
       double precision vav,aux(4*mnpw),ev(1:2*mnpw)
       complex*16 work(1:mngrid),dummy
cccc       complex*16 hd(1:mnset)
       complex*16 hd(mnpw,mnpw)
       complex*16 cc


       nl1=-ng1(is)/2
       nh1=ng1(is)/2
       nl2=-ng2(is)/2
       nh2=ng2(is)/2
       nl3=-ng3(is)/2
       nh3=ng3(is)/2

C==== if inG >= 0, vr is passed, do FFT
       if(inG .ge. 0) then
        do 10 i=1,ngrid(is)
          work(i)=dcmplx(vr(i),0.0d0)
10      continue
        call fft (work,ng1(is),ng2(is),ng3(is),-1)
       endif

       vav=dreal(work(1))

       ih=0
       do 20 jpw=1,npw
          do 30 ipw=1,jpw
             ih=ih+1
             if (ipw.eq.jpw) then
                gkk=(kpc(1)+gpw(ipw,1))**2+
     &          (kpc(2)+gpw(ipw,2))**2+(kpc(3)+gpw(ipw,3))**2
                gkk=gkk*sigma
cc                hd(ih)=vav*wg(ipw)**2+gkk
                hd(ipw,jpw)=vav*wg(ipw)**2+gkk
             else
                nd1=ngpw(jpw,1)-ngpw(ipw,1)
                nd2=ngpw(jpw,2)-ngpw(ipw,2)
                nd3=ngpw(jpw,3)-ngpw(ipw,3)
                if (nd1.le.nh1.and.nd2.le.nh2.and.nd3.le.nh3) then
                 if (nd1.gt.nl1.and.nd2.gt.nl2.and.nd3.gt.nl3) then
cc                   hd(ih)=work(pfft(nd1,nd2,nd3))*wg(ipw)*wg(jpw)
                   hd(ipw,jpw)=work(pfft(nd1,nd2,nd3))*wg(ipw)*wg(jpw)
                 else if (nd1.eq.nl1.or.nd2.eq.nl2.or.nd2.eq.nl3) then
cc                   hd(ih)=dconjg(work(pfft(-nd1,-nd2,-nd3)))*
cc     &              wg(ipw)*wg(jpw)
                   hd(ipw,jpw)=dconjg(work(pfft(-nd1,-nd2,-nd3)))*
     &              wg(ipw)*wg(jpw)
                 else
cc                   hd(ih)=dcmplx(0.0d0,0.0d0)
                   hd(ipw,jpw)=dcmplx(0.0d0,0.0d0)
                 end if
                else
cc                   hd(ih)=dcmplx(0.0d0,0.0d0)
                   hd(ipw,jpw)=dcmplx(0.0d0,0.0d0)
                end if
             end if

            hd(jpw,ipw)=dconjg(hd(ipw,jpw))
30        continue
20     continue

C===== now, diagonalize the Hamiltonian
        naux=4*mnpw
c	call zhpev (20,hd,ev,dummy,mnpw,npw,aux,naux)

       if(iproj.eq.0) then
c       call cheev('N','U',npw,hd,mnpw,ev,workx,lwork,
c     &    workrx,info)
       call zheev('N','U',npw,hd,mnpw,ev,workx,lwork,
     &    workrx,info)
       else
c       call cheev('V','U',npw,hd,mnpw,ev,workx,lwork,
c     &    workrx,info)
       call zheev('V','U',npw,hd,mnpw,ev,workx,lwork,
     &    workrx,info)
       do i=1,npw
       do j=1,npw
       zz(i,j)=hd(i,j)
       enddo
       enddo
       endif


      if(iproj.eq.2) then
      open(10,file="ug.out",form="unformatted")
      rewind(10)
      write(6,*) "npw=",npw
      iflag=1
      write(10) iflag, npw
      write(10) ((zz(j1,j2),j1=1,npw),j2=1,npw)
      close(10)
      write(6,*) "dump wavefunction in wgz.dump, then stop"
      do i=1,16
      write(6,*) i,ev(i)*27.211396/2
      enddo
      stop
      endif

      if(iproj.eq.1) then
      do j1=1,13
      do j2=1,13
      cc=dcmplx(0.d0,0.d0)
      do i=1,npw
      cc=cc+zz(i,j1)*dconjg(zz_st(i,j2))
      enddo
      proj(j2,j1)=cdabs(cc)**2
      enddo
      enddo
      endif



      return
      end
