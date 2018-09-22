C==== check the local variables and the indices of the Matrix
      subroutine solve_so (kpc,vr,work,inG,ev)

       implicit none
       include 'epm.inc'

       integer nl1,nh1,nl2,nh2,nl3,nh,inG,info,ntot2,iarray
       integer i,j,ipw,jpw,ih,nd1,nd2,nd3
       integer igx,jgx,naux,ii,nh3
       double precision kpc(3),vr(1:mngrid)
       double precision vav,ev
       double precision aux(1:4*mnpw)
       double precision pi
       double precision g1x,g1y,g1z,gxi
       double precision g2x,g2y,g2z,gxj
       double precision fi1,fi2,fj1,fj2
       double precision akkx,akky,akkz
       double precision fa1,fa2,fa3,fa4
       double precision vk1k2,gkk
       real*8 fa(mnat),fp(mnat)
       real*8 x(mnat),y(mnat),z(mnat)

       complex*16 hd(mnpw,mnpw)
       complex*16 htemp(mnpw*(2*mnpw+1))
       complex*16 work(1:mngrid),dummy
       complex*16 cai,cc1,cc2,cc3

       naux = 4*mnpw
       ntot2 = 2*2*npw

ccccccccccccccccccccccccccccc

       pi=4*datan(1.d0)
       do i=1,natoms(is)
            do j=1,itotps
               if (psnum(j)==atomnum(i,is)) then
                 fp(i)=-psSO(j)*4*pi/omega*3
                 exit
               end if
            end do          
       enddo


       do i=1,npw
       do j=1,npw
       hd(i,j)=dcmplx(0.d0,0.d0)
       enddo
       enddo

        do i=1,mnpw*(2*mnpw+1)
         htemp(i)=dcmplx(0.d0,0.d0)
        enddo

       do i=1,natoms(is)
cccccc tau(1,i,is),tau(2,i,is),tau(3,i,is) are the Cartesian Coord, not
cccccc the supercell edge coord. 
       x(i)=tau(1,i,is)*alat
       y(i)=tau(2,i,is)*alat
       z(i)=tau(3,i,is)*alat
       enddo

       cai=dcmplx(0.d0,1.d0)

       nl1=-ng1(is)/2
       nh1=ng1(is)/2
       nl2=-ng2(is)/2
       nh2=ng2(is)/2
       nl3=-ng3(is)/2
       nh3=ng3(is)/2

C==============
C  if inG >= 0, vr is passed, and therefore do FFT 
       if(inG .ge. 0) then
        do 10 i=1,ngrid(is)
          work(i)=dcmplx(vr(i),0.0d0)
10      continue
        call fft (work,ng1(is),ng2(is),ng3(is),-1)
       endif

       vav=dreal(work(1))

       do 20 jpw=1,npw
          do 30 ipw=1,npw
             if (ipw.eq.jpw) then
                gkk=(kpc(1)+gpw(ipw,1))**2+
     &          (kpc(2)+gpw(ipw,2))**2+(kpc(3)+gpw(ipw,3))**2
                gkk=gkk*sigma
                hd(ipw,jpw)=vav*wg(ipw)**2+gkk
             else
                nd1=ngpw(jpw,1)-ngpw(ipw,1)
                nd2=ngpw(jpw,2)-ngpw(ipw,2)
                nd3=ngpw(jpw,3)-ngpw(ipw,3)
                if (nd1.le.nh1.and.nd2.le.nh2.and.nd3.le.nh3) then
                 if (nd1.gt.nl1.and.nd2.gt.nl2.and.nd3.gt.nl3) then
                   hd(ipw,jpw)=work(pfft(nd1,nd2,nd3))
     &                           *wg(ipw)*wg(jpw)
                 else if (nd1.eq.nl1.or.nd2.eq.nl2.or.nd2.eq.nl3) then
                   hd(ipw,jpw)=dconjg(work(pfft(-nd1,-nd2,-nd3)))*
     &              wg(ipw)*wg(jpw)
                 else
                   hd(ipw,jpw)=dcmplx(0.0d0,0.0d0)
                 end if
                else
                   hd(ipw,jpw)=dcmplx(0.0d0,0.0d0)
                end if
             end if
30        continue
20     continue

C===== repack the Hamiltonian matrix
      do ipw = 1,npw
        do jpw = ipw,npw
         htemp((ntot2-ipw)*(ipw-1)/2+jpw) = hd(ipw,jpw)
         htemp((ntot2-ipw-npw)*(ipw+npw-1)/2+jpw+npw) = hd(ipw,jpw)
        enddo
      enddo

      do 40 i=1,npw
      g1x=kpc(1)+gpw(i,1)
      g1y=kpc(2)+gpw(i,2)
      g1z=kpc(3)+gpw(i,3)
      gxi=dsqrt(g1x**2+g1y**2+g1z**2)

      igx=1+gxi/gA(100)*99.d0
      fi1=(gA(igx+1)-gxi)/(gA(igx+1)-gA(igx))
      fi2=(gxi-gA(igx))/(gA(igx+1)-gA(igx))

      do 40 j=1,npw
      g2x=kpc(1)+gpw(j,1)
      g2y=kpc(2)+gpw(j,2)
      g2z=kpc(3)+gpw(j,3)
      gxj=dsqrt(g2x**2+g2y**2+g2z**2)

      jgx=1+gxj/gA(100)*99.d0
      fj1=(gA(jgx+1)-gxj)/(gA(jgx+1)-gA(jgx))
      fj2=(gxj-gA(jgx))/(gA(jgx+1)-gA(jgx))

      if(gxi.lt.1.E-8.or.gxj.lt.1.E-8) then
      akkx=0.d0
      akky=0.d0
      akkz=0.d0
      else
      akkx=(g1y*g2z-g1z*g2y)/gxi/gxj
      akky=(g1z*g2x-g1x*g2z)/gxi/gxj
      akkz=(g1x*g2y-g1y*g2x)/gxi/gxj
      endif

      fa1=fi1*fj1
      fa2=fi2*fj1
      fa3=fi1*fj2
      fa4=fi2*fj2

      vk1k2=vso(igx,jgx)*fa1+vso(igx+1,jgx)*fa2
     &  +vso(igx,jgx+1)*fa3+vso(igx+1,jgx+1)*fa4

ccccccccccc this is a very heavy computation, it losts all the
ccccccccccc advantages of using vr -> vg fft. That is why it is so
ccccccccccc slow for SO calculation.
      do ii=1,natoms(is)

cccccc not sure about the sign
      cc1=cdexp(-cai*(g1x*x(ii)+g1y*y(ii)+g1z*z(ii)))
      cc2=cdexp(-cai*(g2x*x(ii)+g2y*y(ii)+g2z*z(ii)))
      cc3=dconjg(cc1)*cc2

      if(j .ge. i) then
         iarray = (ntot2-i)*(i-1)/2+j
         htemp(iarray)=htemp(iarray)+cai*akkz*vk1k2*fp(ii)*cc3*
     $        atweight0(ii,is)

         iarray = (ntot2-i-npw)*(i+npw-1)/2+j+npw
         htemp(iarray)=htemp(iarray)-cai*akkz*vk1k2*fp(ii)*cc3*
     $        atweight0(ii,is)
      endif

      iarray = (ntot2-i)*(i-1)/2+j+npw
      htemp(iarray)=htemp(iarray)+(cai*akkx+akky)*vk1k2*fp(ii)*cc3*
     $        atweight0(ii,is)

      enddo

40    continue

      call zhpev (0,htemp,ev,dummy,2*mnpw,2*npw,aux,naux)

      end
