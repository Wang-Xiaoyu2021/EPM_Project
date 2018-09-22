      subroutine funx (iflag,nprop,nparam,xc,fvecc,iw,liw,w,lw)
      
      implicit none

      include 'epm.inc' 

      integer k,i,j,ips,writeout,iflag,nparam,liw,lw,nprop
      integer iw(liw),icase,ii,ist,itest
      double precision e_xx,e_zz,dvol,hydrofac
      double precision xc(nparam),fvecc(mproperty),w(lw)
      double precision kp(1:3),ev(1:mnpw*2),vr(1:mngrid)
      double precision pi,vbmoffset,temp,temp1,temp2,dltk,sum
      double precision cbm_shift,vbm_shift,psflocal(mnat)
      complex*16 psi,vc(1:mngrid)
      real*8 proj(mnpw*2,mnpw*2)
      complex*16 zz_st(mnpw*2,mnpw*2,mns)
      integer iproj(mns),ipcbm(mns)
      integer ibandcbm,ibandvbm,ibandvbm_1,ibandvbm_2
      real*8 s,smax
      integer j1

ccc added by jwluo
      integer ikpt,nkpt,nkpt_cb1,nkpt_vb1,indx
      integer ig1c(nstruct)
      double precision sscb1,ssvb1,kpr(1:3),evGW(1:10)
      double precision vbm_1offset,vbmstrn,dltkx,dltky,dltkz
      double precision tempel1c,tempex1c
      character*30 filess,fileband,filekpt

      common  /comproj/zz_st,iproj,ipcbm

      e_zz = 0.001d0

      pi=4.0d0*atan(1.0d0)

C======================================================
C==== rather tricky here
C==== psf(cation) = (dV/V)*psf0(ips)*exp(-10*q*q)
C==== using a local variable, psflocal(ips), to get around the problem 
C==== whenever the unit cell changes, set 
C==== psf(cation)= (dV/V)*psflocal(cation)*exp(-10*q*q)
C=======================================================
      k=0
      do 20 ips=1,itotps
         if(moveSO(ips).eq.1) then
            k=k+1
            psSO(ips) = xc(k)*psSO0(ips)
         endif
         if(mvpsf(ips).eq.1) then
            k=k+1
            psfvol(ips) = xc(k)*psf0(ips)
         endif
         if(mvpsbeta(ips).eq.1) then
            k=k+1
            psbeta(ips) = xc(k)*psbeta0(ips)
         endif
         
         do 30 j=1,ngauss(ips)
            if (mvpsa(j,ips).eq.1) then
               k=k+1
               psa(j,ips)=xc(k)*psa0(j,ips)
            end if
            if (mvpsb(j,ips).eq.1) then
               k=k+1
               psb(j,ips)=xc(k)*psb0(j,ips)
            end if
            if (mvpsc(j,ips).eq.1) then
               k=k+1
               psc(j,ips)=xc(k)*psc0(j,ips)
            end if
 30      continue
 20   continue
      
      if (ifit_sig==0) then
         sigma=sigma0
      else
         k=k+1
         sigma=xc(k)*sigma0
      end if

      if (k.ne.nparam) then
         write(6,*)'oh no'
         stop
      end if

c     Loop over each of the structures
      do is=1,nstruct
         
c     Reweight all the atoms in the structure by (1+alpha_4.deltav)
         do i=1,natoms(is)
            do j=1,itotps
               if (psnum(j)==atomnum(i,is)) then
                  atweight(i,is)=atweight0(i,is)
     $                 *(1+psfvol(j)*deltav(i,is))
                  exit
               end if
            end do
         end do
c -----------------------------------------------------
c for the SQS structure, call funxsqs()
c -----------------------------------------------------
         if(natoms(is) .gt. 8) then
         call funxsqs ()
         goto 999
         endif 
c -----------------------------------------------------
c define structure parameters; 
c for the ideal zincblend structure, psf = 0
c -----------------------------------------------------
         psf(1:itotps)=0.0d0
         
         alat=alat0(is)
         a(:,:)=a0(:,:,is)*alat
         tau(:,:,is)=tau0(:,:,is)

         call gridg ()
         call vcell (vr,vc,-1)
c ---------------------------------------------------
c Gamma-point energy levels
c ---------------------------------------------------
         kp(1)=0.0d0    !Gamma point
         kp(2)=0.0d0    !
         kp(3)=0.0d0    !

         call pwk(kp,vr,vc,-1,ev,proj,iproj(is),zz_st(1,1,is))

         if(iproj(is).eq.0) then
         vbmoffset=ev(iband(4,is))
         vbm_1offset=ev(iband(3,is))    ! for structures with crystal field
         ig1c(is)=iband(5,is)
         if(ievbm(is) .gt. 0) calcvalue(ievbm(is),is)
     $        =ev(iband(4,is))
         if(ievbm_1(is) .gt. 0) calcvalue(ievbm_1(is),is)
     $        =ev(iband(3,is))-vbmoffset
         if(ieg1c(is) .gt. 0) calcvalue(ieg1c(is),is)
     $        =ev(iband(5,is))-vbmoffset
         if(ieg15c(is) .gt. 0) calcvalue(ieg15c(is),is)
     $        =ev(iband(6,is))-vbmoffset
         if(ieg1v(is) .gt. 0) calcvalue(ieg1v(is),is)
     $        =ev(iband(1,is))-vbmoffset
         if(idso0(is) .gt. 0) calcvalue(idso0(is),is)
     $        =ev(iband(4,is))-ev(iband(2,is))
         if(idso0p(is) .gt. 0) calcvalue(idso0p(is),is)
     $        =ev(iband(7,is))-ev(iband(6,is))
         if(idcf0(is).gt. 0) calcvalue(idcf0(is),is)
     $        =ev(iband(4,is))-ev(iband(3,is))
c swith ordering between s-type and p-type bands
         if(natoms(is).eq.2 .and. ev(iband(8,is))
     &      -ev(iband(7,is)).gt.1.d-4) then
           if(ieg1c(is) .gt. 0) calcvalue(ieg1c(is),is)
     $         =ev(iband(8,is))-vbmoffset
           if(ieg15c(is) .gt. 0) calcvalue(ieg15c(is),is)
     $         =ev(iband(5,is))-vbmoffset
           ig1c(is)=iband(8,is)
           if(idso0p(is) .gt. 0) calcvalue(idso0p(is),is)
     $         =ev(iband(6,is))-ev(iband(5,is))
         endif
         endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc to decide ibandcbm based on projection
         if(iproj(is).eq.1) then

         if(iSOps.eq.1) then

         smax=-100.d0
         do j1=1,25,2
         s=proj(ipcbm(is),j1)+proj(ipcbm(is)+1,j1)
         if(s.gt.smax) then
         smax=s
         ibandcbm=j1
         endif
         enddo
         if(smax.lt.0.1) then
         write(6,*) "warning, small smax, ibandcbm", smax, ibandcbm
         endif

         if(ibandcbm.gt.iband(4,is)) then
         ibandvbm=iband(4,is)
         else
         ibandvbm=iband(4,is)+2
         endif

         if(ibandcbm.gt.iband(3,is)) then
         ibandvbm_1=iband(3,is)
         else
         ibandvbm_1=iband(3,is)+2
         endif

         if(ibandcbm.gt.iband(2,is)) then
         ibandvbm_2=iband(2,is)
         else
         ibandvbm_2=iband(2,is)+2
         endif

         endif      ! iSOps=1

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if(iSOps.eq.0) then

         smax=-100.d0
         do j1=1,13
         s=proj(ipcbm(is),j1)
         if(s.gt.smax) then
         smax=s
         ibandcbm=j1
         endif
         enddo
         if(smax.lt.0.1) then
         write(6,*) "warning, small smax, ibandcbm", smax, ibandcbm
         endif

         if(ibandcbm.gt.iband(4,is)) then
         ibandvbm=iband(4,is)
         else
         ibandvbm=iband(4,is)+1
         endif

         if(ibandcbm.gt.iband(3,is)) then
         ibandvbm_1=iband(3,is)
         else
         ibandvbm_1=iband(3,is)+1
         endif

         if(ibandcbm.gt.iband(2,is)) then
         ibandvbm_2=iband(2,is)
         else
         ibandvbm_2=iband(2,is)+1
         endif

         endif      ! iSOps=0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         vbmoffset=ev(ibandvbm)

         if(ievbm(is) .gt. 0) calcvalue(ievbm(is),is)
     $        =ev(ibandvbm)
         if(ievbm_1(is) .gt. 0) calcvalue(ievbm_1(is),is)
     $        =ev(ibandvbm_1)-vbmoffset
         if(ieg1c(is) .gt. 0) calcvalue(ieg1c(is),is)
     $        =ev(ibandcbm)-vbmoffset
         if(ieg15c(is) .gt. 0)calcvalue(ieg15c(is),is)
     $        =ev(iband(6,is))-vbmoffset
         if(ieg1v(is) .gt. 0) calcvalue(ieg1v(is),is)
     $        =ev(iband(1,is))-vbmoffset
         if(idso0(is) .gt. 0) calcvalue(idso0(is),is)
     $        =ev(ibandvbm)-ev(ibandvbm_2)

         endif     ! for iproj.eq.1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         itest=iex1v(is)+iex3v(is)+iex5v(is)+iex1c(is)+iex3c(is)
         itest=itest+iel1v(is)+iel2v(is)+iel3v(is)+iel1c(is)+
     &               idso1(is)+idso2(is)
         itest=itest+img1c(is)+imghh100(is)+imglh100(is)+
     &              imghh111(is)
         itest=itest+idpa(is)+idpag1c(is)+idpag15v(is)+
     &               idpaGX(is)+iegp3pc(is)+iegm3pc(is)+idpaGL(is)
         itest=itest+imx1cl(is)+imx1ct(is)+iml1cl(is)+iml1ct(is)
         itest=itest+idso0p(is)+idcf0(is)+isscb1(is)+issvb1(is)+
     &               img1c100(is)+img1c001(is)+imghh001(is)+
     &               imglh001(is)+ilgap(is)+il57(is)+il46(is)
         itest=itest+isscb1(is)+issvb1(is)+iengvb1(is)+
     &               iengvb2(is)+iengvb3(is)+iengvb4(is)+
     &               iengcb1(is)+iengcb2(is)+iengcb3(is)+
     &               iengcb4(is)
         do i=1,mvbo
         itest=itest+ivbo(i,is)
         enddo

         if(itest.le.0) goto 999
c ---------------------------------------------------
c X-point energy levels
c ---------------------------------------------------
         kp(1)=0.5d0     ! X point
         kp(2)=0.5d0     ! kp(i) are in reciprocal coordinates
         kp(3)=0.0d0     !
         call pwk (kp,vr,vc,-1,ev,proj,0,zz_st)
         if(iex1v(is) .gt. 0) calcvalue(iex1v(is),is)=ev(iband(1,is))
     $        -vbmoffset
         if(iex3v(is) .gt. 0) calcvalue(iex3v(is),is)=ev(iband(2,is))
     $        -vbmoffset
         if(iex5v(is) .gt. 0) calcvalue(iex5v(is),is)=ev(iband(4,is))
     $        -vbmoffset
         if(idso2(is).gt.0) calcvalue(idso2(is),is)=ev(iband(4,is))
     $        -ev(iband(3,is))
         if(iex3c(is) .gt. 0) calcvalue(iex3c(is),is)=ev(iband(6,is))
     $        -vbmoffset
         if(idcamel(is) .gt. 0) calcvalue(idcamel(is),is)=
     $        ev(iband(6,is))-ev(iband(5,is))

         kp(1)=0.5*xpt(1,is)+0.5*xpt(2,is)  ! X-point read from fit.d
         kp(2)=0.5*xpt(1,is)+0.5*xpt(3,is)
         kp(3)=0.5*xpt(2,is)+0.5*xpt(3,is)
         call pwk (kp,vr,vc,-1,ev,proj,0,zz_st)
         if(iex1c(is) .gt. 0) calcvalue(iex1c(is),is)=ev(iband(5,is)) 
     $        -vbmoffset 
c ---------------------------------------------------
c L-point energy levels
c ---------------------------------------------------
         kp(1)=0.5d0    ! L point
         kp(2)=0.5d0    ! kp(i) not in reciprocal coordinates
         kp(3)=0.5d0    !
         call pwk (kp,vr,vc,-1,ev,proj,0,zz_st)
         if(iel1v(is) .gt. 0) calcvalue(iel1v(is),is)=ev(iband(1,is))
     $        -vbmoffset
         if(iel2v(is) .gt. 0) calcvalue(iel2v(is),is)=ev(iband(2,is))
     $        -vbmoffset
         if(iel3v(is) .gt. 0) calcvalue(iel3v(is),is)=ev(iband(4,is))
     $        -vbmoffset
         if(iel1c(is) .gt. 0) calcvalue(iel1c(is),is)=ev(iband(5,is))
     $        -vbmoffset
         if(iel3c(is) .gt. 0) calcvalue(iel3c(is),is)=ev(iband(6,is))
     $        -vbmoffset
         if(idso1(is) .gt. 0) calcvalue(idso1(is),is)=ev(iband(4,is))
     $        -ev(iband(3,is))
         if(ilgap(is) .gt. 0) calcvalue(ilgap(is),is)=ev(iband(5,is))
     $        -ev(iband(4,is))
         if(il57(is) .gt. 0) calcvalue(il57(is),is)=ev(iband(6,is))
     $        -ev(iband(4,is))
         if(il46(is) .gt. 0) calcvalue(il46(is),is)=ev(iband(5,is))
     $        -ev(iband(3,is))
     
c ---------------------------------------------------
c effective masses at Gamma-point
c ---------------------------------------------------
c********effective masses along [100] direction 
         kp(1)=0.001d0       ! kp(i) not in reciprocal coordinates
         kp(2)=0.001d0       !
         kp(3)=0.000d0       !
         call pwk (kp,vr,vc,-1,ev,proj,0,zz_st)

         dltk=0.001*4.0d0*pi/alat

         temp=ev(iband(5,is))-vbmoffset
         if(img1c(is) .gt. 0)  
     &        calcvalue(img1c(is),is)=dltk**2/Ryd/
     $        (temp-calcvalue(ieg1c(is),is))
         if(img1c100(is) .gt. 0)
     $        calcvalue(img1c100(is),is)=dltk**2/Ryd/
     $        (temp-calcvalue(ieg1c(is),is))

         temp=ev(iband(4,is))-vbmoffset
         if(imghh100(is) .gt. 0)
     $        calcvalue(imghh100(is),is)=-dltk**2/Ryd/temp

ccccccccc The definition of lh is different for SO and NSO
         if(iSOps.eq.0) temp=ev(iband(2,is))-vbmoffset
         if(iSOps.eq.1) temp=ev(iband(3,is))-vbmoffset
         if(imglh100(is) .gt. 0)
     $        calcvalue(imglh100(is),is)=-dltk**2/Ryd/temp

c*******effective masses along [001] 
c for structure with cyrstal field
         if(img1c001(is)+imghh001(is)+imglh001(is).ge.1)then
         kpr(1)=0.0d0
         kpr(2)=0.0d0
         kpr(3)=0.001d0
         kp(1)=kpr(1)*a0(1,1,is)+kpr(2)*a0(1,2,is)+kpr(3)*a0(1,3,is)
         kp(2)=kpr(1)*a0(2,1,is)+kpr(2)*a0(2,2,is)+kpr(3)*a0(2,3,is)
         kp(3)=kpr(1)*a0(3,1,is)+kpr(2)*a0(3,2,is)+kpr(3)*a0(3,3,is)
         call pwk (kp,vr,vc,-1,ev,proj,0,zz_st)

         dltk=0.001*2.0d0*pi/alat

         temp=(ev(iband(5,is))+ev(iband(5,is)-1))/2-vbmoffset
         if(img1c001(is) .gt. 0)
     &        calcvalue(img1c001(is),is)=dltk**2/Ryd/
     $        (temp-calcvalue(ieg1c(is),is))

         temp=(ev(iband(4,is))+ev(iband(4,is)-1))/2-vbmoffset
         if(imghh001(is) .gt. 0)
     $        calcvalue(imghh001(is),is)=-dltk**2/Ryd/temp

         temp=(ev(iband(3,is))+ev(iband(3,is)-1))/2-vbm_1offset
         if(imglh001(is) .gt. 0)
     $        calcvalue(imglh001(is),is)=-dltk**2/Ryd/temp

         endif  ! along [001]

C******heavy hole along [111]
         kp(1)=0.003d0
         kp(2)=0.003d0
         kp(3)=0.003d0
         call pwk (kp,vr,vc,-1,ev,proj,0,zz_st)

         dltk=0.003*2.0d0*pi*dsqrt(3.0d0)/alat

         if(iSOps.eq.0) then
         temp=ev(iband(4,is))-vbmoffset
         endif

         if(iSOps.eq.1) then
ccccc for (111) direction, the hh is splitted into two lines
ccccc and the lh is not splitted.
         ist=iband(4,is)
         if(dabs(ev(ist-1)-ev(ist-2)).gt.1.D-7) then
cccccc ev(ist-1) is the other hh 
         temp=(ev(ist)+ev(ist-1))/2-vbmoffset
         else
cccccc ev(ist-1) and ev(ist-2) are degenerate, thus, they are the lh
cccccc the other hh is ev(ist-3)
         temp=(ev(ist)+ev(ist-3))/2-vbmoffset
         endif
         endif

         if(imghh111(is) .gt. 0)
     $        calcvalue(imghh111(is),is)=-dltk**2/Ryd/temp

ccccccccc The definition of lh is different for SO and NSO
         if(iSOps.eq.0) temp=ev(iband(2,is))-vbm_1offset
         if(iSOps.eq.1) temp=ev(iband(3,is))-vbm_1offset
         if(imglh111(is) .gt. 0)
     $        calcvalue(imglh111(is),is)=-dltk**2/Ryd/temp

c*******effective masses along [110] 
       if(imglh110(is)+imghh110(is) .gt.0) then
       dltkx=0.003d0
       dltky=0.003d0
       dltkz=0.d0
       kpr(1)=dltkx
       kpr(2)=dltky
       kpr(3)=dltkz
       kp(1)=kpr(1)*a0(1,1,is)+kpr(2)*a0(1,2,is)+kpr(3)*a0(1,3,is)
       kp(2)=kpr(1)*a0(2,1,is)+kpr(2)*a0(2,2,is)+kpr(3)*a0(2,3,is)
       kp(3)=kpr(1)*a0(3,1,is)+kpr(2)*a0(3,2,is)+kpr(3)*a0(3,3,is)

       call pwk (kp,vr,vc,-1,ev,proj,0,zz_st)

       dltk=dsqrt((dltkx**2+dltky**2+dltkz**2))*2.d0*pi/alat

       ist=iband(4,is)
       if(iSOps.eq.0) temp=ev(ist)-vbmoffset
       if(iSOps.eq.1) temp=(ev(ist-1)+ev(ist))/2.d0-vbmoffset
       if(imghh110(is) .gt. 0)
     $       calcvalue(imghh110(is),is)=-dltk**2/Ryd/temp

       ist=iband(3,is)
       if(iSOps.eq.0) temp=ev(ist)-vbm_1offset
       if(iSOps.eq.1) temp=(ev(ist-1)+ev(ist))/2.d0-vbm_1offset
       if(imglh110(is) .gt. 0)
     $       calcvalue(imglh110(is),is)=-dltk**2/Ryd/temp
       endif

c ---------------------------------------------------
c effective masses at X-point
c ---------------------------------------------------
c* longitudinal direction 
       if(iex1c(is).gt.0) then
       if(imx1cl(is).gt.0) then
       dltkx=0.d0
       dltky=0.d0
       dltkz=0.d0
       if(abs(xpt(1,is)).gt.0.1) dltkx=0.003d0
       if(abs(xpt(2,is)).gt.0.1) dltky=0.003d0
       if(abs(xpt(3,is)).gt.0.1) dltkz=0.003d0

       kpr(1)=xpt(1,is)+dltkx
       kpr(2)=xpt(2,is)+dltky
       kpr(3)=xpt(3,is)+dltkz
       kp(1)=kpr(1)*a0(1,1,is)+kpr(2)*a0(1,2,is)+kpr(3)*a0(1,3,is)
       kp(2)=kpr(1)*a0(2,1,is)+kpr(2)*a0(2,2,is)+kpr(3)*a0(2,3,is)
       kp(3)=kpr(1)*a0(3,1,is)+kpr(2)*a0(3,2,is)+kpr(3)*a0(3,3,is)

       call pwk (kp,vr,vc,-1,ev,proj,0,zz_st)
       
       dltk=(dltkx+dltky+dltkz)*2.d0*pi/alat
       temp=ev(iband(5,is))-vbmoffset-calcvalue(iex1c(is),is)
       calcvalue(imx1cl(is),is)=dltk**2/Ryd/temp
       endif

c* transveral direction
       if(imx1ct(is).gt.0)then
       dltkx=0.d0
       dltky=0.d0
       dltkz=0.d0
       if(xpt(1,is).gt.0.1) dltky=0.003d0
       if(xpt(2,is).gt.0.1) dltkz=0.003d0
       if(xpt(3,is).gt.0.1) dltkx=0.003d0
       
       kpr(1)=xpt(1,is)+dltkx
       kpr(2)=xpt(2,is)+dltky
       kpr(3)=xpt(3,is)+dltkz
       kp(1)=kpr(1)*a0(1,1,is)+kpr(2)*a0(1,2,is)+kpr(3)*a0(1,3,is)
       kp(2)=kpr(1)*a0(2,1,is)+kpr(2)*a0(2,2,is)+kpr(3)*a0(2,3,is)
       kp(3)=kpr(1)*a0(3,1,is)+kpr(2)*a0(3,2,is)+kpr(3)*a0(3,3,is)

       call pwk (kp,vr,vc,-1,ev,proj,0,zz_st)

       dltk=(dltkx+dltky+dltkz)*2.d0*pi/alat
       ist=iband(5,is)
       temp=(ev(ist-1)+ev(ist))/2.d0-vbmoffset-calcvalue(iex1c(is),is)
       calcvalue(imx1ct(is),is)=dltk**2/Ryd/temp
       endif
       endif   

c ---------------------------------------------------
c effective masses at L-point
c ---------------------------------------------------
c* longitudinal direction
       if(iml1cl(is)+iml3vl(is) .gt. 0) then
       dltkx=0.003d0
       dltky=dltkx
       dltkz=dltkx
       kpr(1)=0.5+dltkx
       kpr(2)=0.5+dltky
       kpr(3)=0.5+dltkz
       kp(1)=kpr(1)*a0(1,1,is)+kpr(2)*a0(1,2,is)+kpr(3)*a0(1,3,is)
       kp(2)=kpr(1)*a0(2,1,is)+kpr(2)*a0(2,2,is)+kpr(3)*a0(2,3,is)
       kp(3)=kpr(1)*a0(3,1,is)+kpr(2)*a0(3,2,is)+kpr(3)*a0(3,3,is)

       call pwk (kp,vr,vc,-1,ev,proj,0,zz_st)
       dltk=dsqrt(dltkx**2+dltky**2+dltkz**2)*2.0d0*pi/alat

       if(iel1c(is) .eq.0) then
         write(6,*) "Warning when calculate ml1cl: iel1c(is),is=",
     &             iel1c(is),is
       endif
       ist=iband(5,is)
       temp=(ev(ist)+ev(ist-1))/2.d0-vbmoffset-calcvalue(iel1c(is),is)
       calcvalue(iml1cl(is),is)=dltk**2/Ryd/temp

       if(iel3v(is) .eq.0) then
         write(6,*) "Warning when calculate ml3vl: iel3v(is),is=",
     &            iel3v(is),is
       endif
       ist=iband(4,is)
       temp=(ev(ist)+ev(ist-1))/2.d0-vbmoffset-calcvalue(iel3v(is),is)
       calcvalue(iml3vl(is),is)=-dltk**2/Ryd/temp

       endif    !  if(iml1cl(is)+iml3vl(is) .gt. 0)

       if(iml1ct(is)+iml3vt(is) .gt. 0) then
       dltkx=0.005d0
       dltky=dltkx
       dltkz=-2*dltkx
       kpr(1)=0.5+dltkx
       kpr(2)=0.5+dltky
       kpr(3)=0.5+dltkz
       kp(1)=kpr(1)*a0(1,1,is)+kpr(2)*a0(1,2,is)+kpr(3)*a0(1,3,is)
       kp(2)=kpr(1)*a0(2,1,is)+kpr(2)*a0(2,2,is)+kpr(3)*a0(2,3,is)
       kp(3)=kpr(1)*a0(3,1,is)+kpr(2)*a0(3,2,is)+kpr(3)*a0(3,3,is)

       call pwk (kp,vr,vc,-1,ev,proj,0,zz_st)
       dltk=dsqrt(dltkx**2+dltky**2+dltkz**2)*2.0d0*pi/alat

       if(iel1c(is) .eq.0) then
         write(6,*) "Warning when calculate ml1ct: iel1c(is),is=",
     &             iel1c(is),is
       endif
       ist=iband(5,is)
       temp=(ev(ist)+ev(ist-1))/2.d0-vbmoffset-calcvalue(iel1c(is),is)
       calcvalue(iml1ct(is),is)=dltk**2/Ryd/temp

       if(iel3v(is) .eq.0) then
         write(6,*) "Warning when calculate ml3vt: iel3v(is),is=",
     &             iel3v(is),is
       endif
       ist=iband(4,is)
       temp=(ev(ist)+ev(ist-1))/2.d0-vbmoffset-calcvalue(iel3v(is),is)
       calcvalue(iml3vt(is),is)=-dltk**2/Ryd/temp

       endif     ! if(iml1ct(is)+iml3vt(is) .gt. 0)

c ---------------------------------------------------
c spin-splitting at certain k-points
c ---------------------------------------------------
      if(iSOps .eq. 1 .and. (isscb1(is)+issvb1(is)).gt.0) then     
        if(isscb1(is).gt.0) calcvalue(isscb1(is),is)=0.d0
        if(issvb1(is).gt.0) calcvalue(issvb1(is),is)=0.d0
        nkpt_cb1=0
        nkpt_vb1=0
        filess=trim(inputfile(is))//"_SS"
        open (unit=4,file=filess,status='old')
        read(4,*) nkpt
        do 777 ikpt=1,nkpt
        read(4,*) kpr(1:3), sscb1,ssvb1        !k-point units in 2\pi/a
        sscb1=sscb1/Ryd
        ssvb1=ssvb1/Ryd
        kp(1)=kpr(1)*a0(1,1,is)+kpr(2)*a0(1,2,is)+kpr(3)*a0(1,3,is) 
        kp(2)=kpr(1)*a0(2,1,is)+kpr(2)*a0(2,2,is)+kpr(3)*a0(2,3,is)
        kp(3)=kpr(1)*a0(3,1,is)+kpr(2)*a0(3,2,is)+kpr(3)*a0(3,3,is)
        call pwk (kp,vr,vc,-1,ev,proj,0,zz_st)

        if(isscb1(is).gt.0 .and.sscb1.gt.1.d-6) then
          nkpt_cb1=nkpt_cb1+1
          calcvalue(isscb1(is),is)=calcvalue(isscb1(is),is)+
     $    ((ev(iband(5,is))-ev(iband(5,is)-1)-sscb1))**2
        endif

        if(issvb1(is).gt.0 .and. ssvb1.gt.1.d-6) then
          nkpt_vb1=nkpt_vb1+1
          calcvalue(issvb1(is),is)=calcvalue(issvb1(is),is)+
     $    ((ev(iband(4,is))-ev(iband(4,is)-1)-ssvb1))**2
        endif
777     continue 
        close(4)

      if(isscb1(is).gt.0) calcvalue(isscb1(is),is)=
     $                    dsqrt(calcvalue(isscb1(is),is)/nkpt_cb1)
      if(issvb1(is).gt.0) calcvalue(issvb1(is),is)=
     $                    dsqrt(calcvalue(issvb1(is),is)/nkpt_vb1)

      endif
c ---------------------------------------------------
c GW band structure in whole Brillouin zone
c ---------------------------------------------------
      itest=issvb1(is)+iengvb1(is)+iengvb2(is)+iengvb3(is)+
     &      iengvb4(is)+iengcb1(is)+iengcb2(is)+iengcb3(is)+
     &      iengcb4(is)

      if(itest.gt.0) then

        if(iengcb1(is).gt.0) calcvalue(iengcb1(is),is)=0.d0
        if(iengcb2(is).gt.0) calcvalue(iengcb2(is),is)=0.d0
        if(iengcb3(is).gt.0) calcvalue(iengcb3(is),is)=0.d0
        if(iengcb4(is).gt.0) calcvalue(iengcb4(is),is)=0.d0
        if(iengvb1(is).gt.0) calcvalue(iengvb1(is),is)=0.d0
        if(iengvb2(is).gt.0) calcvalue(iengvb2(is),is)=0.d0
        if(iengvb3(is).gt.0) calcvalue(iengvb3(is),is)=0.d0
        if(iengvb4(is).gt.0) calcvalue(iengvb4(is),is)=0.d0

        if(iSOPs .eq. 0) vbm_shift=the_data(iengvb1(is),is)%sovalue  
        if(iSOPs .eq. 1) vbm_shift=the_data(iengvb1(is),is)%nsovalue  
        vbm_shift=vbm_shift/Ryd

        filekpt=trim(inputfile(is))//"_kpt"
        fileband=trim(inputfile(is))//"_GW"
        open (unit=11,file=filekpt,status='old')
        open (unit=12,file=fileband,status='old')
        read(11,*) nkpt
        read(12,*)                 ! comment line

        do 888 ikpt=1,nkpt
        read(11,*) indx,kpr(1:3)
        read(12,*) evGW(1:8)       !k-point units in 2\pi/a
        evGW=evGW/Ryd

        kp(1)=kpr(1)*a0(1,1,is)+kpr(2)*a0(1,2,is)+kpr(3)*a0(1,3,is)
        kp(2)=kpr(1)*a0(2,1,is)+kpr(2)*a0(2,2,is)+kpr(3)*a0(2,3,is)
        kp(3)=kpr(1)*a0(3,1,is)+kpr(2)*a0(3,2,is)+kpr(3)*a0(3,3,is)

c       write(6,'(i3,3f8.3)') ikpt,kp(1:3)

        call pwk (kp,vr,vc,-1,ev,proj,0,zz_st)

        if(iengvb1(is).gt.0) 
     &     calcvalue(iengvb1(is),is)=calcvalue(iengvb1(is),is)+
     $     (ev(iband(4,is))-evGW(4)-vbm_shift)**2
        if(iengvb2(is).gt.0) 
     &     calcvalue(iengvb2(is),is)=calcvalue(iengvb2(is),is)+
     $     (ev(iband(3,is))-evGW(3)-vbm_shift)**2
        if(iengvb3(is).gt.0) 
     &     calcvalue(iengvb3(is),is)=calcvalue(iengvb3(is),is)+
     $     (ev(iband(2,is))-evGW(2)-vbm_shift)**2
        if(iengvb4(is).gt.0) 
     &     calcvalue(iengvb4(is),is)=calcvalue(iengvb4(is),is)+
     $     (ev(iband(1,is))-evGW(1)-vbm_shift)**2
        if(iengcb1(is).gt.0) 
     &     calcvalue(iengcb1(is),is)=calcvalue(iengcb1(is),is)+
     $     (ev(iband(5,is))-evGW(5)-vbm_shift)**2
        if(iengcb2(is).gt.0) 
     &     calcvalue(iengcb2(is),is)=calcvalue(iengcb2(is),is)+
     $     (ev(iband(6,is))-evGW(6)-vbm_shift)**2
        if(iengcb3(is).gt.0) 
     &     calcvalue(iengcb3(is),is)=calcvalue(iengcb3(is),is)+
     $     (ev(iband(7,is))-evGW(7)-vbm_shift)**2
        if(iengcb4(is).gt.0) 
     &     calcvalue(iengcb4(is),is)=calcvalue(iengcb4(is),is)+
     $     (ev(iband(8,is))-evGW(8)-vbm_shift)**2

888     continue
        close(11)
        close(12)

      if(iengcb1(is).gt.0) calcvalue(iengcb1(is),is)=
     &    dsqrt(calcvalue(iengcb1(is),is)/nkpt)
      if(iengcb2(is).gt.0) calcvalue(iengcb2(is),is)=
     &    dsqrt(calcvalue(iengcb2(is),is)/nkpt)
      if(iengcb3(is).gt.0) calcvalue(iengcb3(is),is)=
     &    dsqrt(calcvalue(iengcb3(is),is)/nkpt)
      if(iengcb4(is).gt.0) calcvalue(iengcb4(is),is)=
     &    dsqrt(calcvalue(iengcb4(is),is)/nkpt)
      if(iengvb1(is).gt.0) calcvalue(iengvb1(is),is)=
     &    dsqrt(calcvalue(iengvb1(is),is)/nkpt)
      if(iengvb2(is).gt.0) calcvalue(iengvb2(is),is)=
     &    dsqrt(calcvalue(iengvb2(is),is)/nkpt)
      if(iengvb3(is).gt.0) calcvalue(iengvb3(is),is)=
     &    dsqrt(calcvalue(iengvb3(is),is)/nkpt)
      if(iengvb4(is).gt.0) calcvalue(iengvb4(is),is)=
     &    dsqrt(calcvalue(iengvb4(is),is)/nkpt)

      endif      ! if(itest.gt.0) then
c ---------------------------------------------------------------
c     Now change the lattice constant uniformaly and re-calculate the
c     bandstructure, to obtain the deformation potentials.
c     There is no need to move the atomic positions as vcell will scale them
c     according to alat as we want.
c ---------------------------------------------------------------
C=====now, time to deform the crystal, and change the psf(cation)
C=====
         dvol = 0.001

         alat=alat0(is)*(1.0d0+dvol)
c     dv/v       hydrofac = dvol*(3.0d0+dvol*(3.0d0+dvol))
         hydrofac = dvol*3.0d0

         do ips=1,itotps
            psf(ips)=hydrofac*psfvol(ips)
         enddo

         a(:,:)=a0(:,:,is)*alat

         call gridg ()
         call vcell (vr,vc,-1)

         kp(1)=0.0d0
         kp(2)=0.0d0
         kp(3)=0.0d0
         call pwk (kp,vr,vc,-1,ev,proj,0,zz_st)
         temp1=ev(iband(4,is))
         temp2=ev(iband(5,is))
         vbmstrn=ev(iband(4,is))

         vbm_shift = ev(iband(4,is))-calcvalue(ievbm(is),is)
         cbm_shift = ev(ig1c(is))
     $        -(calcvalue(ievbm(is),is)+calcvalue(ieg1c(is),is))

         if(idpa(is) .gt. 0) calcvalue(idpa(is),is)
     $        =(ev(ig1c(is))-ev(iband(4,is))
     &        -calcvalue(ieg1c(is),is))/hydrofac

         if(idpag1c(is) .gt. 0) 
     $        calcvalue(idpag1c(is),is)=cbm_shift/hydrofac

         if(idpag15v(is) .gt. 0) 
     &        calcvalue(idpag15v(is),is)=vbm_shift/hydrofac

ccc dpaX1c: the hydrostatic deformation potential of X_1c-Gamma_15v gap.
c         kp(1)=0.5d0
c         kp(2)=0.5d0
c         kp(3)=0.0d0
         kp(1)=0.5*xpt(1,is)+0.5*xpt(2,is)  ! X-point read from fit.d
         kp(2)=0.5*xpt(1,is)+0.5*xpt(3,is)
         kp(3)=0.5*xpt(2,is)+0.5*xpt(3,is)
        
         call pwk (kp,vr,vc,-1,ev,proj,0,zz_st)
         temp2=ev(iband(5,is))

         if(idpaGX(is)>0) calcvalue(idpaGX(is),is)
     $        =(temp2-vbmstrn-calcvalue(ieX1c(is),is))/hydrofac

c*jwluo* dpaL: the hydrostatic deformation potential of L_1c-Gamma_15v gap.
         if(idpaGL(is)>0) then
         kp(1)=0.5d0
         kp(2)=0.5d0
         kp(3)=0.5d0
         call pwk (kp,vr,vc,-1,ev,proj,0,zz_st)
         temp2=ev(iband(5,is))

         if(iel1c(is)>0) calcvalue(idpaGL(is),is)
     &       =(temp2-vbmstrn-calcvalue(iel1c(is),is))/hydrofac

         endif

c*********Deform the lattice constant by +- 3 percent

         do icase=1,2
            if((icase.eq.1.and.iegp3pc(is)>0).or.(icase.eq.2.and.
     &           iegm3pc(is)>0)) then 
               if(icase.eq.1) dvol=0.03
               if(icase.eq.2) dvol=-0.03

               alat=alat0(is)*(1.0d0+dvol)
               hydrofac = dvol*(3.0d0+dvol*(3.0d0+dvol))

               do ips=1,itotps
                  psf(ips)=hydrofac*psfvol(ips)
               enddo

               a=a0(:,:,is)*alat

               call gridg ()
               call vcell (vr,vc,-1)

               kp(1)=0.0d0
               kp(2)=0.0d0
               kp(3)=0.0d0
               call pwk (kp,vr,vc,-1,ev,proj,0,zz_st)

               if(icase.eq.1) calcvalue(iegp3pc(is),is)
     $              =ev(iband(5,is))-ev(iband(4,is))
               if(icase.eq.2) calcvalue(iegm3pc(is),is)
     $              =ev(iband(5,is))-ev(iband(4,is))

            end if
         enddo
c -------------------------------------------------------------------------         
c     Now change the lattice constant non-uniformaly and re-calculate the
c     bandstructure, to obtain the epitaxial deformation potentials.
c     In the x-y plane the lattice constant is decreased by 1% and in the
c     z direction the lattice constant is increased by 1%.  The value is fitted
c     to the result of Eq.(8) in PRB 49, p14337.
c     The crystal field splits the 3-fold degeneracy with 2 bands 1/3 above the
c     original value and 1 band 2/3 below the original value.
C --------------------------------------------------------------------------
C==   as you see, e_xx is calculated so as to DV/V = 0
C==   
         alat=alat0(is)
         e_xx = 1.0d0/dsqrt(1.0d0+e_zz)-1.0d0
         hydrofac=2.0d0*e_xx+e_zz
         do ips=1,itotps
            psf(ips)=hydrofac*psfvol(ips)
         enddo

         a(1,1)=a0(1,1,is)*(1+e_xx)*alat
         a(1,2)=a0(1,2,is)*(1+e_xx)*alat
         a(1,3)=a0(1,3,is)*(1+e_zz)*alat
         a(2,1)=a0(2,1,is)*(1+e_xx)*alat
         a(2,2)=a0(2,2,is)*(1+e_xx)*alat
         a(2,3)=a0(2,3,is)*(1+e_zz)*alat
         a(3,1)=a0(3,1,is)*(1+e_xx)*alat
         a(3,2)=a0(3,2,is)*(1+e_xx)*alat
         a(3,3)=a0(3,3,is)*(1+e_zz)*alat

c     As well as rescaling the lattice vectors, have to move atomic positions
c     The In stays at the origin and the P moves according to the strain.
ccc   tau(1,:,is), etc, are the Cartesian Coord in unit of alat, not the
ccc   supercell edge coord.

         tau(1,:,is) = tau0(1,:,is)*(1+e_xx)
         tau(2,:,is) = tau0(2,:,is)*(1+e_xx)
         tau(3,:,is) = tau0(3,:,is)*(1+e_zz)

         call gridg ()
         call vcell (vr,vc,-1)

         kp(1)=0.0d0
         kp(2)=0.0d0
         kp(3)=0.0d0

         call pwk(kp,vr,vc,-1,ev,proj,0,zz_st)

C     = deformation potential under [001] strain PRB 49, 14337 (1994) Eq.25
         if(idpb100(is) .gt. 0) then

            if(iSOps .eq. 0)      
     &           calcvalue(idpb100(is),is)=
     $           -(ev(iband(4,is))-ev(iband(2,is)))/3.0d0/(e_zz-e_xx)

            if(iSOps .eq. 1) then     
               if(idso0(is).le.0) then
               write(6,*) "must calc. edso0 to calc. dpb100, stop"
               stop
               endif
               temp = ev(iband(4,is))-ev(iband(3,is))+
     &              ev(iband(4,is))-ev(iband(2,is))
     $              -calcvalue(idso0(is),is)

               calcvalue(idpb100(is),is)=-temp/3.0d0/(e_zz-e_xx)
            endif               !if(SO)

         endif                  ! b(100)

c ----------------------------------------------------------------
c Fit the VBO if there are any to be fitted
C now, again, add the hydrostatic term        
c ----------------------------------------------------------------
         do i=1,mvbo
            if (ivbo(i,is).gt.0) then 

               alat = alat0(is)
               hydrofac=(1.0d0+eps_paral(i,is))**2
     $              *(1.0d0+eps_perp(i,is))-1

               do ips=1,itotps
                  psf(ips)=hydrofac*psfvol(ips)
               enddo

               a(1,1)=a0(1,1,is)*alat*(1.0d0+eps_paral(i,is))
               a(1,2)=a0(1,2,is)*alat*(1.0d0+eps_paral(i,is))
               a(1,3)=a0(1,3,is)*alat*(1.0d0+eps_perp(i,is))
               a(2,1)=a0(2,1,is)*alat*(1.0d0+eps_paral(i,is))
               a(2,2)=a0(2,2,is)*alat*(1.0d0+eps_paral(i,is))
               a(2,3)=a0(2,3,is)*alat*(1.0d0+eps_perp(i,is))
               a(3,1)=a0(3,1,is)*alat*(1.0d0+eps_paral(i,is))
               a(3,2)=a0(3,2,is)*alat*(1.0d0+eps_paral(i,is))
               a(3,3)=a0(3,3,is)*alat*(1.0d0+eps_perp(i,is))
               
c     Move the atomic positions
               
               tau(1,:,is) = tau0(1,:,is)*(1.0d0+eps_paral(i,is))
               tau(2,:,is) = tau0(2,:,is)*(1.0d0+eps_paral(i,is))
               tau(3,:,is) = tau0(3,:,is)*(1.0d0+eps_perp(i,is))
               
               call gridg ()
               call vcell (vr,vc,-1)
               
               kp(1)=0.0d0
               kp(2)=0.0d0
               kp(3)=0.0d0
               call pwk(kp,vr,vc,-1,ev,proj,0,zz_st)
               calcvalue(ivbo(i,is),is)=ev(iband(4,is))-vbovbm(i,is)/Ryd

            end if
            
         end do


999   continue

      end do                    ! Loop over structures
c ------------------------------------------------------------------      
c     Now set up the residuals whose sum of squares is to be minimised
c     eveything is in eV, now
c ------------------------------------------------------------------      
      calcvalue=calcvalue*Ryd
      
cc      write(6,*) "nstruct=",nstruct
      sum=0.d0
      ii=0
      do is=1,nstruct
         if(iSOps .eq. 0) then
            do i=1,nsprop(is)
               ii=ii+1
               fvecc(ii)=the_data(i,is)%weight*
     &              (calcvalue(i,is)-the_data(i,is)%nsovalue)/
     $              the_data(i,is)%nsovalue
            enddo
         else
            do i=1,nsprop(is)
               ii=ii+1
               fvecc(ii)=the_data(i,is)%weight*
     &              (calcvalue(i,is)-the_data(i,is)%sovalue)/
     $              the_data(i,is)%sovalue
               sum=sum+fvecc(ii)**2
c -----------------------------------------------------------------
cccc test, test
c         write(6,800) is,i,the_data(i,is)%sovalue,calcvalue(i,is),
c     &         fvecc(ii)**2,the_data(i,is)%weight
c800      format("test", 2(i3,1x),f13.8,f13.8,f13.8,f13.8)
cccc test, test
c -----------------------------------------------------------------
            enddo
         endif
      end do

      
      end

