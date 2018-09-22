cc      subroutine funx (iflag,nprop,nparam,xc,fvecc,iw,liw,w,lw)
      program band_str
      
      implicit none

      include 'epm.inc' 
      !include a parameter bank
      integer mparam
      parameter (mparam=24)

      integer k,i,j,ips,writeout,iflag,nparam,nprop
      integer icase,ii,ist,itest,iii,nk
      double precision e_xx,e_zz,dvol,hydrofac
      double precision xc(mparam),fvecc(mproperty)
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

      common  /comproj/zz_st,iproj,ipcbm
      !set a common space with three parameters
      call readdata(nparam)


      e_zz = 0.001d0

      pi=4.0d0*atan(1.0d0)

      k=0
      xc=1.d0
C======================================================
C==== rather tricky here
C==== psf(cation) = (dV/V)*psf0(ips)*exp(-10*q*q)
C==== using a local variable, psflocal(ips), to get around the problem 
C==== whenever the unit cell changes, set 
C==== psf(cation)= (dV/V)*psflocal(cation)*exp(-10*q*q)
C=======================================================
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

      open(17,file="band_str.input") 
      rewind(17)
      open(11,file="band_str.out") 
      rewind(11)
      read(17,*) is,nk  ! index of the structure,num of k-points
         
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

C==== for the ideal zincblend structure, psf = 0        
         psf(1:itotps)=0.0d0
         
         alat=alat0(is)
         a(:,:)=a0(:,:,is)*alat
         tau(:,:,is)=tau0(:,:,is)

         call gridg ()
         call vcell (vr,vc,-1)

         do iii=1,nk
         read(17,*) kp(1),kp(2),kp(3)

         call pwk(kp,vr,vc,-1,ev,proj,iproj(is),zz_st(1,1,is))

        write(6,300) kp(1),kp(2),kp(3),
     &        (27.211396/2*ev(i),i=1,16)
        write(11,300) kp(1),kp(2),kp(3),
     &        (27.211396/2*ev(i),i=1,16)
300     format(3(f10.7,1x),16(f10.5,1x))
         enddo
        close(17)
        close(11)
        


      stop
      end

