      subroutine funxsqs ( )
      
      implicit none

      include 'epm.inc' 

      integer k,i,j,itest
      double precision e_xx,e_zz,dvol,hydrofac
      double precision kp(1:3),ev(1:mnpw*2),vr(1:mngrid)
      double precision pi,vbmoffset,temp,temp1,temp2,dltk,sum
      complex*16 psi,vc(1:mngrid)
      real*8 proj(mnpw*2,mnpw*2)
      complex*16 zz_st(mnpw*2,mnpw*2,mns)
      integer iproj(mns),ipcbm(mns)
      integer ibandcbm,ibandvbm,ibandvbm_1,ibandvbm_2
      real*8 s,smax
      integer j1
ccc added by jwluo
      double precision blat(3,3) ! reciprocal lattice and its inverse
      double precision tmp1(3)

      common  /comproj/zz_st,iproj,ipcbm


      e_zz = 0.001d0

      pi=4.0d0*atan(1.0d0)

      k=0


C==== for the ideal zincblend structure, psf = 0        
         psf(1:itotps)=0.0d0
         
         alat=alat0(is)
         a(:,:)=a0(:,:,is)*alat
         tau(:,:,is)=tau0(:,:,is)

         call gridg ()
         call vcell (vr,vc,-1)

ccccccccc arbitrary
         tmp1(1)=0.1
         tmp1(2)=0.2
         tmp1(3)=0.3

ccccccccc inverse  reciprocal lattices blat after gaussj
         blat(:,1)=b1(:)
         blat(:,2)=b2(:)
         blat(:,3)=b3(:)
         call gaussj(blat,3,3,tmp1,1,1)
         
ccccccccccccc Gamma point cccccccccccccccccccccccccccccccc
         kp(1)=0.0d0    !Gamma point
         kp(2)=0.0d0    !
         kp(3)=0.0d0    !

         call pwk(kp,vr,vc,-1,ev,proj,iproj(is),zz_st(1,1,is))

         if(iproj(is).eq.0) then
         vbmoffset=ev(iband(4,is))
         if(ievbm(is) .gt. 0) calcvalue(ievbm(is),is)
     $        =ev(iband(4,is))
         if(ievbm_1(is) .gt. 0) calcvalue(ievbm_1(is),is)
     $        =ev(iband(3,is))-vbmoffset
         if(ieg1c(is) .gt. 0) calcvalue(ieg1c(is),is)
     $        =ev(iband(5,is))-vbmoffset
         if(ieg15c(is) .gt. 0)calcvalue(ieg15c(is),is)
     $        =ev(iband(6,is))-vbmoffset
         if(ieg1v(is) .gt. 0) calcvalue(ieg1v(is),is)
     $        =ev(iband(1,is))-vbmoffset
         if(idso0(is) .gt. 0) calcvalue(idso0(is),is)
     $        =ev(iband(4,is))-ev(iband(2,is))
         endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         itest=iex1v(is)+iex3v(is)+iex5v(is)+iex1c(is)+iex3c(is)
         itest=itest+iel1v(is)+iel2v(is)+iel3v(is)+iel1c(is)+
     &               idso1(is)+idso2(is)
         itest=itest+img1c(is)+imghh100(is)+imglh100(is)+
     &              imghh111(is)
         itest=itest+idpa(is)+idpag1c(is)+idpag15v(is)+
     &               idpaGX(is)+iegp3pc(is)+iegm3pc(is)+idpaGL(is)
         itest=itest+imx1cl(is)+imx1ct(is)+iml1cl(is)+iml1ct(is)
         do i=1,mvbo
         itest=itest+ivbo(i,is)
         enddo

ccccccccc X-point ccccccccccccccccccccccccccccccccccccccccc
         if(iex1v(is).gt.0 .or. iex3v(is).gt.0 .or.
     1      iex5v(is).gt.0 .or. idso2(is).gt.0 .or.
     2      iex1c(is).gt.0 .or. iex3c(is).gt.0) then

         tmp1(1)=1.d0*2.0*pi/alat  ! X point in Cartesian
         tmp1(2)=0.d0
         tmp1(3)=0.d0
         

         kp(1)=tmp1(1)*blat(1,1)+tmp1(2)*blat(1,2)
     &        +tmp1(3)*blat(1,3)
         kp(2)=tmp1(1)*blat(2,1)+tmp1(2)*blat(2,2)
     &        +tmp1(3)*blat(2,3)
         kp(3)=tmp1(1)*blat(3,1)+tmp1(2)*blat(3,2)
     &        +tmp1(3)*blat(3,3)

         call pwk (kp,vr,vc,-1,ev,proj,0,zz_st)
         if(iex1v(is) .gt. 0) calcvalue(iex1v(is),is)=ev(iband(1,is))
     $        -vbmoffset
         if(iex3v(is) .gt. 0) calcvalue(iex3v(is),is)=ev(iband(2,is))
     $        -vbmoffset
         if(iex5v(is) .gt. 0) calcvalue(iex5v(is),is)=ev(iband(4,is))
     $        -vbmoffset
         if(idso2(is).gt.0) calcvalue(idso2(is),is)=ev(iband(4,is))
     $        -ev(iband(3,is))
         if(iex1c(is) .gt. 0) calcvalue(iex1c(is),is)=ev(iband(5,is))
     $        -vbmoffset
         if(iex3c(is) .gt. 0) calcvalue(iex3c(is),is)=ev(iband(6,is))
     $        -vbmoffset

         endif  !if(iex1v(is).gt.0 .or. iex3v(is).gt.0 .or.

         if((xpt(1,is)+xpt(2,is)+xpt(3,is)).ne.1.d0) then
         if(iex1c(is).gt.0 .or. iex3c(is).gt.0) then
         tmp1(1)=xpt(1,is)*2.0*pi/alat   ! X-point in Cartesian
         tmp1(2)=xpt(2,is)*2.0*pi/alat
         tmp1(3)=xpt(3,is)*2.0*pi/alat
         
c change X-point from cartesian to reciprocal lattices
         kp(1)=tmp1(1)*blat(1,1)+tmp1(2)*blat(1,2)
     &        +tmp1(3)*blat(1,3)
         kp(2)=tmp1(1)*blat(2,1)+tmp1(2)*blat(2,2)
     &        +tmp1(3)*blat(2,3)
         kp(3)=tmp1(1)*blat(3,1)+tmp1(2)*blat(3,2)
     &        +tmp1(3)*blat(3,3)

         call pwk (kp,vr,vc,-1,ev,proj,0,zz_st)
         if(iex1c(is) .gt. 0) calcvalue(iex1c(is),is)=ev(iband(5,is)) 
     $        -vbmoffset 
         if(iex3c(is) .gt. 0) calcvalue(iex3c(is),is)=ev(iband(6,is))
     $        -vbmoffset

         endif  !if(iex1c(is).gt.0 .or. iex3c(is).gt.0)
         endif  !if(xpt(1,is)+xpt(2,is)+xpt(3,is).ne.1.d0)

ccccccccccc L-point  cccccccccccccccccccccccccccccccccccccccc
         if(iel1v(is).gt.0 .or. iel2v(is).gt.0 .or.
     1      iel3v(is).gt.0 .or. iel1c(is).gt.0 .or.
     2      iel3c(is).gt.0 .or. idso1(is).gt.0 ) then
         tmp1(1)=0.5d0*2.0*pi/alat    ! L point
         tmp1(2)=0.5d0*2.0*pi/alat  
         tmp1(3)=0.5d0*2.0*pi/alat    !

c change k-point unit form cartsian to reciprocal lattices
         kp(1)=tmp1(1)*blat(1,1)+tmp1(2)*blat(1,2)
     &        +tmp1(3)*blat(1,3)
         kp(2)=tmp1(1)*blat(2,1)+tmp1(2)*blat(2,2)
     &        +tmp1(3)*blat(2,3)
         kp(3)=tmp1(1)*blat(3,1)+tmp1(2)*blat(3,2)
     &        +tmp1(3)*blat(3,3)

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
     
         endif   ! if(iel1v(is).gt.0 .or. iel2v(is).gt.0 .or.

      return
      end

