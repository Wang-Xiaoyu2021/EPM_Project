      program evaluate_ps
      implicit none 
      include 'epm.inc'

      integer k,i,j,ips,writeout,iflag,nparam,lw,nprop,idpal,idpax
      double precision kp(1:3),kpc(1:3),ev(1:mnpw*2),vr(1:mngrid)
      double precision pi,vbmoffset,temp,temp1,temp2,temp3,dltk
      double precision dpbcheck,bandcenter,vbm_shift,cbm_shift
      double precision e_zz,e_xx,dvol,hydrofac,de_zz,invkpoints
      complex*16 psi,vc(1:mngrid)

      integer,parameter :: nkpoints=100
      
      invkpoints = 1.0d0/dble(nkpoints)
      idpal=0
      idpax=0
      
      call readdata(nparam)

      de_zz = 0.001d0
      e_zz = 0.001d0
      dvol = 0.001d0

      pi=4.0d0*atan(1.0d0)
      
      do 20 ips=1,nps
         psf(ips)=0.0d0*psf0(ips)
         psbeta(ips)=psbeta0(ips)
         psSO(ips) = psSO0(ips)
         do 30 j=1,ngauss(ips)
            psa(j,ips)=psa0(j,ips)
            psb(j,ips)=psb0(j,ips)
            psc(j,ips)=psc0(j,ips)
 30      continue
 20   continue

      
      alat=alat0
      a1(1)=a01(1)*alat
      a1(2)=a01(2)*alat
      a1(3)=a01(3)*alat
      a2(1)=a02(1)*alat
      a2(2)=a02(2)*alat
      a2(3)=a02(3)*alat
      a3(1)=a03(1)*alat
      a3(2)=a03(2)*alat
      a3(3)=a03(3)*alat
      
      tau(1,1) = tau0(1,1)
      tau(2,1) = tau0(2,1)
      tau(3,1) = tau0(3,1)
      tau(1,2) = tau0(1,2)
      tau(2,2) = tau0(2,2)
      tau(3,2) = tau0(3,2)

      call recvec (a1,a2,a3,b1,b2,b3,omega)
      call gridg ()
      call vcell (vr,vc,-1)

      kp(1)=0.0d0
      kp(2)=0.0d0
      kp(3)=0.0d0
      call pwk (kp,vr,vc,-1,ev)
      vbmoffset=ev(band4)
      if(ievbm .gt. 0) calcvalue(ievbm) =ev(band4)
      if(ieg1c .gt. 0) calcvalue(ieg1c) =ev(band5)-vbmoffset
      if(ieg15c .gt. 0) calcvalue(ieg15c)=ev(band6)-vbmoffset
      if(ieg1v .gt. 0) calcvalue(ieg1v) =ev(band1)-vbmoffset
      if(idso0 .gt. 0) calcvalue(idso0) =ev(band4)-ev(band2)	
      if(iSOps .eq. 1) bandcenter=vbmoffset-calcvalue(idso0)/3.0d0

      write(6,119) alat,ev(band2)*Ryd,ev(band3)*Ryd,
     &     ev(band4)*Ryd,ev(band5)*Ryd

      kp(1)=0.5d0
      kp(2)=0.5d0
      kp(3)=0.0d0
      call pwk (kp,vr,vc,-1,ev)
      if(iex1v .gt. 0) calcvalue(iex1v)=ev(band1)-vbmoffset
      if(iex3v .gt. 0) calcvalue(iex3v)=ev(band2)-vbmoffset
      if(iex5v .gt. 0) calcvalue(iex5v)=ev(band4)-vbmoffset
      if(iex1c .gt. 0) calcvalue(iex1c)=ev(band5)-vbmoffset
      if(iex3c .gt. 0) calcvalue(iex3c)=ev(band6)-vbmoffset

      kp(1)=0.5d0
      kp(2)=0.5d0
      kp(3)=0.5d0
      call pwk (kp,vr,vc,-1,ev)
      if(iel1v .gt. 0) calcvalue(iel1v)=ev(band1)-vbmoffset
      if(iel2v .gt. 0) calcvalue(iel2v)=ev(band2)-vbmoffset
      if(iel3v .gt. 0) calcvalue(iel3v)=ev(band4)-vbmoffset
      if(iel1c .gt. 0) calcvalue(iel1c)=ev(band5)-vbmoffset
      if(idso1 .gt. 0) calcvalue(idso1)=ev(band4)-ev(band3)
      
      kp(1)=0.001d0
      kp(2)=0.001d0
      kp(3)=0.000d0
      call pwk (kp,vr,vc,-1,ev)

      dltk=0.001*4.0d0*pi/alat

      temp=ev(band5)-vbmoffset
      if(img1c .gt. 0)  
     &     calcvalue(img1c)=dltk**2/Ryd/(temp-calcvalue(ieg1c))

      temp=ev(band4)-vbmoffset
      if(imghh100 .gt. 0) calcvalue(imghh100)=-dltk**2/Ryd/temp

      temp=ev(band3)-vbmoffset
      if(imglh100 .gt. 0) calcvalue(imglh100)=-dltk**2/Ryd/temp

C******heavy hole along [111]
      kp(1)=0.001d0
      kp(2)=0.001d0
      kp(3)=0.001d0
      call pwk (kp,vr,vc,-1,ev)

      dltk=0.001*2.0d0*pi*dsqrt(3.0d0)/alat

      temp=ev(band4)-vbmoffset
      if(imghh111 .gt. 0) calcvalue(imghh111)=-dltk**2/Ryd/temp

c     Now change the lattice constant uniformaly and re-calculate the
c     bandstructure, to obtain the deformation potentials.
c     There is no need to move the atomic positions as vcell will scale them
c     according to alat as we want.

      alat=alat0*(1.0d0+dvol)
C==== DV/V
C     hydrofac = dvol*(3.0d0+dvol*(3.0d0+dvol))
C==== Tr(e)
      hydrofac = 3.0d0*dvol

      do ips=1,nps
         psf(ips)=hydrofac*psfvol(ips)
      enddo

      a1(1)=a01(1)*alat
      a1(2)=a01(2)*alat
      a1(3)=a01(3)*alat
      a2(1)=a02(1)*alat
      a2(2)=a02(2)*alat
      a2(3)=a02(3)*alat
      a3(1)=a03(1)*alat
      a3(2)=a03(2)*alat
      a3(3)=a03(3)*alat


      call recvec (a1,a2,a3,b1,b2,b3,omega)
      call gridg ()
      call vcell (vr,vc,-1)

      kp(1)=0.0d0
      kp(2)=0.0d0
      kp(3)=0.0d0
      call pwk (kp,vr,vc,-1,ev)
      temp1=ev(band4)
      temp2=ev(band5)

      vbm_shift = ev(band4)-calcvalue(ievbm)
      cbm_shift = ev(band5)-(calcvalue(ievbm)+calcvalue(ieg1c))

      if(idpa .gt. 0) 
     &     calcvalue(idpa)=(temp2-temp1-calcvalue(ieg1c))/hydrofac

C=====dpaX

      kp(1)=0.5d0
      kp(2)=0.5d0
      kp(3)=0.0d0

      call pwk (kp,vr,vc,-1,ev)
      temp2=ev(band5)

      if(idpaX .gt. 0)         
     &     calcvalue(idpaX)=(temp2-temp1-calcvalue(iex1c))/hydrofac

C==== dpxL

      kp(1)=0.5d0
      kp(2)=0.5d0
      kp(3)=0.5d0

      call pwk (kp,vr,vc,-1,ev)
      temp2=ev(band5)

      if(idpaL .gt. 0) 
     &     calcvalue(idpaL)=(temp2-temp1-calcvalue(iel1c))/hydrofac

c     Now change the lattice constant non-uniformaly and re-calculate the
c     bandstructure, to obtain the epitaxial deformation potentials.
c     In the x-y plane the lattice constant is decreased by 1% and in the
c     z direction the lattice constant is increased by 1%.  The value is fitted
c     to the result of Eq.(8) in PRB 49, p14337.
c     The crystal field splits the 3-fold degeneracy with 2 bands 1/3 above the
c     original value and 1 band 2/3 below the original value.

      alat=alat0
      do ips=1,nps
         psf(ips)=0.0d0
      enddo
      e_xx = 1.0d0/dsqrt(1.0d0+e_zz)-1.0d0
      a1(1)=a01(1)*(1+e_xx)*alat
      a1(2)=a01(2)*(1+e_xx)*alat
      a1(3)=a01(3)*(1+e_zz)*alat
      a2(1)=a02(1)*(1+e_xx)*alat
      a2(2)=a02(2)*(1+e_xx)*alat
      a2(3)=a02(3)*(1+e_zz)*alat
      a3(1)=a03(1)*(1+e_xx)*alat
      a3(2)=a03(2)*(1+e_xx)*alat
      a3(3)=a03(3)*(1+e_zz)*alat

c     As well as rescaling the lattice vectors, have to move atomic positions
c     The In stays at the origin and the P moves according to the strain.

      tau(1,1) = tau0(1,1)
      tau(2,1) = tau0(2,1)
      tau(3,1) = tau0(3,1)
      tau(1,2) = tau0(1,2)*(1+e_xx)
      tau(2,2) = tau0(2,2)*(1+e_xx)
      tau(3,2) = tau0(3,2)*(1+e_zz)

      call recvec (a1,a2,a3,b1,b2,b3,omega)
      call gridg ()
      call vcell (vr,vc,-1)

      kp(1)=0.0d0
      kp(2)=0.0d0
      kp(3)=0.0d0

      call pwk(kp,vr,vc,-1,ev)

C     = deformation potential under [001] strain PRB 49, 14337 (1994) Eq.25
      if(idpb100 .gt. 0) then

         if(iSOps .eq. 0)      
     &        calcvalue(idpb100)=
     &        -(ev(band4)-ev(band2))/3.0d0/(e_zz-e_xx)

         if(iSOps .eq. 1) then     

C==== ignore the strain effect on the so
            temp = ev(band4)-ev(band3)+ev(band4)-ev(band2)-
     &           calcvalue(idso0)

C==== include the strain effect on the so
c     temp1=2.0d0*ev(band4)-ev(band3)-ev(band2)
c     temp2=(ev(band3)-ev(band2))**2

c     temp = (temp1-dsqrt((3.0d0*temp2-temp1*temp1)/2.0d0))/2.0d0

            calcvalue(idpb100)=-temp/3.0d0/(e_zz-e_xx)

         endif                  !if(SO)

      endif                     ! b(100)



C==== b[111]

      alat=alat0
      do ips=1,nps
         psf(ips)=0.0d0
      enddo

      a1(1)=(a01(1)+e_zz*(a01(2)+a01(3)))*alat
      a1(2)=(a01(2)+e_zz*(a01(3)+a01(1)))*alat
      a1(3)=(a01(3)+e_zz*(a01(1)+a01(2)))*alat
      a2(1)=(a02(1)+e_zz*(a02(2)+a02(3)))*alat
      a2(2)=(a02(2)+e_zz*(a02(3)+a02(1)))*alat
      a2(3)=(a02(3)+e_zz*(a02(1)+a02(2)))*alat
      a3(1)=(a03(1)+e_zz*(a03(2)+a03(3)))*alat
      a3(2)=(a03(2)+e_zz*(a03(3)+a03(1)))*alat
      a3(3)=(a03(3)+e_zz*(a03(1)+a03(2)))*alat

c     As well as rescaling the lattice vectors, have to move atomic positions
c     The In stays at the origin and the P moves according to the strain.

      tau(1,1) = tau0(1,1)
      tau(2,1) = tau0(2,1)
      tau(3,1) = tau0(3,1)
      tau(1,2) = tau0(1,2)+e_zz*(tau0(2,2)+tau0(3,2))
      tau(2,2) = tau0(2,2)+e_zz*(tau0(3,2)+tau0(1,2))
      tau(3,2) = tau0(3,2)+e_zz*(tau0(1,2)+tau0(2,2))

      call recvec (a1,a2,a3,b1,b2,b3,omega)
      call gridg ()
      call vcell (vr,vc,-1)

      kp(1)=0.0d0
      kp(2)=0.0d0
      kp(3)=0.0d0

      call pwk(kp,vr,vc,-1,ev)

C     = deformation potential under [111] strain PRB 49, 14337 (1994) Eq.25
      if(idpd111 .gt. 0) then

         if(iSOps .eq. 0)      
     &        calcvalue(idpd111)=
     &        -(ev(band4)-ev(band2))/3.0d0/dsqrt(3.0d0)/e_zz

         if(iSOps .eq. 1) then     

C==== ignore the strain effect on the so
c     temp = ev(band4)-ev(band3)+ev(band4)-ev(band2)-
c     &          calcvalue(idso0)

C==== include the strain effect on the so
            temp1=2.0d0*ev(band4)-ev(band3)-ev(band2)
            temp2=(ev(band3)-ev(band2))**2

            temp = (temp1-dsqrt((3.0d0*temp2-temp1*temp1)/2.0d0))/2.0d0

            calcvalue(idpd111)=-temp/3.0d0/dsqrt(3.0d0)/e_zz

         endif                  !if(SO)

      endif                     ! d(111)



c     Fit the VBO if there are any to be fitted

      write(6,*) 'a(||)      VBM  '
      
      do i=1,mvbo
         if (ivbo(i).gt.0) then 
            alat = alat0
C===  DV/V
C     hydrofac=(1.0d0+eps_paral(i))**2*(1.0d0+eps_perp(i))-1.d0
C===  Tr(e)
            hydrofac=2.0d0*eps_paral(i)+eps_perp(i)
C============
            do ips=1,nps
               psf(ips)=hydrofac*psf0(ips)
            enddo

c     alat=vboalat(i)
            a1(1)=a01(1)*alat*(1.0d0+eps_paral(i))
            a1(2)=a01(2)*alat*(1.0d0+eps_paral(i))
            a1(3)=a01(3)*alat*(1.0d0+eps_perp(i))
            a2(1)=a02(1)*alat*(1.0d0+eps_paral(i))
            a2(2)=a02(2)*alat*(1.0d0+eps_paral(i))
            a2(3)=a02(3)*alat*(1.0d0+eps_perp(i))
            a3(1)=a03(1)*alat*(1.0d0+eps_paral(i))
            a3(2)=a03(2)*alat*(1.0d0+eps_paral(i))
            a3(3)=a03(3)*alat*(1.0d0+eps_perp(i))

            
c     Move the atomic positions
            
            tau(1,1) = tau0(1,1)
            tau(2,1) = tau0(2,1)
            tau(3,1) = tau0(3,1)
            tau(1,2) = tau0(1,2)*(1.0d0+eps_paral(i))
            tau(2,2) = tau0(2,2)*(1.0d0+eps_paral(i))
            tau(3,2) = tau0(3,2)*(1.0d0+eps_perp(i))

            call recvec (a1,a2,a3,b1,b2,b3,omega)
            call gridg ()
            call vcell (vr,vc,-1)
            
            kp(1)=0.0d0
            kp(2)=0.0d0
            kp(3)=0.0d0
            call pwk(kp,vr,vc,-1,ev)
            calcvalue(ivbo(i))=ev(band4)-vbovbm(i)/Ryd

            write(6,119) (1.0d0+eps_paral(i))*alat0,
     &           ev(band2)*Ryd,ev(band3)*Ryd,
     &           ev(band4)*Ryd,ev(band5)*Ryd

 119        format(5f10.5)
         end if
      end do

      do i=1,nproperty
         calcvalue(i) = calcvalue(i)*Ryd
      enddo


c     Now set up the residuals whose sum of squares is to be minimised
c     eveyrhing is in eV, now

      write(6,*)'       Property   Expt.  Value   ',
     &     'Calulated Value   Error'
      write(6,*)'       --------   -------------   ',
     &     ' --------------   -------------'
      if(iSOps .eq. 0) then
         do i=1,nproperty
            write(6,100)the_data(i)%property,the_data(i)%nsovalue,
     &           calcvalue(i),
     &           calcvalue(i)-the_data(i)%nsovalue
         end do
      else
         do i=1,nproperty
            write(6,100)the_data(i)%property,the_data(i)%sovalue,
     &           calcvalue(i),
     &           calcvalue(i)-the_data(i)%sovalue
         enddo
      endif
 100  format(a20,f10.5,f10.5,f10.5)


C======band structure L-Gamma

      alat=alat0

      do ips=1,nps
         psf(ips)=0.0d0
      enddo

      a1(1)=a01(1)*alat
      a1(2)=a01(2)*alat
      a1(3)=a01(3)*alat
      a2(1)=a02(1)*alat
      a2(2)=a02(2)*alat
      a2(3)=a02(3)*alat
      a3(1)=a03(1)*alat
      a3(2)=a03(2)*alat
      a3(3)=a03(3)*alat

      tau(1,1) = tau0(1,1)
      tau(2,1) = tau0(2,1)
      tau(3,1) = tau0(3,1)
      tau(1,2) = tau0(1,2)
      tau(2,2) = tau0(2,2)
      tau(3,2) = tau0(3,2)

      call recvec (a1,a2,a3,b1,b2,b3,omega)
      call gridg ()
      call vcell (vr,vc,-1)

 111  format(i5,8f10.5)
 112  format(i5,16f10.5)

      open(unit=8,file='bsLGXW.out',form='formatted',
     &     status='unknown')

c     L to Gamma
      
      do i=0,nkpoints
         kp(1)=(0.5-0.5*invkpoints*dble(i))
         kp(2)=(0.5-0.5*invkpoints*dble(i))
         kp(3)=(0.5-0.5*invkpoints*dble(i))

         kpc(1)=kp(1)*b1(1)+kp(2)*b2(1)+kp(3)*b3(1)
         kpc(2)=kp(1)*b1(2)+kp(2)*b2(2)+kp(3)*b3(2)
         kpc(3)=kp(1)*b1(3)+kp(2)*b2(3)+kp(3)*b3(3)

         call pwset (kpc)

         if(iSOps .eq. 0) then
            call solve_nso(kpc,vr,vc,-1,ev)
            write(8,111) i, (ev(j)*Ryd,j=1,8)
         endif
         if(iSOps .eq. 1) then
            call solve_so(kpc,vr,vc,-1,ev)
            write(8,112) i, (ev(j)*Ryd,j=1,16)
         endif

      enddo
      
c     Gamma to X
      
      do i=1,nkpoints
         kp(1)=0.5*invkpoints*dble(i)
         kp(2)=0.5*invkpoints*dble(i)
         kp(3)=0.0d0

         kpc(1)=kp(1)*b1(1)+kp(2)*b2(1)+kp(3)*b3(1)
         kpc(2)=kp(1)*b1(2)+kp(2)*b2(2)+kp(3)*b3(2)
         kpc(3)=kp(1)*b1(3)+kp(2)*b2(3)+kp(3)*b3(3)

         call pwset (kpc)

         if(iSOps .eq. 0) then
            call solve_nso(kpc,vr,vc,-1,ev)
            write(8,111) i+nkpoints, (ev(j)*Ryd,j=1,8)
         endif
         if(iSOps .eq. 1) then
            call solve_so(kpc,vr,vc,-1,ev)
            write(8,112) i+nkpoints, (ev(j)*Ryd,j=1,16)
         endif

      enddo

c     X to W
      
      do i=1,nkpoints
         kp(1)=0.5 + 0.25*invkpoints*dble(i)
         kp(2)=0.5
         kp(3)=0.25*invkpoints*dble(i)

         kpc(1)=kp(1)*b1(1)+kp(2)*b2(1)+kp(3)*b3(1)
         kpc(2)=kp(1)*b1(2)+kp(2)*b2(2)+kp(3)*b3(2)
         kpc(3)=kp(1)*b1(3)+kp(2)*b2(3)+kp(3)*b3(3)

         call pwset (kpc)

         if(iSOps .eq. 0) then
            call solve_nso(kpc,vr,vc,-1,ev)
            write(8,111) i+2*nkpoints, (ev(j)*Ryd,j=1,8)
         endif
         if(iSOps .eq. 1) then
            call solve_so(kpc,vr,vc,-1,ev)
            write(8,112) i+2*nkpoints, (ev(j)*Ryd,j=1,16)
         endif

      enddo
      
      close(8)
c     band strucure changes under hydrostatic pressure
      write(6,*) ' DV/V     band2    band3   band4   band5'

      do i=-10,10

         dvol = 0.001*dble(i)
C     hydrofac = dvol*(3.0d0+dvol*(3.0d0+dvol))
         hydrofac = 3.0d0*dvol
         do ips=1,nps
            psf(ips)=hydrofac*psf0(ips)
         enddo

c     hydrofac = 1.0d0/hydrofac

         alat=alat0*(1.0d0+dvol)
         a1(1)=a01(1)*alat
         a1(2)=a01(2)*alat
         a1(3)=a01(3)*alat
         a2(1)=a02(1)*alat
         a2(2)=a02(2)*alat
         a2(3)=a02(3)*alat
         a3(1)=a03(1)*alat
         a3(2)=a03(2)*alat
         a3(3)=a03(3)*alat

         call recvec (a1,a2,a3,b1,b2,b3,omega)
         call gridg ()
         call vcell (vr,vc,-1)

         kp(1)=0.0d0
         kp(2)=0.0d0
         kp(3)=0.0d0
         call pwk (kp,vr,vc,-1,ev)

         write(6,200) hydrofac,ev(band2)*Ryd,
     &        ev(band3)*Ryd,ev(band4)*Ryd,ev(band5)*Ryd

      enddo
c     band structure changes under biaxial strain
      write(6,*) '  Tr(e)      band2     band3     band4    band5'

      do i=-10,10
         e_zz = de_zz*dble(i)
         e_xx = 1.0d0/dsqrt(1.0d0+e_zz)-1.0d0

         hydrofac = 2.0d0*e_xx+e_zz
         do ips=1,nps
            psf(ips)=hydrofac*psfvol(ips)
         enddo

         alat=alat0
         a1(1)=a01(1)*(1+e_xx)*alat
         a1(2)=a01(2)*(1+e_xx)*alat
         a1(3)=a01(3)*(1+e_zz)*alat
         a2(1)=a02(1)*(1+e_xx)*alat
         a2(2)=a02(2)*(1+e_xx)*alat
         a2(3)=a02(3)*(1+e_zz)*alat
         a3(1)=a03(1)*(1+e_xx)*alat
         a3(2)=a03(2)*(1+e_xx)*alat
         a3(3)=a03(3)*(1+e_zz)*alat

c     As well as rescaling the lattice vectors, have to move atomic positions
c     The In stays at the origin and the P moves according to the strain.

         tau(1,1) = tau0(1,1)
         tau(2,1) = tau0(2,1)
         tau(3,1) = tau0(3,1)
         tau(1,2) = tau0(1,2)*(1+e_xx)
         tau(2,2) = tau0(2,2)*(1+e_xx)
         tau(3,2) = tau0(3,2)*(1+e_zz)

         call recvec (a1,a2,a3,b1,b2,b3,omega)
         call gridg ()
         call vcell (vr,vc,-1)

         kp(1)=0.0d0
         kp(2)=0.0d0
         kp(3)=0.0d0

         call pwk(kp,vr,vc,-1,ev)
         write(6,200) e_xx,ev(band2)*Ryd,
     &        ev(band3)*Ryd,ev(band4)*Ryd,ev(band5)*Ryd

      enddo
 200  format(5(f10.5))
      return
      end

