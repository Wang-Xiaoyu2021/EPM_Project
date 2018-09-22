      subroutine funx (iflag,nprop,nparam,xc,fvecc,iw,liw,w,lw)
      
      implicit none

      include 'epm.inc' 

      integer k,i,j,ips,writeout,iflag,nparam,liw,lw,nprop
      integer iw(liw),icase,ii,ist,itest
      double precision e_xx,e_zz,dvol,hydrofac
      double precision xc(nparam),fvecc(mproperty),w(lw)
      double precision kp(1:3),ev(1:mnpw*2),vr(1:mngrid)
      double precision pi,vbmoffset,temp,temp1,temp2,dltk
      double precision cbm_shift,vbm_shift,psflocal(mnat)
      complex*16 psi,vc(1:mngrid)

      e_zz = 0.001d0

      pi=4.0d0*atan(1.0d0)

      k=0
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

C==== for the ideal zincblend structure, psf = 0        
         psf(1:itotps)=0.0d0
         
         alat=alat0(is)
         a(:,:)=a0(:,:,is)*alat
         tau(:,:,is)=tau0(:,:,is)

         call gridg ()
         call vcell (vr,vc,-1)

         kp(1)=0.0d0
         kp(2)=0.0d0
         kp(3)=0.0d0
         call pwk (kp,vr,vc,-1,ev)

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



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	 itest=iex1v(is)+iex3v(is)+iex5v(is)+iex1c(is)+iex3c(is)
	 itest=itest+iel1v(is)+iel2v(is)+iel3v(is)+iel1c(is)+
     &               idso1(is)
	 itest=itest+img1c(is)+imghh100(is)+imglh100(is)+
     &	             imghh111(is)
	 itest=itest+idpa(is)+idpag1c(is)+idpag15v(is)+
     &               idpaX1c(is)+iegp3pc(is)+iegm3pc(is)
	 do i=1,mvbo
         itest=itest+ivbo(i,is)
	 enddo

	 if(itest.le.0) goto 999
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


         kp(1)=0.5d0
         kp(2)=0.5d0
         kp(3)=0.0d0
         call pwk (kp,vr,vc,-1,ev)
         if(iex1v(is) .gt. 0) calcvalue(iex1v(is),is)=ev(iband(1,is))
     $        -vbmoffset
         if(iex3v(is) .gt. 0) calcvalue(iex3v(is),is)=ev(iband(2,is))
     $        -vbmoffset
         if(iex5v(is) .gt. 0) calcvalue(iex5v(is),is)=ev(iband(4,is))
     $        -vbmoffset
         if(iex1c(is) .gt. 0) calcvalue(iex1c(is),is)=ev(iband(5,is))
     $        -vbmoffset
         if(iex3c(is) .gt. 0) calcvalue(iex3c(is),is)=ev(iband(6,is))
     $        -vbmoffset

         kp(1)=0.5d0
         kp(2)=0.5d0
         kp(3)=0.5d0
         call pwk (kp,vr,vc,-1,ev)
         if(iel1v(is) .gt. 0) calcvalue(iel1v(is),is)=ev(iband(1,is))
     $        -vbmoffset
         if(iel2v(is) .gt. 0) calcvalue(iel2v(is),is)=ev(iband(2,is))
     $        -vbmoffset
         if(iel3v(is) .gt. 0) calcvalue(iel3v(is),is)=ev(iband(4,is))
     $        -vbmoffset
         if(iel1c(is) .gt. 0) calcvalue(iel1c(is),is)=ev(iband(5,is))
     $        -vbmoffset
         if(idso1(is) .gt. 0) calcvalue(idso1(is),is)=ev(iband(4,is))
     $        -ev(iband(3,is))
	 
         kp(1)=0.001d0
         kp(2)=0.001d0
         kp(3)=0.000d0
         call pwk (kp,vr,vc,-1,ev)

         dltk=0.001*4.0d0*pi/alat

         temp=ev(iband(5,is))-vbmoffset
         if(img1c(is) .gt. 0)  
     &        calcvalue(img1c(is),is)=dltk**2/Ryd/
     $        (temp-calcvalue(ieg1c(is),is))

         temp=ev(iband(4,is))-vbmoffset
         if(imghh100(is) .gt. 0)
     $        calcvalue(imghh100(is),is)=-dltk**2/Ryd/temp

         temp=ev(iband(3,is))-vbmoffset
         if(imglh100(is) .gt. 0)
     $        calcvalue(imglh100(is),is)=-dltk**2/Ryd/temp

C******heavy hole along [111]
         kp(1)=0.003d0
         kp(2)=0.003d0
         kp(3)=0.003d0
         call pwk (kp,vr,vc,-1,ev)

         dltk=0.003*2.0d0*pi*dsqrt(3.0d0)/alat

         if(iSOps.eq.0) then
         temp=ev(iband(4,is))-vbmoffset
	 endif

	 if(iSOps.eq.1) then
ccccc for (110) direction, the hh is splitted into two lines
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

c     Now change the lattice constant uniformaly and re-calculate the
c     bandstructure, to obtain the deformation potentials.
c     There is no need to move the atomic positions as vcell will scale them
c     according to alat as we want.
C=================================================================
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
         call pwk (kp,vr,vc,-1,ev)
         temp1=ev(iband(4,is))
         temp2=ev(iband(5,is))

         vbm_shift = ev(iband(4,is))-calcvalue(ievbm(is),is)
         cbm_shift = ev(iband(5,is))
     $        -(calcvalue(ievbm(is),is)+calcvalue(ieg1c(is),is))

         if(idpa(is) .gt. 0) 
     &        calcvalue(idpa(is),is)
     $        =(temp2-temp1-calcvalue(ieg1c(is),is))/hydrofac

         if(idpag1c(is) .gt. 0) calcvalue(idpag1c(is),is)
     $        =(temp2-vbmoffset-calcvalue(ieg1c(is),is))/hydrofac

         if(idpag15v(is) .gt. 0) 
     &        calcvalue(idpag15v(is),is)=(temp1-vbmoffset)/hydrofac

         kp(1)=0.5d0
         kp(2)=0.5d0
         kp(3)=0.0d0
         call pwk (kp,vr,vc,-1,ev)
         temp1=ev(iband(4,is))
         temp2=ev(iband(5,is))

         vbm_shift = ev(iband(4,is))-calcvalue(ievbm(is),is)
         cbm_shift = ev(iband(5,is))-(calcvalue(ievbm(is),is)
     $        +calcvalue(ieg1c(is),is))

         if(idpaX1c(is)>0) calcvalue(idpaX1c(is),is)
     $        =(temp2-vbmoffset-calcvalue(ieX1c(is),is))/hydrofac

c     Deform the lattice constant by +- 3 percent

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
               call pwk (kp,vr,vc,-1,ev)

               if(icase.eq.1) calcvalue(iegp3pc(is),is)
     $              =ev(iband(5,is))-ev(iband(4,is))
               if(icase.eq.2) calcvalue(iegm3pc(is),is)
     $              =ev(iband(5,is))-ev(iband(4,is))

            end if
         enddo
         
c     Now change the lattice constant non-uniformaly and re-calculate the
c     bandstructure, to obtain the epitaxial deformation potentials.
c     In the x-y plane the lattice constant is decreased by 1% and in the
c     z direction the lattice constant is increased by 1%.  The value is fitted
c     to the result of Eq.(8) in PRB 49, p14337.
c     The crystal field splits the 3-fold degeneracy with 2 bands 1/3 above the
c     original value and 1 band 2/3 below the original value.
C=====================
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

         call pwk(kp,vr,vc,-1,ev)

C     = deformation potential under [001] strain PRB 49, 14337 (1994) Eq.25
         if(idpb100(is) .gt. 0) then

            if(iSOps .eq. 0)      
     &           calcvalue(idpb100(is),is)=
     $           -(ev(iband(4,is))-ev(iband(2,is)))/3.0d0/(e_zz-e_xx)

            if(iSOps .eq. 1) then     
               temp = ev(iband(4,is))-ev(iband(3,is))+
     &              ev(iband(4,is))-ev(iband(2,is))
     $              -calcvalue(idso0(is),is)

               calcvalue(idpb100(is),is)=-temp/3.0d0/(e_zz-e_xx)
            endif               !if(SO)

         endif                  ! b(100)

c     Fit the VBO if there are any to be fitted
C==== now, again, add the hydrostatic term        
C==== 
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
               call pwk(kp,vr,vc,-1,ev)
               calcvalue(ivbo(i,is),is)=ev(iband(4,is))-vbovbm(i,is)/Ryd
            end if
            
         end do


999   continue

      end do                    ! Loop over structures
      
c     Now set up the residuals whose sum of squares is to be minimised
c     eveything is in eV, now
      
      calcvalue=calcvalue*Ryd
      
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
            enddo
         endif
      end do
      
      end

