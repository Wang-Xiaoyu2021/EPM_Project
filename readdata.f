      subroutine readdata(nparam)

      implicit none
      include "epm.inc"

      character*30,allocatable,dimension(:) :: pseudofile

      integer i,ips,jps,nparam,iat,j,tempmove,match,ivalue,iv,ii,
     $     npseudofile,iloop,itemp,isum,j1,j2,npw_t
      integer temp_mvpsf,temp_mvbeta,ierr,temppstype,ifac,pstype0
      double precision tempps,vboweight,vbo,temp
      real*8 x1,x2,x3,ecut_all,smth_all,Ecut_t,Smth_t,scalkin_t
      character*2 tempname
cc      character*70 linestring,file_dum
      character*70 file_dum
      character*120 linestring

      real*8 proj(mnpw*2,mnpw*2)
      complex*16 zz_st(mnpw*2,mnpw*2,mns)
      integer iproj(mns),ipcbm(mns),iflag(mns)
      common  /comproj/zz_st,iproj,ipcbm
      common /EcutSmth/Ecut_t,Smth_t,scalkin_t


      type (fitdata) tempdata
      ivbo(:,:)=0

c -----------------------------------------------
c Read in the fitting data
c -----------------------------------------------
      open (unit=4,file='fit.d',status='old')
      read(4,*)nstruct
      do i=1,nstruct
         read(4,*)structfile(i)
         read(4,*)ng1(i),ng2(i),ng3(i)
         read(4,*)inputfile(i)
	 read(4,*)iproj(i),ipcbm(i),file_dum
         read(4,*)xpt(1,i),xpt(2,i),xpt(3,i)
         read(4,*)
      if(iproj(i).eq.1) then
      open(10,file=file_dum,form="unformatted")
      rewind(10)
      read(10) iflag(i),npw_t
      read(10) ((zz_st(j1,j2,i),j1=1,npw_t),j2=1,npw_t)
      close(10)
      endif
      end do
      read(4,*)npseudofile
      allocate(pseudofile(npseudofile))
      do i=1,npseudofile
         read(4,*)pseudofile(i)
      end do
      read(4,*)
      read(4,*) iSOps

      do i=1,nstruct
      if(iproj(i).eq.1.and.iflag(i).ne.iSOps+1) then
      write(6,*) "input iSO in ug.input does not march iSO in fit, stop"
      stop
      endif
      enddo

      if (iSOps==1) write(6,*)'Spin-Orbit interaction is ON'
      if (iSOps==0) write(6,*)'Spin-Orbit interaction is OFF'
      read(4,*)
      read(4,*)ifit_sig
      read(4,*)sigma0
      write(6,*)'Reweight of KE term ',sigma0
      if (ifit_sig==0) then 
         write(6,*)'Fixed Reweight'
      else
         write(6,*)'Fitting Reweight'
      end if

      read(4,*) ecut_all,smth_all
      do i=1,nstruct
         ecut(i)=ecut_all
         smth(i)=smth_all
      enddo
      close(4)

      write(6,*) "debug: close fit.d file "
      nparam=0
      if (ifit_sig/=0) nparam=1
c -----------------------------------------------
c set Initial bounds for pseudopotential parameters
c -----------------------------------------------
      psSO0L(:)=0.d0
      psSO0U(:)=2.d0
      psvol0L(:)=0.1
      psvol0U(:)=200.0
      psf0L(:)=0.d0
      psf0U(:)=1.d0
      psbeta0L(:)=0.d0
      psbeta0U(:)=10.d0
      psa0L(:,:)=0.d0
      psa0U(:,:)=300.d0
      psb0L(:,:)=0.d0
      psb0U(:,:)=dsqrt(20.d0)
      psc0L(:,:)=0.d0
      psc0U(:,:)=10.d0
c -----------------------------------------------
c     Read in the pseudopotential parameters
c -----------------------------------------------
      itotps=0
      do i=1,npseudofile
         open (unit=4,file=pseudofile(i),status='old')
         read (4,*) nps, pstype
         if(itotps.eq.0) pstype0=pstype
         if(pstype.ne.pstype0) then
         write(6,*) "EPM format is inconsistent, stop"
         stop
         endif

         read (4,*) Ecut_t, Smth_t, scalkin_t

         if(dabs(Ecut_t-ecut_all)+dabs(Smth_t-smth_all)+
     &       dabs(scalkin_t-Sigma0).gt.0.01) then
         write(6,*) "Ecut,Smth and Sigma(scalkin) in fit.d and\
     &     pseudo.xx.fit do not agree, stop"
           stop
           endif
c***** format 5 or 51 *********************************
c* in this special format, 5, some parameters are fixed.
         if(pstype.eq.5.or.pstype.eq.51) then
         do 10 ips=1,nps
            itotps=itotps+1
            ngauss(itotps)=2

            psvol(itotps)=1.d0    ! all these param are fixed in format type 5 
            psbeta(itotps)=0.d0
            mvpsbeta(itotps)=0
            psb(2,itotps)=2.d0
            mvpsb(2,itotps)=0
            psc(2,itotps)=0.d0
            mvpsc(2,itotps)=0

c* set bounds of parameters 
            psSO0L(itotps)=0.d0
            psSO0U(itotps)=2.d0
            psf0L(itotps)=0.d0
            psf0U(itotps)=1.d0
            psa0L(1,itotps)=1.d0                  ! 1.0 < a1 < 500.0
            psa0U(1,itotps)=500.d0                
            psb0L(1,itotps)=1.d0                  ! 1.0 < b1 < 10.0
            psb0U(1,itotps)=10.d0                 
            psc0L(1,itotps)=0.d0                  ! 0.0 < c1 < Ecut
            psc0U(1,itotps)=Ecut_t               
            psa0L(2,itotps)=0.d0                  ! 0.1 < a2 < 10.0
            psa0U(2,itotps)=10.d0               
            psb0L(2,itotps)=100.d0                ! not used 
            psb0U(2,itotps)=100.d0                ! not used 
            psc0L(2,itotps)=100.d0                ! not used
            psc0U(2,itotps)=100.d0                ! not used

            read (4,*)
            read (4,*) psname(itotps),psnum(itotps)
            read (4,*) psSO(itotps),moveSO(itotps)
            read (4,*) psfvol(itotps),mvpsf(itotps)
            read (4,*) psa(1,itotps),mvpsa(1,itotps)
            read (4,*) psb(1,itotps),mvpsb(1,itotps)
            read (4,*) psc(1,itotps),mvpsc(1,itotps)
            read (4,*) psa(2,itotps),mvpsa(2,itotps)

          if(iSOps.eq.0.and.moveSO(itotps).eq.1) then
          write(6,*) "iSOps=0 in fit.d, but moveSO.eq.1 in pseudo, stop"
          stop
          endif

         nparam=nparam+moveSO(itotps)+mvpsf(itotps)+
     &    mvpsa(1,itotps)+mvpsb(1,itotps)+mvpsc(1,itotps)+
     &    mvpsa(2,itotps)
 10      continue
c***** format 1 *********************************
      else if(pstype.eq.1) then
         do 20 ips=1,nps
         itotps=itotps+1
         read(4,*)
         read(4,*) psname(itotps),psnum(itotps)
         read(4,*) psvol(itotps),ngauss(itotps)     
         read(4,*) psSO(itotps),moveSO(itotps)
         read(4,*) psfvol(itotps), mvpsf(itotps)
         read(4,*) psbeta(itotps), mvpsbeta(itotps)
         nparam=nparam+moveSO(itotps)+mvpsf(itotps)+mvpsbeta(itotps)

         do 30 j=1,ngauss(itotps)
         read(4,*) psa(j,itotps),psb(j,itotps),psc(j,itotps),
     &     mvpsa(j,itotps), mvpsb(j,itotps), mvpsc(j,itotps)
         nparam=nparam+mvpsa(j,itotps)+mvpsb(j,itotps)+mvpsc(j,itotps) 
30       continue

c* set bounds of parameters for format 1
         psSO0L(itotps)=0.d0          ! 0<= SO < 1.d0
         psSO0U(itotps)=2.d0
         psvol0L(itotps)=0.1d0        ! 0.1 < vol < 200.0
         psvol0U(itotps)=200.d0
         psf0L(itotps)=0.d0           ! 0.0 < fstr < 1.0
         psf0U(itotps)=1.d0
         psbeta0L(itotps)=10.d0       ! this parameter was not used
         psbeta0U(itotps)=10.d0       ! this parameter was not used
         psa0L(:,itotps)=0.d-3        ! 0.0 < ai < 300.d0
         psa0U(:,itotps)=300.d0
         psb0L(:,itotps)=0.d0         ! 0.0 < bi < dsqrt(Ecut)
         psb0U(:,itotps)=dsqrt(20.d0)
         psc0L(:,itotps)=0.1d0        ! 0.1 < ci < 10.0
         psc0U(:,itotps)=10.d0
20       continue
c***** format 52 *********************************
      else if(pstype.eq.52) then
         do 40 ips=1,nps
            itotps=itotps+1
            ngauss(itotps)=2

            psvol(itotps)=1.d0    
            psbeta(itotps)=0.d0
            mvpsbeta(itotps)=0
            psb(2,itotps)=0.d0
            mvpsb(2,itotps)=0
            psc(2,itotps)=0.d0
            mvpsc(2,itotps)=0

c* set bounds of parameter for format 52
            psSO0L(itotps)=0.d0
            psSO0U(itotps)=2.d0
            psf0L(itotps)=0.d0
            psf0U(itotps)=1.d0
            psa0L(1,itotps)=1.0d0                 ! 1.0 < a1 < 1000.0
            psa0U(1,itotps)=1000.d0              
            psb0L(1,itotps)=1.d0                  ! 1.0 < b1 < 10.0
            psb0U(1,itotps)=10.d0              
            psc0L(1,itotps)=0.d0                  ! 0.0 < c1 < Ecut
            psc0U(1,itotps)=Ecut_t       
            psa0L(2,itotps)=-1.d0                 ! -1.0 < a2 < Ecut/2
            psa0U(2,itotps)=Ecut_t/2.d0
            psb0L(2,itotps)=0.d0                  ! 0.0 < b2 < 10.0
            psb0U(2,itotps)=10.d0       
            psc0L(2,itotps)=0.d0                  ! 0.0 < c2 < 1.0
            psc0U(2,itotps)=1.d0                

            read (4,*)
            read (4,*) psname(itotps),psnum(itotps)
            read (4,*) psSO(itotps),moveSO(itotps)
            read (4,*) psfvol(itotps),mvpsf(itotps)
            read (4,*) psa(1,itotps),mvpsa(1,itotps)
            read (4,*) psb(1,itotps),mvpsb(1,itotps)
            read (4,*) psc(1,itotps),mvpsc(1,itotps)
            read (4,*) psa(2,itotps),mvpsa(2,itotps)
            read (4,*) psb(2,itotps),mvpsb(2,itotps)
            read (4,*) psc(2,itotps),mvpsc(2,itotps)

          if(iSOps.eq.0.and.moveSO(itotps).eq.1) then
          write(6,*) "iSOps=0 in fit.d, but moveSO.eq.1 in pseudo, stop"
          stop
          endif

         nparam=nparam+moveSO(itotps)+mvpsf(itotps)+
     &    mvpsa(1,itotps)+mvpsb(1,itotps)+mvpsc(1,itotps)+
     &    mvpsa(2,itotps)+mvpsb(2,itotps)+mvpsc(2,itotps)
40      continue
c***** format 53 *********************************
      else if(pstype.eq.53) then
         do 50 ips=1,nps
            itotps=itotps+1
            ngauss(itotps)=2

            psvol(itotps)=1.d0
            psbeta(itotps)=0.d0
            mvpsbeta(itotps)=0
            psb(2,itotps)=2.d0
            mvpsb(2,itotps)=0
            psc(2,itotps)=0.d0
            mvpsc(2,itotps)=0.d0

            read (4,*)
            read (4,*) psname(itotps),psnum(itotps)
            read (4,*) psSO(itotps),moveSO(itotps)
            read (4,*) psfvol(itotps),mvpsf(itotps)
            read (4,*) psbeta(itotps),mvpsbeta(itotps)
            read (4,*) psa(1,itotps),mvpsa(1,itotps)
            read (4,*) psb(1,itotps),mvpsb(1,itotps)
            read (4,*) psc(1,itotps),mvpsc(1,itotps)
            read (4,*) psa(2,itotps),mvpsa(2,itotps)

          if(iSOps.eq.0.and.moveSO(itotps).eq.1) then
          write(6,*) "iSOps=0 in fit.d, but moveSO.eq.1 in pseudo, stop"
          stop
          endif

         nparam=nparam+
     &    moveSO(itotps)+mvpsf(itotps)+mvpsbeta(itotps)+
     &    mvpsa(1,itotps)+mvpsb(1,itotps)+mvpsc(1,itotps)+
     &    mvpsa(2,itotps)+mvpsb(2,itotps)+mvpsc(2,itotps)
50      continue
c***** format 54 *********************************
      else if(pstype.eq.54) then
         do 60 ips=1,nps
            itotps=itotps+1
            ngauss(itotps)=2

c* set bounds of parameter for format 54
            psSO0L(itotps)=0.d0
            psSO0U(itotps)=2.d0
            psf0L(itotps)=0.d0
            psf0U(itotps)=1.d0
            psbeta0L(itotps)=0.d0                ! 0.0 < qmax < dsqrt(Ecut)
            psbeta0U(itotps)=dsqrt(Ecut_t)     
            psa0L(1,itotps)=0.1d0                ! 0.1 < V0 < 300.0
            psa0U(1,itotps)=300.d0              
            psb0L(1,itotps)=0.1d0                ! 0.1 < pl < 1.0
            psb0U(1,itotps)=1.d0      
            psc0L(1,itotps)=1.d0                 ! 1.0 < e1 < 10.0
            psc0U(1,itotps)=10.d0     
            psa0L(2,itotps)=0.d0                 ! 0.0 < Vmax < 50.0
            psa0U(2,itotps)=50.d0   
            psb0L(2,itotps)=0.1d0                ! 0.1 < pr < 1.0
            psb0U(2,itotps)=10.d0        
            psc0L(2,itotps)=1.d0                 ! 1.0 < er < 10.0
            psc0U(2,itotps)=10.d0        

            read (4,*)
            read (4,*) psname(itotps),psnum(itotps)
            read (4,*) psSO(itotps),moveSO(itotps)
            read (4,*) psfvol(itotps),mvpsf(itotps)
            read (4,*) psbeta(itotps),mvpsbeta(itotps)    ! qmax
            read (4,*) psa(1,itotps),mvpsa(1,itotps)      ! V0
            read (4,*) psb(1,itotps),mvpsb(1,itotps)      ! pl
            read (4,*) psc(1,itotps),mvpsc(1,itotps)      ! el
            read (4,*) psa(2,itotps),mvpsa(2,itotps)      ! Vmax
            read (4,*) psb(2,itotps),mvpsb(2,itotps)      ! pr
            read (4,*) psc(2,itotps),mvpsc(2,itotps)      ! er

          if(iSOps.eq.0.and.moveSO(itotps).eq.1) then
          write(6,*) "iSOps=0 in fit.d, but moveSO.eq.1 in pseudo, stop"
          stop
          endif

         nparam=nparam+
     &    moveSO(itotps)+mvpsf(itotps)+mvpsbeta(itotps)+
     &    mvpsa(1,itotps)+mvpsb(1,itotps)+mvpsc(1,itotps)+
     &    mvpsa(2,itotps)+mvpsb(2,itotps)+mvpsc(2,itotps)
60      continue
      else
         write(6,*) "undefined EPM format, stop", pstype
         stop
      endif
      close (unit=4)
      end do    ! all the pseudopotetial files
      write(6,*)'Read in a total of ',itotps,' pseudopotentials'
      write(6,*)'There are a total of',nparam,'parameters to be fitted'
      write(6,*)

      
      do 110 ips=1,itotps
         psvol0(ips)=psvol(ips)
         psf0(ips)=psfvol(ips)
         psbeta0(ips)=psbeta(ips)
         psSO0(ips)=psSO(ips)
         do 120 j=1,ngauss(ips)
            psa0(j,ips)=psa(j,ips)
            psb0(j,ips)=psb(j,ips)
            psc0(j,ips)=psc(j,ips)
 120     continue
 110  continue

c -----------------------------------------------
c Read in the experimental data to be fitted for each structure
c -----------------------------------------------
         nsprop = 0
         nproperty=0
         do i=1,nstruct
         iv=0    ! for the counting of evbm_strain, it can have multiple ones
            write(6,*)'Reading expt. data from ',inputfile(i)
            open (unit=7,file=inputfile(i),form='formatted',
     $           status='old')
            write(6,*)'    Input Quantity      Non So.      S.O.',
     &           '    weight'
            write(6,*)'    --------------       -----     ------',
     &        '    ------'

            read(7,*)
            read(7,*)
         data:  do ii=1,100
         call getline(7,linestring,ierr)
         if (ierr==1) goto 2000
         if (ierr==2) goto 2000
         read(linestring,*) tempdata
         nproperty = nproperty + 1
         nsprop(i)=nsprop(i) + 1
         select case (tempdata%property)
            case('eg1v') 
               ieg1v(i) = nsprop(i)
            case('evbm') 
               ievbm(i) = nsprop(i)
            case('evbm_1') 
               ievbm_1(i) = nsprop(i)
            case('eg1c') 
               ieg1c(i) = nsprop(i)
            case('eg15c') 
               ieg15c(i) = nsprop(i)
            case('ex1v') 
               iex1v(i) = nsprop(i)
            case('ex3v')
               iex3v(i) = nsprop(i)
            case('ex5v') 
               iex5v(i) = nsprop(i)
            case('ex1c') 
               iex1c(i) = nsprop(i)
            case('ex3c') 
               iex3c(i) = nsprop(i)
            case('el1v')
               iel1v(i) = nsprop(i)
            case('el2v') 
               iel2v(i) = nsprop(i)
            case('el3v') 
               iel3v(i) = nsprop(i)
            case('el1c') 
               iel1c(i) = nsprop(i)
            case('el3c') 
               iel3c(i) = nsprop(i)
            case('mg1c') 
               img1c(i) = nsprop(i)
            case('mghh100') 
               imghh100(i) = nsprop(i)
            case('mghh111')
               imghh111(i) = nsprop(i)
            case('mghh110')
               imghh110(i) = nsprop(i)
            case('mglh100') 
               imglh100(i) = nsprop(i)
            case('mglh111') 
               imglh111(i) = nsprop(i)
            case('mglh110') 
               imglh110(i) = nsprop(i)
            case('dpa') 
               idpa(i) = nsprop(i)
            case('dpag1c') 
               idpag1c(i) = nsprop(i)
            case('dpag15v') 
               idpag15v(i) = nsprop(i)
            case('dpaGX') 
               idpaGX(i) = nsprop(i)
            case('dpb100') 
               idpb100(i) = nsprop(i)
            case('edso0') 
               idso0(i) = nsprop(i)
            case('edso1')
               idso1(i) = nsprop(i) 
            case('edso2')
               idso2(i) = nsprop(i)
            case('egp3pc')
               iegp3pc(i) = nsprop(i) 
            case('egm3pc')
               iegm3pc(i) = nsprop(i) 
            case('evbm_strain')  ! special ones, replace the old VBO
               iv=iv+1
               ivbo(iv,i) = nsprop(i)
               vbovbm(iv,i) = 0.d0    ! no longer used in this new version
            read(linestring,*) tempdata, eps_paral(iv,i),
     &          eps_perp(iv,i)      ! need strain inform for this
            case('ml1cl')
               iml1cl(i) = nsprop(i)
            case('ml1ct')
               iml1ct(i) = nsprop(i)
            case('ml3vl')
               iml3vl(i) = nsprop(i)
            case('ml3vt')
               iml3vt(i) = nsprop(i)
            case('mx1cl')
               imx1cl(i) = nsprop(i)
            case('mx1ct')
               imx1ct(i) = nsprop(i)
            case('dpaGL')
               idpaGL(i) = nsprop(i)
            case('dcamel')
               idcamel(i)  = nsprop(i)
            case('edso0p') 
               idso0p(i) = nsprop(i)
            case('edcf0')
               idcf0(i) = nsprop(i)
            case('mg1c100')
               img1c100(i) = nsprop(i)
            case('mg1c001')
               img1c001(i) = nsprop(i)
            case('mghh001')
               imghh001(i) = nsprop(i)
            case('mglh001')
               imglh001(i) = nsprop(i)
            case('lgap')
               ilgap(i) = nsprop(i)
            case('l57')
               il57(i) = nsprop(i)
            case('l46')
               il46(i) = nsprop(i)
            case('nvband')    ! for structures with num of val. bands rather than 4
               invband(i) = nsprop(i)
            case('sscb1')     ! following targets are for fit to GW
               isscb1(i) = nsprop(i)
            case('ssvb1')
               issvb1(i) = nsprop(i)
            case('engvb1')
               iengvb1(i) = nsprop(i)
            case('engvb2')
               iengvb2(i) = nsprop(i)
            case('engvb3')
               iengvb3(i) = nsprop(i)
            case('engvb4')
               iengvb4(i) = nsprop(i)
            case('engcb1')
               iengcb1(i) = nsprop(i)
            case('engcb2')
               iengcb2(i) = nsprop(i)
            case('engcb3')
               iengcb3(i) = nsprop(i)
            case('engcb4')
               iengcb4(i) = nsprop(i)
            case default
               write(6,*)'Could not find a match for input ',
     &        tempdata%property
         stop
      end select
      the_data(nsprop(i),i)=tempdata
      write(6,100)tempdata
 100  format(a20,3(2x,f10.5))
      end do data

 2000  continue
      close(7)
cccccccccccccccccccccccccccccccc

      isum=0
      write(6,*) nsprop(i),
     $     ' Properties to be fitted for ',structfile(i)
      write(6,*)
      
      end do                    !Loop over structures
      write(6,*)nproperty,' total expt. properties to be fitted'
      write(6,*)
c ----------------------------------------------------
c Read in the spin orbit data file
c ----------------------------------------------------
      if (iSOps.eq.1) then
c         open(10,file='VG1G2.SO',form='unformatted')
         open(10,file='VG1G2.SO.form')
         rewind(10)
         read(10,*)gA
         read(10,*)vso
         close(10)
      end if
c ----------------------------------------------------
c Read in the positions of the atoms for the structure files
c ----------------------------------------------------
      do i=1,nstruct
         write(6,*)'Reading in structure ',structfile(i)
         open (unit=4,file=structfile(i),status='old')
         read (4,*) natoms(i)
	 read (4,*) alat0(i)     
         read (4,*) a0(1,:,i)
         read (4,*) a0(2,:,i)
         read (4,*) a0(3,:,i)
         do 101 iat=1,natoms(i)
ccccc use the unified atom.config convention, ie. the input coord
ccccc is the supercell coord, not the x,y,z
cccccc tau0 is the x,y,z coord in terms of alat0
         read(4,*) atomnum(iat,i),x1,x2,x3,
     $           deltav(iat,i),atweight0(iat,i)

        tau0(1,iat,i)=a0(1,1,i)*x1+a0(2,1,i)*x2+a0(3,1,i)*x3
        tau0(2,iat,i)=a0(1,2,i)*x1+a0(2,2,i)*x2+a0(3,2,i)*x3
        tau0(3,iat,i)=a0(1,3,i)*x1+a0(2,3,i)*x2+a0(3,3,i)*x3

            match=0
            do jps=1,itotps
               if (atomnum(iat,i)==psnum(jps)) match=1
            enddo
            if (match==0) then
               write(6,*)'Cannot match atom number ',atomnum(iat,i)
               stop
            end if
 101     continue
         close(4)
      end do
c -----------------------------------------------------
c     Assign the naming indicies for the bands.
c     Do this by calculating the total no. of electrons in the supercell.
c -----------------------------------------------------
      do is=1,nstruct
         temp=0.0
         do i=1,natoms(is)
            temp=temp+atweight0(i,is)
         end do
         ifac=int(temp*1.0001)/2
         if(iSOps.eq.1) iband(4,is)=8*ifac
         if(iSOps.eq.0) iband(4,is)=4*ifac
         if(invband(is).gt.0) then
           if(iSOps.eq.1) iband(4,is)
     &       =2*the_data(invband(is),is)%sovalue
           if(iSOps.eq.0) iband(4,is)
     &       =the_data(invband(is),is)%nsovalue
         endif

         if(iSOps.eq.1) iband(5,is)=iband(4,is)+2
         if(iSOps.eq.1) iband(6,is)=iband(4,is)+4
         if(iSOps.eq.1) iband(7,is)=iband(4,is)+6
         if(iSOps.eq.1) iband(8,is)=iband(4,is)+8
c        if(iSOps.eq.1) iband(4,is)=iband(4,is)-1
         if(iSOps.eq.1) iband(3,is)=iband(4,is)-2
         if(iSOps.eq.1) iband(2,is)=iband(4,is)-4
         if(iSOps.eq.1) iband(1,is)=iband(4,is)-6

c         if(iSOps.eq.1) iband(5,is)=8*ifac+1
         if(iSOps.eq.0) iband(5,is)=iband(4,is)+1
         if(iSOps.eq.0) iband(6,is)=iband(4,is)+2
         if(iSOps.eq.0) iband(7,is)=iband(4,is)+3
         if(iSOps.eq.0) iband(8,is)=iband(4,is)+4
         if(iSOps.eq.0) iband(3,is)=iband(4,is)-1
         if(iSOps.eq.0) iband(2,is)=iband(4,is)-2
         if(iSOps.eq.0) iband(1,is)=iband(4,is)-3

ccccccccccc iband(4,is): VBM (even num in SO), iband(5,is): CBM (even num in SO)

      end do


      end
