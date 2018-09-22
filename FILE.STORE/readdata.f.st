      subroutine readdata(nparam)

      implicit none
      include "epm.inc"

      character*30,allocatable,dimension(:) :: inputfile,
     $     pseudofile,controlfile

      integer i,ips,jps,nparam,iat,j,tempmove,match,ivalue,iv,ii,
     $     npseudofile,iloop,itemp,isum,j1,j2,npw_t
      integer temp_mvpsf,temp_mvbeta,ierr,temppstype,ifac
      double precision tempps,vboweight,vbo,temp
      real*8 x1,x2,x3,ecut_all,smth_all
      character*2 tempname
cc      character*70 linestring,file_dum
      character*70 file_dum
      character*120 linestring

      real*8 proj(mnpw*2,mnpw*2)
      complex*16 zz_st(mnpw*2,mnpw*2,mns)
      integer iproj(mns),ipcbm(mns),iflag(mns)
      common  /comproj/zz_st,iproj,ipcbm



      type (fitdata) tempdata
      ivbo(:,:)=0

c     Read in the fitting data

      open (unit=4,file='fit.d',status='old')
      read(4,*)nstruct
      allocate(inputfile(nstruct))
      do i=1,nstruct
         read(4,*)structfile(i)
         read(4,*)ng1(i),ng2(i),ng3(i)
         read(4,*)inputfile(i)
	 read(4,*)iproj(i),ipcbm(i),file_dum
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
      allocate(controlfile(npseudofile))
      do i=1,npseudofile
         read(4,*)pseudofile(i)
         read(4,*)controlfile(i)
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

      nparam=0
      if (ifit_sig/=0) nparam=1

c     Read in the pseudopotential parameters
      
      itotps=0
      do i=1,npseudofile
         open (unit=4,file=pseudofile(i),status='old')
         read (4,*)nps
         do 10 ips=1,nps
            itotps=itotps+1
            read (4,*) psname(itotps),psnum(itotps)
            read (4,*) psvol(itotps),psfvol(itotps),psbeta(itotps),
     $           ngauss(itotps)
            do 20 j=1,ngauss(itotps)
               read (4,*) psa(j,itotps),psb(j,itotps),psc(j,itotps) 
 20         continue
 10      continue
         close (unit=4)
      end do
      write(6,*)'Read in a total of ',itotps,' pseudopotentials'
      
      do 110 ips=1,itotps
         psvol0(ips)=psvol(ips)
         psf0(ips)=psfvol(ips)
         psbeta0(ips)=psbeta(ips)
         do 120 j=1,ngauss(ips)
            psa0(j,ips)=psa(j,ips)
            psb0(j,ips)=psb(j,ips)
            psc0(j,ips)=psc(j,ips)
 120     continue
 110  continue

c     Read in the control data for each of the pseudopotentials

      temppstype=0

      do i=1,npseudofile
         open(unit=7,file=controlfile(i),form='formatted',status='old')
         read(7,*) pstype
         if (i>1 .and. pstype/=temppstype) then
            write(6,*)'Pseudopotentials must have same form!'
            stop
         end if
         temppstype=pstype

         do ips=1,2
            read(7,*)
            read(7,*)itemp
            match=0
            do jps=1,itotps
               if (itemp==psnum(jps)) then
                  read(7,*)
                  read(7,*)psSO(jps)
                  psSO0(jps)=psSO(jps)

                  read(7,*)moveSO(jps),mvpsf(jps),mvpsbeta(jps)
                  if (moveSO(jps)==1)  nparam=nparam+1
                  if (mvpsf(jps)==1)   nparam=nparam+1
                  if (mvpsbeta(jps)==1)nparam=nparam+1
                  do j=1,ngauss(jps)
                     read(7,*)mvpsa(j,jps),mvpsb(j,jps),mvpsc(j,jps)
                     if (mvpsa(j,jps).eq.1) nparam=nparam+1
                     if (mvpsb(j,jps).eq.1) nparam=nparam+1
                     if (mvpsc(j,jps).eq.1) nparam=nparam+1
                  end do     
                  match=1
               end if
            end do
            if (match.eq.0) then
               write(6,*)'Could not match pseudo no. ',itemp
               stop
            end if
         end do
      end do
      write(6,*)'There are a total of',nparam,'parameters to be fitted'
      write(6,*)

c     Read in the experimental data to be fitted for each structure

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
            case('el3v') 
               iel3v(i) = nsprop(i)
            case('el2v') 
               iel2v(i) = nsprop(i)
            case('el1c') 
               iel1c(i) = nsprop(i)
            case('mg1c') 
               img1c(i) = nsprop(i)
            case('mghh100') 
               imghh100(i) = nsprop(i)
            case('mghh111')
               imghh111(i) = nsprop(i)
            case('mglh100') 
               imglh100(i) = nsprop(i)
            case('dpa') 
               idpa(i) = nsprop(i)
            case('dpag1c') 
               idpag1c(i) = nsprop(i)
            case('dpag15v') 
               idpag15v(i) = nsprop(i)
            case('dpaX1c') 
               idpaX1c(i) = nsprop(i)
            case('dpb100') 
               idpb100(i) = nsprop(i)
            case('edso0') 
               idso0(i) = nsprop(i)
            case('edso1')
               idso1(i) = nsprop(i) 
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
            case default
               write(6,*)'Could not find a match for input ',
     &        tempdata%property
         stop
      end select
      the_data(nsprop(i),i)=tempdata
      write(6,100)tempdata
 100  format(a20,3(2x,f10.5))
      end do data

c     Read in the Valence Band Offset data for each structure

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

c     Read in the spin orbit data file

      if (iSOps.eq.1) then
         open(10,file='VG1G2.SO',form='unformatted')
         rewind(10)
         read(10)gA
         read(10)vso
         close(10)
      end if

c     Read in the positions of the atoms for the structure files

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
c            read(4,*) atomnum(iat,i),tau0(:,iat,i),
c     $           deltav(iat,i),atweight0(iat,i)
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

c     Assign the naming indicies for the bands.
c     Do this by calculating the total no. of electrons in the supercell.
      
      do is=1,nstruct
         temp=0.0
         do i=1,natoms(is)
            temp=temp+atweight0(i,is)
         end do
         ifac=int(temp*1.0001)/2
         do i=1,6
            if (iSOps.eq.0) iband(i,is)=i*ifac
            if (iSOps.eq.1) iband(i,is)=2*i*ifac
         end do

	 if(iSOps.eq.1) iband(4,is)=8*ifac
	 if(iSOps.eq.1) iband(5,is)=8*ifac+2
	 if(iSOps.eq.0) iband(4,is)=4*ifac
	 if(iSOps.eq.0) iband(5,is)=4*ifac+1

ccccccccccc iband(4,is): VBM (even num in SO), iband(5,is): CBM (even num in SO)

      end do

      end
