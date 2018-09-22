c  Parameters

      integer mproperty,mnat,mntp,mnps,mng1,mng2,mng3,mngrid
      integer mnw1,mnw2,mnw3,mnpw,mnset,mnkp,mnpar,mngauss,mvbo

c      parameter (mproperty=50,mnat=200,mntp=10,mnps=10)
      parameter (mproperty=60,mnat=200,mntp=10,mnps=10)
      parameter (mng1=16,mng2=20,mng3=40,mngrid=mng1*mng2*mng3)
      parameter (mnw1=mng1/2,mnw2=mng2/2,mnw3=mng3/2)
      parameter (mnpw=700,mnset=mnpw*mnpw,mnkp=10)
      parameter (mngauss=8,mnpar=40,mvbo=4)
      double precision Ryd
      parameter(Ryd=13.605698d0)
      integer,parameter :: mns=5 

c The indicies for the experimental properties

      integer ieg1v(mns),ievbm(mns),ievbm_1(mns),ieg1c(mns),
     &  ieg15c(mns),idso0(mns),iex1v(mns),iex3v(mns),iex5v(mns),
     &  iex1c(mns),iex3c(mns),iel2v(mns),iel3v(mns),iel1c(mns),
     &  img1c(mns),imghh100(mns),imglh100(mns),idpa(mns),
     &  idpb100(mns),idpd111(mns),idso1(mns),iel1v(mns),
     &  imghh111(mns),ivbo(mvbo,mns),idpag1c(mns),idpag15v(mns),
     &  idpaX1c(mns),iegp3pc(mns),iegm3pc(mns)

      common/exptindex/ieg1v,ievbm,ievbm_1,ieg1c,ieg15c,idso0,iex1v,
     &                iex3v,iex5v,iex1c,iex3c,iel2v,iel3v,iel1c,img1c,
     &                 imghh100,imglh100,idpa,idpb100,idpd111,idso1,
     &                 iel1v,ivbo,imghh111,idpag1c,idpag15v,idpaX1c,
     &                 iegp3pc,iegm3pc

c Indicies for the different band
      
      integer iband(6,mns)
      common /bands/iband

c New kinetic energy term
	real*8 sigma,sigma0,delta0
	integer ifit_sig	
	common /sig/sigma,sigma0,delta0,ifit_sig

c Experimental Data

      type fitdata
        sequence
        character*8 property
        double precision nsovalue
        double precision sovalue 
        double precision weight
      end type fitdata

      double precision calcvalue(mproperty,mns)
      type (fitdata) the_data(mproperty,mns)      

      integer nproperty
      common /exptdata/calcvalue,the_data,nproperty

c VBO data

      double precision eps_paral(mvbo,mns),eps_perp(mvbo,mns),
     &                 vbovbm(mvbo,mns)
      common /vbodata/ eps_paral,eps_perp,vbovbm

c Pseudopotential Data

       character*2 psname(1:mnps)
       integer pstype, nps,ngauss(1:mnps)
       double precision psvol(1:mnps),psfvol(1:mnps),
     &                  psf(1:mnps),psbeta(1:mnps),
     &                  psa(1:mngauss,1:mnps),psb(1:mngauss,1:mnps),
     &                  psc(1:mngauss,1:mnps)
       
       common /pseudo/ psvol,psfvol,psf,psbeta,psa,psb,psc,
     &                 pstype,nps,psname,ngauss

c Reference pseudopotential data

       double precision psvol0(1:mnps),psf0(1:mnps),psbeta0(1:mnps),
     &                  psa0(1:mngauss,1:mnps),psb0(1:mngauss,1:mnps),
     &                  psc0(1:mngauss,1:mnps)
	integer itotps,psnum(mnps)
       common /psref/ psvol0,psf0,psbeta0,psa0,psb0,psc0,itotps,psnum

c Spin Orbit Data

      double precision psSO(1:mnat),psSO0(1:mnat),gA(100),vso(100,100)
      integer iSOps
      common /SOpseudo/psSO,psSO0,gA,vso,iSOps

c Data on which paramters to fit

      integer mvpsa,mvpsb,mvpsc,moveSO,mvpsf,mvpsbeta
      common /mvpseudo/mvpsa(1:mngauss,1:mnps),mvpsb(1:mngauss,1:mnps),
     &                 mvpsc(1:mngauss,1:mnps),moveSO(1:mnps),
     &                 mvpsf(1:mnps),mvpsbeta(1:mnps)

c Data defining the cells in each structure file

      double precision alat0(mns),a0(3,3,mns)
      double precision alat,a(3,3),b1(3),b2(3),b3(3),omega
      common /cell/alat0,alat,a0,a,b1,b2,b3,omega

c Data defining the grid

      integer ng1(mns),ng2(mns),ng3(mns),ngrid(mns)
      integer pfft(-mnw1:mnw1,-mnw2:mnw2,-mnw3:mnw3)
      double precision ecut(mns),Smth(mns)
      common /grid/ Smth,ecut,ng1,ng2,ng3,ngrid,pfft

      integer ngpw(1:mnpw,1:3),npw
      double precision gpw(1:mnpw,1:3),wg(1:mnpw)
      common /basis/ gpw,wg,ngpw,npw

c Data defining the atoms

      double precision tau(3,mnat,mns),tau0(3,mnat,mns),
     &                 atweight(mnat,mns),atweight0(mnat,mns),
     &                 deltav(mnat,mns)
      integer natoms(mns),atomnum(mnat,mns) 
      common /atoms/ tau,tau0,natoms,atomnum,atweight,
     &               atweight0,deltav

      integer ntypes,tpnat(1:mntp),tppsd(1:mntp),
     &        tpindx(1:mnat,1:mntp),tpnum(mntp)
      common /atmtp/ntypes,tpnat,tppsd,tpindx,tpnum
 
c Which structure we are currently dealing with

      integer is,nsprop(mns),nstruct
      character*30 structfile(mns)
      common /structure/ is,nsprop,nstruct,structfile
