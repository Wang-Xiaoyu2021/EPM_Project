      subroutine monitor(m,n,xc,fvecc,fjacc,ljc,s,igrade,
     &                   niter,nf,iw,liw,w,lw)

      implicit none
      include 'epm.inc'

       integer liw,ljc,igrade,ips,lw,i,j,k
       integer niter,nf,iw(liw),n,m,ii,iis
       double precision fvecc(m),fjacc(ljc,n),s(n)
       double precision w(lw),xc(n),temp,calcvalue_tmp
       real*8 Ecut_t,Smth_t,scalkin_t
       common /EcutSmth/Ecut_t,Smth_t,scalkin_t

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc xc(k) is the most updated one, the psSO(ips)'s are not
ccc  psSO(ips) is the parameter for the last trial direction run. 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       write(6,*)
       write(6,*)'Iteration no. ',niter
       write(6,*)'Values of the pseudopotential parameters'
       write(6,*)'----------------------------------------'

       

c        do 70 ips=1,itotps
c           write(6,*)psname(ips),psnum(ips)
c           write(6,*) " ---- S0,  af,  a0,   a1,    a2,   a3 -------"
c           write(6,200) psSO(ips),psfvol(ips),psa(1,ips),psb(1,ips)
c     &     ,psc(1,ips), psa(2,ips)
 200       format(6(ES15.7,1x))
c 70     continue

        write(6,*)'Weight of Kinetic Energy',sigma
        write(6,*)

*********************************** write out the param files
        open (unit=4,file="pseudo_all.out")
        rewind(4)
c        if(pstype.ne.5) then
c        write(6,*) "pstype.ne.5, stop, monitor.f"
c        stop
c        endif

        write (4,102) itotps, pstype
c        write(4,101) Ecut_t,Smth_t,scalkin_t
        write(4,101) Ecut_t,Smth_t,sigma
        write(6,101) Ecut_t,Smth_t,sigma
101     format(3(f8.5,1x), "          ! Ecut, Smth, scalkin")
102     format(2(i5,3x), "          ! num_pot, pot. type")
c type5 EPM
         if(pstype.eq.5) then
         do ips=1,itotps
         write(4,*)
         write(6,*)
         write(4,*) psname(ips),psnum(ips), "      ! Atom, atomic_num"
         write(6,*) psname(ips),psnum(ips), "      ! Atom, atomic_num"
         write(4,111) psSO(ips),moveSO(ips)
         write(6,111) psSO(ips),moveSO(ips)
         write(4,112) psfvol(ips),mvpsf(ips)
         write(6,112) psfvol(ips),mvpsf(ips)
         write(4,113) psa(1,ips),mvpsa(1,ips)
         write(6,113) psa(1,ips),mvpsa(1,ips)
         write(4,114) psb(1,ips),mvpsb(1,ips)
         write(6,114) psb(1,ips),mvpsb(1,ips)
         write(4,115) psc(1,ips),mvpsc(1,ips)
         write(6,115) psc(1,ips),mvpsc(1,ips)
         write(4,116) psa(2,ips),mvpsa(2,ips)
         write(6,116) psa(2,ips),mvpsa(2,ips)
c        write(4,124) psbeta(ips),mvpsbeta(ips)
c        write(6,124) psbeta(ips),mvpsbeta(ips)
         enddo

c type1 EPM
         else if(pstype.eq.1) then
         do ips=1,itotps
         write(4,*)
         write(6,*)
         write(4,*) psname(ips),psnum(ips), "      ! Atom, atomic_num"
         write(6,*) psname(ips),psnum(ips), "      ! Atom, atomic_num"
         write(4,118) psvol(ips),ngauss(ips)
         write(6,118) psvol(ips),ngauss(ips)
         write(4,111) psSO(ips),moveSO(ips)
         write(6,111) psSO(ips),moveSO(ips)
         write(4,112) psfvol(ips),mvpsf(ips)
         write(6,112) psfvol(ips),mvpsf(ips)
         write(4,119) psbeta(ips), mvpsbeta(ips)
         write(6,119) psbeta(ips), mvpsbeta(ips)
         do 300 j=1,ngauss(ips)
         write(4,120)  psa(j,ips),psb(j,ips),psc(j,ips),
     &       mvpsa(j,ips), mvpsb(j,ips), mvpsc(j,ips)
         write(6,120)  psa(j,ips),psb(j,ips),psc(j,ips),
     &       mvpsa(j,ips), mvpsb(j,ips), mvpsc(j,ips)
300      continue
         enddo

c type52 EPM
         else if(pstype.eq.52) then
         do ips=1,itotps
         write(4,*)
         write(6,*)
         write(4,*) psname(ips),psnum(ips), "      ! Atom, atomic_num"
         write(6,*) psname(ips),psnum(ips), "      ! Atom, atomic_num"
         write(4,111) psSO(ips),moveSO(ips)
         write(6,111) psSO(ips),moveSO(ips)
         write(4,112) psfvol(ips),mvpsf(ips)
         write(6,112) psfvol(ips),mvpsf(ips)
         write(4,113) psa(1,ips),mvpsa(1,ips)
         write(6,113) psa(1,ips),mvpsa(1,ips)
         write(4,114) psb(1,ips),mvpsb(1,ips)
         write(6,114) psb(1,ips),mvpsb(1,ips)
         write(4,115) psc(1,ips),mvpsc(1,ips)
         write(6,115) psc(1,ips),mvpsc(1,ips)
         write(4,116) psa(2,ips),mvpsa(2,ips)
         write(6,116) psa(2,ips),mvpsa(2,ips)
         write(4,121) psb(2,ips),mvpsb(2,ips)
         write(6,121) psb(2,ips),mvpsb(2,ips)
         write(4,122) psc(2,ips),mvpsc(2,ips)
         write(6,122) psc(2,ips),mvpsc(2,ips)
         enddo

c type53 EPM
         else if(pstype.eq.53) then
         do ips=1,itotps
         write(4,*)
         write(6,*)
         write(4,*) psname(ips),psnum(ips), "      ! Atom, atomic_num"
         write(6,*) psname(ips),psnum(ips), "      ! Atom, atomic_num"
         write(4,111) psSO(ips),moveSO(ips)
         write(6,111) psSO(ips),moveSO(ips)
         write(4,112) psfvol(ips),mvpsf(ips)
         write(6,112) psfvol(ips),mvpsf(ips)
         write(4,123) psbeta(ips),mvpsbeta(ips)
         write(6,123) psbeta(ips),mvpsbeta(ips)
         write(4,113) psa(1,ips),mvpsa(1,ips)
         write(6,113) psa(1,ips),mvpsa(1,ips)
         write(4,114) psb(1,ips),mvpsb(1,ips)
         write(6,114) psb(1,ips),mvpsb(1,ips)
         write(4,115) psc(1,ips),mvpsc(1,ips)
         write(6,115) psc(1,ips),mvpsc(1,ips)
         write(4,116) psa(2,ips),mvpsa(2,ips)
         write(6,116) psa(2,ips),mvpsa(2,ips)
         enddo

c type54 EPM
         else if(pstype.eq.54) then
         do ips=1,itotps
         write(4,*)
         write(6,*)
         write(4,*) psname(ips),psnum(ips), "      ! Atom, atomic_num"
         write(6,*) psname(ips),psnum(ips), "      ! Atom, atomic_num"
         write(4,111) psSO(ips),moveSO(ips)
         write(6,111) psSO(ips),moveSO(ips)
         write(4,112) psfvol(ips),mvpsf(ips)
         write(6,112) psfvol(ips),mvpsf(ips)
         write(4,119) psbeta(ips),mvpsbeta(ips)
         write(6,119) psbeta(ips),mvpsbeta(ips)
         write(4,113) psa(1,ips),mvpsa(1,ips)
         write(6,113) psa(1,ips),mvpsa(1,ips)
         write(4,114) psb(1,ips),mvpsb(1,ips)
         write(6,114) psb(1,ips),mvpsb(1,ips)
         write(4,115) psc(1,ips),mvpsc(1,ips)
         write(6,115) psc(1,ips),mvpsc(1,ips)
         write(4,116) psa(2,ips),mvpsa(2,ips)
         write(6,116) psa(2,ips),mvpsa(2,ips)
         write(4,121) psb(2,ips),mvpsb(2,ips)
         write(6,121) psb(2,ips),mvpsb(2,ips)
         write(4,122) psc(2,ips),mvpsc(2,ips)
         write(6,122) psc(2,ips),mvpsc(2,ips)
         enddo

c not defined EPM format
         else
         write(6,*) "unregonized EPM format, stop"
         stop
         endif
   
111      format(ES20.13,4x,i3,10x, "!SO param")
112      format(ES20.13,4x,i3,10x, "!af param for strain")
113      format(ES20.13,4x,i3,10x, "!a0 param")
114      format(ES20.13,4x,i3,10x, "!a1 param")
115      format(ES20.13,4x,i3,10x, "!a2 param")
116      format(ES20.13,4x,i3,10x, "!a3 param")
117      format(ES20.13,4x,i3,10x, "!kin_scal param")
118      format(ES20.13,4x,i3,10x, "!norm. atomic volume, ngauss")
119      format(ES20.13,4x,i3,10x, "!beta param")
120      format(3(ES20.13,1x),3(1x,i2),2x, "!gaussian's param")
121      format(ES20.13,4x,i3,10x, "!a4 param")
122      format(ES20.13,4x,i3,10x, "!a5 param")
123      format(ES20.13,4x,i3,10x, "!beta param for strain")
124      format(ES20.13,4x,i3,10x, "!psbeta param")
         write(4,117) sigma,ifit_sig
         write(4,*)
         if(iSOps.eq.0)  write(4,*)  "non-SO calc."
         if(iSOps.eq.1)  write(4,*)  "SO calc."
         write(4,*) "**************************"
         

        
        ii=0
        do iis=1,nstruct
           write(6,*)'Results for structure ',structfile(iis)
           write(6,'(A,3(1x,f8.5))')'X-point:', 
     &        xpt(1,iis),xpt(2,iis),xpt(3,iis)
           write(6,'(a,i2,a,a20)') "The ", iband(4,iis), 
     &        "th state is the VBM of ", structfile(iis)
           write(6,*)'  Property       Ideal_Value   ',
     &       'fitted_Value  Error    Weighted_Err**2   Weight'
           write(6,*)'  --------             ------   ',
     &          '-------     -----     ------------'
           write(6,*)
           write(4,*)'Results for structure ',structfile(iis)
           write(4,'(A,3(1x,f8.5))')'X-point:', 
     &        xpt(1,iis),xpt(2,iis),xpt(3,iis)
           write(4,'(a,i2,a,a20)') "The ", iband(4,iis), 
     &        "th state is the VBM of ", structfile(iis)
           write(4,*)'  Property       Ideal_Value   ',
     &        'fitted_Value  Error    Weighted_Err**2   Weight'
           write(4,*)'  --------             ------   ',
     &          '-------     -----     ------------'
           write(4,*)
           if(isOps .eq. 0) then
               do i=1,nsprop(iis)
                 ii=ii+1
                  if (the_data(i,iis)%weight > 0.0d0) then
            calcvalue_tmp=the_data(i,iis)%nsovalue+fvecc(ii)*
     $         the_data(i,iis)%nsovalue/the_data(i,iis)%weight

                    write(6,100)the_data(i,iis)%property,
     $                   the_data(i,iis)%nsovalue,
     $                   calcvalue_tmp,calcvalue_tmp-
     $                   the_data(i,iis)%nsovalue,fvecc(ii)**2,
     $                   the_data(i,iis)%weight
                    write(4,100)the_data(i,iis)%property,
     $                   the_data(i,iis)%nsovalue,
     $                   calcvalue_tmp,calcvalue_tmp-
     $                   the_data(i,iis)%nsovalue,fvecc(ii)**2,
     $                   the_data(i,iis)%weight
                 end if
              end do
           else
              do i=1,nsprop(iis)
                 ii=ii+1
                 if (the_data(i,iis)%weight > 0.0d0) then

            calcvalue_tmp=the_data(i,iis)%sovalue+fvecc(ii)*
     $         the_data(i,iis)%sovalue/the_data(i,iis)%weight

                    write(6,100)the_data(i,iis)%property,
     $                   the_data(i,iis)%sovalue,
     $                   calcvalue_tmp,calcvalue_tmp
     $                   -the_data(i,iis)%sovalue,fvecc(ii)**2
     $                    ,the_data(i,iis)%weight
                    write(4,100)the_data(i,iis)%property,
     $                   the_data(i,iis)%sovalue,
     $                   calcvalue_tmp,calcvalue_tmp
     $                   -the_data(i,iis)%sovalue,fvecc(ii)**2
     $                    ,the_data(i,iis)%weight
                 end if
              end do
           endif
           write(6,*)
        end do
c 100    format(a20,f10.5,f10.5,f10.5,f13.8)
c 100    format(a20,f10.5,f10.5,f10.5,f13.8,f13.8)
 100    format(a16,f10.5,f13.8,f13.8,f17.8,f10.3)
 

        temp=0.0d0
        ii=0
        do iis=1,nstruct
           do i=1,nsprop(iis)
              ii=ii+1
              temp=temp+fvecc(ii)**2
           end do
        end do
        write(6,*)'Total Error = ',temp
        write(4,*)'Total Error = ',temp

        close(4)



        end




