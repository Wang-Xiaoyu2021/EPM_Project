      subroutine monitor(m,n,xc,fvecc,fjacc,ljc,s,igrade,
     &                   niter,nf,iw,liw,w,lw)

      implicit none
       include 'epm.inc'

       integer liw,ljc,igrade,ips,lw,i,j,k
       integer niter,nf,iw(liw),n,m,ii,iis
       double precision fvecc(m),fjacc(ljc,n),s(n)
       double precision w(lw),xc(n),temp

       write(6,*)
       write(6,*)'Iteration no. ',niter
       write(6,*)'Values of the pseudopotential parameters'
       write(6,*)'----------------------------------------'

        do 70 ips=1,itotps
           write(6,*)psname(ips),psnum(ips)
           write(6,200)psvol0(ips),psfvol(ips),psbeta(ips),ngauss(ips)
 200       format(f7.3,2(2x,f20.10),4x,i2)
	   do 80 j=1,ngauss(ips)
	      write (6,*)psa(j,ips),abs(psb(j,ips)),abs(psc(j,ips))
 80        continue
 70     continue

        write(6,*)
        write(6,*)'Spin Orbit Parameters'
        write(6,*)'---------------------'
        if(iSOps .eq.1) then
           do ips=1,itotps
           write(6,101) psnum(ips),psSO(ips)
 101       format(I7,F8.5)
           enddo
        endif

        write(6,*)'Weight of Kinetic Energy',sigma
        write(6,*)
        
        ii=0
        do iis=1,nstruct
           write(6,*)'Results for structure ',structfile(iis)
           write(6,*)'       Property   Current Value   ',
     &          'Error   Weighted Rel Err'
           write(6,*)'       --------   -------------   ',
     &          ' -----   ----------------'
           write(6,*)
           if(isOps .eq. 0) then
               do i=1,nsprop(iis)
                 ii=ii+1
                  if (the_data(i,iis)%weight > 0.0d0) then
                    write(6,100)the_data(i,iis)%property,
     $                   calcvalue(i,iis),calcvalue(i,iis)-
     $                   the_data(i,iis)%nsovalue,fvecc(ii)
                 end if
              end do
           else
              do i=1,nsprop(iis)
                 ii=ii+1
                 if (the_data(i,iis)%weight > 0.0d0) then
                    write(6,100)the_data(i,iis)%property,
     $                   calcvalue(i,iis),calcvalue(i,iis)
     $                   -the_data(i,iis)%sovalue,fvecc(ii)
                 end if
              end do
           endif
           write(6,*)
        end do
 100    format(a20,f10.5,f10.5,f10.5)

        temp=0.0d0
        ii=0
        do iis=1,nstruct
           do i=1,nsprop(iis)
              ii=ii+1
              temp=temp+fvecc(ii)**2
           end do
        end do
        write(6,*)'Total Error = ',temp

        end



