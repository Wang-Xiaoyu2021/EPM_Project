      subroutine readps(fileps)
  
      implicit none
      integer ips,j 
      character*20 fileps
      
      include 'epm.inc'

     open (unit=4,file=fileps,status='old')
      read (4,*) nps
      do 10 ips=1,nps
         read (4,*) psname(ips)
         read (4,*) psvol(ips),psf(ips),psbeta(ips),ngauss(ips)
         do 20 j=1,ngauss(ips)
            read (4,*) psa(j,ips),psb(j,ips),psc(j,ips) 
 20      continue
 10   continue
        close (unit=4)
	
       do 110 ips=1,nps
          psvol0(ips)=psvol(ips)
          psf0(ips)=psf(ips)
          psbeta0(ips)=psbeta(ips)
	  do 120 j=1,ngauss(ips)
             psa0(j,ips)=psa(j,ips)
             psb0(j,ips)=psb(j,ips)
             psc0(j,ips)=psc(j,ips)
120        continue
          psSO0(ips)=psSO(ips)
110     continue

        end
