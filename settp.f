      subroutine settp ()

        implicit none
        include 'epm.inc'

	integer iat,itp,ips,nt,chk

100     format (' Error in subroutine settp: 
     &  pseudopotential not found')

        ntypes=1
        tpnum(1)=atomnum(1,is)
        tpnat(1)=1 
	tpindx(1,1)=1

        do 10 iat=2,natoms(is)
	   chk=0
	   do 20 nt=1,ntypes
	      if (atomnum(iat,is).eq.tpnum(nt)) then 
		 tpnat(nt)=tpnat(nt)+1 
		 tpindx(tpnat(nt),nt)=iat
		 chk=1
              end if
20         continue	 
           if (chk.eq.0) then
	      ntypes=ntypes+1
	      tpnum(ntypes)=atomnum(iat,is)
	      tpnat(ntypes)=1
	      tpindx(1,ntypes)=iat
           end if
10      continue

        do 30 itp=1,ntypes
	   ips=0
40	   ips=ips+1
	   if (psnum(ips).eq.tpnum(itp)) then
	      tppsd(itp)=ips
           else
	      go to 40
	   end if
30      continue

        end
