      
      subroutine setpar (xp,np,lambda)

        implicit none
	integer np,i,j
	double precision xp,lambda,fact

        include 'param.d'

        dimension xp(1:mnpar,1:mnpar+1)

        fact=1.0d0+lambda

        do 10 i=1,np
	   if (xp(i,1).eq.0.0d0) then
	      xp(i,1)=1.0d0
           end if
10      continue

        do 20 j=2,np+1
	   do 30 i=1,np
              xp(i,j)=xp(i,1)
30         continue
20      continue

        do 40 j=2,np+1
           xp(j-1,j)=xp(j-1,j)*fact
40      continue


        return
	end
