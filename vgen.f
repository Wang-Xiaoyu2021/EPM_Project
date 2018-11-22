           subroutine vgen (vc)

C==== declaration needs to be checked
         implicit none
         include 'epm.inc'

	 integer n1,n2,n3
	 integer ig,ig1,ig2,ig3,i,iat,itp,ips
	 double precision volume,twopi,gtau,vtp,x1,x2,x3
	 double precision gn1,gn2,gn3,g
	 complex*16 stp,vc

	 dimension vc(1:mngrid)

	 volume=omega
         twopi=8.0d0*atan(1.0d0)

         ig=0
         do 10 ig3=0,ng3-1
	    n3=ig3
	    if (ig3.gt.ng3/2) n3=ig3-ng3
	    do 20 ig2=0,ng2-1
	       n2=ig2
	       if (ig2.gt.ng2/2) n2=ig2-ng2
	       do 30 ig1=0,ng1-1
		  n1=ig1
		  if (ig1.gt.ng1/2) n1=ig1-ng1

		  ig=ig+1
		  gn1=n1*b1(1)+n2*b2(1)+n3*b3(1)
		  gn2=n1*b1(2)+n2*b2(2)+n3*b3(2)
		  gn3=n1*b1(3)+n2*b2(3)+n3*b3(3)
		  g=sqrt(gn1**2+gn2**2+gn3**2)

	          vc(ig)=dcmplx(0.0d0,0.0d0)
                  do 40 itp=1,ntypes
	             stp=dcmplx(0.0d0,0.0d0)
		
	             do 50 iat=1,tpnat(itp)
		        x1=tau(1,tpindx(iat,itp))
		        x2=tau(2,tpindx(iat,itp))
		        x3=tau(3,tpindx(iat,itp))
		        gtau=alat*(gn1*x1+gn2*x2+gn3*x3)
 		        stp=stp+dcmplx(cos(gtau),-sin(gtau))
50                   continue

                     call formf (tppsd(itp),g,vtp)
                     vc(ig)=vc(ig)+stp*vtp/volume

40                continue
30             continue
20          continue
10       continue

         return
         end                  
