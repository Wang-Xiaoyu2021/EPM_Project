      subroutine vcell (vr,vc,inG)

      implicit none
      include 'epm.inc'

      integer inG,i
      double precision vr(1:mngrid)
      complex*16 vc(1:mngrid)

      integer n1,n2,n3
      integer ig,ig1,ig2,ig3,iat,itp,ips
      double precision volume,twopi,gtau,vtp,x1,x2,x3
      double precision gn1,gn2,gn3,g
      complex*16 stp


      call settp ()
      call recvec (a,b1,b2,b3,omega)
      volume=omega
      twopi=8.0d0*atan(1.0d0)
      
      ig=0
      do 10 ig3=0,ng3(is)-1
         n3=ig3
         if (ig3.gt.ng3(is)/2) n3=ig3-ng3(is)
         do 20 ig2=0,ng2(is)-1
            n2=ig2
            if (ig2.gt.ng2(is)/2) n2=ig2-ng2(is)
            do 30 ig1=0,ng1(is)-1
               n1=ig1
               if (ig1.gt.ng1(is)/2) n1=ig1-ng1(is)

               ig=ig+1
               gn1=n1*b1(1)+n2*b2(1)+n3*b3(1)  ! (n)(b),(vect multply matrix)
               gn2=n1*b1(2)+n2*b2(2)+n3*b3(2)  ! obtain a recip vector g
               gn3=n1*b1(3)+n2*b2(3)+n3*b3(3)
               g=sqrt(gn1**2+gn2**2+gn3**2)

               vc(ig)=dcmplx(0.0d0,0.0d0)
               do 40 itp=1,ntypes
                  stp=dcmplx(0.0d0,0.0d0)
                  
                  do 50 iat=1,tpnat(itp)
                     x1=tau(1,tpindx(iat,itp),is)
                     x2=tau(2,tpindx(iat,itp),is)
                     x3=tau(3,tpindx(iat,itp),is)
                     gtau=alat*(gn1*x1+gn2*x2+gn3*x3)
                     stp=stp+dcmplx(cos(gtau),-sin(gtau))
     $                    *atweight(tpindx(iat,itp),is)   !tpindx(:) has no is

 50               continue

ccccccc this following formf is not quit correct, if the cell has a original
ccccccc local strain and a additional cell deformation. But for all small
ccccccc deformation, it is O.K
ccccccc This is really crazy. It assumes that there is only one atom for
ccccccc each type. Not easy to understand the logic here. But for small 
ccccccc deformations, it works. Who wrote the original crapy code ?! 
ccccccc commented by LWW. 

                  call formf (tppsd(itp),g,vtp)  ! get value of v(g)   
                   vc(ig)=vc(ig)+stp*vtp/volume

 40            continue
 30         continue
 20      continue
 10   continue

C     Since vgen gives pseudopotentials in G space,
C     do FFT if necessary  

      if(inG .ge. 0) then
         call fft(vc,ng1(is),ng2(is),ng3(is),1)

         do 130 i=1,ngrid(is) 
            vr(i)=dreal(vc(i))
 130     continue
       endif

      end



