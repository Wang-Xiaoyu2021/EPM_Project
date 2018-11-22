      subroutine pwset (kpc)
cccccc plane wave basis set

       implicit none

       include 'epm.inc'

       integer nm1,nm2,nm3,n1,n2,n3,ipw
       double precision kpc(3)
       double precision bb1,bb2,bb3
       double precision gn1,gn2,gn3,ggn
       double precision x,pi,onemsmth

       onemsmth = 1.0d0-Smth(is)

       pi=4.d0*datan(1.d0)

       bb1=b1(1)**2+b1(2)**2+b1(3)**2
       bb2=b2(1)**2+b2(2)**2+b2(3)**2
       bb3=b3(1)**2+b3(2)**2+b3(3)**2
       nm1=int(sqrt(ecut(is)/bb1))+1
       nm2=int(sqrt(ecut(is)/bb2))+1
       nm3=int(sqrt(ecut(is)/bb3))+1
     
       if (nm1.gt.ng1(is)/2.or.nm2.gt.ng2(is)/2
     $      .or.nm3.gt.ng3(is)/2) then
          write (6,*) ' Error '
          stop
       end if
       ipw=0
       do 10 n3=-nm3,nm3 
          do 20 n2=-nm2,nm2
            do 30 n1=-nm1,nm1

                gn1=n1*b1(1)+n2*b2(1)+n3*b3(1)
                gn2=n1*b1(2)+n2*b2(2)+n3*b3(2)
                gn3=n1*b1(3)+n2*b2(3)+n3*b3(3)
                ggn=(kpc(1)+gn1)**2+
     &          (kpc(2)+gn2)**2+(kpc(3)+gn3)**2
                if (ggn.le.ecut(is)) then
                   ipw=ipw+1
                   ngpw(ipw,1)=n1
                   ngpw(ipw,2)=n2
                   ngpw(ipw,3)=n3
                   gpw(ipw,1)=gn1
                   gpw(ipw,2)=gn2
                   gpw(ipw,3)=gn3
                   if(ggn.lt.ecut(is)*Smth(is)) then
                   wg(ipw)=1.d0
                   else
                   x=(ggn-Smth(is)*ecut(is))/(onemsmth*ecut(is))*pi
                   wg(ipw)=(dcos(x)+1.d0)*0.5d0
                   endif
                end if
30           continue
20        continue
10     continue

       npw=ipw

       if(npw.gt.mnpw) then
       write(6,*) "npw>mnpw, stop", npw, mnpw
       stop
       endif

       return
       end

