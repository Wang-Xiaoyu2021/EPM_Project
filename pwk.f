      subroutine pwk (kp,vr,vc,inG,ev,proj,iproj,zz_st)

       implicit double precision (a-h,o-z)

       include 'epm.inc'

       integer i,n,ipw,jpw,ik,ierr,inG
       double precision kp(3),vr(1:mngrid),ev(1:mnpw*2),kpc(3)
       complex*16 vc(1:mngrid)
       real*8 proj(mnpw*2,mnpw*2)
       complex*16 zz_st(mnpw*2,mnpw*2)
       integer iproj

300    format (///' k =',3f8.5,'     Number of PWs =',i5/)
310    format ('    Eigenvalue','   Energy (eV)'/) 
320    format (6x,i4,4x,f14.6) 

       kpc(1)=kp(1)*b1(1)+kp(2)*b2(1)+kp(3)*b3(1)
       kpc(2)=kp(1)*b1(2)+kp(2)*b2(2)+kp(3)*b3(2)
       kpc(3)=kp(1)*b1(3)+kp(2)*b2(3)+kp(3)*b3(3)

       call pwset (kpc)


       if(iSOps .eq. 0) call solve_nso(kpc,vr,vc,inG,ev,
     &       proj,iproj,zz_st)
       if(iSOps .eq. 1) call solve_so(kpc,vr,vc,inG,ev,
     &       proj,iproj,zz_st)
 

       end

