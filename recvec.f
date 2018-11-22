       subroutine recvec (a,b1,b2,b3,omega)

       implicit none
       double precision a(3,3),b1(3),b2(3),b3(3),omega
       double precision a12(3),a23(3),a31(3),twopi,tpim,omgs

       twopi=8.0d0*atan(1.0d0)

       a12(1)=a(1,2)*a(2,3)-a(1,3)*a(2,2)  ! a1y*a2z-a1z*a2y
       a12(2)=a(1,3)*a(2,1)-a(1,1)*a(2,3)  ! a1z*a2x-a1x*a2z
       a12(3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)  ! a1x*a2y-a1y*a2x
       a23(1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)  ! a1,a2,a3 is the base vector of
       a23(2)=a(2,3)*a(3,1)-a(2,1)*a(3,3)  ! primary cell not supercell.
       a23(3)=a(2,1)*a(3,2)-a(2,2)*a(3,1)  
       a31(1)=a(3,2)*a(1,3)-a(3,3)*a(1,2)
       a31(2)=a(3,3)*a(1,1)-a(3,1)*a(1,3)
       a31(3)=a(3,1)*a(1,2)-a(3,2)*a(1,1)

       omgs=a(1,1)*a23(1)+a(1,2)*a23(2)+a(1,3)*a23(3) ! a1*(a2xa3)
       tpim=twopi/omgs

       b1(1)=tpim*a23(1)                  ! b1,b2,b3 is the recip base vector of
       b1(2)=tpim*a23(2)                  ! primary cell (i.e. fcc)
       b1(3)=tpim*a23(3)
       b2(1)=tpim*a31(1)
       b2(2)=tpim*a31(2)
       b2(3)=tpim*a31(3)
       b3(1)=tpim*a12(1)
       b3(2)=tpim*a12(2)
       b3(3)=tpim*a12(3)
 
       omega=abs(omgs)

       end
