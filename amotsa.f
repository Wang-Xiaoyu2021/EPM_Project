       FUNCTION amotsa(p,y,psum,mp,np,ndim,pb,yb,funk,ihi,yhi,fac)
       INTEGER ihi,mp,ndim,np,NMAX
       DOUBLE PRECISION amotsa,fac,yb,yhi,
     & p(mp,np),pb(np),psum(np),y(mp),funk, rnd
       PARAMETER (NMAX=200)
       EXTERNAL funk
C USES funk,ran1

       INTEGER idum,j
       DOUBLE PRECISION fac1,fac2,tt,yflu,ytry,ptry(NMAX)
       COMMON /ambsa/ tt,idum
       fac1=(1.-fac)/ndim
       fac2=fac1-fac
       do j=1,ndim
           ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
       enddo

       ytry=funk(ptry)
       if (ytry.le.yb) then
           do j=1,ndim
               pb(j)=ptry(j)
           enddo
           yb=ytry
       endif
       call random_number(rnd)
       yflu=ytry-tt*log(rnd) 
       if (yflu.lt.yhi) then
           y(ihi)=ytry
           yhi=yflu
           do j=1,ndim
               psum(j)=psum(j)-p(ihi,j)+ptry(j)
               p(ihi,j)=ptry(j)
           enddo
       endif
       amotsa=yflu
       return
       END
