       SUBROUTINE amebsa(p,y,mp,np,ndim,pb,yb,ftol,funk,iter,temptr)
       INTEGER iter,mp,ndim,np,NMAX
       DOUBLE PRECISION ftol,temptr,yb,p(mp,np),pb(np),y(mp),funk
       PARAMETER (NMAX=200)
       EXTERNAL funk
! USES amotsa,funk,ran1

       INTEGER i,idum,ihi,ilo,j,m,n
       DOUBLE PRECISION rtol,sum,swap,tt,yhi,ylo,ynhi,ysave,
     & yt,ytry,psum(NMAX), amotsa, rnd
       COMMON /ambsa/ tt,idum
       tt=-temptr
       idum = -5
1      do n=1,ndim 
          sum=0. 
          do m=1,ndim+1
              sum=sum+p(m,n)
          enddo
          psum(n)=sum
       enddo  
2      ilo=1 
       ihi=2
       call random_number(rnd)
       ylo=y(1)+tt*log(rnd) 
       ynhi=ylo
       call random_number(rnd)
       yhi=y(2)+tt*log(rnd)
       if (ylo.gt.yhi) then
           ihi=1
           ilo=2
           ynhi=yhi
           yhi=ylo
           ylo=ynhi
       endif
       do i=3,ndim+1
           call random_number(rnd) 
           yt=y(i)+tt*log(rnd) 
           if(yt.le.ylo) then
               ilo=i
               ylo=yt
           endif
           if(yt.gt.yhi) then
               ynhi=yhi
               ihi=i
               yhi=yt
           else if(yt.gt.ynhi) then
               ynhi=yt
               endif
       enddo
       rtol=2.*abs(yhi-ylo)/(abs(yhi)+abs(ylo))

       if (rtol.lt.ftol.or.iter.lt.0) then 
           swap=y(1)
           y(1)=y(ilo)
           y(ilo)=swap
           do n=1,ndim
               swap=p(1,n)
               p(1,n)=p(ilo,n)
               p(ilo,n)=swap
           enddo
           return
       endif
       iter=iter-2
       ytry=amotsa(p,y,psum,mp,np,ndim,pb,yb,funk,ihi,yhi,-1.0d0)
       if (ytry.le.ylo) then
           ytry=amotsa(p,y,psum,mp,np,ndim,pb,yb,funk,ihi,yhi,
     & 2.0d0)
       else if (ytry.ge.ynhi) then
           ysave=yhi
           ytry=amotsa(p,y,psum,mp,np,ndim,pb,yb,funk,ihi,yhi,
     & 0.5d0)
           if (ytry.ge.ysave) then 
               do i=1,ndim+1
                   if(i.ne.ilo)then
                       do j=1,ndim
                           psum(j)=0.5d0*(p(i,j)+p(ilo,j))
                           p(i,j)=psum(j)
                       enddo
                       y(i)=funk(psum)
                   endif
               enddo
               iter=iter-ndim
               goto 1
           endif
       else
           iter=iter+1
       endif
       goto 2
       END
