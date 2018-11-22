      subroutine penalty(fvecc,nprop)

      implicit none
      include 'epm.inc'

      integer ips,j,i,nprop
      double precision fvecc(mproperty)

      if(pstype .eq. 1) then
c Add an extra cost that penalises values of c greater than 3.
        do ips=1,nps
           do j=1,ngauss(ips)
              if (mvpsc(j,ips) .eq. 1 .and. psc(j,ips).gt.4.0d0) then
                 do i=1,nprop
                    fvecc(i)=fvecc(i)+fvecc(i)*0.1d0*(psc(j,ips)-4)**2
                 end do
              end if
           end do
        end do
        
c Add an extra cost that penalises values of b greater than 10.
        do ips=1,nps
           do j=1,ngauss(ips)
              if (mvpsc(j,ips) .eq. 1 .and. psb(j,ips).gt.10.0d0) then
                 do i=1,nprop
                fvecc(i)=fvecc(i)+fvecc(i)*0.01d0*(psb(j,ips)-10)**2
                 end do
              end if
           end do
        end do
       endif

        return
	end
