        subroutine formf (kps,g,vtp)

        implicit none
        include 'epm.inc'

        integer kps,j
        double precision g,vtp,power
        double precision vlr,vsr,vtp1,vtp2,c
        double precision g2,gmax2

c
c  pstype = 1, gaussians
c
        if(pstype .eq. 1) then    ! EPM with gaussians 

        vlr=psvol(kps)*(1.0d0+psf(kps)*exp(-psbeta(kps)*g*g))
        vsr=0.0d0
c Use the absolute values of the b and c parameters as negative values
c correspond to explosive exponentials and large g components of the potential
c respectively.

        do 10 j=1,ngauss(kps)
c        vsr=vsr+psa(j,kps)*exp(-abs(psb(j,kps))*(g-abs(psc(j,kps)))**2)
c        vsr=vsr+psa(j,kps)*exp(-abs(psc(j,kps))*(g-abs(psb(j,kps)))**2)
        vsr=vsr+psa(j,kps)*exp(-(psc(j,kps))*(g-(psb(j,kps)))**2)
10      continue

        vtp=vlr*vsr
        endif
c
c  pstype = 2,5, exponential
c
        if(pstype .eq. 2 .or. pstype.eq.5) then   ! type 5 used in SLCBB

        power = psb(2,kps) ! to give more freedom if you like to.

        vlr = psvol(kps)*(1.0d0+psf(kps)*exp(-psbeta(kps)*g*g))
        vtp = vlr*psa(1,kps)*(g*g-psb(1,kps))/
     &        (psc(1,kps)*dexp(psa(2,kps)*g*g)-1.0d0)

        endif
c
c  pstype = 51, exponential formula of cohen
c
       if(pstype .eq. 51) then
        vlr = psvol(kps)*(1.0d0+psf(kps)*exp(-psbeta(kps)*g*g))
        vtp = vlr*psa(1,kps)*(g*g-psb(1,kps))/
     &        (dexp(psc(1,kps)*(g*g-psa(2,kps)))+1.0d0)
        endif
c
c  pstype = 52, exponential formula of Schluter from PRB39, 7974(1989).
c  
        if(pstype .eq. 52) then
        vlr = psvol(kps)*(1.0d0+psf(kps)*exp(-psbeta(kps)*g*g))
        vtp = vlr*psa(1,kps)*(g*g-psb(1,kps))/
     &        (dexp(psc(1,kps)*(g*g-psa(2,kps)))+1.0d0)
        vtp = vtp*(0.5d0*tanh((psb(2,kps)-g*g)/psc(2,kps))+0.5d0)
        endif
c
c  pstype = 53, exponential formula of type5 with more freedom
c
       if(pstype.eq.53) then   ! type 5 used in SLCBB
        power = psb(2,kps) ! to give more freedom if you like to.
        vlr = psvol(kps)*(1.0d0+psf(kps)*exp(-psbeta(kps)*g*g))
        vtp = vlr*psa(1,kps)*(g*g-psb(1,kps))/
     &        (psc(1,kps)*dexp(psa(2,kps)*g*g)-1.0d0)
        endif
c
c pstype = 54, exponential formula of Goano from JAP88,6467(2000). 
c
       if(pstype.eq.54) then
         g2=g*g
         gmax2=psbeta(kps)*psbeta(kps)
         if(g2.le.gmax2) then
           vsr = 1.d0-dexp(-(psb(1,kps)*gmax2)**psc(1,kps))
           vsr = 1.d0/vsr
           vlr = 1.d0-dexp(-(psb(1,kps)*(gmax2-g2))**psc(1,kps))
           vtp = -psa(1,kps)+(psa(2,kps)+psa(1,kps))*(1.d0-vsr*vlr)
         else
           vlr = dexp(-(psb(2,kps)*(g2-gmax2))**psc(2,kps))
           vtp = psa(2,kps)*vlr
         endif
       endif
c
c  check
c
        if(pstype.ne.2.and.pstype.ne.1.and.pstype.ne.5.and. 
     &      pstype.ne.51. and. pstype.ne.52.and.pstype.ne.53
     &      .and. pstype.ne.54) then
        write(6,*) "pstype not found in form.f, stop", pstype
        stop
        endif


        return
        end


