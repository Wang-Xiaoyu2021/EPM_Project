               subroutine formf (kps,g,vtp)

        implicit none
        include 'epm.inc'

        integer kps
        double precision g,vtp
	double precision power

	power = psa(2,kps) ! to give more freedom if you like to.

        vtp = psvol(kps)*(g**power-psa(1,kps))/
     &        (psb(1,kps)*dexp(psc(1,kps)*g*g)-1)


        return
        end

