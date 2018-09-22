         subroutine fft (work,ng1,ng2,ng3,sign)

ccccccccccc this subroutine is never called in the main program. 
ccccccccccc So, it is never used in multifit.x

        implicit none
        integer lwrk
        parameter (lwrk=1000)
	integer ng1,ng2,ng3,sign,i,isign
        real*8 scale
	complex*16 work(ng1*ng2*ng3)
        real*8 chdr(ng1*ng2*ng3),chdi(ng1*ng2*ng3)    ! dynamic allocation in F90
        real*8 wrk(lwrk)

        do i=1,ng1*ng2*ng3
        chdr(i)=dreal(work(i))
        chdi(i)=dimag(work(i))
        enddo


        isign=-sign

	call cfft(ng1,ng2,ng3,chdr,chdi,wrk,lwrk,isign)    ! the sign ?

        scale=1.d0
        if(sign.eq.-1) scale=1.d0/(ng1*ng2*ng3)

        do i=1,ng1*ng2*ng3
        work(i)=dcmplx(chdr(i),chdi(i))*scale
        enddo

	return
	end
