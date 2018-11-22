        subroutine hdiag (hd,npw,ev,psi)

        implicit none
	integer npw,naux
	double precision ev,aux
	complex*16 hd,psi

        include 'param.d'

        dimension hd(1:mnset)
	dimension ev(1:mnpw),psi(1:mnpw,1:mnpw)
        dimension aux(1:4*mnpw)

        naux=4*mnpw
	call zhpev (20,hd,ev,psi,mnpw,npw,aux,naux)

	return
	end
