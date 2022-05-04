subroutine lcount(w90basis,autovetor,dcont,pcont)

	implicit none

	integer :: i

	integer :: w90basis,wbasisaux
	double complex,dimension(w90basis) :: autovetor

	double precision :: dcont,pcont

	dcont = 0.0
	pcont = 0.0

	wbasisaux=w90basis/2

	do i=1,5

		dcont=dcont+REAL(autovetor(i)*conjg(autovetor(i)))
		dcont=dcont+REAL(autovetor(wbasisaux+i)*conjg(autovetor(wbasisaux+i)))

	end do

	do i=1,6

		pcont=pcont+REAL(autovetor(5+i)*conjg(autovetor(5+i)))
		pcont=pcont+REAL(autovetor(wbasisaux+5+i)*conjg(autovetor(wbasisaux+5+i)))

	end do



end subroutine lcount



subroutine spincont(w90basis,autovetor,cup,cdown,dft,systype)

	implicit none

	integer :: i,j

	integer :: w90basis,wbasisaux
	double complex,dimension(w90basis) :: autovetor

	double precision :: cup,cdown
	character(len=1) :: dft
	character(len=4) :: systype	

	cup=0d0
	cdown=0d0
	
	wbasisaux=w90basis/2
	
	if (systype .eq. "NP") then 
	
		cup=0d0
		cdown=0d0
	
	else if (systype .eq. "SP") then
	
	 do i=1,wbasisaux

		cup=cup+REAL(autovetor(i)*conjg(autovetor(i)))

		cdown=cdown+REAL(autovetor(wbasisaux+i)*conjg(autovetor(wbasisaux+i)))


	 end do	
	
	
	else if ((systype .eq. "SOC") .and. (dft .eq. "V") ) then
	
	 do i=1,wbasisaux

		cup=cup+REAL(autovetor(i)*conjg(autovetor(i)))

		cdown=cdown+REAL(autovetor(wbasisaux+i)*conjg(autovetor(wbasisaux+i)))


	 end do
	
	else
	
	 do i=1,wbasisaux

		cup=cup+REAL(autovetor((2*i)-1)*conjg(autovetor((2*i)-1)))

		cdown=cdown+REAL(autovetor(2*i)*conjg(autovetor(2*i)))


	 end do
	
	
	end if
	
	

	
end subroutine spincont

subroutine endos(rlat,dimse,n1,n2,n3,en,etb,fosc,sme,ints)


	implicit none

	character(len=2) :: systype

	double precision :: en
	integer :: dimse !dimensão do sistema => nkpts*nbands
	double precision :: sme
	double precision ::  ints
	integer :: n1,n2,n3,i,j

	double precision :: vc

	double complex :: intsaux

	double precision,dimension(3,3) :: rlat

	double precision,dimension(dimse) :: etb
	double precision,dimension(dimse) :: fosc
	
	double precision :: caux,dimk,gaussian

	double precision,parameter :: pi=acos(-1.)
	

	dimk = dble(n1*n2*n3)

	caux=(1.0/dimk)


	ints=0.


	do i=1,dimse

		
		intsaux = fosc(i)*gaussian(en-etb(i),sme)


		ints=ints+intsaux


	end do


	ints=ints*caux!*(1./elux)


end subroutine endos

subroutine layercont(w90basis,wf,orbweight,lc)

	implicit none

	integer :: w90basis
	double complex,dimension(w90basis) :: wf
	double precision,dimension(w90basis) :: orbweight 
	integer :: i,j

	double complex :: lcaux
	double precision :: lc

	lcaux = 0.0


	do i=1,w90basis

		lcaux = wf(i)*orbweight(i)*conjg(wf(i))+lcaux

	end do	


	lc  = real(lcaux)


end subroutine


