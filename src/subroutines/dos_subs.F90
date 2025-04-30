subroutine lcount(w90basis,autovetor,dcont,pcont)

	implicit none

	integer :: i

	integer :: w90basis,wbasisaux
	complex,dimension(w90basis) :: autovetor

	real :: dcont,pcont

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



subroutine spincont(w90basis,autovetor,cup,cdown,dft,systype,ovptb)

	implicit none

	integer :: i,j

	integer :: w90basis,wbasisaux
	
	complex,dimension(w90basis,w90basis)  :: ovptb
	
	complex,dimension(w90basis) :: autovetor
	
	complex,dimension(w90basis) :: eigvup,eigvdown
	
	complex,dimension(w90basis) :: orbup,orbdown
	
	complex   :: cup1,cdown1

	real :: cup,cdown
	character(len=1) :: dft
	character(len=4) :: systype	

	cup=0d0
	cdown=0d0
	
	wbasisaux=w90basis/2
	
	orbup = 0d0
	orbdown = 0d0
	
	do i=1,wbasisaux
	
		orbup(i) = 1.0
		orbdown(wbasisaux+i) = 1.0
	
	end do
		
	
	if (systype .eq. "NP") then 
	
		cup=0d0
		cdown=0d0
	
	else if (systype .eq. "SP") then
	
	
	
		if (dft .eq. "S") then
		
			do i=1,w90basis
			
				eigvup(i)= autovetor(i)*orbup(i)
				eigvdown(i)= autovetor(i)*orbdown(j)
			
			end do

			call sandwich(w90basis,eigvup,ovptb,eigvup,cup1)
			call sandwich(w90basis,eigvdown,ovptb,eigvdown,cdown1)
			
			cup = real(cup1)
			cdown= real(cdown1)
		
		else
	
	 	do i=1,wbasisaux

			cup=cup+REAL(autovetor(i)*conjg(autovetor(i)))

			cdown=cdown+REAL(autovetor(wbasisaux+i)*conjg(autovetor(wbasisaux+i)))


	 	end do	
	
		end if
	
	else if ((systype .eq. "SOC") .and. (dft .eq. "V") ) then
	
	 do i=1,wbasisaux

		cup=cup+REAL(autovetor(i)*conjg(autovetor(i)))

		cdown=cdown+REAL(autovetor(wbasisaux+i)*conjg(autovetor(wbasisaux+i)))


	 end do
	
	else if ((systype .eq. "SOC") .and. (dft .eq. "W") ) then
	
	 do i=1,wbasisaux

		cup=cup+REAL(autovetor((2*i)-1)*conjg(autovetor((2*i)-1)))

		cdown=cdown+REAL(autovetor(2*i)*conjg(autovetor(2*i)))


	 end do
	 
	 else if ((systype .eq. "SOC") .and. (dft .eq. "S") ) then
	 
	 			do i=1,w90basis
			
				eigvup(i)= autovetor(i)*orbup(i)
				eigvdown(i)= autovetor(i)*orbdown(j)
			
			end do

			call sandwich(w90basis,eigvup,ovptb,eigvup,cup1)
			call sandwich(w90basis,eigvdown,ovptb,eigvdown,cdown1)
			
			cup = real(cup1)
			cdown= real(cdown1)
	 
	 else
	 
	  do i=1,wbasisaux

		cup=cup+REAL(autovetor(i)*conjg(autovetor(i)))

		cdown=cdown+REAL(autovetor(wbasisaux+i)*conjg(autovetor(wbasisaux+i)))

	  end do
	
	
	end if
	
	

	
end subroutine spincont

subroutine endos(rlat,dimse,n1,n2,n3,en,etb,fosc,sme,ints)


	implicit none

	character(len=2) :: systype

	real :: en
	integer :: dimse !dimensão do sistema => nkpts*nbands
	real :: sme
	real ::  ints
	integer :: n1,n2,n3,i,j

	real :: vc

	complex :: intsaux

	real,dimension(3,3) :: rlat

	real,dimension(dimse) :: etb
	real,dimension(dimse) :: fosc
	
	real :: caux,dimk,gaussian

	real,parameter :: pi=acos(-1.)
	

	dimk = dble(n1*n2*n3)

	caux=(1.0/dimk)


	ints=0.


	do i=1,dimse

		
		intsaux = fosc(i)*gaussian(en-etb(i),sme)


		ints=ints+intsaux


	end do


	ints=ints*caux!*(1./elux)


end subroutine endos

subroutine layercont(w90basis,dft,wf,orbweight,ovptb,lc)

	implicit none

	integer :: w90basis
	complex,dimension(w90basis) :: wf,wfc
	real,dimension(w90basis) :: orbweight 
	integer :: i,j

	complex :: lcaux
	real :: lc
	character(len=1) :: dft	
	complex,dimension(w90basis,w90basis) :: ovptb


	select case (dft)
	
	
	case ("S")
	
	 do i=1,w90basis
	 
	 	wf(i)=wf(i)*orbweight(i)
	 	
	 
	 end do
	 
	  
	 call sandwich(w90basis,wf,ovptb,wf,lcaux)
	
	 lc  = real(lcaux)
	 
	
	case default
	

	 lcaux = 0.0



	 do i=1,w90basis

		lcaux = wf(i)*orbweight(i)*conjg(wf(i))+lcaux

	 end do	


	 lc  = real(lcaux)

	end select

end subroutine


