subroutine emissionpl(dimbse,omega,refractive,dipole,exctemp,eta,emin,eplm)

	implicit none
	
	integer :: dimbse
	real :: emin,eta,exctemp
	real,dimension(dimbse,4) :: dipole
	real :: omega
	real,dimension(3) :: refractive
	
	real,parameter :: pi=acos(-1.)
	complex,parameter :: imag = cmplx(0.0,1.0)
	real,parameter :: hbar = 6.582119569E-16 !hbar planck's constant
	real,parameter :: keV = 8.617330350E-5  !Boltzmann's constant eV/K
	real,parameter :: hc3 = 1.0 !177.714E8 !ev*angstrom/femtoseconds**2
	
	integer :: i
	complex,dimension(3) :: epl
	real,dimension(3) :: eplm, aux2
	real :: aux1
	complex :: auxc
	
	epl = 0.0

	
	do i=1,dimbse
	
		aux1 = exp(-(dipole(i,1)-emin)/(kev*exctemp))
		
		auxc = 1.0/(omega-dipole(i,1)+(eta*imag))
		!auxc = eta/((omega-dipole(i,1)**2+eta**2))
		
		aux2(1) = (omega**3)*refractive(1)*dipole(i,2)/(pi*pi*hc3)
		aux2(2) = (omega**3)*refractive(2)*dipole(i,3)/(pi*pi*hc3)
		aux2(3) = (omega**3)*refractive(3)*dipole(i,4)/(pi*pi*hc3)
		
		epl(1) = epl(1)+ aux2(1)*aux1*auxc
		epl(2) = epl(2)+ aux2(2)*aux1*auxc
		epl(3) = epl(3)+ aux2(3)*aux1*auxc	
	
	end do


		eplm = -aimag(epl) !unit: angstrom*eV*femtosecond
		!eplm = real(epl) !unit: angstrom*eV*femtosecond

end subroutine emissionpl
