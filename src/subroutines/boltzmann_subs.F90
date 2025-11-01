function ehlft(n,nocpk,elft,hlft) !lifetime (s) in boltzmann equation

	implicit none
	
	integer :: n,nocpk
	real :: elft,hlft
	real :: ehlft
	
	if ( n .gt. nocpk) then
	
		ehlft = elft
	
	else
	
		ehlft = hlft
	
	end if

end function ehlft

function dfermidist(mu,temp,energy) !minus 1 multiplied for the Energy derivative of Fermi-Dirac distribution


	implicit none
	
	real,parameter :: keV = 8.617330350E-5	
	real :: dfermidist
	real :: mu
	real :: energy
	real :: temp
	real :: fdis
	
	real :: aux0,aux1,aux2
	
	aux0 = (energy-mu)/(keV*temp)
	aux1 = exp((energy-mu)/(keV*temp))
	aux2 = (keV*temp)*(aux1+1.0)**2
	
	!fdis = 1.0/(aux1+1.00)
	
	!dfermidist = (1.00/(keV*temp))*(fdis)*(1.00-fdis)
	
	if ( aux2 .gt. 1.0E+37) then
	 dfermidist = 0.0 
	else
	 dfermidist = aux1/aux2
	end if

end function dfermidist


function df1(mu,temp,energy)

	implicit none
	
	real,parameter :: keV = 8.617330350E-5	
	real :: df1
	real :: mu
	real :: energy
	real :: temp
	real :: fdis
	
	real :: aux0,aux1,aux2
	
	aux0 = (energy-mu)/(keV*temp)
	aux1 = exp((energy-mu)/(keV*temp))
	aux2 = (keV*temp)*(aux1+1.0)**2
	
	df1 = (1.0/(4.0*keV*temp))*(1.0/cosh((energy-mu)/(2.0*keV*temp)))**2

end function df1

subroutine eleccond(nsteps,energies,tdfij,mu,temp,sij) !calculate the electrical conductivity in 1/ohm * angstrom


	implicit none
	real,parameter :: keV = 8.617330350E-5  !Boltzmann's constant eV/K
	real,parameter :: echarge = 1.602176620898E-19 !electron charge in coulomb	
	integer :: i,j,nsteps
	real,dimension(nsteps) :: energies
	real :: dfermidist,mu,temp,df1
	real,dimension(nsteps,3) :: tdfij
	real,dimension(3) :: sij,aux0,aux1,aux2
	
	real :: de
	
	de = abs(energies(2)-energies(1))

	sij = 0.0
	
	do i=1,nsteps-1
	
		!aux0 = tdfij(i,1)!*dfermidist(mu,temp,energies(i))
		!aux1 = tdfij(i+1,1)!*dfermidist(mu,temp,energies(i+1))

		aux0(1) = tdfij(i,1)*dfermidist(mu,temp,energies(i))
		aux1(1) = tdfij(i+1,1)*dfermidist(mu,temp,energies(i+1))
		aux2(1) = (aux0(1)+aux1(1))*(de/2.0)	
		sij(1) = sij(1) + aux2(1)

		aux0(2) = tdfij(i,2)*dfermidist(mu,temp,energies(i))
		aux1(2) = tdfij(i+1,2)*dfermidist(mu,temp,energies(i+1))
		aux2(2) = (aux0(2)+aux1(2))*(de/2.0)
		sij(2) = sij(2) + aux2(2)

		aux0(3) = tdfij(i,3)*dfermidist(mu,temp,energies(i))
		aux1(3) = tdfij(i+1,3)*dfermidist(mu,temp,energies(i+1))
		aux2(3) = (aux0(3)+aux1(3))*(de/2.0)		
		sij(3) = sij(3) + aux2(3)

		!aux0(4) = tdfij(i,4)*dfermidist(mu,temp,energies(i))
		!aux1(4) = tdfij(i+1,4)*dfermidist(mu,temp,energies(i+1))
		!aux2(4) = (aux0(4)+aux1(4))*(de/2.0)		
		!sij(4) = sij(4) + aux2(4)

		!aux0(5) = tdfij(i,5)*dfermidist(mu,temp,energies(i))
		!aux1(5) = tdfij(i+1,5)*dfermidist(mu,temp,energies(i+1))
		!aux2(5) = (aux0(5)+aux1(5))*(de/2.0)		
		!sij(5) = sij(5) + aux2(5)

		!aux0(6) = tdfij(i,6)*dfermidist(mu,temp,energies(i))
		!aux1(6) = tdfij(i+1,6)*dfermidist(mu,temp,energies(i+1))
		!aux2(6) = (aux0(6)+aux1(6))*(de/2.0)		
		!sij(6) = sij(6) + aux2(6)										
	end do
	
	       sij = sij*echarge*echarge

end subroutine eleccond

subroutine electermcond(nsteps,energies,tdfij,mu,temp,kij) !calculate electron termal conductivity in watt/kelvin * angstrom


	implicit none
	real,parameter :: keV = 8.617330350E-5  !Boltzmann's constant eV/K
	integer :: i,j,nsteps
	real,dimension(nsteps) :: energies
	real :: dfermidist,mu,temp
	real,dimension(nsteps,3) :: tdfij
	real,dimension(3) :: kij
	
	real :: de,aux0,aux1		
	
	de = abs(energies(2)-energies(1))
	
	kij = 0.0
	
	do i=1,3
	 do j=1,nsteps-1	

		aux0 = dfermidist(mu,temp,energies(j))*(energies(j)-mu)*(energies(j)-mu)*tdfij(j,i)
		aux1 = dfermidist(mu,temp,energies(j+1))*(energies(j+1)-mu)*(energies(j+1)-mu)*tdfij(j+1,i)

		kij(i) = kij(i) + ((aux0+aux1)*(de/2.0))

	 end do
	 
	 	kij(i) = kij(i)/temp
	end do
	
	

end subroutine electermcond

subroutine seebeck(nsteps,energies,tdfij,mu,temp,sij,seij) !calculate seebeck coefficient in volt/kelvin


	implicit none
	real,parameter :: keV = 8.617330350E-5  !Boltzmann's constant eV/K
	real,parameter :: echarge = 1.602176620898E-19 !electron charge in coulomb
	integer :: i,j,nsteps
	real,dimension(nsteps) :: energies
	real :: dfermidist,mu,temp
	real,dimension(nsteps,3) :: tdfij
	real,dimension(3) :: sij,seij
	
	real,dimension(3) :: aux3
	real :: de,aux1,aux0
	
	de = abs(energies(2)-energies(1))	
	
	aux3(1) = echarge/(temp*sij(1))
	aux3(2) = echarge/(temp*sij(2))
	aux3(3) = echarge/(temp*sij(3))
	!aux3(4) = echarge/(temp*sij(4))
	!aux3(5) = echarge/(temp*sij(5))
	!aux3(6) = echarge/(temp*sij(6))					
	
	seij = 0.0
	
	do i=1,nsteps-1
	
		aux0 = dfermidist(mu,temp,energies(i))*(energies(i)-mu)*tdfij(i,1)
		aux1 = dfermidist(mu,temp,energies(i+1))*(energies(i+1)-mu)*tdfij(i+1,1)
		
		seij(1) = seij(1) + ((aux0+aux1)*(de/2.0))
		
		aux0 = dfermidist(mu,temp,energies(i))*(energies(i)-mu)*tdfij(i,2)
		aux1 = dfermidist(mu,temp,energies(i+1))*(energies(i+1)-mu)*tdfij(i+1,2)
		
		seij(2) = seij(2) + ((aux0+aux1)*(de/2.0))
		
		aux0 = dfermidist(mu,temp,energies(i))*(energies(i)-mu)*tdfij(i,3)
		aux1 = dfermidist(mu,temp,energies(i+1))*(energies(i+1)-mu)*tdfij(i+1,3)
		
		seij(3) = seij(3) + ((aux0+aux1)*(de/2.0))
		
		!aux0 = dfermidist(mu,temp,energies(i))*(energies(i)-mu)*tdfij(i,4)
		!aux1 = dfermidist(mu,temp,energies(i+1))*(energies(i+1)-mu)*tdfij(i+1,4)
		
		!seij(4) = seij(4) + ((aux0+aux1)*(de/2.0))
		
		!aux0 = dfermidist(mu,temp,energies(i))*(energies(i)-mu)*tdfij(i,5)
		!aux1 = dfermidist(mu,temp,energies(i+1))*(energies(i+1)-mu)*tdfij(i+1,5)
		
		!seij(5) = seij(5) + ((aux0+aux1)*(de/2.0))
		
		!aux0 = dfermidist(mu,temp,energies(i))*(energies(i)-mu)*tdfij(i,6)
		!aux1 = dfermidist(mu,temp,energies(i+1))*(energies(i+1)-mu)*tdfij(i+1,6)
		
		!seij(6) = seij(6) + ((aux0+aux1)*(de/2.0))										
	
	end do
	
	seij(1) = seij(1)*aux3(1)
	seij(2) = seij(2)*aux3(2)
	seij(3) = seij(3)*aux3(3)
	!seij(4) = seij(4)*aux3(4)
	!seij(5) = seij(5)*aux3(5)
	!seij(6) = seij(6)*aux3(6)							
	

end subroutine seebeck


subroutine sigmas(nsteps,energies,tdfij,mu,temp,sigsij) !calculate sigma*seebecj coefficient in volt/kelvin

	
	implicit none
	real,parameter :: keV = 8.617330350E-5  !Boltzmann's constant eV/K
	real,parameter :: echarge = 1.602176620898E-19 !electron charge in coulomb
	integer :: i,j,nsteps
	real,dimension(nsteps) :: energies
	real :: dfermidist,mu,temp
	real,dimension(nsteps,3) :: tdfij
	real,dimension(3) :: sigsij
	
	real,dimension(3) :: aux3
	real :: de
	real,dimension(3) :: aux1,aux0
	
	de = abs(energies(2)-energies(1))	
	
	aux3(1) = echarge/(temp)
	aux3(2) = echarge/(temp)
	aux3(3) = echarge/(temp)
	!aux3(4) = echarge/(temp*sij(4))
	!aux3(5) = echarge/(temp*sij(5))
	!aux3(6) = echarge/(temp*sij(6))					
	
	sigsij = 0.0
	
	do i=1,nsteps-1
	
		aux0(1) = dfermidist(mu,temp,energies(i))*(mu-energies(i))*tdfij(i,1)
		aux1(1) = dfermidist(mu,temp,energies(i+1))*(mu-energies(i+1))*tdfij(i+1,1)
		
		sigsij(1) = sigsij(1) + ((aux0(1)+aux1(1))*(de/2.0))
		
		aux0(2) = dfermidist(mu,temp,energies(i))*(mu-energies(i))*tdfij(i,2)
		aux1(2) = dfermidist(mu,temp,energies(i+1))*(mu-energies(i+1))*tdfij(i+1,2)
		
		sigsij(2) = sigsij(2) + ((aux0(2)+aux1(2))*(de/2.0))
		
		aux0(3) = dfermidist(mu,temp,energies(i))*(mu-energies(i))*tdfij(i,3)
		aux1(3) = dfermidist(mu,temp,energies(i+1))*(mu-energies(i+1))*tdfij(i+1,3)
		
		sigsij(3) = sigsij(3) + ((aux0(3)+aux1(3))*(de/2.0))
		
		!aux0 = dfermidist(mu,temp,energies(i))*(energies(i)-mu)*tdfij(i,4)
		!aux1 = dfermidist(mu,temp,energies(i+1))*(energies(i+1)-mu)*tdfij(i+1,4)
		
		!seij(4) = seij(4) + ((aux0+aux1)*(de/2.0))
		
		!aux0 = dfermidist(mu,temp,energies(i))*(energies(i)-mu)*tdfij(i,5)
		!aux1 = dfermidist(mu,temp,energies(i+1))*(energies(i+1)-mu)*tdfij(i+1,5)
		
		!seij(5) = seij(5) + ((aux0+aux1)*(de/2.0))
		
		!aux0 = dfermidist(mu,temp,energies(i))*(energies(i)-mu)*tdfij(i,6)
		!aux1 = dfermidist(mu,temp,energies(i+1))*(energies(i+1)-mu)*tdfij(i+1,6)
		
		!seij(6) = seij(6) + ((aux0+aux1)*(de/2.0))										
	
	end do
	
	sigsij(1) = sigsij(1)*aux3(1)
	sigsij(2) = sigsij(2)*aux3(2)
	sigsij(3) = sigsij(3)*aux3(3)
	!seij(4) = seij(4)*aux3(4)
	!seij(5) = seij(5)*aux3(5)
	!seij(6) = seij(6)*aux3(6)							
	

end subroutine sigmas



subroutine tdf(systype,rlat,w90basis,nk,eigenvalues,nocpk,sme,energy,elft,hlft,vx,vy,vz,tij) !transport distribution function in 1/(angstrom * s * eV)

	implicit none 
	
	integer :: i,j
	integer :: w90basis,nk
	real :: sme,energy
	integer,dimension(nk) :: nocpk
	real,dimension(nk,w90basis) :: eigenvalues
	real,dimension(3,3) :: rlat
	real :: ehlft
	real,dimension(3) :: elft,hlft
	real,dimension(w90basis,nk) :: vx,vy,vz
	real,dimension(3) :: tij
	real :: gaussian
	
	real :: vcell,divfactor
	character(len=4) :: systype
	
	call vcell3D(rlat,vcell)
	
	if (systype .eq. 'NP') then
	
	divfactor = 2.0/(vcell*real(nk))
	
	else
	
	divfactor = 1.0/(vcell*real(nk))
	
	end if
	
	tij(1) = 0.0
	tij(2) = 0.0
	tij(3) = 0.0
	!tij(4) = 0.0
	!tij(5) = 0.0
	!tij(6) = 0.0
	
	do i=1,w90basis
	
	   do j=1,nk
	   
	     tij(1) = tij(1) + vx(i,j)*vx(i,j)*gaussian((energy-eigenvalues(j,i)),sme)
	     tij(2) = tij(2) + vy(i,j)*vy(i,j)*gaussian((energy-eigenvalues(j,i)),sme)
	     tij(3) = tij(3) + vz(i,j)*vz(i,j)*gaussian((energy-eigenvalues(j,i)),sme)
	   
	    !if (energy .gt. 0.0) then
	    
	    ! tij(1) = tij(1) + vx(i,j)*vx(i,j)*elft(1)*gaussian((energy-eigenvalues(j,i)),sme)
	    ! tij(2) = tij(2) + vy(i,j)*vy(i,j)*elft(2)*gaussian((energy-eigenvalues(j,i)),sme)
	    ! tij(3) = tij(3) + vz(i,j)*vz(i,j)*elft(3)*gaussian((energy-eigenvalues(j,i)),sme)
	    
	    !else
	    
	    ! tij(1) = tij(1) + vx(i,j)*vx(i,j)*hlft(1)*gaussian((energy-eigenvalues(j,i)),sme)
	    ! tij(2) = tij(2) + vy(i,j)*vy(i,j)*hlft(2)*gaussian((energy-eigenvalues(j,i)),sme)
	    ! tij(3) = tij(3) + vz(i,j)*vz(i,j)*hlft(3)*gaussian((energy-eigenvalues(j,i)),sme)
	    
	    !end if
		
	     !tij(1) = tij(1) + vx(i,j)*vx(i,j)*ehlft(i,nocpk(j),elft(1),hlft(1))*gaussian((energy-eigenvalues(j,i)),sme)
	     !tij(2) = tij(2) + vy(i,j)*vy(i,j)*ehlft(i,nocpk(j),elft(2),hlft(2))*gaussian((energy-eigenvalues(j,i)),sme)
	     !tij(3) = tij(3) + vz(i,j)*vz(i,j)*ehlft(i,nocpk(j),elft(3),hlft(3))*gaussian((energy-eigenvalues(j,i)),sme)			
	     !tij(4) = tij(4) + vx(i,j)*vy(i,j)*ehlft(i,nocpk(j),elft,hlft)*gaussian((energy-eigenvalues(j,i)),sme)
	     !tij(5) = tij(5) + vx(i,j)*vz(i,j)*ehlft(i,nocpk(j),elft,hlft)*gaussian((energy-eigenvalues(j,i)),sme)
	     !tij(6) = tij(6) + vy(i,j)*vz(i,j)*ehlft(i,nocpk(j),elft,hlft)*gaussian((energy-eigenvalues(j,i)),sme)		
		
	   end do
	
	end do
	
	tij(1) = tij(1)*divfactor
	tij(2) = tij(2)*divfactor
	tij(3) = tij(3)*divfactor
	!tij(4) = tij(4)*divfactor
	!tij(5) = tij(5)*divfactor
	!tij(6) = tij(6)*divfactor

end subroutine tdf



subroutine bndvel(kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		    ihopmatrices,wf,vx,vy,vz) !calculate band velocity (angstrom/s)


	implicit none
	real,parameter :: hbar = 6.582119569E-16 !hbar planck's constant
	integer :: i,j,k
	integer :: w90basis,nvec
	real :: kx,ky,kz
	integer,dimension(nvec) :: ffactor
	real,dimension(3,3) :: rlat
	real,dimension(nvec,3) :: rvec
	
	complex,dimension(w90basis,w90basis) :: hx,hy,hz
	complex,dimension(w90basis) :: wf
	
	real,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices

	
	real :: vx,vy,vz	
	
	call  hlm(kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		    ihopmatrices,hx,hy,hz)
		    
	call sandwich(w90basis,wf,hx,wf,vx)
	call sandwich(w90basis,wf,hy,wf,vy)
	call sandwich(w90basis,wf,hz,wf,vz)	
	
	vx = vx/hbar
	vy = vy/hbar
	vz = vz/hbar 	    

end subroutine bndvel

