module constants

	implicit none

	real,parameter :: c =299792458.0 !speed of light, m/s
	real,parameter :: h =6.62607004081E-34 !Planck's constant J*s (W)
	real,parameter :: heV =4.135667516E-15  !Planck's constant eV*s
	real,parameter :: k =1.3806485279E-23  !Boltzmann's constant J/K
	real,parameter :: keV =8.617330350E-5  !Boltzmann's constant eV/K
	real,parameter :: e =1.602176620898E-19 !Coulomb

	real,parameter :: jev=6.24150913E18 !conversion factor joule to eV
	real,parameter :: pi=acos(-1.0)
	real,parameter :: tkv=293.15 !temperatura em kelvin
	real,parameter :: thick=1E-6 !espessura em metros

	!intervalo para o grafico do SLME vc thickness
	real,parameter :: thickmin=0E-7
	real,parameter :: thickmax=1E-6
	

	!energy window solar spectra in eV
	!real,parameter :: e0 = 0.31
	!real,parameter :: ef = 4.42801

	!parametros da subrotina de minimização
	real,parameter :: delta_tol= 10**(-7.0d0)
	real,parameter :: delta= 0.2D+00
	integer,parameter :: k_max = 200 

end module constants
