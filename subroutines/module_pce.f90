module constants

	implicit none

	double precision,parameter :: c =299792458.0 !speed of light, m/s
	double precision,parameter :: h =6.62607004081E-34 !Planck's constant J*s (W)
	double precision,parameter :: heV =4.135667516E-15  !Planck's constant eV*s
	double precision,parameter :: k =1.3806485279E-23  !Boltzmann's constant J/K
	double precision,parameter :: keV =8.617330350E-5  !Boltzmann's constant eV/K
	double precision,parameter :: e =1.602176620898E-19 !Coulomb

	double precision,parameter :: jev=6.24150913E18 !conversion factor joule to eV
	double precision,parameter :: pi=acos(-1.0)
	double precision,parameter :: tkv=293.15 !temperatura em kelvin
	double precision,parameter :: thick=1E-6 !espessura em metros

	!intervalo para o grafico do SLME vc thickness
	double precision,parameter :: thickmin=0E-7
	double precision,parameter :: thickmax=1E-6
	

	!energy window solar spectra in eV
	!double precision,parameter :: e0 = 0.31
	!double precision,parameter :: ef = 4.42801

	!parametros da subrotina de minimização
	double precision,parameter :: delta_tol= 10**(-7.0d0)
	double precision,parameter :: delta= 0.2D+00
	integer,parameter :: k_max = 200 

end module constants
