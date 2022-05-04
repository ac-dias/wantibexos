!potencial 2D Keldysh (Ridolfi) DOI:https://doi.org/10.1103/PhysRevB.97.205409
function v2dk(kpt1,kpt2,ediel,rlat,ngrid,lc,tolr)

	implicit none

	double precision,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	double precision,parameter:: pi=acos(-1.)
	double precision,parameter :: alpha1 = 1.76
	double precision,parameter :: alpha2 = 1.0
	double precision,parameter :: alpha3 = 0d0

	integer,dimension(3) :: ngrid
	double precision,dimension(3) :: kpt1,kpt2
	double precision,dimension(3,3) :: rlat
	double precision :: tolr
	double precision,dimension(3) :: ediel
	double precision :: lc
	double precision :: gridaux1,ed


	double precision :: a0,modk,vc
	double precision :: r0,vbz,auxi
	double precision :: v2dk

	call alat2D(rlat,a0)
	call modvec(kpt1,kpt2,modk)
	call vcell2D(rlat,vc)

	r0= ((ediel(2)-1.0)*lc)/(ediel(1)+ediel(3))

	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))

	gridaux1 = dble(ngrid(1)*ngrid(2))

	auxi = (2.*pi*r0)/(a0*sqrt(gridaux1))

	ed = (ediel(1)+ediel(3))/(2.0)


	if (modk .lt. tolr) then

		v2dk = vbz*(cic/ed)*(a0*sqrt(gridaux1)/(2.*pi))*(alpha1+auxi*alpha2+alpha3*auxi**2)
	else 

		v2dk = vbz*(cic/ed)*(1./(modk*(1+(r0*modk))))
	end if

	


end function 

!potencial 3D tradicional

function vcoul(kpt1,kpt2,rlat,ngrid,tolr)

	implicit none

	double precision,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	double precision,parameter:: pi=acos(-1.)

	integer,dimension(3) :: ngrid
	double precision,dimension(3) :: kpt1,kpt2
	double precision,dimension(3,3) :: rlat
	double precision,dimension(3) :: ediel

	double precision :: modk,ed,tolr,vbz,vc
	double precision :: vcoul

	call modvec(kpt1,kpt2,modk)
	call vcell3D(rlat,vc)

	ed = 1.0

	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))

	if (modk .lt. tolr) then

		vcoul = 0.0
	else 

		vcoul = vbz*(cic/ed)*(1.0/(modk*modk))
	end if

end function 

!potencial 3D considerando ambiente dielétrico

function v3diel(kpt1,kpt2,ediel,rlat,ngrid,tolr)

	implicit none

	double precision,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	double precision,parameter:: pi=acos(-1.)

	integer,dimension(3) :: ngrid
	double precision,dimension(3) :: kpt1,kpt2
	double precision,dimension(3,3) :: rlat
	double precision :: tolr
	double precision,dimension(3) :: ediel

	double precision :: modk,ed,vbz,vc
	double precision :: v3diel

	call modvec(kpt1,kpt2,modk)
	call vcell3D(rlat,vc)

	ed = ediel(2)

	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))

	if (modk .lt. tolr) then

		v3diel = 0.0
	else 

		v3diel = vbz*(cic/ed)*(1.0/(modk*modk))
	end if



end function 

!potencial 2D truncado (DOI: 10.1103/PhysRevB.73.205119)

function v2dt(kpt1,kpt2,ngrid,rlat,tolr)

	implicit none

	double precision,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	double precision,parameter:: pi=acos(-1.)

	double precision,dimension(3) :: kpt1,kpt2,vkpt
	double precision,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid

	double precision :: vc,vbz,modk,factor
	double precision :: gpar,gz,rc
	double precision :: v2dt,tolr

	double precision :: aux1,aux2,aux3,aux4,aux5

	call vcell3D(rlat,vc)
	call modvec(kpt1,kpt2,modk)

	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))

	vkpt=kpt1-kpt2

	gz= abs(vkpt(3))

	gpar= sqrt((vkpt(1)*vkpt(1))+(vkpt(2)*vkpt(2)))

	rc = 0.5*rlat(3,3)

	!factor= 4.0*pi
	factor=1.0


	if ( (gpar .lt. tolr) .and. (gz .lt. tolr) ) then

		!v2dt = (vbz*cic)*(-2.0*pi*rc*rc)
		v2dt = (vbz*cic)*(-0.5*rc*rc)		

	else if ((gpar .lt. tolr) .and. (gz .ge. tolr)) then

		v2dt = (vbz*cic)*((factor)/(modk*modk))*(1.0-cos(gz*rc)-(gz*rc*sin(gz*rc)))

	else 
		aux1= gz/gpar
		aux2= gpar*rc
		aux3= gz*rc
		
		aux4= aux1*sin(aux3)
		aux5= cos(aux3)

		v2dt = (vbz*cic)*((factor)/(modk*modk))*(1.0+(dexp(-aux2)*(aux4-aux5)))

	end if



end function 

!potencial 0D truncado (DOI: 10.1103/PhysRevB.73.205119)

function v0dt(kpt1,kpt2,ngrid,rlat,tolr)

	implicit none

	double precision,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	double precision,parameter:: pi=acos(-1.)

	double precision,dimension(3) :: kpt1,kpt2
	double precision,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid

	double precision :: vc,modk,cr,vbz,factor
	double precision,dimension(3) :: vsize
	double precision :: v0dt,tolr

	call vcell3D(rlat,vc)
	call modvec(kpt1,kpt2,modk)

	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))

	call vecsize(rlat(1,:),vsize(1))
	call vecsize(rlat(2,:),vsize(2))
	call vecsize(rlat(3,:),vsize(3))

	cr = MIN(vsize(1),vsize(2),vsize(3))
	cr = 0.5*cr

	!factor= 4.0*pi
	factor=1.0


	if (modk .lt. tolr) then

		!v0dt = (vbz*cic)*(2.0*pi)*cr*cr
		v0dt = (vbz*cic)*(0.5)*cr*cr		
	else 

		v0dt = (vbz*cic)*((factor)/(modk*modk))*(1.0-cos(cr*modk))
	end if


end function 

!potencial 2D truncado v2 (DOI: 10.1103/PhysRevB.73.233103)


function v2dt2(kpt1,kpt2,ngrid,rlat,lc,tolr)

	implicit none

	double precision,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	double precision,parameter:: pi=acos(-1.)

	double precision,dimension(3) :: kpt1,kpt2,vkpt
	double precision,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid

	double precision :: vc,modk,vbz,qxy,lc
	double precision,dimension(3) :: vsize
	double precision :: v2dt2,factor
	double precision :: tolr

	!factor= 4.0*pi
	factor=1.0

	call vcell3D(rlat,vc)
	call modvec(kpt1,kpt2,modk)

	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))

	vkpt=kpt1-kpt2

	qxy= sqrt((vkpt(1)*vkpt(1))+(vkpt(2)*vkpt(2)))


	if (modk .lt. tolr) then 

		v2dt2 = 0.0
	else 

		v2dt2 = (vbz*cic)*((factor)/(modk*modk))*(1.0-dexp(-0.5*qxy*lc)*cos(0.5*lc*vkpt(3)))
	end if
	


end function

!potencial Rytova-Keldysh (DOI: 10.1103/PhysRevB.98.125308)

function v2drk(kpt1,kpt2,ngrid,rlat,ediel,lc,ez,w,r0,tolr)

	implicit none

	double precision,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	double precision,parameter:: pi=acos(-1.)

	double precision,dimension(3) :: kpt1,kpt2,vkpt
	double precision,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid
	double precision :: v2drk,modk
	double precision,dimension(3) :: ediel
	double precision :: vc,lc,vbz,tolr

	
	!parameters
	double precision :: w,ew
	double precision :: ez,epar
	double precision :: eta,kappa
	double precision :: et,eb,pb,pt
	double precision :: r0
	double precision :: aux1,aux2,aux3
	

	call vcell2D(rlat,vc)
	call modvec(kpt1,kpt2,modk)
	
	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))
	!vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3)))
	
	!r0= ((ediel(2)-1.0)*rlat(3,3))/(ediel(1)+ediel(3))
	!r0 = 40.0
	!lc = 6.5


	epar = ediel(2)
	et = ediel(1)
	eb = ediel(3)
	
	eta = sqrt(epar/ez)
	kappa = sqrt(epar*ez)
	
	pb = (eb-kappa)/(eb+kappa)
	pt = (et-kappa)/(et+kappa)
	

	
	if (modk .lt. tolr) then

		v2drk = 0.0
	else 
	
		aux1 = (1.0-(pb*pt*dexp(-2.0*modk*eta*lc)))*kappa
		aux2 = (1.0-(pt*dexp(-eta*modk*lc)))*(1.0-(pb*dexp(-eta*modk*lc)))
		aux3 = r0*modk*dexp(-modk*w)
	
		ew = (aux1/aux2)+aux3
		
		v2drk = (vbz*cic)*dexp(-modk*w)*(1.0/ew)*(1.0/modk)
	
	end if
	
!novos inputs : w,ez,r0

end function

!potential ohno (DOI: 10.1103/PhysRevB.89.235410)

function v2dohono(kpt1,kpt2,ngrid,rlat,ediel,w,ez,tolr)

	implicit none

	double precision,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	double precision,parameter:: pi=acos(-1.)

	double precision,dimension(3) :: kpt1,kpt2,vkpt
	double precision,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid
	double precision :: v2dohono,modk
	double precision,dimension(3) :: ediel
	double precision :: vc,lc,vbz,w,ez,tolr
	
	
	call vcell2D(rlat,vc)
	call modvec(kpt1,kpt2,modk)
	
	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))	
	!vbz= 1.0/vc !/((ngrid(1)*ngrid(2)*ngrid(3)))	
	
	if (modk .lt. tolr) then

		v2dohono = 0.0
	else 
		
		!v2dohono = (vbz*cic)*(1.0/(4.0*pi*pi))*dexp(-w*modk)*(1.0/(ez*modk))
		v2dohono = (vbz*cic)*dexp(-w*modk)*(1.0/(ez*modk))
	
	end if


!novos inputs : w

end function


