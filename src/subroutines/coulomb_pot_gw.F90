function wk(sysdim,ngrid,rlat)

	implicit none
	character(len=5) :: sysdim
	integer,dimension(3) :: ngrid
	real,parameter:: pi=acos(-1.)
	real,dimension(3,3) :: rlat
	real :: nk,wk
	real :: vbz
	real :: aux
	
	nk = ngrid(1)*ngrid(2)*ngrid(3)
	
	select case (sysdim)
	
	case("3D")
	
	 call vcell3D(rlat,vbz)
	 aux = (2.0*pi)*(2.0*pi)*(2.0*pi)
	 wk = 1.0/(nk*vbz*aux)	
	
	case("2D")
	
	 call vcell2D(rlat,vbz)
	 aux = (2.0*pi)*(2.0*pi)
	 wk = 1.0/(nk*vbz*aux)	
	
	case("1D")
	
	 vbz = rlat(1,1)
	 aux = 2.0*pi
	 wk = 1.0/(nk*vbz*aux)	
	
	case default
	
		write(*,*) "Wrong value for system dimension"
		STOP
	
	end select


end function


!potencial 2D Keldysh (Ridolfi) DOI:https://doi.org/10.1103/PhysRevB.97.205409
function v2dkgw(kpt1,kpt2,ediel,rlat,ngrid,lc,tolr)

	implicit none

	real,parameter :: cic= (0.0904756)*10**(3) !constante da interação coulombiana (e^2/2 e0)
	real,parameter:: pi=acos(-1.)
	real,parameter :: alpha1 = 1.76
	real,parameter :: alpha2 = 1.0
	real,parameter :: alpha3 = 0d0

	integer,dimension(3) :: ngrid
	real,dimension(3) :: kpt1,kpt2
	real,dimension(3,3) :: rlat
	real :: tolr
	real,dimension(3) :: ediel
	real :: lc
	real :: gridaux1,ed


	real :: a0,modk,vc
	real :: r0,vbz,auxi
	real :: v2dkgw

	call alat2D(rlat,a0)
	call modvec(kpt1,kpt2,modk)
	!call vcell2D(rlat,vc)

	r0= ((ediel(2)-1.0)*lc)/(ediel(1)+ediel(3))

	!vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))
	vbz = 1.0
	gridaux1 = dble(ngrid(1)*ngrid(2))

	auxi = (2.*pi*r0)/(a0*sqrt(gridaux1))

	ed = (ediel(1)+ediel(3))/(2.0)


	if (modk .lt. tolr) then

		v2dkgw = vbz*(cic/ed)*(a0*sqrt(gridaux1)/(2.*pi))*(alpha1+auxi*alpha2+alpha3*auxi**2)
	else 

		v2dkgw = vbz*(cic/ed)*(1./(modk*(1+(r0*modk))))
	end if

	


end function 

!potencial 3D tradicional

function vcoulgw(kpt1,kpt2,rlat,ngrid,tolr)

	implicit none

	real,parameter :: cic= (0.0904756)*10**(3) !constante da interação coulombiana (e^2/2 e0)
	real,parameter:: pi=acos(-1.)

	integer,dimension(3) :: ngrid
	real,dimension(3) :: kpt1,kpt2
	real,dimension(3,3) :: rlat
	real,dimension(3) :: ediel

	real :: modk,ed,tolr,vbz,vc
	real :: vcoulgw

	call modvec(kpt1,kpt2,modk)
	!call vcell3D(rlat,vc)

	ed = 1.0

	!vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))
	 vbz = 1.0
	if (modk .lt. tolr) then

		vcoulgw = 0.0
	else 

		vcoulgw = vbz*(cic/ed)*(1.0/(modk*modk))
	end if

end function 

!potencial 3D considerando ambiente dielétrico

function v3dielgw(kpt1,kpt2,ediel,rlat,ngrid,tolr)

	implicit none

	real,parameter :: cic= (0.0904756)*10**(3) !constante da interação coulombiana (e^2/2 e0)
	real,parameter:: pi=acos(-1.)

	integer,dimension(3) :: ngrid
	real,dimension(3) :: kpt1,kpt2
	real,dimension(3,3) :: rlat
	real :: tolr
	real,dimension(3) :: ediel

	real :: modk,ed,vbz,vc
	real :: v3dielgw

	call modvec(kpt1,kpt2,modk)
	!call vcell3D(rlat,vc)

	ed = ediel(2)

	!vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))
	vbz = 1.0
	if (modk .lt. tolr) then

		v3dielgw = 0.0
	else 

		v3dielgw = vbz*(cic/ed)*(1.0/(modk*modk))
	end if



end function 

!potencial 2D tradicional

function v2dgw(kpt1,kpt2,rlat,ngrid,tolr)

	implicit none

	real,parameter :: cic= (0.0904756)*10**(3) !constante da interação coulombiana (e^2/2 e0)
	real,parameter:: pi=acos(-1.)

	integer,dimension(3) :: ngrid
	real,dimension(3) :: kpt1,kpt2
	real,dimension(3,3) :: rlat
	real,dimension(3) :: ediel

	real :: modk,ed,tolr,vbz,vc
	real :: v2dgw

	call modvec(kpt1,kpt2,modk)
	!call vcell2D(rlat,vc)

	ed = 1.0

	!vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))
	vbz = 1.0
	if (modk .lt. tolr) then

		v2dgw = 0.0
	else 

		v2dgw = vbz*(cic/ed)*(1.0/(modk))
	end if

end function 


!potencial 2D considerando ambiente dielétrico

function v2dielgw(kpt1,kpt2,ediel,rlat,ngrid,tolr)

	implicit none

	real,parameter :: cic= (0.0904756)*10**(3) !constante da interação coulombiana (e^2/2 e0)
	real,parameter:: pi=acos(-1.)

	integer,dimension(3) :: ngrid
	real,dimension(3) :: kpt1,kpt2
	real,dimension(3,3) :: rlat
	real :: tolr
	real,dimension(3) :: ediel

	real :: modk,ed,vbz,vc
	real :: v2dielgw

	call modvec(kpt1,kpt2,modk)
	!call vcell2D(rlat,vc)

	ed = ediel(2)

	!vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))
	vbz = 1.0
	if (modk .lt. tolr) then

		v2dielgw = 0.0
	else 

		v2dielgw = vbz*(cic/ed)*(1.0/(modk))
	end if



end function 


!potencial 2D truncado (DOI: 10.1103/PhysRevB.73.205119)

function v2dtgw(kpt1,kpt2,ngrid,rlat,tolr)

	implicit none

	real,parameter :: cic= (0.0904756)*10**(3) !constante da interação coulombiana (e^2/2 e0)
	real,parameter:: pi=acos(-1.)

	real,dimension(3) :: kpt1,kpt2,vkpt
	real,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid

	real :: vc,vbz,modk,factor
	real :: gpar,gz,rc
	real :: v2dtgw,tolr

	real :: aux1,aux2,aux3,aux4,aux5

	!call vcell3D(rlat,vc)
	call modvec(kpt1,kpt2,modk)

	!vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))
	vbz = 1.0

	vkpt=kpt1-kpt2

	gz= abs(vkpt(3))

	gpar= sqrt((vkpt(1)*vkpt(1))+(vkpt(2)*vkpt(2)))

	rc = 0.5*rlat(3,3)

	!factor= 4.0*pi
	factor=1.0


	if ( (gpar .lt. tolr) .and. (gz .lt. tolr) ) then

		!v2dt = (vbz*cic)*(-2.0*pi*rc*rc)
		v2dtgw = (vbz*cic)*(-0.5*rc*rc)		

	else if ((gpar .lt. tolr) .and. (gz .ge. tolr)) then

		v2dtgw = (vbz*cic)*((factor)/(modk*modk))*(1.0-cos(gz*rc)-(gz*rc*sin(gz*rc)))

	else 
		aux1= gz/gpar
		aux2= gpar*rc
		aux3= gz*rc
		
		aux4= aux1*sin(aux3)
		aux5= cos(aux3)

		v2dtgw = (vbz*cic)*((factor)/(modk*modk))*(1.0+(exp(-aux2)*(aux4-aux5)))

	end if



end function 

!potencial 0D truncado (DOI: 10.1103/PhysRevB.73.205119)

function v0dtgw(kpt1,kpt2,ngrid,rlat,tolr)

	implicit none

	real,parameter :: cic= (0.0904756)*10**(3) !constante da interação coulombiana (e^2/2 e0)
	real,parameter:: pi=acos(-1.)

	real,dimension(3) :: kpt1,kpt2
	real,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid

	real :: vc,modk,cr,vbz,factor
	real,dimension(3) :: vsize
	real :: v0dtgw,tolr

	!call vcell3D(rlat,vc)
	call modvec(kpt1,kpt2,modk)

	!vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))
	vbz = 1.0
	
	call vecsize(rlat(1,:),vsize(1))
	call vecsize(rlat(2,:),vsize(2))
	call vecsize(rlat(3,:),vsize(3))

	cr = MIN(vsize(1),vsize(2),vsize(3))
	cr = 0.5*cr

	!factor= 4.0*pi
	factor=1.0


	if (modk .lt. tolr) then

		!v0dt = (vbz*cic)*(2.0*pi)*cr*cr
		v0dtgw = (vbz*cic)*(0.5)*cr*cr		
	else 

		v0dtgw = (vbz*cic)*((factor)/(modk*modk))*(1.0-cos(cr*modk))
	end if


end function 

!potencial 2D truncado v2 (DOI: 10.1103/PhysRevB.73.233103)


function v2dt2gw(kpt1,kpt2,ngrid,rlat,lc,tolr)

	implicit none

	real,parameter :: cic= (0.0904756)*10**(3) !constante da interação coulombiana (e^2/2 e0)
	real,parameter:: pi=acos(-1.)

	real,dimension(3) :: kpt1,kpt2,vkpt
	real,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid

	real :: vc,modk,vbz,qxy,lc
	real,dimension(3) :: vsize
	real :: v2dt2gw,factor
	real :: tolr

	!factor= 4.0*pi
	factor=1.0

	!call vcell3D(rlat,vc)
	call modvec(kpt1,kpt2,modk)

	!vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))
	vbz = 1.0
	
	vkpt=kpt1-kpt2

	qxy= sqrt((vkpt(1)*vkpt(1))+(vkpt(2)*vkpt(2)))


	if (modk .lt. tolr) then 

		v2dt2gw = 0.0
	else 

		v2dt2gw = (vbz*cic)*((factor)/(modk*modk))*(1.0-exp(-0.5*qxy*lc)*cos(0.5*lc*vkpt(3)))
	end if
	


end function

!potencial Rytova-Keldysh (DOI: 10.1103/PhysRevB.98.125308)

function v2drkgw(kpt1,kpt2,ngrid,rlat,ediel,lc,ez,w,r0,tolr)

	implicit none

	real,parameter :: cic= (0.0904756)*10**(3) !constante da interação coulombiana (e^2/2 e0)
	real,parameter:: pi=acos(-1.)

	real,dimension(3) :: kpt1,kpt2,vkpt
	real,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid
	real :: v2drkgw,modk
	real,dimension(3) :: ediel
	real :: vc,lc,vbz,tolr

	
	!parameters
	real :: w,ew
	real :: ez,epar
	real :: eta,kappa
	real :: et,eb,pb,pt
	real :: r0
	real :: aux1,aux2,aux3
	

	!call vcell2D(rlat,vc)
	call modvec(kpt1,kpt2,modk)
	
	!vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))
	vbz = 1.0
	
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
	
		aux1 = (1.0-(pb*pt*exp(-2.0*modk*eta*lc)))*kappa
		aux2 = (1.0-(pt*exp(-eta*modk*lc)))*(1.0-(pb*exp(-eta*modk*lc)))
		aux3 = r0*modk*exp(-modk*w)
	
		ew = (aux1/aux2)+aux3
		
		v2drk = (vbz*cic)*exp(-modk*w)*(1.0/ew)*(1.0/modk)
	
	end if
	
!novos inputs : w,ez,r0

end function

!potential ohno (DOI: 10.1103/PhysRevB.89.235410)

function v2dohonogw(kpt1,kpt2,ngrid,rlat,ediel,w,ez,tolr)

	implicit none

	real,parameter :: cic= (0.0904756)*10**(3) !constante da interação coulombiana (e^2/2 e0)
	real,parameter:: pi=acos(-1.)

	real,dimension(3) :: kpt1,kpt2,vkpt
	real,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid
	real :: v2dohonogw,modk
	real,dimension(3) :: ediel
	real :: vc,lc,vbz,w,ez,tolr
	
	
	!call vcell2D(rlat,vc)
	call modvec(kpt1,kpt2,modk)
	
	!vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))	
	vbz = 1.0	
	
	if (modk .lt. tolr) then

		v2dohonogw = 0.0
	else 
		
		!v2dohono = (vbz*cic)*(1.0/(4.0*pi*pi))*exp(-w*modk)*(1.0/(ez*modk))
		v2dohonogw = (vbz*cic)*exp(-w*modk)*(1.0/(ez*modk))
	
	end if
	
	!novos inputs : w

end function

function v1dgw(kpt1,kpt2,ngrid,rlat,lc,tolr)
	USE Exponential_Integral
	implicit none
	
	real,parameter :: cic= (0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	real,parameter:: pi=acos(-1.)
	real :: aux1,aux2,res
	real,dimension(3) :: kpt1,kpt2
	real,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid
	real :: lc,tolr,modk
	
	double precision :: res1,arg
	
	real :: v1dgw
	
	!exponential integral variables
	INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)
  	COMPLEX(dp) :: z, cw
  	REAL(dp)    :: tol
  	INTEGER     :: n, kode, kerr
  	COMPLEX (dp)   :: arg2
  	
    	z = CMPLX(1.0_dp, 1.0_dp, KIND=dp)  ! argumento complexo do ponto da funcao
  	n = 1                               ! oprdem da funcao que oce vai avaliar
  	kode = 1                            ! KODE diz se tu vai querer a funcao no ponto ou exp(-z) * funcao no ponto
  	tol = 1.0e-12_dp                    ! Aqui eh a tolerancia para as integrais		

	res1 = 0.0	
	
	call modvec(kpt1,kpt2,modk)
	
	if (modk .lt. tolr) then
	
		v1dgw = 0.0
	
	else
	
		arg = (lc*modk)*(lc*modk)
		!call CALCEI(-arg,res1,1)
		arg2 = CMPLX(arg, 0.0_dp, KIND=dp)  ! argumento complexo do ponto da funcao
		CALL cexqad(arg2, n, kode, tol, cw, kerr)
		
		res1 = real(cw)
		
		aux1 = cic/(2.0*pi)
		
		v1dgw = -aux1*exp(arg)*res1
	
	end if
	
	

end function v1dgw


function v1dielgw(kpt1,kpt2,ngrid,rlat,lc,ediel,tolr)
	USE Exponential_Integral
	implicit none
	
	real,parameter :: cic= (0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	real,parameter:: pi=acos(-1.)
	real :: aux1,aux2,res
	real,dimension(3) :: kpt1,kpt2
	real,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid
	real :: lc,tolr,modk
	real,dimension(3) :: ediel
	
	double precision :: res1,arg
	
	real :: v1dielgw
	
	!exponential integral variables
	INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)
  	COMPLEX(dp) :: z, cw
  	REAL(dp)    :: tol
  	INTEGER     :: n, kode, kerr
  	COMPLEX (dp)   :: arg2
  	
    	z = CMPLX(1.0_dp, 1.0_dp, KIND=dp)  ! argumento complexo do ponto da funcao
  	n = 1                               ! oprdem da funcao que oce vai avaliar
  	kode = 1                            ! KODE diz se tu vai querer a funcao no ponto ou exp(-z) * funcao no ponto
  	tol = 1.0e-12_dp                    ! Aqui eh a tolerancia para as integrais		
	
	res1 = 0.0
	
	call modvec(kpt1,kpt2,modk)
	
	if (modk .lt. tolr) then
	
		v1dielgw = 0.0
	
	else
	
		arg = (lc*modk)*(lc*modk)
		!call CALCEI(-arg,res1,1)
		arg2 = CMPLX(arg, 0.0_dp, KIND=dp)  ! argumento complexo do ponto da funcao
		CALL cexqad(arg2, n, kode, tol, cw, kerr)
		
		res1 = real(cw)
		
		aux1 = cic/(2.0*pi*ediel(2))
		
		v1dielgw = -aux1*exp(arg)*res1
	
	end if
	
	

end function v1dielgw


!1D truncated potential  (DOI: 10.1103/PhysRevB.73.205119)

function v1dtgw(kpt1,kpt2,ngrid,rlat,tolr,lc)

	implicit none
	real,parameter :: cic= (0.0904756)*10**(3) !constante da interação coulombiana (e^2/2 e0)
	real,parameter:: pi=acos(-1.)
	
	real,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid
	real,dimension(3) :: kpt1,kpt2,vkpt
	real :: tolr,lc
	real :: v1dtgw,modk,vc,vbz
	
	real :: c1,c2,qyz,res
	real :: j0,j1,k0,k1

	!call vcell3D(rlat,vc)
	call modvec(kpt1,kpt2,modk)
	
	!lc = (1.0+sqrt(2.0))*rlat(1,1)
	
	!vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))
	vbz = 1.0
	vkpt=kpt1-kpt2
	
	qyz=sqrt( (vkpt(2)*vkpt(2)) +(vkpt(3)*vkpt(3)) )
	
	
	if ( (abs(vkpt(1)) .lt. tolr) .and. (qyz .ge. tolr) ) then
	
		call integralJ0pot(lc,qyz,res)
	
		v1dtgw = (vbz*cic)*res
	
	else if  ( (abs(vkpt(1)) .lt. tolr) .and. (qyz .lt. tolr) ) then
	

		v1dtgw = (vbz*cic)*((lc*lc)/16.0)*(2.0*log(lc)-1.0) 
	
	else

		call caljy0( qyz*lc, j0, 0 )
		call caljy1( qyz*lc, j1, 0 )
		call calck0( abs(vkpt(1))*lc, k0, 1 )
		call calck1( abs(vkpt(1))*lc, k1, 1 )
		
		c1 = qyz*lc*j1*k0
		
		c2 = abs(vkpt(1))*lc*j0*k1
		
		v1dtgw =	(vbz*cic)*(1.0+c1+c2)
	
	end if


end function



