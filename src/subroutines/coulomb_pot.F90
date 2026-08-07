!potencial 2D Keldysh (Ridolfi) DOI:https://doi.org/10.1103/PhysRevB.97.205409
function v2dk(kpt1,kpt2,ediel,rlat,ngrid,lc,tolr)

	implicit none

	real,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
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
	real :: v2dk

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

	real,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	real,parameter:: pi=acos(-1.)

	integer,dimension(3) :: ngrid
	real,dimension(3) :: kpt1,kpt2
	real,dimension(3,3) :: rlat
	real,dimension(3) :: ediel

	real :: modk,ed,tolr,vbz,vc
	real :: vcoul

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

	real,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	real,parameter:: pi=acos(-1.)

	integer,dimension(3) :: ngrid
	real,dimension(3) :: kpt1,kpt2
	real,dimension(3,3) :: rlat
	real :: tolr
	real,dimension(3) :: ediel

	real :: modk,ed,vbz,vc
	real :: v3diel

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

!potencial 2D tradicional

function v2d(kpt1,kpt2,rlat,ngrid,tolr)

	implicit none

	real,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	real,parameter:: pi=acos(-1.)

	integer,dimension(3) :: ngrid
	real,dimension(3) :: kpt1,kpt2
	real,dimension(3,3) :: rlat
	real,dimension(3) :: ediel

	real :: modk,ed,tolr,vbz,vc
	real :: v2d

	call modvec(kpt1,kpt2,modk)
	call vcell2D(rlat,vc)
	

	ed = 1.0

	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))

	if (modk .lt. tolr) then

		v2d = 0.0
	else 

		v2d = vbz*(cic/ed)*(1.0/(modk))
	end if

end function 


!potencial 2D considerando ambiente dielétrico

function v2diel(kpt1,kpt2,ediel,rlat,ngrid,tolr)

	implicit none

	real,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	real,parameter:: pi=acos(-1.)

	integer,dimension(3) :: ngrid
	real,dimension(3) :: kpt1,kpt2
	real,dimension(3,3) :: rlat
	real :: tolr
	real,dimension(3) :: ediel

	real :: modk,ed,vbz,vc
	real :: v2diel

	call modvec(kpt1,kpt2,modk)
	call vcell2D(rlat,vc)
	

	ed = ediel(2)

	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))

	if (modk .lt. tolr) then

		v2diel = 0.0
	else 

		v2diel = vbz*(cic/ed)*(1.0/(modk))
	end if



end function 


!potencial 2D truncado (DOI: 10.1103/PhysRevB.73.205119)

function v2dt(kpt1,kpt2,ngrid,rlat,tolr)

	implicit none

	real,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	real,parameter:: pi=acos(-1.)

	real,dimension(3) :: kpt1,kpt2,vkpt
	real,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid

	real :: vc,vbz,modk,factor
	real :: gpar,gz,rc
	real :: v2dt,tolr

	real :: aux1,aux2,aux3,aux4,aux5

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

		v2dt = (vbz*cic)*((factor)/(modk*modk))*(1.0+(exp(-aux2)*(aux4-aux5)))

	end if



end function 


! Slab-truncated 2D potential with the q=0 value averaged over the
! reciprocal-space sampling cell. Away from q=0 this is identical to V2DT.
function v2dtavg(kpt1,kpt2,ngrid,rlat,tolr)

	implicit none

	real,parameter :: cic=-(0.0904756)*10**3
	real,dimension(3) :: kpt1,kpt2,vkpt
	real,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid
	real :: tolr,modk,v2dtavg
	real :: vc,vbz,factor,gpar,gz,rc
	real :: aux1,aux2,aux3,aux4,aux5
	real(kind=8) :: v2dtavg_cell

	logical,save :: cache_valid=.false.
	integer,dimension(3),save :: cached_ngrid=(/0,0,0/)
	real,dimension(3,3),save :: cached_rlat=0.0
	real,save :: cached_value=0.0

	call vcell3D(rlat,vc)
	call modvec(kpt1,kpt2,modk)

	vbz=1.0/((ngrid(1)*ngrid(2)*ngrid(3))*vc)
	vkpt=kpt1-kpt2
	gz=abs(vkpt(3))
	gpar=sqrt(vkpt(1)*vkpt(1)+vkpt(2)*vkpt(2))
	rc=0.5*rlat(3,3)
	factor=1.0

	if ((gpar .lt. tolr) .and. (gz .lt. tolr)) then
		! The BSE construction is OpenMP parallel. Cache the mesh-dependent
		! average so the reciprocal-cell quadrature is evaluated only once.
		!$omp critical (v2dtavg_cache)
		if (.not. cache_valid) then
			cached_value=real(v2dtavg_cell(ngrid,rlat))
			cached_ngrid=ngrid
			cached_rlat=rlat
			cache_valid=.true.
		else if (any(cached_ngrid .ne. ngrid) .or. any(cached_rlat .ne. rlat)) then
			cached_value=real(v2dtavg_cell(ngrid,rlat))
			cached_ngrid=ngrid
			cached_rlat=rlat
		end if
		v2dtavg=cached_value
		!$omp end critical (v2dtavg_cache)

	else if ((gpar .lt. tolr) .and. (gz .ge. tolr)) then
		v2dtavg=(vbz*cic)*(factor/(modk*modk)) &
		         *(1.0-cos(gz*rc)-(gz*rc*sin(gz*rc)))

	else
		aux1=gz/gpar
		aux2=gpar*rc
		aux3=gz*rc
		aux4=aux1*sin(aux3)
		aux5=cos(aux3)
		v2dtavg=(vbz*cic)*(factor/(modk*modk)) &
		         *(1.0+exp(-aux2)*(aux4-aux5))
	end if

end function v2dtavg


! Average [1-exp(-q*Rc)]/q**2 over the parallelogram generated by
! b1/ngrid(1) and b2/ngrid(2). Four Duffy-mapped triangles remove the
! integrable 1/q singularity from the numerical quadrature.
real(kind=8) function v2dtavg_cell(ngrid,rlat)

	implicit none

	integer,dimension(3) :: ngrid
	real,dimension(3,3) :: rlat

	integer,parameter :: nq=16
	real(kind=8),parameter,dimension(nq) :: xg=(/ &
		-0.989400934991650d0,-0.944575023073233d0,-0.865631202387832d0, &
		-0.755404408355003d0,-0.617876244402644d0,-0.458016777657227d0, &
		-0.281603550779259d0,-0.095012509837637d0, 0.095012509837637d0, &
		 0.281603550779259d0, 0.458016777657227d0, 0.617876244402644d0, &
		 0.755404408355003d0, 0.865631202387832d0, 0.944575023073233d0, &
		 0.989400934991650d0 /)
	real(kind=8),parameter,dimension(nq) :: wg=(/ &
		0.027152459411754d0,0.062253523938648d0,0.095158511682493d0, &
		0.124628971255534d0,0.149595988816577d0,0.169156519395003d0, &
		0.182603415044924d0,0.189450610455069d0,0.189450610455069d0, &
		0.182603415044924d0,0.169156519395003d0,0.149595988816577d0, &
		0.124628971255534d0,0.095158511682493d0,0.062253523938648d0, &
		0.027152459411754d0 /)

	real(kind=8),parameter :: pi=acos(-1.0d0)
	real(kind=8),parameter :: cic=-(0.0904756d0)*10.0d0**3
	real(kind=8),dimension(3) :: a1,a2,a3,cross23,cross31
	real(kind=8),dimension(3) :: b1,b2,dq1,dq2
	real(kind=8),dimension(3,4) :: vertex
	real(kind=8),dimension(3) :: va,vb,crossab,direction,qvec
	real(kind=8) :: volume,cell_area,rc
	real(kind=8) :: integral,jacobian,s,t,ws,wt
	real(kind=8) :: qnorm,x,kernel,avg_kernel,grid_size
	integer :: itri,inext,i,j

	a1=dble(rlat(1,:))
	a2=dble(rlat(2,:))
	a3=dble(rlat(3,:))

	cross23(1)=a2(2)*a3(3)-a2(3)*a3(2)
	cross23(2)=a2(3)*a3(1)-a2(1)*a3(3)
	cross23(3)=a2(1)*a3(2)-a2(2)*a3(1)
	volume=a1(1)*cross23(1)+a1(2)*cross23(2)+a1(3)*cross23(3)

	cross31(1)=a3(2)*a1(3)-a3(3)*a1(2)
	cross31(2)=a3(3)*a1(1)-a3(1)*a1(3)
	cross31(3)=a3(1)*a1(2)-a3(2)*a1(1)

	b1=(2.0d0*pi/volume)*cross23
	b2=(2.0d0*pi/volume)*cross31
	dq1=b1/dble(ngrid(1))
	dq2=b2/dble(ngrid(2))

	vertex(:,1)= 0.5d0*(dq1+dq2)
	vertex(:,2)= 0.5d0*(-dq1+dq2)
	vertex(:,3)=-0.5d0*(dq1+dq2)
	vertex(:,4)= 0.5d0*(dq1-dq2)

	crossab(1)=dq1(2)*dq2(3)-dq1(3)*dq2(2)
	crossab(2)=dq1(3)*dq2(1)-dq1(1)*dq2(3)
	crossab(3)=dq1(1)*dq2(2)-dq1(2)*dq2(1)
	cell_area=sqrt(sum(crossab*crossab))

	rc=0.5d0*dble(rlat(3,3))
	integral=0.0d0

	do itri=1,4
		inext=mod(itri,4)+1
		va=vertex(:,itri)
		vb=vertex(:,inext)

		crossab(1)=va(2)*vb(3)-va(3)*vb(2)
		crossab(2)=va(3)*vb(1)-va(1)*vb(3)
		crossab(3)=va(1)*vb(2)-va(2)*vb(1)
		jacobian=sqrt(sum(crossab*crossab))

		do i=1,nq
			s=0.5d0*(xg(i)+1.0d0)
			ws=0.5d0*wg(i)
			do j=1,nq
				t=0.5d0*(xg(j)+1.0d0)
				wt=0.5d0*wg(j)
				direction=(1.0d0-t)*va+t*vb
				qvec=s*direction
				qnorm=sqrt(sum(qvec*qvec))
				x=rc*qnorm

				if (abs(x) .lt. 1.0d-5) then
					kernel=rc/qnorm-0.5d0*rc*rc &
					       +(rc**3)*qnorm/6.0d0-(rc**4)*qnorm*qnorm/24.0d0
				else
					kernel=(1.0d0-exp(-x))/(qnorm*qnorm)
				end if

				integral=integral+ws*wt*jacobian*s*kernel
			end do
		end do
	end do

	avg_kernel=integral/cell_area
	grid_size=dble(ngrid(1))*dble(ngrid(2))*dble(ngrid(3))
	v2dtavg_cell=cic*avg_kernel/(grid_size*abs(volume))

end function v2dtavg_cell


!potencial 0D truncado (DOI: 10.1103/PhysRevB.73.205119)

function v0dt(kpt1,kpt2,ngrid,rlat,tolr)

	implicit none

	real,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	real,parameter:: pi=acos(-1.)

	real,dimension(3) :: kpt1,kpt2
	real,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid

	real :: vc,modk,cr,vbz,factor
	real,dimension(3) :: vsize
	real :: v0dt,tolr

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

	real,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	real,parameter:: pi=acos(-1.)

	real,dimension(3) :: kpt1,kpt2,vkpt
	real,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid

	real :: vc,modk,vbz,qxy,lc
	real,dimension(3) :: vsize
	real :: v2dt2,factor
	real :: tolr

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

		v2dt2 = (vbz*cic)*((factor)/(modk*modk))*(1.0-exp(-0.5*qxy*lc)*cos(0.5*lc*vkpt(3)))
	end if
	


end function

!potencial Rytova-Keldysh (DOI: 10.1103/PhysRevB.98.125308)

function v2drk(kpt1,kpt2,ngrid,rlat,ediel,lc,ez,w,r0,tolr)

	implicit none

	real,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	real,parameter:: pi=acos(-1.)

	real,dimension(3) :: kpt1,kpt2,vkpt
	real,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid
	real :: v2drk,modk
	real,dimension(3) :: ediel
	real :: vc,lc,vbz,tolr

	
	!parameters
	real :: w,ew
	real :: ez,epar
	real :: eta,kappa
	real :: et,eb,pb,pt
	real :: r0
	real :: aux1,aux2,aux3
	

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
	
		aux1 = (1.0-(pb*pt*exp(-2.0*modk*eta*lc)))*kappa
		aux2 = (1.0-(pt*exp(-eta*modk*lc)))*(1.0-(pb*exp(-eta*modk*lc)))
		aux3 = r0*modk*exp(-modk*w)
	
		ew = (aux1/aux2)+aux3
		
		v2drk = (vbz*cic)*exp(-modk*w)*(1.0/ew)*(1.0/modk)
	
	end if
	
!novos inputs : w,ez,r0

end function

!potential ohno (DOI: 10.1103/PhysRevB.89.235410)

function v2dohono(kpt1,kpt2,ngrid,rlat,ediel,w,ez,tolr)

	implicit none

	real,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	real,parameter:: pi=acos(-1.)

	real,dimension(3) :: kpt1,kpt2,vkpt
	real,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid
	real :: v2dohono,modk
	real,dimension(3) :: ediel
	real :: vc,lc,vbz,w,ez,tolr
	
	
	call vcell2D(rlat,vc)
	call modvec(kpt1,kpt2,modk)
	
	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))	
	!vbz= 1.0/vc !/((ngrid(1)*ngrid(2)*ngrid(3)))	
	
	if (modk .lt. tolr) then

		v2dohono = 0.0
	else 
		
		!v2dohono = (vbz*cic)*(1.0/(4.0*pi*pi))*exp(-w*modk)*(1.0/(ez*modk))
		v2dohono = (vbz*cic)*exp(-w*modk)*(1.0/(ez*modk))
	
	end if
	
	!novos inputs : w

end function

function v1d(kpt1,kpt2,ngrid,rlat,lc,tolr)
	USE Exponential_Integral
	implicit none
	
	real,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	real,parameter:: pi=acos(-1.)
	real :: aux1,aux2,res
	real,dimension(3) :: kpt1,kpt2
	real,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid
	real :: lc,tolr,modk
	
	double precision :: res1,arg
	
	real :: v1d,vbz
	
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
	
	vbz = 1./((ngrid(1)*ngrid(2)*ngrid(3))*(rlat(1,1)))
	
	
	call modvec(kpt1,kpt2,modk)
	
	if (modk .lt. tolr) then
	
		v1d = 0.0
	
	else
	
		arg = (lc*modk)*(lc*modk)
		!call CALCEI(-arg,res1,1)
		arg2 = CMPLX(arg, 0.0_dp, KIND=dp)  ! argumento complexo do ponto da funcao
		CALL cexqad(arg2, n, kode, tol, cw, kerr)
		
		res1 = real(cw)
				
		aux1 = cic/(2.0*pi)
		
		v1d = -aux1*exp(arg)*res1*vbz
	
	end if
	
	

end function v1d


function v1diel(kpt1,kpt2,ngrid,rlat,lc,ediel,tolr)

	USE Exponential_Integral
	implicit none
	
	real,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	real,parameter:: pi=acos(-1.)
	real :: aux1,aux2,res
	real,dimension(3) :: kpt1,kpt2
	real,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid
	real :: lc,tolr,modk,vbz
	real,dimension(3) :: ediel
	
	double precision :: res1,arg
	
	real :: v1diel
	
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
	
	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(rlat(1,1)))
	
	res1 = 0.0
	
	
	call modvec(kpt1,kpt2,modk)
	
	if (modk .lt. tolr) then
	
		v1diel = 0.0
	
	else
	
		arg = (lc*modk)*(lc*modk)
		!call CALCEI(-arg,res1,1)
		arg2 = CMPLX(arg, 0.0_dp, KIND=dp)  ! argumento complexo do ponto da funcao
		CALL cexqad(arg2, n, kode, tol, cw, kerr)
		
		res1 = real(cw)
		
		aux1 = cic/(2.0*pi*ediel(2))
		
		v1diel = -aux1*exp(arg)*res1*vbz
	
	end if
	
	

end function v1diel

!1D truncated potential  (DOI: 10.1103/PhysRevB.73.205119)

function v1dt(kpt1,kpt2,ngrid,rlat,tolr,lc)

	implicit none
	real,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	real,parameter:: pi=acos(-1.)
	
	real,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid
	real,dimension(3) :: kpt1,kpt2,vkpt
	real :: tolr,lc
	real :: v1dt,modk,vc,vbz
	
	real :: c1,c2,qyz,res
	real :: j0,j1,k0,k1

	call vcell3D(rlat,vc)
	call modvec(kpt1,kpt2,modk)
	
	!lc = (1.0+sqrt(2.0))*rlat(1,1)
	
	vbz= 1./((ngrid(1)*ngrid(2)*ngrid(3))*(vc))
	vkpt=kpt1-kpt2
	
	qyz=sqrt( (vkpt(2)*vkpt(2)) +(vkpt(3)*vkpt(3)) )
	
	
	if ( (abs(vkpt(1)) .lt. tolr) .and. (qyz .ge. tolr) ) then
	
		call integralJ0pot(lc,qyz,res)
	
		v1dt = (vbz*cic)*res
	
	else if  ( (abs(vkpt(1)) .lt. tolr) .and. (qyz .lt. tolr) ) then
	

		v1dt = (vbz*cic)*((lc*lc)/16.0)*(2.0*log(lc)-1.0) 
	
	else

		call caljy0( qyz*lc, j0, 0 )
		call caljy1( qyz*lc, j1, 0 )
		call calck0( abs(vkpt(1))*lc, k0, 1 )
		call calck1( abs(vkpt(1))*lc, k1, 1 )
		
		c1 = qyz*lc*j1*k0
		
		c2 = abs(vkpt(1))*lc*j0*k1
		
		v1dt =	(vbz*cic)*(1.0+c1+c2)
	
	end if


end function

