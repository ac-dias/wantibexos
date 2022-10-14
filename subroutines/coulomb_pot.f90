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

!1D truncated potential  (DOI: 10.1103/PhysRevB.73.205119)

function v1dt(kpt1,kpt2,ngrid,rlat,tolr,lc)

	implicit none
	double precision,parameter :: cic= -(0.0904756)*10**(3) !constante da interação coulombiana (-e^2/2 e0)
	double precision,parameter:: pi=acos(-1.)
	
	double precision,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid
	double precision,dimension(3) :: kpt1,kpt2,vkpt
	double precision :: tolr,lc
	double precision :: v1dt,modk,vc,vbz
	
	double precision :: c1,c2,qyz,res
	double precision :: j0,j1,k0,k1

	call vcell3D(rlat,vc)
	call modvec(kpt1,kpt2,modk)
	
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


!#############Bessel Functions######################

subroutine integralJ0pot(rmax,gperp,res)

	implicit none
	double precision :: rmax, h, res
	integer :: i
	integer,parameter :: n=100
	double precision :: gperp
	double precision :: j0
	double precision,dimension(n) :: r , f
	double precision,parameter :: eta=1E-04
	
	do i=1,n
	
	 r(i)=((i-1d0)/(n-1d0))*rmax+eta
	 
	 call caljy0(r(i)*gperp, j0, 0 )
	   
	 f(i)= r(i)*j0*log(r(i))
	 
	
	end do

	h = r(2)-r(1)
	res = 0.0

	do i=1,n-1
	
		res= res + (h/2.0)*(f(i)+f(i+1))
	
	end do


end subroutine



subroutine calck0 ( arg, result, jint )

!*****************************************************************************80
!
!! CALCK0 computes various K0 Bessel functions.
!
!  Discussion:
!
!    This routine computes modified Bessel functions of the second kind
!    and order zero, K0(X) and EXP(X)*K0(X), for real
!    arguments X.
!
!    The main computation evaluates slightly modified forms of near
!    minimax rational approximations generated by Russon and Blair,
!    Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
!    1969.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  0 < ARG is
!    always required.  If JINT = 1, then the argument must also be
!    less than XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = K0(x);
!    2, RESULT = exp(x) * K0(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    1, K0(x);
!    2, exp(x) * K0(x);
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) jint
  real ( kind = 8 ) arg
  real ( kind = 8 ) f(4)
  real ( kind = 8 ) g(3)
  real ( kind = 8 ) p(6)
  real ( kind = 8 ) pp(10)
  real ( kind = 8 ) q(2)
  real ( kind = 8 ) qq(10)
  real ( kind = 8 ) result
  real ( kind = 8 ) sumf
  real ( kind = 8 ) sumg
  real ( kind = 8 ) sump
  real ( kind = 8 ) sumq
  real ( kind = 8 ) temp
  real ( kind = 8 ) x
  real ( kind = 8 ) xinf
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xsmall
  real ( kind = 8 ) xx
!
!  Machine-dependent constants
!
  data xsmall /1.11d-16/
  data xinf /1.79d+308/
  data xmax /705.342d0/
!
!  Coefficients for XSMALL <= ARG <= 1.0
!
  data   p/ 5.8599221412826100000d-04, 1.3166052564989571850d-01, &
            1.1999463724910714109d+01, 4.6850901201934832188d+02, &
            5.9169059852270512312d+03, 2.4708152720399552679d+03/
  data   q/-2.4994418972832303646d+02, 2.1312714303849120380d+04/
  data   f/-1.6414452837299064100d+00,-2.9601657892958843866d+02, &
           -1.7733784684952985886d+04,-4.0320340761145482298d+05/
  data   g/-2.5064972445877992730d+02, 2.9865713163054025489d+04, &
           -1.6128136304458193998d+06/
!
!  Coefficients for  1.0 < ARG
!
  data  pp/ 1.1394980557384778174d+02, 3.6832589957340267940d+03, &
            3.1075408980684392399d+04, 1.0577068948034021957d+05, &
            1.7398867902565686251d+05, 1.5097646353289914539d+05, &
            7.1557062783764037541d+04, 1.8321525870183537725d+04, &
            2.3444738764199315021d+03, 1.1600249425076035558d+02/
  data  qq/ 2.0013443064949242491d+02, 4.4329628889746408858d+03, &
            3.1474655750295278825d+04, 9.7418829762268075784d+04, &
            1.5144644673520157801d+05, 1.2689839587977598727d+05, &
            5.8824616785857027752d+04, 1.4847228371802360957d+04, &
            1.8821890840982713696d+03, 9.2556599177304839811d+01/

  x = arg
!
!  0.0 < ARG <= 1.0.
!
  if ( 0.0D+00 < x ) then

    if ( x <= 1.0D+00 ) then

      temp = log ( x )

      if ( x < xsmall ) then
!
!  Return for small ARG.
!
        result = p(6) / q(2) - temp

      else

        xx = x * x

        sump = (((( &
                 p(1) &
          * xx + p(2) ) &
          * xx + p(3) ) &
          * xx + p(4) ) &
          * xx + p(5) ) &
          * xx + p(6)

        sumq = ( xx + q(1) ) * xx + q(2)
        sumf = ( ( &
                 f(1) &
          * xx + f(2) ) &
          * xx + f(3) ) &
          * xx + f(4)

        sumg = ( ( xx + g(1) ) * xx + g(2) ) * xx + g(3)

        result = sump / sumq - xx * sumf * temp / sumg - temp

        if ( jint == 2 ) then
          result = result * exp ( x )
        end if

      end if

    else if ( jint == 1 .and. xmax < x ) then
!
!  Error return for XMAX < ARG.
!
      result = 0.0D+00

    else
!
!  1.0 < ARG.
!
      xx = 1.0D+00 / x
      sump = pp(1)
      do i = 2, 10
        sump = sump * xx + pp(i)
      end do

      sumq = xx
      do i = 1, 9
        sumq = ( sumq + qq(i) ) * xx
      end do
      sumq = sumq + qq(10)
      result = sump / sumq / sqrt ( x )

      if ( jint == 1 ) then
        result = result * exp ( -x )
      end if

    end if

  else
!
!  Error return for ARG <= 0.0.
!
    result = xinf

  end if

  return
end


subroutine caljy0 ( arg, result, jint )

!*****************************************************************************80
!
!! CALJY0 computes various J0 and Y0 Bessel functions.
!
!  Discussion:
!
!    This routine computes zero-order Bessel functions of the first and
!    second kind (J0 and Y0), for real arguments X, where 0 < X <= XMAX
!    for Y0, and |X| <= XMAX for J0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  If JINT = 0, ARG
!    must satisfy
!     -XMAX < ARG < XMAX;
!    If JINT = 1, then ARG must satisfy
!      0 < ARG < XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    0, RESULT = J0(x);
!    1, RESULT = Y0(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    0, J0(x);
!    1, Y0(x);
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) jint
  real ( kind = 8 ) arg
  real ( kind = 8 ) ax
  real ( kind = 8 ) cons
  real ( kind = 8 ) down
  real ( kind = 8 ) eight
  real ( kind = 8 ) five5
  real ( kind = 8 ) oneov8
  real ( kind = 8 ) pi2
  real ( kind = 8 ) pj0(7)
  real ( kind = 8 ) pj1(8)
  real ( kind = 8 ) plg(4)
  real ( kind = 8 ) prod
  real ( kind = 8 ) py0(6)
  real ( kind = 8 ) py1(7)
  real ( kind = 8 ) py2(8)
  real ( kind = 8 ) p0(6)
  real ( kind = 8 ) p1(6)
  real ( kind = 8 ) p17
  real ( kind = 8 ) qj0(5)
  real ( kind = 8 ) qj1(7)
  real ( kind = 8 ) qlg(4)
  real ( kind = 8 ) qy0(5)
  real ( kind = 8 ) qy1(6)
  real ( kind = 8 ) qy2(7)
  real ( kind = 8 ) q0(5)
  real ( kind = 8 ) q1(5)
  real ( kind = 8 ) resj
  real ( kind = 8 ) result
  real ( kind = 8 ) r0
  real ( kind = 8 ) r1
  real ( kind = 8 ) sixty4
  real ( kind = 8 ) three
  real ( kind = 8 ) twopi
  real ( kind = 8 ) twopi1
  real ( kind = 8 ) twopi2
  real ( kind = 8 ) two56
  real ( kind = 8 ) up
  real ( kind = 8 ) w
  real ( kind = 8 ) wsq
  real ( kind = 8 ) xden
  real ( kind = 8 ) xinf
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xnum
  real ( kind = 8 ) xsmall
  real ( kind = 8 ) xj0
  real ( kind = 8 ) xj1
  real ( kind = 8 ) xj01
  real ( kind = 8 ) xj02
  real ( kind = 8 ) xj11
  real ( kind = 8 ) xj12
  real ( kind = 8 ) xy
  real ( kind = 8 ) xy0
  real ( kind = 8 ) xy01
  real ( kind = 8 ) xy02
  real ( kind = 8 ) xy1
  real ( kind = 8 ) xy11
  real ( kind = 8 ) xy12
  real ( kind = 8 ) xy2
  real ( kind = 8 ) xy21
  real ( kind = 8 ) xy22
  real ( kind = 8 ) z
  real ( kind = 8 ) zsq
!
!  Mathematical constants
!  CONS = ln(.5) + Euler's gamma
!
  data three /3.0d0 /
  data eight /8.0d0/
  data five5 / 5.5d0 /
  data sixty4 /64.0d0 /
  data oneov8 /0.125d0 /
  data p17 /1.716d-1/
  data two56 /256.0d0/
  data cons / -1.1593151565841244881d-1/
  data pi2 /6.3661977236758134308d-1/
  data twopi /6.2831853071795864769d0/
  data twopi1 /6.28125d0 /
  data twopi2 / 1.9353071795864769253d-3/
!
!  Machine-dependent constants
!
  data xmax /1.07d+09/
  data xsmall /9.31d-10/
  data xinf /1.7d+38/
!
!  Zeroes of Bessel functions
!
  data xj0 /2.4048255576957727686d+0/
  data xj1 /5.5200781102863106496d+0/
  data xy0 /8.9357696627916752158d-1/
  data xy1 /3.9576784193148578684d+0/
  data xy2 /7.0860510603017726976d+0/
  data xj01 / 616.0d+0/
  data xj02 /-1.4244423042272313784d-03/
  data xj11 /1413.0d+0/
  data xj12 / 5.4686028631064959660d-04/
  data xy01 / 228.0d+0/
  data xy02 / 2.9519662791675215849d-03/
  data xy11 /1013.0d+0/
  data xy12 / 6.4716931485786837568d-04/
  data xy21 /1814.0d+0/
  data xy22 / 1.1356030177269762362d-04/
!
!  Coefficients for rational approximation to ln(x/a)
!
  data plg/-2.4562334077563243311d+01,2.3642701335621505212d+02, &
           -5.4989956895857911039d+02,3.5687548468071500413d+02/
  data qlg/-3.5553900764052419184d+01,1.9400230218539473193d+02, &
           -3.3442903192607538956d+02,1.7843774234035750207d+02/
!
!  Coefficients for rational approximation of
!  J0(X) / (X**2 - XJ0**2),  XSMALL < |X| <= 4.0
!
  data pj0/6.6302997904833794242d+06,-6.2140700423540120665d+08, &
           2.7282507878605942706d+10,-4.1298668500990866786d+11, &
          -1.2117036164593528341d-01, 1.0344222815443188943d+02, &
          -3.6629814655107086448d+04/
  data qj0/4.5612696224219938200d+05, 1.3985097372263433271d+08, &
           2.6328198300859648632d+10, 2.3883787996332290397d+12, &
           9.3614022392337710626d+02/
!
!  Coefficients for rational approximation of
!  J0(X) / (X**2 - XJ1**2), 4.0 < |X| <= 8.0
!
  data pj1/4.4176707025325087628d+03, 1.1725046279757103576d+04, &
           1.0341910641583726701d+04,-7.2879702464464618998d+03, &
          -1.2254078161378989535d+04,-1.8319397969392084011d+03, &
           4.8591703355916499363d+01, 7.4321196680624245801d+02/
  data qj1/3.3307310774649071172d+02,-2.9458766545509337327d+03, &
           1.8680990008359188352d+04,-8.4055062591169562211d+04, &
           2.4599102262586308984d+05,-3.5783478026152301072d+05, &
          -2.5258076240801555057d+01/
!
!  Coefficients for rational approximation of
!  (Y0(X) - 2 LN(X/XY0) J0(X)) / (X**2 - XY0**2),
!  XSMALL < |X| <= 3.0
!
  data py0/1.0102532948020907590d+04,-2.1287548474401797963d+06, &
           2.0422274357376619816d+08,-8.3716255451260504098d+09, &
           1.0723538782003176831d+11,-1.8402381979244993524d+01/
  data qy0/6.6475986689240190091d+02, 2.3889393209447253406d+05, &
           5.5662956624278251596d+07, 8.1617187777290363573d+09, &
           5.8873865738997033405d+11/
!
!  Coefficients for rational approximation of
!  (Y0(X) - 2 LN(X/XY1) J0(X)) / (X**2 - XY1**2),
!  3.0 < |X| <= 5.5
!
  data py1/-1.4566865832663635920d+04, 4.6905288611678631510d+06, &
           -6.9590439394619619534d+08, 4.3600098638603061642d+10, &
           -5.5107435206722644429d+11,-2.2213976967566192242d+13, &
            1.7427031242901594547d+01/
  data qy1/ 8.3030857612070288823d+02, 4.0669982352539552018d+05, &
            1.3960202770986831075d+08, 3.4015103849971240096d+10, &
            5.4266824419412347550d+12, 4.3386146580707264428d+14/
!
!  Coefficients for rational approximation of
!  (Y0(X) - 2 LN(X/XY2) J0(X)) / (X**2 - XY2**2),
!  5.5 < |X| <= 8.0
!
  data py2/ 2.1363534169313901632d+04,-1.0085539923498211426d+07, &
            2.1958827170518100757d+09,-1.9363051266772083678d+11, &
           -1.2829912364088687306d+11, 6.7016641869173237784d+14, &
           -8.0728726905150210443d+15,-1.7439661319197499338d+01/
  data qy2/ 8.7903362168128450017d+02, 5.3924739209768057030d+05, &
            2.4727219475672302327d+08, 8.6926121104209825246d+10, &
            2.2598377924042897629d+13, 3.9272425569640309819d+15, &
            3.4563724628846457519d+17/
!
!  Coefficients for Hart,s approximation, 8.0 < |X|.
!
  data p0/3.4806486443249270347d+03, 2.1170523380864944322d+04, &
          4.1345386639580765797d+04, 2.2779090197304684302d+04, &
          8.8961548424210455236d-01, 1.5376201909008354296d+02/
  data q0/3.5028735138235608207d+03, 2.1215350561880115730d+04, &
          4.1370412495510416640d+04, 2.2779090197304684318d+04, &
          1.5711159858080893649d+02/
  data p1/-2.2300261666214198472d+01,-1.1183429920482737611d+02, &
          -1.8591953644342993800d+02,-8.9226600200800094098d+01, &
          -8.8033303048680751817d-03,-1.2441026745835638459d+00/
  data q1/1.4887231232283756582d+03, 7.2642780169211018836d+03, &
          1.1951131543434613647d+04, 5.7105024128512061905d+03, &
          9.0593769594993125859d+01/
!
!  Check for error conditions.
!
  ax = abs ( arg )

  if ( jint == 1 .and. arg <= 0.0D+00 ) then
    result = -xinf
    return
  else if ( xmax < ax ) then
    result = 0.0D+00
    return
  end if

  if ( eight < ax ) then
    go to 800
  end if

  if ( ax <= xsmall ) then
    if ( jint == 0 ) then
      result = 1.0D+00
    else
      result = pi2 * ( log ( ax ) + cons )
    end if
    return
  end if
!
!  Calculate J0 for appropriate interval, preserving
!  accuracy near the zero of J0.
!
  zsq = ax * ax

  if ( ax <= 4.0D+00 ) then
    xnum = ( pj0(5) * zsq + pj0(6) ) * zsq + pj0(7)
    xden = zsq + qj0(5)
    do i = 1, 4
      xnum = xnum * zsq + pj0(i)
      xden = xden * zsq + qj0(i)
    end do
    prod = ( ( ax - xj01 / two56 ) - xj02 ) * ( ax + xj0 )
  else
    wsq = 1.0D+00 - zsq / sixty4
    xnum = pj1(7) * wsq + pj1(8)
    xden = wsq + qj1(7)
    do i = 1, 6
      xnum = xnum * wsq + pj1(i)
      xden = xden * wsq + qj1(i)
    end do
    prod = ( ax + xj1 ) * ( ( ax - xj11 / two56 ) - xj12 )
  end if

  result = prod * xnum / xden

  if ( jint == 0 ) then
    return
  end if
!
!  Calculate Y0.  First find  RESJ = pi/2 ln(x/xn) J0(x),
!  where xn is a zero of Y0.
!
  if ( ax <= three ) then
    up = ( ax - xy01 / two56 ) - xy02
    xy = xy0
  else if ( ax <= five5 ) then
    up = ( ax - xy11 / two56 ) - xy12
    xy = xy1
  else
    up = ( ax - xy21 / two56 ) - xy22
    xy = xy2
  end if

  down = ax + xy

  if ( abs ( up ) < p17 * down ) then
    w = up / down
    wsq = w * w
    xnum = plg(1)
    xden = wsq + qlg(1)
    do i = 2, 4
      xnum = xnum * wsq + plg(i)
      xden = xden * wsq + qlg(i)
    end do
    resj = pi2 * result * w * xnum / xden
  else
    resj = pi2 * result * log ( ax / xy )
  end if
!
!  Now calculate Y0 for appropriate interval, preserving
!  accuracy near the zero of Y0.
!
  if ( ax <= three ) then
    xnum = py0(6) * zsq + py0(1)
    xden = zsq + qy0(1)
    do i = 2, 5
      xnum = xnum * zsq + py0(i)
      xden = xden * zsq + qy0(i)
    end do
  else if ( ax <= five5 ) then
    xnum = py1(7) * zsq + py1(1)
    xden = zsq + qy1(1)
    do i = 2, 6
      xnum = xnum * zsq + py1(i)
      xden = xden * zsq + qy1(i)
    end do
  else
    xnum = py2(8) * zsq + py2(1)
    xden = zsq + qy2(1)
    do i = 2, 7
      xnum = xnum * zsq + py2(i)
      xden = xden * zsq + qy2(i)
    end do
  end if

  result = resj + up * down * xnum / xden

  return
!
!  Calculate J0 or Y0 for 8.0 < |ARG|.
!
  800 continue

  z = eight / ax
  w = ax / twopi
  w = aint ( w ) + oneov8
  w = ( ax - w * twopi1 ) - w * twopi2
  zsq = z * z
  xnum = p0(5) * zsq + p0(6)
  xden = zsq + q0(5)
  up = p1(5) * zsq + p1(6)
  down = zsq + q1(5)

  do i = 1, 4
    xnum = xnum * zsq + p0(i)
    xden = xden * zsq + q0(i)
    up = up * zsq + p1(i)
    down = down * zsq + q1(i)
  end do

  r0 = xnum / xden
  r1 = up / down

  if ( jint == 0 ) then
    result = sqrt ( pi2 / ax ) &
      * ( r0 * cos ( w ) - z * r1 * sin ( w ) )
  else
    result = sqrt ( pi2 / ax ) &
      * ( r0 * sin ( w ) + z * r1 * cos ( w ) )
  end if

  return
end

subroutine caljy1 ( arg, result, jint )

!*****************************************************************************80
!
!! CALJY1 computes various J1 and Y1 Bessel functions.
!
!  Discussion:
!
!    This routine computes first-order Bessel functions of the first and
!    second kind (J1 and Y1), for real arguments X, where 0 < X <= XMAX
!    for Y1, and |X| <= XMAX for J1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thatcher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  If JINT = 0, ARG
!    must satisfy
!     -XMAX < ARG < XMAX;
!    If JINT = 1, then ARG must satisfy
!      0 < ARG < XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    0, RESULT = J1(x);
!    1, RESULT = Y1(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    0, J1(x);
!    1, Y1(x);
!
  implicit none

  real ( kind = 8 ) arg
  real ( kind = 8 ) ax
  real ( kind = 8 ) down
  real ( kind = 8 ) eight
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jint
  real ( kind = 8 ) pi2
  real ( kind = 8 ) pj0(7)
  real ( kind = 8 ) pj1(8)
  real ( kind = 8 ) plg(4)
  real ( kind = 8 ) prod
  real ( kind = 8 ) py0(7)
  real ( kind = 8 ) py1(9)
  real ( kind = 8 ) p0(6)
  real ( kind = 8 ) p1(6)
  real ( kind = 8 ) p17
  real ( kind = 8 ) qj0(5)
  real ( kind = 8 ) qj1(7)
  real ( kind = 8 ) qlg(4)
  real ( kind = 8 ) qy0(6)
  real ( kind = 8 ) qy1(8)
  real ( kind = 8 ) q0(6)
  real ( kind = 8 ) q1(6)
  real ( kind = 8 ) resj
  real ( kind = 8 ) result
  real ( kind = 8 ) rtpi2
  real ( kind = 8 ) r0
  real ( kind = 8 ) r1
  real ( kind = 8 ) throv8
  real ( kind = 8 ) twopi
  real ( kind = 8 ) twopi1
  real ( kind = 8 ) twopi2
  real ( kind = 8 ) two56
  real ( kind = 8 ) up
  real ( kind = 8 ) w
  real ( kind = 8 ) wsq
  real ( kind = 8 ) xden
  real ( kind = 8 ) xinf
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xnum
  real ( kind = 8 ) xsmall
  real ( kind = 8 ) xj0
  real ( kind = 8 ) xj1
  real ( kind = 8 ) xj01
  real ( kind = 8 ) xj02
  real ( kind = 8 ) xj11
  real ( kind = 8 ) xj12
  real ( kind = 8 ) xy
  real ( kind = 8 ) xy0
  real ( kind = 8 ) xy01
  real ( kind = 8 ) xy02
  real ( kind = 8 ) xy1
  real ( kind = 8 ) xy11
  real ( kind = 8 ) xy12
  real ( kind = 8 ) z
  real ( kind = 8 ) zsq
!
!  Mathematical constants
!
  data eight /8.0d0/
  data throv8 /0.375d0/
  data pi2 /6.3661977236758134308d-1/
  data p17 /1.716d-1/
  data twopi /6.2831853071795864769d+0/
  data twopi1 /6.28125d0/
  data twopi2 /1.9353071795864769253d-03/
  data two56 /256.0d+0/
  data rtpi2 /7.9788456080286535588d-1/
!
!  Machine-dependent constants
!
  data xmax /1.07d+09/
  data xsmall /9.31d-10/
  data xinf /1.7d+38/
!
!  Zeroes of Bessel functions
!
  data xj0 /3.8317059702075123156d+0/
  data xj1 /7.0155866698156187535d+0/
  data xy0 /2.1971413260310170351d+0/
  data xy1 /5.4296810407941351328d+0/
  data xj01 / 981.0d+0/
  data xj02 /-3.2527979248768438556d-04/
  data xj11 /1796.0d+0/
  data xj12 /-3.8330184381246462950d-05/
  data xy01 / 562.0d+0/
  data xy02 / 1.8288260310170351490d-03/
  data xy11 /1390.0d+0/
  data xy12 /-6.4592058648672279948d-06/
!
!  Coefficients for rational approximation to ln(x/a)
!
  data plg/-2.4562334077563243311d+01,2.3642701335621505212d+02, &
           -5.4989956895857911039d+02,3.5687548468071500413d+02/
  data qlg/-3.5553900764052419184d+01,1.9400230218539473193d+02, &
           -3.3442903192607538956d+02,1.7843774234035750207d+02/
!
!  Coefficients for rational approximation of
!  J1(X) / (X * (X**2 - XJ0**2)), XSMALL < |X| <=  4.0
!
  data pj0/9.8062904098958257677d+05,-1.1548696764841276794d+08, &
         6.6781041261492395835d+09,-1.4258509801366645672d+11, &
        -4.4615792982775076130d+03, 1.0650724020080236441d+01, &
        -1.0767857011487300348d-02/
  data qj0/5.9117614494174794095d+05, 2.0228375140097033958d+08, &
         4.2091902282580133541d+10, 4.1868604460820175290d+12, &
         1.0742272239517380498d+03/
!
!  Coefficients for rational approximation of
!  J1(X) / (X * (X**2 - XJ1**2)), 4.0 < |X| <= 8.0
!
  data pj1/4.6179191852758252280d+00,-7.1329006872560947377d+03, &
         4.5039658105749078904d+06,-1.4437717718363239107d+09, &
         2.3569285397217157313d+11,-1.6324168293282543629d+13, &
         1.1357022719979468624d+14, 1.0051899717115285432d+15/
  data qj1/1.1267125065029138050d+06, 6.4872502899596389593d+08, &
         2.7622777286244082666d+11, 8.4899346165481429307d+13, &
         1.7128800897135812012d+16, 1.7253905888447681194d+18, &
         1.3886978985861357615d+03/
!
!  Coefficients for rational approximation of
!  (Y1(X) - 2 LN(X/XY0) J1(X)) / (X**2 - XY0**2),
!  XSMALL < |X| <=  4.0
!
  data py0/2.2157953222280260820d+05,-5.9157479997408395984d+07, &
           7.2144548214502560419d+09,-3.7595974497819597599d+11, &
           5.4708611716525426053d+12, 4.0535726612579544093d+13, &
          -3.1714424660046133456d+02/
  data qy0/8.2079908168393867438d+02, 3.8136470753052572164d+05, &
           1.2250435122182963220d+08, 2.7800352738690585613d+10, &
           4.1272286200406461981d+12, 3.0737873921079286084d+14/
!
!  Coefficients for rational approximation of
!  (Y1(X) - 2 LN(X/XY1) J1(X)) / (X**2 - XY1**2),
!  4.0 < |X| <= 8.0
!
  data py1/ 1.9153806858264202986d+06,-1.1957961912070617006d+09, &
            3.7453673962438488783d+11,-5.9530713129741981618d+13, &
            4.0686275289804744814d+15,-2.3638408497043134724d+16, &
           -5.6808094574724204577d+18, 1.1514276357909013326d+19, &
           -1.2337180442012953128d+03/
  data qy1/ 1.2855164849321609336d+03, 1.0453748201934079734d+06, &
            6.3550318087088919566d+08, 3.0221766852960403645d+11, &
            1.1187010065856971027d+14, 3.0837179548112881950d+16, &
            5.6968198822857178911d+18, 5.3321844313316185697d+20/
!
!  Coefficients for Hart's approximation, 8.0 < |X|.
!
  data p0/-1.0982405543459346727d+05,-1.5235293511811373833d+06, &
           -6.6033732483649391093d+06,-9.9422465050776411957d+06, &
           -4.4357578167941278571d+06,-1.6116166443246101165d+03/
  data q0/-1.0726385991103820119d+05,-1.5118095066341608816d+06, &
           -6.5853394797230870728d+06,-9.9341243899345856590d+06, &
           -4.4357578167941278568d+06,-1.4550094401904961825d+03/
  data p1/ 1.7063754290207680021d+03, 1.8494262873223866797d+04, &
            6.6178836581270835179d+04, 8.5145160675335701966d+04, &
            3.3220913409857223519d+04, 3.5265133846636032186d+01/
  data q1/ 3.7890229745772202641d+04, 4.0029443582266975117d+05, &
            1.4194606696037208929d+06, 1.8194580422439972989d+06, &
            7.0871281941028743574d+05, 8.6383677696049909675d+02/
!
!  Check for error conditions.
!
  ax = abs ( arg )

  if ( jint == 1 .and. ( arg <= 0.0D+00 .or. &
    ( arg < 0.5D+00 .and. ax * xinf < pi2 ) ) ) then
    result = -xinf
    return
  else if ( xmax < ax ) then
    result = 0.0D+00
    return
  end if

  if ( eight < ax ) then
    go to 800
  else if ( ax <= xsmall ) then
    if ( jint == 0 ) then
      result = arg * 0.5D+00
    else
      result = -pi2 / ax
    end if
    return
  end if
!
!  Calculate J1 for appropriate interval, preserving
!  accuracy near the zero of J1.
!
  zsq = ax * ax

  if ( ax <= 4.0D+00 ) then
    xnum = ( pj0(7) * zsq + pj0(6) ) * zsq + pj0(5)
    xden = zsq + qj0(5)
    do i = 1, 4
      xnum = xnum * zsq + pj0(i)
      xden = xden * zsq + qj0(i)
    end do
    prod = arg * ( ( ax - xj01 / two56 ) - xj02 ) * ( ax + xj0 )
  else
    xnum = pj1(1)
    xden = ( zsq + qj1(7) ) * zsq + qj1(1)
    do i = 2, 6
      xnum = xnum * zsq + pj1(i)
      xden = xden * zsq + qj1(i)
    end do
    xnum = xnum * ( ax - eight ) * ( ax + eight ) + pj1(7)
    xnum = xnum * ( ax - 4.0D+00 ) * ( ax + 4.0D+00 ) + pj1(8)
    prod = arg * ( ( ax - xj11 / two56 ) - xj12 ) * ( ax + xj1 )
  end if

  result = prod * ( xnum / xden )

  if ( jint == 0 ) then
    return
  end if
!
!  Calculate Y1.  First find RESJ = pi/2 ln(x/xn) J1(x),
!  where xn is a zero of Y1.
!
  if ( ax <= 4.0D+00 ) then
    up = ( ax - xy01 / two56 ) - xy02
    xy = xy0
  else
    up = ( ax - xy11 / two56 ) - xy12
    xy = xy1
  end if

  down = ax + xy

  if ( abs ( up ) < p17 * down ) then
    w = up / down
    wsq = w * w
    xnum = plg(1)
    xden = wsq + qlg(1)
    do i = 2, 4
      xnum = xnum * wsq + plg(i)
      xden = xden * wsq + qlg(i)
    end do
    resj = pi2 * result * w * xnum / xden
  else
    resj = pi2 * result * log ( ax / xy )
  end if
!
!  Now calculate Y1 for appropriate interval, preserving
!  accuracy near the zero of Y1.
!
  if ( ax <= 4.0D+00 ) then
    xnum = py0(7) * zsq + py0(1)
    xden = zsq + qy0(1)
    do i = 2, 6
      xnum = xnum * zsq + py0(i)
      xden = xden * zsq + qy0(i)
    end do
  else
    xnum = py1(9) * zsq + py1(1)
    xden = zsq + qy1(1)
    do i = 2, 8
      xnum = xnum * zsq + py1(i)
      xden = xden * zsq + qy1(i)
    end do
  end if

  result = resj + ( up * down / ax ) * xnum / xden
  return
!
!  Calculate J1 or Y1 for 8.0 < |ARG|.
!
  800 continue

  z = eight / ax
  w = aint ( ax / twopi ) + throv8
  w = ( ax - w * twopi1 ) - w * twopi2
  zsq = z * z
  xnum = p0(6)
  xden = zsq + q0(6)
  up = p1(6)
  down = zsq + q1(6)

  do i = 1, 5
    xnum = xnum * zsq + p0(i)
    xden = xden * zsq + q0(i)
    up = up * zsq + p1(i)
    down = down * zsq + q1(i)
  end do

  r0 = xnum / xden
  r1 = up / down

  if ( jint == 0 ) then
    result = ( rtpi2 / sqrt ( ax ) ) &
      * ( r0 * cos ( w ) - z * r1 * sin ( w ) )
  else
    result = ( rtpi2 / sqrt ( ax ) ) &
      * ( r0 * sin ( w ) + z * r1 * cos ( w ) )
  end if

  if ( jint == 0 .and. arg < 0.0D+00 ) then
    result = -result
  end if

  return
end

subroutine calck1 ( arg, result, jint )

!*****************************************************************************80
!
!! CALCK1 computes various K1 Bessel functions.
!
!  Discussion:
!
!    This routine computes modified Bessel functions of the second kind
!    and order one, K1(X) and EXP(X)*K1(X), for real arguments X.
!
!    The main computation evaluates slightly modified forms of near
!    minimax rational approximations generated by Russon and Blair,
!    Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
!    1969.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody, Laura Stoltz.
!    FORTRAN90 version by John Burkardt.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ARG, the argument.  XLEAST < ARG is
!    always required.  If JINT = 1, then the argument must also be
!    less than XMAX.
!
!    Output, real ( kind = 8 ) RESULT, the value of the function,
!    which depends on the input value of JINT:
!    1, RESULT = K1(x);
!    2, RESULT = exp(x) * K1(x);
!
!    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
!    1, K1(x);
!    2, exp(x) * K1(x);
!
  implicit none

  real ( kind = 8 ) arg
  real ( kind = 8 ) f(5)
  real ( kind = 8 ) g(3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jint
  real ( kind = 8 ) p(5)
  real ( kind = 8 ) pp(11)
  real ( kind = 8 ) q(3)
  real ( kind = 8 ) qq(9)
  real ( kind = 8 ) result
  real ( kind = 8 ) sumf
  real ( kind = 8 ) sumg
  real ( kind = 8 ) sump
  real ( kind = 8 ) sumq
  real ( kind = 8 ) x
  real ( kind = 8 ) xinf
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xleast
  real ( kind = 8 ) xsmall
  real ( kind = 8 ) xx
!
!  Machine-dependent constants
!
  data xleast /2.23d-308/
  data xsmall /1.11d-16/
  data xinf /1.79d+308/
  data xmax /705.343d+0/
!
!  Coefficients for  XLEAST <=  ARG  <= 1.0
!
  data   p/ 4.8127070456878442310d-1, 9.9991373567429309922d+1, &
            7.1885382604084798576d+3, 1.7733324035147015630d+5, &
            7.1938920065420586101d+5/
  data   q/-2.8143915754538725829d+2, 3.7264298672067697862d+4, &
           -2.2149374878243304548d+6/
  data   f/-2.2795590826955002390d-1,-5.3103913335180275253d+1, &
           -4.5051623763436087023d+3,-1.4758069205414222471d+5, &
           -1.3531161492785421328d+6/
  data   g/-3.0507151578787595807d+2, 4.3117653211351080007d+4, &
           -2.7062322985570842656d+6/
!
!  Coefficients for  1.0 < ARG
!
  data  pp/ 6.4257745859173138767d-2, 7.5584584631176030810d+0, &
            1.3182609918569941308d+2, 8.1094256146537402173d+2, &
            2.3123742209168871550d+3, 3.4540675585544584407d+3, &
            2.8590657697910288226d+3, 1.3319486433183221990d+3, &
            3.4122953486801312910d+2, 4.4137176114230414036d+1, &
            2.2196792496874548962d+0/
  data  qq/ 3.6001069306861518855d+1, 3.3031020088765390854d+2, &
            1.2082692316002348638d+3, 2.1181000487171943810d+3, &
            1.9448440788918006154d+3, 9.6929165726802648634d+2, &
            2.5951223655579051357d+2, 3.4552228452758912848d+1, &
            1.7710478032601086579d+0/

  x = arg
!
!  Error return for ARG < XLEAST.
!
  if ( x < xleast ) then

    result = xinf
!
!  XLEAST <= ARG <= 1.0.
!
  else if ( x <= 1.0D+00 ) then

    if ( x < xsmall ) then
!
!  Return for small ARG.
!
      result = 1.0D+00 / x

    else

      xx = x * x

      sump = (((( &
               p(1) &
        * xx + p(2) ) &
        * xx + p(3) ) &
        * xx + p(4) ) &
        * xx + p(5) ) &
        * xx + q(3)

      sumq = (( &
          xx + q(1) ) &
        * xx + q(2) ) &
        * xx + q(3)

      sumf = ((( &
               f(1) &
        * xx + f(2) ) &
        * xx + f(3) ) &
        * xx + f(4) ) &
        * xx + f(5)

      sumg = (( &
          xx + g(1) ) &
        * xx + g(2) ) &
        * xx + g(3)

      result = ( xx * log ( x ) * sumf / sumg + sump / sumq ) / x

      if ( jint == 2 ) then
        result = result * exp ( x )
      end if

    end if

  else if ( jint == 1 .and. xmax < x ) then
!
!  Error return for XMAX < ARG.
!
    result = 0.0D+00

  else
!
!  1.0 < ARG.
!
    xx = 1.0D+00 / x

    sump = pp(1)
    do i = 2, 11
      sump = sump * xx + pp(i)
    end do

    sumq = xx
    do i = 1, 8
      sumq = ( sumq + qq(i) ) * xx
    end do
    sumq = sumq + qq(9)

    result = sump / sumq / sqrt ( x )

    if ( jint == 1 ) then
      result = result * exp ( -x )
    end if

  end if

  return
end

