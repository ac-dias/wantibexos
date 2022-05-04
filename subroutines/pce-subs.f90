
subroutine slme(ndim,ctemp,trial,bbpflux,sepflux,absc,isolar,pin,fr,pm,jsc1,j01,j1,j2,vmax,voc,ff,pce)

	implicit none

	integer :: ndim,m,kp,i
	double precision,dimension(ndim,2) :: bbpflux,sepflux,isolar
	double precision,dimension(ndim,2) :: absc,auxvec
	double precision :: pce,fr,pin,pm
	double precision :: j0r,j0,jsc,thickness,jsc1,j01,j1,j2
	double precision :: voc, vmax, aux1,trial,ff

	double precision,dimension(1) :: x0,x,y0,y
	double precision :: fx
	double precision,external :: negpower,vocf
	
	double precision,parameter :: c =299792458.0 !speed of light, m/s
	double precision,parameter :: h =6.62607004081E-34 !Planck's constant J*s (W)
	double precision,parameter :: heV =4.135667516E-15  !Planck's constant eV*s
	double precision,parameter :: k =1.3806485279E-23  !Boltzmann's constant J/K
	double precision,parameter :: keV =8.617330350E-5  !Boltzmann's constant eV/K
	double precision,parameter :: e =1.602176620898E-19 !Coulomb

	double precision,parameter :: jev=6.24150913E18 !conversion factor joule to eV
	double precision,parameter :: pi=acos(-1.0)	
	double precision :: tkv,ctemp
	
	double precision,parameter :: delta_tol= 10**(-7.0d0)
	double precision,parameter :: delta= 0.2D+00
	integer,parameter :: k_max = 200 	

	common /dados/ jsc,j0,tkv 
	
	tkv = ctemp
	
	!call absorbance(ndim,iabscoef,thickness,absc)
	!absc=1.0

	!calculando j0r e j0
	auxvec = 0.0
	do i=1,ndim
	auxvec(i,1) = sepflux(i,1)
	auxvec(i,2) = bbpflux(i,2)*absc(i,2) 
	end do


	call integral1D(ndim,auxvec,j0r)
	j0r = e*pi*j0r
	j0 = j0r/fr

	j01= j0

	!write(*,*) "J0r:",j0r
	!write(*,*) "J0:",j0

	!write(*,*) fr

	!calculando jsc
	auxvec = 0.0
	do i=1,ndim
	auxvec(i,1) = isolar(i,1)
	auxvec(i,2) = sepflux(i,2)*absc(i,2) 
	end do
		
	call integral1D(ndim,auxvec,jsc)
	jsc = jsc*e

	!write(*,*) "Jsc:",jsc



	m=1
	x0(1)=trial
	call compass_search ( negpower, m, x0, delta_tol,delta, k_max, x, fx, kp )
	!x(1) = 0.08113602496880366
	pm = -negpower(m, x(1))
	

	vmax = x(1)
	

	aux1= (e*vmax)/(k*tkv)

	j2= j0*(dexp(aux1)-1.0)
	j1= jsc-j0*(dexp(aux1)-1.0)

	jsc1 = jsc

	!write(*,*) 'Voc:',x(1)
	!write(*,*) 'P_in:',pin
	!write(*,*) 'P_m:',pm
	!pm=1.0
	pce = pm/pin

	!write(*,*) pce
	
	m=1
	y0(1)=vmax
	call compass_search ( vocf, m, y0, delta_tol,delta, k_max, y, fx, kp )
	voc = y(1)	

	ff = (j1*vmax)/(jsc1*voc)

end subroutine slme

subroutine sq(ndim,ctemp,bbpflux,sepflux,iabscoef,isolar,pin,fr,pm,jsc1,j01,j1,j2,vmax,voc,ff,pce)


	implicit none

	integer :: ndim,m,kp,i
	double precision,dimension(ndim,2) :: bbpflux,sepflux,iabscoef,isolar
	double precision,dimension(ndim,2) :: auxvec
	double precision :: pce,fr,pin,pm
	double precision :: j0r,j0,jsc,thickness,jsc1,j01,j1,j2
	double precision :: voc,vmax,aux1,ff

	double precision,dimension(1) :: x0,x,y,y0
	double precision :: fx
	double precision,external :: negpower,vocf
	
	double precision,parameter :: c =299792458.0 !speed of light, m/s
	double precision,parameter :: h =6.62607004081E-34 !Planck's constant J*s (W)
	double precision,parameter :: heV =4.135667516E-15  !Planck's constant eV*s
	double precision,parameter :: k =1.3806485279E-23  !Boltzmann's constant J/K
	double precision,parameter :: keV =8.617330350E-5  !Boltzmann's constant eV/K
	double precision,parameter :: e =1.602176620898E-19 !Coulomb

	double precision,parameter :: jev=6.24150913E18 !conversion factor joule to eV
	double precision,parameter :: pi=acos(-1.0)	
	double precision :: tkv,ctemp
	
	double precision,parameter :: delta_tol= 10**(-7.0d0)
	double precision,parameter :: delta= 0.2D+00
	integer,parameter :: k_max = 200 	

	common /dados/ jsc,j0,tkv 
	
	tkv = ctemp	

	!calculando j0r e j0
	auxvec = 0.0
	do i=1,ndim
	auxvec(i,1) = sepflux(i,1)
	auxvec(i,2) = bbpflux(i,2)*iabscoef(i,2) 
	end do

	call integral1D(ndim,auxvec,j0r)
	j0r = e*pi*j0r
	j0 = j0r/fr

	j01= j0

	!write(*,*) "J0r:",j0r
	!write(*,*) "J0:",j0

	!write(*,*) fr

	!calculando jsc
	auxvec = 0.0
	do i=1,ndim
	auxvec(i,1) = isolar(i,1)
	auxvec(i,2) = sepflux(i,2)*iabscoef(i,2) 
	end do
		
	call integral1D(ndim,auxvec,jsc)
	jsc = jsc*e

	!write(*,*) "Jsc:",jsc

	m=1
	x0(1)=0.1
	call compass_search ( negpower, m, x0, delta_tol,delta, k_max, x, fx, kp )
	!x(1) = 0.08113602496880366
	pm = -negpower(m, x(1))
	vmax = x(1)
	
	aux1= (e*vmax)/(k*tkv)

	j2= j0*(dexp(aux1)-1.0)
	j1= jsc-j0*(dexp(aux1)-1.0)

	jsc1 = jsc

	!write(*,*) 'Voc:',x(1)
	!write(*,*) 'P_in:',pin
	!write(*,*) 'P_m:',pm
	!pm=1.0
	pce = pm/pin

	!write(*,*) pce

	m=1
	y0(1)=vmax
	call compass_search ( vocf, m, y0, delta_tol,delta, k_max, y, fx, kp )
	voc = y(1)


	ff = (j1*vmax)/(jsc1*voc)


end subroutine sq


function vocf(m, x)

	implicit none

	integer :: m
	double precision :: x
	double precision :: jsc,j0
	double precision :: vocf,tkv
	
	double precision,parameter :: e =1.602176620898E-19 !Coulomb
	double precision,parameter :: k =1.3806485279E-23  !Boltzmann's constant J/K	

	double precision :: aux1

	common /dados/ jsc,j0,tkv 


	aux1= (e*x)/(k*tkv)

	vocf=abs(jsc-j0*(dexp(aux1)-1.0))
		

end function

function negpower(m, x)

	implicit none

	integer :: m
	double precision :: x
	double precision :: jsc,j0
	double precision :: negpower,tkv

	double precision :: aux1
	
	double precision,parameter :: e =1.602176620898E-19 !Coulomb
	double precision,parameter :: k =1.3806485279E-23  !Boltzmann's constant J/K		

	common /dados/ jsc,j0,tkv


	aux1= (e*x)/(k*tkv)

	negpower=-(jsc-j0*(dexp(aux1)-1.0))*x
		

end function

subroutine blackbody(ndim,tkv,isolarm,bbir,bbpflux)


	implicit none

	integer :: i,ndim
	double precision,dimension(ndim,2) :: isolarm
	double precision,dimension(ndim,2) :: bbir,bbpflux
	double precision,parameter :: c =299792458.0 !speed of light, m/s
	double precision,parameter :: h =6.62607004081E-34 !Planck's constant J*s (W)
	double precision,parameter :: heV =4.135667516E-15  !Planck's constant eV*s
	double precision,parameter :: k =1.3806485279E-23  !Boltzmann's constant J/K
	double precision,parameter :: keV =8.617330350E-5  !Boltzmann's constant eV/K
	double precision,parameter :: e =1.602176620898E-19 !Coulomb

	double precision,parameter :: jev=6.24150913E18 !conversion factor joule to eV
	double precision,parameter :: pi=acos(-1.0)	

	double precision :: aux1,tkv


	aux1=2.0*pi*h*c*c


	do i=1,ndim

		bbir(i,1) = isolarm(i,1)
		bbir(i,2) = (aux1/(isolarm(i,1)**5))*(1.0/(dexp(h*c/(isolarm(i,1)*k*tkv))-1.0))

		bbpflux(i,1) = isolarm(i,1)
		bbpflux(i,2) = bbir(i,2)*(isolarm(i,1)/(h*c))
	
	end do

end subroutine blackbody

subroutine absorbance2x(ndim,input,thickness,absc)


	implicit none
	
	integer :: i
	integer :: ndim
	double precision,dimension(ndim,2) :: input,absc
	double precision :: thickness
	double precision :: aux
	
	double precision,parameter :: c =299792458.0 !speed of light, m/s
	double precision,parameter :: h =6.62607004081E-34 !Planck's constant J*s (W)
	double precision,parameter :: heV =4.135667516E-15  !Planck's constant eV*s
	double precision,parameter :: k =1.3806485279E-23  !Boltzmann's constant J/K
	double precision,parameter :: keV =8.617330350E-5  !Boltzmann's constant eV/K
	double precision,parameter :: e =1.602176620898E-19 !Coulomb

	double precision,parameter :: jev=6.24150913E18 !conversion factor joule to eV
	double precision,parameter :: pi=acos(-1.0)	


	do i=1,ndim

	absc(i,1) = input(i,1)
	absc(i,2) = 1.0-dexp(-2.0*thickness*input(i,2))
	!absc(i,2) = 1.0
	
	aux = 2*absc(i,2)
	if ( aux .gt. 1.0) then
		aux = 1.0
	end if
	
	absc(i,2) = aux

	end do

end subroutine absorbance2x

subroutine absorbance(ndim,input,thickness,absc)

	implicit none
	
	integer :: i
	integer :: ndim
	double precision,dimension(ndim,2) :: input,absc
	double precision :: thickness


	do i=1,ndim

	absc(i,1) = input(i,1)
	absc(i,2) = 1.0-dexp(-2.0*thickness*input(i,2))
	!absc(i,2) = 1.0

	end do

end subroutine absorbance


subroutine interp1Da(ndim,vec,x,res)

	implicit none

	integer :: ndim,i,iaux
	double precision,dimension(ndim,2) :: vec
	double precision :: x,res

	do i=1,ndim

		if (x .le. vec(i,1)) then

			iaux=i

			go to 105
		else
			continue
		end if

	end do

105 continue
        
	res= vec(iaux-1,2) + (x-vec(iaux-1,1))*((vec(iaux,2)-vec(iaux-1,2))/(vec(iaux,1)-vec(iaux-1,1)))

end subroutine interp1Da

subroutine integral1D(ndim,vec,resultado)

	implicit none
	
	integer :: ndim,i
	double precision :: resultado,f,flag1,flag2
	double precision:: dx
	
	double precision, dimension(ndim,2) :: vec
		
	dx= abs(vec(2,1)-vec(1,1))


	flag2=0
	
	do i=1,ndim-1
	
		flag1= (vec(i,2)+vec(i+1,2))*(dx/2.)
		
		flag2=flag2+flag1
	
		!write(*,*) "progresso:",i,"/",npassos

	end do
	
	resultado=flag2


end subroutine integral1D

!subrotina de minimizacao
subroutine compass_search ( function_handle, m , x0, delta_tol, delta_init, &
  k_max, x, fx, k )

!*****************************************************************************80
!
!! COMPASS_SEARCH carries out a direct search minimization algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Tamara Kolda, Robert Michael Lewis, Virginia Torczon,
!    Optimization by Direct Search: New Perspectives on Some Classical 
!    and Modern Methods,
!    SIAM Review,
!    Volume 45, Number 3, 2003, pages 385-482. 
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) FUNCTION_HANDLE, the name of
!    a FORTRAN90 function which evaluates the function to be minimized, of the
!    form FUNCTION FUNCTION_HANDLE ( M, X ).
!
!    Input, integer ( kind = 4 ) M, the number of variables.
!
!    Input, real ( kind = 8 ) X0(M), a starting estimate for the minimizer.
!
!    Input, real ( kind = 8 ) DELTA_TOL, the smallest step size that is allowed.
!
!    Input, real ( kind = 8 ) DELTA_INIT, the starting stepsize.  
!
!    Input, integer ( kind = 4 ) K_MAX, the maximum number of steps allowed.
!
!    Output, real ( kind = 8 ) X(M), the estimated minimizer.
!
!    Output, real ( kind = 8 ) FX, the function value at X.
!
!    Output, integer ( kind = 4 ) K, the number of steps taken.
!
  implicit none

  integer ( kind = 4 ) m
  real ( kind = 8 ) :: j0, jsc

  logical decrease
  real ( kind = 8 ) delta
  real ( kind = 8 ) delta_init
  real ( kind = 8 ) delta_tol
  real ( kind = 8 ), external :: function_handle
  real ( kind = 8 ) fx
  real ( kind = 8 ) fxd
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_max
  real ( kind = 8 ) s
  real ( kind = 8 ) x(m)
  real ( kind = 8 ) x0(m)
  real ( kind = 8 ) xd(m)

  k = 0
  x(1:m) = x0(1:m)
  fx = function_handle ( m, x )

  if ( delta_tol <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMPASS_SEARCH - Fatal error!'
    write ( *, '(a)' ) '  DELTA_TOL <= 0.0.'
    write ( *, '(a,g14.6)' ) '  DELTA_TOL = ', delta_tol
    stop
  end if

  if ( delta_init <= delta_tol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'COMPASS_SEARCH - Fatal error!'
    write ( *, '(a)' ) '  DELTA_INIT < DELTA_TOL.'
    write ( *, '(a,g14.6)' ) '  DELTA_INIT = ', delta_init
    write ( *, '(a,g14.6)' ) '  DELTA_TOL = ', delta_tol
    stop
  end if

  delta = delta_init

  do while ( k < k_max )

    k = k + 1
!
!  For each coordinate direction I, seek a lower function value
!  by increasing or decreasing X(I) by DELTA.
!
    decrease = .false.
    s = + 1.0D+00
    i = 1

    do ii = 1, 2 * m

      xd = x
      xd(i) = xd(i) + s * delta
      fxd = function_handle ( m, xd )
!
!  As soon as a decrease is noticed, accept the new point.
!
      if ( fxd < fx ) then
        x = xd
        fx = fxd
        decrease = .true.
        exit
      end if

      s = - s
      if ( s == + 1.0D+00 ) then
        i = i + 1
      end if

    end do
!
!  If no decrease occurred, reduce DELTA.
!
    if ( .not. decrease ) then
      delta = delta / 2.0D+00
      if ( delta < delta_tol ) then
        exit
      end if
    end if

  end do

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
