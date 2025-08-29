function enocp(fermishift,en)

	implicit none
	
	real :: enocp,fermishift,en
	
	if (en .gt. fermishift) then
	
		enocp = 0.0
	
	else
	
		enocp = 1.0
	
	end if


end function

subroutine ovpsqr(dft,w90basis,wf1,sk1,wf2,sk2,ovpmn)

	implicit none
	character(len=1) :: dft
	integer :: w90basis
	real,dimension(w90basis,w90basis) :: sk1,sk2,skaux
	complex,dimension(w90basis) :: wf1,wf2,wf1c
	real :: ovpmn
	complex :: vc
	
	select case (dft)
	
	case("S")
	
	call vecconjg(wf1,w90basis,wf1c)
	
	skaux = 0.5*(sk1+sk2)
	
	call sandwich(w90basis,wf1c,skaux,wf2,vc)
	
	ovpmn = vc*conjg(vc)
	
	case default
	
	call vecconjg(wf1,w90basis,wf1c)
	call prodintsq(wf1c,wf2,w90basis,vc)
	
	ovpmn = vc*conjg(vc)
	
	end select


end subroutine ovpsqr

subroutine polarization(fermishift,sysdim,ngrid,rlat,w,nomega,w90basis,systype,eta,ek,ekq,ovpmnkkq,pol) !given polarization for a given q for all w-freq

	implicit none
	real :: fermishift
	character(len=5) :: sysdim
	integer,dimension(3) :: ngrid
	integer :: nk = ngrid(1)*ngrid(2)*ngrid(3)
	real,parameter:: pi=acos(-1.)
	real,dimension(3,3) :: rlat	
	
	integer :: nomega
	real,dimension(nomega) :: w
	
	integer :: w90basis
	character(len=4) :: systype
	real :: eta
	real,dimension(w90basis,nk) :: ek,ekq
	real,dimension(w90basis,w90basis,nk) :: ovpmnkkq
	
	complex,dimension(nomega) :: pol
	
	real :: gs,wk,waux,enocp
	
	real :: aux1
	complex :: aux2
	
	integer :: i,j,k,l
	
	waux = wk(sysdim,ngrid,rlat)
	
	select case (systype)
	
	case("NP")
	
		gs = 2.0
	
	case default
	
		gs = 1.0	
	end select
	
	pol = 0.0
	
	do l=1,nomega
	  do i=1,nk
	    do j=1,w90basis
	      do k=1,w90basis
	    
	    	aux1 = enocp(fermishift,ek(j,i)) - enocp(fermishift,ekq(k,i))
	    	
	    	aux2 = w(l) - (ekq(k,i)-ek(j,i)) + eta*cmplx(0.0,1.0)
	    	
	    	pol(l) = pol(l) + ((aux1/aux2)*ovpmnkkq(j,k,i))
	    
	      end do
	    end do
	  end do
	end do
	
		pol = pol*gs*waux
	
end subroutine polarization

subroutine w0coul(qpt,nomega,rlat,ngrid,nk,coultype,ediel,lc,ez,w,r0,tolr,w0)  !given W0 (Coulomb interaction) for a given q for all w-freq

	implicit none
	integer :: nomega,i
	real,dimension(3,3) :: rlat
	integer,dimension(3) :: ngrid
	integer :: nk
	real,dimension(3) :: qpt,v0
	character(len=5) :: coultype
	complex,dimension(nomega) :: pol,w0
	real :: vq
	real :: tolr
	real,dimension(3) :: ediel
	real :: lc,ez,w,r0
	
	real :: v2dkgw,vcoulgw,v3dielgw,v2dgw,v2dielgw,v2dtgw
	real :: v0dtgw,v2dt2gw,v2drkgw,v2dohonogw,v1dtgw,v1dgw,v1dielgw
	
	v0 = 0.0
	
	select case (coultype)

	case("V2DK")

		vq= v2dk(qpt,v0,ediel,rlat,ngrid,lc,tolr)

	case("V3D")

		vq= vcoul(qpt,v0,rlat,ngrid,tolr)

	case("V3DL")

		vq= v3diel(qpt,v0,ediel,rlat,ngrid,tolr)
		
	case("V2D")

		vq= v2d(qpt,v0,rlat,ngrid,tolr)

	case("V2DL")

		vq= v2diel(qpt,v0,ediel,rlat,ngrid,tolr)		

	case("V2DT")

		vq= v2dt(qpt,v0,ngrid,rlat,tolr)

	case("V2DT2")

		vq= v2dt2(qpt,v0,ngrid,rlat,lc,tolr)
		
	case("V2DOH")

		vq= v2dohono(qpt,v0,ngrid,rlat,ediel,w,ez,tolr)
		
	case("V2DRK")

		vq= v2drk(qpt,v0,ngrid,rlat,ediel,lc,ez,w,r0,tolr)
		
	case("V1D")
	
		vq= v1dgw(qpt,v0,ngrid,rlat,lc,tolr)
		
		
	case("V1DL")
	
		vq= v1dielgw(qpt,v0,ngrid,rlat,lc,ediel,tolr)
						
		
	case("V1DT")
	
		vq= v1dt(qpt,v0,ngrid,rlat,tolr,lc)	
				

	case("V0DT")

		vq= v0dt(qpt,v0,ngrid,rlat,tolr)

	case default

		write(*,*) "Wrong Coulomb Potential"
		STOP

	end select
	
	do i=1,nomega
	
		w0(i) = (vq)/(1.0-(vq*pol(i)))
	
	end do


end subroutine w0

subroutine self_en_x(sysdim,rlat,n,kpt,ngrid,qpt,nocpq,coultype,ediel,lc,ez,w,r0,tolr,sxnk)

	implicit none
	character(len=5) :: sysdim
	real,dimension(3,3) :: rlat
	integer :: n
	real,dimension(3) :: kpt
	integer,dimension(3) :: ngrid
	
	integer :: nq=ngrid(1)*ngrid(2)*ngrid(3)
	real,dimension(nk,3) :: qpt
	integer,dimension(nk) :: nocpq
	
	character(len=5) :: coultype
	real,dimension(3) :: ediel
	real :: lc,ez,w,r0,tolr
	
	real :: v2dkgw,vcoulgw,v3dielgw,v2dgw,v2dielgw,v2dtgw
	real :: v0dtgw,v2dt2gw,v2drkgw,v2dohonogw,v1dtgw,v1dgw,v1dielgw
	real :: vq
	
	real,dimension(w90basis,w90basis,nq) :: ovpmnkq
	real :: waux,sxnk
	
	sxnk = 0.0
	
	waux = wk(sysdim,ngrid,rlat)
	
	do i=1,nq
	
	select case (coultype)

	case("V2DK")

		vq= v2dk(kpt,qpt(i,:),ediel,rlat,ngrid,lc,tolr)

	case("V3D")

		vq= vcoul(kpt,qpt(i,:),rlat,ngrid,tolr)

	case("V3DL")

		vq= v3diel(kpt,qpt(i,:),ediel,rlat,ngrid,tolr)
		
	case("V2D")

		vq= v2d(kpt,qpt(i,:),rlat,ngrid,tolr)

	case("V2DL")

		vq= v2diel(kpt,qpt(i,:),ediel,rlat,ngrid,tolr)		

	case("V2DT")

		vq= v2dt(kpt,qpt(i,:),ngrid,rlat,tolr)

	case("V2DT2")

		vq= v2dt2(kpt,qpt(i,:),ngrid,rlat,lc,tolr)
		
	case("V2DOH")

		vq= v2dohono(kpt,qpt(i,:),ngrid,rlat,ediel,w,ez,tolr)
		
	case("V2DRK")

		vq= v2drk(kpt,qpt(i,:),ngrid,rlat,ediel,lc,ez,w,r0,tolr)
		
	case("V1D")
	
		vq= v1dgw(kpt,qpt(i,:),ngrid,rlat,lc,tolr)
		
		
	case("V1DL")
	
		vq= v1dielgw(kpt,qpt(i,:),ngrid,rlat,lc,ediel,tolr)		
		
	case("V1DT")
	
		vq= v1dt(kpt,qpt(i,:),ngrid,rlat,tolr,lc)	
				

	case("V0DT")

		vq= v0dt(kpt,qpt(i,:),ngrid,rlat,tolr)

	case default

		write(*,*) "Wrong Coulomb Potential"
		STOP

	end select
	
	  do j=1,nocpq(i)
	  
	  	sxnk = sxnk - (waux*vq*ovpmnkq(n,j,i))
	  
	  end do
	
	
	end do	
	
	

end subroutine self_en_x

function auxself1(fermishift,en,em,omega,eta)

	implicit none
	complex :: auxself1
	real :: fermishift,en,em,omega,eta
	
	real :: fmkq
	complex :: aux1,aux2
	
	if (em .gt. fermishift) then
	
		fmkq = 0.0
	
	else 
	
		fmkq = 1.0
	end if
	
	aux1 = (1.0 - fmkq)/(en-(em+omega)+eta*cmplx(0.0,1.0))
	aux2 = (fmkq)/(en-(em-omega)+eta*cmplx(0.0,1.0))	

	auxself1 = aux1+aux2

end function

function auxself2(fermishift,en,em,omega,eta)

	implicit none
	real :: auxself2
	real :: fermishift,en,em,omega,eta
	
	real :: fmkq
	real :: xp,xm,etasqr
	real :: aux1,aux2,aux3,aux4
	
	if (em .gt. fermishift) then
	
		fmkq = 0.0
	
	else 
	
		fmkq = 1.0
	end if
	
	etasqr = eta*eta
	xp = en - em + omega		
	xm = en - em - omega
	
	aux1 = ((xm*xm)+(etasqr))*((xm*xm)+(etasqr))
	aux2 = ((xp*xp)+(etasqr))*((xp*xp)+(etasqr))
	
	aux3 = (1.0-fmkq)*(etasqr-(xm*xm))
	aux4 = (fmkq)*(etasqr-(xp*xp))
	
	auxself2 = (aux3/aux1) + (aux4/aux2)	

end function



