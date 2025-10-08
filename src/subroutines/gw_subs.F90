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

subroutine polarization(fermishift,sysdim,ngrid,nk,rlat,w,nomega,w90basis,systype,eta,ek,ekq,ovpmnkkq,pol) !given polarization for a given q for all w-freq

	implicit none
	real :: fermishift
	character(len=5) :: sysdim
	integer,dimension(3) :: ngrid
	integer :: nk 
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
	
	!select case (systype)
	
	!case("NP")
	
	!	gs = 2.0
	
	!case default
	
		gs = 1.0	
	!end select
	
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

		vq= v2dkgw(qpt,v0,ediel,rlat,ngrid,lc,tolr)

	case("V3D")

		vq= vcoulgw(qpt,v0,rlat,ngrid,tolr)

	case("V3DL")

		vq= v3dielgw(qpt,v0,ediel,rlat,ngrid,tolr)
		
	case("V2D")

		vq= v2dgw(qpt,v0,rlat,ngrid,tolr)

	case("V2DL")

		vq= v2dielgw(qpt,v0,ediel,rlat,ngrid,tolr)		

	case("V2DT")

		vq= v2dtgw(qpt,v0,ngrid,rlat,tolr)

	case("V2DT2")

		vq= v2dt2gw(qpt,v0,ngrid,rlat,lc,tolr)
		
	case("V2DOH")

		vq= v2dohonogw(qpt,v0,ngrid,rlat,ediel,w,ez,tolr)
		
	case("V2DRK")

		vq= v2drkgw(qpt,v0,ngrid,rlat,ediel,lc,ez,w,r0,tolr)
		
	case("V1D")
	
		vq= v1dgw(qpt,v0,ngrid,rlat,lc,tolr)
		
		
	case("V1DL")
	
		vq= v1dielgw(qpt,v0,ngrid,rlat,lc,ediel,tolr)
						
		
	case("V1DT")
	
		vq= v1dtgw(qpt,v0,ngrid,rlat,tolr,lc)	
				

	case("V0DT")

		vq= v0dtgw(qpt,v0,ngrid,rlat,tolr)

	case default

		write(*,*) "Wrong Coulomb Potential"
		STOP

	end select
	
	do i=1,nomega
	
		w0(i) = (vq)/(1.0-(vq*pol(i)))
	
	end do


end subroutine w0coul

subroutine self_en_x(w90basis,sysdim,rlat,n,kpt,ngrid,nq,qpt,nocpq,coultype,ediel,lc,ez,w,r0,tolr,ovpmnkq,sxnk)

	implicit none
	
	integer :: i,j
	
	integer :: w90basis
	character(len=5) :: sysdim
	real,dimension(3,3) :: rlat
	integer :: n
	real,dimension(3) :: kpt
	integer,dimension(3) :: ngrid
	
	integer :: nq!=ngrid(1)*ngrid(2)*ngrid(3)
	real,dimension(nq,3) :: qpt
	integer,dimension(nq) :: nocpq
	
	character(len=5) :: coultype
	real,dimension(3) :: ediel
	real :: lc,ez,w,r0,tolr
	
	real :: v2dkgw,vcoulgw,v3dielgw,v2dgw,v2dielgw,v2dtgw
	real :: v0dtgw,v2dt2gw,v2drkgw,v2dohonogw,v1dtgw,v1dgw,v1dielgw
	real :: vq
	
	real,dimension(w90basis,w90basis,nq) :: ovpmnkq
	real :: wk,waux,sxnk
	
	sxnk = 0.0
	
	waux = wk(sysdim,ngrid,rlat)
	
	do i=1,nq
	
	select case (coultype)

	case("V2DK")

		vq= v2dkgw(kpt,qpt(i,:),ediel,rlat,ngrid,lc,tolr)

	case("V3D")

		vq= vcoulgw(kpt,qpt(i,:),rlat,ngrid,tolr)

	case("V3DL")

		vq= v3dielgw(kpt,qpt(i,:),ediel,rlat,ngrid,tolr)
		
	case("V2D")

		vq= v2dgw(kpt,qpt(i,:),rlat,ngrid,tolr)

	case("V2DL")

		vq= v2dielgw(kpt,qpt(i,:),ediel,rlat,ngrid,tolr)		

	case("V2DT")

		vq= v2dtgw(kpt,qpt(i,:),ngrid,rlat,tolr)

	case("V2DT2")

		vq= v2dt2gw(kpt,qpt(i,:),ngrid,rlat,lc,tolr)
		
	case("V2DOH")

		vq= v2dohonogw(kpt,qpt(i,:),ngrid,rlat,ediel,w,ez,tolr)
		
	case("V2DRK")

		vq= v2drkgw(kpt,qpt(i,:),ngrid,rlat,ediel,lc,ez,w,r0,tolr)
		
	case("V1D")
	
		vq= v1dgw(kpt,qpt(i,:),ngrid,rlat,lc,tolr)
		
		
	case("V1DL")
	
		vq= v1dielgw(kpt,qpt(i,:),ngrid,rlat,lc,ediel,tolr)		
		
	case("V1DT")
	
		vq= v1dtgw(kpt,qpt(i,:),ngrid,rlat,tolr,lc)	
				

	case("V0DT")

		vq= v0dtgw(kpt,qpt(i,:),ngrid,rlat,tolr)

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

subroutine self_c_znk(w90basis,sysdim,rlat,ngrid,nq,fermishift,n,enk,emkmq,kpt,omega,nomega,w0,ovpmnkkmq,eta,scnk,znk)

	implicit none
	
	integer :: i,j,k
	real,parameter :: pi=acos(-1.)

	real :: fermishift
	integer :: w90basis
	character(len=5) :: sysdim
	real,dimension(3,3) :: rlat
	real,dimension(3) :: kpt

	integer,dimension(3) :: ngrid
	integer :: nq!=ngrid(1)*ngrid(2)*ngrid(3)
	!real,dimension(nq,3) :: qpt
	
	integer :: n
	real :: enk
	integer :: nomega
	real,dimension(nomega) :: omega
	real :: domega
	
	complex,dimension(nq,nomega) :: w0
	real,dimension(w90basis,w90basis,nq) :: ovpmnkkmq
	real,dimension(nq,w90basis) :: emkmq
	
	real :: wk,waux
	
	real :: eta
	complex :: scnk
	real :: dscnk, znk
	
	real :: auxself2
	complex :: auxself1
	
	real :: integral2,aux3,aux4
	complex :: integral1,aux1,aux2
	
	waux = wk(sysdim,ngrid,rlat)
	
	scnk = 0.0
	dscnk = 0.0
	
	domega = (omega(2)-omega(1))/pi
	
	do i=1,nq
	  do j=1,w90basis
	
	
	    integral1 = 0.0
	    integral2 = 0.0
	    !1/2 simpson integral
	    do k=1,nomega-1
	      
	      aux1 = aimag(w0(i,k))*auxself1(fermishift,enk,emkmq(i,j),omega(k),eta) 
	      aux2 = aimag(w0(i,k+1))*auxself1(fermishift,enk,emkmq(i,j),omega(k+1),eta)
	      
	      integral1 = integral1 + ((aux1+aux2)/2.0)*domega
	      
	      aux3 = aimag(w0(i,k))*auxself2(fermishift,enk,emkmq(i,j),omega(k),eta)
	      aux4 = aimag(w0(i,k+1))*auxself2(fermishift,enk,emkmq(i,j),omega(k+1),eta)
	      
	      integral2 = integral2 + ((aux3+aux4)/2.0)*domega
	
	    end do
	    
	      scnk = scnk + (integral1*ovpmnkkmq(n,j,i)*waux)
	      
	      dscnk = dscnk + (integral2*ovpmnkkmq(n,j,i)*waux)
	    
	  end do
	end do
	
	znk = 1.0/(1.0 - dscnk)


end subroutine self_c_znk




