subroutine transmittance(ni,ns,tmax,optcond,tr)

	implicit none
	
	real :: ni,ns,tmax,tr
	complex :: optcond,traux,aux2
	
	real,parameter:: pi=acos(-1.)
	real,parameter :: clight=3.0E+08
	
	aux2 = (4.0*pi)*tmax*optcond/clight
	
	traux = (2.0)/(ni+ns+aux2)
	
	tr = ns*traux*conjg(traux)

end subroutine transmittance

subroutine reflectance(ni,ns,tmax,optcond,ref)

	real :: ni,ns,tmax,ref
	complex :: optcond,refaux,aux2
	
	real,parameter:: pi=acos(-1.)
	real,parameter :: clight=3.0E+08
	
	aux2 = (4.0*pi)*tmax*optcond/clight
	
	refaux = (ni-ns-aux2)/(1.0+ns+aux2)
	
	ref = refaux*conjg(refaux)

end subroutine reflectance




subroutine optcondcalc(e1,e2,efoton,optcond)

	implicit none
	
	complex,parameter :: imag=cmplx(0.0,-1.0)
	real :: e1,e2
	complex :: edielec,optcond
	real,parameter :: hbar = 6.582119569E-16 !hbar planck's constant
	real,parameter:: pi=acos(-1.)
	real :: efoton,omega,aux
	
	edielec = cmplx(e1,e2)
	
	omega = efoton/hbar
	
	aux = omega/(4.0*pi)
	
	optcond = imag*aux*(edielec-1.0)


end subroutine optcondcalc


subroutine refracao(e1,e2,refrac)

	implicit none

	real :: e1,e2
	real :: aux,refrac

	aux = sqrt((e1*e1)+(e2*e2))

	refrac = sqrt((aux+e1)/2.0)

end subroutine refracao

subroutine extincao(e1,e2,extinc)

	implicit none

	real :: e1,e2
	real :: aux,extinc

	aux = sqrt((e1*e1)+(e2*e2))

	extinc = sqrt((aux-e1)/2.0)

end subroutine extincao

subroutine reflectibilidade(e1,e2,reflec)

	implicit none

	real :: e1,e2
	real :: reflec,refrac,extinc
	real :: aux1,aux2,aux3

	call extincao(e1,e2,extinc)
	call refracao(e1,e2,refrac)

	aux1 = (refrac-1.0)*(refrac-1.0)
	aux2 = (refrac+1.0)*(refrac+1.0)
	aux3 = extinc*extinc

	reflec = (aux1+aux3)/(aux2+aux3)
	

end subroutine reflectibilidade

subroutine abscoef(e1,e2,efoton,absc)

	implicit none

	real :: e1,e2
	real :: absc
	real :: aux,efoton
	real,parameter :: hc=19.746E-06

	aux = sqrt((e1*e1)+(e2*e2))
	absc=(sqrt(2.0)*efoton)*(1.0/hc)*sqrt(aux-e1)

end subroutine abscoef

subroutine enloss(e1,e2,eloss)

	implicit none

	real :: e1,e2
	real :: eloss
	real :: aux
	
	aux = (e1*e1)+(e2*e2)
	
	eloss= e2/aux
	
end subroutine enloss




