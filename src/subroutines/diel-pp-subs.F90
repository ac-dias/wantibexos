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




