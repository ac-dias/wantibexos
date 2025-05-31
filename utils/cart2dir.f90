subroutine matvec(matriz,vetor,n,veout)

	implicit none

	integer :: n !ordem da matriz e tamanho do vetor
	complex, dimension(n,n) :: matriz
	complex, dimension(n):: vetor,veout
	complex :: flag
	integer :: i,j

	veout=0
	
	do i=1,n


		do j=1,n

			flag=matriz(i,j)*vetor(j)

			veout(i)=flag+veout(i)
		

		end do
	

	end do


end subroutine matvec

SUBROUTINE invers(h11inv,n1)

	integer n1
	complex h11inv(n1,n1)
        integer    lda1,lwk1,if1,ipiv1(n1)
        complex wk1(64*n1)
	lda1     =n1
	lwk1     =64*n1

        call zgetrf(n1,n1,h11inv,lda1,ipiv1,if1)
  	if (if1.eq.0) then
	call zgetri(n1,h11inv,lda1,ipiv1,wk1,lwk1,if1)
	else
	write(*,*)'error in routine invers'
	endif

   

END SUBROUTINE invers

subroutine prodvec(v1,v2,vx) !calcula o vetor oriundo do produto vetorial entre dois vetores v1 X v2

	implicit none
	
	real,dimension(3) :: v1,v2,vx

	vx(1) = (v1(2)*v2(3))-(v1(3)*v2(2))
	vx(2) = (v1(3)*v2(1))-(v1(1)*v2(3))
	vx(3) = (v1(1)*v2(2))-(v1(2)*v2(1))


end subroutine prodvec

subroutine recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3) !calcula os vetores da rede recíproca, a partir dos vetores da rede real

	implicit none
	
	real,parameter :: pi=acos(-1.)
	
	real,dimension(3) :: rlat1,rlat2,rlat3
	real,dimension(3) :: blat1,blat2,blat3

	real,dimension(3) :: v23,v31,v12
	real :: vol

	call prodvec(rlat2,rlat3,v23)
	call prodvec(rlat3,rlat1,v31)
	call prodvec(rlat1,rlat2,v12)

	vol= abs((rlat1(1)*v23(1))+(rlat1(2)*v23(2))+(rlat1(3)*v23(3)))

	blat1 = ((2.0*pi)/vol)*v23
	blat2 = ((2.0*pi)/vol)*v31
	blat3 = ((2.0*pi)/vol)*v12


end subroutine recvec

subroutine conv_kpt_d2c(rlat1,rlat2,rlat3,kpt)!convert kpoint coordinates from direct to cartesian

	implicit none
	
	real,parameter :: pi=acos(-1.)
	real,dimension(3) :: rlat1,rlat2,rlat3
	real,dimension(3) :: blat1,blat2,blat3
	real,dimension(3) :: kpt,kptaux
	
	
	call recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3)
	
	kptaux(1) = kpt(1)*blat1(1)+kpt(2)*blat2(1)+kpt(3)*blat3(1)
	kptaux(2) = kpt(1)*blat1(2)+kpt(2)*blat2(2)+kpt(3)*blat3(2)
	kptaux(3) = kpt(1)*blat1(3)+kpt(2)*blat2(3)+kpt(3)*blat3(3)
	
	kpt = kptaux	 	

end subroutine conv_kpt_d2c
