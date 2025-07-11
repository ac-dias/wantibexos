program main

	implicit none 
	
	integer :: i,nk,flagi
	real,dimension(3,3) :: cvec,rvec
	real,allocatable,dimension(:,:) :: kcar,kdir,kcar2
	
	character(len=70) :: flagx
	real :: flag
	
	character(len=70) :: a,b
	character(len=70) :: add
	integer :: erro
	character(len=70) :: emfile
	character(len=70) :: params	
	
	!lendo input wantibexos
	do 

	read(*,*,iostat=erro) a,b
	if (erro/=0) exit

	select case (a)
	
	case ("EM_FILE=")

	read(b,*) emfile 
	
	case ("PARAMS_FILE=")

		params = b	
	
	case default
	 continue

	end select


	end do
	
	OPEN(UNIT=100, FILE= params,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening hamiltonian hr input file"
	OPEN(UNIT=200, FILE= emfile,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening effective mass input file"  
	OPEN(UNIT=300, FILE= "results.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening results output file"  
  	    	
    	
   					
	
	!lendo vetores da rede no hamiltoniano do input
	read(100,*) flagx
	read(100,*) flag
	read(100,*) flag
			
	do i=1,3
	
		read(100,*) cvec(i,1),cvec(i,2),cvec(i,3)
	
	end do
	
	!vetores da rede reciproca
	call recvec(cvec(1,:),cvec(2,:),cvec(3,:),rvec(1,:),rvec(2,:),rvec(3,:))	
	
	!lendo os pontos k para serem transformados
	read(200,*) nk
	read(200,*) flagx
	
	!nk=1
	allocate(kcar(nk,3),kdir(nk,3),kcar2(nk,3))
	
	!kcar(1,1) = -0.27564469E+00 
	!kcar(1,2) = -0.11937930E+00
	!kcar(1,3) = 0.24334151E+00		
	
	
	do i=1,nk
		read(200,*) kcar(i,1),kcar(i,2),kcar(i,3),flagi,flagx
		
		write(300,*) "nk=",i,"kptcar",kcar(i,1),kcar(i,2),kcar(i,3)
	end do

		write(300,*) "##########################"
	
	do i=1,nk
	
		call kptcar2kptdir(kcar(i,:),rvec,kdir(i,:))
		
		write(300,*)"nk=",i,"kptdir",kdir(i,1),kdir(i,2),kdir(i,3)
		
	end do	

		write(300,*) "##########################"

	do i=1,nk
		
		call conv_kpt_d2c(cvec(1,:),cvec(2,:),cvec(3,:),kdir(i,:),kcar2(i,:))
		
		write(300,*)"nk=",i,"kptcar",kcar2(i,1),kcar2(i,2),kcar2(i,3)	
	end do
	
	deallocate(kcar,kdir,kcar2)
	
	close(100)
	close(200)
	close(300)



end program main


subroutine kptcar2kptdir(kcar,rvec,kdir) !converte o ponto k cartesiano para k em coordenadas diretas

	implicit none
	
	real,dimension(3,3) :: rvec,c2dm
	real,dimension(3) :: kcar,kdir
	
	call c2dmatrix(rvec,c2dm)
	
	call matvec(c2dm,kcar,3,kdir)


end subroutine


subroutine c2dmatrix(rvec,c2dm) !matriz que converte coordenadas reais para coordenadas diretas


	implicit none
	
	real,dimension(3,3) :: rvec
	real,dimension(3,3) :: c2dm,c2dmaux
	
	c2dm(1,1) = rvec(1,1)
	c2dm(1,2) = rvec(2,1)
	c2dm(1,3) = rvec(3,1)
	
	c2dm(2,1) = rvec(1,2)
	c2dm(2,2) = rvec(2,2)
	c2dm(2,3) = rvec(3,2)
	
	c2dm(3,1) = rvec(1,3)
	c2dm(3,2) = rvec(2,3)
	c2dm(3,3) = rvec(3,3)
	
	c2dmaux = c2dm
	
	call invers(c2dm,3)		

	!write(*,*) matmul(c2dm,c2dmaux)

end subroutine



subroutine matvec(matriz,vetor,n,veout)

	implicit none

	integer :: n !ordem da matriz e tamanho do vetor
	real, dimension(n,n) :: matriz
	real, dimension(n):: vetor,veout
	real :: flag
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
	real h11inv(n1,n1)
        integer    lda1,lwk1,if1,ipiv1(n1)
        real wk1(64*n1)
	lda1     =n1
	lwk1     =64*n1

        call sgetrf(n1,n1,h11inv,lda1,ipiv1,if1)
  	if (if1.eq.0) then
	call sgetri(n1,h11inv,lda1,ipiv1,wk1,lwk1,if1)
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

subroutine conv_kpt_d2c(rlat1,rlat2,rlat3,kpt,kptout)!convert kpoint coordinates from direct to cartesian

	implicit none
	
	real,parameter :: pi=acos(-1.)
	real,dimension(3) :: rlat1,rlat2,rlat3
	real,dimension(3) :: blat1,blat2,blat3
	real,dimension(3) :: kpt,kptout
	
	
	call recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3)
	
	kptout(1) = kpt(1)*blat1(1)+kpt(2)*blat2(1)+kpt(3)*blat3(1)
	kptout(2) = kpt(1)*blat1(2)+kpt(2)*blat2(2)+kpt(3)*blat3(2)
	kptout(3) = kpt(1)*blat1(3)+kpt(2)*blat2(3)+kpt(3)*blat3(3)
	
	!kpt = kptaux	 	

end subroutine conv_kpt_d2c
