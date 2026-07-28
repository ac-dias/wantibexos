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

subroutine kpath(outputfolder,rlat1,rlat2,rlat3,nks,ks,npts,kpt)

	implicit none

	character(len=70) :: outputfolder    !pasta saida
	integer :: i,j,erro
	integer :: nks,npts
	real,dimension(3) :: rlat1,rlat2,rlat3
	real,dimension(3) :: blat1,blat2,blat3

	real,dimension((nks/2)*npts,4) :: kpt

	real,dimension(nks,3) :: ks
	real,dimension(nks,3) :: ksaux

	real,dimension(nks/2) :: kdis !distancia entre os pontos k do caminho
	real :: kdistot !soma do caminho todo
	real :: kdisaux


	
	OPEN(UNIT=1733, FILE=trim(outputfolder)//'KLABELS.dat',STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening KLABELS output file"


	call recvec(rlat1,rlat2,rlat3,blat1,blat2,blat3)

	do i=1,nks
	
		ksaux(i,1) = ks(i,1)*blat1(1)+ks(i,2)*blat2(1)+ks(i,3)*blat3(1)
		ksaux(i,2) = ks(i,1)*blat1(2)+ks(i,2)*blat2(2)+ks(i,3)*blat3(2)
		ksaux(i,3) = ks(i,1)*blat1(3)+ks(i,2)*blat2(3)+ks(i,3)*blat3(3)
	
	end do



	kdistot= 0.0
	do i=1,nks/2

		call distvec(ksaux(2*i,:),ksaux(2*i-1,:),kdis(i))

		!kdis(i)=sqrt((ksaux(2*i,1)-ksaux(2*i-1,1))**2+(ksaux(2*i,2)-ksaux(2*i-1,2))**2+(ksaux(2*i,3)-ksaux(2*i-1,3))**2)
		kdistot=kdistot+kdis(i)

	end do

	kdis=kdis/kdistot


	kdisaux=0.0
	write(1733,*) real(kdisaux)

	 do j=1,nks/2
		
	  do i=1,npts

		!kpt((j-1)*npts+i) = (kdisaux)+(kdis(j))*dble((i)/(npts))
		kpt((j-1)*npts+i,1) = (kdisaux)+(kdis(j))*(i-1.)/(npts-1.)
		kpt((j-1)*npts+i,2) = ksaux(2*j-1,1)+ (ksaux(2*j,1)-ksaux(2*j-1,1))*(i-1.)/(npts-1.)
		kpt((j-1)*npts+i,3) = ksaux(2*j-1,2)+ (ksaux(2*j,2)-ksaux(2*j-1,2))*(i-1.)/(npts-1.)
		kpt((j-1)*npts+i,4) = ksaux(2*j-1,3)+ (ksaux(2*j,3)-ksaux(2*j-1,3))*(i-1.)/(npts-1.)	
		
		!write(1733,*) kpt((j-1)*npts+i,2),kpt((j-1)*npts+i,3),kpt((j-1)*npts+i,4)
		

	  end do

		kdisaux=kdisaux+kdis(j)
		write(1733,*) real(kdisaux)
	 end do


	close(1733)


end subroutine kpath

subroutine distvec(vec1,vec2,distvecout) !subroutina para calcular a distância entre dois vetores

	implicit none

	real,dimension(3) :: vec1,vec2
	real :: distvecout

	distvecout=sqrt((vec1(1)-vec2(1))**2+(vec1(2)-vec2(2))**2+(vec1(3)-vec2(3))**2)


end subroutine distvec

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

subroutine interpolate_real(qpt,ngrid,nkpt,eqijk,qptmesh,rinterp)

	implicit none
	integer :: i
	integer,dimension(3) :: ngrid
	integer :: nkpt
	
	real :: cardinal_weight
	real,dimension(nkpt) :: eqijk
	real,dimension(3) :: qpt
	real,dimension(nkpt,3) :: qptmesh
	real :: rinterp
	
	real :: lx,ly,lz
	
	rinterp = 0.0
	
	do i=1,nkpt
	
		lx = cardinal_weight(qpt(1),ngrid(1),qptmesh(i,1))
		ly = cardinal_weight(qpt(2),ngrid(2),qptmesh(i,2))
		lz = cardinal_weight(qpt(3),ngrid(3),qptmesh(i,3))
		
		rinterp = rinterp + eqijk(i)*lx*ly*lz
	
	
	end do



end subroutine interpolate_real

function cardinal_weight(qpt,ngrid,qptmesh)

	implicit none
	
	real :: qpt,qptmesh
	integer :: ngrid
	real,parameter :: pi=acos(-1.)
	
	real :: numerator
	real :: denominator
	real :: cardinal_weight
	
	if (ngrid .eq. 1) then
	
		cardinal_weight = 1.000
		
	else
	
		if (abs(qpt-qptmesh) .lt. 1.0E-13 ) then
		
			cardinal_weight = 1.000
		else
		
			numerator = sin(pi*real(ngrid)*(qpt-qptmesh))
			
			if (mod(ngrid,2) .eq. 1) then
			
				denominator = real(ngrid)*sin(pi*(qpt-qptmesh))
			
			else
			
				denominator = real(ngrid)*tan(pi*(qpt-qptmesh))
			
			end if
		
			cardinal_weight = numerator/denominator
		
		
		end if
	
	
	
	end if	



end function cardinal_weight
