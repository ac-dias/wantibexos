subroutine hamiltonian(w90basis,nvec,rvec,hopmatrices,ihopmatrices,ffactor,kx,ky,kz,htb)

	implicit none
	
	integer :: i,j
	integer :: w90basis,nvec
	real :: kx,ky,kz
	real,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices
	real,dimension(nvec,3) :: rvec
	complex,dimension(w90basis,w90basis) :: htb
	complex,dimension(nvec) :: kvecs
	complex,parameter :: imag=cmplx(0.0,1.0)
	integer,dimension(nvec) :: ffactor


		do i=1,nvec

		kvecs(i) = cmplx(0.0,(kx*rvec(i,1))+(ky*rvec(i,2))+(kz*rvec(i,3)))


		end do


		htb = 0.0

		do i=1,nvec

			htb=htb+exp(kvecs(i))*(hopmatrices(i,:,:)+imag*ihopmatrices(i,:,:))*(1.0/real(ffactor(i)))

		end do

end subroutine hamiltonian


subroutine overlap(nbasis,nvec,rvec,s0,kx,ky,kz,stb) !escreve a matriz S(K)

	implicit none
	integer :: nbasis,nvec,i
	real,dimension(nvec,3) :: rvec
	real :: kx,ky,kz
	real,dimension(nvec,nbasis,nbasis) :: s0
	complex,dimension(nbasis,nbasis) :: stb
	complex,parameter :: imag=cmplx(0.0,1.0)
	
	complex,dimension(nvec) :: kvecs
	
	do i=1,nvec

	kvecs(i) = cmplx(0.0,(kx*rvec(i,1))+(ky*rvec(i,2))+(kz*rvec(i,3)))
	

	end do
	
	stb= 0.0
	
	do i=1,nvec

		stb=stb+exp(kvecs(i))*s0(i,:,:)
			
	end do


end subroutine overlap

subroutine fermilvl(w90basis,efermi,fermimt)

	implicit none

	integer :: w90basis,i
	real :: efermi
	complex,dimension(w90basis,w90basis) :: fermimt

	fermimt = cmplx(0.0,0.0)

	do i=1,w90basis

		fermimt(i,i) = cmplx(-efermi,0.0)

	end do


end subroutine fermilvl

subroutine spin_exchange(w90basis,dft,systype,exc,mag,ovptb,hexc)

	implicit none
	
	integer :: w90basis
	character(len=1) :: dft
	character(len=4) :: systype
	
	real :: exc
	real,dimension(3) :: mag
	
	complex,dimension(w90basis,w90basis) :: spmz,spmx,spmy,hexc,ovptb
	
	
	if (systype .eq. "NP") then
	
	hexc = 0.0
	
	else
	
	call spin_z_matrix(w90basis,dft,systype,ovptb,spmz)
	call spin_x_matrix(w90basis,dft,systype,ovptb,spmx)		
	call spin_y_matrix(w90basis,dft,systype,ovptb,spmy)
	
	spmz = spmz*mag(3)*exc
	spmx = spmx*mag(1)*exc
	spmy = spmy*mag(2)*exc
	
	hexc= spmx+spmy+spmz

	!hexc = 0.0

	end if

end subroutine spin_exchange





