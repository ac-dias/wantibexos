subroutine berryct2(nocp,kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		  ihopmatrices,efermi,energias,autovetores,nthread,bxy,bxz,byz)
		   

	implicit none

	integer :: nthread

	real :: kx,ky,kz,ezz
	integer :: nocp,w90basis
	integer :: nvec
	integer :: n,i,j

	real,parameter :: ktol = 0.001

	integer,dimension(nvec) :: ffactor

	real,dimension(3,3) :: rlat
	real,dimension(nvec,3) :: rvec

	integer :: nsite

	complex, dimension(w90basis,w90basis) :: hx,hy,hz
	complex, dimension(w90basis,w90basis) :: autovetores,autovetores2

	real,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices

	real, dimension(w90basis) :: energias,energias2

	!real,dimension(3) :: mag

	real :: efermi

	complex, dimension(w90basis) :: vflagx,vflagy,vflagz
	complex, dimension(w90basis) :: vn,vm,vl

	complex,dimension(w90basis) :: vx,vy,vz

	complex,dimension(w90basis) :: vxa,vya
	complex, dimension(w90basis) :: vflagxa,vflagya
	complex, dimension(w90basis) :: vna,vma

	complex :: bxy,bxz,byz
	complex :: bfxy,bfxz,bfyz

	real:: gammas

	real :: berryaux


	
	call hlm(kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		 ihopmatrices,hx,hy,hz)
	


	bxy=cmplx(0.0,0.0)
	bxz=cmplx(0.0,0.0)
	byz=cmplx(0.0,0.0)
	

	do n=1,nocp,1
	

	do i=nocp+1,w90basis,1

	call vecconjg(autovetores(n,:),w90basis,vn)

	call matvec(hx,autovetores(i,:),w90basis,vflagx)

	call prodintsq(vn,vflagx,w90basis,vx(i))

	!write(*,*) vx(i)


	call vecconjg(autovetores(n,:),w90basis,vm)

	call matvec(hy,autovetores(i,:),w90basis,vflagy)

	call prodintsq(vm,vflagy,w90basis,vy(i))
	

	!write(*,*) vy(i)
	
	call vecconjg(autovetores(n,:),w90basis,vl)

	call matvec(hz,autovetores(i,:),w90basis,vflagz)

	call prodintsq(vl,vflagz,w90basis,vz(i))
	

	if (i == n) then


	bfxy=0.0
	bfxz=0.0
	bfyz=0.0
			

	else if (abs(1.0/(energias(i)-energias(n))) .eq. (1.0/0.0) ) then


	bfxy=0.0
	bfxz=0.0
	bfyz=0.0	

	else
	
	berryaux= (energias(n)-energias(i))**2
	

	bfxy=(((vx(i))*conjg(vy(i)))/(berryaux))
	bfxz=(((vx(i))*conjg(vz(i)))/(berryaux))
	bfyz=(((vy(i))*conjg(vz(i)))/(berryaux))

						

	end if


	bxy=bxy+bfxy
	bxz=bxz+bfxz
	byz=byz+bfyz
					

	end do

	end do



end subroutine berryct2


subroutine berryct(nocp,kx,ky,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		  ihopmatrices,efermi,energias,autovetores,nthread,gammas,berry)
		   

	implicit none

	integer :: nthread

	real :: kx,ky,kz,ezz
	integer :: nocp,w90basis
	integer :: nvec
	integer :: n,i,j

	real,parameter :: ktol = 0.001

	integer,dimension(nvec) :: ffactor

	real,dimension(3,3) :: rlat
	real,dimension(nvec,3) :: rvec

	integer :: nsite

	complex, dimension(w90basis,w90basis) :: hx,hy,hz
	complex, dimension(w90basis,w90basis) :: autovetores,autovetores2

	real,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices

	real, dimension(w90basis) :: energias,energias2

	!real,dimension(3) :: mag

	real :: efermi

	complex, dimension(w90basis) :: vflagx,vflagy
	complex, dimension(w90basis) :: vn,vm

	complex,dimension(w90basis) :: vx,vy

	complex,dimension(w90basis) :: vxa,vya
	complex, dimension(w90basis) :: vflagxa,vflagya
	complex, dimension(w90basis) :: vna,vma

	complex :: berry,berryf

	real:: gammas

	real :: berryaux


	
	call hlm(kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		 ihopmatrices,hx,hy,hz)
	

	berry=cmplx(0.,0.)

	do n=1,nocp,1
	

	do i=nocp+1,w90basis,1

	call vecconjg(autovetores(n,:),w90basis,vn)

	call matvec(hx,autovetores(i,:),w90basis,vflagx)

	call prodintsq(vn,vflagx,w90basis,vx(i))

	!write(*,*) vx(i)


	call vecconjg(autovetores(i,:),w90basis,vm)

	call matvec(hy,autovetores(n,:),w90basis,vflagy)

	call prodintsq(vm,vflagy,w90basis,vy(i))

	!write(*,*) vy(i)

	if (i == n) then

	berryf=0.


	else if (abs(energias(i)-energias(n)) .lt. 0.0001 ) then

	
	!call eigsys(nthread,kx+ktol,ky+ktol,w90basis,nvec,rlat1,rlat2,rlat3,rvec,hopmatrices,&
	!	   ihopmatrices,efermi,energias2,autovetores2)

	!call vecconjg(autovetores2(n,:),w90basis,vna)

	!call matvec(hx,autovetores2(i,:),w90basis,vflagxa)

	!call prodintsq(vna,vflagxa,w90basis,vxa(i))


	!call vecconjg(autovetores2(i,:),w90basis,vma)

	!call matvec(hy,autovetores2(n,:),w90basis,vflagya)

	!call prodintsq(vma,vflagya,w90basis,vya(i))


	!berryf=-2.*(((vxa(i))*vya(i))/((energias2(n)-energias2(i))**2+(gammas)**2))

	berryf=0.

	else
	
	berryaux= (energias(n)-energias(i))**2+(gammas)**2
	berryf=-2.*(((vx(i))*vy(i))/(berryaux))

	end if

	berry=berry+berryf

	end do

	end do



end subroutine berryct


subroutine berryc(n,nocp,kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		  ihopmatrices,efermi,energias,autovetores,nthread,gammas,berry)
		   

	implicit none

	integer :: nthread

	real :: kx,ky,kz,ezz
	integer :: nocp,w90basis
	integer :: nvec
	integer :: n,i,j

	real,parameter :: ktol = 0.001


	real,dimension(3,3) :: rlat
	real,dimension(nvec,3) :: rvec
	integer,dimension(nvec) :: ffactor
	integer :: nsite

	complex, dimension(w90basis,w90basis) :: hx,hy,hz
	complex, dimension(w90basis,w90basis) :: autovetores,autovetores2

	real,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices

	real, dimension(w90basis) :: energias,energias2

	!real,dimension(3) :: mag

	real :: efermi

	complex, dimension(w90basis) :: vflagx,vflagy
	complex, dimension(w90basis) :: vn,vm

	complex,dimension(w90basis) :: vx,vy

	complex,dimension(w90basis) :: vxa,vya
	complex, dimension(w90basis) :: vflagxa,vflagya
	complex, dimension(w90basis) :: vna,vma

	complex :: berry,berryf

	real:: gammas


	
	call hlm(kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		 ihopmatrices,hx,hy,hz)
	

	berry=cmplx(0.,0.)
	

	do i=nocp+1,w90basis,1

	call vecconjg(autovetores(n,:),w90basis,vn)

	call matvec(hx,autovetores(i,:),w90basis,vflagx)

	call prodintsq(vn,vflagx,w90basis,vx(i))

	!write(*,*) vx(i)


	call vecconjg(autovetores(i,:),w90basis,vm)

	call matvec(hy,autovetores(n,:),w90basis,vflagy)

	call prodintsq(vm,vflagy,w90basis,vy(i))

	!write(*,*) vy(i)

	if (i == n) then

	berryf=0.


	else if (abs(energias(i)-energias(n)) .lt. 0.0001 ) then

	
	!call eigsys(nthread,kx+ktol,ky+ktol,w90basis,nvec,rlat1,rlat2,rlat3,rvec,hopmatrices,&
	!	   ihopmatrices,efermi,energias2,autovetores2)

	!call vecconjg(autovetores2(n,:),w90basis,vna)

	!call matvec(hx,autovetores2(i,:),w90basis,vflagxa)

	!call prodintsq(vna,vflagxa,w90basis,vxa(i))


	!call vecconjg(autovetores2(i,:),w90basis,vma)

	!call matvec(hy,autovetores2(n,:),w90basis,vflagya)

	!call prodintsq(vma,vflagya,w90basis,vya(i))


	!berryf=-2.*(((vxa(i))*vya(i))/((energias2(n)-energias2(i))**2+(gammas)**2))

	berryf=0.

	else

	berryf=-2.*(((vx(i))*vy(i))/((energias(n)-energias(i))**2+(gammas)**2))

	end if

	berry=berry+berryf

	end do



end subroutine berryc



	

