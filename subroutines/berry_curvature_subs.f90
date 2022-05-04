subroutine berryct2(nocp,kx,ky,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		  ihopmatrices,efermi,energias,autovetores,nthread,gammas,bxx,bxy,bxz,byy,byz,bzz)
		   

	implicit none

	integer :: nthread

	double precision :: kx,ky,kz,ezz
	integer :: nocp,w90basis
	integer :: nvec
	integer :: n,i,j

	double precision,parameter :: ktol = 0.001

	integer,dimension(w90basis) :: ffactor

	double precision,dimension(3,3) :: rlat
	integer,dimension(nvec,3) :: rvec

	integer :: nsite

	double complex, dimension(w90basis,w90basis) :: hx,hy,hz
	double complex, dimension(w90basis,w90basis) :: autovetores,autovetores2

	double precision,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices

	double precision, dimension(w90basis) :: energias,energias2

	!double precision,dimension(3) :: mag

	double precision :: efermi

	double complex, dimension(w90basis) :: vflagx,vflagy,vflagz
	double complex, dimension(w90basis) :: vn,vm,vl

	double complex,dimension(w90basis) :: vx,vy,vz

	double complex,dimension(w90basis) :: vxa,vya
	double complex, dimension(w90basis) :: vflagxa,vflagya
	double complex, dimension(w90basis) :: vna,vma

	double complex :: bxx,bxy,bxz,byy,byz,bzz
	double complex :: bfxx,bfxy,bfxz,bfyy,bfyz,bfzz

	double precision:: gammas

	double precision :: berryaux


	
	call hlm(kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		 ihopmatrices,hx,hy,hz)
	

	bxx=cmplx(0.0,0.0)
	bxy=cmplx(0.0,0.0)
	bxz=cmplx(0.0,0.0)
	byy=cmplx(0.0,0.0)
	byz=cmplx(0.0,0.0)
	bzz=cmplx(0.0,0.0)		

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

	bfxx=0.
	bfxy=0.
	bfxz=0.
	bfyy=0.
	bfyz=0.
	bfzz=0.			

	else if (abs(energias(i)-energias(n)) .lt. 0.0001 ) then

	bfxx=0.
	bfxy=0.
	bfxz=0.
	bfyy=0.
	bfyz=0.
	bfzz=0.	

	else
	
	berryaux= (energias(n)-energias(i))**2+(gammas)**2
	
	bfxx=-2.*(((vx(i))*conjg(vx(i)))/(berryaux))
	bfxy=-2.*(((vx(i))*conjg(vy(i)))/(berryaux))
	bfxz=-2.*(((vx(i))*conjg(vz(i)))/(berryaux))
	bfyy=-2.*(((vy(i))*conjg(vy(i)))/(berryaux))
	bfyz=-2.*(((vy(i))*conjg(vz(i)))/(berryaux))
	bfzz=-2.*(((vz(i))*conjg(vz(i)))/(berryaux))					

	end if

	bxx=bxx+bfxx
	bxy=bxy+bfxy
	bxz=bxz+bfxz
	byy=byy+bfyy
	byz=byz+bfyz
	bzz=bzz+bfzz					

	end do

	end do



end subroutine berryct2


subroutine berryct(nocp,kx,ky,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		  ihopmatrices,efermi,energias,autovetores,nthread,gammas,berry)
		   

	implicit none

	integer :: nthread

	double precision :: kx,ky,kz,ezz
	integer :: nocp,w90basis
	integer :: nvec
	integer :: n,i,j

	double precision,parameter :: ktol = 0.001

	integer,dimension(w90basis) :: ffactor

	double precision,dimension(3,3) :: rlat
	integer,dimension(nvec,3) :: rvec

	integer :: nsite

	double complex, dimension(w90basis,w90basis) :: hx,hy,hz
	double complex, dimension(w90basis,w90basis) :: autovetores,autovetores2

	double precision,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices

	double precision, dimension(w90basis) :: energias,energias2

	!double precision,dimension(3) :: mag

	double precision :: efermi

	double complex, dimension(w90basis) :: vflagx,vflagy
	double complex, dimension(w90basis) :: vn,vm

	double complex,dimension(w90basis) :: vx,vy

	double complex,dimension(w90basis) :: vxa,vya
	double complex, dimension(w90basis) :: vflagxa,vflagya
	double complex, dimension(w90basis) :: vna,vma

	double complex :: berry,berryf

	double precision:: gammas

	double precision :: berryaux


	
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

	double precision :: kx,ky,kz,ezz
	integer :: nocp,w90basis
	integer :: nvec
	integer :: n,i,j

	double precision,parameter :: ktol = 0.001


	double precision,dimension(3,3) :: rlat
	integer,dimension(nvec,3) :: rvec
	integer,dimension(w90basis) :: ffactor
	integer :: nsite

	double complex, dimension(w90basis,w90basis) :: hx,hy,hz
	double complex, dimension(w90basis,w90basis) :: autovetores,autovetores2

	double precision,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices

	double precision, dimension(w90basis) :: energias,energias2

	!double precision,dimension(3) :: mag

	double precision :: efermi

	double complex, dimension(w90basis) :: vflagx,vflagy
	double complex, dimension(w90basis) :: vn,vm

	double complex,dimension(w90basis) :: vx,vy

	double complex,dimension(w90basis) :: vxa,vya
	double complex, dimension(w90basis) :: vflagxa,vflagya
	double complex, dimension(w90basis) :: vna,vma

	double complex :: berry,berryf

	double precision:: gammas


	
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



	

