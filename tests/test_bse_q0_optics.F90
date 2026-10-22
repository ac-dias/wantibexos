program test_bse_q0_optics

	use bse_q_optics
	implicit none

	type(rmn_data) :: data
	complex :: vv(2),vc(2),hrx,hry,hrz,expected(3),correction(3)
	real :: rlat(3,3),kpoint(3),ev,ec,sme
	integer :: failures
	logical :: ok
	character(len=256) :: message

	failures=0
	rlat=0.0
	rlat(1,1)=1.0
	rlat(2,2)=1.0
	rlat(3,3)=1.0
	kpoint=(/0.37,0.0,0.0/)
	ev=-1.0
	ec=2.0
	sme=0.5

	data%num_wann=2
	data%nrpts=1
	allocate(data%rvec(3,1),data%rmn(3,2,2,1))
	data%rvec=0
	data%rmn=cmplx(0.0,0.0)
	data%rmn(1,1,2,1)=cmplx(1.0,-2.0)
	data%rmn(2,1,2,1)=cmplx(-0.5,0.25)

	vv=(/cmplx(0.0,0.0),cmplx(1.0,0.0)/)
	vc=(/cmplx(1.0,0.0),cmplx(0.0,0.0)/)
	hrx=cmplx(5.0,-1.0)
	hry=cmplx(-3.0,0.5)
	hrz=cmplx(2.0,4.0)

	correction=cmplx(0.0,ec-ev)/cmplx(ec-ev,sme)* &
		(/cmplx(1.0,-2.0),cmplx(-0.5,0.25),cmplx(0.0,0.0)/)
	expected=(/hrx,hry,hrz/)+correction

	call rmn_apply_q0_optical_correction(data,kpoint,rlat,ev,vv,ec,vc,sme,hrx,hry,hrz,ok,message)
	if (.not. ok) then
		write(*,*) trim(message)
		failures=failures+1
	else
		if (abs(hrx-expected(1)) > 1.0e-6) failures=failures+1
		if (abs(hry-expected(2)) > 1.0e-6) failures=failures+1
		if (abs(hrz-expected(3)) > 1.0e-6) failures=failures+1
	end if

	sme=0.0
	hrx=cmplx(0.0,0.0)
	hry=cmplx(0.0,0.0)
	hrz=cmplx(0.0,0.0)
	expected=cmplx(0.0,1.0)* &
		(/cmplx(1.0,-2.0),cmplx(-0.5,0.25),cmplx(0.0,0.0)/)
	call rmn_apply_q0_optical_correction(data,kpoint,rlat,ev,vv,ec,vc,sme,hrx,hry,hrz,ok,message)
	if (.not. ok) then
		write(*,*) trim(message)
		failures=failures+1
	else
		if (abs(hrx-expected(1)) > 1.0e-6) failures=failures+1
		if (abs(hry-expected(2)) > 1.0e-6) failures=failures+1
		if (abs(hrz-expected(3)) > 1.0e-6) failures=failures+1
	end if

	call rmn_destroy(data)
	if (failures == 0) then
		write(*,*) 'bse_q0_optics: PASS'
	else
		write(*,*) 'bse_q0_optics: FAIL; mismatches: ',failures
		stop 1
	end if

end program test_bse_q0_optics
