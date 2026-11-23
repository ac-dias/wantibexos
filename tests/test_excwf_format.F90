! Output of the exciton wavefunctions: excwf (BSE_WF= T, direct BSE) and excwfi (q-path BSE) write exc_wf_*.dat
! and exc_den_*.dat.  The integer edit descriptors of their formats need a width with gfortran; this test
! calls both routines on a two-k-point system with nc = 1, nv = 2, reads the files back and checks the header,
! the number and the 8 columns of every line (the three integer columns exactly, the others to the printed
! precision), the length of the data lines, and the k-space density.  Run it with make -j1 test-excwf.
program test_excwf_format

	implicit none

	integer,parameter :: ngkpt=2,nc=1,nv=2,nstates=ngkpt*nc*nv
	integer,parameter :: excnum=9871,qptnum=2
	character(len=70) :: outputfolder
	integer :: nocp(ngkpt),stt(nstates,4)
	real :: kpt(ngkpt,3),qpt(4),excenergy
	complex :: ewf(nstates)
	integer :: failures

	failures=0
	outputfolder='./'
	excenergy=1.2345
	! occupied states per k point; the valence bands of k are nocp-nv+1..nocp, the conduction band nocp+1
	nocp=(/4,5/)
	! stt(i,:) = (state, valence index 1..nv, conduction index nv+1, k point)
	stt(1,:)=(/1,1,3,1/)
	stt(2,:)=(/2,2,3,1/)
	stt(3,:)=(/3,1,3,2/)
	stt(4,:)=(/4,2,3,2/)
	kpt(1,:)=(/0.0,0.0,0.0/)
	kpt(2,:)=(/0.25,-0.125,0.5/)
	qpt=(/0.0,0.5,0.25,0.0/)
	ewf(1)=cmplx(0.5,0.5)
	ewf(2)=cmplx(-0.25,0.125)
	ewf(3)=cmplx(0.0,-0.5)
	ewf(4)=cmplx(0.75,0.0)

	call excwf(outputfolder,ngkpt,kpt,nc,nv,nocp,stt,excenergy,excnum,ewf)
	call check_files('exc_wf_9871.dat','exc_den_9871.dat',7)

	call excwfi(outputfolder,ngkpt,kpt,qpt,nc,nv,nocp,stt,excenergy,excnum,qptnum,ewf)
	call check_files('exc_wf_9871_2.dat','exc_den_9871_2.dat',8)

	if (failures == 0) then
		write(*,*) 'excwf_format: PASS'
	else
		write(*,*) 'excwf_format: FAIL; mismatches: ',failures
		stop 1
	end if

contains

	subroutine check_files(wffile,denfile,nheader)

		character(len=*) :: wffile,denfile
		integer :: nheader
		integer :: i,k,u,ios,iocc,icond,ival,ncols
		real :: x(3),re,im,den(ngkpt),dens
		character(len=200) :: line

		! the wavefunction file: header, then one line per transition (nocp, k, conduction band, valence band, amplitude)
		u=20
		open(unit=u,file=wffile,status='old',iostat=ios)
		if (ios /= 0) then
			write(*,*) 'cannot open ',wffile
			failures=failures+1
			return
		end if
		do i=1,nheader
			read(u,"(A)",iostat=ios) line
			if (ios /= 0 .or. line(2:2) /= '#') then
				write(*,*) wffile,': header line ',i,' is missing or is not a comment'
				failures=failures+1
			end if
		end do
		do i=1,nstates
			read(u,"(A)",iostat=ios) line
			if (ios /= 0) then
				write(*,*) wffile,': data line ',i,' is missing'
				failures=failures+1
				exit
			end if
			! (I6,3F15.4,2I6,2E15.4): 93 characters
			if (len_trim(line) /= 93) then
				write(*,*) wffile,': data line ',i,' has ',len_trim(line),' characters, expected 93'
				failures=failures+1
			end if
			ncols=count_columns(line)
			if (ncols /= 8) then
				write(*,*) wffile,': data line ',i,' has ',ncols,' columns, expected 8'
				failures=failures+1
			end if
			read(line,*,iostat=ios) iocc,x(1),x(2),x(3),icond,ival,re,im
			if (ios /= 0) then
				write(*,*) wffile,': data line ',i,' cannot be read back: ',trim(line)
				failures=failures+1
				cycle
			end if
			k=stt(i,4)
			if (iocc /= nocp(k)) failures=failures+1
			if (icond /= nocp(k)+stt(i,3)-nv) failures=failures+1
			if (ival /= nocp(k)-nv+stt(i,2)) failures=failures+1
			if (any(abs(x-kpt(k,:)) > 1.0e-4)) failures=failures+1
			if (abs(re-real(ewf(i))) > 1.0e-4*max(1.0,abs(real(ewf(i))))) failures=failures+1
			if (abs(im-aimag(ewf(i))) > 1.0e-4*max(1.0,abs(aimag(ewf(i))))) failures=failures+1
		end do
		read(u,"(A)",iostat=ios) line
		if (ios == 0) then
			write(*,*) wffile,': more lines than states'
			failures=failures+1
		end if
		close(u,status='delete')

		! the density file: header, then |ewf|^2 summed over the nc*nv transitions of every k point
		den=0.0
		do i=1,nstates
			den(stt(i,4))=den(stt(i,4))+real(ewf(i)*conjg(ewf(i)))
		end do
		open(unit=u,file=denfile,status='old',iostat=ios)
		if (ios /= 0) then
			write(*,*) 'cannot open ',denfile
			failures=failures+1
			return
		end if
		do i=1,nheader
			read(u,"(A)",iostat=ios) line
			if (ios /= 0 .or. line(2:2) /= '#') then
				write(*,*) denfile,': header line ',i,' is missing or is not a comment'
				failures=failures+1
			end if
		end do
		do k=1,ngkpt
			read(u,*,iostat=ios) x(1),x(2),x(3),dens
			if (ios /= 0) then
				write(*,*) denfile,': line ',k,' cannot be read'
				failures=failures+1
				cycle
			end if
			if (any(abs(x-kpt(k,:)) > 1.0e-4)) failures=failures+1
			if (abs(dens-den(k)) > 1.0e-3*max(1.0,den(k))) failures=failures+1
		end do
		close(u,status='delete')

	end subroutine check_files

	integer function count_columns(line)

		character(len=*) :: line
		integer :: i
		logical :: inside

		count_columns=0
		inside=.false.
		do i=1,len_trim(line)
			if (line(i:i) /= ' ') then
				if (.not. inside) count_columns=count_columns+1
				inside=.true.
			else
				inside=.false.
			end if
		end do

	end function count_columns

end program test_excwf_format
