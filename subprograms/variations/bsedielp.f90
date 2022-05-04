subroutine dielbsep2(nthread,dimse,excitonvec,hopt1,hopt2,activity) !adaptado por João Paulo, mais lenta e resultado nthreads maior que o previsto

	use omp_lib
	implicit none


	integer :: dimse,nthread
	double precision,dimension(dimse) :: activity, exciton
	double precision :: actaux,actaux2
	double complex,dimension(dimse) :: hopt1,hopt2
	double complex,dimension(dimse,dimse) :: excitonvec

	integer :: i,j,k,no
	

	call OMP_SET_NUM_THREADS(nthread)

	activity=0.0
	!actaux=0.0


	! $OMP DO PRIVATE(actaux)
	! $OMP PARALLEL DO PRIVATE(actaux)
	do i=1,dimse

		actaux=0.0
		

		do j=1,dimse
			!$OMP PARALLEL PRIVATE(actaux2)
			actaux2=0.0
			do k=1,dimse
				actaux2=actaux2+(excitonvec(j,i)*hopt1(j)*conjg(excitonvec(k,i))*conjg(hopt2(k)))
			end do
			
			!$OMP CRITICAL
			actaux=actaux+actaux2
			!$OMP END CRITICAL
			!$OMP END PARALLEL
		end do

		activity(i)=actaux

	end do
	!$OMP END PARALLEL DO



end subroutine dielbsep2

