program test_bse_mpi_blocks

	use mpi
	implicit none

	integer,parameter :: element_count=11
	integer,parameter :: block_elements=3
	integer :: ierr,rank,nodes,i,failures,total_failures,file_handle,complex_bytes
	integer(kind=8) :: file_element_count
	integer(kind=MPI_OFFSET_KIND) :: file_offset
	real :: real_values(element_count),reduced_values(element_count),expected_sum
	complex :: complex_values(element_count),expected_complex
	character(len=*),parameter :: test_filename='test_bse_mpi_blocks.dat'

	call MPI_INIT(ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,nodes,ierr)

	failures=0
	real_values=-1.0
	complex_values=cmplx(-1.0,-1.0)
	if (rank == 0) then
		do i=1,element_count
			real_values(i)=real(i)
			complex_values(i)=cmplx(real(i),-real(i))
		end do
	end if

	call bse_mpi_bcast_real_blocks(real_values,int(element_count,kind=8),0,MPI_COMM_WORLD,ierr,block_elements)
	if (ierr /= MPI_SUCCESS) failures=failures+1
	do i=1,element_count
		if (abs(real_values(i)-real(i)) > 1.0e-6) failures=failures+1
	end do

	call bse_mpi_bcast_complex_blocks(complex_values,int(element_count,kind=8),0,MPI_COMM_WORLD,ierr,block_elements)
	if (ierr /= MPI_SUCCESS) failures=failures+1
	do i=1,element_count
		if (abs(complex_values(i)-cmplx(real(i),-real(i))) > 1.0e-6) failures=failures+1
	end do

	reduced_values=real(rank+1)
	call bse_mpi_reduce_real_sum_blocks(reduced_values,int(element_count,kind=8),0,rank,MPI_COMM_WORLD,ierr,block_elements)
	if (ierr /= MPI_SUCCESS) failures=failures+1
	if (rank == 0) then
		expected_sum=real(nodes*(nodes+1))/2.0
		do i=1,element_count
			if (abs(reduced_values(i)-expected_sum) > 1.0e-6) failures=failures+1
		end do
	end if

	complex_values=cmplx(real(rank+1),-real(rank+1))
	call bse_mpi_allreduce_complex_sum_blocks(complex_values,int(element_count,kind=8),MPI_COMM_WORLD,ierr,block_elements)
	if (ierr /= MPI_SUCCESS) failures=failures+1
	expected_sum=real(nodes*(nodes+1))/2.0
	expected_complex=cmplx(expected_sum,-expected_sum)
	do i=1,element_count
		if (abs(complex_values(i)-expected_complex) > 1.0e-6) failures=failures+1
	end do

	! The BSE checkpoint path uses collective MPI-IO.  Give each rank a
	! contiguous region and force the same three-element transfer blocks.
	call MPI_FILE_OPEN(MPI_COMM_WORLD,test_filename,IOR(MPI_MODE_CREATE,MPI_MODE_RDWR),MPI_INFO_NULL,file_handle,ierr)
	if (ierr /= MPI_SUCCESS) failures=failures+1
	call MPI_FILE_SET_SIZE(file_handle,int(0,kind=MPI_OFFSET_KIND),ierr)
	if (ierr /= MPI_SUCCESS) failures=failures+1
	call MPI_TYPE_SIZE(MPI_COMPLEX,complex_bytes,ierr)
	file_offset=int(rank*element_count*complex_bytes,kind=MPI_OFFSET_KIND)
	call MPI_FILE_SET_VIEW(file_handle,file_offset,MPI_COMPLEX,MPI_COMPLEX,'native',MPI_INFO_NULL,ierr)
	if (ierr /= MPI_SUCCESS) failures=failures+1
	file_element_count=int(max(1,element_count-3*rank),kind=8)
	complex_values=cmplx(real(100*(rank+1)),-real(100*(rank+1)))
	call bse_mpi_file_write_all_complex_blocks(file_handle,complex_values,file_element_count,MPI_COMM_WORLD,block_elements,ierr)
	if (ierr /= MPI_SUCCESS) failures=failures+1
	call MPI_FILE_CLOSE(file_handle,ierr)
	if (ierr /= MPI_SUCCESS) failures=failures+1

	complex_values=cmplx(0.0,0.0)
	call MPI_FILE_OPEN(MPI_COMM_WORLD,test_filename,MPI_MODE_RDONLY,MPI_INFO_NULL,file_handle,ierr)
	if (ierr /= MPI_SUCCESS) failures=failures+1
	call MPI_FILE_SET_VIEW(file_handle,file_offset,MPI_COMPLEX,MPI_COMPLEX,'native',MPI_INFO_NULL,ierr)
	if (ierr /= MPI_SUCCESS) failures=failures+1
	call bse_mpi_file_read_all_complex_blocks(file_handle,complex_values,file_element_count,MPI_COMM_WORLD,block_elements,ierr)
	if (ierr /= MPI_SUCCESS) failures=failures+1
	call MPI_FILE_CLOSE(file_handle,ierr)
	if (ierr /= MPI_SUCCESS) failures=failures+1
	expected_complex=cmplx(real(100*(rank+1)),-real(100*(rank+1)))
	do i=1,int(file_element_count)
		if (abs(complex_values(i)-expected_complex) > 1.0e-6) failures=failures+1
	end do
	do i=int(file_element_count)+1,element_count
		if (abs(complex_values(i)) > 1.0e-6) failures=failures+1
	end do
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	if (rank == 0) then
		call MPI_FILE_DELETE(test_filename,MPI_INFO_NULL,ierr)
		if (ierr /= MPI_SUCCESS) failures=failures+1
	end if
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	call MPI_ALLREDUCE(failures,total_failures,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
	if (rank == 0) then
		if (total_failures == 0) then
			write(*,*) 'bse_mpi_blocks: PASS (11 elements in blocks of 3)'
		else
			write(*,*) 'bse_mpi_blocks: FAIL; mismatches: ',total_failures
		end if
	end if
	call MPI_FINALIZE(ierr)
	if (total_failures /= 0) stop 1

end program test_bse_mpi_blocks
