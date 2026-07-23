subroutine bse_hamiltonian_write(filename,metadata,nrows,ncols,hbse,ok)

	implicit none

	character(len=*),intent(in) :: filename
	integer,dimension(3),intent(in) :: metadata
	integer,intent(in) :: nrows,ncols
	complex,intent(in) :: hbse(nrows,ncols)
	logical,intent(out) :: ok
	character(len=8),parameter :: magic = 'WTBSEH02'
	integer,parameter :: format_version = 2
	integer :: unit,io_status

	ok = .false.
	open(newunit=unit,file=trim(filename),status='replace',access='stream', &
	     form='unformatted',action='write',iostat=io_status)
	if (io_status /= 0) return

	write(unit,iostat=io_status) magic,format_version,metadata
	if (io_status == 0) write(unit,iostat=io_status) hbse
	close(unit)
	ok = (io_status == 0)

end subroutine bse_hamiltonian_write


subroutine bse_hamiltonian_read(filename,metadata,nrows,ncols,hbse,ok)

	implicit none

	character(len=*),intent(in) :: filename
	integer,dimension(3),intent(in) :: metadata
	integer,intent(in) :: nrows,ncols
	complex,intent(out) :: hbse(nrows,ncols)
	logical,intent(out) :: ok
	character(len=8),parameter :: magic = 'WTBSEH02'
	integer,parameter :: format_version = 2
	character(len=8) :: file_magic
	integer :: unit,io_status,file_version,file_metadata(3)

	ok = .false.
	open(newunit=unit,file=trim(filename),status='old',access='stream', &
	     form='unformatted',action='read',iostat=io_status)
	if (io_status /= 0) return

	read(unit,iostat=io_status) file_magic,file_version,file_metadata
	if (io_status /= 0) then
		close(unit)
		return
	end if

	if (file_magic /= magic .or. file_version /= format_version .or. &
	    any(file_metadata /= metadata)) then
		close(unit)
		return
	end if

	read(unit,iostat=io_status) hbse
	close(unit)
	ok = (io_status == 0)

end subroutine bse_hamiltonian_read


#ifdef MPI
subroutine bse_hamiltonian_write_parallel(filename,metadata,dimbse,hbse,locr,locc,lld, &
									  mb,nb,nprow,npcol,myrow,mycol,comm,ok)

	use mpi
	implicit none

	character(len=*),intent(in) :: filename
	integer,dimension(3),intent(in) :: metadata
	integer,intent(in) :: dimbse,locr,locc,lld,mb,nb,nprow,npcol,myrow,mycol,comm
	complex,intent(in) :: hbse(lld,max(1,locc))
	logical,intent(out) :: ok
	character(len=8),parameter :: magic = 'WTBSEH02'
	integer,parameter :: format_version = 2
	integer :: file_handle,filetype,io_status,mpi_ierr,status(MPI_STATUS_SIZE)
	integer :: char_bytes,integer_bytes,complex_bytes,darray_rank
	integer :: gsizes(2),distribs(2),dargs(2),psizes(2)
	integer(kind=MPI_OFFSET_KIND) :: header_bytes,file_bytes

	ok = .false.
	call MPI_Type_size(MPI_CHARACTER,char_bytes,mpi_ierr)
	call MPI_Type_size(MPI_INTEGER,integer_bytes,mpi_ierr)
	call MPI_Type_size(MPI_COMPLEX,complex_bytes,mpi_ierr)
	if (mpi_ierr /= MPI_SUCCESS) return

	header_bytes = int(8*char_bytes + 4*integer_bytes,MPI_OFFSET_KIND)
	file_bytes = header_bytes + int(dimbse,MPI_OFFSET_KIND)*int(dimbse,MPI_OFFSET_KIND)* &
		     int(complex_bytes,MPI_OFFSET_KIND)

	call MPI_File_open(comm,trim(filename),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,file_handle,io_status)
	if (io_status /= MPI_SUCCESS) return
	call MPI_File_set_size(file_handle,file_bytes,io_status)
	if (io_status /= MPI_SUCCESS) then
		call MPI_File_close(file_handle,io_status)
		return
	end if

	if (myrow == 0 .and. mycol == 0) then
		call MPI_File_write_at(file_handle,0_MPI_OFFSET_KIND,magic,8,MPI_CHARACTER,status,io_status)
		if (io_status == MPI_SUCCESS) then
			call MPI_File_write_at(file_handle,int(8*char_bytes,MPI_OFFSET_KIND),format_version,1,MPI_INTEGER,status,io_status)
		end if
		if (io_status == MPI_SUCCESS) then
			call MPI_File_write_at(file_handle,int(8*char_bytes+integer_bytes,MPI_OFFSET_KIND), &
							   metadata,3,MPI_INTEGER,status,io_status)
		end if
	end if
	call MPI_Bcast(io_status,1,MPI_INTEGER,0,comm,mpi_ierr)
	if (io_status /= MPI_SUCCESS .or. mpi_ierr /= MPI_SUCCESS) then
		call MPI_File_close(file_handle,mpi_ierr)
		return
	end if

	gsizes = (/ dimbse,dimbse /)
	distribs = (/ MPI_DISTRIBUTE_CYCLIC,MPI_DISTRIBUTE_CYCLIC /)
	dargs = (/ mb,nb /)
	psizes = (/ nprow,npcol /)
	darray_rank = myrow*npcol + mycol
	filetype = MPI_DATATYPE_NULL
	call MPI_Type_create_darray(nprow*npcol,darray_rank,2,gsizes,distribs,dargs,psizes, &
						   MPI_ORDER_FORTRAN,MPI_COMPLEX,filetype,io_status)
	if (io_status == MPI_SUCCESS) call MPI_Type_commit(filetype,io_status)
	if (io_status == MPI_SUCCESS) then
		call MPI_File_set_view(file_handle,header_bytes,MPI_COMPLEX,filetype,'native',MPI_INFO_NULL,io_status)
	end if
	if (io_status == MPI_SUCCESS) then
		call MPI_File_write_all(file_handle,hbse,locr*locc,MPI_COMPLEX,status,io_status)
	end if
	if (filetype /= MPI_DATATYPE_NULL) then
		call MPI_Type_free(filetype,mpi_ierr)
		if (io_status == MPI_SUCCESS .and. mpi_ierr /= MPI_SUCCESS) io_status = mpi_ierr
	end if
	call MPI_File_close(file_handle,mpi_ierr)
	if (io_status == MPI_SUCCESS .and. mpi_ierr /= MPI_SUCCESS) io_status = mpi_ierr
	ok = (io_status == MPI_SUCCESS)

end subroutine bse_hamiltonian_write_parallel


subroutine bse_hamiltonian_read_parallel(filename,metadata,dimbse,hbse,locr,locc,lld, &
									 mb,nb,nprow,npcol,myrow,mycol,comm,ok)

	use mpi
	implicit none

	character(len=*),intent(in) :: filename
	integer,dimension(3),intent(in) :: metadata
	integer,intent(in) :: dimbse,locr,locc,lld,mb,nb,nprow,npcol,myrow,mycol,comm
	complex,intent(out) :: hbse(lld,max(1,locc))
	logical,intent(out) :: ok
	character(len=8),parameter :: magic = 'WTBSEH02'
	integer,parameter :: format_version = 2
	character(len=8) :: file_magic
	integer :: file_handle,filetype,io_status,mpi_ierr,status(MPI_STATUS_SIZE)
	integer :: char_bytes,integer_bytes,complex_bytes,darray_rank,file_version,file_metadata(3)
	integer :: gsizes(2),distribs(2),dargs(2),psizes(2)
	integer(kind=MPI_OFFSET_KIND) :: header_bytes,file_bytes,actual_file_bytes
	logical :: header_ok

	ok = .false.
	call MPI_Type_size(MPI_CHARACTER,char_bytes,mpi_ierr)
	call MPI_Type_size(MPI_INTEGER,integer_bytes,mpi_ierr)
	call MPI_Type_size(MPI_COMPLEX,complex_bytes,mpi_ierr)
	if (mpi_ierr /= MPI_SUCCESS) return

	header_bytes = int(8*char_bytes + 4*integer_bytes,MPI_OFFSET_KIND)
	file_bytes = header_bytes + int(dimbse,MPI_OFFSET_KIND)*int(dimbse,MPI_OFFSET_KIND)* &
		     int(complex_bytes,MPI_OFFSET_KIND)

	call MPI_File_open(comm,trim(filename),MPI_MODE_RDONLY,MPI_INFO_NULL,file_handle,io_status)
	if (io_status /= MPI_SUCCESS) return

	header_ok = .true.
	if (myrow == 0 .and. mycol == 0) then
		call MPI_File_read_at(file_handle,0_MPI_OFFSET_KIND,file_magic,8,MPI_CHARACTER,status,io_status)
		if (io_status == MPI_SUCCESS) then
			call MPI_File_read_at(file_handle,int(8*char_bytes,MPI_OFFSET_KIND),file_version,1,MPI_INTEGER,status,io_status)
		end if
		if (io_status == MPI_SUCCESS) then
			call MPI_File_read_at(file_handle,int(8*char_bytes+integer_bytes,MPI_OFFSET_KIND), &
							  file_metadata,3,MPI_INTEGER,status,io_status)
		end if
		if (io_status == MPI_SUCCESS) call MPI_File_get_size(file_handle,actual_file_bytes,io_status)
		if (io_status /= MPI_SUCCESS .or. file_magic /= magic .or. file_version /= format_version .or. &
		    any(file_metadata /= metadata) .or. actual_file_bytes < file_bytes) header_ok = .false.
	end if
	call MPI_Bcast(header_ok,1,MPI_LOGICAL,0,comm,mpi_ierr)
	if (.not. header_ok .or. mpi_ierr /= MPI_SUCCESS) then
		call MPI_File_close(file_handle,mpi_ierr)
		return
	end if

	gsizes = (/ dimbse,dimbse /)
	distribs = (/ MPI_DISTRIBUTE_CYCLIC,MPI_DISTRIBUTE_CYCLIC /)
	dargs = (/ mb,nb /)
	psizes = (/ nprow,npcol /)
	darray_rank = myrow*npcol + mycol
	filetype = MPI_DATATYPE_NULL
	call MPI_Type_create_darray(nprow*npcol,darray_rank,2,gsizes,distribs,dargs,psizes, &
						   MPI_ORDER_FORTRAN,MPI_COMPLEX,filetype,io_status)
	if (io_status == MPI_SUCCESS) call MPI_Type_commit(filetype,io_status)
	if (io_status == MPI_SUCCESS) then
		call MPI_File_set_view(file_handle,header_bytes,MPI_COMPLEX,filetype,'native',MPI_INFO_NULL,io_status)
	end if
	if (io_status == MPI_SUCCESS) then
		call MPI_File_read_all(file_handle,hbse,locr*locc,MPI_COMPLEX,status,io_status)
	end if
	if (filetype /= MPI_DATATYPE_NULL) then
		call MPI_Type_free(filetype,mpi_ierr)
		if (io_status == MPI_SUCCESS .and. mpi_ierr /= MPI_SUCCESS) io_status = mpi_ierr
	end if
	call MPI_File_close(file_handle,mpi_ierr)
	if (io_status == MPI_SUCCESS .and. mpi_ierr /= MPI_SUCCESS) io_status = mpi_ierr
	ok = (io_status == MPI_SUCCESS)

end subroutine bse_hamiltonian_read_parallel
#endif
