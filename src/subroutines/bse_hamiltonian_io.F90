subroutine bse_hamiltonian_write(filename,metadata,nrows,ncols,hbse,ok)

	implicit none

	character(len=*),intent(in) :: filename
	integer,dimension(7),intent(in) :: metadata
	integer,intent(in) :: nrows,ncols
	complex,intent(in) :: hbse(nrows,ncols)
	logical,intent(out) :: ok
	character(len=8),parameter :: magic = 'WTBSEH01'
	integer,parameter :: format_version = 1
	integer :: unit,io_status

	ok = .false.
	open(newunit=unit,file=trim(filename),status='replace',access='stream', &
	     form='unformatted',action='write',iostat=io_status)
	if (io_status /= 0) return

	write(unit,iostat=io_status) magic,format_version,metadata,nrows,ncols
	if (io_status == 0) write(unit,iostat=io_status) hbse
	close(unit)
	ok = (io_status == 0)

end subroutine bse_hamiltonian_write


subroutine bse_hamiltonian_read(filename,metadata,nrows,ncols,hbse,ok)

	implicit none

	character(len=*),intent(in) :: filename
	integer,dimension(7),intent(in) :: metadata
	integer,intent(in) :: nrows,ncols
	complex,intent(out) :: hbse(nrows,ncols)
	logical,intent(out) :: ok
	character(len=8),parameter :: magic = 'WTBSEH01'
	integer,parameter :: format_version = 1
	character(len=8) :: file_magic
	integer :: unit,io_status,file_version,file_metadata(7),file_nrows,file_ncols

	ok = .false.
	open(newunit=unit,file=trim(filename),status='old',access='stream', &
	     form='unformatted',action='read',iostat=io_status)
	if (io_status /= 0) return

	read(unit,iostat=io_status) file_magic,file_version,file_metadata,file_nrows,file_ncols
	if (io_status /= 0) then
		close(unit)
		return
	end if

	if (file_magic /= magic .or. file_version /= format_version .or. &
	    any(file_metadata /= metadata) .or. file_nrows /= nrows .or. file_ncols /= ncols) then
		close(unit)
		return
	end if

	read(unit,iostat=io_status) hbse
	close(unit)
	ok = (io_status == 0)

end subroutine bse_hamiltonian_read
