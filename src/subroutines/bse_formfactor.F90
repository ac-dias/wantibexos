module bse_formfactor

	! Orbital form factors of the tight-binding basis (BSE_FF= T, BSE_FF_FILE),
	!
	!   F_mn(R; Q) = <m,0| exp(iQ.r) |n,R>,
	!
	! the *_ff.bin files of utils/siesta2wtb.py, paoflow2wtb.py and
	! wannier2wtb.py with --formfactor (layout in utils/wtb_formfactor.py).
	! With them the pair densities of the Q=0 BSE kernel are those of the
	! orbitals of the basis instead of point charges at their centres,
	!
	!   <c k1| exp(iq.r) |c' k2> = sum_mn c_m(k1)^* c'_n(k2) B_mn(k2; q),
	!   B(k; q) = sum_R exp(ik.R) F(R; q),
	!
	! for DFT=W and DFT=S alike (F(R; 0) is the overlap S(R), and the
	! eigenvectors are those of H(k) = sum_R exp(ik.R) H(R)). The direct term
	! takes q = k1-k2 at each of its shortest images (bse_q_images), which the
	! file holds for every mesh that divides its own. The exchange term
	! (BSE_FF_EXCHANGE= T) takes the G != 0 of the file, up to BSE_FF_ECUT:
	!
	!   K^x(vck, v'c'k') = spinf sum_G v(G) rho_vck(G) rho_v'c'k'(G)^*,
	!   rho_vck(G) = <c k| exp(iG.r) |v k>,
	!
	! v(G) the bare Coulomb interaction e^2/(eps0 G^2) per N_k Omega, cut at
	! half the cell height for SYSDIM 2D as in V2DT, and spinf = 2 (the
	! singlets of an unpolarized model, NP) or 1. The G = 0 term is left out,
	! so the BSE gives the macroscopic dielectric function with local fields.

	implicit none
	private

	type,public :: ff_data
		logical :: direct=.false.                 ! the kernel uses the form factors
		logical :: exchange=.false.               ! and adds the exchange term
		integer :: id=0
		integer :: norb=0,nr=0,nq=0,ndir=0,nexc=0
		integer,dimension(3) :: mesh=0
		real(8) :: ecut=0.0d0                     ! eV, the G of the file
		real(8),dimension(3,3) :: lat=0.0d0       ! lat(:,i) = a_i (A)
		real(8),allocatable :: tau(:,:)           ! (3,norb) orbital centres (A)
		real,allocatable :: rcart(:,:)            ! (3,nr) lattice vectors R (A)
		real(8),allocatable :: qc(:,:)            ! (3,nq) Q in units of b_i
		complex,allocatable :: f(:,:,:,:)         ! (norb,norb,nr,nq) F(m,n,iR,iQ)
		integer,dimension(3) :: klo=0,khi=0
		integer,allocatable :: qkey(:,:,:)        ! nint(Q_i mesh_i) -> iQ, direct Q
		integer :: ng=0                           ! exchange: the G used,
		real :: gmax=0.0                          ! the largest hbar^2 G^2/2m (eV),
		complex,allocatable :: u(:,:)             ! (ng,dimbse) sqrt(spinf v(G)) rho(G)
	end type ff_data

	integer,save :: ff_count=0

	! B(k2; q) of the last element of each thread, for the next one with the
	! same k2 and q (the band pairs of one pair of k-points)
	complex,allocatable,save :: bcache(:,:)
	integer,dimension(3),save :: bkey=-1
	!$omp threadprivate(bcache,bkey)

	public :: ff_read
	public :: ff_destroy
	public :: ff_setup_check
	public :: ff_mesh_check
	public :: ff_overlap_check
	public :: ff_q_index
	public :: ff_direct_vertices
	public :: ff_exchange_setup
	public :: ff_exchange_element

contains


subroutine ff_read(filename,ff,ok,message)

	! Read a form-factor file: little-endian stream, single-precision F

	character(len=*),intent(in) :: filename
	type(ff_data),intent(inout) :: ff
	logical,intent(out) :: ok
	character(len=*),intent(out) :: message

	character(len=8) :: magic
	integer(4),dimension(8) :: hdr
	integer(4),allocatable :: rv(:,:)
	integer :: unit,ios,ir,iq
	integer,dimension(3) :: n

	ok=.false.
	message=''
	call ff_destroy(ff)
	open(newunit=unit,file=trim(filename),access='stream',form='unformatted', &
		status='old',action='read',iostat=ios)
	if (ios /= 0) then
		message='Unable to open the form-factor file '//trim(filename)
		return
	end if
	read(unit,iostat=ios) magic
	if (ios /= 0 .or. magic /= 'WTBFF001') then
		message=trim(filename)//' is not a form-factor file of utils/wtb_formfactor.py'
		close(unit)
		return
	end if
	read(unit,iostat=ios) hdr
	if (ios == 0) read(unit,iostat=ios) ff%ecut
	if (ios == 0) read(unit,iostat=ios) ff%lat
	if (ios /= 0 .or. any(hdr(1:5) < 0) .or. hdr(1) == 0 .or. hdr(2) == 0 .or. &
		hdr(4)+hdr(5) /= hdr(3) .or. any(hdr(6:8) < 1)) then
		message='Corrupt header in the form-factor file '//trim(filename)
		close(unit)
		return
	end if
	ff%norb=hdr(1)
	ff%nr=hdr(2)
	ff%nq=hdr(3)
	ff%ndir=hdr(4)
	ff%nexc=hdr(5)
	ff%mesh=hdr(6:8)
	allocate(ff%tau(3,ff%norb),rv(3,ff%nr),ff%qc(3,ff%nq))
	allocate(ff%f(ff%norb,ff%norb,ff%nr,ff%nq))
	read(unit,iostat=ios) ff%tau
	if (ios == 0) read(unit,iostat=ios) rv
	if (ios == 0) read(unit,iostat=ios) ff%qc
	if (ios == 0) read(unit,iostat=ios) ff%f
	close(unit)
	if (ios /= 0) then
		message='The form-factor file '//trim(filename)//' is shorter than its header says'
		call ff_destroy(ff)
		return
	end if

	allocate(ff%rcart(3,ff%nr))
	do ir=1,ff%nr
		ff%rcart(:,ir)=real(rv(1,ir)*ff%lat(:,1)+rv(2,ir)*ff%lat(:,2)+rv(3,ir)*ff%lat(:,3))
	end do

	! the direct Q by their integer coordinates on the mesh of the file
	ff%klo=huge(1)
	ff%khi=-huge(1)
	do iq=1,ff%ndir
		n=nint(ff%qc(:,iq)*ff%mesh)
		ff%klo=min(ff%klo,n)
		ff%khi=max(ff%khi,n)
	end do
	if (ff%ndir == 0) then
		ff%klo=0
		ff%khi=0
	end if
	allocate(ff%qkey(ff%klo(1):ff%khi(1),ff%klo(2):ff%khi(2),ff%klo(3):ff%khi(3)))
	ff%qkey=0
	do iq=1,ff%ndir
		n=nint(ff%qc(:,iq)*ff%mesh)
		if (ff%qkey(n(1),n(2),n(3)) /= 0) then
			message='Repeated direct Q in the form-factor file '//trim(filename)
			call ff_destroy(ff)
			return
		end if
		ff%qkey(n(1),n(2),n(3))=iq
	end do

	ff_count=ff_count+1
	ff%id=ff_count
	ok=.true.

end subroutine ff_read


subroutine ff_destroy(ff)

	type(ff_data),intent(inout) :: ff

	if (allocated(ff%tau)) deallocate(ff%tau)
	if (allocated(ff%rcart)) deallocate(ff%rcart)
	if (allocated(ff%qc)) deallocate(ff%qc)
	if (allocated(ff%f)) deallocate(ff%f)
	if (allocated(ff%qkey)) deallocate(ff%qkey)
	if (allocated(ff%u)) deallocate(ff%u)
	ff%direct=.false.
	ff%exchange=.false.
	ff%id=0
	ff%norb=0
	ff%nr=0
	ff%nq=0
	ff%ndir=0
	ff%nexc=0
	ff%ng=0

end subroutine ff_destroy


subroutine ff_setup_check(ff,w90basis,rlat,ngrid,ok,message)

	! The file must belong to this model: its basis size, its lattice, and a
	! mesh that the BSE mesh divides (the q = k1-k2 of the BSE mesh are then
	! among the direct Q of the file)

	type(ff_data),intent(in) :: ff
	integer,intent(in) :: w90basis
	real,dimension(3,3),intent(in) :: rlat
	integer,dimension(3),intent(in) :: ngrid
	logical,intent(out) :: ok
	character(len=*),intent(out) :: message

	integer :: i
	real :: dlat

	ok=.false.
	message=''
	if (ff%norb /= w90basis) then
		write(message,'(A,I0,A,I0)') 'Form factors of ',ff%norb, &
			' orbitals for a tight-binding basis of ',w90basis
		return
	end if
	dlat=0.0
	do i=1,3
		dlat=max(dlat,maxval(abs(real(ff%lat(:,i))-rlat(i,:))))
	end do
	if (dlat > 1.0e-3) then
		write(message,'(A,ES9.2,A)') 'The lattice of the form-factor file differs from that of '// &
			'the tight-binding file by ',dlat,' A'
		return
	end if
	if (any(mod(ff%mesh,ngrid) /= 0)) then
		write(message,'(A,3(1X,I0),A,3(1X,I0))') 'The form factors were written for the mesh', &
			ff%mesh,', which the BSE mesh does not divide:',ngrid
		return
	end if
	ok=.true.

end subroutine ff_setup_check


subroutine ff_mesh_check(ff,kpts,rlat,ok,message)

	! Every shortest image of every q = k1-k2 of the BSE mesh (k1 the first
	! k-point, k2 all of them: the classes q+G of a uniform mesh) in the file

	type(ff_data),intent(in) :: ff
	real,dimension(:,:),intent(in) :: kpts          ! (3,nk)
	real,dimension(3,3),intent(in) :: rlat
	logical,intent(out) :: ok
	character(len=*),intent(out) :: message

	real,parameter :: pi=acos(-1.)
	integer :: ik,img,nimg,i
	real,dimension(3,27) :: gshift
	real,dimension(3) :: q

	ok=.false.
	message=''
	do ik=1,size(kpts,2)
		call bse_q_images(kpts(:,1),kpts(:,ik),rlat,nimg,gshift)
		do img=1,nimg
			q=kpts(:,1)-kpts(:,ik)-gshift(:,img)
			if (ff_q_index(ff,q,rlat) == 0) then
				write(message,'(A,3F10.5,A)') 'q =',(dot_product(q,rlat(i,:))/(2.0*pi),i=1,3), &
					' (units of b_i) of the BSE mesh is not in the form-factor file'
				return
			end if
		end do
	end do
	ok=.true.

end subroutine ff_mesh_check


subroutine ff_overlap_check(ff,stt,vector,kpts,rlat,dev)

	! max |<n k|n' k> - delta_nn'| over the bands of the BSE transitions,
	! with the overlaps from F(R; 0): the file and the eigenvectors use the
	! same basis and the same Bloch sums (DFT=S: F(R; 0) = S(R))

	type(ff_data),intent(in) :: ff
	integer,dimension(:,:),intent(in) :: stt          ! (4,dimbse): state, v, c, k
	complex,dimension(:,:,:),intent(in) :: vector     ! (orbital, band, k)
	real,dimension(:,:),intent(in) :: kpts            ! (3,nk)
	real,dimension(3,3),intent(in) :: rlat
	real,intent(out) :: dev

	integer :: ik,it,iq0,ir,nlo,nhi,n1,n2
	real :: ang,d
	complex :: s
	complex,allocatable :: b(:,:)

	dev=0.0
	iq0=ff_q_index(ff,(/ 0.0,0.0,0.0 /),rlat)
	if (iq0 == 0) return
	allocate(b(ff%norb,ff%norb))
	do ik=1,size(kpts,2)
		nlo=huge(1)
		nhi=-huge(1)
		do it=1,size(stt,2)
			if (stt(4,it) /= ik) cycle
			nlo=min(nlo,stt(2,it))
			nhi=max(nhi,stt(3,it))
		end do
		if (nhi < nlo) cycle
		b=cmplx(0.0,0.0)
		do ir=1,ff%nr
			ang=dot_product(kpts(:,ik),ff%rcart(:,ir))
			b=b+cmplx(cos(ang),sin(ang))*ff%f(:,:,ir,iq0)
		end do
		do n2=nlo,nhi
			do n1=nlo,nhi
				s=dot_product(vector(:,n1,ik),matmul(b,vector(:,n2,ik)))
				d=abs(s)
				if (n1 == n2) d=abs(s-cmplx(1.0,0.0))
				dev=max(dev,d)
			end do
		end do
	end do
	deallocate(b)

end subroutine ff_overlap_check


integer function ff_q_index(ff,q,rlat)

	! the direct Q of the file equal to q (Cartesian, 1/A), 0 if none

	type(ff_data),intent(in) :: ff
	real,dimension(3),intent(in) :: q
	real,dimension(3,3),intent(in) :: rlat

	real,parameter :: pi=acos(-1.)
	real,dimension(3) :: x
	integer,dimension(3) :: n
	integer :: i

	ff_q_index=0
	do i=1,3
		x(i)=dot_product(q,rlat(i,:))/(2.0*pi)*real(ff%mesh(i))
	end do
	n=nint(x)
	if (any(abs(x-real(n)) > 1.0e-2)) return
	if (any(n < ff%klo) .or. any(n > ff%khi)) return
	ff_q_index=ff%qkey(n(1),n(2),n(3))

end function ff_q_index


subroutine ff_direct_vertices(ff,q,kpt2,ik2,rlat,c1,c2,v1,v2,vc,vv)

	! vc = <c1 k1| exp(iq.r) |c2 k2> and the hole vertex
	! vv = <v2 k2| exp(-iq.r) |v1 k1> = conjg(<v1 k1| exp(iq.r) |v2 k2>)
	! for the image q = k1-k2-G of the transfer, from B(k2; q)

	type(ff_data),intent(in) :: ff
	real,dimension(3),intent(in) :: q,kpt2
	integer,intent(in) :: ik2
	real,dimension(3,3),intent(in) :: rlat
	complex,dimension(:),intent(in) :: c1,c2,v1,v2
	complex,intent(out) :: vc,vv

	integer :: iq,ir
	real :: ang

	iq=ff_q_index(ff,q,rlat)
	if (iq == 0) stop 'BSE_FF: a q of the BSE mesh is missing from BSE_FF_FILE'
	if (bkey(1) /= ff%id .or. bkey(2) /= ik2 .or. bkey(3) /= iq) then
		if (allocated(bcache)) then
			if (size(bcache,1) /= ff%norb) deallocate(bcache)
		end if
		if (.not. allocated(bcache)) allocate(bcache(ff%norb,ff%norb))
		bcache=cmplx(0.0,0.0)
		do ir=1,ff%nr
			ang=dot_product(kpt2,ff%rcart(:,ir))
			bcache=bcache+cmplx(cos(ang),sin(ang))*ff%f(:,:,ir,iq)
		end do
		bkey=(/ ff%id,ik2,iq /)
	end if
	vc=dot_product(c1,matmul(bcache,c2))
	vv=conjg(dot_product(v1,matmul(bcache,v2)))

end subroutine ff_direct_vertices


subroutine ff_exchange_setup(ff,stt,vector,kpts,rlat,ngrid,sysdim,spinf,ecut,ok,message)

	! U(G,T) = sqrt(spinf v(G)) <c k| exp(iG.r) |v k> for each transition
	! T = (v,c,k) of stt and each G != 0 of the file with hbar^2 G^2/2m up to
	! ecut (eV; all of them for ecut <= 0), so that the exchange term is
	! K^x(T,T') = sum_G U(G,T) U(G,T')^*

	type(ff_data),intent(inout) :: ff
	integer,dimension(:,:),intent(in) :: stt          ! (4,dimbse): state, v, c, k
	complex,dimension(:,:,:),intent(in) :: vector     ! (orbital, band, k)
	real,dimension(:,:),intent(in) :: kpts            ! (3,nk)
	real,dimension(3,3),intent(in) :: rlat
	integer,dimension(3),intent(in) :: ngrid
	character(len=*),intent(in) :: sysdim
	real,intent(in) :: spinf,ecut
	logical,intent(out) :: ok
	character(len=*),intent(out) :: message

	real(8),parameter :: cic=-0.0904756d3       ! -e^2/(2 eps0), eV A, as in coulomb_pot
	real(8),parameter :: hb2m2=3.80998212d0      ! hbar^2/2m, eV A^2
	real(8),dimension(3) :: b1,b2,b3,g
	real :: vcell,ang
	real(8) :: rc,g2,gpar,gz,tz,vg
	integer :: nk,dimbse,ig,ng,iq,ik,it,ir
	integer,allocatable :: gsel(:),first(:),order(:),fill(:)
	real,allocatable :: wg(:)
	complex,allocatable :: b(:,:)

	ok=.false.
	message=''
	if (trim(sysdim) /= '2D' .and. trim(sysdim) /= '3D') then
		message='BSE_FF_EXCHANGE needs SYSDIM 2D or 3D'
		return
	end if
	! the G from the lattice of the file, in double precision: the T(G) of the
	! G on the axis are then exactly 0 or 2
	b1=ff_cross(ff%lat(:,2),ff%lat(:,3))
	b2=ff_cross(ff%lat(:,3),ff%lat(:,1))
	b3=ff_cross(ff%lat(:,1),ff%lat(:,2))
	g2=2.0d0*acos(-1.0d0)/dot_product(ff%lat(:,1),b1)
	b1=g2*b1
	b2=g2*b2
	b3=g2*b3
	call vcell3D(rlat,vcell)
	nk=ngrid(1)*ngrid(2)*ngrid(3)
	dimbse=size(stt,2)
	rc=0.5d0*ff%lat(3,3)

	! the G of the file within the cutoff, and their weights sqrt(spinf v(G))
	allocate(gsel(ff%nexc),wg(ff%nexc))
	ng=0
	ff%gmax=0.0
	do ig=1,ff%nexc
		iq=ff%ndir+ig
		g=ff%qc(1,iq)*b1+ff%qc(2,iq)*b2+ff%qc(3,iq)*b3
		g2=dot_product(g,g)
		if (ecut > 0.0 .and. hb2m2*g2 > ecut*(1.0d0+1.0d-5)) cycle
		if (trim(sysdim) == '2D') then
			! the slab-truncated interaction of V2DT, cut at |z| = c/2
			gpar=sqrt(g(1)*g(1)+g(2)*g(2))
			gz=abs(g(3))
			if (gpar < 1.0d-6) then
				tz=1.0d0-cos(gz*rc)-gz*rc*sin(gz*rc)
			else
				tz=1.0d0+exp(-gpar*rc)*((gz/gpar)*sin(gz*rc)-cos(gz*rc))
			end if
		else
			tz=1.0d0
		end if
		vg=-2.0d0*cic*tz/(dble(nk)*dble(vcell)*g2)
		if (vg <= 1.0d-12) cycle
		ff%gmax=max(ff%gmax,real(hb2m2*g2))
		ng=ng+1
		gsel(ng)=iq
		wg(ng)=real(sqrt(spinf*vg))
	end do
	ff%ng=ng
	if (allocated(ff%u)) deallocate(ff%u)
	allocate(ff%u(max(ng,1),dimbse))
	ff%u=cmplx(0.0,0.0)
	if (ng == 0) then
		ok=.true.
		return
	end if

	! the transitions of each k-point
	allocate(first(nk+1),order(dimbse),fill(nk))
	first=0
	do it=1,dimbse
		first(stt(4,it)+1)=first(stt(4,it)+1)+1
	end do
	first(1)=1
	do ik=1,nk
		first(ik+1)=first(ik+1)+first(ik)
	end do
	fill=0
	do it=1,dimbse
		ik=stt(4,it)
		order(first(ik)+fill(ik))=it
		fill(ik)=fill(ik)+1
	end do

	!$omp parallel do private(ik,ig,iq,ir,it,ang,b) schedule(dynamic)
	do ik=1,nk
		if (first(ik+1) == first(ik)) cycle
		allocate(b(ff%norb,ff%norb))
		do ig=1,ng
			iq=gsel(ig)
			b=cmplx(0.0,0.0)
			do ir=1,ff%nr
				ang=dot_product(kpts(:,ik),ff%rcart(:,ir))
				b=b+cmplx(cos(ang),sin(ang))*ff%f(:,:,ir,iq)
			end do
			do it=first(ik),first(ik+1)-1
				ff%u(ig,order(it))=wg(ig)*dot_product(vector(:,stt(3,order(it)),ik), &
					matmul(b,vector(:,stt(2,order(it)),ik)))
			end do
		end do
		deallocate(b)
	end do
	!$omp end parallel do

	deallocate(gsel,wg,first,order,fill)
	ok=.true.

end subroutine ff_exchange_setup


function ff_cross(a,b)

	real(8),dimension(3),intent(in) :: a,b
	real(8),dimension(3) :: ff_cross

	ff_cross=(/ a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1) /)

end function ff_cross


complex function ff_exchange_element(ff,i,j)

	! K^x(T_i,T_j) = sum_G U(G,T_i) U(G,T_j)^*

	type(ff_data),intent(in) :: ff
	integer,intent(in) :: i,j

	ff_exchange_element=dot_product(ff%u(:,j),ff%u(:,i))

end function ff_exchange_element


end module bse_formfactor
