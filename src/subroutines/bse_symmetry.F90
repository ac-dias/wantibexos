module bse_symmetry

	! Symmetry group of a tight-binding system: the point group and the space-group operations of the structure given by
	! the lattice, the orbital centres (sites) and what the orbitals are, checked against the Hamiltonian, and the part of
	! the group that maps the k mesh onto itself.  One analysis serves the standalone utility (utils/wtb_symmetry.F90,
	! bin/wtb_symmetry.x) and the output of the BSE solver (log_bse_optics.dat).
	!
	! What is searched.  The candidates are the rotations R (proper and improper, up to 48) that map the lattice onto
	! itself, found from the lattice vectors of the length and mutual angles of a_1, a_2, a_3.  For each R the
	! translations t = r_j - R r_0 that map a reference site r_0 onto a site r_j of its own kind are tried, and (R, t) is
	! an operation of the space group when every site goes onto a site of the same kind, R r_s + t = r_s' + L (L a
	! lattice vector).  A site is a group of orbitals with the same centre; two sites are of the same kind when they have
	! the same species, number of orbitals and shells (l of the orbitals, when the basis_set file gives them) and the same
	! on-site energy, summed over the orbitals of the site (that sum does not depend on how the orbitals of a shell
	! are rotated into each other).  The operations (R, t) with the same R are one operation of the point group;
	! R = 1 with t /= 0 are pure translations: the cell is a supercell of a smaller one.
	!
	! Hamiltonian check.  The geometry only says that the structure has the symmetry; the Hamiltonian can have less
	! (a strained or perturbed model with the old positions).  For every operation the norm of the block of H(R) (and of
	! S(R)) between the orbitals of two sites, <a,0|H|b,R_n>, must equal that of the image blocks <a',0|H|b',L_b-L_a+R R_n>.
	! This does not need the rotation matrices of the orbitals.
	!
	! k mesh.  R acts on the fractional coordinates of k by F_ij = (R b_j).a_i / 2 pi; the mesh f = (i/n + s) is mapped
	! onto itself when F (1/n_j) e_j and F s - s are multiples of 1/n.  The operations that keep the mesh are the
	! symmetry group that can be used with that mesh; F = -1 is the time-reversal pairing k -> -k of BSE_REAL.

	implicit none

	integer,parameter :: sym_maxop = 48
	real(kind=8),parameter :: sym_pi = 3.14159265358979323846d0

	type sym_group
		integer :: nop = 0                         ! operations of the point group
		integer :: ntrans = 1                      ! pure translations of the cell, with the identity
		integer :: nsite = 0, ntype = 0
		real(kind=8) :: rot(3,3,sym_maxop) = 0.0d0 ! Cartesian rotation matrices
		real(kind=8) :: tra(3,sym_maxop) = 0.0d0   ! fractional translations of the operations
		integer :: cls(sym_maxop) = 0              ! class: 1 E, 2 C2, 3 C3, 4 C4, 5 C6, 6 i, 7 sigma, 8 S3, 9 S4, 10 S6
		character(len=24) :: label(sym_maxop) = ''
		real(kind=8) :: axis(3,sym_maxop) = 0.0d0
		real(kind=8) :: angle(sym_maxop) = 0.0d0
		logical :: meshok(sym_maxop) = .false.
		real(kind=8) :: hdev(sym_maxop) = -1.0d0   ! largest relative deviation of the H blocks, -1 not checked
		real(kind=8) :: sdev(sym_maxop) = -1.0d0
		character(len=8) :: schoenflies = '?'
		character(len=10) :: hermann = '?'
	end type sym_group

	character(len=*),parameter :: sym_cls_name(10) = (/ 'E    ','C2   ','C3   ','C4   ','C6   ','i    ','sigma','S3   ','S4   ','S6   ' /)

contains

	subroutine sym_inv3(a,ainv)
		real(kind=8),intent(in) :: a(3,3)
		real(kind=8),intent(out) :: ainv(3,3)
		real(kind=8) :: det
		det = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) - a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) + a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
		ainv(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))/det
		ainv(1,2) = (a(1,3)*a(3,2)-a(1,2)*a(3,3))/det
		ainv(1,3) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))/det
		ainv(2,1) = (a(2,3)*a(3,1)-a(2,1)*a(3,3))/det
		ainv(2,2) = (a(1,1)*a(3,3)-a(1,3)*a(3,1))/det
		ainv(2,3) = (a(1,3)*a(2,1)-a(1,1)*a(2,3))/det
		ainv(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))/det
		ainv(3,2) = (a(1,2)*a(3,1)-a(1,1)*a(3,2))/det
		ainv(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))/det
	end subroutine sym_inv3

	function sym_det3(a) result(det)
		real(kind=8),intent(in) :: a(3,3)
		real(kind=8) :: det
		det = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) - a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) + a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
	end function sym_det3

	subroutine sym_lattice_rotations(lat,tol,nrot,rot)

		! the rotations (3 x 3, Cartesian, proper and improper) that map the lattice, rows lat(i,:) = a_i, onto itself

		real(kind=8),intent(in) :: lat(3,3),tol
		integer,intent(out) :: nrot
		real(kind=8),intent(out) :: rot(3,3,sym_maxop)

		integer,parameter :: nn = 6
		real(kind=8),allocatable :: v(:,:),vl(:)
		integer :: nv,n1,n2,n3,i1,i2,i3,m
		real(kind=8) :: maxlen,len3(3),vec(3),g12,g13,g23,atinv(3,3),at(3,3),ap(3,3),r(3,3),rr(3,3),dev
		real(kind=8) :: a1(3),a2(3),a3(3)
		logical :: ok

		a1 = lat(1,:)
		a2 = lat(2,:)
		a3 = lat(3,:)
		len3(1) = sqrt(dot_product(a1,a1))
		len3(2) = sqrt(dot_product(a2,a2))
		len3(3) = sqrt(dot_product(a3,a3))
		maxlen = maxval(len3)
		g12 = dot_product(a1,a2)
		g13 = dot_product(a1,a3)
		g23 = dot_product(a2,a3)

		allocate(v(3,(2*nn+1)**3),vl((2*nn+1)**3))
		nv = 0
		do n1=-nn,nn
			do n2=-nn,nn
				do n3=-nn,nn
					vec = n1*a1+n2*a2+n3*a3
					if (sqrt(dot_product(vec,vec)) <= maxlen+tol .and. sqrt(dot_product(vec,vec)) > tol) then
						nv = nv+1
						v(:,nv) = vec
						vl(nv) = sqrt(dot_product(vec,vec))
					end if
				end do
			end do
		end do

		! columns of at are a_i; R at = ap (columns a_i')
		do m=1,3
			at(:,m) = lat(m,:)
		end do
		call sym_inv3(at,atinv)

		nrot = 0
		do i1=1,nv
			if (abs(vl(i1)-len3(1)) > tol) cycle
			do i2=1,nv
				if (abs(vl(i2)-len3(2)) > tol) cycle
				if (abs(dot_product(v(:,i1),v(:,i2))-g12) > tol*maxlen) cycle
				do i3=1,nv
					if (abs(vl(i3)-len3(3)) > tol) cycle
					if (abs(dot_product(v(:,i1),v(:,i3))-g13) > tol*maxlen) cycle
					if (abs(dot_product(v(:,i2),v(:,i3))-g23) > tol*maxlen) cycle
					ap(:,1) = v(:,i1)
					ap(:,2) = v(:,i2)
					ap(:,3) = v(:,i3)
					r = matmul(ap,atinv)
					rr = matmul(r,transpose(r))
					dev = maxval(abs(rr-reshape((/1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0/),(/3,3/))))
					ok = dev < 1.0d-3 .and. abs(abs(sym_det3(r))-1.0d0) < 1.0d-3 .and. nrot < sym_maxop
					if (ok) then
						nrot = nrot+1
						rot(:,:,nrot) = r
					end if
				end do
			end do
		end do
		deallocate(v,vl)

	end subroutine sym_lattice_rotations

	subroutine sym_classify(r,cls,axis,angle,label)

		! class of a rotation matrix, its axis (rotation axis; mirror normal; axis of the rotoreflection) and angle in degrees

		real(kind=8),intent(in) :: r(3,3)
		integer,intent(out) :: cls
		real(kind=8),intent(out) :: axis(3),angle
		character(len=*),intent(out) :: label

		real(kind=8) :: rp(3,3),tr,det,th,s,m(3,3),nrm
		integer :: i,j,jbest
		real(kind=8) :: best

		det = sym_det3(r)
		rp = r
		if (det < 0.0d0) rp = -r
		tr = rp(1,1)+rp(2,2)+rp(3,3)
		th = acos(max(-1.0d0,min(1.0d0,(tr-1.0d0)/2.0d0)))
		axis = (/0.0d0,0.0d0,1.0d0/)
		angle = 0.0d0
		if (th > 1.0d-3) then
			if (th < sym_pi-1.0d-3) then
				s = 2.0d0*sin(th)
				axis = (/rp(3,2)-rp(2,3),rp(1,3)-rp(3,1),rp(2,1)-rp(1,2)/)/s
			else
				m = 0.5d0*(rp+reshape((/1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0/),(/3,3/)))
				best = -1.0d0
				jbest = 1
				do j=1,3
					nrm = sqrt(dot_product(m(:,j),m(:,j)))
					if (nrm > best) then
						best = nrm
						jbest = j
					end if
				end do
				axis = m(:,jbest)/best
			end if
			nrm = sqrt(dot_product(axis,axis))
			axis = axis/nrm
			angle = th
			! canonical sense of the axis: first non-zero component positive
			do i=1,3
				if (abs(axis(i)) > 1.0d-6) then
					if (axis(i) < 0.0d0) then
						axis = -axis
						angle = 2.0d0*sym_pi-angle
					end if
					exit
				end if
			end do
		end if
		angle = angle*180.0d0/sym_pi

		if (det > 0.0d0) then
			if (th < 1.0d-3) then
				cls = 1
			else if (abs(tr+1.0d0) < 0.1d0) then
				cls = 2
			else if (abs(tr) < 0.1d0) then
				cls = 3
			else if (abs(tr-1.0d0) < 0.1d0) then
				cls = 4
			else
				cls = 5
			end if
		else
			tr = r(1,1)+r(2,2)+r(3,3)
			if (abs(tr+3.0d0) < 0.1d0) then
				cls = 6
			else if (abs(tr-1.0d0) < 0.1d0) then
				cls = 7
			else if (abs(tr+2.0d0) < 0.1d0) then
				cls = 8
			else if (abs(tr+1.0d0) < 0.1d0) then
				cls = 9
			else
				cls = 10
			end if
		end if
		label = trim(sym_cls_name(cls))

	end subroutine sym_classify

	subroutine sym_name_group(g)

		! Schoenflies and Hermann-Mauguin name of the point group from the numbers of operations of each class

		type(sym_group),intent(inout) :: g
		integer :: n(10),i
		character(len=8),parameter :: names(32) = (/ 'C1  ','Ci  ','Cs  ','C2  ','C2v ','C2h ','D2  ','D2h ','C3  ','C3v ','S6  ', &
			'D3  ','D3d ','C3h ','D3h ','C4  ','C4v ','S4  ','C4h ','D2d ','D4  ','D4h ','C6  ','C6v ','C6h ','D6  ','D6h ', &
			'T   ','Th  ','Td  ','O   ','Oh  ' /)
		character(len=10),parameter :: hm(32) = (/ '1         ','-1        ','m         ','2         ','mm2       ','2/m       ', &
			'222       ','mmm       ','3         ','3m        ','-3        ','32        ','-3m       ','-6        ','-6m2      ', &
			'4         ','4mm       ','-4        ','4/m       ','-42m      ','422       ','4/mmm     ','6         ','6mm       ', &
			'6/m       ','622       ','6/mmm     ','23        ','m-3       ','-43m      ','432       ','m-3m      ' /)
		! counts of [C2 C3 C4 C6 i sigma S3 S4 S6] (E not counted) of each group
		integer,parameter :: tab(9,32) = reshape((/ &
			0,0,0,0,0,0,0,0,0,  0,0,0,0,1,0,0,0,0,  0,0,0,0,0,1,0,0,0,  1,0,0,0,0,0,0,0,0,  1,0,0,0,0,2,0,0,0,  1,0,0,0,1,1,0,0,0, &
			3,0,0,0,0,0,0,0,0,  3,0,0,0,1,3,0,0,0,  0,2,0,0,0,0,0,0,0,  0,2,0,0,0,3,0,0,0,  0,2,0,0,1,0,0,0,2,  3,2,0,0,0,0,0,0,0, &
			3,2,0,0,1,3,0,0,2,  0,2,0,0,0,1,2,0,0,  3,2,0,0,0,4,2,0,0,  1,0,2,0,0,0,0,0,0,  1,0,2,0,0,4,0,0,0,  1,0,0,0,0,0,0,2,0, &
			1,0,2,0,1,1,0,2,0,  3,0,0,0,0,2,0,2,0,  5,0,2,0,0,0,0,0,0,  5,0,2,0,1,5,0,2,0,  1,2,0,2,0,0,0,0,0,  1,2,0,2,0,6,0,0,0, &
			1,2,0,2,1,1,2,0,2,  7,2,0,2,0,0,0,0,0,  7,2,0,2,1,7,2,0,2,  3,8,0,0,0,0,0,0,0,  3,8,0,0,1,3,0,0,8,  3,8,0,0,0,6,0,6,0, &
			9,8,6,0,0,0,0,0,0,  9,8,6,0,1,9,0,6,8 /),(/9,32/))

		n = 0
		do i=1,g%nop
			n(g%cls(i)) = n(g%cls(i))+1
		end do
		g%schoenflies = '?'
		g%hermann = '?'
		do i=1,32
			if (all(tab(:,i) == (/n(2),n(3),n(4),n(5),n(6),n(7),n(8),n(9),n(10)/))) then
				g%schoenflies = names(i)
				g%hermann = hm(i)
				exit
			end if
		end do

	end subroutine sym_name_group

	subroutine sym_find(lat,norb,tau,have_pos,spec,lq,have_lab,onsite,have_h,tol,rclu,g,nsite,site_of,perm_all,off_all,nfound,site_pos,site_kind,site_norb,nlatrot)

		! Space-group operations of the structure.  lat(i,:) = a_i; tau(norb,3): orbital centres; spec, lq: species and
		! shell l of each orbital (have_lab); onsite: diagonal of H(R=0) (have_h); without positions only the lattice
		! group is returned (every operation of the lattice, nsite = 0).  perm_all(site,op) and off_all(3,site,op) are
		! the images of the sites for the operations that were found (the first translation of each rotation).

		real(kind=8),intent(in) :: lat(3,3),tol,rclu
		integer,intent(in) :: norb
		real(kind=8),intent(in) :: tau(norb,3)
		logical,intent(in) :: have_pos,have_lab,have_h
		character(len=*),intent(in) :: spec(norb)
		integer,intent(in) :: lq(norb)
		real(kind=8),intent(in) :: onsite(norb)
		type(sym_group),intent(out) :: g
		integer,intent(out) :: nsite
		integer,allocatable,intent(out) :: site_of(:),perm_all(:,:),off_all(:,:,:)
		integer,intent(out) :: nfound,nlatrot
		real(kind=8),allocatable,intent(out) :: site_pos(:,:)
		integer,allocatable,intent(out) :: site_kind(:),site_norb(:)

		real(kind=8) :: rots(3,3,sym_maxop),ainv(3,3),pos(3),fp(3),fq(3),d(3),tvec(3),tfrac(3),rp(3)
		real(kind=8),allocatable :: spos(:,:),sfrac(:,:),strace(:)
		integer,allocatable :: sn(:),snl(:,:),stype(:),first(:),perm(:),offs(:,:),tcount(:)
		character(len=8),allocatable :: sspec(:)
		integer :: nrot,i,j,s,s2,t,m,k,ntype,nref,iref,nt,ir,ie
		logical :: same,ok,found,used(norb)
		integer,allocatable :: permsave(:,:),offsave(:,:,:)
		integer :: ntr_this,nm(3),order(sym_maxop)
		real(kind=8) :: lt0(3,3)
		real(kind=8) :: dd,ax_tmp(3)
		character(len=24) :: lab_tmp

		call sym_lattice_rotations(lat,tol,nrot,rots)
		nlatrot = nrot

		nsite = 0
		g%nop = 0
		g%ntrans = 1
		if (.not. have_pos) then
			g%nop = nrot
			do ir=1,nrot
				g%rot(:,:,ir) = rots(:,:,ir)
				g%tra(:,ir) = 0.0d0
			end do
			allocate(site_of(norb),perm_all(1,max(1,nrot)),off_all(3,1,max(1,nrot)),site_pos(3,1),site_kind(1),site_norb(1))
			site_of = 0
			site_pos = 0.0d0
			site_kind = 0
			site_norb = 0
			perm_all = 0
			off_all = 0
			nfound = nrot
			call sym_finish(g,order)
			return
		end if

		lt0 = transpose(lat)
		call sym_inv3(lt0,ainv)     ! r = sum_i f_i a_i  ->  f = ainv r

		! sites: orbitals with the same centre
		allocate(spos(3,norb),sn(norb),snl(0:4,norb),sspec(norb),strace(norb),site_of(norb))
		sn = 0
		snl = 0
		strace = 0.0d0
		do i=1,norb
			found = .false.
			do s=1,nsite
				if (sqrt(sum((tau(i,:)-spos(:,s)/dble(sn(s)))**2)) < rclu) then
					found = .true.
					exit
				end if
			end do
			if (.not. found) then
				nsite = nsite+1
				s = nsite
				spos(:,s) = 0.0d0
				sspec(s) = ''
				if (have_lab) sspec(s) = spec(i)
			end if
			site_of(i) = s
			spos(:,s) = spos(:,s)+tau(i,:)
			sn(s) = sn(s)+1
			if (have_lab) then
				if (lq(i) >= 0 .and. lq(i) <= 4) snl(lq(i),s) = snl(lq(i),s)+1
			end if
			if (have_h) strace(s) = strace(s)+onsite(i)
		end do

		! spos was the sum of the centres of the orbitals of each site: the centroid is the position of the site
		do s=1,nsite
			spos(:,s) = spos(:,s)/dble(sn(s))
		end do

		! kinds of site; the on-site energy is a good invariant only with the species of the orbitals (a shell is
		! rotated into itself); single Wannier functions are told apart by the Hamiltonian check instead
		allocate(stype(nsite),first(nsite))
		ntype = 0
		do s=1,nsite
			stype(s) = 0
			do t=1,ntype
				s2 = first(t)
				same = (trim(sspec(s)) == trim(sspec(s2))) .and. sn(s) == sn(s2) .and. all(snl(:,s) == snl(:,s2))
				if (same .and. have_h .and. have_lab) same = abs(strace(s)-strace(s2)) <= 0.05d0*sn(s)
				if (same) then
					stype(s) = t
					exit
				end if
			end do
			if (stype(s) == 0) then
				ntype = ntype+1
				first(ntype) = s
				stype(s) = ntype
			end if
		end do
		g%nsite = nsite
		g%ntype = ntype

		allocate(sfrac(3,nsite))
		do s=1,nsite
			sfrac(:,s) = matmul(ainv,spos(:,s))
		end do

		! the reference site: of the kind with the fewest sites
		allocate(tcount(ntype))
		tcount = 0
		do s=1,nsite
			tcount(stype(s)) = tcount(stype(s))+1
		end do
		iref = 1
		do s=1,nsite
			if (tcount(stype(s)) < tcount(stype(iref))) iref = s
		end do

		allocate(site_pos(3,nsite),site_kind(nsite),site_norb(nsite))
		site_pos = spos(:,1:nsite)
		site_kind = stype
		site_norb = sn(1:nsite)
		allocate(perm(nsite),offs(3,nsite))
		allocate(permsave(nsite,sym_maxop),offsave(3,nsite,sym_maxop))
		g%nop = 0
		do ir=1,nrot
			ntr_this = 0
			do j=1,nsite
				if (stype(j) /= stype(iref)) cycle
				! translation that takes the reference site onto site j
				rp = matmul(rots(:,:,ir),spos(:,iref))
				tvec = spos(:,j)-rp
				ok = .true.
				perm = 0
				do s=1,nsite
					pos = matmul(rots(:,:,ir),spos(:,s))+tvec
					fp = matmul(ainv,pos)
					found = .false.
					do s2=1,nsite
						if (stype(s2) /= stype(s)) cycle
						fq = fp-sfrac(:,s2)
						nm = nint(fq)
						d = matmul(transpose(lat),fq-dble(nm))
						if (sqrt(dot_product(d,d)) < tol) then
							perm(s) = s2
							offs(:,s) = nm
							found = .true.
							exit
						end if
					end do
					if (.not. found) then
						ok = .false.
						exit
					end if
				end do
				if (ok) then
					! the images must be distinct sites
					used = .false.
					do s=1,nsite
						if (used(perm(s))) ok = .false.
						used(perm(s)) = .true.
					end do
				end if
				if (.not. ok) cycle
				ntr_this = ntr_this+1
				if (ntr_this == 1) then
					g%nop = g%nop+1
					g%rot(:,:,g%nop) = rots(:,:,ir)
					tfrac = matmul(ainv,tvec)
					tfrac = tfrac-floor(tfrac+1.0d-4)
					where (abs(tfrac) < 1.0d-4) tfrac = 0.0d0
					g%tra(:,g%nop) = tfrac
					permsave(:,g%nop) = perm
					offsave(:,:,g%nop) = offs
				end if
			end do
			if (ntr_this > 0) then
				call sym_classify(rots(:,:,ir),m,ax_tmp,dd,lab_tmp)
				if (m == 1) g%ntrans = ntr_this
			end if
		end do
		nfound = g%nop
		allocate(perm_all(nsite,max(1,g%nop)),off_all(3,nsite,max(1,g%nop)))
		do ir=1,g%nop
			perm_all(:,ir) = permsave(:,ir)
			off_all(:,:,ir) = offsave(:,:,ir)
		end do
		call sym_finish(g,order)
		! the images of the sites follow the order of the operations
		permsave(:,1:g%nop) = perm_all(:,order(1:g%nop))
		offsave(:,:,1:g%nop) = off_all(:,:,order(1:g%nop))
		perm_all(:,1:g%nop) = permsave(:,1:g%nop)
		off_all(:,:,1:g%nop) = offsave(:,:,1:g%nop)
		deallocate(spos,sn,snl,sspec,strace,stype,first,sfrac,tcount,perm,offs,permsave,offsave)

	end subroutine sym_find

	subroutine sym_finish(g,order)

		! classes, axes, labels, order of the operations (identity first), name of the group; order(i) is the
		! original index of the operation that is now the i-th

		type(sym_group),intent(inout) :: g
		integer,intent(out) :: order(sym_maxop)
		integer :: i,j,k
		type(sym_group) :: h

		order = 0

		do i=1,g%nop
			call sym_classify(g%rot(:,:,i),g%cls(i),g%axis(:,i),g%angle(i),g%label(i))
		end do
		! sort by class, then by angle and axis
		do i=1,g%nop
			order(i) = i
		end do
		do i=2,g%nop
			j = i
			do while (j > 1)
				if (sym_before(g,order(j),order(j-1))) then
					k = order(j)
					order(j) = order(j-1)
					order(j-1) = k
					j = j-1
				else
					exit
				end if
			end do
		end do
		h = g
		do i=1,g%nop
			g%rot(:,:,i) = h%rot(:,:,order(i))
			g%tra(:,i) = h%tra(:,order(i))
			g%cls(i) = h%cls(order(i))
			g%label(i) = h%label(order(i))
			g%axis(:,i) = h%axis(:,order(i))
			g%angle(i) = h%angle(order(i))
		end do
		call sym_name_group(g)

	end subroutine sym_finish

	function sym_before(g,i,j) result(b)
		type(sym_group),intent(in) :: g
		integer,intent(in) :: i,j
		logical :: b
		if (g%cls(i) /= g%cls(j)) then
			b = g%cls(i) < g%cls(j)
		else if (abs(g%angle(i)-g%angle(j)) > 1.0d-3) then
			b = g%angle(i) < g%angle(j)
		else
			b = g%axis(1,i)+10.0d0*g%axis(2,i)+100.0d0*g%axis(3,i) > g%axis(1,j)+10.0d0*g%axis(2,j)+100.0d0*g%axis(3,j)
		end if
	end function sym_before

	elemental function sym_clean(x) result(y)
		real(kind=8),intent(in) :: x
		real(kind=8) :: y
		y = x
		if (abs(x) < 5.0d-6) y = 0.0d0
	end function sym_clean

	function sym_frac_map(r,lat) result(f)

		! F = A R A^-1 (A: rows a_i): the rotation acting on the fractional coordinates of k, f' = F f

		real(kind=8),intent(in) :: r(3,3),lat(3,3)
		real(kind=8) :: f(3,3),linv(3,3)
		call sym_inv3(lat,linv)
		f = matmul(matmul(lat,r),linv)
	end function sym_frac_map

	function sym_mesh_ok(f,ngrid,mshift) result(ok)

		! the mesh f = (i/n + s) is mapped onto itself by f -> F f (mod 1)

		real(kind=8),intent(in) :: f(3,3)
		integer,intent(in) :: ngrid(3)
		real(kind=8),intent(in) :: mshift(3)
		logical :: ok
		integer :: i,m
		real(kind=8) :: x

		ok = .true.
		do m=1,3
			do i=1,3
				x = f(i,m)*dble(ngrid(i))/dble(ngrid(m))
				if (abs(x-anint(x)) > 1.0d-3) ok = .false.
			end do
		end do
		do i=1,3
			x = (dot_product(f(i,:),mshift)-mshift(i))*dble(ngrid(i))
			if (abs(x-anint(x)) > 1.0d-3) ok = .false.
		end do

	end function sym_mesh_ok

	subroutine sym_hcheck(lat,g,norb,nsite,site_of,perm,off,nvec,rvec,hop,ihop,ovp,ran)

		! Hamiltonian check of the operations of g: relative deviation of the site-pair norm of H(R) (g%hdev) and S(R)
		! (g%sdev) between each block and its image, as described at the top of the module

		real(kind=8),intent(in) :: lat(3,3)
		type(sym_group),intent(inout) :: g
		integer,intent(in) :: norb,nsite,nvec
		integer,intent(in) :: site_of(norb),perm(nsite,*),off(3,nsite,*)
		real,intent(in) :: rvec(nvec,3)
		real,intent(in) :: hop(nvec,norb,norb),ihop(nvec,norb,norb),ovp(nvec,norb,norb)
		logical,intent(out) :: ran

		real(kind=8),allocatable :: nh(:,:,:),ns(:,:,:)
		integer,allocatable :: ivol(:,:,:),nint3(:,:),osite(:),ofirst(:)
		real(kind=8) :: ainv(3,3),f(3),dh,ds,hmax,smax,v1,v2,x
		integer :: iv,i,j,a,b,ir,n1min,n1max,n2min,n2max,n3min,n3max,m(3),mp(3),ivp
		integer :: rn(3,3),sizev
		integer(kind=8) :: cost
		real(kind=8) :: lt(3,3)

		ran = .false.
		g%hdev = -1.0d0
		g%sdev = -1.0d0
		cost = int(nvec,kind=8)*int(norb,kind=8)*int(norb,kind=8)
		if (cost > 6000000000_8 .or. int(nsite,kind=8)**2*int(nvec,kind=8) > 80000000_8) return

		lt = transpose(lat)
		call sym_inv3(lt,ainv)
		allocate(nint3(3,nvec))
		do iv=1,nvec
			f = matmul(ainv,dble(rvec(iv,:)))
			nint3(:,iv) = nint(f)
			if (maxval(abs(f-dble(nint3(:,iv)))) > 1.0d-2) return        ! R is not a lattice vector: not checked
		end do
		n1min = minval(nint3(1,:)); n1max = maxval(nint3(1,:))
		n2min = minval(nint3(2,:)); n2max = maxval(nint3(2,:))
		n3min = minval(nint3(3,:)); n3max = maxval(nint3(3,:))
		sizev = (n1max-n1min+1)*(n2max-n2min+1)*(n3max-n3min+1)
		if (sizev > 20000000) return
		allocate(ivol(n1min:n1max,n2min:n2max,n3min:n3max))
		ivol = 0
		do iv=1,nvec
			ivol(nint3(1,iv),nint3(2,iv),nint3(3,iv)) = iv
		end do

		! orbitals of each site
		allocate(ofirst(nsite+1),osite(norb))
		ofirst = 0
		do i=1,norb
			ofirst(site_of(i)+1) = ofirst(site_of(i)+1)+1
		end do
		ofirst(1) = 1
		do a=1,nsite
			ofirst(a+1) = ofirst(a)+ofirst(a+1)
		end do
		allocate(nh(nsite,nsite,nvec),ns(nsite,nsite,nvec))
		nh = 0.0d0
		ns = 0.0d0
		do iv=1,nvec
			do j=1,norb
				b = site_of(j)
				do i=1,norb
					a = site_of(i)
					nh(a,b,iv) = nh(a,b,iv)+dble(hop(iv,i,j))**2+dble(ihop(iv,i,j))**2
					ns(a,b,iv) = ns(a,b,iv)+dble(ovp(iv,i,j))**2
				end do
			end do
		end do
		nh = sqrt(nh)
		ns = sqrt(ns)
		hmax = max(maxval(nh),1.0d-30)
		smax = max(maxval(ns),1.0d-30)

		do ir=1,g%nop
			! the rotation on the integer coordinates of the lattice vectors
			rn = nint(matmul(matmul(ainv,g%rot(:,:,ir)),lt))
			dh = 0.0d0
			ds = 0.0d0
			do b=1,nsite
				do a=1,nsite
					do iv=1,nvec
						m = nint3(:,iv)
						mp = off(:,b,ir)-off(:,a,ir)+matmul(rn,m)
						if (mp(1) < n1min .or. mp(1) > n1max .or. mp(2) < n2min .or. mp(2) > n2max .or. &
						    mp(3) < n3min .or. mp(3) > n3max) then
							ivp = 0
						else
							ivp = ivol(mp(1),mp(2),mp(3))
						end if
						v1 = nh(a,b,iv)
						v2 = 0.0d0
						x = 0.0d0
						if (ivp > 0) then
							v2 = nh(perm(a,ir),perm(b,ir),ivp)
							x = ns(perm(a,ir),perm(b,ir),ivp)
						end if
						dh = max(dh,abs(v1-v2))
						ds = max(ds,abs(ns(a,b,iv)-x))
					end do
				end do
			end do
			g%hdev(ir) = dh/hmax
			g%sdev(ir) = ds/smax
		end do
		ran = .true.
		deallocate(nh,ns,ivol,nint3,ofirst,osite)

	end subroutine sym_hcheck

	subroutine sym_analyze(unit,lat,norb,tau,have_pos,spec,lq,have_lab,nvec,rvec,hop,ihop,ovp,have_h,ngrid,mshift,tol,rclu,g)

		! The whole analysis, written to unit (and returned in g).  lat(i,:) = a_i in Angstrom; tau(norb,3) the orbital
		! centres; spec, lq the species and shell l of each orbital (have_lab); hop, ihop, ovp the H(R) and S(R) of the
		! tight-binding file with their Cartesian lattice vectors rvec (have_h); ngrid, mshift the k mesh (the shift is a
		! fraction of the reciprocal lattice vectors, as SHIFT_1..3= of the input).

		integer,intent(in) :: unit,norb,nvec,ngrid(3)
		real(kind=8),intent(in) :: lat(3,3),tol,rclu,mshift(3)
		real(kind=8),intent(in) :: tau(norb,3)
		logical,intent(in) :: have_pos,have_lab,have_h
		character(len=*),intent(in) :: spec(norb)
		integer,intent(in) :: lq(norb)
		real,intent(in) :: rvec(nvec,3)
		real,intent(in) :: hop(nvec,norb,norb),ihop(nvec,norb,norb),ovp(nvec,norb,norb)
		type(sym_group),intent(out) :: g

		real(kind=8),allocatable :: onsite(:)
		integer :: i,j,iv0,nsite,nfound,k,ncls(10),nmesh,nuse,s,t
		integer,allocatable :: site_of(:),perm_all(:,:),off_all(:,:,:),site_kind(:),site_norb(:)
		real(kind=8),allocatable :: site_pos(:,:)
		integer :: nlatrot
		real(kind=8) :: f(3,3),fm(3,3),rmin,vol
		logical :: hran,trmesh,anybad
		type(sym_group) :: gm,gu
		character(len=24) :: lab
		character(len=8) :: sp
		real(kind=8),parameter :: htol = 5.0d-3

		allocate(onsite(norb))
		onsite = 0.0d0
		iv0 = 0
		if (have_h) then
			do i=1,nvec
				if (sqrt(dble(rvec(i,1))**2+dble(rvec(i,2))**2+dble(rvec(i,3))**2) < 1.0d-5) iv0 = i
			end do
			if (iv0 > 0) then
				do j=1,norb
					onsite(j) = dble(hop(iv0,j,j))
				end do
			end if
		end if

		call sym_find(lat,norb,tau,have_pos,spec,lq,have_lab,onsite,have_h .and. iv0 > 0,tol,rclu,g,nsite,site_of,perm_all,off_all,nfound, &
			site_pos,site_kind,site_norb,nlatrot)

		write(unit,'(A)') ' ----------------------------------------------------------------------------'
		write(unit,'(A)') ' SYMMETRY OF THE SYSTEM'
		write(unit,'(A)') ' ----------------------------------------------------------------------------'
		write(unit,'(A)') ' Lattice vectors (Angstrom):'
		do i=1,3
			write(unit,'(3F14.6)') lat(i,:)
		end do
		vol = abs(sym_det3(lat))
		write(unit,'(A,F14.4,A)') ' Cell volume: ',vol,' Angstrom^3'
		if (have_pos) then
			write(unit,'(A,I0,A,I0,A,I0,A,ES8.1,A,ES8.1,A)') ' Orbitals: ',norb,' on ',g%nsite,' sites of ',g%ntype, &
				' kinds (orbitals within ',rclu,' Angstrom form a site, positions compared to ',tol,' Angstrom)'
			if (.not. have_lab) write(unit,'(A)') ' (no species and shells: sites are told apart by their number of orbitals and by the Hamiltonian check)'
			if (g%nsite > 0) then
				! the kinds of site, listed by their first site
				call sym_list_kinds(unit,norb,nsite,site_of,spec,lq,have_lab,onsite,have_h .and. iv0 > 0 .and. have_lab)
				if (g%nsite <= 40) then
					write(unit,'(A)') '   sites (Cartesian, Angstrom):  #  kind  orbitals           x           y           z'
					do s=1,g%nsite
						write(unit,'(A,I4,I6,I9,3F12.5)') '                              ',s,site_kind(s),site_norb(s),site_pos(:,s)
					end do
				end if
			end if
		else
			write(unit,'(A)') ' No orbital centres are available (DFT=W without BSE_CENTER_FILE): only the symmetry of'
			write(unit,'(A)') ' the lattice is reported, which is the largest group the structure can have'
		end if
		if (have_pos) then
			if (g%ntrans > 1) then
				write(unit,'(A,I0,A)') ' Pure translations of the cell: ',g%ntrans, &
					' (the cell is a supercell: it repeats a smaller cell this many times)'
			else
				write(unit,'(A)') ' Pure translations of the cell: 1 (the cell is primitive)'
			end if
		end if
		write(unit,'(A,A,A,A,A,I0)') ' Point group: ',trim(g%schoenflies),'  (Hermann-Mauguin ',trim(g%hermann),'),  order ',g%nop
		if (have_pos .and. g%nop < nlatrot) then
			write(unit,'(A,I0,A,I0,A)') ' The lattice has ',nlatrot,' operations, the structure keeps ',g%nop,' of them.'
			if (g%nop == 1 .and. nlatrot > 1) then
				write(unit,'(A)') ' If the model should have more symmetry, its orbital centres may be symmetric only to a larger tolerance'
				write(unit,'(A)') ' (set WTB_SYM_TOL, in Angstrom, e.g. 0.05), or the orbitals of an atom were not grouped into one site (WTB_SYM_CLUSTER, in Angstrom,'
				write(unit,'(A)') ' the distance within which Wannier centres form one site; try 0.5).'
			end if
		end if
		ncls = 0
		do i=1,g%nop
			ncls(g%cls(i)) = ncls(g%cls(i))+1
		end do
		write(unit,'(A)',advance='no') ' Operations by class:'
		do k=1,10
			if (ncls(k) > 0) write(unit,'(A,I0,A)',advance='no') ' ',ncls(k),trim(sym_cls_name(k))
		end do
		write(unit,*)

		! the mesh
		nmesh = 0
		do i=1,g%nop
			f = sym_frac_map(g%rot(:,:,i),lat)
			g%meshok(i) = sym_mesh_ok(f,ngrid,mshift)
			if (g%meshok(i)) nmesh = nmesh+1
		end do
		fm = 0.0d0
		do i=1,3
			fm(i,i) = -1.0d0
		end do
		trmesh = sym_mesh_ok(fm,ngrid,mshift)

		! the Hamiltonian
		hran = .false.
		if (have_h .and. have_pos .and. g%nsite > 0) then
			call sym_hcheck(lat,g,norb,nsite,site_of,perm_all,off_all,nvec,rvec,hop,ihop,ovp,hran)
		end if

		write(unit,'(A)') ' Operations of the point group (rotation about the axis, normal of the mirror plane, axis of the'
		write(unit,'(A)') ' improper rotation; fractional translation of the first space-group operation with this rotation):'
		write(unit,'(A)') '   #  class   axis / normal (Cartesian)     angle   translation (fractional)   mesh   H dev.   S dev.'
		do i=1,g%nop
			if (g%hdev(i) >= 0.0d0) then
				write(unit,'(I4,2X,A5,2X,3F9.5,F9.2,2X,3F8.4,5X,A3,2ES9.1)') i,trim(g%label(i)),sym_clean(g%axis(:,i)),g%angle(i),g%tra(:,i), &
					merge('yes','no ',g%meshok(i)),g%hdev(i),g%sdev(i)
			else
				write(unit,'(I4,2X,A5,2X,3F9.5,F9.2,2X,3F8.4,5X,A3,A)') i,trim(g%label(i)),sym_clean(g%axis(:,i)),g%angle(i),g%tra(:,i), &
					merge('yes','no ',g%meshok(i)),'      n/a      n/a'
			end if
		end do
		write(unit,'(A)') ' Rotation matrices (Cartesian):'
		do i=1,g%nop
			write(unit,'(I4,2X,A5,3(2X,3F9.5))') i,trim(g%label(i)),(sym_clean(g%rot(1,j,i)),j=1,3),(sym_clean(g%rot(2,j,i)),j=1,3),(sym_clean(g%rot(3,j,i)),j=1,3)
		end do

		write(unit,'(A,3I5,A,3F8.4,A)') ' k mesh ',ngrid,'  shift',mshift,' (fractions of the reciprocal vectors):'
		write(unit,'(A,I0,A,I0,A)') '   operations that map the mesh onto itself: ',nmesh,' of ',g%nop,' (an operation that does not'
		write(unit,'(A)') '   keep the mesh cannot be used to reduce a calculation on it)'
		if (nmesh > 0) then
			gm%nop = 0
			do i=1,g%nop
				if (g%meshok(i)) then
					gm%nop = gm%nop+1
					gm%rot(:,:,gm%nop) = g%rot(:,:,i)
					gm%cls(gm%nop) = g%cls(i)
				end if
			end do
			call sym_name_group(gm)
			write(unit,'(A,A,A,A,A,I0)') '   point group of the mesh and the structure: ',trim(gm%schoenflies),' (', &
				trim(gm%hermann),'), order ',gm%nop
		end if
		if (trmesh) then
			write(unit,'(A)') '   the mesh is symmetric under k -> -k (time reversal pairs its k points, BSE_REAL= T can be used)'
		else
			write(unit,'(A)') '   the mesh is NOT symmetric under k -> -k (BSE_REAL= T cannot pair its k points: use an unshifted mesh'
			write(unit,'(A)') '   or a shift of half a mesh step)'
		end if

		if (hran) then
			anybad = .false.
			do i=1,g%nop
				if (g%hdev(i) > htol .or. g%sdev(i) > htol) anybad = .true.
			end do
			write(unit,'(A,ES9.2,A,ES9.2,A)') ' Hamiltonian check: largest deviation of the site-pair norms of H(R) ', &
				maxval(g%hdev(1:g%nop)),' and of S(R) ',maxval(g%sdev(1:g%nop)),' (relative to the largest block)'
			if (anybad) then
				write(unit,'(A,ES8.1,A)') ' WARNING: the Hamiltonian does not have the full symmetry of the structure (deviation above ',htol,'):'
				do i=1,g%nop
					if (g%hdev(i) > htol .or. g%sdev(i) > htol) write(unit,'(A,I0,A,A)') '   operation ',i,' ',trim(g%label(i))
				end do
				write(unit,'(A)') ' (a strained or perturbed model with the positions of the unperturbed one, or too small a tolerance)'
			else
				write(unit,'(A)') ' The Hamiltonian has the symmetry of the structure (all deviations below 5e-3).'
			end if
			gu%nop = 0
			do i=1,g%nop
				if (g%meshok(i) .and. g%hdev(i) <= htol .and. g%sdev(i) <= htol) then
					gu%nop = gu%nop+1
					gu%rot(:,:,gu%nop) = g%rot(:,:,i)
					gu%cls(gu%nop) = g%cls(i)
				end if
			end do
			if (gu%nop > 0) then
				call sym_name_group(gu)
				if (gu%schoenflies == '?') then
					write(unit,'(A,I0,A)') ' The ',gu%nop,' operations that keep the mesh and the Hamiltonian do not form a group (tolerances too tight)'
				else
					write(unit,'(A,A,A,A,A,I0)') ' Group usable with this mesh and this Hamiltonian: ',trim(gu%schoenflies),' (', &
						trim(gu%hermann),'), order ',gu%nop
				end if
			end if
		else if (have_h .and. have_pos) then
			write(unit,'(A)') ' Hamiltonian check skipped (the model is too large for it, or R is not a lattice vector)'
		end if
		write(unit,'(A)') ' ----------------------------------------------------------------------------'
		call flush(unit)

	end subroutine sym_analyze

	subroutine sym_list_kinds(unit,norb,nsite,site_of,spec,lq,have_lab,onsite,have_h)

		! one line per kind of site: species, number of orbitals, sites of the kind

		integer,intent(in) :: unit,norb,nsite
		integer,intent(in) :: site_of(norb),lq(norb)
		character(len=*),intent(in) :: spec(norb)
		logical,intent(in) :: have_lab,have_h
		real(kind=8),intent(in) :: onsite(norb)

		integer :: s,i,t,nk
		integer :: sn(nsite),stype(nsite),snl(0:4,nsite),first(nsite),cnt(nsite)
		character(len=8) :: sspec(nsite)
		real(kind=8) :: tr(nsite)
		logical :: same

		sn = 0
		snl = 0
		tr = 0.0d0
		sspec = ''
		do i=1,norb
			s = site_of(i)
			if (sn(s) == 0 .and. have_lab) sspec(s) = spec(i)
			sn(s) = sn(s)+1
			if (have_lab) then
				if (lq(i) >= 0 .and. lq(i) <= 4) snl(lq(i),s) = snl(lq(i),s)+1
			end if
			if (have_h) tr(s) = tr(s)+onsite(i)
		end do
		nk = 0
		cnt = 0
		do s=1,nsite
			stype(s) = 0
			do t=1,nk
				same = trim(sspec(s)) == trim(sspec(first(t))) .and. sn(s) == sn(first(t)) .and. all(snl(:,s) == snl(:,first(t)))
				if (same .and. have_h) same = abs(tr(s)-tr(first(t))) <= 0.05d0*sn(s)
				if (same) then
					stype(s) = t
					exit
				end if
			end do
			if (stype(s) == 0) then
				nk = nk+1
				first(nk) = s
				stype(s) = nk
			end if
			cnt(stype(s)) = cnt(stype(s))+1
		end do
		do t=1,nk
			s = first(t)
			if (have_lab) then
				write(unit,'(A,I0,A,A,A,I0,A,5(1X,I0),A,F10.3,A,I0,A)') '   kind ',t,': species ',trim(sspec(s)),', ',sn(s), &
					' orbitals per site (number with l = 0..4:',snl(:,s),'), on-site sum ',tr(s),' eV, ',cnt(t),' sites'
			else
				write(unit,'(A,I0,A,I0,A,F10.3,A,I0,A)') '   kind ',t,': ',sn(s),' orbitals per site, on-site sum ',tr(s),' eV, ',cnt(t),' sites'
			end if
		end do

	end subroutine sym_list_kinds

end module bse_symmetry
