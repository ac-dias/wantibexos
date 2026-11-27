module bse_eph

    implicit none
    private

    type, public :: eph_data
        integer :: nwan   = 0
        integer :: nk     = 0
        integer :: nq     = 0
        integer :: nmodes = 0
        real,    allocatable :: kpts(:,:)        ! (3, nk)  fractional reciprocal coords
        real,    allocatable :: qpts(:,:)        ! (3, nq)  fractional reciprocal coords
        real,    allocatable :: omega(:,:)       ! (nmodes, nq)  phonon frequencies in eV
        complex, allocatable :: g(:,:,:,:,:)     ! (nwan, nwan, nmodes, nk, nq)  Wannier-gauge g^W
    end type eph_data

    public :: eph_read
    public :: eph_destroy
    public :: eph_build_kplusq
    public :: exciton_phonon_linewidth
    public :: exciton_phonon_phassist

contains

!---------------------------------------------------------------------------
subroutine eph_destroy(eph)
    type(eph_data), intent(inout) :: eph
    if (allocated(eph%kpts))  deallocate(eph%kpts)
    if (allocated(eph%qpts))  deallocate(eph%qpts)
    if (allocated(eph%omega)) deallocate(eph%omega)
    if (allocated(eph%g))     deallocate(eph%g)
    eph%nwan = 0;  eph%nk = 0;  eph%nq = 0;  eph%nmodes = 0
end subroutine eph_destroy


!---------------------------------------------------------------------------
subroutine eph_read(filename, eph, ok, message)
    ! Reads Wannier-gauge el-ph matrix elements g^W_{mn,nu}(k,q).
    ! m,n are Wannier function indices 1..nwan.
    ! All MPI ranks read independently (no Node==0 guard), following the
    ! rmn_read convention in bse_q_optics.F90.
    !
    ! File format:
    !   line 1 : comment (ignored)
    !   line 2 : nk  nq  nmodes  nwan
    !   nk lines: kx ky kz     (fractional reciprocal coords)
    !   nq lines: qx qy qz
    !   nk*nq*nmodes*nwan^2 data lines:
    !     iq  nu  omega_eV  ik  m  n  Re(g)  Im(g)
    !   Indices 1-based; m,n are Wannier indices 1..nwan.

    character(len=*), intent(in)    :: filename
    type(eph_data),   intent(inout) :: eph
    logical,          intent(out)   :: ok
    character(len=*), intent(out)   :: message

    integer :: unit, ios, ik, iq, nu, m, n
    integer :: iq_f, nu_f, ik_f, m_f, n_f
    real    :: omega_f, re, im
    character(len=256) :: header

    ok = .false.;  message = ''
    call eph_destroy(eph)

    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        message = 'Cannot open EPH file: '//trim(filename);  return
    end if

    read(unit,'(A)',iostat=ios) header
    if (ios == 0) read(unit,*,iostat=ios) eph%nk, eph%nq, eph%nmodes, eph%nwan
    if (ios /= 0 .or. eph%nk <= 0 .or. eph%nq <= 0 .or. &
        eph%nmodes <= 0 .or. eph%nwan <= 0) then
        message = 'Invalid header in EPH file';  close(unit);  call eph_destroy(eph);  return
    end if

    allocate(eph%kpts(3,eph%nk), eph%qpts(3,eph%nq), &
             eph%omega(eph%nmodes,eph%nq), &
             eph%g(eph%nwan,eph%nwan,eph%nmodes,eph%nk,eph%nq), stat=ios)
    if (ios /= 0) then
        message = 'Allocation failed for EPH arrays';  close(unit);  call eph_destroy(eph);  return
    end if
    eph%g = cmplx(0.0,0.0);  eph%omega = 0.0

    do ik = 1, eph%nk
        read(unit,*,iostat=ios) eph%kpts(1,ik), eph%kpts(2,ik), eph%kpts(3,ik)
        if (ios /= 0) then
            message = 'Error reading k-point coordinates from EPH file'
            close(unit);  call eph_destroy(eph);  return
        end if
    end do
    do iq = 1, eph%nq
        read(unit,*,iostat=ios) eph%qpts(1,iq), eph%qpts(2,iq), eph%qpts(3,iq)
        if (ios /= 0) then
            message = 'Error reading q-point coordinates from EPH file'
            close(unit);  call eph_destroy(eph);  return
        end if
    end do

    do iq = 1, eph%nq
      do nu = 1, eph%nmodes
        do ik = 1, eph%nk
          do m = 1, eph%nwan
            do n = 1, eph%nwan
                read(unit,*,iostat=ios) iq_f, nu_f, omega_f, ik_f, m_f, n_f, re, im
                if (ios /= 0) then
                    message = 'Premature end or error in EPH data records'
                    close(unit);  call eph_destroy(eph);  return
                end if
                if (iq_f < 1 .or. iq_f > eph%nq  .or. &
                    nu_f < 1 .or. nu_f > eph%nmodes .or. &
                    ik_f < 1 .or. ik_f > eph%nk   .or. &
                    m_f  < 1 .or. m_f  > eph%nwan  .or. &
                    n_f  < 1 .or. n_f  > eph%nwan) then
                    message = 'Out-of-range index in EPH data record'
                    close(unit);  call eph_destroy(eph);  return
                end if
                eph%g(m_f,n_f,nu_f,ik_f,iq_f) = cmplx(re,im)
                eph%omega(nu_f,iq_f) = omega_f
            end do
          end do
        end do
      end do
    end do

    close(unit)
    ok = .true.

end subroutine eph_read


!---------------------------------------------------------------------------
subroutine eph_build_kplusq(kpts_cart, nk, qpts_frac, nq, rlat, ktol, kplusq_idx)
    ! Builds lookup: kplusq_idx(ik,iq) = index of k+q in kpts grid.
    ! kpts_cart(3,nk): BSE k-points in Cartesian reciprocal units (Ang^-1).
    ! qpts_frac(3,nq): EPH q-points as fractions of the reciprocal lattice vectors.
    ! rlat(3,3):       real-space lattice vectors (rlat(i,:) = i-th lattice vector in Ang).
    !                  Reciprocal lattice vectors blat are derived internally via recvec.
    ! ktol:            tolerance in fractional units for matching k+q to a grid point.
    ! Returns -1 if k+q is not found within ktol.

    integer, intent(in)  :: nk, nq
    real,    intent(in)  :: kpts_cart(3,nk), qpts_frac(3,nq), rlat(3,3), ktol
    integer, intent(out) :: kplusq_idx(nk,nq)

    integer :: ik, iq, jk
    real    :: q_cart(3), kpq(3), dk(3), dist, bvol
    real    :: blat(3,3)     ! reciprocal lattice vectors (from recvec)
    real    :: inv1(3), inv2(3), inv3(3)  ! rows of blat^{-1} scaled by bvol
    real    :: kpq_frac(3), kj_frac(3)

    ! Compute reciprocal lattice vectors from real-space ones.
    ! blat(i,:) = i-th reciprocal lattice vector.
    call recvec(rlat(1,:), rlat(2,:), rlat(3,:), blat(1,:), blat(2,:), blat(3,:))

    ! Precompute cross products for blat^{-1} (same trick as recvec):
    !   frac_coord = (k_cart . (b_j x b_k)) / (b_i . (b_j x b_k))
    inv1(1) = blat(2,2)*blat(3,3) - blat(2,3)*blat(3,2)
    inv1(2) = blat(2,3)*blat(3,1) - blat(2,1)*blat(3,3)
    inv1(3) = blat(2,1)*blat(3,2) - blat(2,2)*blat(3,1)
    inv2(1) = blat(3,2)*blat(1,3) - blat(3,3)*blat(1,2)
    inv2(2) = blat(3,3)*blat(1,1) - blat(3,1)*blat(1,3)
    inv2(3) = blat(3,1)*blat(1,2) - blat(3,2)*blat(1,1)
    inv3(1) = blat(1,2)*blat(2,3) - blat(1,3)*blat(2,2)
    inv3(2) = blat(1,3)*blat(2,1) - blat(1,1)*blat(2,3)
    inv3(3) = blat(1,1)*blat(2,2) - blat(1,2)*blat(2,1)
    bvol = blat(1,1)*inv1(1) + blat(1,2)*inv1(2) + blat(1,3)*inv1(3)

    do iq = 1, nq
        ! Convert q from fractional (of reciprocal lattice) to Cartesian
        q_cart = qpts_frac(1,iq)*blat(1,:) + qpts_frac(2,iq)*blat(2,:) + qpts_frac(3,iq)*blat(3,:)
        do ik = 1, nk
            kpq = kpts_cart(:,ik) + q_cart
            ! Convert kpq to fractional, fold into [0,1)
            kpq_frac(1) = (kpq(1)*inv1(1)+kpq(2)*inv1(2)+kpq(3)*inv1(3)) / bvol
            kpq_frac(2) = (kpq(1)*inv2(1)+kpq(2)*inv2(2)+kpq(3)*inv2(3)) / bvol
            kpq_frac(3) = (kpq(1)*inv3(1)+kpq(2)*inv3(2)+kpq(3)*inv3(3)) / bvol
            kpq_frac = kpq_frac - floor(kpq_frac)
            kplusq_idx(ik,iq) = -1
            do jk = 1, nk
                kj_frac(1) = (kpts_cart(1,jk)*inv1(1)+kpts_cart(2,jk)*inv1(2)+kpts_cart(3,jk)*inv1(3)) / bvol
                kj_frac(2) = (kpts_cart(1,jk)*inv2(1)+kpts_cart(2,jk)*inv2(2)+kpts_cart(3,jk)*inv2(3)) / bvol
                kj_frac(3) = (kpts_cart(1,jk)*inv3(1)+kpts_cart(2,jk)*inv3(2)+kpts_cart(3,jk)*inv3(3)) / bvol
                dk = kpq_frac - kj_frac
                dk = dk - nint(dk)
                dist = sqrt(dk(1)**2 + dk(2)**2 + dk(3)**2)
                if (dist < ktol) then
                    kplusq_idx(ik,iq) = jk
                    exit
                end if
            end do
        end do
    end do

end subroutine eph_build_kplusq


!---------------------------------------------------------------------------
! Private helpers
!---------------------------------------------------------------------------

real function bose_einstein(omega, temp)
    real, intent(in) :: omega, temp
    real, parameter  :: kb = 8.617333e-5   ! eV/K
    if (temp <= 0.0 .or. omega <= 0.0) then
        bose_einstein = 0.0
    else
        bose_einstein = 1.0 / (exp(omega / (kb*temp)) - 1.0)
    end if
end function bose_einstein


real function lorentzian(x, eta)
    real, intent(in) :: x, eta
    real, parameter  :: pi = 3.14159265358979
    lorentzian = eta / (pi * (x**2 + eta**2))
end function lorentzian


! Build reverse pair-state lookup: pair_idx(iv,ic_local,ik) -> pair state j.
! stt_bse(2,j) = iv  (valence band, 1..nv in vector)
! stt_bse(3,j) = nv+ic  (conduction band in vector, so ic_local = stt_bse(3,j)-nv, 1..nc)
! stt_bse(4,j) = ik
subroutine build_pair_idx(stt_bse, dimbse, nv, nc, ngkpt, pair_idx)
    integer, intent(in)  :: dimbse, nv, nc, ngkpt
    integer, intent(in)  :: stt_bse(4,dimbse)
    integer, intent(out) :: pair_idx(nv,nc,ngkpt)

    integer :: j, iv, ic, ik

    pair_idx = 0
    do j = 1, dimbse
        iv = stt_bse(2,j)
        ic = stt_bse(3,j) - nv
        ik = stt_bse(4,j)
        if (iv >= 1 .and. iv <= nv .and. &
            ic >= 1 .and. ic <= nc .and. &
            ik >= 1 .and. ik <= ngkpt) then
            pair_idx(iv,ic,ik) = j
        end if
    end do

end subroutine build_pair_idx


! Rotate nwan x nwan Wannier-gauge g^W to an nb x nb Bloch-gauge block.
! g_out(b,b') = sum_{mn} conj(U(m,b,ik_l)) * gW(m,n) * U(n,b',ik_r)
! Bands b_start .. b_start+nb-1 are extracted from vector(:,:,ik_l/r).
! Phase cancellation in g^{ex} is guaranteed because the same vector(:)
! arrays appear here and in the BSE Hamiltonian build.
subroutine rotate_gW_block(nwan, nb, b_start, gW, vector, ik_l, ik_r, nc_nv, ngkpt, g_out)
    integer, intent(in)  :: nwan, nb, b_start, ik_l, ik_r, nc_nv, ngkpt
    complex, intent(in)  :: gW(nwan,nwan), vector(nwan,nc_nv,ngkpt)
    complex, intent(out) :: g_out(nb,nb)

    complex :: vl(nwan,nb), vr(nwan,nb), tmp(nwan,nb)
    integer :: ib

    do ib = 1, nb
        vl(:,ib) = vector(:, b_start+ib-1, ik_l)
        vr(:,ib) = vector(:, b_start+ib-1, ik_r)
    end do

    ! tmp = gW * vr   [nwan x nb]
    call cgemm('N','N', nwan, nb, nwan, cmplx(1.0,0.0), &
               gW, nwan, vr, nwan, cmplx(0.0,0.0), tmp, nwan)
    ! g_out = vl^H * tmp   [nb x nb]
    call cgemm('C','N', nb, nb, nwan, cmplx(1.0,0.0), &
               vl, nwan, tmp, nwan, cmplx(0.0,0.0), g_out, nb)

end subroutine rotate_gW_block


! Accumulate gex(dimbse,dimbse) for one (iq,nu) into the provided array.
! gex(lam,lp) = g^{ex}_{lp->lam}(q,nu), flat-dispersion approximation:
!   A^{lam}_{cv,k+q}(Q=q) approx A^{lam}_{cv,k+q}(Q=0) from Q=0 BSE.
!
! Electron (+): for each iv,ik: dressed_c = g_ee * hbse[pairs(iv,ic',ik),:];
!               gex += hbse[pairs(iv,ic,ikq),:]^H * dressed_c
! Hole    (-): for each ic,ik: dressed_v = g_hh * hbse[pairs(iv',ic,ikq),:];
!               gex -= hbse[pairs(iv,ic,ik),:]^H  * dressed_v
subroutine compute_gex(eph, iq, nu, vector, w90basis, nc, nv, ngkpt, nc_nv, &
                       hbse, dimbse, pair_idx, kplusq_idx, &
                       strip_c, dressed_c, ket_c, strip_v, dressed_v, ket_v, gex)
    type(eph_data), intent(in) :: eph
    integer, intent(in)        :: iq, nu, w90basis, nc, nv, ngkpt, nc_nv, dimbse
    complex, intent(in)        :: vector(w90basis,nc_nv,ngkpt)
    complex, intent(in)        :: hbse(dimbse,dimbse)
    integer, intent(in)        :: pair_idx(nv,nc,ngkpt)
    integer, intent(in)        :: kplusq_idx(ngkpt,eph%nq)
    complex, intent(inout)     :: strip_c(nc,dimbse), dressed_c(nc,dimbse), ket_c(nc,dimbse)
    complex, intent(inout)     :: strip_v(nv,dimbse), dressed_v(nv,dimbse), ket_v(nv,dimbse)
    complex, intent(out)       :: gex(dimbse,dimbse)

    integer :: ik, ikq, ic, iv, icp, ivp, j, jp
    complex :: gW(eph%nwan,eph%nwan)
    complex :: g_ee(nc,nc), g_hh(nv,nv)

    gex = cmplx(0.0,0.0)

    do ik = 1, ngkpt
        ikq = kplusq_idx(ik,iq)
        if (ikq < 1) cycle

        gW = eph%g(:,:,nu,ik,iq)

        ! g_ee(ic,ic') = U^H_c(k+q) * gW * U_c(k)   [b_start = nv+1 for conduction]
        call rotate_gW_block(eph%nwan, nc, nv+1, gW, vector, ikq, ik, nc_nv, ngkpt, g_ee)

        ! g_hh(iv,iv') = U^H_v(k) * gW * U_v(k+q)   [b_start = 1 for valence]
        call rotate_gW_block(eph%nwan, nv, 1,    gW, vector, ik, ikq, nc_nv, ngkpt, g_hh)

        ! --- Electron contribution (add) ---
        do iv = 1, nv
            do icp = 1, nc
                j = pair_idx(iv,icp,ik)
                if (j > 0) then
                    strip_c(icp,:) = hbse(j,:)
                else
                    strip_c(icp,:) = cmplx(0.0,0.0)
                end if
            end do
            ! dressed_c = g_ee * strip_c   [nc x dimbse]
            call cgemm('N','N', nc, dimbse, nc, cmplx(1.0,0.0), &
                       g_ee, nc, strip_c, nc, cmplx(0.0,0.0), dressed_c, nc)
            do ic = 1, nc
                j = pair_idx(iv,ic,ikq)
                if (j > 0) then
                    ket_c(ic,:) = hbse(j,:)
                else
                    ket_c(ic,:) = cmplx(0.0,0.0)
                end if
            end do
            ! gex += ket_c^H * dressed_c   [dimbse x dimbse]
            call cgemm('C','N', dimbse, dimbse, nc, cmplx(1.0,0.0), &
                       ket_c, nc, dressed_c, nc, cmplx(1.0,0.0), gex, dimbse)
        end do

        ! --- Hole contribution (subtract) ---
        do ic = 1, nc
            do ivp = 1, nv
                jp = pair_idx(ivp,ic,ikq)
                if (jp > 0) then
                    strip_v(ivp,:) = hbse(jp,:)
                else
                    strip_v(ivp,:) = cmplx(0.0,0.0)
                end if
            end do
            ! dressed_v = g_hh * strip_v   [nv x dimbse]
            call cgemm('N','N', nv, dimbse, nv, cmplx(1.0,0.0), &
                       g_hh, nv, strip_v, nv, cmplx(0.0,0.0), dressed_v, nv)
            do iv = 1, nv
                j = pair_idx(iv,ic,ik)
                if (j > 0) then
                    ket_v(iv,:) = hbse(j,:)
                else
                    ket_v(iv,:) = cmplx(0.0,0.0)
                end if
            end do
            ! gex -= ket_v^H * dressed_v   [dimbse x dimbse]
            call cgemm('C','N', dimbse, dimbse, nv, cmplx(-1.0,0.0), &
                       ket_v, nv, dressed_v, nv, cmplx(1.0,0.0), gex, dimbse)
        end do

    end do  ! ik

end subroutine compute_gex


!---------------------------------------------------------------------------
subroutine exciton_phonon_linewidth(eph, vector, w90basis, nc, nv, ngkpt, &
           hbse, stt_bse, dimbse, W, kplusq_idx, temp, eta, Gamma, Node, Nodes)
    ! Gauge-invariant exciton linewidth (Fermi's golden rule):
    !   Gamma(lam) = pi * sum_{q,nu,lp} |gex(lam,lp)|^2
    !                * [(n+1)*L(E_lam-E_lp-omega,eta) + n*L(E_lam-E_lp+omega,eta)]
    ! Uses the flat-dispersion (Q=0) approximation for BSE amplitudes.
    ! MPI: round-robin over q-points.

#ifdef MPI
    use mpi
#endif

    type(eph_data), intent(in) :: eph
    integer, intent(in)        :: w90basis, nc, nv, ngkpt, dimbse
    complex, intent(in)        :: vector(w90basis,nc+nv,ngkpt)
    complex, intent(in)        :: hbse(dimbse,dimbse)
    integer, intent(in)        :: stt_bse(4,dimbse)
    real,    intent(in)        :: W(dimbse)
    integer, intent(in)        :: kplusq_idx(ngkpt,eph%nq)
    real,    intent(in)        :: temp, eta
    real,    intent(out)       :: Gamma(dimbse)
    integer, intent(in)        :: Node, Nodes

    integer :: iq, nu, lam, lp
    real    :: omega_qnu, n_qnu, dE, wt
    complex, allocatable :: gex(:,:)
    complex, allocatable :: strip_c(:,:), dressed_c(:,:), ket_c(:,:)
    complex, allocatable :: strip_v(:,:), dressed_v(:,:), ket_v(:,:)
    integer, allocatable :: pair_idx(:,:,:)
    real    :: Gamma_loc(dimbse)
    real, parameter :: pi = 3.14159265358979
#ifdef MPI
    integer :: MPIError
    real    :: Gamma_tot(dimbse)
#endif

    Gamma_loc = 0.0
    allocate(gex(dimbse,dimbse))
    allocate(strip_c(nc,dimbse), dressed_c(nc,dimbse), ket_c(nc,dimbse))
    allocate(strip_v(nv,dimbse), dressed_v(nv,dimbse), ket_v(nv,dimbse))
    allocate(pair_idx(nv,nc,ngkpt))
    call build_pair_idx(stt_bse, dimbse, nv, nc, ngkpt, pair_idx)

    do iq = 1 + Node, eph%nq, Nodes   ! round-robin over q

      do nu = 1, eph%nmodes
        omega_qnu = eph%omega(nu,iq)
        n_qnu     = bose_einstein(omega_qnu, temp)

        call compute_gex(eph, iq, nu, vector, w90basis, nc, nv, ngkpt, nc+nv, &
                         hbse, dimbse, pair_idx, kplusq_idx, &
                         strip_c, dressed_c, ket_c, strip_v, dressed_v, ket_v, gex)

        !$omp parallel do private(lp,dE,wt) schedule(static)
        do lam = 1, dimbse
          do lp = 1, dimbse
            dE = W(lam) - W(lp)
            wt = (n_qnu+1.0)*lorentzian(dE-omega_qnu,eta) &
               +  n_qnu      *lorentzian(dE+omega_qnu,eta)
            Gamma_loc(lam) = Gamma_loc(lam) + abs(gex(lam,lp))**2 * wt
          end do
        end do
        !$omp end parallel do

      end do  ! nu
    end do  ! iq

    Gamma_loc = Gamma_loc * pi

    deallocate(gex, strip_c, dressed_c, ket_c, strip_v, dressed_v, ket_v, pair_idx)

#ifdef MPI
    call MPI_ALLREDUCE(Gamma_loc, Gamma, dimbse, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, MPIError)
#else
    Gamma = Gamma_loc
#endif

end subroutine exciton_phonon_linewidth


!---------------------------------------------------------------------------
subroutine exciton_phonon_phassist(eph, vector, w90basis, nc, nv, ngkpt, &
           hbse, stt_bse, dimbse, W, mu_x, mu_y, mu_z, kplusq_idx, &
           temp, eta, omega_grid, nomega, chi2_xx, chi2_yy, chi2_zz, Node, Nodes)
    ! Incoherent phonon-assisted absorption spectrum (gauge-invariant):
    !   chi2_xx(omega) += sum_{lam,lp,q,nu}
    !     |gex(lam,lp)|^2 * |mu_x(lp)|^2
    !     * [(n+1)*L(omega-W(lam)-omega_qnu,eta) + n*L(omega-W(lam)+omega_qnu,eta)]
    !
    ! |mu(lp)|^2 is the optical weight of the initial bright exciton lp.
    ! The incoherent formula is gauge-invariant under global phase rotations
    ! of exciton eigenvectors (unlike the coherent |gex.mu|^2 form).

#ifdef MPI
    use mpi
#endif

    type(eph_data), intent(in) :: eph
    integer, intent(in)        :: w90basis, nc, nv, ngkpt, dimbse, nomega
    complex, intent(in)        :: vector(w90basis,nc+nv,ngkpt)
    complex, intent(in)        :: hbse(dimbse,dimbse)
    integer, intent(in)        :: stt_bse(4,dimbse)
    real,    intent(in)        :: W(dimbse)
    complex, intent(in)        :: mu_x(dimbse), mu_y(dimbse), mu_z(dimbse)
    integer, intent(in)        :: kplusq_idx(ngkpt,eph%nq)
    real,    intent(in)        :: temp, eta
    real,    intent(in)        :: omega_grid(nomega)
    real,    intent(out)       :: chi2_xx(nomega), chi2_yy(nomega), chi2_zz(nomega)
    integer, intent(in)        :: Node, Nodes

    integer :: iq, nu, lam, lp, iomega
    real    :: omega_qnu, n_qnu, gex2, muxx, muyy, muzz, dE_em, dE_ab, Lval
    complex, allocatable :: gex(:,:)
    complex, allocatable :: strip_c(:,:), dressed_c(:,:), ket_c(:,:)
    complex, allocatable :: strip_v(:,:), dressed_v(:,:), ket_v(:,:)
    integer, allocatable :: pair_idx(:,:,:)
    real    :: chi2_xx_loc(nomega), chi2_yy_loc(nomega), chi2_zz_loc(nomega)
    real    :: mu2_x(dimbse), mu2_y(dimbse), mu2_z(dimbse)
#ifdef MPI
    integer :: MPIError
#endif

    chi2_xx_loc = 0.0;  chi2_yy_loc = 0.0;  chi2_zz_loc = 0.0
    allocate(gex(dimbse,dimbse))
    allocate(strip_c(nc,dimbse), dressed_c(nc,dimbse), ket_c(nc,dimbse))
    allocate(strip_v(nv,dimbse), dressed_v(nv,dimbse), ket_v(nv,dimbse))
    allocate(pair_idx(nv,nc,ngkpt))
    call build_pair_idx(stt_bse, dimbse, nv, nc, ngkpt, pair_idx)

    mu2_x = real(mu_x * conjg(mu_x))
    mu2_y = real(mu_y * conjg(mu_y))
    mu2_z = real(mu_z * conjg(mu_z))

    do iq = 1 + Node, eph%nq, Nodes

      do nu = 1, eph%nmodes
        omega_qnu = eph%omega(nu,iq)
        n_qnu     = bose_einstein(omega_qnu, temp)

        call compute_gex(eph, iq, nu, vector, w90basis, nc, nv, ngkpt, nc+nv, &
                         hbse, dimbse, pair_idx, kplusq_idx, &
                         strip_c, dressed_c, ket_c, strip_v, dressed_v, ket_v, gex)

        ! Incoherent accumulation; OMP reduction on fixed-size stack arrays.
        !$omp parallel do private(lp,gex2,muxx,muyy,muzz,iomega,dE_em,dE_ab,Lval) &
        !$omp&            reduction(+:chi2_xx_loc,chi2_yy_loc,chi2_zz_loc) schedule(static)
        do lam = 1, dimbse
          do lp = 1, dimbse
            gex2 = abs(gex(lam,lp))**2
            muxx = gex2 * mu2_x(lp)
            muyy = gex2 * mu2_y(lp)
            muzz = gex2 * mu2_z(lp)
            do iomega = 1, nomega
                dE_em = omega_grid(iomega) - W(lam) - omega_qnu
                dE_ab = omega_grid(iomega) - W(lam) + omega_qnu
                Lval  = (n_qnu+1.0)*lorentzian(dE_em,eta) + n_qnu*lorentzian(dE_ab,eta)
                chi2_xx_loc(iomega) = chi2_xx_loc(iomega) + muxx * Lval
                chi2_yy_loc(iomega) = chi2_yy_loc(iomega) + muyy * Lval
                chi2_zz_loc(iomega) = chi2_zz_loc(iomega) + muzz * Lval
            end do
          end do
        end do
        !$omp end parallel do

      end do  ! nu
    end do  ! iq

    deallocate(gex, strip_c, dressed_c, ket_c, strip_v, dressed_v, ket_v, pair_idx)

#ifdef MPI
    call MPI_ALLREDUCE(chi2_xx_loc, chi2_xx, nomega, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, MPIError)
    call MPI_ALLREDUCE(chi2_yy_loc, chi2_yy, nomega, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, MPIError)
    call MPI_ALLREDUCE(chi2_zz_loc, chi2_zz, nomega, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, MPIError)
#else
    chi2_xx = chi2_xx_loc
    chi2_yy = chi2_yy_loc
    chi2_zz = chi2_zz_loc
#endif

end subroutine exciton_phonon_phassist


end module bse_eph
