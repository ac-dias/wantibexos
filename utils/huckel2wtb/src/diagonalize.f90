! -----------------------------------------------------------------------------
! MODULE DIAGONALIZE: In this module we have the subroutines that diagonalizes 
! The system hamiltonian. The subroutine were taken from the numerical recipes
!
! (1) ENERGIES: calculates the band structure of a material.
!
! (2) KPOINTS: Generates the k-points according with the bravais lattice for 
!              Band-Structure Calculations
!
! LAST MODIFICATION: 06/07/2025
! -----------------------------------------------------------------------------

module diagonalize
use overlaps_jc
implicit none

public:: energies, kpoints

CONTAINS

!--------------------------------------------------------------------------------
! SUBROUTINE ENERGIES: Calculates the hamiltonian eigenvalues in a previously
!                      readen k-path from de BANDS_LINES.IN input file
!--------------------------------------------------------------------------------

subroutine energies(ncell,numtp,nat,tnao,norbv,itype,rx,ry,rz,pv1,pv2,pv3,oe,esit,keht,cf,zeta,&
                 sip,ldm,dbov,lptb,nzt,nval,lval,eshift,nelect,nkpt,kpoint,hamk,ovlk,eband,rcut,&
                 noso,maxnn,nbnn,lisnn,posnn,iref,flp,nsybz,kpcoord,hspbz)
implicit none

! ----------------- Variables and Parameters of the Subroutine ----------------

  ! Externals:

    integer, parameter  :: lmax = 2, nmax = 7 
    integer, intent(in) :: ncell(3), numtp, nat, tnao, norbv(numtp), itype(nat)
    integer, intent(in) :: sip(nat), ldm(0:lmax), dbov(numtp), lptb(numtp,0:lmax)
    integer, intent(in) :: nzt(numtp,3), nval(numtp,3), lval(numtp,3)
    integer, intent(in) :: nelect, nkpt, noso, flp, nsybz
    integer, intent(in) :: maxnn, nbnn(nat), lisnn(nat,maxnn), iref(nat)
    real*8 , intent(in) :: oe(numtp,9), esit(numtp,3), cf(numtp,3,2), zeta(numtp,3,2)
    real*8 , intent(in) :: keht(numtp,numtp), rx(nat), ry(nat), rz(nat)
    real*8 , intent(in) :: pv1(3), pv2(3), pv3(3), kpoint(nkpt,3), rcut
    real*8 , intent(in) :: posnn(nat,maxnn,3), kpcoord(nkpt)
    real*8 , intent(inout) :: eband(nkpt,noso)
    complex*16, intent(inout) :: hamk(nkpt,noso,noso), ovlk(nkpt,noso,noso)
    character*5, intent(in)   :: hspbz(nsybz)

  ! Internals:

    integer, parameter   :: nospc = 9
    real*8, parameter    :: pi = 4.d0*atan(1.d0)
    integer              :: i, j, k, jj, oi, oj, ct, ic, ix, iy, iz
    integer              :: jo, jd, io, id, li, lj, idmax
    integer              :: isoc, jsoc, imax(1), imin(1)
    integer              :: izt, izoi, izoj, ik, idv, ikcount, mdcount, testdim
    integer              :: llsel(0:lmax,0:lmax), nnsel(nmax,nmax)
    real*8               :: ovout(25), ui(3), uj(3), rij(3), kp(nkpt,3), eshift
    real*8               :: kpb(3), kdotr, kh, delt, emax, emin, egap, ucvec(3)
    complex*16           :: hk(tnao,tnao), sk(tnao,tnao)

  ! Parametros para a subrotina de diagonalizacao (LAPACK, ZHEGVD):

    integer                  :: itp, ndm, lda, ldb, lwork, lrwork, liwork, info
    complex*16, allocatable  :: A(:,:), B(:,:), work(:)
    real*8, allocatable      :: w(:), rwork(:), hop(:,:,:,:), sov(:,:,:,:), orben(:,:)
    integer, allocatable     :: iwork(:)
    character*1              :: jobz, uplo
! ------------------------------------------------------------------------

  ! ========================================================================     
  !    Setting the parameters for the lapack diagonalization subroutine:
  ! ========================================================================     

    ndm = tnao
    itp = 1
     jobz = 'V'
      uplo = 'U'
       lda = ndm 
       ldb = ndm
      lwork = 2*ndm + ndm*ndm
     lrwork = 1 + 5*ndm + 2*ndm*ndm
    liwork = 3 + 5*ndm
    allocate(A(ndm,ndm),B(ndm,ndm),W(ndm),work(lwork),rwork(lrwork),iwork(liwork))

  ! ========================================================================     
  !              BUILDING THE K-INDEPENDENT HOP and SOV MATRICES:
  ! ========================================================================     

    idmax = maxval(dbov(:))
    allocate(hop(nat,idmax,maxnn,idmax), sov(nat,idmax,maxnn,idmax))
    allocate(orben(numtp,idmax))

  ! Loading the orbital energy array:

    do i = 1, numtp
     do j = 1, dbov(i)
      orben(i,j) = oe(i,j)
     enddo
    enddo

  ! Building the HOP matrix, which saves all the hoppings:

    do i = 1, nat                ! Loop over the atoms in the central cell:

      do jj = 1, nbnn(i)        ! Loop over the J atoms

          j = lisnn(i,jj) 

          if(jj == iref(i)) then
            do k = 1, dbov(itype(i))
             hop(i,k,jj,k) = orben(itype(i),k)
             sov(i,k,jj,k) = 1.d0
            enddo

          else          

       rij(1) = posnn(i,jj,1)
       rij(2) = posnn(i,jj,2)
       rij(3) = posnn(i,jj,3)

  ! Loop over the I and J's orbitals:

        do oi = 1, norbv(itype(i))   ! Loop over the I's orbitals:
          do oj = 1, norbv(itype(j))   ! Loop over the J's orbitals:

            li = lval(itype(i),oi)
            lj = lval(itype(j),oj)

  ! Overlap Calculation:

       call calcov_jc(i,j,oi,oj,nat,numtp,nzt,nval,lval,itype,zeta,cf,rij,ovout)

       ct = 0
       do jo = 1, ldm(lj)
        jd = lptb(itype(j),lj) + jo - 1

         do io = 1, ldm(li)
          id = lptb(itype(i),li) + io - 1

           ct = ct + 1

           kh = keht(itype(i),itype(j))
           kdotr = dot_product(kpb,rij)
           hop(i,id,jj,jd) = 0.5d0*kh*(esit(itype(i),oi) + esit(itype(j),oj))*ovout(ct)
           sov(i,id,jj,jd) = ovout(ct)

         enddo
       enddo

          enddo         ! End of the loop over the J's orbitals
        enddo         ! End of the loop over the I's orbitals

          endif

      enddo         ! End of the Loop over the J atoms
    enddo         ! End of the loop over the I atoms

  ! Apllying the shift to the hamiltoniana matrix elements:

    hop(:,:,:,:) = hop(:,:,:,:) + eshift*sov(:,:,:,:)

  ! Writing the Hamiltonian for the WantiBexos code:

    call writeham(nat,ncell,tnao,numtp,dbov,itype,sip,maxnn,idmax,&
                  nbnn,lisnn,posnn,hop,sov,pv1,pv2,pv3)

! ------------------------------------------------------------------------
!                    DIAGONALIZATION OF THE HAMILTONIAN
! ------------------------------------------------------------------------

    ikcount = 1

    do ik = 1, nkpt     

    kpb(1) = kpoint(ik,1)
    kpb(2) = kpoint(ik,2)
    kpb(3) = kpoint(ik,3)

  ! Initializing the arrays SK and HK:

    sk(:,:) = cmplx(0.0d0,0.0d0)
    hk(:,:) = cmplx(0.0d0,0.0d0)

    do i = 1, nat                ! Loop over the atoms in the central cell:
      do jj = 1, nbnn(i)            ! Loop over the J atoms

           j = lisnn(i,jj)

           ! MAIN DIAGONAL of the HK and SK matrices:

            if(jj == iref(i)) then 
              do k = 1, dbov(itype(i))
               ic = sip(i) + k
               hk(ic,ic) = hk(ic,ic) + cmplx(hop(i,k,jj,k))
               sk(ic,ic) = sk(ic,ic) + cmplx(sov(i,k,jj,k))
              enddo

            else

       rij(1) = posnn(i,jj,1)
       rij(2) = posnn(i,jj,2)
       rij(3) = posnn(i,jj,3)

       

  ! Loop over the I and J's orbitals:

        do oi = 1, dbov(itype(i))   ! Loop over the I's orbitals:
          do oj = 1, dbov(itype(j))   ! Loop over the J's orbitals:

  ! Hamiltonian and overlap matrices:

       id = sip(i) + oi
       jd = sip(j) + oj
       kdotr = dot_product(kpb,rij)
       hk(id,jd) = hk(id,jd) + hop(i,oi,jj,oj)*cmplx(cos(kdotr),sin(kdotr))
       sk(id,jd) = sk(id,jd) + sov(i,oi,jj,oj)*cmplx(cos(kdotr),sin(kdotr))

          enddo         ! End of the loop over the J's orbitals
        enddo         ! End of the loop over the I's orbitals

            endif

      enddo         ! End of the Loop over the J atoms
    enddo         ! End of the loop over the I atoms

    A(:,:) = hk(:,:)
    B(:,:) = sk(:,:)

     ovlk(ikcount,:,:) = B(:,:)   ! CARE HERE: The values are saved before the diagonalization

    call zhegvd(itp,jobz,uplo,ndm,A,lda,B,ldb,W,work,lwork,rwork,lrwork,IWORK,LIWORK,INFO)

  ! Storing the bands in EBAND, HK(k) and OVLK(k) to calculate DOS:

    eband(ikcount,:)   = W(:)
     hamk(ikcount,:,:) = A(:,:)
    ikcount = ikcount + 1

    W = 0.0d0
    A = 0.0d0
    B = 0.0d0

    enddo         ! End of the loop over the k-points

    deallocate(A,B)

  ! ---------------------------------------------------------------------------
  !  POST PROCESSING: writing the band structure
  ! ---------------------------------------------------------------------------

  ! Writing the bands:

    open(9, file='BandsHuckel.dat', status ='unknown')

    do ik = 1, nkpt
      write(9,'(f5.2,10f10.4)') kpcoord(ik), (eband(ik,j), j = flp - 3, flp + 4)
    enddo

  ! Printing the GAP and the K-point at the minimum:

    imax = maxloc(eband(:,flp))
    imin = minloc(eband(:,flp+1))
    emin = maxval(eband(:,flp))
    emax = minval(eband(:,flp+1))
    egap = emax - emin

  write(*,'(a,1f10.4)') '  Gap Energy                 =', egap
  write(*,'(a,1f10.4)') '  Fermi Energy               =', eband(imax,flp)
  write(*,'(a,1i6   )') '  IMAX =                     =', imax
  write(*,'(a,3f10.4)') '  K at the Maximum of the VB =', kpoint(imax,1), kpoint(imax,2), kpoint(imax,3)
  write(*,'(a,3f10.4)') '  K at the Minimum of the CB =', kpoint(imin,1), kpoint(imin,2), kpoint(imin,3)
  print* 

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

! ------------------------------------------------------------------------------
!                          End of the subroutine
! ------------------------------------------------------------------------------

                         end subroutine energies

! ------------------------------------------------------------------------------
subroutine kpoints(nkpt,kcoord,pi,kp,pv1,pv2,pv3)
implicit none

! ----------------- Variables and Parameters of the Subroutine ----------------
! Externals:
integer, intent(in) :: nkpt
real*8,  intent(in) :: kcoord(nkpt,3), pi, pv1(3), pv2(3), pv3(3)
real*8,  intent(out):: kp(nkpt,3)    
! Internals
real*8              :: b1(3), b2(3), b3(3), vuc(3)
integer             :: i, j
!-----------------------------------------------------------------------

! Volum of the primitive cell

vuc(1) = pv1(1)*(pv2(2)*pv3(3) - pv2(3)*pv3(2)) + pv1(2)*(pv2(3)*pv3(1) - pv2(1)*pv3(3)) + &
         pv1(3)*(pv2(1)*pv3(2) - pv2(2)*pv3(1))

vuc(2) = pv2(1)*(pv3(2)*pv1(3) - pv3(3)*pv1(2)) + pv2(2)*(pv3(3)*pv1(1) - pv3(1)*pv1(3)) + &
         pv2(3)*(pv3(1)*pv1(2) - pv3(2)*pv1(1))

vuc(3) = pv3(1)*(pv1(2)*pv2(3) - pv1(3)*pv2(2)) + pv3(2)*(pv1(3)*pv2(1) - pv1(1)*pv2(3)) + &
         pv3(3)*(pv1(1)*pv2(2) - pv1(2)*pv2(1))

    ! Primitive vectors of the reciprocal lattice

  b1(1) = 2.*pi*(pv2(2)*pv3(3) - pv2(3)*pv3(2))/vuc(1)
  b1(2) = 2.*pi*(pv2(3)*pv3(1) - pv2(1)*pv3(3))/vuc(1)
  b1(3) = 2.*pi*(pv2(1)*pv3(2) - pv2(2)*pv3(1))/vuc(1)

  b2(1) = 2.*pi*(pv3(2)*pv1(3) - pv3(3)*pv1(2))/vuc(2)
  b2(2) = 2.*pi*(pv3(3)*pv1(1) - pv3(1)*pv1(3))/vuc(2)
  b2(3) = 2.*pi*(pv3(1)*pv1(2) - pv3(2)*pv1(1))/vuc(2)

  b3(1) = 2.*pi*(pv1(2)*pv2(3) - pv1(3)*pv2(2))/vuc(3)
  b3(2) = 2.*pi*(pv1(3)*pv2(1) - pv1(1)*pv2(3))/vuc(3)
  b3(3) = 2.*pi*(pv1(1)*pv2(2) - pv1(2)*pv2(1))/vuc(3)

    do i = 1, nkpt
     kp(i,1) = kcoord(i,1)*b1(1) + kcoord(i,2)*b2(1) + kcoord(i,3)*b3(1)
      kp(i,2) = kcoord(i,1)*b1(2) + kcoord(i,2)*b2(2) + kcoord(i,3)*b3(2)
     kp(i,3) = kcoord(i,1)*b1(3) + kcoord(i,2)*b2(3) + kcoord(i,3)*b3(3)
    enddo

  ! End of the subroutine

    end subroutine kpoints

  !====================================================================
  ! SUBROUTINE WRITEHAM: writes the the Huckel hamiltonian in a format 
  !                      which can be readen for the WANTIBEXOS program                         
  !====================================================================

    subroutine writeham(nat,ncell,tnao,numtp,dbov,itype,sip,maxnn,idmax,&
                        nbnn,lisnn,posnn,hop,sov,pv1,pv2,pv3)
    implicit none

  ! ------------- Variables and Parameters of the Subroutine -------------

    ! Externals:

    integer, intent(in) :: nat, ncell(3), tnao, numtp, dbov(numtp), itype(nat) 
    integer, intent(in) :: sip(nat), maxnn, idmax
    integer, intent(in) :: nbnn(nat), lisnn(nat,maxnn)
    real*8 , intent(in) :: posnn(nat,maxnn,3), hop(nat,idmax,maxnn,idmax)
    real*8 , intent(in) :: sov(nat,idmax,maxnn,idmax), pv1(3), pv2(3), pv3(3)

    ! Internals:

    integer              :: i, j, jj, io, jo, k, kk, ko, id, jd, kd
    integer              :: icc, ica, m, nctot, nlines
    integer, allocatable :: atlab(:)
    real*8               :: reh, imh, ovl, rij(3), zero
    real*8, allocatable  :: cvec(:,:), hamwtbx(:,:,:), ovlwtbx(:,:,:)

  !--------------------------------------------------------------------

   open(71,file='tb_hr.dat', status='unknown')
   zero = 0.0d0

               ! **************************************
               !       Generating the Cell Vectors
               ! **************************************

        nctot = (2*ncell(1) + 1)*(2*ncell(2) + 1)*(2*ncell(3) + 1)
        allocate(cvec(nctot,3), atlab(maxnn))

        icc = 1
        ica = 1
        do i = -ncell(1), ncell(1)
          do j = -ncell(2), ncell(2)
            do k = -ncell(3), ncell(3)
              cvec(icc,1) = i*pv1(1) + j*pv2(1) + k*pv3(1)
              cvec(icc,2) = i*pv1(2) + j*pv2(2) + k*pv3(2)
              cvec(icc,3) = i*pv1(3) + j*pv2(3) + k*pv3(3)
                do m = 1, nat
                  atlab(ica) = m
                  ica = ica + 1
                enddo
              icc = icc + 1
            enddo
          enddo
        enddo

               ! **************************************
               !  Huckel to WantiBexos data Structure
               ! **************************************

       allocate(hamwtbx(nctot,tnao,tnao), ovlwtbx(nctot,tnao,tnao))

       do i = 1, nctot
         do j = 1, nat
           do jo = 1, dbov(itype(j))
            jd = sip(j) + jo
             do k = 1, nat
               do ko = 1, dbov(itype(k))

                 kd = sip(k) + ko
                 kk = nat*(i - 1) + k
                 hamwtbx(i,jd,kd) = hop(j,jo,kk,ko)
                 ovlwtbx(i,jd,kd) = sov(j,jo,kk,ko)

               enddo
             enddo
           enddo
         enddo
       enddo

               ! **************************************
               !         Writing the Hamiltonian
               ! **************************************

   nlines = nctot*tnao*tnao

   write(71,*) 'NP'
   write(71,*) zero
   write(71,*) zero
   write(71,'(3f14.8)') (pv1(j), j = 1, 3)
   write(71,'(3f14.8)') (pv2(j), j = 1, 3)
   write(71,'(3f14.8)') (pv3(j), j = 1, 3)
   write(71,*) tnao
   write(71,*) nctot
   write(71,'(a)') '#rcell x   rcell y   rcell z   i   j   ReH   ImH   S'

        do i = 1, nctot
          do j = 1, tnao
            do k = 1, tnao
              write(71,'(3f10.4,2i8,3f14.8)') cvec(i,1), cvec(i,2), cvec(i,3), j, k, hamwtbx(i,j,k), zero, ovlwtbx(i,j,k)
            enddo
          enddo
        enddo

  !====================================================================
                     end subroutine writeham
  !====================================================================

! -----------------------------------------------------------------------------
                        end module diagonalize
! -----------------------------------------------------------------------------
