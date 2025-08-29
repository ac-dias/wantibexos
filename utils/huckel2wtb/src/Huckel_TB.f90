! -----------------------------------------------------------------------------
!                               Program HUCKEL_TB
! -----------------------------------------------------------------------------
! Simplified version of the Huckel code from A. S. Martins, which calculates the
! Band structure of a give material and generates the file HAMIL.DAT for the
! WantiBexos code. Some remarks:
! 
! (1) Uses double zeta basis for representing each atomic orbital
! 
! (2) The overlaps among the orbitals are calculated numerically
!
! (3) All the input data is readen from a just single file
!
! (4) SOC was not included in this version
!
! -----------------------------------------------------------------------------
!  Date: 24/08/2025
! -----------------------------------------------------------------------------

    program huckel_tb
    use diagonalize

    implicit none

  !----------------------- Variables and Arrays ---------------------------------
  ! (i) Parameters:

   integer, parameter       :: nmax = 7, lmax = 2
   real*8, parameter        :: pi = 4.d0*atan(1.0d0) 

  ! (ii) From the input files:

   integer                  :: nat, nauc, numtp, ncell(3)
   real*8                   :: lpar(3), rcut
   integer, allocatable     :: zat(:), itp(:)
   integer, allocatable     :: norbv(:), nval(:,:), lval(:,:), nzt(:,:), neao(:,:)
   real*8 , allocatable     :: zeta(:,:,:), cf(:,:,:), esit(:,:), keht(:,:)
   real*8 , allocatable     :: kcoord(:,:)
   character*5, allocatable :: zsimb(:), hspbz(:)
   character*40             :: cltit
 
  ! (iii) Internals

   real*8,  allocatable     :: rx(:), ry(:), rz(:), oe(:,:), hr(:,:), sr(:,:), W(:)
   real*8                   :: eshift, vuc,  version
   real*8                   :: pv1(3), pv2(3), pv3(3)
   real*8                   :: ui(3), uj(3), rij(3), rdij
   integer, allocatable     :: itype(:), nngh(:), lisngh(:,:), dbov(:)
   integer, allocatable     :: ldm(:), sip(:), lptb(:,:), nspc(:), nbnn(:)
   integer, allocatable     :: numnn(:), lisnn(:,:), iref(:)
   real*8, allocatable      :: posnn(:,:,:)
   integer                  :: ni, nj, li, lj, llp, i, j, k, m, n, id, jd, ik, kk, today(3)
   integer                  :: ix, iy, iz, ib, iatom, noso, maxnn, icount, icnn
   integer                  :: nnn, tnao, nctot, ErrorFlag, nelect, imin(1), imax(1)

  ! (iv) Internals - Reciprocal space

   integer                  :: nsybz, ndiv, nlines, nkpt, idv, ikcount, flp
   real*8                   :: emin, emax, egap 
   real*8                   :: up, ur, us, b1(3), b2(3), b3(3)
   real*8, allocatable      :: eband(:,:)
   real*8, allocatable      :: kp(:,:), kpmod(:), kpmodacc(:), kpoint(:,:), kpcoord(:)
   complex*16, allocatable  :: hamk(:,:,:), ovlk(:,:,:)

  ! (v) NAMELISTS:

   namelist /globalvars/ nauc, numtp, nsybz, ndiv, cltit, eshift
   namelist /unitcell/ ncell, lpar, pv1, pv2, pv3, itype, rx, ry, rz
   namelist /orbitals/ zat, zsimb, norbv, nval, lval, neao, rcut, keht, zeta, cf, esit
   namelist /kpath/ kcoord, hspbz

  ! -----------------------------------------------------------------------------

    version = 1.0

  print*
  print*, "--------------------- PROGRAM HUCKEL2WTB -----------------------------"
  write(*,'(a,f3.1)') '     VERSION: ', version
  print*, "----------------------------------------------------------------------"
  print*
  print*, "Calculates the electronic structure of a crystalline material in the  "
  print*, "Tight-Binding aproximation. The radial part of the atomic orbitals are" 
  print*, "Described  as linear combination  of Slater-Type Functions (STO).  The"
  print*, "Elements of the Hamiltonian matrix follows the Huckel prescription: the" 
  print*, "Overlaps among the orbitals are explicitely computed."
  print* 
  print*, "Unit System:"
  print*
  print*, "DISTANCES: Angstrons (A)"
  print*, " ENERGIES: electron volts (eV)"
  print*, "   FORCES: eV/Angs"
  print*
  print*, "In the current version  of the program, the maximum allowed values for"
  print*, "The L quantum numbers is LMAX = 2 (s, p and d)"

  ! -----------------------------------------------------------------------------
  !                (I) READING THE DATA FROM THE INPUT FILE                                      
  ! -----------------------------------------------------------------------------
  ! All data is stored in just one input file, grouped in different namelists
  ! -----------------------------------------------------------------------------

   ! (a) Global Variables:

     read(*,NML=globalvars)
     nat = nauc

   ! (b) Unit Cell Related variables:

     allocate(itype(nat), rx(nat), ry(nat), rz(nat))
     read(*,NML=unitcell) 

   ! Scaling the lattice Primitive Vectors:

     pv1(:) = pv1(:)*lpar(1)
     pv2(:) = pv2(:)*lpar(2)
     pv3(:) = pv3(:)*lpar(3)

   ! Volum of the primitive cell

     vuc = pv1(1)*(pv2(2)*pv3(3) - pv2(3)*pv3(2)) + pv1(2)*(pv2(3)*pv3(1) - pv2(1)*pv3(3)) + &
           pv1(3)*(pv2(1)*pv3(2) - pv2(2)*pv3(1))

     vuc = dabs(vuc)   ! For non-orthogonal primitive vectors

   ! (c) Orbital's and Huckel's parameters:

     allocate(zat(numtp), zsimb(numtp), norbv(numtp))
     allocate(nval(numtp,3), lval(numtp,3),neao(numtp,3))
     allocate(keht(numtp,numtp))
     allocate(nzt(numtp,3))        ! 1 = s, 2 = p and 3 = d
     allocate(zeta(numtp,3,2), cf(numtp,3,2))
     allocate(esit(numtp,3))
     nzt(:,:) = 2                  ! Double Zeta basis

     read(*,NML=orbitals)  

  ! Off Diagonal values of KEHT:

    do i = 1, numtp - 1
      do j = i + 1, numtp
        keht(i,j) = dsqrt(keht(i,i)*keht(j,j))    !0.50d0*(keht(i,i) + keht(j,j))
        keht(j,i) = keht(i,j)
      enddo
    enddo

    allocate(kcoord(nsybz,3), kp(nsybz,3), kpmod(nsybz), kpmodacc(nsybz))
    allocate(hspbz(nsybz))

   ! (d) Kpath:

    read(*,NML=kpath)  

  ! ----------------------------------------------------------------------------
  !                (II) SETTING THE INITIAL VARIABLES                           
  ! ----------------------------------------------------------------------------
  ! Defining the vectors related to the basis orbitals: LDM(0:lmax) gives the
  ! Number of orbitals for each L value and DBOV(NUMTP) gives the total number of
  ! Orbitals for each specie of the system:
  ! ----------------------------------------------------------------------------

     allocate(ldm(0:lmax), dbov(numtp))

     ldm(0) = 1
     ldm(1) = 3
     ldm(2) = 5

  ! Assigning the total number of valence orbitals according to the basis dimension:

     dbov = 0
     do i = 1, numtp
      do j = 1, norbv(i)
       dbov(i) = dbov(i) + ldm(lval(i,j))
      enddo
     enddo

  ! Pointer for each L value inside the basis: for each specie, this pointer
  ! Gives the initial position in the basis for the orbital of L quantum number.

     allocate(lptb(numtp,0:lmax))

     do i = 1, numtp

     llp = norbv(i)
     select case (llp)
     case(1) ! Basis with only S orbital
     lptb(i,0) = 1
     case(2) ! Basis with S and P orbitals
     lptb(i,0) = 1
     lptb(i,1) = 2
     case(3) ! Basis with S, P and D orbitals
     lptb(i,0) = 1
     lptb(i,1) = 2
     lptb(i,2) = 5
     end select

     enddo

  ! ------------------------------------------------------------------------------
  !                 (III) GENERATION OF NEIGHBOHR's RELATED ARRAYS                
  ! ------------------------------------------------------------------------------
  ! These arrays are needed for a efficient construction of the matrix elements
  ! Of the hamiltonian.
  ! ------------------------------------------------------------------------------

    allocate(nbnn(nat), iref(nat))
    nctot = (2*ncell(1) + 1)*(2*ncell(2) + 1)*(2*ncell(3) + 1)
    nbnn(:) = nat*nctot

            ! ---------------------------------------------------
            ! (iii.a) CALCULATING MAXNN         
            ! ---------------------------------------------------

    do i = 1, nat       ! Loop over the atoms in the central cell:
    icount = 0

    ! Loop over the unit cells:

        do ix = -ncell(1), ncell(1)
          do iy = -ncell(2), ncell(2)
            do iz = -ncell(3), ncell(3)

            do j = 1, nat                ! Loop over the J atoms

            icount = icount + 1
            if((j == i).and.(ix == 0).and.(iy == 0).and.(iz == 0)) iref(i) = icount

      ! Relative position

            enddo   ! End of the J atoms

            enddo
          enddo
        enddo

    enddo   ! End of the I atoms

            ! ---------------------------------------------------
            ! (iii.b) LOADING LISNN and POSNN
            ! ---------------------------------------------------

    maxnn = nat*nctot
    allocate(lisnn(nat,maxnn), posnn(nat,maxnn,3))

    do i = 1, nat       ! Loop over the atoms in the central cell:

    icount = 1
    ui(1) = rx(i)
    ui(2) = ry(i)
    ui(3) = rz(i)

    ! Loop over the unit cells:

        do ix = -ncell(1), ncell(1)
          do iy = -ncell(2), ncell(2)
            do iz = -ncell(3), ncell(3)
 
              do j = 1, nat                ! Loop over the J atoms

              uj(1) = rx(j) + ix*pv1(1) + iy*pv2(1) + iz*pv3(1)
              uj(2) = ry(j) + ix*pv1(2) + iy*pv2(2) + iz*pv3(2)
              uj(3) = rz(j) + ix*pv1(3) + iy*pv2(3) + iz*pv3(3)

      ! Relative position

              rij(1) = uj(1) - ui(1)
              rij(2) = uj(2) - ui(2)
              rij(3) = uj(3) - ui(3)

              rdij = dsqrt(dot_product(rij,rij))
                  lisnn(i,icount) = j
                posnn(i,icount,1) = rij(1)
                posnn(i,icount,2) = rij(2)
                posnn(i,icount,3) = rij(3)
                 icount = icount + 1

              enddo   ! End of the J atoms

            enddo
          enddo
        enddo

    enddo   ! End of the I atoms

! ------------------------------------------------------------------------------
!              Calculating the total number of electrons NELECT:
! ------------------------------------------------------------------------------

    nelect = 0

    do i = 1, nat
      do k = 1, norbv(itype(i))
        nelect = nelect + neao(itype(i),k)
      enddo
    enddo

! -----------------------------------------------------------------------------
!             Printing the information readen from the input file:
! -----------------------------------------------------------------------------

print*
print*, "----------------------------------------------------------------------"
print*, "          SUMMARY OF THE DATA READ FROM THE INPUT FILE                "
print*, "----------------------------------------------------------------------"
print*
print*, " SYSTEM PARAMETERS:"
print*
write(*,'(a,a)')    ' System                                    =  ',cltit
write(*,'(a,i8)')   ' Number of atoms in the unit cell (NAUC)   =', nauc
write(*,'(a,i8)')   ' Number of atoms (NAT)                     =', nat
write(*,'(a,i8)')   ' Number of atom types (NUMTP)              =', numtp
write(*,'(a,f8.3)') ' Lattice Parameter (LPAR)                  =', lpar(1)
write(*,'(a,f8.3)') ' C/A Ratio (COVA)                          =', lpar(3)/lpar(1)
write(*,'(a,f12.3)') ' Unit Cell Volum (VUC)                    =', vuc  
print*

print*, "  INFORMATION ABOUT THE ATOM TYPES:"
print*
write(*,'(a)')   '(a) Atomic numbers (ZAT) and Symbols:          '
print*
do k = 1, numtp
write(*,'(a,a,i1,a,a,i2,a,a,a)') ' ZAT','(',k,')',' = ',zat(k),' (',zsimb(k),')'
enddo
print*
write(*,'(a)')   '(b) Valence orbitals for each atom type:       '
print*
do k = 1, numtp
write(*,'(a,a,i1,a,a,i2,a,a,a)') ' NORBV','(',k,')',' = ', norbv(k)
enddo
print*
write(*,'(a)')   '(c) Number of Electrons of the Atomic Orbitals:       '
print*
do k = 1, numtp
write(*,'(a,a,a,a,a,a,a,a,a)') ' | ' , 'Specie' ,   ' | ' ,  '  (n,l)', '  | ', 'NEAO', ' |  ',  'ESIT', '   |  '
print*
  do j = 1, norbv(k)
write(*,'(i7,a,i1,a,i1,a,i7,f10.3)') k, '       (',nval(k,j),',',lval(k,j),')', neao(k,j), esit(k,j)
     enddo
print*
enddo
write(*,'(a,i6)') ' Total Number of Electrons: ', nelect
print*
write(*,'(a)')      '(d) Informations about the STOs for each specie:'
do i = 1, numtp
print*
write(*,'(a,a,a)') ' (',zsimb(i),'):'
print*
write(*,'(a,a,a,a,a,a,a,a)') ' | ' , 'Specie' ,   ' | ' ,  ' (n,l)', ' | ', 'ZETA', '  | ', ' CJ '
print*
  do j = 1, norbv(i)
     do k = 1, nzt(i,j)
write(*,'(i7,a,i1,a,i1,a,2f8.3)') i, '      (',nval(i,j),',',lval(i,j),')', zeta(i,j,k), cf(i,j,k)
     enddo
  enddo
enddo
print*

! -----------------------------------------------------------------------------
! Counting the total number of atoms of each specie:
! -----------------------------------------------------------------------------

  allocate(nspc(numtp))
  nspc(:) = 0

    do i = 1, nat
    j = itype(i)
      nspc(j) = nspc(j) + 1
    enddo

  ! Total number of atomic orbitals:

    tnao = 0
    do i = 1, numtp
      tnao = tnao + nspc(i)*dbov(i)
    enddo

    write(*,*) ' Dimension of the Basis (TNAO) =', tnao
  
  ! Creating the SIP pointer, that gives the number of atomic orbitals up the I
  ! Site.

  allocate(sip(nat))

    sip(1) = 0   

    do i = 2, nat
      sip(i) = sip(i-1) + dbov(itype(i-1))
    enddo

  ! Allocating the orbital energies according the specie:

    allocate(oe(numtp,9))

    do i = 1, numtp
     do j = 1, norbv(i)
      do k = 1, ldm(lval(i,j))
       id = lptb(i,lval(i,j)) + k - 1
       oe(i,id) = esit(i,j)
      enddo
     enddo
    enddo

! -----------------------------------------------------------------------------
!                             K-POINT GENERATION                               
! -----------------------------------------------------------------------------

    nlines = nsybz - 1

  ! Mapping the KCOORD vectors, in reduced units, in the KP arrays, in units
  ! Of Angs^(-1):

    call kpoints(nsybz,kcoord,pi,kp,pv1,pv2,pv3)

  ! Modules of the difference (kp(i+1) - kp(i)):

    kpmodacc(1) = 0.0d0

    do i = 1, nlines   
     kpmod(i) = dsqrt((kp(i+1,1) - kp(i,1))**2 + (kp(i+1,2) - kp(i,2))**2 + (kp(i+1,3) - kp(i,3))**2)
      kpmodacc(i+1) = kpmodacc(i) + kpmod(i)
    enddo

  ! Defining the total number of K-points:

    nkpt = nlines*ndiv + 1

  ! Allocating and generating the KPOINTS:

    allocate(kpoint(nkpt,3),kpcoord(nkpt))

    ikcount = 1 
    do ik = 1, nlines 
      do idv = 0, ndiv - 1
        kpoint(ikcount,1) = kp(ik,1) + idv*(kp(ik+1,1) - kp(ik,1))/dfloat(ndiv)
        kpoint(ikcount,2) = kp(ik,2) + idv*(kp(ik+1,2) - kp(ik,2))/dfloat(ndiv)
        kpoint(ikcount,3) = kp(ik,3) + idv*(kp(ik+1,3) - kp(ik,3))/dfloat(ndiv)
        kpcoord(ikcount) = kpmodacc(ik) + idv*kpmod(ik)/dfloat(ndiv)
        ikcount = ikcount + 1
      enddo
    enddo

    kpoint(nkpt,:) = kp(nsybz,:)
    kpcoord(nkpt) = kpmodacc(nlines) + kpmod(nlines)

  ! ---------------------- CALCULATION OF THE EIGENSTATES ------------------------
  ! The difference between KCALC = 1 and KCALC = 3,4 options is the k-point 
  ! Generation: in the former case we generate the k-points along the lines that
  ! Connect the high-symmetry points of the Brillouin zone and, in the later
  ! Case, the k-kpoints are generated according the Monkhost-Pack recipe. The
  ! Case KCALC = 2 just calculates the energies at specific K points (similar to
  ! Supercell calculations
  ! -----------------------------------------------------------------------------

  ! Allocating EBAND, HAMK and OVLK:

   noso = tnao
   allocate(eband(nkpt,tnao),hamk(nkpt,tnao,tnao),ovlk(nkpt,tnao,tnao))

     eband(:,:) = 0.0d0
    hamk(:,:,:) = cmplx(0.0d0,0.0d0)
    ovlk(:,:,:) = cmplx(0.0d0,0.0d0)

  ! File for the Band Structure and determining the HOMO level:          

    if(mod(nelect,2) == 0) then
     flp = nelect/2
      else
     flp = (nelect + 1)/2
    endif

  ! Band structure Calculation:    

  call energies(ncell,numtp,nat,tnao,norbv,itype,rx,ry,rz,pv1,pv2,pv3,oe,esit,keht,cf,zeta,&
                sip,ldm,dbov,lptb,nzt,nval,lval,eshift,nelect,nkpt,kpoint,hamk,ovlk,eband,rcut,&
                noso,maxnn,nbnn,lisnn,posnn,iref,flp,nsybz,kpcoord,hspbz)

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

  ! Deallocations:

    deallocate(zat, zsimb, norbv, nval, lval, neao, nzt)
    deallocate(zeta, cf, esit, oe, ldm, dbov, lptb)
    deallocate(rx, ry, rz, itype)

! -----------------------------------------------------------------------------
                          end program huckel_tb
! -----------------------------------------------------------------------------
