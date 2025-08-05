! -----------------------------------------------------------------------------
!                             MODULE OVERLAP
!
! This module contains all subroutines needed to calculate the overlap between
! The basis orbitals (S, P and D) and to build the Huckel Hamiltonian. The 
! Overlap subroutines were taken from J. Cerda's Green Code.
!
! LAST MODIFICATION: 27/01/2019
! -----------------------------------------------------------------------------

module overlaps_jc
implicit none

! Common Variables to the module:

public :: calcov_jc, mov, abfns, lovlap

contains

! -----------------------------------------------------------------------------
! Compute the overlap matrix among the orbitals of the atoms. In the symmetric
! Case, that is, for systems where the atoms are arranjed in a ordered array,
! The overlap matrix includes only the atoms that belong to the basis of the
! Unit cell. WARNING: double zeta basis and S, P and D orbitals only!!
! -----------------------------------------------------------------------------

subroutine calcov_jc(i,j,k,n,nat,numtp,nzt,nval,lval,itype,zeta,cf,rij,ovout)

!----------------------- Variables and Arrays ---------------------------------
! (I) Extern and parameters:
integer, parameter  :: lmax = 2
integer, intent(in) :: i, j, k, n, nat, numtp
integer, intent(in) :: nzt(numtp,3), nval(numtp,3), lval(numtp,3), itype(nat)
real*8,  intent(in) :: zeta(numtp,3,2), cf(numtp,3,2), rij(3)
! (II) Intern:
integer, parameter  :: nzzi = 2, nzzj = 2, nn = 4
integer             :: ni, nj, li, lj, iz, jz, llopt, nnopt, kj, lp
real*8, intent(out) :: ovout(25)
real*8              :: rho, ovlp(0:lmax), ovin(nn), fc
real*8              :: cl, cm, cn, cl2, cm2, cn2
real*8              :: zti(nzzi), ztj(nzzj), cofi(nzzi), cofj(nzzj)
! -----------------------------------------------------------------------------

  ! Initializing OVOUT:

    ovout(:) = 0.0d0

  ! Convertion factor for the distances:

    fc = 0.529177d0

  ! Assigning the ni and li quantum numbers of (I,K) orbital:

    ni  = nval(itype(i),k)
    li  = lval(itype(i),k)

  ! Assigning the nj and lj quantum numbers of (J,N) orbital:
  
    nj  = nval(itype(j),n)
    lj  = lval(itype(j),n)

! Distance |rj - ri| with PBC (in Bohrs):

        rho = dsqrt(rij(1)**2 + rij(2)**2 + rij(3)**2)/fc

  zti(1) = zeta(itype(i),k,1) 
  zti(2) = zeta(itype(i),k,2) 
  ztj(1) = zeta(itype(j),n,1) 
  ztj(2) = zeta(itype(j),n,2) 
 cofi(1) = cf(itype(i),k,1)
 cofi(2) = cf(itype(i),k,2)
 cofj(1) = cf(itype(j),n,1)
 cofj(2) = cf(itype(j),n,2)

! Overlap <ni,li|nj,lj>:

  ovin = 0.0d0
  call mov(ovin, rho, ni, nj, li, lj, nn, nzzi, zti, cofi, nzzj, ztj, cofj)

! Applying the Two-Center Slater and Koster Rules (SIGNALS TO BE VERYFIED):
! Direction cosines:
!
!   cl  : (xj-xi)/|rj-ri| = rij(1)/|rj - ri|
!   cm  : (yj-yi)/|rj-ri| = rij(2)/|rj - ri|
!   cn  : (zj-zi)/|rj-ri| = rij(3)/|rj - ri|
!   cl2,cm2,cn2: squares of the direction cosines cl,cm,cn

  cl = rij(1)/(rho*fc)
  cm = rij(2)/(rho*fc)
  cn = rij(3)/(rho*fc)
  cl2 = cl*cl
  cm2 = cm*cm
  cn2 = cn*cn

call twocenter(li,lj,ovin,ovout,cl,cl2,cm,cm2,cn,cn2)

! End of the subroutine

end subroutine calcov_jc


!  ******************************************************************** 
!  *                                                                 
!  *     SUBROUTINE MOV          CALLED FROM MOVLAP                   *
!  *                                                                  *
!  *                                                                  *
!  *       MOV      SUBROUTINE TO CALCULATE THE OVERLAP               *
!  *                COMPONENT INDEPENDENT OF THE ANGLE BETWEEN        *
!  *                THE ATOMS                                         *
!  *                                                                  *
!  *       SUBROUTINES USED:                                          *
!  *                                                                  *
!  *             LOVLAP AND ABFNS                                     *
!  *                                                                  *
!  *                                                                  *
!  *  WRITTEN BY CHARLES WILKER AUGUST 1982                           *
!  *    sigma= ov(1)
!  *    pi   = ov(2)
!  *    delta= ov(3)
!  *    phi  = ov(4)
!  *
!  *                                                                  *
!  ********************************************************************
subroutine mov( ov, r, na, nb, la, lb, nn, nzza, chia, cofa, nzzb, chib, cofb )
implicit none

integer , intent(in) :: na, nb, la, lb, nn, nzza, nzzb
real*8  , intent(in) :: r
real*8  , intent(in) :: cofa(nzza), cofb(nzzb), chia(nzza), chib(nzzb)
real*8  , intent(out):: ov(nn)
integer :: ij, maxn, izza, izzb, m
real*8  :: rll(4)

!     THIS ONLY WORKS FOR PRINCIPAL QUANTUM L < OR = 7
      real*8  :: a(30), b(30), ska, skb, zz

!---- factors required to conform SK tables
      real*8  :: fac_ll(3,0:2,0:2)
      data fac_ll /  1.0d0 , 1.0d0 , 1.0d0 ,& ! s-s
                    -1.0d0 , 1.0d0 , 1.0d0 ,& ! p-s
                     1.0d0 , 1.0d0 , 1.0d0 ,& ! d-s
                    -1.0d0 , 1.0d0 , 1.0d0 ,& ! s-p
                    -1.0d0 , 1.0d0 , 1.0d0 ,& ! p-p
                     1.0d0 ,-1.0d0 , 1.0d0 ,& ! d-p
                     1.0d0 , 1.0d0 , 1.0d0 ,& ! s-d
                     1.0d0 ,-1.0d0 , 1.0d0 ,& ! p-d
                     1.0d0 ,-1.0d0 , 1.0d0 /  ! d-d

!!    COMMON/LOCLAP/SK1,SK2,R,L1,L2,M,N1,N2,MAXN

!-----------------------------------------------------------------------

ov = 0.0 
maxn=na+nb

!  loop over orbitals

do izza = 1, nzza
 do izzb = 1, nzzb
  ska = chia(izza)
  skb = chib(izzb)

  call abfns( a, b, ska, skb, r, maxn )

  zz = cofa(izza) * cofb(izzb)
  do ij = 1 , nn
   m = ij - 1
   call lovlap( rll(ij) ,a, b, ska, skb, r, la, lb, m, na, nb  )
   ov( ij ) = ov( ij ) + zz*rll(ij)
  enddo
!!      write(iout,'(a,2i3,a,2f7.2,2(a,f7.2),a,20e9.2)') 
!!   &                   ' izza/b=',izza,izzb,
!!   &                   ' ska/b=',ska,skb,' zz=',zz,' r=',r,
!!   &                   ' rll=',rll(1:nn)
!!      write(iout,'(a,30e9.2)') ' a=',a(1:5)
!!      write(iout,'(a,30e9.2)') ' b=',b(1:5)

 enddo
enddo

!     required to match SK formulas
do ij = 1 , nn
 ov(ij) = ov(ij) * fac_ll(ij,la,lb)
enddo

return
end subroutine mov
!  ********************************************************************
!  *                                                                  *
!  *     SUBROUTINE ABFNS        CALLED FROM MOVLAP                   *
!  *                                                                  *
!  *                                                                  *
!  *       ABFNS    SUBROUTINE TO CALCULATE THE AB FUNCTIONS          *
!  *                                                                  *
!  *       SUBROUTINES USED:                                          *
!  *                                                                  *
!  *             NONE                                                 *
!  *                                                                  *
!  *                                                                  *
!  *       ORIGIN LOST IN ANTIQUITY                                   *
!  *                                                                  *
!  ********************************************************************
SUBROUTINE ABFNS(A,B, sk1, sk2, rr, maxcal )
implicit none

integer , intent(in)  :: maxcal
real*8  , intent(in)  :: rr, sk1, sk2
real*8  , intent(out) :: a(30), b(30)

integer  :: j, ix, ir, is, il, k, in, i
real*8   :: rho1, rho2, c, d, h, r, ra, rho22, t, tr

!!    COMMON/LOCLAP/SK1,SK2,RR,L1,L2,M,N1,N2,MAXCAL

a(:) = 0.0
b(:) = 0.0

J=MAXCAL+1
RHO1=0.5d0*(SK1+SK2)*RR
RHO2=0.5d0*(SK1-SK2)*RR
if ( abs(rho1).gt.165.0 .or. abs(rho2).gt.165.0 ) return

C=DEXP(-RHO1)
A(1)=C/RHO1
do i=2,j
 A(I)=(DBLE(FLOAT(I-1))*A(I-1)+C)/RHO1
enddo
IX=J
      IR=nint(DABS(2.0d0*RHO2))
      IS=MIN0(IR+1,19)
      IF(RHO2) 25,35,25
 25   D=DEXP(RHO2)
      H=1.0d0/D

!  IF THE VALUE OF RHO2 IS TOO SMALL THE SINH MUST BE OBTAINED
!  BY SUMMING THE INFINITE SERIES RATHER THAN BY ADDITION OF
!  TWO EXPONENTIALS.
      R=D-H
      IF(DABS(R)-0.1d0) 26,28,28
 26   RA=RHO2
      RHO22=RHO2*RHO2
      T=RHO2
      do i=2,50,2
       T=T*RHO22/DBLE(FLOAT(I*I+I))
       RA=RA+T
       IF(T.LT.1.D-30) exit
      enddo
      R=RA+RA
!
!  AS MANY SUCCESSIVE B-FUNCTIONS ARE GENERATED FROM B-0 BY THE
!  RECURRENCE FORMULAE AS ACCURACY WILL PERMIT.
 28   B(1)=R/RHO2
      DO 51 I=2,IX,IS
      IF(IR.EQ.0) GO TO 40
      IL=IS-1
      IF(1.GT.IL) GO TO 9050
      DO 31 J=1,IL
      K=I+J-1
      IF((-1)**K) 29,29,30
 29   B(K)=(R+DBLE(FLOAT(K-1))*B(K-1))/RHO2
      GO TO 31
 30   B(K)=-(D+H-DBLE(FLOAT(K-1))*B(K-1))/RHO2
 31   CONTINUE
 9050 CONTINUE
 40   IN=I+IS-1

!  AFTER THE RECURRENCE FORMULAE HAVE BEEN APPLIED AN APPROPRIATE
!  NUMBER OF TIMES THE NEXT B-FUNCTION IS OBTAINED BY SUMMATION
!  OF THE INFINITE SERIES.
      IF(IN-IX) 39,39,38
 39   IF((-1)**IN) 44,44,42
 42   TR=RHO2
      B(IN)=-2.0d0*TR/DBLE(FLOAT(IN+1))
      DO 43 J=1,500
      TR=TR*RHO2**2/DBLE(FLOAT((2*J)*(2*J+1)))
      IF(DABS(TR/B(IN))-1.0D-7 ) 51,51,43
 43   B(IN)=B(IN)-2.0d0*TR/DBLE(FLOAT(IN+1+2*J))
 44   TR=1.0d0
      B(IN)=2.0d0*TR/DBLE(FLOAT(IN))
      DO 46 J=1,500
      TR=TR*RHO2**2/DBLE(FLOAT((2*J)*(2*J-1)))
      IF(DABS(TR/B(IN))-1.0D-7 ) 51,51,46
 46   B(IN)=B(IN)+2.0d0*TR/DBLE(FLOAT(IN+2*J))
 51   CONTINUE

!  IF THE ARGUMENT IS ZERO A SEPARATE FORMULA MUST BE USED.
      GO TO 38
 35   continue
      do I=1,IX,2
       B(I)=2.0/DBLE(FLOAT(I))
       B(I+1)=0.0d0
      enddo
 38   RETURN
      END subroutine abfns
!  ********************************************************************
!  *                                                                  *
!  *     SUBROUTINE LOVLAP       CALLED FROM MOVLAP                   *
!  *                                                                  *
!  *                                                                  *
!  *       LOVLAP   SUBROUTINE TO CALCULATE THE OVERLAP               *
!  *                COMPONENT INDEPENDENT OF THE ANGLE BETWEEN        *
!  *                THE ATOMS                                         *
!  *                                                                  *
!  *       SUBROUTINES USED:                                          *
!  *                                                                  *
!  *             NONE                                                 *
!  *                                                                  *
!  *                                                                  *
!  *       ORIGIN LOST IN ANTIQUITY                                   *
!  *                                                                  *
!  ********************************************************************
      subroutine lovlap(STRAD,A,B, sk1, sk2, r, l1, l2, m1, n1, n2  )
      implicit none

real*8  , intent(in) :: a(30), b(30), sk1, sk2, r
integer , intent(in) :: l1, l2, m1, n1, n2
real*8  , intent(out):: strad

integer  :: m2, ir, ip, jend, kend, ieb, j, ju, iab, icb
integer  :: k, ku, iev, ibb, idb, i6, i5, i4, i3, i2, i1, i
real*8   :: rhoa, rhob, rhoap, rhoab, rhopo, term, terma
real*8   :: value0, value1, value2, value3, value4
real*8   :: con1, con12

      logical ,save :: ftime = .true.
      real*8  ,save :: fact(0:25)            
      real*8        :: BINCOE(8,8),BINCO(64)

!    COMMON/LOCLAP/SK1,SK2,R,L1,L2,M1,N1,N2,MAXN                      
     EQUIVALENCE (BINCOE(1,1),BINCO(1))                               

      DATA BINCO/8*1.0d0,0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,&
           7.0d0,2*0.0d0,1.0d0,3.0d0,6.0d0,10.0d0,15.0d0,21.0d0,&
           3*0.0d0,1.0d0,4.0d0,10.0d0,20.0d0,35.0d0,4*0.0d0,1.0d0,&
           5.0d0,15.0d0,35.0d0,5*0.0d0,1.0d0,6.0d0,21.0d0,6*0.0d0,&
           1.0d0,7.0d0,7*0.0d0,1.0d0/

!-----------------------------------------------------------------------
      if ( ftime ) then
       fact(0)=1.0d0
       fact(1)=1.0d0
       do i=2, 25                                                     
        fact(i)=fact(i-1)*dfloat(i-1)
       enddo
       ftime = .false.
      endif

      M2=M1                                                            
      STRAD=0.0d0
      RHOA=R*SK1                                                       
      RHOB=R*SK2                                                       
      RHOAP=RHOA**N1                                                   
      RHOAP=RHOAP*RHOAP                                                
      RHOAP=RHOAP*RHOA                                                 
      RHOAB=RHOB**N2                                                   
      RHOAB=RHOAB*RHOAB                                                
      RHOAB=RHOAB*RHOB                                                 
      RHOPO=RHOAP*RHOAB                                                
      TERMA=0.5d0**(L1+L2+1)*dsqrt(dfloat((L1+L1+1)*(L2+L2+1))*&
     fact(L1-M1+1)*fact(L2-M1+1)/(fact(N1+N1+1)*fact(N2+N2+1)*&
     fact(L1+M1+1)*fact(L2+M1+1))*RHOPO)                             

      JEND=1+((L1-M1)/2)                                               
      KEND=1+((L2-M2)/2)                                               
      IEB=M1+1                                                         
      do J=1,JEND                                                   
       JU=J-1                                                           
       IAB=N1-L1+JU+JU+1                                                
       ICB=L1-M1-JU-JU+1                                                
       CON1=fact(L1+L1-JU-JU+1)/(fact(L1-M1-JU-JU+1)*fact(JU+1)*fact(L1-JU+1))
       do K=1,KEND
        KU=K-1
        CON12=CON1*fact(L2+L2-KU-KU+1)/(fact(L2-M2-KU-KU+1)*fact(KU+1)*fact(L2-KU+1))
        IEV=JU+KU+L2
        IF(2*(IEV/2).NE.IEV) CON12=-CON12
        IBB=N2-L2+KU+KU+1
        IDB=L2-M2-KU-KU+1

        value0=0.0d0
        do I6=1,IEB
         do I5=1,IEB
          VALUE1=BINCOE(IEB,I6)*BINCOE(IEB,I5)
          IEV=I5+I6
          IF(2*(IEV/2).NE.IEV) VALUE1=-VALUE1
          do  I4=1,IDB
           VALUE1=-VALUE1
           VALUE2=BINCOE(IDB,I4)*VALUE1
           do I3=1,ICB
            VALUE3=BINCOE(ICB,I3)*VALUE2
            do I2=1,IBB
             VALUE3=-VALUE3
             VALUE4=BINCOE(IBB,I2)*VALUE3
             do I1=1,IAB
              TERM=VALUE4*BINCOE(IAB,I1)
              IR=I1+I2+IEB+IEB-I6-I6-I3+IDB-I4+ICB-1
              IP=IAB-I1+IBB-I2+IEB+IEB-I5-I5+ICB-I3+IDB-I4+1
              value0=value0+A(IP)*B(IR)*TERM
             enddo ! i1
            enddo ! i2
           enddo ! i3
          enddo ! i4
         enddo ! i5
        enddo ! i6
        STRAD=STRAD+value0*CON12
       enddo
      enddo
      STRAD=STRAD*TERMA                                                

      RETURN                                                           
      END subroutine lovlap

!-----------------------------------------------------------------------
!
! Given two l-shells li,lj, centered at atoms ri and rj, this subroutine
! generates all the < ljm | lim > by applying the two-center SK formulas
! 
! Input:
!
!   ovin: sigma, pi & delta components for current <li|lj> interaction
!         (they may correspond to the overlap or the Hamiltonian)
!   cl  : (xj-xi)/|rj-ri|
!   cm  : (yj-yi)/|rj-ri|
!   cn  : (zj-zi)/|rj-ri|
!   cl2,cm2,cn2: squares of the direction cosines cl,cm,cn
!   
!
! Output:
!   ovout: overlap/Hamiltonian elements are given in ovout
!
! Care with signs:
!
! ----------------------------------------------------------------------

subroutine twocenter(li,lj,ovin,ovout,cl,cl2,cm,cm2,cn,cn2)

      implicit none

      integer , intent(in) :: li,lj
      real*8,  intent(in)  :: ovin(3) , cl , cl2 , cm , cm2 , cn ,cn2
      real*8,  intent(out) :: ovout(25)

      integer  :: itype
      integer  :: kao2(0:2,0:2)
      real*8   :: sqr3
      data kao2 / 0 , 3 , 6 , 1 , 4 , 7 , 2 , 5 , 8 /

! ----------------------------------------------------------------------

      sqr3 = sqrt(3.)
      itype = kao2(li,lj)

      select case(itype)

! --  itype = 0 => (S(I):S(J))
      case(0)
        ovout(1)= ovin(1)

! --  itype = 1 => (S(I):P(J))
      case(1)
        ovout(1) = cl*ovin(1)                           ! s,px
        ovout(2) = cm*ovin(1)                           ! s,py
        ovout(3) = cn*ovin(1)                           ! s,pz

! --  itype = 2 => (S(I):D(J))
      case(2)
        ovout(1) = 0.5d0*sqr3*(cl2-cm2)   * ovin(1)     ! s,dx2-y2
        ovout(2) = (cn2-0.5d0*(cl2+cm2))  * ovin(1)     ! s,dz2
        ovout(3) = sqr3*cl*cm             * ovin(1)     ! s,dxy
        ovout(4) = sqr3*cn*cl             * ovin(1)     ! s,dzx
        ovout(5) = sqr3*cm*cn             * ovin(1)     ! s,dyz

! --  itype = 3 => (P(I):S(J)) 
      case(3)
        ovout(1) =-cl*ovin(1)                           ! px,s
        ovout(2) =-cm*ovin(1)                           ! py,s
        ovout(3) =-cn*ovin(1)                           ! pz,s

! --  itype = 4 => (P(I):P(J))
      case(4)
        ovout(1) =   cl2*ovin(1) + (1.0d0-cl2)*ovin(2)  ! px,px
        ovout(2) = cm*cl*ovin(1) - cm*cl* ovin(2)       ! py,px
        ovout(3) = cn*cl*ovin(1) - cn*cl* ovin(2)       ! pz,px
        ovout(4) = ovout(2)                             ! px,py
        ovout(5) =   cm2*ovin(1) + (1.0d0-cm2)*ovin(2)  ! py,py
        ovout(6) = cn*cm*ovin(1) - cn*cm*ovin(2)        ! pz,py
        ovout(7) = ovout(3)                             ! px,pz
        ovout(8) = ovout(6)                             ! py,pz
        ovout(9) =   cn2*ovin(1) + (1.0d0-cn2)*ovin(2)  ! pz,pz

! --  itype = 5 => (P(I):D(J))

      case(5)
        ovout( 1)= 0.5d0*sqr3*cl*(cl2-cm2)*ovin(1)  + cl*(1.0d0-cl2+cm2)*ovin(2)  ! px,dx2
        ovout( 2)= 0.5d0*sqr3*cm*(cl2-cm2)*ovin(1)  - cm*(1.0d0+cl2-cm2)*ovin(2)  ! py,dx2
        ovout( 3)= 0.5d0*sqr3*cn*(cl2-cm2)*ovin(1)  - cn*(cl2-cm2)*ovin(2)        ! pz,dx2
        ovout( 4)= cl*(cn2-0.5d0*(cl2+cm2))*ovin(1) - sqr3*cl*cn2*ovin(2)         ! px,dz2
        ovout( 5)= cm*(cn2-0.5d0*(cl2+cm2))*ovin(1) - sqr3*cm*cn2* ovin(2)        ! py,dz2
        ovout( 6)= cn*(cn2-0.5d0*(cl2+cm2))*ovin(1) + sqr3*cn*(cl2+cm2)* ovin(2)  ! pz,dz2
        ovout( 7)= sqr3*cl2*cm*ovin(1) + cm*(1.0d0-2.0d0*cl2)* ovin(2)            ! pxdxy
        ovout( 8)= sqr3*cm2*cl*ovin(1) + cl*(1.0d0-2.0d0*cm2)* ovin(2)            ! pydxy
        ovout( 9)= sqr3*cl*cm*cn*ovin(1) - 2.0d0*cl*cm*cn* ovin(2)                ! pzdxy
        ovout(10)= sqr3*cl2*cn*ovin(1) + cn*(1.0d0-2.0d0*cl2)* ovin(2)            ! px,dzx
        ovout(11)= ovout(9)                                                       ! py,dzx
        ovout(12)= sqr3*cn2*cl*ovin(1) + cl*(1.0d0-2.0d0*cn2)* ovin(2)            ! pz,dzx
        ovout(13)= ovout(9)                                                       ! px,dyz
        ovout(14)= sqr3*cm2*cn*ovin(1) + cn*(1.0d0-2.0d0*cm2)* ovin(2)            ! py,dyz
        ovout(15)= sqr3*cn2*cm*ovin(1) + cm*(1.0d0-2.0d0*cn2)* ovin(2)            ! py,dyz

! --  itype = 6 => (D(I):S(J))
      case(6)
        ovout(1) =  0.5d0*sqr3*(cl2-cm2)  * ovin(1)    ! s,dx2
        ovout(2) = (cn2-0.5d0*(cl2+cm2))  * ovin(1)    ! s,dz2
        ovout(3) = sqr3*cl*cm             * ovin(1)    ! s,dxy
        ovout(4) = sqr3*cn*cl             * ovin(1)    ! s,dzx
        ovout(5) = sqr3*cm*cn             * ovin(1)    ! s,dyz

! --  itype = 7 => (D(I):P(J))   ! Signals changed here too
      case(7)
        ovout( 1) =-0.5d0*sqr3*cl*(cl2-cm2)*ovin(1) - cl*(1.0d0-cl2+cm2)*ovin(2)   ! dx2,px
        ovout( 2) =-cl*(cn2-0.5d0*(cl2+cm2))*ovin(1) + sqr3*cl*cn2*ovin(2)         ! dz2,px
        ovout( 3) =-sqr3*cl2*cm*ovin(1) - cm*(1.0d0-2.0d0*cl2)*ovin(2)             ! dxy,px
        ovout( 4) =-sqr3*cl2*cn*ovin(1) - cn*(1.0d0-2.0d0*cl2)*ovin(2)             ! dzx,px
        ovout( 5) =-sqr3*cl*cm*cn*ovin(1) + 2.0d0*cl*cm*cn*ovin(2)                 ! dyz,px
        ovout( 6) =-0.5d0*sqr3*cm*(cl2-cm2)*ovin(1) + cm*(1.0d0+cl2-cm2)*ovin(2)   ! dx2,py
        ovout( 7) =-cm*(cn2-0.5d0*(cl2+cm2))*ovin(1) + sqr3*cm*cn2*ovin(2)         ! dz2,py
        ovout( 8) =-sqr3*cm2*cl*ovin(1) - cl*(1.0d0-2.0d0*cm2)*ovin(2)             ! dxy,py
        ovout( 9) = ovout(5)                                                       ! dzx,py
        ovout(10) =-sqr3*cm2*cn*ovin(1) - cn*(1.0d0-2.0d0*cm2)*ovin(2)             ! dyz,py
        ovout(11) =-0.5d0*sqr3*cn*(cl2-cm2)*ovin(1) + cn*(cl2-cm2)*ovin(2)         ! dx2,pz
        ovout(12) =-cn*(cn2-0.5d0*(cl2+cm2))*ovin(1) - sqr3*cn*(cl2+cm2)*ovin(2)   ! dz2,pz
        ovout(13) = ovout(5)                                                       ! dxy,pz
        ovout(14) =-sqr3*cn2*cl*ovin(1) - cl*(1.0d0-2.0d0*cn2)*ovin(2)             ! dzx,pz
        ovout(15) =-sqr3*cn2*cm*ovin(1) - cm*(1.0d0-2.0d0*cn2)*ovin(2)             ! dyz,pz

! --  itype = 8 =>  (D(I):D(J))
      case(8)
        ovout( 1)= 0.75d0*((cl2-cm2)**2)*ovin(1) + (cl2+cm2-(cl2-cm2)**2)*ovin(2) &
                  + (cn2+0.25*(cl2-cm2)**2)*ovin(3)                                      ! dx2,dx2
        ovout( 2)= 0.5d0*sqr3*(cl2-cm2)*(cn2-0.5d0*(cl2+cm2))*ovin(1) &
                  + sqr3*cn2*(cm2-cl2)*ovin(2) + (0.25*sqr3)*(1.+cn2)*(cl2-cm2)*ovin(3)  ! dz2,dx2
        ovout( 3)= 1.5d0*cl*cm*(cl2-cm2)*ovin(1) + 2.0d0*cl*cm*(cm2-cl2)*ovin(2) &
                  + 0.5d0*cl*cm*(cl2-cm2)*ovin(3)                                        ! dxy,dx2
        ovout( 4)= 1.5d0*cn*cl*(cl2-cm2)*ovin(1) + cn*cl*(1. - 2.*(cl2-cm2))*ovin(2) &
                  - cn*cl*(1.0d0-0.5d0*(cl2-cm2))*ovin(3)                                ! dzx,dx2
        ovout( 5)= 1.5d0*cm*cn*(cl2-cm2)*ovin(1) - cm*cn*(1. + 2.*(cl2-cm2))*ovin(2) &
                  + cm*cn*(1.0d0+0.5d0*(cl2-cm2))*ovin(3)                                ! dyz,dx2
        ovout( 6)= ovout(2)                                                              ! dx2,dz2 = dz2,dx2
        ovout( 7)= (cn2-0.5d0*(cl2+cm2))**2*ovin(1) + 3.0d0*cn2*(cl2+cm2)*ovin(2) &
                  +0.75d0*(cl2+cm2)**2*ovin(3)                                          ! dz2,dz2

        ovout( 8)= sqr3*cl*cm*(cn2-0.5*(cl2+cm2))*ovin(1) - sqr3*2.*cl*cm*cn2*ovin(2) &
                  +0.5d0*sqr3*cl*cm*(1.0d0+cn2)*ovin(3)                                 ! dxy,dz2
        ovout( 9)= sqr3*cl*cn*(cn2-0.5*(cl2+cm2))*ovin(1) &
                  +sqr3*cl*cn*(cl2+cm2-cn2)*ovin(2) - 0.5*sqr3*cl*cn*(cl2+cm2)*ovin(3)  ! dzx,dz2
        ovout(10) = sqr3*cm*cn*(cn2-0.5d0*(cl2+cm2))*ovin(1) &
                   +sqr3*cm*cn*(cl2+cm2-cn2)*ovin(2)- 0.5*sqr3*cm*cn*(cl2+cm2)*ovin(3)  ! dyz,dz2
        ovout(11)= ovout(3)                                                             ! dx2,dxy = dxy,dx2
        ovout(12)= ovout(8)                                                             ! dz2,dxy = dxy,dz2
        ovout(13)= 3.0d0*cl2*cm2*ovin(1) + (cl2+cm2-4.0d0*cl2*cm2)*ovin(2) &
                  + (cn2+cl2*cm2)*ovin(3)                                               ! dxy,dxy
        ovout(14)= 3.0d0*cn*cl2*cm*ovin(1) + cn*cm*(1.0d0-4.0d0*cl2)*ovin(2) &
                  + cn*cm*(cl2-1.0d0)*ovin(3)                                           ! dzx,dxy
        ovout(15)= 3.0d0*cm2*cn*cl*ovin(1) + cn*cl*(1.0d0-4.0d0*cm2)*ovin(2) &
                  + cn*cl*(cm2-1.0d0)*ovin(3)                                           ! dyz,dxy
        ovout(16)= ovout(4)                                                             ! dx2,dzx = dzx,dx2
        ovout(17)= ovout(9)                                                             ! dz2,dzx = dzx,dz2
        ovout(18)= ovout(14)                                                            ! dxy,dzx = dzx,dxy
        ovout(19)= 3.0d0*cn2*cl2*ovin(1) + (cn2+cl2-4.0d0*cn2*cl2)*ovin(2) &
                  + (cm2+cn2*cl2)*ovin(3)                                               ! dzx,dzy
        ovout(20)= 3.0d0*cm*cn2*cl*ovin(1) + cm*cl*(1.0d0-4.0d0*cn2)*ovin(2) &
                  + cm*cl*(cn2-1.0d0)*ovin(3)                                           ! dyz,dzx
        ovout(21)= ovout(5)                                                             ! dx2,dyz = dyz,dx2
        ovout(22)= ovout(10)                                                            ! dz2,dyz = dyz,dz2
        ovout(23)= ovout(15)                                                            ! dxy,dyz = dyz,dxy
        ovout(24)= ovout(20)                                                            ! dzx,dyz = dyz,dzx
        ovout(25)= 3.0d0*cm2*cn2*ovin(1) + (cm2+cn2-4.0d0*cm2*cn2)*ovin(2) &
                  + (cl2+cm2*cn2)*ovin(3)                                               !dyz,dyz

      case default
       write(*,'(a,i3,/a)') 'twocenter: parameter ITYPE had an ilegal value =', itype
       stop

      end select

      return
      end subroutine twocenter

!************************************************************************
  end module overlaps_jc
!  ******************************************************************** 
