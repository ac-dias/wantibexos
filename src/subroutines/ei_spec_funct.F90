!https://jblevins.org/mirror/amiller/toms683.f90
MODULE Exponential_Integral

!   ALGORITHM 683, COLLECTED ALGORITHMS FROM ACM.
!   THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!   VOL. 16, NO. 2, PP. 178-182.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

! COMMON /qad/ fnm, gln, x, y, iflag
REAL (dp), SAVE  :: fnm, gln, x, y
INTEGER, SAVE    :: iflag
PRIVATE
PUBLIC  :: cexint, cexqad

CONTAINS


SUBROUTINE cexqad(z, n, kode, tol, cw, kerr)

!  CEXQAD COMPUTES EXPONENTIAL INTEGRALS E(N,Z) OF A COMPLEX (dp) ARGUMENT
!  Z BY QUADRATURE ON  (T**(N-1)*EXP(-T)/(Z+T))/(N-1)!  FROM T=0 TO
!  T=INFINITY.  KODE=1 RETURNS CW=E(N,Z) WHILE KODE=2 RETURNS
!  CW=E(N,Z)*CEXP(Z).  TOL IS THE REQUESTED RELATIVE ERROR AND KERR.NE.0 IS
!  AN ERROR FLAG INDICATING A PREMATURE TRUNCATION OF AN INTEGRAL IN CEXQAD.
!  THE QUADRATURE IS DONE AS TWO REAL INTEGRALS SINCE Z IS COMPLEX (dp).

COMPLEX (dp), INTENT(IN)   :: z
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: kode
REAL (dp), INTENT(IN)      :: tol
COMPLEX (dp), INTENT(OUT)  :: cw
INTEGER, INTENT(OUT)       :: kerr

! COMMON /qad/ fnm, gln, x, y, iflag
! EXTERNAL fqcex
REAL (dp)  :: a, ans, b, fn, gm, s1, s2, xtol
INTEGER    :: i, ierr

fn = n
fnm = fn - 1.0
gm = 1.0_dp
IF (n /= 1) THEN
  DO  i = 2, n
    gm = gm * (i-1)
  END DO
END IF
gln = LOG(gm)
kerr = 0
x = REAL(z)
y = AIMAG(z)
!-----------------------------------------------------------------------
!     REAL PART OF E(N,Z)
!-----------------------------------------------------------------------
iflag = 1
s1 = 0.0_dp
b = 0.0_dp
DO  i = 1, 100
  a = b
  b = b + 5.0_dp
  xtol = tol
  CALL gaus8(fqcex, a, b, xtol, ans, ierr)
  s1 = s1 + ans
  IF (ABS(ans) < ABS(s1)*tol) GO TO 30
END DO
kerr = 1
!-----------------------------------------------------------------------
!     IMAGINARY PART OF E(N,Z)
!-----------------------------------------------------------------------
30 iflag = 2
s2 = 0.0_dp
IF (y /= 0.0_dp) THEN
  b = 0.0_dp
  DO  i = 1, 100
    a = b
    b = b + 5.0_dp
    xtol = tol
    CALL gaus8(fqcex, a, b, xtol, ans, ierr)
    s2 = s2 + ans
    IF (ABS(ans) < ABS(s2)*tol) GO TO 50
  END DO
  kerr = 2
END IF

50 cw = CMPLX(s1, -s2, KIND=dp)
IF (kode == 1) cw = cw * EXP(-z)
RETURN
END SUBROUTINE cexqad



FUNCTION fqcex(t) RESULT(fn_val)

REAL (dp), INTENT(IN)  :: t
REAL (dp)              :: fn_val

! COMMON /qad/ fnm, gln, x, y, iflag

REAL (dp)  :: a, b

a = -t + fnm * LOG(t) - gln
a = EXP(a)
IF (iflag /= 2) THEN
  b = (x+t) / ((x+t)**2 + y**2)
ELSE
  b = y / ((x+t)**2 + y**2)
END IF
fn_val = a * b
RETURN
END FUNCTION fqcex



FUNCTION g8(fun, x, h) RESULT(fn_val)
! Replaces the statement function g8 which was part of GAUS8

REAL (dp), INTENT(IN)  :: x, h
REAL (dp)              :: fn_val

INTERFACE
  FUNCTION fun(x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    REAL (dp)              :: fn_val
  END FUNCTION fun
END INTERFACE

REAL (dp), PARAMETER  :: x1 = 1.83434642495649805D-01, x2 = 5.25532409916328986D-01, &
                         x3 = 7.96666477413626740D-01, x4 = 9.60289856497536232D-01
REAL (dp), PARAMETER  :: w1 = 3.62683783378361983D-01, w2 = 3.13706645877887287D-01, &
                         w3 = 2.22381034453374471D-01, w4 = 1.01228536290376259D-01

fn_val = h * ((w1*(fun(x-x1*h) + fun(x+x1*h)) + w2*(fun(x-x2*h) +  &
         fun(x+x2*h))) + (w3*(fun(x-x3*h) + fun(x+x3*h)) + w4*(fun(x-x4*h) + &
         fun(x+x4*h))))

RETURN
END FUNCTION g8



SUBROUTINE gaus8(fun, a, b, ERR, ans, ierr)

!  WRITTEN BY R.E. JONES

!  ABSTRACT
!     GAUS8 INTEGRATES REAL FUNCTIONS OF ONE VARIABLE OVER FINITE INTERVALS
!     USING AN ADAPTIVE 8-POINT LEGENDRE-GAUSS ALGORITHM.
!     GAUS8 IS INTENDED PRIMARILY FOR HIGH ACCURACY INTEGRATION OR
!     INTEGRATION OF SMOOTH FUNCTIONS.

!     GAUS8 CALLS I1MACH, R1MACH, XERROR

!  DESCRIPTION OF ARGUMENTS

!     INPUT--
!     FUN - NAME OF EXTERNAL FUNCTION TO BE INTEGRATED.  THIS NAME MUST BE
!           IN AN EXTERNAL STATEMENT IN THE CALLING PROGRAM.
!           FUN MUST BE A FUNCTION OF ONE REAL ARGUMENT.  THE VALUE OF THE
!           ARGUMENT TO FUN IS THE VARIABLE OF INTEGRATION WHICH RANGES
!           FROM A TO B.
!     A   - LOWER LIMIT OF INTEGRAL
!     B   - UPPER LIMIT OF INTEGRAL (MAY BE LESS THAN A)
!     ERR - IS A REQUESTED PSEUDORELATIVE ERROR TOLERANCE.  NORMALLY PICK A
!           VALUE OF ABS(ERR) SO THAT STOL < ABS(ERR) <= 1.0E-3 WHERE STOL
!           IS THE DOUBLE PRECISION UNIT ROUNDOFF = EPSILON(0.0_dp).
!           ANS WILL NORMALLY HAVE NO MORE ERROR THAN ABS(ERR) TIMES THE
!           INTEGRAL OF THE ABSOLUTE VALUE OF FUN(X).
!           USUALLY, SMALLER VALUES FOR ERR YIELD MORE ACCURACY AND
!           REQUIRE MORE FUNCTION EVALUATIONS.

!           A NEGATIVE VALUE FOR ERR CAUSES AN ESTIMATE OF THE
!           ABSOLUTE ERROR IN ANS TO BE RETURNED IN ERR.  NOTE THAT
!           ERR MUST BE A VARIABLE (NOT A CONSTANT) IN THIS CASE.
!           NOTE ALSO THAT THE USER MUST RESET THE VALUE OF ERR
!           BEFORE MAKING ANY MORE CALLS THAT USE THE VARIABLE ERR.

!     OUTPUT--
!     ERR - WILL BE AN ESTIMATE OF THE ABSOLUTE ERROR IN ANS IF THE
!           INPUT VALUE OF ERR WAS NEGATIVE.  (ERR IS UNCHANGED IF
!           THE INPUT VALUE OF ERR WAS NONNEGATIVE.)  THE ESTIMATED
!           ERROR IS SOLELY FOR INFORMATION TO THE USER AND SHOULD
!           NOT BE USED AS A CORRECTION TO THE COMPUTED INTEGRAL.
!     ANS - COMPUTED VALUE OF INTEGRAL
!     IERR- A STATUS CODE
!         --NORMAL CODES
!            1 ANS MOST LIKELY MEETS REQUESTED ERROR TOLERANCE, OR A=B.
!           -1 A AND B ARE TOO NEARLY EQUAL TO ALLOW NORMAL INTEGRATION.
!              ANS IS SET TO ZERO.
!         --ABNORMAL CODE
!            2 ANS PROBABLY DOES NOT MEET REQUESTED ERROR TOLERANCE.
!***END PROLOGUE

REAL (dp), INTENT(IN)      :: a
REAL (dp), INTENT(IN)      :: b
REAL (dp), INTENT(IN OUT)  :: ERR
REAL (dp), INTENT(OUT)     :: ans
INTEGER, INTENT(OUT)       :: ierr

! EXTERNAL fun
INTERFACE
  FUNCTION fun(x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    REAL (dp)              :: fn_val
  END FUNCTION fun
END INTERFACE

INTEGER    :: k, l, lmn, lmx, lr(30), mxl, nbits, nib, nlmx
REAL (dp)  :: aa(30), ae, anib, area, c, ce, ee, ef, eps, est, gl, glr,  &
              gr(30), hh(30), tol, vl(30), vr
INTEGER, PARAMETER    :: nlmn = 1, kmx = 5000, kml = 6
INTEGER, SAVE         :: icall = 0
REAL (dp), PARAMETER  :: sq2 = 1.41421356_dp

!     INITIALIZE

IF (icall /= 0) THEN
  WRITE(*, *) 'GAUS8- GAUS8 CALLED RECURSIVELY; NOT ALLOWED HERE'
  RETURN
END IF

icall = 1
k = DIGITS(0.0_dp)
anib = LOG10(DBLE(RADIX(0.0_dp))) * k / 0.30102000_dp
nbits = INT(anib)
nlmx = (nbits*5) / 8
ans = 0.0_dp
ierr = 1
ce = 0.0_dp
IF (a /= b) THEN
  lmx = nlmx
  lmn = nlmn
  IF (b /= 0.0_dp) THEN
    IF (SIGN(1.0_dp,b)*a > 0.0_dp) THEN
      c = ABS(1.0_dp-a/b)
      IF (c <= 0.1_dp) THEN
        IF (c <= 0.0_dp) GO TO 100
        anib = 0.5_dp - LOG(c) / 0.69314718_dp
        nib = anib
        lmx = MIN(nlmx, nbits-nib-7)
        IF (lmx < 1) GO TO 90
        lmn = MIN(lmn,lmx)
      END IF
    END IF
  END IF
  tol = MAX(ABS(ERR), 2.0_dp**(5-nbits)) / 2.0_dp
  IF (ERR == 0.0_dp) tol = SQRT(EPSILON(0.0_dp))
  eps = tol
  hh(1) = (b-a) / 4.0_dp
  aa(1) = a
  lr(1) = 1
  l = 1
  est = g8(fun, aa(l)+2.0_dp*hh(l), 2.0_dp*hh(l))
  k = 8
  area = ABS(est)
  ef = 0.5_dp
  mxl = 0
  
!     COMPUTE REFINED ESTIMATES, ESTIMATE THE ERROR, ETC.
  
  10   gl = g8(fun, aa(l)+hh(l), hh(l))
  gr(l) = g8(fun, aa(l)+3.0_dp*hh(l), hh(l))
  k = k + 16
  area = area + (ABS(gl) + ABS(gr(l)) - ABS(est))
!     IF (L < LMN) GO TO 11
  glr = gl + gr(l)
  ee = ABS(est-glr) * ef
  ae = MAX(eps*area, tol*ABS(glr))
  IF (ee-ae > 0.0) THEN
    GO TO 40
  ELSE
    GO TO 30
  END IF
  20 mxl = 1
  30 ce = ce + (est-glr)
  IF (lr(l) > 0) THEN
    GO TO 70
  ELSE
    GO TO 50
  END IF
  
!     CONSIDER THE LEFT HALF OF THIS LEVEL
  
  40 IF (k > kmx) lmx = kml
  IF (l >= lmx) GO TO 20
  l = l + 1
  eps = eps * 0.5_dp
  ef = ef / sq2
  hh(l) = hh(l-1) * 0.5_dp
  lr(l) = -1
  aa(l) = aa(l-1)
  est = gl
  GO TO 10
  
!     PROCEED TO RIGHT HALF AT THIS LEVEL
  
  50 vl(l) = glr
  60 est = gr(l-1)
  lr(l) = 1
  aa(l) = aa(l) + 4.0_dp * hh(l)
  GO TO 10
  
!     RETURN ONE LEVEL
  
  70 vr = glr
  80 IF (l > 1) THEN
    l = l - 1
    eps = eps * 2.0_dp
    ef = ef * sq2
    IF (lr(l) <= 0) THEN
      vl(l) = vl(l+1) + vr
      GO TO 60
    END IF
    vr = vl(l+1) + vr
    GO TO 80
  END IF
  
!      EXIT
  
  ans = vr
  IF (mxl == 0 .OR. ABS(ce) <= 2.0_dp*tol*area) GO TO 100
  ierr = 2
  WRITE(*, *) 'GAUS8- ANS IS PROBABLY INSUFFICIENTLY ACCURATE.'
  GO TO 100

  90 ierr = -1
  WRITE(*, *) 'GAUS8- THE FOLLOWING TEMPORARY DIAGNOSTIC WILL APPEAR ONLY ONCE.'
  WRITE(*, *) 'A AND B ARE TOO NEARLY EQUAL TO ALLOW NORMAL INTEGRATION.'
  WRITE(*, *) 'ANS IS SET TO ZERO, AND IERR=-1.'
END IF
100 icall = 0
IF (ERR < 0.0_dp) ERR = ce
RETURN
END SUBROUTINE gaus8



SUBROUTINE cexint(z, n, kode, tol, m, cy, ierr)
!***BEGIN PROLOGUE  CEXINT
!***DATE WRITTEN   870515   (YYMMDD)
!***REVISION DATE  870515   (YYMMDD)
!***CATEGORY NO.  B5E
!***KEYWORDS  EXPONENTIAL INTEGRALS, SINE INTEGRAL, COSINE INTEGRAL
!***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
!***PURPOSE  TO COMPUTE EXPONENTIAL INTEGRALS OF A COMPLEX (dp) ARGUMENT
!***DESCRIPTION

!   ON KODE=1, CEXINT COMPUTES AN M MEMBER SEQUENCE OF COMPLEX (dp)
!   EXPONENTIAL INTEGRALS CY(J)=E(N+J-1,Z), J=1,...,M, FOR
!   POSITIVE ORDERS N,...,N+M-1 AND COMPLEX (dp) Z IN THE CUT PLANE
!   -PI < ARG(Z) <= PI (N=1 AND Z=CMPLX(0.0,0.0) CANNOT HOLD AT
!   THE SAME TIME).  ON KODE=2, CEXINT COMPUTES SCALED FUNCTIONS

!                    CY(J)=E(N+J-1,Z)*CEXP(Z),      J=1,...,M,

!   WHICH REMOVES THE EXPONENTIAL BEHAVIOR IN BOTH THE LEFT AND
!   RIGHT HALF PLANES.  DEFINITIONS AND NOTATION ARE FOUND IN THE
!   NBS HANDBOOK OF MATHEMATICAL FUNCTIONS (REF. 1).

!   INPUT
!     Z      - Z=CMPLX(X,Y), -PI < ARG(Z) <= PI
!     N      - INTEGER ORDER OF INITIAL E FUNCTION, N=1,2,...
!              (N=1 AND Z=CMPLX(0.0,0.0) IS AN ERROR)
!     KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
!              KODE= 1  RETURNS
!                       CY(J)=E(N+J-1,Z),          J=1,...,M
!                  = 2  RETURNS
!                       CY(J)=E(N+J-1,Z)*CEXP(Z),  J=1,...,M
!     TOL    - PRECISION (ACCURACY) DESIRED FOR THE SEQUENCE,
!              URND <= TOL <= 1.0E-3, WHERE URND IS LIMITED BY
!              URND = MAX(UNIT ROUNDOFF,1.0E-18) AND UNIT
!              ROUNDOFF = EPSILON(0.0_dp)
!     M      - NUMBER OF E FUNCTIONS IN THE SEQUENCE, M >= 1

!   OUTPUT
!     CY     - A COMPLEX (dp) VECTOR WHOSE FIRST M COMPONENTS CONTAIN
!              VALUES FOR THE SEQUENCE
!              CY(J)=E(N+J-1,Z)  OR
!              CY(J)=E(N+J-1,Z)*CEXP(Z), J=1,...,M
!              DEPENDING ON KODE.
!     IERR   - ERROR FLAG
!              IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
!              IERR=1, INPUT ERROR   - NO COMPUTATION
!              IERR=2, UNDERFLOW     - FIRST M COMPONENTS OF CY
!                      SET TO ZERO, CY(J)=CMPLX(0.0,0.0), J=1,M,
!                      REAL(Z) > 0.0 TOO LARGE ON KODE=1
!              IERR=3, OVERFLOW      - NO COMPUTATION,
!                      REAL(Z) < 0.0 TOO SMALL ON KODE=1
!              IERR=4, CABS(Z) OR N+M-1 LARGE - COMPUTATION DONE
!                      BUT LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION
!                      PRODUCE LESS THAN HALF OF MACHINE ACCURACY
!              IERR=5, CABS(Z) OR N+M-1 TOO LARGE - NO COMPUTATION BECAUSE
!                      OF COMPLETE LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION
!              IERR=6, ERROR         - NO COMPUTATION,
!                      ALGORITHM TERMINATION CONDITION NOT MET.
!                      SEE LONG DESCRIPTION ABOUT PARAMETER ICMAX.
!              IERR=7, ERROR         - NO COMPUTATION,
!                      DISCRIMINATION ERROR.  THIS CONDITION SHOULD NEVER OCCUR.

!***LONG DESCRIPTION

!   CEXINT USES A COMBINATION OF POWER SERIES AND BACKWARD RECURRENCE
!   DESCRIBED IN REF. 2 FOR THE COMPLEX (dp) Z PLANE EXCEPT FOR A STRIP
!   2*YB WIDE ABOUT THE NEGATIVE REAL AXIS, WHERE ANALYTIC CONTINUATION IS
!   CARRIED OUT BY LIMITED USE OF TAYLOR SERIES.
!   THE SWITCH FROM BACKWARD RECURRENCE TO TAYLOR SERIES IS NECESSARY BECAUSE
!   BACKWARD RECURRENCE IS SLOWLY CONVERGENT NEAR THE NEGATIVE REAL AXIS.
!   THE BOUNDARIES Y=-YB AND Y=YB WERE DETERMINED SO THAT BACKWARD RECURRENCE
!   WOULD CONVERGE EASILY WITH N AS LARGE AS 100 AND TOL AS SMALL AS 1.0D-18.
!   SUBROUTINE CEXENZ DOES THE BACKWARD RECURRENCE AND SUBROUTINE CACEXI DOES
!   THE ANALYTIC CONTINUATION.  TO START THE CONTINUATION, CACEXI CALLS CEXENZ
!   WITH ZB=CMPLX(X,YB).
!   IF CEXENZ RETURNS IERR=6, THEN YB IS INCREASED BY 0.5 UNTIL CEXENZ RETURNS
!   IERR=0 OR 10 TRIES, WHICHEVER COMES FIRST.
!   WHEN IERR=0, THEN THE ANALYTIC CONTINUATION PROCEEDS VERTICALLY DOWN FROM
!   ZB=CMPLX(X,YB) TO Z=CMPLX(X,Y), 0 <= Y < YB.
!   CONJUGATION IS USED FOR Y < 0.  YB INCREASES AS TOL DECREASES TO KEEP
!   CONVERGENCE RATES UP AND RECURRENCE DOWN.

!   PARAMETER ICDIM=250 ALLOCATES STORAGE FOR THE COEFFICIENTS OF THE BACKWARD
!   RECURRENCE ALGORITHM.  IF THE ALGORITHM TERMINATION CONDITION IS NOT MET
!   IN ICDIM STEPS, THEN RECURRENCE PROCEEDS WITH NO ADDITIONAL STORAGE UNTIL
!   THE TERMINATION CONDITION IS MET OR THE LIMIT ICMAX=2000 IS EXCEEDED.
!   THE PURPOSE OF STORAGE IS TO MAKE THE ALGORITHM MORE EFFICIENT.
!   THE TERMINATION CONDITION IS MET IN LESS THAN 250 STEPS OVER MOST OF
!   THE COMPLEX (dp) PLANE EXCLUDING THE STRIP ABS(Y) < YB, X < 0.
!   EXCEPTIONS TO THIS RULE ARE GENERATED NEAR STRIP BOUNDARIES WHEN N+M-1
!   AND ABS(Z) ARE LARGE AND NEARLY EQUAL.  IN THESE CASES, THE CONVERGENCE
!   IS VERY SLOW AND ADDITIONAL RECURRENCE (UP TO ICMAX) MUST BE USED.
!   ON THE OTHERHAND, THESE REGIONS OF SLOW CONVERGENCE ARE KEPT SMALL BY
!   ADJUSTING YB AS A FUNCTION OF TOL.  THESE REGIONS COULD BE ELIMINATED
!   ENTIRELY BY MAKING YB SUFFICIENTLY LARGE, BUT THE EXPENSE AND INSTABILITY
!   OF CONTINUATION BY TAYLOR SERIES NOT ONLY GOES UP, BUT THE COMPUTATIONAL
!   EXPENSE BECOMES EXCESSIVELY LARGE IN OTHER PARTS OF THE LEFT HALF PLANE
!   (Y < YB) WHERE THE BACKWARD RECURRENCE ALGORITHM WOULD CONVERGE RAPIDLY.

!   DERIVATIVES FOR SUCCESSIVE POWER SERIES ARE NOT COMPUTED BY EVALUATING
!   DERIVATIVES OF A PREVIOUS POWER SERIES.  BECAUSE OF THE RELATION

!     (1)           DE(N,Z)/DZ = - E(N-1,Z),

!   SUCCESSIVE DERIVATIVES AT Z ARE GIVEN BY LOWER ORDER FUNCTIONS AND CAN BE
!   COMPUTED IN A STABLE FASHION BY BACKWARD RECURRENCE USING (2) PROVIDED
!   THAT THE BEGINNING ORDER NUB IS SMALLER THAN THE ARGUMENT.
!   TO ACHIEVE THIS FOR ALL INTERMEDIATE VALUES ZZ BETWEEN ZB AND Z, WE TAKE
!   NUB=MINO(N+M-1,INT(CABS(Z)+0.5)).
!   TO START, E(NUB,ZB) IS EVALUATED BY THE BACKWARD RECURRENCE ALGORITHM OF
!   REF. 3.  TO CONTINUE THE FUNCTION FROM ZB TO Z VIA INTERMEDIATE VALUES ZZ,
!   DERIVATIVES OF E(NUB,ZB) ARE COMPUTED BY BACKWARD RECURRENCE ON (2).
!   THIS ALLOWS A STEP (NO LARGER THAN 0.5) TO ZZ FOR E(NUB,ZZ) USING THE
!   TAYLOR SERIES ABOUT ZB.  NOW, WE APPLY (2) AGAIN STARTING AT E(NUB,ZZ) FOR
!   THE DERIVATIVES AT ZZ AND TAKE ANOTHER STEP, ETC.  NOTICE THAT THE
!   STABILITY CONDITION FOR BACKWARD RECURRENCE, NUB <= ABS(Z) <= ABS(ZZ)
!   <= ABS(ZB), IS SATISFIED FOR ALL INTERMEDIATE VALUES ZZ.  THE FINAL
!   SEQUENCE FOR ORDERS N,...,N+M-1 IS GENERATED FROM (2) BY BACKWARD
!   RECURRENCE, FORWARD RECURRENCE OR BOTH ONCE E(NUB,Z) HAS BEEN COMPUTED.

!   RECURRENCE WITH THE RELATION

!      (2)     N*E(N+1,Z) + Z*E(N,Z) = CEXP(-Z)

!   IN A DIRECTION AWAY FROM THE INTEGER CLOSEST TO ABS(Z) IS STABLE.
!   FOR NEGATIVE ORDERS, THE RECURRENCE

!         E( 0,Z) = CEXP(-Z)/Z
!         E(-N,Z) = ( CEXP(-Z)+N*E(-N+1,Z) )/Z   ,N=1,2,...

!   IS NUMERICALLY STABLE FOR ALL Z.

!   THE (CAPITAL) SINE AND COSINE INTEGRALS CAN BE COMPUTED FROM

!           SI(Z) =  (E(1,I*Z)-E(1,-I*Z))/(2*I) + PI/2
!           CI(Z) = -(E(1,I*Z)+E(1,-I*Z))/2

!   IN -PI/2 < ARG(Z) <= PI/2, (I**2=-1), WHILE THE PRINCIPAL
!   VALUED EXPONENTIAL INTEGRAL EI(X) CAN BE COMPUTED FROM

!       EI( X) = -(E(1,-X+I*0)+E(1,-X-I*0))/2 = -REAL(E(1,-X))
!       EI(-X) = -REAL(E(1,X))

!   FOR X > 0.0 TO AN ACCURACY TOL.  IF Z = X > 0 THEN THE REAL SINE AND
!   COSINE INTEGRALS ARE GIVEN BY

!           SI(X) = AIMAG(E(1,I*X)) + PI/2
!           CI(X) = -REAL(E(1,I*X)) .

!   THE ANALYTIC CONTINUATION TO OTHER SHEETS CAN BE DONE BY THE RELATIONS

!   E(N,Z*CEXP(2*PI*M*I)) = E(N,Z) - 2*PI*M*I*(-Z)**(N-1)/(N-1)!

!   WHERE M=+1 OR M=-1 AND I**2=-1.

!***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND
!           I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF COMMERCE, 1955.

!         COMPUTATION OF EXPONENTIAL INTEGRALS OF COMPLEX (dp) ARGUMENT
!           BY D. E. AMOS, ACM TRANS. MATH. SOFTWARE

!         COMPUTATION OF EXPONENTIAL INTEGRALS BY D. E. AMOS, ACM TRANS.
!           MATH. SOFTWARE, VOL 6, NO. 3 SEPTEMBER 1980, PP. 365-377;
!           ALGORITHM 556, EXPONENTIAL INTEGRALS, PP. 420-428.

!         REMARK ON ALGORITHM 556
!           BY D. E. AMOS, ACM TRANS. MATH. SOFTWARE, VOL 9, NO. 4
!           DECEMBER 1983, P. 525.

!         UNIFORM ASYMPTOTIC EXPANSIONS FOR EXPONENTIAL INTEGRALS E(N,X)
!           AND BICKLEY FUNCTIONS KI(N,X) BY D. E. AMOS, ACM TRANS. MATH.
!           SOFTWARE, VOL 9, NO. 4, DECEMBER. 1983, PP. 467-479;
!           ALGORITHM 609, A PORTABLE FORTRAN SUBROUTINE FOR BICKLEY
!           FUNCTIONS KI(N,X), PP. 480-493.

!***ROUTINES CALLED  CACEXI,CEXENZ,I1MACH,R1MACH
!***END PROLOGUE  CEXINT

COMPLEX (dp), INTENT(IN)   :: z
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: kode
REAL (dp), INTENT(IN)      :: tol
INTEGER, INTENT(IN)        :: m
COMPLEX (dp), INTENT(OUT)  :: cy(m)
INTEGER, INTENT(OUT)       :: ierr

INTEGER      :: i, k, k1, k2
REAL (dp)    :: aa, alim, az, bb, ca(250), d, elim, fn, rbry, r1m5, urnd, x,  &
                y, yb, htol
COMPLEX (dp) :: cb(250)
!-----------------------------------------------------------------------
!     DIMENSION CA(ICDIM),CB(ICDIM)
!-----------------------------------------------------------------------
INTEGER, PARAMETER  :: icdim = 250

!***FIRST EXECUTABLE STATEMENT  CEXINT
ierr = 0
x = REAL(z, KIND=dp)
y = AIMAG(z)
IF (x == 0.0_dp .AND. y == 0.0_dp .AND. n == 1) ierr = 1
IF (n < 1) ierr = 1
IF (kode < 1 .OR. kode > 2) ierr = 1
IF (m < 1) ierr = 1
urnd = MAX(EPSILON(0.0_dp), 1.0E-18)
IF (tol < urnd .OR. tol > 1.0E-3) ierr = 1
IF (ierr /= 0) RETURN
IF (x /= 0.0_dp .OR. y /= 0.0_dp .OR. n <= 1) THEN
!-----------------------------------------------------------------------
!     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
!     URND IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
!     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
!     EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/URND    AND
!     EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*URND       ARE INTERVALS NEAR
!     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
!-----------------------------------------------------------------------
  k1 = MINEXPONENT(0.0_dp)
  k2 = MAXEXPONENT(0.0_dp)
  r1m5 = LOG10(DBLE(RADIX(0.0_dp)))
  k = MIN(ABS(k1),ABS(k2))
  elim = 2.303_dp * (k*r1m5 - 3.0_dp)
  k1 = DIGITS(0.0_dp) - 1
  aa = r1m5 * k1
  aa = aa * 2.303
  alim = elim + MAX(-aa, -41.45_dp)
  rbry = 2.0
  IF (urnd > 1.0E-8) rbry = 1.0_dp
!-----------------------------------------------------------------------
!     TEST VARIABLES FOR RANGE. ABS(Z) CANNOT BE LARGER THAN THE ARGUMENT
!     OF THE INT( ) FUNCTION.
!-----------------------------------------------------------------------
  az = ABS(z)
  fn = n+m-1
  aa = 0.5_dp / urnd
  bb = HUGE(0.0_dp) * 0.5_dp
  aa = MIN(aa,bb)
  IF (az > aa) GO TO 20
  IF (fn > aa) GO TO 20
  aa = SQRT(aa)
  IF (az > aa) ierr = 4
  IF (fn > aa) ierr = 4
  IF (x >= 0.0_dp) THEN
!-----------------------------------------------------------------------
!     BACKWARD RECURRENCE FOR THE RIGHT HALF PLANE, X >= 0.0E0
!-----------------------------------------------------------------------
    CALL cexenz(z, n, kode, m, cy, ierr, rbry, tol, elim, alim, icdim, ca, cb)
    RETURN
  END IF
  IF (az <= rbry) THEN
!-----------------------------------------------------------------------
!     POWER SERIES FOR ABS(Z) <= RBRY AND X < 0.0E0
!-----------------------------------------------------------------------
    CALL cexenz(z, n, kode, m, cy, ierr, rbry, tol, elim, alim, icdim, ca, cb)
    RETURN
  END IF
  d = -0.4342945_dp * LOG(tol)
  yb = 10.5_dp - 0.538460_dp * (18.0_dp-d)
  IF (ABS(y) >= yb) THEN
!-----------------------------------------------------------------------
!     BACKWARD RECURRENCE EXTERIOR TO THE STRIP ABS(Y) < YB, X < 0.0
!-----------------------------------------------------------------------
    htol = 0.125_dp * tol
    CALL cexenz(z, n, kode, m, cy, ierr, rbry, htol, elim, alim, icdim, ca, cb )
    RETURN
  END IF
!-----------------------------------------------------------------------
!     TAYLOR SERIES IN CACEXI FOR ANALYTIC CONTINUATION
!-----------------------------------------------------------------------
  CALL cacexi(z, n, kode, m, cy, ierr, yb, rbry, tol, elim, alim, icdim, ca, cb )
  RETURN
END IF
DO  i = 1, m
  cy(i) = 1.0_dp / (n+i-2)
END DO
RETURN

20 ierr = 5
RETURN
END SUBROUTINE cexint



SUBROUTINE cexenz(z, n, kode, m, cy, ierr, rbry, tol, elim, alim, icdim, ca, cb)
!***BEGIN PROLOGUE  CEXENZ
!***REFER TO  CEXINT

!   CEXENZ COMPUTES THE EXPONENTIAL INTEGRAL BY MEANS OF POWER SERIES AND A
!   BACKWARD RECURRENCE ALGORITHM FOR THE CONFLUENT HYPERGEOMETRIC
!   REPRESENTATION.

!***ROUTINES CALLED PSIXN
!***END PROLOGUE  CEXENZ

COMPLEX (dp), INTENT(IN)   :: z
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN)        :: m
COMPLEX (dp), INTENT(OUT)  :: cy(m)
INTEGER, INTENT(OUT)       :: ierr
REAL (dp), INTENT(IN)      :: rbry
REAL (dp), INTENT(IN)      :: tol
REAL (dp), INTENT(IN)      :: elim
REAL (dp), INTENT(IN)      :: alim
INTEGER, INTENT(IN)        :: icdim
REAL (dp), INTENT(OUT)     :: ca(icdim)
COMPLEX (dp), INTENT(OUT)  :: cb(icdim)

INTEGER       :: i, ic, icase, ict, ik, ind, iz, jset, k, kflag, kn, ks, ml, &
                 mu, nd, nm
REAL (dp)     :: aa, aam, aams, aem, ah, ak, ap1, at, az, bk, bt, dk, ERR, &
                 fc, fnm, rtola, tola, x, xtol, y, ck
COMPLEX (dp)  :: cp1, cp2, cpt, cat, cbt, cy1, cy2, cyy(2), cnorm, cs, cak, &
                 emz, caa, tz, fz, ct, scle, rscle

REAL (dp), PARAMETER     :: euler = -5.77215664901532861D-01
COMPLEX (dp), PARAMETER  :: czero = (0.0_dp, 0.0_dp), cone = (1.0_dp, 0.0_dp)
INTEGER, PARAMETER       :: icmax = 2000

ierr = 0
scle = cone
rscle = cone
x = REAL(z, KIND=dp)
y = AIMAG(z)
az = ABS(z)
IF (az <= rbry) THEN
!-----------------------------------------------------------------------
!     SERIES FOR E(N,Z) FOR ABS(Z) <= RBRY
!-----------------------------------------------------------------------
  iz = az + 0.5_dp
!-----------------------------------------------------------------------
!     ICASE=1 MEANS INTEGER CLOSEST TO ABS(Z) IS 2 AND N=1
!     ICASE=2 MEANS INTEGER CLOSEST TO ABS(Z) IS 0,1, OR 2 AND N >= 2
!-----------------------------------------------------------------------
  icase = 2
  IF (iz > n) icase = 1
  nm = n - icase + 1
  nd = nm + 1
  ind = 3 - icase
  mu = m - ind
  ml = 1
  ks = nd
  fnm = nm
  cs = czero
  xtol = 0.3333_dp * tol
  aam = 1.0_dp
  IF (nd /= 1) THEN
    aam = 1.0_dp / fnm
    cs = aam
  END IF
  caa = cone
  aa = 1.0_dp
  ak = 1.0_dp
!-----------------------------------------------------------------------
!     LIMIT INDEX I TO IK TO PREVENT UNDERFLOW ON SMALL VALUES OF Z
!-----------------------------------------------------------------------
  ik = 35
  IF (az < xtol*aam) ik = 1
  DO  i = 1, ik
    at = 1.0_dp / ak
    caa = -caa * z * at
    aa = aa * az * at
    IF (i /= nm) THEN
      cs = cs - caa / (ak-fnm)
    ELSE
      cs = cs + caa * (-LOG(z) + psixn(nd))
    END IF
    IF (aa <= xtol*ABS(cs)) EXIT
    ak = ak + 1.0_dp
  END DO

  IF (nd == 1) cs = cs + (-LOG(z) + euler)
  IF (kode == 2) cs = cs * EXP(z)
  cy(1) = cs
  ct = cs
  emz = cone
  IF (m /= 1) THEN
    cy(ind) = cs
    ak = ks
    IF (kode == 1) emz = EXP(-z)
    SELECT CASE ( icase )
      CASE (    1)
        GO TO 140
      CASE (    2)
        GO TO 160
    END SELECT
  END IF
  IF (icase == 2) RETURN
  IF (kode == 1) emz = EXP(-z)
  cy(1) = (emz-cs) / z
  RETURN
END IF
!-----------------------------------------------------------------------
!     BACKWARD RECURSIVE MILLER ALGORITHM FOR
!              E(N,Z)=EXP(-Z)*(Z**(N-1))*U(N,N,Z)
!     WITH RECURSION AWAY FROM N=INTEGER CLOSEST TO ABS(Z)
!     U(A,B,Z) IS THE SECOND CONFLUENT HYPERGEOMETRIC FUNCTION
!-----------------------------------------------------------------------
emz = cone
IF (kode /= 2) THEN
!-----------------------------------------------------------------------
!     SCALE NEAR EXPONENT EXTREMES ON KODE=1
!-----------------------------------------------------------------------
  IF (x >= 0.0_dp) THEN
    at = x + (n+m-1)
    ct = CMPLX(at, y, KIND=dp)
    aa = ABS(ct)
    at = x + LOG(aa)
    IF (at > elim) GO TO 30
    kflag = 1
  ELSE
    at = x
    IF (at < (-elim)) GO TO 180
    kflag = 2
  END IF
  IF (ABS(at) >= alim) THEN
    tola = EXP(alim-elim)
    rtola = 1.0_dp / tola
    IF (kflag /= 2) THEN
      scle  = rtola
      rscle = tola
    ELSE
      scle  = tola
      rscle = rtola
    END IF
    emz = scle
  END IF
  emz = emz * EXP(-z)
END IF
iz = az + 0.5_dp
kn = n + m - 1
IF (kn <= iz) GO TO 50
IF (n < iz .AND. iz < kn) GO TO 80
IF (n >= iz) GO TO 70
ierr = 7
RETURN

30 ierr = 2
cy(1:m) = czero
RETURN

50 icase = 1
ks = kn
ml = m - 1
mu = -1
ind = m
IF (kn > 1) GO TO 90

60 ks = 2
icase = 3
GO TO 90

70 icase = 2
ind = 1
ks = n
mu = m - 1
IF (n > 1) GO TO 90
IF (kn == 1) GO TO 60
iz = 2

80 icase = 1
ks = iz
ml = iz - n
ind = ml + 1
mu = kn - iz

90 ik = ks / 2
ah = ik
jset = 1 + ks - 2 * ik
!-----------------------------------------------------------------------
!     START COMPUTATION FOR
!              CYY(1) = C*U( A , A ,Z)    JSET=1
!              CYY(1) = C*U(A+1,A+1,Z)    JSET=2
!     FOR AN EVEN INTEGER A.
!-----------------------------------------------------------------------
ic = 0
aa = ah + ah
caa = aa
aam = aa - 1.0_dp
aams = aam * aam
tz = z + z
fz = tz + tz
ak = ah
xtol = tol
ct = aams + fz * ah
cak = z + caa
aem = (ak+1.0_dp) / xtol
aem = aem / ABS(cak)
bk = aa
ck = ah * ah
!-----------------------------------------------------------------------
!     FORWARD RECURSION FOR P(IC),P(IC+1) AND INDEX IC FOR BACKWARD
!     RECURSION
!-----------------------------------------------------------------------
cp1 = czero
cp2 = cone

100 ic = ic + 1
IF (ic <= icdim) THEN
  ak = ak + 1.0_dp
  ck = ck + 1.0_dp
  at = bk / (bk+ak+ck)
  bk = bk + ak + ak
  cat = at
  ca(ic) = at
  bt = 1.0_dp / (ak + 1.0_dp)
  cbt = (ak + ak + z) * bt
  cb(ic) = cbt
  cpt = cp2
  cp2 = cbt * cp2 - cat * cp1
  cp1 = cpt
  ct = ct + fz
  aem = aem * at
  bt = ABS(ct)
  ERR = aem / SQRT(bt)
  ap1 = ABS(cp1)
  IF (ERR*(ak+1.0_dp)/ap1 > ap1) GO TO 100
ELSE
!-----------------------------------------------------------------------
!     CONTINUE FORWARD RECURRENCE UNINDEXED WHEN IC EXCEEDS ICDIM
!-----------------------------------------------------------------------
  ic = ic - 1

  110 ic = ic + 1
  IF (ic > icmax) GO TO 190
  ak = ak + 1.0_dp
  ck = ck + 1.0_dp
  at = bk / (bk+ak+ck)
  bk = bk + ak + ak
  cat = at
  bt = 1.0_dp / (ak+1.0_dp)
  cbt = (ak+ak + z) * bt
  cpt = cp2
  cp2 = cbt * cp2 - cat * cp1
  cp1 = cpt
  ct = ct + fz
  aem = aem * at
  bt = ABS(ct)
  ERR = aem / SQRT(bt)
  ap1 = ABS(cp1)
  IF (ERR*(ak+1.0_dp)/ap1 > ap1) GO TO 110
END IF
fc = ic
at = ((fc+1.0_dp)/(ak+1.0_dp)) * ((ak+ah)/(ak+1.0_dp))
cat = at * SQRT(ct/(ct+fz))
cy2 = cat * (cp1/cp2)
cy1 = cone
!-----------------------------------------------------------------------
!     BACKWARD RECURRENCE FOR
!             CY1=             C*U( A ,A,Z)
!             CY2= C*(A/(1+A/2))*U(A+1,A,Z)
!-----------------------------------------------------------------------
IF (ic > icdim) THEN
!-----------------------------------------------------------------------
!     BACKWARD RECUR UNINDEXED WHEN IC EXCEEDS ICDIM
!-----------------------------------------------------------------------
  bt = aa + x

  120 ak = ah + fc
  bk = ak + 1.0_dp
  ck = aam + fc
  dk = bt + fc + fc
  at = (ak/fc) * (bk/ck)
  cbt = CMPLX(dk/bk, y/bk, KIND=dp)
  cpt = cy1
  cy1 = (cbt*cy1 - cy2) * at
  cy2 = cpt
  fc = fc - 1.0_dp
  ic = ic - 1
  IF (ic > icdim) GO TO 120
END IF
ict = ic
DO  k = 1, ict
  at = 1.0_dp / ca(ic)
  cpt = cy1
  cy1 = (cb(ic)*cy1 - cy2) * at
  cy2 = cpt
  ic = ic - 1
END DO
!-----------------------------------------------------------------------
!     THE CONTIGUOUS RELATION
!              Z*U(B,C+1,Z)=(C-B)*U(B,C,Z)+U(B-1,C,Z)
!     WITH  B=A+1 , C=A IS USED FOR
!             CYY(2) = C * U(A+1,A+1,Z)
!     Z IS INCORPORATED INTO THE NORMALIZING RELATION FOR CNORM.
!-----------------------------------------------------------------------
cpt = cy2 / cy1
cnorm = cone - cpt * (ah + 1.0_dp) / aa
cyy(1) = cone / (cnorm*caa+z)
cyy(2) = cnorm * cyy(1)
IF (icase /= 3) THEN
  ct = emz * cyy(jset)
  cy(ind) = ct * rscle
  cs = ct
  IF (m == 1) RETURN
  ak = ks
  SELECT CASE ( icase )
    CASE (    1)
      GO TO 140
    CASE (    2)
      GO TO 160
  END SELECT
END IF
!-----------------------------------------------------------------------
!     RECURSION SECTION  N*E(N+1,Z) + Z*E(N,Z)=EMZ
!-----------------------------------------------------------------------
ct = emz * (cone - cyy(1)) / z
cy(1) = ct * rscle
RETURN

140 caa = ak
tz = cone / z
k = ind - 1
DO  i = 1, ml
  caa = caa - cone
  ct = (emz-caa*ct) * tz
  cy(k) = ct * rscle
  k = k - 1
END DO
IF (mu <= 0) RETURN
ak = ks

160 k = ind
DO  i = 1, mu
  cs = (emz - z*cs) / ak
  cy(k+1) = cs * rscle
  ak = ak + 1.0_dp
  k = k + 1
END DO
RETURN

180 ierr = 3
RETURN

190 ierr = 6
RETURN
END SUBROUTINE cexenz



SUBROUTINE cacexi(z, nu, kode, n, y, ierr, yb, rbry, tol, elim, alim, icdim, &
                  ca, cb)
!***BEGIN PROLOGUE  CACEXI
!***REFER TO  CEXINT

!  CACEXI COMPUTES THE ANALYTIC CONTINUATION OF THE EXPONENTIAL INTEGRAL
!  FOR X < 0 AND -YB < Y < YB BY TAYLOR SERIES IN INCREMENTS OF HALF A UNIT.
!  THE CONTINUATION PROCEEDS VERTICALLY DOWN FROM ZB=CMPLX(X,YB) TO
!  Z=CMPLX(X,Y) FOR 0.0 <= Y < YB.  CONJUGATION IS USED FOR Y < 0.0E0.

!***ROUTINES CALLED  CEXENZ
!***END PROLOGUE  CACEXI

COMPLEX (dp), INTENT(IN)      :: z
INTEGER, INTENT(IN)           :: nu
INTEGER, INTENT(IN)           :: kode
INTEGER, INTENT(IN)           :: n
COMPLEX (dp), INTENT(OUT)     :: y(n)
INTEGER, INTENT(OUT)          :: ierr
REAL (dp), INTENT(IN OUT)     :: yb
REAL (dp), INTENT(IN OUT)     :: rbry
REAL (dp), INTENT(IN)         :: tol
REAL (dp), INTENT(IN)         :: elim
REAL (dp), INTENT(IN)         :: alim
INTEGER, INTENT(IN)           :: icdim
REAL (dp), INTENT(IN OUT)     :: ca(icdim)
COMPLEX (dp), INTENT(IN OUT)  :: cb(icdim)

INTEGER       :: i, iaz, il, is, iy, k, kb, kl, kmax, ks, kyb, nb, nflg, nub
REAL (dp)     :: az, del, fk, rtola, tola, yt,  &
                 zi, zid, zr, xtol, rzi, atrm, fj, rw, rq(64), asum, htol
COMPLEX (dp)  :: yy(1), cex, sum, trm, zz, zp(64), zt, scle,  &
                 rscle, trms, zw, cezt
COMPLEX (dp), PARAMETER  :: cone = (1.0_dp, 0.0_dp)

scle = cone
rscle = cone
zr = REAL(z, KIND=dp)
zi = AIMAG(z)
zid = zi
kyb = 0
IF (zi < 0.0_dp) zid = -zid
az = ABS(z)
iaz = az + 0.5_dp
!-----------------------------------------------------------------------
!     SET ORDER NUB=MIN(N+M-1,INT(ABS(Z)+0.5)), GENERATE REMAINDER
!     OF THE SEQUENCE AFTER THE COMPUTATION OF E(NUB,Z)
!-----------------------------------------------------------------------
IF (nu < iaz) THEN
  IF (nu+n-1 <= iaz) GO TO 10
  nub = iaz
  nb = 0
  nflg = 3
  kb = nub - nu + 1
  ks = kb + 1
  kl = n
  is = 1
  il = kb - 1
  GO TO 30
END IF
nub = MAX(iaz, 1)
nb = nu - nub
nflg = 1
kb = 1
ks = 2
kl = n
GO TO 30

10 nub = nu + n - 1
nb = 0
nflg = 2
is = 2
il = n
kb = n
GO TO 30
!-----------------------------------------------------------------------
!     SET PARAMETERS FOR ANALYTIC CONTINUATION FROM Y=YB INTO THE REGION
!     0 <= ZID < YB.
!-----------------------------------------------------------------------
20 yb = yb + 0.5_dp
kyb = kyb + 1
IF (kyb > 10) RETURN
ierr = 0

30 del = yb - zid
!-----------------------------------------------------------------------
!     MAKE DEL LARGE ENOUGH TO AVOID UNDERFLOW IN GENERATION OF POWERS
!-----------------------------------------------------------------------
IF (ABS(del) <= 1.0E-4) THEN
  yb = yb + 1.0E-4
  del = yb - zid
END IF
htol = 0.125_dp * tol
zz = CMPLX(zr, yb, KIND=dp)
CALL cexenz(zz, nub, 2, 1, yy, ierr, rbry, htol, elim, alim, icdim, ca, cb)
IF (ierr == 6) GO TO 20
!-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION VIA TAYLOR SERIES FOR ORDER NUB
!-----------------------------------------------------------------------
iy = del + del
yt = del / (iy+1)
sum = yy(1)
htol = 0.25_dp * tol
trm = sum
cezt = CMPLX(COS(yt), -SIN(yt), KIND=dp)
zw = cone
zp(1) = cone
fk = 1.0_dp
fj = nub - 1
zt = cone / zz
!-----------------------------------------------------------------------
!     TERMS OF THE SERIES TRM=E(NUB-K,ZZ)*(YT**K)/K!,  K=0,1,... ARE
!     GENERATED BY BACKWARD RECURRENCE.  E IS SCALED BY CEXP(ZZ).
!-----------------------------------------------------------------------
DO  k = 2, 64
  rw = yt / fk
  zw = zw * CMPLX(0.0_dp, rw, KIND=dp)
  rw = fj * rw
  trm = (zw - CMPLX(0.0_dp, rw, KIND=dp)*trm) * zt
  sum = sum + trm
  zp(k) = zw
  rq(k) = rw
  asum = ABS(sum)
  atrm = ABS(trm)
  IF (atrm < htol*asum) GO TO 50
  fk = fk + 1.0_dp
  fj = fj - 1.0_dp
END DO
k = 64

50 kmax = k
sum = sum * cezt
IF (iy /= 0) THEN
  DO  i = 1, iy
    rzi = (iy-i+1) * yt + zid
    zz = CMPLX(zr, rzi, KIND=dp)
    zt = cone / zz
    trm = sum
    DO  k = 2, kmax
      trm = (zp(k) - CMPLX(0.0_dp, rq(k), KIND=dp)*trm) * zt
      sum = sum + trm
    END DO
    atrm = ABS(trm)
    asum = ABS(sum)
    xtol = htol * asum
    IF (atrm >= xtol) THEN
      IF (kmax < 64) THEN
        kmax = kmax + 1
        DO  k = kmax, 64
          rw = yt / fk
          zw = zw * CMPLX(0.0_dp, rw, KIND=dp)
          rw = fj * rw
          trm = (zw - CMPLX(0.0_dp, rw, KIND=dp)*trm) * zt
          sum = sum + trm
          zp(k) = zw
          rq(k) = rw
          atrm = ABS(trm)
          IF (atrm < xtol) GO TO 80
          fk = fk + 1.0_dp
          fj = fj - 1.0_dp
        END DO
        k = 64

        80 kmax = k
      END IF
    END IF
    sum = sum * cezt
  END DO
END IF
!-----------------------------------------------------------------------
!     FORM SEQUENCE UP OR DOWN FROM ORDER NUB
!-----------------------------------------------------------------------
IF (zi < 0.0_dp) sum = CONJG(sum)
cex = cone
!-----------------------------------------------------------------------
!     SCALE NEAR OVERFLOW LIMIT ON KODE=1
!-----------------------------------------------------------------------
IF (kode /= 2) THEN
  IF (ABS(zr) >= alim) THEN
    IF (ABS(zr) > elim) GO TO 130
    tola = EXP(alim-elim)
    rtola = 1.0_dp / tola
    scle = tola
    rscle = rtola
  END IF
  cex = scle * EXP(-z)
END IF
trm = sum * cex
y(kb) = trm * rscle
trms = trm
IF (n == 1 .AND. nflg /= 1) RETURN
IF (nflg /= 2) THEN
  fk = nub
  IF (nflg == 1 .AND. nb /= 0) THEN
    DO  k = 1, nb
      trm = (cex - z*trm) / fk
      fk = fk + 1.0_dp
    END DO
  END IF
  y(kb) = trm * rscle
  trms = trm
  IF (n == 1) RETURN
  DO  k = ks, kl
    trm = (cex - z*trm) / fk
    y(k) = trm * rscle
    fk = fk + 1.0_dp
  END DO
  IF (nflg == 1) RETURN
END IF
k = kb - 1
fk = nub - 1
zt = cone / z
DO  i = is, il
  trms = (cex - fk*trms) * zt
  y(k) = trms * rscle
  k = k - 1
  fk = fk - 1.0_dp
END DO
RETURN

130 ierr = 3
RETURN
END SUBROUTINE cacexi



FUNCTION psixn(n) RESULT(fn_val)
!***BEGIN PROLOGUE  PSIXN
!***REFER TO  EXINT,BSKIN

!  THIS SUBROUTINE RETURNS VALUES OF PSI(X)=DERIVATIVE OF LOG GAMMA(X),
!  X > 0.0 AT INTEGER ARGUMENTS.  A TABLE LOOK-UP IS PERFORMED FOR N <= 100,
!  AND THE ASYMPTOTIC EXPANSION IS EVALUATED FOR N > 100.

!***END PROLOGUE  PSIXN

INTEGER, INTENT(IN)  :: n
REAL (dp)            :: fn_val

INTEGER    :: k
REAL (dp)  :: ax, fn, rfn2, trm, s, wdtol
!-----------------------------------------------------------------------
!             PSIXN(N), N = 1,100
!-----------------------------------------------------------------------
REAL (dp), PARAMETER  :: c(100) = (/ -5.77215664901532861D-01,  &
  4.22784335098467139D-01, 9.22784335098467139D-01, 1.25611766843180047_dp, &
  1.50611766843180047_dp, 1.70611766843180047_dp, 1.87278433509846714_dp,  &
  2.01564147795561000_dp, 2.14064147795561000_dp, 2.25175258906672111_dp,  &
  2.35175258906672111_dp, 2.44266167997581202_dp, 2.52599501330914535_dp,  &
  2.60291809023222227_dp, 2.67434666166079370_dp, 2.74101332832746037_dp,  &
  2.80351332832746037_dp, 2.86233685773922507_dp, 2.91789241329478063_dp,  &
  2.97052399224214905_dp, 3.02052399224214905_dp, 3.06814303986119667_dp,  &
  3.11359758531574212_dp, 3.15707584618530734_dp, 3.19874251285197401_dp,  &
  3.23874251285197401_dp, 3.27720405131351247_dp, 3.31424108835054951_dp,  &
  3.34995537406483522_dp, 3.38443813268552488_dp, 3.41777146601885821_dp,  &
  3.45002953053498724_dp, 3.48127953053498724_dp, 3.51158256083801755_dp,  &
  3.54099432554389990_dp, 3.56956575411532847_dp, 3.59734353189310625_dp,  &
  3.62437055892013327_dp, 3.65068634839381748_dp, 3.67632737403484313_dp,  &
  3.70132737403484313_dp, 3.72571761793728215_dp, 3.74952714174680596_dp,  &
  3.77278295570029433_dp, 3.79551022842756706_dp, 3.81773245064978928_dp,  &
  3.83947158108457189_dp, 3.86074817682925274_dp, 3.88158151016258607_dp,  &
  3.90198967342789220_dp, 3.92198967342789220_dp, 3.94159751656514710_dp,  &
  3.96082828579591633_dp, 3.97969621032421822_dp, 3.99821472884273674_dp,  &
  4.01639654702455492_dp, 4.03425368988169777_dp, 4.05179754953082058_dp,  &
  4.06903892884116541_dp, 4.08598808138353829_dp, 4.10265474805020496_dp,  &
  4.11904819067315578_dp, 4.13517722293122029_dp, 4.15105023880423617_dp,  &
  4.16667523880423617_dp, 4.18205985418885155_dp, 4.19721136934036670_dp,  &
  4.21213674247469506_dp, 4.22684262482763624_dp, 4.24133537845082464_dp,  &
  4.25562109273653893_dp, 4.26970559977879245_dp, 4.28359448866768134_dp,  &
  4.29729311880466764_dp, 4.31080663231818115_dp, 4.32413996565151449_dp,  &
  4.33729786038835659_dp, 4.35028487337536958_dp, 4.36310538619588240_dp,  &
  4.37576361404398366_dp, 4.38826361404398366_dp, 4.40060929305632934_dp,  &
  4.41280441500754886_dp, 4.42485260777863319_dp, 4.43675736968339510_dp,  &
  4.44852207556574804_dp, 4.46014998254249223_dp, 4.47164423541605544_dp,  &
  4.48300787177969181_dp, 4.49424382683587158_dp, 4.50535493794698269_dp,  &
  4.51634394893599368_dp, 4.52721351415338499_dp, 4.53796620232542800_dp,  &
  4.54860450019776842_dp, 4.55913081598724211_dp, 4.56954748265390877_dp,  &
  4.57985676100442424_dp, 4.59006084263707730_dp, 4.60016185273808740_dp /)
!-----------------------------------------------------------------------
!             COEFFICIENTS OF ASYMPTOTIC EXPANSION
!-----------------------------------------------------------------------
REAL (dp), PARAMETER  :: b(6) =  &
  (/ 8.33333333333333333D-02, -8.33333333333333333D-03, 3.96825396825396825D-03,   &
    -4.16666666666666666D-03, 7.57575757575757576D-03, -2.10927960927960928D-02 /)

IF (n <= 100) THEN
  fn_val = c(n)
  RETURN
END IF
wdtol = MAX(EPSILON(0.0_dp), 1.0D-18)
fn = n
ax = 1.0_dp
s = -0.5_dp / fn
IF (ABS(s) > wdtol) THEN
  rfn2 = 1.0_dp / (fn*fn)
  DO  k = 1, 6
    ax = ax * rfn2
    trm = -b(k) * ax
    IF (ABS(trm) < wdtol) EXIT
    s = s + trm
  END DO
END IF

fn_val = s + LOG(fn)
RETURN
END FUNCTION psixn

END MODULE Exponential_Integral





