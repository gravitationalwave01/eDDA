      SUBROUTINE CAXPY(N,CA,CX,INCX,CY,INCY)
      USE DDPRECISION,ONLY : WP

!     constant times a vector plus a vector:
! given:
!       vectors CX,CY
!       scalar CA
! computes
!     CY = CY + CA*CX
!
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)

      COMPLEX (WP) :: CX(*),CY(*),CA
      INTEGER :: I,INCX,INCY,IX,IY,N

      IF (N<=0) RETURN
      IF (ABS(REAL(CA))+ABS(AIMAG(CA))==0.0_WP) RETURN
      IF (INCX==1 .AND. INCY==1) GO TO 20

!        code for unequal increments or equal increments
!          not equal to 1

      IX = 1
      IY = 1
      IF (INCX<0) IX = (-N+1)*INCX+1
      IF (INCY<0) IY = (-N+1)*INCY+1
      DO I = 1,N
        CY(IY) = CY(IY)+CA*CX(IX)
        IX = IX+INCX
        IY = IY+INCY
      END DO
      RETURN

!        code for both increments equal to 1

20    DO I = 1,N
        CY(I) = CY(I)+CA*CX(I)
      END DO
      RETURN
    END SUBROUTINE CAXPY
    SUBROUTINE CSCAL(N,CA,CX,INCX)
      USE DDPRECISION,ONLY : WP

!     scales a vector by a constant.
!     jack dongarra, 3/11/78.
!     modified to correct problem with negative increment, 8/21/90.

      COMPLEX (WP) :: CA,CX(1)
      INTEGER :: I,INCX,IX,N

      IF (N<=0) RETURN
      IF (INCX==1) GO TO 20

!        code for increment not equal to 1

      IX = 1
      IF (INCX<0) IX = (-N+1)*INCX+1
      DO I = 1,N
        CX(IX) = CA*CX(IX)
        IX = IX+INCX
      END DO
      RETURN

!        code for increment equal to 1

20    DO I = 1,N
        CX(I) = CA*CX(I)
      END DO
      RETURN
    END SUBROUTINE CSCAL
    SUBROUTINE CSWAP(N,CX,INCX,CY,INCY)
      USE DDPRECISION,ONLY : WP

!     interchanges two vectors.
!     jack dongarra, 3/11/78.

      COMPLEX (WP) :: CX(1),CY(1),CTEMP
      INTEGER :: I,IX,IY,N,INCX,INCY

      IF (N<=0) RETURN
      IF (INCX==1 .AND. INCY==1) GO TO 20

!       code for unequal increments or equal increments not equal
!         to 1

      IX = 1
      IY = 1
      IF (INCX<0) IX = (-N+1)*INCX+1
      IF (INCY<0) IY = (-N+1)*INCY+1
      DO I = 1,N
        CTEMP = CX(IX)
        CX(IX) = CY(IY)
        CY(IY) = CTEMP
        IX = IX+INCX
        IY = IY+INCY
      END DO
      RETURN

!       code for both increments equal to 1
20    DO I = 1,N
        CTEMP = CX(I)
        CX(I) = CY(I)
        CY(I) = CTEMP
      END DO
      RETURN
      END SUBROUTINE CSWAP

      INTEGER FUNCTION ICAMAX(N,CX,INCX)
      USE DDPRECISION,ONLY : WP

!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)

      COMPLEX (WP) :: CX(*)
      REAL (WP) :: SMAX
      INTEGER :: I,INCX,IX,N
      COMPLEX (WP) :: ZDUM
      REAL (WP) :: CABS1

      CABS1(ZDUM) = ABS(REAL(ZDUM))+ABS(AIMAG(ZDUM))

      ICAMAX = 0
      IF (N<1 .OR. INCX<=0) RETURN
      ICAMAX = 1
      IF (N==1) RETURN
      IF (INCX==1) GO TO 20

!        code for increment not equal to 1

      IX = 1
      SMAX = CABS1(CX(1))
      IX = IX+INCX
      DO I = 2,N
        IF (CABS1(CX(IX))<=SMAX) GO TO 5
        ICAMAX = I
        SMAX = CABS1(CX(IX))
5       IX = IX+INCX
      END DO
      RETURN

!        code for increment equal to 1

20    SMAX = CABS1(CX(1))
      DO I = 2,N
        IF (CABS1(CX(I))<=SMAX) GO TO 30
        ICAMAX = I
        SMAX = CABS1(CX(I))
30    END DO
      RETURN
      END FUNCTION ICAMAX

!---------------------------------------------------------------------------

      INTEGER FUNCTION ISAMAX(N,SX,INCX)
      USE DDPRECISION,ONLY : WP

!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)

      REAL (WP) :: SX(*),SMAX
      INTEGER :: I,INCX,IX,N

      ISAMAX = 0
      IF (N<1 .OR. INCX<=0) RETURN
      ISAMAX = 1
      IF (N==1) RETURN
      IF (INCX==1) GO TO 20

!        code for increment not equal to 1

      IX = 1
      SMAX = ABS(SX(1))
      IX = IX+INCX
      DO I = 2,N
        IF (ABS(SX(IX))<=SMAX) GO TO 5
        ISAMAX = I
        SMAX = ABS(SX(IX))
5       IX = IX+INCX
      END DO
      RETURN

!        code for increment equal to 1

20    SMAX = ABS(SX(1))
      DO I = 2,N
        IF (ABS(SX(I))<=SMAX) GO TO 30
        ISAMAX = I
        SMAX = ABS(SX(I))
30    END DO
      RETURN
    END FUNCTION ISAMAX
    LOGICAL FUNCTION LSAME(CA,CB)

!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994

!     .. Scalar Arguments ..
      CHARACTER :: CA,CB
!     ..

!  Purpose
!  =======

!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.

!  Arguments
!  =========

!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.

! =====================================================================

!     .. Intrinsic Functions ..
      INTRINSIC ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER :: INTA,INTB,ZCODE
!     ..
!     .. Executable Statements ..

!     Test if the characters are equal

      LSAME = CA == CB
      IF (LSAME) RETURN

!     Now test for equivalence if both characters are alphabetic.

      ZCODE = ICHAR('Z')

!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.

      INTA = ICHAR(CA)
      INTB = ICHAR(CB)

      IF (ZCODE==90 .OR. ZCODE==122) THEN

!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.

        IF (INTA>=97 .AND. INTA<=122) INTA = INTA - 32
        IF (INTB>=97 .AND. INTB<=122) INTB = INTB - 32

      ELSE IF (ZCODE==233 .OR. ZCODE==169) THEN

!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.

        IF (INTA>=129 .AND. INTA<=137 .OR. INTA>=145 .AND. INTA<=153 .OR. &
          INTA>=162 .AND. INTA<=169) INTA = INTA+64
        IF (INTB>=129 .AND. INTB<=137 .OR. INTB>=145 .AND. INTB<=153 .OR. &
          INTB>=162 .AND. INTB<=169) INTB = INTB+64

      ELSE IF (ZCODE==218 .OR. ZCODE==250) THEN

!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.

        IF (INTA>=225 .AND. INTA<=250) INTA = INTA - 32
        IF (INTB>=225 .AND. INTB<=250) INTB = INTB - 32
      END IF
      LSAME = INTA == INTB

!     RETURN

!     End of LSAME

    END FUNCTION LSAME
    SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
      USE DDPRECISION,ONLY : WP

!     constant times a vector plus a vector.
!     uses unrolled loop for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)

      REAL (WP) :: SX(*),SY(*),SA
      INTEGER :: I,INCX,INCY,IX,IY,M,MP1,N

      IF (N<=0) RETURN
      IF (SA==0.0_WP) RETURN
      IF (INCX==1 .AND. INCY==1) GO TO 20

!        code for unequal increments or equal increments
!          not equal to 1

      IX = 1
      IY = 1
      IF (INCX<0) IX = (-N+1)*INCX+1
      IF (INCY<0) IY = (-N+1)*INCY+1
      DO I = 1,N
        SY(IY) = SY(IY)+SA*SX(IX)
        IX = IX+INCX
        IY = IY+INCY
      END DO
      RETURN

!        code for both increments equal to 1


!        clean-up loop

20    M = MOD(N,4)
      IF (M==0) GO TO 40
      DO I = 1,M
        SY(I) = SY(I)+SA*SX(I)
      END DO
      IF (N<4) RETURN
40    MP1 = M+1
      DO I = MP1,N,4
        SY(I) = SY(I)+SA*SX(I)
        SY(I+1) = SY(I+1)+SA*SX(I+1)
        SY(I+2) = SY(I+2)+SA*SX(I+2)
        SY(I+3) = SY(I+3)+SA*SX(I+3)
      END DO
      RETURN
    END SUBROUTINE SAXPY
    SUBROUTINE SCOPY(N,SX,INCX,SY,INCY)
      USE DDPRECISION,ONLY : WP

!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to 1.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)

      REAL (WP) :: SX(*),SY(*)
      INTEGER :: I,INCX,INCY,IX,IY,M,MP1,N

      IF (N<=0) RETURN
      IF (INCX==1 .AND. INCY==1) GO TO 20

!        code for unequal increments or equal increments
!          not equal to 1

      IX = 1
      IY = 1
      IF (INCX<0) IX = (-N+1)*INCX+1
      IF (INCY<0) IY = (-N+1)*INCY+1
      DO I = 1,N
        SY(IY) = SX(IX)
        IX = IX+INCX
        IY = IY+INCY
      END DO
      RETURN

!        code for both increments equal to 1


!        clean-up loop

20    M = MOD(N,7)
      IF (M==0) GO TO 40
      DO I = 1,M
        SY(I) = SX(I)
      END DO
      IF (N<7) RETURN
40    MP1 = M+1
      DO I = MP1,N,7
        SY(I) = SX(I)
        SY(I+1) = SX(I+1)
        SY(I+2) = SX(I+2)
        SY(I+3) = SX(I+3)
        SY(I+4) = SX(I+4)
        SY(I+5) = SX(I+5)
        SY(I+6) = SX(I+6)
      END DO
      RETURN
    END SUBROUTINE SCOPY
    FUNCTION SDOT(N,SX,INCX,SY,INCY)
      USE DDPRECISION,ONLY : WP
      REAL (WP) :: SDOT

!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)

      REAL (WP) :: SX(*),SY(*),STEMP
      INTEGER :: I,INCX,INCY,IX,IY,M,MP1,N

      STEMP = 0.0E0_WP
      SDOT = 0.0E0_WP
      IF (N<=0) RETURN
      IF (INCX==1 .AND. INCY==1) GO TO 20

!        code for unequal increments or equal increments
!          not equal to 1

      IX = 1
      IY = 1
      IF (INCX<0) IX = (-N+1)*INCX+1
      IF (INCY<0) IY = (-N+1)*INCY+1
      DO I = 1,N
        STEMP = STEMP+SX(IX)*SY(IY)
        IX = IX+INCX
        IY = IY+INCY
      END DO
      SDOT = STEMP
      RETURN

!        code for both increments equal to 1


!        clean-up loop

20    M = MOD(N,5)
      IF (M==0) GO TO 40
      DO I = 1,M
        STEMP = STEMP+SX(I)*SY(I)
      END DO
      IF (N<5) GO TO 60
40    MP1 = M+1
      DO I = MP1,N,5
        STEMP = STEMP+SX(I)*SY(I)+SX(I+1)*SY(I+1)+SX(I+2)*SY(I+2)+ &
          SX(I+3)*SY(I+3)+SX(I+4)*SY(I+4)
      END DO
60    SDOT = STEMP
      RETURN
    END FUNCTION SDOT
    SUBROUTINE SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      USE DDPRECISION,ONLY : WP
!     .. Scalar Arguments ..
      CHARACTER (1) :: TRANSA,TRANSB
      INTEGER :: M,N,K,LDA,LDB,LDC
      REAL (WP) :: ALPHA,BETA
!     .. Array Arguments ..
      REAL (WP) :: A(LDA,*),B(LDB,*),C(LDC,*)
!     ..

!  Purpose
!  =======

!  SGEMM  performs one of the matrix-matrix operations

!     C := alpha*op( A )*op( B ) + beta*C,

!  where  op( X ) is one of

!     op( X ) = X   or   op( X ) = X',

!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.

!  Parameters
!  ==========

!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:

!              TRANSA = 'N' or 'n',  op( A ) = A.

!              TRANSA = 'T' or 't',  op( A ) = A'.

!              TRANSA = 'C' or 'c',  op( A ) = A'.

!           Unchanged on exit.

!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:

!              TRANSB = 'N' or 'n',  op( B ) = B.

!              TRANSB = 'T' or 't',  op( B ) = B'.

!              TRANSB = 'C' or 'c',  op( B ) = B'.

!           Unchanged on exit.

!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.

!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.

!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.

!  ALPHA  - REAL            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.

!  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.

!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.

!  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.

!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.

!  BETA   - REAL            .
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.

!  C      - REAL             array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).

!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.


!  Level 3 Blas routine.

!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.


!     .. External Functions ..
      LOGICAL :: LSAME
      EXTERNAL LSAME
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     .. Local Scalars ..
      LOGICAL :: NOTA,NOTB
      INTEGER :: I,INFO,J,L,NCOLA,NROWA,NROWB
      REAL (WP) :: TEMP
!     .. Parameters ..
      REAL (WP) :: ONE,ZERO
      PARAMETER (ONE=1.0E+0_WP,ZERO=0.0E+0_WP)
!     ..
!     .. Executable Statements ..

!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.

      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      IF (NOTA) THEN
        NROWA = M
        NCOLA = K
      ELSE
        NROWA = K
        NCOLA = M
      END IF
      IF (NOTB) THEN
        NROWB = K
      ELSE
        NROWB = N
      END IF

!     Test the input parameters.

      INFO = 0
      IF (( .NOT. NOTA) .AND. ( .NOT. LSAME(TRANSA,'C')) .AND. ( .NOT. LSAME( &
          TRANSA,'T'))) THEN
        INFO = 1
      ELSE IF (( .NOT. NOTB) .AND. ( .NOT. LSAME(TRANSB, &
          'C')) .AND. ( .NOT. LSAME(TRANSB,'T'))) THEN
        INFO = 2
      ELSE IF (M<0) THEN
        INFO = 3
      ELSE IF (N<0) THEN
        INFO = 4
      ELSE IF (K<0) THEN
        INFO = 5
      ELSE IF (LDA<MAX(1,NROWA)) THEN
        INFO = 8
      ELSE IF (LDB<MAX(1,NROWB)) THEN
        INFO = 10
      ELSE IF (LDC<MAX(1,M)) THEN
        INFO = 13
      END IF
      IF (INFO/=0) THEN
        CALL XERBLA('SGEMM ',INFO)
        RETURN
      END IF

!     Quick return if possible.

      IF ((M==0) .OR. (N==0) .OR. (((ALPHA==ZERO) .OR. (K==0)) .AND. (BETA== &
        ONE))) RETURN

!     And if  alpha.eq.zero.

      IF (ALPHA==ZERO) THEN
        IF (BETA==ZERO) THEN
          DO J = 1,N
            DO I = 1,M
              C(I,J) = ZERO
            END DO
          END DO
        ELSE
          DO J = 1,N
            DO I = 1,M
              C(I,J) = BETA*C(I,J)
            END DO
          END DO
        END IF
        RETURN
      END IF

!     Start the operations.

      IF (NOTB) THEN
        IF (NOTA) THEN

!           Form  C := alpha*A*B + beta*C.

          DO J = 1,N
            IF (BETA==ZERO) THEN
              DO I = 1,M
                C(I,J) = ZERO
              END DO
            ELSE IF (BETA/=ONE) THEN
              DO I = 1,M
                C(I,J) = BETA*C(I,J)
              END DO
            END IF
            DO L = 1,K
              IF (B(L,J)/=ZERO) THEN
                TEMP = ALPHA*B(L,J)
                DO I = 1,M
                  C(I,J) = C(I,J)+TEMP*A(I,L)
                END DO
              END IF
            END DO
          END DO
        ELSE

!           Form  C := alpha*A'*B + beta*C

          DO J = 1,N
            DO I = 1,M
              TEMP = ZERO
              DO L = 1,K
                TEMP = TEMP+A(L,I)*B(L,J)
              END DO
              IF (BETA==ZERO) THEN
                C(I,J) = ALPHA*TEMP
              ELSE
                C(I,J) = ALPHA*TEMP+BETA*C(I,J)
              END IF
            END DO
          END DO
        END IF
      ELSE
        IF (NOTA) THEN

!           Form  C := alpha*A*B' + beta*C

          DO J = 1,N
            IF (BETA==ZERO) THEN
              DO I = 1,M
                C(I,J) = ZERO
              END DO
            ELSE IF (BETA/=ONE) THEN
              DO I = 1,M
                C(I,J) = BETA*C(I,J)
              END DO
            END IF
            DO L = 1,K
              IF (B(J,L)/=ZERO) THEN
                TEMP = ALPHA*B(J,L)
                DO I = 1,M
                  C(I,J) = C(I,J)+TEMP*A(I,L)
                END DO
              END IF
            END DO
          END DO
        ELSE

!           Form  C := alpha*A'*B' + beta*C

          DO J = 1,N
            DO I = 1,M
              TEMP = ZERO
              DO L = 1,K
                TEMP = TEMP+A(L,I)*B(J,L)
              END DO
              IF (BETA==ZERO) THEN
                C(I,J) = ALPHA*TEMP
              ELSE
                C(I,J) = ALPHA*TEMP+BETA*C(I,J)
              END IF
            END DO
          END DO
        END IF
      END IF

      RETURN

!     End of SGEMM .

    END SUBROUTINE SGEMM
    SUBROUTINE SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      USE DDPRECISION,ONLY : WP
!     .. Scalar Arguments ..
      REAL (WP) :: ALPHA,BETA
      INTEGER :: INCX,INCY,LDA,M,N
      CHARACTER (1) :: TRANS
!     .. Array Arguments ..
      REAL (WP) :: A(LDA,*),X(*),Y(*)
!     ..

!  Purpose
!  =======

!  SGEMV  performs one of the matrix-vector operations

!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,

!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.

!  Parameters
!  ==========

!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:

!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.

!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.

!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.

!           Unchanged on exit.

!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.

!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.

!  ALPHA  - REAL            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.

!  A      - REAL             array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.

!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.

!  X      - REAL             array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.

!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.

!  BETA   - REAL            .
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.

!  Y      - REAL             array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.

!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.


!  Level 2 Blas routine.

!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.


!     .. Parameters ..
      REAL (WP) :: ONE,ZERO
      PARAMETER (ONE=1.0E+0_WP,ZERO=0.0E+0_WP)
!     .. Local Scalars ..
      REAL (WP) :: TEMP
      INTEGER :: I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
!     .. External Functions ..
      LOGICAL :: LSAME
      EXTERNAL LSAME
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

      INFO = 0
      IF ( .NOT. LSAME(TRANS,'N') .AND. .NOT. LSAME(TRANS,'T') .AND. &
          .NOT. LSAME(TRANS,'C')) THEN
        INFO = 1
      ELSE IF (M<0) THEN
        INFO = 2
      ELSE IF (N<0) THEN
        INFO = 3
      ELSE IF (LDA<MAX(1,M)) THEN
        INFO = 6
      ELSE IF (INCX==0) THEN
        INFO = 8
      ELSE IF (INCY==0) THEN
        INFO = 11
      END IF
      IF (INFO/=0) THEN
        CALL XERBLA('SGEMV ',INFO)
        RETURN
      END IF

!     Quick return if possible.

      IF ((M==0) .OR. (N==0) .OR. ((ALPHA==ZERO) .AND. (BETA==ONE))) RETURN

!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.

      IF (LSAME(TRANS,'N')) THEN
        LENX = N
        LENY = M
      ELSE
        LENX = M
        LENY = N
      END IF
      IF (INCX>0) THEN
        KX = 1
      ELSE
        KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY>0) THEN
        KY = 1
      ELSE
        KY = 1 - (LENY-1)*INCY
      END IF

!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.

!     First form  y := beta*y.

      IF (BETA/=ONE) THEN
        IF (INCY==1) THEN
          IF (BETA==ZERO) THEN
            DO I = 1,LENY
              Y(I) = ZERO
            END DO
          ELSE
            DO I = 1,LENY
              Y(I) = BETA*Y(I)
            END DO
          END IF
        ELSE
          IY = KY
          IF (BETA==ZERO) THEN
            DO I = 1,LENY
              Y(IY) = ZERO
              IY = IY+INCY
            END DO
          ELSE
            DO I = 1,LENY
              Y(IY) = BETA*Y(IY)
              IY = IY+INCY
            END DO
          END IF
        END IF
      END IF
      IF (ALPHA==ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN

!        Form  y := alpha*A*x + y.

        JX = KX
        IF (INCY==1) THEN
          DO J = 1,N
            IF (X(JX)/=ZERO) THEN
              TEMP = ALPHA*X(JX)
              DO I = 1,M
                Y(I) = Y(I)+TEMP*A(I,J)
              END DO
            END IF
            JX = JX+INCX
          END DO
        ELSE
          DO J = 1,N
            IF (X(JX)/=ZERO) THEN
              TEMP = ALPHA*X(JX)
              IY = KY
              DO I = 1,M
                Y(IY) = Y(IY)+TEMP*A(I,J)
                IY = IY+INCY
              END DO
            END IF
            JX = JX+INCX
          END DO
        END IF
      ELSE

!        Form  y := alpha*A'*x + y.

        JY = KY
        IF (INCX==1) THEN
          DO J = 1,N
            TEMP = ZERO
            DO I = 1,M
              TEMP = TEMP+A(I,J)*X(I)
            END DO
            Y(JY) = Y(JY)+ALPHA*TEMP
            JY = JY+INCY
          END DO
        ELSE
          DO J = 1,N
            TEMP = ZERO
            IX = KX
            DO I = 1,M
              TEMP = TEMP+A(I,J)*X(IX)
              IX = IX+INCX
            END DO
            Y(JY) = Y(JY)+ALPHA*TEMP
            JY = JY+INCY
          END DO
        END IF
      END IF

      RETURN

!     End of SGEMV .

    END SUBROUTINE SGEMV
    SUBROUTINE SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      USE DDPRECISION,ONLY : WP
!     .. Scalar Arguments ..
      REAL (WP) :: ALPHA
      INTEGER :: INCX,INCY,LDA,M,N
!     .. Array Arguments ..
      REAL (WP) :: A(LDA,*),X(*),Y(*)
!     ..

!  Purpose
!  =======

!  SGER   performs the rank 1 operation

!     A := alpha*x*y' + A,

!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.

!  Parameters
!  ==========

!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.

!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.

!  ALPHA  - REAL            .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.

!  X      - REAL             array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.

!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.

!  Y      - REAL             array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.

!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.

!  A      - REAL             array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.

!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.


!  Level 2 Blas routine.

!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.


!     .. Parameters ..
      REAL (WP) :: ZERO
      PARAMETER (ZERO=0.0E+0_WP)
!     .. Local Scalars ..
      REAL (WP) :: TEMP
      INTEGER :: I,INFO,IX,J,JY,KX
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

      INFO = 0
      IF (M<0) THEN
        INFO = 1
      ELSE IF (N<0) THEN
        INFO = 2
      ELSE IF (INCX==0) THEN
        INFO = 5
      ELSE IF (INCY==0) THEN
        INFO = 7
      ELSE IF (LDA<MAX(1,M)) THEN
        INFO = 9
      END IF
      IF (INFO/=0) THEN
        CALL XERBLA('SGER  ',INFO)
        RETURN
      END IF

!     Quick return if possible.

      IF ((M==0) .OR. (N==0) .OR. (ALPHA==ZERO)) RETURN

!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.

      IF (INCY>0) THEN
        JY = 1
      ELSE
        JY = 1 - (N-1)*INCY
      END IF
      IF (INCX==1) THEN
        DO J = 1,N
          IF (Y(JY)/=ZERO) THEN
            TEMP = ALPHA*Y(JY)
            DO I = 1,M
              A(I,J) = A(I,J)+X(I)*TEMP
            END DO
          END IF
          JY = JY+INCY
        END DO
      ELSE
        IF (INCX>0) THEN
          KX = 1
        ELSE
          KX = 1 - (M-1)*INCX
        END IF
        DO J = 1,N
          IF (Y(JY)/=ZERO) THEN
            TEMP = ALPHA*Y(JY)
            IX = KX
            DO I = 1,M
              A(I,J) = A(I,J)+X(IX)*TEMP
              IX = IX+INCX
            END DO
          END IF
          JY = JY+INCY
        END DO
      END IF

      RETURN

!     End of SGER  .

    END SUBROUTINE SGER
    FUNCTION SNRM2(N,X,INCX)
      USE DDPRECISION,ONLY : WP
      REAL (WP) :: SNRM2
!     .. Scalar Arguments ..
      INTEGER :: INCX,N
!     .. Array Arguments ..
      REAL (WP) :: X(*)
!     ..

!  SNRM2 returns the euclidean norm of a vector via the function
!  name, so that

!     SNRM2 := sqrt( x'*x )



!  -- This version written on 25-October-1982.
!     Modified on 14-October-1993 to inline the call to SLASSQ.
!     Sven Hammarling, Nag Ltd.


!     .. Parameters ..
      REAL (WP) :: ONE,ZERO
      PARAMETER (ONE=1.0E+0_WP,ZERO=0.0E+0_WP)
!     .. Local Scalars ..
      INTEGER :: IX
      REAL (WP) :: ABSXI,NORM,SCALE,SSQ
!     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
!     ..
!     .. Executable Statements ..
      IF (N<1 .OR. INCX<1) THEN
        NORM = ZERO
      ELSE IF (N==1) THEN
        NORM = ABS(X(1))
      ELSE
        SCALE = ZERO
        SSQ = ONE
!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL SLASSQ( N,X,INCX,SCALE,SSQ )

        DO IX = 1,1+(N-1)*INCX,INCX
          IF (X(IX)/=ZERO) THEN
            ABSXI = ABS(X(IX))
            IF (SCALE<ABSXI) THEN
              SSQ = ONE+SSQ*(SCALE/ABSXI)**2
              SCALE = ABSXI
            ELSE
              SSQ = SSQ+(ABSXI/SCALE)**2
            END IF
          END IF
        END DO
        NORM = SCALE*SQRT(SSQ)
      END IF

      SNRM2 = NORM
      RETURN

!     End of SNRM2.

    END FUNCTION SNRM2
    SUBROUTINE SROT(N,SX,INCX,SY,INCY,C,S)
      USE DDPRECISION,ONLY : WP

!     applies a plane rotation.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)

      REAL (WP) :: SX(*),SY(*),STEMP,C,S
      INTEGER :: I,INCX,INCY,IX,IY,N

      IF (N<=0) RETURN
      IF (INCX==1 .AND. INCY==1) GO TO 20

!       code for unequal increments or equal increments not equal
!         to 1

      IX = 1
      IY = 1
      IF (INCX<0) IX = (-N+1)*INCX+1
      IF (INCY<0) IY = (-N+1)*INCY+1
      DO I = 1,N
        STEMP = C*SX(IX)+S*SY(IY)
        SY(IY) = C*SY(IY) - S*SX(IX)
        SX(IX) = STEMP
        IX = IX+INCX
        IY = IY+INCY
      END DO
      RETURN

!       code for both increments equal to 1

20    DO I = 1,N
        STEMP = C*SX(I)+S*SY(I)
        SY(I) = C*SY(I) - S*SX(I)
        SX(I) = STEMP
      END DO
      RETURN
    END SUBROUTINE SROT
    SUBROUTINE SSCAL(N,SA,SX,INCX)
      USE DDPRECISION,ONLY : WP

!     scales a vector by a constant.
!     uses unrolled loops for increment equal to 1.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)

      REAL (WP) :: SA,SX(*)
      INTEGER :: I,INCX,M,MP1,N,NINCX

      IF (N<=0 .OR. INCX<=0) RETURN
      IF (INCX==1) GO TO 20

!        code for increment not equal to 1

      NINCX = N*INCX
      DO I = 1,NINCX,INCX
        SX(I) = SA*SX(I)
      END DO
      RETURN

!        code for increment equal to 1


!        clean-up loop

20    M = MOD(N,5)
      IF (M==0) GO TO 40
      DO I = 1,M
        SX(I) = SA*SX(I)
      END DO
      IF (N<5) RETURN
40    MP1 = M+1
      DO I = MP1,N,5
        SX(I) = SA*SX(I)
        SX(I+1) = SA*SX(I+1)
        SX(I+2) = SA*SX(I+2)
        SX(I+3) = SA*SX(I+3)
        SX(I+4) = SA*SX(I+4)
      END DO
      RETURN
    END SUBROUTINE SSCAL
    SUBROUTINE SSWAP(N,SX,INCX,SY,INCY)
      USE DDPRECISION,ONLY : WP

!     interchanges two vectors.
!     uses unrolled loops for increments equal to 1.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)

      REAL (WP) :: SX(*),SY(*),STEMP
      INTEGER :: I,INCX,INCY,IX,IY,M,MP1,N

      IF (N<=0) RETURN
      IF (INCX==1 .AND. INCY==1) GO TO 20

!       code for unequal increments or equal increments not equal
!         to 1

      IX = 1
      IY = 1
      IF (INCX<0) IX = (-N+1)*INCX+1
      IF (INCY<0) IY = (-N+1)*INCY+1
      DO I = 1,N
        STEMP = SX(IX)
        SX(IX) = SY(IY)
        SY(IY) = STEMP
        IX = IX+INCX
        IY = IY+INCY
      END DO
      RETURN

!       code for both increments equal to 1


!       clean-up loop

20    M = MOD(N,3)
      IF (M==0) GO TO 40
      DO I = 1,M
        STEMP = SX(I)
        SX(I) = SY(I)
        SY(I) = STEMP
      END DO
      IF (N<3) RETURN
40    MP1 = M+1
      DO I = MP1,N,3
        STEMP = SX(I)
        SX(I) = SY(I)
        SY(I) = STEMP
        STEMP = SX(I+1)
        SX(I+1) = SY(I+1)
        SY(I+1) = STEMP
        STEMP = SX(I+2)
        SX(I+2) = SY(I+2)
        SY(I+2) = STEMP
      END DO
      RETURN
    END SUBROUTINE SSWAP
    SUBROUTINE STRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
      USE DDPRECISION,ONLY : WP
!     .. Scalar Arguments ..
      CHARACTER (1) :: SIDE,UPLO,TRANSA,DIAG
      INTEGER :: M,N,LDA,LDB
      REAL (WP) :: ALPHA
!     .. Array Arguments ..
      REAL (WP) :: A(LDA,*),B(LDB,*)
!     ..

!  Purpose
!  =======

!  STRMM  performs one of the matrix-matrix operations

!     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),

!  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of

!     op( A ) = A   or   op( A ) = A'.

!  Parameters
!  ==========

!  SIDE   - CHARACTER*1.
!           On entry,  SIDE specifies whether  op( A ) multiplies B from
!           the left or right as follows:

!              SIDE = 'L' or 'l'   B := alpha*op( A )*B.

!              SIDE = 'R' or 'r'   B := alpha*B*op( A ).

!           Unchanged on exit.

!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:

!              UPLO = 'U' or 'u'   A is an upper triangular matrix.

!              UPLO = 'L' or 'l'   A is a lower triangular matrix.

!           Unchanged on exit.

!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:

!              TRANSA = 'N' or 'n'   op( A ) = A.

!              TRANSA = 'T' or 't'   op( A ) = A'.

!              TRANSA = 'C' or 'c'   op( A ) = A'.

!           Unchanged on exit.

!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:

!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.

!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.

!           Unchanged on exit.

!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.

!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.

!  ALPHA  - REAL            .
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.

!  A      - REAL             array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.

!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.

!  B      - REAL             array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain the matrix  B,  and  on exit  is overwritten  by the
!           transformed matrix.

!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.


!  Level 3 Blas routine.

!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.


!     .. External Functions ..
      LOGICAL :: LSAME
      EXTERNAL LSAME
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     .. Local Scalars ..
      LOGICAL :: LSIDE,NOUNIT,UPPER
      INTEGER :: I,INFO,J,K,NROWA
      REAL (WP) :: TEMP
!     .. Parameters ..
      REAL (WP) :: ONE,ZERO
      PARAMETER (ONE=1.0E+0_WP,ZERO=0.0E+0_WP)
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
        NROWA = M
      ELSE
        NROWA = N
      END IF
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')

      INFO = 0
      IF (( .NOT. LSIDE) .AND. ( .NOT. LSAME(SIDE,'R'))) THEN
        INFO = 1
      ELSE IF (( .NOT. UPPER) .AND. ( .NOT. LSAME(UPLO,'L'))) THEN
        INFO = 2
      ELSE IF (( .NOT. LSAME(TRANSA,'N')) .AND. ( .NOT. LSAME(TRANSA, &
          'T')) .AND. ( .NOT. LSAME(TRANSA,'C'))) THEN
        INFO = 3
      ELSE IF (( .NOT. LSAME(DIAG,'U')) .AND. ( .NOT. LSAME(DIAG,'N'))) THEN
        INFO = 4
      ELSE IF (M<0) THEN
        INFO = 5
      ELSE IF (N<0) THEN
        INFO = 6
      ELSE IF (LDA<MAX(1,NROWA)) THEN
        INFO = 9
      ELSE IF (LDB<MAX(1,M)) THEN
        INFO = 11
      END IF
      IF (INFO/=0) THEN
        CALL XERBLA('STRMM ',INFO)
        RETURN
      END IF

!     Quick return if possible.

      IF (N==0) RETURN

!     And when  alpha.eq.zero.

      IF (ALPHA==ZERO) THEN
        DO J = 1,N
          DO I = 1,M
            B(I,J) = ZERO
          END DO
        END DO
        RETURN
      END IF

!     Start the operations.

      IF (LSIDE) THEN
        IF (LSAME(TRANSA,'N')) THEN

!           Form  B := alpha*A*B.

          IF (UPPER) THEN
            DO J = 1,N
              DO K = 1,M
                IF (B(K,J)/=ZERO) THEN
                  TEMP = ALPHA*B(K,J)
                  DO I = 1,K - 1
                    B(I,J) = B(I,J)+TEMP*A(I,K)
                  END DO
                  IF (NOUNIT) TEMP = TEMP*A(K,K)
                  B(K,J) = TEMP
                END IF
              END DO
            END DO
          ELSE
            DO J = 1,N
              DO K = M,1,-1
                IF (B(K,J)/=ZERO) THEN
                  TEMP = ALPHA*B(K,J)
                  B(K,J) = TEMP
                  IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                  DO I = K+1,M
                    B(I,J) = B(I,J)+TEMP*A(I,K)
                  END DO
                END IF
              END DO
            END DO
          END IF
        ELSE

!           Form  B := alpha*A'*B.

          IF (UPPER) THEN
            DO J = 1,N
              DO I = M,1,-1
                TEMP = B(I,J)
                IF (NOUNIT) TEMP = TEMP*A(I,I)
                DO K = 1,I - 1
                  TEMP = TEMP+A(K,I)*B(K,J)
                END DO
                B(I,J) = ALPHA*TEMP
              END DO
            END DO
          ELSE
            DO J = 1,N
              DO I = 1,M
                TEMP = B(I,J)
                IF (NOUNIT) TEMP = TEMP*A(I,I)
                DO K = I+1,M
                  TEMP = TEMP+A(K,I)*B(K,J)
                END DO
                B(I,J) = ALPHA*TEMP
              END DO
            END DO
          END IF
        END IF
      ELSE
        IF (LSAME(TRANSA,'N')) THEN

!           Form  B := alpha*B*A.

          IF (UPPER) THEN
            DO J = N,1,-1
              TEMP = ALPHA
              IF (NOUNIT) TEMP = TEMP*A(J,J)
              DO I = 1,M
                B(I,J) = TEMP*B(I,J)
              END DO
              DO K = 1,J - 1
                IF (A(K,J)/=ZERO) THEN
                  TEMP = ALPHA*A(K,J)
                  DO I = 1,M
                    B(I,J) = B(I,J)+TEMP*B(I,K)
                  END DO
                END IF
              END DO
            END DO
          ELSE
            DO J = 1,N
              TEMP = ALPHA
              IF (NOUNIT) TEMP = TEMP*A(J,J)
              DO I = 1,M
                B(I,J) = TEMP*B(I,J)
              END DO
              DO K = J+1,N
                IF (A(K,J)/=ZERO) THEN
                  TEMP = ALPHA*A(K,J)
                  DO I = 1,M
                    B(I,J) = B(I,J)+TEMP*B(I,K)
                  END DO
                END IF
              END DO
            END DO
          END IF
        ELSE

!           Form  B := alpha*B*A'.

          IF (UPPER) THEN
            DO K = 1,N
              DO J = 1,K - 1
                IF (A(J,K)/=ZERO) THEN
                  TEMP = ALPHA*A(J,K)
                  DO I = 1,M
                    B(I,J) = B(I,J)+TEMP*B(I,K)
                  END DO
                END IF
              END DO
              TEMP = ALPHA
              IF (NOUNIT) TEMP = TEMP*A(K,K)
              IF (TEMP/=ONE) THEN
                DO I = 1,M
                  B(I,K) = TEMP*B(I,K)
                END DO
              END IF
            END DO
          ELSE
            DO K = N,1,-1
              DO J = K+1,N
                IF (A(J,K)/=ZERO) THEN
                  TEMP = ALPHA*A(J,K)
                  DO I = 1,M
                    B(I,J) = B(I,J)+TEMP*B(I,K)
                  END DO
                END IF
              END DO
              TEMP = ALPHA
              IF (NOUNIT) TEMP = TEMP*A(K,K)
              IF (TEMP/=ONE) THEN
                DO I = 1,M
                  B(I,K) = TEMP*B(I,K)
                END DO
              END IF
            END DO
          END IF
        END IF
      END IF

      RETURN

!     End of STRMM .

    END SUBROUTINE STRMM
    SUBROUTINE STRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      USE DDPRECISION,ONLY : WP
!     .. Scalar Arguments ..
      INTEGER :: INCX,LDA,N
      CHARACTER (1) :: DIAG,TRANS,UPLO
!     .. Array Arguments ..
      REAL (WP) :: A(LDA,*),X(*)
!     ..

!  Purpose
!  =======

!  STRMV  performs one of the matrix-vector operations

!     x := A*x,   or   x := A'*x,

!  where x is an n element vector and  A is an n by n unit, or non-unit,
!  upper or lower triangular matrix.

!  Parameters
!  ==========

!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:

!              UPLO = 'U' or 'u'   A is an upper triangular matrix.

!              UPLO = 'L' or 'l'   A is a lower triangular matrix.

!           Unchanged on exit.

!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:

!              TRANS = 'N' or 'n'   x := A*x.

!              TRANS = 'T' or 't'   x := A'*x.

!              TRANS = 'C' or 'c'   x := A'*x.

!           Unchanged on exit.

!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:

!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.

!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.

!           Unchanged on exit.

!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.

!  A      - REAL             array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular matrix and the strictly lower triangular part of
!           A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular matrix and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced either, but are assumed to be unity.
!           Unchanged on exit.

!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.

!  X      - REAL             array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.

!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.


!  Level 2 Blas routine.

!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.


!     .. Parameters ..
      REAL (WP) :: ZERO
      PARAMETER (ZERO=0.0E+0_WP)
!     .. Local Scalars ..
      REAL (WP) :: TEMP
      INTEGER :: I,INFO,IX,J,JX,KX
      LOGICAL :: NOUNIT
!     .. External Functions ..
      LOGICAL :: LSAME
      EXTERNAL LSAME
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

      INFO = 0
      IF ( .NOT. LSAME(UPLO,'U') .AND. .NOT. LSAME(UPLO,'L')) THEN
        INFO = 1
      ELSE IF ( .NOT. LSAME(TRANS,'N') .AND. .NOT. LSAME(TRANS,'T') .AND. &
          .NOT. LSAME(TRANS,'C')) THEN
        INFO = 2
      ELSE IF ( .NOT. LSAME(DIAG,'U') .AND. .NOT. LSAME(DIAG,'N')) THEN
        INFO = 3
      ELSE IF (N<0) THEN
        INFO = 4
      ELSE IF (LDA<MAX(1,N)) THEN
        INFO = 6
      ELSE IF (INCX==0) THEN
        INFO = 8
      END IF
      IF (INFO/=0) THEN
        CALL XERBLA('STRMV ',INFO)
        RETURN
      END IF

!     Quick return if possible.

      IF (N==0) RETURN

      NOUNIT = LSAME(DIAG,'N')

!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.

      IF (INCX<=0) THEN
        KX = 1 - (N-1)*INCX
      ELSE IF (INCX/=1) THEN
        KX = 1
      END IF

!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.

      IF (LSAME(TRANS,'N')) THEN

!        Form  x := A*x.

        IF (LSAME(UPLO,'U')) THEN
          IF (INCX==1) THEN
            DO J = 1,N
              IF (X(J)/=ZERO) THEN
                TEMP = X(J)
                DO I = 1,J - 1
                  X(I) = X(I)+TEMP*A(I,J)
                END DO
                IF (NOUNIT) X(J) = X(J)*A(J,J)
              END IF
            END DO
          ELSE
            JX = KX
            DO J = 1,N
              IF (X(JX)/=ZERO) THEN
                TEMP = X(JX)
                IX = KX
                DO I = 1,J - 1
                  X(IX) = X(IX)+TEMP*A(I,J)
                  IX = IX+INCX
                END DO
                IF (NOUNIT) X(JX) = X(JX)*A(J,J)
              END IF
              JX = JX+INCX
            END DO
          END IF
        ELSE
          IF (INCX==1) THEN
            DO J = N,1,-1
              IF (X(J)/=ZERO) THEN
                TEMP = X(J)
                DO I = N,J+1,-1
                  X(I) = X(I)+TEMP*A(I,J)
                END DO
                IF (NOUNIT) X(J) = X(J)*A(J,J)
              END IF
            END DO
          ELSE
            KX = KX+(N-1)*INCX
            JX = KX
            DO J = N,1,-1
              IF (X(JX)/=ZERO) THEN
                TEMP = X(JX)
                IX = KX
                DO I = N,J+1,-1
                  X(IX) = X(IX)+TEMP*A(I,J)
                  IX = IX - INCX
                END DO
                IF (NOUNIT) X(JX) = X(JX)*A(J,J)
              END IF
              JX = JX - INCX
            END DO
          END IF
        END IF
      ELSE

!        Form  x := A'*x.

        IF (LSAME(UPLO,'U')) THEN
          IF (INCX==1) THEN
            DO J = N,1,-1
              TEMP = X(J)
              IF (NOUNIT) TEMP = TEMP*A(J,J)
              DO I = J - 1,1,-1
                TEMP = TEMP+A(I,J)*X(I)
              END DO
              X(J) = TEMP
            END DO
          ELSE
            JX = KX+(N-1)*INCX
            DO J = N,1,-1
              TEMP = X(JX)
              IX = JX
              IF (NOUNIT) TEMP = TEMP*A(J,J)
              DO I = J - 1,1,-1
                IX = IX - INCX
                TEMP = TEMP+A(I,J)*X(IX)
              END DO
              X(JX) = TEMP
              JX = JX - INCX
            END DO
          END IF
        ELSE
          IF (INCX==1) THEN
            DO J = 1,N
              TEMP = X(J)
              IF (NOUNIT) TEMP = TEMP*A(J,J)
              DO I = J+1,N
                TEMP = TEMP+A(I,J)*X(I)
              END DO
              X(J) = TEMP
            END DO
          ELSE
            JX = KX
            DO J = 1,N
              TEMP = X(JX)
              IX = JX
              IF (NOUNIT) TEMP = TEMP*A(J,J)
              DO I = J+1,N
                IX = IX+INCX
                TEMP = TEMP+A(I,J)*X(IX)
              END DO
              X(JX) = TEMP
              JX = JX+INCX
            END DO
          END IF
        END IF
      END IF

      RETURN

!     End of STRMV .

    END SUBROUTINE STRMV
    SUBROUTINE XERBLA(SRNAME,INFO)

!  -- LAPACK auxiliary routine (preliminary version) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

!     .. Scalar Arguments ..
      CHARACTER (6) :: SRNAME
      INTEGER :: INFO
!     ..

!  Purpose
!  =======

!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.

!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.

!  Arguments
!  =========

!  SRNAME  (input) CHARACTER*6
!          The name of the routine which called XERBLA.

!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.


      WRITE (*,FMT=9999) SRNAME,INFO

      STOP

9999  FORMAT (' ** On entry to ',A6,' parameter number ',I2,' had ', &
        'an illegal value')

!     End of XERBLA

    END SUBROUTINE XERBLA
