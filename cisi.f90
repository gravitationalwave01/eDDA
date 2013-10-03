    SUBROUTINE CISI(X,CI,SI)
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE

! Given:
!     X = real argument
! Returns

!     CI = gamma + ln(x) + \int_0^x [(cos(t)-1)/t] dt = "cosine integral"

!     SI = \int_0^x [sin(t)/t] dt = "sine integral"
! 
! Code adapted from Numerical Recipes in FORTRAN (2e) 
! by Press, Teukolsky, Vetterling, and Flannery (1994)

! Arguments:

      REAL(WP) :: CI,SI,X

! Local variables:

      INTEGER :: MAXIT
      REAL(WP) :: EPMIN,EPS,EULER,FPMIN,PIBY2,TMIN
      PARAMETER(EPS=6.E-8,EULER=0.57721566_WP,MAXIT=100,  &
                PIBY2=1.5707963_WP,FPMIN=1.E-30,TMIN=2._WP)
      INTEGER :: I,K
      REAL(WP) :: A,ERR,FACT,SGN,SU,SUMC,SUMS,T,TERM,ABSC
      COMPLEX(WP) :: B,C,D,DEL,H
      LOGICAL ODD
      ABSC(H)=ABS(REAL(H))+ABS(AIMAG(H))
      T=ABS(X)
      IF(T.EQ.0.)THEN
	 SI=0.
	 CI=-1._WP/FPMIN
	 RETURN
      ENDIF
      IF(T.GT.TMIN)THEN
	 B=CMPLX(1._WP,T)
	 C=1._WP/FPMIN
	 D=1._WP/B
	 H=D
	 DO I=2,MAXIT
            A=-(I-1)**2
            B=B+2.
            D=1._WP/(A*D+B)
            C=B+A/C
            DEL=C*D
            H=H*DEL
            IF(ABSC(DEL-1._WP).LT.EPS)GOTO 1
         ENDDO
         CALL ERRMSG('FATAL','CISI',              &
                     'Fatal error: failed in cisi')
1	 CONTINUE
	 H=CMPLX(COS(T),-SIN(T))*H
	 CI=-REAL(H)
	 SI=PIBY2+AIMAG(H)
      ELSE
	 IF(T.LT.SQRT(FPMIN))THEN
	    SUMC=0.
            SUMS=T
	 ELSE
	    SU=0.
	    SUMS=0.
	    SUMC=0.
	    SGN=1._WP
	    FACT=1._WP
	    ODD=.TRUE.
	    DO K=1,MAXIT
               FACT=FACT*T/K
               TERM=FACT/K
               SU=SU+SGN*TERM
               ERR=TERM/ABS(SU)
               IF(ODD)THEN
                  SGN=-SGN
                  SUMS=SU
                  SU=SUMC
               ELSE
                  SUMC=SU
                  SU=SUMS
               ENDIF
               IF(ERR.LT.EPS)GOTO 2
               ODD=.NOT.ODD
            ENDDO
            CALL ERRMSG('FATAL','CISI',          &
                        'maxits exceeded in CISI')
	 ENDIF
2	 SI=SUMS
	 CI=SUMC+LOG(T)+EULER
      ENDIF
      IF(X.LT.0._WP)SI=-SI
      RETURN
    END SUBROUTINE CISI

