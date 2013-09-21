    SUBROUTINE PETR(X,B,WRK,LDA,IPAR,SPAR,MATVEC,CMATVEC)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE
! solves Ax=b

! Input:
! B(LOCLEN) - RHS (doesn't change on output)
! WRK(LDA, LOCLEN) - scratch array
! IPAR - defines integer paramteres (see documentation of PIM)
! SPAR - define real parameters
! LDA - leading dimension of arrays
! MATVEC -- name of procedure for computing y=Ax
! CMATVEC -- name of procedure for computing y= conj(A') x

! Output:
! X(LOCLEN)

! Reference:
!    Petravic,M and Kuo-Petravic,G. 1979, JCP, 32,263
! History:
!    Late 80's  (Coded by Draine)
!    Changes to include "cprod" (PJF+BTD)
!    95.06.01 (PJF) re-written to conform to PIM
!    95.08.11 (BTD) minor editing (comments, etc.)
!    08.05.12 (BTD) change CDUMMY -> CDUMMY(1) following suggestion
!                   by Art Lazanoff, NASA Ames
! end history
!=======================================================================
! Parameters:
      INTEGER :: IACE,IAXI,IGI,IPI,IQI,IR
      PARAMETER (IACE=1,IGI=2,IPI=3,IQI=4,IAXI=5,IR=6)
! Arguments:
      INTEGER :: LDA
      INTEGER :: IPAR(*)
      REAL(WP) :: SPAR(*)
      COMPLEX(WP) :: B(*),WRK(LDA,*),X(*)
      EXTERNAL CMATVEC,MATVEC

! Local variables
      COMPLEX(WP) :: CDUMMY(1)
      INTEGER :: IDUMMY,ITNO,L,LOCLEN,MAXIT,STATUS
      REAL(WP) :: ALPHAI,BETAI1,BNORM,EPSILON,EXITNORM,GIGI, &
                  GI1GI1,QIQI,RHSSTOP
      REAL(WP) :: SCSETRHSSTOP
      EXTERNAL DUMMY,CDUMMY_1,PSCNRM2

!***********************************************************************
      LOCLEN=IPAR(4)
      MAXIT=IPAR(10)
      EPSILON=SPAR(1)

! Compute RHSSTOP=EPSILON*|B|

      RHSSTOP=SCSETRHSSTOP(B,CDUMMY,EPSILON,IPAR,CDUMMY_1,PSCNRM2)

! Compute conjg(A')*B

      CALL CMATVEC(B,WRK(1,IACE),IDUMMY)
      BNORM=0._WP
      DO L=1,LOCLEN
         BNORM=BNORM+REAL(B(L)*CONJG(B(L)))
         WRK(L,IGI)=WRK(L,IACE)
         WRK(L,IPI)=WRK(L,IGI)
      ENDDO

! Compute |QI>=A|PI>

      CALL MATVEC(WRK(1,IPI),WRK(1,IQI),IDUMMY)
      QIQI=0._WP
      GIGI=0._WP
      DO L=1,LOCLEN
         GIGI=GIGI+REAL(CONJG(WRK(L,IGI))*WRK(L,IGI))
         QIQI=QIQI+REAL(CONJG(WRK(L,IQI))*WRK(L,IQI))
      ENDDO
      ALPHAI=GIGI/QIQI
      DO L=1,LOCLEN
         X(L)=X(L)+ALPHAI*WRK(L,IPI)
      ENDDO

! compute |AX1>:

      CALL MATVEC(X,WRK(1,IAXI),IDUMMY)

      DO ITNO=2,MAXIT
         IPAR(11)=ITNO

! Transfer <GI|GI> -> <GI-1|GI-1>:

         GI1GI1=GIGI

! compute |GI>=AC|E>-AC|AXI>:

         CALL CMATVEC(WRK(1,IAXI),WRK(1,IGI),IDUMMY)

! Compute GIGI=<GI|GI>:
         GIGI=0._WP
         DO L=1,LOCLEN
            WRK(L,IGI)=WRK(L,IACE)-WRK(L,IGI)
         ENDDO
         DO L=1,LOCLEN
            GIGI=GIGI+REAL(CONJG(WRK(L,IGI))*WRK(L,IGI))
         ENDDO

! Compute BETAI-1=<GI|GI>/<GI-1|GI-1>:

         BETAI1=GIGI/GI1GI1

! Compute |PI>=|GI>+BETAI-1|PI-1>:

         DO L=1,LOCLEN
            WRK(L,IPI)=WRK(L,IGI)+BETAI1*WRK(L,IPI)
         ENDDO

! Compute |QI>=A|PI>:

         CALL MATVEC(WRK(1,IPI),WRK(1,IQI),IDUMMY)

! Compute <QI|QI>:

         QIQI=0._WP
         DO L=1,LOCLEN
            QIQI=QIQI+REAL(CONJG(WRK(L,IQI))*WRK(L,IQI))
         ENDDO

! Compute ALPHAI=<GI|GI>/<QI|QI>:

         ALPHAI=GIGI/QIQI

! Compute |XI+1>=|XI>+ALPHAI*|PI>:

         DO L=1,LOCLEN
            X(L)=X(L)+ALPHAI*WRK(L,IPI)
         ENDDO

! Except every 10TH iteration, compute |AXI+1>=|AXI>+ALPHAI*|QI>
! Draine: warning this part perhaps not checked

         IF(ITNO/=10*(ITNO/10))THEN
            DO L=1,LOCLEN
               WRK(L,IAXI)=WRK(L,IAXI)+ALPHAI*WRK(L,IQI)
            ENDDO
         ELSE
            CALL MATVEC(X,WRK(1,IAXI),IDUMMY)
         ENDIF

! Compute residual vector |RI>=A|XI>-|B>

         DO L=1,LOCLEN
            WRK(L,IR)=WRK(L,IAXI)-B(L)
         ENDDO

! Call STOPCRIT to check whether to terminate iteration.
! Criterion: terminate if <RI|RI> < TOL**2 * <B|B>
! Return with STATUS=0 if this condition is satisfied.

         CALL STOPCRIT(B,WRK(1,IR),CDUMMY,CDUMMY,CDUMMY,WRK,RHSSTOP,IDUMMY, &
                       EXITNORM,STATUS,IPAR,DUMMY,DUMMY,DUMMY,DUMMY,PSCNRM2)

! Call PROGRESS to report "error" = sqrt(<RI|RI>)/sqrt(<B|B>)

         CALL PROGRESS(LOCLEN,ITNO,EXITNORM,CDUMMY,CDUMMY,CDUMMY)

         IF(STATUS==0)GOTO 500
      ENDDO
500   CONTINUE
      IPAR(12)=0
      RETURN
    END SUBROUTINE PETR
