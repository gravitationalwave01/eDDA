    SUBROUTINE DIVIDE(CDIVID,X1,X2,NX,MXARY,ARY)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE
! Arguments:
      REAL (WP) :: X1, X2
      INTEGER :: MXARY, NX
      CHARACTER :: CDIVID*3
      REAL (WP) :: ARY(MXARY)

! Local variables:
      INTEGER :: I
      REAL (WP) :: DELTA

! External Subroutines:
      EXTERNAL ERRMSG, WRIMSG

! Intrinsic Functions:
      INTRINSIC REAL
!***********************************************************************
! Given:
!       CDIVID='LIN','INV', or 'LOG'
!       X1 = lower limit to interval
!       X2 = upper limit to interval
!       NX = number of elements
!       MXARY=dimensioning information for array ARY

! Returns:
!       ARY(1-NX)=vector of points with
!                 ARY(1)=X1
!                 ARY(NX)=X2
!                 ARY(J) values spaced either
!                        linearly (if CDIVID.EQ.'LIN')
!                        linearly in 1/X (if CDIVID.EQ.'INV')
!                        logarithmically (if CDIVID.EQ.'LOG')

! Copyright (C) 1993, B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
      IF (NX>MXARY) THEN
        CALL ERRMSG('FATAL','DIVIDE','  NX .GT. MXARY ')
      END IF
!      IF(X1.GT.X2)THEN
!         CALL ERRMSG('FATAL','DIVIDE','  X1 .GT. X2 ')
!      ENDIF
      IF (NX<1) THEN
        CALL ERRMSG('FATAL','DIVIDE','  NX .LT. 1 ')
      END IF
      IF (NX==1) THEN
        CALL WRIMSG('DIVIDE',' Only one element initialized ')
        ARY(1) = X1
      ELSE
        IF (CDIVID=='LIN') THEN
          DELTA = (X2-X1)/REAL(NX-1,KIND=WP)
          DO I = 1, NX
            ARY(I) = X1 + DELTA*REAL(I-1,KIND=WP)
          END DO
        ELSE IF (CDIVID=='INV') THEN
          DELTA = (1._WP/X2-1._WP/X1)/REAL(NX-1,KIND=WP)
          DO I = 1, NX
            ARY(I) = 1._WP/(1._WP/X1+DELTA*REAL(I-1,KIND=WP))
          END DO
        ELSE IF (CDIVID=='LOG') THEN
          DELTA = LOG(X2/X1)/REAL(NX-1,KIND=WP)
          DO I = 1, NX
            ARY(I) = EXP(LOG(X1)+DELTA*REAL(I-1,KIND=WP))
          END DO
        END IF
      END IF
      RETURN
    END SUBROUTINE DIVIDE
