    SUBROUTINE INTERP(X,Y,Z,XA,YA,ZA,IDVOUT,MXTAB,NTAB)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE
! Arguments:
      INTEGER :: IDVOUT, MXTAB, NTAB
      REAL (WP) :: X, Y, Z
      REAL (WP) :: XA(MXTAB), YA(MXTAB), ZA(MXTAB)
! Local variables:
      INTEGER :: I1, I2, I3, I4, INC
      REAL (WP) :: SGN, X1, X2, X3, X4, Y1, Y2, Y3, Y4
!***********************************************************************
! Given array of tabulated values XA(1-NTAB), YA(1-NTAB), ZA(1-NTAB),
! and given independent variable X, this routine interpolates for two
! dependent variables: Y and Z.
! Uses only parabolic interpolation for safety.
! When called with X outside range of tabulated values XA, returns
! extrapolated values Y and Z but issues warning statement.
! Note: XA values must be either monotonically increasing OR
!       monotonically decreasing.
! B.T.Draine, Princeton Univ. Observatory
! 89/12/19 (BTD): Revised.
! 91/05/02 (BTD): Added IDVOUT to argument list.
!                 Changed WRITE(0 -> WRITE(IDVOUT
! 92/03/24 (BTD): Added test IF(I2.GE.NTAB) in case INTERP is used for
!                 data sets with different lengths.

! Copyright (C) 1993, B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
      DATA I2/1/

      IF (I2>=NTAB) I2 = NTAB - 1
      INC = 0
      SGN = 1.E0_WP
!*** Check whether X is increasing or decreasing:
      IF (XA(1)>=XA(NTAB)) SGN = -1.E0_WP
!*** Check whether outside table limits
      IF ((SGN*(X-XA(1)))<0.E0_WP) THEN
        I2 = 2
        WRITE (IDVOUT,FMT=6990) X
        GO TO 4600
      END IF
      IF (SGN*(X-XA(NTAB))>0.E0_WP) THEN
        I2 = NTAB - 1
        WRITE (IDVOUT,FMT=6990) X
        GO TO 4600
      END IF
!*** X is within table limits.  Find I2
1000  IF (SGN*(XA(I2)-X)) 2000, 3000, 4000
2000  IF (INC) 2100, 2200, 2200
2100  IF (I2+2-NTAB) 4700, 4700, 2500
2200  INC = 1
      I2 = I2 + 1
      IF (I2+1<=NTAB) GO TO 1000
2500  I2 = NTAB - 1
      GO TO 4600
3000  Y = YA(I2)
      Z = ZA(I2)
      RETURN
4000  IF (INC) 4200, 4200, 4100
4100  I2 = I2 - 1
      IF (I2-2) 4500, 4700, 4700
4200  INC = -1
      I2 = I2 - 1
      IF (I2>=2) GO TO 1000
4500  I2 = 2
4600  I1 = I2 - 1
      I3 = I2 + 1
      X1 = XA(I1)
      X2 = XA(I2)
      X3 = XA(I3)
      Y1 = YA(I1)
      Y2 = YA(I2)
      Y3 = YA(I3)
      CALL PARAB3(X,Y,X1,X2,X3,Y1,Y2,Y3)
      Y1 = ZA(I1)
      Y2 = ZA(I2)
      Y3 = ZA(I3)
      CALL PARAB3(X,Z,X1,X2,X3,Y1,Y2,Y3)
      RETURN
4700  I1 = I2 - 1
      I3 = I2 + 1
      I4 = I2 + 2
      X1 = XA(I1)
      X2 = XA(I2)
      X3 = XA(I3)
      X4 = XA(I4)
      Y1 = YA(I1)
      Y2 = YA(I2)
      Y3 = YA(I3)
      Y4 = YA(I4)
      CALL PARAB4(X,Y,X1,X2,X3,X4,Y1,Y2,Y3,Y4)
      Y1 = ZA(I1)
      Y2 = ZA(I2)
      Y3 = ZA(I3)
      Y4 = ZA(I4)
      CALL PARAB4(X,Z,X1,X2,X3,X4,Y1,Y2,Y3,Y4)
      RETURN
6990  FORMAT ('Warning from INTERP: outside table limits for X=',1P,D12.5)
    END SUBROUTINE INTERP

    SUBROUTINE PARAB3(X,Y,X1,X2,X3,Y1,Y2,Y3)
      USE DDPRECISION, ONLY : WP
      REAL (WP) :: A, B, X, X1, X2, X3, Y, Y1, Y2, Y3
!***********************************************************************
! Subroutine PARAB3 does parabolic interpolation, with parabola
! constrained to fit (x1,y1),(x2,y2),(x3,y3) exactly.
! B.T.Draine, Institute for Advanced Study, March 1980.

! Copyright (C) 1993, B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
      A = (Y3-Y2-(X3-X2)*(Y1-Y2)/(X1-X2))/((X3-X2)*(X3-X1))
      B = (Y1-Y2)/(X1-X2) - A*(X1-X2)
      Y = (A*(X-X2)+B)*(X-X2) + Y2
      RETURN
    END SUBROUTINE PARAB3

    SUBROUTINE PARAB4(X,Y,X1,X2,X3,X4,Y1,Y2,Y3,Y4)
      USE DDPRECISION, ONLY : WP
      REAL (WP) :: A, B, X, X1, X2, X3, X4, Y, Y1, Y2, Y3, Y4
!***********************************************************************
! Subroutine PARAB4 is designed to do parabolic interpolation, with
! parabola constrained to match (x2,y2) and (x3,y3) exactly, and to
! minimize sum of squared deviations from (x1,y1) and (x4,y4).
! It is assumed that x1.lt.x2.le.x.lt.x3.lt.x4
! Fit parabola Y=A*(X-X2)**2+B*(X-X2)+Y2
! B.T.Draine, Institute for Advanced Study, March 1980.

! Copyright (C) 1993, B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!**********************************************************************
      A = ((X1-X2)*(Y3-Y2)/(X3-X2)+Y2-Y1)*(X1-X2)*(X1-X3) + &
        ((X4-X2)*(Y3-Y2)/(X3-X2)+Y2-Y4)*(X4-X2)*(X4-X3)
      A = -A/(((X1-X2)*(X1-X3))**2+((X4-X2)*(X4-X3))**2)
      B = (Y3-Y2)/(X3-X2) - A*(X3-X2)
      Y = (A*(X-X2)+B)*(X-X2) + Y2
      RETURN
    END SUBROUTINE PARAB4
