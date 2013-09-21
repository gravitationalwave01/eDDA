    SUBROUTINE ROT2(A,THETA,RM)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE
!***********************************************************************
! This subroutine finds the (3 by 3) rotation matrix for the
! given axis (vector) --- a and the rotation angle --- theta.
! On input:
! a(3)    --- axis of rotation  (vector, any length)
! theta   --- angle of rotation (in radians) around vector a
! On output:
! rm(3,3) --- rotation matrix such that new vector v2 is obtained
!             from the old vector v1 by the transformation:
!             v2 = rm \dot v1
! (If more rotations are needed: generate matrices rm with the
! help of this subroutine and (matrix) multiply them before
! performing v2 = rm \dot v1.)

! Reference:
! Leubner, C., 1977,  Coordinate-free  rotation operator,
! Am. J. Phys. 47, 727---729.

! History:
! Written by PJF, May 1989, for use in the rotation-averaged
! DDA calculations.

! Copyright (C) 1993, B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
! pre-calculate some stuff:
!     .. Scalar Arguments ..
      REAL (WP) :: THETA
!     ..
!     .. Array Arguments ..
      REAL (WP) :: A(3), RM(3,3)
!     ..
!     .. Local Scalars ..
      REAL (WP) :: ANORM, CT, ST
      INTEGER :: I, J
!     ..
!     .. Local Arrays ..
      REAL (WP) :: AHAT(3)
!     ..
!     .. External Subroutines ..
      EXTERNAL XNORM3
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC COS, SIN
!     ..
      CT = COS(THETA)
      ST = SIN(THETA)
      CALL XNORM3(A,ANORM)
      AHAT(1) = A(1)/ANORM
      AHAT(2) = A(2)/ANORM
      AHAT(3) = A(3)/ANORM
!  cos(\theta) {\bf 1}  term:
      RM(1,1) = CT
      RM(2,2) = CT
      RM(3,3) = CT
! Skew-term  a \times {\bf 1}:
      RM(1,2) = -ST*AHAT(3)
      RM(1,3) = ST*AHAT(2)
      RM(2,1) = ST*AHAT(3)
      RM(2,3) = -ST*AHAT(1)
      RM(3,1) = -ST*AHAT(2)
      RM(3,2) = ST*AHAT(1)
! aa-dyadic (outer-product):
      DO I = 1, 3
        DO J = 1, 3
          RM(I,J) = RM(I,J) + (1.E0_WP-CT)*AHAT(J)*AHAT(I)
        END DO
      END DO
      RETURN
    END SUBROUTINE ROT2

    SUBROUTINE XNORM3(A,XN)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

! Arguments

      REAL (WP) :: XN
      REAL (WP) :: A(3)

! Intrinsic Functions ..

      INTRINSIC SQRT

! auxilary for "rotation" routines

! Copyright (C) 1993, B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

      XN = SQRT(A(1)**2+A(2)**2+A(3)**2)
      RETURN
    END SUBROUTINE XNORM3

    SUBROUTINE MULT3(A,B,C)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

! Arguments:

      REAL (WP) :: A(3,3), B(3,3), C(3,3)

! Local variables:

      REAL (WP) :: S
      INTEGER :: I, J, K
!***********************************************************************
! Purpose:
! Given:  A(I,J)= 3x3 matrix
!         B(J,K)= 3x3 matrix
! Returns:C(I,K)= matrix product of A and B
! This is auxilary for "rotation" routines

! Copyright (C) 1993, B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!**********************************************************************

      DO J = 1, 3
        DO I = 1, 3
          S = 0.E0_WP
          DO K = 1, 3
            S = S + A(I,K)*B(K,J)
          END DO
          C(I,J) = S
        END DO
      END DO
      RETURN
    END SUBROUTINE MULT3

    SUBROUTINE PROD3(RM,VIN,VOUT)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

! Arguments;

      REAL (WP) :: RM(3,3), VIN(3), VOUT(3)

! Local variables:

      INTEGER :: I, J

! Given:  RM(I,J)=3x3 matrix
!         VIN(J)=3-vector
! Returns:VOUT(I)=RM*VIN

! Auxilary for "rotation" routines

! Copyright (C) 1993, B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

      DO I = 1, 3
        VOUT(I) = 0.E0_WP
        DO J = 1, 3
          VOUT(I) = VOUT(I) + RM(I,J)*VIN(J)
        END DO
      END DO
      RETURN
    END SUBROUTINE PROD3

    SUBROUTINE PROD3C(RM,CVIN,CVOUT)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

! Arguments:

      REAL (WP) :: RM(3,3)
      COMPLEX (WP) :: CVIN(3), CVOUT(3)

! Local variables:

      INTEGER :: I, J

! Given:
!       RM=3x3 rotation matrix
!       CVIN=complex 3-vector
! Computes:
!       CVOUT=RM \dot CVIN

! Copyright (C) 1993, B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
      DO I = 1, 3
        CVOUT(I) = (0.E0_WP,0.E0_WP)
        DO J = 1, 3
          CVOUT(I) = CVOUT(I) + RM(I,J)*CVIN(J)
        END DO
      END DO
      RETURN
    END SUBROUTINE PROD3C

    SUBROUTINE PROD3V(RM,VIN,VOUT,NV,NVMAX)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

! Arguments:

      INTEGER :: NV, NVMAX
      REAL (WP) :: RM(3,3), VIN(3,NVMAX), VOUT(3,NVMAX)

! Local variables:
      INTEGER :: I, J, N
! Given:
!       RM = 3x3 rotation matrix
!       VIN = NV 3-vectors
! Computes:
!       VOUT=RM \dot VIN for each of the NV 3-vectors

! Copyright (C) 1993,1996 B.T. Draine and P.J. Flatau
! History
! 96.08.09 (BTD) trivial modification to remove obsolescent construct.
! end history
! This code is covered by the GNU General Public License.
!***********************************************************************

      DO N = 1, NV
        DO I = 1, 3
          VOUT(I,N) = 0.E0_WP
        END DO
      END DO

      DO N = 1, NV
        DO J = 1, 3
          DO I = 1, 3
            VOUT(I,N) = VOUT(I,N) + RM(I,J)*VIN(J,N)
          END DO
        END DO
      END DO

      RETURN
    END SUBROUTINE PROD3V
