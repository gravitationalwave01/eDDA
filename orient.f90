    SUBROUTINE ORIENT(BETAMI,BETAMX,THETMI,THETMX,PHIMIN,PHIMAX,MXBETA,MXTHET, &
        MXPHI,NBETA,NTHETA,NPHI,BETA,THETA,PHI,WGTA,WGTB)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE
! Arguments:
      INTEGER :: MXBETA, MXTHET, MXPHI, NBETA, NTHETA, NPHI
      REAL (WP) :: BETAMI, BETAMX, THETMI, THETMX, PHIMIN, PHIMAX
      REAL (WP) :: BETA(MXBETA), THETA(MXTHET), PHI(MXPHI), &
        WGTA(MXTHET,MXPHI), WGTB(MXBETA)
! Local variables:
      INTEGER :: I, J
      REAL (WP) :: DELTA, WG
!***********************************************************************
! Given:
!        BETAMI=minimum value of beta (radians)
!        BETAMX=maximum value of beta (radians)
!        THETMI=minimum value of theta (radians)
!        THETMX=maximum value of theta (radians)
!        PHIMIN=minimum value of phi (radians)
!        PHIMAX=maximum value of phi (radians)
!        MXBETA,MXTHET,MXPHI=parameters for dimensioning of arrays
!                             BETA,THETA,PHI
!        NBETA=desired number of values of beta
!        NTHETA=desired number of values of theta
!        NPHI=desired number of values of PHI
! Returns:
!        BETA(1-NBETA)=beta values (radians)
!        THETA(1-NTHETA)=theta values (radians)
!        PHI(1-NPHI)=phi values (radians)
!        WGTA(1-NTHETA,1-NPHI)=weighting of each orientation of target
!                              axis a1 in Lab Frame
!                              (sum of weights = 1)
!        WGTB(1-NBETA)=weighting of each rotation of target around a1
!                              (sum of weights = 1)
!            Note: it is assumed that target orientation weight function
!                  can be factored into WGTA*WGTB -- i.e., that rotation
!                  around a1 are decoupled from orientation of a1.

! Purpose: to generate list of desired target orientations

! Definitions: beta=angle of rotation of target around target axis a1
!                   (beta=0 is defined to be such that a2 lies in a1-x
!                    plane, with a2 in direction of increasing theta)
!              theta=angle between target axis a1 and x axis
!                   (0.le.theta.le.pi)
!              phi=angle of rotation of axis a1 around axis
!              (phi=0 is defined to be such that a1 lies in xy plane)
! Present version assumes:
!        beta to be uniformly distributed between BETAMI and BETAMX
!        cos(theta) to be uniformly distributed between cos(THETMI) and
!                   cos(THETMX)
!        phi to be uniformly distributed between PHIMIN and PHIMAX

! Values assigned to BETA and PHI are midpoints of uniform intervals in
! beta and phi
! Note that if NBETA=1, BETA is set to midpoint of range in phi
! Likewise, if NPHI=1, PHI is set to midpoint of range in phi

! If NTHETA=1: specify a single value of THETA
!        cos(theta)=0.5*(cos(thetami)+cos(thetamx))

! If NTHETA>1:
!    If NTHETA is even:
!       divide range of cos(theta) into NTHETA equal intervals
!       set THETA to midpoints of these intervals
!       give intervals equal weights
!    If NTHETA is odd:
!       divide range of cos(theta) into NTHETA-1 equal intervals
!       set THETA to endpoints of these intervals; first value of
!       THETA is THETAMI and last value of THETA is THETAMX

! Note: values chosen for beta, phi are at midpoints of
!        uniform intervals

! B.T.Draine, Princeton Univ. Obs., 89.11.20
! History:
! 92.04.01 (BTD) Modify to allow Simpson's rule or trapezoidal rule
!                for integration over cos(theta).
! 99.03.01 (BTD) Modified to properly assign THETA values when NPHI>1
!                Code up to this time was not consistent with
!                description in UserGuide.  This problem was brought
!                to my attention by Miroslav Kocifaj
!                (Miroslav.Kocifaj@swh.sk).
! 99.03.05 (BTD) Corrected errors in assignment of THETA values.
!                Thanks to Miroslav Kocifaj for detecting these.
! 99.04.30 (BTD) Corrected errors in assignment of BETA and
!                PHI values when BETAMI.ne.0 or PHIMIN.ne.0
!                These errors (inadvertent omission of BETAMI and
!                PHIMIN from expressions assigning BETA(J) and PHI(J))
!                were introduced in 99.03.05 revision.
!                Thanks to Henriette Lemke for calling attention to
!                these errors.
! Copyright (C) 1993,1999 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

! Assign values to BETA:

      DELTA = (BETAMX-BETAMI)/REAL(NBETA,KIND=WP)
      DO J = 1, NBETA
        BETA(J) = BETAMI + DELTA*(REAL(J,KIND=WP)-0.5_WP)
      END DO

! Assign values to THETA:

      IF (NTHETA==2*(NTHETA/2)) THEN

! THETA is even:

        DELTA = (COS(THETMX)-COS(THETMI))/REAL(NTHETA,KIND=WP)
        DO J = 1, NTHETA
          THETA(J) = ACOS(COS(THETMI)+DELTA*(REAL(J,KIND=WP)-0.5_WP))
        END DO
      ELSE

! THETA is odd:

        IF (NTHETA==1) THEN
          THETA(1) = ACOS(0.5_WP*(COS(THETMI)+COS(THETMX)))
        ELSE
          DELTA = (COS(THETMX)-COS(THETMI))/REAL(NTHETA-1,KIND=WP)
          THETA(1) = THETMI
          THETA(NTHETA) = THETMX
          DO J = 2, NTHETA - 1
            THETA(J) = ACOS(COS(THETMI)+DELTA*REAL(J-1,KIND=WP))
          END DO
        END IF
      END IF

! Assign values to PHI:

      DELTA = (PHIMAX-PHIMIN)/REAL(NPHI,KIND=WP)
      DO J = 1, NPHI
        PHI(J) = PHIMIN + DELTA*(REAL(J,KIND=WP)-0.5_WP)
      END DO

!** Specify weight factors WGTA, WGTB
!   (Note: weight function WGTA = 4*pi*P/(NTHETA*NPHI), where
!    P=(probability/solid angle) of orientation in direction THETA,PHI .
!    The orientational averaging program automatically samples uniformly
!    in cos(theta) to allow for d(solid angle)=sin(theta)d(theta)d(phi).

!   Present version assumes random orientations.

!   When NTHETA is >1 and even, we use trapezoidal integration
!   When NTHETA is >1 and odd, we use Simpson's rule integration.

      IF (NTHETA==1 .OR. NTHETA==2*(NTHETA/2)) THEN

! either 1 or even number of theta values: midpoints of intervals

        WG = 1._WP/REAL(NTHETA*NPHI,KIND=WP)
        DO J = 1, NPHI
          DO I = 1, NTHETA
            WGTA(I,J) = WG
          END DO
        END DO
      ELSE

! odd number >1 of theta values: use Simpson's rule weighting

        WG = 1._WP/REAL(3*(NTHETA-1)*NPHI,KIND=WP)
        DO J = 1, NPHI
          WGTA(1,J) = WG
          WGTA(NTHETA,J) = WG
        END DO
        WG = 4._WP/REAL(3*(NTHETA-1)*NPHI,KIND=WP)
        DO I = 2, NTHETA, 2
          DO J = 1, NPHI
            WGTA(I,J) = WG
          END DO
        END DO
        IF (NTHETA>=5) THEN
          WG = 2._WP/REAL(3*(NTHETA-1)*NPHI,KIND=WP)
          DO I = 3, NTHETA - 1, 2
            DO J = 1, NPHI
              WGTA(I,J) = WG
            END DO
          END DO
        END IF
      END IF

! Assign weights for rotation angle beta:

      WG = 1._WP/REAL(NBETA,KIND=WP)
      DO I = 1, NBETA
        WGTB(I) = WG
      END DO

      RETURN
    END SUBROUTINE ORIENT
