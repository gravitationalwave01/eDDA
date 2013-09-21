    SUBROUTINE ROTATE(A1,A2,AK1,CXE01,CXE02,BETA,THETA,PHI,EN0R,CXE01R,CXE02R, &
        MXSCA,NSCAT,ENSC,EM1,EM2,AKSR,EM1R,EM2R)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

! Arguments:

      INTEGER :: MXSCA, NSCAT
      REAL (WP) :: AK1, BETA, THETA, PHI
      REAL (WP) :: A1(3), A2(3), AKSR(3,MXSCA), EM1(3,MXSCA), EM1R(3,MXSCA), &
        EM2(3,MXSCA), EM2R(3,MXSCA), EN0R(3), ENSC(3,MXSCA)
      COMPLEX (WP) :: CXE01(3), CXE02(3), CXE01R(3), CXE02R(3)

! Local variables

      INTEGER :: I
      REAL (WP) :: BETA0, PHI0, PI, SINTHE, TERM, THETA0
      REAL (WP) :: RM1(3,3), RM2(3,3), RM3(3,3), VEC(3)
!***********************************************************************
! Given:
!       A1(1-3)=axis 1 of (original, unrotated) target
!       A2(1-3)=axis 2 of (original, unrotated) target
!       AK1=magnitude of k vector
!       CXE01(1-3)=original incident polarization vector 1
!       CXE02(1-3)=original incident polarization vector 2
!       BETA=angle (radians) through which target is to be rotated
!            around target axis A1
!       THETA=angle (radians)which target axis A1 is to make with
!            Lab x-axis
!       PHI=angle (radians) between Lab x,A1 plane and Lab x,y plane
!       ENSC(1-3,1-NSCAT)=scattering directions in Lab Frame
!       EM1(1-3,1-NSCAT)=scattering polarization vector 1 in Lab Frame
!       EM2(1-3,1-NSCAT)=scattering polarization vector 2 in Lab Frame

! Returns:
!       EN0R(1-3)=incident propagation vector in Target Frame
!       CXE01R(1-3)=incident polariz. vector 1 in Target Frame
!       CXE02R(1-3)=incident polariz. vector 2 in Target Frame
!       AKSR(1-3,1-NSCAT)=scattering k vectors in Target Frame
!       EM1R(1-3,1-NSCAT)=scattering pol. vector 1 in Target Frame
!       EM2R(1-3,1-NSCAT)=scattering pol. vector 2 in Target Frame

! Purpose:
!       In the Lab Frame, we hold the propagation direction (x-axis)
!       and incident polarization vectors (CXE01 and CXE02) fixed,
!       and rotate target to an orientation specified by the three
!       angles BETA,THETA,PHI.
!       However, computations in the main program DDSCAT are carried out
!       in the Target Frame, in which the target is held fixed and the
!       propagation vectors and polarization vectors are rotated.
!       Hence, this routine computes the incident propagation vector,
!       scattering propagation vectors, and associated polarization
!       vectors in the Target Frame.

! B. T. Draine, Princeton Univ. Obs., 89.11.20
! History:
! 91.04.22 (BTD) Eliminated vector A1R and removed superfluous line
!                computing A1R
! 96.11.03 (BTD) Corrected comments
! 03.06.06 (BTD) Corrected calculation of BETA0,PHI0 -- it appears that
!                these were not correct when a_1z < 0
!                                     and/or a_3x < 0
! 05.03.29 (BTD) Add a few lines of code to guard against possibility
!                that numerical roundoff will cause argument of
!                ACOS to have absolute value > 1.
!                Declare new veriable TERM
! end history

! Copyright (C) 1996,2003,2005 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

      PI = 4._WP*ATAN(1._WP)

! Determine initial orientation THETA0,PHI0,BETA0 of target:
! theory:

!     theta = arccos(a_1x)

! if sin(theta) .ne. 0 , then

!     cos(phi)=a_1y/sin(theta)      sin(phi)=a_1z/sin(theta)
!     if(a_1z/sin(theta) > 0) then phi = arccos(a_1y/sin(theta))
!     if(a_1z/sin(theta) < 0) then phi = 2*pi-arccos(a_1y/sin(theta))

!     cos(beta)=-a_2x/sin(theta)    sin(beta)=a_3x/sin(theta)
!     if(a_3x/sin(theta) > 0) then beta = arccos(-a_2x/sin(theta))
!     if(a_3x/sin(theta) < 0) then beta = 2*pi-arccos(-a_2x/sin(theta))

! note: a_3x = a_1y*a_2z - a_1z*a_2y

! if sin(theta) = 0 , then

!     phi = 0

!     cos(beta)=a_2y                sin(beta)=a_2z
!     if(a_2z) > 0)      then      beta = arccos(a_2y)
!     if(a_2z) < 0)      then      beta = 2*pi-arccos(a_2y)

      THETA0 = ACOS(A1(1))
      SINTHE = SIN(THETA0)
      IF (SINTHE/=0._WP) THEN

!-------------------------------
! guard against roundoff errors:
        TERM = A1(2)/SINTHE
        IF (TERM>1._WP) TERM = 1._WP
        IF (TERM<-1._WP) TERM = -1._WP
!-------------------------------

        IF (A1(3)>=0._WP) THEN
          PHI0 = ACOS(TERM)
        ELSE
          PHI0 = 2._WP*PI - ACOS(TERM)
        END IF

!------------------------------
! guard against roundoff errors:
        TERM = -A2(1)/SINTHE
        IF (TERM>1._WP) TERM = 1._WP
        IF (TERM<-1._WP) TERM = -1
!-------------------------------

        IF (A1(2)*A2(3)-A1(3)*A2(2)>=0._WP) THEN
          BETA0 = ACOS(TERM)
        ELSE
          BETA0 = 2._WP*PI - ACOS(TERM)
        END IF
      ELSE
        PHI0 = 0._WP
        IF (A2(3)>=0._WP) THEN
          BETA0 = ACOS(A2(2))
        ELSE
          BETA0 = 2._WP*PI - ACOS(A2(2))
        END IF
      END IF

!**** First rotate the target through angle PHI-PHI0 around x axis
!    (i.e., rotate the lab through angle PHI0-PHI around x axis)
!** obtain rotation matrix:

      VEC(1) = 1._WP
      VEC(2) = 0._WP
      VEC(3) = 0._WP

      CALL ROT2(VEC,PHI0-PHI,RM1)

!**** Now rotate the target through angle THETA-THETA0 around axis
!     VEC = xaxis cross a1
!     (i.e., rotate the lab through angle THETA0-THETA around
!     xaxis cross a1 )
!** first compute rotation axis

      VEC(1) = 0._WP
      VEC(2) = -A1(3)
      VEC(3) = A1(2)

!*** For special case where a1=xaxis, we want rotation axis=a1 cross a2

      IF ((VEC(2)**2+VEC(3)**2)<1.E-10_WP) THEN
        VEC(2) = A1(3)*A2(1) - A1(1)*A2(3)
        VEC(3) = A1(1)*A2(2) - A1(2)*A2(1)
      END IF

      CALL ROT2(VEC,THETA0-THETA,RM2)
      CALL MULT3(RM2,RM1,RM3)

!**** Now rotate the target through angle BETA-BETA0 around a1
!     (i.e., rotate the lab through angle BETA0-BETA around a1)

      CALL ROT2(A1,BETA0-BETA,RM1)
      CALL MULT3(RM1,RM3,RM2)

!**** RM2 is now the full rotation matrix

!**** Now rotate all required vectors:
! Incident propagation vector:

      VEC(1) = 1._WP
      VEC(2) = 0._WP
      VEC(3) = 0._WP

      CALL PROD3(RM2,VEC,EN0R)

! Incident polarization vectors:

      CALL PROD3C(RM2,CXE01,CXE01R)
      CALL PROD3C(RM2,CXE02,CXE02R)

! Scattering vectors:

      CALL PROD3V(RM2,ENSC,AKSR,NSCAT,MXSCA)

      DO I = 1, NSCAT
        AKSR(1,I) = AK1*AKSR(1,I)
        AKSR(2,I) = AK1*AKSR(2,I)
        AKSR(3,I) = AK1*AKSR(3,I)
      END DO

! Scattering polarization vectors:

      CALL PROD3V(RM2,EM1,EM1R,NSCAT,MXSCA)
      CALL PROD3V(RM2,EM2,EM2R,NSCAT,MXSCA)

      RETURN
    END SUBROUTINE ROTATE
