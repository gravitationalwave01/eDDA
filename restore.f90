    SUBROUTINE RESTORE(CXV,CXUNOC,IOCC,MXN3,MXNAT,NAT,NAT0)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE
!***********************************************************************
! Given:
!      CXV(1-3*NAT0) defined for NAT0 occupied lattice sites
!      IOCC(1-NAT) = 0 or 1 depending on whether lattice site
!                      is vacant or occupied
!      CXUNOC      = value to use for unoccupied sites
!      MXN3,MXNAT  = dimensioning information
!      NAT         = number of lattice sites
!      NAT0        = number of occupied lattice sites
! Returns
!      CXV(1-3*NAT) defined for NAT lattice sites, occupied and
!                   unoccupied, with CXV(J)=CXUNOC for unoccupied sites
! Arguments:
      INTEGER :: MXN3, MXNAT, NAT, NAT0
      INTEGER*2 :: IOCC(MXNAT)
      COMPLEX (WP) :: CXUNOC
      COMPLEX (WP) :: CXV(MXN3)
! Local variables:
      INTEGER :: J, JUN, M

! B.T.Draine, Princeton Univ. Observatory
! History:
! 98.01.01 (BTD): Created for use in DDSCAT.6
! 06.07.07 (BTD): corrected typo in comments
! end history

! Copyright (C) 1998,2006 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
      DO M = 2, 0, -1
        JUN = NAT - NAT0
        DO J = NAT, 1, -1
          IF (IOCC(J)==0) THEN
            JUN = JUN - 1
            CXV(J+M*NAT) = CXUNOC
          ELSE
            CXV(J+M*NAT) = CXV(M*NAT0+J-JUN)
          END IF
        END DO
      END DO
      RETURN
    END SUBROUTINE RESTORE
