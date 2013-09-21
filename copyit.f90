    SUBROUTINE COPYIT(CXA,CXB,NPTS)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE
!*** Arguments:
      INTEGER :: NPTS
      COMPLEX (WP) :: CXA(*), CXB(*)

!*** Local variables:
      INTEGER :: IPTS
!***********************************************************************
! Purpose: To copy a complex vector from one memory location to another.
!***********************************************************************
      DO IPTS = 1, NPTS
        CXB(IPTS) = CXA(IPTS)
      END DO
      RETURN
    END SUBROUTINE COPYIT
