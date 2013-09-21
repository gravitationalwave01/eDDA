    SUBROUTINE NAMID(MYID,CFLLOG)
      IMPLICIT NONE
!***********************************************************************
! Purpose: to generate unique file names for log file written
!          each MPI process

! Present version allows up to 1000 MPI processes (000-999)

! Given:
!         MYID=integer (0-999) MPI ID number

! Returns:
!          CFLLOG=name for output file containing running output
!                 from DDSCAT
!                 = ddscat.log_nnn
!                   where nnn = MYID

! B.T.Draine, Princeton Univ. Observatory, 2003.04.12
! History:
! 03.04.12 (BTD): Created using NAMER as an example
! end history

! Copyright (C) 2003, B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
! Arguments:

      INTEGER :: MYID
      CHARACTER :: CFLLOG*14

! Local variables:

      CHARACTER :: A1*11, ZLOG*14
      CHARACTER :: N(0:9)*1, V(1:3)*1, W(1:14)*1

! Equivalences:
      EQUIVALENCE (A1,W(1)), (V,W(12))
      EQUIVALENCE (W,ZLOG)

! DATA statements:

      DATA N/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
!***********************************************************************
!*** Filename CFLLOG will be of the form ddscat.log_NNN
!     where NNN = MYID

      A1 = 'ddscat.log_'

      V(1) = N(MYID/100)
      V(2) = N((MYID-100*(MYID/100))/10)
      V(3) = N(MYID-10*(MYID/10))

      CFLLOG = ZLOG

      RETURN

    END SUBROUTINE NAMID
