    SUBROUTINE ERRMSG(CSTATS,CSUBRT,CMSGNM)
      IMPLICIT NONE
! Arguments:
      CHARACTER :: CMSGNM*(*), CSTATS*(*), CSUBRT*(*)

! Local variables:
      INTEGER :: IOERR

!***********************************************************************
! Given:
!       CSTATS = 'WARNING' or 'FATAL'
!       CSUBRT = name of subroutine
!       CMESGN = message
! Prints a warning message in a standardized way,
! and STOPs if CSTATS='FATAL'

! Copyright (C) 1993, B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.

! History:
! 96.11.14 (PJF) Remove "getset" and hardwire "ioerr"
! 04.05.23 (BTD) cleanup
! end history
!***********************************************************************
      DATA IOERR/6/

      IF (CSTATS=='FATAL') THEN
        WRITE (IOERR,FMT=9000) CSUBRT
        WRITE (IOERR,FMT=9010) CMSGNM
        WRITE (IOERR,FMT=9020)
        STOP
      ELSE IF (CSTATS=='WARNING') THEN
        WRITE (IOERR,FMT=9030) CSUBRT
        WRITE (IOERR,FMT=9010) CMSGNM
      END IF
      RETURN
9000  FORMAT (/' >>>>> FATAL ERROR IN PROCEDURE: ',A)
9010  FORMAT (1X,A)
9020  FORMAT (' >>>>> EXECUTION ABORTED ')
9030  FORMAT (/' >>>>> WARNING IN PROCEDURE: ',A)
    END SUBROUTINE ERRMSG
