    SUBROUTINE WRIMSG(CSUBRT,CMSGNM)
      IMPLICIT NONE
! Standard procedure for writing messages
! History:
! 96.11.14 (PJF) Remove "getset" and hardwire" ioout"
! 96.11.20 (BTD) change IOOUT to IDVOUT

! Copyright (C) 1993,1996 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
! arguments
      CHARACTER :: CMSGNM*(*), CSUBRT*(*)

! local variables:

      INTEGER :: IDVOUT
      SAVE IDVOUT

! Note: on some systems (e.g., Solaris) IDVOUT=0 generates unbuffered
! output to "standard output".  On other systems IDVOUT=0 may not be
! valid; then set IDVOUT=6 to get "standard output", probably buffered.

      DATA IDVOUT/0/

      WRITE(IDVOUT,FMT=9000)CSUBRT,CMSGNM
9000  FORMAT(' >',A,' ',A)
      RETURN
    END SUBROUTINE WRIMSG
