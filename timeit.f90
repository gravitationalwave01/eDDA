      SUBROUTINE TIMEIT(CMSGTM,DTIME)
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE
      CHARACTER :: CMSGTM*(*)
      REAL(WP) :: DTIME

!***********************************************************************
! Subroutine TIMEIT
!
! Given:
!       CMSGTM = string
!
! Returns:
!       If odd-numbered call (first, third, etc.): no input/output:
!       
!       If even-numbered call:
!         * Call WRIMSG to print out CMSGTM and elapsed cpu time since
!           previous call
!         * Return
!           DTIME = elapsed cputime (sec) on
!   
! This version of timeit uses the system call "etime" available under
! Linux, Solaris, and some other bsd-like Unix systems (e.g., Convex unix).
!
!***********************************************************************
! history
! 94.06.20 (PJF) add dtime to the formal parameters
! 95.06.19 (PJF) modified to suppress printing in CMSGTM='noprint'
! 04.04.05 (BTD) changed IF(CMSGTM(1:7).NE.. to IF(CMSGTM.NE..
! 07.07.31 (BTD) rewritten to use f95 intrinsic routine CPU_TIME
! 08.03.14 (BTD) add USE DDPRECISION,ONLY: WP
!                REAL(WP) DTIME
! end history
! Copyright (C) 1993,1994,1995,2004,2007 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!=======================================================================
! External system calls:

!     EXTERNAL CPU_TIME

! Local variables:

      CHARACTER :: CSTA*3,CMSGNM*70
      REAL :: T1,T2

! External subroutines:

      EXTERNAL WRIMSG

! Local variables to be saved:

      SAVE CSTA,T1

! Data statements:

      DATA CSTA/'ONE'/

      IF(CSTA.EQ.'ONE')THEN
         CSTA='TWO'
         CALL CPU_TIME(T1)
      ELSEIF(CSTA.EQ.'TWO')THEN
         CSTA='ONE'
         CALL CPU_TIME(T2)
         DTIME=T2-T1
         IF(CMSGTM.NE.'NOPRINT'.AND.CMSGTM.NE.'noprint')THEN
            WRITE(CMSGNM,FMT='(A,A)')' Timing results for: ',CMSGTM
            CALL WRIMSG('TIMEIT',CMSGNM)
            IF(DTIME.GT.1.E4)THEN
               WRITE(CMSGNM,FMT='(F9.0, A)')DTIME,' = CPU time (sec)'
            ELSE
               WRITE(CMSGNM,FMT='(F9.3, A)')DTIME,' = CPU time (sec)'
            ENDIF
            CALL WRIMSG('TIMEIT',CMSGNM)
         ENDIF
      ENDIF
      RETURN
      END
