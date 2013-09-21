    SUBROUTINE CXFFTW(CX,MX,MY,MZ,ISIGN)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

! Arguments:

      INTEGER :: ISIGN, MX, MY, MZ
      COMPLEX (WP) :: CX(MX,MY,MZ)

! Local variables:

      CHARACTER :: CMSGNM*70

!***********************************************************************

! Purpose: This is a dummy routine to substitute for CXFFTW to allow
!          DDSCAT to be used on systems where the FFTW package has not
!          been installed.  If called, it will generate a fatal error
!          with an error message reporting DDSCAT has been compiled
!          with this dummy routine and option FFTW is unavailable.

! History:
! 00.07.05 (BTD): created
! 03.07.13 (BTD): changed option FFTWFJ to option FFTW21
! 04.03.31 (BTD): replaced WRITE(0,  with CALL WRIMSG(...
! end history
! Copyright (C) 2000, 2003
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.

!***********************************************************************

      WRITE (CMSGNM,FMT='(A,A)') &
        'FATAL ERROR: DDSCAT has been compiled with dummy version of', &
        ' CXFFTW'
      CALL WRIMSG('CXFFTW',CMSGNM)
      WRITE (CMSGNM,FMT='(A)') &
        ' *** option FFTW21 cannot be used with dummy routine'
      CALL WRIMSG('CXFFTW',CMSGNM)
      WRITE (CMSGNM,FMT='(A)') &
        ' *** to enable option FFTW21 it is necessary to:'
      CALL WRIMSG('CXFFTW',CMSGNM)
      WRITE (CMSGNM,FMT='(A)') &
        ' --- have fftw 2.1.x library installed on system'
      CALL WRIMSG('CXFFTW',CMSGNM)
      WRITE (CMSGNM,FMT='(A,A)') &
        ' --- in Makefile, variable LIBFFTW should include ', &
        'correct path to libsfftw.a'
      CALL WRIMSG('CXFFTW',CMSGNM)
      WRITE (CMSGNM,FMT='(A,A)') ' --- in Makefile, define cxfftw = cxfftw', &
        ' (rather than dummycxfftw.f)'
      CALL WRIMSG('CXFFTW',CMSGNM)
      WRITE (CMSGNM,FMT='(A)') ' --- then make ddscat'
      CALL WRIMSG('CXFFTW',CMSGNM)
      STOP
    END SUBROUTINE CXFFTW
