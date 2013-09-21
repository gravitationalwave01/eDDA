    SUBROUTINE CXFFT3_MKL(CX,MX,MY,MZ,ISIGN)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

! Arguments

      INTEGER :: ISIGN,MX,MY,MZ
      COMPLEX(WP) :: CX(MX*MY*MZ)

! Local variables

      CHARACTER :: CMSGNM*70

!=======================================================================

! Purpose: This is a dummy routine to substitute for CXFFT3_MKL to
!          allow DDSCAT to be used on systems where the Intel MKL
!          library is not available.  If called, it will generate
!          a fatal error with an explanatory error message.
! History
! 08.06.05 (BTD): created
! end history
! Copyright (C) 2008
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!=======================================================================

      WRITE(CMSGNM,FMT='(A,A)') &
        'FATAL ERROR: DDSCAT compiled with cxfft3_mkl_fake'
      CALL WRIMSG('CXFFT3_MKL',CMSGNM)
      WRITE(CMSGNM,FMT='(A)') &
        ' *** option FFTMKL cannot be used with dummy routine'
      CALL WRIMSG('CXFFT3_MKL',CMSGNM)
      WRITE(CMSGNM,FMT='(A)') &
        ' *** to enable option FFTMKL it is necessary to:'
      CALL WRIMSG('CXFFT3_MKL',CMSGNM)
      WRITE(CMSGNM,FMT='(A)') &
        '     have Intel Math Kernel Library (MKL) installed on system'
      CALL WRIMSG('CXFFT3_MKL',CMSGNM)
      WRITE(CMSGNM,FMT='(A)') &
        '     and edit Makefile to use cxfft3_mkl.f90 and mkl_dfti.f90'
      CALL WRIMSG('CXFFT3_MKL',CMSGNM)

      STOP
    END SUBROUTINE CXFFT3_MKL
