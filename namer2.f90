    SUBROUTINE NAMER2(IWAV,IRAD,IDIR,CFLFML)
      IMPLICIT NONE
!***********************************************************************
! Purpose: to generate file names for specific (wavelength,
!                                               size,
!                                               incident direction)
! Present version allows up to 1000 wavelengths, (000-999)
!                              1000 sizes        (000-999)
!                              1000 directions   (000-999)
! Given:
!         IWAV=integer (0-999) identifying wavelength
!         IRAD=integer (0-999) identifying size
!         IDIR=integer (0-999) identifying direction
! Returns:
!         CFLFML =name for output file containing complex scattering
!                 amplitudes f_ml for selected directions
!
! B.T.Draine, Princeton Univ. Observatory, 2007
! History:
! 07.08.31 (BTD): Adapted from NAMER
! Copyright (C) 2007 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
! Arguments:

      INTEGER :: IDIR, IDIR0, IRAD, IRAD0, IWAV, IWAV0
      CHARACTER :: CFLFML*16

! Local variables:

      CHARACTER :: N(10)*1, V1(3)*1, V2(3)*1, V3(3)*1

! DATA statements:

      DATA N/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
      SAVE N

!***********************************************************************
!*** Filename CFLSCA will be of the form w001r01k001
!     (for IWAV=2,IRAD=2,IDIR=2)

!*** Set IWAV0=IWAV-1 to run over range 0-999
!    Set IRAD0=IRAD-1 to run over range 0-99
!    Set IDIR0=IDIR-1 to run over range 0-999

      IWAV0=IWAV - 1
      IRAD0=IRAD - 1
      IDIR0=IDIR - 1
      V1(1)=N(IWAV0/100+1)
      V1(2)=N((IWAV0-100*(IWAV0/100))/10+1)
      V1(3)=N(IWAV0-10*(IWAV0/10)+1)
      V2(1)=N(IRAD0/100+1)
      V2(2)=N((IRAD0-100*(IRAD0/100))/10+1)
      V2(3)=N(IRAD0-10*(IRAD0/10)+1)
      V3(1)=N(IDIR0/100+1)
      V3(2)=N((IDIR0-100*(IDIR0/100))/10+1)
      V3(3)=N(IDIR0-10*(IDIR0/10)+1)

      CFLFML='w'//V1(1)//V1(2)//V1(3)//'r'//V2(1)//V2(2)//V2(3)//&
        'k'//V3(1)//V3(2)//V3(3)//'.fml'

      RETURN

    END SUBROUTINE NAMER2
