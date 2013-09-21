    SUBROUTINE NAMER(IWAV,IRAD,IDIR,CFLPOL1,CFLPOL2,CFLSCA,CFLAVG)
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
!         CFLSCA =name for output file containing Qext,Qabs,Qpha,Qsca,g
!                 and fml values for selected directions
!         CFLAVG =name for output file for given wavelength and size,
!                 containing scattering etc. averaged over inc.dir.
!         CFLPOL1=name for output file for given wavelength and size,
!                 containing polarization vector for input pol 1
!         CFLPOL2=name for output file for given wavelength and size,
!                 containing polarization vector for input pol 2

! B.T.Draine, Princeton Univ. Observatory, 1988
! History:
! 90.11.21 (BTD): Modified to allow 3 digits to identify direction
!                 and to change indices to begin from w00r00k000.sca
! 05.10.11 (BTD): Modified to allow 3 digits to identify wavelength
! end history
! 06.04.08 (BTD): Rewritten and extended to CFLPOL1,CFLPOL2
! 07.08.31 (BTD): Added third digit to size identifier IDIR
! Copyright (C) 1993,2005,2006,2007 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
! Arguments:

      INTEGER :: IDIR,IDIR0,IRAD,IRAD0,IWAV,IWAV0
      CHARACTER :: CFLAVG*15,CFLPOL1*17,CFLPOL2*17,CFLSCA*16

! Local variables:

      CHARACTER :: N(10)*1,V1(3)*1,V2(3)*1,V3(3)*1

! DATA statements:

      DATA N/'0','1','2','3','4','5','6','7','8','9'/
      SAVE N

!***********************************************************************
!*** Filename CFLSCA will be of the form w001r001k001
!     (for IWAV=2,IRAD=2,IDIR=2)

!*** Set IWAV0=IWAV-1 to run over range 0-999
!    Set IRAD0=IRAD-1 to run over range 0-999
!    Set IDIR0=IDIR-1 to run over range 0-999

      IWAV0=IWAV-1
      IRAD0=IRAD-1
      IDIR0=IDIR-1
      V1(1)=N(IWAV0/100+1)
      V1(2)=N((IWAV0-100*(IWAV0/100))/10+1)
      V1(3)=N(IWAV0-10*(IWAV0/10)+1)
      V2(1)=N(IRAD0/100+1)
      V2(2)=N((IRAD0-100*(IRAD0/100))/10+1)
      V2(3)=N(IRAD0-10*(IRAD0/10)+1)
      V3(1)=N(IDIR0/100+1)
      V3(2)=N((IDIR0-100*(IDIR0/100))/10+1)
      V3(3)=N(IDIR0-10*(IDIR0/10)+1)

      CFLAVG='w'//V1(1)//V1(2)//V1(3)//'r'//V2(1)//V2(2)//V2(3)//'.avg'
      CFLSCA='w'//V1(1)//V1(2)//V1(3)//'r'//V2(1)//V2(2)//V2(3)// &
        'k'//V3(1)//V3(2)//V3(3)//'.sca'
      CFLPOL1='w'//V1(1)//V1(2)//V1(3)//'r'//V2(1)//V2(2)//V2(3)// &
        'k'//V3(1)//V3(2)//V3(3)//'.pol1'
      CFLPOL2='w'//V1(1)//V1(2)//V1(3)//'r'//V2(1)//V2(2)//V2(3)// &
        'k'//V3(1)//V3(2)//V3(3)//'.pol2'

      RETURN

    END SUBROUTINE NAMER
