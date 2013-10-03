    SUBROUTINE NAMER2(IWAV0,IRAD0,IDIR0,NORICHAR,CFLFML)
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
!         IDIR=integer (0-...) identifying direction
! Returns:
!         CFLFML =name for output file containing complex scattering
!                 amplitudes f_ml for selected directions
!
! B.T.Draine, Princeton Univ. Observatory, 2007
! History:
! 07.08.31 (BTD): Adapted from NAMER
! 13.03.22 (BTD): namer2_v2
!                 * add support for up to 1e6 orientations
! end history
! Copyright (C) 2007,2013 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
! Arguments:

      INTEGER :: IDIR,IDIR0,IRAD,IRAD0,IWAV,IWAV0,NORICHAR
      CHARACTER ::         &
         CFLFML(13+NORICHAR)

! Local variables:

      CHARACTER :: N(0:9)*1,V1(3)*1,V2(3)*1,V3(1:NORICHAR)*1

! DATA statements:

      DATA N/'0','1','2','3','4','5','6','7','8','9'/
      SAVE N

!***********************************************************************
!*** Filename CFLSCA will be of the form w001r01k001
!     (for IWAV=2,IRAD=2,IDIR=2)

!*** diagnostic
!      write(0,*)'namer2 ckpt 0, norichar=',norichar
!      write(0,*)'       iwav0,irad0,idir0=',iwav0,irad0,idir0
!***
      V1(1)=N(IWAV0/100)
      V1(2)=N((IWAV0-100*(IWAV0/100))/10)
      V1(3)=N(IWAV0-10*(IWAV0/10))
      V2(1)=N(IRAD0/100)
      V2(2)=N((IRAD0-100*(IRAD0/100))/10)
      V2(3)=N(IRAD0-10*(IRAD0/10))
      IF(NORICHAR.EQ.1)THEN
         V3(1)=N(IDIR0)
      ELSEIF(NORICHAR.EQ.2)THEN
         V3(1)=N(IDIR0/10)
         V3(2)=N(IDIR0-10*(IDIR0/10))
      ELSEIF(NORICHAR.EQ.3)THEN
         V3(1)=N(IDIR0/100)
         V3(2)=N((IDIR0-100*(IDIR0/100))/10)
         V3(3)=N(IDIR0-10*(IDIR0/10))
      ELSEIF(NORICHAR.EQ.4)THEN
         V3(1)=N(IDIR0/1000)
         V3(2)=N((IDIR0-1000*(IDIR0/1000))/100)
         V3(3)=N((IDIR0-100*(IDIR0/100))/10)
         V3(4)=N(IDIR0-10*(IDIR0/10))
      ELSEIF(NORICHAR.EQ.5)THEN
         V3(1)=N(IDIR0/10000)
         V3(2)=N((IDIR0-10000*(IDIR0/10000))/1000)
         V3(3)=N((IDIR0-1000*(IDIR0/1000))/100)
         V3(4)=N((IDIR0-100*(IDIR0/100))/10)
         V3(5)=N(IDIR0-10*(IDIR0/10))
      ELSEIF(NORICHAR.EQ.6)THEN
         V3(1)=N(IDIR0/100000)
         V3(2)=N((IDIR0-100000*(IDIR0/100000))/10000)
         V3(3)=N((IDIR0-10000*(IDIR0/10000))/1000)
         V3(4)=N((IDIR0-1000*(IDIR0/1000))/100)
         V3(5)=N(IDIR0-100*(IDIR0/100))
         V3(6)=N(IDIR0-10*(IDIR0/10))
      ELSE
         WRITE(0,*)'>NAMER2 fatal error: ',                 &
                   'cannot handle more than 1e6 orientations'
         STOP
      ENDIF

      CFLFML(1)='w'
      CFLFML(2:4)=V1(1:3)
      CFLFML(5)='r'
      CFLFML(6:8)=V2(1:3)
      CFLFML(9)='k'
      CFLFML(10:9+NORICHAR)=V3(1:NORICHAR)
      CFLFML(10+NORICHAR)='.'
      CFLFML(11+NORICHAR)='f'
      CFLFML(12+NORICHAR)='m'
      CFLFML(13+NORICHAR)='l'

!*** diagnostic
!      write(0,*)' namer2 ckpt 9, cflfml=',cflfml
!***
      RETURN

    END SUBROUTINE NAMER2
