    SUBROUTINE NAMER(IWAV0,IRAD0,IDIR0,NORICHAR,CFLPOL1,CFLPOL2, &
                     CFLSCA,CFLAVG,CFLE1,CFLE2,CFLEB1,CFLEB2)    !
      IMPLICIT NONE
!**************************** namer_v4 *********************************

!***********************************************************************
! Purpose: to generate file names for specific (wavelength,
!                                               size,
!                                               orientation)
! Present version allows up to 1000 wavelengths, (000-999)
!                              1000 sizes        (000-999)
!                              arbitrary number of orientations (000-...)
! Given:
!         IWAV0=integer (0-999) identifying wavelength
!         IRAD0=integer (0-999) identifying size
!         IDIR0=integer (0-...) identifying orientation
!         NORICHAR=integer = 1 if          NORI <     11
!                            2 if     10 < NORI <    101
!                            3 if    100 < NORI <   1001
!                            4 if   1000 < NORI <  10001
!                            5 if  10000 < NORI < 100000
!                            6 if 100000 < NORI <1000000 etc
! Returns:
!         CFLSCA =name for output file containing Qext,Qabs,Qpha,Qsca,g
!                 and fml values for selected directions
!         CFLAVG =name for output file for given wavelength and size,
!                 containing scattering etc. averaged over inc.dir.
!         CFLPOL1=name for output file for given wavelength and size,
!                 containing polarization vector for input pol 1
!         CFLPOL2=name for output file for given wavelength and size,
!                 containing polarization vector for input pol 2
!         CFLE1  =name for "nearfield" output file with E in rectangular
!                 volume for input pol 1
!         CFLE2  =name for "nearfield" output file with E in rectangular
!                 volume for input pol 2
!         CFLEB1  =name for "nearfield" output file with E and B in rectangular
!                 volume for input pol 1
!         CFLEB2  =name for "nearfield" output file with E and B in rectangular
!                 volume for input pol 2

! B.T.Draine, Princeton Univ. Observatory, 1988
! History:
! 90.11.21 (BTD): Modified to allow 3 digits to identify direction
!                 and to change indices to begin from w00r00k000.sca
! 05.10.11 (BTD): Modified to allow 3 digits to identify wavelength
! 06.04.08 (BTD): Rewritten and extended to CFLPOL1,CFLPOL2
! 07.08.31 (BTD): Added third digit to size identifier IDIR
! 11.08.11 (BTD): v7.2.1
!                 * Added CFLE1 and CFLE2 to argument list and code
! 12.07.10 (BTD): v7.3 and namer_v2
!                 * Added CFLB1 and CFLB2 to argument list and code
! 12.12.25 (BTD): namer_v3
!                 * change CFLB1 and CFLB2 to CFLEB1 and CFLEB2
!                 * make CFLEB1,CFLEB2 character*16
! 13.03.22 (BTD): namer_v4
!                 * support for up to 1e6 orientations
!                 * add NORICHAR to argument list
! end history
! Copyright (C) 1993,2005,2006,2007,2011,2012,2013 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
! Arguments:

      INTEGER :: IDIR,IDIR0,IRAD,IRAD0,IWAV,IWAV0,NORICHAR
      CHARACTER ::             &
         CFLAVG(12),           &
         CFLE1(12+NORICHAR),   &
         CFLE2(12+NORICHAR),   &
         CFLEB1(13+NORICHAR),  &
         CFLEB2(13+NORICHAR),  &
         CFLSCA(13+NORICHAR),  &
         CFLPOL1(14+NORICHAR), &
         CFLPOL2(14+NORICHAR)  !

! Local variables:

      CHARACTER :: N(0:9)*1,V1(3)*1,V2(3)*1,V3(1:NORICHAR)*1

! DATA statements:

      DATA N/'0','1','2','3','4','5','6','7','8','9'/
      SAVE N

!***********************************************************************
!*** Filename CFLSCA will be of the form w001r001k001
!     (for IWAV=2,IRAD=2,IDIR=2)

!*** diagnostic
!      write(0,*)'namer ckpt 0'
!      write(0,*)' norichar=',norichar
!      write(0,*)'iwav0,irad0,idir0=',iwav0,irad0,idir0
!***
      V1(1)=N(IWAV0/100)
      V1(2)=N((IWAV0-100*(IWAV0/100))/10)
      V1(3)=N(IWAV0-10*(IWAV0/10))
      V2(1)=N(IRAD0/100)
      V2(2)=N((IRAD0-100*(IRAD0/100))/10)
      V2(3)=N(IRAD0-10*(IRAD0/10))
      IF(NORICHAR.EQ.1)THEN
         V3(1)=N(IDIR0)
!*** diagnostic
!         write(0,*)'namer ckpt 2, norichar=',norichar
!***     
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
         V3(5)=N((IDIR0-100*(IDIR0/100))/10)
         V3(6)=N(IDIR0-10*(IDIR0/10))
      ELSE
         WRITE(0,*)'>NAMER fatal error: ',                  &
                   'cannot handle more than 1e6 orientations'
         STOP
      ENDIF
      CFLAVG(1)='w'
      CFLAVG(2:4)=V1(1:3)
      CFLAVG(5)='r'
      CFLAVG(6:8)=V2(1:3)
      CFLAVG(9)='.'
      CFLAVG(10)='a'
      CFLAVG(11)='v'
      CFLAVG(12)='g'

      CFLSCA(1:8)=CFLAVG(1:8)
      CFLSCA(9)='k'
      CFLSCA(10:9+NORICHAR)=V3(1:NORICHAR)
      CFLSCA(10+NORICHAR)='.'
      CFLSCA(11+NORICHAR)='s'
      CFLSCA(12+NORICHAR)='c'
      CFLSCA(13+NORICHAR)='a'

      CFLPOL1(1:10+NORICHAR)=CFLSCA(1:10+NORICHAR)
      CFLPOL1(11+NORICHAR)='p'
      CFLPOL1(12+NORICHAR)='o'
      CFLPOL1(13+NORICHAR)='l'
      CFLPOL1(14+NORICHAR)='1'

      CFLPOL2(1:13+NORICHAR)=CFLPOL1(1:13+NORICHAR)
      CFLPOL2(14+NORICHAR)='2'

      CFLE1(1:10+NORICHAR)=CFLPOL1(1:10+NORICHAR)
      CFLE1(11+NORICHAR)='E'
      CFLE1(12+NORICHAR)='1'

      CFLE2(1:11+NORICHAR)=CFLE1(1:11+NORICHAR)
      CFLE2(12+NORICHAR)='2'

      CFLEB1(1:11+NORICHAR)=CFLE1(1:11+NORICHAR)
      CFLEB1(12+NORICHAR)='B'
      CFLEB1(13+NORICHAR)='1'

      CFLEB2(1:12+NORICHAR)=CFLEB1(1:12+NORICHAR)
      CFLEB2(13+NORICHAR)='2'
!*** diagnostic
!      write(0,*)'namer ckpt 9: cflavg=',cflavg
!      write(0,*)'cflsca=',cflsca
!      write(0,*)'cflpol1=',cflpol1
!      write(0,*)'cflpol2=',cflpol2
!      write(0,*)'cfle1=',cfle1
!      write(0,*)'cfle2=',cfle2
!      write(0,*)'cfleb1=',cfleb1
!      write(0,*)'cfleb2=',cfleb2
!***
      RETURN

    END SUBROUTINE NAMER
