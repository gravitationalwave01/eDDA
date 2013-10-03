    SUBROUTINE WRITEFML(JPBC,IORI,IORTH,IRAD,IWAV,MXCOMP,MXSCA, &
                        NAT0,NAVG,NCOMP,NSCAT,NORICHAR,CALPHA,  &
                        CDESCR,CMDFFT,CMDSOL,CMDFRM,CSHAPE,     &
                        CSTAMP,AEFF,AK1,AKR,NAMBIENT,BETAD,     &
                        PHID,THETAD,TOL,WAVE,XX,A1,A2,PHIN,     &
                        THETAN,CXE01,CXE01R,CXE02,              &
                        CXE02R,CXEPS,CXRFR,CXF11,CXF21,CXF12,   &
                        CXF22,PYD,PZD)                          !
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

! Scalar Arguments ..

      REAL(WP) :: AEFF,AK1,BETAD,NAMBIENT,PHID,PYD,PZD,THETAD,TOL,WAVE,XX
      INTEGER :: IORI,IORTH,IRAD,IWAV,JPBC,MXCOMP,MXSCA, &
                 NAT0,NAVG,NCOMP,NORICHAR,NSCAT          !

      CHARACTER :: CALPHA*6,CMDFFT*6,CMDFRM*6, CMDSOL*6, &
                   CSHAPE*9,CFLFML*19,CSTAMP*26,CDESCR*67

! Array Arguments ..

      COMPLEX(WP) ::   &
        CXE01(3),      &
        CXE01R(3),     &
        CXE02(3),      &
        CXE02R(3),     &
        CXEPS(MXCOMP), &
        CXF11(MXSCA),  &
        CXF12(MXSCA),  &
        CXF21(MXSCA),  &
        CXF22(MXSCA),  &
        CXRFR(MXCOMP)
      REAL(WP) ::    &
        A1(3),       &
        A2(3),       &
        AKR(3),      &
        PHIN(MXSCA), &
        THETAN(MXSCA)

! Local variables:

      CHARACTER :: CFRAME*12
      INTEGER :: J,ND
      REAL(WP) :: DEGRAD,MKD,PHIND,PI,RVAR1,RVAR2,RVAR3,SINALPHA,THETND
      COMPLEX(WP) :: CXFAC,CXI

! External Subroutines:

      EXTERNAL NAMER2

! Intrinsic Functions:

      INTRINSIC CONJG,SQRT
      DATA CFLFML/'                   '/
      SAVE CFLFML

!***********************************************************************
! Subroutine WRITEFML

! Purpose of this module is to write out the complex scattering
! amplitudes f_ml

! History:
! 07.08.31 (BTD) created using WRITESCA as template
! 08.06.30 (BTD) changed normalization factor CXFAC for case JPBC=3
!                [target periodic in y and z directions]
! 10.05.09 (BTD) added CMDSOL to argument list
! 11.08.30 (BTD) v7.2.2 added NAMBIENT to argument list
! 13.03.22 (BTD) v7.3.0
!                * CFLFML: initialize and SAVE
!                * extend to CHARACTER*19 to allow up to 1e6 orientations
! end history

! Copyright (C) 2007,2008,2010,2011,2013
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

! diagnostic
!      write(0,*)'writefml ckpt 0, NAT0=',NAT0
!

! If IORI > 0, then write out scattering properties for this specific
! target orientation.

      PI=4._WP*ATAN(1._WP)
      DEGRAD=180._WP/PI
      CXI=(0._WP,1._WP)

      IF(CMDFRM=='LFRAME' .OR. CMDFRM=='lframe')THEN
        CFRAME='Lab Frame'
      ELSEIF(CMDFRM=='TFRAME' .OR. CMDFRM=='tframe')THEN
        CFRAME='Target Frame'
      ELSE
        CFRAME='????????'
      ENDIF
!*** diagnostic
!      write(0,*)'writefml ckpt 1' 
!      write(0,*)'   iwav,irad,iori=',iwav,irad,iori
!      write(0,*)'   norichar=',norichar
!***
      CALL NAMER2(IWAV-1,IRAD-1,IORI-1,NORICHAR,CFLFML)
!*** diagnostic
!      write(0,*)'writefml ckpt 1, cflfml=',cflfml
!***
      OPEN (UNIT=8,FILE=CFLFML,STATUS='UNKNOWN')
      WRITE(8,FMT=9030)CSTAMP,CDESCR,CMDSOL,CALPHA,CSHAPE,NAT0

      WRITE(8,FMT=9032)AEFF,WAVE,XX,NAMBIENT
      DO J=1,NCOMP
         MKD=(4._WP*PI/(3._WP*NAT0))**(1._WP/3._WP)* &
             SQRT(REAL(CXRFR(J)*CONJG(CXRFR(J))))*XX
         WRITE(8,FMT=9031)CXRFR(J),CXEPS(J),MKD,J
      ENDDO
      WRITE(8,FMT=9033)TOL,NAVG,A1,A2
      WRITE(8,FMT=9035)AKR,CXE01R,CXE02R
      RVAR1=AK1
      RVAR2=0._WP
      RVAR3=0._WP
      WRITE(8,FMT=9037)RVAR1,RVAR2,RVAR3,CXE01,CXE02
      WRITE(8,FMT=9040)BETAD,THETAD,PHID

      IF(JPBC==0)THEN
         CXFAC=1._WP
         WRITE(8,FMT=9050)

! for periodic structures (JPC=1,2, or 3), it is assumed that 
! target axis a_1 = x_TF
!             a_2 = y_TF
!             a_3 = z_TF
! JPBC=1 : target is periodic in y direction
!      2                         z direction
!      3                         y and z directions

      ELSEIF(JPBC==1)THEN

         SINALPHA=SQRT(1._WP-(SIN(THETAD/DEGRAD)*COS(BETAD/DEGRAD))**2)
         CXFAC=SQRT(2._WP*PI*CXI/SINALPHA)/(AK1*PYD)
         RVAR1=AK1*PYD
         WRITE(8,FMT=9051)PYD,RVAR1
      ELSEIF(JPBC==2)THEN
         SINALPHA=SQRT(1._WP-(SIN(THETAD/DEGRAD)*SIN(BETAD/DEGRAD))**2)
         CXFAC=SQRT(2._WP*PI*CXI/SINALPHA)/(AK1*PZD)
         RVAR1=AK1*PZD
         WRITE(8,FMT=9052)PZD,RVAR1
      ELSEIF(JPBC==3)THEN
! 080630 BTD not certain how we want to set normalization factor CXFAC
!            THETAD = angle relative to normal
!            ALPHA  = 90-THETAD
!                   = angle between scattered vector and surface
         SINALPHA=COS(THETAD/DEGRAD)
         CXFAC=2._WP*PI*CXI/(AK1*PYD*AK1*PZD*SINALPHA)
         RVAR1=AK1*PYD
         RVAR2=AK1*PZD
         WRITE(8,FMT=9053)PYD,PZD,RVAR1,RVAR2
      ENDIF

      DO ND=1,NSCAT

! convert scattering angles to degrees

         PHIND=DEGRAD*PHIN(ND)
         THETND=DEGRAD*THETAN(ND)
         IF(IORTH==1)THEN
            WRITE(8,FMT=9070)THETND,PHIND,                        &
                             (CXFAC*CXF11(ND)),(CXFAC*CXF21(ND))
         ELSE
            WRITE(8,FMT=9071)THETND,PHIND,                        &
                             (CXFAC*CXF11(ND)),(CXFAC*CXF21(ND)), &
                             (CXFAC*CXF12(ND)),(CXFAC*CXF22(ND))
         ENDIF
      ENDDO

9030  FORMAT (' DDSCAT --- ',A/' TARGET ---',A/' ',A,                  &
        ' --- method of solution '/' ',A,                              &
        ' --- prescription for polarizabilies'/' ',A,' --- shape '/I7, &
        ' = NAT0 = number of dipoles')
! following code to be enabled for noncubic treatment
! 9030 FORMAT(' DDSCAT --- ',A,/,' TARGET ---',A,/,' ',A,
!     &       ' --- method of solution ',/,' ',A,
!     &       ' --- prescription for polarizabilies',/,' ',A,
!     &       ' --- shape ',/,I7,' = NAT0 = number of dipoles',
!     &       3F8.5,' = normalized lattice spacings dx,dy,dz')
9031  FORMAT ('n= (',F7.4,' , ',F7.4,'),  eps.= (',F8.4,' , ',F7.4, &
        ')  |m|kd=',0P,F8.4,' for subs.',I2)
9032  FORMAT(                                                         &
         '  AEFF=',F11.5,' =',' effective radius (physical units)',/, &
         '  WAVE=',F11.5,' = vacuum wavelength (physical units)',/,   &
         'K*AEFF=',F11.5,' = 2*pi*aeff/lambda',/,                     &
         'NAMBIENT=',F9.5,' = refractive index of ambient medium')
9033  FORMAT ('   TOL=',1P,E10.3,' = error tolerance for CCG method'/ &
        '  NAVG=',I6,' = (theta,phi) values used in comp. of Qsca,g'/ &
        0P,'(',F8.5,2F9.5,') = target axis A1 in Target Frame'/       &
        1X,'(',F8.5,2F9.5,') = target axis A2 in Target Frame')
9035  FORMAT (1X,'(',F8.5,2F9.5,') = k vector (latt. units) in TF',/,1X, &
        3('(',F8.5,',',F8.5,')'),'=inc.pol.vec. 1 in TF',/,1X,           &
        3('(',F8.5,',',F8.5,')'),'=inc.pol.vec. 2 in TF')
9037  FORMAT (1X,'(',F8.5,2F9.5,') = k vector (latt. units) in Lab Frame',/, &
        1X,3('(',F8.5,',',F8.5,')'),'=inc.pol.vec. 1 in LF',/,1X,            &
        3('(',F8.5,',',F8.5,')'),'=inc.pol.vec. 2 in LF')
9040  FORMAT (' BETA =',F7.3,' = rotation of target around A1',/,          &
        ' THETA=',F7.3,' = angle between A1 and k',/,                      &
        '  PHI =',F7.3,' = rotation of A1 around k')
9050  FORMAT(5X,'Finite target:',/,                                        &
        5X,'e_m dot E(r) = i*exp(ikr)*f_ml*E_inc(0)/(kr)',/,               &
        5X,'m=1 in scatt. plane, m=2 perp to scatt. plane',/,/,            &
        ' theta   phi',2X,'Re(f_11)',3X,'Im(f_11)',3X,'Re(f_21)',3X,       &
        'Im(f_21)',3X,'Re(f_12)',3X,'Im(f_12)',3X,'Re(f_22)',3X,'Im(f_22)')
9051  FORMAT(                                                              &
        5X,'Target periodic in y direction with L_y/d=',1P,E10.3,2X,       &
        'kL_y = ',E10.3,/,                                                 &
        5X,'e_m dot E(r) = exp(ikr)*f_ml*[e_l dot E_inc(0)]/sqrt(kR)',     &
        5X,'[R = dist. from target]',/,                                    &
        5X,'m=1 in scatt. plane, m=2 perp to scatt. plane',/,/,            &
        ' alpha  zeta',2X,'Re(f_11)',3X,'Im(f_11)',3X,'Re(f_21)',3X,       &
        'Im(f_21)',3X,'Re(f_12)',3X,'Im(f_12)',3X,'Re(f_22)',3X,'Im(f_22)')
9052  FORMAT(                                                              &
        5X,'Target periodic in z direction with L_z/d=',1P,E10.3,2X,       &
        'kL_z = ',E10.3,/,                                                 &
        5X,'e_m dot E(r) = exp(ikr)*f_ml*[e_l dot E_inc(0)]/sqrt(kR)',     &
        5X,'[R = dist. from target]',/,                                    & 
        5X,'m=1 in scatt. plane, m=2 perp to scatt. plane',/,/,            &
        ' alpha  zeta',2X,'Re(f_11)',3X,'Im(f_11)',3X,'Re(f_21)',3X,       &
        'Im(f_21)',3X,'Re(f_12)',3X,'Im(f_12)',3X,'Re(f_22)',3X,'Im(f_22)')
9053  FORMAT(5X,                                                           &
        'Target periodic in y,z directions with L_y/d,L_z/d=',1P,2E10.3,   &
        2X,'kL_y=',E10.3,' kL_z=',E10.3,/,                                 &
        5X,' e_m dot E(r) = exp(ikr)*f_ml*[e_l dot E_inc(0)]',/,           &
        5X,'m=1 in scatt. plane, m=2 perp to scatt. plane',/,/,            &
        ' alpha   phi',2X,'Re(f_11)',3X,'Im(f_11)',3X,'Re(f_21)',3X,       &
        'Im(f_21)',3X,'Re(f_12)',3X,'Im(f_12)',3X,'Re(f_22)',3X,'Im(f_22)')
9070  FORMAT(2F6.1,1P,4E11.3)
9071  FORMAT(2F6.1,1P,8E11.3)

    END SUBROUTINE WRITEFML
