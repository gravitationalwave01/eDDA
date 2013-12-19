    SUBROUTINE WRITESCA(MXNX,MXNY,MXNZ,NX,NY,NZ,DX,ICOMP,NAMBIENT,IXYZ0,JPBC,  &
                        MXNAT,MXN3,MYID,NAT,NAT3,WAVEA,MXWAV,NWAV,AEFFA,MXRAD, &
                        NRAD,ITHETA,IBETA,IPHI,THETA,MXTHET,NTHETA,BETA,       &
                        MXBETA,NBETA,PHI,MXPHI,NPHI,NSMELTS,TIMERS,MXTIMERS,   &
                        NTIMERS,CBINFILE,IOBIN,CBINFLAG,IDNC,ILIN10,ILIN12,    &
                        IORI,IWRKSC,IORTH,IRAD,IWAV,MXCOMP,MXSCA,NAT0,NAVG,    &
                        ITNUM,MXITER,NCOMP,NORI,NORICHAR,NSCAT,CALPHA,CDESCR,  &
                        CFLAVG,CFLSCA,CMDFFT,CMDSOL,CMDFRM,CMDTRQ,CSHAPE,      &
                        CSTAMP,AEFF,AK1,AK_TF,BETAD,BETMID,BETMXD,ETASCA,PHID, &
                        PHIMID,PHIMXD,THETAD,THTMID,THTMXD,TOL,WAVE,XX,A1_TF,  &
                        A2_TF,PHIN,QABS,QABSUM,QBKSCA,QBKSUM,QEXSUM,QEXT,QPHA, &
                        QPHSUM,QSCAG,QSCAG2,QSCAT,QSCGSUM,QSCG2SUM,QSCSUM,     &
                        QTRQAB,QTRQABSUM,QTRQSC,QTRQSCSUM,S1111,S2121,SM,      &
                        SMIND1,SMIND2,SMORI,THETAN,CX1121,CXE01_LF,CXE01_TF,   &
                        CXE02_LF,CXE02_TF,CXEPS,CXRFR,CXF11,CXF21,CXF12,CXF22, &
                        PYDDX,PZDDX,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX)             !
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE
!--------------------- writesca_v7 ------------------------------------------
! Scalar Arguments ..

      INTEGER :: IBETA,IDNC,ILIN10,ILIN12,IOBIN,IORI,IORTH,IPHI,IRAD,ITHETA, &
         IWAV,IWRKSC,JPBC,MXBETA,MXCOMP,MXITER,MXN3,MXNAT,MXNX,MXNY,MXNZ,    &
         MXPHI,MXRAD,MXSCA,MXTHET,MXTIMERS,MXWAV,MYID,NAT,NAT0,NAT3,NAVG,    &
         NBETA,NCOMP,NORI,NORICHAR,NPHI,NRAD,NSCAT,NSMELTS,NTHETA,NTIMERS,   &
         NWAV,NX,NY,NZ
      REAL(WP) :: AEFF,AK1,BETAD,BETMID,BETMXD,ETASCA,NAMBIENT,        &
         PHID,PHIMID,PHIMXD,PYDDX,PZDDX,THETAD,THTMID,THTMXD,TOL,WAVE, &
         XMAX,XMIN,XX,YMAX,YMIN,ZMAX,ZMIN
      INTEGER :: ITNUM(2)
      CHARACTER :: CALPHA*6,CBINFLAG*6,CMDFFT*6,CMDFRM*6,CMDSOL*6,CMDTRQ*6, &
         CSHAPE*9,CFLAVG*12,CBINFILE*14,CFLSCA*19,CSTAMP*26,CDESCR*67

! Array Arguments ..

      COMPLEX(WP) ::    &
         CX1121(MXSCA), &
         CXE01_LF(3),   &
         CXE01_TF(3),   &
         CXE02_LF(3),   &
         CXE02_TF(3),   &
         CXEPS(MXCOMP), &
         CXF11(MXSCA),  &
         CXF12(MXSCA),  &
         CXF21(MXSCA),  &
         CXF22(MXSCA),  &
         CXRFR(MXCOMP)
      REAL(WP) ::          &
         A1_TF(3),         &
         A2_TF(3),         &
         AEFFA(MXRAD),     &
         AK_TF(3),         &
         BETA(MXBETA),     &
         DX(3),            &
         PHI(MXPHI),       &
         PHIN(MXSCA),      &
         QABS(2),          &
         QABSUM(2),        &
         QBKSCA(2),        &
         QBKSUM(2),        &
         QEXSUM(2),        &
         QEXT(2),          &
         QPHA(2),          &
         QPHSUM(2),        &
         QSCAG(3,2),       &
         QSCAG2(2),        &
         QSCAT(2),         &
         QSCG2SUM(2),      &
         QSCGSUM(3,2),     &
         QSCSUM(2),        &
         QTRQAB(3,2),      &
         QTRQABSUM(3,2),   &
         QTRQSC(3,2),      &
         QTRQSCSUM(3,2),   &
         S1111(MXSCA),     &
         S2121(MXSCA),     &
         SM(4,4,MXSCA),    &
         SMORI(4,4,MXSCA), &
         THETA(MXTHET),    &
         THETAN(MXSCA),    &
         TIMERS(MXTIMERS), &
         WAVEA(MXWAV)
      INTEGER*2 :: &
         ICOMP(MXN3)
      INTEGER ::        &
         IXYZ0(NAT0,3), &
         SMIND1(9),     &
         SMIND2(9)

! Local variables:

      CHARACTER :: CFLE1*18,CFLE2*18,CFLPOL1*20,CFLPOL2*20,CFRAME*12, &
                   CFLEB1*19,CFLEB2*19
      COMPLEX(WP) :: CXVAR1
      REAL(WP) :: DAEFF,DEGRAD,DPHYS,MKD,PHIND,PI,PPOL,RVAR1,RVAR2,RVAR3,THETND,ENERGY 
                  !Energy added by SMC 8.5.13 following NWB 7/12/12
      INTEGER :: I,J,ND
      REAL(WP) ::  &
         A1_LF(3), &
         A2_LF(3), &
         ABSCO(2), &
         G(2),     &
         G2(2),    &
         QAV(18)

! External Subroutines:

      EXTERNAL NAMER,WRITEBIN

! Intrinsic Functions:

      INTRINSIC CONJG,SQRT

! Save statement:

      SAVE DEGRAD

! Data statements:

      DATA DEGRAD/57.29578_WP/
!***********************************************************************
! Subroutine WRITESCA
! Purpose of WRITESCA is to collect a lot of the output code into
! a single module.
! Given:
!    JPBC   = 0 for isolated finite target
!           = 1 for target periodic in y direction in Target Frame
!           = 2 for target periodic in z direction in Target Frame
!           = 3 for target periodic in y and z directions in TF
!
! Returns: CFLE1,CFLE2 = file names for stored P and E
!          CFLEB1,CFLEB2 = file names for stored P,E, and B

! If IORI = 0, then print orientational average
! If IORI > 0, print output for specific orientation

! History:
! 96.11.14 (PJF) added binary write options (see also info in
!                "versions.f" about IDL postprocessing)
!                possible "checkpointing" could be done here "close(iobi
!                open(iobin,file=cbinfile,form='unformatted',access='app
!                (which is non-standard Fortran 77 (but standard in F90)
!                open(iobin,file=cbinfile,form='unformatted',status='old
! 96.11.21 (BTD) added call to WRITENET if called with CNETFLAG='NCCLOS'
!                in order to get call to NCCLOS out of DDSCAT.f
! 98.01.11 (BTD) added DX(3) to argument list to pass information on
!                nonuniform lattice spacing [development of DDSCAT.6.0]
! 98.11.16 (BTD) corrected CX01R -> CXE01R in argument list of WRITEBIN
! 98.12.12 (BTD) Introduce NSMELTS,SMIND1,SMIND2 to argument list to
!                allow selection of Muller matrix elements to be
!                printed out.  Add code to allow from 1 to 9 matrix
!                elements to be printed.
! 98.12.13 (BTD) Modify to print out validity criterion |m|kd
!                for each composition
! 03.01.30 (BTD) disabled output line to print out dx,dy,dz
! 03.02.13 (BTD) added one more digit of precision to output Q values
!                with appropriate modifications to format statements
! 03.04.14 (BTD) added ITNUM(2) and MXITER to argument list
!                added ITNUM(2) = number of iterations
!                and MXITER = max number of iterations allowed
!                to output statements 9150 and 9155
! 03.09.01 (BTD) cosmetic change to output statement: g -> g(1)
! 03.10.23 (BTD) removed ICTHM and IPHIM from argument list
!                removed ICTHM and IPHIM from argument list of
!                WRITENET and WRITEBIN
!                removed ICTHM and IPHIM from WRITE(..,9033)
!                added NAVG to argument list
!                added NAVG to WRITE(..,9033)
!                added NAVG to argument list of WRITENET and WRITEBIN
!                write NAVG to qtable,wxxrxxkxxx.sca, wxxrxxori.avg
! 04.02.22 (BTD) added IWRKSC to argument list of WRITESCA and use to
!                control writing to ascii output files
! 04.05.21 (BTD) reordered QSCGSUM and QSCG2SUM in argument list
! 04.12.29 (BTD) corrected error in writing Mueller matrix for the
!                "sca" files when IWRKSC=1
! 05.10.11 (BTD) changed CFLAVG*13 to CFLAVG*14
!                changed CFLSCA*14 to CFLSCA*15
! 05.10.11 (BTD) added CMDFRM to argument list
!                added variable CFRAME
!                modified 9080 FORMAT to print out string CFRAME
!                specifying whether (theta,phi) are in Lab Frame or
!                Target Frame
! 06.10.07 (BTD) added JPBC to argument list
! 06.10.08 (BTD) modified output statements when JPBC=1,2,or 3
! 06.11.28 (BTD) corrected error in calculation of degree of linear
!                polarization of scattered light (two places below)
! 07.06.22 (BTD) original calculation of PPOL was correct -- return
!                to original formula
! 07.08.13 (BTD) modified output format statements to report scattering
!                angles to 0.01 degree
! 07.09.11 (BTD) changed IXYZ0 from INTEGER*2 to INTEGER
! 08.01.12 (BTD) v2
!                * added PYDDX,PZDDX to argument list
!                * for JPBC=3, compute ABSCO(J) = absorption coefficient
!                * for JPBC=3, write out ABSCO(J) and average over two
!                  incident polarizations
! 08.03.15 (BTD) v7.0.5
!                changed IXYZ0(MXNAT,3) -> IXYZ0(NAT0,3)
! 08.04.13 (BTD) v7.0.5
!                * modified format so the first two S_ij elements are
!                  written out with 5 sig figs instead of 4 
! 08.05.09 (BTD) * added MYID to argument list
! 08.05.12 (BTD) * following suggestion from Art Lazanoff, initialize
!                  DAEFF=0._WP and DPHYS=0._WP [WHY?]
! 08.05.29 (BTD) v7.0.6
!                * changed declarations
!                  SM(MXSCA,4,4) -> SM(4,4,MXSCA)
!                  SMORI(MXSCA,4,4) -> SMORI(4,4,MXSCA)
!                * made corresponding changes to write statements
! 08.07.16 (BTD) * change format of w000r000k000.sca and w000r000.avg
!                  files when JPBC > 0, since in these cases we do not
!                  carry out sums over scattering directions.
! 08.07.17 (BTD) * added A1L,A2L,AK to w000r000k000.sca
! 08.08.29 (BTD) * Remove calls to WRITENET -- no longer supported
!                * Removed CNETFLAG and CNETFILE from argument list.
! 11.08.16 (BTD) v5 / v7.2.1
! 11.08.30 (BTD) v5 / v7.2.2
!                * add NAMBIENT to arg list
!                * write NAMBIENT to .sca output file
! 11.11.18 (BTD) * Change notation
!                  A1     -> A1_TF
!                  A2     -> A2_TF
!                  A1L    -> A1_LF
!                  A2L    -> A2_LF
!                  AKR    -> AK_TF
!                  CXE01  -> CXE01_LF
!                  CXE02  -> CXE02_LF
!                  CXE01R -> CXE01_TF
!                  CXE02R -> CXE02_TF
! 12.01.30 (BTD) * correct typo in WRITE(8,FMT=9020)...
!                  CMDFFT -> CMDSOL
!                * add CMDSOL to argument list (here and in DDSCAT.f90)
! 12.12.23 (BTD) v6 / v7.3.0
!                * minor changes to formatting -- put space between
!                  end of number and ")"
! 13.03.22 (BTD) * changed 
!                  CFLAVG*15  -> CFLAVG*12
!                  CFLE1*15   -> CFLE1*18
!                  CFLE2*15   -> CFLE2*18
!                  CFLB1*15   -> CFLEB1*19
!                  CFLB2*15   -> CFLEB2*19
!                  CFLPOL1*17 -> CFLPOL1*20 
!                  CFLPOL2*17 -> CFLPOL2*20 
!                * added NORICHAR to argument list, and to argument list
!                  in call to NAMER
! end history
! Copyright (C) 1996,1998,2003,2004,2005,2006,2007,2008,2011,2012,2013
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

!*** diagnostic
!      write(0,*)'writesca_v6 ckpt 1, myid=',myid
!      write(0,*)'   cflsca=',cflsca
!      write(0,*)'   cflavg=',cflavg
!***
      PI=4._WP*ATAN(1._WP)

! compute DAEFF=d/aeff and DPHYS=d (physical units)

      DAEFF=(4._WP*PI/(3._WP*REAL(NAT0,KIND=WP)))**(1._WP/3._WP)
      DPHYS=WAVE*AK1/(2._WP*PI)

! If IORI > 0, then write out scattering properties for this specific
! target orientation.

      IF(CMDFRM=='LFRAME'.OR.CMDFRM=='lframe')THEN
         CFRAME='Lab Frame'
      ELSEIF(CMDFRM=='TFRAME'.OR.CMDFRM=='tframe')THEN
         CFRAME='Target Frame'
      ELSE
         CFRAME='????????'
      ENDIF

      IF(IORI>0)THEN

!*** Compute G=<cos(theta)> and G2=<cos^2(theta)>

         DO J=1,IORTH
            IF(JPBC==0)THEN
               G(J)=QSCAG(1,J)/QSCAT(J)
               G2(J)=QSCAG2(J)/QSCAT(J)
            ELSE

! QSCAT and QSCAG are not available when JPBC > 0

               G(J)=0._WP
               G2(J)=0._WP
            ENDIF
         ENDDO

!*** diagnostic
!         write(0,*)'writesca_v6 ckpt 2, myid=',myid
!***
         IF(IWRKSC==1)THEN

            CALL NAMER(IWAV-1,IRAD-1,IORI-1,NORICHAR,CFLPOL1,CFLPOL2, &
                       CFLSCA,CFLAVG,CFLE1,CFLE2,CFLEB1,CFLEB2)
!*** diagnostic
!            write(0,*)'writesca_v6 ckpt 2.1, myid=',myid,' cflsca=',cflsca
!            write(0,*)'                   itheta=',itheta
!            write(0,*)'                   iphi=',iphi
!***
            OPEN(UNIT=8,FILE=CFLSCA,STATUS='UNKNOWN')

            WRITE(8,FMT=9020)CSTAMP,CDESCR,CALPHA,CMDSOL,CSHAPE,NAT0,DAEFF, &
                             DPHYS
            WRITE(8,FMT=9030)(XMIN*DPHYS),(XMAX*DPHYS),(YMIN*DPHYS), &
                             (YMAX*DPHYS),(ZMIN*DPHYS),(ZMAX*DPHYS)

! following code to be enabled for noncubic treatment:
!     &         ,DX

            WRITE(8,FMT=9032)AEFF,WAVE,XX,NAMBIENT
            DO J=1,NCOMP
               MKD=DAEFF*SQRT(REAL(CXRFR(J)*CONJG(CXRFR(J))))*XX
               WRITE(8,FMT=9031)CXRFR(J),CXEPS(J),MKD,J
            ENDDO
            WRITE(8,FMT=9033)TOL,A1_TF,A2_TF
            IF(JPBC==0)WRITE(8,FMT=9034)NAVG
            WRITE(8,FMT=9035)AK_TF,CXE01_TF,CXE02_TF

! calculate A1_LF = target axis A1 in Lab Frame
!           A2_LF = target axis A2 in Lab Frame

            A1_LF(1)=COS(THETA(ITHETA))
            A1_LF(2)=SIN(THETA(ITHETA))*COS(PHI(IPHI))
            A1_LF(3)=SIN(THETA(ITHETA))*SIN(PHI(IPHI))
            A2_LF(1)=-SIN(THETA(ITHETA))*COS(BETA(IBETA))
            A2_LF(2)=COS(THETA(ITHETA))*COS(BETA(IBETA))*COS(PHI(IPHI))- &
                   SIN(BETA(IBETA))*SIN(PHI(IPHI))
            A2_LF(3)=COS(THETA(ITHETA))*COS(BETA(IBETA))*SIN(PHI(IPHI))+ &
                   SIN(BETA(IBETA))*COS(PHI(IPHI))
            WRITE(8,FMT=9036)A1_LF,A2_LF
            RVAR1=AK1
            RVAR2=0._WP
            RVAR3=0._WP
            WRITE(8,FMT=9037)RVAR1,RVAR2,RVAR3,CXE01_LF,CXE02_LF
            WRITE(8,FMT=9040)BETAD,THETAD,PHID
            IF(JPBC==0)THEN
               WRITE(8,FMT=9041)ETASCA
               WRITE(8,FMT=9050)QEXT(1),QABS(1),QSCAT(1),G(1),G2(1), &
                                QBKSCA(1),QPHA(1)
            ELSEIF(JPBC==1.OR.JPBC==2)THEN
               WRITE(8,9350)
            ELSEIF(JPBC==3)THEN

! abs.coeff = 1-R-T
!           = [Qabs*pi/cos(theta_i)]*(d/L_y)*(d/L_z)*(3*N/4*pi)^{2/3}

               ABSCO(1)=PI*(AK1/ABS(AK_TF(1)))*QABS(1)*                    &
                        (0.75_WP*REAL(NAT0)/PI)**(2._WP/3._WP)/(PYDDX*PZDDX)
               WRITE(8,9550)ABSCO(1),ITNUM(1),MXITER
            ENDIF
         ENDIF

         IF(IORTH==1)THEN

! Assign "missing value" value if there is only one polarization

            DO I=1,18
               QAV(I)=-999._WP
            ENDDO
            QEXT(2)=-999._WP
            QABS(2)=-999._WP
            QSCAT(2)=-999._WP
            QBKSCA(2)=-999._WP
            QPHA(2)=-999._WP
            G(2)=-999._WP
            G2(2)=-999._WP
            DO J=1,3
               QSCAG(J,2)=-999._WP
            ENDDO
            QSCAG2(2)=-999._WP
         ENDIF

         IF(IORTH==2)THEN
            IF(JPBC==0)THEN
               QAV(1)=.5_WP*(QEXT(1)+QEXT(2))
               QAV(2)=.5_WP*(QABS(1)+QABS(2))
               QAV(3)=.5_WP*(QSCAT(1)+QSCAT(2))
               QAV(4)=.5_WP*(QSCAG(1,1)+QSCAG(1,2))/QAV(3)
               QAV(5)=.5_WP*(QSCAG2(1)+QSCAG2(2))/QAV(3)
               QAV(6)=.5_WP*(QBKSCA(1)+QBKSCA(2))
               QAV(7)=.5_WP*(QPHA(1)+QPHA(2))
               QAV(8)=QEXT(1)-QEXT(2)
               QAV(9)=QPHA(1)-QPHA(2)
               DO J=1,3
                  QAV(9+J)=.5_WP*(QSCAG(J,1)+QSCAG(J,2))
               ENDDO
               IF(CMDTRQ=='DOTORQ')THEN
                  DO J=1,3
                     QAV(12+J)=.5_WP*(QTRQAB(J,1)+QTRQAB(J,2))
                     QAV(15+J)=.5_WP*(QTRQSC(J,1)+QTRQSC(J,2))
                  ENDDO
               ENDIF
            ELSE
               QAV(1)=.5_WP*(QEXT(1)+QEXT(2))
               QAV(2)=.5_WP*(QABS(1)+QABS(2))
               DO J=3,6
                  QAV(J)=0._WP
               ENDDO
               QAV(7)=.5_WP*(QPHA(1)+QPHA(2))
               QAV(8)=QEXT(1)-QEXT(2)
               QAV(9)=QPHA(1)-QPHA(2)
               DO J=10,18
                  QAV(J)=0._WP
               ENDDO
            ENDIF
            IF(IWRKSC==1)THEN
               IF(JPBC==0)THEN
                  WRITE(8,FMT=9055)QEXT(2),QABS(2),QSCAT(2),G(2),G2(2), &
                                   QBKSCA(2),QPHA(2),(QAV(J),J=1,9)
               ELSEIF(JPBC==1.OR.JPBC==2)THEN
                  WRITE(8,FMT=9355)
               ELSEIF(JPBC==3)THEN
                  ABSCO(2)=PI*(AK1/ABS(AK_TF(1)))*QABS(2)*                    &
                           (0.75_WP*REAL(NAT0)/PI)**(2._WP/3._WP)/(PYDDX*PZDDX)
                  QAV(2)=0.5_WP*(ABSCO(1)+ABSCO(2))
                  WRITE(8,9551)ABSCO(2),ITNUM(2),MXITER,QAV(2)
               ENDIF
            ENDIF
         ENDIF

         IF(IWRKSC==1)THEN
            IF(JPBC==0)THEN
               WRITE(8,FMT=9150)(QSCAG(J,1),J=1,3),ITNUM(1),MXITER,NAVG
               IF(IORTH==2)THEN
                  WRITE(8,FMT=9155)(QSCAG(J,2),J=1,3),ITNUM(2),MXITER,NAVG, &
                                   (QAV(J),J=10,12)
               ENDIF
            ENDIF
         ENDIF

         IF(CMDTRQ=='DOTORQ')THEN
            IF(IWRKSC==1)THEN
               WRITE(8,FMT=9250)(QTRQAB(J,1),J=1,3),(QTRQSC(J,1),J=1,3)
               IF(IORTH==2)WRITE(8,FMT=9255)(QTRQAB(J,2),J=1,3), &
                           (QTRQSC(J,2),J=1,3),(QAV(J),J=13,18)
            ENDIF
         ELSE

! add missing value code (set to -999.)

            DO I=1,3
               DO J=1,2
                  QTRQAB(I,J)=-999._WP
                  QTRQSC(I,J)=-999._WP
               ENDDO
            ENDDO

         ENDIF

         IF(IORTH==1)THEN

            IF(IWRKSC==1)WRITE(8,FMT=9060)

            DO ND=1,NSCAT

! convert scattering angles to degrees

               PHIND=DEGRAD*PHIN(ND)
               THETND=DEGRAD*THETAN(ND)
               RVAR1=REAL(CONJG(CXF11(ND))*CXF11(ND))
               RVAR2=REAL(CONJG(CXF21(ND))*CXF21(ND))
               CXVAR1=CONJG(CXF11(ND))*CXF21(ND)
               IF(IWRKSC==1)WRITE(8,FMT=9070)ND,THETND,PHIND,RVAR1,RVAR2, &
                                             CXVAR1
            ENDDO

         ELSEIF(IORTH==2)THEN

            IF(IWRKSC==1)THEN

! write file header information

               WRITE(8,FMT=9080)CFRAME
               IF(JPBC==0.OR.JPBC==3)THEN
                  IF(NSMELTS==1)THEN
                     WRITE(8,FMT=9081)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
                  ELSEIF(NSMELTS==2)THEN
                     WRITE(8,FMT=9082)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
                  ELSEIF(NSMELTS==3)THEN
                     WRITE(8,FMT=9083)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
                  ELSEIF(NSMELTS==4)THEN
                     WRITE(8,FMT=9084)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
                  ELSEIF(NSMELTS==5)THEN
                     WRITE(8,FMT=9085)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
                  ELSEIF(NSMELTS==6)THEN
                     WRITE(8,FMT=9086)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
                  ELSEIF(NSMELTS==7)THEN
                     WRITE(8,FMT=9087)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
                  ELSEIF(NSMELTS==8)THEN
                     WRITE(8,FMT=9088)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
                  ELSEIF(NSMELTS==9)THEN
                     WRITE(8,FMT=9089)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
                  ENDIF
               ELSEIF(JPBC==1.OR.JPBC==2)THEN
                  IF(NSMELTS==1)THEN
                     WRITE(8,FMT=9181)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
                  ELSEIF(NSMELTS==2)THEN
                     WRITE(8,FMT=9182)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
                  ELSEIF(NSMELTS==3)THEN
                     WRITE(8,FMT=9183)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
                  ELSEIF(NSMELTS==4)THEN
                     WRITE(8,FMT=9184)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
                  ELSEIF(NSMELTS==5)THEN
                     WRITE(8,FMT=9185)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
                  ELSEIF(NSMELTS==6)THEN
                     WRITE(8,FMT=9186)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
                  ELSEIF(NSMELTS==7)THEN
                     WRITE(8,FMT=9187)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
                  ELSEIF(NSMELTS==8)THEN
                     WRITE(8,FMT=9188)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
                  ELSEIF(NSMELTS==9)THEN
                     WRITE(8,FMT=9189)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
                  ENDIF
               ENDIF

               DO ND=1,NSCAT
  
! convert scattering angles to degrees

                  PHIND=DEGRAD*PHIN(ND)
                  THETND=DEGRAD*THETAN(ND)

!*** diagnostic
!                  write(0,*)'writesca_v6 ckpt 3, nd=',nd
!                  write(0,*)'            phin(nd),thetan(nd)=',
!     &                      phin(nd),thetan(nd)
!***

! Compute PPOL = degree of linear polarization of scattered light
! for incident unpolarized light
!*** diagnostic
!                  write(0,*)'writesca_v6 ckpt 4, sm(nd,1,1)=',sm(nd,1,1)
!***
                  PPOL=SQRT(SM(2,1,ND)**2+SM(3,1,ND)**2)/SM(1,1,ND)

                  WRITE(8,FMT=9090)THETND,PHIND,PPOL,                     &
                                   (SM(SMIND1(J),SMIND2(J),ND),J=1,NSMELTS)
               ENDDO
            ENDIF
         ENDIF

         IF(IWRKSC==1)CLOSE(8)

      ELSE

! This module writes out results from orientational averaging.
! (arrive here when called with IORI=0)

         OPEN(UNIT=8,FILE=CFLAVG,STATUS='UNKNOWN')
         WRITE(8,FMT=9020)CSTAMP,CDESCR,CALPHA,CMDSOL,CSHAPE,NAT0,DAEFF, &
                          DPHYS

! following code to be enabled for noncubic treatment:
!         WRITE(8,FMT=9020)CSTAMP,CDESCR,CMDSOL,CALPHA,CSHAPE,NAT0,DX

         WRITE(8,FMT=9032)AEFF,WAVE,XX,NAMBIENT
         DO J=1,NCOMP
            MKD =(4._WP*PI/(3._WP*NAT0))**(1._WP/3._WP)* &
                 SQRT(REAL(CXRFR(J)*CONJG(CXRFR(J))))*XX
            WRITE(8,FMT=9031)CXRFR(J),CXEPS(J),MKD,J
         ENDDO
         WRITE(8,FMT=9033)TOL,A1_TF,A2_TF
         IF(JPBC==0)WRITE(8,FMT=9034)NAVG
         RVAR1=AK1
         RVAR2=0._WP
         RVAR3=0._WP
         WRITE(8,FMT=9037)RVAR1,RVAR2,RVAR3,CXE01_LF,CXE02_LF
         WRITE(8,FMT=9042)BETMID,BETMXD,NBETA,THTMID,THTMXD,NTHETA, &
                          PHIMID,PHIMXD,NPHI
         IF(JPBC==0)WRITE(8,FMT=9041)ETASCA
         WRITE(8,FMT=9043)NORI,IORTH
         IF(JPBC==0)THEN
            RVAR1=QSCGSUM(1,1)/QSCSUM(1)
            RVAR2=QSCG2SUM(1)/QSCSUM(1)
            WRITE(8,FMT=9050)QEXSUM(1),QABSUM(1),QSCSUM(1),RVAR1,RVAR2, &
                             QBKSUM(1),QPHSUM(1)
         ELSE
            QSCSUM(1)=QEXSUM(1)-QABSUM(1)
            WRITE(8,FMT=9051)QEXSUM(1),QABSUM(1),QSCSUM(1)
         ENDIF
         IF(IORTH==1)THEN
            IF(JPBC==0)THEN
               WRITE(8,FMT=9150)(QSCGSUM(J,1),J=1,3),ITNUM(1),MXITER,NAVG
               IF(CMDTRQ=='DOTORQ')WRITE(8,FMT=9250)(QTRQABSUM(J,1),J=1,3), &
                                   (QTRQSCSUM(J,1),J=1,3)
            ELSE
               WRITE(8,FMT=9151)
            ENDIF
            WRITE(8,FMT=9060)
            DO ND=1,NSCAT

! Convert scattering angles to degrees

               PHIND=DEGRAD*PHIN(ND)
               THETND=DEGRAD*THETAN(ND)
               WRITE(8,FMT=9070)ND,THETND,PHIND,S1111(ND),S2121(ND), &
                                CX1121(ND)
            ENDDO

! Write orientationally-averaged Q values (except Qpha) to 'qtable':

!*** diagnostic
!            write(0,*)'writesca_v6 ckpt 5, about to open qtable...'
!***
            OPEN(UNIT=10,FILE='gammatable',STATUS='OLD') 
            !Name changed to 'gammatable' by SMC 8.5.13 following NWB 7/12/12
            DO ND=1,ILIN10
               READ(10,FMT=*)
            ENDDO
            IF(JPBC==0)THEN
               RVAR1=QSCGSUM(1,1)/QSCSUM(1)
               RVAR2=QSCG2SUM(1)/QSCSUM(1)
            ENDIF
!*** diagnostic
!            write(0,*)'writesca_v6 ckpt 6, about to write to qtable...'
!***
            !Inserted by SMC 8.5.13 following NWB 7/12/12
            !***Compute energy loss NWB 7/12/12
            ENERGY = 1.2398419301E+0_WP/WAVE
            !Output modified by NWB 7/12/12 to only have wavelength and Gamma (as QEXSUM(1))
            WRITE(10,FMT=9056)WAVE,ENERGY,QEXSUM(1)
            ILIN10=ILIN10+1
            CLOSE(10)

            !Original Code commented out by SMC 8.5.13 following NWB 7/12/12
            !IF(JPBC==0)THEN
            !   WRITE(10,FMT=9056)AEFF,WAVE,QEXSUM(1),QABSUM(1),QSCSUM(1), &
            !                     RVAR1,RVAR2,QBKSUM(1),NAVG
            !ELSE
            !   WRITE(10,FMT=9057)AEFF,WAVE,QEXSUM(1),QABSUM(1),QSCSUM(1)
            !ENDIF
            !ILIN10=ILIN10+1
            !CLOSE(10)
            !OPEN(UNIT=12,FILE='qtable2',STATUS='OLD')
            !DO ND=1,ILIN12
            !   READ(12,FMT=*)
            !ENDDO
            !WRITE(12,FMT=9057)AEFF,WAVE,QPHSUM(1)
            !ILIN12=ILIN12+1
            !CLOSE(12)

         ELSEIF(IORTH==2)THEN
            IF(JPBC==0)THEN
               RVAR1=QSCGSUM(1,2)/QSCSUM(2)
               QAV(1)=.5_WP*(QEXSUM(1)+QEXSUM(2))
               QAV(2)=.5_WP*(QABSUM(1)+QABSUM(2))
               QAV(3)=.5_WP*(QSCSUM(1)+QSCSUM(2))
               QAV(4)=.5_WP*(QSCGSUM(1,1)+QSCGSUM(1,2))/QAV(3)
               QAV(5)=.5_WP*(QSCG2SUM(1)+QSCG2SUM(2))/QAV(3)
               QAV(6)=.5_WP*(QBKSUM(1)+QBKSUM(2))
               QAV(7)=.5_WP*(QPHSUM(1)+QPHSUM(2))
               QAV(8)=QEXSUM(1)-QEXSUM(2)
               QAV(9)=QPHSUM(1)-QPHSUM(2)
               DO J=1,3
                  QAV(9+J)=.5_WP*(QSCGSUM(J,1)+QSCGSUM(J,2))
               ENDDO
               IF(CMDTRQ=='DOTORQ')THEN
                  DO J=1,3
                     QAV(12+J)=.5_WP*(QTRQABSUM(J,1)+QTRQABSUM(J,2))
                     QAV(15+J)=.5_WP*(QTRQSCSUM(J,1)+QTRQSCSUM(J,2))
                  ENDDO
               ENDIF
               RVAR2=.5_WP*(QSCG2SUM(1)+QSCG2SUM(2))/QAV(3)
            ELSE
               QSCSUM(2)=QEXSUM(2)-QABSUM(2)
               RVAR1=0._WP
               RVAR2=0._WP
               QAV(1)=.5_WP*(QEXSUM(1)+QEXSUM(2))
               QAV(2)=.5_WP*(QABSUM(1)+QABSUM(2))
               QAV(3)=.5_WP*(QSCSUM(1)+QSCSUM(2))
               DO J=4,6
                  QAV(J)=0._WP
               ENDDO
               QAV(7)=.5_WP*(QPHSUM(1)+QPHSUM(2))
               QAV(8)=QEXSUM(1)-QEXSUM(2)
               QAV(9)=QPHSUM(1)-QPHSUM(2)
               DO J=10,12
                  QAV(J)=0._WP
               ENDDO
            ENDIF
            IF(JPBC==0)THEN
               WRITE(8,FMT=9055)QEXSUM(2),QABSUM(2),QSCSUM(2),RVAR1,RVAR2, &
                                QBKSUM(2),QPHSUM(2),(QAV(J),J=1,9)
               WRITE(8,FMT=9150)(QSCGSUM(J,1),J=1,3),ITNUM(1),MXITER,NAVG
               WRITE(8,FMT=9155)(QSCGSUM(J,2),J=1,3),ITNUM(2),MXITER,NAVG, &
                                (QAV(J),J=10,12)
               IF(CMDTRQ=='DOTORQ')THEN
                  WRITE(8,FMT=9250)(QTRQABSUM(J,1),J=1,3), &
                                   (QTRQSCSUM(J,1),J=1,3)
                  WRITE(8,FMT=9255)(QTRQABSUM(J,2),J=1,3),               &
                                   (QTRQSCSUM(J,2),J=1,3),(QAV(J),J=13,18)
               ENDIF
            ELSE
               WRITE(8,FMT=9054)QEXSUM(2),QABSUM(2),QSCSUM(2), &
                                (QAV(J),J=1,3),QAV(8)
            ENDIF

! Write orientationally-averaged Q values (except Qpha) to 'qtable':

            OPEN(UNIT=10,FILE='qtable',STATUS='OLD')
            DO ND=1,ILIN10
               READ(10,FMT=*)
            ENDDO
            WRITE(10,FMT=9056)AEFF,WAVE,QAV(1),QAV(2),QAV(3),QAV(4), &
                              QAV(5),QAV(6),NAVG
            ILIN10=ILIN10+1
            CLOSE(10)
            OPEN(UNIT=12,FILE='qtable2',STATUS='OLD')
            DO ND=1,ILIN12
               READ(12,FMT=*)
            ENDDO
            WRITE(12,FMT=9057)AEFF,WAVE,QAV(7),QAV(8),QAV(9)
            ILIN12=ILIN12+1
            CLOSE(12)

            WRITE(8,FMT=9080)CFRAME
            IF(NSMELTS==1)THEN
               WRITE(8,FMT=9081)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
            ELSEIF(NSMELTS==2)THEN
               WRITE(8,FMT=9082)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
            ELSEIF(NSMELTS==3)THEN
               WRITE(8,FMT=9083)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
            ELSEIF(NSMELTS==4)THEN
               WRITE(8,FMT=9084)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
            ELSEIF(NSMELTS==5)THEN
               WRITE(8,FMT=9085)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
            ELSEIF(NSMELTS==6)THEN
               WRITE(8,FMT=9086)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
            ELSEIF(NSMELTS==7)THEN
               WRITE(8,FMT=9087)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
            ELSEIF(NSMELTS==8)THEN
               WRITE(8,FMT=9088)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
            ELSEIF(NSMELTS==9)THEN
               WRITE(8,FMT=9089)(SMIND1(J),SMIND2(J),J=1,NSMELTS)
            ENDIF

            DO ND=1,NSCAT

! Convert scattering angles to degrees

               PHIND=DEGRAD*PHIN(ND)
               THETND=DEGRAD*THETAN(ND)

! Compute PPOL = degree of linear polarization of scattered light
! for incident unpolarized light (Bohren & Huffman p. 382)

               PPOL=SQRT(SMORI(2,1,ND)**2+SMORI(3,1,ND)**2)/SMORI(1,1,ND)

               WRITE(8,FMT=9090)THETND,PHIND,PPOL,                        &
                                (SMORI(SMIND1(J),SMIND2(J),ND),J=1,NSMELTS)
            ENDDO

         ENDIF

         CLOSE(8)
      ENDIF

! write FORTRAN unformatted file

      IF(CBINFLAG/='NOTBIN')THEN

         CALL WRITEBIN(CBINFLAG,CBINFILE,IOBIN,CSTAMP,CDESCR,CMDFFT,CALPHA,  &
                       CSHAPE,CMDTRQ,MXWAV,MXRAD,MXTHET,MXBETA,MXPHI,MXN3,   &
                       MXNAT,MXSCA,MXTIMERS,IORTH,NX,NY,NZ,NAT0,NAT,NAT3,    &
                       NAVG,NSCAT,NCOMP,NORI,NWAV,NRAD,NTHETA,NPHI,NBETA,    &
                       NTIMERS,MXNX,MXNY,MXNZ,MXCOMP,BETMID,BETMXD,THTMID,   &
                       THTMXD,PHIMID,PHIMXD,TOL,DX,WAVEA,AEFFA,THETA,BETA,   &
                       PHI,A1_TF,A2_TF,ICOMP,IXYZ0,IWAV,IRAD,IORI,AEFF,WAVE, &
                       XX,AK1,BETAD,THETAD,PHID,AK_TF,CXE01_TF,CXE02_TF,     &
                       CXRFR,CXEPS,QEXT,QABS,QSCAT,G,QBKSCA,QPHA,QAV,QTRQAB, &
                       QTRQSC,SM,PHIN,THETAN,TIMERS)

      ENDIF

      RETURN

9020  FORMAT(' DDSCAT --- ',A,/,' TARGET ---',A,/,' ',A,             &
         ' --- DDA method',/,' ',A,                         &
         ' --- CCG method',/,                   &
         ' ',A,' --- shape ',/,I8,5X,'= NAT0 = number of dipoles',/, &
         F12.8,' = d/aeff for this target [d=dipole spacing]',/,     &
         F12.6,' = d (physical units)')
9030  FORMAT(                                                               &
         '----- physical extent of target volume in Target Frame ------',/, &
         2F14.6,' = xmin,xmax (physical units)',/,                          &
         2F14.6,' = ymin,ymax (physical units)',/,                          &
         2F14.6,' = zmin,zmax (physical units)')
! following code to be enabled for noncubic treatment
! 9030 FORMAT(' DDSCAT --- ',A,/,' TARGET ---',A,/,' ',A,
!     &       ' --- method of solution ',/,' ',A,
!     &       ' --- prescription for polarizabilies',/,' ',A,
!     &       ' --- shape ',/,I7,' = NAT0 = number of dipoles',
!     &       3F8.5,' = normalized lattice spacings dx,dy,dz')
9031  FORMAT(                                                  &
         'n= (',F7.4,' , ',F7.4,'),  eps.= (',F8.4,' , ',F7.4, &
         ')  |m|kd=',0P,F8.4,' for subs.',I2)
9032  FORMAT(                                                          &
         '  AEFF=',F14.6,' =',' effective radius (physical units)',/,  &
         '  WAVE=',F14.6,' = wavelength (in vacuo, physical units)',/, &
         'K*AEFF=',F14.6,' = 2*pi*aeff/lambda',/,                      &
         'NAMBIENT=',F12.6,' = refractive index of ambient medium')
9033  FORMAT(                                                       &
         '   TOL=',1P,E10.3,' = error tolerance for CCG method',/,  &
         0P,'(',F8.5,2F9.5,' ) = target axis A1 in Target Frame',/, &
         '(',F8.5,2F9.5,' ) = target axis A2 in Target Frame')
9034  FORMAT('  NAVG=',I6,' = (theta,phi) values used in comp. of Qsca,g')
9035  FORMAT(                                                  &
         '(',F8.5,2F9.5,' ) = k vector (latt. units) in TF',/, &
         3('(',F8.5,',',F8.5,' )'),'=inc.pol.vec. 1 in TF',/,  &
         3('(',F8.5,',',F8.5,' )'),'=inc.pol.vec. 2 in TF')
9036  FORMAT(                                                 &
         '(',F8.5,2F9.5,' ) = target axis A1 in Lab Frame',/, &
         '(',F8.5,2F9.5,' ) = target axis A2 in Lab Frame')
9037  FORMAT(                                                         &
         '(',F8.5,2F9.5,' ) = k vector (latt. units) in Lab Frame',/, &
         3('(',F8.5,',',F8.5,' )'),'=inc.pol.vec. 1 in LF',/,         &
         3('(',F8.5,',',F8.5,' )'),'=inc.pol.vec. 2 in LF')
9040  FORMAT(                                                &
         ' BETA =',F7.3,' = rotation of target around A1',/, &
         ' THETA=',F7.3,' = angle between A1 and k',/,       &
         '  PHI =',F7.3,' = rotation of A1 around k')
9041  FORMAT(                                                            &
         F7.4,' = ETASCA = param. controlling # of scatt. dirs used to', &
         ' calculate <cos> etc.')
9042  FORMAT(                                                 &
         2F8.3,' = beta_min, beta_max ;  NBETA =',I2,/,2F8.3, &
         ' = theta_min, theta_max; NTHETA=',I2,/,2F8.3,       &
         ' = phi_min, phi_max   ;   NPHI =',I2,/,F7.4)
9043  FORMAT(                                                    &
         ' Results averaged over ',I4,' target orientations',/,  &
         '                   and ',I4,' incident polarizations')
9050  FORMAT(                                                         &
         10X,'Qext',8X,'Qabs',8X,'Qsca',6X,'g(1)=<cos>',2X,'<cos^2>', &
         5X,'Qbk',7X,'Qpha',/,1X,'JO=1:',1P,4E12.4,2E11.4,E12.4)
9051  FORMAT(10X,'Qext',7X,'Qabs',7X,'Qsca',/,1X,'JO=1:',1P,3E12.4)
9054  FORMAT(1X,'JO=2:',1P,3E12.4,/,1X,'mean:',3E12.4,/,1X,'Qpol= ',E11.4)
9055  FORMAT(                                                            &
         1X,'JO=2:',1P,3E12.4,E12.4,2E11.4,E12.4,/,1X,'mean:',         &
         4E12.4,2E11.4,E12.4,/,1X,'Qpol=',E12.4,50X,'dQpha=',E12.4)
9056  FORMAT(1P,E10.4,4E11.4,E12.4,E11.4,E11.4,I6)
9057  FORMAT(1P,E10.4,E11.4,3E12.4)
9060  FORMAT(                                                            &
         '**** Selected scattering directions ',                         &
         ' [note: incident pol state 1 is FIXED!]',/,                    &
         ' ND THETA   PHI <|f11|^2> <|f21|^2> Re<f11*f21>',' Im<f11*f21>')
9061  FORMAT(                                                            &
         '**** Selected scattering directions ',                         &
         ' [note: incident pol state 1 is FIXED!]',/,                    &
         ' ND orderM zeta <|f11|^2> <|f21|^2> Re<f11*f21>',' Im<f11*f21>')
9063  FORMAT(                                                              &
         '**** Selected scattering directions ',                           &
         ' [note: incident pol state 1 is FIXED!]',/,                      &
         ' ND orderM orderN <|f11|^2> <|f21|^2> Re<f11*f21>',' Im<f11*f21>')
9070  FORMAT(I3,2F6.1,1P,2E10.3,2E11.3)
9080  FORMAT(                                                    &
         12X,'Mueller matrix elements for selected scattering ', &
         'directions in ',A12)
9081  FORMAT(' theta    phi    Pol.    S_',I1,I1)
9082  FORMAT(' theta    phi    Pol.    S_',I1,I1,8X,'S_',I1,I1)
9083  FORMAT(' theta    phi    Pol.    S_',I1,I1,8X,'S_',I1,I1,8X,'S_',I1,I1)
9084  FORMAT(' theta    phi    Pol.    S_',I1,I1,8X,'S_',I1,I1,8X,'S_',I1,I1, &
        7X,'S_',I1,I1)
9085  FORMAT(' theta    phi    Pol.    S_',I1,I1,8X,'S_',I1,I1,8X,'S_',I1,I1, &
        7X,'S_',I1,I1,7X,'S_',I1,I1)
9086  FORMAT(' theta    phi    Pol.    S_',I1,I1,8X,'S_',I1,I1,8X,'S_',I1,I1, &
        7X,'S_',I1,I1,7X,'S_',I1,I1,7X,'S_',I1,I1)
9087  FORMAT(' theta    phi    Pol.    S_',I1,I1,8X,'S_',I1,I1,8X,'S_',I1,I1, &
        7X,'S_',I1,I1,7X,'S_',I1,I1,7X,'S_',I1,I1,7X,'S_',I1,I1)
9088  FORMAT(' theta    phi    Pol.    S_',I1,I1,8X,'S_',I1,I1,8X,'S_',I1,I1, &
        7X,'S_',I1,I1,7X,'S_',I1,I1,7X,'S_',I1,I1,7X,'S_',I1,I1,7X,'S_',I1,I1)
9089  FORMAT(' theta    phi    Pol.    S_',I1,I1,8X,'S_',I1,I1,8X,'S_',I1,I1, &
        7X,'S_',I1,I1,7X,'S_',I1,I1,7X,'S_',I1,I1,7X,'S_',I1,I1,7X,'S_',I1,   &
        I1,7X,'S_',I1,I1)
!-----------------------------------------------------------------------
9181  FORMAT(' alpha   zeta    Pol.    S_',I1,I1)
9182  FORMAT(' alpha   zeta    Pol.    S_',I1,I1,8X,'S_',I1,I1)
9183  FORMAT(' alpha   zeta    Pol.    S_',I1,I1,8X,'S_',I1,I1,8X,'S_',I1,I1)
9184  FORMAT(' alpha   zeta    Pol.    S_',I1,I1,8X,'S_',I1,I1,8X,'S_',I1,I1, &
        7X,'S_',I1,I1)
9185  FORMAT(' alpha   zeta    Pol.    S_',I1,I1,8X,'S_',I1,I1,8X,'S_',I1,I1, &
        7X,'S_',I1,I1,7X,'S_',I1,I1,7X,'S_',I1,I1)
9186  FORMAT(' alpha   zeta    Pol.    S_',I1,I1,8X,'S_',I1,I1,8X,'S_',I1,I1, &
        7X,'S_',I1,I1,7X,'S_',I1,I1,7X,'S_',I1,I1)
9187  FORMAT(' alpha   zeta    Pol.    S_',I1,I1,8X,'S_',I1,I1,8X,'S_',I1,I1, &
        7X,'S_',I1,I1,7X,'S_',I1,I1,7X,'S_',I1,I1,7X,'S_',I1,I1)
9188  FORMAT(' alpha   zeta    Pol.    S_',I1,I1,8X,'S_',I1,I1,8X,'S_',I1,I1, &
        7X,'S_',I1,I1,7X,'S_',I1,I1,7X,'S_',I1,I1,7X,'S_',I1,I1,7X,'S_',I1,I1)
9189  FORMAT(' alpha   zeta    Pol.    S_',I1,I1,8X,'S_',I1,I1,8X,'S_',I1,I1, &
        7X,'S_',I1,I1,7X,'S_',I1,I1,7X,'S_',I1,I1,7X,'S_',I1,I1,7X,'S_',I1,   &
        I1,7X,'S_',I1,I1)
!-----------------------------------------------------------------------
9090  FORMAT(F6.2,F7.2,F9.5,1P,2E12.4,7E11.3)
9150  FORMAT(9X,'Qsca*g(1)   Qsca*g(2)   Qsca*g(3)   iter  mxiter',2X, &
        'Nsca',/,1X,'JO=1:',1P,3E12.4,3I7)
9151  FORMAT(43X,' iter  mxiter')
9155  FORMAT(1X,'JO=2:',1P,3E12.4,3I7,/,1X,'mean:',1P,3E12.4)
9250  FORMAT(8X,'Qtrqab(1)  Qtrqab(2)  Qtrqab(3)  ', &
        'Qtrqsc(1)  Qtrqsc(2)  Qtrqsc(3)',/,1X,'JO=1:',1P,6E11.3)
9255  FORMAT(1X,'JO=2:',1P,6E11.3,/,1X,'mean:',1P,6E11.3)
9350  FORMAT('*** need to figure out what we want to print for total cross section/unit length for JPBC=1,2 ***')
9355  FORMAT('*** need to figure out what we want to print for total cross section/unit length for JPBC=1,2 ***')
9550  FORMAT(' absorption coeff.  iter  mxiter',/,1X,'JO=1: ',1PE11.4,2I6)
9551  FORMAT(1X,'JO=2: ',1PE11.4,2I6,/,1X,'mean: ',1PE11.4)
    END SUBROUTINE WRITESCA
