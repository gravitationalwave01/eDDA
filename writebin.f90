    SUBROUTINE WRITEBIN(CBINFLAG,CBINFILE,IOBIN,CSTAMP,CDESCR,CMDFFT,CALPHA, &
        CSHAPE,CMDTRQ,MXWAV,MXRAD,MXTHET,MXBETA,MXPHI,MXN3,MXNAT,MXSCA,      &
        MXTIMERS,IORTH,NX,NY,NZ,NAT0,NAT,NAT3,NAVG,NSCAT,NCOMP,NORI,NWAV,    &
        NRAD,NTHETA,NPHI,NBETA,NTIMERS,MXNX,MXNY,MXNZ,MXCOMP,BETMID,BETMXD,  &
        THTMID,THTMXD,PHIMID,PHIMXD,TOL,DX,WAVEA,AEFFA,THETA,BETA,PHI,A1,A2, &
        ICOMP,IXYZ0,IWAV,IRAD,IORI,AEFF,WAVE,XX,AK1,BETAD,THETAD,PHID,AKR,   &
        CXE01R,CXE02R,CXRFR,CXEPS,QEXT,QABS,QSCAT,G,QBKSCA,QPHA,QAV,QTRQAB,  &
        QTRQSC,SM,PHIN,THETAN,TIMERS)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE
!     .. Scalar Arguments ..
      REAL (WP) :: AEFF,AK1,BETAD,BETMID,BETMXD,PHID,PHIMID,PHIMXD, &
        THETAD,THTMID,THTMXD,TOL,WAVE,XX
      INTEGER :: IOBIN,IORI,IORTH,IRAD,IWAV,MXBETA,MXCOMP,MXN3,MXNAT, &
        MXNX,MXNY,MXNZ,MXPHI,MXRAD,MXSCA,MXTHET,MXTIMERS,MXWAV,NAT,   &
        NAT0,NAT3,NAVG,NBETA,NCOMP,NORI,NPHI,NRAD,NSCAT,NTHETA,       &
        NTIMERS,NWAV,NX,NY,NZ
      CHARACTER :: CALPHA*6,CBINFLAG*6,CMDFFT*6,CMDTRQ*6,CSHAPE*9, &
        CBINFILE*14,CSTAMP*26,CDESCR*67
!     .. Array Arguments ..
      COMPLEX (WP) :: CXE01R(3),CXE02R(3),CXEPS(MXCOMP),CXRFR(MXCOMP)
      REAL (WP) :: A1(3),A2(3),AEFFA(MXRAD),AKR(3),BETA(MXBETA),DX(3), &
        G(2),PHI(MXPHI),PHIN(MXSCA),QABS(2),QAV(17),QBKSCA(2),QEXT(2), &
        QPHA(2),QSCAT(2),QTRQAB(3,2),QTRQSC(3,2),SM(MXSCA,4,4),        &
        THETA(MXTHET),THETAN(MXSCA),TIMERS(MXTIMERS),WAVEA(MXWAV)
      INTEGER*2 :: ICOMP(MXN3)
      INTEGER :: IXYZ0(MXNAT,3)
!     .. Local Scalars ..
      INTEGER :: I,IENTRY,IS,J
!     .. External Subroutines ..
      EXTERNAL WRIMSG
!     .. Data statements ..
      DATA IENTRY/0/
      SAVE IENTRY

!***********************************************************************
! Subroutine  "writebin" writes FORTRAN unformatted binary file
! to file cbinfile ('dd.bin')

! Original version created by P.J.Flatau, University of California, SD
! History:
! 96.11.15 (PJF): First version.
! 98.01.11 (BTD): Added DX(3) to argument list, and write DX to
!                 binary output upon entry (DX doesn't change).
! 99.06.30 (BTD): Added SAVE IENTRY statement.
! 03.10.23 (BTD): Removed ICTHM and IPHIM from argument list
!                 Removed ICTHM and IPHIM from WRITE statement
!                 Added NAVG to argument list.
!                 Note: manner in which data is written could be
!                 improved along lines followed in WRITENET.
!                 However, does not appear to be demand, so such
!                 changes are postponed indefinitely.
! 07.09.11 (BTD): changed IXYZ0 from INTEGER*2 to INTEGER
! 13.05.02 (BTD): corrected error found by Vasyl Choliy
!                 IXYZ0(I,J) should only be written for I=1-NAT0,
!                 not I=1-NAT
! end history
! Copyright (C) 1996,1998,1999,2003,2007,2013 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

      IF(IENTRY==0)THEN
         CALL WRIMSG('WRITEBIN','First call to binary write')
         IENTRY = 1
! General run information:
         WRITE(IOBIN)CSTAMP,CDESCR,CMDFFT,CALPHA,CSHAPE,CMDTRQ
! mostly dimensions:
         WRITE(IOBIN)IORTH,NX,NY,NZ,NAT0,NAT,NAT3,NSCAT,NCOMP,NORI, &
            NWAV,NRAD,NTHETA,NPHI,NBETA,NTIMERS
! maximum dimension (not really needed)
         WRITE(IOBIN)MXNX,MXNY,MXNZ,MXCOMP,MXSCA
         WRITE(IOBIN)BETMID,BETMXD,THTMID,THTMXD,PHIMID,PHIMXD
! dump scattering angles
         WRITE(IOBIN)(PHIN(I),I=1,NSCAT),(THETAN(I),I=1,NSCAT)
! misc run information
         WRITE(IOBIN) TOL,NAVG
         WRITE(IOBIN) (WAVEA(I),I=1,NWAV),(AEFFA(I),I=1,NRAD),                 &
            (THETA(I),I=1,NTHETA),(BETA(I),I=1,NBETA),(PHI(I),I=1,NPHI),       &
            (A1(I),I=1,3),(A2(I),I=1,3),(ICOMP(I),I=1,NAT3),                   &
            (IXYZ0(I,1),I=1,NAT0),(IXYZ0(I,2),I=1,NAT0),(IXYZ0(I,3),I=1,NAT0), &
            (DX(I),I=1,3)

! endif  "if(cbinflag.eq.'USEBIN') then"
      ENDIF

! binary write of  for particular "iwav, irad, itheta, ibeta iphi"
! Let us try to keep only  one version of the binary file. Thus, in
! case of IORTH=1 or 'DOTORQ' some values are "missing"
! flatau -->draine  do we need smori for "ppol" do we need cxf_ij ?

      IF(IORI/=0.AND.CBINFLAG=='ALLBIN')THEN
         WRITE(IOBIN)IWAV,IRAD,IORI,(TIMERS(I),I=1,NTIMERS)
         WRITE(IOBIN)AEFF,WAVE,XX,AK1,BETAD,THETAD,PHID,AKR,CXE01R, &
            CXE02R,(CXRFR(I),I=1,NCOMP),(CXEPS(I),I=1,NCOMP)
         WRITE(IOBIN)QEXT,QABS,QSCAT,G,QBKSCA,QPHA,QAV,QTRQAB,QTRQSC
! write sm as 16 independent matrices (IDL handles 7 subscripts only)
         DO I=1,4
            DO J=1,4
               WRITE(IOBIN)(SM(IS,I,J),IS=1,NSCAT)
            ENDDO
         ENDDO
! endif "if(iori. ne. 0 .and. cbinflag.eq.'ALLCDF') then"
      ENDIF

      IF(IORI==0)THEN

! BTD 98.05.01 add CONTINUE statement in lieu of actual code

         CONTINUE

! iori=0 case
!flatau ---> note. I need to add code for orientational averages
! but may be I can recalculate it from sm for iwav, irad, itheta,
! ---> dump weights

      ENDIF

      RETURN

    END SUBROUTINE WRITEBIN
