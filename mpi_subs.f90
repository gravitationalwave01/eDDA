!=======================================================================
! Module mpi_subs.f
! This module contains three main subroutines 
!    SUBROUTINE COLSUM
!    SUBROUTINE SHARE0
!    SUBROUTINE SHARE1
!    SUBROUTINE SHARE2
! and auxiliary routines
!    SUBROUTINE MPI_BCAST_CHAR
!    SUBROUTINE MPI_BCAST_INT
!    SUBROUTINE MPI_BCAST_INT2
!    SUBROUTINE MPI_BCAST_REAL
!    SUBROUTINE MPI_BCAST_CPLX

! COLSUM, SHARE1, and SHARE2 are used by DDSCAT to manage information 
! in all the various processes when MPI is being used. 

! MPI_BCAST_CHAR, MPI_BCAST_INT, MPI_BCAST_INT2, MPI_BCAST_REAL, and
! MPI_BCAST_CPLX are used here to provide "jackets" for MPI_BCAST for
! different variable types, to eliminate warning messages from fortran
! compilers.
!
! COLSUM, SHARE1, and SHARE2
! were originally written by Matthew Collinge,
! Princeton University Observatory, 2002.11.15

! Code has been modified to provide "jackets" for MPI_BCAST for
! different variable types, to eliminate warning messages from fortran
! compilers.

! history: 
! 02.11.15 (BTD) Working version developed by Matthew Collinge
! 03.07.13 (BTD) cosmetic changes during preparation of official
!                release of version 6.0
! 04.04.09 (BTD) modifications to introduce jacket routines
!                MPI_BCAST_CHAR 
!                MPI_BCAST_INT
!                MPI_BCAST_INT2
!                MPI_BCAST_REAL
!                MPI_BCAST_CPLX
! 08.05.01 (BTD) added MYID to arg list of SHARE0
! 08.05.09 (BTD) added LACE,LAXI,LCLM,LGI,LPI,LQI,LSC0 to arg list of SHARE0
! end history
!
! Copyright (C) 2002,2003,2004,2007,2008 B.T. Draine and M. Collinge
!
! This code is covered by the GNU General Public License
!=======================================================================

    SUBROUTINE COLSUM(IORTH,MYID,MXSCA,NSCAT,QSCSUM,QSCSUM_1,QABSUM,     &
                      QABSUM_1,QEXSUM,QEXSUM_1,QBKSUM,QBKSUM_1,QPHSUM,   &
                      QPHSUM_1,QSCG2SUM,QSCG2SUM_1,QSCGSUM,QSCGSUM_1,    &
                      QTRQABSUM,QTRQABSUM_1,QTRQSCSUM,QTRQSCSUM_1,S1111, &
                      S1111_1,S2121,S2121_1,CX1121,CX1121_1,SMORI,SMORI_1)

      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

!-----------------------------------------------------------------------

      INCLUDE 'mpif.h'

!-----------------------------------------------------------------------

! arguments:

      INTEGER :: IORTH,MXSCA,NSCAT,MYID
      REAL(WP) ::              &
         QABSUM(IORTH),        &
         QABSUM_1(IORTH),      &
         QBKSUM(IORTH),        &
         QBKSUM_1(IORTH),      &
         QEXSUM(IORTH),        &
         QEXSUM_1(IORTH),      &
         QPHSUM(IORTH),        &
         QPHSUM_1(IORTH),      &
         QSCG2SUM(IORTH),      &
         QSCG2SUM_1(IORTH),    &
         QSCGSUM(3,IORTH),     &
         QSCGSUM_1(3,IORTH),   &
         QSCSUM(IORTH),        &
         QSCSUM_1(IORTH),      &
         QTRQABSUM(3,IORTH),   &
         QTRQABSUM_1(3,IORTH), &
         QTRQSCSUM(3,IORTH),   &
         QTRQSCSUM_1(3,IORTH), &
         S1111(MXSCA),         &
         S1111_1(MXSCA),       &
         S2121(MXSCA),         &
         S2121_1(MXSCA),       &
         SMORI(4,4,MXSCA),     &
         SMORI_1(4,4,MXSCA)
      COMPLEX(WP) ::    &
         CX1121(MXSCA), &
         CX1121_1(MXSCA)

! local variables:

      LOGICAL :: SINGLE
      INTEGER :: IERR,NDAT

!=======================================================================
! Subroutine COLSUM
! Given:
!    IORTH                   = 1 for single incident polarization
!                              2 for both incident polarizations
!    MYID                    = process rank
!    MXSCA                   = max number of scattering directions
!    NSCAT                   = actual number of scattering directions
!    QSCSUM_1(1-IORTH)       = sum of Q_sca for orientations done by MYID
!    QABSUM_1(1-IORTH)       = sum of Q_abs "
!    QEXSUM_1(1-IORTH)       = sum of Q_ext "
!    QBKSUM_1(1-IORTH)       = sum of Q_bk  "
!    QPHSUM_1(1-IORTH)       = sum of Q_ph  "
!    QSCGSUM_1(1-3,1-IORTH)  = sum of Q_sca*<cos> for orientation done by MYID
!    QSCG2SUM_1(1-IORTH)     = sum of Q_sca*<cos^2> "
!    QTRQABSUM_1(1-3,1-IORTH)= sum of Q_trq,abs for orientations done by MYID
!    QTRQSCSUM_1(1-3,1-IORTH)= sum of Q_trq,sca "
!    S1111_1(1-NSCAT)        = sum of S1111 for orientations done by MYID
!    S2121_1(1-NSCAT)        = sum of S2121 for orientations done by MYID
!    SMORI_1(1-NSCAT,1-4,1-4)= sum of Mueller matrix for orientations done
!                              by MYID 
!    QSCSUM(1-IORTH)         = running sum of Q_sca
!    QABSUM(1-IORTH)         = running sum of Q_abs
!    QEXSUM(1-IORTH)         = running sum of Q_ext
!    QBKSUM(1-IORTH)         = running sum of Q_bk
!    QPHSUM(1-IORTH)         = running sum of Q_ph
!    QSCGSUM(1-3,1-IORTH)    = running sum of Q_sca*<cos>
!    QSCG2SUM(1-IORTH)       = running sum of Q_sca*<cos^2>
!    QTRQABSUM(1-3,1-IORTH)  = running sum of Q_trq,abs
!    QTRQSCSUM(1-3,1-IORTH)  = running sum of Q_trq,sca
!    S1111(1-NSCAT)          = running sum of S1111
!    S2121(1-NSCAT)          = running sum of S2121
!    SMORI(1-NSCAT,1-4,1-4)  = running sum of Mueller matrix for orientations
!                              done by MYID

! Returns
!    QSCSUM(1-IORTH)         = updated Q_sca sum over orientations
!    QABSUM(1-IORTH)         = updated Q_abs sum over orientations
!    QEXSUM(1-IORTH)         = updated Q_ext sum over orientations
!    QBKSUM(1-IORTH)         = updated Q_bk sum over orientations
!    QPHSUM(1-IORTH)         = updated Q_pha sum over orientations
!    QSCGSUM(1-3,1-IORTH)    = updated Q_sca*<cos> sum over orientations
!    QSCG2SUM(1-IORTH)       = updated Q_sca*<cos^2> sum over orientations
!    QTRQABSUM(1-3,1-IORTH)  = updated Q_trq,abs sum over orientations
!    QTRQSCSUM(1-3,1-IORTH)  = updated Q_trq,sca sum over orientations
!    S1111(1-NSCAT)          = updated S1111 sum over orientations
!    S2121(1-NSCAT)          = updated S2121 sum over orientations
!    SMORI(1-NSCAT,1-4,1-4)  = updated Mueller matrix sum over orientations
!
! Originally written by Matthew Collinge
! history
! 03.10.23 (BTD) Added variable QSCG2SUM(1-IORTH)
! 08.05.10 (BTD) ver7.0.5
!                * major revision: add to argument list
!                   CX1121_1
!                   QABSUM_1
!                   QBKSUM_1
!                   QEXTSUM_1
!                   QPHSUM_1
!                   QSCAGSUM_1
!                   QSCG2SUM_1
!                   QSCASUM_1
!                   QTRQABSUM_1
!                   QTRQSCSUM_1
!                   S1111_1
!                   S2121_1
!                   SMORI_1
! 08.05.11 (BTD) rewrite to transfer NDAT values per call to MPI_REDUCE
!                where NDAT=IORTH or IORTH*3 or NSCAT or 4*4*NSCAT
! 08.05.29 (BTD) v7.0.6
!                * change declarations
!                  SMORI(MXSCA,4,4) -> SMORI(4,4,MXSCA)
!                  in order to accomplish update of SMORI using SMORI_1
!                  with a single call to MPI_REDUCE
! end history
!=======================================================================
!*** diagnostic
!      write(0,*)'COLSUM ckpt 1, myid=',myid
!***
      IF(WP==KIND(0.E0))THEN
         SINGLE=.TRUE.
      ELSEIF(WP==KIND(0.D0))THEN
         SINGLE=.FALSE.
      ELSE
         WRITE(0,*)'fatal error in mpi_subs: unable to determine precision'
         STOP
      ENDIF
!*** diagnostic
!      write(0,*)'COLSUM ckpt 2, myid=',myid
!***
!*** diagnostic
!      write(0,*)'COLSUM ckpt 3, myid=',myid
      IF(SINGLE)THEN
!*** diagnostic
!         write(0,*)'COLSUM ckpt 4, myid=',myid
!***
         CALL MPI_REDUCE(QSCSUM_1(1),QSCSUM(1),IORTH, &
                         MPI_REAL,MPI_SUM,0,          &
                         MPI_COMM_WORLD,IERR)
      ELSE
!*** diagnostic
!         write(0,*)'COLSUM ckpt 5, myid=',myid
!***
         CALL MPI_REDUCE(QSCSUM_1(1),QSCSUM(1),IORTH,    &
                         MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                         MPI_COMM_WORLD,IERR)
!*** diagnostic
!         write(0,*)'COLSUM ckpt 5.1, myid=',myid
!***
      ENDIF
!*** diagnostic
!      write(0,*)'COLSUM ckpt 6, myid=',myid
!***
      IF(SINGLE)THEN
!*** diagnostic
!      write(0,*)'COLSUM ckpt 7, myid=',myid
!***
         CALL MPI_REDUCE(QABSUM_1(1),QABSUM(1),IORTH, &
                         MPI_REAL,MPI_SUM,0,          &
                         MPI_COMM_WORLD,IERR)
      ELSE
!*** diagnostic
!      write(0,*)'COLSUM ckpt 8, myid=',myid
!***
         CALL MPI_REDUCE(QABSUM_1(1),QABSUM(1),IORTH,    &
                         MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                         MPI_COMM_WORLD,IERR)
      ENDIF
!*** diagnostic
!      write(0,*)'COLSUM ckpt 9, myid=',myid
!***
      IF(SINGLE)THEN
         CALL MPI_REDUCE(QEXSUM_1(1),QEXSUM(1),IORTH, &
                         MPI_REAL,MPI_SUM,0,          &
                         MPI_COMM_WORLD,IERR)
      ELSE
         CALL MPI_REDUCE(QEXSUM_1(1),QEXSUM(1),IORTH,    &
                         MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                         MPI_COMM_WORLD,IERR)
      ENDIF
!*** diagnostic
!      write(0,*)'COLSUM ckpt 10, myid=',myid
!***
      IF(SINGLE)THEN
         CALL MPI_REDUCE(QBKSUM_1(1),QBKSUM(1),IORTH, &
                         MPI_REAL,MPI_SUM,0,          &
                         MPI_COMM_WORLD,IERR)
      ELSE
         CALL MPI_REDUCE(QBKSUM_1(1),QBKSUM(1),IORTH,    &
                         MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                         MPI_COMM_WORLD,IERR)
      ENDIF
!*** diagnostic
!      write(0,*)'COLSUM ckpt 11, myid=',myid
!***
      IF(SINGLE)THEN
         CALL MPI_REDUCE(QPHSUM_1(1),QPHSUM(1),IORTH, &
                         MPI_REAL,MPI_SUM,0,          &
                         MPI_COMM_WORLD,IERR)
      ELSE
         CALL MPI_REDUCE(QPHSUM_1(1),QPHSUM(1),IORTH,    &
                         MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                         MPI_COMM_WORLD,IERR)
      ENDIF
!*** diagnostic
!      write(0,*)'COLSUM ckpt 12, myid=',myid
!***
      IF(SINGLE)THEN
         CALL MPI_REDUCE(QSCG2SUM_1(1),QSCG2SUM(1),IORTH, &
                         MPI_REAL,MPI_SUM,0,              &
                         MPI_COMM_WORLD,IERR)
      ELSE
         CALL MPI_REDUCE(QSCG2SUM_1(1),QSCG2SUM(1),IORTH, &
                         MPI_DOUBLE_PRECISION,MPI_SUM,0,  &
                         MPI_COMM_WORLD,IERR)
      ENDIF
!*** diagnostic
!      write(0,*)'COLSUM ckpt 13, myid=',myid
!***
      NDAT=3*IORTH
      IF(SINGLE)THEN
         CALL MPI_REDUCE(QSCGSUM_1(1,1),QSCGSUM(1,1),NDAT, &
                         MPI_REAL,MPI_SUM,0,               &
                         MPI_COMM_WORLD,IERR)
      ELSE
         CALL MPI_REDUCE(QSCGSUM_1(1,1),QSCGSUM(1,1),NDAT, &
                         MPI_DOUBLE_PRECISION,MPI_SUM,0,   &
                         MPI_COMM_WORLD,IERR)
      ENDIF
!*** diagnostic
!         write(0,*)'COLSUM ckpt 15, myid=',myid
!***
      IF(SINGLE)THEN
         CALL MPI_REDUCE(QTRQABSUM_1(1,1),QTRQABSUM(1,1),NDAT, &
                         MPI_REAL,MPI_SUM,0,                   &
                         MPI_COMM_WORLD,IERR)
      ELSE
         CALL MPI_REDUCE(QTRQABSUM_1(1,1),QTRQABSUM(1,1),NDAT, &
                         MPI_DOUBLE_PRECISION,MPI_SUM,0,       &
                         MPI_COMM_WORLD,IERR)
      ENDIF
!*** diagnostic
!      write(0,*)'COLSUM ckpt 16, myid=',myid
!***
      IF(SINGLE)THEN
         CALL MPI_REDUCE(QTRQSCSUM_1(1,1),QTRQSCSUM(1,1),NDAT, &
                         MPI_REAL,MPI_SUM,0,                   &
                         MPI_COMM_WORLD,IERR)
      ELSE
         CALL MPI_REDUCE(QTRQSCSUM_1(1,1),QTRQSCSUM(1,1),NDAT, &
                        MPI_DOUBLE_PRECISION,MPI_SUM,0,       &
                        MPI_COMM_WORLD,IERR)
      ENDIF
!*** diagnostic
!      write(0,*)'COLSUM ckpt 17, myid=',myid
!***
      IF(SINGLE)THEN
         CALL MPI_REDUCE(S1111_1(1),S1111(1),NSCAT, &
                         MPI_REAL,MPI_SUM,0,          &
                         MPI_COMM_WORLD,IERR)
      ELSE
         CALL MPI_REDUCE(S1111_1(1),S1111(1),NSCAT,    &
                         MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                         MPI_COMM_WORLD,IERR)
      ENDIF
      IF(SINGLE)THEN
         CALL MPI_REDUCE(S2121_1(1),S2121(1),NSCAT, &
                         MPI_REAL,MPI_SUM,0,          &
                         MPI_COMM_WORLD,IERR)
      ELSE
         CALL MPI_REDUCE(S2121_1(1),S2121(1),NSCAT,    &
                         MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                         MPI_COMM_WORLD,IERR)
      ENDIF
      IF(SINGLE)THEN
         CALL MPI_REDUCE(CX1121_1(1),CX1121(1),NSCAT, &
                         MPI_COMPLEX,MPI_SUM,0,         &
                         MPI_COMM_WORLD,IERR)
      ELSE
         CALL MPI_REDUCE(CX1121_1(1),CX1121(1),NSCAT, &
                         MPI_DOUBLE_COMPLEX,MPI_SUM,0,  &
                         MPI_COMM_WORLD,IERR)
      ENDIF
!*** diagnostic
!   write(0,*)'COLSUM ckpt 19, myid=',myid
!***
      NDAT=4*4*NSCAT

      IF(SINGLE)THEN
         CALL MPI_REDUCE(SMORI_1(1,1,1),SMORI(1,1,1),NDAT, &
                         MPI_REAL,MPI_SUM,0,               &
                         MPI_COMM_WORLD,IERR)
      ELSE
         CALL MPI_REDUCE(SMORI_1(1,1,1),SMORI(1,1,1),NDAT, &
                         MPI_DOUBLE_PRECISION,MPI_SUM,0,   &
                         MPI_COMM_WORLD,IERR)
      ENDIF
!*** diagnostic
!      write(0,*)'COLSUM ckpt 99, myid=',myid
!***
      RETURN
    END SUBROUTINE COLSUM

!=============================================================================

    SUBROUTINE SHARE0(LACE,LAXI,LCLM,LGI,LPI,LQI,                     &
                      MXN3,MXNAT,MXNX,MXNY,MXNZ,MXPBC,MXCXSC,MYID,NAT0)

      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

!-----------------------------------------------------------------------

      INCLUDE 'mpif.h'

!-----------------------------------------------------------------------
! arguments

      INTEGER :: LACE,LAXI,LCLM,LGI,LPI,LQI,                    &
                 MXN3,MXNAT,MXNX,MXNY,MXNZ,MXPBC,MXCXSC,MYID,NAT0

! local variables:

      INTEGER :: IERR

! History:
! 08.05.01 (BTD) v7.0.5 
!                * added MYID to argument list of SHARE0
! end history

      CALL MPI_BCAST_INT(LACE,1,0,IERR)
      CALL MPI_BCAST_INT(LAXI,1,0,IERR)
      CALL MPI_BCAST_INT(LCLM,1,0,IERR)
      CALL MPI_BCAST_INT(LGI,1,0,IERR)
      CALL MPI_BCAST_INT(LPI,1,0,IERR)
      CALL MPI_BCAST_INT(LQI,1,0,IERR)

!*** diagnostic
!      write(0,*)'share0 ckpt 1, myid=',myid,' mxn3=',mxn3
!***
      CALL MPI_BCAST_INT(MXN3,1,0,IERR)
!*** diagnostic
!      write(0,*)'share0 ckpt 2, myid=',myid,' mxn3=',mxn3,' mxnat=',mxnat
!***
      CALL MPI_BCAST_INT(MXNAT,1,0,IERR)
!*** diagnostic
!      write(0,*)'share0 ckpt 3, myid=',myid,' mxnat=',mxnat,' mxnx=',mxnx
!***
      CALL MPI_BCAST_INT(MXNX,1,0,IERR)
!*** diagnostic
!      write(0,*)'share0 ckpt 4, myid=',myid,' mxnx=',mxnx,' mxny=',mxny
!***
      CALL MPI_BCAST_INT(MXNY,1,0,IERR)
!*** diagnostic
!      write(0,*)'share0 ckpt 5, myid=',myid,' mxny=',mxny,' mxnz=',mxnz
!***
      CALL MPI_BCAST_INT(MXNZ,1,0,IERR)
!*** diagnostic
!      write(0,*)'share0 ckpt 6, myid=',myid,' mxnz=',mxnz,' mxpbc=',mxpbc
!***
      CALL MPI_BCAST_INT(MXPBC,1,0,IERR)
!*** diagnostic
!      write(0,*)'share0 ckpt 7, myid=',myid,' mxpbc=',mxpbc,' mxcxsc=',mxcxsc
!***
      CALL MPI_BCAST_INT(MXCXSC,1,0,IERR)
!*** diagnostic
!      write(0,*)'share0 ckpt 8, myid=',myid,' mxcxsc=',mxcxsc,' nat0=',nat0
!***
      CALL MPI_BCAST_INT(NAT0,1,0,IERR)
!*** diagnostic
!      write(0,*)'share0 ckpt 9, myid=',myid,' nat0=',nat0,' return...'
!***
      RETURN

    END SUBROUTINE SHARE0

!=============================================================================

    SUBROUTINE SHARE1(A1,A2,A3,AEFFA,BETA,BETADF,BETMID,BETMXD,CALPHA,       &
                      CBINFLAG,CDESCR,CFLEPS,CMDFFT,CMDFRM,CMDSOL,CMDTRQ,    &
                      CSHAPE,CXE01,CXE02,DAEFF,DX,ENSC,EM1,EM2,ETASCA,GAMMA, &
                      IANISO,ICOMP,IDVERR,IDVOUT,ILIN10,ILIN12,INIT,IOBIN,   &
                      IOCC,IORI0,IORTH,IPBC,IRAD0,IWAV0,IWRKSC,IWRPOL,IXYZ0, &
                      JPBC,MXBETA,MXCOMP,MXN3,MXNAT,MXNX,MXNY,MXNZ,MXPBC,    &
                      MXPHI,MXRAD,MXSCA,MXTHET,MXWAV,MYID,NAT,NAT0,NAT3,     &
                      NBETA,NBETH,NCOMP,NORI,NPHI,NRAD,NSCAT,NSMELTS,NTHETA, &
                      NWAV,NX,NY,NZ,ORDERM,ORDERN,PHI,PHIDF,PHIMID,PHIMXD,   &
                      PHIN,PYD,PYDDX,PZD,PZDDX,SHPAR,SMIND1,SMIND2,THETA,    &
                      THETADF,THETAN,THTMID,THTMXD,TOL,WAVEA,WGTA,WGTB,X0,   &
                      XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN)

      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

!-----------------------------------------------------------------------

      INCLUDE 'mpif.h'

!-----------------------------------------------------------------------
! arguments
! dimensioning information:

      INTEGER :: MXNAT,MXN3,MXBETA,MXCOMP,MXNX,MXNY,MXNZ,MXPBC, &
                 MXPHI,MXRAD,MXSCA,MXTHET,MXWAV,MYID

      CHARACTER :: CBINFLAG*(*),CALPHA*6,CDESCR*67,CMDFFT*6, &
                   CMDFRM*6,CMDSOL*6,CSHAPE*9,CMDTRQ*6

      CHARACTER*60 :: CFLEPS(MXCOMP)

      INTEGER :: IANISO,IDVERR,IDVOUT,ILIN10,ILIN12,INIT,IOBIN,        &
                 IORTH,IORI0,IPBC,IRAD0,IWAV0,IWRKSC,IWRPOL,JPBC,      &
                 NAT,NAT0,NAT3,NBETA,NBETH,NCOMP,NORI,NPHI,NRAD,NSCAT, &
                 NSMELTS,NTHETA,NWAV,NX,NY,NZ

      INTEGER*2 ::    &
         ICOMP(MXN3), &
         IOCC(MXNAT)

      INTEGER ::        &
         IXYZ0(NAT0,3), &
         SMIND1(9),     &
         SMIND2(9)

      REAL(WP) :: BETMID,BETMXD,DAEFF,ETASCA,GAMMA,PHIMID,PHIMXD,PYD,PYDDX, &
                  PZD,PZDDX,THTMID,THTMXD,TOL,XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN

      REAL(WP) ::            &
         A1(3),              &
         A2(3),              &
         A3(3),              &
         AEFFA(MXRAD),       &
         BETA(MXBETA),       &
         BETADF(NAT0),       &
         DX(3),              &
         EM1(3,MXSCA),       &
         EM2(3,MXSCA),       &
         ENSC(3,MXSCA),      &
         ORDERM(MXSCA),      &
         ORDERN(MXSCA),      &
         PHI(MXPHI),         &
         PHIDF(NAT0),        &
         PHIN(MXSCA),        &
         SHPAR(12),          &
         THETA(MXTHET),      &
         THETADF(NAT0),      &
         THETAN(MXSCA),      &
         WAVEA(MXWAV),       &
         WGTA(MXTHET,MXPHI), &
         WGTB(MXBETA),       &
         X0(3)

      COMPLEX(WP) :: &
         CXE01(3),   &
         CXE02(3)

! local variables:
      LOGICAL :: SINGLE
      INTEGER :: IERR
!=======================================================================
! Subroutine SHARE1
! Given:

! Returns:


! This routine uses the MPI data type MPI_INTEGER2, which may 
! not be supported in all implementations of MPI.

! The format of the call to MPI_BCAST is as follows:
!
!     CALL MPI_BCAST(SENDPAR,NUMPAR,MPI_TYPEPAR,SOURCE_PROCESS,
!    &     MPI_COMM_WORLD,IERR)
!
! SENDPAR        = information to communicate with other processes
! NUMPAR         = number of units of information in SENDPAR
! MPI_TYPEPAR    = MPI data type of SENDPAR
! SOURCE_PROCESS = identifier of the process sending the broadcast.
!                  Always 0 (master) in this code.
! MPI_COMM_WORLD = name of MPI communicator
! IERR           = Error code. Not used by this code.
!
! collinge:
! what are the sizes of cbinflag and cnetflag?
! apparently 6, or maaaaybe 7 in fortran convention?
! I am not sure if this means of sending characters actually 
! works or not. Needs to be tested.
!
! Originally written by Matthew Collinge,
! Princeton University Observatory, 2002.11.15

! history: 
! 02.11.15 (BTD) Working version developed by Matthew Collinge
! 03.07.13 (BTD) cosmetic changes during preparation of official
!                release of version 6.0
! 03.10.23 (BTD) removed ICTHM and IPHIM from argument list: no longer
!                used.
!                removed ICTHM and IPHIM from calls to SHARE1
!                removed MPI_BCAST calls for ICTHM and IPHIM
! 03.10.24 (BTD) added ETASCA to argument list
!                added MPI_BCAST call for ETASCA
! 04.09.14 (BTD) added BETADF,PHIDF,THETADF to argument list
!                add MPI_BCAST calls for BETADF,PHIDF,THETADF
! 05.08.04 (BTD) added CMDFRM to argument list
!                added MPI_BCAST_CHAR call for CMDFRM
! 07.01.18 (BTD) added A3(3),IWRPOL,IPBC,JPBC,MXNX,MXNY,MXNZ,MXPBC,
!                ORDERM(MXSCA),ORDERN(MXSCA),
!                PYD,PYDDX,PZD,PZDDX 
!                to argument list of SHARE1, with appropriate
!                MPI_BCAST... calls below
! 07.06.21 (BTD) add X0 to argument list of SHARE1
!                add call to MPI_BCAST_REAL
! 08.01.21 (BTD) add IANISO to argument list of SHARE1
!                add call to MPI_BCAST_INT
!                changed MPI_BCAST_REAL call to share SHPAR(12)
! 08.02.04 (BTD) corrected typos found by Y.Shen
!                add declaration of IXYZ0
! 08.02.10 (BTD) corrected declaration of SHPAR(6)->SHPAR(12) 
! 08.02.12 (BTD) changed CSHAPE*6 to CSHAPE*9, and changed
!                   CALL MPI_BCAST_CHAR(CSHAPE,6,0,IERR)
!                to CALL MPI_BCAST_CHAR(CSHAPE,9,0,IERR)
!                removed CDIEL from argument list
!                removed CALL MPI_BCAST_CHAR(CDIEL,6,0,IERR)
! 08.02.13 (BTD) added MXNAT0 to argument list
!                added CALL_MPI_BCAST_INT calls for
!                * MXBETA
!                * MXNAT
!                * MXNAT0
!                * MXPHI
!                * MXSCA
!                * MXTHET
!                * MXWAV
!                changed dimensioning of IOCC,IXYZ0,BETADF,THETADF,PHIDF from
!                MXNAT to MXNAT0
! 08.02.14 (BTD) * corrected error in dimension of IOCC:
!                  was MXN3, changed to MXNAT
!                also corrected call to MPI_BCAST_INT2
! 08.02.14 (BTD) * changed allocation of IOCC to IOCC(MXNAT) [not MXNAT0]
!                * changed allocation from IXYZ0(MXNAT,3) to IXYZ0(NAT0,3)
!                * changed allocation from BETADF(MXNAT) to BETADF(NAT0)
!                * changed allocation from PHIDF(MXNAT) to PHIDF(NAT0)
!                * changed allocation from THETADF(MXNAT) to THETADF(NAT0)
! 08.02.15 (BTD) * deleted CALL MPI_BCAST_INT(MXNAT0,1,0,IERR)
!                * changed CALL MPI_BCAST_INT(IXYZ0,3*MXNAT0,0,IERR)
!                  to CALL MPI_BCAST_INT(IXYZ0,3*NAT0,0,IERR)
! 08.02.17 (BTD) * add DAEFF to argument list
!                * add CALL_MPI_BCAST_REAL(DAEFF)
! 08.03.11 (BTD) ver7.0.5
!                * add ALPHA to argument list
!                * add CALL_MPI_BCAST_REAL(ALPHA)
! 08.04.19 (BTD) * changed notation: ALPHA -> GAMMA
!                  reordered argument list
! 08.05.01 (BTD) * added MYID to argument list
! 08.07.22 (BTD) * added XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN to argument list
! end history
! Copyright (C) 2002,2003,2004,2005,2007,2008 
!               B.T. Draine, M. Collinge, and P.J. Flatau
!
! This code is covered by the GNU General Public License
!=======================================================================

! character variables:
!*** diagnostic
!      write(0,*)'share1 ckpt 1, myid=',myid
!***
      CALL MPI_BCAST_CHAR(CALPHA,6,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 2, myid=',myid
!***
      CALL MPI_BCAST_CHAR(CBINFLAG,6,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 3, myid=',myid
!***
      CALL MPI_BCAST_CHAR(CDESCR,67,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 4, myid=',myid
!***
      CALL MPI_BCAST_CHAR(CFLEPS,60*MXCOMP,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 5, myid=',myid
!***
      CALL MPI_BCAST_CHAR(CMDFFT,6,0,IERR)
      CALL MPI_BCAST_CHAR(CMDFRM,6,0,IERR)
      CALL MPI_BCAST_CHAR(CMDSOL,6,0,IERR)
      CALL MPI_BCAST_CHAR(CMDTRQ,6,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 10, myid=',myid
!***
      CALL MPI_BCAST_CHAR(CSHAPE,9,0,IERR)


! integer scalars:

      CALL MPI_BCAST_INT(IANISO,1,0,IERR)
      CALL MPI_BCAST_INT(IDVERR,1,0,IERR)
      CALL MPI_BCAST_INT(IDVOUT,1,0,IERR)
      CALL MPI_BCAST_INT(ILIN10,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 15, myid=',myid
!***
      CALL MPI_BCAST_INT(ILIN12,1,0,IERR)
      CALL MPI_BCAST_INT(INIT,1,0,IERR)
      CALL MPI_BCAST_INT(IOBIN,1,0,IERR)
      CALL MPI_BCAST_INT(IORTH,1,0,IERR)
      CALL MPI_BCAST_INT(IORI0,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 20, myid=',myid
!***
      CALL MPI_BCAST_INT(IPBC,1,0,IERR)
      CALL MPI_BCAST_INT(IRAD0,1,0,IERR)
      CALL MPI_BCAST_INT(IWAV0,1,0,IERR)
      CALL MPI_BCAST_INT(IWRKSC,1,0,IERR)
      CALL MPI_BCAST_INT(IWRPOL,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 25, myid=',myid
!***
      CALL MPI_BCAST_INT(JPBC,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 26, myid=',myid
      call mpi_barrier(mpi_comm_world,ierr)
!***
      CALL MPI_BCAST_INT(MXBETA,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 27, myid=',myid
      call mpi_barrier(mpi_comm_world,ierr)
!***
      CALL MPI_BCAST_INT(MXNAT,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 28, myid=',myid
      call mpi_barrier(mpi_comm_world,ierr)
!***
      CALL MPI_BCAST_INT(MXNX,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 29, myid=',myid
      call mpi_barrier(mpi_comm_world,ierr)
!***
      CALL MPI_BCAST_INT(MXNY,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 30, myid=',myid
      call mpi_barrier(mpi_comm_world,ierr)
!***
      CALL MPI_BCAST_INT(MXNZ,1,0,IERR)
      CALL MPI_BCAST_INT(MXPBC,1,0,IERR)
      CALL MPI_BCAST_INT(MXPHI,1,0,IERR)
      CALL MPI_BCAST_INT(MXSCA,1,0,IERR)
      CALL MPI_BCAST_INT(MXTHET,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 35, myid=',myid
      call mpi_barrier(mpi_comm_world,ierr)
!***
      CALL MPI_BCAST_INT(MXWAV,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 36, myid=',myid
      call mpi_barrier(mpi_comm_world,ierr)
!***
      CALL MPI_BCAST_INT(NAT,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 37, myid=',myid
      call mpi_barrier(mpi_comm_world,ierr)
!***
      CALL MPI_BCAST_INT(NAT0,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 38, myid=',myid
      call mpi_barrier(mpi_comm_world,ierr)
!***
      CALL MPI_BCAST_INT(NAT3,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 39, myid=',myid
      call mpi_barrier(mpi_comm_world,ierr)
!***
      CALL MPI_BCAST_INT(NBETA,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 40, myid=',myid
      call mpi_barrier(mpi_comm_world,ierr)
!***
      CALL MPI_BCAST_INT(NBETH,1,0,IERR)
      CALL MPI_BCAST_INT(NCOMP,1,0,IERR)
      CALL MPI_BCAST_INT(NORI,1,0,IERR)
      CALL MPI_BCAST_INT(NPHI,1,0,IERR)
      CALL MPI_BCAST_INT(NRAD,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 45, myid=',myid
!***
      CALL MPI_BCAST_INT(NSCAT,1,0,IERR)
      CALL MPI_BCAST_INT(NSMELTS,1,0,IERR)
      CALL MPI_BCAST_INT(NTHETA,1,0,IERR)
      CALL MPI_BCAST_INT(NWAV,1,0,IERR)
      CALL MPI_BCAST_INT(NX,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 50, myid=',myid
!***
      CALL MPI_BCAST_INT(NY,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 51, myid=',myid,' nz=',nz
      call mpi_barrier(mpi_comm_world,ierr)
!***
      CALL MPI_BCAST_INT(NZ,1,0,IERR)

! integer*2 arrays:
!*** diagnostic
!      write(0,*)'share1 ckpt 52, myid=',myid
!      write(0,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!      write(0,*)'about to call mpi_bcast_int2(icomp,...) with mxn3=',mxn3
      call mpi_barrier(mpi_comm_world,ierr)
!***
      CALL MPI_BCAST_INT2(ICOMP,MXN3,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 53, myid=',myid
!      write(0,*)'about to call mpi_bcast_int2(iocc,...) with mxnat=',mxnat
!***
      CALL MPI_BCAST_INT2(IOCC,MXNAT,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 54, myid=',myid
!      write(0,*)'returned from mpi_bcast_int2'
!***
! integer arrays:

      CALL MPI_BCAST_INT(IXYZ0,3*NAT0,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 55, myid=',myid
!***
      CALL MPI_BCAST_INT(SMIND1,NSMELTS,0,IERR)
      CALL MPI_BCAST_INT(SMIND2,NSMELTS,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 57, myid=',myid
!***

! determine precision

      IF(WP==KIND(0.E0))THEN
         SINGLE=.TRUE.
      ELSEIF(WP==KIND(0.D0))THEN
         SINGLE=.FALSE.
      ELSE
         WRITE(0,*)'fatal error in mpi_subs: unable to determine precision'
         STOP
      ENDIF

! Real scalars:

!*** diagnostic
!      write(0,*)'share1 ckpt 58, myid=',myid
!***
      CALL MPI_BCAST_REAL(SINGLE,BETMID,1,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,BETMXD,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 60, myid=',myid
!***
      CALL MPI_BCAST_REAL(SINGLE,DAEFF,1,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,ETASCA,1,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,GAMMA,1,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,PHIMID,1,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,PHIMXD,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 65, myid=',myid
!***
      CALL MPI_BCAST_REAL(SINGLE,PYD,1,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,PYDDX,1,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,PZD,1,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,PZDDX,1,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,THTMID,1,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 70, myid=',myid
!***
      CALL MPI_BCAST_REAL(SINGLE,THTMXD,1,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,TOL,1,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,XMAX,1,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,XMIN,1,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,YMAX,1,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,YMIN,1,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,ZMAX,1,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,ZMIN,1,0,IERR)

! Real arrays:

      CALL MPI_BCAST_REAL(SINGLE,A1,3,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,A2,3,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,A3,3,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 75, myid=',myid
!***
      CALL MPI_BCAST_REAL(SINGLE,AEFFA,MXRAD,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,BETA,MXBETA,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,BETADF,NAT0,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,DX,3,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,EM1,3*MXSCA,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 80, myid=',myid
!***
      CALL MPI_BCAST_REAL(SINGLE,EM2,3*MXSCA,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,ENSC,3*MXSCA,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,ORDERM,MXSCA,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,ORDERN,MXSCA,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,PHI,MXPHI,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 85, myid=',myid
!***
      CALL MPI_BCAST_REAL(SINGLE,PHIDF,NAT0,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,PHIN,MXSCA,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,SHPAR,12,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,THETA,MXTHET,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,THETADF,NAT0,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 90, myid=',myid
!***
      CALL MPI_BCAST_REAL(SINGLE,THETAN,MXSCA,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,WAVEA,MXWAV,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,WGTA,MXTHET*MXPHI,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,WGTB,MXBETA,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,X0,3,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 95, myid=',myid
!***

! Complex arrays:

      CALL MPI_BCAST_CPLX(SINGLE,CXE01,3,0,IERR)
      CALL MPI_BCAST_CPLX(SINGLE,CXE02,3,0,IERR)
!*** diagnostic
!      write(0,*)'share1 ckpt 97, myid=',myid,' returning from share1'
!***

      RETURN

    END SUBROUTINE SHARE1

!=======================================================================

    SUBROUTINE SHARE2(MXCOMP,MXWAVT,CXEPS,E1A,E2A,WVA)

      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

!-----------------------------------------------------------------------

      INCLUDE 'mpif.h'

!-----------------------------------------------------------------------
! arguments:

      INTEGER :: MXCOMP,MXWAVT,IERR

      REAL(WP) ::     &
         E1A(MXWAVT), &
         E2A(MXWAVT), &
         WVA(MXWAVT)
      COMPLEX(WP) ::  &
         CXEPS(MXCOMP)

! local variable:

      LOGICAL SINGLE

!=======================================================================
! Subroutine SHARE2
! Given:
!    MXCOMP
!    MXWAVT
!    CXEPS
!    E1A
!    E2A
!    WVA

! history:
! Originally written by Matthew Collinge,
! Princeton University Observatory, 2002.11.15

! history: 
! 02.11.15 (BTD) Working version developed by Matthew Collinge
! 03.07.13 (BTD) cosmetic changes during preparation of official
!                release of version 6.0
! end history
!
! Copyright (C) 2002,2003 B.T. Draine and M. Collinge

!=======================================================================

      IF(WP==KIND(0.E0))THEN
         SINGLE=.TRUE.
      ELSEIF(WP==KIND(0.D0))THEN
         SINGLE=.FALSE.
      ELSE
         WRITE(0,*)'fatal error in SHARE2: unable to determine precision'
         STOP
      ENDIF

      CALL MPI_BCAST_REAL(SINGLE,E1A,MXWAVT,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,E2A,MXWAVT,0,IERR)
      CALL MPI_BCAST_REAL(SINGLE,WVA,MXWAVT,0,IERR)
      CALL MPI_BCAST_CPLX(SINGLE,CXEPS,MXCOMP,0,IERR)

      RETURN
    END SUBROUTINE SHARE2

