! mpi_fake.f contains dummy subroutines that take the place
! of
!    SUBROUTINE COLSUM
!    SUBROUTINE SHARE0
!    SUBROUTINE SHARE1
!    SUBROUTINE SHARE2
!    SUBROUTINE MPI_INIT
!    SUBROUTINE MPI_COMM_RANK
!    SUBROUTINE MPI_COMM_SIZE
!    SUBROUTINE MPI_FINALIZE
!    SUBROUTINE MPI_BARRIER
! which are called when MPI is being used.
! These dummy routines do nothing when MPI is not being used.

    SUBROUTINE COLSUM(IORTH,MYID,MXSCA,NSCAT,QSCSUM,QSCSUM_1,QABSUM,     &
                      QABSUM_1,QEXSUM,QEXSUM_1,QBKSUM,QBKSUM_1,QPHSUM,   &
                      QPHSUM_1,QSCG2SUM,QSCG2SUM_1,QSCGSUM,QSCGSUM_1,    &
                      QTRQABSUM,QTRQABSUM_1,QTRQSCSUM,QTRQSCSUM_1,S1111, &
                      S1111_1,S2121,S2121_1,CX1121,CX1121_1,SMORI,SMORI_1)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE
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
! local variables
      INTEGER J,K,ND

      DO ND=1,IORTH
         QSCSUM(ND)=QSCSUM(ND)+QSCSUM_1(ND)
         QABSUM(ND)=QABSUM(ND)+QABSUM_1(ND)
         QEXSUM(ND)=QEXSUM(ND)+QEXSUM_1(ND)
         QBKSUM(ND)=QBKSUM(ND)+QBKSUM_1(ND)
         QPHSUM(ND)=QPHSUM(ND)+QPHSUM_1(ND)
         QSCG2SUM(ND)=QSCG2SUM(ND)+QSCG2SUM_1(ND)
      ENDDO
      DO J=1,IORTH
         DO ND=1,3
            QSCGSUM(ND,J)=QSCGSUM(ND,J)+QSCGSUM_1(ND,J)
            QTRQABSUM(ND,J)=QTRQABSUM(ND,J)+QTRQABSUM_1(ND,J)
            QTRQSCSUM(ND,J)=QTRQSCSUM(ND,J)+QTRQSCSUM_1(ND,J)
         ENDDO
      ENDDO
      DO ND=1,NSCAT
         CX1121(ND)=CX1121(ND)+CX1121_1(ND)
         S1111(ND)=S1111(ND)+S1111_1(ND)
         S2121(ND)=S2121(ND)+S2121_1(ND)
      ENDDO
      DO ND=1,NSCAT
         DO K=1,4
            DO J=1,4
               SMORI(J,K,ND)=SMORI(J,K,ND)+SMORI_1(J,K,ND)
            ENDDO
         ENDDO
      ENDDO
         
      RETURN
    END SUBROUTINE COLSUM

    SUBROUTINE SHARE0(LACE,LAXI,LCLM,LGI,LPI,LQI,                     &
                      MXN3,MXNAT,MXNX,MXNY,MXNZ,MXPBC,MXCXSC,MYID,NAT0)

      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE
! arguments
      INTEGER :: LACE,LAXI,LCLM,LGI,LPI,LQI,                    &
                 MXN3,MXNAT,MXNX,MXNY,MXNZ,MXPBC,MXCXSC,MYID,NAT0
      RETURN
    END SUBROUTINE SHARE0
    
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

! arguments

      INTEGER :: MXNAT,MXN3,MXBETA,MXCOMP,MXNX,MXNY,MXNZ,MXPBC, &
                 MXPHI,MXRAD,MXSCA,MXTHET,MXWAV
      CHARACTER :: CALPHA*6,CBINFLAG*(*),CDESCR*67,CMDFFT*6, &
                   CMDFRM*6,CSHAPE*9,CMDSOL*6,CMDTRQ*6
      CHARACTER(60) :: &
         CFLEPS(MXCOMP)
      INTEGER :: IANISO,IDVERR,IDVOUT,ILIN10,ILIN12,INIT,IOBIN,IORTH,IORI0,  &
                 IPBC,IRAD0,IWAV0,IWRKSC,IWRPOL,JPBC,MYID,NAT,NAT0,NAT3,     &
                 NBETA,NBETH,NCOMP,NORI,NPHI,NRAD,NSCAT,NSMELTS,NTHETA,NWAV, &
                 NX,NY,NZ
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
         SHPAR(10),          &
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
!-----------------------------------------------------------------------
      RETURN
    END SUBROUTINE SHARE1

!=======================================================================

    SUBROUTINE SHARE2(MXCOMP,MXWAVT,CXEPS,E1A,E2A,WVA)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE
! arguments:
      INTEGER :: MXCOMP,MXWAVT,IERR
      REAL(WP) ::     &
         E1A(MXWAVT), &
         E2A(MXWAVT), &
         WVA(MXWAVT)
      COMPLEX(WP) :: CXEPS(MXCOMP)
!-----------------------------------------------------------------------
      RETURN
    END SUBROUTINE SHARE2

!=======================================================================

    SUBROUTINE MPI_INIT(IERR)
      IMPLICIT NONE
! arguments:
      INTEGER :: IERR
!-----------------------------------------------------------------------
      IERR=0
      RETURN
    END SUBROUTINE MPI_INIT

!=======================================================================

    SUBROUTINE MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      IMPLICIT NONE
! arguments:
      INTEGER :: IERR,MPI_COMM_WORLD,MYID
!-----------------------------------------------------------------------
      MYID=0
      RETURN
    END SUBROUTINE MPI_COMM_RANK

!=======================================================================

    SUBROUTINE MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
      IMPLICIT NONE
! arguments:
      INTEGER :: IERR,MPI_COMM_WORLD,NUMPROCS
!-----------------------------------------------------------------------
      NUMPROCS=1
      RETURN
    END SUBROUTINE MPI_COMM_SIZE

!=======================================================================

    SUBROUTINE MPI_FINALIZE(IERR)
      IMPLICIT NONE
! arguments:
      INTEGER :: IERR
!-----------------------------------------------------------------------
      IERR=0
      RETURN
    END SUBROUTINE MPI_FINALIZE

!=======================================================================

    SUBROUTINE MPI_BARRIER(MPI_COMM_WORLD,IERR)
      IMPLICIT NONE
! arguments:
      INTEGER :: IERR,MPI_COMM_WORLD
!-----------------------------------------------------------------------
      IERR=0
      RETURN
    END SUBROUTINE MPI_BARRIER
