    SUBROUTINE WRITEPOL(NRWORD,MXNAT,NX,NY,NZ,NAT0,IANISO,  &
                        ICOMP,IXYZ0,PYD,PZD,AKD,DX,X0,WAVE, &
                        BETADF,THETADF,PHIDF,CXE0,CXADIA,   &
                        CXAOFF,CXPOL,CFLPOL)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

! arguments:

      INTEGER :: IANISO,MXNAT,NAT0,NRWORD
      INTEGER*2 ::    &
         ICOMP(MXNAT,3)
      INTEGER ::      &
         IXYZ0(NAT0,3)
      INTEGER :: NX,NY,NZ
      REAL(WP) :: PYD,PZD,WAVE
      REAL(WP) ::       &
         AKD(3),        &
         BETADF(NAT0),  &
         DX(3),         &
         PHIDF(NAT0),   &
         THETADF(NAT0), &
         X0(3)
      COMPLEX(WP) ::      &
         CXADIA(MXNAT,3), &
         CXAOFF(MXNAT,3), &
         CXE0(3),         &
         CXPOL(NAT0,3)
      CHARACTER :: CFLPOL*17

! local variables:

      INTEGER :: J,NA

!***********************************************************************
! subroutine WRITEPOL
! given:
!    NRWORD   = length (bytes) of real word [integer]
!    MXNAT    = dimensioning parameter
!    NX       = (1/d)*(maximum extent of target in x-direction) (TF)
!    NY       = (1/d)*(maximum extent of target in y-direction) (TF)
!    NZ       = (1/d)*(maximum extent of target in z-direction) (TF)
!    NAT0     = number of occupied dipole sites
!    IANISO   = 0 for isotropic materials
!               1 for anisotropic materials with optical axes aligned
!                 with Target Frame (i.e., DF=TF)
!               2 for general anisotropic material, with arbitrary
!                 orientation of DF at each site relative to TF
!    ICOMP(IA,JD)=composition of dipole # IA
!                 for principal axis JD=1,2,3 in "dielectric frame"
!                 where dielectric tensor is diagonalized
!    IXYZ0(IA,JD)=x/d,y/d,z/d for dipole # IA, JD=1-3
!                 where (x,y,z) are in Target Frame
!    PYD      = periodicity/d of target replication in y_TF direction
!             = 0 if no replication (single target)
!    PZD      = periodicity/d of target replication in z_TF direction
!             = 0 if no replication (single target)
!    AKD(1-3)=k_x*d,k_y*d,k_z*d where (k_x,k_y,k_z)=k in ambient medium
!             in Target Frame
!    DX(1-3)  = lattice spacing/d in x,y,z directions, with
!               DX(1)*DX(2)*DX(3)=1
!               normally have DX(1)=DX(2)=DX(3)=1
!    X0(1-3)  = location/(DX*d) in TF of lattice site IX=0,IY=0,IZ=0
!    WAVE     = wavelength in ambient medium (physical units)
!    CXE0(1-3)=E_x,E_y,E_z for incident wave in ambient medium
!              in Target Frame
!    CXADIA(3i-2) =alphainv_xx at location i in TF
!          (3i-1) =alphainv_yy at location i in TF
!          (3i)   =alphainv_zz at location i in TF
!          where alphainv = inverse of (complex polarizability tensor)/d
!    CXAOFF(3i-2) =alpha_yz=alphainv_zy at location i in TF
!          (3i-1) =alpha_zx=alphainv_xz at location i in TF
!          (3i)   =alpha_xy=alphainv_yx at location i in TF
!          where alphainv = inverse of (complex polarizability matrix/d^
!    CXPOL(i)        =P_x at location i in TF
!         (NAT0+i)   =P_y at location i
!         (2*NAT0+i) =P_z at location i
!    CFLPOL = filename

! writes NAT0,IXYZ,CXPOL to unformatted file CFLPOL

! Note: IXYZ0, CXPOL, CXADIA, CXAOFF are "reduced" arrays
!       IXYZ0 is produced by subroutine EXTEND
!       CXPOL is reduced by subroutine REDUCE called from GETFML

!       CXADIA is reduced by REDUCE called from GETFLM
!       CXAOFF is reduced by REDUCE called from GETFLM

!       where "reduced" arrays are limited to values for occupied sites

! Note: if IANISO=0 or 1: do not store BETADF,THETADF,PHIDF
!       if IANISO=2       do store BETADF,THETADF,PHIDF

! B.T. Draine, Princeton Univ. Observatory, 2006.04.08
! history
! 06.04.08 (BTD) first written
! 06.04.10 (BTD) revised
! 06.09.21 (BTD) add PYD,PZD, and DX to argument list and to output file
! 07.06.21 (BTD) modify for use with DDSCAT 7.0.2:
!                * add X0 to argument list
!                * add X0 to WRITE statement
! 07.09.11 (BTD) changed IXYZ from INTEGER*2 to INTEGER
! 08.01.12 (BTD) add WAVE to argument list
!                add WAVE to WRITE statement
! 08.01.13 (BTD) add NRWORD to argument list and to WRITE statement
! 08.01.17 (BTD) * add IANISO to argument list, 
!                * write IANISO to file
!                * use IANISO to determine whether BETADF,THETADF,PHIDF
!                  need to be written out
!                * addition cleanup
! 08.03.15 (BTD) v7.0.5
!                * corrected error in dimensioning
!                  IXYZ0(MXNAT,3) -> IXYZ0(NAT0,3)
!                  BETADF(MXNAT) -> BETADF(NAT0)
!                  PHIDF(MXNAT) -> PHIDF(NAT0)
!                  THETADF(MXNAT) -> THETADF(NAT0)
! 08.05.21 (BTD) * CXPOL is reduced in DDSCAT before calling WRITEPOL:
!                  change dimensioning 
!                  CXPOL(MXNAT,3) -> CXPOL(NAT0,3)
!
! end history

! Copyright (C) 2006,2007,2008
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

      NA=NAT0
      OPEN(UNIT=22,FILE=CFLPOL,FORM='UNFORMATTED')
      WRITE(22)NRWORD,NAT0,IANISO
      WRITE(22)NX,NY,NZ,(IXYZ0(J,1),J=1,NA),(IXYZ0(J,2),J=1,NA),            &
               (IXYZ0(J,3),J=1,NA),PYD,PZD,AKD,DX,X0,WAVE,CXE0,             &
               (CXADIA(J,1),J=1,NA),(CXADIA(J,2),J=1,NA),                   &
               (CXADIA(J,3),J=1,NA),(CXAOFF(J,1),J=1,NA),                   &
               (CXAOFF(J,2),J=1,NA),(CXAOFF(J,3),J=1,NA),                   &
               (CXPOL(J,1),J=1,NA),(CXPOL(J,2),J=1,NA),(CXPOL(J,3),J=1,NA), &
               (ICOMP(J,1),J=1,NA),(ICOMP(J,2),J=1,NA),(ICOMP(J,3),J=1,NA)
      IF(IANISO==2)THEN
         WRITE(22)(BETADF(J),J=1,NA),(THETADF(J),J=1,NA),(PHIDF(J),J=1,NA)
      ENDIF
      CLOSE(22)
      RETURN
    END SUBROUTINE WRITEPOL
