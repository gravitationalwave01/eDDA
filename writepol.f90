    SUBROUTINE WRITEPOL(NRWORD,MXNAT,NX,NY,NZ,NAT0,IANISO,ICOMP,ISCR1,    &
                        IXYZ0,JPBC,PYD,PZD,GAMMA,NAMBIENT,A1,A2,AEFF,     &
                        AKR,DX,X0,WAVE,BETADF,THETADF,PHIDF,CXE0R,CXADIA, &
                        CXAOFF,CXPOL,CFLPOL)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE
!---------------------- writepol_v2 -------------------------------
! arguments:

      INTEGER :: IANISO,MXNAT,NAT0,NRWORD
      INTEGER*2 ::       &
         ICOMP(MXNAT,3), &
         ISCR1(3*MXNAT)
      INTEGER ::      &
         IXYZ0(NAT0,3)
      INTEGER :: JPBC,NX,NY,NZ
      REAL(WP) :: AEFF,GAMMA,NAMBIENT,PYD,PZD,WAVE
      REAL(WP) ::       &
         A1(3),         &
         A2(3),         &
         AKR(3),        &
         BETADF(NAT0),  &
         DX(3),         &
         PHIDF(NAT0),   &
         THETADF(NAT0), &
         X0(3)
      COMPLEX(WP) ::      &
         CXADIA(MXNAT*3), &
         CXAOFF(MXNAT,3), &
         CXE0R(3),        &
         CXPOL(MXNAT*3)
      CHARACTER :: CFLPOL*17

! local variables:

      INTEGER :: IOPOL,IXMAX,IXMIN,IYMAX,IYMIN,IZMAX,IZMIN,J,JA,NA

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
!    ICOMP(IA,JD)=composition of site # IA
!                 for principal axis JD=1,2,3 in "dielectric frame"
!                 where dielectric tensor is diagonalized
!    IXYZ0(IA,JD)=x/d,y/d,z/d for physical dipole # IA, JD=1-3
!                 where (x,y,z) are in Target Frame
!    JPBC
!    PYD      = periodicity/d of target replication in y_TF direction
!             = 0 if no replication (single target)
!    PZD      = periodicity/d of target replication in z_TF direction
!             = 0 if no replication (single target)
!    GAMMA
!    NAMBIENT = (real) refractive index of ambient medium
!    A1
!    A2
!    AKR(1-3) = k_x*d,k_y*d,k_z*d where (k_x,k_y,k_z)=k in ambient medium
!               in Target Frame
!    DX(1-3)  = lattice spacing/d in x,y,z directions, with
!               DX(1)*DX(2)*DX(3)=1
!               normally have DX(1)=DX(2)=DX(3)=1
!    X0(1-3)  = location/(DX*d) in TF of lattice site IX=0,IY=0,IZ=0
!    WAVE     = wavelength in ambient medium (physical units)
!    CXE0R(1-3)=E_x,E_y,E_z for incident wave in ambient medium
!              in Target Frame
!??? BTD 110819 I don't think this is correct 
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
! 11.07.16 (BTD) v7.2.1
!                * modified write statements to give more flexibilitiy
!                  in reading parts of file before memory allocation
! 11.07.18 (BTD) * add JPBC to argument list
!                * write JPBC to stored file
!                * add target axis vectors A1,A2 to argument list
!                * write A1,A2 to stored file
!                * add GAMMA to argument list
!                * write GAMMA to stored file
! 11.08.01 (BTD) * modified structure of output file for convenience
!                  in subsequent reading
! 11.08.16 (BTD) * disable diagnostic write statements
! 11.08.30 (BTD) v7.2.2
!                * add NAMBIENT to arg list
!                * write NAMBIENT to pol file
! 11.08.31 (BTD) * changed from FORM='UNFORMATTED' to ACCESS='STREAM'
! 12.06.03 (BTD) v2
! end history

! Copyright (C) 2006,2007,2008,2011
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
!*** diagnostic
!      write(0,*)' writepol_v2 ckpt 1'
!      write(0,*)' icomp(1,1)=',icomp(1,1),' writepol_v2 ckpt 1'
!***
      NA=NAT0
      IXMIN=IXYZ0(1,1)
      IYMIN=IXYZ0(1,2)
      IZMIN=IXYZ0(1,3)
      IXMAX=IXMIN
      IYMAX=IYMIN
      IZMAX=IZMIN
      DO J=2,NAT0
         IF(IXYZ0(J,1)<IXMIN)IXMIN=IXYZ0(J,1)
         IF(IXYZ0(J,1)>IXMAX)IXMAX=IXYZ0(J,1)
         IF(IXYZ0(J,2)<IYMIN)IYMIN=IXYZ0(J,2)
         IF(IXYZ0(J,2)>IYMAX)IYMAX=IXYZ0(J,2)
         IF(IXYZ0(J,3)<IZMIN)IZMIN=IXYZ0(J,3)
         IF(IXYZ0(J,3)>IZMAX)IZMAX=IXYZ0(J,3)
      ENDDO
      IOPOL=32
!      OPEN(UNIT=IOPOL,FILE=CFLPOL,FORM='UNFORMATTED')
      OPEN(UNIT=IOPOL,FILE=CFLPOL,ACCESS='STREAM')

!*** diagnostic
!      write(0,*)'writepol_v2 ckpt 1'
!      write(0,*)' icomp(1,1)=',icomp(1,1),' writepol_v2 ckpt 1'
!      write(0,*)'   about to write file=',cflpol
!***
      WRITE(IOPOL)NRWORD,NAT0,IANISO,                        &
                  IXMIN,IXMAX,IYMIN,IYMAX,IZMIN,IZMAX,NX,NY,NZ
!*** diagnostic
!      write(0,*)'writepol ckpt 2'
!      write(0,*)' nx,ny,nz=',nx,ny,nz
!      write(0,*)' icomp(1,1)=',icomp(1,1),' writepol_v2 ckpt 2'
!***
      WRITE(IOPOL)(IXYZ0(J,1),J=1,NAT0),(IXYZ0(J,2),J=1,NAT0), &
                  (IXYZ0(J,3),J=1,NAT0),JPBC,PYD,PZD,GAMMA,    &
                  A1,A2,DX,X0

! reduce composition vector to only 3*NAT0 elements
      JA=0
      DO J=1,NX*NY*NZ
         IF(ICOMP(J,1)>0)THEN
            JA=JA+1
            ISCR1(JA)=ICOMP(J,1)
            ISCR1(JA+NAT0)=ICOMP(J,2)
            ISCR1(JA+2*NAT0)=ICOMP(J,3)
         ENDIF
      ENDDO

!*** diagnostic
!      write(0,*)'writepol_v2 ckpt 3'
!      write(0,*)' icomp(1,1)=',icomp(1,1),' writepol_v2 ckpt 3'
!      write(0,*)'           ja=',ja
!      write(0,*)'         nat0=',nat0
!      write(0,*)'   j iscr1(j) = reduced ICOMP'
!      do j=1,3*nat0
!         write(0,fmt='(i5,i3)')j,iscr1(j)
!      enddo
!***

! write reduced composition vector ICOMP
      WRITE(IOPOL)(ISCR1(J),J=1,NAT0),(ISCR1(NAT0+J),J=1,NAT0), &
                  (ISCR1(2*NAT0+J),J=1,NAT0)

! write refractive index of ambient medium
      WRITE(IOPOL)NAMBIENT
!***
! polarization vector CXPOL has already been reduced

!*** diagnostic
!      write(0,*)'writepol ckpt 4: cxpol being written to disk'
!      write(0,*)'   j  ----- cxpol(j) -----  [reduced]'
!      do ja=1,3
!         do j=1,nat0
!            write(0,fmt='(i5,1p2e11.3)')(nat0*(ja-1)+j),cxpol(j,ja)
!         enddo
!      enddo
!*** diagnostic

      WRITE(IOPOL)WAVE
      WRITE(IOPOL)AEFF
      WRITE(IOPOL)AKR
      WRITE(IOPOL)CXE0R
      WRITE(IOPOL)(CXPOL(J),J=1,NAT0),(CXPOL(NAT0+J),J=1,NAT0), &
                  (CXPOL(2*NAT0+J),J=1,NAT0)

!*** diagnostic
!            write(0,*)'writepol ckpt 5'
!            write(0,*)'      nx,ny,nz=',nx,ny,nz
!            write(0,*)'          nat0=',nat0
!            if(nat<301), write out A matrix
!            if(nx*ny*nz<301)then
!               write(0,*)'diagnostic for nat < 301'
!               write(0,*)'                    reduced cxadia'
!               write(0,*)'   j ix iy iz ',       &
!                         ' -------  A_x -------', &
!                         ' -------  A_y -------', &
!                         ' -------  A_z -------'
!               do j=1,nat0
!                  write(0,fmt='(i5,3i3,1p6e11.3)')j,ixyz0(j,1),ixyz0(j,2), &
!                        ixyz0(j,3),cxadia(j),cxadia(j+nat0),cxadia(j+2*nat0) 
!               enddo
!               write(0,*)'                    reduced P'
!               write(0,*)'   j ix iy iz ',       &
!                         ' -------  P_x -------', &
!                         ' -------  P_y -------', &
!                         ' -------  P_z -------'
!               do j=1,nat0
!                  write(0,fmt='(i5,3i3,1p6e11.3)')j,ixyz0(j,1),ixyz0(j,2), &
!                     ixyz0(j,3),cxpol(j),cxpol(nat0+j),cxpol(2*nat0+j) 
!               enddo
!            endif ! endif(nat<301)
!***

! store elements of A matrix for future use
      WRITE(IOPOL)(CXADIA(J),J=1,NAT0),(CXADIA(J+NAT0),J=1,NAT0), &
                  (CXADIA(J+2*NAT0),J=1,NAT0),(CXAOFF(J,1),J=1,NAT0), &
                  (CXAOFF(J,2),J=1,NAT0),(CXAOFF(J,3),J=1,NAT0)
      IF(IANISO==2)THEN
         WRITE(IOPOL)(BETADF(J),J=1,NAT0),(THETADF(J),J=1,NAT0), &
                     (PHIDF(J),J=1,NAT0)
      ENDIF
      CLOSE(IOPOL)
!*** diagnostic
!      write(0,*)'writepol_v2 ckpt 99 icomp(1,1)=',icomp(1,1)
!**
      RETURN
    END SUBROUTINE WRITEPOL
