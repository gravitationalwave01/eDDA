    SUBROUTINE READPOL(MODE,MXNAT,NX,NY,NZ,NAT0,IANISO,ICOMP,IXYZ0,PYD,PZD, &
                       AKD,DX,X0,WAVE,BETADF,THETADF,PHIDF,CXE0,CXADIA,     &
                       CXAOFF,CXPOL,CFLPOL)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

!------------------------- readpol_v2 ---------------------------------------
! arguments:

      INTEGER :: IANISO,MODE,MXNAT,NAT0,NX,NY,NZ

      INTEGER*2 :: &
         ICOMP(MXNAT,3)
      INTEGER ::      &
         IXYZ0(MXNAT,3)

      REAL(WP) :: PYD,PZD,WAVE
      REAL(WP) ::        &
         AKD(3),         &
         BETADF(MXNAT),  &
         DX(3),          &
         PHIDF(MXNAT),   &
         THETADF(MXNAT), &
         X0(3)
      COMPLEX(WP) ::      &
         CXADIA(MXNAT,3), &
         CXAOFF(MXNAT,3), &
         CXE0(3),         &
         CXPOL(MXNAT,3)
      CHARACTER :: CFLPOL*80

! local variables:

      INTEGER :: J,NA,NRWORD

!***********************************************************************
! subroutine READPOL
! given:
!    CFLPOL= file name
!    MXNAT = dimensioning parameter

! when called with MODE=0:
! returns:
!    NAT0   = number of dipoles in target
!    NRWORD = real word length (byte)
!           = 4 if data was stored in single precision
!           = 8 if data was stored in double precision
!
! when called with MODE=1
! returns:
!    NAT0   = number of occupied dipole sites
!    NX     = (1/d)*(maximum extent of target in x-direction) (TF)
!    NY     = (1/d)*(maximum extent of target in y-direction) (TF)
!    NZ     = (1/d)*(maximum extent of target in z-direction) (TF)
!    IXYZ0(IA,JD)=x/(d*DX(1)),y/(d*DX(2)),z/(d*DX(3))
!                 for dipole # IA, JD=1-3
!                 where (x,y,z) are in Target Frame
!    DX(1-3) =lattice spacing/d in x,y,z directions
!            =(1,1,1) for cubic lattice
!    X0(1-3) = vector specifying location (x,y,z) in TF of
!              lattice site (0,0,0):
!              X0(1)=x/(d*DX(1))
!              X0(2)=y/(d*DX(2))
!              X0(3)=z/(d*DX(3))
!    PYD     = periodicity/d of target replication in y_TF direction
!            = 0 if no replication in y direction
!    PZD     = periodicity/d of target replication in z_TF direction
!            = 0 if no replication in z direction
!    AKD(1-3)=k_x*d,k_y*d,k_z*d where (k_x,k_y,k_z)=k in ambient medium
!             in Target Frame
!    WAVE  = wavelength in ambient medium (physical units)
!    CXE0(1-3)=E_x,E_y,E_z for incident wave in ambient medium
!              in Target Frame
!    CXADIA(3i-2) =alphainv_xx at location i in TF
!          (3i-1) =alphainv_yy at location i in TF
!          (3i)   =alphainv_zz at location i in TF
!          where alphainv = inverse of (complex polarizability tensor)/d
!    CXAOFF(3i-2) =alphainv_zy at location i in TF
!          (3i-1) =alphainv_xz at location i in TF
!          (3i)   =alphainv_yx at location i in TF
!          where alphainv = inverse of (complex polarizability matrix/d^
!    CXPOL(i)        =P_x at location i in TF
!         (NAT0+i)   =P_y at location i
!         (2*NAT0+i) =P_z at location i

! Note: IXYZ0, CXPOL, CXADIA, CXAOFF are "reduced" arrays
!       IXYZ0 is produced by subroutine EXTEND
!       CXPOL is reduced by subroutine REDUCE called from GETFML

!       CXADIA is reduced by REDUCE called from GETFLM
!       CXAOFF is reduced by REDUCE called from GETFLM

!       where "reduced" arrays are limited to values for occupied sites

!       If one wishes to reconstitute "full" arrays

!       CXP(JX,JY,JZ,1-3)   polarization at (JX,JY,JZ)
!       CXAD(1-3,JX,JY,JZ)  diagonal elements of alphainv =
!                           inverse of polarizability
!                           tensor alpha/d^3 at (JX,JY,JZ)
!       CXAO(1-3,JX,JY,JZ)  off-diag elements yz,xz,xy of alphainv
!                           at (JX,JY,JZ)


!       one will need to add code like the following

!      COMPLEX
!     &   CXAD(3,NX,NY,NZ),
!     &   CXAO(3,NX,NY,NZ),
!     &   CXP(NX,NY,NZ,3)
!      DO JC=1,3
!         DO JZ=1,NZ
!            DO JY=1,NY
!               DO JX=1,NX
!                  CXP(JX,JY,JZ,JC)=0.
!                  CXAD(JC,JX,JY,JZ)=0.
!                  CXAO(JC,JX,JY,JZ)=0.
!               ENDDO
!            ENDDO
!         ENDDO
!      ENDDO
!      DO JA=1,NAT0
!         JX=IXYZ0(JA,1)
!         JY=IXYZ0(JA,2)
!         JZ=IXYZ0(JA,3)
!         DO JC=1,3
!            CXP(JX,JY,JZ,JC)=CXPOL(JA+(JC-1)*NAT0)
!            CXAD(JC,JX,JY,JZ)=CXADIA(3*(JA-1)+JC)
!            CXAO(JC,JX,JY,JZ)=CXAOFF(3*(JA-1)+JC)
!         ENDDO
!      ENDDO

! B.T. Draine, Princeton Univ. Observatory, 2006.04.08
! history
! 06.04.08 (BTD) first written
! 06.04.10 (BTD) revised
! 06.07.07 (BTD) revised comments
! 06.09.21 (BTD) add PYD,PZD,DX to argument list
!                add PYD,PZD,DX to list of variable read from input file
! 07.06.21 (BTD) modify to use X0(1-3) written out by DDSCAT v7.0.2
!                * add X0 to argument
!                * add X0 to READ statement
! 07.09.11 (BTD) changed IXYZ0 from INTEGER*2 to INTEGER
! 08.01.17 (BTD) * add IANISO to argument list
!                * read IANISO from file
!                * use IANISO to determine whether BETADF,THETADF,
!                  PHIDF need to be read from file; if IANISO<2 then
!                  set BETADF=0,THETADF=0,PHIDF=0
! end history

! Copyright (C) 2006,2007,2008
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
!*** diagnostic
!      write(0,*)'entered readpol...'
!      write(0,*)'cflpol=',cflpol
!***
      OPEN(UNIT=22,FILE=CFLPOL,FORM='UNFORMATTED')
!*** diagnostic
!      write(0,*)'readpol ckpt alpha'
!***
      READ(22)NRWORD,NAT0,IANISO
!*** diagnostic
!      write(0,*)'readpol ckpt beta'
!***
      IF(WP==KIND(0.E0).AND.NRWORD==8)THEN
         WRITE(0,*)'Fatal error in READPOL:'
         WRITE(0,*)'   DDfield was compiled in single precision,'
         WRITE(0,*)'   but ',CFLPOL,' contains double precision data'
         STOP
      ELSEIF(WP==KIND(0.D0).AND.NRWORD==4)THEN
         WRITE(0,*)'Fatal error in READPOL:'
         WRITE(0,*)'   DDfield was compiled in double precision,'
         WRITE(0,*)'   but ',CFLPOL,' contains single precision data'
         STOP
      ENDIF
      IF(MODE==0)THEN
         CLOSE(22)
         RETURN
      ENDIF

!*** sanity check

      IF(NAT0>MXNAT)THEN
         WRITE(0,*)'Fatal error in READPOL:'
         WRITE(0,*)'   attempting to read file, but NAT0 > MXNAT'
         CLOSE(22)
         STOP
      ENDIF

!*** diagnostic
!      write(0,*)'readpol ckpt gamma'
!***
      NA=NAT0
      READ(22)NX,NY,NZ,(IXYZ0(J,1),J=1,NA),(IXYZ0(J,2),J=1,NA),            &
              (IXYZ0(J,3),J=1,NA),PYD,PZD,AKD,DX,X0,WAVE,CXE0,             &
              (CXADIA(J,1),J=1,NA),(CXADIA(J,2),J=1,NA),                   &
              (CXADIA(J,3),J=1,NA),(CXAOFF(J,1),J=1,NA),                   &
              (CXAOFF(J,2),J=1,NA),(CXAOFF(J,3),J=1,NA),                   &
              (CXPOL(J,1),J=1,NA),(CXPOL(J,2),J=1,NA),(CXPOL(J,3),J=1,NA), &
              (ICOMP(J,1),J=1,NA),(ICOMP(J,2),J=1,NA),(ICOMP(J,3),J=1,NA)
      IF(IANISO==2)THEN
         READ(22)(BETADF(J),J=1,NA),(THETADF(J),J=1,NA),(PHIDF(J),J=1,NA)
      ELSE
         DO J=1,NA
            BETADF(J)=0._WP
            THETADF(J)=0._WP
            PHIDF(J)=0._WP
         ENDDO
      ENDIF
      CLOSE(22)

      RETURN
    END SUBROUTINE READPOL
