    PROGRAM DDPOL
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

! arguments of READPOL
      INTEGER :: IANISO,MODE,MXNAT,NAT0
      INTEGER*2,ALLOCATABLE :: ICOMP(:,:)
      INTEGER,ALLOCATABLE :: IXYZ0(:,:)
      INTEGER :: NX,NY,NZ
      REAL(WP) :: PYD,PZD,WAVE
      REAL(WP) :: &
         AKD(3),  &
         DX(3),   &
         X0(3)
      REAL(WP),ALLOCATABLE :: &
         BETADF(:),           &
         PHIDF(:),            &
         THETADF(:)
      COMPLEX(WP) :: &
         CXE0(3)
      COMPLEX(WP),ALLOCATABLE :: &
         CXADIA(:,:),            &
         CXAOFF(:,:),            &
         CXPOL(:,:)
      CHARACTER :: CFLPOL*80

! local variables

      INTEGER ::                       &
         IAT,JA,JD,JX,JXMAX,JXMIN,     &
         JY,JYMAX,JYMIN,JZ,JZMAX,JZMIN
      REAL(WP) ::                          &
         CWORD,DSTORAGE,EX0,EY0,EZ0,KD,MB, &
         PI,PSIX,PSIY,PSIZ,                &
         RWORD,STORAGE,STORAGE0,X,Y,Z

      COMPLEX(WP) :: &
         CXB0(3)

!=======================================================================
!                          DDpol
!
! Program DDpol takes the target polarization array output by DDSCAT,
! and reports polarization at selected points in the target.
!
! Input: file DDpol.in with list of coordinates, given as
!        x/dx(1), y/dx(2), z/dx(3)
!
!        where dx(j)=lattice spacing in direction j
!        with dx(1)*dx(2)*dx(3)=d**3
!
! history:
! 06.09.22 (BTD) created from DDfield as template
!                write PYD,PZD to output files
! 08.01.17 (BTD) modified to use new version of READPOL, with additional
!                arguments MODE,IANISO,ICOMP,BETADF,THETADF,PHIDF
!                modified to write out information concerning target
!                geometry
! 08.10.30 (BTD) corrected 
!                * typos in memory allocation
!                  (errors pointed out by Michel Devel)
!                * must allocate BETADF,THETADF,PHIDF even when
!                  IANISO=1, since readpol writes zeros.
! 09.08.12 (BTD) bug reported by Rodrigo Alcazar de la Osa (Univ. of
!                Cantabria, Spain)
!                corrected declaration of CFLPOL*16 -> CFLPOL*80
!                as required by READPOL
! 09.08.13 (BTD) several declaration statements 
!                arrays were inconsistent with READPOL:
!                changes:
!                INTEGER -> INTEGER*2 for ICOMP
!                REAL -> REAL(WP) for real variables
!                COMPLEX -> COMPLEX(WP) for complex variables
! end history

! Copyright (C) 2006,2007,2008,2009
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!=======================================================================

      PI=4.*ATAN(1.)

! for storage computations:

      MB=REAL(1024**2)
      IF(WP==KIND(0.E0))THEN
         RWORD=4._WP
         CWORD=2._WP*RWORD
         STORAGE0=6.79
      ELSEIF(WP==KIND(0.D0))THEN
         RWORD=8._WP
         CWORD=2._WP*RWORD
         STORAGE0=6.794
      ELSE
         WRITE(0,*)'Fatal error in DDfield: unable to determine word length'
         STOP
      ENDIF
      STORAGE=STORAGE0

! Input control file:

      OPEN(UNIT=3,FILE='DDpol.in')

! Output files:

      OPEN(UNIT=7,FILE='DDpol.out')

! Read name of file containing stored polarization information:

! temporary allocation

      MXNAT=1
      ALLOCATE(ICOMP(MXNAT,3))
      ALLOCATE(IXYZ0(MXNAT,3))
      ALLOCATE(BETADF(MXNAT))
      ALLOCATE(THETADF(MXNAT))
      ALLOCATE(PHIDF(MXNAT))
      ALLOCATE(CXADIA(MXNAT,3))
      ALLOCATE(CXAOFF(MXNAT,3))
      ALLOCATE(CXPOL(MXNAT,3))

      READ(3,*)CFLPOL

! skip one line so that one can use the same file structure as DDfield.in
      READ(3,*)

! Use routine READPOL to read in the polarization array from file CFLPOL.
! call first with MODE=0 to obtain size information

      MODE=0
      CALL READPOL(MODE,MXNAT,NX,NY,NZ,NAT0,IANISO,ICOMP,IXYZ0,PYD,PZD,AKD,  &
                   DX,X0,WAVE,BETADF,THETADF,PHIDF,CXE0,CXADIA,CXAOFF,CXPOL, &
                   CFLPOL)

! Now allocate necessary storage:

      DEALLOCATE(ICOMP)
      DEALLOCATE(IXYZ0)
      DEALLOCATE(CXADIA)
      DEALLOCATE(CXAOFF)
      DEALLOCATE(CXPOL)

      MXNAT=NAT0
      STORAGE=STORAGE0

      DSTORAGE=REAL(2*3*MXNAT)/MB
      STORAGE=STORAGE+DSTORAGE
      WRITE(0,6605)DSTORAGE,STORAGE
      ALLOCATE(ICOMP(MXNAT,3))

      DSTORAGE=REAL(4*3*MXNAT)/MB
      STORAGE=STORAGE+DSTORAGE
      WRITE(0,6610)DSTORAGE,STORAGE
      ALLOCATE(IXYZ0(MXNAT,3))

      DSTORAGE=CWORD*REAL(3*MXNAT)/MB
      STORAGE=STORAGE+DSTORAGE
      WRITE(0,6620)DSTORAGE,STORAGE
      ALLOCATE(CXADIA(MXNAT,3))

      DSTORAGE=CWORD*REAL(3*MXNAT)/MB
      STORAGE=STORAGE+DSTORAGE
      WRITE(0,6630)DSTORAGE,STORAGE
      ALLOCATE(CXAOFF(MXNAT,3))

      DSTORAGE=CWORD*REAL(3*MXNAT)/MB
      STORAGE=STORAGE+DSTORAGE
      WRITE(0,6640)DSTORAGE,STORAGE
      ALLOCATE(CXPOL(MXNAT,3))

      DEALLOCATE(BETADF)
      DEALLOCATE(THETADF)
      DEALLOCATE(PHIDF)

      DSTORAGE=RWORD*REAL(MXNAT)/MB
      STORAGE=STORAGE+DSTORAGE
      WRITE(0,6660)DSTORAGE,STORAGE
      ALLOCATE(BETADF(MXNAT))

      DSTORAGE=RWORD*REAL(MXNAT)/MB
      STORAGE=STORAGE+DSTORAGE
      WRITE(0,6670)DSTORAGE,STORAGE
      ALLOCATE(THETADF(MXNAT))

      DSTORAGE=RWORD*REAL(MXNAT)/MB
      STORAGE=STORAGE+DSTORAGE
      WRITE(0,6680)DSTORAGE,STORAGE
      ALLOCATE(PHIDF(MXNAT))

      MODE=1
      CALL READPOL(MODE,MXNAT,NX,NY,NZ,NAT0,IANISO,ICOMP,IXYZ0,PYD,PZD,AKD,  &
                   DX,X0,WAVE,BETADF,THETADF,PHIDF,CXE0,CXADIA,CXAOFF,CXPOL, &
                   CFLPOL)

! Check target dimensions

      JXMIN=IXYZ0(1,1)
      JYMIN=IXYZ0(1,2)
      JZMIN=IXYZ0(1,3)
      JXMAX=JXMIN
      JYMAX=JYMIN
      JZMAX=JZMIN
      DO JA=2,NAT0
         JX=IXYZ0(JA,1)
         JY=IXYZ0(JA,2)
         JZ=IXYZ0(JA,3)
         IF(JX.LT.JXMIN)JXMIN=JX
         IF(JX.GT.JXMAX)JXMAX=JX
         IF(JY.LT.JYMIN)JYMIN=JY
         IF(JY.GT.JYMAX)JYMAX=JY
         IF(JZ.LT.JZMIN)JZMIN=JZ
         IF(JZ.GT.JZMAX)JZMAX=JZ
      ENDDO

      KD=0.
      DO JD=1,3
         KD=KD+AKD(JD)**2
      ENDDO
      KD=SQRT(KD)

! KD = k*d

! Calculate incident B field:

      CXB0(1)=(AKD(2)*CXE0(3)-AKD(3)*CXE0(2))/KD
      CXB0(2)=(AKD(3)*CXE0(1)-AKD(1)*CXE0(3))/KD
      CXB0(3)=(AKD(1)*CXE0(2)-AKD(2)*CXE0(1))/KD

      WRITE(7,7001)NAT0,JXMIN,JXMAX,JYMIN,JYMAX,JZMIN,JZMAX,PYD,PZD,AKD,CXE0

 1000 READ(3,*,END=9000,ERR=9000)X,Y,Z
      WRITE(0,6100)X,Y,Z

      DO JA=1,NAT0
         IF(IXYZ0(JA,1).EQ.INT(X))THEN
            IF(IXYZ0(JA,2).EQ.NINT(Y))THEN
               IF(IXYZ0(JA,3).EQ.NINT(Z))THEN
                  IAT=JA
               ENDIF
            ENDIF
         ENDIF
      ENDDO
            
      EX0=ABS(CXPOL(IAT,1))
      PSIX=0.
      IF(EX0.GT.0.)THEN
         PSIX=ACOS(REAL(CXPOL(IAT,1))/EX0)
         IF(AIMAG(CXPOL(IAT,1)).LT.0.)PSIX=2.*PI-PSIX
      ENDIF
      EX0=4.*PI*EX0

      EY0=ABS(CXPOL(IAT,2))
      PSIY=0.
      IF(EY0.GT.0.)THEN
         PSIY=ACOS(REAL(CXPOL(IAT,2))/EY0)
         IF(AIMAG(CXPOL(IAT,2)).LT.0.)PSIY=2*PI-PSIY
      ENDIF
      EY0=4.*PI*EY0

      EZ0=ABS(CXPOL(IAT,3))
      PSIZ=0.
      IF(EZ0.GT.0.)THEN
         PSIZ=ACOS(REAL(CXPOL(IAT,3))/EZ0)
         IF(AIMAG(CXPOL(IAT,3)).LT.0.)PSIZ=2*PI-PSIZ
      ENDIF
      EZ0=4.*PI*EZ0

      PSIX=180.*PSIX/PI
      PSIY=180.*PSIY/PI
      PSIZ=180.*PSIZ/PI

      WRITE(7,7100)X,Y,Z,EX0,PSIX,EY0,PSIY,EZ0,PSIZ

      GOTO 1000
!-----------------------------------------------------------------------
 9000 CLOSE(7)
      STOP
 6100 FORMAT('calc E and B at (x/dx,y/dy,z/dz)=',3F10.3)
 6605 FORMAT('allocating',F8.3,' MB for ICOMP ; total=',F10.3,' MB')
 6610 FORMAT('allocating',F8.3,' MB for IXYZ0 ; total=',F10.3,' MB')
 6620 FORMAT('allocating',F8.3,' MB for CXADIA; total=',F10.3,' MB')
 6630 FORMAT('allocating',F8.3,' MB for CXAOFF; total=',F10.3,' MB')
 6640 FORMAT('allocating',F8.3,' MB for CXPOL ; total=',F10.3,' MB')
 6650 FORMAT('allocating',F8.3,' MB for DCXE  ; total=',F10.3,' MB')
 6660 FORMAT('allocating',F8.3,' MB for BETADF; total=',F10.3,' MB')
 6670 FORMAT('allocating',F8.3,' MB for THETADF;total=',F10.3,' MB')
 6680 FORMAT('allocating',F8.3,' MB for PHIDF;  total=',F10.3,' MB')
 7001 FORMAT(I10,' = number of dipoles in Target',/,              &
             'Target dipole positions in Target Frame:',/,        &
             2I6,' = min,max values of x/d',/,                    &
             2I6,' = min,max values of y/d',/,                    &
             2I6,' = min,max values of z/d',/,/,                  &
             2F10.4,' = PYD, PZD = period_y/dy, period_z/dz',/,/, &
             1PE12.5,' = k_x*d for incident wave',/,              &
             1PE12.5,' = k_y*d for incident wave',/,              &
             1PE12.5,' = k_z*d for incident wave',/,              &
             '(',0PF10.6,' , ',0PF10.6,' ) = E_inc,x(0,0,0)',/,   &
             '(',0PF10.6,' , ',0PF10.6,' ) = E_inc,y(0,0,0)',/,   &
             '(',0PF10.6,' , ',0PF10.6,' ) = E_inc,z(0,0,0)',/,   &
             4X,'x/d',6X,'y/d',6X,'z/d',2X,'4pi|P_x| psi_x',1X,   &
             '4pi|P_y| psi_y',1X,'4pi|P_z| psi_z')
 7100 FORMAT(3F9.4,F9.5,F6.1,F8.5,F6.1,F8.5,F6.1)
      END
