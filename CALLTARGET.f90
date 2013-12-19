    PROGRAM CALLTARGET

      USE DDPRECISION,ONLY : WP
      USE DDCOMMON_1,ONLY : DX
      USE DDCOMMON_6,ONLY : GAMMA,PYD,PZD,MXNATF,MXNXF,MXNYF,MXNZF,NAT,NAT3, &
                            NAT0,NX,NY,NZ,MXN3F,IDVOUT,IPBC
      IMPLICIT NONE

! Program CALLTARGET is used to generate target arrays using target
! generation routines employed by DDSCAT.

! B.T.Draine, Princeton Univ. Obs., 91/1/1
! History:
! 91.05.24 (BTD): Modified to allow entry target names in lower case
! 91.07.16 (BTD): Replaced DO...ENDDO with DO #...# CONTINUE
! 92.12.16 (BTD): Added IDVSHP to argument list of TARGET
! 93.01.07 (BTD): Added options TWOELL and TWOAEL
! 94.05.15 (BTD): Added option CONELL
! 95.12.11 (BTD): Added option BLOCKS
! 96.01.04 (BTD): Added option DW1996
! 96.01.25 (BTD): Added option TWOSPH
! 96.01.26 (BTD): Replaced IX,IY,IZ by IXYZ here and in TARGET
!                 Added calls to SIZER
! 96.01.29 (BTD): Support selection of principal axis calculation
!                 in TAR2SP via SHPAR(6)
! 97.04.25 (BTD): Added option BLK_AN
! 97.04.30 (BTD): Removed code offering choice on target axis for
!                 HEXGON option (target axis now hardwired within
!                 TARHEX)
! 00.06.12 (BTD): Added option NSPHER
! 02.02.13 (BTD): Added option PRISM3
! 02.11.12 (BTD): Added option SPHARM
! 03.11.06 (BTD): Modified to input parameters to generate random
!                 SPHARM target
!                 Added DX(3) to argument list of TARGET
! 03.11.08 (BTD): Added option GSPHER
! 04.03.31 (BTD): Replaced all WRITE(0, with WRITE(IDVOUT,
! 04.04.01 (BTD): Added NPY,NPZ to argument list of TARGET,
!                 so we can use new version of TARGET that supports
!                 periodic boundary conditions
! 04.10.02 (BTD): Declared BETADF,PHIDF,THETADF
!                 Added BETADF,PHIDF,THETADF to argument list of
!                 subroutine TARGET
! 04.10.19 (BTD): Modified to allow repeated trials of different
!                 SHPAR(1) for NSPHER and NANSPH options
! 05.06.16 (BTD): Replaced integer NPY,NPZ by real PYD,PZD in
!                 argument list of TARGET
! 06.09.13 (BTD): Added option CYLCAP
! 06.09.15 (BTD): Modified option HEXGON to allow orientation to
!                 be changed.
! 06.12.08 (BTD): Modified to support PBXPBC option
! 08.08.30 (BTD): Added IANISO to argument list of TARGET
! 08.08.31 (BTD): Added options including BISLINPBC
! 08.09.12 (BTD): Correcte typos DSKRCTGNL -> DSKRCTNGL
! 09.09.22 (BTD): Added NCOMP_NEED to argument list in each call to
!                 TARGET
! 10.02.06 (BTD): Added support for options 
!                 * CYLNDRPBC
!                 * HEXGONPBC
!                 * LYRSLBPBC
!                 * RCTGL_PBC
!                 * SLAB_HOLE
!                 * SLBHOLPBC
!                 * SPHRN_PBC
!                 * TRILYRPBC
! end history
!
! Copyright (C) 1993,1994,1995,1996,1997,2000,2002,2003,2004,2005
!               2006,2008,2009,2010 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
!
! Adjustable parameters:
!    MXNAT=max.number of dipoles in target
!
      INTEGER :: MXN3,MXNAT

! Arguments of TARGET:

      CHARACTER :: CDESCR*67,CFLSHP*80,CSHAPE*9

      INTEGER*2,ALLOCATABLE :: &
     &   ICOMP(:,:)

      INTEGER,ALLOCATABLE :: &
     &   IXYZ(:,:)

      INTEGER :: IANISO,IDVSHP,IOSHP
      REAL(WP) ::   &
         A1(3),     &
         A2(3),     &
         SHPAR(12), &
         X0(3)
      REAL(WP),ALLOCATABLE :: &
         BETADF(:),           &
         PHIDF(:),            &
         THETADF(:)

! Local variables:

      INTEGER :: J,JJ,NCOMP_NEED
      INTEGER :: LXYZ(3)
      REAL(WP) :: ERR,F,PI,RGYR2,RXY,RXZ,RYX,RYZ,RZX,RZY, &
              X,XCM,Y,YCM,Z,ZCM
      REAL(WP) :: DDX,DDY,DDZ,X2,Y2,Z2
!***********************************************************************
      DATA IOSHP/15/

      MXNAT=2000000
      MXN3=3*MXNAT

      ALLOCATE(IXYZ(MXNAT,3))
      ALLOCATE(ICOMP(MXNAT,3))
      ALLOCATE(BETADF(MXNAT))
      ALLOCATE(PHIDF(MXNAT))
      ALLOCATE(THETADF(MXNAT))

      PI=4.*ATAN(1.)

! Define variables CFLSHP and IDVSHP:

      CFLSHP='none'
      IDVSHP=10
 1000 WRITE(IDVOUT,*)'What shape? (Enter choice in quotes) Choices:'
      WRITE(IDVOUT,*)'  ANIELLIPS = anisotropic ellipsoid'
      WRITE(IDVOUT,*)'  ANI_ELL_2 = 2 touching identical anisotropic ', &
                     'ellipsoids'
      WRITE(IDVOUT,*)'  ANI_ELL_3 = 3 touching identical anisotropic ', &
                     'ellipsoids'
      WRITE(IDVOUT,*)'  ANIRCTNGL = anisotropic rectangular solid'
      WRITE(IDVOUT,*)'  CONELLIPS = 2 concentric ellipsoids'
      WRITE(IDVOUT,*)'  CYLINDER1 = homogenous, isotropic cylinder'
      WRITE(IDVOUT,*)'  CYLNDRCAP = homogenous, isotropic cylinder with ', &
                     'hemispherical endcaps'
      WRITE(IDVOUT,*)'  DSKRCTNGL = disk on homogeneous rectangular block'
      WRITE(IDVOUT,*)'  DW1996TAR = 13 cube target used by Draine & ', &
                     'Weingartner (1996)'
      WRITE(IDVOUT,*)'  ELLIPSOID = homogenous, isotropic ellipsoid'
      WRITE(IDVOUT,FMT='(A,A)')'  ELLIPSO_2 = 2 touching ellipsoids; 2nd ', &
                     'ellipsoid displaced from 1st in x-direction'
      WRITE(IDVOUT,FMT='(A,A)')'  ELLIPSO_3 = 3 touching ellipsoids; 2nd ', &
                     'ellipsoid displaced from 1st in +x-direction,'
      WRITE(IDVOUT,*)'              3rd displaced from 2nd in +x-direction'
!      WRITE(IDVOUT,*)'  GAUSS_SPH = gaussian sphere target'
      WRITE(IDVOUT,*)'  HEX_PRISM = hexagonal prism'
      WRITE(IDVOUT,FMT='(A,A)')'  MLTBLOCKS = construct from cubic blocks', &
                     ' input from file "blocks.par"'
      WRITE(IDVOUT,*)'  RCTGLPRSM = rectangular prism'
      WRITE(IDVOUT,*)'  RCTGLBLK3 = stack of 3 rectangular blocks'
      WRITE(IDVOUT,*)'  SLAB_HOLE = rect. block with cylindrical hole'
      WRITE(IDVOUT,*)'  SPHERES_N = union of N spheres'
      WRITE(IDVOUT,*)'  SPHROID_2 = 2 touching spheroids (with angle ', &
                     'phi between symm. axes)'
      WRITE(IDVOUT,*)'  SPH_ANI_N = N-spheres of anisotropic materials'
      WRITE(IDVOUT,*)'  TETRAHDRN = regular tetrahedron'
      WRITE(IDVOUT,*)'  TRNGLPRSM = triangular prism'
!---------------------- PBC target options start here ---------------------
      WRITE(IDVOUT,*)'  BISLINPBC = bilayer slab with parallel lines (one TUC)'
      WRITE(IDVOUT,*)'  CYLNDRPBC = circular slice (TUC for infinite cylinder)'
      WRITE(IDVOUT,*)'  DSKBLYPBC = disk on bilayer slab (one TUC)'
      WRITE(IDVOUT,*)'  DSKRCTPBC = disk on slab (one TUC)'
      WRITE(IDVOUT,*)'  HEXGONPBC = hexagonal slab (one TUC)'
      WRITE(IDVOUT,*)'  LYRSLBPBC = rect. block with 2,3, or 4 layers (one TUC)'
      WRITE(IDVOUT,*)'  RECRECPBC = 2 stacked blocks (one TUC)'
      WRITE(IDVOUT,*)'  RCTGL_PBC = 1 block (one TUC)'
      WRITE(IDVOUT,*)'  SLBHOLPBC = block with cylindrical hole (one TUC)'
      WRITE(IDVOUT,*)'  SPHRN_PBC = N possibly anisotropic spheres (one TUC)'
      WRITE(IDVOUT,*)'  TRILYRPBC = 3 stacked rect. blocks (one TUC)'
! at this time, restrict to cubic lattice:

      DX(1)=1.
      DX(2)=1.
      DX(3)=1.

      READ(*,*)CSHAPE
!-----------------------------------------------------------------------
      IF(CSHAPE=='ANIELLIPS'.OR.CSHAPE=='aniellips')THEN
         CSHAPE='ANIELLIPS'
         DO J=1,100
            WRITE(IDVOUT,*)'Enter x-length/d, y-length/d, z-length/d ', &
                           '(0 0 0 to stop)'
            READ(*,*)SHPAR(1),SHPAR(2),SHPAR(3)
            IF(SHPAR(1).LE.0.)STOP
            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT,  &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                 &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,      &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
            WRITE(IDVOUT,7003)SHPAR(1),SHPAR(2),SHPAR(3),NAT
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='ANI_ELL_2'.OR.CSHAPE=='ani_ell_2')THEN
         CSHAPE='ANI_ELL_2'
         DO J=1,100
            WRITE(IDVOUT,*)'Enter length of one ellipsoid in x, y, z', &
                           ' directions (0 0 0 to stop)'
            READ(*,*)SHPAR(1),SHPAR(2),SHPAR(3)
            IF(SHPAR(1).LE.0)STOP
            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
            WRITE(IDVOUT,7003)(SHPAR(JJ),JJ=1,3),NAT
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='ANI_ELL_3'.OR.CSHAPE=='ani_ell_3')THEN
         CSHAPE='ANI_ELL_3'
         DO J=1,100
            WRITE(IDVOUT,*)'Enter length of one ellipsoid in x, y, z', &
                           ' directions (0 0 0 to stop)'
            READ(*,*)SHPAR(1),SHPAR(2),SHPAR(3)
            IF(SHPAR(1).LE.0)STOP
            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
            WRITE(IDVOUT,7003)(SHPAR(JJ),JJ=1,3),NAT
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='ANIRCTNGL'.OR.CSHAPE=='anirctngl')THEN
         CSHAPE='ANIRCTNGL'
         DO J=1,100
            WRITE(IDVOUT,*)'Enter length/d in x, y, z', &
                           ' directions (0 0 0 to stop)'
            READ(*,*)SHPAR(1),SHPAR(2),SHPAR(3)
            IF(SHPAR(1).LE.0)STOP
            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
            WRITE(IDVOUT,7003)(SHPAR(JJ),JJ=1,3),NAT
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='BISLINPBC'.OR.CSHAPE=='bislinpbc')THEN
         CSHAPE='BISLINPBC'
         WRITE(IDVOUT,*)'TUC for bilayer slab with parallel line on top'
         DO J=1,100
            WRITE(IDVOUT,*)'Enter x-thickness/d of line (0 to stop)'
            READ(*,*)SHPAR(1)
            IF(SHPAR(1).LE.0)STOP
            WRITE(IDVOUT,*)'Enter y-width/d of line'
            READ(*,*)SHPAR(2)
            WRITE(IDVOUT,*)'Enter x-thickness/d of upper layer'
            READ(*,*)SHPAR(3)
            WRITE(IDVOUT,*)'Enter x-thickness/d of lower layer'
            READ(*,*)SHPAR(4)
            WRITE(IDVOUT,*)'Enter L_y/d dimension of slab'
            READ(*,*)SHPAR(5)
            SHPAR(6)=0   ! P_y/d
            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
            WRITE(IDVOUT,7003)(SHPAR(JJ),JJ=1,3),NAT
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='CONELLIPS'.OR.CSHAPE=='conellips')THEN
         CSHAPE='CONELLIPS'
         DO J=1,100
            WRITE(IDVOUT,*)'Enter length of outer ellipsoid (icomp=2) ', &
                           'in x,y,z directions (lattice units) (0 0 0', &
                           ' to stop)'
            READ(*,*)SHPAR(1),SHPAR(2),SHPAR(3)
            IF(SHPAR(1).LE.0)STOP
            WRITE(IDVOUT,*)'Enter length of inner ellipsoid (icomp=1) ', &
                           'in x,y,z directions (lattice units)'
            READ(*,*)SHPAR(4),SHPAR(5),SHPAR(6)
            IF(SHPAR(4).LE.SHPAR(1).AND. &
               SHPAR(5).LE.SHPAR(2).AND. &
               SHPAR(6).LE.SHPAR(3))THEN
               CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,   &
                           MXNAT,SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,      &
                           MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0, &
                           IANISO,NCOMP_NEED)
               CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
               WRITE(IDVOUT,7006)(SHPAR(JJ),JJ=1,6),NAT
            ELSE
               WRITE(IDVOUT,*)'Input error: inner ellipsoid is not ', &
                              'contained within outer ellipsoid'
            ENDIF
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='CYLINDER1'.OR.CSHAPE=='cylinder1')THEN
         CSHAPE='CYLINDER1'
         DO J=1,100
            WRITE(IDVOUT,*)'Enter cylinder length (0 to stop)'
            READ(*,*)SHPAR(1)
            IF(SHPAR(1).LE.0.)STOP
            WRITE(IDVOUT,*)'Enter cylinder diameter'
            READ(*,*)SHPAR(2)
            SHPAR(3)=1     ! 1 for cylinder axis in x direction
            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
            WRITE(IDVOUT,7002)SHPAR(1),SHPAR(2),NAT
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='CYLNDRCAP'.OR.CSHAPE=='cylndrcap')THEN
         CSHAPE='CYLNDRCAP'
         WRITE(IDVOUT,*)' CYLNDRCAP: cylinder with hemispherical endcaps'
         DO J=1,100
            WRITE(IDVOUT,*)'Enter cylinder length (not including ', &
                           'end-caps) (0 to stop)'
            READ(*,*)SHPAR(1)
            IF(SHPAR(1).LE.0.)STOP
            WRITE(IDVOUT,*)'Enter cylinder diameter'
            READ(*,*)SHPAR(2)
            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
            WRITE(IDVOUT,7002)SHPAR(1),SHPAR(2),NAT
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='CYLNDRPBC'.OR.CSHAPE=='cylndrpbc')THEN
         CSHAPE='CYLNDRPBC'
         WRITE(IDVOUT,*)' CYLNDRPBC: TUC for infinite cylinder' 
         DO J=1,100
            WRITE(IDVOUT,*)'Enter cylinder length/d (normally 1)', &
                           '(0 to stop)'
            READ(*,*)SHPAR(1)
            IF(SHPAR(1).LE.0.)STOP
            WRITE(IDVOUT,*)'Enter cylinder diameter/d'
            READ(*,*)SHPAR(2)
            WRITE(0,*)'Enter 1 for cylinder axis || x_TF'
            WRITE(0,*)'      2 for cylinder axis || y_TF'
            WRITE(0,*)'      3 for cylinder axis || z_TF'
            READ(*,*)SHPAR(3)
            SHPAR(4)=0   ! P_y/d
            SHPAR(5)=0   ! P_z/d
            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
            WRITE(IDVOUT,7005)SHPAR(1),SHPAR(2),SHPAR(3), &
                              SHPAR(4),SHPAR(5),NAT
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='DSKBLYPBC'.OR.CSHAPE=='dskblypbc')THEN
         CSHAPE='DSKBLYPBC'
         WRITE(IDVOUT,*)'TUC for array of (disk on top of bilayer slab)'
         WRITE(IDVOUT,*)'With disk (comp. 1) located on x=0 surface'
         WRITE(IDVOUT,*)'Enter disk thickness/d and diameter/d (0 0 to stop)'
         READ(*,*)SHPAR(1),SHPAR(2)
         WRITE(IDVOUT,*)'Enter thickness/d of upper layer of slab (comp. 2)'
         READ(*,*)SHPAR(3)
         WRITE(IDVOUT,*)'Enter thickness/d of lower layer of slab (comp. 3)'
         READ(*,*)SHPAR(4)
         WRITE(IDVOUT,*)'Enter slab extent Ly/d  Lz/d'
         READ(*,*)SHPAR(5),SHPAR(6)
         SHPAR(7)=0._WP   ! P_y/d
         SHPAR(8)=0._WP   ! P_z/d
         CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                     SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                     MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                     IANISO,NCOMP_NEED)
         CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='DSKRCTNGL'.OR.CSHAPE=='dskrctngl')THEN
         CSHAPE='DSKRCTNGL'
         WRITE(IDVOUT,*)'Rectangular slab with dimension Lx x Ly x Lz'
         WRITE(IDVOUT,*)'With disk located on x=0 surface'
         WRITE(IDVOUT,*)'Enter disk thickness/d and diameter/d'
         READ(*,*)SHPAR(1),SHPAR(2)
         WRITE(IDVOUT,*)'Enter Lx/d  Ly/d  Lz/d'
         READ(*,*)SHPAR(3),SHPAR(4),SHPAR(5)
         CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                     SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                     MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                     IANISO,NCOMP_NEED)
         CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='DSKRCTPBC'.OR.CSHAPE=='dskrctpbc')THEN
         CSHAPE='DSKRCTPBC'
         WRITE(IDVOUT,*)'TUC for array of (disk on slab)'
         WRITE(IDVOUT,*)'With disk located on x=0 surface'
         WRITE(IDVOUT,*)'Enter disk thickness/d and diameter/d (0 0 to stop)'
         READ(*,*)SHPAR(1),SHPAR(2)
         WRITE(IDVOUT,*)'Enter slab dimensions Lx/d  Ly/d  Lz/d'
         READ(*,*)SHPAR(3),SHPAR(4),SHPAR(5)
         SHPAR(6)=0._WP   ! P_y/d
         SHPAR(7)=0._WP   ! P_z/d
         CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                     SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                     MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                     IANISO,NCOMP_NEED)
         CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='DW1996TAR'.OR.CSHAPE=='dw1996tar')THEN
         CSHAPE='DW1996TAR'
         WRITE(IDVOUT,*)'Enter block size (lattice units)'
         READ(*,*)SHPAR(1)
         CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                     SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                     MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                     IANISO,NCOMP_NEED)
         CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
         WRITE(IDVOUT,7001)SHPAR(1),NAT
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='ELLIPSOID'.OR.CSHAPE=='ellipsoid')THEN
         CSHAPE='ELLIPSOID'
         DO J=1,100
            WRITE(IDVOUT,*)'Enter ellipsoid diameters xv,yv,zv', &
                           ' (0,0,0 to stop)'
            READ(*,*)SHPAR(1),SHPAR(2),SHPAR(3)
            IF(SHPAR(1).LE.0.)stop
            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
            XCM=0.
            YCM=0.
            ZCM=0.
            DO JJ=1,NAT
               X=IXYZ(JJ,1)
               Y=IXYZ(JJ,2)
               Z=IXYZ(JJ,3)
               XCM=XCM+X
               YCM=YCM+Y
               ZCM=ZCM+Z
            ENDDO
            XCM=XCM/REAL(NAT)
            YCM=YCM/REAL(NAT)
            ZCM=ZCM/REAL(NAT)
            X2=0.
            Y2=0.
            Z2=0.
            DO JJ=1,NAT
               DDX=IXYZ(JJ,1)-XCM
               DDY=IXYZ(JJ,2)-YCM
               DDZ=IXYZ(JJ,3)-ZCM
               X2=X2+DDX*DDX
               Y2=Y2+DDY*DDY
               Z2=Z2+DDZ*DDZ
            ENDDO
            X2=X2/NAT+1./12.
            Y2=Y2/NAT+1./12.
            Z2=Z2/NAT+1./12.
            ERR=20.*PI*SQRT(5.*REAL(X2*Y2*Z2))/(3.*REAL(NAT))-1.
            RXY=REAL(SQRT(X2/Y2))
            RYZ=REAL(SQRT(Y2/Z2))
            RZX=REAL(SQRT(Z2/X2))
            RYX=1./RXY
            RZY=1./RYZ
            RXZ=1./RZX
            X2=2.*SQRT(5.*X2)
            Y2=2.*SQRT(5.*Y2)
            Z2=2.*SQRT(5.*Z2)
            WRITE(*,6200)SHPAR(1),SHPAR(2),SHPAR(3),NAT,X2,Y2,Z2,ERR, &
                         RYX,RZX,RYZ,RXY,RXZ,RZY
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='ELLIPSO_2'.OR.CSHAPE=='ellipso_2')THEN
         CSHAPE='ELLIPSO_2'
         DO J=1,100
            WRITE(IDVOUT,*)'Enter x,y,z length/d for one ', &
                           'ellipsoid (0 0 0 to stop)'
            READ(*,*)SHPAR(1),SHPAR(2),SHPAR(3)
            IF(SHPAR(1).LE.0.)STOP
            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
            WRITE(IDVOUT,7003)SHPAR(1),SHPAR(2),SHPAR(3),NAT
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='ELLIPSO_3'.OR.CSHAPE=='ellipso_3')THEN
         CSHAPE='ELLIPSO_3'
         DO J=1,100
            WRITE(IDVOUT,*)'Enter x,y,z length/d for one ', &
                           'ellipsoid (0 0 0 to stop)'
            READ(*,*)SHPAR(1),SHPAR(2),SHPAR(3)
            IF(SHPAR(1).LE.0.)STOP
            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
            WRITE(IDVOUT,7003)SHPAR(1),SHPAR(2),SHPAR(3),NAT
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='HEX_PRISM'.OR.CSHAPE=='hex_prism')THEN
         CSHAPE='HEX_PRISM'
         DO J=1,100
            WRITE(IDVOUT,*)'Enter length (along symmetry axis)', &
                           ' (0 to stop)'
            READ(*,*)SHPAR(1)
            IF(SHPAR(1).LE.0.)STOP
            WRITE(IDVOUT,*)'Enter length of 1 hexagon side'
            READ(*,*)SHPAR(2)
            WRITE(IDVOUT,*)'Specify orientation of axes a1 (hex axis) ', &
                           'and a2'
            WRITE(IDVOUT,*)' 1 for a1=x,a2=y'
            WRITE(IDVOUT,*)' 2 for a1=x,a2=z'
            WRITE(IDVOUT,*)' 3 for a1=y,a2=x'
            WRITE(IDVOUT,*)' 4 for a1=y,a2=z'
            WRITE(IDVOUT,*)' 5 for a1=z,a2=x'
            WRITE(IDVOUT,*)' 6 for a1=z,a2=y'
            READ(*,*)SHPAR(3)
            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
            WRITE(IDVOUT,7003)(SHPAR(JJ),JJ=1,3),NAT
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='HEXGONPBC'.OR.CSHAPE=='hexgonpbc')THEN
         CSHAPE='HEXGONPBC'
         DO J=1,100
            WRITE(IDVOUT,*)'Enter length (along symmetry axis)', &
                           ' (0 to stop)'
            READ(*,*)SHPAR(1)
            IF(SHPAR(1).LE.0.)STOP
            WRITE(IDVOUT,FMT='(A,A)')'Enter 2*length of 1 hexagon side/d', &
                                     ' = (max diameter)/d'
            READ(*,*)SHPAR(2)
            WRITE(IDVOUT,*)'Specify orientation of axes a1 (hex axis) ', &
                           'and a2 (normal to one hex edge)'
            WRITE(IDVOUT,*)' 1 for a1=x,a2=y'
            WRITE(IDVOUT,*)' 2 for a1=x,a2=z'
            WRITE(IDVOUT,*)' 3 for a1=y,a2=x'
            WRITE(IDVOUT,*)' 4 for a1=y,a2=z'
            WRITE(IDVOUT,*)' 5 for a1=z,a2=x'
            WRITE(IDVOUT,*)' 6 for a1=z,a2=y'
            READ(*,*)SHPAR(3)
            SHPAR(4)=0._WP   !  P_y/d
            SHPAR(5)=0._WP   !  P_z/d
            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
            WRITE(IDVOUT,7003)(SHPAR(JJ),JJ=1,3),NAT
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='LAYRDSLAB'.OR.CSHAPE=='layrdslab')THEN
         CSHAPE='LAYRDSLAB'
         WRITE(IDVOUT,*)'LAYRDSLAB = block with 2,4, or 4 layers'
         DO J=1,100
            WRITE(IDVOUT,*)'Enter x,y,z length/d ', &
                           ' (0 to stop)'
            READ(*,*)SHPAR(1),SHPAR(2),SHPAR(3)
            IF(SHPAR(1).LE.0.)STOP
            WRITE(IDVOUT,*)'Enter fractions f1,f2,f3 (f1+f2+f3.le.1)'
            WRITE(IDVOUT,*)'      (for bilayer slab set f1+f2=1, f3=0)'
            READ(*,*)SHPAR(4),SHPAR(5),SHPAR(6)
            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
            WRITE(IDVOUT,7001)SHPAR(1),NAT
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='LYRSLBPBC'.OR.CSHAPE=='lyrslbpbc')THEN
         CSHAPE='LYRSLBPBC'
         WRITE(IDVOUT,*)'LYRSLBPBC TUC=block with 4 layers'
         DO J=1,100
            WRITE(IDVOUT,*)'Enter x,y,z length/d ', &
                           ' (0 0 0 to stop)'
            READ(*,*)SHPAR(1),SHPAR(2),SHPAR(3)
            IF(SHPAR(1).LE.0.)STOP
            WRITE(IDVOUT,*)'Enter fractions f1,f2,f3,f4 (f1+f2+f3+f4=1)'
            READ(*,*)SHPAR(4),SHPAR(5),SHPAR(6),SHPAR(7)
            SHPAR(8)=0._WP   ! P_y/d
            SHPAR(9)=0._WP   ! P_z/d
            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
            WRITE(IDVOUT,7001)SHPAR(1),NAT
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='MLTBLOCKS'.OR.CSHAPE=='mltblocks')THEN
         CSHAPE='MLTBLOCKS'
         WRITE(IDVOUT,*)'MLTBLOCKS : target defined by file blocks.par'
         CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                     SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                     MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                     IANISO,NCOMP_NEED)
         CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
         WRITE(IDVOUT,7000)NAT
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='RCTGL_PBC'.OR.CSHAPE=='rctgl_pbc')THEN
         CSHAPE='RCTGL_PBC'
         WRITE(IDVOUT,*)'Enter block dimensions X/d, Y/d, Z/d'
         READ(*,*)SHPAR(1),SHPAR(2),SHPAR(3)
         SHPAR(4)=0._WP   ! P_y/d
         SHPAR(5)=0._WP   ! P_z/d
         CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                     SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                     MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                     IANISO,NCOMP_NEED)
         CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
         WRITE(IDVOUT,7000)NAT
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='RCTGLPRSM'.OR.CSHAPE=='rctglprsm')THEN
         CSHAPE='RCTGLPRSM'
         DO J=1,100
            WRITE(IDVOUT,*)'Enter x-length (0 to stop)'
            READ(*,*)SHPAR(1)
            IF(SHPAR(1).LE.0.)STOP
            WRITE(IDVOUT,*)'Enter y-length'
            READ(*,*)SHPAR(2)
            WRITE(IDVOUT,*)'Enter z-length'
            READ(*,*)SHPAR(3)
            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
            WRITE(IDVOUT,7003)(SHPAR(JJ),JJ=1,3),NAT
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='RCTGLBLK3'.OR.CSHAPE=='rctglblk3')THEN
         CSHAPE='RCTGLBLK3'
         DO J=1,100
            WRITE(IDVOUT,*)'Enter x,y,z length/d for first block ', &
                           ' (0 0 0 to stop)'
            READ(*,*)SHPAR(1),SHPAR(2),SHPAR(3)
            IF(SHPAR(1).LE.0.)STOP
            WRITE(IDVOUT,*)'Enter x,y,z length/d for second block '
            READ(*,*)SHPAR(4),SHPAR(5),SHPAR(6)
            WRITE(IDVOUT,*)'Enter x,y,z length/d for third block '
            READ(*,*)SHPAR(7),SHPAR(8),SHPAR(9)
            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
            WRITE(IDVOUT,7001)SHPAR(1),NAT
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='RECRECPBC'.OR.CSHAPE=='recrecpbc')THEN
         CSHAPE='RECRECPBC'
         DO J=1,100
            WRITE(IDVOUT,*)'Enter x,y,z length/d for first block ', &
                           ' (0 0 0 to stop)'
            READ(*,*)SHPAR(1),SHPAR(2),SHPAR(3)
            IF(SHPAR(1).LE.0.)STOP
            WRITE(IDVOUT,*)'Enter x,y,z length/d for second block '
            READ(*,*)SHPAR(4),SHPAR(5),SHPAR(6)
            SHPAR(7)=0._WP   ! PY/d
            SHPAR(8)=0._WP   ! PZ/d
            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
            WRITE(IDVOUT,7006)SHPAR(1),SHPAR(2),SHPAR(3),SHPAR(4), &
                              SHPAR(5),SHPAR(6),NAT
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='SLAB_HOLE'.OR.CSHAPE=='slab_hole')THEN
         CSHAPE='SLAB_HOLE'
         DO J=1,100
            WRITE(IDVOUT,*)'Slab dimensions = a  b  c ; ', &
                           'cylindrical hole radius = r'
            WRITE(IDVOUT,*)'Enter a/d for block (0 to stop)'
            READ(*,*)SHPAR(1)
            IF(SHPAR(1)<0.5_WP)STOP
            WRITE(IDVOUT,*)'Enter aspect ratios b/a  c/a'
            READ(*,*)SHPAR(2),SHPAR(3)
            WRITE(IDVOUT,*)'Enter r/a for cylindrical hole'
            READ(*,*)SHPAR(4)
            IF(SHPAR(4)<0.5_WP*SQRT(SHPAR(2)**2+SHPAR(3)**2))THEN
               CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                           SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                           MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                           IANISO,NCOMP_NEED)
               CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
               WRITE(IDVOUT,7004)SHPAR(1),SHPAR(2),SHPAR(3),SHPAR(4),NAT
            ELSE
               WRITE(IDVOUT,*)'hole is larger than y-z dimensions of block'
            ENDIF
         ENDDO
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='SLBHOLPBC'.OR.CSHAPE=='slbholpbc')THEN
         CSHAPE='SLBHOLPBC'
         DO J=1,100
            WRITE(IDVOUT,*)'TUC dimensions = a  b  c ; ', &
                           'cylindrical hole radius = r'
            WRITE(IDVOUT,*)'Enter a/d for TUC (0 to stop)'
            READ(*,*)SHPAR(1)
            WRITE(IDVOUT,*)'Enter aspect ratios b/a  c/a'
            READ(*,*)SHPAR(2),SHPAR(3)
            IF(SHPAR(1)<0.5_WP)STOP
            WRITE(IDVOUT,*)'Enter r/a for cylindrical hole'
            READ(*,*)SHPAR(4)
            IF(SHPAR(4)<0.5_WP*SQRT(SHPAR(2)**2+SHPAR(3)**2))THEN
               SHPAR(5)=0._WP   ! P_y/d
               SHPAR(6)=0._WP   ! P_z/d
               CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                           SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                           MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                           IANISO,NCOMP_NEED)
               CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
               WRITE(IDVOUT,7006)SHPAR(1),SHPAR(2),SHPAR(3),  &
                                 SHPAR(4),SHPAR(5),SHPAR(6),NAT
            ELSE
               WRITE(IDVOUT,*)'hole is larger than y-z dimensions of block'
            ENDIF
         ENDDO
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='SPHERES_N'.OR.CSHAPE=='spheres_n')THEN
         CSHAPE='SPHERES_N'
         WRITE(IDVOUT,*)'Enter DIAMX = target extent in x-direction'
         READ(*,*)SHPAR(1)
         WRITE(IDVOUT,*)'Enter PRINAX=0 for a_1=(1,0,0),a_2=(0,1,0) ', &
                        'in TF'
         WRITE(IDVOUT,*)'            =1 for a_1,a_2 = principal axes'
         READ(*,*)SHPAR(2)
         WRITE(IDVOUT,*)'Enter name of sphere parameter file (e.g., ', &
                   'tarnsp.par [in quotes])'
         READ(*,*)CFLSHP
         CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                     SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                     MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                     IANISO,NCOMP_NEED)
         CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='SPHROID_2'.OR.CSHAPE=='sphroid_2')THEN
         CSHAPE='SPHROID_2'
         WRITE(IDVOUT,*)'Enter length a_1 and diameter b_1 of 1st ', &
                        'spheroid'
         READ(*,*)SHPAR(1),SHPAR(2)
         WRITE(IDVOUT,*)'Enter length a_2 and diameter b_2 of 2nd ', &
                        'spheroid'
         READ(*,*)SHPAR(3),SHPAR(4)
         WRITE(IDVOUT,*)'Enter angle phi (deg) between axes a_1 and ', &
                        'a_2'
         READ(*,*)SHPAR(5)
         WRITE(IDVOUT,*)'Enter PRINAX=0 for a_1=(1,0,0),a_2=(0,1,0) ', &
                        'in TF'
         WRITE(IDVOUT,*)'            =1 for a_1,a_2 = principal axes'
         READ(*,*)SHPAR(6)
         CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                     SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                     MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                     IANISO,NCOMP_NEED)
         CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
         WRITE(IDVOUT,7006)
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='SPHRN_PBC'.OR.CSHAPE=='sphrn_pbc')THEN
         CSHAPE='SPHRN_PBC'
         WRITE(IDVOUT,*)'Enter diam_x/d = target extent in x-direction '
         READ(*,*)SHPAR(1)
         WRITE(IDVOUT,*)'Enter filename for sphere locations (in quotes)'
         READ(*,*)CFLSHP
         SHPAR(2)=0._WP   ! Py/aeff
         SHPAR(3)=0._WP   ! PZ/aeff
         CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                     SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                     MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                     IANISO,NCOMP_NEED)
         CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
         WRITE(IDVOUT,7006)
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='SPH_ANI_N'.OR.CSHAPE=='sph_ani_n')THEN
         CSHAPE='SPH_ANI_N'
         WRITE(IDVOUT,*)'Enter PRINAX=0 for a_1=(1,0,0),a_2=(0,1,0) ', &
                        'in TF'
         WRITE(IDVOUT,*)'            =1 for a_1,a_2 = principal axes'
         READ(*,*)SHPAR(2)
         WRITE(IDVOUT,*)'Enter name of sphere parameter file (e.g., ', &
                   'tarnsp.par [in quotes])'
         READ(*,*)CFLSHP
         DO J=1,100
            WRITE(IDVOUT,*)'Enter DIAMX = target extent in x-direction', &
                           ' (0 to stop)'
            READ(*,*)SHPAR(1)
            IF(SHPAR(1).LE.0.)STOP

            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='TETRAHDRN'.OR.CSHAPE=='tetrahdrn')THEN
         CSHAPE='TETRAHDRN'
         DO J=1,100
            WRITE(IDVOUT,*)'Enter tetrahedron side length ', &
                           ' (0 to stop)'
            READ(*,*)SHPAR(1)
            SHPAR(2)=SHPAR(1)
            SHPAR(3)=SHPAR(1)
            IF(SHPAR(1).LE.0.)STOP
            CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                        SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                        MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                        IANISO,NCOMP_NEED)
            CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
            WRITE(IDVOUT,7001)SHPAR(1),NAT
         ENDDO
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='TRILYRPBC'.OR.CSHAPE=='trilyrpbc')THEN
         CSHAPE='TRILYRPBC'
         WRITE(IDVOUT,*)'TRILYRPBC TUC=three stacked rectangular blocks'
         WRITE(IDVOUT,*)'Enter size x/d, y/d, z/d for upper block'
         READ(*,*)SHPAR(1),SHPAR(2),SHPAR(3)
         WRITE(IDVOUT,*)'Enter size x/d, y/d, z/d for middle block'
         READ(*,*)SHPAR(4),SHPAR(5),SHPAR(6)
         WRITE(IDVOUT,*)'Enter size x/d, y/d, z/d for bottom block'
         READ(*,*)SHPAR(7),SHPAR(8),SHPAR(9)
         SHPAR(10)=0._WP   ! P_y/d
         SHPAR(11)=0._WP   ! P_z/d
         CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                     SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                     MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                     IANISO,NCOMP_NEED)
         CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
         STOP
!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='TRNGLPRSM'.OR.CSHAPE=='trnglprsm')THEN
         CSHAPE='TRNGLPRSM'
         WRITE(IDVOUT,*)'Enter a/d = first triangle side/d (0 to stop)'
         READ(*,*)SHPAR(1)
         IF(SHPAR(1).LE.0.)STOP
         WRITE(IDVOUT,*)'Enter b/a , c/a, L/a'
         READ(*,*)SHPAR(2),SHPAR(3),SHPAR(4)
         CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                     SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                     MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                     IANISO,NCOMP_NEED)
         CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
         STOP
!-----------------------------------------------------------------------
! Note: this target option is not described in the DDSCAT UserGuide.
!       if you are interested in using it, please contact 
!       draine@astro.princeton.edu

      ELSEIF(CSHAPE=='GAUSS_SPH'.OR.CSHAPE=='gauss_sph')THEN
         CSHAPE='GAUSS_SPH'
         WRITE(IDVOUT,*)'Enter a_eff/d'
         READ(*,*)SHPAR(1)
         WRITE(IDVOUT,*)'Enter power law index beta > 1 (0 to input', &
                   ' a_lm from file targspher.alm)'
         READ(*,*)SHPAR(3)
         IF(SHPAR(3)==0)THEN
            WRITE(IDVOUT,*)'OK -- will read input params from targspher.alm'
            CFLSHP='targspher.alm'
         ELSEIF(SHPAR(3).GT.1.)THEN
            WRITE(IDVOUT,*)'Enter Lmax > 0 (e.g., 5)'
            READ(*,*)SHPAR(4)
            WRITE(IDVOUT,*)'Enter <s^2> (e.g., 0.05) [this includes ', &
                      'all a_lm up to NMAX (=100)'
            READ(*,*)SHPAR(5)
            WRITE(IDVOUT,*)'Enter integer seed > 0'
            READ(*,*)SHPAR(6)
         ELSE
            WRITE(IDVOUT,*)'Invalid value for BETA'
         ENDIF
         WRITE(IDVOUT,*)'Enter PRINAX=0 for a_1=(1,0,0),a_2=(0,1,0) ', &
                        'in TF'
         WRITE(IDVOUT,*)'             1 for a_1,a_2 = principal axes'
         READ(*,*)SHPAR(2)
!*** diagnostic
         WRITE(IDVOUT,*)'about to call target with ioshp=',ioshp
!***
         CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT, &
                     SHPAR,DX,NAT,IXYZ,ICOMP,IDVOUT,                &
                     MXN3,NAT3,PYD,PZD,BETADF,PHIDF,THETADF,X0,     &
                     IANISO,NCOMP_NEED)
         CALL SIZER(MXNAT,NAT,IXYZ,LXYZ)
         STOP
!-----------------------------------------------------------------------
      ELSE
         WRITE(IDVOUT,*)'CALLTARGET does not know how to make shape = ',CSHAPE
         GOTO 1000
      ENDIF
 6000 FORMAT(' AX=',F7.4,' NAT=',I7)
 6100 FORMAT(' AX=',F7.4,' AY=',F7.4,' AZ=',F7.4,' NAT=',I7,' F=',F10.6)
 6200 FORMAT('    XV =',F7.4,'    YV =',F7.4,'    ZV =',F7.4,   &
             ' NAT=',I8,/,                                      &
             ' 2a_eff=',F7.4,' 2b_eff=',F7.4,' 2c_eff=',F7.4,   &
             ' err=',F8.7,/,                                    &
             '   y:x =',F7.4,'   z:x =',F7.4,'   y:z =',F7.4,/, &
             '   x:y =',F7.4,'   x:z =',F7.4,'   z:y =',F7.4)
 7000 FORMAT('NAT=',I8)
 7001 FORMAT('SHPAR(1)=',F10.5,' NAT=',I8)
 7002 FORMAT('SHPAR(1-2)=',2F10.5,' NAT=',I8)
 7003 FORMAT('SHPAR(1-3)=',3F10.5,' NAT=',I8)
 7004 FORMAT('SHPAR(1-4)=',4F10.5,' NAT=',I8)
 7005 FORMAT('SHPAR(1-5)=',5F10.5,' NAT=',I8)
 7006 FORMAT('SHPAR(1-6_=',6F10.5,' NAT=',I8)
    END PROGRAM CALLTARGET



