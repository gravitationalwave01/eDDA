    SUBROUTINE TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT,SHPAR,DX, &
                      NAT,IXYZ,ICOMP,IDVOUT,MXN3,NAT3,PYD,PZD,BETADF,PHIDF,   &
                      THETADF,X0,IANISO,NCOMP_NEED)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

!------------------------------ target_v3 ------------------------------
! Arguments:

      INTEGER :: MXNAT,MXN3,NCOMP_NEED

      CHARACTER :: CDESCR*67,CFLSHP*80,CSHAPE*9
      INTEGER*2 ::   &
         ICOMP(MXNAT,3)
      INTEGER ::     &
         IXYZ(MXNAT,3)
      INTEGER :: IANISO,IDVOUT,IDVSHP,IOSHP,NAT,NAT3
      REAL(WP) :: PYD,PZD
      REAL(WP) ::        &
         A1(3),          &
         A2(3),          &
         BETADF(MXNAT),  &
         DX(3),          &
         PHIDF(MXNAT),   &
         SHPAR(12),      &
         THETADF(MXNAT), &
         X0(3)
      EXTERNAL ERRMSG,REASHP,TAR2EL,TAR2SP,TAR3EL,TARANIREC,TARBLOCKS,TARCEL, &
         TARCYL,TARCYLCAP,TARELL,TARGSPHER,TARHEX,TARLYRSLAB,TARNAS,TARNSP,   &
         TARPBXN,TARPRSM,TARRCTBLCK3,TARREC,TARSLABHOLE,TARTET

! Local variables:

      INTEGER :: BLKSIZ,IPRINAX,J,J1,NBLOCKS
      INTEGER ::  &
         XYZB(3,100)
      REAL(WP) :: PI

!***********************************************************************
! Given:
!    CSHAPE = string describing target geometry; will determine which
!             "shape" routine is invoked.  Current options are:
!             
! ANIELLIPS, 3 args
! ANIFRMFIL, 0
! ANI_ELL_2, 3
! ANI_ELL_3, 3
! ANIRCTNGL, 3
! BISLINPBC, 6
! CONELLIPS, 6
! CYLINDER1, 3
! CYLNDRCAP, 2
! CYLNDRPBC, 5
! DSKBLYPBC, 8
! DSKRCTNGL, 5
! DSKRCTPBC, 7
! DW1996TAR, 1
! ELLIPSOID, 3
! ELLIPSO_2, 3
! ELLIPSO_3, 3
! FROM_FILE, 0
! GAUSS_SPH, 6
! HEXGONPBC, 5
! HEX_PRISM, 3
! LAYRDSLAB, 7
! LYRSLBPBC, 9
! MLTBLOCKS, 0 [read data from blocks.par]
! RCTGLPRSM, 3
! RCTGLBLK3, 9
! RCTGL_PBC, 5
! RCTG_RCTG, 6
! RECRECPBC, 8
! SLAB_HOLE, 4
! SLBHOLPBC, 6
! SPHERES_N, 2 [also need file name]
! SPHRN_PBC, 3 [also need file name]
! SPHROID_2, 6
! SPH_ANI_N, 2 [also need file name]
! TETRAHDRN, 1
! TRNGLPRSM, 4
! UNIAXICYL, 2


!    IOSHP = device number for "target.out" file
!          = -1 to suppress printing of "target.out"

!    SHPAR(1-10) = up to 10 parameters defining target geometry

!    DX(1-3)= d_x/d, d_y/d, d_z/d for lattice, where d=(d_x*d_y*d_z)**(1
!             is the effective lattice spacing, and d_x,d_y,d_z are latt
!             spacings in x,y,z directions

!    MXNAT = limit on largest allowed value of NAT
!    MXN3  = 3*MXNAT
!    IDVOUT = 1 to print out information on target shape
!            -1 to suppress printing information on target shape
!    CFLSHP = name of target shape file if CSHAPE='FROM_FILE' or
!                                          CSHAPE='ANIFRMFIL'
!    IDVSHP = device number to use to read target shape if
!             CSHAPE='FROM_FILE' or CSHAPE='ANIFRMFIL'

! Returns:

!    IANISO  = 0 for targets composed of isotropic materials
!            = 1 for targets composed of aniosotropic materials,
!                but with optical axes || x_TF,y_TF,z_TF
!            = 2 for general anisotropic targets
!    A1,A2   = two 3-vectors defining initial target orientation
!    CDESCR  = string describing target
!    NAT     = number of dipoles in target
!    NAT3    = 3*NAT
!    IXYZ(1-NAT,1-3)=integers fixing x,y,z coordinates of
!              occupied sites in target
!    ICOMP(1-NAT,1-3)=composition for E in x,y,z direction at each site
!    BETADF(1-NAT)=angle specifying orientation of "Dielectric Frame"
!                  for locations 1-NAT
!    PHIDF(1-NAT)=angle specifying orientation of "Dielectric Frame"
!                 for locations 1-NAT
!    THETADF(1-NAT)=angle specifying orientation of "Dielectric Frame"
!                   for locations 1-NAT
!    X0(1-3) = offset for IXYZ array.
!              Defined so that dipole with IXYZ(J,1-3)=(IX,IY,IZ)
!              is located at x_TF/(d*DX(1)) = IX+X0(1)
!                            y_TF/(d*DX(2)) = IY+X0(2)
!                            z_TF/(d*DX(3)) = IZ+X0(3)
!              where (x,y,z)_TF = (0,0,0) corresponds to some
!              well-defined location in the target, as specified in
!              the UserGuide for each target geometry.
!              i.e., lattice site IXYZ=(0,0,0) is located at
!                    (x,y,z)_TF = (X0(1)*DX(1), X0(2)*DX(2), X0(3)*DX(3))
!
!    PYD = 0 if isolated target
!        = (periodicity in target frame y direction)/DX(2) if periodic
!           boundary conditions are to be used
!    PZD = 0 if isolated target
!        = (periodicity in target frame z direction)/DX(3) if periodic
!           boundary conditions are to be used

! P.J.Flatau, Colorado State University
! B.T.Draine, Princeton Univ. Observatory
! History:
! 90.11.08 (BTD): Changed DO...ENDDO to DO #...# CONTINUE
! 90.12.02 (BTD): Changed ordering of data in ICOMP
! 90.12.14 (BTD): Changed option MKRECT -> RCTGLPRSM
! 91.05.02 (BTD): Added IDVOUT to argument list
!                 Added IDVOUT to argument list for subr. REASHP
! 91.05.23 (BTD): Added A1,A2 to arg. list for subr. REASHP
! 91.11.12 (BTD): Added IDVSHP to arg. list for TARGET and REASHP
! 92.09.21 (BTD): Added option UNIELL (uniaxial ellipsoid)
! 93.01.07 (BTD): Added option ELLIPSO_2 (two touching ellipsoids,
!                 consisting of materials 1 and 2)
! 92.01.07 (BTD): Added option ANI_ELL_2 (two touching anisotropic
!                 ellipsoids, 1st with diel.fnc.1,2,3 along x,y,z
!                 2nd with diel.fnc.4,5,6 along x,y,z
! 92.01.19 (BTD): Added option ELLIPSO_3 (three touching anisotropic
!                 ellipsoids, 1st with diel.fnc.1,2,3 along x,y,z
!                 2nd with diel.fnc.4,5,6 along x,y,z
!                 3rd with diel.fnc.7,8,9 along x,y,z
! 93.03.12 (BTD): Changed CDESCR*60 -> CDESCR*67
!                 Moved WRITE(IDVOUT,9010) from DDSCAT to TARGET
! 93.12.14 (BTD): Corrected handling of ELLIPSO_3 (had assigned
!                 ICOMP=1 to 50% of sites, ICOMP=2 to other 50%)
! 94.01.27 (BTD): Replaced SHPAR1,SHPAR2,SHPAR3 by SHPAR(1-6) to
!                 allow up to 6 shape parameters
!                 Added call to TARCEL to generate two concentric
!                 ellipsoids
! 94.03.21 (BTD): Changed "CONCEL" to "CONELLIPS" for consistency with
!                 REAPAR
! 95.12.11 (BTD): Modified to handle "MLTBLOCKS" as target option, using
!                 subroutine TARBLOCKS
! 96.01.04 (BTD): Modified to handle "DW1996TAR" as target option
! 96.01.25 (BTD): Modified to handle "SPHROID_2" as target option
! 96.01.26 (BTD): Replaced IX,IY,IZ by IXYZ
! 96.01.29 (BTD): Added variable IPRINAX to control choice of axes A1,A2
!                 in option MLTBLOCKS
! 97.06.04 (BTD): Corrected bug in handling of option "UNIAXICYL" -- was
!                 not setting composition correctly after changing REAPA
!                 so that no longer had choice of cylinder axis.
!                 Cylinder axis is always in x direction in target frame
!                 Thanks for Mike Wolff for calling attention to the
!                 bug.
! 97.11.01 (BTD): add DX to argument list to allow use of noncubic latti
!                 add DX to argument list of TARELL
! 97.12.26 (BTD): add DX to argument list of TARREC, TARCYL, TARHEX,
!                 TARTET, TAR2EL, TAR3EL, TARCEL, TAR2SP
! 98.04.27 (BTD): restored missing comma in argument list of TAR2EL
! 98.12.29 (BTD): added new option 'LAYRDSLAB' for layered slab
! 00.06.12 (BTD): modified to support option TARNSP for multisphere
!                 target, including variable lattice spacing DX
!                 all data type conversions now explicit
! 00.07.18 (BTD): corrected bug in calls to tarblocks
! 00.11.02 (BTD): added ICOMP to arg. list for TAR2EL,TAR3EL,TARCYL,
!                 TARELL,TARHEX,TARREC,RCTGLPRSM,TARTET
! 02.02.12 (BTD): added option TRNGLPRSM and call to routine TARPRSM
! 03.05.21 (BTD): added option SPHARM = irregular target
!                 generated using spherical harmonics
!                 target generated by TARSPHARM
! 03.05.22 (BTD): bug fix: define CFLSHP='spharm.par' when calling
!                 TARSPHARM
! 03.11.06 (BTD): modify to use new version of TARSPHARM
! 04.03.19 (BTD): change CFLSHP*13 to CFLSHP*80 to accomodate long file
!                 names
! 04.04.01 (BTD): added arguments NPY,NPZ to support periodic b.c.
!                 option (at present time, used only for target option
!                 SPHRN_PBC)
! 04.04.01 (BTD): Add DX to argument list of TARGSPHER
!                 Remove call to TARSPHARM -- this is superseded by
!                 TARGSPHER
! 04.04.29 (BTD): Add option ANIRCTNGL for anisotropic rectangular target
! Begin DDSCAT 6.2:
! 04.09.14 (BTD): Add material orientation angles BETADF,PHIDF,THETADF
!                 to argument list.
!                 Add target option ANIFRMFIL to read properties of general
!                 anisotropic target from file.
!                 Add target option SPH_ANI_N to read properties of target
!                 consisting of N anisotropic spheres
! 05.06.16 (BTD): Change from integers NYD,NZD to real PYD,PZD, because
!                 there is no reason for periods to be integer multiples
!                 of lattice spacing d.
! 05.08.03 (BTD): Add new target option RCTGL_PBC
! 06.09.13 (BTD): Add new target option CYLNDRCAP
! 06.09.15 (BTD): Modify to pass additional variable to TARHEX to allow
!                 orientation in TF to be selected.
! 06.09.15 (BTD): Add new target option HEXGONPBC
! 06.09.15 (BTD): Modify to pass additional variable to TARCYL to allow
!                 orientation in TF to be selected.
! 06.09.15 (BTD): Add new target option CYLNDRPBC
! 06.10.23 (BTD): Set a1=(1,0,0),a2=(0,1,0) for options CYLNDRPBC,HEXGONPBC
! 06.12.08 (BTD): Added support for target option TARPBX
!                 (pillbox on slab)
! 07.01.18 (BTD): Fixed bug in assignment of value of PYD,PZD
! 07.01.20 (BTD): Added support for option DSKRCTGNL
! 07.02.23 (BTD): Added support for option SLBPBC
! 07.02.25 (BTD): Fixed bugs.
! 07.06.19 (BTD): Add X0 (target offset info) to argument list
!                 Modified to use new TARPBX with X0
! 07.06.20 (BTD): Add X0 to argument list for TARCYL
! 07.09.09 (BTD): Cosmetic changes.
!                 Make sure that all PBC targets have
!                 a1=x_TF, a_2=y_TF
! 07.09.11 (BTD): Changed IXYZ from INTEGER*2 to INTEGER
! 07.10.26 (BTD): Replace call to TARPBX by call to TARPBXN
!                 Add support for target option DSKBLYPBC
! 07.10.27 (BTD): * Changed SHPAR(6) -> SHPAR(10)
!                 * Changed all target names to new set of 9 character
!                   names.
!                 * Rorganized calls to target generating routines.
! 08.01.13 (BTD): Modified to add X0 to argument list of almost all
!                 target generation routines, as well as REASHP
! 08.01.17 (BTD): Added IANISO to argument list
!                 add code to set IANISO for each target option
! 08.02.01 (BTD): Changed SHPAR(10) -> SHPAR(12)
!                 Added target options 
!                 * RCTGLBLK3  (isolated target: 3 rectangular blocks)
!                 * BISLINPBC  (periodic target: bilayer slab with 
!                               parallel lines on top)
!                 * TRILYRPBC  (periodic target: 3 layer rectangular structure)
! 08.03.13 (BTD): v7.0.5
!                 corrected bug in handling of BISLINPBC
! 08.04.24 (BTD): Add new target option
!                 * TARHOLSLB  (slab with circular hole)
!                   using routine TARHOLESLAB
! 08.08.08 (BTD): Added X0 to argument list of TARBLOCKS
!                 Added X0 to argument list of TAR2SP
!                 Added X0 to argument list of TARHEX
!                 Added X0 to argument list of TARNAS
! 08.08.30 (BTD): Added X0 to argument list of TARPRSM
!                 Added X0 to argument list of TARTET
! 08.09.12 (BTD): Corrected typos DSKRCTGNL -> DSKRCTNGL
! 08.10.28 (BTD): Added X0 to argument list of TARBLOCKS for target
!                 option MLTBLOCKS
! 09.02.02 (BTD): ver7.0.8
!                 For CSHAPE=='FROM_FILE', added code to set IANISO=2
!                 if any of the angles BETADF,PHIDF,THETADF are nonzero
! 09.09.10 (BTD): ver7.1.0
!                 Removed previous change: should use ANIFRMFIL if
!                 dielectric tensor may be nondiagonal in TF
!                 Added support for options
!                 ANIFILPBC (TUC from file, dielectric tensor may be
!                            nondiagonal in TF -- read angles
!                            BETADF,THETADF,PHIDF)
!                 FRMFILPBC (TUC from file, dielectric tensor assumed
!                            to be diagonal in TF)
! 09.09.11 (BTD): Added NCOMP_NEED to argument list
!                 and in calls to REASHP (to allow checking to see if
!                 enough refractive index files hav been provided)
! 09.09.17 (BTD): initialize NCOMP_NEED for when REASHP is not called
! 10.02.01 (BTD): target_v3
!                 add new target options
!                 * SLAB_HOLE
!                 * SLBHOLPBC
! 10.02.06 (BTD): * corrected handling of option BISLINPBC: had incorrect
!                   order of aguments of tarrctblck3
!                 * discontinue target option TARSLBLIN, as this target
!                   geometry can already be generated with RCTGLBLK3
! 10.03.03 (BTD): * corrected handling of option FRMFILPBC
!                   (had failed to set PYD and PZD)
!                 * corrected handling of option ANIFILPBC
!                   (had failed to set PYD and PZD)
! end history

! Copyright (C) 1993,1994,1995,1996,1997,1998,2000,2002,2003,2004,2005
!               2006,2007,2008,2009,2010 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
!*** diagnostic
!      write(0,*)'entered target, CSHAPE=',CSHAPE
!***
! initialize X0 to zero, because some target routines have not yet
! been upgraded to initialize X0

      DO J=1,3
         X0(J)=0._WP
      ENDDO

! initialize NCOMP_NEED to 1 (will be reset if REASHP is called)

      NCOMP_NEED=1

! Most targets are isolated, so set PYD=PZD=0 as the default, to be
! changed only in case of targets with periodic b.c.

      PYD=0._WP
      PZD=0._WP

! Most targets have BETADF=0,PHIDF=0,THETADF=0 so initialize these
! before proceeding.

      DO J=1,MXNAT
         BETADF(J)=0._WP
      ENDDO
      DO J=1,MXNAT
         PHIDF(J)=0._WP
      ENDDO
      DO J=1,MXNAT
         THETADF(J)=0._WP
      ENDDO

! Most targets use isotropic materials, so initialize IANISO=0
! before proceeding

      IANISO=0

!*** diagnostic
!      write(0,*)'target ckpt 1'
!***
!-----------------------------------------------------------------------
! ANIELLIPS -> anisotropic ellipsoid

      IF(CSHAPE=='ANIELLIPS')THEN
         CALL TARELL(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),DX,X0,CDESCR,IOSHP, &
                     MXNAT,NAT,IXYZ,ICOMP)

! Dielectric constants 1,2,3 for E along directions x,y,z in Target Fram

         IANISO=1
         DO J=1,NAT
            ICOMP(J,1)=1
            ICOMP(J,2)=2
            ICOMP(J,3)=3
         ENDDO

! ANIFILPBC -> for each dipole j, read one line from CFLSHP with
!              IXYZ(j,1-3),ICOMP(j,1-3),BETADF(j),THETADF(j),PHIDF(j)
!              where 
!              IXYZ(j,1-3) specifies location of dipole j in lattice units
!              ICOMP(j,1-3) specifies compositions corresponding to
!                 3 principal axes of dielectric tensor at j
!              BETADF(j),THETADF(j),PHIDF(j) specifies orientation of
!                 dielectric tensor at location j

      ELSEIF(CSHAPE=='ANIFILPBC')THEN

         CALL REASHP(CSHAPE,CFLSHP,IDVOUT,IDVSHP,A1,A2,DX,X0,BETADF,PHIDF,   &
                     THETADF,CDESCR,MXNAT,NAT,IXYZ,ICOMP,MXN3,NAT3,NCOMP_NEED)

         IANISO=1
         DO J=1,NAT
            IF(BETADF(J)/=0.)IANISO=2
            IF(PHIDF(J)/=0.)IANISO=2
            IF(THETADF(J)/=0.)IANISO=2            
         ENDDO
         PYD=SHPAR(1)
         PZD=SHPAR(2)

!-----------------------------------------------------------------------
! ANIFRMFIL -> read list of occupied sites and composition and material
!           orientation from file
!           file CFLSHP must specify, for each dipole j=1-NAT
!           IXYZ(j,1-3),ICOMP(j,1-3),betadf(j),phidf(j),thetadf(j)

      ELSEIF(CSHAPE=='ANIFRMFIL')THEN

         CALL REASHP(CSHAPE,CFLSHP,IDVOUT,IDVSHP,A1,A2,DX,X0,BETADF,PHIDF,   &
                     THETADF,CDESCR,MXNAT,NAT,IXYZ,ICOMP,MXN3,NAT3,NCOMP_NEED)
         IANISO=1
         DO J=1,NAT
            IF(BETADF(J)/=0.)IANISO=2
            IF(PHIDF(J)/=0.)IANISO=2
            IF(THETADF(J)/=0.)IANISO=2
         ENDDO

!-----------------------------------------------------------------------
! ANI_ELL_2 -> two touching identical ellipsoids with anisotropic
!           functions.
!           second ellipsoid is displaced from first in x-direction
!           1st ellipsoid has dielectric fcns 1,2,3 in x,y,z dirs.
!           2nd ellipsoid has dielectric fcns 4,5,6 in x,y,z dirs.

      ELSEIF(CSHAPE=='ANI_ELL_2')THEN
         CALL TAR2EL(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),DX,X0,CDESCR,IOSHP, &
                     MXNAT,NAT,IXYZ,ICOMP)
         IANISO=1
         J1=NAT/2
         DO J=1,J1
            ICOMP(J,1)=1
            ICOMP(J,2)=2
            ICOMP(J,3)=3
         ENDDO
         DO J=J1+1,NAT
            ICOMP(J,1)=4
            ICOMP(J,2)=5
            ICOMP(J,3)=6
         ENDDO

!-----------------------------------------------------------------------
! ANI_ELL_3 -> three touching identical ellipsoids with anisotropic
!           functions.  Need to override TAR3EL.
!           second ellipsoid is displaced from first in x-direction
!           1st ellipsoid has dielectric fcns 1,2,3 in x,y,z dirs.
!           2nd ellipsoid has dielectric fcns 4,5,6 in x,y,z dirs.
!           3rd ellipsoid has dielectric fcns 7,8.9 in x,y,z dirs.

      ELSEIF(CSHAPE=='ANI_ELL_3')THEN
         CALL TAR3EL(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),DX,X0,CDESCR,IOSHP, &
                     MXNAT,NAT,IXYZ,ICOMP)
         IANISO=1
         J1=NAT/3
         DO J=1,J1
            ICOMP(J,1)=1
            ICOMP(J,2)=2
            ICOMP(J,3)=3
         ENDDO
         DO J=J1+1,2*J1
            ICOMP(J,1)=4
            ICOMP(J,2)=5
            ICOMP(J,3)=6
         ENDDO
         DO J=2*J1+1,NAT
            ICOMP(J,1)=7
            ICOMP(J,2)=8
            ICOMP(J,3)=9
         ENDDO

!-----------------------------------------------------------------------
! ANIRCTNGL -> homogeneous, anisotropic, rectangular target
!              optical axes of material are parallel to x_TF,y_TF,z_TF

      ELSEIF(CSHAPE=='ANIRCTNGL')THEN
         CALL TARANIREC(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),DX,X0,CDESCR,IOSHP, &
                        MXNAT,NAT,IXYZ,ICOMP)
         IANISO=1

!-----------------------------------------------------------------------
! BISLINPBC -> bilayer slab with a rectangular-cross-section line on top,
!              extending infinitely in z-direction
!              if SHPAR(5)=0, then finite in y direction
!              SHPAR(1) = x-thickness/d of line  [material 1]
!              SHPAR(2) = y-width/d of line
!              SHPAR(3) = x-thickness/d of upper slab layer [material 2]
!              SHPAR(4) = x-thickness/d of lower slab layer [material 3]
!              SHPAR(5) = y-width/d of both lower slab layers
!              SHPAR(6) = periodicity/d of structure in y direction

      ELSEIF(CSHAPE=='BISLINPBC')THEN

! 2010.02.06 (BTD) correction
!         CALL TARRCTBLK3(A1,A2,SHPAR(1),SHPAR(3),SHPAR(4),SHPAR(2),SHPAR(5), &
!                         SHPAR(5),1._WP,1._WP,1._WP,DX,X0,CDESCR,IOSHP,      &
!                         MXNAT,NAT,IXYZ,ICOMP)
         CALL TARRCTBLK3(A1,A2,SHPAR(1),SHPAR(2),1._WP,SHPAR(3),SHPAR(5),  &
                         1._WP,SHPAR(4),SHPAR(5),1._WP,DX,X0,CDESCR,IOSHP, &
                         MXNAT,NAT,IXYZ,ICOMP)
         IANISO=0
         PYD=SHPAR(6)
         PZD=1._WP

!-----------------------------------------------------------------------
! CONELLIPS -> two concentric ellipsoids with isotropic dielectric function
!           1 (inner) and 2 (outer)

      ELSEIF(CSHAPE=='CONELLIPS')THEN
         CALL TARCEL(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),SHPAR(4),SHPAR(5), &
                     SHPAR(6),DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
         IANISO=0

!-----------------------------------------------------------------------
! CYLINDER1 -> homogeneous, isotropic cylinder

      ELSEIF(CSHAPE=='CYLINDER1')THEN
         CALL TARCYL(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),DX,X0,CDESCR,IOSHP, &
                     MXNAT,NAT,IXYZ,ICOMP)
         IANISO=0

!-----------------------------------------------------------------------
! CYLNDRCAP -> homogeneous, isotropic cylinder with hemispherical caps

      ELSEIF(CSHAPE=='CYLNDRCAP')THEN
         CALL TARCYLCAP(A1,A2,SHPAR(1),SHPAR(2),DX,X0,CDESCR,IOSHP,MXNAT,NAT, &
                        IXYZ,ICOMP)
         IANISO=0

!-----------------------------------------------------------------------
! CYLNDRPBC -> homogeneous, isotropic cylinder, repeated with
!           periodicity SHPAR(4) in y direction,
!                       SHPAR(5) in z direction

      ELSEIF(CSHAPE=='CYLNDRPBC')THEN
         CALL TARCYL(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),DX,X0,CDESCR,IOSHP, &
                     MXNAT,NAT,IXYZ,ICOMP)

! set target axis a1 = (1,0,0) in TF
!                 a2 = (0,1,0) in TF

         DO J=1,3
            A1(J)=0._WP
            A2(J)=0._WP
         ENDDO
         A1(1)=1._WP
         A2(2)=1._WP
         PYD=SHPAR(4)
         PZD=SHPAR(5)
         IANISO=0

!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='DSKBLYPBC')THEN

! SLPBX2 -> disk on bilayer slab (isolated target)
!           SHPAR(1) = (disk thickness)/d  [material 1]
!           SHPAR(2) = (disk diameter)/d               
!           SHPAR(3) = (x thickness of slab material 2)/d
!           SHPAR(4) = (x thickness of slab material 3)/d
!           SHPAR(5) = (y length of slab)/d
!           SHPAR(6) = (z length of slab)/d
!           SHPAR(7) = (y-periodicity)/d   [SHPAR(8).GE.SHPAR(5)]
!           SHPAR(8) = (z-periodicity)/d   [SHPAR(8).GE.SPHAR(6)]

! set target axis a1 = (1,0,0) in TF
!                 a2 = (0,1,0) in TF

         CALL TARPBXN(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),SHPAR(4),SHPAR(5), &
                     SHPAR(6),DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
         IANISO=0
         PYD=SHPAR(7)
         PZD=SHPAR(8)

!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='DSKRCTNGL')THEN

! DSKRCTNGL -> disk on homogeneous slab (isolated target)
!           SHPAR(1) = (disk thickness)/d  [material 1]
!           SHPAR(2) = (disk diameter)/d               
!           SHPAR(3) = (x thickness of slab)/d [material 2]
!           SHPAR(4) = (y length of slab)/d
!           SHPAR(5) = (z length of slab)/d

! set target axis a1 = (1,0,0) in TF
!                 a2 = (0,1,0) in TF

         CALL TARPBXN(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),0._WP,SHPAR(4), &
                     SHPAR(5),DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
         IANISO=0

!!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='DSKRCTPBC')THEN

! DSKRCTPBC -> pillbox on slab with periodic boundary conditions
!           SHPAR(1) = (disk thickness)/d    [material 1
!           SHPAR(2) = (disk diameter)/d
!           SHPAR(3) = (x thickness of slab)/d
!           SHPAR(4) = (y length of slab)/d
!           SHPAR(5) = (z length of slab)/d
!           SHPAR(6) = (repeat length in y direction)/d
!           SHPAR(7) = (repeat length in z direction)/d

! set target axis a1 = (1,0,0) in TF
!                 a2 = (0,1,0) in TF

         CALL TARPBXN(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),0._WP,SHPAR(4), &
                      SHPAR(5),DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
         IANISO=0
         PYD=SHPAR(6)
         PZD=SHPAR(7)

!-----------------------------------------------------------------------
! DW1996TAR -> the 13-cube target studied by Draine & Weingartner 1996
!           SHPAR(1)=width per block in lattice units

      ELSEIF(CSHAPE=='DW1996TAR')THEN
         IPRINAX=1
         NBLOCKS=13
         BLKSIZ=NINT(SHPAR(1))
         XYZB(1,1)=0
         XYZB(2,1)=1
         XYZB(3,1)=0
         XYZB(1,2)=1
         XYZB(2,2)=1
         XYZB(3,2)=0
         XYZB(1,3)=0
         XYZB(2,3)=2
         XYZB(3,3)=0
         XYZB(1,4)=1
         XYZB(2,4)=2
         XYZB(3,4)=0
         XYZB(1,5)=0
         XYZB(2,5)=1
         XYZB(3,5)=1
         XYZB(1,6)=1
         XYZB(2,6)=1
         XYZB(3,6)=1
         XYZB(1,7)=0
         XYZB(2,7)=2
         XYZB(3,7)=1
         XYZB(1,8)=1
         XYZB(2,8)=2
         XYZB(3,8)=1
         XYZB(1,9)=2
         XYZB(2,9)=1
         XYZB(3,9)=0
         XYZB(1,10)=2
         XYZB(2,10)=2
         XYZB(3,10)=0
         XYZB(1,11)=0
         XYZB(2,11)=0
         XYZB(3,11)=1
         XYZB(1,12)=0
         XYZB(2,12)=0
         XYZB(3,12)=2
         XYZB(1,13)=0
         XYZB(2,13)=1
         XYZB(3,13)=2
! diagnostic
!         write(0,*)'target ckpt 10, blksiz=',blksiz
!
         CALL TARBLOCKS(A1,A2,DX,NBLOCKS,BLKSIZ,XYZB,X0,IPRINAX,IOSHP, &
                        CDESCR,MXNAT,NAT,IXYZ,ICOMP)
! diagnostic
!         write(0,*)'target ckpt 11, nat=',nat
!
         IANISO=0
         CDESCR='13-cube target used by Draine & Weingartner 1996'

!-----------------------------------------------------------------------
! ELLIPSOID -> homogeneous, isotropic, ellipsoidal target

      ELSEIF(CSHAPE=='ELLIPSOID')THEN
!*** diagnostic
!         write(0,*)'target ckpt 20'
!***
         CALL TARELL(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),DX,X0,CDESCR,IOSHP, &
                     MXNAT,NAT,IXYZ,ICOMP)
!*** diagnostic
!         write(0,*)'target ckpt 21'
!***
         IANISO=0

!-----------------------------------------------------------------------
! ELLIPSO_2 -> two touching identical ellipsoids, of materials 1 and 2
!           second ellipsoid is displaced from first in x-direction

      ELSEIF(CSHAPE=='ELLIPSO_2')THEN
         CALL TAR2EL(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),DX,X0,CDESCR,IOSHP, &
                     MXNAT,NAT,IXYZ,ICOMP)
         IANISO=0

!-----------------------------------------------------------------------
! ELLIPSO_3 -> three touching identical ellipsoids, of materials 1,2,3
!           second ellipsoid is displaced from first in +x-direction
!           third ellipsoid is displaced from second in +x-direction
!           1st,2nd,3rd ellipsoid has ICOMP=1,2,3 (set in TAR3EL)

      ELSEIF(CSHAPE=='ELLIPSO_3')THEN
         CALL TAR3EL(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),DX,X0,CDESCR,IOSHP, &
                     MXNAT,NAT,IXYZ,ICOMP)
         IANISO=0

!-----------------------------------------------------------------------
! FROM_FILE -> read list of locations IXYZ(J,1-3) and compositions 
!              ICOMP(J,1-3) for occupied sites J=1-NAT
!              from file CFLSHP
!              REASHP should leave BETADF,PHIDF,THETADF untouched
!              If target is anisotropic, principal axes of dielectric tensor
!              must be aligned with x,y,z axes in target frame.

      ELSEIF(CSHAPE=='FROM_FILE')THEN
         CALL REASHP(CSHAPE,CFLSHP,IDVOUT,IDVSHP,A1,A2,DX,X0,BETADF,PHIDF,   &
                     THETADF,CDESCR,MXNAT,NAT,IXYZ,ICOMP,MXN3,NAT3,NCOMP_NEED)
         IANISO=1

!-----------------------------------------------------------------------
! FRMFILPBC -> read list of locations IXYZ(J,1-3) and compositions
!              ICOMP(J,1-3) for occupied sites J=1-NAT in TUC
!              from file CFLSHP
!              REASHP should leave BETADF,PHIDF,THETADF unchanged
!              If target is anisotropic, principal axes of dielectric
!              tensor must be aligned with x,y,z axes in target frame.

      ELSEIF(CSHAPE=='FRMFILPBC')THEN
         CALL REASHP(CSHAPE,CFLSHP,IDVOUT,IDVSHP,A1,A2,DX,X0,BETADF,PHIDF,   &
                     THETADF,CDESCR,MXNAT,NAT,IXYZ,ICOMP,MXN3,NAT3,NCOMP_NEED)
         IANISO=1
         PYD=SHPAR(1)
         PZD=SHPAR(2)

!-----------------------------------------------------------------------
! GAUSS_SPH -> gaussian sphere target
! SHPAR(1)= a_eff/d
! SHPAR(2)= 0. to use A1=(1,0,0),A2=(0,1,0) in TF
!           1. to use A1,A2 = principal axes with largest,2nd largest
!              moment of inertia
! SHPAR(3)= BETA : for L>1, contribution of different L values to
!                  <s^2> varies as L^{-BETA}
! SHPAR(4)= LMAX : maximum value of L to use
! SHPAR(5)= <s^2>
! SHPAR(6)= RSEED

      ELSEIF(CSHAPE=='GAUSS_SPH')THEN
!*** diagnostic
!         WRITE(0,*)'about to call targspher with shpar(1)=',SHPAR(1)
!***
         CALL TARGSPHER(SHPAR(1),SHPAR(2),SHPAR(3),SHPAR(4),SHPAR(5)   ,    &
                        SHPAR(6),A1,A2,DX,X0,CFLSHP,CDESCR,IOSHP,MXNAT,NAT, &
                        IXYZ,ICOMP)
         IANISO=0

!-----------------------------------------------------------------------
! HEXGONPBC -> homogeneous, isotropic hexagonal prism, repeated with
!           periodicity SHPAR(4) in y direction,
!                       SHPAR(5) in z direction

      ELSEIF(CSHAPE=='HEXGONPBC')THEN
         CALL TARHEX(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),DX,X0,CDESCR,IOSHP, &
                     MXNAT,NAT,IXYZ,ICOMP)

! set target axis a1 = (1,0,0) in TF
!                 a2 = (0,1,0) in TF

         DO J=1,3
            A1(J)=0._WP
            A2(J)=0._WP
         ENDDO
         A1(1)=1._WP
         A2(2)=1._WP

         IANISO=0
         PYD=SHPAR(4)
         PZD=SHPAR(5)

!-----------------------------------------------------------------------
! HEX_PRISM -> homogeneous, isotropic hexagonal prism

      ELSEIF(CSHAPE=='HEX_PRISM')THEN
         CALL TARHEX(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),DX,X0,CDESCR,IOSHP, &
                     MXNAT,NAT,IXYZ,ICOMP)
         IANISO=0

!-----------------------------------------------------------------------
! LAYRDSLAB -> Layered slab:

! SHPAR(1) = (x-thickness of multilayer slab)/[d*DX(1)]
! SHPAR(2) = (y-extent of slab)/[d*DX(2)]
! SHPAR(3) = (z-extent of slab)/[d*DX(3)]
! SHPAR(4) = fraction f1 of slab occupied by material 1
! SHPAR(5) = fraction f2 of slab occupied by material 2
! SHPAR(6) = fraction f3 of slab occupied by material 3
! SHPAR(7) = fraction f4 of slab occupied by material 4

      ELSEIF(CSHAPE=='LAYRDSLAB')THEN
         CALL TARLYRSLAB(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),SHPAR(4),SHPAR(5),  &
                         SHPAR(6),SHPAR(7),DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ, &
                         ICOMP)
         IANISO=0

! slab has dimensions a,b,c in x,y,z dimensions
! current version allows for up to 4 layers
! slab is layered in x direction
! top surface is composition 1, with interface at x_TF=0
! lower surface has boundary at x_TF=-SHPAR(1)*d

! for one layer only, set f2=f3=f4=0
! for two layers only, set f3=f4=0
! for three layers only, set f4=0

! SHPAR(1-3)= a/d , b/d , c/d  where d=(dx*dy*dz)**{1/3}
! SHPAR(4)  = x1/a
! SHPAR(5)  = x2/a
! SHPAR(6)  = x3/a
!
!-----------------------------------------------------------------------
! LYRSLBPBC -> Multilayer brick repeated in y and z directions

! SHPAR(1) = (x-thickness of multilayer slab)/[d*DX(1)]
! SHPAR(2) = (y-extent of slab)/[d*DX(2)]
! SHPAR(3) = (z-extent of slab)/[d*DX(3)]
! SHPAR(4) = fraction of slab occupied by material 1 (top layer)
! SHPAR(5) = fraction of slab occupied by material 2
! SHPAR(6) = fraction of slab occupied by material 3
! SHPAR(7) = fraction of slab occupied by material 4
! SHPAR(8) = (repeat length in y direction)/[d*DX(2)]
! SHPAR(9) = (repeat length in z direction)/[d*DX(3)]

      ELSEIF(CSHAPE=='LYRSLBPBC')THEN
         CALL TARLYRSLAB(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),SHPAR(4),SHPAR(5),  &
                         SHPAR(6),SHPAR(7),DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ, &
                         ICOMP)

! sets target axis a1 = (1,0,0) in TF
!                  a2 = (0,1,0) in TF
! top surface of brick has x_TF = 0 
! X0(1-3) is at center of top surface.
! lower surface is at x_TF=-SHPAR1*d*DX(1)

         IANISO=0
         PYD=SHPAR(8)
         PZD=SHPAR(9)

!-----------------------------------------------------------------------
! MLTBLOCKS -> target is constructed from identical cubic blocks
!           locations of blocks, and size of blocks (lattice units)
!           is specified in file 'blocks.par'
!           file structure:
!              target description
!              IPRINAX = 0 (a1=100,a2=010) or 1 (a1,a2=principal axes)
!              N = number of blocks
!              SIZE = width of one block, in lattice units
!              X Y Z = location of block 1 (in units of block width)
!              ...
!              X Y Z =    "     "    "   N  "   "    "    "     "

      ELSEIF(CSHAPE=='MLTBLOCKS')THEN
         OPEN(UNIT=IDVSHP,FILE='blocks.par',STATUS='OLD',FORM='FORMATTED')
         READ(IDVSHP,FMT='(A67)')CDESCR
         READ(IDVSHP,*)IPRINAX
         READ(IDVSHP,*)NBLOCKS
         READ(IDVSHP,*)BLKSIZ
         DO J=1,NBLOCKS
            READ(IDVSHP,*)XYZB(1,J),XYZB(2,J),XYZB(3,J)
         ENDDO
         CLOSE(IDVSHP)
         CALL TARBLOCKS(A1,A2,DX,NBLOCKS,BLKSIZ,XYZB,X0,IPRINAX,IOSHP,CDESCR, &
                        MXNAT,NAT,IXYZ,ICOMP)
         IANISO=0

!-----------------------------------------------------------------------
! RCTGLPRSM -> homogeneous, isotropic, rectangular target

      ELSEIF(CSHAPE=='RCTGLPRSM')THEN
         CALL TARREC(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),DX,X0,CDESCR,IOSHP, &
                     MXNAT,NAT,IXYZ,ICOMP)
         IANISO=0

!-----------------------------------------------------------------------
! RCTGL_PBC -> rectangular brick with periodic boundary conditions
!           with A1=(1,0,0), A2=(0,1,0) in target frame
! SHPAR(1) = extent of brick in target frame x direction/d
! SHPAR(2) = extent of brick in target frame y direction)/d
! SHPAR(3) = extent of brick in target frame z direction)/d
! SHPAR(4) = periodicity in target frame y direction/d
!            [SHPAR(4) .GE. SHPAR(2)]
! SHPAR(5) = periodicity in target frame z direction/d
!            [SHPAR(5) .GE. SHPAR(3)]

      ELSEIF(CSHAPE=='RCTGL_PBC')THEN
         CALL TARREC(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),DX,X0,CDESCR,IOSHP, &
                    MXNAT,NAT,IXYZ,ICOMP)

! note: TARREC sets target axis a1 = (1,0,0) in TF
!                               a2 = (0,1,0) in TF
         IANISO=0
         PYD=SHPAR(4)/DX(2)
         PZD=SHPAR(5)/DX(3)

!-----------------------------------------------------------------------
! RCTGLBLK3 -> stack of 3 rectangular solids (isolated target)
!              x axis passes through centers
!              each slab consists of isotropic material
!           SHPAR(1) = (upper solid thickness in x direction)/d [material 1]
!           SHPAR(2) = (upper solid width in y direction)/d
!           SHPAR(3) = (upper solid width in z direction)/d
!           SHPAR(4) = (middle solid thickness in x direction)/d [material 2]
!           SHPAR(5) = (middle solid width in y direction)/d
!           SHPAR(6) = (middle solid width in z direction)/d
!           SHPAR(7) = (lower solid thickness in x direction)/d [material 3]
!           SHPAR(8) = (lower solid width in y direction)/d
!           SHPAR(9) = (lower solid width in z direction)/d

         ELSEIF(CSHAPE=='RCTGLBLK3')THEN
            CALL TARRCTBLK3(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),SHPAR(4),       &
                            SHPAR(5),SHPAR(6),SHPAR(7),SHPAR(8),SHPAR(9),DX, &
                            X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
            IANISO=0
         
!-----------------------------------------------------------------------
! RCTG_RCTG -> rectangular solid on top resting on another rectangular solid
!              (isolated target)
!           SHPAR(1) = (upper solid thickness in x direction)/d  [material 1]
!           SHPAR(2) = (upper solid width in y direction)/d
!           SHPAR(3) = (upper solid width in z direction)/d
!           SHPAR(4) = (x thickness of lower solid)/d [material 2]
!           SHPAR(5) = (y length of lower solid)/d
!           SHPAR(6) = (z length of lower solid)/d

!     target axis a1 = (1,0,0) in TF
!                 a2 = (0,1,0) in TF

      ELSEIF(CSHAPE=='RCTG_RCTG')THEN

         CALL TARRECREC(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),SHPAR(4),SHPAR(5), &
                        SHPAR(6),DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
         IANISO=0

!-----------------------------------------------------------------------
      ELSEIF(CSHAPE=='RECRECPBC')THEN

! RECRECPBC -> rectangular solid on top resting on another rectangular solid
!              repeated periodically in y and z directions
!           SHPAR(1) = (upper solid thickness in x direction)/d  [material 1]
!           SHPAR(2) = (upper solid width in y direction)/d
!           SHPAR(3) = (upper solid width in z direction)/d
!           SHPAR(4) = (x thickness of lower solid)/d [material 2]
!           SHPAR(5) = (y length of lower solid)/d
!           SHPAR(6) = (z length of lower solid)/d
!           SHPAR(7) = (repeat period in y direction)/d 
!           SHPAR(8) = (repeat period in z direction)/d
!                      NB: must have SHPAR(7).GE.MAX(SHPAR(2),SHPAR(5))
!                                    SHPAR(8).GE.MAX(SHPAR(3),SHPAR(6))
!     target axis a1 = (1,0,0) in TF
!                 a2 = (0,1,0) in TF
!*** diagnostic
!         write(0,*)'in target, about to call tarrecrec...'
!***
         CALL TARRECREC(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),SHPAR(4),SHPAR(5), &
                        SHPAR(6),DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
         IANISO=0
         PYD=SHPAR(7)
         PZD=SHPAR(8)

!-----------------------------------------------------------------------
! SLAB_HOLE -> rectangular slab with a cylindrical hole through center
!              axis of cylindrical hole is in x direction
!              slab extent is (a,b,c) in (x,y,z)_TF directions
!              cylindrical hole radius is r
!            SHPAR(1) = a = (thickness of slab in x direction)/d
!            SHPAR(2) = b/a 
!            SHPAR(3) = c/a
!            SHPAR(4) = r/a

      ELSEIF(CSHAPE=='SLAB_HOLE')THEN
!*** diagnostic
!         write(0,*)'target ckpt 30, target option SLAB_HOLE'
!
         CALL TARSLABHOLE(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),SHPAR(4), &
                          DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
         IANISO=0

!-----------------------------------------------------------------------
! SLBHOLPBC -> 1-d or 2-d array of rectangular blocks with cylindrical hole
!              through center of each block
!              axis of cylindrical hole is in x-direction
!              blocks are organized in 2-d lattice with
!              periodicity P_y and P_z in y- and z-directions
!              each block has extent (a,b,c) in (x,y,z)_TF directions
!              cylindrical hole radius = r
!            SHPAR(1) = a/d
!            SHPAR(2) = b/a
!            SHPAR(3) = c/a
!            SHPAR(4) = r/a
!            SHPAR(5) = P_y/d
!            SHPAR(6) = P_z/d
!            if P_y=0 or P_z=0, array is 1-d
!            if both P_y and P_z are nonzero, array is 2-d
!            if P_y>0, it is required that P_y>a
!            if P_z>0, it is required that P_z>b

      ELSEIF(CSHAPE=='SLBHOLPBC')THEN
! diagnostic
!         write(0,*)'target ckpt 40 SLBHOLPBC'
!
         CALL TARSLABHOLE(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),SHPAR(4), &
                          DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
         IANISO=0
         PYD=SHPAR(5)
         PZD=SHPAR(6)

! require that either P_y=0 or P_y > SHPAR(2)*SHPAR(1)

         IF(PYD/(SHPAR(1)*SHPAR(2))<=0.9999999_WP)THEN
            CALL ERRMSG('FATAL','TARGET','0 < P_y < a')
         ELSEIF(PZD/(SHPAR(1)*SHPAR(3))<=0.9999999_WP)THEN
            CALL ERRMSG('FATAL','TARGET','0 < P_z < b')
         ENDIF

!-----------------------------------------------------------------------
! SPHERES_N -> Multisphere target (uniform composition):

      ELSEIF(CSHAPE=='SPHERES_N')THEN
! diagnostic
!         write(0,*)'target ckpt 50 spheres_n'
!
         CALL TARNSP(A1,A2,SHPAR(1),SHPAR(2),DX,X0,CFLSHP,CDESCR,IDVSHP, &
                     IOSHP,MXNAT,NAT,IXYZ,ICOMP)
         IANISO=0

!-----------------------------------------------------------------------
! SPHRN_PBC -> cluster of N spheres with periodic boundary conditions
!           with A1=(1,0,0), A2=(0,1,0) in target frame
!           spheres may be of differing, possibly anisotropic,
!           compositions.
! SHPAR(1) = DIAMX = max extent in target frame x direction/d
! SHPAR(2) = PYAEFF = (periodicity in target y direction)/AEFF
! SHPAR(3) = PZAEFF = (periodicity in target z direction)/AEFF

      ELSEIF(CSHAPE=='SPHRN_PBC')THEN
         CALL TARNAS(A1,A2,SHPAR(1),0._WP,DX,X0,BETADF,PHIDF,THETADF,CFLSHP, &
                    CDESCR,IDVSHP,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
! diagnostic
!         write(0,*)'target ckpt 60 - sphrn_pbc'
!         write(0,*)'                x0(1-3)=',x0
!
         IANISO=0.
         PI=4._WP*ATAN(1._WP)
         PYD=SHPAR(2)*(3._WP*REAL(NAT,KIND=WP)/(4._WP*PI))**(1._WP/3._WP)/DX(2)
         PZD=SHPAR(3)*(3._WP*REAL(NAT,KIND=WP)/(4._WP*PI))**(1._WP/3._WP)/DX(3)

!-----------------------------------------------------------------------
! SPHROID_2 -> two touching spheroids

      ELSEIF(CSHAPE=='SPHROID_2')THEN
         CALL TAR2SP(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),SHPAR(4),SHPAR(5), &
                     SHPAR(6),DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)

!-----------------------------------------------------------------------
! SPH_ANI_N -> Multisphere target of general anisotropic composition.
!           Different spheres can have different composition, and
!           constituent material can have arbitrary orientation

      ELSEIF(CSHAPE=='SPH_ANI_N')THEN

         CALL TARNAS(A1,A2,SHPAR(1),SHPAR(2),DX,X0,BETADF,PHIDF,THETADF, &
                     CFLSHP,CDESCR,IDVSHP,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
         IANISO=1
         DO J=1,NAT
            IF(BETADF(J)/=0.)IANISO=2
            IF(PHIDF(J)/=0.)IANISO=2
            IF(THETADF(J)/=0.)IANISO=2
         ENDDO

!-----------------------------------------------------------------------
! 2010.02.06 (BTD) Eliminate this option -- this target geometry can already
!                  be created with option RCTGLBLK3

! TARSLBLIN -> rectangular slab with rectangular "line" on top
!
!      ELSEIF(CSHAPE=='TARSLBLIN')THEN

!         CALL TARSLBLIN(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),SHPAR(4),SHPAR(5), &
!                        SHPAR(6),DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)

!-----------------------------------------------------------------------
! TETRAHDRN -> homogeneous, isotropic tetrahedron

      ELSEIF(CSHAPE=='TETRAHDRN')THEN
         CALL TARTET(A1,A2,SHPAR(1),DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
         IANISO=0

!-----------------------------------------------------------------------
! TRILYRPBC -> trilayer rectangular structure, periodic in y and z directions
!              SHPAR(1) = x-thickness of upper layer [material 1]
!              SHPAR(2) = y-width/d of upper layer
!              SHPAR(3) = z-width/d of upper layer
!              SHPAR(4) = x-thickness/d of middle layer [material 2]
!              SHPAR(5) = y-width/d of middle layer
!              SHPAR(6) = z-width/d of middle layer
!              SHPAR(7) = x-width/d of lower layer
!              SHPAR(8) = y-width/d of lower layer
!              SHPAR(9) = z-width/d of lower layer
!              SHPAR(10)= period/d in y direction
!              SHPAR(11)= period/d in z direction

      ELSEIF(CSHAPE=='TRILYRPBC')THEN
         CALL TARRCTBLK3(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),SHPAR(4),SHPAR(5), &
                         SHPAR(6),SHPAR(7),SHPAR(8),SHPAR(9),DX,X0,CDESCR,   &
                         IOSHP,MXNAT,NAT,IXYZ,ICOMP)
         IANISO=0
         PYD=SHPAR(10)
         PZD=SHPAR(11)

!-----------------------------------------------------------------------
! TRNGLPRSM -> triangular prism:

      ELSEIF(CSHAPE=='TRNGLPRSM')THEN
         CALL TARPRSM(A1,A2,SHPAR(1),SHPAR(2),SHPAR(3),SHPAR(4),DX,X0,CDESCR, &
                      IOSHP,MXNAT,NAT,IXYZ,ICOMP)
         IANISO=0

!-----------------------------------------------------------------------
! UNIAXICYL -> uniaxial cylinder
! Cylinder axis is in x direction in Target Frame.
! Assume c axis to be parallel to cylinder axis
! Assume component 1 to be dielectric const. for E par. c axis
! Assume component 2 to be dielectric const. for E perp. c axis

      ELSEIF(CSHAPE=='UNIAXICYL')THEN
         CALL TARCYL(A1,A2,SHPAR(1),SHPAR(2),1._WP,DX,X0,CDESCR,IOSHP,MXNAT, &
                     NAT,IXYZ,ICOMP)

! Set composition

         DO J=1,NAT
            ICOMP(J,1)=1
            ICOMP(J,2)=2
            ICOMP(J,3)=2
         ENDDO
         IANISO=1

!-----------------------------------------------------------------------

! add any other specific shape routines here:

!-----------------------------------------------------------------------
! Reach here only in event of error:

      ELSE
         CALL ERRMSG('FATAL','TARGET',                           &
                     ' Unknown shape requested: check ddscat.par')
      ENDIF

! Check to see that arrays are large enough

      IF(NAT>MXNAT)CALL ERRMSG('FATAL','TARGET',     &
        ' MXNAT < NAT; must increase MXNAT in DDSCAT')

! Write out target description

      WRITE(IDVOUT,FMT=9010)CDESCR,NAT

! Before returning, compute:

      NAT3=3*NAT
      RETURN
9010  FORMAT(' >TARGET:',A,/,6X,I9,' = NAT0 = number of dipoles in target')
    END SUBROUTINE TARGET
