    SUBROUTINE REAPAR(CBINFLAG,AEFFA,BETAMI,BETAMX,CALPHA,CFLEPS,CFLPAR,     &
                      CFLSHP,CMDSOL,CMDFFT,CMDFRM,CSHAPE,CMDTRQ,CXE01,CXE02, &
                      DEGRAD,ETASCA,GAMMA,IDVERR,IDVOUT,INIT,IOPAR,IORTH,    &
                      JPBC,IWRKSC,IWRPOL,MXNX,MXNY,MXNZ,MXBETA,MXCOMP,MXPHI, &
                      MXRAD,MXSCA,MXTHET,MXWAV,NBETA,NCOMP,NPHI,NRAD,NSCA,   &
                      NTHETA,NWAV,IWAV0,IRAD0,IORI0,NSMELTS,PHIMIN,PHIMAX,   &
                      PHIN,THETAN,THETMI,THETMX,TOL,WAVEA,SHPAR,SMIND1,      &
                      SMIND2,DX,CENTER,ELENERGY) !ARGUMENT 'CENTER' and 'ELENERGY' ADDED BY NWB 3/13/12 and 7/12/12
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE

! Arguments:

      CHARACTER :: CALPHA*(*),CBINFLAG*(*),CFLPAR*(*),CFLSHP*80,CMDSOL*(*), &
                   CMDFFT*(*),CMDFRM*(*),CMDTRQ*(*),CSHAPE*(*)
      INTEGER :: IDVERR,IDVOUT,INIT,IOPAR,IORI0,IORTH,JPBC,IRAD0,IWAV0,      &
                 IWRKSC,IWRPOL,MXBETA,MXCOMP,MXNX,MXNY,MXNZ,MXPHI,MXRAD,     &
                 MXSCA,MXTHET,MXWAV,NBETA,NCOMP,NPHI,NRAD,NSCA,NSMELTS,NWAV, &
                 NTHETA
      INTEGER ::    &
         SMIND1(9), &
         SMIND2(9)
      CHARACTER(60) :: &
         CFLEPS(MXCOMP)
      REAL(WP) :: BETAMI,BETAMX,DEGRAD,ETASCA,GAMMA,PHIMIN,PHIMAX, &
                  THETMI,THETMX,TOL,ELENERGY !ELENERGY ADDED BY NWB 7/12/12
      REAL(WP) ::       &
         AEFFA(MXRAD),  &
         DX(3),         &
         PHIN(MXSCA),   &
         SHPAR(12),     &
         THETAN(MXSCA), &
         WAVEA(MXWAV),  &
         CENTER(3) !CENTER added by NWB 3/13/12
      COMPLEX(WP) :: &
         CXE01(3),   &
         CXE02(3)

! Local variables:

      CHARACTER :: CDIVID*3,CLINE*70,CMSGNM*70,CWHERE*70
      INTEGER :: J,JPLANE,NPLANES,NSCA0
      REAL(WP) :: AEFEND,AEFINI,DELTA,DTHETA,E1,PHI1,THETA1,THETA2, &
         WAVEND,WAVINI

! External routines:

      EXTERNAL DIVIDE,ERRMSG,WRIMSG

!***********************************************************************
! Subroutine REAPAR handles the reading of input parameters from the
! "ddscat.par" file, as well as elementary processing with those input
! parameters to generate arrays.

! Given:
!       CFLPAR = name of file (normally 'ddscat.par')

! Returns:
!       CMDTRQ   = 'NOTORQ' or 'DOTORQ'
!       CMDSOL   = 'PBCGST' or 'PBCGS2' or 'PETRKP'
!       CMDFFT   = 'GPFAFT' or 'FFTW21' or 'FFTMKL'
!       CALPHA   = 'LATTDR' or 'GKDLDR'
!       CBINFLAG = 'ALLBIN' or 'ORIBIN' or 'NOTBIN'
!       CSHAPE   = 'RCTGLPRSM' or 'ELLIPSOID' or ... other valid shape
!       JPBC     = 0 :isolated finite targets
!                = 1 :target periodic in y direction, finite in z
!                = 2 :target finite in y direction, periodic in z
!                = 3 :target periodic in both y and z directions
!       NCOMP    = number of compositions
!       CFLEPS(1-NCOMP) = filenames for dielectric functions
!       INIT     = 0 or 1 or 2
!       TOL      = error tolerance
!       GAMMA    = parameter controlling integration limit for PBC
!       ETASCA   = parameter controlling accuracy of calculation
!                  of <cos(theta)>
!                  1 is OK, 0.5 gives high accuracy
!       NWAV     = number of wavelengths
!       WAVEA(1-NWAV) = wavelengths (physical units)
!       NRAD     = number of radii
!       AEFFA(1-NRAD) = target effective radii (physical units)
!       CXE01(1-3)= (complex) polarization state 1
!       CXE02(1-3)= complex polarization state 2 (orthogonal to 1)
!       IORTH    = 1 or 2 (number of incident polarization states to use
!       IWRKSC   = 0 or 1 (not write/write ".sca" files)
!       IWRPOL   = 0 or 1 (not write/write ".pol" files)
!       NBETA    = number of beta values for target rotation
!       BETAMI   = minimum beta value (rad) [input from file in deg]
!       BETAMX   = maximum beta value (rad) [input from file in deg]
!       NTHETA   = number of theta values for target rotation
!       THETAMI  = minimum theta value (rad) [input from file in deg]
!       THETAMX  = maximum theta value (rad) [input from file in deg]
!       NPHI     = number of PHI values for target rotation
!       PHIMIN   = minimum PHI value (rad) [input from file in deg]
!       PHIMAX   = maximum PHI value (rad) [input from file in deg]
!       NSMELTS  = number of elements of S matrix to calculate
!                  (less than or equal to 9)
!       SMINDI(J) = indices of S matrix elements to be calculated
!                   (e.g., 11 , 21 , 31 , 41 , 12 , 13)
!       CMDFRM   = 'LFRAME' or 'TFRAME' (Lab frame or Target frame)
!       NSCA     = number of scattering directions

! If JPBC=0:
!       THETAN(1-NSCA) = theta for scattering directions (rad)
!       PHIN(1-NSCA) = phi for scattering directions (rad)

! If JPBC=1 or 2:
!       THETAN(1-NSCA) = order_y or order_z for scattering directions
!       PHIN(1-NSCA)   = psi (radians) for scattering directions,
!                        where psi = rotation around axis of
!                        target periodicity (y_TF if JPBC=1, z_TF if JPB

! If JPBC=3: periodic in both y and z directions
!       THETAN(1-NSCA) = order_y for diffraction
!       PHIN(1-NSCA)   = order_z for diffraction

! Original version created by P.J.Flatau, Colorado State Univ.
! Modified by B.T.Draine, Princeton Univ. Obs.
! History:
! 90.11.09 (BTD): Modified so as to remove sensitivity to errors
!                 in sign of entered value of DTHETA.
! 90.11.30 (BTD): Added test to compare CMDFFT and MXMEM
! 90.12.14 (BTD): Changed option MKRECT to RCTGLPRSM
! 91.01.03 (BTD): Added MXBETA,MXPHI,MXTHET to argument list
!                 and added checking of NBETA,NPHI,NTHETA
! 91.05.08 (BTD): Added CALPHA to allow choice of prescription
!                 for polarizabilities
! 91.05.08 (BTD): Some rewriting done.
! 91.05.14 (BTD): Added separate options LDRXYZ and LDRAVG
! 91.05.15 (BTD): Added new option LDR000
! 91.05.23 (BTD): Added variable IWRKSC to argument list and
!                 added code to read IWRKSC
! 91.05.30 (BTD): Changed LDRXYZ->LDR100
!                 Added option LDR111
! 91.11.13 (BTD): Now three options: DRAI88, GOBR88, LATTDR
!                 (subroutine ALPHA now allows for directional
!                 dependence for LATTDR option)
! 92.05.14 (BTD): Allow for experimental option LDRISO
! 92.09.21 (BTD): Allow for option ANIELLIPS (ellipsoid of anisotropic
!                 material); change option GRPHIT -> UNIAXICYL (cylinder
!                 of uniaxial material)
! 93.01.07 (BTD): Allow for options ELLIPSO_2 and ANI_ELL_2 (two touching
!                 ellipsoids, each either isotropic or anisotropic)
! 93.03.11 (BTD): Deleted all code associated with variable CMACHN
!                 which previously was used to identify machine/OS
!                 This information was not used or needed.
! 93.12.15 (BTD): Added code to check for compatibility of NCOMP
!                 when CSHAPE=ELLIPSO_2,ELLIPSO_3,ANI_ELL_2,or ANI_ELL_3
! 94.01.27 (BTD): Replaced SHPAR1,SHPAR2,SHPAR3 by SHPAR(1-6)
!                 Added CONELLIPS (two concentric ellipsoids) as target
!                 option.
! 94.06.20 (PJF): Added code to handle NEWTMP
! 94.12.02 (BTD): Added options
!                    LB1991 (Dungey-Bohren 1991 polarizability)
!                    LR1994 (Lumme-Rahola 1994 polarizability)
! 95.05.26 (BTD): Corrected construction of orthogonal polarization
!                 state for case where E01 is complex (e.g.,
!                 circular polarization)
! 95.06.14 (BTD): Added code to check that user is not requesting
!                 f_ml when illuminating target with other then
!                 linearly polarized light
! 95.06.19 (BTD): Added argument CMDSOL to allow selection of
!                 different CCG methods
!                 Added argument CMDTRQ to allow torque calculations
!                 to be skipped.
! 96.01.05 (BTD): Added shape options MLTBLOCKS and DW1996TAR
! 96.01.25 (BTD): Added shape option SPHROID_2
! 96.01.25 (BTD): Disable code added 95.06.14 so that user can
!                 calculate f_ml even when illuminating target with
!                 elliptically polarized light.
! 96.10.18 (BTD): Changed option NEWTMP to GPFAFT (Generalized Prime
!                 Factor Algorithm for Fourier Transform).
! 96.11.06 (BTD): Added provision to
!                 (1) produce fatal error if CXE01(1).NE.(0.,0.)
!                 (2) normalize CXE01 if necessary
! 96.11.14 (PJF): cbinflag for binary file
! 96.11.15 (PJF): netCDF netflag
! 96.11.21 (BTD): modified handling of CBINFLAG and CNETFLAG to
!                 be consistent with other option flags
! 97.11.01 (BTD): add DX to argument list
!                 read DX from ddscat.par
!                 renormalize DX (input DX not required to be
!                 normalized with DX(1)*DX(2)*DX(3)=1.)
! 98.12.06 (BTD): modified to read NSMELTS,SMELTS from ddscat.par
!                 and pass NSMELTS,SMIND1,SMIND2 as arguments to calling
!                 program DDSCAT
! 98.12.21 (BTD): changed dimension of CFLPAR from CHARACTER*40 to
!                 CHARACTER*60 to allow longer file names
!                 (also changed in DDSCAT.f and dielec.f)
! 98.12.29 (BTD): modified to support new option LYRD_SLAB
! 99.01.26 (BTD): modified to allow first IWAV,IRAD,IORI to be
!                 specified in ddscat.par (useful for restarting
!                 jobs which did not complete)
! 99.06.29 (BTD): corrected input of SHPAR for option SPHROID_2
!                 previously had input only SHPAR(1-5),
!                 now input SHPAR(1-6)
! 00.06.12 (BTD): modified to support SPHERES_N option
! 00.06.23 (BTD): modified to support FFTWFJ option
!                 remove MXMEM from argument list
!                 since no longer supporting TMPRTN FFT option which
!                 required additional memory
! 00.07.05 (BTD): modified to print more information about CMDFFT
!                 option
! 02.02.12 (BTD): modified to support TRNGLPRSM option
! 03.01.29 (BTD): modified to force DX(1-3)=1.
! 03.05.21 (BTD): modified to support SPHARM option
! 03.05.22 (BTD): bug fix re SPHARM option
! 03.07.13 (BTD): changed FFTWFJ to FFTW21 in anticipation of future
!                 option FFTW3x (FFTW 2.1.x and FFTW 3.0.1 have
!                 different interfaces, so they must be distinguished).
! 03.10.23 (BTD): removed ICTHM and IPHIM from argument list,
!                 removed code to read these values from ddscat.par
! 03.10.24 (BTD): Add parameter ETASCA to allow user to control
!                 angular resolution used for scattering calculations.
!                 Add checks against invalid input
! 04.02.26 (BTD): Modify to accept polarizability option 'GKDLDR'
!                 for Gutkowicz-Krusin & Draine (2004) lattice
!                 dispersion relation.
! 04.04.01 (BTD): Modify to support option NPSPBC (multi-sphere cluster
!                 with periodic boundary conditions)
! 04.04.05 (BTD): Remove SPHARM option: GAUSS_SPH is sufficient
! 04.04.29 (BTD): Added ANIRCTNGL option for homogeneous anisotropic
!                 rectangular target
! 04.05.23 (BTD): Added CWHERE string to provide error information
!                 in case of malformed ddscat.par file.
! 04.09.14 (BTD): Permit new target option ANIFRMFIL, to read target data
!                 from file, with target allowed to be anistropic
!                 with arbitrary orientations of local "Dielectric Frame
!                 (in which local dielectric tensor is diagonalized)
!                 relative to Target Frame.
! 04.10.12 (BTD): Modified to support new target option SPH_ANI_N, for
!                 a target composed of N possibly anisotropic spheres.
! 04.10.14 (BTD): Modified to support option of reading wavelengths
!                 from file 'wave.tab' if CDIVID='TAB' in ddscat.par
! 05.03.19 (BTD): Added 2 lines to set CLFEPS(1) in case using options
!                 H2OLIQ or H2OICE (previously were left undefined,
!                 which cause problems with output statement)
! 05.08.04 (BTD): Added requirement to read string CMDFRM where
!                 CMDFRM = 'LFRAME' to specify scattering directions in
!                          Lab Frame (rel. to xlab,ylab,zlab,
!                          where xlab = direction of incident beam
!                 CMDFRM = 'TFRAME' to specify scattering directions in
!                          Target Frame (rel. to a1,a2,a3)
! 05.10.18 (BTD): Added check to ensure that NBETA*NTHETA*NPHI < 1001
!                 if IWRKSC=1 (since current version of NAMER is limited
!                 to maximum of 1000 directions)
! 06.04.10 (BTD): Added IWRPOL to argument list, and modified REAPAR
!                 to read IWRPOL from ddscat.par in order to allow
!                 user to elect to write or to not write out polarizatio
!                 array for subsequent processing.
! 06.09.13 (BTD): Added handling of option CYLNDRCAP
! 06.09.15 (BTD): Modified to support new version of target routine
!                 TARHEX, that now requires SHPAR(3) to select prism
!                 orientation.
! 06.09.15 (BTD): Modified to support new target option HEXGONPBC
! 06.10.04 (BTD): Extension to support new scattering directions
!                 required for far-field scattering when PBC is used.
! 06.10.11 (BTD): Corrected errors in handling of JPBC=2,3,
! 06.11.29 (BTD): Added code to catch error in input CMDFRM
! 06.12.08 (BTD): Added support for target option DSKRCTPBC
! 07.01.20 (BTD): Added support for target option DSKRCTNGL
! 07.02.23 (BTD): Added support for target option LYRSLBPBC
! 07.02.25 (BTD): Corrected bug (typo)
! 07.08.05 (BTD): Version 7.0.3
!                 Add code to read MXNX,MXNY,MXNZ from ddscat.par
! 07.08.12 (BTD): Allow CMDSOL='PBCGS2'
! 07.10.26 (BTD): Add support for target option DSKBLYPBC
!                 reorder input shape paramers SHPAR for DSKRCTPBC
!                 changed SHPAR(6) -> SHPAR(10)
! 07.10.27 (BTD): Further improve of reapar
!                 Eliminate CDIEL -- option 'TABLES' is now the only
!                 one used, so no need to input and carry around CDIEL
! 08.02.01 (BTD): Changed SHPAR(10) -> SHPAR(12) to allow up to 12
!                 shape parameters to be input
!                 Added support for target option BLSLINPBC
!                 Added support for target option RCTGLBLK3
!                 Added support for target option TRILYRPBC
! 08.03.11 (BTD): corrected typo: BLSLINPBC -> BISLINPBC
! 08.03.11 (BTD): ver7.0.5
!                 Added ALPHA to argument list.
!                 Modified to read parameter ALPHAA from ddscat.par
! 08.03.12 (BTD): Modified to read NPLANES from ddscat.par
!                 so that we no longer rely on END condition when reading
!                 scattering directions or scattering planes from ddscat.par
! 08.03.12 (BTD): Modified to use new form of ddscat.par
!                 eliminate reading of INIT (since we always set INIT=0)
!                 add new descriptive lines separating ddscat.par into sections
! 08.04.19 (BTD): Changed notation: ALPHA -> GAMMA
!                 reordered argument list
! 08.06.07 (BTD): Added support for FFT option FFTMKL
! 08.07.22 (BTD): Added sanity check to catch users who call for SPHERES_N but
!                 specify multiple dielectric functions
! 08.08.29 (BTD): Removed CNETFLAG from argument list
!                 Removed code to read CNETFLAG from ddscat.par
! 08.09.12 (BTD): Corrected typo DSKRCTGNL -> DSKRCTNGL
! 08.09.16 (BTD): Corrected bug reported by Wang Lin, Nankai Univ.:
!                 typo SPHERN_PBC -> SPHRN_PBC caused SPHERN_PBC to be 
!                 assigned JPBC=0 regardless of values of SHPAR(2) or SHPAR(3).
! 08.09.23 (BTD): Corrected typo DSKBLYPRC -> DSKLBYLRPBC line 613
!                 [typo reported by Georges Levi, Univ-Paris-Diderot]
! 09.08.27 (BTD): Added code to test for insane value of TOL
!                 (which could result if other than NCOMP diel. fns.
!                  are listed in ddscat.par)
! 09.09.10 (BTD): ver7.1.0
!                 Added support for new target options
!                   ANIFILPBC
!                   FRMFILPBC
!                 Added ability to recognize lowercase target options
! 09.12.11 (BTD): Corrected coding error which prevented use of
!                 target options DW1996TAR and MLTBLOCKS
! 10.01.02 (BTD): Add code to write out SHPAR info to running output
! 10.02.05 (BTD): Corrected error that caused DW1996TAR option to fail
! end history
! Copyright (C) 1993,1994,1995,1996,1997,1998,1999,2000,2002,2003,2004
!               2005,2006,2007,2008,2009,2010 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.

!***********************************************************************

      OPEN(UNIT=IOPAR,FILE=CFLPAR,STATUS='OLD')
      CWHERE='error reading line 1 of ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)CLINE
      CALL WRIMSG('REAPAR',CLINE)
      CWHERE='error reading line 2 of ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)CLINE
      CALL WRIMSG(' ',CLINE)

!=======================================================================
!    Specify whether torques are to be calculated:

!    CMDTRQ*6 = 'DOTORQ' -  calculate torque on grain
!             = 'NOTORQ' -  do not calculate torque on grain

      CWHERE='error reading CMDTRQ from ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)CMDTRQ
      IF(CMDTRQ=='dotorq')CMDTRQ='DOTORQ'
      IF(CMDTRQ=='notorq')CMDTRQ='NOTORQ'
      IF(CMDTRQ=='DOTORQ')THEN
         WRITE(CMSGNM,FMT='(A,A)')CMDTRQ,' - compute torques '
         CALL WRIMSG('REAPAR',CMSGNM)
      ELSEIF(CMDTRQ=='NOTORQ')THEN
         WRITE(CMSGNM,FMT='(A,A)')CMDTRQ,' - do not compute torques '
         CALL WRIMSG('REAPAR',CMSGNM)
      ELSE
! diagnostic
         write(0,*)'reapar ckpt 0.1 : CMDTRQ=',CMDTRQ
!
         CALL ERRMSG('FATAL','REAPAR',' wrong definition of CMDTRQ')
      ENDIF

!=======================================================================
!    Define method used for iterative solution of complex linear equations

!    CMDSOL*6 = 'PETRKP' -  Petravic & Kuo-Petravic method
!             = 'PBCGST' -  PIM BiConjugate Gradient with Stabilization
!             = 'PBCGS2' -  M.A Botcheve implementation of BiCGstab 
!                           enhanced to improve convergence properties
!                           with finite precision arithmetic

      CWHERE='error reading CMDSOL from ddscat.par'
      READ(IOPAR,FMT=*,ERR=99) CMDSOL
      IF(CMDSOL=='petrkp')CMDSOL='PETRKP'
      IF(CMDSOL=='pbcgst')CMDSOL='PBCGST'
      IF(CMDSOL=='pbcgs2')CMDSOL='PBCGS2'
      IF(CMDSOL=='PETRKP'.OR. &
         CMDSOL=='PBCGST'.OR. &
         CMDSOL=='PBCGS2')THEN
         WRITE(CMSGNM,FMT='(A,A)')CMDSOL,' - CCG Method  '
         CALL WRIMSG('REAPAR',CMSGNM)
      ELSE
         CALL ERRMSG('FATAL','REAPAR',' wrong definition of CMDSOL')
      ENDIF

!=======================================================================
!    Define FFT method:

!    CMDFFT*6 = 'GPFAFT' -  GPFA code of Temperton
!             = 'FFTW21' -  FFTW 2.1.x code of Frigo & Johnson
!             = 'FFTMKL' -  Use DFTI from Intel MKL

      CWHERE='error reading CMDFFT from ddscat.par'
      READ(IOPAR,FMT=*,ERR=99) CMDFFT
      IF(CMDFFT=='FFTW21')THEN
         WRITE(CMSGNM,FMT='(A,A,A)')CMDFFT,' - using FFTW 2.1.x ', &
                                    'package from Frigo & Johnson'
         CALL WRIMSG('REAPAR',CMSGNM)
      ELSEIF(CMDFFT=='GPFAFT')THEN
         WRITE(CMSGNM,FMT='(A,A,A)')CMDFFT,' - using GPFA package ', &
                                    'from Clive Temperton'
         CALL WRIMSG('REAPAR',CMSGNM)
      ELSEIF(CMDFFT=='FFTMKL')THEN
         WRITE(CMSGNM,FMT='(A,A,A)')CMDFFT,' - using DFTI from Intel ', &
                                    'Math Kernel Library (MKL)'
         CALL WRIMSG('REAPAR',CMSGNM)
      ELSE
         WRITE(CMSGNM,FMT='(A,A)')'DDSCAT 7.0 only supports FFT options', &
                                  ' FFTW21, GPFAFT, and FFTMKL'
         CALL WRIMSG('REAPAR',CMSGNM)
         CALL ERRMSG('FATAL','REAPAR',' wrong definition of CMDFFT')
      ENDIF

!=======================================================================

!    Define prescription for computing polarizabilities:

!    CALPHA*6 = 'LATTDR' - Lattice Dispersion Relation of
!                          Draine & Goodman (1993)
!               'GKDLDR' - Lattice Dispersion Relation of
!                          Gutkowicz-Krusin & Draine (2004)

      CWHERE='error reading CALPHA from ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)CALPHA
      IF(CALPHA=='LATTDR')THEN
         CALL WRIMSG('REAPAR','LATTDR - Draine & Goodman (1993) LDR for alpha')
      ELSEIF(CALPHA=='GKDLDR')THEN
         CALL WRIMSG('REAPAR', &
                    'GKDLDR - Gutkowicz-Krusin & Draine (2004) LDR for alpha')
      ELSE
         CALL ERRMSG('FATAL','REAPAR',' wrong definition of CALPHA')
      ENDIF

!=======================================================================
! binary file flag

      CWHERE='error reading CBINFLAG from ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)CBINFLAG
      IF(CBINFLAG=='allbin')CBINFLAG='ALLBIN'
      IF(CBINFLAG=='oribin')CBINFLAG='ORIBIN'
      IF(CBINFLAG=='notbin')CBINFLAG='NOTBIN'
      IF(CBINFLAG=='ALLBIN'.OR.CBINFLAG=='ORIBIN'.OR.CBINFLAG=='NOTBIN')THEN
         WRITE(CMSGNM,FMT='(A,A)')CBINFLAG,' - Unformatted binary dump option'
         CALL WRIMSG('REAPAR',CMSGNM)
      ELSE
         CALL ERRMSG('FATAL','REAPAR',' Wrong definition of CBINFLAG')
      ENDIF

!=======================================================================
! Read upper bound on target extent
! skip line:
      READ(IOPAR,FMT=*,ERR=99)CLINE
      CALL WRIMSG('REAPAR',CLINE)
      CWHERE='error reading target size bounds MXNX,MXNY,MXNZ from ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)MXNX,MXNY,MXNZ
      WRITE(CMSGNM,FMT='(A,3I6,A)')'allow MXNX,MXNY,MXNZ=',MXNX,MXNY,MXNZ, &
                                   ' for target generation'
      CALL WRIMSG('REAPAR',CMSGNM)
!=======================================================================

! e-BEAM CENTROID READER
! Added by NWB 3/13/12

      READ(IOPAR,FMT=*,ERR=99)CENTER

!======================================================================= 

! ELECTRON RELATIVISTIC ENERGY READER
! Added by NWB 7/12/12

      READ(IOPAR,FMT=*,ERR=99)ELENERGY

!=======================================================================
!    Define shape:

! skip line:

      READ(IOPAR,FMT=*,ERR=99)CLINE
      CALL WRIMSG('REAPAR',CLINE)

!    CSHAPE*9 = 'FROM_FILE'  shape and composition data will later
!                            be read from file CFLSHP
!                            for dielectric tensors that are diagonal
!                            in Target Frame
!             = 'ANIFRMFIL'  read shape and composition data from file,
!                            for general anisotropic dielectric tensors
!                            with orientation angles THETADF,PHIDF,BETAD
!                            relative to Target Frame.
!             = 'FRMFILPBC'  shape and composition data for TUC will later
!                            be read from file CFLSHP
!                            for dielectric tensors that are diagonal
!                            in Target Frame
!                            PYD and PZD are input via ddscat.par
!             = 'ANIFILPBC'  shape and composition data for TUC will later
!                            be read from file CFLSHP
!                            for general anisotropic dielectricc tensors
!                            with orientation angles THETADF,PHIDF,BETAD
!                            relative to Target Frame.
!                            PYD and PZD are input via ddscat.par
!             = 'ANIELLIPS'  ellipsoid of anisotropic material
!             = 'ANIRCTNGL'  homogeneous anisotropic rectangular target
!             = 'ANI_ELL_2'  two touching anisotropic ellipsoids of
!                            materials 1-6
!             = 'ANI_ELL_3'  three touching anisotropic ellipsoids of
!                            materials 1-9
!             = 'BISLINPBC'  bilayer slab with periodic grid of lines
!                            parallel to z on top, with y-period/d=PYD
!                            [or, if PYD=0, a single line]
!             = 'CONELLIPS'  two concentric ellipsoids of materials 1,2
!             = 'CYLINDER1'  homogeneous finite cylinder
!             = 'CYLNDRCAP'  homogeneous cylinder with hemispherical endcaps
!             = 'CYLNDRPBC'  1-d or 2-d array of finite cylinders
!             = 'DSKBLYPBC'  1-d or 2-d array of disk on bilayer rect. slab
!             = 'DSKRCTNGL'  single disk on rectangular slab
!             = 'DSKRCTPBC'  1-d or 2-d array of disk on rectangular slab
!             = 'DW1996TAR'  13-cube target used by Draine & Weingartner 1996
!             = 'ELLIPSOID'  ellipsoid (homogeneous and isotropic)
!             = 'ELLIPSPBC'  1-d or 2-d array of ellipsoids
!             = 'ELLIPSO_2'  two touching isotropic ellipsoids of 
!                            materials 1 and 2
!             = 'ELLIPSO_3'  three touching isotropic ellipsoids of
!                            materials 1,2,3
!             = 'GAUSS_SPH'  gaussian sphere target
!             = 'HEXGONPBC'  1-d or 2-d array of hexagonal prisms
!             = 'HEX_PRISM'  homogeneous hexagonal prism
!             = 'LYRD_SLAB'  layered slab target, with up to 4 separate
!                            material layers
!             = 'LYRSLBPBC'  1-d or 2-d array of layered rect. slab targets,
!                            with up to 4 material layers
!             = 'MLTBLOCKS'  collection of cubic blocks defined by
!                            data in file 'blocks.par'
!             = 'RCTGLBLK3'  isolated target: 3 rectangular blocks with
!                            centers on x-axis
!             = 'RCTGLPRSM'  homogeneous rectangular prism
!             = 'RCTGL_PBC'  1-d or 2-d array of rectangular targets
!             = 'SLAB_HOLE'  rectangular block with cylindrical hole
!             = 'SLBHOLPBC'  1-d or 2-d array of rectangular blocks
!                            with cylindrical hole
!             = 'SPHERES_N'  multisphere target = union of N spheres
!             = 'SPHRN_PBC'  1-d or 2-d array of multisphere target
!             = 'SPHROID_2'  two touching spheroids with symmetry axes at
!                            specified angle!
!             = 'SPH_ANI_N'  multisphere target, with spheres that
!                            can have different, anisotropic, composition
!             = 'TETRAHDRN'  regular tetrahedron
!             = 'TRILYRPBC'  periodic target: 3 layer rectangular structure
!             = 'TRNGLPRSM'  triangular prism (homogeneous and isotropic)
!             = 'UNIAXICYL'  cylinder of unixaxial material

      CWHERE='error reading CSHAPE from ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)CSHAPE
      IF(CSHAPE=='aniellips')CSHAPE='ANIELLIPS'
      IF(CSHAPE=='ani_ell_2')CSHAPE='ANI_ELL_2'
      IF(CSHAPE=='anifilpbc')CSHAPE='ANIFILPBC'
      IF(CSHAPE=='anifrmfil')CSHAPE='ANIFRMFIL'
      IF(CSHAPE=='anirctngl')CSHAPE='ANIRCTNGL'
      IF(CSHAPE=='bislinpbc')CSHAPE='BISLINPBC'
      IF(CSHAPE=='conellips')CSHAPE='CONELLIPS'
      IF(CSHAPE=='cylinder1')CSHAPE='CYLINDER1'
      IF(CSHAPE=='cylndrcap')CSHAPE='CYLNDRCAP'
      IF(CSHAPE=='cylndrpbc')CSHAPE='CYLNDRPBC'
      IF(CSHAPE=='dskblypbc')CSHAPE='DSKBLYPBC'
      IF(CSHAPE=='dskrctpbc')CSHAPE='DSKRCTPBC'
      IF(CSHAPE=='dw1996tar')CSHAPE='DW1996TAR'
      IF(CSHAPE=='ellipsoid')CSHAPE='ELLIPSOID'
      IF(CSHAPE=='ellipso_2')CSHAPE='ELLIPSO_2'
      IF(CSHAPE=='ellipso_3')CSHAPE='ELLIPSO_3'
      IF(CSHAPE=='ellipspbc')CSHAPE='ELLIPSPBC'
      IF(CSHAPE=='frmfilpbc')CSHAPE='FRMFILPBC'
      IF(CSHAPE=='from_file')CSHAPE='FROM_FILE'
      IF(CSHAPE=='gauss_sph')CSHAPE='GAUSS_SPH'
      IF(CSHAPE=='hexgonpbc')CSHAPE='HEXGONPBC'
      IF(CSHAPE=='hex_prism')CSHAPE='HEX_PRISM'
      IF(CSHAPE=='layrdslab')CSHAPE='LAYRDSLAB'
      IF(CSHAPE=='lyrslbpbc')CSHAPE='LYRSLBPBC'
      IF(CSHAPE=='mltblocks')CSHAPE='MLTBLOCKS'
      IF(CSHAPE=='rctglblk3')CSHAPE='RCTGLBLK3'
      IF(CSHAPE=='rctglprsm')CSHAPE='RCTGLPRSM'
      IF(CSHAPE=='rctgl_pbc')CSHAPE='RCTGL_PBC'
      IF(CSHAPE=='rctg_rctg')CSHAPE='RCTG_RCTG'
      IF(CSHAPE=='recrecpbc')CSHAPE='RECRECPBC'
      IF(CSHAPE=='slab_hole')CSHAPE='SLAB_HOLE'
      IF(CSHAPE=='slbholpbc')CSHAPE='SLBHOLPBC'
      IF(CSHAPE=='spheres_n')CSHAPE='SPHERES_N'
      IF(CSHAPE=='sphrn_pbc')CSHAPE='SPHRN_PBC'
      IF(CSHAPE=='sphroid_2')CSHAPE='SPHROID_2'
      IF(CSHAPE=='sph_ani_n')CSHAPE='SPH_ANI_N'
      IF(CSHAPE=='tetrahdrn')CSHAPE='TETRAHDRN'
      IF(CSHAPE=='trilyrpbc')CSHAPE='TRILYRPBC'
      IF(CSHAPE=='trnglprsm')CSHAPE='TRNGLPRSM'
      IF(CSHAPE=='uniaxicyl')CSHAPE='UNIAXICYL'
 
    IF(CSHAPE=='ANIELLIPS'.OR. &
         CSHAPE=='ANI_ELL_2'.OR. &
         CSHAPE=='ANI_ELL_3'.OR. &
         CSHAPE=='ANIFRMFIL'.OR. &
         CSHAPE=='ANIRCTNGL'.OR. &
         CSHAPE=='BISLINPBC'.OR. &
         CSHAPE=='CONELLIPS'.OR. &
         CSHAPE=='CYLINDER1'.OR. &
         CSHAPE=='CYLNDRCAP'.OR. &
         CSHAPE=='CYLNDRPBC'.OR. &
         CSHAPE=='DSKBLYPBC'.OR. &
         CSHAPE=='DSKRCTNGL'.OR. &
         CSHAPE=='DSKRCTPBC'.OR. &
         CSHAPE=='DW1996TAR'.OR. &
         CSHAPE=='ELLIPSOID'.OR. &
         CSHAPE=='ELLIPSO_2'.OR. &
         CSHAPE=='ELLIPSO_3'.OR. &
         CSHAPE=='ELLIPSPBC'.OR. &
         CSHAPE=='FRMFILPBC'.OR. &
         CSHAPE=='FROM_FILE'.OR. &
         CSHAPE=='GAUSS_SPH'.OR. &
         CSHAPE=='HEXGONPBC'.OR. &
         CSHAPE=='HEX_PRISM'.OR. &
         CSHAPE=='LAYRDSLAB'.OR. &
         CSHAPE=='LYRSLBPBC'.OR. &
         CSHAPE=='MLTBLOCKS'.OR. &
         CSHAPE=='RCTGLBLK3'.OR. &
         CSHAPE=='RCTGLPRSM'.OR. &
         CSHAPE=='RCTGL_PBC'.OR. &
         CSHAPE=='RCTG_RCTG'.OR. &
         CSHAPE=='RECRECPBC'.OR. &
         CSHAPE=='SLAB_HOLE'.OR. &
         CSHAPE=='SLBHOLPBC'.OR. &
         CSHAPE=='SPHERES_N'.OR. &
         CSHAPE=='SPHRN_PBC'.OR. &
         CSHAPE=='SPHROID_2'.OR. &
         CSHAPE=='SPH_ANI_N'.OR. &
         CSHAPE=='TETRAHDRN'.OR. &
         CSHAPE=='TRILYRPBC'.OR. &
         CSHAPE=='TRNGLPRSM'.OR. &
         CSHAPE=='UNIAXICYL')THEN
         WRITE(CMSGNM,FMT='(A,A)')CSHAPE,' - Shape definition '
         CALL WRIMSG('REAPAR',CMSGNM)
      ELSE
         WRITE(CMSGNM,FMT='(A,A)')'CSHAPE=',CSHAPE
         CALL WRIMSG('REAPAR',CMSGNM)
         CALL ERRMSG('FATAL','REAPAR',                            &
                     ' Unrecognized shape directive')
      ENDIF

! Read shape parameters

! targets needing 0 shape parameters (skip one line):

      IF(CSHAPE=='FROM_FILE'.OR.CSHAPE=='ANIFRMFIL'.OR. &
         CSHAPE=='MLTBLOCKS')THEN
         READ(IOPAR,FMT=*,ERR=99)

! targets needing 1 shape parameter:

      ELSEIF(CSHAPE=='DW1996TAR'.OR.CSHAPE=='TETRAHDRN')THEN
         CWHERE='error reading 1 shape parameter in ddscat.par'
         READ(IOPAR,FMT=*,ERR=99)SHPAR(1)
         WRITE(CMSGNM,FMT='(0PF10.5,A)')SHPAR(1),' = SHPAR(1)'
         CALL WRIMSG('TARGET',CMSGNM)

! targets needing 2 shape parameters:

      ELSEIF(CSHAPE=='CYLNDRCAP'.OR.CSHAPE=='SPHROID_2'.OR. &
             CSHAPE=='UNIAXICYL')THEN
         CWHERE='error reading 2 shape parameters in ddscat.par'
         READ(IOPAR,FMT=*,ERR=99)SHPAR(1),SHPAR(2)
         DO J=1,2
            WRITE(CMSGNM,FMT='(0PF10.5,A,I1,A)')SHPAR(J),' = SHPAR(',J,')'
            CALL WRIMSG('TARGET',CMSGNM)
         ENDDO

! targets needing 2 shape parameters and file name:

      ELSEIF(CSHAPE=='FRMFILPBC'.OR.CSHAPE=='SPHERES_N'.OR. &
             CSHAPE=='SPH_ANI_N')THEN
         CWHERE='error reading 2 shape parameters and filename in ddscat.par'
         READ(IOPAR,FMT=*,ERR=99)SHPAR(1),SHPAR(2),CFLSHP
         CALL WRIMSG('TARGET','shape file=')
         CALL WRIMSG('TARGET',CFLSHP)
         DO J=1,2
            WRITE(CMSGNM,FMT='(0PF10.5,A,I1,A)')SHPAR(J),' = SHPAR(',J,')'
            CALL WRIMSG('TARGET',CMSGNM)
         ENDDO

! targets needing 3 shape parameters:

      ELSEIF(CSHAPE=='ANIELLIPS'.OR.CSHAPE=='ANI_ELL_2'.OR. &
             CSHAPE=='ANI_ELL_3'.OR.CSHAPE=='ANIRCTNGL'.OR. &
             CSHAPE=='CYLINDER1'.OR.CSHAPE=='ELLIPSOID'.OR. &
             CSHAPE=='ELLIPSO_2'.OR.CSHAPE=='ELLIPSO_3'.OR. &
             CSHAPE=='HEX_PRISM'.OR.CSHAPE=='RCTGLPRSM')THEN
         CWHERE='error reading 2 shape parameters in ddscat.par'
         READ(IOPAR,FMT=*,ERR=99)SHPAR(1),SHPAR(2),SHPAR(3)
         DO J=1,3
            WRITE(CMSGNM,FMT='(0PF10.5,A,I1,A)')SHPAR(J),' = SHPAR(',J,')'
            CALL WRIMSG('TARGET',CMSGNM)
         ENDDO

! targets needing 3 shape parameters and file name:

      ELSEIF(CSHAPE=='SPHRN_PBC')THEN
         CWHERE='error reading 3 shape parameters and filename in ddscat.par'
         READ(IOPAR,FMT=*,ERR=99)SHPAR(1),SHPAR(2),SHPAR(3),CFLSHP
         CALL WRIMSG('TARGET','shape file=')
         CALL WRIMSG('TARGET',CFLSHP)
         DO J=1,3
            WRITE(CMSGNM,FMT='(0PF10.5,A,I1,A)')SHPAR(J),' = SHPAR(',J,')'
            CALL WRIMSG('TARGET',CMSGNM)
         ENDDO
         
! targets needing 4 shape parameters

      ELSEIF(CSHAPE=='SLAB_HOLE'.OR.CSHAPE=='TRNGLPRSM')THEN
         CWHERE='error reading 4 shape parameters in ddscat.par'
         READ(IOPAR,FMT=*,ERR=99)SHPAR(1),SHPAR(2),SHPAR(3),SHPAR(4)
         DO J=1,4
            WRITE(CMSGNM,FMT='(0PF10.5,A,I1,A)')SHPAR(J),' = SHPAR(',J,')'
            CALL WRIMSG('TARGET',CMSGNM)
         ENDDO

! targets needing 5 shape parameters

      ELSEIF(CSHAPE=='CYLNDRPBC'.OR.CSHAPE=='DSKRCTNGL'.OR. &
             CSHAPE=='HEXGONPBC'.OR.CSHAPE=='RCTGL_PBC')THEN
         CWHERE='error reading 5 shape parameters in ddscat.par'
         READ(IOPAR,FMT=*,ERR=99)(SHPAR(J),J=1,5)
         DO J=1,5
            WRITE(CMSGNM,FMT='(0PF10.5,A,I1,A)')SHPAR(J),' = SHPAR(',J,')'
            CALL WRIMSG('TARGET',CMSGNM)
         ENDDO

! targets needing 6 shape parameters:

      ELSEIF(CSHAPE=='BISLINPBC'.OR.CSHAPE=='CONELLIPS'.OR. &
             CSHAPE=='GAUSS_SPH'.OR.CSHAPE=='RCTG_RCTG'.OR. &
             CSHAPE=='SLBHOLPBC'.OR.CSHAPE=='SPHROID_2')THEN
         CWHERE='error reading 6 shape parameters in ddscat.par'
         READ(IOPAR,FMT=*,ERR=99)(SHPAR(J),J=1,6)
         DO J=1,6
            WRITE(CMSGNM,FMT='(0PF10.5,A,I1,A)')SHPAR(J),' = SHPAR(',J,')'
            CALL WRIMSG('TARGET',CMSGNM)
         ENDDO

! targets needing 7 shape parameters:

      ELSEIF(CSHAPE=='DSKRCTPBC'.OR.CSHAPE=='LAYRDSLAB')THEN
         CWHERE='error reading 7 shape parameters in ddscat.par'
         READ(IOPAR,FMT=*,ERR=99)(SHPAR(J),J=1,7)
         DO J=1,7
            WRITE(CMSGNM,FMT='(0PF10.5,A,I1,A)')SHPAR(J),' = SHPAR(',J,')'
            CALL WRIMSG('TARGET',CMSGNM)
         ENDDO

! targets needing 8 shape parameters:

      ELSEIF(CSHAPE=='DSKBLYPBC'.OR.CSHAPE=='RECRECPBC')THEN
         CWHERE='error reading 8 shape parameters in ddscat.par'
         READ(IOPAR,FMT=*,ERR=99)(SHPAR(J),J=1,8)
         DO J=1,8
            WRITE(CMSGNM,FMT='(0PF10.5,A,I1,A)')SHPAR(J),' = SHPAR(',J,')'
            CALL WRIMSG('TARGET',CMSGNM)
         ENDDO

! targets needing 9 shape parameters:

      ELSEIF(CSHAPE=='LYRSLBPBC'.OR.CSHAPE=='RCTGLBLK3')THEN
         CWHERE='error reading 9 shape parameters in ddscat.par'
         READ(IOPAR,FMT=*,ERR=99)(SHPAR(J),J=1,9)
         DO J=1,9
            WRITE(CMSGNM,FMT='(0PF10.5,A,I1,A)')SHPAR(J),' = SHPAR(',J,')'
            CALL WRIMSG('TARGET',CMSGNM)
         ENDDO

! targets needing 10 shape parameters: none

! targets needing 11 shape parameters:

      ELSEIF(CSHAPE=='TRILYRPBC')THEN
         CWHERE='error reading 11 shape parameters in ddscat.par'
         READ(IOPAR,FMT=*,ERR=99)(SHPAR(J),J=1,11)
         DO J=1,11
            WRITE(CMSGNM,FMT='(0PF10.5,A,I2,A)')SHPAR(J),' = SHPAR(',J,')'
            CALL WRIMSG('TARGET',CMSGNM)
         ENDDO

      ELSE
         CALL ERRMSG('FATAL','REAPAR',                                   &
            ' No instructions for reading shape parameters for this shape')
      ENDIF

!  JPBC=0 if PBC are not used
!       1 if PBC in y direction only
!       2 if PBC in z direction only
!       3 if PBC in y and z directions

      JPBC=0

! need to set JPBC for PBC targets

      IF(CSHAPE=='BISLINPBC')THEN
         IF(SHPAR(6)>0._WP)JPBC=1
         JPBC=JPBC+2
      ENDIF
      IF(CSHAPE=='ANIFILPBC')THEN
         IF(SHPAR(1)>0._WP)JPBC=1
         IF(SHPAR(2)>0._WP)JPBC=JPBC+2
      ENDIF
      IF(CSHAPE=='CYLNDRPBC')THEN
         IF(SHPAR(4)>0._WP)JPBC=1
         IF(SHPAR(5)>0._WP)JPBC=JPBC+2
      ENDIF
      IF(CSHAPE=='DSKBLYPBC')THEN
         IF(SHPAR(7)>0._WP)JPBC=1
         IF(SHPAR(8)>0._WP)JPBC=JPBC+2
      ENDIF
      IF(CSHAPE=='DSKRCTPBC')THEN
         IF(SHPAR(6)>0._WP)JPBC=1
         IF(SHPAR(7)>0._WP)JPBC=JPBC+2
      ENDIF
      IF(CSHAPE=='FRMFILPBC')THEN
         IF(SHPAR(1)>0._WP)JPBC=1
         IF(SHPAR(2)>0._WP)JPBC=JPBC+2
      ENDIF
      IF(CSHAPE=='HEXGONPBC')THEN
         IF(SHPAR(4)>0._WP)JPBC=1
         IF(SHPAR(5)>0._WP)JPBC=JPBC+2
      ENDIF
      IF(CSHAPE=='LYRSLBPBC')THEN
         IF(SHPAR(8)>0._WP)JPBC=1
         IF(SHPAR(9)>0._WP)JPBC=JPBC+2
      ENDIF
      IF(CSHAPE=='RCTGL_PBC')THEN
         IF(SHPAR(4)>0._WP)JPBC=1
         IF(SHPAR(5)>0._WP)JPBC=JPBC+2
      ENDIF
      IF(CSHAPE=='RECRECPBC')THEN
         IF(SHPAR(7)>0._WP)JPBC=1
         IF(SHPAR(8)>0._WP)JPBC=JPBC+2
      ENDIF
      IF(CSHAPE=='SLBHOLPBC')THEN
         IF(SHPAR(5)>0._WP)JPBC=1
         IF(SHPAR(6)>0._WP)JPBC=JPBC+2
      ENDIF
      IF(CSHAPE=='SPHRN_PBC')THEN
         IF(SHPAR(2)>0._WP)JPBC=1
         IF(SHPAR(3)>0._WP)JPBC=JPBC+2
      ENDIF
      IF(CSHAPE=='TRILYRPBC')THEN
         IF(SHPAR(10)>0._WP)JPBC=1
         IF(SHPAR(11)>0._WP)JPBC=JPBC+2
      ENDIF

!=======================================================================
! following code disabled 03.01.29
! (retain for use in future noncubic version)

!   Obtain lattice anisotropy parameters DX(1-3)
!   For cubic lattice, DX(1)=DX(2)=DX(3)
!   Note that we do not require here that DX(1)*DX(2)*DX(3)=1 upon
!   input; DX is renormalized here before being returned to DDSCAT


!      READ(IOPAR,FMT=*,ERR=99)DX(1),DX(2),DX(3)
!      WRITE(CMSGNM,FMT='(F8.3,F8.3,F8.3,A)')DX,
!     &   ' = relative lattice spacings dx,dy,dz'
!      CALL WRIMSG('REAPAR',CMSGNM)

!      DELTA=(DX(1)*DX(2)*DX(3))**(1./3.)
!      DX(1)=DX(1)/DELTA
!      DX(2)=DX(2)/DELTA
!      DX(3)=DX(3)/DELTA
!      WRITE(CMSGNM,FMT='(F8.3,F8.3,F8.3,A)')DX,
!     &   ' = normalized lattice spacings dx,dy,dz'
!      CALL WRIMSG('REAPAR',CMSGNM)

! and replaced by following:

      DX(1)=1._WP
      DX(2)=1._WP
      DX(3)=1._WP

!=======================================================================

!   Obtain names of file(s) containing dielectric function(s)
!   NCOMP = number of different dielectric functions
!   CFLEPS(1-NCOMP) = names of files containing dielectric functions

      CWHERE='error reading NCOMP in ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)NCOMP
      WRITE(CMSGNM,FMT='(A,I4)')' NCOMP=',NCOMP
      CALL WRIMSG('REAPAR',CMSGNM)

! Check that MXCOMP is large enough

      IF(NCOMP>MXCOMP)CALL ERRMSG('FATAL','REAPAR', &
                      'NCOMP > MXCOMP: must increase MXCOMP in DDSCAT')

!*** Check that NCOMP=2 if CSHAPE=UNIAXICYL *******************************
!                     3           ANIELLIPS
!                     2           ELLIPSO_2
!                     3           ELLIPSO_3
!                     6           ANI_ELL_2
!                     9           ANI_ELL_3
!                     2           CONELLIPS
!                     2           SPHROID_2
!                     2-4         LYRD_SLAB
!                     1-4         LYRSLBPBC

      IF(CSHAPE=='UNIAXICYL'.AND.NCOMP/=2)CALL ERRMSG('FATAL','REAPAR', &
         'NCOMP must be 2 for option UNIAXICYL')
      IF(CSHAPE=='ANIELLIPS'.AND.NCOMP/=3)CALL ERRMSG('FATAL','REAPAR', &
         'NCOMP must be 3 for option ANIELLIPS')
      IF(CSHAPE=='ELLIPSO_2'.AND.NCOMP/=2)CALL ERRMSG('FATAL','REAPAR', &
         'NCOMP must be 2 for option ELLIPSO_2')
      IF(CSHAPE=='ELLIPSO_3'.AND.NCOMP/=3)CALL ERRMSG('FATAL','REAPAR', &
         'NCOMP must be 3 for option ELLIPSO_3')
      IF(CSHAPE=='ANI_ELL_2'.AND.NCOMP/=6)CALL ERRMSG('FATAL','REAPAR', &
         'NCOMP must be 6 for option ANI_ELL_2')
      IF(CSHAPE=='ANI_ELL_3'.AND.NCOMP/=9)CALL ERRMSG('FATAL','REAPAR', &
         'NCOMP must be 9 for option ANI_ELL_3')
      IF(CSHAPE=='CONELLIPS'.AND.NCOMP/=2)CALL ERRMSG('FATAL','REAPAR', &
         'NCOMP must be 2 for option CONELLIPS')
      IF(CSHAPE=='SPHROID_2'.AND.NCOMP/=2)CALL ERRMSG('FATAL','REAPAR', &
         'NCOMP must be 2 for option SPHROID_2')
      IF(CSHAPE=='LYRSLBPBC')THEN
         J=4
         IF(SHPAR(7)<=0._WP)J=3
         IF(SHPAR(7)<=0._WP.AND.SHPAR(6)<=0._WP)J=2
         IF(SHPAR(7)<=0._WP.AND.SHPAR(6)<=0._WP.AND.SHPAR(5)<=0._WP)J=1
         WRITE(CMSGNM,FMT='(A,I2)')'LYRSLBPBC: need NCOMP=',J
         IF(NCOMP<J)CALL ERRMSG('FATAL','REAPAR',CMSGNM)
      ENDIF
      IF(CSHAPE=='LYRD_SLAB')THEN
         J=4
         IF(SHPAR(7)<=0._WP)J=3
         IF(SHPAR(7)<=0._WP.AND.SHPAR(6)<=0._WP)J=2
         IF(SHPAR(7)<=0._WP.AND.SHPAR(6)<=0._WP.AND.SHPAR(6)<=0._WP)J=1
         WRITE(CMSGNM,FMT='(A,I2)')'LYRD_SLAB: need NCOMP=',J
         IF(NCOMP<J)CALL ERRMSG('FATAL','REAPAR',CMSGNM)
      ENDIF
      IF(CSHAPE=='SPHERES_N'.AND.NCOMP>1)THEN
         WRITE(CMSGNM,FMT='(A,A)')'SPHERES_N is for single isotropic ', &
               ' material -- try option SPH_ANI_N ?'
         CALL ERRMSG('FATAL','REAPAR',CMSGNM)
      ENDIF

!=======================================================================

      CWHERE='error reading filename(s) for diel.fn. in ddscat.par'
      DO J=1,NCOMP
         READ(IOPAR,FMT=*,ERR=99)CFLEPS(J)
         WRITE(CMSGNM,FMT='(I2,A,A)')J,' ',CFLEPS(J)
         CALL WRIMSG('REAPAR',CMSGNM)
      ENDDO

!=======================================================================

!    Define INITIALIZATION:

!    INIT = 0 to start with |X0>=0
!           1 to obtain |X0> from 4th-order expansion in polarizability
!           2 to read |X0> from file solvp.in
! disabled 08.03.12 v7.0.5 since we always use INIT=0
! skip line:
!      READ(IOPAR,FMT=*,ERR=99)CLINE
!      CALL WRIMSG(' ',CLINE)
!      CWHERE='error reading INIT in ddscat.par'
!      READ(IOPAR,FMT=*,ERR=99)INIT
!      IF(INIT==0)THEN
!         CALL WRIMSG('REAPAR','INIT=0 to start with |X0>=0 (CCG method)')
!      ELSEIF(INIT==1)THEN
!         CALL WRIMSG('REAPAR', &
!                     'INIT=1 to obtain |X0> from 4th-order expansion (CCG)')
!      ELSEIF(INIT==2)THEN
!         CALL WRIMSG('REAPAR','INIT=2 to read |X0> from file solvp.in')
!      ELSE
!         CALL ERRMSG('FATAL','REAPAR',' Wrong value of INIT')
!      ENDIF

!***********************************************************************

!    Define error tolerance:
!    TOL= maximum acceptable value of |Ax-E|/|E|

! skip line:

      READ(IOPAR,FMT=*,ERR=99)CLINE
      CALL WRIMSG('REAPAR',CLINE)
 
      CWHERE='error reading TOL in ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)TOL
! 2009.08.27 (BTD) add code to catch problem with
! we exoect TOL to be between 1e-10 and <1.
! if outside this range, there has probably been an error reading TOL
      IF(TOL.GE.1.E-10_WP.AND.TOL.LT.1.E0_WP)THEN
         GOTO 2000
      ELSE

!    Note: if the number of diel.fn. files is less than NCOMP or
!          greater than NCOMP, this will lead to subsequent errors
!          in reading ddscat.par, which will probably affect reading
!          of TOL

         WRITE(CMSGNM,FMT='(1PE10.3,A)')TOL,' = TOL'
         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(A)')                                 &
               'Appears that there has been an error reading TOL'
         CALL WRIMSG('REAPAR',CMSGNM)
         WRITE(CMSGNM,FMT='(A)')                             &
               'Check whether there are NCOMP diel.fn. files'
         CALL WRIMSG('REAPAR',CMSGNM)
         GOTO 99
      ENDIF
 2000 WRITE(CMSGNM,FMT='(1PE10.3,A)')TOL,                          &
            ' = TOL = max. acceptable normalized residue |Ax-E|/|E|'
      CALL WRIMSG('REAPAR',CMSGNM)
      IF(TOL<1.E-5_WP.OR.TOL>1.E-1_WP)THEN
         CALL ERRMSG('WARNING','REAPAR',' Strange value of TOL ')
      ENDIF

!***********************************************************************

!    Define summation limit parameter GAMMA for PBC calculations
!    summations will be carried out to distance r=2/(k*alpha)
!    with suppression factor exp[-(alpha*kr)^4]=exp[-16]
!    GAMMA is used only when JPBC=1,2, or 3

! skip line:
      READ(IOPAR,FMT=*,ERR=99)CLINE
      CALL WRIMSG(' ',CLINE)
      CWHERE='error reading GAMMA in ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)GAMMA
      IF(JPBC>0)THEN
         WRITE(CMSGNM,FMT='(1PE10.3,A)')GAMMA,                     &
               ' = GAMMA = replica dipole summation limiter for PBC'
         CALL WRIMSG('REAPAR',CMSGNM)
         IF(GAMMA>1.E-1_WP.OR.GAMMA<1.E-4_WP)THEN
            CALL ERRMSG('WARNING','REAPAR',' Strange value of GAMMA ')
         ENDIF
      ELSE
         WRITE(CMSGNM,FMT='(A)')'[GAMMA is not used in present calculation]'
         CALL WRIMSG('REAPAR',CMSGNM)
      ENDIF

!=======================================================================

!    Define angular resolution used for calculating radiation force,
!    <cos>, <cos^2>, and radiation torque (if CMTORQ='DOTORQ')

!     ETASCA = parameter controlling number of scattering angles
!            = 1 is a good choice, but can reduce this if higher
!                accuracy is required for angular averages

!     number of scattering angles will be approximately
!            6*[(3+x)/ETASCA]^2

      READ(IOPAR,FMT=*,ERR=99)CLINE
      CALL WRIMSG(' ',CLINE)
      CWHERE='error reading ETASCA in ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)ETASCA
      IF(ETASCA>1.E-3_WP)THEN
         WRITE(CMSGNM,FMT='(F7.4,A)')ETASCA, &
               ' = ETASCA (parameter controlling number of scattering angles'
         CALL WRIMSG('REAPAR',CMSGNM)
      ELSE
         WRITE(CMSGNM,FMT='(1PE10.3,A)')ETASCA, &
              ' is not an appropriate value for ETASCA: check ddscat.par'
         CALL ERRMSG('FATAL','REAPAR',CMSGNM)
      ENDIF

!=======================================================================

!    Define WAVELENGTH :

!    WAVINI = first wavelength (physical units)
!    WAVEND = last wavelength (physical units)
!    NWAV   = number of wavelengths
!    CDIVID = 'LIN', 'LOG', or 'INV' for equal increments in
!              wavelength, log(wavelength), or frequency
!           = 'TAB' to read wavelengths from file 'wave.tab'

      READ(IOPAR,FMT=*,ERR=99)CLINE
      CALL WRIMSG(' ',CLINE)
      CWHERE='error reading firstwave, lastwave, number, how in ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)WAVINI,WAVEND,NWAV,CDIVID
      IF(CDIVID/='LIN'.AND.CDIVID/='LOG'.AND.CDIVID/='INV'.AND. &
         CDIVID/='TAB')THEN
         CALL ERRMSG('FATAL','REAPAR', &
                     ' CDIVID for wavelengths must be LIN,LOG,INV, or TAB')
      ENDIF
      IF(CDIVID/='TAB')THEN
         WRITE(CMSGNM,FMT='(I3,A,F7.4,A,F7.4)')NWAV,' wavelengths from ', &
                                               WAVINI,' to ',WAVEND
         CALL WRIMSG('REAPAR',CMSGNM)
         CALL DIVIDE(CDIVID,WAVINI,WAVEND,NWAV,MXWAV,WAVEA)
      ELSE
         OPEN(UNIT=28,FILE='wave.tab',STATUS='OLD')

! skip one header line

         READ(28,*)
         NWAV=0
900      READ(28,*,END=1000,ERR=1000)WAVINI
         NWAV=NWAV+1
         IF(NWAV>MXWAV)THEN
            WRITE(CMSGNM,FMT='(A,I4,A,I4,A)')'NWAV=',NWAV,' > MXWAV=', &
                  MXWAV,' need to recompile with larger MXWAV'
            CALL WRIMSG('REAPAR',CMSGNM)
            CALL ERRMSG('FATAL','REAPAR',' NWAV > MXWAV')
         ENDIF
         WAVEA(NWAV)=WAVINI
!*** diagnostic
!         write(0,*)'reapar ckpt 1: NWAV=',NWAV,' WAVEA(NWAV)=',WAVEA(NWAV)
!****
         GOTO 900
1000     CLOSE(28)
         WRITE(CMSGNM,FMT='(I3,A,1PE11.5,A,1PE11.5)')NWAV,    &
               ' wavelengths from ',WAVEA(1),' to ',WAVEA(NWAV)
         CALL WRIMSG('REAPAR',CMSGNM)
      ENDIF

!=======================================================================

!    Define AEFF = radius of equal-volume sphere (physical units)

      READ(IOPAR,FMT=*,ERR=99)CLINE
      CALL WRIMSG(' ',CLINE)
      CWHERE='error reading firstaeff, lastaeff, number, how in ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)AEFINI,AEFEND,NRAD,CDIVID
      IF(CDIVID/='LIN'.AND.CDIVID/='LOG'.AND.CDIVID/='INV')THEN
         CALL ERRMSG('FATAL','REAPAR', &
                     ' CDIVID for radii must be LIN,LOG, or INV')
      ENDIF
      WRITE(CMSGNM,FMT='(I3,A,F7.4,A,F7.4)')NRAD,' eff. radii from ', &
                                            AEFINI,' to ',AEFEND
      CALL WRIMSG('REAPAR',CMSGNM)
      CALL DIVIDE(CDIVID,AEFINI,AEFEND,NRAD,MXRAD,AEFFA)

!=======================================================================

!    Define incident polarizations (in Lab Frame)
!    It is assumed that incident radiation is along x axis in Lab Frame

! skip one line:

      READ(IOPAR,FMT=*,ERR=99)CLINE
      CALL WRIMSG(' ',CLINE)

! Read complex polarization vector CXE01=e01 (normalize if necessary).

      CWHERE='error reading CXE01 in ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)CXE01
      E1=REAL(SQRT(CXE01(1)*CONJG(CXE01(1))))
      IF(E1/=0._WP)THEN
         CALL ERRMSG('FATAL','REAPAR',' CXE01(1) must be zero!')
      ENDIF

! Normalize:

      E1=SQRT(REAL(CXE01(2)*CONJG(CXE01(2))+CXE01(3)*CONJG(CXE01(3))))
      CXE01(2)=CXE01(2)/E1
      CXE01(3)=CXE01(3)/E1

! Construct orthogonal normalized polarization vector CXE02=e02
! using xhat cross e01* = e02

      CXE02(1)=(0._WP,0._WP)
      CXE02(2)=-CONJG(CXE01(3))
      CXE02(3)=CONJG(CXE01(2))

! IORTH=1 to calculate only for single polarization
!      =2 to also calculate for orthogonal polarization

      READ(IOPAR,FMT=*,ERR=99)IORTH
      WRITE(CMSGNM,FMT='(A,I2)')'IORTH=',IORTH
      CALL WRIMSG('REAPAR',CMSGNM)

!***********************************************************************

!    Specify whether or not to write ".sca" files
! skip line:
      READ(IOPAR,FMT=*,ERR=99)CLINE
      CALL WRIMSG('REAPAR',CLINE)

! IWRKSC=0 to NOT write ".sca" file for each target orientation
!       =1 to write ".sca" file for each target orientation

      CWHERE='error reading IWRKSC in ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)IWRKSC
      WRITE(CMSGNM,FMT='(A,I2)')'IWRKSC=',IWRKSC
      CALL WRIMSG('REAPAR',CMSGNM)

!***********************************************************************

!    Specify whether or not to write ".pol" files

! IWRPOL=0 to NOT write ".pol1" and ".pol2" file for each
!          target orientation (beta,theta)
!       =1 to write ".pol1" and ".pol2" file for each
!          target orientation (beta,theta)

      CWHERE='error reading IWRPOL in ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)IWRPOL
      WRITE(CMSGNM,FMT='(A,I2)')'IWRPOL=',IWRPOL
      CALL WRIMSG('REAPAR',CMSGNM)

!***********************************************************************

!    Read information determining target rotations

! skip line
      READ(IOPAR,FMT=*,ERR=99)CLINE
      CALL WRIMSG(' ',CLINE)
      CWHERE='error reading BETAMI,BETAMX,NBETA in ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)BETAMI,BETAMX,NBETA
      CWHERE='error reading THETMI,THETMX,NTHETA in ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)THETMI,THETMX,NTHETA
      CWHERE='error reading PHIMIN,PHIMAX,NPHI in ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)PHIMIN,PHIMAX,NPHI

! check that user has not requested more than 1000 orientations and
! IWRKSC=1

      IF(IWRKSC>0.AND.(NBETA*NTHETA*NPHI)>1000)THEN
         CALL ERRMSG('FATAL','REAPAR',                                        &
                     'error: if IWRKSC=1, NBETA*NTHETA*NPHI must be .le. 1000')
      ENDIF

      READ(IOPAR,FMT=*,ERR=99)CLINE
      CALL WRIMSG(' ',CLINE)
      CWHERE='error reading IWAV0,IRAD0,IORI0 in ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)IWAV0,IRAD0,IORI0

!*** Check that values of NBETA,NTHETA,NPHI are OK:

      IF(NBETA>MXBETA)CALL ERRMSG('FATAL','REAPAR',              &
                                  ' NBETA>MXBETA: must recompile')
      IF(NPHI>MXPHI)CALL ERRMSG('FATAL','REAPAR',            &
                                ' NPHI>MXPHI: must recompile')
      IF(NTHETA>MXTHET)CALL ERRMSG('FATAL','REAPAR',               &
                                   ' NTHETA>MXTHET: must recompile')

!***  Check that IWAV0,IRAD0,IORI0 are OK:

      IF(IWAV0+1>NWAV)CALL ERRMSG('FATAL','REAPAR','IWAV0+1 > NWAV')
      IF(IRAD0+1>NRAD)CALL ERRMSG('FATAL','REAPAR','IRAD0+1 > NRAD')
      IF(IORI0+1>NBETA*NTHETA*NPHI)CALL ERRMSG('FATAL','REAPAR',            &
                                               'IORI0+1 > NBETA*NTHETA*NPHI')

! If NPHI>1, then set IORTH=2 regardless of value input.

      IF(IORTH==1.AND.NPHI>1)THEN
         IORTH=2
         CALL WRIMSG('REAPAR','Set IORTH=2 since NPHI>1 ')
      ENDIF

      IF(IORTH==1)THEN
         CALL WRIMSG('REAPAR','Calculate only for single polarization ')
      ELSEIF(IORTH==2)THEN
         CALL WRIMSG('REAPAR','Do orthogonal polarization for each case ')
      ELSE
         CALL ERRMSG('FATAL','REAPAR',' WRONG VALUE OF IORTH ')
      ENDIF
      WRITE(CMSGNM,FMT='(2F7.2,A,I4)')BETAMI,BETAMX, &
         ' Range of BETA values ; NBETA =',NBETA
      CALL WRIMSG('REAPAR',CMSGNM)
      WRITE(CMSGNM,FMT='(2F7.2,A,I4)')THETMI,THETMX, &
         ' Range of THETA values; NTHETA=',NTHETA
      CALL WRIMSG('REAPAR',CMSGNM)
      WRITE(CMSGNM,FMT='(2F7.2,A,I4)')PHIMIN,PHIMAX, &
         ' Range of PHI values ;   NPHI =',NPHI
      CALL WRIMSG('REAPAR',CMSGNM)

! Convert from degrees to radians

      BETAMI=BETAMI/DEGRAD
      BETAMX=BETAMX/DEGRAD
      THETMI=THETMI/DEGRAD
      THETMX=THETMX/DEGRAD
      PHIMIN=PHIMIN/DEGRAD
      PHIMAX=PHIMAX/DEGRAD

!=======================================================================

!    Specify elements of scattering matrix to be printed out

      READ(IOPAR,FMT=*,ERR=99)CLINE
      CALL WRIMSG(' ',CLINE)
      CWHERE='error reading NSMELTS in ddscat.par'
      READ(IOPAR,FMT=*,ERR=99)NSMELTS

! NSMELTS = number of elements of scattering matrix element to be
!           printed out.  Must be no greater than 9
!           if NSMELTS is zero, then 6 "default" elements are output

      IF(NSMELTS>9)THEN
         WRITE(IDVOUT,FMT='(A,I2,A)')'Error:',NSMELTS,             &
                          ' > 9 is not permitted by output routines'
        CALL ERRMSG('FATAL','REAPAR',                         &
                    'NSMELTS > 9 is not allowed in ddscat.par')
      ENDIF
      IF(NSMELTS<=0)THEN
         SMIND1(1)=11
         SMIND1(2)=21
         SMIND1(3)=31
         SMIND1(4)=41
         SMIND1(5)=12
         SMIND1(6)=13
         NSMELTS=6
         READ(IOPAR,FMT=*)
      ELSE
         CWHERE='error reading indices ij of S_ij elements in ddscat.par'
         READ(IOPAR,FMT=*,ERR=99)(SMIND1(J),J=1,NSMELTS)
      ENDIF
      DO J=1,NSMELTS
         SMIND2(J)=SMIND1(J) - 10*(SMIND1(J)/10)
         SMIND1(J)=SMIND1(J)/10
      ENDDO

!=======================================================================

!    Specify scattering directions to be calculated
!    Two options:
!    CMDFRM = 'LFRAME': specify scattering directions n in Lab Frame
!                       (where incident beam is in x-direction)
!         THETAN,PHIN = Direction of n from n0 in calculation of the
!                       scattering matrix; cos(THETAN) is n0 \dot n
!                       and PHIN is azimuthal angle of n from Lab xy pla
!                       Note: THETA1, THETA2, and PHI for each scatterin
!                       plane are entered in degrees, and immediately
!                       converted to radians.  Arbitrary number
!                       of scattering planes may be considered.
!    CMDFRM = 'TFRAME': specify scattering directions in Target Frame,
!                       defined by target axes a1,a2,a3
! If JPBC=0:
!         THETAN,PHIN = Direction of n relative to a1,a2,a3:
!                       THETAN = angle between n and a1
!                       PHIN = angle between a1-n plane and a1-a2 plane
!    JPBC=1:
!         THETAN      = Diffraction order along y_TF axis
!         PHIN        = azimuthal angle around y_TF
!    JPBC=2:
!         THETAN      = Diffraction order along z_TF axis
!         PHIN        = azimuthal angle around z_TF
!    JPBC=3
!         THETAN      = Diffraction order along y_TF
!         PHIN        = Diffraction order along z_TF
!      where we first run through transmitted waves 1 -> NSCAT/2
!            and then run through reflected waves NSCAT/2+1 -> NSCAT

! Three cases:
! JPBC=0: single isolated target.
!         specify
!            phi for scattering plane
!            thetamin, thetamax, dtheta for scattering plane

! JPBC=1,2: periodic in one dimension.
!         specify
!            diffraction order in direction of target periodicity
!            phimin, phimax, dphi for scattering cone

! JPBC=3: periodic in two dimensions
!         specify
!            difraction order in y direction and order in z direction

      NSCA=0
      IF(JPBC==0)THEN
         READ(IOPAR,FMT=*,ERR=99)CLINE
         CWHERE='error reading CMDFRM (only LFRAME or TFRAME allowed)'
         READ(IOPAR,FMT=*,ERR=99)CMDFRM
         IF(CMDFRM/='LFRAME'.AND.CMDFRM/='lframe'.AND. &
            CMDFRM/='TFRAME'.AND.CMDFRM/='tframe')GOTO 99
         CALL WRIMSG('REAPAR',CLINE)
         CWHERE='JPBC=0, error reading NPLANES = number of scat. planes'
         READ(IOPAR,FMT=*,ERR=99)NPLANES
         WRITE(CMSGNM,FMT='(I4,A)')NPLANES,                       &
                                   ' = number of scattering planes'
         CWHERE='JPBC=0, error reading scattering directions in ddscat.par'
         IF(NPLANES>0)THEN
            DO JPLANE=1,NPLANES
               READ(IOPAR,FMT=*,ERR=99,END=99)PHI1,THETA1,THETA2,DTHETA
               WRITE(CMSGNM,FMT='(3F7.1,A,I4)')PHI1,THETA1,THETA2,    &
                     ' = phi, theta_min, theta_max for scatt. plane', &
                     JPLANE
               CALL WRIMSG('REAPAR',CMSGNM)
               IF(ABS(THETA1-THETA2)>0._WP.AND.DTHETA==0._WP)THEN
                  CALL ERRMSG('FATAL','REAPAR','DTHETA=0 in ddscat.par!')
               ENDIF

! Convert to radians:

               PHI1=PHI1/DEGRAD
               THETA1=THETA1/DEGRAD
               THETA2=THETA2/DEGRAD
               DTHETA=DTHETA/DEGRAD

! Allow for possibility that user did not enter correct sign for
! DTHETA (if I did it, others will too...)

               IF(THETA2<THETA1)DTHETA=-ABS(DTHETA)

! compute theta values for this scattering plane

               NSCA0=1
               DELTA=THETA2-THETA1
               IF(DELTA/=0)NSCA0=1+NINT(DELTA/DTHETA)

! Check that there is enough storage space:

               IF((NSCA+NSCA0)>MXSCA)THEN
                  WRITE(IDVOUT,FMT='(I6,I6,A)')MXSCA,NSCA0,       &
                        ' =MXSCA,NSCA0 (dirs in this scatt. plane)'
                  CALL ERRMSG('FATAL','REAPAR','MXSCA too small in DDSCAT')
               ENDIF

               THETAN(NSCA+1)=THETA1
               IF(NSCA0>2)THEN
                  DO J=2,NSCA0-1
                     THETAN(NSCA+J)=THETA1+(J-1)*DTHETA
                  ENDDO
               ENDIF

               THETAN(NSCA+NSCA0)=THETA2
               DO J=1,NSCA0
                  PHIN(NSCA+J)=PHI1
               ENDDO

               NSCA=NSCA+NSCA0
               WRITE(CMSGNM,FMT='(I4,A)')NSCA0,                           &
                  ' = number of scattering angles in this scattering plane'
               CALL WRIMSG('REAPAR',CMSGNM)
            ENDDO
         ENDIF

!=======================================================================

      ELSEIF(JPBC==1.OR.JPBC==2)THEN

! enter diffraction order, phimin, phimax, dphi for cone
! we will store order in PHIN
!               phi   in THETAN

         READ(IOPAR,FMT=*,ERR=99)CLINE
         CWHERE='error reading CMDFRM (only LFRAME or TFRAME allowed)'
         READ(IOPAR,FMT=*,ERR=99)CMDFRM
         IF(CMDFRM/='LFRAME'.AND.CMDFRM/='lframe'.AND.  &
            CMDFRM/='TFRAME'.AND.CMDFRM/='tframe')GOTO 99
         IF(CMDFRM=='LFRAME'.OR.CMDFRM=='lframe')CALL ERRMSG('FATAL', &
            'REAPAR',' cannot use LFRAME when JPBC=1 or 2')
         CALL WRIMSG(' ',CLINE)
         CWHERE='JPBC=1 or 2, error reading number of scattering cones'
         READ(IOPAR,FMT=*,ERR=99)NPLANES
         WRITE(CMSGNM,FMT='(I4,A)')NPLANES,                      &
                                   ' = number of scattering cones'
         CALL WRIMSG('REAPAR',CMSGNM)
         CWHERE='JPBC=1 or 2, error reading scattering cones in ddscat.par'
         IF(NPLANES>0)THEN
            DO JPLANE=1,NPLANES
32             READ(IOPAR,FMT=*,ERR=99,END=99)PHI1,THETA1,THETA2,DTHETA
               WRITE(CMSGNM,FMT='(3F7.1,A)')PHI1,THETA1,THETA2,       &
                     ' = order, zeta_min, zeta_max for scattering cone'
               CALL WRIMSG('REAPAR',CMSGNM)
! Convert to radians
               THETA1=THETA1/DEGRAD
               THETA2=THETA2/DEGRAD
               DTHETA=DTHETA/DEGRAD
               IF(THETA2<THETA1)DTHETA=-ABS(DTHETA)

! compute zeta values for this scattering cone

               NSCA0=1
               DELTA=THETA2 - THETA1
               IF(DELTA/=0)NSCA0=1+NINT(DELTA/DTHETA)

! Check that there is sufficient storage

               IF((NSCA+NSCA0)>MXSCA)THEN
                  WRITE(IDVOUT,FMT='(I6,I6,A)')MXSCA,NSCA0,       &
                        ' =MXSCA,NSCA0 (dirs in this scatt. plane)'
                  CALL ERRMSG('FATAL','REAPAR','MXSCA too small in DDSCAT')
               ENDIF
               THETAN(NSCA+1)=THETA1
               IF(NSCA0>2)THEN
                  DO J=2,NSCA0 - 1
                     THETAN(NSCA+J)=THETA1+(J-1)*DTHETA
                  ENDDO
               ENDIF
               THETAN(NSCA+NSCA0)=THETA2
               DO J=1,NSCA0
                  PHIN(NSCA+J)=PHI1
               ENDDO
               NSCA=NSCA+NSCA0
               WRITE(CMSGNM,FMT='(I4,A)')NSCA0,                             &
                     ' = number of scattering angles in this scattering cone'
               CALL WRIMSG('REAPAR',CMSGNM)
            ENDDO
         ENDIF

!=======================================================================

      ELSEIF(JPBC==3)THEN

! Target is doubly periodic: scattering only in discrete directions
! enter diffraction order_y, order_z
! we will store order_y in PHIN
!               order_z in THETAN
! every order_y,order_z case will first be used for transmission, then
! again for reflection

         READ(IOPAR,FMT=*,ERR=99)CLINE
         CWHERE='error reading CMDFRM (only LFRAME or TFRAME allowed)'
         READ(IOPAR,FMT=*,ERR=99)CMDFRM
         IF(CMDFRM=='LFRAME'.OR.CMDFRM=='lframe')CALL ERRMSG('FATAL', &
            'REAPAR',' cannot use LFRAME when JPBC=3')
         CALL WRIMSG(' ',CLINE)
         NSCA=0
         CWHERE='error reading number of transmission diffr. orders for JPBC=3'
         READ(IOPAR,FMT=*,ERR=99)NPLANES
         WRITE(CMSGNM,FMT='(I4,A)')NPLANES,                     &
               ' = number of diffraction orders for transmission'
         CALL WRIMSG('REAPAR',CMSGNM)
         CWHERE='error reading diffraction orders in ddscat.par for JPBC=3'
         IF(NPLANES>0)THEN
            DO JPLANE=1,NPLANES  
               READ(IOPAR,FMT=*,ERR=99,END=99)PHIN(NSCA+1),THETAN(NSCA+1)
               NSCA=NSCA+1
            ENDDO
            DO J=1,NSCA
               PHIN(J+NSCA)=PHIN(J)
               THETAN(J+NSCA)=THETAN(J)
            ENDDO
            NSCA=2*NSCA
         ENDIF
      ENDIF
!=======================================================================
      CLOSE(IOPAR)
!*** diagnostic
!      write(0,*)'reapar ckpt 9: MXNX,MXNY,MXNZ=',MXNX,MXNY,MXNZ
!***
      RETURN
99    CONTINUE
      CALL WRIMSG('REAPAR',CWHERE)
      CALL ERRMSG('FATAL','REAPAR',' Error reading ddscat.par file')
      RETURN
    END SUBROUTINE REAPAR
