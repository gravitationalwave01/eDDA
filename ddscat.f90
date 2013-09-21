    PROGRAM DDSCAT
!-------------------------------- v7.1.0 -------------------------------
      USE DDPRECISION,ONLY: WP
      USE DDCOMMON_1,ONLY: AKR,DX
      USE DDCOMMON_2,ONLY: CXADIA
      USE DDCOMMON_3,ONLY: CXZC
      USE DDCOMMON_4,ONLY: CXZW
      USE DDCOMMON_5,ONLY: IOCC
      USE DDCOMMON_6,ONLY: GAMMA,PYD,PZD,MXNATF,MXNXF,MXNYF,MXNZF,NAT,NAT3, &
                           NAT0,NX,NY,NZ,MXN3F,IDVOUT,IPBC
      USE DDCOMMON_7,ONLY: CXAOFF
      USE DDCOMMON_8,ONLY: CMDFFT
      USE DDCOMMON_10,ONLY: MYID
      IMPLICIT NONE

!-----------------------------------------------------------------------

!           Instructions for Enabling or Disabling MPI
!           ==========================================

! To enable MPI use (this requires that you have installed the "full"
! version of DDSCAT):
!                      uncomment INCLUDE 'mpif.h' statement
!                      comment out INTEGER MPI_COMM_WORLD statement
!                      compile mpi_subs.f

! For non-MPI version (either the "plain" version of DDSCAT, or the
! "full" version with MPI support disabled:

!                      comment out INCLUDE 'mpif.h' statement
!                      uncomment INTEGER MPI_COMM_WORLD statement
!                      compile mpi_fake.f

#ifdef mpi
      INCLUDE 'mpif.h'
#endif
#ifndef mpi
      INTEGER :: MPI_COMM_WORLD
#endif

!-----------------------------------------------------------------------

!**********************************************************************

!    DDSCAT is a program to use the Discrete Dipole Approximation (DDA)
!    to calculate the absorption and scattering properties of targets
!    of arbitrary geometry and dielectric properties.  The only
!    approximation is that the target may be approximated by an array
!    of point dipoles on a cubic lattice.
!    DDA theory and validity criteria are discussed in
!    Draine, B.T. 1988,
!       "The Discrete-Dipole Approximation and its Application to
!       Interstellar Graphite Grains",
!       Astrophys. J., 333, pp. 848-872
!    Draine, B.T., & Flatau, P.J. 1994,
!       "Discrete dipole approximation for scattering calculations",
!       J. Opt. Soc. Am. A, 11, 1491-1499

!    The original version of DDSCAT was developed by B. T. Draine,
!    Princeton University Observatory.

!    Subsequent versions  of DDSCAT were jointly developed by
!    B. T. Draine (Princeton University Observatory) and
!    P. J. Flatau (University of California, San Diego, Scripps Inst.
!                  of Oceanography, La Jolla)

!    DDSCAT.7.0.4 includes the following features:
!    1. Optional generation of target arrays for to approximate several
!       standard target shapes, including
!       a. ellipsoid (sphere is special case) ("ELLIPSOID")
!       b. cylindrical prism ("CYLINDER1")
!       c. rectangular prism ("RCTGLPRSM")
!       d. hexagonal prism ("HEX_PRISM")
!       e. regular tetrahedron ("TETRAHDRN")
!       f. cylindrical prism with anisotropic diel. tensor ("UNIAXICYL")
!          (uniaxial material with crystal axis = cylinder axis)
!       g. ellipsoid with anisotropic diel. tensor ("ANIELLIPS")
!       h. two touching ellipsoids with either isotropic or aniostropic
!          dielectric tensors ("ELLIPSO_2" or "ANI_ELL_2")
!       i. three thouching ellipsoids with either isotropic or aniostrop
!          dielectric constants ("ELLIPSO_3" or "ANI_ELL_3")
!       j. two concentric ellipsoids ("CONELLIPS")
!       k. target consisting of cubic blocks ("MULTBLOKS")
!       l. target consisting of union of N spheres ("SPHERES_N")
!       m. triangular prism ("TRNGLPRSM")
!       n. layered slab ("LAYRDSLAB")
!       ... and additional target options (see UserGuide)
!    2. Alternatively, program may read in target array properties from
!       ascii file ("FROM_FILE")
!    3. User may easily specify target orientation with respect to
!       incident plane wave.
!    4. Solution is found iteratively using Complex Conjugate Gradient
!       methods.
!    5. User may select FFT implementation:
!       a. GPFA code of Temperton (recommended for general use).
!       b. FFTW code of Frigo and Johnson
!       c. Convex vector library implementation (if on Convex);
!    6. User may select prescription for determining dipole
!       polarizabilities:
!       a. 'LATTDR' = Lattice Dispersion Relation of 
!                     Draine & Goodman 1993
!       b. 'GKDLDR' = Lattice Dispersion Relation of
!                     Gutkowicz-Krusin & Draine 2004
!    7. Automatic computation of cross sections for:
!       a. absorption;
!       b. extinction;
!       c. scattering;
!       d. vector radiation pressure;
!       e. vector torque.
!    8. User may easily specify scattering directions for which
!       scattering matrixs is to be computed.
!    9. Automatic computation of orientational averages, with
!       orientation weights provided by simple subroutine.
!   10. Multicomponent targets are allowed.
!   11. Anisotropic dielectric tensors are allowed, although with
!       the restriction that the principal axes must be aligned with
!       the xyz directions in the "target frame".
!   12. Wavelength-dependent dielectric functions are read from user-
!       provided files, with automatic interpolation if necessary.
!   13. Automatic looping over radii and wavelengths if desired.
!   14. User may select up to nine Mueller matrix elements for
!       computation.
!   15. DDSCAT.6.0 introduces support for MPI (Message Passing Interface
!       for parallel computations of different scattering orientations.
!   16. New with DDSCAT.6.1:
!       Automatic selection of scattering directions for computation
!       of angular averages: radiation pressure force, <cos>, <cos^2>,
!       and radiation torque.  User sets a single parameter ETASCA
!       which controls angular resolution.  Reducing ETASCA will lead
!       to greater accuracy in calculation of angular averages.
!       ETASCA=1 should provide accuracy of better than 1% for angular
!       averages.
!   17. New with DDSCAT 6.2:
!       gaussian sphere targets
!       arbitrary orientation of anisotropic constituent material
!   18. New with DDSCAT 7.0:
!       support for infinite periodic targets
! History:
! Note: Overall history of significant changes to the DDSCAT package may
!       be found in the comments in subroutine version.f .
! Here we limit comments to changes to this program module, DDSCAT.f:
! 90.11.30 (BTD): Removed LSC3 portion of CXSC complex scratch space.
! 90.12.03 (BTD): Reordered vectors.
!                 Removed LSC2 and LSC1 portions of CXSC scratch space.
! 90.12.21 (BTD): Added code to create output files
!                 'qtable' (summary of Q values)
!                 'mtable' (summary of diel. func. for material 1)
! 91.01.03 (BTD): Added MXBETA,MXPHI,MXTHET to argument list for REAPAR
! 91.01.05 (BTD): Changed I4 -> I7 in format statements 9010,9011.
!                 Set IOSHP=-1 to suppress printing of target.out file.
! 91.05.08 (BTD): Added CALPHA to argument lists for REAPAR and ALPHA
!                 and added CALPHA to various WRITE statements.
! 91.05.23 (BTD): Added IWRKSC to arg. list for REAPAR and use it to
!                 control creation of wAArBBkCCC.sca files
!                 Move lines calculating G.
!                 Move code writing to qtable and mtable.
!                 Added IDVOUT to argument list for TARGET
! 91.08.14 (BTD): Add QBKSCA to arg. list for GETFML
!                 Modify code to print out Q_bk
! 91.08.15 (BTD): Provide different headings for qtable depending on
!                 whether IORTH=1 or 2 (add statement 9043).
! 91.09.17 (BTD): Remove calls to ALPHA and EVALA (needed to move these
!                 to subroutine GETFML since alpha from Lattice
!                 Dispersion Relation depends on propagation direction
!                 and polarization state)
!                 Add CALPHA,CXEPS,ICOMP,MXCOMP to argument list for
!                 GETFML (information required by ALPHA)
! 91.11.12 (BTD): Added variable IDVSHP to argument list for TARGET
!                 (device number for REASHP to use in reading shape file
! 93.01.15 (BTD): Divided output file qtable into 2 files:
!                 qtable (containing Q_ext,Q_abs,Q_sca,g,Q_bk) and
!                 qtable2 (containing Q_pha, Q_pol, Q_cpol)
! 93.01.16 (BTD): Modify so that qtable, qtable2, and mtable are closed
!                 after writing, and reopened for each new write
! 93.01.20 (BTD): Add MXNX, MXNY, MXNZ to argument list for subroutine
!                 EXTEND to permit checking of target size against
!                 maximum allowed dimensions
! 93.03.11 (BTD): Deleted all code associated with unused variable
!                 CMACHN (originally included to identify machine/OS)
! 93.03.12 (BTD): Changed CDESCR*60 -> CDESCR*67
!                 Moved WRITE(IDVOUT,9010) to subr. TARGET
! 93.06.02 (BTD): Corrected error in FORMAT statement 9045
! 93.06.03 (BTD): Corrected error in computation of angle-averaged
!                 backscattering cross section QBKSUM(JO) (had neglected
!                 to reset sum to zero for each new target).
! 93.09.28 (BTD): Add "fix" to work around Sun compiler/OS bug (see
!                 description below.
! 93.12.15 (BTD): Corrected calculation of PPOL.  Added comments on
!                 correspondence between our notation and elements
!                 of 2x2 amplitude scattering matrix and 4x4 Mueller
!                 scattering matrix.
! 94.01.27 (BTD): Replaced SHPAR1,SHPAR2,SHPAR3 by SHPAR(1-6) to
!                 allow up to 6 shape parameters
! 94.06.20 (PJF): Comment about NEWTMP FFT and corresponding changes
!                 in reapar, eself, extend.  Add FFT timing code
!                 TSTFFT to distribution and changes to cxfft3 and
!                 timeit routines.
! 94.12.20 (BTD): Introduce subroutine VERSION to set value of
!                 variable CSTAMP.  version.f will now serve as
!                 location to maintain log of significant code
!                 changes.
! 95.06.15 (BTD): Changed CSTAMP from character*11 to character*22
!                 to include a date.
! 95.06.14 (BTD): Create version 5x.1 for internal use.
!                 Version 5x.1 computes vector force and torque
!                 on grain due to incident radiation.
!                 Changed QSCATG(2) to QSCAG(3,2) to allow for
!                 asymmetric scattering by target.
!                 Added variable QTRQSC(1-3,1-2) and added to
!                 argument list of GETFML
!                 Replaced variable GSUM(1-2) to QSCGSUM(1-3,1-2)
!                 Added variables QTRQAB(1-3,1-2) and
!                 QTRQABSUM(1-3,1-2)
! 95.06.20 (PJF+: Introduce changes into DDSCAT in order to support
!           BTD): use of modular iterative solvers (e.g., CCG).
!               : add CMDSOL to argument list to specify method
!                 of solution
!               : add CMDTRQ to argument list of REAPAR and GETFML
! 95.07.10 (BTD): deleted redundant declarations of NAT,NAT0,NAT3,
!                 NX,NY,NZ (already declared in COMMON/M6/
! 95.07.20 (BTD): added new scratch array SCRRS2(MXNAT) for use in
!                 torque calculations in SCAT
! 95.07.27 (BTD): Added new scratch variable CXSCR1 for GETFML to
!                 use for storing first polarization solution while
!                 solving for second.
! 96.01.26 (BTD): Replaced IX0,IY0,IZ0 by IXYZ0
!                 (with simultaneous change in TARGET)
! 96.02.23 (BTD): Modified format statements to get additional digits
!                 for a_1, a_2, k vector, and pol. vectors.
! 96.11.05 (BTD): Modified to remove various computations to
!                 subroutines GETMUELLER and WRITESCA
! 96.11.14 (PJF): Add TIMERS, MXTIMERS, NTIMERS
!                 Add cbinflag, cbinfile, IOBIN (binary file)
!                 Add cbinflag to formal parameters of REAPAR routine
!                 Remove GETSET routine calls
! 96.11.15 (PJF): Add cnetflag, cnetfile to reapar call. Add netCDF clos
! 96.11.21 (BTD): Declared IDNC,NTIMERS,etc.
!                 Removed call to NCCLOS
!                 Added call to WRITESCA in order to call WRITENET in
!                 order to call NCCLOS in order to close netCDF file
! 96.12.02 (BTD): initialized smori to zero before looping over
!                 orientations
!                 eliminated variables S1212,S222,CX1112,CX1122,
!                 CX1221,CX1222,CX2122
! 97.07.24 (BTD): Changed argument list for GETMUELLER, replacing
!                 PHI by PHIN, in order to correct error in computation
!                 of Mueller matrix elements and amplitude scattering
!                 matrix elements (see GETMUELLER).
! 97.11.01 (BTD): Modified to allow noncubic rectangular lattice.
!                 Information on lattice anisotropy is contained
!                 in variable DX(1-3)=(dx/d,dy/d,dz/d) where
!                 dx,dy,dz=lattice spacings in x,y,z directions
!                 d=(dx*dy*dz)**(1/3)
!                 Note that by definition, DX(1)*DX(2)*DX(3)=1
!                 Pass DX through COMMON/M1/ for use by subroutines
!                 MATVEC and CMATVEC
! 97.12.25 (BTD): Added CXAOFF to COMMON/M7/
! 98.01.13 (BTD): Added DX to argument list in last 2 calls to WRITESCA
! 98.05.01 (BTD): Removed DATA/IDVOUT/0 since IDVOUT is in COMMON
!                 Reduced number of continuation lines in CALL WRITESCA
!                 statements
! 98.08.09 (BTD): Slight change to FORMAT statement 9020
! 98.12.12 (BTD): Add NSMELTS,SMIND1,SMIND2 to argument lists of REAPAR
!                 and WRITESCA, to allow selection of Muller matrix
!                 elements
! 98.12.21 (BTD): changed dimension of CFLPAR from CHARACTER*40 to
!                 CHARACTER*60 to allow longer file names
!                 (also changed in reapar.f and dielec.f)
! 99.01.26 (BTD): modify so we can restart calculation for specified
!                 IWAV0,IRAD0,IORI0
!                 now need to input IWAV0,IRAD0,IORI0 through ddscat.par
! 00.06.12 (BTD): added CFLSHP to argument list of REAPAR to support
!                 option SPHERES_N (need to pass CFLSPH from REAPAR to
!                 TARGET)
! 03.02.13 (BTD): added 1 more digit of precision to qtable and qtable2
!                 outputs; modified header line in stmnts 9044,9046
! 03.04.13 (BTD): added character variable CFLLOG to contain name of
!                 running output file
!                 added call to NAMID to generate a unique name for outp
!                 file, to support MPI use.
!                 Unfortunately, when this is done, output to device 0
!                 is no longer unbuffered.
!                 If it is necessary to reenable unbuffered output for
!                 debugging purposes,
!                 OPEN(UNIT=IDVOUT,FILE=CFLLOG)
!                 should be commented out
!                 Added variables ITNUM(2) and MXITER to argument lists
!                 for GETFML and WRITESCA
! 03.07.13 (BTD): Added new variable ITNUMMX(2) to store maximum number
!                 of iterations taken for any of the orientations for
!                 each size and wavelength, for output in orientational
!                 average files waarbbkccc.sca
! Develop version 6.1
! 03.10.23 (BTD): Removed ICTHM and IPHIM from argument lists of
!                 REAPAR, GETFML, SHARE1, and WRITESCA -- these
!                 variables are no longer used, since SCAT now
!                 determines scattering directions automatically
!                 Add ETASCA to argument list of REAPAR, SHARE1,
!                    GETFML, WRITESCA
!                 Added NAVG to argument list of GETFML and WRITESCA
!                 Add IWRKSC to argument list of SHARE1
! 04.02.22 (BTD): Add IWRKSC to argument list of WRITESCA and SHARE1
! 04.02.25 (BTD): Added final call to WRIMSG reporting normal
!                 termination
! 04.04.01 (BTD): Added NPY,NPZ to COMMON/M6/ to support periodic
!                 boundary condition option
! 04.04.04 (BTD): Added NPY,NPZ to argument list of TARGET to be
!                 compatible with new version of TARGET
! 04.05.21 (BTD): One of the calls to WRITESCA had QSCGSUM and
!                 QSCG2SUM in incorrect order. Fixed.
! 04.05.22 (BTD): Added IF(CBINFLAG.NE.'NOTBIN') conditional
!                 before closing dd.bin (file never opened otherwise)
! 04.09.14 (BTD): Add new variables THETADF,PHIDF,BETADF to specify
!                 orientation of "Dielectric Frame" for each dipole.
! 05.06.16 (BTD): Replaced integer variables NPY,NPZ by
!                 real variables PYD,PZD in
!                 COMMON/M6/
!                 argument list of TARGET
! 05.08.04 (BTD): Added new character variable CMDFRM to be read by
!                 subroutine REAPAR from input file ddscat.par.
!                 CMDFRM=LFRAME : angles THETAN,PHIN are relative to
!                                 Lab Frame (xlab,ylab,zlab)
!                 CMDFRM=TFRAME : angles THETAN,PHIN are relative to
!                                 Target Frame (a1,a2,a3)
!                 Modified DDSCAT.f to define ENSC,EM1,EM2 accordingly
! 05.09.26 (BTD): Corrected bug in calculation of scattering for
!                 option CMDFRM = TFRAME
! 05.10.11 (BTD): Modified to support up to 1000 wavelengths:
!                 changed MXWAV from 100 to 1000
!                 changed CFLAVG*13 to CFLAVG*14
!                 changed CFLSCA*14 to CFLSCA*15
!                 (also changes in NAMER and WRITESCA)
! 05.10.11 (BTD): Added CMDFRM to argument list of WRITESCA
! 05.10.18 (BTD): Modified to allow use of more than 1000 orientations
!                 if IWRKSC=0 (previously would exceed array limits
!                 in subroutine NAMER with IORI>1000).
! 06.04.14 (BTD): Modified to call routine WRITEPOL to write files
!                 wxxxryykzzz.pol1 and wxxxryykzzz.pol2 with
!                 polarization array and other information
! 06.09.21 (BTD): Modified to add PYD,PZD, and DX to argument list of
!                 WRITEPOL
! 06.09.28 (BTD): *** Version 7.0.0 ***
! 06.10.05 (BTD): Added call to PBCSCAVEC to prepare scattering
!                 vectors for PBC targets
! 06.12.24 (BTD): Added ORDERM,ORDERN to argument list of subroutine
!                 GETMUELLER
! 06.12.28 (BTD): Added arrays XLR(3),YLR(3),ZLR(3),A3(3)
!                 Added XLR,YLR,ZLR to argument list of PBCSCAVEC
! 07.01.18 (BTD): Added A3(3),IWRPOL,IPBC,JPBC,MXNX,MXNY,MXNZ,MXPBC,
!                 ORDERM(MXSCA),ORDERN(MXSCA),
!                 PYD,PYDDX,PZD,PZDDX to argument list of SHARE1
! 07.06.20 (BTD): v7.0.2:
!                 * define X0(3)
!                 * add X0 to argument list of TARGET
!                 * add X0 to argument list of EXTEND
!                 * add X0 to argument list of GETFML
!                 * add X0 to argument list of WRITEPOL
!                 * add X0 to argument list of SHARE1
! 07.06.22 (BTD): removed THETAN and PHIN from argument list of GETFML
!                 (not needed, because scattering directions are
!                  specified through vector AKSR)
! 07.06.30 (BTD): moved CMDFFT from COMMON/M6/... CMDFFT
!                 to COMMON/M8/CMDFFT
!                 moved PYD,PZD to beginning of COMMON/M6/
! 07.07.03 (BTD): add ENSC to argument list of GETMUELLER
!                 add ENSC to argument list of PBCSCAVEC
! 07.07.08 (BTD): add EM1,EM2 to argument list of PBCSCAVEC
!                 add EM1,EM2 to argument list of GETMUELLER
! 07.08.04 (BTD): Version 7.0.3
!                 * replaced COMMON/M1/ with USE MODULE DDCOMMON_1
!                 * replaced COMMON/M2/ with USE MODULE DDCOMMON_2
!                   add dynamic allocation of CXADIA
!                 * replaced COMMON/M3/ with USE MODULE DDCOMMON_3
!                   add dynamic allocation of CXZC
!                 * replaced COMMON/M4/ with USE MODULE DDCOMMON_4
!                   add dynamic allocation of CXZW
!                 * replaced COMMON/M5/ with USE MODULE DDCOMMON_5
!                   add dynamic allocation of IOCC
!                 * replaced COMMON/M6/ with USE MODULE DDCOMMON_6
!                 * replaced COMMON/M7/ with USE MODULE DDCOMMON_7
!                   add dynamic allocation of CXAOFF
!                 * replaced COMMON/M8/ with USE MODULE DDCOMMON_8
! 07.08.05 (BTD): * Modified to do dynamic memory allocation based
!                   on actual target size
!                 * Read size upper bounds MXNX,MXNY,MXNZ from
!                   ddscat.par
!                 * Add MXNX,MXNY,MXNZ to argument list of REAPAR
!                 * Add output line reporting whether single- or
!                   double-precision version used
!                 * Add output to report memory allocation, so that
!                   if system memory limits prevent execution, there
!                   will be an indication of this in log file
! 07.08.31 (BTD): * Added call to new routines NAMER2 and WRITEFML
!                   to write out scattering amplitudes f_ml
! 07.10.07 (BTD): * Added THETA, BETA to argument list of PBCSCAVEC
! 07.10.09 (BTD): * Added EM1R,EM2R to argument list of GETMUELLER
! 07.10.24 (BTD)  * Defined new array ENSCR
!                   Added ENSCR to argument list of GETMUELLER
! 07.10.27 (BTD)  v7.0.4
!                 * changed CSHAPE*6 to CSHAPE*9
!                 * changed shape names to 9 character strings
!                 * changed SHPAR(6) -> SHPAR(10)
! 07.10.28 (BTD)  * eliminated CDIEL -- reading from tables is now standard
! 08.01.06 (BTD)  * added comments
!                 * minor streamlining since CONVEX and TMPRTN are no
!                   longer options
!                 * cosmetic changes
! 08.01.12 (BTD)  * add PYDDX,PZDDX to argument list of WRITESCA
! 08.01.13 (BTD)  * introduce variable NRWORD = length (bytes) of real word
!                 * add to argument list of WRITEPOL, so that this can
!                   be stored in file written by WRITEPOL for sanity check
! 08.01.17 (BTD)  * add new variable IANISO to indicate whether target
!                   is isotropic (IANISO=0), anisotropic with optical
!                   axes || (x,y,z)_TF (IANISO=1), or generally anisotropic
!                   (IANISO=2)
!                 * add IANISO to argument list of TARGET
!                 * modify TARGET to set value of IANISO
!                 * modifications to support changes to WRITEPOL
! 08.01.21 (BTD)  * add IANISO to argument list of SHARE1
! 08.02.01 (BTD)  * changed SHPAR(10) -> SHPAR(12) to allow up to 12 
!                   shape parameters to be passed to TARGET
! 08.02.17 (BTD)  * add call to routine SHARE0 to share dimensioning
!                   information to all processes
!                 * change so ALLOCATION is done by all processes
!                 * add DAEFF to argument list of SHARE1
!                 * add CLOSE(11) to close mtable at end
! 08.03.11 (BTD)  ver7.0.5
!                 * added ALPHA to argument list of REAPAR
!                 * added ALPHA to argument list of SHARE1
!                 * added ALPHA to DDCOMMON_6 (to communicate with CPROD)
! 08.03.15 (BTD)  * added DDCOMMON_10 to transfer MYID to DIRECT_CALC
! 08.04.17 (BTD)  * changed MXRAD from 100 to 1000
! 08.04.19 (BTD)  * changed notation: ALPHA -> GAMMA
!                 * changed order of arguments in argument list of
!                   * REAPAR
!                   * SHARE1 (in mpi_subs.f90)
!                   * GETFML
! 08.05.01 (BTD)  * added MYID to argument list of SHARE0
!                 * added MYID to argument list of SHARE1
! 08.05.09 (BTD)  * added LACE,LAXI,LCLM,LGI,LPI,LQI,LSC0 to argument list
!                   of SHARE0
!                 * added MYID to argument list of WRITESCA
!                 * moved block of code calculating scratch array positions
!                   LACE, etc. so that it is executed for all MYID values
! 08.05.10 (BTD)  * Introduce new arrays
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
! 08.05.12 (BTD) ver7.0.6
!                 * changes introduced by Art Lazanoff, NASA Ames:
!                   added loops initializing IXYZ0,CXZC,CXZW,CXSCR1,SCRRS1
!                   to zero, but BTD does not see why this is helpful...
! 08.05.29 (BTD)  * changed declarations to allow single call to MPI_REDUCE
!                   in routine COLSUM:
!                   SM(MXSCA,4,4)      -> SM(4,4,MXSCA)
!                   SMORI(MXSCA,4,4)   -> SMORI(4,4,MXSCA)
!                   SMORI_1(MXSCA,4,4) -> SMORI_1(4,4,MXSCA)
! 08.07.22 (BTD)  * added XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN to argument list
!                   of SHARE1
! 08.07.23 (BTD)  * removed some MPI_BARRIER calls, added others just
!                   before SHARE1 and SHARE2
! 08.07.27 (BTD)  * changed final allocation of BETADF,PHIDF,THETADF
!                   from NAT0 -> MXNAT
! 08.08.29 (BTD) v7.0.7
!                 * removed variable CNETFLAG and CNETFILE
!                 * removed CNETFLAG from argument list of REAPAR
!                 * removed CNETFLAG from argument list of SHARE1
!                 * removed CNETFLAG and CNETFILE from argument list of
!                   WRITESCA
! 09.09.11 (BTD) v7.0.8
!                 * added variable NCOMP_NEED to communicate with REASHP
!                   via subroutine TARGET
!                 * added check that NCOMP_NEED does not exceed NCOMP
! 10.01.28 (BTD) v7.1.0
!                 * when IWRKSC=0, change IORI from 0 to 1 in call to
!                   NAMER
! 10.01.30 (BTD)  * eliminate restriction against 
!                   CMDFRM='TFRAME' and JPBC=0
! end history
! Copyright (C) 1993,1994,1995,1996,1997,1998,1999,2000,2003,2004,2005,
!               2006,2007,2008,2009,2010 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

!    Adjustable Parameters:

!    MXNX   = max. extent of target in x direction
!    MXNY   =                          y
!    MXNZ   =                          z
!    MXPBC  = 0 if not planning to use PBC option
!           = 1 for best memory use with PBC option
!    MXCOMP = max. number of different dielectric functions
!    MXTHET = max. number of target  rotation angles THETA
!    MXBETA = max. number of target rotation angles BETA
!    MXPHI  = max. number of target rotation angles PHI
!    MXSCA  = max. number of scattered directions
!    MXRAD  = max. number of radii
!    MXWAV  = max. number of wavelengths
!    MXWAVT = max. number of wavelengths in dielectric tables

!*** Set parameter MXPBC
      INTEGER :: MXPBC,MXPBC_SH
      PARAMETER(MXPBC=1)

      INTEGER :: MXNX,MXNY,MXNZ

! experiment with temporary storage for ixyz0 to support first call to target

      INTEGER :: MXNAT0,MXN03
      
!*** Set parameters MXCOMP,MXTHET,MXBETA,MXPHI,MXRAD,MXWAV,MXWAVT,MXSCA

      INTEGER :: MXCOMP,MXCOMP_SH
      PARAMETER(MXCOMP=9)
      INTEGER :: MXTHET,MXTHET_SH,MXBETA,MXBETA_SH,MXPHI,MXPHI_SH
      PARAMETER(MXTHET=100,MXBETA=100,MXPHI=100)
      INTEGER :: MXRAD,MXRAD_SH,MXWAV,MXWAV_SH,MXWAVT
      PARAMETER(MXRAD=1000,MXWAV=1000,MXWAVT=1500)
      INTEGER :: MXSCA,MXSCA_SH
      PARAMETER(MXSCA=10000)
      INTEGER :: MXTIMERS,NTIMERS
      PARAMETER(MXTIMERS=20,NTIMERS=12)

!*** Derived Parameters: **********************************************

      INTEGER :: MXNAT,MXN3
      INTEGER :: MXCXSC

! Local variables:

      CHARACTER :: CMDFRM*6,CMDSOL*6,CMDTRQ*6,CMSGNM*70

      COMPLEX(WP) :: CXEN

      INTEGER :: IANISO,JX,JXMAX,JXMIN,JY,JYMAX,JYMIN,JZ,JZMAX,JZMIN, &
                 NAVG,NRWORD,NCOMP_NEED

      REAL(WP) :: AEFF,AK1,AK3,BETAD,BETAMI,BETAMX,BETMID,BETMXD,COSBET,      &
                  COSPHI,COSTHE,CWORD,DAEFF,DEGRAD,DSTORAGE,ETASCA,FREQ,MB,   &
                  PHID,PHIMAX,PHIMID,PHIMIN,PHIMXD,PI,PIA2,PYDDX,PZDDX,RWORD, &
                  SINBET,SINPHI,SINTHE,STORAGE,STORAGE0,THETAD,THETMI,THETMX, &
                  THTMID,THTMXD,TOL,WAVE,XMAX,XMIN,XX,YMAX,YMIN,ZMAX,ZMIN

      !e-Beam center variable, added by NWB 3/13/12
      REAL(WP) :: CENTER(3)
      !Electron energy variable, added by NWB 7/12/12
      REAL(WP) :: ELENERGY

! Target Properties:
! AEFF=aeff=effective radius of target (physical units)
! PIA2=pi*(aeff/d)**2=pi*(3*NAT/4*pi)**(2/3)

! Incident Wave Properties:
! AK1=(k*d), where k=2*pi/wave=propagation vector and d=dipole spacing
! AK3=(k*d)**3
! AKE2=(k*d \dot E0), where E0 is incident polarization vector
! WAVE=wavelength (in vacuo) in physical units
! XX=k*aeff=2*pi*aeff/wave

! Scattering Properties:
! G=<cos(theta)> for scattered radiation
! QSCGSUM(1-3,1-2)=running (weighted) sum of g(1-3)*Q_sca
!           over target orientations, for incident polarizations 1-2,
!           where
!           g(1)=<cos(theta)> for scattered radiation
!           g(2)=<sin(theta)cos(phi)> for scattered radiation
!           g(3)=<sin(theta)sin(phi)> for scattered radiation
!           theta is measured relative to incident direction
!           phi is measured relative to x,y plane (Lab Frame)
! QPHA(1-2)=(phaseshift cross section)/(pi*aeff**2)
! QEXSUM(1-2)=running sum of QEXT over target orientations
! QABSUM(1-2)=running sum of QABS over target orientations
! QBKSUM(1-2)=running sum of QBKSCA over target orientations
! QSCSUM(1-2)=running sum of QSCA over target orientations
! QSCG2SUM(1-2)=running sum of cos^2*Q_sca
! QTRQSCSUM(1-3,1-2)=running sum of QTRQSC(1-3,1-2) over target orientat
! THETND=scattering angle theta in degrees
! PHIND=scattering angle phi in degrees
! SINPHI=sin(phi)
! DEPOLR=depolarization ratio for scattered of initially polarized
!        light into direction theta,phi
! PPOL=percentage polarization for scattering of initially unpolarized
!      light into direction theta,phi

      INTEGER :: IBETA,IBETA1,IBETH,IBETH1,IDIR,IDNC,IDVERR,IDVSHP,IERR,INIT, &
                 ILIN10,ILIN12,IOBIN,IOPAR,IORI,IORI0,IORI1,IORTH,IOSHP,IPHI, &
                 IPHI1,IRAD,IRAD0,IRAD1,ITASK,ITHETA,ITHETA1,IWAV,IWAV0,      &
                 IWAV1,IWRKSC,IWRPOL,J,JJ,JO,JPBC,LACE,LAXI,LCLM,LGI,LPI,LQI, &
                 LSC0,NBETA,NBETH,NCOMP,NORI,NPHI,NRAD,NSCAT,NSMELTS,NTHETA,  &
                 NUMPROCS,NWAV

      INTEGER :: MXITER
      INTEGER ::     &
         ITNUM(2),   &
         ITNUMMX(2), &
         SMIND1(9),  &
         SMIND2(9)

! MYID     = MPI process identifier that runs from 0-(NUMPROCS-1). 0 is
!            master process.
! NUMPROCS = number of MPI parallel processes (including the master)
! IERR     = MPI error code, which we ignore
! NAT=number of dipoles in target
! NAT3=3*NAT
! NBETA=number of different beta values for target orientation
! NBETH=NBETA*NTHETA=number of target orientations for outer orientation
! NTHETA=number of different theta values for target orientation
! NPHI=number of different phi values for target orientation
! NORI=NBETA*NTHETA*NPHI=total number of target orientations
! NSCAT=number of different scattering directions for each orientation
! NWAV=number of different wavelengths
! NRAD=number of different target radii
! NSMELTS=number of scattering matrix elements to print out
! INIT=0,1,2 for choice of initial vector |x0> in complex conjugate grad
!      method
! IORTH=1 or 2 to do 1 or 2 incident polarizations
! IWRKSC=0 or 1 to suppress or generate "wAArBBkCCC.sca" files for each
!        target orientation
! NX=x-length of extended target
! NY=y-length of extended target
! NZ=z-length of extended target

! SMIND1(1-NSMELTS) = index 1 of scattering matrix elements to
!                     be printed out (e.g., 1 for element S_{13})
! SMIND2(1-NSMELTS) = index 2 of scattering matrix elements to
!                     be printed out (e.g., 3 for element S_{13})

! Variables for fast eDDA NWB 7/11/12
      REAL(WP) :: c, h_bar, h_bar2, e_charge, m_e, velocity, &
                  DielectricConst

      CHARACTER :: CALPHA*6,CSHAPE*9,CSTAMP*26,CFLAVG*15,CFLFML*16, &
        CFLPAR*13,CFLPOL1*17,CFLPOL2*17,CFLSCA*16,CFLSHP*80,CDESCR*67,      &
        CBINFLAG*6,CBINFILE*14
      CHARACTER :: CFLLOG*14
      CHARACTER(60) :: &
        CFLEPS(MXCOMP)
!     ..
!     .. Local Arrays ..

      COMPLEX(WP),ALLOCATABLE :: &
        CXALOF(:),               &
        CXALOS(:),               &
        CXALPH(:),               &
        CXE(:),                  &
        CXSC(:),                 &
        CXSCR1(:),               &
        CXXI(:)

      COMPLEX(WP) ::          &
        CX1121(MXSCA),        &
        CX1121_1(MXSCA),      &
        CXE01(3),             &
        CXE02(3),             &
        CXE01R(3),            &
        CXE02R(3),            &
        CXEPS(MXCOMP),        &
        CXF11(MXSCA),         &
        CXF12(MXSCA),         &
        CXF21(MXSCA),         &
        CXF22(MXSCA),         &
        CXRLOC(MXCOMP+1,3,3), &
        CXRFR(MXCOMP),        &
        CXS1(MXSCA),          &
        CXS2(MXSCA),          &
        CXS3(MXSCA),          &
        CXS4(MXSCA)

! Arrays describing target properties:
! CXADIA(1-3*NAT)=diagonal elements of A matrix
! CXEPS(1-MXCOMP)=dielectric const. for each of MXCOMP materials
! CXREFR(1-3)=complex refractive index at current wavelength
! CXALPH(1-3*NAT)=diagonal elements of dipole polarizability tensor alpha
! CXALOF(1-3*NAT)=off-diagonal elements of dipole polarizability tensor
! CXALOS(1-3*NAT)=off-diagonal elements of dipole polarizability tensor

! Arrays describing incident radiation:
! CXE01(1-3)=incident polarization vector e01 prior to target rotation
! CXE02(1-3)=incident polarization vector e02 prior to target rotation
! CXE01R(1-3)=rotated incident polarization vector e01 in Target Frame
! CXE02R(1-3)=rotated incident polarization vector e02 in Target Frame
! CXE0R(1-3)=incident polarization vector (=CXE01R or CXE02R) being used
! CXE(1-3*NAT)=incident electric field vector in Target Frame

! Array describing polarization of target for one incident wave:
! CXXI(1-3*NAT)=polarization vector at each dipole

! Arrays describing scattering properties of target:
! CXFF11(1-NSCAT)=scattering matrix element f11 for each of NSCAT directions
! CXFF12(1-NSCAT)=                          f12
! CXFF21(1-NSCAT)=                          f21
! CXFF22(1-NSCAT)=                          f22
! CX1121(1-NSCAT)=sum of f11* \times f21 over orientations

! Scratch arrays:
! CXSC(1-MXCXSC)=scratch array used by DDACCG

      INTEGER*2,ALLOCATABLE :: &
        ICOMP(:),              &
        ISCR1(:)

      INTEGER,ALLOCATABLE :: &
        IXYZ0(:,:)

! IXYZ0(1-NAT0,1-3)=[x-X0(1)]/d,
!                   [y-X0(2)]/d,
!                   [z-X0(3)]/d for dipoles 1-NAT0 in real target
!                   where offset vector X0(1-3) is set by routine TARGET
!                   X0(1-3)*DX(1-3) = location of lattice site (0,0,0)

! Scratch arrays:
! ISCR1,ISCR2,ISCR3,ISCR4 (used by EXTEND)


      REAL(WP) ::           &
        A1(3),              &
        A2(3),              &
        A3(3),              &
        AEFFA(MXRAD),       &
        AKSR(3,MXSCA),      &
        BETA(MXBETA),       &
        E1A(MXWAVT),        &
        E2A(MXWAVT),        &
        EM1(3,MXSCA),       &
        EM1R(3,MXSCA),      &
        EM2(3,MXSCA),       &
        EM2R(3,MXSCA),      &
        EN0R(3),            &
        ENSC(3,MXSCA),      &
        ENSCR(3,MXSCA),     &
        ORDERM(MXSCA),      &
        ORDERN(MXSCA),      &
        PHI(MXPHI),         &
        PHIN(MXSCA)
      REAL(WP) ::           &
        QABS(2),            &
        QABSUM(2),          &
        QABSUM_1(2),        &
        QBKSCA(2),          &
        QBKSUM(2),          &
        QBKSUM_1(2),        &
        QEXT(2),            &
        QEXSUM(2),          &
        QEXSUM_1(2),        &
        QPHA(2),            &
        QPHSUM(2),          &
        QPHSUM_1(2),        &
        QSCAT(2),           &
        QSCAG(3,2),         &
        QSCAG2(2),          &
        QSCG2SUM(2),        &
        QSCG2SUM_1(2),      &
        QSCGSUM(3,2),       &
        QSCGSUM_1(3,2),     &
        QSCSUM(2),          &
        QSCSUM_1(2),        &
        QTRQAB(3,2),        &
        QTRQSC(3,2),        &
        QTRQABSUM(3,2),     &
        QTRQABSUM_1(3,2),   &
        QTRQSCSUM(3,2),     &
        QTRQSCSUM_1(3,2)
      REAL(WP) ::           &
        S1111(MXSCA),       &
        S1111_1(MXSCA),     &
        S2121(MXSCA),       &
        S2121_1(MXSCA),     &
        SHPAR(12),          &
        SM(4,4,MXSCA),      &
        SMORI(4,4,MXSCA),   &
        SMORI_1(4,4,MXSCA), &
        THETA(MXTHET),      &
        THETAN(MXSCA),      &
        TIMERS(MXTIMERS),   &
        WAVEA(MXWAV),       &
        WGTA(MXTHET,MXPHI), &
        WGTB(MXBETA),       &
        WVA(MXWAVT),        &
        X0(3),              &
        XLR(3),             &
        YLR(3),             &
        ZLR(3)

      REAL(WP),ALLOCATABLE :: &
        BETADF(:),            &
        PHIDF(:),             &
        SCRRS1(:,:),          &
        SCRRS2(:),            &
        THETADF(:)

! IANISO=0 for isotropic target
!        1 for anisotropic target, but with principal axes of
!              dielectric tensor || x_TF,y_TF,z_TF
!              (i.e., for every dipole, Dielectric Frame = Target Frame)
!        2 for general anisotropic target, where orientation of
!              Dielectric Frame is specified for every dipole
!                                  
! Arrays describing target properties:
! A1(1-3)=direction of target axis a1 in Target Frame
!        =initial direction of target axis a1 in Lab Frame
! A2(1-3)=direction of target axis a2 in Target Frame
!        =initial direction of target axis a2 in Lab Frame
! A3(1-3)=direction of target axis a3=a1xa2 in Target Frame
!        =initial direction of target axis a3 in Lab Frame
! AEFFA(1-NRAD)=values of effective radius of target (physical units)
! SHPAR(1-10)=up to 10 parameters for defining target geometry
! X0(1-3)=offset for IXYZ array, such that
!         x/d = IXYZ(J,1) - X0(1)
!         y/d = IXYZ(J,2) - X0(2)
!         z/d = IXYZ(J,3) - X0(3)
! Arrays describing orientation of "Dielectric Frame" for anisotropic
! material relative to Target Frame, for each dipole location:
! BETADF(1-NAT) = angle BETADF (radians)
! PHIDF(1-NAT) = angle PHIDF (radians)
! THETADF(1-NAT) = angle THETADF (radians)

! Arrays describing target rotation:
! THETA(1-NTHETA)=desired values of target rotation angle theta
! BETA(1-NBETA)=desired values of target rotation angle beta
! PHI(1-NPHI)=desired values of target rotation angle phi
! WGTA(1-NTHETA,1-NPHI)=weight factor for each orientation of target
!                       axis a1 in Lab Frame
! WGTB(1-NBETA)=weight factor for rotation of target around a1
! Arrays describing incident wave:
! WAVEA(1-NWAV)=values of wavelength (physical units)
! EN0R(3)=unit vector in direction of incidence in (rotated) Target
!         Frame
! AKR(3)=k*d vector in (rotated) Target Frame
! XLR(3)=unit vector xlab in target frame
! YLR(3)=unit vector ylab in target frame
! ZLR(3)=unit vector zlab in target frame
! Arrays describing scattering properties:
! ENSC(1-3,1-NSCAT)=unit scattering vectors in Lab Frame
! ENSCR(1-3,1-NSCAT)=unit scattering vectors in Target Frame
! AKSR(1-3,1-NSCAT)=scattering k-vectors in Target Frame
! EM1(1-3,1-NSCAT)=scattering polarization state 1 in Lab Frame
! EM2(1-3,1-NSCAT)=scattering polarization state 2 in Lab Frame
! EM1R(1-3,1-NSCAT)=scattering polarization state 1 in Target Frame
! EM1R(1-3,1-NSCAT)=scattering polarization state 2 in Target Frame
! THETAN(1-NSCAT)=values of theta for NSCAT scattering directions nhat
!                 at which scattering matrix fml is to be calculated
!                 For finite target (JPBC=0), THETAN is specified
!                 in the Lab Frame
! PHIN(1-NSCAT)=values of phi for scattering NSCAT directions nhat at
!               which scattering matrix fml is to be calculated
!               For finite target (JPBC=0), PHIN is specified
!               in the Lab Frame.
! S1111(1-NSCAT)=sum of |f11|**2 over target orientations
! S2121(1-NSCAT)=       |f21|**2
! QABS(1-2)=(absorption cross section)/(pi*aeff**2) for 2 polarizations
! QBKSCA(1-2)=4*pi*(diff.scatt.cross section for theta=pi)/(pi*aeff**2)
!             for 2 polarizations
! QEXT(1-2)=(extinction cross section)/(pi*aeff**2)
! QSCA(1-2)=(scattering cross section)/(pi*aeff**2)=QEXT-QABS
! QSCAT(1-2)=(scattering cross section)/(pi*aeff**2) estimated by SCAT

! Scratch array:
! SCRRS1(1-NAT,3)=scratch array used by DDACCG and SCAT
! SCRRS2(1-NAT)=scratch array used by SCAT

! External subroutines:

      EXTERNAL DIELEC,NAMER,REAPAR

! Intrinsic functions:

      INTRINSIC ATAN,SQRT,REAL

      DATA CFLPAR/'ddscat.par'/
      DATA CFLSHP/'shape.dat'/

! Set device numbers:
! IDVERR = device number for error messages
! IOPAR = device number for reading ddscat.par file
! IOSHP = device number for writing target.out (IOSHP<0 to suppress)
! IDVSHP = device number for reading shape file

      DATA IDVERR/0/
      DATA IDVSHP/10/
      DATA IOPAR/1/
      DATA IOSHP/ -1/

! IDVOUT = device number for running I/O
!          on Solaris IDVOUT=0 produces unbuffered output to standard
!                              output
!                     IDVOUT=6 produces buffered output to standard
!                              output
! Initialize IDVOUT (cannot be done in DATA statement because IDVOUT
! is in COMMON/M6/

      IDVOUT=0
!*** diagnostic
!      write(0,*)'ddscat ckpt 0.9'
!***
!***********************************************************************

! Define variables for eDDA NWB 7/11/12
      c = 3.E8_WP                     !Speed of light, m/s
      e_charge = 1.60217646E-19_WP    !Charge of electron, Coulombs
!      e_charge = 1._WP                !Charge of electron, fundamental
      m_e = 9.10938188E-31_WP         !Mass of electron, kg
!      velocity = 0.5_WP * c           !Velocity of electron, m/s
      DielectricConst = 1._WP         !Dielectric constant -- 1 for vacuum
      h_bar = 1.054572E-34_WP         !h_bar, j-s
      h_bar2 = 6.582119E-19_WP        !h_bar, eV-s

!*********************************************************************** 
! MPI environment setup:

      CALL MPI_INIT(IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      IF(MYID==0)THEN
         WRITE(CMSGNM,FMT='(A,I4)')'NUMPROC=',NUMPROCS
         CALL WRIMSG('DDSCAT',CMSGNM)
      ENDIF

! MPI_INIT is called by each process to initialize the MPI environment
! MPI_COMM_RANK is called by each process to determine its identifier
!                which is stored in MYID
! MPI_COMM_SIZE is called by each process to determine the number of
!                parallel processes, stored in NUMPROCS
! If not using MPI, these calls can be turned into dummies by including
! subroutines in mpi_fake.f (also need to comment out INCLUDE 'mpif.h'
! and uncomment declaration of integer variable MPI_COMM_WORLD as noted
! above)

! Initialize variables to pass to SHARE1:

      MXBETA_SH=MXBETA
      MXCOMP_SH=MXCOMP
      MXPBC_SH=MXPBC
      MXPHI_SH=MXPHI
      MXRAD_SH=MXRAD
      MXSCA_SH=MXSCA
      MXTHET_SH=MXTHET
      MXWAV_SH=MXWAV

! Call subroutine NAMID to obtain unique filename for output file
! written by this process.

!*** diagnostic
!      write(0,*)'ddscat ckpt 1.1, myid=',myid
!      call mpi_barrier(mpi_comm_world,ierr)
!***      

      CALL NAMID(MYID,CFLLOG)
!*** diagnostic
!      write(0,*)'ddscat ckpt 1.2, myid=',myid
!***
!=======================================================================

! When following statement is enabled, output is buffered.
! If using MPI, each MPI process will write to a unique file named
! ddscat.log_xxx , where xxx=000,001,002, etc.
! If it is necessary to have unbuffered output for debugging purposes,
! comment out the following OPEN statement, in which case running
! diagnostic information will be written to standard output.
! Note that in this case, if using MPI, all MPI processes will write
! to a single output stream.

!      OPEN(UNIT=IDVOUT,FILE=CFLLOG)

!=======================================================================
      CALL VERSION(CSTAMP)

      IF(MYID==0)THEN   ! begin if(myid==0)... #1
         WRITE(IDVOUT,FMT=9000)CSTAMP

#ifdef openmp
         WRITE(IDVOUT,FMT=9000)'compiled with OpenMP enabled'
#endif

         IF(WP==KIND(0.E0))THEN
            WRITE(CMSGNM,FMT='(A)')'    Single-precision version'

! for storage computations:
! NRWORD,RWORD = length (bytes) of REAL word
! CWORD = length (bytes) of COMPLEX word

            NRWORD=4
            RWORD=4._WP
            CWORD=2._WP*RWORD
         ELSEIF(WP==KIND(0.D0))THEN
            WRITE(CMSGNM,FMT='(A)')'    Double-precision version'
            NRWORD=8.
            RWORD=8._WP
            CWORD=2._WP*RWORD
         ELSE
            CALL ERRMSG('FATAL','DDSCAT', &
                         'Fatal error determining word length in DDSCAT')
         ENDIF
         CALL WRIMSG('DDSCAT',CMSGNM)
!*** diagnostic
!         write(0,*)'ddscat ckpt 1.3, myid=',myid
!***
      ENDIF   !--- endif(myid==0)... #1

!*** diagnostic
!      write(0,*)'ddscat ckpt 1.4, myid=',myid
!***
      PI=4._WP*ATAN(1._WP)
      DEGRAD=180._WP/PI

!**** Read program control parameters
! CMDFFT=choice of solution methods
! CMDFRM=specify scattering directions in either Lab Frame ('LFRAME')
!        or target frame ('TFRAME')
! CSHAPE=description of shape
! SHPAR(1-10)= up to 10 parameters for defining target geometry
! INIT=selection of |x0> in CCG method
! TOL=error tolerance in CCG method
! MXWAV=dimensioning information for array WAVEA
! NWAV=number of wavelengths to be used
! WAVEA(1-NWAV)=wavelengths (physical units)
! MXRAD=dimensioning information for array AEFFA
! NRAD=number of radii to be used
! AEFFA(1-NRAD)=effective radii (physical units)
! IORTH=1 or 2 for 1 or 2 incident polarization states to be used
! NBETA=number of beta values to be used for target orientation
! BETAMI=minimum value of beta (radians)
! BETAMX=maximum value of beta (radians)
! NTHETA=number of theta values to be used for target orientation
! THETMI=minimum value of theta (radians)
! THETMX=maximum value of theta (radians)
! NPHI=number of phi values to be used for target orientation
! PHIMIN=minimum value of phi (radians)
! PHIMAX=maximum value of phi (radians)
! DEGRAD=degrees/radian (conversion factor)
! MXSCA=dimensioning information for arrays PHIN, THETAN
! NSCAT=number of scattered directions for computation of fml
! THETAN(1-NSCAT)=theta values (radians) for scattered directions
!                 for finite targets (JPBC=0), THETAN is in Lab Frame
! PHIN(1-NSCAT)=phi values (radians) for scattered directions
!                 for finite targets (JPBC=0), PHIN is in Lab Frame
! for JPBC=1 or 2 (infinite target with periodicity in 1-d):
!                 PHIN(1-NSCAT) specifies scattering order
!                 THETAN(1-NSCAT) specifies azimuthal angle of scatterin
!                 direction around periodicity axis
! for JPBC=3 (infinite target with periodicity in y_TF and z_TF)
!                 PHIN(1-NSCAT) specifies scattering order M
!                 THETAN(1-NSCAT) specifies scattering order N

      IF(MYID==0)THEN   ! begin if(myid==0)... #2

!*** diagnostic
!         write(0,*)'ddscat ckpt 1.5, myid=',myid
!***
         CALL REAPAR(CBINFLAG,AEFFA,BETAMI,BETAMX,CALPHA,CFLEPS,CFLPAR,     &
                     CFLSHP,CMDSOL,CMDFFT,CMDFRM,CSHAPE,CMDTRQ,CXE01,CXE02, &
                     DEGRAD,ETASCA,GAMMA,IDVERR,IDVOUT,INIT,IOPAR,IORTH,    &
                     JPBC,IWRKSC,IWRPOL,MXNX,MXNY,MXNZ,MXBETA,MXCOMP,MXPHI, &
                     MXRAD,MXSCA,MXTHET,MXWAV,NBETA,NCOMP,NPHI,NRAD,NSCAT,  &
                     NTHETA,NWAV,IWAV0,IRAD0,IORI0,NSMELTS,PHIMIN,PHIMAX,   &
                     PHIN,THETAN,THETMI,THETMX,TOL,WAVEA,SHPAR,SMIND1,      &
                     SMIND2,DX,CENTER,ELENERGY) !ARGUMENTS'CENTER'AND 'ELENERGY' ADDED BY NWB 3/13/12, 7/12/12
!*** diagnostic
!         write(0,*)'ddscat ckpt 1.6, myid=',myid
!         write(0,*)'nscat=',nscat
!         write(0,*)'thetan(1)=',thetan(1)
!         write(0,*)'thetan(2)=',thetan(2)
!***

!***********************************************************************

! Convert electron energy from keV to electron velocity in m/s
! NWB 7/12/12

         ELENERGY = ELENERGY*1.60217646E-16_WP !Convert to Joules
         velocity = c * (1._WP - ((ELENERGY / (m_e * c**2._WP)) + 1 ) &
                    **(-2._WP))**0.5_WP

!***********************************************************************
! Allocate memory required by TARGET and EXTEND

         MXNAT0=MXNX*MXNY*MXNZ
         MXN03=3*MXNAT0
         MB=REAL(1024**2)

! STORAGE0 = estimate of storage (bytes) that is independent of target size
        
         IF(KIND(STORAGE)==KIND(0.E0))THEN
!***diagnostic
!            write(0,*)'ddscat ckpt 1.61: single-precision'
!***
            STORAGE0=37.6E6/MB     !35.0E6 MB 
         ELSE
!***diagnostic
!            write(0,*)'ddscat ckpt 1.62: double-precision'
!***
            STORAGE0=41.9E6/MB    !40.0D6 MB
         ENDIF

!*** diagnostic
!         write(0,*)'ddscat ckpt 1.7'
!***
         STORAGE=STORAGE0
         WRITE(CMSGNM,FMT='(A)')'Allocate memory for first call to TARGET...'
         CALL WRIMSG('DDSCAT',CMSGNM)

! allocate IXYZ0

         DSTORAGE=REAL(2*MXNAT0*3)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
            ' MB for IXYZ0;  total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
         ALLOCATE(IXYZ0(MXNAT0,3))
!-----------------------------------------------------------------------
! following introduced by Art Lazanoff to touch the array CXZW
! which supposedly helps with allocation and distribution of this memory
! to the different processors...
! BTD: I can't really see why the following is beneficial...

!$omp parallel do
         DO J=1,MXNAT0
            IXYZ0(J,1)=0
            IXYZ0(J,2)=0
            IXYZ0(J,3)=0
         ENDDO
!$omp end parallel do
!-----------------------------------------------------------------------


! allocate ICOMP

         DSTORAGE=REAL(2*MXN03)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
            ' MB for ICOMP;  total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
         ALLOCATE(ICOMP(MXN03))

! allocate IOCC

         DSTORAGE=REAL(2*MXNAT0)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
            ' MB for IOCC;   total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
         ALLOCATE(IOCC(MXNAT0))

! allocate ISCR1

         DSTORAGE=REAL(2*MXN03)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
            ' MB for ISCR1;  total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
         ALLOCATE(ISCR1(MXN03))

! allocate BETADF

         DSTORAGE=RWORD*REAL(MXNAT0)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
            ' MB for BETADF; total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
         ALLOCATE(BETADF(MXNAT0))

! allocate PHIDF

         DSTORAGE=RWORD*REAL(MXNAT0)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
            ' MB for PHIDF;  total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
         ALLOCATE(PHIDF(MXNAT0))

! allocate THETADF

         DSTORAGE=RWORD*REAL(MXNAT0)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
            ' MB for THETADF;total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
         ALLOCATE(THETADF(MXNAT0))

! allocate SCRRS2

         DSTORAGE=RWORD*REAL(MXNAT0)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
            ' MB for SCRRS2; total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
         ALLOCATE(SCRRS2(MXNAT0))

!*** diagnostic
!        write(0,*)'completed preliminary allocation of storage'
!        write(0,*)'about to call target...'
!***
!-----------------------------------------------------------------------------
! Specify target properties

! A1(1-3)=target axis 1 in Target Frame
! A2(1-3)=target axis 2 in Target Frame

         CALL TARGET(A1,A2,CSHAPE,CFLSHP,IDVSHP,IOSHP,CDESCR,MXNAT0,SHPAR,DX, &
                     NAT0,IXYZ0,ICOMP,IDVOUT,MXN03,NAT3,PYD,PZD,BETADF,PHIDF, &
                     THETADF,X0,IANISO,NCOMP_NEED)
!*** diagnostic
!         write(0,*)'ddscat ckpt 1.75'
!         write(0,*)'     pyd=',pyd
!         write(0,*)'     pzd=',pzd
!***
         

! check to see if NCOMP_NEED is compatible with NCOMP

         IF(NCOMP_NEED.GT.NCOMP)THEN
            WRITE(CMSGNM,'(A,I4,A,I4,A)')'Error:NCOMP_NEED=',NCOMP_NEED,' >', &
                                       NCOMP,' = NCOMP'
            CALL WRIMSG('DDSCAT',CMSGNM)
            WRITE(CMSGNM,'(A,A)')'ddscat.par does not list enough ',  &
                                 'different refractive index filenames'
            CALL ERRMSG('FATAL','DDSCAT','NCOMP is too small')
         ENDIF


! write out some information on target geometry

         JXMIN=IXYZ0(1,1)
         JYMIN=IXYZ0(1,2)
         JZMIN=IXYZ0(1,3)
         JXMAX=JXMIN
         JYMAX=JYMIN
         JZMAX=JZMIN
         DO J=2,NAT0
            JX=IXYZ0(J,1)
            JY=IXYZ0(J,2)
            JZ=IXYZ0(J,3)
            IF(JX<JXMIN)JXMIN=JX
            IF(JX>JXMAX)JXMAX=JX
            IF(JY<JYMIN)JYMIN=JY
            IF(JY>JYMAX)JYMAX=JY
            IF(JZ<JZMIN)JZMIN=JZ
            IF(JZ>JZMAX)JZMAX=JZ
         ENDDO
         XMIN=REAL(JXMIN)-0.5_WP+X0(1)
         XMAX=REAL(JXMAX)+0.5_WP+X0(1)
         YMIN=REAL(JYMIN)-0.5_WP+X0(2)
         YMAX=REAL(JYMAX)+0.5_WP+X0(2)
         ZMIN=REAL(JZMIN)-0.5_WP+X0(3)
         ZMAX=REAL(JZMAX)+0.5_WP+X0(3)
         WRITE(CMSGNM,FMT='(A)')'returned from TARGET:'
         CALL WRIMSG('DDSCAT',CMSGNM)
         WRITE(CMSGNM,FMT='(2I6,A)')JXMIN,JXMAX,' = min, max values of JX'
         CALL WRIMSG('DDSCAT',CMSGNM)
         WRITE(CMSGNM,FMT='(2I6,A)')JYMIN,JYMAX,' = min, max values of JY'
         CALL WRIMSG('DDSCAT',CMSGNM)
         WRITE(CMSGNM,FMT='(2I6,A)')JZMIN,JZMAX,' = min, max values of JZ'
         CALL WRIMSG('DDSCAT',CMSGNM)
         WRITE(CMSGNM,FMT='(3F10.5,A)')X0,' = X0(1-3)'
         CALL WRIMSG('DDSCAT',CMSGNM)
         WRITE(CMSGNM,FMT='(A)')                                 &
               ' --- physical extent of target volume (occupied sites) ---'
         CALL WRIMSG('DDSCAT',CMSGNM)
         WRITE(CMSGNM,FMT='(2F12.5,A)')XMIN,XMAX,' = min, max values of x/d'
         CALL WRIMSG('DDSCAT',CMSGNM)
         WRITE(CMSGNM,FMT='(2F12.5,A)')YMIN,YMAX,' = min, max values of y/d'
         CALL WRIMSG('DDSCAT',CMSGNM)
         WRITE(CMSGNM,FMT='(2F12.5,A)')ZMIN,ZMAX,' = min, max values of z/d'
         CALL WRIMSG('DDSCAT',CMSGNM)

! Return from TARGET with following arrays defined:
!      NAT0   = number of occupied lattice sites
! IXYZ0(J,1-3)= (IX,IY,IZ) for occupied sites J=1-NAT0
!     X0(1-3) = specifies location in TF of lattice site (IX,IY,IZ)=(0,0,0) :
!       x_TF  = (IX+X0(1))*d*DX(1)
!       y_TF  = (IY+X0(2))*d*DX(2)
!       z_TF  = (IZ+X0(3))*d*DX(3)
!     A1(1-3) = target axis A1 in TF
!     A2(1-3) = target axis A2 in TF
!       PYD   = periodicity/[d*DX(2)] in y_TF direction (for JPBC=1 or 3)
!       PZD   = periodicity/[d*DX(3)] in z_TF direction (for JPBC=2 or 3)
! ICOMP(J,1-3)= composition at lattice site J, for directions 1-3
!               in local "Dielectric Frame" where dielectric tensor is
!               diagonalized
!
! for targets composed of general anisotropic materials (IANISO=2), also define:
!   BETADF(J) = rotation angle beta defining orientation in TF of
!               "dielectric frame" (DF) for site J
!    PHIDF(J) = rot. angle phi defining orientation in TF of DF for site J
!  THETADF(J) = rot. angle theta defining orientation in TF of DF for site J

! Calculate target axis A3 in Target Frame

         A3(1)=A1(2)*A2(3)-A1(3)*A2(2)
         A3(2)=A1(3)*A2(1)-A1(1)*A2(3)
         A3(3)=A1(1)*A2(2)-A1(2)*A2(1)

         PYDDX=PYD*DX(2)
         PZDDX=PZD*DX(3)

! Calculate DAEFF = d/aeff for this target (d=lattice spacing)

         DAEFF=(4._WP*PI/(3._WP*REAL(NAT0,KIND=WP)))**(1._WP/3._WP)

! NAT0=number of dipoles in "real" target
! Extend target to rectangular volume suitable for FFT

         CALL EXTEND(CMDFFT,ICOMP,ISCR1,IDVERR,IDVOUT,IOCC,IXYZ0(1,1),     &
                     IXYZ0(1,2),IXYZ0(1,3),BETADF,PHIDF,THETADF,SCRRS2,X0, &
                     MXNX,MXNY,MXNZ,MXNAT0,MXN03,NAT,NAT0,NAT3,NX,NY,NZ)

! write out some information on target geometry
! following call to EXTEND, as a check

!*** diagnostic
!         write(0,*)'ddscat ckpt 1.8, myid=',myid
!***
         JXMIN=IXYZ0(1,1)
         JYMIN=IXYZ0(1,2)
         JZMIN=IXYZ0(1,3)
         JXMAX=JXMIN
         JYMAX=JYMIN
         JZMAX=JZMIN
         DO J=2,NAT0
            JX=IXYZ0(J,1)
            JY=IXYZ0(J,2)
            JZ=IXYZ0(J,3)
            IF(JX<JXMIN)JXMIN=JX
            IF(JX>JXMAX)JXMAX=JX
            IF(JY<JYMIN)JYMIN=JY
            IF(JY>JYMAX)JYMAX=JY
            IF(JZ<JZMIN)JZMIN=JZ
            IF(JZ>JZMAX)JZMAX=JZ
         ENDDO
         XMIN=REAL(JXMIN)-0.5_WP+X0(1)
         XMAX=REAL(JXMAX)+0.5_WP+X0(1)
         YMIN=REAL(JYMIN)-0.5_WP+X0(2)
         YMAX=REAL(JYMAX)+0.5_WP+X0(2)
         ZMIN=REAL(JZMIN)-0.5_WP+X0(3)
         ZMAX=REAL(JZMAX)+0.5_WP+X0(3)
         WRITE(CMSGNM,FMT='(A)')'returned from EXTEND:'
         CALL WRIMSG('DDSCAT',CMSGNM)
         WRITE(CMSGNM,FMT='(2I6,A)')JXMIN,JXMAX,' = min, max values of JX'
         CALL WRIMSG('DDSCAT',CMSGNM)
         WRITE(CMSGNM,FMT='(2I6,A)')JYMIN,JYMAX,' = min, max values of JY'
         CALL WRIMSG('DDSCAT',CMSGNM)
         WRITE(CMSGNM,FMT='(2I6,A)')JZMIN,JZMAX,' = min, max values of JZ'
         CALL WRIMSG('DDSCAT',CMSGNM)
         WRITE(CMSGNM,FMT='(A)')                               &
               ' ------- physical extent of target volume -------'
         CALL WRIMSG('DDSCAT',CMSGNM)
         WRITE(CMSGNM,FMT='(2F12.5,A)')XMIN,XMAX,' = min, max values of x/d'
         CALL WRIMSG('DDSCAT',CMSGNM)
         WRITE(CMSGNM,FMT='(2F12.5,A)')YMIN,YMAX,' = min, max values of y/d'
         CALL WRIMSG('DDSCAT',CMSGNM)
         WRITE(CMSGNM,FMT='(2F12.5,A)')ZMIN,ZMAX,' = min, max values of z/d'
         CALL WRIMSG('DDSCAT',CMSGNM)
!***

! EXTEND reorders the occupied lattice sites,
!    returns the extents NX,NY,NZ for the extended target.

!    if necessary, shifts JX,JY,JZ for each site to that JX,JY,JZ
!    run from 1-NX,1-NY,1-NZ

!    IXYZ0(JA,1-3) = coordinates JX,JY,JZ of occupied lattice sites JA=1-NAT0

! location JX,JY,JZ in extended target corresponds to 
! JE=NX*NY*(JZ-1)+NX*(JY-1)+JX, where JE runs from 1 to NX*NY*NZ

!    ICOMP(JE,1-3)  = composition for each site in extended target
!                   = 0 if site is unoccupied
!       IOCC(JE)    = 1 if occupied, 0 if not
! BETADF(JA=1-NAT0) = rotation angle beta for DF in reordered physical target
!                    for JA=1-NAT0
!  PHIDF(JA=1-NAT0) = rotation angle phi for DF in reordered physical target
! THETADF(JA=1-NAT0)= rotation angle theta for DF in reordered physical target
!    X0(1-3)        = location in TF of lattice site JX,JY,JZ=0,0,0

! With actual size of extended target determined, 
! now set required values of MXNX,MXNY,MXNZ,MXNAT,MXN3 and reallocate arrays
! IXYZ0,ICOMP,IOCC,ISCR1,BETADF,PHIDF,THETADF,SCRRS2

         MXNX=NX
         MXNY=NY
         MXNZ=NZ
         MXNAT=NX*NY*NZ
         MXN3=3*MXNAT

!*** diagnostic
!         write(0,*)'ddscat ckpt 1.85, myid=',myid,' mxn3=',mxn3
!***


! MXCXSC specifies required complex scratch array CXSC
! increase from 7*MXN3 required by PETR
!           to 10*MXN3 required by PIMCBICGSTAB method
!                 (PIM Complex BiConjugate Gradient with Stabilization)

         MXCXSC=10*MXN3

! Initialize various variables in common blocks which are derived
! from parameters:

         MXNATF=MXNAT
         MXNXF=MXNX
         MXNYF=MXNY
         MXNZF=MXNZ
         MXN3F=MXN3

! Pseudo names of scratch array  positions:

         LACE=1
         LAXI=1*MXN3+1
         LCLM=2*MXN3+1
         LGI=3*MXN3+1
         LPI=4*MXN3+1
         LQI=5*MXN3+1
         LSC0=6*MXN3+1

!*** diagnostic
!         write(0,*)'ddscat ckpt 1.9, myid=',myid
!***
! Reallocate ISCR1:

         WRITE(CMSGNM,FMT='(A)')'Allocate memory required for this problem...'
         CALL WRIMSG('DDSCAT',CMSGNM)

         STORAGE=STORAGE0

         DSTORAGE=REAL(2*MXN3)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
            ' MB for ISCR1;  total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
         DEALLOCATE(ISCR1)
         ALLOCATE(ISCR1(MXN3))

      ENDIF   !--- endif(myid==)... #2

! need to share dimensioning parameters:
!    MXNAT
!    NAT0
!    MXN3
!    MXNX,MXNY,MXNZ
!    MXCXSC
!    MXPBC

!*** diagnostic
!      write(0,*)'ddscat ckpt 1.94, myid=',myid,' mxn3=',mxn3
!***

! this mpi_barrier is required to ensure that LACE, etc have been
! evaluated by MYID=0 before being SHAREd by other threads

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

!*** diagnostic
!      write(0,*)'ddscat ckpt 1.95, myid=',myid
!***      

      CALL SHARE0(LACE,LAXI,LCLM,LGI,LPI,LQI,                        &
                  MXN3,MXNAT,MXNX,MXNY,MXNZ,MXPBC_SH,MXCXSC,MYID,NAT0)

!*** diagnostic
!      write(0,*)'ddscat ckpt 1.96, myid=',myid,' mxn3=',mxn3
!      call mpi_barrier(mpi_comm_world,ierr)
!***      

! Reallocate SCRRS2:

      IF(MYID==0)THEN   !--- begin if(myid==0)... #3
         DSTORAGE=RWORD*REAL(MXNAT)/MB
         STORAGE=STORAGE+DSTORAGE
            WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
               ' MB for SCRRS2; total=',STORAGE,' MB'
            CALL WRIMSG('DDSCAT',CMSGNM)
!*** diagnostic
!      write(0,*)'ddscat ckpt 2, myid=',myid
!***
         DEALLOCATE(SCRRS2)
!*** diagnostic
!      write(0,*)'ddscat ckpt 3, myid=',myid
!***
      ENDIF   !--- end if(myid==0)... #3

!*** diagnostic
!      call mpi_barrier(mpi_comm_world,ierr)
!      write(0,*)'ddscat ckpt 4, myid=',myid
!***
      ALLOCATE(SCRRS2(MXNAT))
!*** diagnostic
!      write(0,*)'ddscat ckpt 5, myid=',myid
!      call mpi_barrier(mpi_comm_world,ierr)
!***

! Reallocate IXYZ0:

      IF(MYID==0)THEN   !--- begin if(myid==0)... #4
         DSTORAGE=REAL(2*NAT0*3)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
            ' MB for IXYZ0;  total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
         DO J=1,NAT0
            ISCR1(J)=IXYZ0(J,1)
            ISCR1(J+NAT0)=IXYZ0(J,2)
            ISCR1(J+2*NAT0)=IXYZ0(J,3)
         ENDDO
!*** diagnostic
!      write(0,*)'ddscat ckpt 6, myid=',myid
!***
         DEALLOCATE(IXYZ0)
!*** diagnostic
!      write(0,*)'ddscat ckpt 7, myid=',myid
!***
      ENDIF   !--- end if(myid==0)... #4
!*** diagnostic
!      call mpi_barrier(mpi_comm_world,ierr)
!      write(0,*)'ddscat ckpt 8, myid=',myid
!***
      ALLOCATE(IXYZ0(NAT0,3))
!*** diagnostic
!      write(0,*)'ddscat ckpt 9, myid=',myid
!***
      IF(MYID==0)THEN   !--- begin if(myid==0)... #5
!$omp parallel do
         DO J=1,NAT0
            DO JJ=1,3
               IXYZ0(J,JJ)=ISCR1((JJ-1)*NAT0+J)
            ENDDO
         ENDDO
!$omp end parallel do

! Note: IXYZ0(1-NAT0,1-3) contain coordinates of occupied lattice sites
!       [reordered to standard order used for arrays produced by REDUCE]
!       BETADF(1-NAT0),PHIDF(1-NAT0),THETADF(1-NAT0) contain orientation
!                         angles for anisotropic material at occupied
!                         lattice sites (not used for isotropic material

! Reallocate ICOMP:
!*** diagnostic
!         write(0,*)'ddscat ckpt 9.5, myid=',myid
!***
         DSTORAGE=REAL(2*MXN3)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
            ' MB for ICOMP;  total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
         DO J=1,NAT3
            ISCR1(J)=ICOMP(J)
         ENDDO
!*** diagnostic
!      write(0,*)'ddscat ckpt 10, myid=',myid
!***
         DEALLOCATE(ICOMP)
!*** diagnostic
!      write(0,*)'ddscat ckpt 11, myid=',myid
!***
      ENDIF   !--- end if(myid==0)... #5
!*** diagnostic
!      write(0,*)'ddscat ckpt 12, myid=',myid
!      call mpi_barrier(mpi_comm_world,ierr)
!***
      ALLOCATE(ICOMP(MXN3))
!*** diagnostic
!      write(0,*)'ddscat ckpt 13, myid=',myid
!***
      IF(MYID==0)THEN   !--- begin if(myid==0)... #6
         DO J=1,NAT3
            ICOMP(J)=ISCR1(J)
         ENDDO

! Reallocate IOCC:

         DSTORAGE=REAL(2*MXNAT)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
            ' MB for IOCC;   total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
         DO J=1,NAT
            ISCR1(J)=IOCC(J)
         ENDDO 
!*** diagnostic
!      write(0,*)'ddscat ckpt 14, myid=',myid
!***
         DEALLOCATE(IOCC)
!*** diagnostic
!      write(0,*)'ddscat ckpt 15, myid=',myid
!***
      ENDIF   !--- end if(myid==0)... #6
!*** diagnostic
!      write(0,*)'ddscat ckpt 16, myid=',myid
!***
      ALLOCATE(IOCC(MXNAT))
!*** diagnostic
!      write(0,*)'ddscat ckpt 17, myid=',myid
!***
      IF(MYID==0)THEN   !--- begin if(myid==0)... #7
         DO J=1,NAT
            IOCC(J)=ISCR1(J)
         ENDDO

! Reallocate BETADF:

         DSTORAGE=RWORD*REAL(NAT0)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
            ' MB for BETADF; total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
         DO J=1,NAT0
            SCRRS2(J)=BETADF(J)
         ENDDO
!*** diagnostic
!      write(0,*)'ddscat ckpt 18, myid=',myid
!***
         DEALLOCATE(BETADF)
!*** diagnostic
!      write(0,*)'ddscat ckpt 19, myid=',myid
!***
      ENDIF   !--- end if(myid==0)... #7
!*** diagnostic
!      write(0,*)'ddscat ckpt 20, myid=',myid
!***
      ALLOCATE(BETADF(MXNAT))
!*** diagnostic
!      write(0,*)'ddscat ckpt 21, myid=',myid
!***
      IF(MYID==0)THEN   !--- begin if(myid==0)... #8
         DO J=1,NAT0
            BETADF(J)=SCRRS2(J)
         ENDDO

! Reallocate PHIDF:

         DSTORAGE=RWORD*REAL(NAT0)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
            ' MB for PHIDF;  total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
         DO J=1,NAT0
            SCRRS2(J)=PHIDF(J)
         ENDDO
!*** diagnostic
!      write(0,*)'ddscat ckpt 22, myid=',myid
!***
         DEALLOCATE(PHIDF)
!*** diagnostic
!      write(0,*)'ddscat ckpt 23, myid=',myid
!***
      ENDIF   !--- end if(myid==0)... #8
!*** diagnostic
!      write(0,*)'ddscat ckpt 24, myid=',myid
!***
      ALLOCATE(PHIDF(MXNAT))
!*** diagnostic
!      write(0,*)'ddscat ckpt 25, myid=',myid
!***
      IF(MYID==0)THEN   !--- begin if(myid==0)... #9
         DO J=1,NAT0
            PHIDF(J)=SCRRS2(J)
         ENDDO

! Reallocate THETADF:

         DSTORAGE=RWORD*REAL(NAT0)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
            ' MB for THETADF;total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
         DO J=1,NAT0
            SCRRS2(J)=THETADF(J)
         ENDDO
!*** diagnostic
!      write(0,*)'ddscat ckpt 26, myid=',myid
!***
         DEALLOCATE(THETADF)
!*** diagnostic
!      write(0,*)'ddscat ckpt 27, myid=',myid
!***
      ENDIF   !--- end if(myid==0)... #9
!*** diagnostic
!      write(0,*)'ddscat ckpt 28, myid=',myid
!***
      ALLOCATE(THETADF(MXNAT))
!*** diagnostic
!      write(0,*)'ddscat ckpt 29, myid=',myid
!***
      IF(MYID==0)THEN   !--- begin if(myid==0)... #10
         DO J=1,NAT0
            THETADF(J)=SCRRS2(J)
         ENDDO

!------------------------------------------------------------------------------

! allocate CXADIA

         DSTORAGE=CWORD*REAL(MXNAT*3)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
            ' MB for CXADIA; total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
      ENDIF   !--- end if(myid==0)... #10
!*** diagnostic
!      write(0,*)'ddscat ckpt 30, myid=',myid
!***
      ALLOCATE(CXADIA(MXNAT,3))
!*** diagnostic
!      write(0,*)'ddscat ckpt 31, myid=',myid
!***
!---------------------------------------------------------

! allocate CXZC
      IF(MYID==0)THEN   !--- begin if(myid==0)... #11
         DSTORAGE=CWORD*REAL((MXNX+1+MXPBC*(MXNX-1))*   &
                             (MXNY+1+MXPBC*(MXNY-1))*   &
                             (MXNZ+1+MXPBC*(MXNZ-1))*6)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
            ' MB for CXZC;   total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
      ENDIF   !--- end if(myid==0)... #11
!*** diagnostic
!      write(0,*)'ddscat ckpt 32, myid=',myid
!***
      ALLOCATE(CXZC(MXNX+1+MXPBC*(MXNX-1), &
                    MXNY+1+MXPBC*(MXNY-1), &
                    MXNZ+1+MXPBC*(MXNZ-1),6))

!-----------------------------------------------------------------------
! following introduced by Art Lazanoff to touch the array CXZW
! which supposedly helps with allocation and distribution of this memory
! to the different processors...
! BTD: I can't really see why the following is beneficial...

      DO JJ=1,6
!$omp parallel do
         DO JZ=1,MXNZ+1+MXPBC*(MXNZ-1)
            DO JY=1,MXNY+1+MXPBC*(MXNY-1)
               DO JX=1,MXNX+1+MXPBC*(MXNX-1)
                  CXZC(JX,JY,JZ,JJ)=(0._WP,0._WP)
               ENDDO
            ENDDO
         ENDDO
!$omp end parallel do
      ENDDO
!-----------------------------------------------------------------------

!*** diagnostic
!      write(0,*)'ddscat ckpt 33, myid=',myid
!***
! allocate CXZW

      IF(MYID==0)THEN   !--- begin if(myid==0)... #12
         DSTORAGE=CWORD*REAL((2*MXNX)*(2*MXNY)*(2*MXNZ)*3)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
                          ' MB for CXZW;   total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
      ENDIF   !--- end if(myid==0)... #12
!*** diagnostic
!      write(0,*)'ddscat ckpt 34, myid=',myid
!***
      ALLOCATE(CXZW(2*MXNX,2*MXNY,2*MXNZ,3))

!-----------------------------------------------------------------------
! following introduced by Art Lazanoff to touch the array CXZW
! which supposedly helps with allocation and distribution of this memory
! to the different processors...
! BTD: I can't really see why the following is beneficial...

      DO JJ=1,3
!$omp parallel do
         DO JZ=1,2*MXNZ
            DO JY=1,2*MXNY
               DO JX=1,2*MXNX 
                  CXZW(JX,JY,JZ,JJ)=(0._WP,0._WP)
               ENDDO
            ENDDO
         ENDDO
!$omp end parallel do
      ENDDO
!-----------------------------------------------------------------------

!*** diagnostic
!      write(0,*)'ddscat ckpt 35, myid=',myid
!***

! allocate CXAOFF

      IF(MYID==0)THEN   !--- begin if(myid==0)... #13
         DSTORAGE=CWORD*REAL(3*MXNAT)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
                          ' MB for CXAOFF; total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
      ENDIF   !--- end if(myid==0)... #13
!*** diagnostic
!      write(0,*)'ddscat ckpt 36, myid=',myid
!***
      ALLOCATE(CXAOFF(MXNAT,3))
!*** diagnostic
!      write(0,*)'ddscat ckpt 37, myid=',myid
!***

! allocate CXALOF

      IF(MYID==0)THEN   !--- begin if(myid==0)... #14
         DSTORAGE=CWORD*REAL(MXN3)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
                          ' MB for CXALOF; total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
      ENDIF   !--- end if(myid==0)... #14
!*** diagnostic
!      write(0,*)'ddscat ckpt 38, myid=',myid
!***
      ALLOCATE(CXALOF(MXN3))
!*** diagnostic
!      write(0,*)'ddscat ckpt 39, myid=',myid
!***

! allocate CXALOS

      IF(MYID==0)THEN   !--- begin if(myid==0)... #15
         DSTORAGE=CWORD*REAL(MXN3)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
                          ' MB for CXALOS; total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
      ENDIF   !--- end if(myid==0)... #15
!*** diagnostic
!      write(0,*)'ddscat ckpt 40, myid=',myid
!***
      ALLOCATE(CXALOS(MXN3))
!*** diagnostic
!      write(0,*)'ddscat ckpt 41, myid=',myid
!***

! allocate CXALPH

      IF(MYID==0)THEN   !--- begin if(myid==0)... #16
         DSTORAGE=CWORD*REAL(MXN3)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
                          ' MB for CXALPH; total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
      ENDIF   !--- end if(myid==0)... #16
!*** diagnostic
!      write(0,*)'ddscat ckpt 42, myid=',myid
!***
      ALLOCATE(CXALPH(MXN3))
!*** diagnostic
!      write(0,*)'ddscat ckpt 44, myid=',myid
!***

! allocate CXE

      IF(MYID==0)THEN   !--- begin if(myid==0)... #18
         DSTORAGE=CWORD*REAL(MXN3)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
                          ' MB for CXE   ; total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
      ENDIF   !--- end if(myid==0)... #18
!*** diagnostic
!      write(0,*)'ddscat ckpt 46, myid=',myid
!***
      ALLOCATE(CXE(MXN3))
!*** diagnostic
!      write(0,*)'ddscat ckpt 47, myid=',myid
!***

! allocate CXSC

      IF(MYID==0)THEN   !--- begin if(myid==0)... #19
         DSTORAGE=CWORD*REAL(MXCXSC)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
                          ' MB for CXSC  ; total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
      ENDIF   !--- end if(myid==0)... #19
!*** diagnostic
!      write(0,*)'ddscat ckpt 48, myid=',myid
!***
      ALLOCATE(CXSC(MXCXSC))
!*** diagnostic
!      write(0,*)'ddscat ckpt 49, myid=',myid
!***

! allocate CXSCR1

      IF(MYID==0)THEN   !--- begin if(myid==0)... #20
         DSTORAGE=CWORD*REAL(MXN3)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
                          ' MB for CXSCR1; total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
      ENDIF   !--- end if(myid==0)... #20
!*** diagnostic
!      write(0,*)'ddscat ckpt 50, myid=',myid
!***
      ALLOCATE(CXSCR1(MXN3))

!-----------------------------------------------------------------------
! following introduced by Art Lazanoff to touch the array CXZW
! which supposedly helps with allocation and distribution of this memory
! to the different processors...
! BTD: I can't really see why the following is beneficial...

!$omp parallel do
      DO J=1,MXN3
         CXSCR1(J)=(0._WP,0._WP)
      ENDDO
!$omp end parallel do
!-----------------------------------------------------------------------

!*** diagnostic
!      write(0,*)'ddscat ckpt 51, myid=',myid
!***

! allocate CXXI

      IF(MYID==0)THEN   !--- begin if(myid==0)... #21
         DSTORAGE=CWORD*REAL(MXN3)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
                          ' MB for CXXI;   total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
      ENDIF   !--- end if(myid==0)... #21
!*** diagnostic
!      write(0,*)'ddscat ckpt 52, myid=',myid
!***
      ALLOCATE(CXXI(MXN3))
!*** diagnostic
!      write(0,*)'ddscat ckpt 53, myid=',myid
!***

! allocate SCRRS1

      IF(MYID==0)THEN   !--- begin if(myid==0)... #22
         DSTORAGE=RWORD*REAL(MXN3)/MB
         STORAGE=STORAGE+DSTORAGE
         WRITE(CMSGNM,FMT='(A,F9.2,A,F9.2,A)')'allocating',DSTORAGE, &
                          ' MB for SCRRS1; total=',STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
         WRITE(CMSGNM,FMT='(16X,A,F9.2,A)')'Total memory allocation =', &
                          STORAGE,' MB'
         CALL WRIMSG('DDSCAT',CMSGNM)
      ENDIF   !--- end if(myid==0)... #22
!*** diagnostic
!      write(0,*)'ddscat ckpt 54, myid=',myid
!***
!***
!btd 08.05.08
!      ALLOCATE(SCRRS1(MXN3))
      ALLOCATE(SCRRS1(MXNAT,3))
!-----------------------------------------------------------------------
! following introduced by Art Lazanoff to touch the array CXZW
! which supposedly helps with allocation and distribution of this memory
! to the different processors...
! BTD: I can't really see why the following is beneficial...

      DO JJ=1,3
!$omp parallel do
         DO J=1,MXNAT
            SCRRS1(J,JJ)=0._WP
         ENDDO
!$omp end parallel do
      ENDDO
!-----------------------------------------------------------------------

!*** diagnostic
!      write(0,*)'ddscat ckpt 55, myid=',myid
!      call mpi_barrier(mpi_comm_world,ierr)
!***

!---------------- memory allocation is complete -----------------------

!080217 BTD change
      IF(MYID==0)THEN   !--- begin if(myid==0)... #23

         IF(JPBC>0)THEN
            DO J=1,NSCAT
               ORDERM(J)=PHIN(J)
               ORDERN(J)=THETAN(J)
            ENDDO
         ENDIF

!-----------------------------------------------------------------------
! JPBC = 0 if PBC not used
!      = 1 for PBC in y direction only
!      = 2 for PBC in z direction only
!      = 3 for PBC in y and z directions

! IPBC controls memory use.
! At this time use same approach for JPBC=1 and 2 as for JPBC=3
! May change this in future.

         IF(JPBC==0)THEN
!----------------------------------------------------------------------
! 10.01.30 BTD eliminate this restriction: not needed
!            IF(CMDFRM=='TFRAME')CALL ERRMSG('FATAL','DDSCAT'  ,      &
!               ' TFRAME = invalid option for finite targets (JPBC=0)')
!----------------------------------------------------------------------
            IPBC=0
         ELSEIF(JPBC==1.OR.JPBC==2)THEN
            IF(CMDFRM=='LFRAME')CALL ERRMSG('FATAL','DDSCAT',         &
               ' LFRAME = invalid option for periodic target (JPBC>0)')
            IPBC=1
         ELSEIF(JPBC==3)THEN
            IF(CMDFRM=='LFRAME')CALL ERRMSG('FATAL','DDSCAT',         &
               ' LFRAME = invalid option for periodic target (JPBC>0)')
            IPBC=1
         ENDIF

! Check that sufficient memory has been allocated:

         IF(IPBC==1.AND.MXPBC==0)THEN
            IF(8*NX*NY*NZ>(MXNX+1)*(MXNY+1)*(MXNZ+1))THEN
               CALL ERRMSG('FATAL','DDSCAT',                      &
                 ' CXSC too small for PBC: recompile with MXPBC=1')
            ENDIF
         ENDIF

         NORI=NTHETA*NBETA*NPHI
         NBETH=NTHETA*NBETA

! Compute angles in degrees for printout

         BETMID=DEGRAD*BETAMI
         BETMXD=DEGRAD*BETAMX
         THTMID=DEGRAD*THETMI
         THTMXD=DEGRAD*THETMX
         PHIMID=DEGRAD*PHIMIN
         PHIMXD=DEGRAD*PHIMAX

! If JPBC=0, obtain scattering vectors and scattering pol. vectors in 
! Lab Frame
!*** diagnostic
!         write(0,*)'ddscat ckpt 56, myid=',myid
!***

         IF(JPBC==0)THEN

!*** diagnostic
!            write(0,*)'in ddscat with jpbc=',jpbc,' about to call SCAVEC'
!***
            CALL SCAVEC(MXSCA,NSCAT,THETAN,PHIN,CXE01,ENSC,EM1,EM2)
!*** diagnostic
!            write(0,*)'returned to ddscat from scavec'
!***

! if JPBC = 0:
! ENSC(1-3,1-NSCAT) = scattering vectors in Lab Frame
! EM1(1-3,1-NSCAT)  = scattering polarization state 1 in Lab Frame
! EM2(1-3,1-NSCAT)  = scattering polarization state 2 in Lab Frame

         ELSE

! if JPBC > 0:
! initialize ENSC,EM1,EM2 to ensure that there are no divisions by zero
! in ROTATE.  Note that for JPBC > 0 these vectors are not actually used
! AKSR,EM1R,EM2R will be calculated later by subroutine PBCSCAVEC

            DO J=1,NSCAT
               DO JJ=1,3
                  ENSC(JJ,J)=0._WP
                  EM1(JJ,J)=0._WP
                  EM2(JJ,J)=0._WP
               ENDDO
               ENSC(1,J)=1._WP
               EM1(2,J)=1._WP
               EM2(3,J)=1._WP
            ENDDO
         ENDIF

! Note: IXYZ0(1-NAT0,1-3) contain coordinates of occupied lattice sites
!       [reordered to standard order used for arrays produced by REDUCE]
!       BETADF(1-NAT0),PHIDF(1-NAT0),THETADF(1-NAT0) contain orientation
!                         angles for anisotropic material at occupied
!                         lattice sites (not used for isotropic material

         WRITE(IDVOUT,FMT=9011)NAT,NX,NY,NZ

! Specify target orientations

! MXBETA=dimensioning information for array BETA of beta values
! MXTHET=dimensioning information for array THETA of theta values
! MXPHI=dimensioning information for array PHI of phi values
! BETAMI,BETAMX=minimum,maximum beta values (radians)
! THETMI,THETMX=minimum,maximum theta values (radians)
! PHIMIN,PHIMAX=minimum,maximum phi values (radians)
! BETA(1-NBETA)=beta values for target orientations (radians)
! THETA(1-NTHETA)=theta values for target orientations (radians)
! PHI(1-NPHI)=phi values for target orientations (radians)
! WGTA(1-NTHETA,1-NPHI)=weight for each orientation of axis a1 in LF
! WGTB(1-NBETA)=weight for each rotation of target around a1

!*** diagnostic
!         write(0,*)'in DDSCAT, about to call ORIENT'
!***
         CALL ORIENT(BETAMI,BETAMX,THETMI,THETMX,PHIMIN,PHIMAX,MXBETA,MXTHET, &
                    MXPHI,NBETA,NTHETA,NPHI,BETA,THETA,PHI,WGTA,WGTB)

!*** diagnostic
!         write(0,*)'returned to DDSCAT from ORIENT'
!***

! Open file 'qtable' for summary of Q values and
!           'mtable' for summary of refractive indices

         OPEN(UNIT=10,FILE='gammatable',STATUS='UNKNOWN') !Name changed to 'gammatable' NWB 7/12/12
!         OPEN(UNIT=11,FILE='mtable',STATUS='UNKNOWN') !qtable2 and mtable edited out NWB 7/12/12
!         OPEN(UNIT=12,FILE='qtable2',STATUS='UNKNOWN')
!*** diagnostic
!         write(0,*)'in DDSCAT, after opening qtable,mtable,qtable2'
!***

! why not use cbinfile ? and open it in writesca ?

         IF(CBINFLAG=='ALLBIN'.OR.CBINFLAG=='ORIBIN')THEN
            IOBIN=51
            OPEN(UNIT=IOBIN,FILE='dd.bin',FORM='unformatted')
         ENDIF

! Write out header for 'qtable', 'qtable2', and 'mtable' files

         WRITE(10,FMT=9020)CSTAMP,CDESCR,CMDFFT,CALPHA,CSHAPE,NAT0
         ILIN10=6
         WRITE(10,FMT=9021)(CFLEPS(J),J=1,NCOMP)
         ILIN10=ILIN10+NCOMP
         !NWB 7/12/12
!         WRITE(12,FMT=9020)CSTAMP,CDESCR,CMDFFT,CALPHA,CSHAPE,NAT0
!         ILIN12=6
!         WRITE(12,FMT=9021)(CFLEPS(J),J=1,NCOMP)
!         ILIN12=ILIN12+NCOMP
         IF(IORTH==1)THEN
            WRITE(10,FMT=9043)BETMID,BETMXD,NBETA,THTMID,THTMXD,NTHETA, &
                              PHIMID,PHIMXD,NPHI,ETASCA,NORI,IORTH
            ILIN10=ILIN10+7
            !NWB 7/12/12
!            WRITE(12,FMT=9045)BETMID,BETMXD,NBETA,THTMID,THTMXD,NTHETA, &
!                              PHIMID,PHIMXD,NPHI,NORI,IORTH
!            ILIN12=ILIN12+6
         ELSE
            WRITE(10,FMT=9043)BETMID,BETMXD,NBETA,THTMID,THTMXD,NTHETA, &
                              PHIMID,PHIMXD,NPHI,ETASCA,NORI,IORTH
            ILIN10=ILIN10+7
            !NWB 7/12/12
!            WRITE(12,FMT=9046)BETMID,BETMXD,NBETA,THTMID,THTMXD,NTHETA, &
!                              PHIMID,PHIMXD,NPHI,NORI,IORTH
!            ILIN12=ILIN12+6
         ENDIF
         CLOSE(10)
         CLOSE(12)
         !NWB 7/12/12
!         WRITE(11,FMT=9020)CSTAMP,CDESCR,CMDFFT,CALPHA,CSHAPE,NAT0
!         WRITE(11,FMT=9021)(CFLEPS(J),J=1,NCOMP)
!         WRITE(11,FMT=9048)
      ENDIF             !--- end if(myid==0)... #23

!*** The master process needs to share information from the above section
!    with the slave processes. All processes call this subroutine, which
!    does nothing if MPI is not being used.

!*** diagnostic
!      write(0,*)'DDSCAT ckpt 56.1, myid=',myid,' ICOMP(1)  =',ICOMP(1)
!      write(0,*)'DDSCAT ckpt 56.2, myid=',myid,' ICOMP(MXN3)=', &
!                ICOMP(MXN3),' MXN3=',MXN3
!      write(0,*)'DDSCAT ckpt 56.3, myid=',myid,' about to call SHARE1...'
!      write(0,*)'                  nz=',nz
!***
! this mpi_barrier is required to ensure that NX, etc computed by MYID=0
! will be evaluated before being SHAREd to the other threads

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

      CALL SHARE1(A1,A2,A3,AEFFA,BETA,BETADF,BETMID,BETMXD,CALPHA,CBINFLAG,   &
                  CDESCR,CFLEPS,CMDFFT,CMDFRM,CMDSOL,CMDTRQ,CSHAPE,           &
                  CXE01,CXE02,DAEFF,DX,ENSC,EM1,EM2,ETASCA,GAMMA,IANISO,      &
                  ICOMP,IDVERR,IDVOUT,ILIN10,ILIN12,INIT,IOBIN,IOCC,IORI0,    &
                  IORTH,IPBC,IRAD0,IWAV0,IWRKSC,IWRPOL,IXYZ0,JPBC,MXBETA_SH,  &
                  MXCOMP_SH,MXN3,MXNAT,MXNX,MXNY,MXNZ,MXPBC_SH,MXPHI_SH,      &
                  MXRAD_SH,MXSCA_SH,MXTHET_SH,MXWAV_SH,MYID,NAT,NAT0,NAT3,    &
                  NBETA,NBETH,NCOMP,NORI,NPHI,NRAD,NSCAT,NSMELTS,NTHETA,NWAV, &
                  NX,NY,NZ,ORDERM,ORDERN,PHI,PHIDF,PHIMID,PHIMXD,PHIN,PYD,    &
                  PYDDX,PZD,PZDDX,SHPAR,SMIND1,SMIND2,THETA,THETADF,THETAN,   &
                  THTMID,THTMXD,TOL,WAVEA,WGTA,WGTB,X0,XMAX,XMIN,YMAX,YMIN,   &
                  ZMAX,ZMIN)

!*** diagnostic
!      write(0,*)'ddscat ckpt 57, myid=',myid
!***
! IWAV0,IRAD0,IORI0 were obtained from REAPAR

!*** diagnostic
!      write(0,*)'returned to DDSCAT from SHARE1'
!***
      IWAV1=IWAV0+1
      IRAD1=IRAD0+1
      IORI1=IORI0+1
      ITHETA1=1+IORI0/(NBETA*NPHI)
      IBETA1=1+IORI0/NPHI-(ITHETA1-1)*NBETA
      IBETH1=1+IORI0/NPHI
      IPHI1=IORI1-(ITHETA1-1)*NBETA*NPHI-(IBETA1-1)*NPHI

!*** Loop over wavelengths:

      DO IWAV=IWAV1,NWAV                         !----- loop over IWAV -------
         WAVE=WAVEA(IWAV)

!*** Determine dielectric properties of target

         IF(MYID==0)THEN

!*** diagnostic
!            write(0,*)'DDSCAT ckpt 57.1, about to call DIELEC'
!***
            CALL DIELEC(WAVE,IDVOUT,CFLEPS,CXEPS,MXCOMP,MXWAVT,NCOMP, &
                        E1A,E2A,WVA)
!*** diagnostic
!            write(0,*)'DDSCAT ckpt 57.2 returned to DDDSCAT from DIELEC'
!***
         ENDIF

!*** Again the master process needs to share information with the slaves

!*** diagnostic
!         write(0,*)'DDSCAT ckpt 57.3, about to call SHARE2'
!***
! 080723 BTD add mpi_barrier to ensure that E1A,E2A have been read by MYID=0
!            before being SHAREd by other threads.

         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

         CALL SHARE2(MXCOMP,MXWAVT,CXEPS,E1A,E2A,WVA)

!*** diagnostic
!         write(0,*)'DDSCAT ckpt 57.4, returned to DDSCAT from SHARE2'
!***

!*** Obtain complex polarizabilities

         DO J=1,NCOMP
            CXRFR(J)=SQRT(CXEPS(J))
            IF(MYID==0)THEN
               WRITE(IDVOUT,FMT=9031)CXRFR(J),CXEPS(J),J
            ENDIF
         ENDDO

         DO IRAD=IRAD1,NRAD                     !**** Loop over IRAD ----------

!*** diagnostic
!            write(0,*)'in DDSCAT, IRAD=',IRAD
!***

            AEFF=AEFFA(IRAD)
            XX=2._WP*PI*AEFF/WAVE

! Compute AK1=length of k vector in "natural" units=k*d
! (natural unit of length = d = lattice spacing)
! Remember that NAT = number of dipoles including "vacuum" sites in
!                     extended target
!               NAT0= number of dipoles in "real" target

            AK1=XX*DAEFF

            IF(MYID==0)THEN
               WRITE(IDVOUT,FMT=9032)AEFF,DAEFF,WAVE,XX,AK1
            ENDIF

            AK3=AK1**3

! PIA2=pi*(aeff/d)**2

            PIA2=PI*(.75_WP*NAT0/PI)**(2._WP/3._WP)

! Initialize various sums over target orientation.

            DO JO=1,IORTH
               QSCSUM(JO)=0._WP
               QABSUM(JO)=0._WP
               QEXSUM(JO)=0._WP
               QBKSUM(JO)=0._WP
               QPHSUM(JO)=0._WP
               QSCSUM_1(JO)=0._WP
               QABSUM_1(JO)=0._WP
               QEXSUM_1(JO)=0._WP
               QBKSUM_1(JO)=0._WP
               QPHSUM_1(JO)=0._WP
               QSCG2SUM(JO)=0._WP
               QSCG2SUM_1(JO)=0._WP
               DO J=1,3
                  QSCGSUM(J,JO)=0._WP
                  QTRQABSUM(J,JO)=0._WP
                  QTRQSCSUM(J,JO)=0._WP
                  QSCGSUM_1(J,JO)=0._WP
                  QTRQABSUM_1(J,JO)=0._WP
                  QTRQSCSUM_1(J,JO)=0._WP
               ENDDO
            ENDDO
            DO IDIR=1,NSCAT
               S1111(IDIR)=0._WP
               S2121(IDIR)=0._WP
               S1111_1(IDIR)=0._WP
               S2121_1(IDIR)=0._WP
               CX1121(IDIR)=(0._WP,0._WP)
               CX1121_1(IDIR)=(0._WP,0._WP)
            ENDDO
            DO IDIR=1,NSCAT
               DO J=1,4
                  DO JO=1,4
                     SMORI(JO,J,IDIR)=0._WP
                     SMORI_1(JO,J,IDIR)=0._WP
                  ENDDO
               ENDDO
            ENDDO

! Loop over target orientations (target rotations in Lab Frame)
! (actually, loop over lab rotations in Target Frame)

!**** The outer orientation loop over IBETH=1,NBETH is divided up among
!     the parallel processes.

            IORI=IORI0

! ITNUMMX will store the maximum number of iterations used for any of
! the orientations for current size and wavelength, so that this can
! be written to the waarbbkccc.sca output file

            DO J=1,2
               ITNUMMX(J)=0
            ENDDO

            DO IBETH=MYID+IBETH1,NBETH,NUMPROCS  !---- loop over IBETH -------
!*** diagnostic
!               write(0,*)'DDSCAT ckpt 57.5, IBETH=',IBETH
!***
               ITHETA=INT((IBETH-1)/NBETA)+1
               THETAD=THETA(ITHETA)*DEGRAD
               IBETA=MOD(IBETH-1,NBETA)+1
               BETAD=BETA(IBETA)*DEGRAD

               DO IPHI=IPHI1,NPHI                 !---- loop over PHI --------

!*** diagnostic
!                  write(0,*)'DDSCAT ckpt 57.6, IPHI=',IPHI
!***
                  PHID=PHI(IPHI)*DEGRAD
                  IORI=(ITHETA-1)*NBETA*NPHI+(IBETA-1)*NPHI+IPHI

!-----------------------------------------------------------------------
!**** For specified target rotations BETA,THETA,PHI in the Lab Frame,
!     and predefined 

!     ENSC(1-3,1-NSCAT) scattering directions in the Lab Frame 
!     EM1(1-3,1-NSCAT) scattered pol. vectors 1 in the Lab Frame
!     EM2(1-3,1-NSCAT) scattered pol. vectors 2 in the Lab Frame

!     subroutine ROTATE computes

!     AKSR(1-3,1-NSCAT) propagation vectors in the Target Frame
!     EM1R(1-3,1-NSCAT) scattered pol vectors 1 in Target Frame
!     EM2R(1-3,1-NSCAT) scattered pol vectors 2 in Target Frame

!     Note: subroutine ROTATE assumes that the angles THETAN,PHIN define
!     scattering directions in the Lab Frame -- if these values are
!     actually for the Target Frame, values of AKSR,EM1R,EM2R returned
!     by ROTATE will later be replaced by correct values.

!*** diagnostic
!                  write(0,*)'DDSCAT ckpt 57.65, about to call ROTATE'
!***
                  CALL ROTATE(A1,A2,AK1,CXE01,CXE02,BETA(IBETA),          &
                              THETA(ITHETA),PHI(IPHI),EN0R,CXE01R,CXE02R, &
                              MXSCA,NSCAT,ENSC,EM1,EM2,AKSR,EM1R,EM2R)
!*** diagnostic
!                  write(0,*)'DDSCAT ckpt 57.66, myid=',myid,              &
!                            ' BETA,THETA,PHI=',BETA(IBETA),THETA(ITHETA), &
!                            PHI(IPHI)
!                  write(0,*)'myid=',myid,' en0r(1-3)=',en0r
!*** end diagnostic
!                  write(0,*)' returned from ROTATE'
!                  write(0,9113)cxe01,cxe02,cxe01r,cxe02r
! 9113 format('cxe01 =',2f8.5,2x,2f8.5,2x,2f8.5,/, &
!             'cxe02 =',2f8.5,2x,2f8.5,2x,2f8.5,/, &
!             'cxe01r=',2f8.5,2x,2f8.5,2x,2f8.5,/, &
!             'cxe02r=',2f8.5,2x,2f8.5,2x,2f8.5)
!*** end diagnostic

! calculate ENSCR = unit scattering vector in Target Frame

                  DO IDIR=1,NSCAT
                     DO J=1,3
                        ENSCR(J,IDIR)=AKSR(J,IDIR)/AK1
                     ENDDO
                  ENDDO

! calculate XLR,YLR,ZLR = lab unit vectors XL,YL,ZL in Target Frame

                  DO J=1,3
                     SINTHE=SIN(THETA(ITHETA))
                     COSTHE=COS(THETA(ITHETA))
                     SINPHI=SIN(PHI(IPHI))
                     COSPHI=COS(PHI(IPHI))
                     SINBET=SIN(BETA(IBETA))
                     COSBET=COS(BETA(IBETA))
                     XLR(J)=A1(J)*COSTHE-A2(J)*SINTHE*COSBET+                &
                            A3(J)*SINTHE*SINBET
                     YLR(J)=A1(J)*SINTHE*COSPHI+A2(J)*(COSTHE*COSBET*COSPHI- &
                            SINBET*SINPHI)-A3(J)*(COSTHE*SINBET*COSPHI+      &
                            COSBET*SINPHI)
                     ZLR(J)=A1(J)*SINTHE*SINPHI+A2(J)*(COSTHE*COSBET*SINPHI+ &
                            SINBET*COSPHI)-A3(J)*(COSTHE*SINBET*SINPHI-      &
                            COSBET*COSPHI)
                  ENDDO

! If CMDFRM = 'TFRAME' then angles THETAN,PHIN define scattering
!                      directions in the Target Frame, so ignore
!                      vectors AKSR,EM1R,EM2R calculated by ROTATE,
!                      and instead simply use ENSC,EM1,EM2

                  IF(CMDFRM=='TFRAME')THEN
                     IF(JPBC==0)THEN
                        DO IDIR=1,NSCAT
                           DO J=1,3
                              AKSR(J,IDIR)=AK1*ENSC(J,IDIR)
                              EM1R(J,IDIR)=EM1(J,IDIR)
                              EM2R(J,IDIR)=EM2(J,IDIR)
                           ENDDO
                        ENDDO
                     ELSE

! JPBC= 1, 2, or 3: call PBCSCAVEC to calculate scattering vectors 
!                        AKSR,EM1R,EM2R in TF,
!                        ENSC,EM1,EM2 in LF

!*** diagnostic
!                        write(0,*)'DDSCAT ckpt 57.67, ',        & 
!                                  'about to call pbcscavec with'
!                        write(0,*)'pyddx=',pyddx
!                        write(0,*)'pzddx=',pzddx
!                        write(0,*)'          phin(1)=',phin(1)
!                        write(0,*)'          phin(2)=',phin(2)
!                        write(0,9393)xlr,ylr,zlr
!  9393 format('xlr(1-3)=',3f8.5,/,'ylr(1-3)=',3f8.5,/,'zlr(1-3)=',3f8.5)
!***
                        CALL PBCSCAVEC(MXSCA,JPBC,NSCAT,PYDDX,PZDDX,A1,A2, &
                                       THETA(ITHETA),BETA(IBETA),          &
                                       XLR,YLR,ZLR,AK1,EN0R,CXE01R,ORDERM, &
                                       ORDERN,THETAN,PHIN,AKSR,EM1R,EM2R,  &
                                       ENSC,EM1,EM2)
!*** diagnostic
!                        write(0,*)'DDSCAT ckpt 57.68: returned from pbcscavec'
!                        write(0,*)'                   phin(1)=',phin(1)
!                        write(0,*)'                   phin(2)=',phin(2)
!                        do idir=1,2
!                           write(0,9411)idir,                         &
!                              ensc(1,idir),ensc(2,idir),ensc(3,idir), &
!                              aksr(1,idir),aksr(2,idir),aksr(3,idir), &
!                              em1(1,idir),em1(2,idir),em1(3,idir),    &
!                              em1r(1,idir),em1r(2,idir),em1r(3,idir), &
!                              em2(1,idir),em2(2,idir),em2(3,idir),    &
!                              em2r(1,idir),em2r(2,idir),em2r(3,idir)
!                        enddo
! 9411 format('nd=',i2,'ensc=',3f9.5,' aksr=',3f9.5,/, &
!                   5x,' em1=',3f9.5,' em1r=',3f9.5,/, &
!                   5x,' em2=',3f9.5,' em2r=',3f9.5)
!*** end diagnostic
                     ENDIF
                  ENDIF

!-----------------------------------------------------------------------

                  DO J=1,3
                     AKR(J)=AK1*EN0R(J)
                  ENDDO

!*** diagnostic
!                  write(0,*)'ddscat ckpt 58, myid=',myid, &
!                            ' ak1=',ak1,' en0r(1-3)=',en0r
!***
! note: if called with IORI > 999, NAMER will have trouble when
!       generating CFLSCA.  If IWRKSC=0, this is immaterial, because
!       CFLSCA will not be used.
!       We will check in REAPAR to halt if user requests IWRKSC=1 and
!       more than 1000 orientations.

!*** diagnostic
!                  write(0,*)'in DDSCAT, about to call NAMER'
!***
                  IF(IWRKSC==1)THEN
                     CALL NAMER(IWAV,IRAD,IORI,CFLPOL1,CFLPOL2,CFLSCA,CFLAVG)
                  ELSE

! 100128 BTD change IORI from 0 to 1 in call to NAMER when IWRKSC==0

                     CALL NAMER(IWAV,IRAD,1,CFLPOL1,CFLPOL2,CFLSCA,CFLAVG)
!*** diagnostic
!                     write(0,*)'ddscat ckpt 58.5, CFLPOL1=',CFLPOL1
!                     write(0,*)'                  CFLPOL2=',CFLPOL2
!                     write(0,*)'                  CFLAVG =',CFLAVG
!***
                  ENDIF
!*** diagnostic
!                  write(0,*)'returned to DDSCAT from NAMER'
!***
!*** diagnostic
!                  write(0,*)'about to call getfml at ckpt alpha'
!                  write(0,9220)(theta(itheta)*180./pi), &
!                        (beta(ibeta)*180./pi),          &
!                        (phi(iphi)*180./pi),            &
!                        akr,cxe01r,cxe02r,              &
!                        aksr(1,1),aksr(2,1),aksr(3,1),  &
!                        em1r(1,1),em1r(2,1),em1r(3,1),  &
!                        em2r(1,1),em2r(2,1),em2r(3,1)
! 9220 format('=========================================',/, &
!             'theta=',f10.5,' beta=',f10.5,' phi=',f10.5,/, &
!             '  akr =',3f10.6,/,                            &
!             'cxe01r=',2f9.5,' ,',2f9.5,' ,',2f9.5,/,       &
!             'cxe01r=',2f9.5,' ,',2f9.5,' ,',2f9.5,/,       &
!             'aksr(j,1)=',3f10.6,/,                         &
!             'em1r(j,1)=',3f10.6,/,                         &
!             'em2r(j,1)=',3f10.6)
!***
!*** diagnostic
!                  write(0,FMT='(A,I3,A,F10.6,A,3F10.6)')                     &
!                     'ddscat ckpt 59, myid=',myid,' ak1=',ak1,' akr(1-3)=',akr
!                  write(0,*)'mxnat=',mxnat
!                  write(0,*)'about to call GETFML...'
!***
                  CALL GETFML(AKR,AK3,AKSR,BETADF,PHIDF,THETADF,DX,X0,CALPHA, &
                              CMDSOL,CMDFFT,CMDTRQ,CSHAPE,CXADIA,CXAOFF,      &
                              CXALPH,CXALOF,CXALOS,CXE,CXE01R,CXE02R,CXEPS,   &
                              CXF11,CXF12,CXF21,CXF22,CXRLOC,CXXI,CXSC,       &
                              CXSCR1,CXZC,CXZW,EM1R,EM2R,ETASCA,GAMMA,IBETH,  &
                              IBETH1,ICOMP,IDVOUT,INIT,IOCC,IORTH,IPBC,IPHI,  &
                              IPHI1,ITASK,IXYZ0,ITNUM,JPBC,MXITER,LACE,LAXI,  &
                              LCLM,LGI,LPI,LQI,LSC0,MXCOMP,MXCXSC,MXN3,MXNAT, &
                              MXNX,MXNY,MXNZ,MXPBC,MXPHI,MXSCA,MYID,NAT,NAT0, &
                              NAT3,NAVG,NCOMP,NSCAT,NX,NY,NZ,PHI,PIA2,QABS,   &
                              QBKSCA,QEXT,QPHA,QSCAT,QSCAG,QSCAG2,QTRQAB,     &
                              QTRQSC,SCRRS1,SCRRS2,SHPAR,TOL,TIMERS,MXTIMERS, &
                              NTIMERS,AEFFA,WAVEA,MXRAD,MXWAV,CENTER,IWRPOL,  &
                              c,h_bar,h_bar2,velocity,e_charge,DielectricConst) 
                              !Arguments AEFFA and after added by NWB 3/8/12

!*** diagnostic
!                  write(0,*)'DDSCAT ckpt 60, myid=',myid
!***
!                  write(0,9221)cxf11(1),cxf12(1),cxf21(1),cxf22(1)
! 9221 format('cxf11(1)=',2f10.6,' cxf21(1)=',2f10.6,/, &
!             'cxf12(1)=',2f10.6,' cxf22(1)=',2f10.6)
!***
! if IORTH=1, GETFML returns reduced polarization array for orientation
!             IPHI=IPHI1 and JO=1 (incident wave CXE01R) in CXSCR1,
!             and scattering amplitude factors CXF11,CXF21 for NSCAT
!             directions
! if IORTH=2, GETFML also returns reduced polarization array for
!             JO=2 (incident wave CXE02R) in CXSC(LACE)
!             and scattering amplitude factors CXF12,CXF22

                  IF(IWRKSC==1)THEN
!*** diagnostic
!                     write(0,*)'DDSCAT ckpt 61, myid=',myid
!***
                     CALL NAMER2(IWAV,IRAD,IORI,CFLFML)
!*** diagnostic
!                     write(0,*)'DDSCAT ckpt 62, myid=',myid
!***
                     CALL WRITEFML(JPBC,IORI,IORTH,IRAD,IWAV,MXCOMP,MXSCA,   &
                                   NAT0,NAVG,NCOMP,NSCAT,CALPHA,CDESCR,      &
                                   CMDFFT,CMDFRM,CSHAPE,CSTAMP,AEFF,AK1,AKR, &
                                   BETAD,PHID,THETAD,TOL,WAVE,XX,A1,A2,PHIN, &
                                   THETAN,CXE01,CXE01R,CXE02,CXE02R,CXEPS,   &
                                   CXRFR,CXF11,CXF21,CXF12,CXF22,PYD,PZD)
!*** diagnostic
!                     write(0,*)'DDSCAT ckpt 63, myid=',myid
!***
                  ENDIF

                  DO J=1,2
                     IF(ITNUM(J)>ITNUMMX(J))ITNUMMX(J)=ITNUM(J)
                  ENDDO

! Compute Mueller scattering matrix

!*** diagnostic
!                     write(0,*)'DDSCAT ckpt 64, myid=',myid
!***
                  CALL GETMUELLER(IBETA,IORTH,IPHI,ITHETA,JPBC,MXBETA,MXSCA,  &
                                  MXPHI,MXTHET,NSCAT,ORDERM,ORDERN,CMDTRQ,    &
                                  AK1,AKSR,ENSC,ENSCR,PYDDX,PZDDX,PHIN,SM,    &
                                  SMORI_1,S1111_1,S2121_1,CX1121_1,CXE01,     &
                                  CXE02,CXF11,CXF12,CXF21,CXF22,CXS1,CXS2,    &
                                  CXS3,CXS4,QABS,QABSUM_1,QBKSCA,QBKSUM_1,    &
                                  QEXSUM_1,QEXT,QPHA,QPHSUM_1,QSCAG,QSCAG2,   &
                                  QSCAT,QSCG2SUM_1,QSCGSUM_1,QSCSUM_1,QTRQAB, &
                                  QTRQABSUM_1,QTRQSC,QTRQSCSUM_1,WGTA,WGTB,   &
                                  EM1,EM2,EM1R,EM2R)
!*** diagnostic
!                  write(0,*)'DDSCAT ckpt 65, myid=',myid, &
!                            'IORTH=',IORTH,' QSCSUM_1(1-2)=',QSCSUM_1
!***
!----------------------------------------------------------------------
!**** Choose whether or not to write out scattering properties for
!     current orientation; conditional may be replaced if desired.

!**** If IWRKSC=1, write scattering properties of current orientation:

                  IF(IWRKSC==1)THEN

!*** diagnostic
!                     write(0,*)'DDSCAT ckpt 66, myid=',myid
!***
                       CALL WRITESCA(MXNX,MXNY,MXNZ,NX,NY,NZ,DX,ICOMP,IXYZ0,  &
                          JPBC,MXNAT,MXN3,MYID,NAT,NAT3,WAVEA,MXWAV,NWAV,     &
                          AEFFA,MXRAD,NRAD,ITHETA,IBETA,IPHI,THETA,MXTHET,    &
                          NTHETA,BETA,MXBETA,NBETA,PHI,MXPHI,NPHI,NSMELTS,    &
                          TIMERS,MXTIMERS,NTIMERS,CBINFILE,IOBIN,CBINFLAG,    &
                          IDNC,ILIN10,ILIN12,IORI,IWRKSC,   &
                          IORTH,IRAD,IWAV,MXCOMP,MXSCA,NAT0,NAVG,ITNUM,       &
                          MXITER,NCOMP,NORI,NSCAT,CALPHA,CDESCR,CFLAVG,       &
                          CFLSCA,CMDFFT,CMDFRM,CMDTRQ,CSHAPE,CSTAMP,AEFF,AK1, &
                          AKR,BETAD,BETMID,BETMXD,ETASCA,PHID,PHIMID,PHIMXD,  &
                          THETAD,THTMID,THTMXD,TOL,WAVE,XX,A1,A2,PHIN,QABS,   &
                          QABSUM_1,QBKSCA,QBKSUM_1,QEXSUM_1,QEXT,QPHA,        &
                          QPHSUM_1,QSCAG,QSCAG2,QSCAT,QSCGSUM_1,QSCG2SUM_1,   &
                          QSCSUM_1,QTRQAB,QTRQABSUM_1,QTRQSC,QTRQSCSUM_1,     &
                          S1111_1,S2121_1,SM,SMIND1,SMIND2,SMORI_1,THETAN,    &
                          CX1121,CXE01,CXE01R,CXE02,CXE02R,CXEPS,CXRFR,CXF11, &
                          CXF21,CXF12,CXF22,PYDDX,PZDDX,XMIN,XMAX,YMIN,YMAX,  &
                          ZMIN,ZMAX)
!*** diagnostic
!                       write(0,*)'DDSCAT ckpt 67, myid=',myid
!***
                  ENDIF

! only write out polarization array if IPHI=1 (for IPHI>1, can
! reconstruct polarization array by linear combination of CFPOL1 and
! CFPOL2)

                  IF(IWRPOL==1.AND.IPHI==1)THEN
!*** diagnostic
!                     write(0,*)'DDSCAT ckpt 68, myid=',myid
!***
                     CALL WRITEPOL(NRWORD,MXNAT,NX,NY,NZ,NAT0,IANISO,  &
                                   ICOMP,IXYZ0,PYD,PZD,AKR,DX,X0,WAVE, &
                                   BETADF,THETADF,PHIDF,CXE01R,CXADIA, &
                                   CXAOFF,CXSCR1,CFLPOL1)
!*** diagnostic
!                     write(0,*)'DDSCAT ckpt 69, myid=',myid
!***
                     IF(IORTH==2)THEN
                        CALL WRITEPOL(NRWORD,MXNAT,NX,NY,NZ,NAT0,IANISO,  &
                                      ICOMP,IXYZ0,PYD,PZD,AKR,DX,X0,WAVE, &
                                      BETADF,THETADF,PHIDF,CXE02R,CXADIA, &
                                      CXAOFF,CXSC(LACE),CFLPOL2)
                     ENDIF
                  ENDIF
               ENDDO ! end loop over IPHI

               IPHI1=1

            ENDDO ! end loop over IBETH

            IBETH1=1

!*** Collect the partial sums over target orientation
!     from the parallel processes. This routine is a dummy
!     if MPI is not being used.

!*** diagnostic
!            write(0,*)'DDSCAT ckpt 70, myid=',myid
!***
            CALL COLSUM(IORTH,MYID,MXSCA,NSCAT,QSCSUM,QSCSUM_1,QABSUM,     &
                        QABSUM_1,QEXSUM,QEXSUM_1,QBKSUM,QBKSUM_1,QPHSUM,   &
                        QPHSUM_1,QSCG2SUM,QSCG2SUM_1,QSCGSUM,QSCGSUM_1,    &
                        QTRQABSUM,QTRQABSUM_1,QTRQSCSUM,QTRQSCSUM_1,S1111, &
                        S1111_1,S2121,S2121_1,CX1121,CX1121_1,SMORI,SMORI_1)
!*** diagnostic
!            write(0,*)'DDSCAT ckpt 71, myid=',myid
!***
!*** Write dielectric information to 'mtable'

            FREQ=1.E4_WP/WAVE
            CXEN=SQRT(CXEPS(1))
            IF(MYID==0)THEN
                !WRITE(11,FMT=9058)WAVE,FREQ,CXEN,CXEPS(1) !Edited out by NWB 7/12/12
            ENDIF

! Set IORI=0 prior to calling WRITESCA in order to print orientational
! averages

            IORI=0

!******** Now print out orientational averages **********

            IF(MYID==0)THEN

               CALL WRITESCA(MXNX,MXNY,MXNZ,NX,NY,NZ,DX,ICOMP,IXYZ0,JPBC,     &
                             MXNAT,MXN3,MYID,NAT,NAT3,WAVEA,MXWAV,NWAV,AEFFA, &
                             MXRAD,NRAD,ITHETA,IBETA,IPHI,THETA,MXTHET,       &
                             NTHETA,BETA,MXBETA,NBETA,PHI,MXPHI,NPHI,NSMELTS, &
                             TIMERS,MXTIMERS,NTIMERS,CBINFILE,IOBIN,CBINFLAG, &
                             IDNC,ILIN10,ILIN12,IORI,       &
                             IWRKSC,IORTH,IRAD,IWAV,MXCOMP,MXSCA,NAT0,NAVG,   &
                             ITNUMMX,MXITER,NCOMP,NORI,NSCAT,CALPHA,CDESCR,   &
                             CFLAVG,CFLSCA,CMDFFT,CMDFRM,CMDTRQ,CSHAPE,       &
                             CSTAMP,AEFF,AK1,AKR,BETAD,BETMID,BETMXD,ETASCA,  &
                             PHID,PHIMID,PHIMXD,THETAD,THTMID,THTMXD,TOL,     &
                             WAVE,XX,A1,A2,PHIN,QABS,QABSUM,QBKSCA,QBKSUM,    &
                             QEXSUM,QEXT,QPHA,QPHSUM,QSCAG,QSCAG2,QSCAT,      &
                             QSCGSUM,QSCG2SUM,QSCSUM,QTRQAB,QTRQABSUM,QTRQSC, &
                             QTRQSCSUM,S1111,S2121,SM,SMIND1,SMIND2,SMORI,    &
                             THETAN,CX1121,CXE01,CXE01R,CXE02,CXE02R,CXEPS,   &
                             CXRFR,CXF11,CXF21,CXF12,CXF22,PYDDX,PZDDX,XMIN,  &
                             XMAX,YMIN,YMAX,ZMIN,ZMAX)

            ENDIF

        ENDDO ! end loop over IRAD
        IRAD1=1

      ENDDO ! end loop over IWAV

! close binary file (this perhaps should be moved to writesca)

      IF(MYID==0)THEN
         IF(CBINFLAG/='NOTBIN')THEN
            CALL WRIMSG('DDSCAT',' close dd.bin')
            CLOSE(IOBIN)
         ENDIF
!         CLOSE(11) ! close mtable
      ENDIF

      CALL WRIMSG('DDSCAT',' terminates normally')

! MPI environment shutdown:

      CALL MPI_FINALIZE(IERR)

      STOP

9000  FORMAT(' >DDSCAT --- ',A)
9011  FORMAT(7X,I7,' = NAT  = number of dipoles in extended target'/7X,3I4, &
        ' = x,y,z length of extended target (Targ. Frame)')
9020  FORMAT(' DDSCAT --- ',A/' TARGET ---',A/' ',A,                      &
        ' --- method of solution '/' ',A,                                 &
        ' --- prescription for polarizabilities'/' ',A,' --- shape '/I7,  &
        ' = NAT0 = number of dipoles')
9021  FORMAT(' ',A)
9031  FORMAT(' m= (',F7.4,' ,',F7.4,'),  epsilon= (',F8.4,' , ',F7.4,     &
         ')  for material',I2)
9032  FORMAT(F12.6,' = AEFF = effective radius (physical units)',/, &
         F12.6,' = d/aeff for this target',/,                       &
         F12.6,' = WAVE = wavelength (physical units)',/,           &
         F12.6,' = k*aeff = 2*pi*aeff/lambda',/,                    &
         F12.6,' = k*d')
9043  FORMAT(2F8.3,' = beta_min, beta_max ;  NBETA =',I2/2F8.3,               &
        ' = theta_min, theta_max; NTHETA=',I2/2F8.3,                          &
        ' = phi_min, phi_max   ;   NPHI =',I2/F7.4,                           &
        ' = ETASCA (param. controlling # of ',                                &
        'scatt. dirs used to calculate <cos> etc.'/' Results averaged over ', &
        I4,' target orientations'/'                   and ',I2,               &
        ' incident polarizations'/1X,'eq. wave',2X,'Energy loss',2X,'Gamma') !Output edited by NWB 7/12/12
9045  FORMAT(2F8.3,' = beta_min, beta_max ;  NBETA =',I2/2F8.3,             &
        ' = theta_min, theta_max; NTHETA=',I2/2F8.3,                        &
        ' = phi_min, phi_max   ;   NPHI =',I2/' Results averaged over ',I4, &
        ' target orientations'/'                   and ',I2,                &
        ' incident polarizations'/'   aeff       wave      Q_pha')
9046  FORMAT(2F8.3,' = beta_min, beta_max ;  NBETA =',I2/2F8.3,              &
        ' = theta_min, theta_max; NTHETA=',I2/2F8.3,                         &
        ' = phi_min, phi_max   ;   NPHI =',I2/' Results averaged over ',I4,  &
        ' target orientations'/'                   and ',I2,                 &
        ' incident polarizations'/'   aeff',7X,'wave',7X,'Q_pha',7X,'Q_pol', &
        7X,'Q_cpol')
9048  FORMAT(' wave(um)   f(cm-1)     Re(m)     Im(m)    Re(eps)   ', &
        'Im(eps)')
9058  FORMAT(1P,E10.4,E11.4,4E10.3)

    END PROGRAM DDSCAT
