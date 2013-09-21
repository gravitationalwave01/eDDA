    SUBROUTINE VERSION(CSTAMP)
      IMPLICIT NONE
! Arguments:
      CHARACTER :: CSTAMP*26

      CSTAMP = 'DDSCAT 7.1.0 [10.03.03]'
!***********************************************************************
! History of major changes to DDSCAT package:
! 90.12.21 (BTD): Added code to DDSCAT to create output files
!                 'qtable' (summary of Q values)
!                 'mtable' (summary of diel. func. for material 1)
! 91.08.14 (BTD): Add QBKSCA to arg. list for GETFML
!                 Modify code to print out Q_bk
! 91.08.15 (BTD): Provide different headings for qtable depending on
!                 whether IORTH=1 or 2 (add statement 9043 to DDSCAT).
! 91.09.17 (BTD): Remove calls to ALPHA and EVALA in DDSCAT (needed to
!                 move these to subroutine GETFML since alpha from Latti
!                 Dispersion Relation depends on propagation direction
!                 and polarization state)
!                 Add CALPHA,CXEPS,ICOMP,MXCOMP to argument list for
!                 GETFML (information required by ALPHA)
! 91.11.12 (BTD): Added variable IDVSHP to argument list for TARGET
!                 (device number for REASHP to use in reading shape file
! 93.01.15 (BTD): Modified DDSCAT to divided output file qtable into 2 f
!                 qtable (containing Q_ext,Q_abs,Q_sca,g,Q_bk) and
!                 qtable2 (containing Q_pha, Q_pol, Q_cpol)
! 93.01.16 (BTD): Modify DDSCAT so that qtable, qtable2, and mtable are
!                 after writing, and reopened for each new write
! 93.01.20 (BTD): Add MXNX, MXNY, MXNZ to argument list for subroutine
!                 EXTEND to permit checking of target size against
!                 maximum allowed dimensions
! 93.03.11 (BTD): Deleted all code associated with unused variable
!                 CMACHN (originally included to identify machine/OS)
! 93.03.12 (BTD): Changed CDESCR*60 -> CDESCR*67
!                 Moved WRITE(IDVOUT,9010) to subr. TARGET
! 93.06.02 (BTD): Corrected error in FORMAT statement 9045 in DDSCAT
! 93.06.03 (BTD): Corrected error in computation of angle-averaged
!                 backscattering cross section QBKSUM(JO) in DDSCAT (had
!                 neglected to reset sum to zero for each new target).
! 93.09.28 (BTD): Add "fix" in DDSCAT.f to work around Sun compiler/OS
!                 bug (see description below.
! 93.12.15 (BTD): Corrected calculation of PPOL in DDSCAT.f  Added
!                 comments on correspondence between our notation and
!                 elements of 2x2 amplitude scattering matrix and 4x4
!                 Mueller scattering matrix.
! 94.01.27 (BTD): Replaced SHPAR1,SHPAR2,SHPAR3 by SHPAR(1-6) to
!                 allow up to 6 shape parameters in DDSCAT, REAPAR,
!                 TARGET.
! 94.06.20 (PJF): Version 4c.1:
!                 Modifications to various routines to support use of
!                 GPFA FFT code of Temperton by specifying option
!                 NEWTMP in ddscat.par.  Changes made in
!                 cxfft3, eself, extend, reapar
!                 Add FFT timing code TSTFFT to distribution.
!                 Change made to timeit_xxx routines for compatibility
!                 with TSTFFT
! 94.12.20 (BTD): Add subroutine VERSION to consolidate log of
!                 significant changes to DDSCAT package.
! 95.04.07 (BTD): REASHP: Changed CDESCR*60->CDESCR*67 for compatibility
!                         with calling routine TARGET
!                 TARCYL: corrected error in code
! 95.05.15 (BTD): REASHP: Corrected READ statement
!                         added module for inhomogeneous/anisotropic
!                         targets.
! 95.05.26 (BTD): REAPAR: corrected construction of orthogonal
!                         polarization when E01 is complex (e.g.,
!                         circular polarization).
! 95.06.14 (BTD): REAPAR: added code to check that user is not requestin
!                         computation of f_ml when using other than
!                         linearly polarized E01
! 95.06.14 (BTD): Version 5a.1: first development
!                 Differs from version 4d in computing:
!                 full radiation pressure force (transverse components,
!                 in additional to longitudinal component computed
!                 previously);
!                 Torque due to scattering and absorption of radiation
!                 as measured by "vector torque cross section".
!                 Changes made in:
!                 DDSCAT
!                 GETFML (including argument list)
!                 SCAT (including argument list)
!                 VERSION (CSTAMP is now CHARACTER*22, so that a date
!                    can be attached to version number)
! 95.06.15 (BTD): Added CXE to argument list of SCAT
!                 Added new control parameters
!                 CMDSOL (to choose solution method)
!                 CMDTRQ (to choose whether or not to compute torques)
! 95.06.20 (PJF+: Major modification of DDSCAT and GETFML to support
!           BTD): use of modular iterative solvers (e.g., CCGPACK and
!                 PIM routines).
!               : add CMDSOL to argument list to specify method
!                 of solution, with "supported" options
!                 PETRKP = Petravic&Kuo-Petravic CCG routine (from
!                          CCGPACK, coded by P.J. Flatau based on
!                          earlier code by B.T. Draine)
!                 PBCGST = Preconditioned Biconjugate Gradient with
!                          Stabilization (from PIM package, coded by
!                          R. Dias da Cunha and T. Hopkins)
!               : add CMDTRQ to argument list of REAPAR and GETFML
! 95.07.21 (BTD): Major revision to subroutine SCAT
!                 Added additional scratch arrays to DDSCAT,GETFML,SCAT
! 95.07.24 (BTD): Corrected errors in SCAT
! 95.08.14 (BTD): Modified GETFML and PROGRESS to include
!                 COMMON/NORMRES/ERRSCAL to allow renormalization of
!                 error |Ax-b| relative to |b|.
! 96.01.05 (BTD): Version 5a.2:
!                 Differs from version 5a.1 in including
!                 shape options BLOCKS and DW1996, including
!                 automatic calculation of principal axes A1 and A2
!                 with largest and second-largest moment of inertia
!                 eigenvalue and eigenvector problem is solved using
!                 routines BALANC, BALBAK, EIGV, ELMHS0, ELTRN0, HQR2,
!                 IPMPAR, and SPMPAR taken from the Naval Surface
!                 Warfare Center Library of Mathematics Subroutines.
! 96.01.25 (BTD): Added shape option TWOSPH (routine TAR2SP, and
!                 modifications to TARGET,REAPAR).
! 96.02.23 (BTD): Subroutine PRINAXIS now carries out calculation of
!                 principal axes and eigenvalues of moment of inertia
!                 for arbitary grain shape (if invoked by appropriate
!                 "target generation" subroutine), and prints out
!                 3 principal axes and eigenvalues
! 96.10.17 (BTD): Version 5a.3:
!                 Added NAT0 to argument list for routine ALPHA.
!                 This required modification of ALPHA and GETFML.
! 96.10.18 (BTD): Changed NEWTMP to GPFAFT (for Generalized Prime
!                 Factor Algorithm for Fourier Transform).
! 96.11.06 (BTD): Version 5a.4:
!                 Removed code from DDSCAT and moved to subroutines
!                 GETMUELLER and WRITESCA
!                 Changed definition of scattering angle phi.
!                 Previously, phi was measured from plane containing
!                 incident k vector and Re(CXE01)
!                 Henceforth, phi is measured from Lab x,y plane
!                 (incident radiation is along x axis).
!                 This involved changes to SCAVEC
!                 Modified GETMUELLER so that for IORTH=2 it
!                 automatically computes the complete Mueller scattering
!                 matrix for each target orientation, and for the
!                 orientational average.
!                 Modified WRITESCA so that it writes out selected
!                 elements of the Mueller matrix.
! 96.11.14 (PJF): Add timers(mxtimers), ntimers to DDSCAT and getfml
!                 timers contains information about CPU times of
!                 various parts of the code and informations needed
!                 for timming (e.g. number of iterations)
!                 Further changes of "WRITESCA" to include binary
!                 file option. This should be of use, for example,
!                 for users of IDL.
!                 Added "cbinflag" to DDSCA, reapar. Added
!                 "cbinfile", "iobin" and "open" unformatted file in DDS
!                 "close(iobin) in DDSCAT.
!                 Removed "getset" routine. Hardwired "ioerr" in
!                 "errmsg", and "ioout" in "wrimsg".
!                 Combined all "rotation" code in rot.2
!                 Added dd.pro package for DDSCAT output postprocessing.
!                 The package includes IDL reading code of binary
!                 DDSCAT file and Bohren and Huffman Mie routine coded i
!                 Added 10,000 lines of LAPACK code to solve 3x3
!                 eigenvalue problem. Removed NSWC l(eispack)
!                 driver because (a) it contains ald r1mach
!                 tedious to support code (b) because LAPACK is
!                 more modern.
!                 Add several *.ps files to doc and BibTeX
!                 bibliography of recent DDA references.
! 96.11.15 (PJF): Added NetCDF binary option to "writesca.f"
!                 One needs "include" in writesca.f which is non-standar
!                 in FORTRAN77 but will probably work on all architectur
!                 Also -lnetcdf flag is needed during the loading step
!                 if NetCDF is being used. Otherwise  "empty" NetCDF rou
!                 should be substittuted to satisfy loader (c.f. dummyne
!                 "ncclose" is executed from DDSCAT.f instead of writesc
!                 NetCDF is available from http://www.unidata.ucar.edu/
!                 for most popular computers.
!                 Added "writebin.f" and "writecdf.f" for 5a7
! 96.11.21 (BTD): moved NCCLOS to WRITENET so that all dependence on
!                 netCDF is through subroutine WRITENET
!                 This required changes to DDSCAT,WRITESCA, and
!                 WRITENET, and adding additional call to WRITESCA
!                 from DDSCAT, and additional call to WRITENET from
!                 WRITESCA
!                 added variable IDVOUT2 to COMMON/NORMERR to permit
!                 device number for standard output to be passed to
!                 subroutine PROGRESS
! 96.12.02 (BTD): Corrected error in DDSCAT.f : SMORI was not
!                 initialized, and hence was incorrect except (possibly)
!                 for first wavelength and radius.
!                 Also eliminated some superfluous variables from
!                 DDSCAT.f
! 97.04.28 (BTD): Corrected bug in TARHEX: NAT was not initialized to 0
! 97.06.04 (BTD): Corrected bug in TARGET, affecting only target option
!                 "UNICYL" (target composition was not being set
!                 correctly).
! 97.07.24 (BTD): Corrected bug in GETMUELLER: was using PHI rather
!                 than PHIN in computation of amplitude scattering
!                 matrix and Mueller matrix from f_ml computed by
!                 GETFML.  Required change in GETMUELLER argument,
!                 list, with corresponding change required in DDSCAT.f
!                 Also corrected bugs in GETMUELLER pertaining to
!                 computation of Mueller matrix elements S_14 and
!                 S_31.  Note that because S_31 was incorrect, PPOL
!                 (evaluated in WRITESCA) was numerically incorrect.

! --------------- begin creation of version 6.0 -----------------------
!                 continuing to support "public" version 5a7 while
!                 carrying on private development of version 6.0
! ---------------------------------------------------------------------
! 98.01.01 (BTD): Version 6.0, with capability of using noncubic
!                 rectangular lattices.  Extensive changes required,
!                 because polarizability tensor will now be nondiagonal
!                 even when dielectric tensor is isotropic.
!                 Changes include following modules:
!                    DDSCAT
!                    ALPHA
!                    CMATVEC (in cprod.f)
!                    CPROD
!                    ESELF
!                    EVALA
!                    EVALE
!                    EVALQ
!                    GETFML
!                    MATVEC (in cprod.f)
!                    PRINAXIS
!                    REASHP
!                    SCAT
!                    TARELL
!                    TARGET
!                    TARBLOCKS
!                    TARCEL
!                    TARCYL
!                    TARHEX
!                    TARREC
!                    TARTET
!                    TAR2EL
!                    TAR3EL
!                    TAR2SP
!                 and additional module
!                    NONCUBIC (called by subroutine ALPHA)
! 98.01.20 (BTD): Modified subroutine ALPHA to use (and print out)
!                 constant C4 in evaluation of first-order corrections
!                 to alpha.
! 98.03.10 (BTD): Changes to:
!                    ERRMSG
!                    REAPAR
!                    TAR2EL
!                    TAR2SP
!                    TAR3EL
!                    TARBLOCKS
!                    TARCEL
!                    TARELL
!                    TARHEX
!                    TARTET
! 98.04.27 (BTD): Minor changes made in
!                    ALPHA
!                    CGCOMMON
!                    EVALQ
!                    GETFML
!                    NONCUBIC
!                    PRINAXIS
!                    RESTORE
!                    TARCYL
!                    TARGET
! 98.10.07 (BTD) Rewrote cgcommon.f, function SCNRM2,
!                because g77 optimization generates bad code
!                complicated code replaced with simple code
! -------------- [separate creation of version 5a9 at this time] ----
! 98.10.08 (BTD) modified cgcommon.f, subroutine PROGRESS,
!                to keep track of minimum error, and to
!                note when new error minimum is attained.  This is
!                to observe convergence behavior of PBCGST.
! 98.10.26 (BTD) modified subroutine PROGRESS in cgcommon.f
!                to also print information on
!                rate of convergence (variable RATE).
! 98.10.27 (BTD) modified subroutine SMACHCONS in cgcommon.f
!                because Solaris compiler complained about value of
!                OVERFLOW parameter.  Reduced by factor 2
! 98.11.16 (BTD) Corrected error in value of CXE01R passed to WRITEBIN
!                for binary write option.
! 98.12.07 (BTD) Found and corrected bug reported by Timo Nousiainen
!                (tpnousia@lumi.meteo.helsinki.fi) which caused Mueller
!                matrix elements S_12 and S_13 to be computed
!                due to inconsistency in definition of f_ml.  Removed
!                code in GETFML which had been inserted to change
!                f_ml to refer to "usual" incident pol states l=1,2
!                whereas GETFML assumed that incident pol states
!                l=1,2 were those specified by the user in ddscat.par
!                Stick with this latter convention.
!                Modified GETMUELLER so that Mueller matrix elements
!                are computed correctly even when user-specified
!                incident polarization state CXE01 is not linearly
!                polarized.
! 98.12.07 (BTD) Edited UserGuide to give more detailed (and now correct
!                discussion of computation of Mueller matrix elements
!                from the f_ml.
! 98.12.07 (BTD) Inserted (commented out) code into getfml to allow
!                verification of final error in approximate solution
!                to confirm reliability of error reported by iteration
!                machinery.  After confirming that these errors appear
!                to be reasonable, code was commented out, but can be
!                reactivated at future date if repetition of such tests
!                is desired.
! 98.12.16 (BTD) In getfml.f, increase upper limit on iterations from
!                300 to 10000
! 99.01.26 (BTD) Modified DDSCAT.f and reapar.f to allow user to specify
!                first IWAV,IRAD,IORI (normally 0 0 0).  This allows
!                one to restart calculational runs which abort after
!                doing only some of the desired wavelengths, radii,
!                and orientations.  Note: structure of ddscat.par is
!                changed as a result!
! 99.03.01 (BTD) rewrote subroutine orient.f to make it handle odd
!                NTHETA as described in the UserGuide
! 99.03.05 (BTD) corrected errors in new orient.f
! 99.04.30 (BTD) corrected new errors in orient.f (introduced in
!                99.03.05 revision)
! 99.06.30 (BTD) modifications to reapar (was not reading correct
!                number of SHPAR values for case TWOSPH) and tar2sp
!                in order to deal correctly with two spheroid case.
!                (Although have not yet implemented changes necessary
!                to do TWOSPH with noncubic lattice.)
!-----------------------------------------------------------------------
! 00.06.12 (BTD) [Separate release of version 5a10 at this time]
!-----------------------------------------------------------------------
! 00.06.12 (BTD) New target option NSPHER supported by routine TARNSP
!                Replaced old LAPACK source code with latest version
!                of routines from www.netlib.org
!                New LAPACK routines required by PRINAXIS are now
!                collected in two files, lapackblas.f and lapacksubs.f
!                with new LAPACK routines, PRINAXIS now executes properl
!                when compiled with g77
!                To support TARNSP, it was necessary to modify argument
!                lists for REAPAR and TARGET to include character
!                variable CFLSHP to pass name of file containing
!                description of multisphere target
!                Minor "tidying" of a number of files to make
!                all data type conversion explicit.
!                All variables in CXFFT3 now declared explicitly.
! 01.04.21 (BTD) copied over changes made to REASHP in version 5a10
!                modified all target routines so that in all cases
!                the target.out file created when using CALLTARGET
!                will contain composition information, so that there
!                is now a standard format for shape.dat input files
!                Modified as-delivered REASHP to read ICOMP.
! 01.04.21 (BTD) copied over changes made to TARNSP to correct bug
!                in version 5a10
! 02.02.12 (BTD) added support for target option PRISM3, with
!                new routine TARPRSM and concomitant changes to
!                TARGET and REAPAR
! 02.11.15 (BTD) working with Matthew Collinge to add MPI capability to
!                code.
!                DDSCAT.f now calls
!                       MPI_INIT
!                       MPI_COMM_RANK
!                       MPI_COMM_SIZE
!                       MPI_FINALIZE
!                created new module mpi_subs.f with following routines:
!                    SUBROUTINE COLSUM, which calls MPI_REDUCE
!                    SUBROUTINE SHARE1, which calls MPI_BCAST
!                    SUBROUTINE SHARE2, which calls MPI_BCAST
!                for use on systems where MPI is not installed,
!                created new module mpi_fake.f with dummy routines
!                    SUBROUTINE MPI_INIT
!                    SUBROUTINE MPI_COMM_RANK
!                    SUBROUTINE MPI_COMM_SIZE
!                    SUBROUTINE MPI_BCAST
!                    SUBROUTINE MPI_REDUCE
!                    SUBROUTINE MPI_FINALIZE
!                so that calls to these subroutines can remain in
!                the code.
! 03.02.13 (BTD) added 1 more digit of precision to qtable and qtable2
!                outputs
! 03.04.13 (BTD) working with Matthew Collinge to add MPI capability to
!                added character variable CFFLOG to contain name of
!                running output file
!                added routine NAMID to generate a unique name for outpu
!                file printed by each MPI-spawned process.
!                output to logfile is no longer unbuffered (if required
!                for debugging purposes, this can be reactivated -- see
!                comments in DDSCAT.f
!                added variables ITNUM(2) and MXITER to argument lists
!                for GETFLM and WRITESCA
! 03.07.13 (BTD) added new variable ITNUMMX(2) to DDSCAT so that
!                waarbbori.avg file will contain *maximum* number of
!                iterations required by any of the orientations being
!                averaged over
! 03.08.01 (BTD) version 6.0 released
! 03.10.24 (BTD) Modifications for version 6.1:
!                Changed method used for choosing scattering
!                directions for calculation of radiation force and
!                torque.
!                * Eliminate use of ICTHM and IPHIM -- these are
!                  removed from many argument lists, and are no longer
!                  read from file ddscat.par
!                * Add variable ETASCA to control angular resolution
!                  used for angular averaging.  ETASCA is now input
!                  from file ddscat.par, and passed to subroutines
!                  GETFML, SCAT, and WRITESCA
!                * Add variable NAVG to record number of scattering
!                  directions used.  Pass NAVG from SCAT to GETFML to
!                  DDSCAT to WRITESCA to WRITENET and WRITEBIN
!                * Add code to calculate <cos^2> in SCAT, and pass
!                  variable CSCAG2(2) from SCAT to DDSCAT to WRITESCA
!                * Modify output format of qtable, wxxrxxkxxx.sca,
!                  and wxxrxxori.avg to report NAVG
! 03.10.30 (BTD) corrected typo in WRITENET
! 04.02 22 (BTD):Modifications to writenet:
!                 * add QEXSUM,QABSUM,QSCSUM,QBKSUM,QPHSUM,QSCGSUM,
!                       QTRQABSUM,QTRQSCSUM
!                   to argument list
!                 * add various quantities (e.g. Qext) averaged over
!                   polarization to netcdf output when IORI>0
!                 * add various quantities (e.g. Qext) averaged over
!                   orientation (and polarization)
!                   to netcdf output when IORI=0
!                 * reorganize netCDF output following experimental
!                   version ver6.01
! 04.02.22 (BTD) added IWRKSC to argument list of WRITESCA and use to
!                control writing to ascii output files
! 04.02.25 (BTD) corrected typo in TARPRSM
! 04.03.31 (BTD) modified LAPACK routine SLAMC1 so that it behaves
!                as intended even after optimization by pgf77.  Note
!                that fix presumes that arithmetic is base 2.
! 04.04.01 (BTD) * added calls to TARGSPHER in subroutine TARGET
!                * added DX to argument list of TARGSPHER
!                * added arguments NPY,NPZ to support periodic b.c.
!                  option (at present time, used only for target option
!                  NSPPBC)
!                * added code to REAPAR to support option NSPPBC
! 04.04.04 (BTD) * added NPY,NPZ to argument list of TARGET
!                   and in call to TARGET from DDSCAT
!                * added variables to argument list of fake version of
!                  WRITENET in dummywritenet.f
!                * reorganize file mpi_subs.f so that mpi_subs includes
!                   SUBROUTINE COLSUM
!                   SUBROUTINE SHARE1
!                   SUBROUTINE SHARE2
!                  and mpi_fake.f includes dummy versions of
!                   SUBROUTINE COLSUM
!                   SUBROUTINE SHARE1
!                   SUBROUTINE SHARE2
!                   SUBROUTINE MPI_INIT
!                   SUBROUTINE MPI_COMM_RANK
!                   SUBROUTINE MPI_COMM_SIZE
!                   SUBROUTINE MPI_FINALIZE
!                  Subroutine mpi_subs is no longer used in "plain"
!                  version -- all the ersatz MPI calls used by DDSCAT
!                  are supplied in mpi_fake.f
! 04.04.29 (BTD) Added target option ANIREC (homogeneous anisotropic
!                rectangular target
!                To accomplish this, modified subroutines
!                  REAPAR
!                  TARGET
!                and added new routine
!                  TARANIREC
! 04.05.21 (BTD) Fixed bug in DDSCAT.f: one of the calls to WRITESCA
!                had QSCGSUM and QSCG2SUM in incorrect order.
! 04.05.22 (BTD) Fixed bug in computation of true error in getfml.f
! 04.05.23 (BTD) Added CWHERE string in reapar.f to provide helpful
!                information to user in event of malformed ddscat.par
!                file

! Begin DDSCAT 6.2:
! 04.09.14 (BTD) Added THETADF,PHIDF,BETADF to argument list of
!                subroutine TARGET
!                define new local variable R(3,3),RI(3,3),
!                COSBE,COSPH,COSTH,SINBE,SINPH,SINTH,...
!                add code to recalculate CXALPH(1-3) and CXALOF(1-3)
!                if there is nonzero rotation of dielectric frame
!                at local lattice site, AND dielectric tensor is
!                anisotropic.
! 04.09.14 (BTD) Added CSHAPE,BETADF,PHIDF,THETADF to argument list
!                of subroutine REASHP
!                Add support for CSHAPE='ANIFIL' to read dipole
!                locations, composition, and orientation angles
!                BETADF, PHIDF, THETADF for Dielectric Frame relative
!                to Target Frame at each occupied site.
!                Add new shape option CSHAPE=NANSPH to generate
!                target from list of spheres with 3 composition
!                indices plus 3 angles giving orientation of Dielectric
!                Frame relative to Target Frame
! 04.10.13 (BTD) Modified REAPAR to support new shape option NANSPH
!                Modified DIELEC to check that user does not specify
!                refractive index = 1
! 04.10.14 (BTD) Modified REAPAR to support reading wavelengths from
!                file "wave.tab".
! 04.12.29 (BTD) Corrected error in WRITESCA: when IWRKSC=1, had not
!                been writing out correct values of Mueller matrix to
!                wxxrxxkxxx.sca files. (Output error only: wrong
!                variable in WRITE statement. No computational error
!                involved.)
! 05.03.19 (BTD) REAPAR modified:
!                Added 2 lines to set CLFEPS(1) in case using options
!                H2OLIQ or H2OICE (previously were left undefined,
!                which cause problems with output statement)
! 05.03.29 (BTD) ROTATE modified:
!                Added code to guard against roundoff errors leading to
!                argument of ACOS with abs(arg)>1
! 05.07.08 (BTD) Corrected error in routine eself related to option
!                NSPPBC
! 05.08.03 (BTD) Added code to REAPAR and TARGET to support option
!                RCTPBC = rectangular brick with periodic b.c.
! 05.08.04 (BTD) Corrected error in ESELF
! 05.08.04 (BTD) Modified modules
!                   DDSCAT.f
!                   reapar.f
!                   mpisubs.f
!                   mpi_fake.f
!                Added new character variable CMDFRM to be read by
!                subroutine REAPAR from input file ddscat.par.
!                CMDFRM=LFRAME : angles THETAN,PHIN are relative to
!                                Lab Frame (xlab,ylab,zlab)
!                CMDFRM=TFRAME : angles THETAN,PHIN are relative to
!                                Target Frame (a1,a2,a3)
!                Modified DDSCAT.f to define ENSC,EM1,EM2 accordingly
! 05.09.26 (BTD) Corrected bug in DDSCAT regarding calculation of
!                scattering intensities for option CMDFRM = TFRAME
! 05.10.11 (BTD) Modified to support up to 1000 wavelengths:
!                changed MXWAV from 100 to 1000 in DDSCAT
!                changes to DDSCAT,NAMER,WRITESCA:
!                changed CFLAVG*13 to CFLAVG*14
!                changed CFLSCA*14 to CFLSCA*15
!                also other changes to NAMER
! 05.10.11 (BTD) Added CMDFRM to argument list of WRITESCA
! 05.10.18 (BTD) Modified DDSCAT.f to allow use of more than 1000
!                orientations if IWRKSC=0 (previously would exceed array
!                limits in subroutine NAMER with IORI>1000).
!                Modified REAPAR to ensure that user does not request
!                more than 1000 orientations with IWRKSC=1
! 06.04.10 (BTD) * Modified REAPAR to read IWRPOL from ddscat.par
!                  and pass to DDSCAT
!                * Created routine WRITEPOL to write out complex
!                  polarization array as well as additional information
!                  defining target geometry and incident wave
!                * Modified DDSCAT to use IWRPOL and (if IWRPOL=1) to
!                  call subroutine WRITEPOL to write out unformatted
!                  file of reduced polarization array as well as
!                  additional information defining target geometry and
!                  incident wave
!                * Modified NAMER to generate unique names
!                  waaarbbkccc.pol1 and waaarbbkccc.pol2 for
!                  polarization arrays.  Rewrote NAMER to remove
!                  EQUIVALENCE statements and use another method to
!                  construct strings.
! 06.09.13 (BTD) Modified REAPAR and TARGET to handle new target option
!                CYLCAP
!                Added new routine TARCYLCAP to generate homogeneous
!                cylinder with hemispherical caps.
! 06.09.15 (BTD) * Modified TARHEX to allow hex prism orientation in TF
!                  to be selected.
!                  TARHEX now able to produce hexagonal prism in 6
!                  different orientations in TF.
!                * Modified TARGET and REAPAR to support this change.
! 06.09.15 (BTD) Modified TARGET and REAPAR to support new
!                target option HEXPBC (uses routine TARHEX).
! 06.09.15 (BTD) * Modified TARCYL to allow cylinder orientation in TF
!                  to be selected
!                * Modified TARGET and REAPAR to support this change.
! 06.09.21 (BTD) Modified DDSCAT and WRITEPOL to add PYD,PZD,DX to
!                quantities written by WRITEPOL
! 06.09.28 (BTD) *** version 6.2.3 ***
!                Modified
!                * ESELF to do direct calculation of A matrix elements
!                  including contribution form replica dipoles using
!                  new routine DIRECT_CALC
!                  (cannot employ symmetry tricks
!                  because they don't apply when kx or kz .ne.0)
! 06.09.29 (BTD) Modified
!                * ESELF to change way calculation is done when IPBC=1
!                  - calculate and store full Fourier transform of A
!                    matrix in CXZC
!                  - added variable IPBC to argument list to increase
!                    dimensioning of CXZC when PBC is used.
!                  - store full CXZC when IPBC=1
!                * DDSCAT to change dimensioning of IPBC when new
!                  parameter MXPBC=1
!                * GETFML
! 06.10.05 (BTD) Modified
!                * REAPAR to improve handling of input for PBC targets,
!                  specifically selection of scattering directions
!                * PBCSCAVEC = new routine added to calculate
!                  scattering directions for PBC targets
! 06.10.05 (BTD) fixed minor bugs in WRITENET
! 06.10.07 (BTD) Extend to handle computation of S matrix for 1-d
!                and 2-d periodic targets.
!                Modified:
!                * DDSCAT
!                * GETMUELLER
!                * PBCSCAVEC
!                * WRITESCA (added variables to arg list, but no
!                  change in operation)
! 06.10.23 (BTD) Set a1=(1,0,0),a2=(0,1,0) for target options
!                RCTPBC,CYLPBC,HEXPBC
! 06.11.28 (BTD) Modified WRITESCA to correct error in calculation of
!                degree of linear polarization of scattered light
! 06.11.29 (BTD) Modified REAPAR:
!                added code to catch error in input CMDFRM
! 06.12.24 (BTD) Modified DDSCAT
!                Added ORDERM,ORDERN to argument list of GETMUELLER
!                Modified GETMUELLER
!                * Added ORDERM,ORDERN to argument list.
!                * When JPBC=3,(target periodic in two dimensions) use
!                  ORDERM,ORDERN to identify (0,0) forward scattering
!                  direction for special-case calculation of Mueller
!                  matrix
!                Modified EXTEND
!                * Added 1 as first element of NF235, so that we can
!                  treat 1d targets that are 1 layer thick, or
!                  plane-parallel 2d target using single dipole line.
! 06.12.28 (BTD) Modified DDSCAT
!                * Added arrays A3,XLR,YLR,ZLR
!                * Added XLR,YLR,ZLR to argument list of PBCSCAVEC
!                Modified PBCSCAVEC
!                * corrected error in calculation of unit vector for
!                  scattered wave in special case
!                * added XLR,YLR,ZLR to argument list
!                  use XLR,YLR,ZLR in defining unit pol vectors in speci
!                  case of scattering angle theta = 0 or pi (in which
!                  case scattering plane is not uniquely defined).
! 07.01.18 (BTD) Modified DDSCAT,SHARE1
!                * added A3(3),IWRPOL,IPBC,JPBC,MXNX,MXNY,MXNZ,MXPBC,
!                  ORDERM(MXSCA),ORDERN(MXSCA),
!                  PYD,PYDDX,PZD,PZDDX to argument list of SHARE1
!                * added appropriate MPI_BCAST... calls in SHARE1
! 07.01.19 (BTD) * Corrected bug in PBCSCAVEC
! 07.01.20 (BTD) Modified REAPAR,TARGET
!                * added support for target option SLBPBX
! 07.02.23 (BTD) Modified REAPAR,TARGET
!                * added support for target option SLBPBC
!                Create new routine TARSLBPBC for 4-layer slab of
!                infinite extent in target y and z directions
!                (using periodic boundary conditions)
! 07.06.21 (BTD) *** Version 7.0.2.  
!                Modified to handle vector X0(3) definin
!                physical location in TF of lattice site (0,0,0).
!                This facilitates subsequent use of DDfield for
!                near-field calculations (simplifies specification of
!                locations).
!                Also modified to change way azimuthal scattering angle
!                is defined for targets with 1-d periodicity.
!                Modified:
!                * TARPBX
!                * TARCYL
!                * TARGET
!                * EXTEND
!                * GETFML
!                * SCAT
!                * EVALE
!                * WRITEPOL
!                * READPOL
!                * MPI_FAKE
!                * MPI_SUBS
!                * PBCSCAVEC
!                * DDSCAT
!                * DDfield
! 07.06.22 (BTD) * SCAT: corrected error in computation of PPOL (returne
!                  to formula used prior to 06.11.28)
!                * GETMUELLER: corrected multiplicative error
!                  in calculation of Mueller matrix elements for
!                  targets with 1-d periodicity
!                * PBCSCAVEC: changed definition of azimuthal angle
!                  zeta for scattered directions -- changed so that
!                  forward scattering has zeta=0
! 07.06.30 (BTD) DDSCAT, DIAGL (in cprod.f), MATVEC (in cprod.f),
!                CMATVEC (in cprod.f):
!                * moved CMDFFT from COMMON/M6/... CMDFFT
!                  to COMMON/M8/CMDFFT
!                * moved PYD,PZD to beginning of COMMON/M6/
!                ALPHA:
!                * eliminated COMMON/M6/
!                * added COMMON/M8/CMDFFT
! 07.07.03 (BTD) modified
!                * DDSCAT
!                * GETMUELLER
!                * PBCSCAVEC
! 07.07.08 (BTD) modified
!                * DDSCAT: added EM1,EM2 to argument list for PBCSCAVEC
!                  and GETMUELLER
!                * PBCSCAVEC: added EM1,EM2 to argument list
!                  added code to compute EM1,EM2 when JPBC > 0
!                * GETMUELLER: added EM1, EM2 to argument list
!                  use EM1,EM2 to compute S matrix in case of forward
!                  scattering
! 07.08.03 (PJF+BTD) **** Version 7.0.3
!                * PJF converted bulk of code from f77 to f90 using
!                  NAGware
!                * DDPRECISION: module to specify whether DDSCAT uses
!                  single or double precision storage and arithmetic
!                  DDPRECISION module is used by PROGRAM DDSCAT and
!                  all subroutines.
!                * TIMEIT: now uses f95 intrinsic function CPU_TIME
!                * PRINAXIS: modified to use DSYEVJ3 to find 
!                  eigenvalues and eigenvectors of moment of inertia
!                  tensor
!                * LAPACK subroutines have been removed from the code
! 07.08.04 (BTD) * DDSCAT
!                  * replaced COMMON/M1/ with USE MODULE DDCOMMON_1
!                  * replaced COMMON/M2/ with USE MODULE DDCOMMON_2
!                    add dynamic allocation of CXADIA
!                  * replaced COMMON/M3/ with USE MODULE DDCOMMON_3
!                    add dynamic allocation of CXZC
!                  * replaced COMMON/M4/ with USE MODULE DDCOMMON_4
!                    add dynamic allocation of CXZW
!                  * replaced COMMON/M5/ with USE MODULE DDCOMMON_5
!                    add dynamic allocation of IOCC
!                  * replaced COMMON/M6/ with USE MODULE DDCOMMON_6
!                  * replaced COMMON/M7/ with USE MODULE DDCOMMON_7
!                    add dynamic allocation of CXAOFF
!                  * replaced COMMON/M8/ with USE MODULE DDCOMMON_8
!                * ALPHA (in alpha.f90)
!                  * replaced COMMON/M8/ with USE MODULE DDCOMMON_8
!                * DIAGL (in cprod.f90)
!                  * replaced COMMON/M2/ with USE MODULE DDCOMMON_2
!                  * replaced COMMON/M6/ with USE MODULE DDCOMMON_6
!                  * removed COMMON/M8/ (was not used)
!                * CMATVEC (in cprod.f90)
!                  * renamed AK(3)->AKR(3) for consistency with
!                    COMMON/M1/ in DDSCAT
!                  * replaced COMMON/M1/ with USE MODULE DDCOMMON_1
!                  * replaced COMMON/M2/ with USE MODULE DDCOMMON_2
!                  * replaced COMMON/M3/ with USE MODULE DDCOMMON_3
!                  * replaced COMMON/M4/ with USE MODULE DDCOMMON_4
!                  * replaced COMMON/M5/ with USE MODULE DDCOMMON_5
!                  * replaced COMMON/M6/ with USE MODULE DDCOMMON_6
!                  * replaced COMMON/M7/ with USE MODULE DDCOMMON_7
!                  * replaced COMMON/M8/ with USE MODULE DDCOMMON_8
!                * MATVEC (in cprod.f90)
!                  * renamed AK(3)->AKR(3) for consistency with COMMON/M1/
!                  * in DDSCAT
!                  * replaced COMMON/M1/ with USE MODULE DDCOMMON_1
!                  * replaced COMMON/M2/ with USE MODULE DDCOMMON_2
!                  * replaced COMMON/M3/ with USE MODULE DDCOMMON_3
!                  * replaced COMMON/M4/ with USE MODULE DDCOMMON_4
!                  * replaced COMMON/M5/ with USE MODULE DDCOMMON_5
!                  * replaced COMMON/M6/ with USE MODULE DDCOMMON_6
!                  * replaced COMMON/M7/ with USE MODULE DDCOMMON_7
!                  * replaced COMMON/M8/ with USE MODULE DDCOMMON_8
!                * PROGRESS (in cgcommon.f90)
!                  * replaced COMMON/NORMERR/ with USE MODULE DDCOMMON_9
!                * GETFML (in getfml.f90)
!                  * replaced COMMON/NORMERR/ with USE MODULE DDCOMMON_9
! 07.08.05 (BTD) * DDSCAT
!                  * extensive modification to use dynamic memory allocation
!                  * add MXNX,MXNY,MXNZ to CALL REAPAR argument list in
!                    use this MXNX,MXNY,MXNZ on first call to TARGET
!                    after actual size of extended target is set by EXTEND,
!                    then reallocate minimum memory required for problem
!                * REAPAR (in reapar.f90)
!                  * modified to read MXNX,MXNY,MXNZ from ddscat.par
!                  * add MXNX,MXNY,MXNZ to argument list
!                * DMACHCONS (in cgcommon.f90)
!                  * modified setting of MACHEPS,UNDERFLOW, and OVERFLOW
!                    based on whether single- or double-precision
!                    version is used 
!                    Question: do changes slow execution??
! 07.10.27 (BTD) * REAPAR
!                  * allow new target options PB2PBC, SLBLIN, SLNPBC
!                * TARGET
!                  * support new target options PB2PBC, SLBBOX, SLNPBC
! 08.01.12 (BTD) * DDSCAT
!                  * add PYDDX,PZDDX to argument list of WRITESCA
! 08.01.13 (BTD) * DDSCAT
!                  * introduce variable NRWORD = length (bytes) of real word
!                  * add to argument list of WRITEPOL, so that this can
!                    be stored in file written by WRITEPOL for sanity check
!                * WRITEPOL
!                  * provide option for more efficient storage when
!                    BETADF,THETADF,PHIDF are trivial because target
!                    material is either (1) isotropic or (2) has
!                    optical axes aligned with x_TF,y_TF,z_TF
! 08.01.17 (BTD) * DDSCAT
!                  * add new variable IANISO to record whether target
!                    is isotropic (IANISO=0), anisotropic but with
!                    principal axes everywhere aligned with x_TF,y_TF,z_TF
!                    (IANISO=1) or generally anisotropic (IANISO=2)
!                * TARGET
!                  * modify TARGET to set value of IANISO
!                * WRITEPOL
!                  * modifications to use IANISO to determine whether
!                    anisotropy info in BETADF,THETADF,PHIDF needs to
!                    be written out to file
! 08.01.21 (BTD) * DDSCAT
!                  * add IANISO to argument list of SHARE1
! 08.02.01 (BTD) * DDSCAT
!                  * changed SHPAR(10) to SHPAR(12)
!                * REAPAR
!                  * changed SHPAR(10) to SHPAR(12)
!                  * added support for target options 
!                    BISLINPBC
!                    RCTGLBLK3
!                    TRILYRPBC
!                * TARGET
!                  * changed SHPAR(10) to SHPAR(12)
!                    BISLINPBC
!                    RCTGLBLK3
!                    TRILYRLKPBC
!                * added new routine TARRCTBLK3 to support
!                  target shapes BISLINPBC, RCTGLBLK3, TRILYRLKPBC
! 08.02.09 (BTD) Cosmetic changes (affecting appearance of running output)
!                in following routines:
!                * DDSCAT
!                * REAPAR
!                * DSYEVJ3
! 08.02.13 (BTD) * Added MXNAT0 to argument list of SHARE1
!                * A number of changes in SHARE1 to correct dimensioning
!                  of IOCC,IXYZ0,BETADF,THETADF,PHIDF from MXNAT to MXNAT0
!                * Added CALL_MPI_BCAST_INT calls for 
!                  * MXBETA
!                  * MXNAT
!                  * MXNAT0
!                  * MXPHI
!                  * MXSCA
!                  * MXTHET
!                  * MXWAV
! 08.02.17 (BTD) * DDSCAT: allocation of all arrays needed
!                  by parallel processes now allocated by each process
!                  rather than just master process
!                * DDSCAT: call new routine SHARE0 to distribute information
!                  needed for memory allocation
!                * DDSCAT: added CLOSE(11) to close mtable
!                * Created new routine SHARE0
!                  in mpi_subs.f90 and mpi_fake.f90
!                * Added DAEFF to argument list of SHARE1
!                DDSCAT 7.0.4 now appears to run properly under MPI !
! 08.03.11 (BTD) * reapar.f90: corrected typo BLSLINPBC -> BISLINPBC
! 08.04.19 (BTD) change in notation: ALPHA -> GAMMA
!                This affects routines
!                * DDSCAT
!                * REAPAR (changed order in argument list)
!                * GETFML (changed order in argument list)
!                * SHARE1 in mpi_subs.f (changed order in argument list)
!                * SHARE1 in mpi_fake.f (changed order in argument list)
!                * DIAGL (DDCOMMON_6)
!                * MATVEC (DDCOMMON_6)
!                * CMATVEC (DDCOMMON_6)
!                * CPROD (DDCOMMON_6) and call to ESELF
!                * ESELF
!                * DDCOMMON_6 (changed ALPHA -> GAMMA)
!                and program
!                * DDfield
! 08.05.12 (BTD) ver7.0.6:
!                * DDSCAT : added code mods suggested by Art Lazanoff
!                  - initialize IXYZ0 to 0            [WHY?]
!                    with !$omp parallel do / !$omp end parallel do
!                  - initialize CXZC to (0._WP,0._WP) [WHY?]
!                    with !$omp parallel do / !$omp end parallel do
!                  - initialize CXZW to (0._WP,0._WP) [WHY?]
!                    with !$omp parallel do / !$omp end parallel do
!                  - initialize CXSCR1 to (0._WP,0._WP)
!                    with !$omp parallel do / !$omp end parallel do
!                  - iniitialize SCRRS1 to 0._WP
!                    with !$omp parallel do / !$omp end parallel do
!                * CPROD  : added calls to TIMEIT to time call to ESELF
!                * dummy.f90 : add this module to hold dummy routines
!                  - DUMMY
!                  - CDUMMY
!                  - CDUMMY_1
!                * ESELF  : added code mods suggested by Art Lazanoff
!                           to support OpenMP if compiled with preprocessor
!                           flag -fpp -Deself_omp
!                * GETFML : changed declaration of DUMMY from REAL to EXTERNAL
!                * PETR   : (in ccgpack.f90)
!                  - change CDUMMY -> CDUMMY(1)
!                  - declare DUMMY,CDUMMY_1 as EXTERNAL
!                * pim.f90 package :added code mods suggested by Art Lazanoff
!                  changed PARAMETER(SPARSIZ=2) -> PARAMETER(SPARSIZ=6) in
!                  - PIMCBICG
!                  - PIMCBICGSTAB
!                  - PIMCCG
!                  - PIMCCGNE
!                  - PIMCCGNR
!                  - PIMCCHEBYSHEV
!                  - PIMCQMR
!                  - PIMCRBICGSTAB
!                  - PIMCRGCR
!                * SCAT   : added code mods suggested by Art Lazanoff
!                           to support OpenMP if compiled with preprocessor
!                           flag -fpp -Dscat_omp
!                * TARGSPHER:
!                  - added X0 to argument list
!                  - added code to calculate X0
!                * WRITESCA: following suggestion from Lazanoff, initialize
!                  DAEFF=0._WP and DPHYS=0._WP [WHY?]
!                * ZBCG2
!                  - changed WORK(1:N,J) in arguments to WORK(1,J) [10 places]
! 08.05.21 (BTD) * WRITEPOL : change dimensioning 
!                  CXPOL(MXNAT,3) -> CXPOL(NAT0,3)
! 08.05.30 (BTD) * DDSCAT : changed declarations                 
!                  SM(MXSCA,4,4)      -> SM(4,4,MXSCA)
!                  SMORI(MXSCA,4,4)   -> SMORI(4,4,MXSCA)
!                  SMORI_1(MXSCA,4,4) -> SMORI_1(4,4,MXSCA)
!                * GETMUELLER : changed declarations
!                  SM(MXSCA,4,4)    -> SM(4,4,MXSCA)
!                  SMORI(MXSCA,4,4) -> SMORI(4,4,MXSCA)
!                  and corresponding reordering of indices in lines
!                  where SM and SMORI are evaluated
!                * COLSUM in mpi_fake.f90: changed declarations
!                  SMORI(MXSCA,4,4)   -> SMORI(4,4,MXSCA)
!                  SMORI_1(MXSCA,4,4) -> SMORI_1(4,4,MXSCA)
!                  added code to use 
!                     QSCSUM_1 to update QSCSUM
!                     QABSUM_1 to update QABSUM
!                     QEXSUM_1 to update QEXSUM
!                     QBKSUM_1 to update QBKSUM
!                     QPHSUM_1 to update QPHSUM
!                     QSCG2SUM_1 to update QSCG2SUM
!                     QSCGSUM_1 to update QSCGSUM
!                     QTRQABSUM_1 to update QTRQABSUM
!                     QTRQSCSUM_1 to update QTRQSCSUM
!                     CX1121_1 to update CX1121
!                     S1111_1 to update S1111
!                     S2121_1 to update S2121
!                     SMORI_1 to update SMORI
!                * COLSUM in mpi_subs.f90: changed declarations
!                  SMORI(MXSCA,4,4)   -> SMORI(4,4,MXSCA)
!                  SMORI_1(MXSCA,4,4) -> SMORI_1(4,4,MXSCA)
!                * WRITESCA: changed declarations
!                  SM(MXSCA,4,4)    -> SM(4,4,MXSCA)
!                  SMORI(MXSCA,4,4) -> SMORI(4,4,MXSCA)
!                  and corresponding changes in code
! 08.06.07 (BTD) * CXFFT3MKL routine added (written by Art Lazanoff)
!                * REAPAR : added support for option FFTMKL
!                * ESELF : added calls to CXFFT3_MKL 
! 08.06.23 (BTD) * TARCYLCAP : fixed bug 
! 08.07.22 (BTD) * DDSCAT: added XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN to argument list
!                  of SHARE1
!                * SHARE1: added XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN to argument list
!                  of SHARE1
!                * REAPAR: added sanity check to catch users who select
!                  target option SPHERES_N but also specify multiple materials
! 08.07.23 (BTD) * DDSCAT: removed some MPI_BARRIER calls,
!                  added others just before SHARE1 and SHARE2
! 08.07.27 (BTD) * DDSCAT: 
!                  * eliminated variable CXALDS
!                    (removed from argument list of GETFML)
!                  * changed final allocation of BETADF,PHIDF,THETADF
!                    from NAT0 -> MXNAT
!                * ALPHADIAG: 
!                  * changed dimensioning of BETADF,PHIDF,THETADF
!                    from (NAT0) to (MXNAT)
!                  * eliminated variable CXALDS from argument list
!                * GETFML: 
!                  * changed dimensioning of BETADF,PHIDF,THETADF
!                    from (NAT0) to (MXNAT)
!                  * eliminated variable CXALDS from argument list
!                  * eliminated CXALDS from arg list of ALPHADIAG
! 08.08.07 (BTD) v7.0.7
!                * DDSCAT
!                  * changed #ifdef debug -> #ifdef openmp
!                  * add #ifdef mpi , #ifndef mpi to switch between 
!                    INCLUDE 'mpif.h' 
!                    and
!                    INTEGER :: MPI_COMM_WORLD
!                * eself
!                  * changed #ifdef debug -> #ifdef openmp
!                * scat
!                  * changed #ifdef debug -> #ifdef openmp                
! 08.08.08 (BTD) * REASHP corrected 
!                  FRMFIL -> FROM_FILE
!                  ANIFIL -> ANIFRMFIL
!                * TARGET: added X0 to argument lists of
!                  * TARBLOCKS
!                  * TAR2SP
!                  * TARHEX
!                  * TARNAS
!                  * TARNSP
!                * TARBLOCKS: added X0 to arg list, and code to set X0
!                * TAR2SP: added X0 to arg list, and code to set X0
!                * TARHEX: added X0 to arg list, and code to set X0
!                * TARNAS: added X0 to arg list, and code to set X0
!                * TARNSP: added X0 to arg list, and code to set X0
! 08.08.31 (BTD) * DDSCAT: eliminate netCDF-related code
!                          removed CNETFLAG, CNETFILE 
!                * REAPAR: removed CNETFLAG
!                * SHARE1: removed CNETFLAG
!                * WRITESCA: removed CNETFLAG, CNETFILE
!                * TARGET: made corrections
!                edited various routines to standardize target.out format:
!                * REASHP
!                * TAR2EL
!                * TAR2SP
!                * TAR3EL
!                * TARANIREC
!                * TARBLOCKS
!                * TARCEL
!                * TARCYL
!                * TARCYLCAP
!                * TARELL
!                * TARGSPHER
!                * TARHEX
!                * TARLYRSLAB
!                * TARNAS
!                * TARNSP
!                * TARPBXN
!                * TARPRSM
!                * TARRCTBLK3
!                * TARREC
!                * TARRECREC
!                * TARSLBLIN
!                * TARTET
! 08.09.12 (BTD) * Corrected typo DSKRCTGNL -> DSKRCTNGL in
!                  * REAPAR
!                  * TARGET
!                  * CALLTARGET
!                * Removed extra dummy read statement from REASHP
! 08.09.15 (BTD) * Corrected typo SPHERN_PBC -> SPHRN_PBC which
!                  caused SPHERN_PBC to be assigned JPBC=0 regardless of
!                  values of SHPAR(2) or SHPAR(3)
! 08.09.17 (BTD) * Corrected typos in routines TARNAS and TARNSP that caused
!                  target centroid and therefore X0 to be computed incorrectly
!                  [this only affected value of X0(1-3); it did not affect 
!                  calculations of absorption and scattering ]
! 08.09.23 (BTD) * Corrected typo in routine REAPAR
! 08.10.28 (BTD) * Corrected error in TARGET: added X0 to argument list of
!                  TARBLOCKS for target option MLTBLOCKS
! 09.07.08 (BTD) * Corrected typo in DDfield.f90 affecting evaluation of
!                  XMIN and XMINPHYS (output describing spatial extent of
!                  target)
! 09.08.12 (BTD) * Corrected declaration of CFLPOL in DDpol.f90 to agree
!                  with declaration in readpol.f90 and DDfield.f90. 
! 09.09.10 (BTD) ver7.0.8
!                * Added support in REAPAR and TARGET for new target options
!                  * FRMFILPBC
!                  * ANIFILPBC
!                * Added variable NCOMP_NEED for DDSCAT to communicate
!                  with REASHP via subroutine TARGET
!                * Added check in DDSCAT that NCOMP_NEED does not
!                  exceed NCOMP
! 09.09.17 (BTD) * Added allocation test in REASHP
!                * Added various output warnings in REASHP
! 09.12.11 (BTD) v7.0.8
!                * corrected error in REAPAR that prevented use of
!                  target options DW1996TAR and MLTBLOCKS
! 10.01.06 (BTD) ver 7.1.0
!                * modified tarnsp.f90 and tarnas.f90 to expect input file
!                  with structure
!                    N = number of spheres
!                    comment line
!                    comment line
!                    comment line
!                    comment line
!                    x1 y1 z1 a1
!                    x2 y2 z2 a2
!                       ...
!                    xN yN zN aN
!                * added sanity checks to tarnsp.f90 and tarnas.f90
!                  to report file incompatibility to user
! 10.01.24 (BTD) Corrected inconsistency between UserGuide and code
!                regarding shape parameters for target options
!                RCTGLBLK3 and TRILYRPBC.
!                * modified tarrctblk3.f90 (reordered arg. list).
! 10.01.28 (BTD) in DDSCAT:
!                * when IWRKSC=0, change IORI from 0 to 1 in call to
!                  NAMER
! 10.01.30 (BTD) in DDSCAT
!                * allow use of CMDFRM='TFRAME' for isolated targets
!                  (had been inadvertently forbidden)
!                in reashp.f90
!                 * modified to make more robust and informative
!                  regarding errors in input shape file
! 10.02.02 (BTD) added support for new target options
!                SLAB_HOLE and SLBHOLPBC
!                * added module tarslbhol.f90
!                * modified target.f90
!                * modified reapar.f90
! 10.02.06 (BTD) target_v3.f90
!                * discontinue target option TARSLBLIN
!                * correct handling of target option BISLINPBC
! 10.03.03 (BTD) target_v3.f90
!                * corrected error handling target option FRMFILPBC
! end history

! Copyright (C) 1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,
!               2005,2006,2007,2008,2009,2010 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
      RETURN
    END SUBROUTINE VERSION
