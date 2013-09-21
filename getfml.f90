      SUBROUTINE GETFML(AKR,AK3,AKSR,BETADF,PHIDF,THETADF,DX,X0,CALPHA,  &
                        CMDSOL,CMDFFT,CMDTRQ,CSHAPE,CXADIA,CXAOFF,       &
                        CXALPH,CXALOF,CXALOS,CXE,CXE01R,CXE02R,          &
                        CXEPS,CXF11,CXF12,CXF21,CXF22,CXRLOC,CXPOL,CXSC, &
                        CXSCR1,CXZC,CXZW,EM1R,EM2R,ETASCA,GAMMA,IBETH,   &
                        IBETH1,ICOMP,IDVOUT,INIT,IOCC,IORTH,IPBC,IPHI,   &
                        IPHI1,ITASK,IXYZ0,ITNUM,JPBC,MXITER,LACE,LAXI,   &
                        LCLM,LGI,LPI,LQI,LSC0,MXCOMP,MXCXSC,MXN3,MXNAT,  &
                        MXNX,MXNY,MXNZ,MXPBC,MXPHI,MXSCA,MYID,NAT,NAT0,  &
                        NAT3,NAVG,NCOMP,NSCAT,NX,NY,NZ,PHI,PIA2,QABS,    &
                        QBKSCA,QEXT,QPHA,QSCA,QSCAG,QSCAG2,QTRQAB,       &
                        QTRQSC,SCRRS1,SCRRS2,SHPAR,TOL,TIMERS,MXTIMERS,  &
                        NTIMERS,AEFFA,WAVEA,MXRAD,MXWAV,CENTER,IWRPOL,c, &
                        h_bar,h_bar2,velocity,e_charge,DielectricConst)
                        !Arguments AEFFA and after added by NWB 3/8/12

!----------------------------- v2 ----------------------------------------

      USE DDPRECISION,ONLY: WP
      USE DDCOMMON_9,ONLY: ERRSCAL,IDVOUT2,ITERMX,ITERN
      IMPLICIT NONE

! Arguments:

      CHARACTER :: CALPHA*6,CMDSOL*6,CMDFFT*6,CMDTRQ*6,CSHAPE*9

      INTEGER :: IBETH,IBETH1,IDVOUT,INIT,IORTH,IPBC,IPHI,IPHI1,ITASK, &
         JPBC,LACE,LAXI,LCLM,LGI,LPI,LSC0,LQI,MXCOMP,MXCXSC,MXN3,      &
         MXNAT,MXNX,MXNY,MXNZ,MXPBC,MXPHI,MXSCA,MXTIMERS,MYID,NAT,     &
         NAT0,NAT3,NAVG,NCOMP,NSCAT,NTIMERS,NX,NY,NZ,MXRAD,MXWAV, JJ
         !MXRAD, MXWAV, JJ added by NWB 3/8/12

      INTEGER :: MXITER,IWRPOL !IWRPOL added by NWB 7/12/12
      INTEGER :: & 
         ITNUM(2)

      INTEGER*2 ::    &
         ICOMP(MXN3), &
         IOCC(MXNAT)

      INTEGER ::      &
         IXYZ0(NAT0,3)

      REAL(WP) :: AK3,ETASCA,GAMMA,PIA2,TOL,TOLR

      REAL(WP) ::          &
         AKR(3),           &
         AKSR(3,MXSCA),    &
         BETADF(MXNAT),    &
         DX(3),            &
         EM1R(3,MXSCA),    &
         EM2R(3,MXSCA),    &
         PHI(MXPHI),       &
         PHIDF(MXNAT),     &
         QABS(2),          &
         QBKSCA(2),        &
         QEXT(2),          &
         QPHA(2),          &
         QSCA(2),          &
         QSCAG(3,2),       &
         QSCAG2(2),        &
         QTRQAB(3,2),      &
         QTRQSC(3,2),      &
         SCRRS1(MXNAT,3),  &
         SCRRS2(MXNAT),    &
         SHPAR(12),        &
         THETADF(MXNAT),   &
         TIMERS(MXTIMERS), &
         X0(3),            &
         AEFFA(MXRAD),     &
         WAVEA(MXWAV),     &
         CENTER(3),        &
         c,                &
         h_bar,            &
         h_bar2,           &
         velocity,         &
         e_charge,         &
         DielectricConst   !AEFFA and after by NWB 3/8/12

      COMPLEX(WP) ::                    &
         CXADIA(MXN3),                  &
         CXALPH(MXN3),                  &
         CXALOF(MXN3),                  &
         CXAOFF(MXN3),                  &
         CXALOS(NAT,3),                 &
         CXE(MXN3),                     &
         CXE01R(3),                     &
         CXE02R(3),                     &
         CXEPS(MXCOMP),                 &
         CXF11(MXSCA),                  &
         CXF12(MXSCA),                  &
         CXF21(MXSCA),                  &
         CXF22(MXSCA),                  &
         CXRLOC(MXCOMP+1,3,3),          &
         CXPOL(MXN3),                   &
         CXSC(MXCXSC),                  &
         CXSCR1(MXN3),                  &
         CXZC(MXNX+1+MXPBC*(MXNX-1),    &
              MXNY+1+MXPBC*(MXNY-1),    &
              MXNZ+1+MXPBC*(MXNZ-1),6), &
         CXZW(MXNAT,24)

!-----------------------------------------------------------------------

! Local variables

      COMPLEX(WP) :: &
         CXE0R(3)

      INTEGER :: I,ITER,J,JO,MXITER2,NAT03
      INTEGER :: NO_CG_RESTART, NO_CG_ITER

      INTEGER :: &
         IPAR(13)

      REAL(WP) :: CABS,CBKSCA,CEXT,CPHA,CSCA,CSCAG2,COSPHI,DTIME, &
         E02,FALB,SINPHI

      REAL(WP) ::   &
         CSCAG(3),  &
         CTRQAB(3), &
         CTRQSC(3), &
         SPAR(6)

      SAVE E02

      EXTERNAL CMATVEC,DIAGL,DUMMY,MATVEC,PCSUM,PRECOND,PROGRESS,PSCNRM2

!***********************************************************************
! Function of GETFML is to 
! (1) obtain a solution to the scattering problem for given target 
!     orientation, and 
! (2) return the scattering function fml in preselected scattering 
!     directions.

! Given:
!       AKR(1-3)=(k_x,k_y,k_z)*d for incident wave [in Target Frame]
!       AK3     =(kd)**3
!       AKSR(1-3,1-NSCAT)=(k_x,k_y,k_z)*d for NSCAT scattering direction
!       GAMMA = parameter controlling summation over replica dipoles
!               (smaller alpha -> longer summation)
!       BETADF(1-NAT)=orientation angle beta (radians) describing
!                orientation of "Dielectric Frame" relative to Target Frame
!                for dipoles 1-NAT0
!       PHIDF(1-NAT)=orientation angle phi (radians) describing orientation
!                of Dielectric Frame relative to Target Frame for dipole
!                1-NAT0
!       THETADF(1-NAT)=orientation angle theta (radians) describing
!                orientation of Dielectric Frame relative to Target Frame
!                for dipoles 1-NAT0
!       DX(1-3) =(d_x/d,d_y/d,d_z/d) where d_x,d_y,d_z=x,y,z lattice spacing
!                and d=(d_x*d_y*d_z)**(1/3)=effective lattice spacing
!       CALPHA  =descriptor of method used for assigning polarizabilitie
!               =LATTDR or DRAI88 or GOBR88
!       CMDSOL  =descriptor of method used for iterative solution
!       CMDFFT  =descriptor of method used for FFTs
!       CMDTRQ  =descriptor of whether or not to compute torques
!       CSHAPE  =descriptor of target shape, needed by subroutine ALPHA
!       CXE01R(1-3)=incident polarization state 1 at origin [in TF]
!       CXE02R(1-3)=incident polarization state 2 at origin [in TF]
!       CXEPS(1-NCOMP)=dielectric constant for compositions 1-NCOMP
!       EM1R(1-3,1-NSCAT)=unit scat. polarization vectors 1 in TF
!       EM2R(1-3,1-NSCAT)=unit scat. polarization vectors 2 in TF
!       ETASCA  =parameter controlling number of scattering angles used
!                for calculation of radiation force, <cos>, <cos^2>, and
!                radiation torque
!       IBETH   =MYID+IBETH1 if first time through combined BETA/THETA
!                orientation loop
!       IBETH1  =starting value of IBETH (see above)
!       ICOMP(1-NAT3)=x,y,z "composition" for sites 1-NAT
!       IDVOUT  =device for running output
!       INIT    =0,1,2 for choice of |x0> for CCG method
!       IOCC(1-NAT)=0,1 if site in extended target is vacant,occupied
!       IORTH   =1,2 to do 1,2 orthogonal pol states
!       IPHI    =which phi value for target orientation
!       IPHI1   =starting value of IPHI (see above)
!       IPBC    = 0 if PBC not used
!               = 1 if PBC used
!       JPBC    = 0 if PBC not used
!               = 1 if PBC used in y direction only
!               = 2 if PBC used in z direction only
!               = 3 if PBC used in both y and z directions
!               **NOTE** JPBC is not used by GETFML
!                        When JPBC > 0, GETFML computes f_ml for the
!                        Target Unit Cell.
!       ITASK    not used
!       IXYZ0(1-MXNAT,1-3)=lattice coordinates for sites 1-NAT (in TF)
!                with sites 1-NAT0 in list being the occupied sites.
!       DX(1-3) = location/DX(1-3) in TF corresponding to lattice site
!                 (IX,IY,IZ)=(0,0,0)
!       LACE,LAXI,LCLM,LGI,LPI,LQI,LSC0=integers used to assign location
!                in scratch array
!       MXCOMP  =dimensioning info: maximum number of allowed compositions
!       MXCXSC  =dimensioning info for complex scratch space
!       MXN3    =dimensioning info (3*max number of sites)
!       MXNAT   =dimensioning info (max number of sites)
!       MXNX    =dimensioning info (max x extent of target)
!       MXNY    =dimensioning info (max y extent of target)
!       MXNZ    =dimensioning info (max z extent of target)
!       MXPHI   =dimensioning info (max number of phi vals for targ. ori
!       MXSCA   =dimensioning info (max number of scat. dirs for f_ml)
!       MYID    =parallel process identifier (=0 if only 1 process)
!       NAT     =number of sites in extended target
!       NAT0    =number of sites in original target
!       NAT3    =3*NAT
!       NCOMP   =number of different dielectric tensor elements in target
!       NSCAT   =number of scat. dirs for f_ml
!       NX      =x extent of extended target
!       NY      =y extent of extended target
!       NZ      =z extent of extended target
!       PHI(1-NPHI)=phi values for target orientation
!       PIA2    =\pi*(3*NAT/4\pi)^{2/3} (NAT=number of sites in original
!       SHPAR(1-10)=target shape parameters, needed by subroutine ALPHA
!       TOL     =tolerance for terminating conjugate gradient iteration

! and scratch space:
!       CXRLOC(1-MXCOMP+1,3,3)
!       CXSC(1-MXCXSC)=scratch space
!       CXSCR1(1-MXN3)=scratch space
!       CXZC(1->MXNX+1,1->MXNY+1,1->MXNZ+1,1-6)=scratch space
!       CXZW(1-MXNAT,1-24)=scratch space
!       SCRRS1(1-MXNAT,3)=scratch space
!       SCRRS2(1-MXNAT)=scratch space

! Returns:

!       CXADIA(1-3,1-NAT)=diagonal elements of "A matrix"
!       CXAOFF(1-3,1-NAT)=off-diagonal elements of 3x3 blocks on
!                         diagonal of "A matrix"
!       CXALPH(1-3,1-NAT)=diagonal elements of 3x3 polarizability
!                         tensor for dipoles 1-NAT (in TF)
!       CXALOF(1-3,1-NAT)=off-diagonal elements of 3x3 polarizability
!                         tensor for dipoles 1-NAT
!       CXALOS(J,1-3)=(alpha_23,alpha_31,alpha_12)/d^3 for dipole J=1-NA
!                     Off-diagonal polarizability elements
!       CXE(1-NAT3)=incident x,y,z E field at dipoles 1-NAT (in TF)
!       CXF11(1-NSCAT)=scattering matrix element f_11 for NSCAT dirs
!       CXF12(1-NSCAT)=                          f_12
!       CXF21(1-NSCAT)=                          f_21
!       CXF22(1-NSCAT)=                          f_22
!       CXPOL(1-NAT3)=x,y,z polarization of dipoles 1-NAT in TF
!       QABS(1-2)=Q_abs=C_abs/(PIA2*d**2) for incident pols 1,2
!       QBKSCA(1-2)=diff.scatt.cross section/(PIA2*d^2) for backscat
!                   for inc.pols.1,2
!       QEXT(1-2)=Q_ext=C_ext/(PIA2*d**2) for incident pols 1,2
!       QPHA(1-2)=phase shift cross section for incident pols 1,2
!       QSCA(1-2)=Q_sca=C_sca/(PIA2*d**2) for incident pols 1,2
!       QSCAG(1-3,1-2)=<cos(theta)>*Q_sca for incident pols 1,2
!                      <sin(theta)*cos(phi)>*Q_sca for incident pols 1,2
!                      <sin(theta)*sin(phi)>*Q_sca for incident pols 1,2
!                      [QSCAG(1-3,1-2) is in Lab Frame]
!       QSCAG2(1-2)   =<cos^2(theta)>*Q_sca for inciden tpols 1,2
!       QTRQAB(1-3,1-2)=vector torque cross section/(PIA2*d^2) for
!                      torque due to absorption, along axes 1-3, for
!                      incident pols 1,2
!                      [QTRQAB(1-3,1-2) is in Lab Frame]
!       QTRQSC(1-3,1-2)=vector torque cross section/(PIA2*d^2) for
!                      torque due to scattering, along axes 1-3, for
!                      incident pols 1,2
!                      [QTRQSC(1-3,1-2) is in Lab Frame]
!       TIMERS(1-12)=timing information
!       ITNUM(1-2)= number of iterations taken for polarizations 1 and 2
!       MXITER    = maximum number of iterations allowed.
! if JO=1:
!       CXSCR1(1->3*NAT0)=reduced polarization vector for incident pol. 1
!       CXSC(LACE->LACE-1+3*NAT0)=reduced polarization vector for incident
!                                 pol. 2
! Requires:
!	MATVEC = external routine to compute A*x
!	CMATVEC = external routine to compute conjg(A')*x
!	where A is matrix such that A*pol=E
!	and x is arbitrary vector

! B.T.Draine, Princeton Univ. Observatory, 1990.
! History:
! 90.11.06 (BTD): modified to change argument list and pass array
!                 dimensions.
! 90.11.29 (BTD): further modifications, including arguments
!                 changed calls to EVALE
! 90.11.29 (BTD): changes to use REDUCE so that computations in
!                 EVALQ and SCAT are restricted to occupied sites.
! 90.11.30 (BTD): multiple changes to reduce storage requirements
! 90.12.03 (BTD): change ordering of XYZ0
! 90.12.10 (BTD): remove XYZ0 and use IXYZ0 instead
! 90.12.13 (BTD): add new argument to calls to EVALQ
! 91.01.03 (BTD): modified so that when IORTH=1 no attempt is made
!                 to compute f_{ml} for "standard" definition of
!                 incident pol. states defined by scattering plane.
! 91.05.09 (BTD): after call to SCAT, adjust QSCA and QSCAG to make
!                 QSCA=QEXT-QABS
! 91.05.09 (BTD): eliminate redundant write statement
! 91.08.14 (BTD): add QBKSCA to argument list for GETFML
!                 add CBKSCA to argument list for SCAT
!                 compute QBKSCA from CBKSCA
! 91.09.17 (BTD): move call to subroutine ALPHA from DDSCAT to GETFML
!                 (since polarizability from Lattice Dipersion Relation
!                 depends on direction of propagation and polarization)
!                 add CALPHA, CXEPS, ICOMP, MXCOMP to argument list for
!                 GETFML
! 93.01.14 (BTD): change method for computing QSCA to use result from
!                 SCAT when albedo << .01 , to use (QEXT-QABS) when
!                 albedo >> .01, and to smoothly join these limits.
! 93.03.11 (BTD): deleted all code associated with unused variable
!                 CMACHN (originally included to identify machine/OS)
! 93.07.09 (BTD): modify method for computing QSCA to use .03 as
!                 the dividing point rather than .01 .
! 93.11.23 (BTD): added SAVE E02 to preserve value of E02 between calls
!                 [to eliminate problem encountered on Silicon Graphics
!                 machines when IPHI>1]
! 94.06.20 (PJF): Added DTIME to argument list for TIMEIT
! 95.06.14 (BTD): Changed CSCAG from scalar to CSCAG(3) vector
!                 where CSCAG(1)=C_sca*<cos(theta)>
!                 Changed QSCAG(1-2) to QSCAG(1-3,1-2)
!                 Added CTRQSC(1-3) to argument list of SCAT
!                 Added QTRQSC(1-3,1-2) to argument list of GETFML
! 95.06.15 (BTD): Added variable CXE to argument list of SCAT
!                 Added CTRQAB(1-3) to argument list of SCAT
!                 Added QTRQAB(1-3,1-2) to argument list of GETFML
! 95.06.16 (PJF+: Modified to eliminate call to DDACCG and replace
!           BTD): with code to call PETR (a modularized implementation
!                 of the Petravic-Kuo CCG algorithm).  With this change
!                 it becomes easier to substitute other implementations
!                 of CCG (or other iterative methods, such as QMR).
! 95.06.19 (BTD): Added CMDSOL to argument list to select CCG method
!               : Added CMDTRQ to argument list to allow torque
!                 calculation to be skipped
! 95.07.10 (BTD): Changed calls to TIMEIT (was not doing timing
!                 properly)
! 95.07.17 (BTD): Relocated call to EVALQ.
!                 Call REDUCE to "reduce" CXE as well as CXP, prior
!                 to calling EVALQ and SCAT.
! 95.07.20 (BTD): Added scratch space SCRRS2 for use by SCAT
!                 Added scratch array CXZW(1-MXNAT,10-12) to SCAT
!                 argument list
! 95.07.27 (BTD): Revised way in which CXALPH is "reduced" prior to
!                 calls to EVALQ and SCAT.
!                 Revised to compute CXALPH for each incident E field
!                 [Only reason this is needed is because of slight
!                 dependence of polarizabilities on incident pol state]
!                 Added scratch space CXSCR1 to argument list, to store
!                 first polarization solution while solving for second
! 95.08.11 (BTD): Relocated CALL CINIT to precede IF statements for
!                 selection of solution method.
! 95.08.14 (BTD): Added COMMON/NORMERR/ERRSCAL to support desired
!                 normalization of error printed out by subroutine
!                 PROGRESS (part of cgcommon.f)
! 96.11.14 (PJF): ADD TIMERS, MXTIMERS, NTIMERS  to formal parameters
!                 modify code to return timer results
!                 TIMERS contains information needed to time the
!                 code (for example CPU time or number of iterations)
!                 timers(1) --- CG CPU time for all iterations
!                 timers(2) --- number of iterations needed
!                 timers(3) --- not used
!                 timers(4) --- scat CPU time
!                 ... see code for timers(5-12)
! 96.11.21 (BTD): declared MXTIMERS,NTIMERS,
! 97.02.20 (BTD): removed SNORM2 from EXTERNAL statement
! 97.11.01 (BTD): added DX(1-3) to argument list of GETFML
!                 and to argument lists for EVALE, ALPHA, and SCAT
!                 [development of DDSCAT.6.0]
! 97.12.26 (BTD): added CXALOF to argument list of GETFML
!                 and to argument list of ALPHA, EVALA
!                 [development of DDSCAT.6.0]
! 98.01.01 (BTD): added 2 calls to new routine RESTORE to expand
!                 CXADIA and CXAOFF back to 3*NAT elements after they
!                 have been "reduced" to 3*NAT0 for efficient
!                 execution of EVALQ
! 98.04.27 (BTD): added declarations for arguments CXAOFF and CXALOF
!                 [development of DDSCAT.6.0]
!                 only call CXALPH once for each orientation, since
!                 in new form of Lattice Dispersion Relation the
!                 polarizability does not depend on polarization,
!                 although it does depend on direction of propagation
! 98.11.18 (BTD): corrected misleading comments apropos EM1R,EM2R
! 98.12.04 (BTD): located error which resulted in subsequent incorrect
!                 evaluation of Mueller scattering matrix.
!                 Problem arose because at some point it was decided to
!                 convert to usual convention for incident polarization
!                 states 1 and 2 for computation of f_ml (CXF11,...)
!                 Unfortunately, this was forgotten, and subroutine
!                 GETMUELLER was written assuming that incident pol.
!                 states l=1 and 2 referred to by f_ml were the ones
!                 specified by the user in ddscat.par.  Because this in
!                 fact does seem preferable, we conform to this
!                 prescription and have disabled the portion of GETFML
!                 which converted to the "usual" convention.  Note
!                 that since f_ml is now entirely for internal use,
!                 (except when IORTH=1, but then we do not compute
!                 Muller matrix) there is no motivation to adopt usual
!                 convention.
! 98.12.07 (BTD): Experiment to verify claimed accuracy from PIM
!                 New module has been inserted into code below but
!                 commented out after verifying that reported errors
!                 appear to be close to "true" errors.
!                 However, module can be reactivated by uncommenting
!                 if such verification is desired in future.  Search
!                 below for string "98.12.07" to find the 2 modules.
! 98.12.16 (BTD): Increase upper limit on number of iterations from 300
!                 to 10000
! 03.01.29 (BTD): returned CXE0R to argument list for ALPHA
!                 reactived calls to ALPHA for different polarizations
!                 since now using LATTDR again instead of GKDLDR
! 03.04.12 (BTD): added variables ITNUM(2) and MXITER to argument list
!                 and modified code to assign values to them
!                 added veriables ITERN and ITERMX
! 03.10.23 (BTD): modified argument list for SCAT:
!                    eliminated ICTHM and IPHIM
!                    added CSCAG2 and NAVG
!                 modified argument list of GETFML
!                    eliminated ICTHM and IPHIM
! 03.10.24 (BTD): added ETASCA to argument list
!                 and to argument list of SCAT
! 04.03.31 (BTD): changed WRITE(0 to WRITE(IDVOUT
! 04.05.22 (BTD): cleanup: comment out declarations of unused variables
!                 CXVAR1,CXVAR1,CXONE,CXZERO
!                 remove IF(IPHIM.GT.1) test, and always store
!                 first polarization (even if only doing one
!                 target orientation, overhead is not significant)
! 04.05.22 (BTD): fixed bug in computation of true error
!                 was summing from 3*I+1 to 3*I+3
!                 now sum from 3*I-2 to 3*I
! 04.09.14 (BTD): added new variables BETADF(1-NAT0),PHIDF(1-NAT0),
!                 THETADF(1-NAT0) describing orientation of anisotropic
!                 target material relative to Target Frame, for each
!                 dipole location.
! 06.09.28 (BTD): version 6.2.3:
! 06.09.29 (BTD): added IPBC to argument list of ALPHA (3 places)
!                 [note that CXZC and CXSC are not used in present
!                 version of ALPHA, but we keep CXZC and CXSC in arg.
!                 list in case some future version of ALPHA might
!                 require complex scratch space
! 07.06.20 (BTD): added X0 to argument list
!                 added X0 to argument list in 3 calls to EVALE
!                 added X0 to argument list in 4 calls to SCAT
! 07.06.22 (BTD): removed THETAN and PHIN from argument list
!                 (not needed, since scattering directions are specified
!                 through vectors AKS)
!                 removed THETAN and PHIN from arg list of subroutine
!                 SCAT
! 07.08.04 (BTD): Version 7.0.3
!                 Replaced COMMON/NORMERR/ with USE DDCOMMON_9
! 07.08.10 (PJF): Added call to ZBCG2
! 07.08.15 (BTD): Use TOLR to communicate with ZBCG2
! 07.09.11 (BTD): Changed IXYZ0 from INTEGER*2 to INTEGER
! 07.10.27 (BTD): Changed SHPAR(6) -> SHPAR(10)
! 08.01.06 (BTD): Cosmetic changes
! 08.02.17 (BTD): Changed SHPAR(10) -> SHPAR(12)
! 08.03.11 (BTD): v2 / ver7.0.5
!                 added ALPHA to argument list
! 08.03.14 (BTD): corrected dimensioning of
!                 IXYZ0(MXNAT,3) -> IXYZ0(NAT0,3)
!                 BETADF(MXNAT) -> BETADF(NAT0)
!                 PHIDF(MXNAT) -> PHIDF(NAT0)
!                 THETADF(MXNAT) -> THETADF(NAT0)
! 08.04.19 (BTD): changed notation: ALPHA -> GAMMA
!                 reordered arguments
! 08.05.08 (BTD): changed SCRRS1(MXNAT3) -> SCRRS1(MXNAT,3)
! 08.05.12 (BTD): changed declaration of DUMMY from REAL to EXTERNAL
!                 as suggested by Lazanoff
! 08.07.16 (BTD): changed 
!                 ITERMX=IPAR(10) 
!                 -> 
!                 ITERMX=MXITER
!                 IPAR(10)=MXITER
! end history

! Copyright (C) 1993,1994,1995,1996,1997,1998,2003,2004,2006,2007,2008
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
! diagnostic
!      write(0,*)'getfml ckpt 1, mxnat=',mxnat,' myid=',myid
!***
      IDVOUT2=IDVOUT
      NAT03=3*NAT0

!*** Compute fml directly only if IPHI=1
!    Otherwise use previously computed fml values to obtain new fml

      IF(IPHI==1)THEN

         DO JO=1,IORTH
            IF(JO==1)THEN
               CXE0R(1)=CXE01R(1)
               CXE0R(2)=CXE01R(2)
               CXE0R(3)=CXE01R(3)
            ELSE
               CXE0R(1)=CXE02R(1)
               CXE0R(2)=CXE02R(2)
               CXE0R(3)=CXE02R(3)
            ENDIF
            E02=0._WP
            DO I=1,3
               E02=E02+REAL(CXE0R(I)*CONJG(CXE0R(I)))
            ENDDO

! Compute quantity used for "normalizing" error:

            ERRSCAL=SQRT(NAT0*E02)

!*** call EVALE for NAT sites
!*** diagnostic
!            write(0,*)'getfml ckpt 2, myid=',myid
!***
            CALL EVALE(CXE0R,AKR,DX,X0,IXYZ0,MXNAT,MXN3,NAT,NAT0,NX,NY,NZ,CXE, &
                 AEFFA,WAVEA,MXRAD,MXWAV,CENTER,c,velocity,e_charge, &
                 DielectricConst)!Arguments AEFFA and after added by NWB 3/8/12
!*** diagnostic
!            write(0,*)'getfml ckpt 3, myid=',myid
!***

!***call ALPHA to determine polarizabilities at NAT sites
!    (note: this has to be within the loop over directions and
!     polarizations because Lattice Dispersion Relation polarizabilities
!     depend on direction of propagation and on polarization state)

            CALL ALPHADIAG(AKR,BETADF,PHIDF,THETADF,CALPHA,CXALPH,    &
                           CXALOF,CXALOS,CXE0R,CXEPS,CXSC,CXSCR1,CXZC,CXZW,  &
                           DX,IBETH,IBETH1,ICOMP,IOCC,IPBC,IPHI,IPHI1,IXYZ0, &
                           JO,MYID,MXCOMP,MXNAT,MXN3,NAT,NAT0,NAT3,NCOMP,NX, &
                           NY,NZ,CXRLOC,CSHAPE,SHPAR)

!*** diagnostic
!            write(0,*)'getfml ckpt 4, myid=',myid
!***

!*** call EVALA to prepare A matrix elements

            CALL EVALA(CXADIA,CXAOFF,CXALPH,CXALOF,MXN3,NAT)

! PIMSETPAR sets parameters used by PIM package
! ipar(1) LDA         LDA (Leading dimension of a)
! ipar(2) N           N   (Number of rows/columns of a)
! ipar(3) BLKSZ       N   (Size of block of data; used when data is
!                          partitioned using cyclic mode)
! ipar(4) LOCLEN      N   (Number of elements stored locally; for
!                          sequential=n)
! ipar(5) BASISDIM    C=basis=10 (Dimension of orthogonal basis, used in
!                                 GMRES)
! ipar(6) NPROCS      -1 (Number of processors)
! ipar(7) PROCID      -1 (Processor identification)
! ipar(8) PRECONTYPE PRET (=1 Type of preconditioning
!                     0-no preconditioning, 1-left, 2-right, 3-symmetric
! ipar(9) STOPTYPE   STOPT (Type of stopping criteria used)
! ipar(10) MAXIT     MAXIT (=int(n/2) Maximum number of iterations allowed

! on return from PIM, following are defined:
!      IPAR(11) = itno   (Number of iterations executed)
!      IPAR(12) = exit status 
!                 0: converged
!                -1: no convergence has been achieved
!                -2: soft breakdown, solution may have been found
!                -3: hard breakdown, no solution
!                -4 to -11 : other conditions
!      IPAR(13) = if IPAR(12) = -2 or -3, gives the step number in the
!                 algorithm where a breakdown has occurred.

! Set upper limit on iterations to smaller of 1e4 and NAT0

            MXITER=MIN(10000,NAT03)

! pass information on max.no. of iterations allowed to DDCOMMON_9

!BTD 080716: changed
!            ITERMX=IPAR(10)
            ITERMX=MXITER
            IPAR(10)=MXITER
!-------------------
            SPAR(1)=TOL

! Clean CXE:

            IF(NAT0<NAT)CALL NULLER(CXE,IOCC,MXNAT,MXN3,NAT)

! Iterate to improve CXPOL

!*** diagnostic
!            write(0,*)'getfml ckpt 4.1, myid=',myid
!***
            IF(CMDSOL=='PETRKP')THEN
               CALL PIMSSETPAR(IPAR,SPAR,MXN3,NAT3,NAT3,NAT3, &
                               10,-1,-1,1,5,MXITER,TOL)
               CALL CINIT(NAT3,CMPLX(0._WP,0._WP,KIND=WP),CXPOL,1)
               CALL TIMEIT('PETRKP',DTIME)
               CALL PETR(CXPOL,CXE,CXSC,MXN3,IPAR,SPAR,MATVEC,CMATVEC)
               CALL TIMEIT('PETRKP',DTIME)
               TIMERS(1)=DTIME
               TIMERS(2)=REAL(IPAR(11),KIND=WP)

            ELSEIF(CMDSOL=='PBCGS2')THEN
!*** diagnostic
!               write(0,*)'getfml ckpt 4.2, myid=',myid
!***
               CALL CINIT(NAT3,CMPLX(0._WP,0._WP,KIND=WP),CXPOL,1)
!*** diagnostic
!               write(0,*)'getfml ckpt 4.3, myid=',myid
!***
               CALL TIMEIT('PBCGS2',DTIME)
!*** diagnostic
!               write(0,*)'getfml ckpt 4.4, myid=',myid
!***

!BTD 07.08.15 Use TOLR to input TOL with ZBCG2 (ZBCG2 returns
!             actual achieved tolerance in TOLR)

               TOLR=TOL
!*** diagnostic
!               write(0,*)'getfml ckpt 4.5, myid=',myid
!***

               CALL ZBCG2(.TRUE.,2,NAT3,CXPOL,.FALSE.,CXE,MATVEC,PRECOND, &
                           TOLR,MXITER,CXSC,IPAR(12))

!*** diagnostic
!               write(0,*)'getfml ckpt 4.6, myid=',myid
!***
               CALL TIMEIT('PBCGS2',DTIME)
							   
            ELSEIF(CMDSOL=='PBCGST')THEN

! CALL PIMSGETPAR(IPAR,SPAR,LDA,N,BLKSZ,LOCLEN,BASISDIM,NPROCS,PROCID, 
!                 PRECONTYPE,STOPTYPE,MAXIT,ITNO,STATUS,STEPERR,EPSILON,EXITNORM)
!flatau warning something is wring with timers because they also time eself - check this
!we need to write better timming routine

!flatau I am just hardwiring number of iterations (mxiter is not used anymore)
!BTD 07.08.12 minor revision

               NO_CG_RESTART=4
               IPAR(12)=-1
	       NO_CG_ITER=0

               CALL PIMSSETPAR(IPAR,SPAR,MXN3,NAT3,NAT3,NAT3, &
                               10,-1,-1,1,2,NO_CG_RESTART,TOL)

               CALL CINIT(NAT3,CMPLX(0._WP,0._WP,KIND=WP),CXPOL,1)
               CALL TIMEIT('PBCGST',DTIME)
!	       DO WHILE(IPAR(12)==-1 .AND. NO_CG_ITER < MX_CG_ITER)

               DO ITER=0,(MXITER/NO_CG_RESTART)
                  IF(IPAR(12).NE.0)THEN
                     IF(ITER.GT.0)WRITE(0,*)'restart PIMZBICGSTAG: ',ITER
                     CALL PIMCBICGSTAB(CXPOL,CXE,CXSC,IPAR,SPAR,MATVEC,DIAGL, &
                                       DUMMY,PCSUM,PSCNRM2,PROGRESS)
                  ENDIF
!                 NO_CG_ITER=NO_CG_ITER+IPAR(11)
               ENDDO


               CALL TIMEIT('PBCGST',DTIME)
               TIMERS(1)=DTIME
               TIMERS(2)=REAL(NO_CG_ITER,KIND=WP)
	    ELSE
               STOP 'Error -- invalid CMDSOL in getfml'
            ENDIF

!*** Reduce the vectors CXADIA, CXAOFF, CXPOL and CXE to retain only
!    occupied sites, to eliminate unnecessary calculations for
!    unoccupied sites when calling EVALQ and SCAT

!*** diagnostic
!            write(0,*)'getfml ckpt 4.7, myid=',myid
!***
            CALL REDUCE(CXADIA,IOCC,MXN3,MXNAT,NAT,NAT0)
            CALL REDUCE(CXAOFF,IOCC,MXN3,MXNAT,NAT,NAT0)
            CALL REDUCE(CXE,IOCC,MXN3,MXNAT,NAT,NAT0)
            CALL REDUCE(CXPOL,IOCC,MXN3,MXNAT,NAT,NAT0)

!*** diagnostic
!            write(0,*)'getfml ckpt 4.8, myid=',myid
!***

! ITERN is obtained from DDCOMMON_9

            ITNUM(JO)=ITERN

!*** Store reduced polarization vector P
!    for use in subsequent computations for other phi values
!    CXSCR1     for JO=1
!    CXSC(LACE) for JO=2 (this scratch space no longer needed by solver)

            IF(JO==1)THEN
!*** diagnostic
!            write(0,*)'getfml ckpt 4.81, myid=',myid
!***
               CALL COPYIT(CXPOL,CXSCR1,NAT03)
            ELSEIF(JO==2)THEN
!*** diagnostic
!            write(0,*)'getfml ckpt 4.82, myid=',myid,' lace=',lace
!***
               CALL COPYIT(CXPOL,CXSC(LACE),NAT03)
            ENDIF
!*** diagnostic
!            write(0,*)'getfml ckpt 4.9, myid=',myid
!***

!            CALL EVALQ(CXADIA,CXAOFF,AKR,NAT03,E02,CXE,CXPOL,CABS,CEXT,CPHA, &
!                       MXN3,1,h_bar,h_bar2) !All parameters h_bar and after added by NWB 7/11/12
            CALL EVALQ(NAT03,CXE,CXPOL,CABS,CEXT,CPHA,1,MXN3,h_bar,h_bar2)

            QABS(JO)=CABS/PIA2
            QEXT(JO)=CEXT !NEW CODE NWB 3/8/12
            QPHA(JO)=CPHA/PIA2
            WRITE(IDVOUT,9010)QABS(JO),QEXT(JO)

!*** first call to timeit:

            CALL TIMEIT(' SCAT ',DTIME)

            IF(JO==1)THEN

! call SCAT to compute CXF11,CXF21

!*** diagnostic
!            write(0,*)'getfml ckpt 5, myid=',myid
!            write(0,*)'       NAT0=',nat0
!            write(0,*)'       X0=',X0
!***

               CALL SCAT(AKR,AKSR,DX,EM1R,EM2R,E02,ETASCA,CMDTRQ,CBKSCA,CSCA, &
                         CSCAG,CSCAG2,CTRQAB,CTRQSC,CXE,CXE01R,CXF11,CXF21,   &
                         CXPOL,CXZW(1,1),CXZW(1,4),CXZW(1,7),CXZW(1,10),MXN3, &
                         MXNAT,MXSCA,MYID,JPBC,NAT0,NAT03,NAVG,NSCAT,SCRRS1,  &
                         SCRRS2,IXYZ0,X0)
!*** diagnostic
!            write(0,*)'getfml ckpt 6, myid=',myid
!***

            ELSEIF(JO==2)THEN

! call SCAT to compute CXF12,CXF22

!*** diagnostic
!            write(0,*)'getfml ckpt 7, myid=',myid
!***
               CALL SCAT(AKR,AKSR,DX,EM1R,EM2R,E02,ETASCA,CMDTRQ,CBKSCA,CSCA, &
                         CSCAG,CSCAG2,CTRQAB,CTRQSC,CXE,CXE01R,CXF12,CXF22,   &
                         CXPOL,CXZW(1,1),CXZW(1,4),CXZW(1,7),CXZW(1,10),MXN3, &
                         MXNAT,MXSCA,MYID,JPBC,NAT0,NAT03,NAVG,NSCAT,SCRRS1,  &
                         SCRRS2,IXYZ0,X0)
!*** diagnostic
!            write(0,*)'getfml ckpt 8, myid=',myid
!***

            ENDIF

!*** second call to timeit
!*** diagnostic
!            write(0,*)'getfml ckpt 8.1, myid=',myid
!***

            CALL TIMEIT(' SCAT ',DTIME)
!*** diagnostic
!            write(0,*)'getfml ckpt 8.2, myid=',myid
!            write(0,*)'                dtime=',dtime
!            write(0,*)'                pia2=',pia2
!            write(0,*)'                jo=',jo,' qext(jo)=',qext(jo)
!            write(0,*)'                cbksca=',cbksca
!            write(0,*)'                csca=',csca
!***

            TIMERS(4)=DTIME

            IF(JPBC==0)THEN

! Note: at present, SCAT only calculates CBKSCA and CSCA when JPBC=0
!       therefore only calculate QBKSCA(JO), etc, when JPBC=0

               QBKSCA(JO)=CBKSCA/PIA2
               QSCA(JO)=CSCA/PIA2

!*** Since Q_ext and Q_abs are calculated "exactly", when particles
!    are "large" we obtain Q_sca from (Q_ext-Q_abs) rather than
!    scattering cross section computed by SCAT.  This is particularly
!    advantageous when the particles are large with complex scattering
!    patterns, so that angular integration of the differential scatterin
!    cross section would be either very time-consuming or inaccurate.
!    However, when particles are "small" the scattering cross section
!    is small compared to Q_ext and Q_abs, so that the difference
!    (Q_ext-Q_abs) may be inaccurate.  Accordingly, when particles
!    are "small" we instead accept the scattering cross section
!    computed by SCAT (probably quite accurate when albedo is small beca
!    the scattering pattern tends to be dipolar and therefore well-resol
!    by a reasonable number of scattering directions ICTHM and IPHI).
!    We use an ad-hoc procedure to smoothly make the transitions
!    between these regimes, based on the albedo.
!    The "dividing" line of albedo=.03 is chosen on basis of
!    numerical experiments.  If Q_ext and Q_abs are each calculated to
!    accuracy of 1 times 10^{-4} (say), then the
!    difference between Q_ext and Q_abs would be accurate to
!    0.0001*Q_ext*sqrt(2) .  When albedo=.03 (Q_ext-Q_abs)/Q_sca would
!    then be accurate to .0001*sqrt(2)/.03=.005

               FALB=(QSCA(JO)/(.03_WP*QEXT(JO)))**2
               QSCA(JO)=(QSCA(JO)+(QEXT(JO)-QABS(JO))*FALB)/(1._WP+FALB)

!*** diagnostic
!               write(0,*)'getfml ckpt 8.21, myid=',myid
!***
!------------------------------------

               ! Shuzhou Li added this on 26 NOV 2008
               ! Modified by NWB 29 Mar 2012
               IF (IWRPOL == 1) THEN
                  OPEN(UNIT=1093, FILE='Pol.dat', STATUS='UNKNOWN')
                  OPEN(UNIT=1095, FILE='incE.dat', STATUS='UNKNOWN')
                  DO i = 1, 3*NAT0
                     WRITE(1093, *) REAL(CXPOL(i)), IMAG(CXPOL(i))
                     WRITE(1095, *) REAL(CXE(i)), IMAG(CXE(i))
                  ENDDO
               ENDIF

!-------------------------------------

!    Compute QSCAG(1-3,JO), QTRQSC(1-3,JO), QSCAG2

               DO I=1,3
                  QSCAG(I,JO)=QSCA(JO)*CSCAG(I)/CSCA
               ENDDO
               QSCAG2(JO)=QSCA(JO)*CSCAG2/CSCA
               IF(CMDTRQ=='DOTORQ')THEN
                  DO I=1,3
                     QTRQAB(I,JO)=CTRQAB(I)/PIA2
                     QTRQSC(I,JO)=CTRQSC(I)/PIA2
                  ENDDO
               ENDIF

!*** diagnostic
!               write(0,*)'getfml ckpt 8.22, myid=',myid 
!***
            ENDIF !--- end IF(JPBC==0)
         ENDDO !--- end DO J=1,IORTH for IPHI=1

      ELSEIF(IPHI>=2)THEN

!*** Have previously computed dipole polarization vectors for
!    first phi value (and reduced them to NAT03 elements).
!    Use these two solutions to obtain Qabs,Qext, Qpha, Qsca, g*Qsca,
!    and fml for any additional values of target rotation phi.
!    CXSCR1     contains polarization vector for IPHI=1 and JO=1
!    CXSC(LACE) contains polarization vector for IPHI=1 and JO=2

! Incident polarization state 1 = CXE01R
! Call EVALE to obtain E at NAT0 occupied lattice sites
! first call to timing routine

!*** diagnostic
!         write(0,*)'getfml ckpt 8.3, myid=',myid
!***

         CALL TIMEIT(' EVALE',DTIME)

!*** Call EVALE for NAT0 occupied sites

!*** diagnostic
!            write(0,*)'getfml ckpt 9, myid=',myid
!***

         CALL EVALE(CXE01R,AKR,DX,X0,IXYZ0,MXNAT,MXN3,NAT0,NAT0,NX,NY,NZ,CXE, &
                    AEFFA,WAVEA,MXRAD,MXWAV,CENTER,c,velocity,e_charge, &
                    DielectricConst)!Arguments AEFFA and after added by NWB 3/8/12
!*** diagnostic
!            write(0,*)'getfml ckpt 10, myid=',myid
!***

! second call to timing routine

         CALL TIMEIT(' EVALE',DTIME)
         TIMERS(5)=DTIME

! Compute polarizabilities at NAT sites, then reduce to NAT0 occupied
! sites:

         CALL TIMEIT(' ALPHA',DTIME)

!*** diagnostic
!            write(0,*)'getfml ckpt 11, myid=',myid
!***

         CALL ALPHADIAG(AKR,BETADF,PHIDF,THETADF,CALPHA,CXALPH,CXALOF,      &
                        CXALOS,CXE01R,CXEPS,CXSC,CXSCR1,CXZC,CXZW,DX,IBETH, &
                        IBETH1,ICOMP,IOCC,IPBC,IPHI,IPHI1,IXYZ0,1,MYID,     &
                        MXCOMP,MXNAT,MXN3,NAT,NAT0,NAT3,NCOMP,NX,NY,NZ,     &
                        CXRLOC,CSHAPE,SHPAR)
!*** diagnostic
!            write(0,*)'getfml ckpt 12, myid=',myid
!***

         CALL REDUCE(CXALPH,IOCC,MXN3,MXNAT,NAT,NAT0)
         CALL REDUCE(CXALOF,IOCC,MXN3,MXNAT,NAT,NAT0)
         CALL TIMEIT(' ALPHA',DTIME)
         TIMERS(6)=DTIME

!*** Construct polarization for incident pol. state 1

         COSPHI=COS(PHI(IPHI)-PHI(1))
         SINPHI=SIN(PHI(IPHI)-PHI(1))
         DO I=1,NAT03
            CXPOL(I)=COSPHI*CXSCR1(I) - SINPHI*CXSC(LACE+I-1)
         ENDDO

!*** first call to timing routine

         CALL TIMEIT(' EVALQ',DTIME)

!         CALL EVALQ(CXADIA,CXAOFF,AKR,NAT03,E02,CXE,CXPOL,CABS,CEXT,CPHA, &
!                    MXN3,1,h_bar,h_bar2) !All parameters h_bar and after added by NWB 7/11/12
         CALL EVALQ(NAT03,CXE,CXPOL,CABS,CEXT,CPHA,1,MXN3,h_bar,h_bar2)

!*** second call to timing routine

         CALL TIMEIT(' EVALQ',DTIME)
         TIMERS(7)=DTIME
 
         QABS(1)=CABS/PIA2
         !QEXT(1)=CEXT/PIA2 !ORIGINAL CODE NWB 3/8/12
         QEXT(1)=CEXT !NEW CODE NWB 3/8/12
         QPHA(1)=CPHA/PIA2

!*** first call to timing routine

         CALL TIMEIT(' SCAT',DTIME)

!*** Now call SCAT to compute CXF11,CXF21
!    CXPOL has been reduced to NAT03 elements
!    first NAT03 elements of IXYZ0 correspond to physical sites

!*** diagnostic
!            write(0,*)'getfml ckpt 13, myid=',myid
!***

         CALL SCAT(AKR,AKSR,DX,EM1R,EM2R,E02,ETASCA,CMDTRQ,CBKSCA,CSCA,CSCAG, &
                   CSCAG2,CTRQAB,CTRQSC,CXE,CXE01R,CXF11,CXF21,CXPOL,         &
                   CXZW(1,1),CXZW(1,4),CXZW(1,7),CXZW(1,10),MXN3,MXNAT,MXSCA, &
                   MYID,JPBC,NAT0,NAT03,NAVG,NSCAT,SCRRS1,SCRRS2,IXYZ0,X0)

!*** diagnostic
!            write(0,*)'getfml ckpt 14, myid=',myid
!***
!*** second call to timing routine

         CALL TIMEIT(' SCAT',DTIME)
         TIMERS(8)=DTIME

         QBKSCA(1)=CBKSCA/PIA2
         QSCA(1)=CSCA/PIA2

!*** As above, compute Q_sca from either SCAT or (Q_ext-Q_abs)
!    depending on whether albedo << .03 or albedo >> .03

         FALB=(QSCA(1)/(.03_WP*QEXT(1)))**2
         QSCA(1)=(QSCA(1)+(QEXT(1)-QABS(1))*FALB)/(1._WP+FALB)
         DO I=1,3
            QSCAG(I,1)=QSCA(1)*CSCAG(I)/CSCA
            QTRQAB(I,1)=CTRQAB(I)/PIA2
            QTRQSC(I,1)=CTRQSC(I)/PIA2
         ENDDO
         QSCAG2(1)=QSCA(1)*CSCAG2/CSCA

!*** Polarization state 2:

!*** first call to timing routine

         CALL TIMEIT(' EVALE',DTIME)

!*** Call EVALE for NAT0 occupied sites
!    to obtain appropriate E vector for pol. state 2

!*** diagnostic
!            write(0,*)'getfml ckpt 15, myid=',myid
!***

         CALL EVALE(CXE02R,AKR,DX,X0,IXYZ0,MXNAT,MXN3,NAT0,NAT0,NX,NY,NZ,CXE, &
                    AEFFA,WAVEA,MXRAD,MXWAV,CENTER,c,velocity,e_charge, &
                    DielectricConst)!Arguments AEFFA and after added by NWB 3/8/12

!*** diagnostic
!            write(0,*)'getfml ckpt 16, myid=',myid
!***
!*** second call to timing routine

         CALL TIMEIT(' EVALE',DTIME)
         TIMERS(9)=DTIME

! Compute polarizabilities at NAT sites, then reduce to NAT0 occupied
! sites:

         CALL TIMEIT(' ALPHA',DTIME)

!*** diagnostic
!            write(0,*)'getfml ckpt 17, myid=',myid
!***
         CALL ALPHADIAG(AKR,BETADF,PHIDF,THETADF,CALPHA,CXALPH,CXALOF, &
                        CXALOS,CXE02R,CXEPS,CXSC,CXSCR1,CXZC,CXZW,DX,IBETH,   &
                    IBETH1,ICOMP,IOCC,IPBC,IPHI,IPHI1,IXYZ0,2,MYID,MXCOMP,    &
                    MXNAT,MXN3,NAT,NAT0,NAT3,NCOMP,NX,NY,NZ,CXRLOC,CSHAPE,    &
                    SHPAR)

!*** diagnostic
!            write(0,*)'getfml ckpt 18, myid=',myid
!***
         CALL REDUCE(CXALPH,IOCC,MXN3,MXNAT,NAT,NAT0)
         CALL REDUCE(CXALOF,IOCC,MXN3,MXNAT,NAT,NAT0)
         CALL TIMEIT(' ALPHA',DTIME)
         TIMERS(10)=DTIME

!*** Construct incident polarization vector for pol. state 2

         DO I=1,NAT03
            CXPOL(I)=COSPHI*CXSC(LACE+I-1)+SINPHI*CXSCR1(I)
         ENDDO

!*** first call to timing routine

         CALL TIMEIT(' EVALQ ',DTIME)

!         CALL EVALQ(CXADIA,CXAOFF,AKR,NAT03,E02,CXE,CXPOL,CABS,CEXT,CPHA, &
!                    MXN3,1,h_bar,h_bar2) !All parameters h_bar and after added by NWB 7/11/12
         CALL EVALQ(NAT03,CXE,CXPOL,CABS,CEXT,CPHA,1,MXN3,h_bar,h_bar2)

!*** second call to timing routine

         CALL TIMEIT(' EVALQ ',DTIME)
         TIMERS(11)=DTIME

         QABS(2)=CABS/PIA2
         !QEXT(2)=CEXT/PIA2 !ORIGINAL CODE NWB 3/8/12
         QEXT(2)=CEXT !NEW CODE NWB 3/8/12
         QPHA(2)=CPHA/PIA2

!*** Call SCAT to compute CXF12,CXF22
!    CXPOL has NAT03 elements
!    IXYZ0 was initially given NAT03 elements

!*** first call to timing routine:

         CALL TIMEIT(' SCAT',DTIME)
!*** diagnostic
!            write(0,*)'getfml ckpt 19, myid=',myid
!***

         CALL SCAT(AKR,AKSR,DX,EM1R,EM2R,E02,ETASCA,CMDTRQ,CBKSCA,CSCA,CSCAG, &
                   CSCAG2,CTRQAB,CTRQSC,CXE,CXE01R,CXF12,CXF22,CXPOL,         &
                   CXZW(1,1),CXZW(1,4),CXZW(1,7),CXZW(1,10),MXN3,MXNAT,MXSCA, &
                   MYID,JPBC,NAT0,NAT03,NAVG,NSCAT,SCRRS1,SCRRS2,IXYZ0,X0)

!*** diagnostic
!            write(0,*)'getfml ckpt 20, myid=',myid
!***
!*** second call to timing routine

         CALL TIMEIT(' SCAT',DTIME)
         TIMERS(12)=DTIME

         QBKSCA(2)=CBKSCA/PIA2
         QSCA(2)=CSCA/PIA2

!*** As above, compute Q_sca from either SCAT or (Q_ext-Q_abs)
!    depending on whether albedo << .03 or albedo >> .03

         FALB=(QSCA(2)/(.03_WP*QEXT(2)))**2
         QSCA(2)=(QSCA(2)+(QEXT(2)-QABS(2))*FALB)/(1._WP+FALB)
         DO I=1,3
            QSCAG(I,2)=QSCA(2)*CSCAG(I)/CSCA
            QTRQAB(I,2)=CTRQAB(I)/PIA2
            QTRQSC(I,2)=CTRQSC(I)/PIA2
         ENDDO
         QSCAG2(2)=QSCA(2)*CSCAG2/CSCA

      ENDIF
!*** diagnostic
!      write(0,*)'getfml ckpt 21, myid=',myid,' returning from getfml'
!***
      RETURN
9010  FORMAT (1X,'Q_abs =',1P,E10.3,' Q_ext= ',E10.3)
9200  FORMAT (' final true frac.err=',F12.7)
    END SUBROUTINE GETFML
