    SUBROUTINE ESELF(CMETHD,CXZP,NX,NY,NZ,IPBC,GAMMA,PYD,PZD,AK,AKD, &
                     DX,CXZC,CXZW,CXZE)
      USE DDPRECISION,ONLY: WP
      USE DDCOMMON_0,ONLY: AK2OLD,AK3OLD,IDIPINT,NGRID,WOLD
      IMPLICIT NONE

!----------------------- eself v7 --------------------------------
! Arguments:

      CHARACTER(6) :: CMETHD
      INTEGER :: IPBC,NX,NY,NZ
      REAL(WP) :: AKD,GAMMA,PYD,PZD
      REAL(WP) :: AK(3),DX(3)
      COMPLEX(WP) ::                                                 &
         CXZC(NX+1+IPBC*(NX-1),NY+1+IPBC*(NY-1),NZ+1+IPBC*(NZ-1),6), &
         CXZE(NX,NY,NZ,3),                                           &
         CXZP(NX,NY,NZ,3),                                           &
         CXZW(2*NX,2*NY,2*NZ,*)

! NB: module DDCOMMON_0 must have previously set values of
!       AK2OLD,AK3OLD,WOLD
!    to be used by ESELF

! Local scalars:

      CHARACTER :: CMSGNM*70
      INTEGER :: I,IR,ISGN,J,JR,JX,JY,JZ,JSGN,K,KR,KSGN,M
      REAL(WP) :: AKD2,DTIME,FAC,PYDDX,PZDDX
      COMPLEX(WP) :: CXEX,CXEY,CXEZ,CXXX,CXXY,CXXZ,CXYY,CXYZ,CXZZ

! Local arrays:

      INTEGER :: ISYM(3)

#ifdef openmp
      INTEGER NTHREADS,TID
#endif

      EXTERNAL CXFFTW,DIRECT_CALC,EXTND,PAD,TIMEIT,TRIM
      INTRINSIC MIN,NINT,SIGN
 
!-----------------------------------------------------------------------
! Parameter GAMMA determines the range of the sums when periodic
! boundary conditions are employed.  Dipole-dipole interactions
! are screened by a factor exp[-(gamma*kr)^4]
! The effective
! range/d = 1/(gamma*k*d) = 400 if gamma=5e-3 and kd=0.5
! range/lambda = 1/(2*pi*gamma) = 31.8
 
! The sums are actually continued out to
! r/d = 2*/(gamma*kd) = 800 if gamma=5e-3 , kd=0.5
! [screening factor = exp(-16)=1.1e-7]
!
! 
!-----------------------------------------------------------------------
! subroutine ESELF

!       Given the dipole moments, CXZP, at all points on
!       a rectangular grid, oscillating at frequency AKD,
!       compute the electric field amplitude, CXZE,
!       at each point produced by all the other dipoles except the one
!       at that point.
!       The relationship between the dipoles and the field
!       values at the grid points can be expressed as a convolution,
!       and the convolution is efficiently evaluated using 3D FFTs.

!       options for computation of 3-dimensional FFT:

!    if CMETHD='GPFAFT':
!          Use CXFFT3N interface to GPFA code of Temperton.
!          Good points:
!             -On CRAY the code is on average 30-40% faster in comparison
!              to TMPRTN (and 10-15 faster in comparison to BRENNR)
!             -On scalar machines the code is 2-8 faster in comparison
!              to BRENNR or TMPRTN
!             -Doesn't require additional storage
!          Limitations:
!              -Requires that NX,NY,NZ be of form (2**I)*(3**J)*(5**K);
!               subroutine EXTEND takes care of choosing suitable NX,NY,
!              -The choice of "lvr" variable in gpfa2f, gpfa3f, gpfa5f
!               depends on machine. WARNING: on C90 use lvr=128, on
!               all other CRAY's use lvr=64; wrong lvr will produce WRONG
!               RESULTS. On scalar machines optimal lvr depends on cache
!               length; sub-optimal choice degrades performance but still
!               produces correct results.
!   if CMETHD='FFTW21'
!      Use CXFFTW interface to FFTW (Fastest Fourier Transform in the
!               West) version 2.1.x from Frigo and Johnson.
!
!   if CMETHD='FFTMKL':
!      Use CXFFT3_MKL interface to Intel Math Kernel Library (MKL) FFT
!
! INPUT:

!       CXZP(I,J,K,L)   Lth cartesian component of the dipole
!                       moment at the grid point (I,J,K);
!                       the DIMENSIONed length of CXZP in
!                       the calling routine is CXZP(NX,NY,NZ,3)
!                       [or CXZP(NX*NY*NZ,3) or CXZP(3*NX*NY*NZ)]

!       NX,NY,NZ        Size of grid in x,y,z directions (INTEGER).

!       IPBC          = 0 for isolated target
!                     = 1 to use periodic boundary conditions
!                         (in either y direction, z direction, or both)
!       PYD             (Period of lattice in y direction)/DX(2)
!       PZD             (Period of lattice in z direction)/DX(3)

!       GAMMA         = coefficient used to assist convergence in
!                       lattice sums by suppressing long-range contributions
!                       with factor exp(-(gamma*k*r)^4)

!       DX(1-3)         Lattice spacing in x,y,z directions, in units of
!                       n**(-1./3.) .  Note that with this normalization
!                       we have DX(1)*DX(2)*DX(3)=1.

!       AK(1-3)         k(1-3)*d, where k = k vector in vacuo, and
!                       d = effective lattice spacing = (dx*dy*dz)**(1/3)

!       AKD           = (omega/c)*d = k*d (dimensionless)

!       CXZC            (NX+1)*(NY+1)*(NZ+1)*6 array of Green
!                       function coefficients used
!                       internally by ESELF and
!                       *************NOT TO BE OVERWRITTEN***********
!                       between calls, because these coefficients are
!                       recomputed only if W has changed since the last
!                       call to ESELF.

!       CXZW            Complex, scratch-space vector of length:
!                       2*NX*2*NY*2*NY*3
!                       See comment about FFT usage and CMETHD flag.
!                       Can be overwritten between calls to ESELF
!
! OUTPUT:

!       CXZE(I,J,K,L)   Lth component of dipole-generated electric field
!                       at grid point (I,J,K);
!                       the DECLARED length of CXZE in the calling
!                       program is CXZE(NX,NY,NZ,3)
!                       [or CXZE(NX*NY*NZ,3) or CXZE(3*NX*NY*NZ)]

! Originally written by Jeremy Goodman, 
! Princeton Univ. Observatory, 90.09.22
! History:
! 90.11.29 (BTD): Modified to set untransformed ZC(1,1,1,1-6)=0.
! 90.11.29 (PJF): Modified to use FOURX and CXFFT99
! 90.12.05 (BTD): Modified for new ordering of elements of polarization
!                 and electric field vectors in calling program.  
!                 Modified ESELF and PAD to remove distinction between 
!                 NX,NY,NZ and dimensions of CXZE and CXZP, 
!                 since our new ordering always assumes this.
! 90.12.13 (BTD): Modified to include CONVEX option.
! 92.04.20 (BTD): removed ISYM from argument list of TRIM (was not used)
! 94.06.20 (PJF): modified to call CXFFT3N when CMETHD='NEWTMP'
! 96.10.18 (BTD): changed NEWTMP to GPFAFT
! 97.10.16 (BTD): added DX(1-3) to argument list to support use for
!                 noncubic rectangular lattices.
!                 Added DX to 3 lines computing X(1-3)
! 99.04.26 (BTD): changed notation: CXY -> CXZE
! 00.06.22 (BTD): modified to support option FFTWFJ, with calls to CXFFT
! 00.06.25 (BTD): modified to eliminate calls to FOURX (CMETHD=BRENNR)
!                 and CXFFT3 (CMETHD=TMPRTN), as these options are no
!                 longer supported
! 00.07.05 (BTD): further cleanup
! 03.07.13 (BTD): changed FFTWFJ to FFTW21 to allow future distinction
!                 between FFTW 2.1.x and FFTW 3.0.x
! 04.03.05 (BTD): modified to handle a periodic lattice of scatterers
!                 with periodicity NPY*D(2) in y direction and
!                 NPZ*D(3) in z direction
!                 Add parameter BETA to assist convergence
! 05.06.16 (BTD): Replaced integer NPY,NPZ by real PYD,PZD in
!                 argument list and in calculation
!                 corrected error in calculation of JPZM
! 05.07.08 (BTD): corrected error in calculation of JPZM
! 05.08.01 (BTD): changes to increase efficiency of summations
!                 required for PBC option (new variables AKD2,
!                 PYDDX,PZDDX,X0,Y0,Z0,X2,Y2,X2Y2,CXTERM)
! 05.08.03 (BTD): corrected typo.
! 05.08.04 (BTD): added AK(1-3) to argument list,
!                 added local variables PHASY and PHASYZ
!                 added phase shift exp(i*PHASYZ) for replica dipole
! 06.09.15 (BTD): added comments, reduced BETA to 1e-12 for improved
!                 accuracy (at expense of speed).
! 06.09.23 (BTD): corrected error: sign error in calculation of
!                 diplacements X(2) and X(3) for replica dipoles
! 06.09.23 (BTD): corrected error in computation of factor CXFAC
!                 appearing in summation for CXSUM
! 06.09.28 (BTD): eself v2.0 and DDSCAT 6.2.3:
!                 * put part of calculation of A matrix into subroutine
!                   DIRECT_CALC
!                 * modified to support PBC option
!                 * added IPBC to argument list to support PBC option
!                 * changed dimensioning of CXZC when IPBC=1
!                 * when IPBC=1, store full CXZC rather than only first
!                   octant
! 07.08.06 (BTD): Version 7.0.3
!                 eliminate option CONVEX -> CXC3DFFT
!                 No longer appear to be any sites where this would be
!                 useful.
! 07.10.04 (BTD): Recompute Green function coefficients only if
!                 frequency changes by more than 1 part in 1e6
!                 (previous condition of requiring frequency to be
!                 unchanged was evidently being confused by roundoff
!                 error, leading to unnecessary recalculation).
! 07.10.25 (BTD): Modified so that for IPBC>0, Green functions are
!                 recalculated when direction of incidence is
!                 changed.
! 08.03.11 (BTD): v7.0.5
!                 * added ALPHA to argument list
!                 * eliminated variable BETA
!                 * replaced exp(-beta*(k*r)^4) with exp(-(alpha*k*r)^4)
! 08.04.20 (BTD): * change notation: ALPHA -> GAMMA
! 08.05.12 (BTD): v7.0.6
!                 * added OpenMP directives as suggested by Art Lazanoff, eg.
!                      #ifdef eself_omp
!                      !$omp parallel do
!                      #endif
!                   which will be compiled when the preprocessor flag
!                      -fpp -Deself_omp
!                   is used.
! 08.06.05 (ASL,BTD): eself_v3
!                 * added call to new routine CXFFT3_MKL to use Intel MKL
!                   library for FFT.
!                   This will now be invoked with CMETHD='FFTMKL'
! 08.06.27 (BTD) : eself_v4
!                 * BTD learning the ropes...
! 08.08.07 (BTD): removed #ifdef debug / #endif statements inserted
!                 around !omp directives by Art Lazanoff for testing.
!                 omp directives appear to work properly, they can
!                 continued to be tested as a whole by compiling with or 
!                 without the -openmp flag, and it is desirable for
!                 export code to be able to compile without cpp insofar
!                 as possible
! 11.07.03 (PFJ,BTD): eself_v6
!                 * removed local initialization of AK2OLD,AK3OLD,WOLD
!                   so that these values are now passed from calling
!                   program through module DDCOMMON_0
!                 * NGRID is now passed through module DDCOMMON_0 so
!                   that NGRID may be provided to other routines
! 12.07.06 (BTD): edited comments
! 12.08.02 (IYW): v7.3 : eself_v7
!                 * added DIPINT to arg list
!		  * introduced possibility of calculating Green's coefficients
!		    using "Filtered Coupled Dipole" (FCD) method (cf. ADDA)
!		  * added subroutine CISI for calculation of sine and cosine 
!                   integrals
! 12.08.11 (BTD): * removed DIPINT from arg list of ESELF
!                 * removed DIPINT from arg list of DIRECT_CALC
!                 * added new variable IDIPINT from module DDCOMMON_0
! end history
! Copyright (C) 1993,1994,1996,1997,1999,2000,2003,2004,2005,2006,2007,
!               2008,2011,2012 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

!--------------------------------------------

!*** diagnostic
!      write(0,*)'eself_v7 ckpt 1 with'
!      write(0,*)'       cmethd=',cmethd
!      write(0,*)'     nx,ny,nz=',nx,ny,nz
!      write(0,*)'      idipint=',idipint
!      write(0,*)'         ipbc=',ipbc
!      write(0,*)'        gamma=',gamma
!      write(0,*)'      pyd,pzd=',pyd,pzd
!      write(0,*)'      ak(1-3)=',ak
!      write(0,*)'          akd=',akd
!      write(0,*)'      dx(1-3)=',dx
!      write(0,*)'         wold=',wold
!      write(0,*)'       ak2old=',ak2old
!      write(0,*)'       ak3old=',ak3old
!      write(0,*)'  j jx jy jz k   cxzp(jx,jy,jz,k)'
!      do k=1,3
!         do jz=1,nz
!            do jy=1,ny
!               do jx=1,nx
!                  i=nz*ny*nx*(k-1)+ny*nx*(jz-1)+nx*(jy-1)+jx
!                  write(0,fmt='(i4,i3,i3,i3,i2,1p2e11.3)') &
!                               i,jx,jy,jz,k,cxzp(jx,jy,jz,k)
!               enddo
!            enddo
!         enddo
!      enddo
!***

! check if we can skip recomputation of Green-function coefficients

      IF(PYD.EQ.0._WP.AND.PZD.EQ.0._WP)THEN
         IF(ABS(WOLD-AKD)<1.E-6_WP*AKD)GOTO 70
      ELSEIF(PYD.NE.0._WP.AND.PZD.EQ.0._WP)THEN
         IF(ABS(WOLD-AKD)<1.E-6_WP*AKD.AND.   &
            ABS(AK(2)-AK2OLD)<1.E-6_WP*AKD)GOTO 70
      ELSEIF(PYD.EQ.0._WP.AND.PZD.NE.0._WP)THEN
         IF(ABS(WOLD-AKD)<1.E-6_WP*AKD.AND.   &
            ABS(AK(3)-AK3OLD)<1.E-6_WP*AKD)GOTO 70
      ELSEIF(PYD.NE.0._WP.AND.PZD.NE.0._WP)THEN
         IF(ABS(WOLD-AKD)<1.E-6_WP*AKD.AND.   &
            ABS(AK(2)-AK2OLD)<1.E-6_WP*AKD.AND. &
            ABS(AK(3)-AK3OLD)<1.E-6_WP*AKD)GOTO 70
      ENDIF

!*** diagnostic
!      write(0,*)'eself_v7 ckpt 2: recompute Green-function coefficients'
!***
! AKD.NE.WOLD :

! We have to recompute the Green-function coefficients giving
! components of field strength at a each grid point R
! produced by unit-valued component of dipole moment at
! point R', and then Fourier transform these components.

      WOLD=AKD
      AK2OLD=AK(2)
      AK3OLD=AK(3)
      NGRID=8*NX*NY*NZ
      AKD2=AKD*AKD

! We assume screening function exp(-(gamma*kr)^4) so
! range/d = 1/(gamma*kd) = 400 if gamma=5e-3 and kd=0.5
! although the sums are actually continued out to
! r/d = 2/(gamma*kd) = 800 if gamma=5e-3, kd=0.5
! [screening factor = exp(-16)=1.1e-7]

! PYDDX=PYD*DX(2) = periodicity in Y direction
! PZDDX=PZD*DX(3) = periodicity in Z direction

      IF(PYD>0._WP.OR.PZD>0._WP)THEN
         WRITE(CMSGNM,FMT='(A,2F8.2,A,1PE9.2)')'PBC with PYD, PZD=',PYD, &
                                               PZD,', GAMMA=',GAMMA
         CALL WRIMSG('ESELF ',CMSGNM)
      ENDIF

      PYDDX=PYD*DX(2)
      PZDDX=PZD*DX(3)

!               Compute the actual coefficients:

! Compute 6 independent elements of 3x3 symmetric matrix A_jk,where
! A_jk*P_k = -electric field at location j due to dipole P at location k

! A_jk = (a_1  a_2  a_3)
!        (a_2  a_4  a_5)
!        (a_3  a_5  a_6)_jk

      IF(IPBC==0)THEN

! initialize CXZC(I,J,K,M) = a_M for electric field at (I,J,K)
!                            produced by a dipole at (1,1,1)
!                            and replica dipoles (if PYD or PYZ are
!                            nonzero).

! need to calculate this for all (I,J,K) for one octant:
! I running from 1 to NX, J from 1 to NY, K from 1 to NZ

! Later obtain A_jk values for other octants by using symmetry

!*** diagnostic
!         write(0,*)'eself_v7 ckpt 3: call DIRECT_CALC'
!***
        CALL DIRECT_CALC(1,1,1,NX,NY,NZ,IPBC,DX,AK,AKD,AKD2,GAMMA, &
                         PYDDX,PZDDX,CXZC(1,1,1,1))
!*** diagnostic
!        write(0,*)'eself_v7 ckpt 4'
!        write(0,*)'   returned from direct_calc'
!        write(0,*)'   check for NaN...'
!        jr=0
!        do i=1,nx
!           do j=1,ny
!              do k=1,nz
!                 do m=1,6
!                    if(.not.(abs(cxzc(i,j,k,m))>=0.d0).or. &
!                        abs(cxzc(i,j,k,m))>=1.d100)then
!                       write(0,*)'i,j,k,m=',i,j,k,m,'cxzc=',cxzc(i,j,k,m)
!                       jr=jr+1
!                    endif
!                 enddo
!              enddo
!           enddo
!        enddo
!        write(0,*)'cxy checked for NaN or overflow',jr,' instances found'
!***

! At this point, CXZC(I,J,K,1-6) contains the upper triangular part of t
! symmetric 3 x 3 matrix giving the electric field at grid point (i,j,k)
! produced by a dipole at (1,1,1)

! Fill out CXZC to twice the size in each grid dimension to accomodate
! negative lags [periodicity in each dimension is assumed, so (e.g.)
! nx < i <= 2*nx is equivalent to -nx < i <= 0], exploiting symmetries,
! and then Fourier transform.

! If PYDDX=0 and PZDDX=0 , need only do direct calculation of A matrix
! for first octant, since remaining octants can be obtained by symmetry.
! After calculating A matrix, store only the first octant of the
! transform, since the others can be obtained by symmetry.

!-----------------------------------------------------------------------
! extend a_1 = a_xx(x,y,z) : a -> +a for x -> -x
!                                 +a     y -> -y
!                                 +a     z -> -z
        ISYM(1)=1
        ISYM(2)=1
        ISYM(3)=1
!*** diagnostic
!        write(0,*)'eself_v7 ckpt 5, about to call EXTEND'
!***
        CALL EXTND(CXZC(1,1,1,1),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
!*** diagnostic
!        write(0,*)'eself_v7 ckpt 6, returned from EXTEND'
!***
        IF(CMETHD=='GPFAFT')THEN
!*** diagnostic
!           write(0,*)'eself_v7 ckpt 7'
!***
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
!*** diagnostic
!           write(0,*)'eself_v7 ckpt 8'
!**
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF
!*** diagnostic
!        write(0,*)'eself_v7 ckpt 9, about to call TRIM'
!***
        CALL TRIM(CXZW(1,1,1,1),NX,NY,NZ,CXZC(1,1,1,1))
!*** diagnostic
!        write(0,*)'returned from TRIM'
!***
!-----------------------------------------------------------------------
! extend a_2 = a_xy(x,y,z) : a -> -a for x -> -x
!                                 -a     y -> -y
!                                 +a     z -> -z

        ISYM(1)=-1
        ISYM(2)=-1
        ISYM(3)=1
!*** diagnostic
!        write(0,*)'eself_v7 ckpt 10, about to call EXTND'
!***
        CALL EXTND(CXZC(1,1,1,2),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
        IF(CMETHD=='GPFAFT')THEN
!*** diagnostic
!           write(0,*)'eself_v7 ckpt 11, about to call cxfft3n'
!***
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
!*** diagnostic
!          write(0,*)'eself_v7 ckpt 12, returned from cxfft3n'
!***
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF
!*** diagnostic
!        write(0,*)'eself_v7 ckpt 13, about to call TRIM'
!***
        CALL TRIM(CXZW(1,1,1,1),NX,NY,NZ,CXZC(1,1,1,2))
!*** diagnostic
!        write(0,*)'returned from TRIM'
!***
!-----------------------------------------------------------------------
! extend a_3 = a_xz(x,y,z) : a -> -a for x -> -x
!                                 +a     y -> -y
!                                 -a     z -> -z

        ISYM(1)=-1
        ISYM(2)=1
        ISYM(3)=-1
        CALL EXTND(CXZC(1,1,1,3),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
        IF(CMETHD=='GPFAFT')THEN
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF

        CALL TRIM(CXZW(1,1,1,1),NX,NY,NZ,CXZC(1,1,1,3))

!-----------------------------------------------------------------------
! extend a_4 = a_yy(x,y,z) : a -> +a for x -> -x
!                                 +a     y -> -y
!                                 +a     z -> -z

        ISYM(1)=1
        ISYM(2)=1
        ISYM(3)=1
        CALL EXTND(CXZC(1,1,1,4),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
        IF(CMETHD=='GPFAFT')THEN
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF

        CALL TRIM(CXZW(1,1,1,1),NX,NY,NZ,CXZC(1,1,1,4))

!-----------------------------------------------------------------------
! extend a_5 = a_yz(x,y,z) : a -> +a for x -> -x
!                                 -a     y -> -y
!                                 -a     z -> -z
        ISYM(1)=1
        ISYM(2)=-1
        ISYM(3)=-1
        CALL EXTND(CXZC(1,1,1,5),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
        IF(CMETHD=='GPFAFT')THEN
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF

        CALL TRIM(CXZW(1,1,1,1),NX,NY,NZ,CXZC(1,1,1,5))

!-----------------------------------------------------------------------
! extend a_6 = a_zz(x,y,z) : a -> +a for x -> -x
!                                 +a     y -> -y
!                                 +a     z -> -z
        ISYM(1)=1
        ISYM(2)=1
        ISYM(3)=1
        CALL EXTND(CXZC(1,1,1,6),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
        IF(CMETHD=='GPFAFT')THEN
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF

        CALL TRIM(CXZW(1,1,1,1),NX,NY,NZ,CXZC(1,1,1,6))

      ELSEIF(IPBC==1)THEN

! This point is reached when PYDDX or PZDDX are nonzero.
! When PBC are used for general direction of incident wave,
! all octants of A matrix require direct calculation: symmetries valid
! for single target no longer apply because of position-dependent phases
! of replica dipoles.

!*** diagnostic
!         write(0,*)'eself_v7 ckpt 14'
!*** 
        CALL DIRECT_CALC(-1,-1,-1,NX,NY,NZ,IPBC,DX,AK,AKD,AKD2,GAMMA, &
                         PYDDX,PZDDX,CXZC(1,1,1,1))

!*** diagnostic
!        write(0,*)'eself_v7 ckpt 15'
!***
! The array CXZC(1-2*NX,1-2*NY,1-2*NZ,1-6) of A matrix coefficients
! now covers all octants.

! Fourier transform the A matrix CXZC:

        DO M=1,6
           IF(CMETHD=='GPFAFT')THEN
              CALL CXFFT3N(CXZC(1,1,1,M),2*NX,2*NY,2*NZ,+1)
           ELSEIF(CMETHD=='FFTW21')THEN
              CALL CXFFTW(CXZC(1,1,1,M),2*NX,2*NY,2*NZ,+1)
           ELSEIF(CMETHD.EQ.'FFTMKL')THEN
              CALL CXFFT3_MKL(CXZC(1,1,1,M),2*NX,2*NY,2*NZ,+1)
           ENDIF
        ENDDO

!*** diagnostic
!        write(0,*)'eself_v7 ckpt 16'
!***

! CXZC now contains the full Fourier transform of the A convolution
! and should not be overwritten between calls to ESELF

      ENDIF
!      CALL TIMEIT('ESELF (first call)',DTIME)
!      CALL TIMEIT('ESELF',DTIME)
!-----------------------------------------------------------------------

! End of recomputation of Green-function coefficients

70    CONTINUE
!*** diagnostic
!      write(0,*)'eself_v7 ckpt 17'
!****
! Fourier transform the polarizations:

      DO M=1,3
         CALL PAD(CXZP(1,1,1,M),NX,NY,NZ,CXZW(1,1,1,M))

!*** diagnostic
!        write(0,*)'eself_v7 ckpt 18: returned from PAD: ', &
!                  'check cxzw for NaN or overflow...'
!        jr=0
!        do i=1,2*nx
!           do j=1,2*ny
!              do k=1,2*nz
!                 if(.not.(abs(cxzw(i,j,k,m))>=0.d0).or. &
!                     abs(cxzw(i,j,k,m))>=1.d100)then
!                    write(0,*)'i,j,k,m=',i,j,k,m,'cxzw=',cxzw(i,j,k,m)
!                    jr=jr+1
!                 endif
!              enddo
!           enddo
!        enddo
!        write(0,*)'eself_v7 ckpt 19: cxzw checked for NaN or overflow: ', &
!                  jr,' instances found'
!        if(jr>0)stop
!***

        IF(CMETHD=='GPFAFT')THEN
           CALL CXFFT3N(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,+1)
!*** diagnostic
!        write(0,*)'eself_v7 ckpt 20: returned from CXFFT3N: ', &
!                  'check cxzw for NaN or overflow...'
!        jr=0
!        do i=1,2*nx
!           do j=1,2*ny
!              do k=1,2*nz
!                 if(.not.(abs(cxzw(i,j,k,m))>=0.d0).or. &
!                     abs(cxzw(i,j,k,m))>=1.d100)then
!                    write(0,*)'i,j,k,m=',i,j,k,m,'cxzw=',cxzw(i,j,k,m)
!                    jr=jr+1
!                 endif
!              enddo
!           enddo
!        enddo
!        write(0,*)'eself_v7 ckpt 21: cxzw checked for NaN or overflow',jr, &
!                  ' instances found'
!        if(jr>0)stop
!***
         ELSEIF(CMETHD=='FFTW21')THEN
            CALL CXFFTW(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,+1)
         ELSEIF(CMETHD.EQ.'FFTMKL')THEN
            CALL CXFFT3_MKL(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,+1)
         ENDIF
      ENDDO
!*** diagnostic
!      write(0,*)'eself_v7 ckpt 22'
!***
!***********************************************************************

! Multiply by F.t. of Green function.

      IF(IPBC==0)THEN

!*** diagnostic
!         write(0,*)'eself_v7 ckpt 23, check cxzc(i,j,k,m=1,6)'
!        jr=0
!        do i=1,2*nx
!           do j=1,2*ny
!              do k=1,2*nz
!                 do m=1,6
!                    if(.not.(abs(cxzc(i,j,k,m))>=0.d0).or. &
!                        abs(cxzc(i,j,k,m))>=1.d100)then
!                       write(0,*)'i,j,k,m=',i,j,k,m,'cxzc=',cxzc(i,j,k,m)
!                       jr=jr+1
!                    endif
!                 enddo
!              enddo
!           enddo
!        enddo
!        write(0,*)'eself_v7 ckpt 24: cxzc(i,j,k,m=1-6) checked for NaN ', &
!                  'or overflow: ',jr,' instances found'
!        if(jr>0)stop
!        write(0,*)'eself_v7 ckpt 25: check cxzw(i,j,k,m=1-3)'
!        jr=0
!        do i=1,2*nx
!           do j=1,2*ny
!              do k=1,2*nz
!                 do m=1,3
!                    if(.not.(abs(cxzw(i,j,k,m))>=0.d0).or. &
!                        abs(cxzw(i,j,k,m))>=1.d100)then
!                       write(0,*)'i,j,k,m=',i,j,k,m,'cxzw=',cxzw(i,j,k,m)
!                       jr=jr+1
!                    endif
!                 enddo
!              enddo
!           enddo
!        enddo
!        write(0,*)'eself_v7 ckpt 26: cxzw(i,j,k,m=1-3) checked for NaN or ', &
!                  'overflow: ',jr,' instances found'
!        if(jr>0)stop
!***
!***

! If IPBC=0, then only one octant of F.t. of Green function has been
!            stored, but can recover others using symmetry.

#ifdef openmp
!$OMP PARALLEL DO                                            &
!$OMP&   PRIVATE(K,J,I,KSGN,KR,JSGN,JR,ISGN,IR)              &
!$OMP&   PRIVATE(CXXX,CXXY,CXXZ,CXYY,CXYZ,CXZZ,CXEX,CXEY,CXEZ)
#endif

        DO K=1,2*NZ
           KSGN=NINT(SIGN(1._WP,NZ+1.5_WP-K))
           KR=MIN(K,2*NZ+2-K)
!*** diagnostic
!         write(0,*)'eself_v7 ckpt 27'
!***
           DO J=1,2*NY
              JSGN=NINT(SIGN(1._WP,NY+1.5_WP-J))
              JR=MIN(J,2*NY+2-J)
!*** diagnostic
!         write(0,*)'eself_v7 ckpt 28'
!***
              DO I=1,2*NX
                 ISGN=NINT(SIGN(1._WP,NX+1.5_WP-I))
                 IR=MIN(I,2*NX+2-I)
!*** diagnostic
!         write(0,*)'eself_v7 ckpt 29'
!***
                 CXXX=CXZC(IR,JR,KR,1)
                 CXXY=CXZC(IR,JR,KR,2)*(ISGN*JSGN)
                 CXXZ=CXZC(IR,JR,KR,3)*(ISGN*KSGN)
                 CXYY=CXZC(IR,JR,KR,4)
                 CXYZ=CXZC(IR,JR,KR,5)*(JSGN*KSGN)
                 CXZZ=CXZC(IR,JR,KR,6)
!*** diagnostic
!              if(.not.(abs(cxxx)>=0.d0))write(0,*) &
!                'ir,jr,kr,cxzc(ir,jr,kr,1)=',ir,jr,kr,cxzc(ir,jr,kr,1)
!              if(.not.(abs(cxxy)>=0.d0))write(0,*) &
!                'ir,jr,kr,cxzc(ir,jr,kr,2)=',ir,jr,kr,cxzc(ir,jr,kr,2)
!              if(.not.(abs(cxxz)>=0.d0))write(0,*) &
!                 'ir,jr,kr,cxzc(ir,jr,kr,3)=',ir,jr,kr,cxzc(ir,jr,kr,3)
!              if(.not.(abs(cxyy)>=0.d0))write(0,*) &
!                  'ir,jr,kr,cxzc(ir,jr,kr,4)=',ir,jr,kr,cxzc(ir,jr,kr,4)
!              if(.not.(abs(cxyz)>=0.d0))write(0,*) &
!                  'ir,jr,kr,cxzc(ir,jr,kr,5)=',ir,jr,kr,cxzc(ir,jr,kr,5)
!              if(.not.(abs(cxzz)>=0.d0))write(0,*) &
!                  'ir,jr,kr,cxzc(ir,jr,kr,6)=',ir,jr,kr,cxzc(ir,jr,kr,6)
!***
                 CXEX=CXXX*CXZW(I,J,K,1)+CXXY*CXZW(I,J,K,2)+CXXZ*CXZW(I,J,K,3)
                 CXEY=CXXY*CXZW(I,J,K,1)+CXYY*CXZW(I,J,K,2)+CXYZ*CXZW(I,J,K,3)
                 CXEZ=CXXZ*CXZW(I,J,K,1)+CXYZ*CXZW(I,J,K,2)+CXZZ*CXZW(I,J,K,3)
                      
!*** diagnostic
!              if(.not.(abs(cxex+cxey+cxez)>=0.d0))then
!                 write(0,*)'i,j,j,ir,jr,kr=',i,j,k,ir,jr,kr
!                 write(0,*)'cxzc(ir,jr,kr,1)=',cxzc(ir,jr,kr,1)
!                 write(0,*)'cxzc(ir,jr,kr,2)=',cxzc(ir,jr,kr,2)
!                 write(0,*)'cxzc(ir,jr,kr,3)=',cxzc(ir,jr,kr,3)
!                 write(0,*)'cxzc(ir,jr,kr,4)=',cxzc(ir,jr,kr,4)
!                 write(0,*)'cxzc(ir,jr,kr,5)=',cxzc(ir,jr,kr,5)
!                 write(0,*)'cxzc(ir,jr,kr,6)=',cxzc(ir,jr,kr,6)
!                 write(0,*)'      cxex=',cxex
!                 write(0,*)'      cxey=',cxey
!                 write(0,*)'      cxez=',cxez
!                 write(0,*)'      cxxx=',cxxx
!                 write(0,*)'      cxxy=',cxxy
!                 write(0,*)'      cxxz=',cxxz
!                 write(0,*)'      cxyy=',cxyy
!                 write(0,*)'      cxyz=',cxyz
!                 write(0,*)'      cxzz=',cxzz
!                 write(0,*)'      cxzw(i,j,k,1)=',cxzw(i,j,k,1)
!                 write(0,*)'      cxzw(i,j,k,2)=',cxzw(i,j,k,2)
!                 write(0,*)'      cxzw(i,j,k,3)=',cxzw(i,j,k,3)
!                 stop
!              endif
!***
                 CXZW(I,J,K,1)=CXEX
                 CXZW(I,J,K,2)=CXEY
                 CXZW(I,J,K,3)=CXEZ
!*** diagnostic
!              if(.not.(abs(cxzw(i,j,k,1))>=0.d0))write(0,*) &
!                 'i,j,k,cxzw(i,j,k,1)=',i,j,k,cxzw(i,j,k,1)
!              if(.not.(abs(cxzw(i,j,k,2))>=0.d0))write(0,*) &
!                 'i,j,k,cxzw(i,j,k,2)=',i,j,k,cxzw(i,j,k,2)
!              if(.not.(abs(cxzw(i,j,k,3))>=0.d0))write(0,*) &
!                 'i,j,k,cxzw(i,j,k,3)=',i,j,k,cxzw(i,j,k,3)
!***
               ENDDO
            ENDDO
         ENDDO

!*** diagnostic
!         write(0,*)'eself_v7 ckpt 30'
!***

#ifdef openmp
!$OMP END PARALLEL DO
#endif

      ELSEIF(IPBC==1)THEN

!*** diagnostic
!         write(0,*)'eself_v7 ckpt 31'
!***

! If IPBC=1, then the full F.t. of the Green function has been stored.

#ifdef openmp
!$OMP PARALLEL DO                                                  &
!$OMP&   PRIVATE(K,J,I,CXXX,CXXY,CXXZ,CXYY,CXYZ,CXZZ,CXEX,CXEY,CXEZ)
#endif

         DO K=1,2*NZ
            DO J=1,2*NY
               DO I=1,2*NX
                  CXXX=CXZC(I,J,K,1)
                  CXXY=CXZC(I,J,K,2)
                  CXXZ=CXZC(I,J,K,3)
                  CXYY=CXZC(I,J,K,4)
                  CXYZ=CXZC(I,J,K,5)
                  CXZZ=CXZC(I,J,K,6)
                  CXEX=CXXX*CXZW(I,J,K,1)+CXXY*CXZW(I,J,K,2)+CXXZ*CXZW(I,J,K,3)
                  CXEY=CXXY*CXZW(I,J,K,1)+CXYY*CXZW(I,J,K,2)+CXYZ*CXZW(I,J,K,3)
                  CXEZ=CXXZ*CXZW(I,J,K,1)+CXYZ*CXZW(I,J,K,2)+CXZZ*CXZW(I,J,K,3)
                  CXZW(I,J,K,1)=CXEX
                  CXZW(I,J,K,2)=CXEY
                  CXZW(I,J,K,3)=CXEZ
               ENDDO
            ENDDO
         ENDDO

#ifdef openmp
!$OMP END PARALLEL DO
#endif

      ENDIF

! Inverse Fourier transform to obtain electric field:

      DO M=1,3
         IF(CMETHD=='GPFAFT')THEN
            CALL CXFFT3N(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,-1)
         ELSEIF(CMETHD=='FFTW21')THEN
            CALL CXFFTW(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,-1)
         ELSEIF(CMETHD.EQ.'FFTMKL')THEN
            CALL CXFFT3_MKL(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,-1)
         ENDIF

!*** diagnostic
!         write(0,*)'eself_v7 ckpt 32'
!***
!***********************************************************************

! Note: the Convex FFT routine already normalizes result.
!       For other FFT routines need to divide result by NGRID

         IF(CMETHD=='CONVEX')THEN
            DO K=1,NZ
               DO J=1,NY
                  DO I=1,NX
                     CXZE(I,J,K,M)=CXZW(I,J,K,M)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            FAC=1._WP/REAL(NGRID,KIND=WP)

#ifdef openmp
!$OMP PARALLEL DO     &
!$OMP&   PRIVATE(I,J,K)
#endif

            DO K=1,NZ
               DO J=1,NY
                  DO I=1,NX
                     CXZE(I,J,K,M)=FAC*CXZW(I,J,K,M)
                  ENDDO
               ENDDO


            ENDDO

#ifdef openmp
!$OMP END PARALLEL DO
#endif

         ENDIF
      ENDDO
      RETURN
    END SUBROUTINE ESELF

!***********************************************************************

    SUBROUTINE EXTND(CXA,NX,NY,NZ,ISYM,CXB)
      USE DDPRECISION,ONLY: WP

#ifdef openmp
      USE OMP_LIB   !Art for OpenMP function declarations
#endif

      IMPLICIT NONE

!       Using symmetries, extend first octant of coefficient matrix to
!       other octants.

! Originally written by Jeremy Goodman
! 080620 (ASL) modified to use OpenMP
! 080627 (BTD) eself_v4: further mods related to OpenMP
!              reordered several nested do loops
!              introduced variables J2 and K2 to speed computations
! end history
! Copyright (C) 1993, B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
! Arguments:

      INTEGER :: NX,NY,NZ
      COMPLEX(WP) :: CXA(NX+1,NY+1,NZ+1),CXB(2*NX,2*NY,2*NZ)
      INTEGER :: ISYM(3)

! Local variables:

      INTEGER :: I,J,J2,K,K2
      COMPLEX(WP) :: CXZERO

! SAVE statement:

      SAVE CXZERO

! DATA statement:

      DATA CXZERO/(0._WP,0._WP)/
!***********************************************************************

#ifdef openmp
!$OMP PARALLEL              & 
!$OMP&   PRIVATE(I,J,J2,K,K2)
!$OMP DO
#endif

      DO K=1,NZ

         DO J=1,NY
            DO I=1,NX
               CXB(I,J,K)=CXA(I,J,K)
            ENDDO
         ENDDO

      ENDDO

#ifdef openmp
!$OMP ENDDO
#endif

!-----------------------------------------------------------------------
! x -> -x

!btd 080627 moved I=NX+1 out of loop
! the SINGLE directive specifies that enclosed code is to be executed
! by only one thread in the team

#ifdef openmp
!$OMP SINGLE
#endif

      CXB(NX+1,1:NY,1:NZ)=CXZERO

#ifdef openmp
!$OMP END SINGLE
!$OMP DO 
#endif

      DO K=1,NZ

         DO J=1,NY
            DO I=NX+2,2*NX
               CXB(I,J,K)=CXA(2*NX+2-I,J,K)*ISYM(1)
            ENDDO
         ENDDO

      ENDDO

#ifdef openmp
!$OMP END DO
#endif

!-----------------------------------------------------------------------
! y -> -y

!btd 080627 moved J=NY+1 out of loop, switched order of loops I and J
! the SINGLE directive specifies that enclosed code is to be executed
! by only one thread in the team

#ifdef openmp
!$OMP SINGLE
#endif

      CXB(1:2*NX,NY+1,1:NZ)=CXZERO

#ifdef openmp
!$OMP END SINGLE
!$OMP DO 
#endif

      DO K=1,NZ

         DO J=NY+2,2*NY
            J2=2*NY+2-J
            DO I=1,2*NX
               CXB(I,J,K)=CXB(I,J2,K)*ISYM(2)
            ENDDO
         ENDDO

      ENDDO

#ifdef openmp
!$OMP END DO
#endif

!-----------------------------------------------------------------------
! z -> -z

!Art pulling the 3rd dimension to the outer most loop.
!Art we'll do this expression in only the thread that has NZ+1
!btd 080627 reorder loops: J,I,K -> K,J,I
! the SINGLE directive specifies that enclosed code is to be executed
! by only one thread in the team

#ifdef openmp
!$OMP SINGLE
#endif

      CXB(1:2*NX,1:2*NY,NZ+1)=CXZERO

#ifdef openmp
!$OMP END SINGLE
!$OMP DO 
#endif

      DO K=NZ+2,2*NZ

         K2=2*NZ+2-K
         DO J=1,2*NY
            DO I=1,2*NX
               CXB(I,J,K)=CXB(I,J,K2)*ISYM(3)
            ENDDO
         ENDDO

      ENDDO

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif

      RETURN
    END SUBROUTINE EXTND
!***********************************************************************

    SUBROUTINE TRIM(CXB,NX,NY,NZ,CXA)
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE

! Copy the first octant of CXB into CXA

! Originally written by Jeremy Goodman
! History:
! 92.04.20 (BTD) removed ISYM from argument list of TRIM:
! 06.09.28 (BTD) eself v2.0
! 08.06.27 (BTD) eself_v4
!                
! Copyright (C) 1993,2006 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
! Arguments:

      INTEGER :: NX,NY,NZ
      COMPLEX(WP) ::          &
         CXA(NX+1,NY+1,NZ+1), &
         CXB(2*NX,2*NY,2*NZ)

! Local variables:

      INTEGER :: I,J,K

!***********************************************************************

#ifdef openmp
!$OMP PARALLEL             &
!$OMP&    PRIVATE(I,J,K)   &
!$OMP&    SHARED(NX,NY,NZ) &
!$OMP&    SHARED(CXA,CXB)
!$OMP DO
#endif

      DO K=1,NZ+1

         DO J=1,NY+1
            DO I=1,NX+1
               CXA(I,J,K)=CXB(I,J,K)
            ENDDO
         ENDDO

      ENDDO

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif

      RETURN
    END SUBROUTINE TRIM

!***********************************************************************

    SUBROUTINE PAD(CXA,NX,NY,NZ,CXB)
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE

! Pad the array CXA with zeros out to twice its length in each
! dimension, and put the result in CXB.

! Originally written by Jeremy Goodman

! Copyright (C) 1993, B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
! Arguments:
      INTEGER :: NX,NY,NZ
      COMPLEX(WP) ::       &
         CXA(NX,NY,NZ),    &
         CXB(2*NX,2*NY,2*NZ)
! Local variables:
      COMPLEX(WP) :: CXZERO
      INTEGER :: I,J,K
      SAVE CXZERO
      DATA CXZERO/(0._WP,0._WP)/
!***********************************************************************

#ifdef openmp
!$OMP PARALLEL        &
!$OMP&   PRIVATE(I,J,K)
!$OMP DO
#endif

      DO K=1,2*NZ

         DO J=1,2*NY
            DO I=1,2*NX
               CXB(I,J,K)=CXZERO
            ENDDO

         ENDDO
      ENDDO

#ifdef openmp
!$OMP END DO
!$OMP DO
#endif

      DO K=1,NZ

         DO J=1,NY
            DO I=1,NX
               CXB(I,J,K)=CXA(I,J,K)
            ENDDO
         ENDDO

      ENDDO

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif

      RETURN
    END SUBROUTINE PAD

!-----------------------------------------------------------------------

    SUBROUTINE DIRECT_CALC(IX,IY,IZ,NX,NY,NZ,IPBC,DX,AK,AKD,AKD2,GAMMA, &
                           PYDDX,PZDDX,CXZC)
      USE DDPRECISION,ONLY: WP
      USE DDCOMMON_0,ONLY: IDIPINT
      USE DDCOMMON_10,ONLY: MYID
#ifdef openmp
      USE OMP_LIB   !! Art omp v2.0 supplies variable and function definitions
#endif
      IMPLICIT NONE

! arguments:

      INTEGER :: IPBC,IX,IY,IZ,NX,NY,NZ
      REAL(WP) :: AKD,AKD2,GAMMA,PYDDX,PZDDX
      REAL(WP) :: &
         AK(3),   &
         DX(3)
      COMPLEX (WP) ::                                             &
         CXZC(NX+1+IPBC*(NX-1),NY+1+IPBC*(NY-1),NZ+1+IPBC*(NZ-1),6)

! local variables

      CHARACTER :: CMSGNM*66
      INTEGER :: I,ICXZC,II,IMIN,J,JCXZC,JJ,JMIN,JPY,JPYM,JPZ,JPZM, &
                 K,KCXZC,KMIN,M
      REAL(WP) :: CIMINUS,CIPLUS,COSKFR,COSKR,GAMMAKD4,KF,KFR,KR,           &
                  PHASY,PHASYZ,PI,R,R2,R3,RANGE,RANGE2,RJPY,RJPZ,           &
                  SIMINUS,SIPLUS,SINKFR,SINKR,T0,T1,T2,X0,X2,X2Y2,XX,Y0,Y2,Z0
      REAL(WP) :: X(3)
      COMPLEX(WP) :: CX2PIH,CXA0,CXA1,CXA2,CXEXPIKR,CXFAC, &
                     CXG0,CXG1,CXG2,CXI,CXIKR,CXPHAS,CXZERO
                     
      COMPLEX(WP) :: DCXSUM(6)

#ifdef openmp
      INTEGER NTHREADS,TID
#endif

      SAVE CXZERO,CXI
      DATA CXI/(0._WP,1._WP)/,CXZERO/(0._WP,0._WP)/

!-----------------------------------------------------------------------
! subroutine DIRECT_CALC

! calculates electric field at (I,J,K) due to dipole at (1,1,1)
!                                      plus its replicas
! given:
!        IX,IY,IZ = 1, 1, 1 to do octant I>0, J>0, K>0 only
!                   1, 1,-1              I>0, J>0, all K
!                   1,-1, 1              I>0, K>0, all J
!                   1,-1,-1              I>0, all J, all K
!                  -1, 1, 1              all I, J>0, K>0
!                  -1, 1,-1              all I, J>0, all K
!                  -1,-1, 1              all I, all J, K>0
!                  -1,-1,-1              all I, all J, all K

!        NX,NY,NZ = size of first octant (I = 1 -> NX,
!                                         J = 1 -> NY,
!                                         K = 1 -> NZ)

!        IPBC     = 0 for isolated target
!                 = 1 to use periodic boundary conditions
!                     (y direction, z direction, or both y and z directions)
!                   N.B.: IPBC affects dimensioning of CXZC

!        GAMMA     = real coefficient used to assist convergence of sums
!                    over replica dipoles by suppressing long-range
!                    contributions with factor exp(-gamma*(kr)^2)
!                    typical value gamma = 0.005
!                    The effective
!                    range/d = 1/(gamma*k*d) = 400 if gamma=5e-3 and kd=0.5
!                    range/lambda = 1/(2*pi*gamma) = 31.8
!                    The sums are actually continued out to
!                    r/d = 2*/(gamma*kd) = 800 if gamma=5e-3 , kd=0.5
!                    [screening factor = exp(-16)=1.1e-7]
!
!        PYDDX     = period of lattice in y direction/d
!        PZDDX     = period of lattice in z direction/d
!        DX(1-3)   = lattice spacing in x,y,z direction, in units of
!                    d = n**(-1/3).
!        AK(1-3)   = k(1-3)*d, where k = k vector in vacuo
!        AKD       = |k|d
!        CXZC      = array with dimensions
!                    (NX+1)*(NY+1)*(NZ+1)*6 if IPBC=0
!                    (2*NX)*(2*NY)*(2*NZ)*6 if IPBC=1
!
! and through module DDCOMMON_0:
!
!	 IDIPINT   = 0 for point dipole interaction method
!		   = 1 for filtered coupled dipole (FCD) interaction method 

! returns:
!        CXZC(I,J,K,M) = a_M to calculate E field at (I,J,K)
!                        produced by dipole P at (1,1,1) and replica dipoles
!                        if IPBC > 0
!                        -E_x = a_1*P_x + a_2*P_y + a_3*P_z
!                        -E_y = a_2*P_x + a_4*P_y + a_5*P_z
!                        -E_z = a_3*P_x + a_5*P_y + a_6*P_z 

! JPYM = maximum number of periods in Y direction
!        rmax = JPYM*PYDDX
! note that even for short-range interaction, need to extend sums
! from JPY=-1 to JPY=+1 to include "edge" interactions with next TUC

! PYD=0. for nonperiodic case
!    >0. for periodic boundary conditions in y direction

! Adapted from original subroutine ESELF written by Jeremy Goodman
! history:
! 06.09.28 (BTD) DIRECT_CALC works for isolated target
! 06.09.30 (BTD) further modifications to DIRECT_CALC
!                change definition of BETA, so that given BETA 
!                determines range/wavelength
! 08.03.11 (BTD) v7.0.5
!                * eliminate BETA, introduce ALPHA=BETA**(0.25)
! 08.03.15 (BTD) * added code to estimate time to completion
!                  may need to change way this is done before running under MPI
! 08.03.15 (BTD) * added DDCOMMON_10 to communicate value of MYID.
! 08.04.20 (BTD) * changed notation: ALPHA -> GAMMA
! 08.06.05 (ASL) v7.0.6: 
!                * Arthur S. Lazanoff added OpenMP directives
! 08.07.02 (BTD) * Added call to CPU_TIME(T0) outside of PARALLEL region
! 12.04.16 (BTD) v7.2:
!                * add check to verify number of OpenMP threads               
! 12.07.06 (BTD) v7.2.3 edited comments
! 12.08.02 (IYW) added DIPINT to arg list and FCD Green's coefficients
! 12.08.10 (BTD) **** need to examine how these changes affect OpenMP ***
! 12.08.11 (BTD) v7.3 and eself_v7
!                * removed DIPINT from argument list
!                * added IDIPINT from module DDCOMMON_0
! 12.12.21 (BTD) * cosmetic changes, added comments
! 12.12.27 (BTD) * corrected errors in case IDPINT=1 (FILTRDD)
!                * added comments relating code to notation of
!                  Gay-Balmaz & Martin (2002) [Comp. Phys. Comm. 144, 111] 
!                  and Yurkin, Min & Hoekstra (2010) [PRE 82, 036703]
! 12.12.29 (BTD) * added phase correction for DIPINT=1
!                  (needed for PBC calculations)
! 13.01.05 (BTD) * OpenMP bug fix: add ICXZC to list of PRIVATE variables
!                  change
!                  PRIVATE(JCXZC,JPZM,KCXZC)               &
!                  to
!                  PRIVATE(ICXZC,JCXZC,JPZM,KCXZC)         &
! 13.01.07 (BTD) * corrected typo PRIVATE(ICXCZ -> PRIVATE(ICXZC
!                  (noted 13.01.07 by Zhenpeng Qin)
! end history
! Copyright (C) 1993,1994,1996,1997,1999,2000,2003,2004,2005,2006,2007,
!               2008,2011,2012,2013 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!-----------------------------------------------------------------------
! PI and KF are used if IDIPINT=1 (filtered coupled dipole)

      PI=4.D0*DATAN(1.D0)
      KF=PI

      GAMMAKD4=(GAMMA*AKD)**4
      RANGE=2._WP/(GAMMA*AKD)
      RANGE2=RANGE*RANGE
      IF(PYDDX<=0._WP)THEN
         JPYM=0
      ELSE
         JPYM=1+NINT(RANGE/PYDDX)
      ENDIF

!*** diagnostic
!      write(0,*)'eself_v7 direct_calc ckpt 1'
!      write(0,*)'jpym=',jpym
!***
! Compute 6 independent elements of 3x3 symmetric matrix A_jk, where
! A_jk*P_k = -electric field at location j due to dipole P_k at location

! A_jk = (a_1  a_2  a_3)
!        (a_2  a_4  a_5)
!        (a_3  a_5  a_6)_jk

! initialize CX(I,J,K,M) = a_M for electric field at (I,J,K)
!                          produced by a dipole at (1,1,1)
!                          and replica dipoles (if PYD or PYZ are
!                          nonzero).

! need to calculate this for all (I,J,K) for one octant:
! I running from 1 to NX, J from 1 to NY, K from 1 to NZ

! X0,Y0,Z0 = X(I,J,K) - X(1,1,1) = vector from dipole location (1,1,1)
!                                  to point where E is to be calculated
! IX = +1 -> IMIN = 1
! IX = -1 -> IMIN = 2-NX   (1-IMIN = NX-1)

      IMIN=1+(1-NX)*(1-IX)/2
      JMIN=1+(1-NY)*(1-IY)/2
      KMIN=1+(1-NZ)*(1-IZ)/2

!*** diagnostic
!      write(0,*)'eself_v7 direct_calc ckpt 2'
!***

! Determine elapsed cpu time for these sums

      CALL CPU_TIME(T0)

#ifdef openmp
! fork a team of threads 
!$OMP PARALLEL               &
!$OMP&   PRIVATE(NTHREADS,TID)

      TID=OMP_GET_THREAD_NUM()

!*** diagnostic
!      write(0,*)'eself_v7 direct_calc ckpt 3: hello world from thread = ',TID
!***
! only master thread does this:

      IF(TID.EQ.0)THEN
         NTHREADS=OMP_GET_NUM_THREADS()
         WRITE(CMSGNM,FMT='(A,I3)')'number of OpenMP threads =',NTHREADS
         CALL WRIMSG('DIRECT_CALC',CMSGNM)
      ENDIF
!$OMP END PARALLEL
#endif

! 2013.01.04 (BTD) Zhenpen Qin reports error within following parallelization
! 2013.01.05 (BTD) ICXZC was inadvertently omitted from PRIVATE variables list
!                  add it: change
!                  PRIVATE(JCXZC,JPZM,KCXZC) to
!                  PRIVATE(ICXZC,JCXZC,JPZM,KCXZC)
#ifdef openmp
! 2012.04.27 (BTD) added R to private variables
!                  removed FLUSH directive (should not be needed)
!$OMP PARALLEL                                       &
!$OMP&   PRIVATE(I,II,J,JJ,JPY,JPZ,K,M)              &
!$OMP&   PRIVATE(ICXZC,JCXZC,JPZM,KCXZC)             &
!$OMP&   PRIVATE(CIMINUS,CIPLUS,COSKFR,COSKR,KFR,KR) &
!$OMP&   PRIVATE(PHASY,PHASYZ,R,R2,R3,RJPY,RJPZ)     &
!$OMP&   PRIVATE(SIMINUS,SIPLUS,SINKFR,SINKR)        &
!$omp&   PRIVATE(X,X0,X2,X2Y2,XX,Y0,Y2,Z0)           &
!$OMP&   PRIVATE(CX2PIH,CXA0,CXA1,CXA2)              &
!$OMP&   PRIVATE(CXEXPIKR,CXFAC,CXG0,CXG1,CXG2)      &
!$OMP&   PRIVATE(CXIKR,CXPHAS,DCXSUM)
!$OMP DO
#endif

      DO K=KMIN,NZ              !- loop over K

         Z0=REAL(K-1,KIND=WP)*DX(3)
         IF(K>0)THEN
            KCXZC=K
         ELSE
            KCXZC=2*NZ+K
         ENDIF
         DO J=JMIN,NY              !-- loop over J
            Y0=REAL(J-1,KIND=WP)*DX(2)
            IF(J>0)THEN
               JCXZC=J
            ELSE
               JCXZC=2*NY+J
            ENDIF
            DO I=IMIN,NX              !--- loop over I

! for first dipole, determine time to sum over replicas

               IF(I.EQ.IMIN.AND.J.EQ.JMIN.AND.          &
                  K.EQ.KMIN.AND.MYID==0)CALL CPU_TIME(T1)

               X0=REAL(I-1,KIND=WP)*DX(1)
               IF(I>0)THEN
                  ICXZC=I
               ELSE
                  ICXZC=2*NX+I
               ENDIF
               X(1)=X0
               X2=X0**2
               DO M=1,6
                  DCXSUM(M)=CXZERO
               ENDDO

!*** diagnostic
!               write(0,*)'eself_v7 ckpt 4'
!               write(0,*)'   x2=',x2
!               write(0,*)'   jpym=',jpym
!***

! JPY=0, JPZ=0 gives E field from dipole at (1,1,1)
! general JPY,JPZ gives E field from dipole at
! (1,JPY*NPY+1,JPZ*NPZ+1)
! replica dipoles have same magnitude of polarization as dipole (1,1,1),
! but different phase.
! PHASYZ = phase of replica dipole - phase of dipole (1,1,1)

               DO JPY=-JPYM,JPYM         !---- loop over JPY
                  RJPY=REAL(JPY,KIND=WP)*PYDDX
                  X(2)=Y0-RJPY
                  Y2=X(2)*X(2)
                  X2Y2=X2+Y2
                  PHASY=AK(2)*RJPY

! PZDDX=0. for nonperiodic case
!      >0. for periodic boundary conditions in z direction

                  IF(PZDDX<=0._WP)THEN
                     JPZM=0
                  ELSE
                     JPZM=1+NINT(SQRT(MAX(RANGE2-RJPY**2,0._WP))/PZDDX)
                  ENDIF

!*** diagnostic
!      write(0,*)'jpy,jpzm=',jpy,jpzm
!***
                  DO JPZ=-JPZM,JPZM         !----- loop over JPZ
!*** diagnostic
!      write(0,*)'jpz=',jpz
!***

                     RJPZ=REAL(JPZ,KIND=WP)*PZDDX
                     X(3)=Z0-RJPZ
                     R2=X2Y2+X(3)*X(3)

! skip the self-interaction (R=0) case (I=J=K=1 and JPY=JPZ=0)

                     IF(R2>1.E-6_WP)THEN

                        R=SQRT(R2)
                        R3=R*R2

! PHASYZ = phase at (1,JPY*NPY+1,JPZ*NPZ+1) - phase at (1,1,1)

                        PHASYZ=PHASY+AK(3)*RJPZ
                        KR=AKD*R
                        CXIKR=CXI*KR

! include artificial factor exp[-(gamma*kr)^4] to assist convergence

			IF(IDIPINT==0)THEN

                           CXPHAS=EXP(CXIKR+CXI*PHASYZ-GAMMAKD4*R2*R2)/R3
                           CXFAC=(1._WP-CXIKR)/R2

!*** diagnostic
!                           write(0,*)'eself_v7 direct_calc ckpt 5'
!                           write(0,*)'  cxterm=',cxterm
!                           write(0,*)'  cxphas=',cxphas
!                           write(0,*)'   cxfac=',cxfac
!***
                           
! II=1      -> M=1   a_1 (xx)
!      JJ=2      2   a_2 (xy)
!         3      3   a_3 (xz)
!    2           4   a_4 (yy)
!         3      5   a_5 (yz)
!    3           6   a_6 (zz)

                           M=0
                           DO II=1,3              !------- loop over II
                              M=M+1
                              XX=X(II)**2
                              DCXSUM(M)=DCXSUM(M)-CXPHAS*(AKD2*(XX-R2)+ &
                                        CXFAC*(R2-3._WP*XX))
                              IF(II<3)THEN
                                 DO JJ=II+1,3        !------- loop over JJ
                                    M=M+1
                                    DCXSUM(M)=DCXSUM(M)-                          &
                                              CXPHAS*X(II)*X(JJ)*(AKD2-3._WP*CXFAC)
                                 ENDDO               !------- end loop over JJ
                              ENDIF
                           ENDDO                  !------ end loop over II

			ELSEIF(IDIPINT==1)THEN

! Expressions for the Green coefficients correspond to those in 
! "A library for computing the filtered and non-filtered 3D Green's tensor 
!  associated with infinite homogeneous space and surfaces" 
! by P. Gay-Balmaz and O. Martin (2002; Computer Physics Communications,
! 144, 111-120) apart from a factor of 4*pi.

                           CXEXPIKR=EXP(CXIKR)
                           CXPHAS=EXP(CXI*PHASYZ-GAMMAKD4*R2*R2)
                           KFR=KF*R
                           CALL CISI(KFR+KR,CIPLUS,SIPLUS)
                           CALL CISI(KFR-KR,CIMINUS,SIMINUS)
!*** diagnostic
!                           write(57,fmt='(1pe10.3,1p2e11.3)')kfr+kr,     &
!                                                             ciplus,siplus
!                           write(57,fmt='(1pe10.3,1p2e11.3)')kfr-kr,       &
!                                                             ciminus,siminus
!***
                           COSKR=COS(KR)
                           SINKR=SIN(KR)
                           COSKFR=COS(KFR)
                           SINKFR=SIN(KFR)

! CXA0 = alpha(R) of Gay-Balmaz & Martin 2002
! CXA1 = (d/dR) alpha   = alpha'
! CXA2 = (d2/dR2) alpha = alpha"

                           CXA0=SINKR*(CIPLUS-CIMINUS)+ &
                                COSKR*(PI-SIPLUS-SIMINUS)
                           CXA1=AKD*SINKR*(-PI+SIPLUS+SIMINUS)+         &
                                AKD*COSKR*(CIPLUS-CIMINUS)-2._WP*SINKFR/R
                           CXA2=AKD2*(SINKR*(CIMINUS-CIPLUS)+     &
                                      COSKR*(SIPLUS+SIMINUS-PI))+ &
                                2._WP*(SINKFR-KFR*COSKFR)/R2

! Yurkin, Min & Hoekstra 2010 define G_ij such that G_ij*P_j = E field at i
! so that elements of G_ij are same as elements of our matrix CXZC
! notation:
! CXG0   = g_F(R)              of YMH2010
! CXG1   = (d/dR) g_F   = g_F' of YMH2010
! CXG2   = (d2/dR2) g_F = g_F" of YMH2010  
! CX2PIH = 2*pi*h(R)           of YMH2010

! CXG0 = exp(ikR)/R - alpha/(pi*R)
! CXG1 = exp(ikR)*(-1+ikR)/R^2 - alpha'/(pi*R) + alpha/(pi*R^2)
! CXG2 = exp(ikR)*(2-2ikR-(kR)^2)/R^3 - (2*alpha-2R*alpha'+R^2*alpha")/(pi*R^3)

                           CXG0=(CXEXPIKR-CXA0/PI)/R
                           CXG1=(CXEXPIKR*(CXIKR-1._WP)+(CXA0-CXA1*R)/PI)/R2
                           CXG2=(CXEXPIKR*(2._WP-2._WP*CXIKR-KR**2)- &
                                 (2._WP*(CXA0-R*CXA1)+R2*CXA2)/PI)/R3
                           CX2PIH=(SINKFR-KFR*COSKFR)/(PI*R3)
                           M=0
                           DO II=1,3   ! do II
                              M=M+1
! diagonal elements 
! II=1: M=1 xx
! II=2: M=4 yy
! II=3: M=6 zz
                              DCXSUM(M)=DCXSUM(M)+CXPHAS*(AKD2*CXG0+ &
                                        CXG1/R+(2._WP/3._WP)*CX2PIH+ &
                                        (CXG2/R2-CXG1/R3)*X(II)**2)
                              IF(II<3)THEN
                                 DO JJ=II+1,3  ! do JJ
                                    M=M+1
! off-diagonal elements
! II=1: M=2 xy
! II=1: M=3 xz
! II=2: M=5 yz
                                    DCXSUM(M)=DCXSUM(M)+CXPHAS*           &
                                              (CXG2/R2-CXG1/R3)*X(II)*X(JJ)
                                 ENDDO  ! enddo JJ
                              ENDIF
                           ENDDO  ! enddo II
			ENDIF  ! end elseif(IDIPINT==1)

                     ENDIF  ! endif(R2 > 1e-6)
                  ENDDO                     !----- end loop over JPZ
               ENDDO                     !---- end loop over JPY
!*** diagnostic
!               write(0,*)'eself_v7 direct_calc ckpt 6'
!               write(0,*)'   icxzc,jcxzc,kcxzc=',icxzc,jcxzc,kcxzc
!***
               DO M=1,6
                  CXZC(ICXZC,JCXZC,KCXZC,M)=DCXSUM(M)
!*** diagnostic
!                  write(0,fmt='(a,i1,a,2f10.4)')'   dcxsum(',m,')=',dcxsum(m)
!***
               ENDDO
               IF(I==IMIN.AND.J==JMIN.AND.K==KMIN.AND.MYID==0)THEN

! NB: this call to CPU_TIME occurs within the first thread
!     (I==IMIN, J==JMIN, K==KMIN)

                  CALL CPU_TIME(T2)
                  T2=REAL((NX-IMIN+1)*(NY-JMIN+1)*(NZ-KMIN+1),KIND=WP)*(T2-T1)
                  IF(T2<600._WP)THEN
                     WRITE(CMSGNM,FMT='(A,F8.2,A)')                         &
                        'Estimated total cputime required by DIRECT_CALC=', &
                        T2,' cpu-sec'
                  ELSEIF(T2>=600._WP.AND.T2<3600._WP)THEN
                     WRITE(CMSGNM,FMT='(A,F8.2,A)')                         &
                        'Estimated total cputime required by DIRECT_CALC=', &
                        (T2/60._WP),' cpu-min'
                  ELSE
                     WRITE(CMSGNM,FMT='(A,F8.2,A)')                         &
                        'Estimated total cputime required by DIRECT_CALC=', &
                        (T2/3600._WP),' cpu-hr'
                  ENDIF
                  CALL WRIMSG('DIRECT_CALC',CMSGNM)
                  IF(PYDDX*PZDDX>0)THEN
                     WRITE(CMSGNM,FMT='(A)')          &
                          'cputime scales as 1/gamma^2'
                  ELSE
                     IF(PYDDX+PZDDX==0)THEN
                        WRITE(CMSGNM,FMT='(A)')               &
                              'cputime is independent of gamma'
                     ELSE
                        WRITE(CMSGNM,FMT='(A)')           &
                              '[cputime scales as 1/gamma]'
                     ENDIF
                  ENDIF
                  CALL WRIMSG('DIRECT_CALC',CMSGNM)
               ENDIF
            ENDDO                     !--- end loop over I
         ENDDO                     !-- end loop over J

! 2012.04.27 (BTD) it is not necessary to have a flush(cxzc) here
!                  flush is implied by both $omp END DO and 
!                  $omp END PARALLEL directives

      ENDDO                     !- end loop over K

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif

!*** diagnostic
!      write(0,*)'eself_v7 direct_calc ckpt 7'
!***
      CALL CPU_TIME(T2)
!*** diagnostic
!      write(0,*)'eself_v7 direct_calc ckpt 8'
!***
      T2=T2-T0
      IF(MYID==0)THEN
         IF(T2<600.D0)THEN
            WRITE(CMSGNM,FMT='(A,F8.2,A)')                   &
                  'Actual cputime to complete DIRECT_CALC=', &
                  T2,' cpu-sec'
         ELSEIF(T2>=600._WP.AND.T2<3600._WP)THEN
            WRITE(CMSGNM,FMT='(A,F8.2,A)')                   &
                  'Actual cputime to complete DIRECT_CALC=', &
                  (T2/60._WP),' cpu-min'
         ELSE
            WRITE(CMSGNM,FMT='(A,F8.2,A)')                   &
                  'Actual cputime to complete DIRECT_CALC=', &
                  (T2/3600._WP),' cpu-hr'
         ENDIF
         CALL WRIMSG('DIRECT_CALC',CMSGNM)
      ENDIF

!*** diagnostic
!      write(0,*)'eself_v7 direct_calc ckpt 9'
!***

! If IPBC=1: set the elements with ICXZC=NX+1 or JCXZC=NY+1
!            or KCXZC=NZ+1 to zero

      IF(IPBC==1)THEN

         DO M=1,6

#ifdef openmp
!$OMP PARALLEL                    &
!$OMP&   PRIVATE(ICXZC,JCXZC,KCXZC)
!$OMP DO
#endif

            DO KCXZC=1,2*NZ
               DO JCXZC=1,2*NY
                  CXZC(NX+1,JCXZC,KCXZC,M)=CXZERO
               ENDDO
            ENDDO

#ifdef openmp
!$OMP END DO
!$OMP DO
#endif

            DO KCXZC=1,2*NZ

               DO ICXZC=1,2*NX
                  CXZC(ICXZC,NY+1,KCXZC,M)=CXZERO
               ENDDO

            ENDDO

#ifdef openmp
!$OMP END DO
!$OMP DO
#endif

            DO JCXZC=1,2*NY

               DO ICXZC=1,2*NX
                  CXZC(ICXZC,JCXZC,NZ+1,M)=CXZERO
               ENDDO

            ENDDO

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif

         ENDDO
 
      ENDIF
      RETURN
 6700 FORMAT(F8.2,' cpu-sec')
 6701 FORMAT(F8.2,' cpu-min')
 6702 FORMAT(F8.2,' cpu-hr')
    END SUBROUTINE DIRECT_CALC
