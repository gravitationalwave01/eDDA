    SUBROUTINE ESELF(CMETHD,CXZP,NX,NY,NZ,IPBC,GAMMA,PYD,PZD,AK,AKD,DX, &
                     CXZC,CXZW,CXZE)
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE

!----------------------- eself v4.0 --------------------------------
! Arguments:

      CHARACTER(6) :: CMETHD
      INTEGER :: IPBC,NX,NY,NZ
      REAL(WP) :: AKD,GAMMA,PYD,PZD
      REAL(WP) :: AK(3),DX(3)
      COMPLEX(WP) :: &
         CXZC(NX+1+IPBC*(NX-1),NY+1+IPBC*(NY-1),NZ+1+IPBC*(NZ-1),6), &
         CXZE(NX,NY,NZ,3),    &
         CXZP(NX,NY,NZ,3),    &
         CXZW(2*NX,2*NY,2*NZ,*)

! Local scalars:

      CHARACTER :: CMSGNM*70
      INTEGER :: I,IR,ISGN,J,JR,JSGN,K,KR,KSGN,M,NGRID
      REAL(WP) :: AK2OLD,AK3OLD,AKD2,DTIME,FAC,PYDDX,PZDDX,WOLD
      COMPLEX(WP) :: CXEX,CXEY,CXEZ,CXXX,CXXY,CXXZ,CXYY,CXYZ,CXZZ

! Local arrays:

      INTEGER :: ISYM(3)

      EXTERNAL CXFFTW,DIRECT_CALC,EXTND,PAD,TIMEIT,TRIM
      INTRINSIC MIN,NINT,SIGN
      SAVE AK2OLD,AK3OLD,NGRID,WOLD
      DATA AK2OLD/-999._WP/
      DATA AK3OLD/-999._WP/
      DATA WOLD/-999._WP/

!-----------------------------------------------------------------------
! Parameter GAMMA determines the range of the sums when periodic
! boundary conditions are employed.  Dipole-dipole interactions
! are screened by a factor exp[-(gamma*kr)^4]
! The effective
! range/d = 1/(gamma*k*d) = 2000 if gamma=1e-3 and kd=0.5
!
! The sums are actually continued out to
! r/d = 2*/(gamma*kd) = 4000 if gamma=1e-3 , kd=0.5
! [screening factor = exp(-16)=1.1e-7]

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

!       PYD             (Period of lattice in y direction)/D(2)
!       PZD             (Period of lattice in z direction)/D(3)

!       DX(1-3)         Lattice spacing in x,y,z directions, in units of
!                       n**(-1./3.) .  Note that with this normalization
!                       we have DX(1)*DX(2)*DX(3)=1.

!       AK(1-3)         k(1-3)*d, where k = k vector in vacuo, and
!                       d = effective lattice spacing = (dx*dy*dz)**(1/3

!       AKD             Frequency of oscillation of dipoles and electric
!                       field; also absolute value of wave vector of
!                       incident wave [c=1]  (REAL*4).

!       CXZC            (NX+1)*(NY+1)*(NZ+1)*6 array of Green
!                       function coefficients used
!                       internally by ESELF and
!                       *************NOT TO BE OVERWRITTEN***********
!                       between calls, because these coefficients are
!                       recomputed only if W has changed since the last
!                       call to ESELF.

!       CXZW            Complex, scratch-space vector of length:
!                       2*NX*2*NY*2*NY*3
!                       See comment about FFT usage and CMETHD
!                       flag.
!                       Can be overwritten between calls to ESELF

! OUTPUT:

!       CXZE(I,J,K,L)   Lth component of dipole-generated electric field
!                       at grid point (I,J,K);
!                       the DECLARED length of ZE in the calling
!                       program is CXZE(NX,NY,NZ,3)
!                       [or CXZE(NX*NY*NZ,3) or CXZE(3*NX*NY*NZ)]

! Originally written by Jeremy Goodman, Princeton Univ. Observatory, 90.
! History:
! 90.11.29 (BTD): Modified to set untransformed ZC(1,1,1,1-6)=0.
! 90.11.29 (PJF): Modified to use FOURX and CXFFT99
! 90.12.05 (BTD): Modified for new ordering of elements of polarization
!                 and electric field vectors in calling program.  Modifi
!                 ESELF and PAD to remove distinction between NX,NY,NZ
!                 and dimensions CXZE and CXZP, since our new ordering
!                 always assumes this.
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
! end history
! Copyright (C) 1993,1994,1996,1997,1999,2000,2003,2004,2005,2006,2007,2008
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

!*** diagnostic
!      write(0,*)'eself ckpt 1 with nx,ny,nz=',nx,ny,nz
!      write(0,*)'                 pyd,pzd=',pyd,pzd
!      write(0,*)'                    ipbc=',ipbc
!      write(0,*)'                     akd=',akd
!      write(0,*)'                 dx(1-3)=',dx
!***

!*** diagnostic
!      write(0,*)'entered eself with akd=',akd
!      write(0,*)'                  wold=',wold
!      write(0,*)'entered eself with ak(2)=',ak(2)
!      write(0,*)'                  ak2old=',ak2old
!      write(0,*)'entered eself with ak(3)=',ak(3)
!      write(0,*)'                  ak3old=',ak3old
!      write(0,*)'                  pyd=',pyd
!      write(0,*)'                  pzd=',pzd
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
!      write(0,*)'recompute Green-function coefficients'
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
! range/d = 1/(gamma*kd) = 2000 if gamma=1e-3 and kd=0.5
! although the sums are actually continued out to
! r/d = 2/(gamma*kd) = 4000 if gamma=1e-3, kd=0.5
! [screening factor = exp(-16)=1.1e-7]

! PYD*d(2) = periodicity in Y direction
! PZD*d(3) = periodicity in Z direction

      IF(PYD>0._WP.OR.PZD>0._WP)THEN
         WRITE(CMSGNM,FMT='(A,2F8.2,A,1PE9.2)')'PBC with PYD, PZD=',PYD, &
                                               PZD,', GAMMA=',GAMMA
         CALL WRIMSG('ESELF ',CMSGNM)
      ENDIF

      PYDDX=PYD*DX(2)
      PZDDX=PZD*DX(3)

!               Compute the actual coefficients:

! Compute 6 independent elements of 3x3 symmetric matrix A_jk,where
! A_jk*P_k = -electric field at location j due to dipole P_k at location

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
!         write(0,*)'eself ckpt 2'
!***
        CALL DIRECT_CALC(1,1,1,NX,NY,NZ,IPBC,DX,AK,AKD,AKD2,GAMMA, &
                         PYDDX,PZDDX,CXZC(1,1,1,1))
!*** diagnostic
!         write(0,*)'eself ckpt 3'
!***

!*** diagnostic
!        write(0,*)'returned from direct_calc'
!        write(0,*)'check for NaN...'
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
!        write(0,*)'eself ckpt 4, about to call EXTEND'
!***
        CALL EXTND(CXZC(1,1,1,1),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
!*** diagnostic
!        write(0,*)'eself ckpt 5, returned from EXTEND'
!***
        IF(CMETHD=='GPFAFT')THEN
!*** diagnostic
!           write(0,*)'eself ckpt 6'
!***
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
!*** diagnostic
!           write(0,*)'eself ckpt 7'
!**
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF
!*** diagnostic
!        write(0,*)'eself ckpt 7.1, about to call TRIM'
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
!        write(0,*)'eself ckpt 7.2, about to call EXTND'
!***
        CALL EXTND(CXZC(1,1,1,2),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
        IF(CMETHD=='GPFAFT')THEN
!*** diagnostic
!           write(0,*)'eself ckpt 7.3, about to call cxfft3n'
!***
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
!*** diagnostic
!          write(0,*)'eself ckpt 7.4, returned from cxfft3n'
!***
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF
!*** diagnostic
!        write(0,*)'eself ckpt 7.5, about to call TRIM'
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

        CALL DIRECT_CALC(-1,-1,-1,NX,NY,NZ,IPBC,DX,AK,AKD,AKD2,GAMMA, &
                         PYDDX,PZDDX,CXZC(1,1,1,1))

!*** diagnostic
!        write(0,*)'eself ckpt 8'
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
!        write(0,*)'eself ckpt 9'
!***

! CXZC now contains the full Fourier transform of the A convolution
! and should not be overwritten between calls to ESELF

      ENDIF
      CALL TIMEIT('ESELF (first call)',DTIME)
      CALL TIMEIT('ESELF',DTIME)
!-----------------------------------------------------------------------

! End of recomputation of Green-function coefficients

70    CONTINUE
!*** diagnostic
!      write(0,*)'eself ckpt 10'
!****
! Fourier transform the polarizations:

      DO M=1,3
         CALL PAD(CXZP(1,1,1,M),NX,NY,NZ,CXZW(1,1,1,M))

!*** diagnostic
!        write(0,*)'eself ckpt 11: returned from PAD: ', &
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
!        write(0,*)'eself ckpt 12: cxzw checked for NaN or overflow: ', &
!                  jr,' instances found'
!        if(jr>0)stop
!***

        IF(CMETHD=='GPFAFT')THEN
           CALL CXFFT3N(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,+1)
!*** diagnostic
!        write(0,*)'eself ckpt 13: returned from CXFFT3N: ', &
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
!        write(0,*)'eself ckpt 14: cxzw checked for NaN or overflow',jr, &
!                  ' instances found'
!        if(jr>0)stop
!***
         ELSEIF(CMETHD=='FFTW21')THEN
            CALL CXFFTW(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,+1)
         ELSEIF(CMETHD.EQ.'FFTMKL')THEN
            CALL CXFFT3_MKL(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,+1)
         ENDIF
      ENDDO

!***********************************************************************

! Multiply by F.t. of Green function.

      IF(IPBC==0)THEN

!*** diagnostic
!         write(0,*)'eself ckpt 15, check cxzc(i,j,k,m=1,6)'
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
!        write(0,*)'eself ckpt 16: cxzc(i,j,k,m=1-6) checked for NaN ', &
!                  'or overflow: ',jr,' instances found'
!        if(jr>0)stop
!        write(0,*)'eself ckpt 17: check cxzw(i,j,k,m=1-3)'
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
!        write(0,*)'eself ckpt 18: cxzw(i,j,k,m=1-3) checked for NaN or ', &
!                  'overflow: ',jr,' instances found'
!        if(jr>0)stop
!***
!***

! If IPBC=0, then only one octant of F.t. of Green function has been
!            stored, but can recover others using symmetry.

#ifdef openmp
!$omp parallel do                                          &
!$omp& private(K,J,I,KSGN,KR,JSGN,JR,ISGN,IR)              &
!$omp& private(CXXX,CXXY,CXXZ,CXYY,CXYZ,CXZZ,CXEX,CXEY,CXEZ)
#endif

        DO K=1,2*NZ
           KSGN=NINT(SIGN(1._WP,NZ+1.5_WP-K))
           KR=MIN(K,2*NZ+2-K)
!*** diagnostic
!          write(0,*)'K,KSGN,KR=',K,KSGN,KR
!***
           DO J=1,2*NY
              JSGN=NINT(SIGN(1._WP,NY+1.5_WP-J))
              JR=MIN(J,2*NY+2-J)
!*** diagnostic
!            write(0,*)'   J,JSGN,JR=',J,JSGN,JR
!***
              DO I=1,2*NX
                 ISGN=NINT(SIGN(1._WP,NX+1.5_WP-I))
                 IR=MIN(I,2*NX+2-I)
!*** diagnostic
!              write(0,*)'       I,ISGN,IR=',I,ISGN,IR
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

#ifdef openmp
!$omp end parallel do
#endif

      ELSEIF(IPBC==1)THEN

! If IPBC=1, then the full F.t. of the Green function has been stored.

#ifdef openmp
!$omp parallel do                                                &
!$omp& private(K,J,I,CXXX,CXXY,CXXZ,CXYY,CXYZ,CXZZ,CXEX,CXEY,CXEZ)
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
!$omp end parallel do
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
!$omp parallel do private(k,j,i)
#endif

            DO K=1,NZ
               DO J=1,NY
                  DO I=1,NX
                     CXZE(I,J,K,M)=FAC*CXZW(I,J,K,M)
                  ENDDO
               ENDDO
            ENDDO

#ifdef openmp
!$omp end parallel do
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
!$omp parallel
!$omp do private(k,j,i)
#endif
      DO K=1,NZ
         DO J=1,NY
            DO I=1,NX
               CXB(I,J,K)=CXA(I,J,K)
            ENDDO
         ENDDO
      ENDDO
#ifdef openmp
!$omp enddo
#endif

!-----------------------------------------------------------------------
! x -> -x

#ifdef openmp
!btd 080627 moved I=NX+1 out of loop
!$omp single
#endif
      CXB(NX+1,1:NY,1:NZ)=CXZERO
#ifdef openmp
!$omp end single
!$omp do private(k,j,i)
#endif
      DO K=1,NZ
         DO J=1,NY
            DO I=NX+2,2*NX
               CXB(I,J,K)=CXA(2*NX+2-I,J,K)*ISYM(1)
            ENDDO
         ENDDO
      ENDDO
#ifdef openmp
!$omp end do
#endif

!-----------------------------------------------------------------------
! y -> -y

!btd 080627 moved J=NY+1 out of loop, switched order of loops I and J
#ifdef openmp
!$omp single
#endif
      CXB(1:2*NX,NY+1,1:NZ)=CXZERO
#ifdef openmp
!$omp end single
!$omp do private(k,j,i) private(j2)
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
!$omp end do
#endif

!-----------------------------------------------------------------------
! z -> -z

#ifdef openmp
!Art pulling the 3rd dimension to the outer most loop.
!Art we'll do this expression in only the thread that has NZ+1
!btd 080627 reorder loops: J,I,K -> K,J,I
!
!$omp single
#endif
      CXB(1:2*NX,1:2*NY,NZ+1)=CXZERO
#ifdef openmp
!$omp end single
!$omp do private(k,j,i) private(k2)
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
!$omp end do
!$omp end parallel
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
!$omp parallel do          &
!$omp&    private(k,j,i)   &
!$omp&    shared(NX,NY,NZ) &
!$omp&    shared(CXA,CXB)
#endif
      DO K=1,NZ+1
         DO J=1,NY+1
            DO I=1,NX+1
               CXA(I,J,K)=CXB(I,J,K)
            ENDDO
         ENDDO
      ENDDO
#ifdef openmp
!$omp end parallel do
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
!$omp parallel
!$omp do private(k,j,i)
#endif
      DO K=1,2*NZ
         DO J=1,2*NY
            DO I=1,2*NX
               CXB(I,J,K)=CXZERO
            ENDDO
         ENDDO
      ENDDO
#ifdef openmp
!$omp end do
!$omp do private(k,j,i)
#endif
      DO K=1,NZ
         DO J=1,NY
            DO I=1,NX
               CXB(I,J,K)=CXA(I,J,K)
            ENDDO
         ENDDO
      ENDDO
#ifdef openmp
!$omp end do
!$omp end parallel
#endif

      RETURN
    END SUBROUTINE PAD

!-----------------------------------------------------------------------

    SUBROUTINE DIRECT_CALC(IX,IY,IZ,NX,NY,NZ,IPBC,DX,AK,AKD,AKD2,GAMMA, &
                           PYDDX,PZDDX,CXZC)
      USE DDPRECISION,ONLY: WP
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

      CHARACTER :: CMSGNM*70
      INTEGER :: I,ICXZC,II,IMIN,J,JCXZC,JJ,JMIN,JPY,JPYM,JPZ,JPZM, &
                 K,KCXZC,KMIN,M
      REAL(WP) :: GAMMAKD4,PHASY,PHASYZ,R,R2,R3,RANGE,RANGE2,RJPY, &
                  RJPZ,T0,T1,T2,X0,X2,X2Y2,XX,Y0,Y2,Z0
      REAL(WP) :: X(3)
      COMPLEX(WP) :: CXFAC,CXI,CXPHAS,CXTERM,CXZERO
      COMPLEX(WP) :: DCXSUM(6)
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

!        IPBC     = 0 if only doing first octant (IX=IY=IZ=1)
!                 = 1 otherwise
!                   N.B.: IPBC affects dimensioning of CXZC

!        NX,NY,NZ = size of first octant (I = 1 -> NX,
!                                         J = 1 -> NY,
!                                         K = 1 -> NZ)
!        PYDDX     = period of lattice in y direction/d
!        PZDDX     = period of lattice in z direction/d
!        DX(1-3)   = lattice spacing in x,y,z direction, in units of
!                    d = n**(-1/3).
!        AK(1-3)   = k(1-3)*d, where k = k vector in vacuo
!        AKD       = |k|d
!        CXZC      = array with dimensions
!                    (NX+1)*(NY+1)*(NZ+1)*6 if IPBC=0
!                    (2*NX)*(2*NY)*(2*NZ)*6 if IPBC=1

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
!                
! end history
!-----------------------------------------------------------------------
!*** diagnostic
!      write(0,*)'direct_calc ckpt 1'
!***
      GAMMAKD4=(GAMMA*AKD)**4
      RANGE=2._WP/(GAMMA*AKD)
      RANGE2=RANGE*RANGE
      IF(PYDDX<=0._WP)THEN
         JPYM=0
      ELSE
         JPYM=1+NINT(RANGE/PYDDX)
      ENDIF

!*** diagnostic
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
#ifdef openmp
      if(myid==0.and.omp_get_thread_num()==0)then
!Art omp_get_num_threads is an integer function defined in omp_lib
         write(0,*)'direct_calc ckpt 1: omp_get_num_threads = ', &
                   omp_get_num_threads()
!         write(0,*)'DIRECT_CALC: kmin, nz = ',kmin,nz,nz-kmin+1
!         write(0,*)'DIRECT_CALC: jmin, ny = ',jmin,ny,ny-jmin+1
!         write(0,*)'DIRECT_CALC: imin, nx = ',imin,nx,nx-imin+1
      endif
#endif

! Determine elapsed cpu time for these sums

      CALL CPU_TIME(T0)

#ifdef openmp
!$omp parallel
!$omp do                                                 &
!$omp&   private(K,J,I,M,JPY,JPZ,II)                     &
!$omp&   private(JCXZC,JPZM,KCXZC)                       &
!$omp&   private(PHASY,PHASYZ,R2,R3,RJPY,RJPZ)           &
!$omp&   private(X,X0,X2,X2Y2,XX,Y0,Y2,Z0)               &
!$omp&   private(CXFAC,CXTERM,CXPHAS,DCXSUM)
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

! include artificial factor exp[-(gamma*kr)^4] to assist convergence

                        CXTERM=CXI*AKD*R
                        CXPHAS=EXP(CXTERM+CXI*PHASYZ-GAMMAKD4*R2*R2)/R3
                        CXFAC=(1._WP-CXTERM)/R2

! II=1      -> M=1   a_1 (xx)
!      JJ=2      2   a_2 (xy)
!         3      3   a_3 (xz)
!    2           4   a_4 (yy)
!         3      5   a_5 (yz)
!    3           6   a_6 (zz)

                        M=0
                        DO II=1,3              !------ loop over II
                           M=M+1
                           DCXSUM(M)=DCXSUM(M)-CXPHAS*(AKD2*(X(II)**2-R2)+ &
                                     CXFAC*(R2-3._WP*X(II)**2))
                           IF(II<3)THEN
                              DO JJ=II+1,3        !------- loop over JJ
                                 M=M+1
                                 XX=X(II)*X(JJ)
                                 DCXSUM(M)=DCXSUM(M)-                      &
                                           CXPHAS*(AKD2*XX-CXFAC*(3._WP*XX))
                              ENDDO               !------- end loop over JJ
                           ENDIF
                        ENDDO                  !------ end loop over II
                     ENDIF
                  ENDDO                     !----- end loop over JPZ
               ENDDO                     !---- end loop over JPY
               DO M=1,6
                  CXZC(ICXZC,JCXZC,KCXZC,M)=DCXSUM(M)
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
      ENDDO                     !- end loop over K
#ifdef openmp
!$omp end do
!$omp end parallel
#endif
!*** diagnostic
!      write(0,*)'direct_calc ckpt 2'
!***
      CALL CPU_TIME(T2)
!*** diagnostic
!      write(0,*)'direct_calc ckpt 3'
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
!      write(0,*)'direct_calc ckpt 4'
!***

! If IPBC=1: set the elements with ICXZC=NX+1 or JCXZC=NY+1
!            or KCXZC=NZ+1 to zero

      IF(IPBC==1)THEN
#ifdef openmp
!$omp parallel
#endif 
         DO M=1,6
#ifdef openmp
!$omp do                    &
!$omp&   private(KCXZC,JCXZC)
#endif
            DO KCXZC=1,2*NZ
               DO JCXZC=1,2*NY
                  CXZC(NX+1,JCXZC,KCXZC,M)=CXZERO
               ENDDO
            ENDDO
#ifdef openmp
!$omp end do
!$omp do                    &
!$omp&   private(KCXZC,ICXZC)
#endif
            DO KCXZC=1,2*NZ
               DO ICXZC=1,2*NX
                  CXZC(ICXZC,NY+1,KCXZC,M)=CXZERO
               ENDDO
            ENDDO
#ifdef openmp
!$omp end do
!$omp do                    &
!$omp&   private(JCXZC,ICXZC)
#endif
            DO JCXZC=1,2*NY
               DO ICXZC=1,2*NX
                  CXZC(ICXZC,JCXZC,NZ+1,M)=CXZERO
               ENDDO
            ENDDO
#ifdef openmp
!$omp end do
#endif
         ENDDO
#ifdef openmp
!$omp end parallel
#endif
 
      ENDIF
      RETURN
 6700 FORMAT(F8.2,' cpu-sec')
 6701 FORMAT(F8.2,' cpu-min')
 6702 FORMAT(F8.2,' cpu-hr')
    END SUBROUTINE DIRECT_CALC
