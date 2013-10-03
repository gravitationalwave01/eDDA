    SUBROUTINE BSELF(CMETHD,CXZP,NX,NY,NZ,IPBC,GAMMA,PYD,PZD,AK,AKD,DX, &
                     CXZG,CXZW,CXZB)
      USE DDPRECISION,ONLY: WP
      USE DDCOMMON_0,ONLY: AK2OLD_B,AK3OLD_B,NGRID,WOLD_B
      IMPLICIT NONE

!----------------------- bself v4 --------------------------------
! Arguments:

      CHARACTER(6) :: CMETHD
      INTEGER :: IPBC,NX,NY,NZ
      REAL(WP) :: AKD,GAMMA,PYD,PZD
      REAL(WP) :: AK(3),DX(3)
      COMPLEX(WP) ::                                                 &
         CXZB(NX,NY,NZ,3),                                           &
         CXZG(NX+1+IPBC*(NX-1),NY+1+IPBC*(NY-1),NZ+1+IPBC*(NZ-1),3), &
         CXZP(NX,NY,NZ,3),                                           &
         CXZW(2*NX,2*NY,2*NZ,*)

! NB: module DDCOMMON_0 must have previously set values of
!       AK2OLD_B,AK3OLD_B,WOLD_B
!    to be used by BSELF

! Local scalars:

      CHARACTER :: CMSGNM*70
      INTEGER :: I,IR,ISGN,J,JR,JX,JY,JZ,JSGN,K,KR,KSGN,M
      REAL(WP) :: AKD2,DTIME,FAC,PYDDX,PZDDX
      COMPLEX(WP) :: CXBX,CXBY,CXBZ,CXXY,CXXZ,CXYZ

! Local arrays:

      INTEGER :: ISYM(3)

#ifdef openmp
      INTEGER NTHREADS,TID
#endif

      EXTERNAL CXFFTW,DIRECT_CALCB,EXTND,PAD,TIMEIT,TRIM
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

!-----------------------------------------------------------------------
! subroutine BSELF

!       Given the dipole moments, CXZP, at all points on
!       a rectangular grid, oscillating at frequency AKD,
!       compute the magnetic field amplitude, CXZB,
!       at each point produced by all the other dipoles except the one
!       at that point.
!       The relationship between the dipoles and the field
!       values at the grid points can be expressed as a convolution,
!       and the convolution is efficiently evaluated using 3D FFTs.
!
!       options for computation of 3-dimensional FFT:
!
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
!
!   if CMETHD='FFTW21'
!      Use CXFFTW interface to FFTW (Fastest Fourier Transform in the
!               West) version 2.1.x from Frigo and Johnson.
!
!   if CMETHD='FFTMKL':
!      Use CXFFT3_MKL interface to Intel Math Kernel Library (MKL) FFT
!
! Input:
!
!    CXZP(I,J,K,L)   Lth cartesian component of the dipole
!                    moment at the grid point (I,J,K);
!                    the DIMENSIONed length of CXZP in
!                    the calling routine is CXZP(NX,NY,NZ,3)
!                    [or CXZP(NX*NY*NZ,3) or CXZP(3*NX*NY*NZ)]

!    NX,NY,NZ        Size of grid in x,y,z directions (INTEGER).

!    IPBC          = 0 for isolated target
!                    1 for periodic target

!    GAMMA         = coefficient used to assist convergence of sums
!                    over replica dipoles by suppressing long-range
!                    contributions with factor exp(-gamma*(kr)^2)
!                    typical value gamma = 0.005
!                    The effective
!                    range/d = 1/(gamma*k*d) = 400 if gamma=5e-3 and kd=0.5
!                    range/lambda = 1/(2*pi*gamma) = 31.8
!                    The sums are actually continued out to
!                    r/d = 2*/(gamma*kd) = 800 if gamma=5e-3 , kd=0.5
!                    [screening factor = exp(-16)=1.1e-7]

!       PYD          (Period of lattice in y direction)/DX(2)
!       PZD          (Period of lattice in z direction)/DX(3)

!       DX(1-3)      Lattice spacing in x,y,z directions, in units of
!                    n**(-1./3.) .  Note that with this normalization
!                       we have DX(1)*DX(2)*DX(3)=1.

!       AK(1-3)         k(1-3)*d, where k = k vector in vacuo, and
!                       d = effective lattice spacing = (dx*dy*dz)**(1/3)

!       AKD             = (omega/c)*d = k*d (dimensionless)

!       CXZG            (NX+1)*(NY+1)*(NZ+1)*3 array of Green
!                       function coefficients used
!                       internally by BSELF and
!                       *************NOT TO BE OVERWRITTEN***********
!                       between calls, because these coefficients are
!                       recomputed only if W has changed since the last
!                       call to BSELF.

!       CXZW            Complex, scratch-space vector of length:
!                       2*NX*2*NY*2*NY*3
!                       See comment about FFT usage and CMETHD flag.
!                       Can be overwritten between calls to BSELF

! OUTPUT:

!       CXZB(I,J,K,L)   Lth component of dipole-generated magnetic field
!                       at grid point (I,J,K);
!                       the DECLARED length of ZB in the calling
!                       program is CXZB(NX,NY,NZ,3)
!                       [or CXZB(NX*NY*NZ,3) or CXZB(3*NX*NY*NZ)]
!
!=======================================================================
! subroutine BSELF
! history
! based on subroutine ESELF written originally by Jeremy Goodman
! BSELF created by Ian Wong, Princeton University, July 2012
! 12.07.07 (IW)  v1 written
! 12.07.11 (BTD) v2 created from v1
!                * corrected typos (HXZB -> CXZB)
!                * changed notation
!                  HXZC -> CXZG  (Green function)
!                  HXZW -> CXZW
!                  HXZB -> CXZB
!                * a few changes to comments
! 12.12.21 (BTD) v3
!                * added comments
!                * corrected error for periodic targets
! 12.12.22 (BTD) * revised comments, minor cleanup
! 13.01.03 (BTD) v4
!                * added AK2OLD_B,AK3OLD_B,WOLD_B from DDCOMMON_0
!                * modified to skip recomputation of Green-function
!                  coefficients on second call to BSELF
!                  NB: this requires that CXZG *not* be deallocated
!                  in subroutine NEARFIELD after first call to BSELF
! end history
! Copyright (C) 2013 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!=======================================================================        

! check if we can skip recomputation of Green-function coefficients

      IF(PYD.EQ.0._WP.AND.PZD.EQ.0._WP)THEN
         IF(ABS(WOLD_B-AKD)<1.E-6_WP*AKD)GOTO 70
      ELSEIF(PYD.NE.0._WP.AND.PZD.EQ.0._WP)THEN
         IF(ABS(WOLD_B-AKD)<1.E-6_WP*AKD.AND.   &
            ABS(AK(2)-AK2OLD_B)<1.E-6_WP*AKD)GOTO 70
      ELSEIF(PYD.EQ.0._WP.AND.PZD.NE.0._WP)THEN
         IF(ABS(WOLD_B-AKD)<1.E-6_WP*AKD.AND.   &
            ABS(AK(3)-AK3OLD_B)<1.E-6_WP*AKD)GOTO 70
      ELSEIF(PYD.NE.0._WP.AND.PZD.NE.0._WP)THEN
         IF(ABS(WOLD_B-AKD)<1.E-6_WP*AKD.AND.   &
            ABS(AK(2)-AK2OLD_B)<1.E-6_WP*AKD.AND. &
            ABS(AK(3)-AK3OLD_B)<1.E-6_WP*AKD)GOTO 70
      ENDIF

! Compute Green function coefficients 

! We have to compute the Green-function coefficients giving
! components of magnetic field strength at a each grid point R
! produced by unit-valued component of dipole moment at
! point R', and then Fourier transform these components.

      WOLD_B=AKD
      AK2OLD_B=AK(2)
      AK3OLD_B=AK(3)
      NGRID=8*NX*NY*NZ
      AKD2=AKD*AKD

! We assume screening function exp(-(gamma*kr)^4) so
! range/d = 1/(gamma*kd) = 2000 if gamma=1e-3 and kd=0.5
! although the sums are actually continued out to
! r/d = 2/(gamma*kd) = 4000 if gamma=1e-3, kd=0.5
! [screening factor = exp(-16)=1.1e-7]

! PYD*DX(2) = periodicity in Y direction
! PZD*DX(3) = periodicity in Z direction

      IF(PYD>0._WP.OR.PZD>0._WP)THEN
         WRITE(CMSGNM,FMT='(A,2F8.2,A,1PE9.2)')'PBC with PYD, PZD=',PYD, &
                                               PZD,', GAMMA=',GAMMA
         CALL WRIMSG('BSELF',CMSGNM)
      ENDIF

      PYDDX=PYD*DX(2)
      PZDDX=PZD*DX(3)


! Compute 3 independent elements of 3x3 anti-symmetric matrix C_jk,where
! C_jk*P_k = magnetic field at location j due to dipole P_k at location

! C_jk = ( 0    c_1  c_2)
!        (-c_1  0    c_3)
!        (-c_2 -c_3   0 )_jk

      IF(IPBC==0)THEN

! initialize CXZG(I,J,K,M) = c_M for magnetic field at (I,J,K)
!                            produced by a dipole at (1,1,1)
!                            and replica dipoles (if PYD or PYZ are
!                            nonzero).

! need to calculate this for all (I,J,K) for one octant:
! I running from 1 to NX, J from 1 to NY, K from 1 to NZ

! Later obtain C_jk values for other octants by using symmetry

!*** diagnostic
!         write(0,*)'bself_v4 ckpt 1: call DIRECT_CALCB'
!***
        CALL DIRECT_CALCB(1,1,1,NX,NY,NZ,IPBC,DX,AK,AKD,AKD2,GAMMA, &
                         PYDDX,PZDDX,CXZG(1,1,1,1))
!*** diagnostic
!        write(0,*)'bself_v4 ckpt 2'
!        write(0,*)'   returned from direct_calcb'
!        write(0,*)'   check for NaN...'
!        jr=0
!        do i=1,nx
!           do j=1,ny
!              do k=1,nz
!                 do m=1,3
!                    if(.not.(abs(hxzc(i,j,k,m))>=0.d0).or. &
!                        abs(hxzc(i,j,k,m))>=1.d100)then
!                       write(0,*)'i,j,k,m=',i,j,k,m,'hxzc=',hxzc(i,j,k,m)
!                       jr=jr+1
!                    endif
!                 enddo
!              enddo
!           enddo
!        enddo
!        write(0,*)'bself_v4 ckpt 3: cxy checked for NaN or overflow',jr, &
!                   ' instances found'
!***

! At this point, CXZG(I,J,K,1-3) contains the upper triangular part of the
! anti-symmetric 3 x 3 matrix giving the magnetic field at grid point (i,j,k)
! produced by a dipole at (1,1,1)

! Fill out CXZG to twice the size in each grid dimension to accomodate
! negative lags [periodicity in each dimension is assumed, so (e.g.)
! nx < i <= 2*nx is equivalent to -nx < i <= 0], exploiting symmetries,
! and then Fourier transform.

! If PYDDX=0 and PZDDX=0 , need only do direct calculation of C matrix
! for first octant, since remaining octants can be obtained by symmetry.
! After calculating C matrix, store only the first octant of the
! transform, since the others can be obtained by symmetry.

!-----------------------------------------------------------------------
! extend c_1 = c_xy(x,y,z) : c -> +c for x -> -x
!                                 +c     y -> -y
!                                 -c     z -> -z

        ISYM(1)=1
        ISYM(2)=1
        ISYM(3)=-1
!*** diagnostic
!        write(0,*)'bself_v4 ckpt 4, about to call EXTND'
!***
        CALL EXTND(CXZG(1,1,1,1),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
        IF(CMETHD=='GPFAFT')THEN
!*** diagnostic
!           write(0,*)'bself_v4 ckpt 5, about to call cxfft3n'
!***
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
!*** diagnostic
!          write(0,*)'bself_v4 ckpt 6, returned from cxfft3n'
!***
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF
!*** diagnostic
!        write(0,*)'bself_v4 ckpt 7, about to call TRIM'
!***
        CALL TRIM(CXZW(1,1,1,1),NX,NY,NZ,CXZG(1,1,1,1))
!*** diagnostic
!        write(0,*)'returned from TRIM'
!***
!-----------------------------------------------------------------------
! extend c_2 = c_xz(x,y,z) : c -> +c for x -> -x
!                                 -c     y -> -y
!                                 +c     z -> -z

        ISYM(1)=1
        ISYM(2)=-1
        ISYM(3)=1
        CALL EXTND(CXZG(1,1,1,2),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
        IF(CMETHD=='GPFAFT')THEN
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF

        CALL TRIM(CXZW(1,1,1,1),NX,NY,NZ,CXZG(1,1,1,2))

!-----------------------------------------------------------------------
! extend c_3 = c_yz(x,y,z) : c -> -c for x -> -x
!                                 +c     y -> -y
!                                 +c     z -> -z
        ISYM(1)=-1
        ISYM(2)=1
        ISYM(3)=1
        CALL EXTND(CXZG(1,1,1,3),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
        IF(CMETHD=='GPFAFT')THEN
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF

        CALL TRIM(CXZW(1,1,1,1),NX,NY,NZ,CXZG(1,1,1,3))


      ELSEIF(IPBC==1)THEN

! This point is reached when PYDDX or PZDDX are nonzero.
! When PBC are used for general direction of incident wave,
! all octants of C matrix require direct calculation: symmetries valid
! for single target no longer apply because of position-dependent phases
! of replica dipoles.

! DIRECT_CALCB computes arrays (c_1)_jk , (c_2)_jk  (c_3)_jk
!
!               (   0   c_1  c_2 )
! where B(r_j)= ( -c_1   0   c_3 ) * P(r_k)
!               ( -c_2 -c_3   0  )
!
! CXZG(I,J,K,M) = c_M for (r_j - r_k)/d = (I-1)*xhat + (J-1)*yhat + (K-1)*zhat
!
! when IPBC=1, CXZG includes contribution to B from replica dipoles

!*** diagnostic
!         write(0,fmt='(a,a,4I4)')'bself_v4 ckpt 7.5',    &
!              ' call direct_calcb with NX,NY,NZ=',NX,NY,NZ
!***

        CALL DIRECT_CALCB(-1,-1,-1,NX,NY,NZ,IPBC,DX,AK,AKD,AKD2,GAMMA, &
                         PYDDX,PZDDX,CXZG(1,1,1,1))

!*** diagnostic
!        write(0,*)'bself_v4 ckpt 8'
!***
! The array CXZG(1-2*NX,1-2*NY,1-2*NZ,1-3) of C matrix coefficients
! now covers all octants.

! Fourier transform the C matrix CXZG:

        DO M=1,3
           IF(CMETHD=='GPFAFT')THEN
              CALL CXFFT3N(CXZG(1,1,1,M),2*NX,2*NY,2*NZ,+1)
           ELSEIF(CMETHD=='FFTW21')THEN
              CALL CXFFTW(CXZG(1,1,1,M),2*NX,2*NY,2*NZ,+1)
           ELSEIF(CMETHD.EQ.'FFTMKL')THEN
              CALL CXFFT3_MKL(CXZG(1,1,1,M),2*NX,2*NY,2*NZ,+1)
           ENDIF
        ENDDO

!*** diagnostic
!        write(0,*)'bself_v4 ckpt 9'
!***

! CXZG now contains the full Fourier transform of the C convolution
! and should not be overwritten between calls to BSELF

      ENDIF
!      CALL TIMEIT('BSELF (first call)',DTIME)
!      CALL TIMEIT('BSELF',DTIME)
!-----------------------------------------------------------------------

! End of computation of Green-function coefficients

70    CONTINUE

!*** diagnostic
!      write(0,*)'bself_v4 ckpt 10'
!****
! Fourier transform the polarizations

    DO M=1,3
       CALL PAD(CXZP(1,1,1,M),NX,NY,NZ,CXZW(1,1,1,M))
!*** diagnostic
!        write(0,*)'bself_v4 ckpt 11: returned from PAD: ', &
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
!        write(0,*)'bself_v4 ckpt 12: cxzw checked for NaN or overflow: ', &
!                  jr,' instances found'
!        if(jr>0)stop
!***

        IF(CMETHD=='GPFAFT')THEN
           CALL CXFFT3N(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,+1)
!*** diagnostic
!        write(0,*)'bself_v4 ckpt 13: returned from CXFFT3N: ', &
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
!        write(0,*)'bself_v4 ckpt 14: cxzw checked for NaN or overflow',jr, &
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
!         write(0,*)'bself_v4 ckpt 15, check hxzc(i,j,k,m=1,6)'
!        jr=0
!        do i=1,2*nx
!           do j=1,2*ny
!              do k=1,2*nz
!                 do m=1,3
!                    if(.not.(abs(hxzc(i,j,k,m))>=0.d0).or. &
!                        abs(hxzc(i,j,k,m))>=1.d100)then
!                       write(0,*)'i,j,k,m=',i,j,k,m,'hxzc=',hxzc(i,j,k,m)
!                       jr=jr+1
!                    endif
!                 enddo
!              enddo
!           enddo
!        enddo
!        write(0,*)'bself_v4 ckpt 16: hxzc(i,j,k,m=1-3) checked for NaN ', &
!                  'or overflow: ',jr,' instances found'
!        if(jr>0)stop
!        write(0,*)'bself_v4 ckpt 17: check cxzw(i,j,k,m=1-3)'
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
!        write(0,*)'bself_v4 ckpt 18: cxzw(i,j,k,m=1-3) checked for NaN or ', &
!                  'overflow: ',jr,' instances found'
!        if(jr>0)stop
!***
!***

! If IPBC=0, then only one octant of F.t. of Green function has been
!            stored, but can recover others using symmetry.

#ifdef openmp
!$OMP PARALLEL DO                                            &
!$OMP&   PRIVATE(K,J,I,KSGN,KR,JSGN,JR,ISGN,IR)              &
!$OMP&   PRIVATE(CXXY,CXXZ,CXYZ,CXBX,CXBY,CXBZ)
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
                 CXXY=CXZG(IR,JR,KR,1)*(ISGN*JSGN)
                 CXXZ=CXZG(IR,JR,KR,2)*(ISGN*KSGN)
                 CXYZ=CXZG(IR,JR,KR,3)*(JSGN*KSGN)
                 
!*** diagnostic
!              if(.not.(abs(cxxy)>=0.d0))write(0,*) &
!                'ir,jr,kr,hxzc(ir,jr,kr,1)=',ir,jr,kr,hxzc(ir,jr,kr,1)
!              if(.not.(abs(cxxz)>=0.d0))write(0,*) &
!                 'ir,jr,kr,hxzc(ir,jr,kr,2)=',ir,jr,kr,hxzc(ir,jr,kr,2)
!              if(.not.(abs(cxyz)>=0.d0))write(0,*) &
!                  'ir,jr,kr,hxzc(ir,jr,kr,3)=',ir,jr,kr,hxzc(ir,jr,kr,3)

!***
                 CXBX=CXXY*CXZW(I,J,K,2)+CXXZ*CXZW(I,J,K,3)
                 CXBY=-CXXY*CXZW(I,J,K,1)+CXYZ*CXZW(I,J,K,3)
                 CXBZ=-CXXZ*CXZW(I,J,K,1)-CXYZ*CXZW(I,J,K,2)
                      
!*** diagnostic
!              if(.not.(abs(cxbx+cxby+cxbz)>=0.d0))then
!                 write(0,*)'i,j,j,ir,jr,kr=',i,j,k,ir,jr,kr
!                 write(0,*)'hxzc(ir,jr,kr,1)=',hxzc(ir,jr,kr,1)
!                 write(0,*)'hxzc(ir,jr,kr,2)=',hxzc(ir,jr,kr,2)
!                 write(0,*)'hxzc(ir,jr,kr,3)=',hxzc(ir,jr,kr,3)
!                 write(0,*)'      cxbx=',cxbx
!                 write(0,*)'      cxby=',cxby
!                 write(0,*)'      cxbz=',cxbz
!                 write(0,*)'      cxxy=',cxxy
!                 write(0,*)'      cxxz=',cxxz
!                 write(0,*)'      cxyz=',cxyz
!                 write(0,*)'      cxzw(i,j,k,1)=',cxzw(i,j,k,1)
!                 write(0,*)'      cxzw(i,j,k,2)=',cxzw(i,j,k,2)
!                 write(0,*)'      cxzw(i,j,k,3)=',cxzw(i,j,k,3)
!                 stop
!              endif
!***

! now overwrite CXZW [which contained the Fourier transform of C matrix]
!                     with F.t. of B

                 CXZW(I,J,K,1)=CXBX
                 CXZW(I,J,K,2)=CXBY
                 CXZW(I,J,K,3)=CXBZ

! CXZW is now the Fourier transform of B

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
!$OMP END PARALLEL DO
#endif

      ELSEIF(IPBC==1)THEN

! If IPBC=1, then the full F.t. of the Green function has been stored.

#ifdef openmp
!$OMP PARALLEL DO                                                  &
!$OMP&   PRIVATE(K,J,I,CXXY,CXXZ,CXYZ,CXBX,CXBY,CXBZ)
#endif

         DO K=1,2*NZ
            DO J=1,2*NY
               DO I=1,2*NX
                  CXXY=CXZG(I,J,K,1)
                  CXXZ=CXZG(I,J,K,2)
                  CXYZ=CXZG(I,J,K,3)
                  CXBX=CXXY*CXZW(I,J,K,2)+CXXZ*CXZW(I,J,K,3)
                  CXBY=-CXXY*CXZW(I,J,K,1)+CXYZ*CXZW(I,J,K,3)
                  CXBZ=-CXXZ*CXZW(I,J,K,1)-CXYZ*CXZW(I,J,K,2)

! now overwrite CXZW [which contained the F.t. of C matrix]
!                     with Fourier transform of B

                  CXZW(I,J,K,1)=CXBX
                  CXZW(I,J,K,2)=CXBY
                  CXZW(I,J,K,3)=CXBZ

! CXZW is now the Fourier transform of B

               ENDDO
            ENDDO
         ENDDO

#ifdef openmp
!$OMP END PARALLEL DO
#endif

      ENDIF

! Inverse Fourier transform to obtain magnetic field:

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
                     CXZB(I,J,K,M)=CXZW(I,J,K,M)
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
                     CXZB(I,J,K,M)=FAC*CXZW(I,J,K,M)
                  ENDDO
               ENDDO
            ENDDO

#ifdef openmp
!$OMP END PARALLEL DO
#endif

         ENDIF
      ENDDO
      RETURN
    END SUBROUTINE BSELF

!***********************************************************************

    SUBROUTINE DIRECT_CALCB(IX,IY,IZ,NX,NY,NZ,IPBC,DX,AK,AKD,AKD2,GAMMA, &
                           PYDDX,PZDDX,CXZG)
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
         CXZG(NX+1+IPBC*(NX-1),NY+1+IPBC*(NY-1),NZ+1+IPBC*(NZ-1),3)

! local variables

      CHARACTER :: CMSGNM*66
      INTEGER :: I,ICXZG,II,IMIN,J,JCXZG,JJ,JMIN,JPY,JPYM,JPZ,JPZM, &
                 K,KCXZG,KMIN,M
      REAL(WP) :: GAMMAKD4,PHASY,PHASYZ,R,R2,RANGE,RANGE2,RJPY, &
                  RJPZ,T0,T1,T2,X0,X2,X2Y2,XX,Y0,Z0
      REAL(WP) :: X(3)
      COMPLEX(WP) :: CXCOEFF,CXFAC,CXI,CXIKR,CXPHAS,CXZERO
      COMPLEX(WP) :: DCXSUM(6)

#ifdef openmp
      INTEGER NTHREADS,TID
#endif

      SAVE CXZERO,CXI
      DATA CXI/(0._WP,1._WP)/,CXZERO/(0._WP,0._WP)/

!-----------------------------------------------------------------------
! subroutine DIRECT_CALCB

! calculates magnetic field at (I,J,K) due to dipole at (1,1,1)
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
!                   N.B.: IPBC affects dimensioning of CXZG

!        GAMMA    = factor used to assist convergence of sums over
!                   replica dipoles when IPBC > 0
!                   contribution from replica dipoles is suppressed
!                   by factor exp(-gamma*(kr)^4)
!                   typical value gamma = 0.005
!                   The effective
!                   range/d = 1/(gamma*k*d) = 400 if gamma=5e-3 and kd=0.5
!                   range/lambda = 1/(2*pi*gamma) = 31.8
!                   The sums are actually continued out to
!                   r/d = 2*/(gamma*kd) = 800 if gamma=5e-3 , kd=0.5
!                   [screening factor = exp(-16)=1.1e-7]

!        NX,NY,NZ = size of first octant (I = 1 -> NX,
!                                         J = 1 -> NY,
!                                         K = 1 -> NZ)
!        PYDDX     = period of lattice in y direction/d
!        PZDDX     = period of lattice in z direction/d
!        DX(1-3)   = lattice spacing in x,y,z direction, in units of
!                    d = n**(-1/3).
!        AK(1-3)   = k(1-3)*d, where k = k vector in vacuo
!        AKD       = |k|d
!        AKD2      = |kd|^2
!        CXZG      = array with dimensions
!                    (NX+1)*(NY+1)*(NZ+1)*3 if IPBC=0
!                    (2*NX)*(2*NY)*(2*NZ)*3 if IPBC=1

! returns:
!        CXZG(I,J,K,M) = c_M to calculate magnetic field B at (I,J,K)
!                        contributed by dipole at (1,1,1) and
!                        replica dipoles (if IPBC=1)
!                        B_x =           c_1*P_y + c_2*P_z
!                        B_y = -c_1*P_x          + c_3*P_z
!                        B_z = -c_2*P_x -c_3*P_y

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
! 12.07.07 (IW)  created DIRECT_CALCB from DIRECT_CALC
! 12.07.11 (BTD) v7.3 cosmetic changes to code written by Ian Wong
! 12.12.21 (BTD) bself_v3
!                * added comments
!                * corrected error when used for IPBC>0
! end history
!-----------------------------------------------------------------------

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
! Compute 3 independent elements of 3x3 anti-symmetric matrix C_jk, where
! C_jk*P_k = magnetic field at location j due to dipole P_k at location

! C_jk = ( 0   c_1  c_2)
!        (-c_1  0   c_3)
!        (-c_2 -c_3  0 )_jk

! initialize CXG(I,J,K,M) = c_M for magnetic field at (I,J,K)
!                           produced by a dipole at (1,1,1)
!                           and replica dipoles (if PYD or PYZ are
!                           nonzero).

! need to calculate this for all (I,J,K) for one octant:
! I running from 1 to NX, J from 1 to NY, K from 1 to NZ

! X0,Y0,Z0 = X(I,J,K) - X(1,1,1) = vector from dipole location (1,1,1)
!                                  to point where B is to be calculated
! IX = +1 -> IMIN = 1
! IX = -1 -> IMIN = 2-NX   (1-IMIN = NX-1)
!
! similarly for JMIN and KMIN

      IMIN=1+(1-NX)*(1-IX)/2
      JMIN=1+(1-NY)*(1-IY)/2
      KMIN=1+(1-NZ)*(1-IZ)/2

! Determine elapsed cpu time for these sums

      CALL CPU_TIME(T0)

#ifdef openmp
! fork a team of threads 
!$OMP PARALLEL               &
!$OMP&   PRIVATE(NTHREADS,TID)

      TID=OMP_GET_THREAD_NUM()
!      WRITE(0,*)'bself_v4 direct_calcb ckpt 1: hello world from thread = ',TID

! only master thread does this:

      IF(TID.EQ.0)THEN
         NTHREADS=OMP_GET_NUM_THREADS()
         WRITE(CMSGNM,FMT='(A,I3)')'number of OpenMP threads =',NTHREADS
         CALL WRIMSG('DIRECT_CALCB',CMSGNM)
      ENDIF
!$OMP END PARALLEL
#endif

#ifdef openmp
! 2012.04.27 (BTD) added R to private variables
!                  removed FLUSH directive (should not be needed)
!$OMP PARALLEL                                   &
!$OMP&   PRIVATE(I,II,J,JPY,JPZ,K,M)             &
!$OMP&   PRIVATE(JCXZG,JPZM,KCXZG)               &
!$OMP&   PRIVATE(PHASY,PHASYZ,R,R2,RJPY,RJPZ) &
!$OMP&   PRIVATE(X,X0,X2,X2Y2,XX,Y0,Z0)       &
!$OMP&   PRIVATE(CXFAC,CXIKR,CXPHAS,DCXSUM)
!$OMP DO
#endif

      DO K=KMIN,NZ              !- loop over K

         Z0=REAL(K-1,KIND=WP)*DX(3)
         IF(K>0)THEN
            KCXZG=K
         ELSE
            KCXZG=2*NZ+K
         ENDIF
         DO J=JMIN,NY              !-- loop over J
            Y0=REAL(J-1,KIND=WP)*DX(2)
            IF(J>0)THEN
               JCXZG=J
            ELSE
               JCXZG=2*NY+J
            ENDIF
            DO I=IMIN,NX              !--- loop over I

! for first dipole, determine time to sum over replicas

               IF(I.EQ.IMIN.AND.J.EQ.JMIN.AND.          &
                  K.EQ.KMIN.AND.MYID==0)CALL CPU_TIME(T1)

               X0=REAL(I-1,KIND=WP)*DX(1)
               IF(I>0)THEN
                  ICXZG=I
               ELSE
                  ICXZG=2*NX+I
               ENDIF
               X(1)=X0
               X2=X(1)*X(1)
               DO M=1,3
                  DCXSUM(M)=CXZERO
               ENDDO

! JPY=0, JPZ=0 gives B field from dipole at (1,1,1)
! general JPY,JPZ gives B field from dipole at
! (1,JPY*NPY+1,JPZ*NPZ+1)
! replica dipoles have same magnitude of polarization as dipole (1,1,1),
! but different phase.
! PHASYZ = phase of replica dipole - phase of dipole (1,1,1)

               DO JPY=-JPYM,JPYM         !---- loop over JPY
                  RJPY=REAL(JPY,KIND=WP)*PYDDX*DX(2)
                  X(2)=Y0-RJPY
                  X2Y2=X2+X(2)*X(2)
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

                     RJPZ=REAL(JPZ,KIND=WP)*PZDDX*DX(3)
                     X(3)=Z0-RJPZ
                     R2=X2Y2+X(3)*X(3)

! skip the self-interaction (R=0) case (I=J=K=1 and JPY=JPZ=0)

                     IF(R2>1.E-6_WP)THEN

                        R=SQRT(R2)

! PHASYZ = phase at (1,JPY*NPY+1,JPZ*NPZ+1) - phase at (1,1,1)

                        PHASYZ=PHASY+AK(3)*RJPZ

!     k^2                    1
! B = --- * exp(ikr) * (1 - --- ) * r x p
!     r^2                   ikr 

! p = exp(i*phasyz) * p(1,1,1)

!     k^2                                    1
!   = --- * exp(ikr) * exp(i*phasyz) * (1 - --- ) * r x p(1,1,1)
!     r^2                                   ikr

! f(r) = (k/r)^2 * exp(ikr+i*phasyz) * [1 - 1/(ikr)]

! B = f(r) * r x p(1,1,1)     

! B_x =  c_1*p_y + c_2*p_z
! B_y = -c_1*p_x + c_3*p_z
! B_z = -c_2*p_x - c_3*p_y

! c_1 = -f * r_z
! c_2 =  f * r_y
! c_3 = -f * r_x

! include artificial factor exp[-(gamma*kr)^4] to assist convergence

                        CXIKR=CXI*AKD*R
                        CXPHAS=EXP(CXIKR+CXI*PHASYZ-GAMMAKD4*R2*R2)
                        CXFAC=(1._WP-(1._WP/CXIKR))
			CXCOEFF=AKD2*CXFAC*CXPHAS/R2

                        DCXSUM(1)=DCXSUM(1)-CXCOEFF*X(3)
			DCXSUM(2)=DCXSUM(2)+CXCOEFF*X(2)
			DCXSUM(3)=DCXSUM(3)-CXCOEFF*X(1)
                      ENDIF   !----- endif (R2 > 1e-6)
                  ENDDO                     !----- end loop over JPZ
               ENDDO                     !---- end loop over JPY
               DO M=1,3
                  CXZG(ICXZG,JCXZG,KCXZG,M)=DCXSUM(M)
               ENDDO
               IF(I==IMIN.AND.J==JMIN.AND.K==KMIN.AND.MYID==0)THEN

! NB: this call to CPU_TIME occurs within the first thread
!     (I==IMIN, J==JMIN, K==KMIN)

                  CALL CPU_TIME(T2)
                  T2=REAL((NX-IMIN+1)*(NY-JMIN+1)*(NZ-KMIN+1),KIND=WP)*(T2-T1)
                  IF(T2<600._WP)THEN
                     WRITE(CMSGNM,FMT='(A,F8.2,A)')                          &
                        'Estimated total cputime required by DIRECT_CALCB=', &
                        T2,' cpu-sec'
                  ELSEIF(T2>=600._WP.AND.T2<3600._WP)THEN
                     WRITE(CMSGNM,FMT='(A,F8.2,A)')                          &
                        'Estimated total cputime required by DIRECT_CALCB=', &
                        (T2/60._WP),' cpu-min'
                  ELSE
                     WRITE(CMSGNM,FMT='(A,F8.2,A)')                          &
                        'Estimated total cputime required by DIRECT_CALCB=', &
                        (T2/3600._WP),' cpu-hr'
                  ENDIF
                  CALL WRIMSG('DIRECT_CALCB',CMSGNM)
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
                  CALL WRIMSG('DIRECT_CALCB',CMSGNM)
               ENDIF
            ENDDO                     !--- end loop over I
         ENDDO                     !-- end loop over J

! 2012.04.27 (BTD) it is not necessary to have a flush(hxzc) here
!                  flush is implied by both $omp END DO and 
!                  $omp END PARALLEL directives

      ENDDO                     !- end loop over K

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif

!*** diagnostic
!      write(0,*)'bself_v4 direct_calcb ckpt 2'
!***
      CALL CPU_TIME(T2)
!*** diagnostic
!      write(0,*)'bself_v4 direct_calcb ckpt 3'
!***
      T2=T2-T0
      IF(MYID==0)THEN
         IF(T2<600.D0)THEN
            WRITE(CMSGNM,FMT='(A,F8.2,A)')                    &
                  'Actual cputime to complete DIRECT_CALCB=', &
                  T2,' cpu-sec'
         ELSEIF(T2>=600._WP.AND.T2<3600._WP)THEN
            WRITE(CMSGNM,FMT='(A,F8.2,A)')                    &
                  'Actual cputime to complete DIRECT_CALCB=', &
                  (T2/60._WP),' cpu-min'
         ELSE
            WRITE(CMSGNM,FMT='(A,F8.2,A)')                    &
                  'Actual cputime to complete DIRECT_CALCB=', &
                  (T2/3600._WP),' cpu-hr'
         ENDIF
         CALL WRIMSG('DIRECT_CALCB',CMSGNM)
      ENDIF

!*** diagnostic
!      write(0,*)'bself_v4 direct_calcb ckpt 4'
!***

! If IPBC=1: set the elements with ICXZG=NX+1 or JCXZG=NY+1
!            or KCXZG=NZ+1 to zero

      IF(IPBC==1)THEN

         DO M=1,3

#ifdef openmp
!$OMP PARALLEL                    &
!$OMP&   PRIVATE(ICXZG,JCXZG,KCXZG)
!$OMP DO
#endif

            DO KCXZG=1,2*NZ
               DO JCXZG=1,2*NY
                  CXZG(NX+1,JCXZG,KCXZG,M)=CXZERO
               ENDDO
            ENDDO

#ifdef openmp
!$OMP END DO
!$OMP DO
#endif

            DO KCXZG=1,2*NZ

               DO ICXZG=1,2*NX
                  CXZG(ICXZG,NY+1,KCXZG,M)=CXZERO
               ENDDO

            ENDDO

#ifdef openmp
!$OMP END DO
!$OMP DO
#endif

            DO JCXZG=1,2*NY

               DO ICXZG=1,2*NX
                  CXZG(ICXZG,JCXZG,NZ+1,M)=CXZERO
               ENDDO

            ENDDO

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif

         ENDDO   ! enddo M=1,3
 
      ENDIF   ! endif (IPBC==1)
      RETURN
 6700 FORMAT(F8.2,' cpu-sec')
 6701 FORMAT(F8.2,' cpu-min')
 6702 FORMAT(F8.2,' cpu-hr')
    END SUBROUTINE DIRECT_CALCB

      
