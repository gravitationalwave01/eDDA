    SUBROUTINE DIAGL(U,V,IPAR)
      USE DDPRECISION,ONLY : WP
      USE DDCOMMON_2,ONLY : CXADIA
      USE DDCOMMON_6,ONLY : GAMMA,PYD,PZD,MXNATF,MXNXF,MXNYF,MXNZF,NAT,NAT3, &
                            NAT0,NX,NY,NZ,MXN3F,IDVOUT,IPBC
      IMPLICIT NONE

! Arguments:
      COMPLEX (WP) :: U(*),V(*)
      INTEGER :: IPAR(*)

! cg package
! Left preconditioning subroutine - division by diagonal
! History:
! Originally written by P.J. Flatau, 1993
! 00.06.26 (BTD) Removed MXMEMF from COMMON/M6/
! 04.03.04 (BTD) Added NPY,NPZ to COMMON/M6/ to support periodic
!                boundary conditions
! 05.06.16 (BTD) Replaced integer NPY,NPZ by real PYD,PZD in /M6/
! 06.09.28 (BTD) Added IPBC to COMMON/M6/
! 07.06.30 (BTD) moved CMDFFT from COMMON/M6/... CMDFFT
!                to COMMON/M8/CMDFFT
!                moved PYD,PZD to beginning of COMMON/M6/
! 07.08.04 (BTD) Version 7.0.3
!                * replaced COMMON/M2/ with MODULE DDCOMMON_2
!                * replaced COMMON/M6/ with MODULE DDCOMMON_6
!                * removed COMMON/M8/CMDFFT (CMDFFT not used)
! 08.03.11 (BTD) ver7.0.5
!                * added ALPHA to DDCOMMON_6
! 08.04.20 (BTD) * change notation: ALPHA -> GAMMA
! end history
! Copyright (C) 1993,2000,2004,2005,2006,2007,2008
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License
!-----------------------------------------------------------------------

!*** diagnostic
!      write(0,*)'entered DIAGL with IPBC=',IPBC
!      write(0,*)'about to call CDIV'
!***
      CALL CDIV(U,V,CXADIA,NAT3)
      RETURN
    END SUBROUTINE DIAGL

!=======================================================================

    SUBROUTINE CDIV(U,V,CXA,N)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

! Arguments:

      INTEGER :: N
      COMPLEX (WP) :: U(*),V(*),CXA(*)

! Local variables:
      INTEGER :: I

!-----------------------------------------------------------------------
! Part of conjugate gradient package
! needed by left preconditioning subroutine
! Copyright (C) 2003 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License
!-----------------------------------------------------------------------

!*** diagnostic
!      write(0,*)'entered CDIV with N=',N,' **inefficient computation**'
!***
      DO I=1,N
!*** diagnostic --- BE SURE TO ELIMINATE THIS TEST!!!
!         if(.not.(abs(cxa(i))>=0.d0))then
!            write(0,*)'division by zero for i=',i,' cxa(i)=',cxa(i)
!         endif
!         if(.not.(abs(u(i))>=0.d0))then
!            write(0,*)'problem for i=',i,' u(i)=',u(i)
!         endif
!***
         V(I)=U(I)/CXA(I)
      ENDDO
      RETURN
    END SUBROUTINE CDIV

!=======================================================================

    SUBROUTINE MATVEC(CXX,CXY,IPAR)
      USE DDPRECISION,ONLY : WP
      USE DDCOMMON_1,ONLY : AKR,DX
      USE DDCOMMON_2,ONLY : CXADIA
      USE DDCOMMON_3,ONLY : CXZC
      USE DDCOMMON_4,ONLY : CXZW
      USE DDCOMMON_5,ONLY : IOCC
      USE DDCOMMON_6,ONLY : GAMMA,PYD,PZD,MXNATF,MXNXF,MXNYF,MXNZF,NAT,NAT3, &
                            NAT0,NX,NY,NZ,MXN3F,IDVOUT,IPBC
      USE DDCOMMON_7,ONLY : CXAOFF
      USE DDCOMMON_8,ONLY : CMDFFT
      IMPLICIT NONE

! Arguments:

      COMPLEX(WP) :: CXX(*),CXY(*)
      INTEGER :: IPAR(*)

! Subroutine MATVEC

! matvec  --> cxy = A cxx          'N'
! tmatvec --> cxy = A' cxx         'T'
! cmatvec --> cxy = conjg(A') cxx  'C'


! Notice that DDSCAT has cases 'N' and 'C' implemented in CPROD
! Notice, however, that in DDSCAT 'T'='N'.
! Some CG methods use cases N, C (like Petravic).
! Some CG methods use cases N, T (like PIM)
! Some CG use all three products N, C, T  (parts of CCGPAK)

! Information about matrix A is transfered via commons /M1/-/M7/ with
! the main program. Arrays cxadia, cxzc, cxzwm,  iocc are properly
! dimensioned inside "cprod". You may get warning during the compilation
! about "misaligned common". Ignore these. Also notice mxnatf,
! mxn3f,mxnxf, mxnyf, mxnzf. These are defined by "parameter
! statement" in main code and this is the reason why they have "f" suffi

! ipar - information array, not used, defined for compatibility

! History:
! 97.11.01 (BTD): modified to allow for noncubic lattice, introduced
!                 new vector DX(1-3) via COMMON/M1/
!                 added DX to argument list of CPROD
! 97.12.25 (BTD): added COMMON/M7/CXAOFF to communicate off-diagonal
!                 elements of inverse polarizability tensors.
!                 Added CXAOFF to argument list of CPROD
! 00.06.13 (BTD): Explicit declaration of all variables.
! 00.06.25 (BTD): Removed MXMEMF from COMMON/M6/
! 04.03.04 (BTD): Added NPY,NPZ to COMMON/M6/ to support periodic
!                 boundary conditions
!                 Added NPY,NPZ to argument list of CPROD
! 05.06.16 (BTD): Replaced integer NPY,NPZ by real PYD,PZD
! 06.09.28 (BTD): Version 6.2.3
!                 * Added IPBC to COMMON/M6/
!                 * Added IPBC to argument list of CPROD
! 07.06.30 (BTD): moved CMDFFT from COMMON/M6/... CMDFFT
!                 to COMMON/M8/CMDFFT
!                 moved PYD,PZD to beginning of COMMON/M6/
! 07.08.04 (BTD): Version 7.0.3
!                 * renamed AK(3)->AKR(3) for consistency with COMMON/M1/
!                   in DDSCAT
!                 * replaced COMMON/M1/ with USE MODULE DDCOMMON_1
!                 * replaced COMMON/M2/ with USE MODULE DDCOMMON_2
!                 * replaced COMMON/M3/ with USE MODULE DDCOMMON_3
!                 * replaced COMMON/M4/ with USE MODULE DDCOMMON_4
!                 * replaced COMMON/M5/ with USE MODULE DDCOMMON_5
!                 * replaced COMMON/M6/ with USE MODULE DDCOMMON_6
!                 * replaced COMMON/M7/ with USE MODULE DDCOMMON_7
!                 * replaced COMMON/M8/ with USE MODULE DDCOMMON_8
! 08.03.11 (BTD): v7.0.5
!                 * added ALPHA to DDCOMMON_6
! 08.03.14 (BTD): * changed declaration of argument IPAR
!                   from COMPLEX(WP) :: IPAR(*)
!                   to INTEGER :: IPAR(*)
!                   NB: dummy argument IPAR is not used
! 08.04.20 (BTD): * change notation: ALPHA -> GAMMA
! end history
! Copyright (C) 1993,1997,2000,2004,2005,2006,2007,2008
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
!*** diagnostic
!      write(0,*)'matvec ckpt 1, IPBC=',IPBC
!      write(0,*)'               about to call CPROD with CWHAT=N'
!***
      CALL CPROD(AKR,GAMMA,DX,CMDFFT,'N',CXADIA,CXAOFF,CXX,CXY,CXZC,CXZW,     &
                 IDVOUT,IOCC,MXN3F,MXNATF,MXNXF,MXNYF,MXNZF,NAT,NAT0,NAT3,NX, &
                 NY,NZ,IPBC,PYD,PZD)
!*** diagnostic
!      write(0,*)'matvec ckpt 2'
!***
      RETURN
    END SUBROUTINE MATVEC

!=======================================================================

    SUBROUTINE CMATVEC(CXX,CXY,IPAR)
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
      IMPLICIT NONE

! Arguments:

      COMPLEX(WP) :: CXX(*),CXY(*),IPAR(*)

!-----------------------------------------------------------------------
! Subroutine CMATVEC

! matvec  --> cxy = A cxx
! tmatvec --> cxy = A' cxx
! cmatvec --> cxy = conjg(A') cxx
! History:
! 97.11.01 (BTD): modified to allow for noncubic lattice, introduced
!                 new vector DX(1-3) via COMMON/M1/
!                 added DX to argument list of CPROD
! 97.12.25 (BTD): added COMMON/M7/CXAOFF to communicate off-diagonal
!                 elements of inverse polarizability tensors.
!                 Added CXAOFF to argument list of CPROD
! 00.06.13 (BTD): Explicit declaration of all variables.
! 00.06.26 (BTD): Removed MXMEMF from COMMON/M6/
! 04.03.04 (BTD): Added NPY,NPZ to COMMON/M6/ to support periodic
!                 boundary conditions
!                 Added NPY,NPZ to argument list of CPROD
! 05.06.16 (BTD): Replaced integer NPY,NPZ by real PYD,PZD in
!                 COMMON/M6/ and argument list of CPROD
! 06.09.29 (BTD): Version 6.2.3:
!                 * Added IPBC to COMMON/M6/
!                 * Added IPBC to argument list of CPROD
! 07.06.30 (BTD): moved CMDFFT from COMMON/M6/... CMDFFT
!                 to COMMON/M8/CMDFFT
!                 moved PYD,PZD to beginning of COMMON/M6/
! 07.08.04 (BTD): Version 7.0.3:
!                 * renamed AK(3)->AKR(3) for consistency with
!                   COMMON/M1/ in DDSCAT
!                 * replaced COMMON/M1/ with USE MODULE DDCOMMON_1
!                 * replaced COMMON/M2/ with USE MODULE DDCOMMON_2
!                 * replaced COMMON/M3/ with USE MODULE DDCOMMON_3
!                 * replaced COMMON/M4/ with USE MODULE DDCOMMON_4
!                 * replaced COMMON/M5/ with USE MODULE DDCOMMON_5
!                 * replaced COMMON/M6/ with USE MODULE DDCOMMON_6
!                 * replaced COMMON/M7/ with USE MODULE DDCOMMON_7
!                 * replaced COMMON/M8/ with USE MODULE DDCOMMON_8
! 08.03.11 (BTD): v7.0.5:
!                 * added ALPHA to DDCOMMON_6
! 08.04.20 (BTD): * changed notation: ALPHA -> GAMMA
! end history
! Copyright (C) 1993,1997,2000,2004,2005,2006,2007,2008
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
!*** diagnostic
!      write(0,*)'entered cmatvec with IPBC=',IPBC
!***
      CALL CPROD(AKR,GAMMA,DX,CMDFFT,'C',CXADIA,CXAOFF,CXX,CXY,CXZC,CXZW,     &
                 IDVOUT,IOCC,MXN3F,MXNATF,MXNXF,MXNYF,MXNZF,NAT,NAT0,NAT3,NX, &
                 NY,NZ,IPBC,PYD,PZD)
      RETURN
    END SUBROUTINE CMATVEC

!=======================================================================

      SUBROUTINE CPROD(AKR,GAMMA,DX,CMETHD,CWHAT,CXADIA,CXAOFF,CXX,CXY,CXZC, &
                       CXZW,IDVOUT,IOCC,MXN3,MXNAT,MXNX,MXNY,MXNZ,NAT,NAT0,  &
                       NAT3,NX,NY,NZ,IPBC,PYD,PZD)
      USE DDPRECISION,ONLY : WP
       
      IMPLICIT NONE

! Arguments:

      INTEGER :: IDVOUT,IPBC,MXN3,MXNAT,MXNX,MXNY,MXNZ,NAT,NAT0,NAT3, &
         NX,NY,NZ
      CHARACTER :: CMETHD*6,CWHAT*1
      INTEGER*2 :: &
         IOCC(MXNAT)
      REAL(WP) :: GAMMA,PYD,PZD
      REAL(WP) :: &
         AKR(3),  &
         DX(3)
      COMPLEX(WP) ::                                                 &
         CXADIA(MXN3),                                               &
         CXAOFF(MXN3),                                               &
         CXX(MXN3),                                                  &
         CXY(MXN3),                                                  &
         CXZC(NX+1+IPBC*(NX-1),NY+1+IPBC*(NY-1),NZ+1+IPBC*(NZ-1),6), &
         CXZW(MXNX*MXNY*MXNZ,*)

! Local variables:

      INTEGER :: J1,J2,J3
      REAL(WP) :: AKD,DTIME

! Intrinsic functions:

      INTRINSIC CONJG

!***********************************************************************
! Subroutine CPROD

! Given:
!   AKR(1-3)= k(1-3)*d , where k=k vector in vacuo
!                       d=effective lattice spacing=(dx*dy*dz)**(1/3)
!   GAMMA  = parameter controlling PBC summations over replica dipoles
!            interaction is suppressed by factor exp(-(gamma*k*r)^4)
!            and summations are carried out to r=2/(gamma*k)
!   DX(1-3)=(dx/d, dy/d, dz/d)
!                 where dx,dy,dz=lattice spacing in x,y,z directions.
!                 By definition of effective lattice spacing d, must
!                 have DX(1)*DX(2)*DX(3)=1
!   CMETHD         = character variable specifying method used for FFT
!                    evaluation (see ESELF)
!   CWHAT          = 'N' or 'C' determines what is to be done (see below
!   CXADIA(1-NAT3) = diagonal elements of A(j,k)
!                    (diagonal elements of inverse of polarizability
!                     tensor for dipoles 1-NAT)
!   CXAOFF(1-NAT3) = off-diagonal elements of 3x3 diagonal blocks of
!                    A(j,k): CXAOFF(J,1-3)=a_{23},a_{31},a_{12} from 3x3
!                    matrix a_{ij} for dipole J.
!                    Recall that a_{ij} is inverse of 3x3 polarizability
!                    tensor alpha_{ij} for dipole J.
!   CXX(1-NAT3)    = complex vector
!   CXZC           = array of Fourier transformed Green function
!                    coefficients used internally by ESELF, and generate
!                    by ESELF each time called with new k vector
!                    If IPBC=0: size = (NX+1)*(NY+1)*(NZ+1)*6
!                    If IPBC=1: size = (2*NX)*(2*NY)*(2*NZ)*6
!                    *****Not to be overwritten between calls to ESELF**
!   CXZW           = complex workspace required by ESELF
!   IDVOUT      = device number for output
!   IOCC(1-NAT) = 0 or 1 if site is unoccupied or occupied
!   MXN3,MXNX,MXNY,MXNZ = dimensioning information (see DDSCAT)
!   NAT         = number of sites in extended target
!   NAT0        = number of occupied sites in extended target
!   NAT3        = 3*NAT
!   NX,NY,NZ    = extended target size is (NX*DX(1)) by (NY*DX(2)) by
!                 (NZ*DX(3))
!   PYD         = 0. to do isolated target
!               = (periodicity in y direction)/d(2) for periodic b.c.
!   PZD         = 0. to do isolated target
!               = (periodicity in z direction)/d(3) for periodic b.c.

! Returns:
!   CXY(1-NAT3)
!       If CWHAT = 'N' : CXY(j)=A(j,k)*CXX(k) (sum over k)
!                = 'C' : CXY(j)=conjg(A(k,j))*CXX(k) (sum over k)

! It is assumed that matrix A(j,k) is symmetric
! Diagonal elements of A(j,k) are in vector CXADIA(j)
! Off-diagonal terms coupling dipole component x with direction y, etc
! are stored in vector CXAOFF(j), CXAOFF(j+NAT), CXAOFF(j+2*NAT)
! Other off-diagonal elements of A(j,k) are produced internally by
!   subroutine ESELF, which uses 3d-FFT to accomplish convolution
!   of CXX with A-matrix terms coupling dipole j with other dipoles
!   (i.e., E field at j due to other dipoles)

! Original version of CPROD created by P.J.Flatau, Colorado State Univ.
! Subsequently modified by B.T.Draine, Princeton Univ. Obs.
! History:
! 90.11.04 (BTD): Major modification to incorporate FFT method (FFT
!                 computations are handled by ESELF).
! 90.11.08 (BTD): Remove DO...ENDDO structures; renumber DO...CONTINUEs
! 90.11.21 (BTD): Remove "CCGMVM" option: only FFT stuff retained.
! 90.12.03 (BTD): No longer necessary to change ordering of CXX and CXY
! 90.12.05 (BTD): Changed call to ESELF to eliminate MXNX,MXNY,MXNZ
! 97.11.01 (BTD): Added comments and modified to accomodate noncubic
!                 lattice. Added DX to argument list for CPROD and to
!                 argument list in call to ESELF
! 97.12.25 (BTD): Added CXAOFF to argument list.
! 97.12.28 (BTD): Modified to use CXAOFF in calculation of CXY
! 98.01.01 (BTD): Modified coding to use indices J1,J2,J3 in loops
!                 computing with CXAOFF
! 00.06.13 (BTD): Explicit declaration of all variables.
! 04.03.05 (BTD): Added NPY,NPZ to argument list to support periodic
!                 boundary conditions
! 05.06.16 (BTD): Replaced integer NPY,NPZ by real PYD,PZD
!                 in argument list and in argument list of CPROD
! 05.08.04 (BTD): Added AKR(1-3) to argument list of ESELF
! 06.09.23 (BTD): Added comments
! 06.09.29 (BTD): Version 6.2.3:
!                 * Added IPBC to argument list of CPROD
!                 * Modified dimensioning of array CXZC to use larger
!                   dimension when IPBC=1
!                 * Added IPBC to argument list of ESELF (two places)
! 07.08.04 (BTD): Version 7.0.3:
!                 * renamed AK(3)->AKR(3) for consistency with
!                   COMMON/M1/ in DDSCAT
! 08.01.13 (BTD): cosmetic changes
! 08.03.11 (BTD): v7.0.5
!                 * added ALPHA to argument list
!                 * added ALPHA to argument list of ESELF
! 08.04.20 (BTD): * changed notation: ALPHA -> GAMMA
! 08.05.12 (BTD): v7.0.6
!                 * added calls to TIMEIT to time ESELF
! end history
! Copyright (C) 1993,1997,1998,2000,2004,2005,2006,2007,2008
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!-----------------------------------------------------------------------
!*** diagnostic
!      write(0,*)'cprod ckpt 1 with nx,ny,nz,ipbc=',nx,ny,nz,ipbc
!      write(0,*)'              cwhat=',cwhat
!***

      AKD=SQRT(AKR(1)*AKR(1)+AKR(2)*AKR(2)+AKR(3)*AKR(3))

!*** diagnostic
!      write(0,*)' in cprod, akd=',akd
!***
! We assume that input vector CXX has zeroes for elements corresponding
! to vacuum sites.

      IF(CWHAT=='N')THEN
!***********************************************************************

!        CXX(1-NAT3) = polarizations

!*** diagnostic
!         write(0,*)'cprod ckpt 2: about to call eself from cprod'
!***
! 08.10.02 suppress this...
!         CALL TIMEIT('ESELFCP',DTIME)
         CALL ESELF(CMETHD,CXX,NX,NY,NZ,IPBC,GAMMA,PYD,PZD,AKR,AKD,DX, &
                    CXZC,CXZW,CXY)
! 08.10.02 suppress this...
!         CALL TIMEIT('ESELFCP',DTIME)

!*** diagnostic
!         write(0,*)'cprod ckpt 3'
!***
!        CXY(1-NAT3) is now Efield produced by other dipoles (including
!                    replica dipoles if PYZ or PZD are nonzero)

!*** diagnostic
!        write(0,*)'returned to cprod from eself'
!        write(0,*)'check cxy for NaN...'
!        j2=0
!        do j1=1,nat3
!           if(.not.(abs(cxy(j1))>=0.d0))then
!              write(0,*)'j1=',j1,' cxy(j1)=',cxy(j1)
!              j2=j2+1
!           endif
!        enddo
!        write(0,*)'cxy checked for NaN'
!        write(0,*)'j2=',j2,' instances'
!***
!***********************************************************************
! ESELF omitted contribution from diagonal terms: add these.
! Also remember that ESELF computes -A_jk*x_k

         DO J1=1,NAT3
            CXY(J1)=CXADIA(J1)*CXX(J1)-CXY(J1)
         ENDDO

!*** diagnostic
!        write(0,*)'in cprod ckpt 700, check cxy for NaN...'
!        do j1=1,nat3
!           if(.not.(abs(cxy(j1))>=0.d0))then
!              write(0,*)'j1=',j1,' cxy(j1)=',cxy(j1)
!           endif
!        enddo
!        write(0,*)'cxy checked for NaN'
!***        
! Add contribution from off-diagonal elements of 3x3 diagonal blocks
! Version below assumes that off-diagonal elements a_ij are stored as
! CXAOFF(J,1-3)=(a_23,a_31,a_12) for dipole J.
! 97.12.28 (BTD):

         DO J1=1,NAT
            J2=J1+NAT
            J3=J2+NAT
            CXY(J1)=CXY(J1)+CXAOFF(J2)*CXX(J3)+CXAOFF(J3)*CXX(J2)
            CXY(J2)=CXY(J2)+CXAOFF(J3)*CXX(J1)+CXAOFF(J1)*CXX(J3)
            CXY(J3)=CXY(J3)+CXAOFF(J1)*CXX(J2)+CXAOFF(J2)*CXX(J1)
         ENDDO
!*** diagnostic
!        write(0,*)'in cprod ckpt 701, check cxy for NaN...'
!        do j1=1,nat3
!           if(.not.(abs(cxy(j1))>=0.d0))then
!              write(0,*)'j1=',j1,' cxy(j1)=',cxy(j1)
!           endif
!        enddo
!        write(0,*)'cxy checked for NaN'
!***        

      ELSEIF(CWHAT=='C')THEN

! Need to temporarily replace CXX by conjg(CXX)

         DO J1=1,NAT3
            CXX(J1)=CONJG(CXX(J1))
         ENDDO

!***********************************************************************

         CALL ESELF(CMETHD,CXX,NX,NY,NZ,IPBC,GAMMA,PYD,PZD,AKR,AKD,DX, &
                    CXZC,CXZW,CXY)

!***********************************************************************

! ESELF omitted contribution from diagonal terms: add these.
! Also remember that ESELF computed -A_jk*conjg(x_k), while
! we seek to compute conjg(A_kj)* x_k.  Since A_kj=A_jk, we have
! conjg(A_kj)*x_k=conjg(A_jk*conjg(x_k))

         DO J1=1,NAT3
            CXY(J1)=CXADIA(J1)*CXX(J1)-CXY(J1)
         ENDDO

! Add contribution from off-diagonal elements of 3x3 diagonal blocks
! It is assumed that we have stored off-diagonal elements a_ij as
! CXAOFF(J,1-3)=(a_23,a_31,a_12) for element J
! Note that we could *probably* speed up the code by
! (1) "reducing" product vector CXY immediately after calling ESELF
! (2) restricting CXADIA to occupied sites
! (3) restricting CXAOFF to occupied sites
! but this would require modifications outside this routine, and
! therefore requires further study...
! 97.12.28 (BTD)

         DO J1=1,NAT
            J2=J1+NAT
            J3=J2+NAT
            CXY(J1)=CXY(J1)+CXAOFF(J2)*CXX(J3)+CXAOFF(J3)*CXX(J2)
            CXY(J2)=CXY(J2)+CXAOFF(J3)*CXX(J1)+CXAOFF(J1)*CXX(J3)
            CXY(J3)=CXY(J3)+CXAOFF(J1)*CXX(J2)+CXAOFF(J2)*CXX(J1)
         ENDDO

! Now compute conjugate of CXY, and restore CXX:

         DO J1=1,NAT3
            CXY(J1)=CONJG(CXY(J1))
            CXX(J1)=CONJG(CXX(J1))
         ENDDO
      ELSE
         WRITE (IDVOUT,FMT=*) ' Error in CPROD: CWHAT= ',CWHAT
      ENDIF

!*** Now must "clean" product vector CXY:
!    (zero out elements corresponding to vacuum sites)

!*** diagnostic
!      write(0,*)'about to call nuller to clean cxy...'
!      write(0,*)'mxnat=',mxnat,' mxn3=',mxn3
!      write(0,*)' nat=',nat,' nat0=',nat0
!***

      IF(NAT0<NAT)CALL NULLER(CXY,IOCC,MXNAT,MXN3,NAT)

!*** diagnostic
!      write(0,*)'in cprod, ckpt 707, check xy for NaN'
!      j2=0
!      do j1=1,nat3
!         if(.not.(abs(cxy(j1))>=0.d0))then
!            j2=j2+1
!            write(0,*)'j1=',j1,' cxy(j1)=',cxy(j1)
!         endif
!      enddo
!      write(0,*)'cxy checked for NaN: j2=',j2,' instances'
!***
       RETURN
    END SUBROUTINE CPROD
