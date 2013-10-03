    SUBROUTINE SCAT(AK_TF,AKS_TF,DX,EM1_TF,EM2_TF,E02,ETASCA,CMDTRQ,CBKSCA, &
                    CSCA,CSCAG,CSCAG2,CTRQIN,CTRQSC,CXE_TF,CXE01_TF,CXF1L,  &
                    CXF2L,CXP_TF,CXSCR1,CXSCR2,CXSCR3,CXSCR4,MXN3,MXNAT,    &
                    MXSCA,MYID,JPBC,NAT,NAT3,NAVG,NDIR,SCRRS1,SCRRS2,IXYZ,X0)
!-------------------------------- v5 -----------------------------------------
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

! Arguments:

      CHARACTER :: CMDTRQ*6
      REAL(WP) :: CBKSCA,CSCA,CSCAG2,E02,ETASCA
      INTEGER :: JPBC,MXN3,MXNAT,MXSCA,MYID,NAT,NAT3,NAVG,NDIR
      INTEGER ::   &
         IXYZ(NAT,3)
      COMPLEX(WP) ::     &
         CXE_TF(NAT,3),  &
         CXE01_TF(3),    &
         CXF1L(MXSCA),   &
         CXF2L(MXSCA),   &
         CXP_TF(NAT,3),  &
         CXSCR1(MXN3),   &
         CXSCR2(MXN3),   &
         CXSCR3(MXN3),   &
         CXSCR4(MXNAT,3)
      REAL(WP) ::         &
         AK_TF(3),        &
         AKS_TF(3,MXSCA), &
         CSCAG(3),        &
         CTRQIN(3),       &
         CTRQSC(3),       &
         DX(3),           &
         EM1_TF(3,MXSCA), &
         EM2_TF(3,MXSCA), &
         SCRRS1(NAT,3),   &
         SCRRS2(MXNAT),   &
         X0(3)

! Local variables:

      CHARACTER :: CMSGNM*70
      COMPLEX(WP) :: CXI,CXSCL1,CXSCL2,CXSCL3,CXZERO
      COMPLEX(WP) :: CXSCL1_L,CXSCL2_L,CXSCL3_L,CXF1L_L,CXF2L_L, &
                     CXTRM_K  !Art for local private vars
      REAL(WP) :: AFAC,AK2,AK3,AKK,COSPHI,COSTH,DOMEGA,DPHI,EINC,ES2,PI,PHI, &
                  RRR,SI,SINPHI,SINTH,TERM,THETA,THETA0,THETAL,THETAU,W,X
      INTEGER :: ICOSTH,IPHI,J,K,ND,NPHI,NTHETA
      COMPLEX(WP) ::  &
         CXBS(3),     &
         CXES(3),     &
         CXES_L(3),   &  !Art for private use
         CXTRM(3)
      REAL(WP) ::     &
         AKS0(3),     &
         AKS1(3),     &
         AKS2(3),     &
         AKSN(3),     &
         CTRQSCTF(3), &
         TEMP(MXNAT), &   !Art for local temp
         XYZCM(3)

! Intrinsic functions:

      INTRINSIC COS,EXP,SIN,SQRT,CONJG,REAL

!***********************************************************************

! SCAT computes energy radiated by array of oscillating dipoles
! and corresponding scattering properties if oscillations are in
! response to incident E field.

! Note: *** All vectors are given in the same frame ***
!           (e.g., the Target Frame)

! Given:
!     AK_TF(1-3) = (kx,ky,kz)*d = where (kx,ky,kz)=incident k vector in TF
!               d = (dx*dy*dz)**(1/3) = effective lattice spacing
!     AKS_TF(3,MXSCA)=scattered k vectors in TF
!     DX(1-3) = (dx,dy,dz)/d where dx,dy,dz=x,y,z lattice spacing
!     CMDTRQ = 'DOTORQ' to calculate torques
!            = 'NOTORQ' to skip calculation of torques
!     EM1_TF(3,MXSCA)=pol.vector 1 in TF for each scattering direction
!     EM2_TF(3,MXSCA)=pol.vector 2 in TF for each scattering direction
!     E02 = |E_0|^2 , where E_0 = incident complex E field
!     NAT = number of dipoles
!     NAT3 = 3*NAT
!     IXYZ(NAT,1-3) = [x-X0(1)]/d, [y-X0(2)]/d, [z-X0(3)]/d (integers) 
!                     for dipoles 1-NAT
!     X0(1-3) = location/d in TF of lattice site with IXYZ=(0,0,0)
!     CXE_TF(1-NAT,1-3) = components in TF of E field at each dipole at t=0
!     CXE01_TF(1-3) = (complex) reference polarization vector in TF
!     CXP_TF(1-NAT,1-3) = components in TF of polarization of each dipole at t=0
!     NDIR = number of directions at which scattering matrix elements
!            are to be computed
!     JPBC = 0 for finite target
!            1 for target periodic in y_TF direction
!            2 for target periodic in z_TF direction
!            3 for target periodic in y_TF and z_TF directions
!     MYID = id of this thread 

! Returns:

!     CBKSCA   =differential scattering cross section for theta=pi
!               (lattice units)
!     CSCA     = scattering cross section (lattice units)
!     CSCAG(1) = CSCA*<cos(theta)> , where theta=scattering angle
!     CSCAG(2) = CSCA*<sin(theta)cos(phi)>
!     CSCAG(3) = CSCA*<sin(theta)sin(phi)>
!     CSCAG2   = CSCA*<cos^2(theta)>
!     NAVG     = number of directions used for calculating CSCAG(1-3)
!                and CSCAG2
!   CTRQIN(1-3)=cross section for torque on grain due to incident
!               E and B fields,
!               for torque in directions x,y,z assuming incident
!               radiation is along x direction,
!               and "reference pol state" is in y direction.
!   CTRQSC(1-3)=cross section for torque on grain due to scattered
!               E and B fields,
!               for torque in directions x,y,z assuming incident
!               radiation is along x direction,
!               and "reference pol state" is in y direction
!               Here "torque cross section" is defined so that
!               torque(1-3)=(momentum flux)*(1/k)*C_torque(1-3)
!               Total torque cross section=CTRQIN+CTRQSC
!   CXF1L(NDIR)=f_{1L} for direction NDIR
!   CXF2L(NDIR)=f_{2L} for direction NDIR
!               where f_{1L} connects incident polarization L to
!               outgoing theta polarization (E in scattering plane)
!               f_{2L} connects incident polarization L to outgoing phi
!               polarization (E perpendicular to scattering plane)
!               This is same f_{ml} notation as used by Draine (1988;
!               Ap.J., 333, 848).

! Notes:
!     AKS0(1-3)=k_s*d=scattered k vector (lattice units)
!     AKSN(1-3)=k_s/|k_s|=unit vector in scattering direction
! Scratch vectors introduced to permit vectorization:
!     SCRRS1 = real vector of length.GE.3*NAT
!     SCRRS2 = real vector of length.GE.NAT
!     CXSCR1 = complex scratch array of length 3*MXNAT
!     CXSCR2 = complex scratch array of length 3*MXNAT
!     CXSCR3 = complex scratch array of length 3*MXNAT
!     CXSCR4 = complex scratch array of length 3*MXNAT

! B.T.Draine, Princeton Univ. Obs., 87.01.04
! History:
! 88.06.29 (BTD): modified
! 89.11.28 (BTD): modified
! 90.11.02 (BTD): modified to enable use with vacuum sites.
! 90.11.20 (BTD): corrected error in calculation of f_ml (missing
!                 factor of k).
! 90.12.03 (BTD): changed ordering of XYZ0 and CXP
!                 eliminated cxe0 from argument list
! 90.12.10 (BTD): replaced XYZ0 with IXYZ to reduce memory use
! 91.08.14 (BTD): add CBKSCA to argument list for SCAT
!                 add code to eliminate sum over phi when theta=0 or pi
!                 add code to calculate CBKSCA
! 94.04.04 (BTD): removed superfluous CMPLX( ) operation from 2 lines
!                 in DO 200 loop. Also added new variable TERM to move
!                 some computation out of last DO loop.
! 95.06.14 (BTD): modified to compute moments <sin(theta)cos(phi)>
!                 and <sin(theta)sin(phi)> of scattered radiation
!                 added code to compute torque on grain relative to
!                 grain centroid XYZCM(1-3)
!                 NOTE: handling of grain centroid is not as efficient
!                 as it should be.  Can be recoded to reduce number
!                 of computations.
! 95.06.15 (BTD): modified to compute absorption contribution to
!                 torque on grain.
!                 added CXE and CTRQIN to argument list
! 95.06.15 (BTD): corrected errors found by Joseph Weingartner
! 95.06.16 (BTD): more corrections
! 95.06.20 (BTD): added argument CMDTRQ to allow torque calculations
!                 to be skipped if not desired.
! 95.06.22 (BTD): added factor CXI in terms appearing in evaluation
!                 of scattered torque
! 95.06.29 (BTD): Added loop to compute AKS0(1-3) for special cases
!                 theta-0 and 180. Error found by J. Weingartner.
! 95.07.20 (BTD): Add scratch vector SCRRS2
!                 Add scratch array CXSCR4 to argument list
!                 Rewrite to make torque computation more efficient:
!                 move as many calculations as possible out of loops
!                 over scattering directions:
!                 SCRRS1 now used to store coords relative to centroid
!                 SCRRS2 now used to store (k_s dot r_j)
!                 CXSCR4 is needed to store (p_j cross r_j)
! 95.07.21 (BTD): Rewrite (again!) to make the correspondence
!                 to notation of Draine & Weingartner more apparent.
!                 This involved changing various definitions by
!                 factors of k, and change of sign of CXBS
! 95.07.24 (BTD): Corrected errors in evaluation of both scattering
!                 torque and incident torque
! 95.07.28 (BTD): Reenabled computation of CBKSCA
! 97.11.02 (BTD): Modified to allow use of noncubic lattice:
!                 added DX(1-3) to argument list and used DX in
!                 evaluation of phase factors
! 98.08.10 (BTD): Added DX(1-3) in calculation of scattering amplitudes
!                 FML for selected scattering directions.
! 03.10.21 (BTD): Changed selection of angles used in calculation of
!                 radiation force and torque.
! 03.10.23 (BTD): New method works well.  Number of scattering angles
!                 used in calculation of radiation force (and torque)
!                 will now be set automatically.
!                 Local variable alpha can be used to adjust resolution
!                 if required (alpha now set to 0.5).
!                 Eliminate ICTHM and IPHIM from argument list -- no longer
!                 to be used.
!                 Add CSCAG2 and NAVG to argument list.
! 07.06.20 (BTD): Add X0 to argument list
!                 Modify to use X0 in evaluating phase factors
! 07.06.22 (BTD): Removed THETAN and PHIN from argument list
!                 They are not needed, because scattering directions
!                 are specified by vector AKS
! 07.09.11 (BTD): Changed IXYZ from INTEGER*2 to INTEGER
! 08.01.06 (BTD): Cosmetic changes
! 08.03.14 (BTD): v7.0.5
!                 changed declaration
!                 IXYZ(MXNAT,3) -> IXYZ(NAT,3)
! 08.05.12 (BTD): v7.0.6 changes introduced by Art Lazanoff, NASA Ames:
!                 introduced local variables CXSCL1_L,CXSCL2_L,CXSCL3_L,
!                 and TEMP(MXNAT) 
!                 added numerous OpenMP directives such as
!                    #ifdef scat_omp
!                    !$omp parallel
!                    #endif
!                 which will be compiled only if compiled with
!                 preprocessor flag -fpp -Dscat_omp 
! 08.07.01 (BTD): v3: this is a copy of scat_asl2.f90 as of 08.06.20
! 08.08.06 (BTD): change #ifdef debug -> #ifdef openmp
!                 so that all omp directives are enabled by single
!                 preprocessor flag
! 11.11.15 (BTD): v5: change notation
!                 AK    -> AK_TF
!                 AKS   -> AKS_TF
!                 CXE   -> CXE_TF
!                 CXE00 -> CXE01_TF
!                 CXP   -> CXP_TF
!                 EM1   -> EM1_TF
!                 EM2   -> EM2_TF
! 12.04.27 (BTD): modifications to omp directives
! end history
! Copyright (C) 1993,1994,1995,1997,1998,2003,2007,2008,2011,2012
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
!*** diagnostic
!      write(0,*)'scat ckpt 0: entered SCAT with JPBC=',JPBC
!***
      CXZERO=(0._WP,0._WP)
      CXI=(0._WP,1._WP)
      PI=4._WP*ATAN(1._WP)

      AK2=AK_TF(1)*AK_TF(1)+AK_TF(2)*AK_TF(2)+AK_TF(3)*AK_TF(3)
      AKK=SQRT(AK2)
      AK3=AKK*AK2

!=======================================================================

!         Calculate scattering and possibly torque for finite target

!    Method used to choose scattering directions for calculation
!    of radiation force and torque:

!    Define function
!                        A*(theta/theta0)
!       s(theta)=theta + ----------------
!                        1 + theta/theta0

!    We choose angles theta to be uniformly-distributed in s.
!    With A > 0 this gives increased resolution at small
!    scattering angles, particularly for theta < theta0

!    We take
!       A=1
!       theta0=2*pi/(1.+x)

!    Parameter etasca determines the numbers of angles, chosen so that
!    the largest interval is
!                                 pi/2
!    [Delta theta]_max = etasca * ----
!                                 3+x

!    A reasonable choice is etasca = 1

!    With etasca=1:
!         x=0 --> [Delta theta]_max = 30 deg
!           1                         22.5 deg
!           10                         6.92
!           15                         5 deg

!    With this requirement for [Delta theta]_max,
!    the number of theta values (including 0 and pi) becomes

!                 2*(3+x)       1 + A/(pi+theta0)
!    NTHETA = 1 + ------- * --------------------------
!                 etasca    1 + A*theta0/(pi+theta0)^2

!    Thus, with A=1 and etasca=1
!    x=0: theta0 = 2*pi   NTHETA=7.19818 -> 7
!      1            pi           9.48970 -> 9
!      2          2*pi/3         12.0646 -> 12
!     10          2*pi/11        32.6897 -> 33

! etasca = 1 appears to work well.
! Use smaller value of ETASCA for improved accuracy

! Incident k vector AK_TF defines one axis (for spherical integration)
! Use (real part of) reference polarization vector CXE01_TF to define
! second coordinate axis.
! Construct third coordinate axis by taking cross product of first
! two axes.
! AKS1, AKS2 = vectors of length equal to AK_TF lying along second and
!              third coordinate axes.

!*** diagnostic
!      write(0,*)'scat ckpt 1'
!***
      IF(JPBC==0)THEN

         AKS1(1)=REAL(CXE01_TF(1))
         AKS1(2)=REAL(CXE01_TF(2))
         AKS1(3)=REAL(CXE01_TF(3))

! EINC = magnitude of (real part of) reference electric field

         EINC=SQRT(AKS1(1)**2+AKS1(2)**2+AKS1(3)**2)
         IF(EINC<=0._WP)CALL ERRMSG('FATAL','SCAT',              &
                                    'CXE01_TF IS PURE IMAGINARY!')
         RRR=AKK/EINC
         AKS1(1)=RRR*AKS1(1)
         AKS1(2)=RRR*AKS1(2)
         AKS1(3)=RRR*AKS1(3)
         AKS2(1)=(AK_TF(2)*AKS1(3)-AK_TF(3)*AKS1(2))/AKK
         AKS2(2)=(AK_TF(3)*AKS1(1)-AK_TF(1)*AKS1(3))/AKK
         AKS2(3)=(AK_TF(1)*AKS1(2)-AK_TF(2)*AKS1(1))/AKK

!*** diagnostic
!      write(0,*)'scat ckpt 2'
!***
      
         AFAC=1._WP
         X=AKK*(3._WP*NAT/(4._WP*PI))**(1._WP/3._WP)
         THETA0=2._WP*PI/(1._WP+X)
         NTHETA=1+NINT((2._WP*(3._WP+X)/ETASCA)*(1._WP+             &
                AFAC/(PI+THETA0))/(1._WP+AFAC*THETA0/(PI+THETA0)**2))

! prepare for first interval in theta

         THETA=0._WP
         THETAU=0._WP

         CSCA=0._WP
         DO K=1,3
            CSCAG(K)=0._WP
            CTRQSCTF(K)=0._WP
         ENDDO
         CSCAG2=0._WP

         IF(CMDTRQ=='DOTORQ')THEN

!***********************************************************************

! Additional quantities required for torque computations
! XYZCM(1-3) = coordinates of grain centroid
! SCRRS1(J,K)= r_j = coordinates of dipole J relative to centroid.
! CXSCR3(J) = r_j dot p_j
! CXSCR4(J,1-3) = p_j cross r_j

#ifdef openmp
!$OMP PARALLEL 
#endif

            DO K=1,3

#ifdef openmp
!$OMP SINGLE
#endif


               XYZCM(K)=0._WP

#ifdef openmp
!$OMP END SINGLE
!$OMP DO                  &
!$OMP    PRIVATE(J)       &
!$OMP    REDUCTION(+:XYZCM)
#endif

               DO J=1,NAT
                  XYZCM(K)=XYZCM(K)+(REAL(IXYZ(J,K),KIND=WP)+X0(K))*DX(K)
               ENDDO

#ifdef openmp
!$OMP END DO
#endif

               XYZCM(K)=XYZCM(K)/REAL(NAT,KIND=WP)
            ENDDO
            DO K=1,3

#ifdef openmp
!$OMP DO          &
!$OMP&   PRIVATE(J)
#endif

               DO J=1,NAT
                  SCRRS1(J,K)=(REAL(IXYZ(J,K),KIND=WP)+X0(K))*DX(K)-XYZCM(K)
               ENDDO

#ifdef openmp
!$OMP END DO
#endif

            ENDDO

#ifdef openmp
!$OMP DO            &
!$OMP&   PRIVATE(J,K)
#endif
            DO J=1,NAT
               CXSCR3(J)=CXZERO
               DO K=1,3
                  CXSCR3(J)=CXSCR3(J)+SCRRS1(J,K)*CXP_TF(J,K)
               ENDDO
            ENDDO

#ifdef openmp
!$OMP END DO
!$OMP DO          &
!$OMP&   PRIVATE(J)
#endif

            DO J=1,NAT
               CXSCR4(J,1)=CXP_TF(J,2)*SCRRS1(J,3)-CXP_TF(J,3)*SCRRS1(J,2)
               CXSCR4(J,2)=CXP_TF(J,3)*SCRRS1(J,1)-CXP_TF(J,1)*SCRRS1(J,3)
               CXSCR4(J,3)=CXP_TF(J,1)*SCRRS1(J,2)-CXP_TF(J,2)*SCRRS1(J,1)
            ENDDO

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif

!***********************************************************************

         ELSE   ! NOTORQ case
!*** diagnostic
!           write(0,*)'scat ckpt 3, nat, mxnat = ',nat,mxnat
!***

#ifdef openmp
!$OMP PARALLEL
#endif

            DO K=1,3

#ifdef openmp
!$OMP DO          &
!$OMP&   PRIVATE(J)
#endif

               DO J=1,NAT
!*** diagnostic
!                  write(0,*)'scat ckpt 3.5 myid=',myid,' J,K=',J,K, &
!                            ' ixyz(j,k)=',ixyz(j,k)
!***
                  
                  SCRRS1(J,K)=(REAL(IXYZ(J,K),KIND=WP)+X0(K))*DX(K)
!*** diagnostic
!                  write(0,*)'scat ckpt 3.6'
!***
               ENDDO

#ifdef openmp
!$OMP END DO
#endif

            ENDDO

#ifdef openmp
!$OMP END PARALLEL
#endif

         ENDIF !--- end IF(CMDTRQ=='DOTORQ') elseif

! Proceed to sum over scattering angles:
!*** diagnostic
!      write(0,*)'scat ckpt 4, NTHETA=',NTHETA
!***
         NAVG=0

           DO ICOSTH=1,NTHETA

!*** diagnostic
!      write(0,*)'scat ckpt 5, ICOSTH=',ICOSTH
!***

! set the scattering angle THETA
! and calculate THETAU = scattering angle for next value of ICOSTH

            THETAL=THETA
            THETA=THETAU
            IF(ICOSTH<NTHETA)THEN

! THETAU will be the *next* value of theta:
! SI = value of S for *next* value of ICOSTH

               SI=ICOSTH*(PI+AFAC*PI/(THETA0+PI))/(NTHETA-1)
               TERM=SI-AFAC-THETA0
               THETAU=0.5_WP*(TERM+SQRT(TERM**2+4._WP*THETA0*SI))

            ELSE
               THETAU=PI
            ENDIF

            DOMEGA=2._WP*PI*(COS(0.5_WP*(THETAL+THETA))-  &
                   COS(0.5_WP*(THETA+THETAU)))

            COSTH=COS(THETA)
            SINTH=SIN(THETA)

            IF(ICOSTH==1.OR.ICOSTH==NTHETA)THEN
               NPHI=1
            ELSE
               NPHI=NINT(4._WP*PI*SINTH/(THETAU-THETAL))
               NPHI=MAX(NPHI,3)
            ENDIF

! evaluate weight factor W

            W=DOMEGA/REAL(NPHI,KIND=WP)

            DPHI=2._WP*PI/REAL(NPHI,KIND=WP)
            NAVG=NAVG+NPHI

            DO IPHI=1,NPHI

! offset phi values by 0.5*DPHI each time

               PHI=DPHI*(REAL(IPHI,KIND=WP)-0.5_WP*MOD(ICOSTH,2))
               SINPHI=SIN(PHI)
               COSPHI=COS(PHI)

! Compute scattered k vector AKS0(1-3)
! Initialize CXES(1-3) = scattered electric field

               DO K=1,3
                  AKS0(K)=COSTH*AK_TF(K)+SINTH*(COSPHI*AKS1(K)+SINPHI*AKS2(K))
                  CXES(K)=CXZERO
               ENDDO

! Compute normalized scattering vector = nhat

               DO K=1,3
                  AKSN(K)=AKS0(K)/AKK
               ENDDO

! Now compute CXES(1-3) = sum_j[p_j-nhat(nhat dot p_j)]exp[-ik_s dot x_j
!      =(r/k*k)*exp(-ikr+iwt)*E(k_s,r,t)
!      where E(k_s,r,t) is complex electric field vector propagating
!      with wave vector k_s, at distance r from origin
! Notation:
!    If CMDTRQ='DOTORQ', then r_j = position relative to centroid XYZCM
!    If CMDTRQ='NOTORQ', then r_j = original lattice coords, and XYZCM=0
! Note that since we are not concerned with overall phase shifts for
! evaluation of forces and torques, we do NOT need to worry about change
! in phase resulting from using two different definitions of r_j when
! computing phase factors exp[-i k_s dot r_j]

#ifdef openmp
!$OMP PARALLEL                  & 
!$OMP&   PRIVATE(CXES_L,CXSCL2_L)
#endif

              CXES_L=CXZERO

#ifdef openmp
!$OMP DO          &
!$OMP& PRIVATE(J,K)
#endif

               DO J=1,NAT

! CXSCR1(J) = nhat dot p(j)

                  CXSCR1(J)=CXZERO
                  DO K=1,3
                     CXSCR1(J)=CXSCR1(J)+AKSN(K)*CXP_TF(J,K)
                  ENDDO

! SCRRS2(J) = nhat dot r(j)

                  SCRRS2(J)=0._WP
                  DO K=1,3
                     SCRRS2(J)=SCRRS2(J)+SCRRS1(J,K)*AKSN(K)
                  ENDDO

! CXSCR2(J) = exp(-i*k_s dot r_j) factor for atom j

                  CXSCR2(J)=EXP(-CXI*AKK*SCRRS2(J))

! CXES(1-3) = sum_j[p_j-nhat(nhat dot p_j)]exp[-ik_s dot x_j]

                  DO K=1,3
!Art                     CXES(K)=CXES(K)+                                 &
!                                 (CXP_TF(J,K)-AKSN(K)*CXSCR1(J))*CXSCR2(J)
                     CXES_L(K)=CXES_L(K)+                              &
                               (CXP_TF(J,K)-AKSN(K)*CXSCR1(J))*CXSCR2(J)
                  ENDDO
               ENDDO

#ifdef openmp
!$OMP END DO
!$OMP ATOMIC
#endif

               CXES(1)=CXES(1)+CXES_L(1)

#ifdef openmp
!$OMP ATOMIC
#endif

               CXES(2)=CXES(2)+CXES_L(2)

#ifdef openmp
!$OMP ATOMIC
#endif

               CXES(3)=CXES(3)+CXES_L(3)

               IF(CMDTRQ=='DOTORQ')THEN

!***********************************************************************
! Additional terms required to compute torque on grain
! Have already computed
! CXES(1-3)=sum_j[p_j - k_s*(k_s dot p_j)/|k|^2]exp[-ik_s dot x_j]
!             where p_j = polarization of dipole j
!                   x_j = position of dipole j rel. to CM
!                   k_s = scattered k vector

! CXBS = -nhat cross CXES  = -(r/k*k)*exp(-ikr+iwt)*B(k_s,r,t)

                  CXBS(1)=AKSN(3)*CXES(2)-AKSN(2)*CXES(3)
                  CXBS(2)=AKSN(1)*CXES(3)-AKSN(3)*CXES(1)
                  CXBS(3)=AKSN(2)*CXES(1)-AKSN(1)*CXES(2)

! CXSCL1 = sum_j p_j dot (r_j cross nhat)*exp(-i*k_s dot r_j)
!        = nhat dot sum_j (p_j cross r_j)*exp(-i*k_s dot r_j)
! recall that have previously evaluated
! CXSCR2(J) = exp(-i*k_s dot r_j) factor for atom j
! CXSCR4(J,K)=p_j cross r_j

#ifdef openmp
!$OMP SINGLE
#endif

                  CXSCL1=CXZERO

#ifdef openmp
!$OMP END SINGLE
#endif


!Art                  DO K=1,3
!Art                     CXTRM(K)=CXZERO
!Art                  ENDDO

                  DO K=1,3

#ifdef openmp
!$OMP SINGLE
#endif

                     CXTRM_K=CXZERO

#ifdef openmp
!$OMP END SINGLE
!$OMP DO                    &     
!$OMP&   PRIVATE(J)         &
!$OMP&   REDUCTION(+:CXTRM_K)
#endif

                     DO J=1,NAT
!Art                        CXTRM(K)=CXTRM(K)+CXSCR2(J)*CXSCR4(J,K)
                        CXTRM_K=CXTRM_K+CXSCR2(J)*CXSCR4(J,K)
                     ENDDO

#ifdef openmp
!$OMP END DO
!$OMP ATOMIC
#endif

!Art                     CXSCL1=CXSCL1+AKSN(K)*CXTRM(K)
                     CXSCL1=CXSCL1+AKSN(K)*CXTRM_K
                  ENDDO

! CXSCL2 = sum_j[(r_j dot p_j) -
!                 (nhat dot p_j)*((nhat dot r_j) + 2i/k_s)]*
!                 exp(-i*k_s dot r_j)
! Have previously computed
! CXSCR1(J) = nhat dot p(j)
! CXSCR2(J) = exp(-i*k_s dot r_j) factor for atom j
! CXSCR3(J) = r_j dot p_j
! SCRRS2(J) = nhat dot r(j)

#ifdef openmp
!$OMP SINGLE
#endif

                  CXSCL2=CXZERO
                  CXSCL3=2._WP*CXI/AKK

#ifdef openmp
!$OMP END SINGLE
!$OMP DO                    &
!&OMP    PRIVATE(J)         &       
!$OMP&   REDUCTION (+:CXSCL2)
#endif

                  DO J=1,NAT
                     CXSCL2=CXSCL2+CXSCR2(J)*                      &
                            (CXSCR3(J)-CXSCR1(J)*(SCRRS2(J)+CXSCL3))
                  ENDDO
#ifdef openmp
!$OMP END DO
#endif

! Note following notational correspondence:
!    present    Draine & Weingartner (1996)
!    -------    ---------------------------
!     CXES   =        V_E
!     CXBS   =        V_B
!    CXSCL1  =        S_B
!    CXSCL2  =        S_E

!***********************************************************************
               ENDIF
#ifdef openmp
!$OMP END PARALLEL
#endif

! Have completed computation of vector CXES(1-3) and terms to calculate
! angular momentum flux for scattering direction (THETA,PHI).
! CTRQSCTF(K) is "torque cross section" due to "scattered radiation"
! resulting in torque in lattice direction k (i.e., in "target
! frame".  After CTRQSCTF has been evaluated, we will compute CTRQSC in
! the "scattering frame", where radiation is incident along the x direct
! and "reference polarization state" is along the y direction.

               ES2=0.0_WP
               DO K=1,3
                  ES2=ES2+REAL(CXES(K)*CONJG(CXES(K)))
               ENDDO
               CSCA=CSCA+W*ES2
               CSCAG(1)=CSCAG(1)+W*COSTH*ES2
               CSCAG(2)=CSCAG(2)+W*SINTH*COSPHI*ES2
               CSCAG(3)=CSCAG(3)+W*SINTH*SINPHI*ES2
               CSCAG2=CSCAG2+W*COSTH**2*ES2
               IF(CMDTRQ=='DOTORQ')THEN
                  DO K=1,3
                     CTRQSCTF(K)=CTRQSCTF(K)-W*REAL(CONJG(CXBS(K))*CXSCL2+ &
                                 CONJG(CXES(K))*CXSCL1)
                  ENDDO
               ENDIF
            ENDDO !--- end DO IPHI=1,NPHI

! Save backscattering cross section:

            IF(ICOSTH==NTHETA)CBKSCA=AK2*AK2*ES2/E02

         ENDDO !--- end DO ICOSTH=1,NTHETA
!***********************************************************************

         CSCA=AK2*AK2*CSCA/E02

! Convert CSCAG,CSCAG2, and CTRQSCTF to cross sections measured in
! lattice units:

         DO K=1,3
            CSCAG(K)=AK2*AK2*CSCAG(K)/E02
         ENDDO
         CSCAG2=AK2*AK2*CSCAG2/E02
         IF(CMDTRQ=='DOTORQ')THEN

!***********************************************************************
! Compute the two contributions to the torque on the target:
! (1) torque due to E and B fields of incident radiation (CTRQIN)
! (2) torque due to E and B fields of scattered radiation (CTRQSC)
! These are first computed in the target frame, then transformed to
! the lab frame (frame where incident radiation propagates along
! the x axis, and "reference pol. state 1" is along the y axis).

            DO K=1,3
               CTRQSCTF(K)=AK2*AK2*CTRQSCTF(K)/E02
            ENDDO

! CTRQSCTF(1-3) is the vector torque cross section/|k| measured in the
! target frame.  We now compute CTRQSC(1-3)=the vector torque cross
! section in the "scattering frame", where the radiation is incident alo
! the x axis, and the "reference polarization state" (at t=0) is
! linearly polarized along the y axis.
! Note that the following computation introduces an additional factor
! of |k|.

            CTRQSC(1)=AK_TF(1)*CTRQSCTF(1)+AK_TF(2)*CTRQSCTF(2)+ &
                      AK_TF(3)*CTRQSCTF(3)
            CTRQSC(2)=AKS1(1)*CTRQSCTF(1)+AKS1(2)*CTRQSCTF(2) + &
                      AKS1(3)*CTRQSCTF(3)
            CTRQSC(3)=AKS2(1)*CTRQSCTF(1)+AKS2(2)*CTRQSCTF(2) + &
                      AKS2(3)*CTRQSCTF(3)

! Now compute vector torque cross section due to incident radiation:

            DO K=1,3
               CTRQIN(K)=0._WP
            ENDDO

! Redefine scratch vector CXSCR1:
! CXSCR1(J)=conjg(p_j) dot E_0j at t=0

#ifdef openmp
!$OMP PARALLEL 
!$OMP DO            & 
!$OMP&   PRIVATE(J,K)
#endif

            DO J=1,NAT
               CXSCR1(J)=CXZERO
               DO K=1,3
                  CXSCR1(J)=CXSCR1(J)+CONJG(CXP_TF(J,K))*CXE_TF(J,K)
               ENDDO
            ENDDO

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif

#ifdef openmp
!$OMP SINGLE
#endif

! 12.04.28 (BTD) I cannot see what is accomplished by this
!                isolated OMP SINGLE section

            CXSCL1=CXZERO
            CXSCL2=CXZERO
            CXSCL3=CXZERO

#ifdef openmp
!$OMP END SINGLE
#endif

! Compute CXSCL1,CXSCL2,CXSCL3=components of the torque due to
! incident E and B fields in the target frame.
! 95.06.22 (BTD): added factor CXI in r cross k term appearing in
!                 evaluations of CXSCL1,CXSCL2,CXSCL2
! 95.07.21 (BTD): corrected sign mistake (+CXI -> -CXI)
! 95.07.24 (BTD): changed def. of CXSCR1 to agree with paper,
!                 changed calc. of CXSCLn to more closely mirror paper,
!                 and changed AKS0 to AK_TF

#ifdef openmp
!$OMP PARALLEL 
!$OMP DO                                 &
!$OMP&   PRIVATE(J)                      &
!$OMP&   REDUCTION(+:CXSCL1,CXSCL2,CXSCL3)
#endif
            DO J=1,NAT
               CXSCL1=CXSCL1+CONJG(CXP_TF(J,2))*CXE_TF(J,3)- &
                         CONJG(CXP_TF(J,3))*CXE_TF(J,2)-     &
                         CXI*(AK_TF(2)*SCRRS1(J,3)-          &
                         AK_TF(3)*SCRRS1(J,2))*CXSCR1(J)
               CXSCL2=CXSCL2+CONJG(CXP_TF(J,3))*CXE_TF(J,1)- &
                         CONJG(CXP_TF(J,1))*CXE_TF(J,3)-     &
                         CXI*(AK_TF(3)*SCRRS1(J,1)-          &
                         AK_TF(1)*SCRRS1(J,3))*CXSCR1(J)
               CXSCL3=CXSCL3+CONJG(CXP_TF(J,1))*CXE_TF(J,2)- &
                         CONJG(CXP_TF(J,2))*CXE_TF(J,1)-     &
                         CXI*(AK_TF(1)*SCRRS1(J,2)-          &
                         AK_TF(2)*SCRRS1(J,1))*CXSCR1(J)
            ENDDO

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif

! Compute CTRQIN(1-3)=2*components of torque due to incident E and B
! fields (in lab frame), expressed as a cross section:

            CTRQIN(1)=AK_TF(1)*REAL(CXSCL1)+AK_TF(2)*REAL(CXSCL2)+ &
                      AK_TF(3)*REAL(CXSCL3)
            CTRQIN(2)=AKS1(1)*REAL(CXSCL1)+AKS1(2)*REAL(CXSCL2)+ &
                      AKS1(3)*REAL(CXSCL3)
            CTRQIN(3)=AKS2(1)*REAL(CXSCL1)+AKS2(2)*REAL(CXSCL2)+ &
                      AKS2(3)*REAL(CXSCL3)
            DO K=1,3
               CTRQIN(K)=4._WP*PI*CTRQIN(K)/E02
            ENDDO
         ENDIF

!***********************************************************************
         WRITE(CMSGNM,FMT='(I8,A)')NAVG,                           &
               ' scattering directions used to calculate <cos>, etc.'
         CALL WRIMSG('SCAT',CMSGNM)

         IF(NDIR<=0)RETURN
      ENDIF ! end IF(JPBC==0)

!============== Following code is used to calculate ====================
!               CXF1L(J=1-NDIR) and CXF2L(J=1-NDIR)
!               for JPBC=0,1,2,3 (i.e., all cases)

!*** diagnostic
!      write(0,*)'in scat, ckpt alpha'
!***
!*** Calculate the 2 matrix elements FM1 for selected directions
! F1L(NDIR) is to the polarization state e1 (E in scattering plane)
! F2L(NDIR) is to the polarization state e2 (E perp. to scatt. plane)

! IMPORTANT NOTE: If the incident radiation is not linearly polarized,
! with zero imaginary component to the polarization vector, then the
! F1L and F2L computed here are NOT the "usual" f_ml, which
! describe scattering from linearly polarized states l to linearly
! polarized states m.
! If one wishes to reconstruct the "usual" f_ml even when computing
! with elliptically polarized light, then one must add additional
! code to subroutine SCAT to reconstruct the f_ml from the F1L and
! F2L computed here.

!*** diagnostic
!      write(0,7431)
!      do j=1,nat
!         write(0,7432)j,cxp(j,1),cxp(j,2),cxp(j,3)
!      enddo
! 7431 format(' j  ---------------- cxp(j,1-3) -------------------------')
! 7432 format(i4,2f10.5,1x,2f10.5,1x,2f10.5)
!*** end diagnostic
!*** diagnostic
!      write(0,*)'scat ckpt 6'
!***
      TERM=AK3/SQRT(E02)
      DO ND=1,NDIR
         CXF1L(ND)=CXZERO
         CXF2L(ND)=CXZERO

#ifdef openmp
!$OMP PARALLEL 
!$OMP SINGLE
#endif

         CXF1L_L=CXZERO
         CXF2L_L=CXZERO

#ifdef openmp
!$OMP END SINGLE
!$OMP DO                              &
!$OMP&   PRIVATE(K)                   &
!$OMP&   REDUCTION(+: CXF1L_L, CXF2L_L)
#endif

         DO J=1,NAT
            CXSCR2(J)=CXZERO
            CXSCR3(J)=CXZERO
            CXSCR1(J)=EXP(-CXI*(                                          &
                      AKS_TF(1,ND)*(REAL(IXYZ(J,1),KIND=WP)+X0(1))*DX(1)+ &
                      AKS_TF(2,ND)*(REAL(IXYZ(J,2),KIND=WP)+X0(2))*DX(2)+ &
                      AKS_TF(3,ND)*(REAL(IXYZ(J,3),KIND=WP)+X0(3))*DX(3)))
            DO K=1,3
               CXSCR2(J)=CXSCR2(J)+CXP_TF(J,K)*EM1_TF(K,ND)
               CXSCR3(J)=CXSCR3(J)+CXP_TF(J,K)*EM2_TF(K,ND)
            ENDDO
            CXF1L_L=CXF1L_L+CXSCR2(J)*CXSCR1(J)
            CXF2L_L=CXF2L_L+CXSCR3(J)*CXSCR1(J)
         ENDDO

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif

!*** diagnostic
!         write(0,9900)nd,aks(1,nd),aks(2,nd),aks(3,nd), &
!                         em1(1,nd),em1(2,nd),em1(3,nd), &
!                         em2(1,nd),em2(2,nd),em2(3,nd)
! 9900 format('nd=',i2,' aksr=',3f10.5,/,           &
!             6x,'em1r=',3f10.5,/,6x,'em2r=',3f10.5)
!
!         do j=1,nat
!            write(0,9910)j,cxscr2(j),cxscr3(j)
! 9910       format(i2,' cxscr2=',2f8.4,' cxscr3=',2f8.4)
!         enddo
!***
!Art         DO J=1,NAT
!Art            CXSCR1(J)=EXP(-CXI*(                                          &
!Art                      AKS_TF(1,ND)*(REAL(IXYZ(J,1),KIND=WP)+X0(1))*DX(1)+ &
!Art                      AKS_TF(2,ND)*(REAL(IXYZ(J,2),KIND=WP)+X0(2))*DX(2)+ &
!Art                      AKS_TF(3,ND)*(REAL(IXYZ(J,3),KIND=WP)+X0(3))*DX(3)))
!Art            CXF1L(ND)=CXF1L(ND)+CXSCR2(J)*CXSCR1(J)
!Art            CXF2L(ND)=CXF2L(ND)+CXSCR3(J)*CXSCR1(J)
!Art            CXF1L_L=CXF1L_L+CXSCR2(J)*CXSCR1(J)
!Art            CXF2L_L=CXF2L_L+CXSCR3(J)*CXSCR1(J)
!Art         ENDDO
!Art!$omp end  do

!   Units: Dipole polarizations are in units of [E]*d**3, where d=lattic
!          spacing, and [E]=dimensions of E field
!          Therefore, up to this point CXF1L, CXF2L are in units
!          of [E]*d**3
!          Now multiply by AK3/SQRT(E02) (units of 1/([E]d**3)) to get
!          dimensionless result.

         CXF1L(ND)=CXF1L_L*TERM
         CXF2L(ND)=CXF2L_L*TERM
!Art!$omp end parallel

! diagnostic
!        write(0,9700)nd,cxf1l(nd),cxf2l(nd)
! 9700   format('nd=',i3,' cxf1L=',1p,2e10.3,' cxf2L=',2e10.3)
! end diagnostic

      ENDDO !--- end DO ND=1,NDIR

!*** diagnostic
!      write(0,*)'returning from SCAT'
!***
      RETURN
    END SUBROUTINE SCAT
