    SUBROUTINE ALPHADIAG(AKR,BETADF,PHIDF,THETADF,CALPHA,CXALPH,CXALOF,      &
                         CXALOS,CXE0R,CXEPS,CXSC,CXSCR1,CXZC,CXZW,DX,IBETH,  &
                         IBETH1,ICOMP,IOCC,IPBC,IPHI,IPHI1,IXYZ0,JORTH,MYID, &
                         MXCOMP,MXNAT,MXN3,NAT,NAT0,NAT3,NCOMP,NX,NY,NZ,     &
                         CXRLOC,CSHAPE,SHPAR)
      USE DDPRECISION,ONLY : WP
      USE DDCOMMON_8,ONLY : CMDFFT
      IMPLICIT NONE

!*** Arguments:

      CHARACTER(6) :: CALPHA
      CHARACTER(9) :: CSHAPE
      INTEGER :: IBETH,IBETH1,IPBC,IPHI,IPHI1,JORTH,MXCOMP,MXNAT,MXN3, &
         MYID,NAT,NAT0,NAT3,NCOMP,NX,NY,NZ
      COMPLEX(WP) ::                                                 &
         CXALOF(NAT,3),                                              &
         CXALPH(NAT,3),                                              &
         CXE0R(3),                                                   &
         CXEPS(MXCOMP),                                              &
         CXALOS(NAT,3),                                              &
         CXSC(NAT,3,3),                                              &
         CXSCR1(NAT,3),                                              &
         CXZC(NX+1+IPBC*(NX-1),NY+1+IPBC*(NY-1),NZ+1+IPBC*(NZ-1),6), &
         CXZW(2*NX,2*NY,2*NZ,3),                                     &
         CXRLOC(MXCOMP+1,3,3)
      INTEGER*2 ::     & 
         ICOMP(NAT,3), &
         IOCC(NAT)
      INTEGER :: &
         IXYZ0(NAT0,3)
      REAL(WP) ::       &
         AKR(3),        &
         BETADF(MXNAT), &
         DX(3),         &
         PHIDF(MXNAT),  &
         SHPAR(12),     &
         THETADF(MXNAT)

!*** Local Variables:

      LOGICAL :: INIT
      CHARACTER :: CMSGNM*70
      INTEGER :: IA,IC,IC2,IC3,L
      REAL(WP) :: AK1,AK2,B1,B2,B3,B3L,COSBE,COSPH,COSTH,EMKD,PI, &
                  SINBE,SINPH,SINTH,SUM,LAMBDA
      REAL(WP) :: &
         R(3,3),  &
         RI(3,3)
      COMPLEX(WP) :: CXI,CXRR,CXSUM1,CXSUM2,CXSUM3,CXSUM4,CXSUM5, &
         CXSUM6,CXTERM

!*** Common:

!      CHARACTER(6) :: CMETHD
!-----------------------------------------------------------------------
!      COMMON /M8/CMETHD
!-----------------------------------------------------------------------

      SAVE INIT,CXI,PI
      DATA CXI/(0._WP,1._WP)/
      DATA INIT/.TRUE./

!***********************************************************************
! Given:
!       AKR(1-3)=(kx,ky,kz)*d, where d=effective lattice spacing
!       BETADF(1-NAT)
!       PHIDF(1-NAT)
!       THETADF(1-NAT): orientation angles beta,phi,theta (radians)
!                       specifying orientation of "Dielectric Frame" (DF
!                       relative to the "Target Frame" (TF)

!       CALPHA = 'LATTDR' or 'SCLDR' = polarizability prescription
!       CMETHD = determines 3-d FFT routine used by ESELF
!       CSHAPE = descriptor of target shape
!       CXE0R(1-3) = polarization vector in lattice coordinates
!                    (assumed to be normalized)
!       CXEPS(1-NCOMP)=distinct values of dielectric constant
!       CXSC(1-NAT,1-3,1-3) = complex scratch space for SCLDR calculatio
!       CXSCR1(1-NAT,1-3)   = complex scratch space for SCLDR calculatio
!       CXZC = complex scratch space needed by ESELF
!       CXZW = complex scratch space needed by ESELF
!       DX(1-3)=(dx/d,dy/d,dz/d), where dx,dy,dz=lattice spacings in
!               x,y,z directions, and d=(dx*dy*dz)**(1/3)
!       IBETH= MYID+IBETH1 if first time through combined BETA/THETA
!              orientation loop
!       IBETH1=starting value of IBETH (see above)
!       ICOMP(1-NAT3)=composition identifier for each lattice site and
!                     direction (ICOMP=0 if lattice site is unoccupied)
!                storage scheme:
!                1x,2x,...,NATx,1y,2y,...,NATy,1z,2z,...,NATz
!       IOCC(1-NAT) == 0 if site unoccupied
!                   == 1 if site occupied
!       IPHI = IPHI1 if first time through PHI orientation loop
!       IPHI1= starting value of IPHI (see above)
!       JORTH= 1 if first incident polarization state
!       MXCOMP = dimensioning information
!       MXN3 = dimensioning information
!       MYID = parallel process identifier (=0 if only 1 process)
!       NAT  = number of sites in extended target (incl. vacuum sites)
!       NAT3 = 3*number of sites in extended target (incl. vacuum sites)
!       NCOMP= number of different dielectric tensor elements in target
!       NX   = x-dimension of extended target
!       NY   = y-dimension of extended target
!       NZ   = z-dimension of extended target
!       SHPAR(1-10)=target shape parameters

! Returns:
!       CXALPH(J,1-3)=(alpha_11,alpha_22,alpha_33)/d^3 for dipole J=1-NA
!                     where alpha_ij=complex polarizability tensor.
!       CXALOF(J,1-3)=(alpha_23,alpha_31,alpha_12)/d^3 for dipole J=1-NA

!***
! If CALPHA = LATTDR:
!    Compute dipole polarizability using "Lattice Dispersion Relation"
!    of Draine & Goodman (1993,ApJ,March 10).  It is required that
!    polarizability be such that an infinite lattice of such dipoles
!    reproduce the continuum dispersion relation for radiation
!    propagating with direction and polarization of radiation incident
!    on the DDA target.

! If CALPHA = GKDLDR:
!    Compute dipole polarizability using "Lattice Dispersion Relation"
!    of Gutkowicz-Krusin and Draine (2004)..  It is required that
!    polarizability be such that an infinite lattice of such dipoles
!    reproduce the continuum dispersion relation for radiation
!    propagating with direction and polarization of radiation incident
!    on the DDA target.  This is the recommended option.  It is nearly
!    but not exactly identical to LATTDR.

!***
!    Note: CXALPH = polarizability/d^3
!    In the event that ICOMP=0, then we set CXALPH=1.
!***********************************************************************
! B.T.Draine, Princeton Univ. Obs.
! History:
! 90.09.13 (BTD): Corrected error in DO loop limit.
! 90.11.01 (BTD): Special treatment for "vacuum" sites in order
!                 to allow FFT treatment.
! 90.11.02 (BTD): Set CXALPH=(1.,0.) at vacuum sites.
! 91.04.30 (BTD): Modified to include both O[(kd)^2] correction term
!                 from Goedecke & Obrien (1988) in addition to
!                 radiative reaction correction.
! 91.05.07 (BTD): Modified to allow easy choice among three methods
!                 for computing polarizability.
! 91.05.07 (BTD): Added call to WRIMSG to record which method in use
! 91.05.08 (BTD): Added CALPHA to argument list.
! 91.05.09 (BTD): Added printing of kd and magnitude of correction
!                 terms.
! 91.05.13 (BTD): Experiment with use of approximate directional
!                 average for coefficient of CXEPS in Lattice
!                 Disperision Relation theory
! 91.05.13 (BTD): Corrected error in radiative reaction correction
!                 (error presumably introduced during last week)
! 91.05.14 (BTD): Added separate options LDRXYZ and LDRAVG
! 91.05.15 (BTD): Added option LDR000 to omit CXEPS-dependent
!                 correction term
! 91.05.30 (BTD): Changed LDRXYZ -> LDR100
!                 Added option LDR111
! 91.09.12 (BTD): Corrected error in numerical coefficient for
!                 LDRAVG case
! 91.09.17 (BTD): Eliminate options LDR000,LDRXYZ,LDRAVG
!                 Introduce option LATTDR (LATTice Dispersion Relation)
!                 which includes explicit dependence
!                 on direction and polarization state
!                 remove AK1 from argument list
!                 add AKR, CXE0R to argument list
! 92.05.14 (BTD): Introduce option LDRISO to experiment with
!                 using average directional correction term
!                 rather than using correction for incident direction
!                 and polarization
! 92.05.16 (BTD): Correct error in LDRISO option: average SUM should be
!                 3/15 rather than 4/15
! 97.11.02 (BTD): Add DX to argument list to allow use with noncubic
!                 lattices (additional modifications still required!!)
! 97.12.26 (BTD): Add CXALOF to argument list to allow use with noncubic
!                 lattices (additional modifications still required).
! 97.12.30 (BTD): Add code to evaluate alpha for noncubic case,
!                 using subroutine NONCUBIC to evaluate the lattice
!                 sums R0,R1,R2,R3
! 98.01.14 (BTD): Change sign of radiative reaction correction term.
! 98.01.19 (BTD): Change treatment of R_1,C3, and C4
!                 and output value of C4
!                 temporary modification to allow value of C4 to be
!                 read in from file 'alpha.par'
! 98.01.22 (BTD): Minor change in computation of radiative reaction
!                 correction.
! 98.03.10 (BTD): set C4=-C3, disable reading C4 from 'alpha.par'
! 98.04.27 (BTD): declared IC1,IC2,IC3 as integers.
! 99.02.09 (BTD): experiment with change in expression used to calculate
!                 off-diagonal elements of polarizability tensor
!                 experiment with C4=-C3-R1 to "cancel" the contribution
!                 of R1 to the off-diagonal terms
! 03.01.29 (BTD): remove code to calculate alpha for noncubic lattice
!                 (hold for release in future version)
! 04.02.26 (BTD): add new option 'GKDLDR' to use Gutkowicz-Krusin
!                 and Draine (2004) lattice dispersion relation
! 04.03.31 (BTD): replace WRITE(0,... with CALL WRIMSG(...
! 04.04.01 (BTD): added DUM13,DUM14 to COMMON/M6/ to accomodate use of
!                 COMMON/M6/ in communicating NPY,NPZ to subroutineS
!                 MATVEC and CMATVEC
! 04.05.21 (BTD): cleanup -- eliminated a number of unused local
!                 variables
! 04.09.14 (BTD): modify to allow target to be made of anisotropic
!                 material with arbitrary orientation, with microcrystal
!                 orientation at each lattice site specified by angles
!                 BETADF,PHIDF,THETADF giving orientation of Dielectric
!                 Frame relative to Target Frame
! 05.06.16 (BTD): Changed DUM13,DUM14 in COMMON/M6/ from integer to real
! 06.09.28 (BTD): *** Version 6.2.3 ***
!                 Note: CXZC is not used at all by alpha and could be
!                 deleted from the argument list!
! 06.09.29 (BTD): Added IPBC to argument list
!                 Modified dimensioning of CXZC when IPBC=1
!                 Note: CXZC and CXSC do not appear to be used at
!                 present. Nevertheless, keep CXZC and CXSC in argument
!                 list in case complex scratch space is needed in some
!                 future version of ALPHA
! 07.06.30 (BTD): moved CMDFFT from COMMON/M6/... CMDFFT
!                 to COMMON/M8/CMDFFT
!                 COMMON/M6/ deleted -- no longer needed by alpha
! 07.08.04 (BTD): Version 7.0.3
!                 * replaced COMMON/M8/ with USE MODULE DDCOMMON_8
! 07.09.11 (BTD): Changed IXYZ0 from INTEGER*2 to INTEGER
! 07.10.27 (BTD): Changed SHPAR(6) -> SHPAR(10)
! 08.02.17 (BTD): Changed SHPAR(10) -> SHPAR(12)
! 08.03.11 (BTD): v7.0.5
!                 renamed SUBROUTINE ALPHA -> SUBROUTINE ALPHADIAG
! 08.03.14 (BTD): corrected dimensioning
!                 IXYZ0(MXNAT,3) -> IXYZ0(NAT0,3)
!                 BETADF(MXNAT)  -> BETADF(NAT0)
!                 PHIDF(MXNAT)   -> PHIDF(NAT0)
!                 THETADF(MXNAT) -> THETADF(NAT0)
! 08.07.27 (BTD): corrected dimensioning
!                 BETADF(NAT0)  -> BETADF(MXNAT)
!                 PHIDF(NAT0)   -> PHIDF(MXNAT)
!                 THETADF(NAT0) -> THETADF(MXNAT)
!                 eliminated variable CXALDS from argument list
! end history
! Copyright (C) 1993,1996,1997,1998,1999,2003,2004,2006,2007,2008
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

!*** diagnostic
!      write(0,*)'alphadiag ckpt 1'
!***

      PI=4._WP*ATAN(1._WP)
      AK2=AKR(1)*AKR(1)+AKR(2)*AKR(2)+AKR(3)*AKR(3)
      AK1=SQRT(AK2)
      CXRR=-(AK1*AK2/1.5_WP)*CXI

!*** EMKD = |m|*k_0*d , where |m|=refractive index
!                             k_0 = wave vector in vacuo
!                             d = lattice spacing

      EMKD=SQRT(ABS(CXEPS(1)))*AK1

!***********************************************************************

!*** Lattice dispersion relation (Draine & Goodman 1993)

      IF(CALPHA=='LATTDR')THEN
         WRITE(CMSGNM,FMT='(A,F6.4)')                     &
            'Lattice dispersion relation for |m|k_0d=',EMKD
         CALL WRIMSG('ALPHA ',CMSGNM)

!*** Compute sum (a_j*e_j)^2 , where a_j=unit propagation vector
!                                    e_j=unit polarization vector

         SUM=0.0_WP
         DO L=1,3
            SUM=SUM+(AKR(L)*ABS(CXE0R(L)))**2
         ENDDO
         SUM=SUM/AK2
         B1=-1.8915316_WP*AK2
         B2=(0.1648469_WP-1.7700004_WP*SUM)*AK2
         !NWB 7/3/12
!         PRINT *, 'B1:', B1
!         PRINT *, 'B2:', B2
!         PRINT *, 'AK1:', AK1
!         LAMBDA = 2*3.1416/AK1
!         PRINT *, 'Lambda:', LAMBDA
         DO L=1,3
            DO IA=1,NAT
               IC=ICOMP(IA,L)
               IF(IC>0)THEN

!*** First compute Clausius-Mossotti polarizability:

                  CXTERM=(.75_WP/PI)*(CXEPS(IC)-1._WP)/(CXEPS(IC)+2._WP)

!*** Determine polarizability by requiring that infinite lattice of
!    dipoles have dipersion relation of continuum.

                  CXTERM=CXTERM/(1._WP+CXTERM*(B1+CXEPS(IC)*B2))

!*** Radiative-reaction correction:

                  CXALPH(IA,L)=CXTERM/(1._WP+CXTERM*CXRR)

! set off-diagonal terms to zero

                  CXALOF(IA,L)=0._WP

                  !NWB Diagnostics 7/3/12
!                  IF ( L .EQ. 3 ) THEN
!                     PRINT *, 'CXEPS:', CXEPS(IC)
!                     PRINT *, 'CXALPH(1):', CXALPH(IA,1), CXALPH(IA,2), CXALPH(IA,3)
!                  ENDIF

               ELSEIF(IC==0)THEN

! To avoid divisions by zero, etc., set CXALPH=1 for vacuum sites.

                  CXALPH(IA,L)=1._WP
                  CXALOF(IA,L)=0._WP
               ENDIF
            ENDDO
         ENDDO
         RETURN

!*** Lattice dispersion relation: modified
!    (Gutkowicz-Krusin & Draine 2004)

      ELSEIF(CALPHA=='GKDLDR')THEN

         IF(MYID==0)THEN
            WRITE(CMSGNM,FMT='(A,F6.4)')                               &
                  'GKDLDR Lattice dispersion relation for |m|k_0d=',EMKD
            CALL WRIMSG('ALPHA ',CMSGNM)
         ENDIF

! B1 = (c_1/pi)*ak2 = (-5.9424219/pi)*ak2 = -1.8915316*ak2
! B2 = (c_2/pi)*ak2 = (0.5178819/pi)*ak2 = 0.1648469*ak2
! B3 = -[(3c_2+c_3)/pi]*ak2 = -[(3*0.5178819+4.0069747)/pi]*ak2
!                           = -1.7700004*ak2
! B3L = B3*A(I)**2 where a_i = unit vector in direction of propagation

!*** diagnostic
!         write(0,*)'alphadiag ckpt 2'
!***
         B1=-1.8915316_WP*AK2
         B2=0.1648469_WP*AK2
         B3=-1.7700004_WP*AK2
         DO L=1,3
            B3L=B3*AKR(L)*AKR(L)/AK2
            DO IA=1,NAT
               IC=ICOMP(IA,L)
               IF(IC>0)THEN

!*** First compute Clausius-Mossotti polarizability:

                  CXTERM=(.75_WP/PI)*(CXEPS(IC)-1._WP)/(CXEPS(IC)+2._WP)

!*** Determine polarizability by requiring that infinite lattice of
!    dipoles have dipersion relation of continuum.


                  CXTERM=CXTERM/(1._WP+CXTERM*(B1+CXEPS(IC)*(B2+B3L)))

!*** Radiative-reaction correction:

                  CXALPH(IA,L)=CXTERM/(1._WP+CXTERM*CXRR)

! set off-diagonal terms to zero

                  CXALOF(IA,L)=0._WP

               ELSEIF(IC==0)THEN

! To avoid divisions by zero, etc., set CXALPH=1 for vacuum sites.

                  CXALPH(IA,L)=1._WP
                  CXALOF(IA,L)=0._WP
               ENDIF
            ENDDO
         ENDDO

!*** diagnostic
!         write(0,*)'alphadiag ckpt 3'
!***

! Now enter module to handle possible microcrystal rotation at each
! lattice site.  BETADF,PHIDF,THETADF = 3 rotation angles specifying
! orientation of "Dielectric Frame" (in which dielectric tensor is
! diagonal) relative to Target Frame.

         DO IA=1,NAT
            IC=ICOMP(IA,1)
!*** diagnostic
!            write(0,*)'alphadiag ckpt 3.1, ia=',ia,' ic=',ic
!***
            IF(IC>0)THEN
               IC2=ICOMP(IA,2)
               IC3=ICOMP(IA,3)
               IF(IC/=IC2.OR.IC/=IC3)THEN
!*** diagnostic
!                  write(0,*)'alphadiag ckpt 3.2, ic,ic2,ic3=',ic,ic2,ic3
!***
                  COSTH=COS(THETADF(IA))
                  COSPH=COS(PHIDF(IA))
                  COSBE=COS(BETADF(IA))
                  IF(COSTH*COSBE*COSPH<1._WP)THEN

! if nonzero rotation, recalculate CXALPH and CXALOF

                     SINTH=SIN(THETADF(IA))
                     SINPH=SIN(PHIDF(IA))
                     SINBE=SIN(BETADF(IA))

! Define R = rotation matrix

                     R(1,1)=COSTH
                     R(1,2)=SINTH*COSPH
                     R(1,3)=SINTH*SINPH
                     R(2,1)=-SINTH*COSBE
                     R(2,2)=COSTH*COSBE*COSPH-SINBE*SINPH
                     R(2,3)=COSTH*COSBE*SINPH+SINBE*COSPH
                     R(3,1)=SINTH*SINBE
                     R(3,2)=-COSTH*SINBE*COSPH-COSBE*SINPH
                     R(3,3)=-COSTH*SINBE*SINPH+COSBE*COSPH

! Define RI = inverse of R

                     RI(1,1)=COSTH
                     RI(1,2)=-SINTH*COSBE
                     RI(1,3)=SINTH*SINBE
                     RI(2,1)=SINTH*COSPH
                     RI(2,2)=COSTH*COSBE*COSPH-SINBE*SINPH
                     RI(2,3)=-COSTH*SINBE*COSPH-COSBE*SINPH
                     RI(3,1)=SINTH*SINPH
                     RI(3,2)=COSTH*COSBE*SINPH+SINBE*COSPH
                     RI(3,3)=-COSTH*SINBE*SINPH+COSBE*COSPH

! calculate diagonal elements:

                     CXSUM1=0._WP
                     CXSUM2=0._WP
                     CXSUM3=0._WP
                     CXSUM4=0._WP
                     CXSUM5=0._WP
                     CXSUM6=0._WP
                     DO L=1,3
                        CXSUM1=CXSUM1+R(1,L)*CXALPH(IA,L)*RI(L,1)
                        CXSUM2=CXSUM2+R(2,L)*CXALPH(IA,L)*RI(L,2)
                        CXSUM3=CXSUM3+R(3,L)*CXALPH(IA,L)*RI(L,3)
                        CXSUM4=CXSUM4+R(2,L)*CXALPH(IA,L)*RI(L,3)
                        CXSUM5=CXSUM5+R(3,L)*CXALPH(IA,L)*RI(L,1)
                        CXSUM6=CXSUM6+R(1,L)*CXALPH(IA,L)*RI(L,2)
                     ENDDO
                     CXALPH(IA,1)=CXSUM1
                     CXALPH(IA,2)=CXSUM2
                     CXALPH(IA,3)=CXSUM3
                     CXALOF(IA,1)=CXSUM4
                     CXALOF(IA,2)=CXSUM5
                     CXALOF(IA,3)=CXSUM6
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
!*** diagnostic
!         write(0,*)'alphadiag ckpt 4'
!***
         RETURN
      ELSE
         WRITE(CMSGNM,FMT='(A)') 'Error: invalid option for subroutine ALPHA'
         CALL WRIMSG('ALPHA ',CMSGNM)
         STOP
      ENDIF
    END SUBROUTINE ALPHADIAG
