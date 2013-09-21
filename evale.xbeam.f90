!*****************************Alex Vaschillo 2/9/2012*****************************
!Incorporated code that models a fast electron's interaction with a group of dipoles as in 
!"Optical Excitations in electron microscopy", Rev. Mod. Phys. v. 82 p. 213 equations (4) and (5)
    SUBROUTINE EVALE(CXE00,AKD,DX,X0,IXYZ0,MXNAT,MXN3,NAT,NAT0,NX,NY,NZ,CXE,AEFFA, &
                     WAVEA,MXRAD,MXWAV,Center)
      !Arguments AEFFA and after added by NWB 3/8/12
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

!*** Arguments:

      INTEGER :: MXN3,MXNAT,NAT,NAT0,NX,NY,NZ,MXRAD,MXWAV
      !MXRAD, MXWAV added by NWB 3/8/12
      INTEGER ::     &
         IXYZ0(NAT0,3)
      REAL(WP) :: &
         AKD(3),  &
         DX(3),   &
         X0(3),   &
         AEFFA(MXRAD), &
         WAVEA(MXWAV), &
         Center(3)
         !AEFFA and after added by NWB 3/8/12

! Note: CXE should be dimensioned to CXE(NAT,3) in this routine
!       so that first 3*NAT elements of CXE are employed.
!       XYZ0 should be dimensioned to

      COMPLEX(WP) :: & 
         CXE(NAT,3), &
         CXE00(3)

!***  Local variables:

      COMPLEX(WP) :: CXFAC,CXI
      REAL(WP) :: X,X1,X2,CenterX0(3)
      INTEGER :: IA,IX,IY,IZ,M, JJ

!*** Variables added by Alex Vaschillo:
      REAL(WP) :: c, e_charge, EFieldConstant, omega, gamma, k_mag, DS, PI
      REAL(WP) :: BesselArg
      REAL(WP) :: DielectricConst, velocity
      REAL, DIMENSION(:), ALLOCATABLE :: R 
      REAL(WP) :: Radius !This serves the exact same purpose as the array R() but is used in a different scope (only in the else statement, see below)
      REAL(WP) :: besselk0, besselk1 !Values of bessel functions K0 and K1 (these will be evaluated at BesselArg)
      Real(WP) :: besseli1, besseli0

      Integer :: i !For indexing do loops

!*** Intrinsic functions:

      INTRINSIC EXP, REAL
      INTRINSIC ATAN2, DCOS, DSIN, AIMAG

!*** Fast electron constants

!      c = 3.E8_WP              !Speed of light, m/s
      c = 3.E10_WP             !Speed of light, cm/s (cgs)
!      e_charge = 1.60217646E-19_WP !Charge of electron, Coulombs
      e_charge = 4.8032068E-8_WP !Charge of electron, esu (cgs)
      velocity = 0.5_WP * c   !Velocity of electron
      DielectricConst = 1._WP !Dielectric constant
      PI=4._WP*ATAN(1._WP)

      CenterX0(1) = Center(1) - X0(1) !Center of e-beam
      CenterX0(2) = Center(2) - X0(2)
      CenterX0(3) = Center(3) - X0(3)

!***********************************************************************
! subroutine EVALE

! Given:   CXE00(1-3)=Incident E field at origin (complex) at t=0
!          AKD(1-3)=(kx,ky,kz)*d for incident wave (d=effective
!                    lattice spacing)
!          DX(1-3)=(dx/d,dy/d,dz/d) for lattice (dx,dy,dz=lattice
!                   spacings in x,y,z directions, d=(dx*dy*dz)**(1/3)
!          X0(1-3)=(x,y,z)location/(d*DX(1-3)) in TF of lattice site
!                  with IX=0,IY=0,IZ=0
!          IXYZ0(1-NAT0,3)=[x-x0(1)]/dx, [y-x0(2)]/dy, [z-x0(3)]/dz
!                  for each of NAT0 physical dipoles
!          MXNAT,MXN3=dimensioning information
!          NAT0=number of dipoles in physical target
!          NAT=number of locations at which to calculate CXE

! Returns: CXE(1-NAT,3)=incident E field at NAT locations at t=0

! B.T.Draine, Princeton Univ. Obs., 88.05.09
! History:
! 90.11.06 (BTD): Modified to pass array dimension.
! 90.11.29 (BTD): Modified to allow option for either
!                   physical locations only (NAT=NAT0), or
!                   extended dipole array (NAT>NAT0)
! 90.11.30 (BTD): Corrected error for case NAT>NAT0
! 90.11.30 (BTD): Corrected another error for case NAT>NAT0
! 90.12.03 (BTD): Change ordering of XYZ0 and CXE
! 90.12.05 (BTD): Corrected error in dimensioning of CXE
! 90.12.10 (BTD): Remove XYZ0, replace with IXYZ0
! 97.11.02 (BTD): Add DX to argument list to allow use with
!                 noncubic lattices.
! 07.06.20 (BTD): Add X0 to the argument list to specify location
!                 in TF corresponding to IX=0,IY=0,IZ=0
! 07.09.11 (BTD): Changed IXYZ0 from INTEGER*2 to INTEGER
! 08.03.14 (BTD): v7.05
!                 corrected dimensioning
!                 IXYZ0(MXNAT,3) -> IXYZO(NAT0,3)
! Copyright (C) 1993,1997,2007 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

      CXI=(0._WP,1._WP)

! Evaluate electric field vector at each dipole location.

! If NAT=NAT0, then evaluate E only at occupied sites.
! If NAT>NAT0, then evaluate E at all sites.


!*** Allocate memory to arrays

      ALLOCATE( R(NAT) )

!*** Compute dipole spacing

      DS = AEFFA(1) * (4._WP * PI / (3._WP * NAT0) )**(1._WP/3._WP)
      DS = DS * 1E-6_WP
!      PRINT *, 'Dipole spacing (m):', DS

!*** Compute omega

      omega = 2._WP * PI * c / (WAVEA(1) * 1E-6_WP)
!      PRINT *, 'omega:', omega

!*** Calculate EFieldConstant - the constant that g(r) is multiplied by

      gamma = (1._WP - DielectricConst * (velocity / c)**2._WP)**(-0.5_WP)
      EFieldConstant = 2._WP * e_charge * omega / (velocity**2._WP * gamma &
                       * DielectricConst)
      
      PRINT *, 'X0:', X0
      PRINT *, 'Relavitve coordinates of beam:', Center
      PRINT *, 'Internal coordinates of beam:', Center(1) - X0(1), &
           Center(2) - X0(2), Center(3) - X0(3)
!      PRINT *, 'gamma:', gamma
      PRINT *, 'electron charge:', e_charge
!      PRINT *, 'electron speed:', velocity
!      PRINT *, 'Dielectric constant:', DielectricConst
!      PRINT *, 'epsilon:', DielectricConst
!      PRINT *, 'DX:', DX(1), DX(2), DX(3)
!      PRINT *, 'Aeff:', AEFFA(1)
      PRINT *, 'Wavelength:', WAVEA(1)

      IF (NAT == NAT0) THEN
         
         PRINT *, '*** IF STATEMENT: EField for dipole points ***'
!         PRINT *, 'MXNAT:', MXNAT
!         PRINT *, 'MXN3:', MXN3
!         PRINT *, 'NAT:', NAT, " NAT0: ", NAT0
!         PRINT *,         

         !*** Calculate radius
        
         DO i = 1, NAT0
            R(i) = (IXYZ0(i, 1) - CenterX0(1))**2._WP + &
                   (IXYZ0(i, 2) - CenterX0(2))**2._WP
            R(i) = SQRT(R(i)) * DS * 100._WP !(cgs)
            IF (R(i) .EQ. 0._WP) THEN
               R(i) = 0.01_WP*DS
            END IF
         END DO
     
         !*** Calculate g(r)
         
         DO i = 1, NAT0

            BesselArg = omega * R(i) / (velocity * gamma) !The argument of the Bessel functions
            !PRINT *, 'The Bessel argument:', BesselArg

            CXE(i, 3) = EXP(CXI * omega * DS * (IXYZ0(i, 3) - CenterX0(3)*100._WP) / velocity) ! This is the prefactor that each component of CXE is multiplied by (cgs)
            !PRINT *, 'The prefactor:', CXE(i, 1)
            !PRINT *, 'bessk0:', besselk0(BesselArg)
            !PRINT *, 'bessk1:', besselk1(BesselArg)

            !*** Calculate electric field components at point i                                                                                                                             

            CXE(i, 1) = EFieldConstant * CXE(i, 3) * (CXI * besselk0(BesselArg) &
                        / gamma)
            CXE(i, 2) = EFieldConstant * CXE(i, 3) * (-1._WP * besselk1(BesselArg)) * &
                        DSIN(ATAN2(DS * (IXYZ0(i, 2) - CenterX0(2)) , DS * (IXYZ0(i, 1) &
                        - CenterX0(1))))
            CXE(i, 3) = EFieldConstant * CXE(i, 3) * (-1._WP * besselk1(BesselArg)) * &
                        DCOS(ATAN2(DS * (IXYZ0(i, 2) - CenterX0(2)) , DS * (IXYZ0(i, 1) &
                        - CenterX0(1))))

            !PRINT *, 'CXE of X:', CXE(i, 1)
            !PRINT *, 'CXE of Y:', CXE(i, 2)
            !PRINT *, 'CXE of Z:', CXE(i, 3)
            !PRINT *,

         END DO
     
         !*** Print out the results to a file!
         
!         OPEN(UNIT = 10000, FILE = 'EFieldOutputAlex.txt', STATUS = 'replace')
!         i = 0
!         DO i = 0, NAT0
!
!            WRITE(10000, *), (IXYZ0(i, 1) * DX(1)), (IXYZ0(i, 2) * DX(2)), &
!                             (IXYZ0(i, 3) * DX(3)), REAL(CXE(i, 1)), &
!                             REAL(CXE(i, 2)), REAL(CXE(i, 3))
!         END DO
     
         !Below is the original code that calculated the EField due to an incident wave:
!         DO IA=1,NAT0
!            X=0._WP
!            DO M=1,3
!               X=X+AKD(M)*DX(M)*(REAL(IXYZ0(IA,M),KIND=WP)+X0(M))
!            ENDDO
!            CXFAC=EXP(CXI*X)
!            DO M=1,3
!               CXE(IA,M)=CXE00(M)*CXFAC
!            ENDDO
!         ENDDO

      ELSE

         PRINT *, '*** ELSE STATEMENT: EField for non-dipole points ***'
!         PRINT *, 'MXNAT:', MXNAT
!         PRINT *, 'MXN3:', MXN3
!         PRINT *, 'NAT:', NAT, " NAT0: ", NAT0
!         PRINT *, 
     
         IA=0 !Index that labels each unique point at which the field is calculated
         DO IZ=1,NZ
            DO IY=1,NY
               DO IX=1,NX
              
                  IA = IA+1
                  
                  !*** Calculate Radius

!                  WRITE(10000, *) IX, IY, IZ
                  
                  Radius = (IX - CenterX0(1))**2._WP + &
                           (IY - CenterX0(2))**2._WP
                  Radius = SQRT(Radius) * DS * 100._WP !(cgs)
                  IF (Radius .EQ. 0._WP) THEN
                     Radius = 0.01_WP*DS
                     PRINT *, 'WARNING: RADIUS = 0! Re-set to O.01*DS!'
                     PRINT *, 'IX,IY,IZ:',IX,IY,IZ
                  END IF
!                  PRINT *, 'radius:', Radius
!                  PRINT *, 'The x coordinate:', IX * DS                  
!                  PRINT *, 'The y coordinate:', IY * DS
!                  PRINT *, 'The z coordinate:', IZ * DS
              
                  !*** Calculate g(r)
                  
                  BesselArg = omega * Radius / (velocity * gamma) !The argument of the Bessel functions
!                  PRINT *, 'The Bessel argument:', BesselArg
              
                  CXE(IA, 3) = EXP(CXI * omega * DS * (IZ - CenterX0(3) * 100._WP) / velocity) ! This is the prefactor that each component of CXE is multiplied by (cgs)
!                  PRINT *, 'The prefactor:', CXE(IA, 1)
!                  PRINT *, 'bessk0:', besselk0(BesselArg)
!                  PRINT *, 'bessk1:', besselk1(BesselArg)
              
                  !*** Calculate electric field components at point IA
                  
                  CXE(IA, 1) = EFieldConstant * CXE(IA, 3) * (CXI * besselk0(BesselArg) &
                               / gamma)
                  CXE(IA, 2) = EFieldConstant * CXE(IA, 3) * (-1._WP * besselk1(BesselArg)) * &
                               DSIN(ATAN2(DS * (IY - CenterX0(2)) , DS * (IX - CenterX0(1))))
                  CXE(IA, 3) = EFieldConstant * CXE(IA, 3) * (-1._WP * besselk1(BesselArg)) * &
                               DCOS(ATAN2(DS * (IY - CenterX0(2)) , DS * (IX - CenterX0(1))))
!                  PRINT *, 'CXE of X:', CXE(IA, 1)
!                  PRINT *, 'CXE of Y:', CXE(IA, 2)
!                  PRINT *, 'CXE of Z:', CXE(IA, 3) 
!                  PRINT *,
         
                  !Check for NAN components and print out for error checking
                  DO JJ = 1,3
                     IF(ISNAN(REAL(CXE(IA,JJ)))) THEN

                        PRINT *, 'radius:', Radius
                        PRINT *, 'The x coordinate:', IX * DS
                        PRINT *, 'The y coordinate:', IY * DS
                        PRINT *, 'The z coordinate:', IZ * DS
                        PRINT *, 'The prefactor:', CXE(IA, 1)
                        PRINT *, 'bessk0:', besselk0(BesselArg)
                        PRINT *, 'bessk1:', besselk1(BesselArg)
                        PRINT *, 'CXE of X:', CXE(IA, 1)
                        PRINT *, 'CXE of Y:', CXE(IA, 2)
                        PRINT *, 'CXE of Z:', CXE(IA, 3)
                        PRINT *,
                        
                     ENDIF
                  ENDDO
     
               END DO
            END DO
         END DO
         
         !*** Print out the results to a file
!         OPEN (unit = 10000, file = 'EFieldOutputAlex.txt', status = 'replace')
!         IA = 0
!         DO IZ = 0,NZ
!            DO IY = 0,NY
!               DO IX = 0,NX
!                
!                  IA = IA + 1
!                  
!                  WRITE (10000, *), (IX * DX(1)), (IY * DX(2)), &
!                                    (IZ * DX(3)), REAL(CXE(IA, 1)), AIMAG(CXE(IA, 1)), & 
!                                    REAL(CXE(IA, 2)), AIMAG(CXE(IA, 2)), REAL(CXE(IA, 3)), &
!                                    AIMAG(CXE(IA, 3))
!               END DO
!            END DO
!         END DO
         PRINT *, "IA is: ", IA
     
!Below is the old code that calculates the EField at non-dipole points due to an incident wave
!         DO IZ=1,NZ
!            X1=AKD(3)*DX(3)*(REAL(IZ,KIND=WP)+X0(3))
!            DO IY=1,NY
!               X2=X1+AKD(2)*DX(2)*(REAL(IY,KIND=WP)+X0(2))
!               DO IX=1,NX
!                  IA=IA+1
!                  X=X2+AKD(1)*DX(1)*(REAL(IX,KIND=WP)+X0(1))
!                  CXFAC=EXP(CXI*X)
!                  DO M=1,3
!                     CXE(IA,M)=CXE00(M)*CXFAC
!                  ENDDO
!               ENDDO
!            ENDDO
!         ENDDO

      ENDIF


!      CLOSE (unit = 10000)
      DEALLOCATE( R )
      RETURN
    END SUBROUTINE EVALE
