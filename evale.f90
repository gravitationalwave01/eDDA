!*************************Alex Vaschillo and Nicholas Bigelow 2012*************************
!Incorporated code that models a fast electron's interaction with a group of dipoles as in 
!"Optical Excitations in electron microscopy", Rev. Mod. Phys. v. 82 p. 213 equations (4) and (5)
    SUBROUTINE EVALE(CXE00,AKD,DX,X0,IXYZ0,MXNAT,MXN3,NAT,NAT0,NX,NY,NZ,CXE,AEFFA, &
                     WAVEA,MXRAD,MXWAV,Center,CenterX0,CenterX0R, c,velocity,      &
                     e_charge,DielectricConst,XLR,YLR,ZLR,RM)
      !Arguments AEFFA and after added by NWB 3/8/12
      !CenterX0 added by SMC 14.5.13
      !XLR, YLR, ZLR, CenterX0R added by SMC 15.5.13
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

!*** Arguments:
      INTEGER :: MXN3, MXNAT, NAT, NAT0, NX, NY, NZ, MXRAD, MXWAV
      !MXRAD, MXWAV added by NWB 3/8/12
      INTEGER :: IXYZ0(NAT0,3)
      REAL(WP) :: AKD(3), DX(3), X0(3), AEFFA(MXRAD), WAVEA(MXWAV), Center(3), CenterX0(3), &
                  CenterX0R(3), XLR(3), YLR(3), ZLR(3), RM(3,3)
         !AEFFA and after added by NWB 3/8/12
         !CenterX0 added by SMC 14.5.13
         !CenterX0R, XLR, YLR, ZLR added by SMC 15.5.13

! Note: CXE should be dimensioned to CXE(NAT,3) in this routine
!       so that first 3*NAT elements of CXE are employed.
!       XYZ0 should be dimensioned to
      COMPLEX(WP) :: CXE(NAT,3), CXE00(3)

!***  Local variables:
      COMPLEX(WP) :: CXFAC, CXI, CXE_temp(3)
      REAL(WP) :: X, X1, X2, DVEC(3), R, XP, YP, ZP !, CenterX0(3) !Edited out SMC, !DIST added by SMC 15.5.13
      INTEGER :: IA, IX, IY, IZ, M, JJ

!*** Variables added by Alex Vaschillo:
      REAL(WP) :: c, e_charge, EFieldConstant, omega, gamma, k_mag, DS, PI, &
                  BesselArg, DielectricConst, velocity
      REAL(WP) :: Radius, XPe, YPe, ZPe !This serves the exact same purpose as the array R() but is used in a different scope (only in the else statement, see below)
!Added XPe, YPe, ZPe for else statement SMC 15.5.13
      REAL(WP) :: besselk0, besselk1 !Values of Bessel functions K0 and K1
      REAL(WP) :: besseli1, besseli0 !Values of Bessel functions I0 and I0, necessary for K0 and K1 routines
      INTEGER :: i !For indexing do loops

!*** Intrinsic functions and constants:
      INTRINSIC EXP, REAL
      INTRINSIC ATAN2, DCOS, DSIN, AIMAG, DOT_PRODUCT
      PI = 4._WP*ATAN(1._WP)          !Pi

!*** Define center in internal coordinates
!      CenterX0(1) = Center(1) - X0(1) !Center of e-beam in relative coordinates
!      CenterX0(2) = Center(2) - X0(2)
!      CenterX0(3) = Center(3) - X0(3)
!Edited out SMC 14.5.13

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


!*** Compute dipole spacing in meters
      DS = 1E-6_WP * AEFFA(1) * (4._WP * PI / (3._WP * NAT0) )**(1._WP/3._WP)

!*** Compute omega
      omega = 2._WP * PI * c / (WAVEA(1) * 1E-6_WP) !Conversion from microns to meters

!*** Calculate EFieldConstant - the constant that g(r) is multiplied by
      gamma = (1._WP - DielectricConst * (velocity / c) ** 2._WP) ** (-0.5_WP)
      EFieldConstant = 2._WP * e_charge * omega / (velocity ** 2._WP * gamma &
                       * DielectricConst)
      
      PRINT *, 'Relative coordinates of beam:', Center !Spelling correction SMC 8.5.13
      PRINT *, 'Target coordinates of beam:', CenterX0 !Diagnostic added by SMC 14.5.13
      PRINT *, 'Rotated target coordinates of beam:' , CenterX0R !Diagnostic added by SMC 15.5.13
      PRINT *, 'Target lattice offset:', X0 !Diagnostic added by SMC 14.5.13
      PRINT *, 'Electron speed:', velocity
      PRINT *, 'XLR:', XLR !Diagnostic added by SMC 15.5.13
      PRINT *, 'YLR:', YLR
      PRINT *, 'ZLR:', ZLR

      IF (NAT == NAT0) THEN
         !*** Calculate radius and prevent divide by zero errors
         DO i = 1, NAT0
            DVEC(1) = IXYZ0(i, 1) - CenterX0R(1)
            DVEC(2) = IXYZ0(i, 2) - CenterX0R(2)
            DVEC(3) = IXYZ0(i, 3) - CenterX0R(3)
            !R = (IXYZ0(i, 3) - CenterX0(3)) ** 2._WP + &    !Edited out old code
            !       (IXYZ0(i, 2) - CenterX0(2)) ** 2._WP
            !R = SQRT(R) * DS                             !Edited out old code
            !R(1) = (IXYZ0(1, 3) - CenterX0(3)) ** 2._WP + &    
            !       (IXYZ0(1, 2) - CenterX0(2)) ** 2._WP
            !PRINT *, 'R(', i, ') OLD:', R(1)
            !R(1) = SQRT(R(1)) * DS
            XP = DOT_PRODUCT(DVEC,XLR)
            YP = DOT_PRODUCT(DVEC,YLR)
            ZP = DOT_PRODUCT(DVEC,ZLR)
            R = (ZP ** 2._WP) + (YP ** 2._WP)
            R = SQRT(R) * DS
            IF (R .EQ. 0._WP) THEN !If the radius is zero, set to a small, but finite distance
               R = 0.01_WP * DS
               PRINT *, 'WARNING: RADIUS = 0! Re-set to O.01*DS!'
               PRINT *, 'IX, IY, IZ:', IX, IY, IZ
            END IF
     
            !*** Calculate g(r)
            BesselArg = omega * R / (velocity * gamma) !The argument of the Bessel functions
            !CXE(i, 3) = EXP(CXI * omega * DS * (IXYZ0(i, 1) - CenterX0(1)) / velocity) ! This is the prefactor that each component of CXE is multiplied by
            CXE_temp(3) = EXP(CXI * omega * DS * (XP) / velocity)
            
            !*** Calculate electric field components at point i                                                                                                                             
            CXE_temp(1) = EFieldConstant * CXE_temp(3) * (CXI * besselk0(BesselArg) &
                        / gamma)
            !CXE(i, 2) = EFieldConstant * CXE(i, 3) * (-1._WP * besselk1(BesselArg)) * &
            !            DSIN(ATAN2(DS * (DBLE(IXYZ0(i, 2)) - CenterX0(2)) , DS * &
            !            (DBLE(IXYZ0(i, 3)) - CenterX0(3))))
            !CXE(i, 3) = EFieldConstant * CXE(i, 3) * (-1._WP * besselk1(BesselArg)) * &
            !            DCOS(ATAN2(DS * (DBLE(IXYZ0(i, 2)) - CenterX0(2)) , DS * &
            !            (DBLE(IXYZ0(i, 3)) - CenterX0(3))))
            CXE_temp(2) = EFieldConstant * CXE_temp(3) * (-1._WP * besselk1(BesselArg)) * &
                        DSIN(ATAN2((DBLE(YP)) , (DBLE(ZP))))
            CXE_temp(3) = EFieldConstant * CXE_temp(3) * (-1._WP * besselk1(BesselArg)) * &
                        DCOS(ATAN2((DBLE(YP)) , (DBLE(ZP))))
            !Modified 15.5.13 by SMC
            CALL PROD3C(RM,CXE_temp,CXE(i,:))
         END DO
      ELSE
         IA=0 !Index that labels each unique point at which the field is calculated
         DO IZ=1,NZ
            DO IY=1,NY
               DO IX=1,NX
              
                  IA = IA + 1 !Advance IA
                  
                  !*** Calculate Radius and prevent divide by zero errors
                  DVEC(1) = IX - CenterX0R(1)
                  DVEC(2) = IY- CenterX0R(2)
                  DVEC(3) = IZ - CenterX0R(3)
                  !Radius = (IZ - CenterX0(3)) ** 2._WP + &        !Old code
                  !         (IY - CenterX0(2)) ** 2._WP
                  !Radius = SQRT(Radius) * DS
                  !IF (IA == 4000) THEN
                  !Radius = (IZ - CenterX0(3)) ** 2._WP + &        !Diagnostic
                  !         (IY - CenterX0(2)) ** 2._WP
                  !Radius = SQRT(Radius) * DS
                  !PRINT *, 'Rad(1) OLD:', Radius
		  !ENDIF
                  XPe = DOT_PRODUCT(DVEC,XLR)
                  YPe = DOT_PRODUCT(DVEC,YLR)
                  ZPe = DOT_PRODUCT(DVEC,ZLR)
                  Radius = (ZPe) ** 2._WP + &
                           (YPe) ** 2._WP
                  Radius = SQRT(Radius) * DS
                  IF (IA == 1) THEN                       !Diagnostic
                  !PRINT *, 'Rad(1) NEW:', Radius
                  PRINT *, 'XPe:', XPe
                  PRINT *, 'Check:', DVEC(1)
                  PRINT *, 'YPe:', YPe
                  PRINT *, 'Check:', DVEC(2)
                  PRINT *, 'ZPe:', ZPe
                  PRINT *, 'Check:', DVEC(3)
                  ENDIF
                  IF (Radius .EQ. 0._WP) THEN !If the radius is zero, set to a small, but finite distance
                     Radius = 0.01_WP * DS
                     PRINT *, 'WARNING: RADIUS = 0! Re-set to O.01*DS!'
                     PRINT *, 'IX, IY, IZ:', IX, IY, IZ
                  END IF

                  !*** Calculate g(r)
                  BesselArg = omega * Radius / (velocity * gamma) !The argument of the Bessel functions
                  !CXE(IA, 3) = EXP(CXI * omega * DS * (DBLE(IX) - CenterX0(1)) / velocity) ! This is the prefactor that each component of CXE is multiplied by
                  CXE_temp(3) = EXP(CXI * omega * DS * (XPe) / velocity)
                                  
                  !*** Calculate electric field components at point IA
                  CXE_temp(1) = EFieldConstant * CXE_temp(3) * (CXI * besselk0(BesselArg) &
                               / gamma)
                  !CXE(IA, 2) = EFieldConstant * CXE(IA, 3) * (-1._WP * besselk1(BesselArg)) * &
                  !             DSIN(ATAN2(DS * (DBLE(IY) - CenterX0(2)) , DS * (DBLE(IZ) - &
                  !             CenterX0(3))))
                  !CXE(IA, 3) = EFieldConstant * CXE(IA, 3) * (-1._WP * besselk1(BesselArg)) * &
                  !             DCOS(ATAN2(DS * (DBLE(IY) - CenterX0(2)) , DS * (DBLE(IZ) - &
                  !             CenterX0(3))))
                  CXE_temp(2) = EFieldConstant * CXE_temp(3) * (-1._WP * besselk1(BesselArg)) * &
                               DSIN(ATAN2((DBLE(YPe)) , (DBLE((ZPe)))))
                  CXE_temp(3) = EFieldConstant * CXE_temp(3) * (-1._WP * besselk1(BesselArg)) * &
                               DCOS(ATAN2((DBLE(YPe)) , (DBLE((ZPe)))))
                  
                  CALL PROD3C(RM,CXE_temp,CXE(IA,:))
               END DO
            END DO
         END DO
         
         PRINT *, "IA is: ", IA
      ENDIF

      RETURN
    END SUBROUTINE EVALE
