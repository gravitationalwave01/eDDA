    SUBROUTINE TARCYLCAP(A1,A2,A,B,DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

! Scalar arguments:

      CHARACTER :: CDESCR*67
      INTEGER :: MXNAT,NAT
      REAL(WP) :: A,B

! Array arguments:

      INTEGER*2 ::    &
         ICOMP(MXNAT,3)
      INTEGER ::     & 
         IXYZ(MXNAT,3)
      REAL(WP) :: & 
         A1(3),   &
         A2(3),   &
         DX(3),   &
         X0(3)

! Local scalars:

      INTEGER :: IOSHP,JX,JXL,JXU,JY,JZ,NB,NFAC,NLAY,NX2,NY1,NY2,NZ1,NZ2
      REAL(WP) :: ASPR,PI,R2,REFF2,X2,XCM,Y2,Y2M,YCM,Z2,ZCM
!***********************************************************************
! Purpose: to construct cylindrical target with hemispherical caps
! from "atoms", with cylinder axis along x axis.
! cylinder diameter is b,
! length of cylinder itself is a
! length of cylinder with caps is (a+b)

! Input:
!        A = cylinder length (in units of lattice spacing d)
!            *NOT* including caps
!        B = cylinder diameter (in units of lattice spacing d)
!        DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z
!                directions, and d=(dx*dy*dz)**(1/3)=effective lattice
!                spacing
!        IOSHP=device number for "target.out" file
!             =-1 to suppress printing of "target.out"
!        MXNAT = dimensioning information
! Returns:
!        A1(1-3) = unit vector along cylinder axis
!        A2(1-3) = unit vector perpendicular to cylinder axis
!        CDESCR = string describing target (up to 67 chars.)
!        NAT = number of atoms in target
!        IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms in target
!        ICOMP(1-NAT,1-3)=1 (composition identifier)
!        X0(1-3)=location/d) in Target Frame corresponding to dipole with
!                IXYZ=(0,0,0).  This will be treated at the origin of physical
!                coordinates in the TF.
!                Here set to be centroid of capped cylinder

! B.T.Draine, Princeton Univ. Obs., 2006.09.13
! History:
! 06.09.13 (BTD): Created from tarcyl.f
! 06.12.09 (BTD): initialized NAT=0
! 07.09.11 (BTD): Changed IXYZ from INTEGER*2 to INTEGER
! 08.01.13 (BTD): Cosmetic changes
!                 Added X0 to argument list
!                 Added code setting X0 so that centroid will be origin
! 08.06.23 (BTD): Fixed bug discovered by Prescott and Mulvaney
!                 ICOMP was not being initialized for hemispherical caps.
! 08.08.29 (BTD): Modified format 9020
! end history
! Copyright (C) 2006,2007,2008
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

      PI=4._WP*ATAN(1._WP)

! Cylinder axis is along x axis in target frame.

      DO JX=1,3
        A1(JX)=0._WP
        A2(JX)=0._WP
      ENDDO
      A1(1)=1._WP
      A2(2)=1._WP

! A = cylinder length/d
! B = cylinder diameter/d

! Determine centroid XCM,YCM,ZCM

! If A/d_x is near an odd integer, XCM=0
! If A/d_x is near an even integer, XCM=0.5

! If B/d_y is near an odd or even integer, YCM=0 or 0.5
! If B/d_z is near an odd or even integer, ZCM=0 or 0.5

! Determine limits for testing x,y,z values
! In case of disk axis, run from I=1 to I=INT(A+0.5)
! In radial directions, place atoms as follows:
!     If B is close to even number
!          y = j + 0.5
!          with j running from -int((b+.5)/2) to int((b+.5)/2)-1
!     If B is close to odd number
!          x = j
!          with j running from -int((b+.5)/2) to int((b+.5)/2)

      NX2=INT(A/DX(1)+0.5_WP)

      XCM=0.5_WP*REAL(NX2+1,KIND=WP)

      NB=INT(B/DX(2)+0.5_WP)
      IF(2*(NB/2)<NB)THEN
         YCM=0._WP
         NY1=-NB/2
         NY2=NB/2
      ELSE
         YCM=0.5_WP
         NY1=-NB/2+1
         NY2=NB/2
      ENDIF
      NB=INT(B/DX(3)+0.5_WP)
      IF(2*(NB/2)<NB)THEN
         ZCM=0._WP
         NZ1=-NB/2
         NZ2=NB/2
      ELSE
         ZCM=0.5_WP
         NZ1=-NB/2 + 1
         NZ2=NB/2
      ENDIF

      NAT=0
      DO JZ=NZ1,NZ2
         Z2=((JZ-ZCM)*DX(3))**2
         Y2M=((0.5_WP*B)**2-Z2)/DX(2)**2
         DO JY=NY1,NY2
            IF((JY-YCM)**2<Y2M)THEN
               DO JX=1,NX2
                  NAT=NAT + 1
                  IF(NAT>MXNAT)THEN
                     CALL ERRMSG('FATAL','CYLCAP',' NAT.GT.MXNAT')
                  ENDIF
                  IXYZ(NAT,1)=JX
                  IXYZ(NAT,2)=JY
                  IXYZ(NAT,3)=JZ
               ENDDO
            ENDIF
         ENDDO
      ENDDO

! NLAY = number of layers in cylinder
! NFAC = number of atoms in slice

      NLAY=NX2
      NFAC=NAT/NLAY

! Now add hemispherical caps
! XCM = 0.5*REAL(NX2+1)
! XU  = XCM+0.5*(A+B)
! JXU = INT(XU)
! JXL = NX2+1-JXU

      R2=(0.5_WP*B)**2
      JXU=INT(0.5_WP*(NX2+1+A+B))
      JXL=NX2 + 1 - JXU
      DO JZ=NZ1,NZ2
         Z2=((JZ-ZCM)*DX(3))**2
         DO JY=NY1,NY2
            Y2=((JY-YCM)*DX(2))**2
            DO JX=JXL,0
               X2=((JX-XCM+0.5_WP*A)*DX(1))**2
               IF(X2+Y2+Z2<R2) THEN
                  NAT=NAT + 1
                  IF(NAT>MXNAT) THEN
                     CALL ERRMSG('FATAL','CYLCAP',' NAT.GT.MXNAT')
                  ENDIF
                  IXYZ(NAT,1)=JX
                  IXYZ(NAT,2)=JY
                  IXYZ(NAT,3)=JZ
               ENDIF
            ENDDO
            DO JX=NX2+1,JXU
               X2=((JX-XCM-0.5_WP*A)*DX(1))**2
               IF(X2+Y2+Z2<R2)THEN
                  NAT=NAT+1
                  IF(NAT>MXNAT)THEN
                     CALL ERRMSG('FATAL','CYLCAP',' NAT.GT.MXNAT')
                  ENDIF
                  IXYZ(NAT,1)=JX
                  IXYZ(NAT,2)=JY
                  IXYZ(NAT,3)=JZ
               ENDIF
            ENDDO
         ENDDO
      ENDDO

! Set composition

      DO JZ=1,3
         DO JX=1,NAT
            ICOMP(JX,JZ)=1
         ENDDO
      ENDDO

! Set X0 so that origin is at centroid

      DO JY=1,3
         X0(JY)=0.
         DO JX=1,NAT
            X0(JY)=X0(JY)+IXYZ(JX,JY)
         ENDDO
         X0(JY)=-X0(JY)/REAL(NAT)
      ENDDO

! REFF2=effective radius**2 of disk

      REFF2=REAL(NFAC,KIND=WP)/PI

! ASPR=aspect ratio (length/diameter)

      ASPR=0.5_WP*REAL(NLAY,KIND=WP)/SQRT(REFF2)
      WRITE(CDESCR,FMT='(A,I7,A,I4,A,I4,A,F7.4)')' Cyl.prism, NAT=',NAT, &
            ' NFAC=',NFAC,' NLAY=',NLAY,' asp.ratio=',ASPR

      IF(IOSHP>=0)THEN
         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')
         WRITE(IOSHP,FMT=9020)A,B,NAT,A1,A2,DX,X0
         DO JX=1,NAT
            WRITE(IOSHP,FMT=9030)JX,IXYZ(JX,1),IXYZ(JX,2),IXYZ(JX,3), &
                                 ICOMP(JX,1),ICOMP(JX,2),ICOMP(JX,3)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF
      RETURN
9020  FORMAT(' >TARCYLCAP  cylinder with end-caps: Length=',       &
             F8.4,' Diameter=',F8.4,/,                             &
         I10,' = NAT',/,                                           &
         3F10.6,' = A_1 vector',/,                                 &
         3F10.6,' = A_2 vector',/,                                 &
         3F10.6,' = lattice spacings (d_x,d_y,d_z)/d',/,           &
         3F10.5,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ', &
               'for dipole 0 0 0',/,                               &
         '     JA  IX  IY  IZ ICOMP(x,y,z)')

9030  FORMAT(I7,3I4,3I2)
    END SUBROUTINE TARCYLCAP
