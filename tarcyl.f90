    SUBROUTINE TARCYL(A1,A2,A,B,ORI,DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

! Scalar arguments:

      CHARACTER :: CDESCR*67
      INTEGER :: MXNAT,NAT
      REAL(WP) :: A,B,ORI

! Array arguments:

      INTEGER*2 :: ICOMP(MXNAT,3)
      INTEGER :: IXYZ(MXNAT,3)
      REAL(WP) :: A1(3),A2(3),DX(3),X0(3)

! Local scalars:

      INTEGER :: IORI,IOSHP,JA,JLO,JHI,JX,JY,JZ,NB,NFAC,NLAY,NX1, &
        NX2,NY1,NY2,NZ1,NZ2
      REAL(WP) :: ASPR,PI,R2,REFF2,RX2,RY2,RZ2,XCM,Y2M,YCM,Z2,ZCM

!***********************************************************************
! Subroutine TARCYL
! Purpose: to construct cylindrical target by populating sites on a
! rectangular lattice.

! Input:
!        A = cylinder length (in units of lattice spacing d)
!        B = cylinder diameter (in units of lattice spacing d)
!        ORI = 1. for cylinder axis in x_TF direction
!                 a1=(1,0,0),a2=(0,1,0)
!              2. for cylinder axis in y_TF direction
!                 a1=(0,1,0),a1=(0,0,1)
!              3. for cylinder axis in z_TF direction
!                 a1=(0,0,1),a2=(1,0,0)
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
!        IXYZ(1-NAT,1-3)=[x-X0(1)]/d,[y-X0(2)]/d,[z-X0(3)]/d for 
!                        dipoles in target
!                        where x=y=z=0 at target centroid
!        X0(1-3) = offset vector defined by above condition
!        ICOMP(1-NAT,1-3)=1 (composition identifier)

! B.T.Draine, Princeton Univ. Obs., 87.04.03
! History:
! 89.12.21 (BTD): Modified for use by ddscat 89.12.21
! 90.12.03 (BTD): Modified for conformity with other target routines
! 91.01.05 (BTD): Changed I4 -> I7 when printing NAT
! 91.04.22 (BTD): Added XV,YV,ZV and A1,A2 to WRITE(IOSHP,FMT=9020)
! 93.03.12 (BTD): Changed CDESCR*(*) -> CDESCR*67
! 95.04.07 (BTD): Corrected error in code which affected target
!                 construction for 4.5<B<5.5, 6.5<B<7.5, ...
!                 In these cases code was construction target with
!                 diameter 4, 6, ... dipoles
!                 This error was discovered by Tao Du
! 96.01.25 (BTD): Modified for conformity with DDSCAT.5a
! 96.01.26 (BTD): Replaced IX,IY,IZ by IXYZ
! 97.12.26 (BTD): Added DX(1-3) to argument list to support noncubic
!                 lattice.
! 98.02.11 (BTD): Make changes to generate target for noncubic lattice.
!                 Take opportunity to simplify and streamline.
! 98.04.27 (BTD): Delete unused variables IN,Y2,ROFF,X,Z,IAXIS
!                 Declare variabels Y2M
!                 Delete IAXIS from WRITE(IOSHP,FMT=9020) statement
! 00.11.02 (BTD): Add ICOMP to argument list, return with ICOMP=1
!                 for occupied sites.
!                 write ICOMP to target.out
! 06.09.15 (BTD): Revised to allow selection of cylinder orientation
!                 in TF via argument ORI
! 06.12.09 (BTD): Initialized NAT=0
! 07.06.21 (BTD): Major rewrite
!                 * added X0 to argument list
!                 * revised to work properly for noncubic lattice
! 07.06.26 (BTD): corrected calculation of NLAY and NFAC for all IORI
! 07.09.10 (BTD): changed definition of X0(1-3) to ensure that target
!                 centroid is at x=y=z=0
!                 where dipole J is at (x,y,z)=[IXYZ(J,1-3)+X0(1-3)]*d
! 07.09.11 (BTD): Changed IXYZ from INTEGER*2 to INTEGER
! 08.08.30 (BTD): Modified format 9020
! 08.08.31 (BTD): added code to set A1,A2
! end history
! Copyright (C) 1993,1995,1996,1997,1998,2000,2006,2007,2008
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

      PI=4._WP*ATAN(1._WP)
      IORI=NINT(ORI)

      DO JX=1,3
         A1(JX)=0._WP
         A2(JX)=0._WP
      ENDDO

! A=cylinder length/d
! B=cylinder diameter/d

      R2=(B/2._WP)**2

      IF(ORI==1)THEN

         A1(1)=1._WP
         A2(2)=1._WP

! cylinder axis in x direction
! determine centroid XCM,YCM,ZCM

! If A/d_x is near an odd integer, XCM=0
!                    even             0.5

! If B/d_y is near an odd integer, YCM=0
!                    even             0.5

! If B/d_z is near an odd integer, ZCM=0
!                    even             0.5

         JX=NINT(A/DX(1))
         IF(MOD(JX,2)==0)THEN
            XCM=0.5_WP
         ELSE
            XCM=0._WP
         ENDIF
         JY=NINT(B/DX(2))
         IF(MOD(JY,2)==0)THEN
            YCM=0.5_WP
         ELSE
            YCM=0._WP
         ENDIF
         JZ=NINT(B/DX(3))
         IF(MOD(JZ,2)==0)THEN
            ZCM=0.5_WP
         ELSE
            ZCM=0._WP
         ENDIF
         NX1=-JX
         NX2=JX+1
         NY1=-JY
         NY2=JY+1
         NZ1=-JZ
         NZ2=JZ+1

         JLO=NX2
         JHI=NX1
         DO JZ=NZ1,NZ2
            RZ2=(REAL(JZ,KIND=WP)*DX(3)-ZCM)**2
            DO JY=NY1,NY2
               RY2=(REAL(JY,KIND=WP)*DX(2)-YCM)**2
               IF(RZ2+RY2<=R2)THEN
                  DO JX=NX1,NX2
                     RX2=2._WP*ABS(REAL(JX,KIND=WP)*DX(1)-XCM)
                     IF(RX2<=A)THEN
                        NAT=NAT+1
                        IF(NAT>MXNAT)THEN
                           CALL ERRMSG('FATAL','TARCYL',' NAT.GT.MXNAT')
                        ENDIF
                        IXYZ(NAT,1)=JX
                        IXYZ(NAT,2)=JY
                        IXYZ(NAT,3)=JZ
                        IF(JX<JLO)JLO=JX
                        IF(JX>JHI)JHI=JX
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         NLAY=JHI-JLO+1

      ELSEIF(ORI==2)THEN

         A1(2)=1._WP
         A2(3)=1._WP

! cylinder axis in y direction

         JY=NINT(A/DX(2))
         IF(MOD(JY,2)==0)THEN
            YCM=0.5_WP
         ELSE
            YCM=0._WP
         ENDIF
         JX=NINT(B/DX(1))
         IF(MOD(JX,2)==0)THEN
            XCM=0.5_WP
         ELSE
            XCM=0._WP
         ENDIF
         JZ=NINT(B/DX(3))
         IF(MOD(JZ,2)==0)THEN
            ZCM=0.5_WP
         ELSE
            ZCM=0._WP
         ENDIF
         NX1=-JX
         NX2=JX+1
         NY1=-JY
         NY2=JY+1
         NZ1=-JZ
         NZ2=JZ

         JLO=NY2
         JHI=NY1
         DO JZ=NZ1,NZ2
            RZ2=(REAL(JZ,KIND=WP)*DX(3)-ZCM)**2
            DO JX=NX1,NX2
               RX2=(REAL(JX,KIND=WP)*DX(1)-XCM)**2
               IF(RZ2+RX2<=R2)THEN
                  DO JY=NY1,NY2
                     RY2=2*ABS(REAL(JY,KIND=WP)*DX(2)-YCM)
                     IF(RY2<=A)THEN
                        NAT=NAT+1
                        IF(NAT>MXNAT)THEN
                           CALL ERRMSG('FATAL','TARCYL',' NAT.GT.MXNAT')
                        ENDIF
                        IXYZ(NAT,1)=JX
                        IXYZ(NAT,2)=JY
                        IXYZ(NAT,3)=JZ
                        IF(JY<JLO)JLO=JY
                        IF(JY>JHI)JHI=JY
                     ENDIF
                  ENDDO
              ENDIF
            ENDDO
         ENDDO
         NLAY=JHI-JLO+1

      ELSEIF(ORI==3)THEN

         A1(3)=1._WP
         A2(1)=1._WP

! cylinder axis in z direction
! determine centroid XCM,YCM,ZCM

! If A/d_z is near an odd integer, ZCM=0
!                    even             0.5

! If B/d_x is near an odd integer, XCM=0
!                    even             0.5

! If B/d_y is near an odd integer, YCM=0
!                    even             0.5

         JZ=NINT(A/DX(3))
          IF(MOD(JZ,2)==0)THEN
            ZCM=0.5_WP
         ELSE
            ZCM=0._WP
         ENDIF
         JX=NINT(B/DX(1))
         IF(MOD(JX,2)==0)THEN
            XCM=0.5_WP
         ELSE
            XCM=0._WP
         ENDIF
         JY=NINT(B/DX(2))
         IF(MOD(JY,2)==0)THEN
            YCM=0.5_WP
         ELSE
            YCM=0._WP
         ENDIF
         NX1=-JX
         NX2=JX+1
         NY1=-JY
         NY2=JY+1
         NZ1=-JZ
         NZ2=JZ+1

         JLO=NZ2
         JHI=NZ1
         DO JX=NX1,NX2
            RX2=(REAL(JX,KIND=WP)*DX(1)-XCM)**2
            DO JY=NY1,NY2
               RY2=(REAL(JY,KIND=WP)*DX(2)-YCM)**2
               IF(RX2+RY2<=R2)THEN
                  DO JZ=NZ1,NZ2
                     RZ2=2*ABS(REAL(JZ,KIND=WP)*DX(3)-ZCM)
                     IF(RZ2<=A)THEN
                        NAT=NAT+1
                        IF(NAT>MXNAT)THEN
                           CALL ERRMSG('FATAL','TARCYL',' NAT.GT.MXNAT')
                        ENDIF
                        IXYZ(NAT,1)=JX
                        IXYZ(NAT,2)=JY
                        IXYZ(NAT,3)=JZ
                        IF(JZ<JLO)JLO=JZ
                        IF(JZ>JHI)JHI=JZ
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         NLAY=JHI-JLO+1
      ELSE
         CALL ERRMSG('FATAL','TARCYL',' INVALID ORI')
      ENDIF

! now redetermine location of centroid

      XCM=0._WP
      YCM=0._WP
      ZCM=0._WP
      DO JA=1,NAT
         XCM=XCM+REAL(IXYZ(JA,1),KIND=WP)
         YCM=YCM+REAL(IXYZ(JA,2),KIND=WP)
         ZCM=ZCM+REAL(IXYZ(JA,3),KIND=WP)
      ENDDO
      X0(1)=-XCM/REAL(NAT,KIND=WP)
      X0(2)=-YCM/REAL(NAT,KIND=WP)
      X0(3)=-ZCM/REAL(NAT,KIND=WP)
! diagnostic
!      write(0,*)'diagnostic: in tarcyl, calculated x0=',x0
! end diagnostic

! Set composition

      DO JZ=1,3
         DO JX=1,NAT
            ICOMP(JX,JZ)=1
         ENDDO
      ENDDO

! NLAY = number of layers in cylinder
! NFAC = number of atoms in slice

      NFAC=NAT/NLAY

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
9020  FORMAT(' >TARCYL  cylindrical target: Length=',F8.4,         &
             ' Diameter=',F8.4,/,                                  &
         I10,' = NAT',/,                                           &
         3F10.6,' = A_1 vector',/,                                 &
         3F10.6,' = A_2 vector',/,                                 &
         3F10.6,' = lattice spacings (d_x,d_y,d_z)/d',/,           &
         3F10.5,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ', &
               'for dipole 0 0 0',/,                               &
         '     JA  IX  IY  IZ ICOMP(x,y,z)')
9030  FORMAT(I7,3I4,3I2)
    END SUBROUTINE TARCYL
