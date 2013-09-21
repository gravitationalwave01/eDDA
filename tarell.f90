    SUBROUTINE TARELL(A1,A2,AX,AY,AZ,DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

!** Arguments:

      CHARACTER :: CDESCR*67
      INTEGER :: IOSHP, MXNAT, NAT
      INTEGER*2 :: ICOMP(MXNAT,3)
      INTEGER ::     &
         IXYZ(MXNAT,3)
      REAL(WP) :: AX, AY, AZ
      REAL(WP) :: &
         A1(3),   &
         A2(3),   &
         DX(3),   &
         X0(3)

!** Local variables:

      INTEGER :: JX,JY,JZ,LMX1,LMX2,LMY1,LMY2,LMZ1,LMZ2
      REAL(WP) :: AX2,AY2,AZ2,R,RYZ2,RZ2,X,XOFF,Y,YOFF,Z,ZOFF

!***********************************************************************
! Routine to construct ellipsoid from "atoms"
! Input:
!        AX=(x-length)/d    (d=lattice spacing)
!        AY=(y-length)/d
!        AZ=(z-length)/d
!        DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z
!                directions, and d=(dx*dy*dz)**(1/3)=effective lattice
!                spacing
!        IOSHP=device number for "target.out" file
!             =-1 to suppress printing of "target.out"
!        MXNAT=dimensioning information (max number of atoms)
! Output:
!        A1(1-3)=(1,0,0)=unit vector defining target axis 1 in Target Fr
!        A2(1-3)=(0,1,0)=unit vector defining target axis 2 in Target Fr
!        NAT=number of atoms in target
!        IXYZ(1-NAT,1-3)=(x-x0(1))/d,(y-x0(2))/d,(z-x0(3))/d
!                        for atoms of target
!        CDESCR=description of target (up to 67 characters)
!        ICOMP(1-NAT,1-3)=1 (composition identifier)
!        X0(1-3)=(location/d) in Target Frame corresponding to dipole with
!                IXYZ=(0,0,0).  This will be treated at the origin of physical
!                coordinates in the TF.
!                Here origin is set to be centroid of ellipsoid.

! B.T.Draine, Princeton Univ. Obs.
! History:
! 91.01.05 (BTD): Changed I4 -> I7 when printing NAT
! 91.04.22 (BTD): Added AX,AY,AZ and A1,A2 to WRITE(IOSHP,FMT=9020)
! 91.05.01 (BTD): TEMPORARY hack to be able to reproduce all arrays
!                 used by Draine 88 (in which XOFF=YOFF=ZOFF=0.5 always)
! 92.10.26 (BTD): Added code to set string CDESCR
!                 Added one call to WRIMSG
! 93.01.01 (BTD): Changed method for determining target axes A1, A2
!                 [previous criteria -- A1 = longest dimension,
!                 A2 = intermediate dimension -- created confusion]
!                 new rule: A1=(1,0,0), A2=(0,1,0)
! 93.03.12 (BTD): Changed CDESCR*(*) -> CDESCR*67
! 96.01.26 (BTD): Replaced IX,IY,IZ by IXYZ
! 97.11.02 (BTD): Added DX to argument list, and modified to allow use
!                 of noncubic lattice
! 98.02.11 (BTD): Modified to write DX to file "target.out"
! 00.10.19 (BTD): Corrected to reduce length of string written to
!                 CDESCR for ellipsoidal grain
! 00.11.02 (BTD): Add ICOMP to argument list, return ICOMP=1
!                 Write ICOMP to target.out
! 07.09.11 (BTD): Change IXYZ from INTEGER*2 to INTEGER
! 08.01.13 (BTD): Cosmetic changes
!                 Added X0 to argument list
!                 Add code specifying (x,y,z)_TF = 0 to be at centroid
!                 of ellipsoid
! 08.02.07 (BTD): Corrected typo in evaluation of X0
! 08.08.29 (BTD): Modified format 9020
! end history

! Copyright (C) 1993,1996,1997,1998,2000,2007,2008 
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

! Routine to construct pseudo-ellipsoidal target aray.
! With occupied array sites contained within ellipsoidal surface
! defined by (X/AX*d)**2+(Y/AY*d)**2+(Z/AZ*d)**2=0.25
! Ideal volume V=(pi/6)*AX*AY*AZ*d**3
! where d = effective lattice spacing

! Dipoles are located on lattice at sites
! (x,y,z)=(I+XOFF,J+YOFF,Z+KOFF), I,J,K=integers
!                                 XOFF,YOFF,ZOFF=constants

! For sphere: call with AX=AY=AZ
! For spheroid: call with AX=AY (or AY=AZ or AX=AZ)
! B.T.Draine, Princeton Univ. Obs., 88.08.12

! Criterion for choosing XOFF:
! If AX is close to an even integer, take XOFF=1/2
! If AX is close to an odd integer, take XOFF=0
!-----------------------------------------------------------------------
!*** diagnostic
!      write(0,*)'entered tarell,ax,ay,az=',ax,ay,az
!***
      JX=INT(AX/DX(1)+.5_WP)
      IF(2*(JX/2)==JX)THEN
         XOFF=0.5_WP
      ELSE
         XOFF=0._WP
      ENDIF

! Similar criterion for YOFF:

      JY=INT(AY/DX(2)+.5_WP)
      IF(2*(JY/2)==JY)THEN
         YOFF=0.5_WP
      ELSE
         YOFF=0._WP
      ENDIF

! Similar criterion for ZOFF:

      JZ=INT(AZ/DX(3)+.5_WP)
      IF(2*(JZ/2)==JZ)THEN
         ZOFF=0.5_WP
      ELSE
         ZOFF=0._WP
      ENDIF

      LMX1=-INT(.5_WP*AX/DX(1)+.5_WP)
      LMX2=INT(.5_WP*AX/DX(1)-.25_WP)
      LMY1=-INT(.5_WP*AY/DX(2)+.5_WP)
      LMY2=INT(.5_WP*AY/DX(2)-.25_WP)
      LMZ1=-INT(.5_WP*AZ/DX(3)+.5_WP)
      LMZ2=INT(.5_WP*AZ/DX(3)-.25_WP)
      AX2=AX*AX
      AY2=AY*AY
      AZ2=AZ*AZ
      NAT=0

! Specify target axes A1 and A2
! Convention: A1=(1,0,0) in target frame
!             A2=(0,1,0) in target frame

      DO JX=1,3
         A1(JX)=0._WP
         A2(JX)=0._WP
      ENDDO
      A1(1)=1._WP
      A2(2)=1._WP

!*** diagnostic
!      write(0,*)'tarell, ckpt 2'
!***

! Determine list of occupied sites.

      DO JZ=LMZ1, LMZ2
         Z=(REAL(JZ,KIND=WP)+ZOFF)*DX(3)
         RZ2=Z*Z/AZ2
         IF(RZ2<0.25_WP)THEN
            DO JY=LMY1, LMY2
               Y=(REAL(JY,KIND=WP)+YOFF)*DX(2)
               RYZ2=RZ2 + Y*Y/AY2
               IF(RYZ2<0.25_WP)THEN
                  DO JX=LMX1, LMX2
                     X=(REAL(JX,KIND=WP)+XOFF)*DX(1)
                     R=RYZ2+X*X/AX2
                     IF(R<0.25_WP)THEN
! Site is occupied:
                        NAT=NAT+1
                        IXYZ(NAT,1)=JX
                        IXYZ(NAT,2)=JY
                        IXYZ(NAT,3)=JZ
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO

!*** diagnostic
!      write(0,*)'tarell, ckpt 3'
!***
      IF(NAT>MXNAT)THEN
         CALL ERRMSG('FATAL','TARELL',' NAT.GT.MXNAT ')
      ENDIF

! Set composition

      DO JX=1,3
         DO JZ=1,NAT
            ICOMP(JZ,JX)=1
         ENDDO
      ENDDO

!*** diagnostic
!      write(0,*)'tarell, ckpt 4'
!***
! Set X0 so that origin is at centroid

      DO JY=1,3
         X0(JY)=0._WP
         DO JX=1,NAT
            X0(JY)=X0(JY)+REAL(IXYZ(JX,JY))
         ENDDO
         X0(JY)=-X0(JY)/REAL(NAT)
      ENDDO
!*** diagnostic
!      write(0,*)'tarell, ckpt 5'
!***
!***********************************************************************
! Write target description into string CDESCR
      IF(AX==AY.AND.AX==AZ)THEN
         WRITE(CDESCR,FMT='(A,I7,A,F6.3,F6.3,F6.3,A)')' Sphere,',NAT, &
            ' dipoles,',DX(1),DX(2),DX(3),'=x,y,z lattice spacing'
      ELSE
         WRITE(CDESCR,FMT='(A,I7,A,F6.3,F6.3,F6.3,A)')' Ellipsoid,',NAT, &
            ' dipoles,',DX(1),DX(2),DX(3),'=x,y,z lattice spacing'
      ENDIF
!***********************************************************************
!*** diagnostic
!      write(0,*)'tarell,ckpt 6, ioshp=',ioshp
!***
      IF(IOSHP>=0)THEN
         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')
         WRITE(IOSHP,FMT=9020)AX,AY,AZ,NAT,A1,A2,DX,X0
         DO JX=1,NAT
            WRITE(IOSHP,FMT=9030)JX,IXYZ(JX,1),IXYZ(JX,2),IXYZ(JX,3), &
                                 ICOMP(JX,1),ICOMP(JX,2),ICOMP(JX,3)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF
      RETURN
9020  FORMAT(' >TARELL  ellipsoidal grain; AX,AY,AZ=',3F8.4,/,     &
         I10,' = NAT',/,                                           &
         3F10.6,' = A_1 vector',/,                                 &
         3F10.6,' = A_2 vector',/,                                 &
         3F10.6,' = lattice spacings (d_x,d_y,d_z)/d',/,           &
         3F10.5,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ', &
               'for dipole 0 0 0',/,                               &
         '     JA  IX  IY  IZ ICOMP(x,y,z)')
9030  FORMAT (I7,3I5,3I2)
    END SUBROUTINE TARELL
