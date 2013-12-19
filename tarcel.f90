    SUBROUTINE TARCEL(A1,A2,AX,AY,AZ,BX,BY,BZ,DX,X0,CDESCR,IOSHP,MXNAT,NAT, &
                      IXYZ,ICOMP)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

!** Arguments:

      CHARACTER :: CDESCR*67
      INTEGER :: IOSHP,MXNAT,NAT
      INTEGER*2 ::    &
         ICOMP(MXNAT,3)
      INTEGER ::     &
         IXYZ(MXNAT,3)
      REAL(WP) :: AX,AY,AZ,BX,BY,BZ
      REAL(WP) :: & 
         A1(3),   &
         A2(3),   &
         DX(3),   &
         X0(3)

!** Local variables:

      INTEGER :: JX,JY,JZ,LMX1,LMX2,LMY1,LMY2,LMZ1,LMZ2,NIN
      REAL(WP) :: AX2,AY2,AZ2,BX2,BY2,BZ2,R,RR,RYZ2,RZ2,X,XOFF,Y,YOFF,Z,ZOFF
!***********************************************************************

! Routine to construct target consisting of two materials, with outer
! surface an ellipsoid of dimensions AX,AY,AZ,
! and core/mantle interface a concentric ellipsoid of dimensions BX,BY,BZ
! Input:
!        AX=(x-length)/d of outer ellipsoid  (d=lattice spacing)
!        AY=(y-length)/d "   "     "
!        AZ=(z-length)/d "   "     "
!        BX=(x-length)/d of inner ellipsoid
!        BY=(y-length)/d "   "     "
!        BZ=(z-length)/d "   "     "
!        DX(1-3)=lattice spacing (dx/d,dy/d,dz/d), d=(dx*dy*dz)**(1/3)
!        IOSHP=device number for "target.out" file
!             =-1 to suppress printing of "target.out"
!        MXNAT=dimensioning information (max number of atoms)
! Output:
!        A1(1-3)=(1,0,0)=unit vector defining target axis 1 in Target Frame
!        A2(1-3)=(0,1,0)=unit vector defining target axis 2 in Target Frame
!        NAT=number of atoms in target
!        IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
!        ICOMP(1-NAT,1-3)=1 for sites within inner ellipsoid,
!                        =2 for sites between inner and outer ellipsoids
!        CDESCR=description of target (up to 67 characters)
!        X0(1-3)=(location/d) in Target Frame corresponding to dipole with
!                IXYZ=(0,0,0).  This will be treated at the origin of physical
!                coordinates in the TF.
!                Here set to be centroid of the ellipsoids.

! B.T.Draine, Princeton Univ. Obs.
! History:
! 94.01.26 (AP) : Adapted from TARELL by Antonio Peimbert, Princeton U.
! 94.01.27 (BTD): Modified for compatibility with TARGET, using
!                 AX,AY,AZ,BX,BY,BZ to transfer shape parameters read by
!                 subroutine REAPAR
! 94.03.21 (BTD): Corrected format statement
! 96.01.26 (BTD): Replaced IX,IY,IZ by IXYZ
! 97.12.26 (BTD): Added DX(3) to argument list to support nonuniform
!                 lattice spacing.
!                 *** Note: additional code modifications required
!                     to support noncubic lattice!
! 98.03.06 (BTD): If called with noncubic lattice (not yet supported)
!                 call ERRMSG and halt with fatal error.
! 98.03.07 (BTD): Modify to write DX to file "target.out" for
!                 compatibility with REASHP
! 00.11.02 (BTD): Remove XOFF,YOFF,ZOFF from target.out
!                 cosmetic changes
! 07.09.11 (BTD): Changed IXYZ from INTEGER*2 to INTEGER
! 08.01.13 (BTD): Cosmetic changes
!                 Added X0 to argument list
!                 Added code setting X0 to be centroid of ellipsoid
! 08.08.29 (BTD): Modified format 9020
! end history

! Copyright (C) 1994,1996,1997,1998,2007,2008 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

! Current version of TARCEL is restricted to cubic lattices

      IF(DX(1)/=1._WP.OR.DX(2)/=1._WP)THEN
         CALL ERRMSG('FATAL','TARCEL',                          &
                     ' tarcel does not support noncubic lattice')
      ENDIF

! Inner ellipsoid is defined by
!      (X/AX)**2+(Y/AY)**2+(Z/AZ)**2 = 0.25
! Outer ellipsoid is defined by
!      (X/BX)**2+(Y/BY)**2+(Z/BZ)**2 = 0.25
! Material within inner ellipsoid is of composition 1
! Material between inner and outer ellipsoids is of composition 2

! Ideal volume V=(pi/6)*BX*BY*BZ

! Dipoles are located at sites
! (x,y,z)=(I+XOFF,J+YOFF,Z+KOFF), I,J,K=integers
!                                 XOFF,YOFF,ZOFF=constants

! For concentric spheres: call with AX=AY=AZ, BX=BY=BZ
! For spheroids: call with AX=AY (or AY=AZ or AX=AZ) and BX=BY (or BY=BZ
!   or BX=BZ)

!***********************************************************************

! Check that input parameters have "inner" ellipsoid smaller than
! "outer" ellipsoid:
      IF(AX<BX)CALL ERRMSG('FATAL','TARCEL',' AX < BX ')
      IF(AY<BY)CALL ERRMSG('FATAL','TARCEL',' AY < BY ')
      IF(AZ<BZ)CALL ERRMSG('FATAL','TARCEL',' AZ < BZ ')
! Criterion for choosing XOFF: try to optimize outer ellipsoid surface
! If AX is close to an even integer, take XOFF=1/2
! If AX is close to an odd integer, take XOFF=0
      JX=INT(AX+.5_WP)
      IF(2*(JX/2)==JX)THEN
         XOFF=0.5_WP
      ELSE
         XOFF=0._WP
      ENDIF
! Similar criterion for YOFF:
      JY=INT(AY+.5_WP)
      IF(2*(JY/2)==JY)THEN
         YOFF=0.5_WP
      ELSE
         YOFF=0._WP
      ENDIF
! Similar criterion for ZOFF:
      JZ=INT(AZ+.5_WP)
      IF(2*(JZ/2)==JZ)THEN
         ZOFF=0.5_WP
      ELSE
         ZOFF=0._WP
      ENDIF

      LMX1=-INT(.5_WP*AX+.5_WP)
      LMX2=INT(.5_WP*AX-.25_WP)
      LMY1=-INT(.5_WP*AY+.5_WP)
      LMY2=INT(.5_WP*AY-.25_WP)
      LMZ1=-INT(.5_WP*AZ+.5_WP)
      LMZ2=INT(.5_WP*AZ-.25_WP)

      AX2=AX*AX
      AY2=AY*AY
      AZ2=AZ*AZ
      BX2=BX*BX
      BY2=BY*BY
      BZ2=BZ*BZ

      NAT=0
      NIN=0

! Specify target axes A1 and A2
!   A1=(1,0,0) in target frame
!   A2=(0,1,0) in target frame

      DO JX=1, 3
         A1(JX)=0._WP
         A2(JX)=0._WP
      ENDDO
      A1(1)=1._WP
      A2(2)=1._WP

! Determine list of occupied sites.

      DO JZ=LMZ1,LMZ2
         Z=REAL(JZ,KIND=WP)+ZOFF
         RZ2=Z*Z/AZ2
         IF(RZ2<0.25_WP)THEN
            DO JY=LMY1, LMY2
               Y=REAL(JY,KIND=WP)+YOFF
               RYZ2=RZ2+Y*Y/AY2
               IF(RYZ2<0.25_WP)THEN
                  DO JX=LMX1,LMX2
                     X=REAL(JX,KIND=WP)+XOFF
                     R=RYZ2+X*X/AX2
                     IF(R<0.25_WP)THEN
! Site is occupied:
                        NAT=NAT+1
                        RR=Z*Z/BZ2+Y*Y/BY2+X*X/BX2
                        IF(RR<0.25_WP)THEN
! within inner ellipsoid:
                           ICOMP(NAT,1)=1
                           ICOMP(NAT,2)=1
                           ICOMP(NAT,3)=1
                           NIN=NIN+1
                        ELSE
! between inner and outer ellipsoids:
                           ICOMP(NAT,1)=2
                           ICOMP(NAT,2)=2
                           ICOMP(NAT,3)=2
                        ENDIF
                        IXYZ(NAT,1)=JX
                        IXYZ(NAT,2)=JY
                        IXYZ(NAT,3)=JZ
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO

! locate centroid = origin of coordinates in TF

      DO JX=1,3
         X0(JX)=0._WP
      ENDDO
      DO JX=1,3
         DO JY=1,NAT
            X0(JX)=X0(JX)+REAL(IXYZ(JY,JX))
         ENDDO
         X0(JX)=X0(JX)/REAL(NAT)
      ENDDO
      
      IF(NAT>MXNAT)THEN
         CALL ERRMSG('FATAL','TARELL',' NAT.GT.MXNAT ')
      ENDIF

!***********************************************************************
! Write target description into string CDESCR

      IF(AX==AY .AND. AX==AZ)THEN
         WRITE(CDESCR,FMT='(A,I7,A)')' Spherical target containing',NAT, &
                                     ' dipoles'
      ELSE
         WRITE(CDESCR,FMT='(A,I6,A,I7,A)')' Concentric ellipsoids with',NIN, &
                                          ' inner ', NAT, ' total dipoles'
      ENDIF

!***********************************************************************

      IF(IOSHP>=0)THEN
         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')
         WRITE(IOSHP,FMT=9020)AX,AY,AZ,BX,BY,BZ,NAT,NIN,A1,A2,DX,X0
         DO JX=1,NAT
            WRITE(IOSHP,FMT=9030)JX,IXYZ(JX,1),IXYZ(JX,2),IXYZ(JX,3), &
                                 ICOMP(JX,1),ICOMP(JX,2),ICOMP(JX,3)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF
      RETURN
9020  FORMAT(' >TARCEL: concentric ellipsoids; AX,AY,AZ=',3F8.4,   &
             ' BX,BY,BZ=',3F8.4,/,                                 &
         2I10,' = NAT, NIN,',/,                                    &
         3F10.6,' = A_1 vector',/,                                 &
         3F10.6,' = A_2 vector',/,                                 &
         3F10.6,' = lattice spacings (d_x,d_y,d_z)/d',/,           &
         3F10.5,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ', &
               'for dipole 0 0 0',/,                               &
         '     JA  IX  IY  IZ ICOMP(x,y,z)')
9030  FORMAT(I7,3I4,3I2)
    END SUBROUTINE TARCEL
