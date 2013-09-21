    SUBROUTINE TAR3EL(A1,A2,AX,AY,AZ,DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE
!** Arguments:
      CHARACTER :: CDESCR*67
      INTEGER :: IOSHP,MXNAT,NAT
      INTEGER*2 ::    &
         ICOMP(MXNAT,3)
      INTEGER ::     &
         IXYZ(MXNAT,3)
      REAL(WP) :: AX,AY,AZ
      REAL(WP) :: & 
         A1(3),   &
         A2(3),   &
         DX(3),   &
         X0(3)

!** Local variables:
      CHARACTER :: CMSGNM*70
      INTEGER :: JA,JX,JXMAX,JXMIN,JY,JZ,LMX1,LMX2,LMY1,LMY2,LMZ1, &
        LMZ2,NAT2
      REAL(WP) :: AX2,AY2,AZ2,R,RYZ2,RZ2,X,XOFF,Y,YOFF,Z,ZOFF
!***********************************************************************

! Routine to construct three collinear touching ellipsoids from "atoms"
! Input:
!        AX=(x-length of one ellipsoid)/d    (d=lattice spacing)
!        AY=(y-length of one ellipsoid)/d
!        AZ=(z-length of one ellipsoid)/d
!        DX(1-3)=lattice spacing(dx,dy,dz)/d [d=(dx*dy*dz)**(1/3)]
!        IOSHP=device number for "target.out" file
!             =-1 to suppress printing of "target.out"
!        MXNAT=dimensioning information (max number of atoms)
! Output:
!        A1(1-3)=(1,0,0)=unit vector defining target axis 1 in Target Fr
!        A2(1-3)=(0,1,0)=unit vector defining target axis 2 in Target Fr
!        NAT=number of dipoles in target
!        IXYZ(1-NAT,1-3)=x/d,y/d,z/d for target dipole locations
!        ICOMP(1-NAT,1-3)=x,y,z composition for dipole locations
!                        = 1 for locations in first ellipsoid
!                          2                  second
!                          3                  third
!        CDESCR=description of target (up to 67 characters)
! Note: atoms 1       - NAT/3   are in first ellipsoid
!             NAT/3+1 - 2*NAT/3 are in second ellipsoid
!             2*NAT/3+1 - NAT   are in third ellipsoid
! Ellipsoids are displaced from one another in x-direction
! First ellipsoid is at smallest x values, third at largest
! B.T.Draine, Princeton Univ. Obs.

! History:
! 93.01.20 (BTD): adapted from tar2el.f
! 93.03.12 (BTD): changed CDESCR*(*) -> CDESCR*67
! 96.01.26 (BTD): Replaced IX,IY,IZ by IXYZ
! 97.12.26 (BTD): Added DX(3) to argument list to support nonuniform
!                 lattice spacing.
!                 *** Note: additional code modifications required
!                     to support noncubic lattice!
! 98.03.06 (BTD): If called with noncubic lattice (not yet supported)
!                 call ERRMSG and halt with fatal error.
! 98.03.07 (BTD): Modify to write DX to file "target.out" for
!                 compatibility with REASHP
! 00.11.02 (BTD): Add ICOMP to argument list
!                 Set ICOMP=1,2,3 for first,second,third ellipsoid
!                 Write ICOMP to target.out
! 07.09.11 (BTD): Changed IXYZ from INTEGER*2 to INTEGER
! 08.01.13 (BTD): Cosmetic changes
!                 Added X0 to argument list
!                 Added code to evaluate X0, assuing TF origin to be
!                 located at centroid of 3 ellipses.
! 08.08.30 (BTD): Modified format 9020
! end history

! Copyright (C) 1993,1996,1997,1998,2000,2007,2008 
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

! Current version of TAR3EL is restricted to cubic lattices

      IF(DX(1)/=1._WP.OR.DX(2)/=1._WP)THEN
         CALL ERRMSG('FATAL','TAR3EL', &
                     ' tar3el does not support noncubic lattice')
      ENDIF

! Routine to construct pseudo-ellipsoidal target aray.
! With occupied array sites contained within ellipsoidal surface
! defined by (X/AX)**2+(Y/AY)**2+(Z/AZ)**2=0.25
! Ideal volume V=(pi/6)*AX*AY*AZ

! Dipoles are located at sites
! (x,y,z)=(I+XOFF,J+YOFF,Z+KOFF), I,J,K=integers
!                                 XOFF,YOFF,ZOFF=constants

! For sphere: call with AX=AY=AZ
! For spheroid: call with AX=AY (or AY=AZ or AX=AZ)
! B.T.Draine, Princeton Univ. Obs., 88/8/12

! Criterion for choosing XOFF:
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
      NAT=0

! Specify target axes A1 and A2
! Previous convention: A1 along longest dimension
!                      A2 along intermediate dimension
! New convention: A1=(1,0,0) in target frame
!                 A2=(0,1,0) in target frame

      DO JX=1,3
         A1(JX)=0._WP
         A2(JX)=0._WP
      ENDDO
      A1(1)=1._WP
      A2(2)=1._WP

! Determine list of occupied sites in first ellipsoid

      JXMAX=LMX1
      JXMIN=LMX2
      DO JZ=LMZ1,LMZ2
         Z=REAL(JZ,KIND=WP)+ZOFF
         RZ2=Z*Z/AZ2
         IF(RZ2<0.25_WP)THEN
            DO JY=LMY1,LMY2
               Y=REAL(JY,KIND=WP)+YOFF
               RYZ2=RZ2+Y*Y/AY2
               IF(RYZ2<0.25_WP)THEN
                  DO JX=LMX1,LMX2
                     X=REAL(JX,KIND=WP)+XOFF
                     R=RYZ2+X*X/AX2
                     IF(R<0.25_WP)THEN
! Site is occupied:
                        NAT=NAT+1
                        IXYZ(NAT,1)=JX
                        IXYZ(NAT,2)=JY
                        IXYZ(NAT,3)=JZ
                        IF(NAT>1)THEN
                           IF(JX>JXMAX) JXMAX=JX
                           IF(JX<JXMIN) JXMIN=JX
                        ELSE
                           JXMAX=JX
                           JXMIN=JX
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      DO JX=1,3
         DO JA=1,NAT
            ICOMP(JA,JX)=1
         ENDDO
      ENDDO

      IF(3*NAT>MXNAT)THEN
         CALL ERRMSG('FATAL','TARELL',' NAT.GT.MXNAT ')
      ENDIF

! Now create duplicate second and third ellipsoids

      JXMAX=JXMAX-JXMIN+1
      DO JA=1,NAT
         IXYZ(NAT+JA,1)=IXYZ(JA,1)+JXMAX
         IXYZ(NAT+JA,2)=IXYZ(JA,2)
         IXYZ(NAT+JA,3)=IXYZ(JA,3)
         IXYZ(2*NAT+JA,1)=IXYZ(JA,1)+2*JXMAX
         IXYZ(2*NAT+JA,2)=IXYZ(JA,2)
         IXYZ(2*NAT+JA,3)=IXYZ(JA,3)
      ENDDO
      NAT=3*NAT
      DO JX=1,3
         DO JA=NAT/3+1,2*NAT/3
            ICOMP(JA,JX)=2
         ENDDO
         DO JA=2*NAT/3+1,NAT
            ICOMP(JA,JX)=3
         ENDDO
      ENDDO

! Now find X0, assuming TF origin to be at centroid

      DO JY=1,3
         X0(JY)=0._WP
         DO JA=1,NAT
            X0(JY)=X0(JY)+REAL(IXYZ(JA,JY))
         ENDDO
         X0(JY)=-X0(JY)/REAL(NAT)
      ENDDO

!***********************************************************************
! Write target description into string CDESCR

      NAT2=NAT/3
      IF(AX==AY.AND.AX==AZ)THEN
         WRITE(CDESCR,FMT='(A,I7,A)')' Three spheres, each containing', &
                                     NAT2,' dipoles'
      ELSE
         WRITE(CDESCR,FMT='(A,I7,A)')' Three ellipsoids, each containing', &
                                     NAT2,' dipoles'
      ENDIF
!***********************************************************************
      CMSGNM=CDESCR
      CALL WRIMSG('TARELL',CMSGNM)
      IF(IOSHP>=0)THEN
         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')
         WRITE(IOSHP,FMT=9020) AX,AY,AZ,NAT,A1,A2,DX,X0
         DO JX=1,NAT
            WRITE(IOSHP,FMT=9030)JX,IXYZ(JX,1),IXYZ(JX,2),IXYZ(JX,3), &
                                 ICOMP(JX,1),ICOMP(JX,2),ICOMP(JX,3)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF
      RETURN
9020  FORMAT(' >TARELL  ellipsoidal grain; AX,AY,AZ=',3F8.4,/,    &
         I10,' = NAT',/,                                           &
         3F10.6,' = A_1 vector',/,                                 &
         3F10.6,' = A_2 vector',/,                                 &
         3F10.6,' = lattice spacings (d_x,d_y,d_z)/d',/,           &
         3F10.5,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ', &
               'for dipole 0 0 0',/,                              &
         '     JA  IX  IY  IZ ICOMP(x,y,z)')
9030  FORMAT(I7,3I4,3I2)
    END SUBROUTINE TAR3EL
