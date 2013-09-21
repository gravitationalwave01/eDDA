    SUBROUTINE TARHEX(A1,A2,A,B,ORI,DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

! Scalar arguments:

      CHARACTER :: CDESCR*67
      INTEGER :: IOSHP,MXNAT,NAT
      INTEGER*2 :: ICOMP(MXNAT,3)
      INTEGER :: IXYZ(MXNAT,3)
      REAL(WP) :: A,B,ORI

! Array arguments:

      REAL(WP) :: A1(3),A2(3),DX(3),X0(3)

! Local scalars:

      CHARACTER :: CMSGNM*70
      INTEGER :: IORI,JA,JX,JY,JZ,NA,NB,NC,NFAC,NLAY,NMX,NMY,NMZ
      LOGICAL :: IN
      REAL(WP) :: ACM,ASPR,BCM,BEFF,CCM,X,Y,Z

!***********************************************************************
! Routine to construct hexagonal prism from "atoms"
! Input:
!        A = prism length/d
!        B = (vertex-vertex diameter of hexagon face)/d
!          = 2*(one hexagon side)/d
!        ORI = 1. for a1 in x direction, a2 in y direction
!            = 2. for a1 in x direction, a2 in z direction
!            = 3. for a1 in y direction, a2 in x direction
!            = 4. for a1 in y direction, a2 in z direction
!            = 5. for a1 in z direction, a2 in x direction
!            = 6. for a1 in z direction, a2 in y direction
!              (a1 is parallel to prism axis [normal to hexagonal face]
!               a2 is normal to rectangular face
!        DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z
!                directions, and d=(dx*dy*dz)**(1/3)=effective lattice
!                spacing
!        IOSHP=device number for "target.out" file
!             =-1 to suppress printing of "target.out"
!        MXNAT=dimensioning information (max number of atoms)
! Output:
!        A1(1-3)=unit vector (1,0,0) defining target axis 1 (prism axis)
!        A2(1-3)=unit vector (0,1,0) defining target axis 2 (normal to
!                one of the rectangular faces.
!        X0(1-3)=location/d in TF of dipole with IXYZ=0 0 0
!        NAT=number of atoms in target
!        IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
!        CDESCR=string describing target (up to 67 chars.)

! B.T.Draine, Princeton Univ. Obs., 87.04.03

! History:
! 90.12.02 (BTD): IOSHP is now output device for target.out
!                 (and is now closed before returning)
! 90.12.14 (BTD): Rewritten so that atoms are always on integral sites
!                 but CM coords may be either integral or half-integral
! 91.01.05 (BTD): Change I4 -> I7 when printing NAT.
! 91.04.22 (BTD): Added XV,YV,ZV and A1,A2 to WRITE(IOSHP,FMT=9020)
! 92.09.09 (BTD+PJF): in output, NLAY was too large by factor 2 and NFAC
!                 too small by factor of 2: now fixed.
! 93.03.12 (BTD): Changed CDESCR*(*) -> CDESCR*67
! 95.07.13 (BTD): Simplify: restrict to case where prism axis is in
!                 x direction and y axis is normal to a rectangular
!                 face.
! 96.01.26 (BTD): Replaced IX,IY,IZ by IXYZ
! 96.02.23 (BTD): Delete PI (unused)
! 97.04.28 (BTD): Initialize NAT=0 (failure to do so reported by
!                 Michael Wolff).
! 97.12.26 (BTD): Added DX(3) to argument list to support nonuniform
!                 lattice spacing.
!                 *** Note: additional code modifications required
!                     to support noncubic lattice!
! 98.03.06 (BTD): If called with noncubic lattice (not yet supported)
!                 call ERRMSG and halt with fatal error.
! 98.03.07 (BTD): Modify to write DX to file "target.out" for
!                 compatibility with REASHP
! 00.11.02 (BTD): Added ICOMP to argument list
!                 set ICOMP=1 for occupied sites
!                 write ICOMP to target.out
! 06.09.15 (BTD): Generalized to support 6 different orientations in TF.
!                 Now uses ORI=SHPAR(3) to select orientation.
! 07.09.11 (BTD): Changed IXYZ from INTEGER*2 to INTEGER
! 08.08.08 (BTD): Added X0(3) to argument list
!                 Added code to set X0(1-3)
! 08.08.29 (BTD): Modified format 9020
! end history

! Copyright (C) 1993,1995,1996,1997,1998,2000,2006,2007,2008
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

      IORI=NINT(ORI)

! Current version of TARHEX is restricted to cubic lattices

      IF(DX(1)/=1._WP.OR.DX(2)/=1._WP)THEN
         CALL ERRMSG('FATAL','TARHEX',                          &
                     ' tarhex does not support noncubic lattice')
      ENDIF
      DO JA=1,3
         A1(JA)=0._WP
         A2(JA)=0._WP
      ENDDO

! A = prism length
! B = prism "diameter"
!       B/2 = length of one side of hexagon
!       B= distance between opposite vertices=max.diameter of hexagon
!       for hexagon, area=(3/8)*SQRT(3)*B**2
!       for hex prism, volume=A*area=(3/8)*SQRT(3)*A*B**2
! Now determine limits for testing x,y,z values
! Along axis, run from 1 to NA=INT(A+0.5)
! Perp. to axis, run from 1 to NB=INT(B+0.5)
      NA=INT(A+0.5_WP)
      NB=INT(B+0.5_WP)
      NC=INT(.5_WP*SQRT(3._WP)*B+0.5_WP)
! ACM="center" of figure in axial direction
! BCM="center" of figure in B direction
! CCM="center" of figure in C direction
      ACM=REAL(NA,KIND=WP)/2._WP+0.5_WP
      BCM=REAL(NB,KIND=WP)/2._WP+0.5_WP
      CCM=REAL(NC,KIND=WP)/2._WP+0.5_WP
      IF(IORI==1)THEN
         NMX=NA
         NMY=NC
         NMZ=NB
         A1(1)=1._WP
         A2(2)=1._WP
      ELSEIF(IORI==2)THEN
         NMX=NA
         NMZ=NC
         NMY=NB
         A1(1)=1._WP
         A2(3)=1._WP
      ELSEIF(IORI==3)THEN
         NMY=NA
         NMX=NC
         NMZ=NB
         A1(2)=1._WP
         A2(1)=1._WP
      ELSEIF(IORI==4)THEN
         NMY=NA
         NMZ=NC
         NMX=NB
         A1(2)=1._WP
         A2(3)=1._WP
      ELSEIF(IORI==5)THEN
         NMZ=NA
         NMX=NC
         NMY=NB
         A1(3)=1._WP
         A2(1)=1._WP
      ELSEIF(IORI==6)THEN
         NMZ=NA
         NMY=NC
         NMX=NB
         A1(3)=1._WP
         A2(2)=1._WP
      ELSE
         CALL ERRMSG('FATAL','TARHEX',' invalid value for SHPAR(3)')
      ENDIF
      NAT=0
      DO JZ=1,NMZ
         Z=REAL(JZ,KIND=WP)
         DO JY=1,NMY
            Y=REAL(JY,KIND=WP)
            DO JX=1,NMX
               X=REAL(JX,KIND=WP)
               IF(IORI==1)THEN
                  CALL TESTHEX(ACM,BCM,CCM,A,B,X,Y,Z,IN)
               ELSEIF(IORI==2)THEN
                  CALL TESTHEX(ACM,BCM,CCM,A,B,X,Z,Y,IN)
               ELSEIF(IORI==3)THEN
                  CALL TESTHEX(ACM,BCM,CCM,A,B,Y,X,Z,IN)
               ELSEIF(IORI==4)THEN
                  CALL TESTHEX(ACM,BCM,CCM,A,B,Y,Z,X,IN)
               ELSEIF(IORI==5)THEN
                  CALL TESTHEX(ACM,BCM,CCM,A,B,Z,X,Y,IN)
               ELSEIF(IORI==6)THEN
                  CALL TESTHEX(ACM,BCM,CCM,A,B,Z,Y,X,IN)
               ENDIF
               IF(IN)THEN
                  NAT=NAT+1
                  IXYZ(NAT,1)=JX
                  IXYZ(NAT,2)=JY
                  IXYZ(NAT,3)=JZ
               ENDIF
            ENDDO
         ENDDO
      ENDDO

! Set X0 = location in TF of IXYZ=0 0 0
! We assume that origin of TF is located at centroid of hex prism
! Set composition

      DO JZ=1,3
         X0(JZ)=0._WP
         DO JX=1,NAT
            ICOMP(JX,JZ)=1
            X0(JZ)=X0(JZ)+REAL(IXYZ(JX,JZ),KIND=WP)
         ENDDO
         X0(JZ)=-X0(JZ)/REAL(NAT,KIND=WP)
      ENDDO

! Now compute some geometric properties of pseudo-hexagonal prism
!   NLAY = number of layers = effective length of prism/d
!   NFAC = number of atoms in one hexagonal face

      NLAY=INT(A)
      NFAC=NAT/NLAY

! BEFF = effective length of hexagonal size/d
!        computed for hexagon of area NFAC*d**2
!                                     =(3/2)*SQRT(3)*BEFF**2
! Note: for hexagon, vertex-vertex diameter = 2*BEFF

      BEFF=.6204032_WP*SQRT(REAL(NFAC,KIND=WP))

! ASPR = effective aspect ratio of target = length/(2*BEFF)

      ASPR=.5_WP*REAL(NLAY,KIND=WP)/BEFF

! Note: string CDESCR will be printed by subr. TARGET

      WRITE(CDESCR,FMT='(A,I7,A)')' Hexagonal prism of NAT=',NAT,' dipoles'

! Here print any additional target info which is desired by using
! subroutine WRIMSG

      WRITE(CMSGNM,FMT='(A,3F7.4,A)')' A1 = (',A1,') = hexagon symmetry axis'
      CALL WRIMSG('TARHEX',CMSGNM)
      WRITE(CMSGNM,FMT='(A,3F7.4,A)')' A2 = (',A2, ') = center-vertex axis'
      CALL WRIMSG('TARHEX',CMSGNM)
      WRITE(CMSGNM,FMT='(A,I7,A,I4,A,I4,A,F7.4)')' NAT=',NAT,' NFAC=', &
                                NFAC,' NLAY=',NLAY,' aspect ratio=',ASPR
      CALL WRIMSG('TARHEX',CMSGNM)
      WRITE(CDESCR,FMT='(A,I7,A,I4,A,I4,A,F7.4)')' Hex prism,NAT=',NAT, &
                           ' NFAC=',NFAC,' NLAY=',NLAY,' asp.ratio=',ASPR

! Now store atom coordinates in file

      IF(IOSHP>0)THEN
         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')
         WRITE(IOSHP,9020)A,B,NAT,A1,A2,DX,X0
         DO JX=1,NAT
            WRITE(IOSHP,FMT=9030)JX,IXYZ(JX,1),IXYZ(JX,2),IXYZ(JX,3), &
                                 ICOMP(JX,1),ICOMP(JX,2),ICOMP(JX,3)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF
      RETURN
9020  FORMAT(' Hexagonal Prism:  A=',F7.4,' B=',F7.4,/,            &
         I10,' = NAT ',/,                                          &
         3F10.6,' = A_1 vector',/,                                 &
         3F10.6,' = A_2 vector',/,                                 &
         3F10.6,' = lattice spacings (d_x,d_y,d_z)/d',/,           &
         3F10.5,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ', &
               'for dipole 0 0 0',/,                               &
         '     JA  IX  IY  IZ ICOMP(x,y,z)')

9030  FORMAT(I7,3I4,3I2)
    END SUBROUTINE TARHEX

    SUBROUTINE TESTHEX(ACM,BCM,CCM,A,B,X,Y,Z,IN)
      USE DDPRECISION,ONLY : WP
! Arguments:
      REAL(WP) :: ACM,BCM,CCM,A,B,X,Y,Z
      LOGICAL :: IN
! Local variables:
      REAL(WP) :: AX,Q,U,V

!*** ACM,BCM,CCM=location of center in A,B,C directions
!    A direction is along axis
!    B direction is vertex to vertex of hexagon
!    C direction is face to face of hexagon

      AX=X
      U=Y-CCM
      V=Z-BCM
      IN=.FALSE.
!*** Test along axis:
      IF(2._WP*ABS(AX-ACM)>A) RETURN
!*** Test along C direction
! .4330127=SQRT(3)/4
      IF(ABS(U)>0.4330127_WP*B) RETURN
!*** Now test whether closer to CM than line bounding edge
! .5773503=1/SQRT(3)
      Q=0.5_WP*ABS(U)+0.8660254_WP*ABS(V)
      IF(Q>0.4330127_WP*B) RETURN
      IN=.TRUE.
      RETURN
    END SUBROUTINE TESTHEX
