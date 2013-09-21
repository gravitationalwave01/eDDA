    SUBROUTINE TARTET(A1,A2,AX,DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

!** Arguments:
      CHARACTER :: CDESCR*67
      INTEGER :: IOSHP,MXNAT,NAT
      INTEGER*2 :: ICOMP(MXNAT,3)
      INTEGER :: IXYZ(MXNAT,3)
      REAL(WP) :: AX
      REAL(WP) :: A1(3),A2(3),DX(3),X0(3)
!** Local variables:
      INTEGER :: I,IMAX,IMIN,J,JMAX,JMIN,JX,K,KMAX,KMIN
      REAL(WP) :: FY,S,X,XOFF,Y,YMAX,YMIN,YOFF,Z,ZMAX,ZMAX0,ZOFF
!***********************************************************************

! Routine to construct regular tetrahedral target array
! one of the faces is parallel to the y-z plane
! one of the edges of this face is parallel to x-y plane

! Input:
!        AX = length of one edge of tetrahedron
!        DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z
!                directions, and d=(dx*dy*dz)**(1/3)=effective lattice
!                spacing
!        IOSHP=device number for "target.out" file
!             =-1 to suppress printing of "target.out"
!        MXNAT = dimensioning information
! Returns:
!        A1(1-3) = unit vector (1,0,0) (along one axis of tetrahedron)
!        A2(1-3) = unit vector (0,1,0)
!        CDESCR = string describing target (up to 67 chars.)
!        NAT = number of atoms in target
!        IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms in target
!        ICOMP(1-NAT,1-3)=1 (composition identifier)
!        X0(1-3) = location/d in TF of site with IXYZ=0 0 0
! Occupied array sites are those within tetrahedral surface
! Size:
!    S=length of each side of tetrahedron
!    Volume = S**3/(6*sqrt(2))
! Orientation:
!    Center of mass is at origin.
!    One face is parallel to yz plane.
!    Projection onto yz plane has vertex on y axis.
!    S=length of one side (in lattice units)
!    Vertices are at
!    A=( sqrt(3/8),         0,    0)*S
!    B=(-sqrt(1/24), sqrt(1/3),    0)*S
!    C=(-sqrt(1/24),-sqrt(1/12),-1/2)*S
!    D=(-sqrt(1/24),-sqrt(1/12), 1/2)*S
! Length in x direction = S*sqrt(2/3)
!           y           = S*sqrt(3)/2
!           z           = S
! Angle(AOB)=arccos(-1/3)=109.4712 deg.
! OA=OB=OC=OD=sqrt(3/8)*S
! Occupied sites are assumed to be located at
! (X,Y,Z)=(I+XOFF,J+YOFF,K+ZOFF)
! where I,J,K are integers, and XOFF,YOFF,ZOFF are constants.
! Program sets XOFF,YOFF,ZOFF depending on choice of parameter S.

! Criterion for choosing XOFF:
!    Base of tetrahedron is located at x = -sqrt(1/24)*S
!    Let IMIN be value of I for this plane
!    Choose XOFF so that IMIN+XOFF = -sqrt(1/24)*S + 0.5
!    with -0.5 < XOFF < 0.5
! Criterion for choosing YOFF:
!    One edge of tetrahedron is located at y= -S/(2*sqrt(3))
!    Let JMIN be value of J for this line
!    Choose YOFF so that JMIN+YOFF = -S/sqrt(12) + 0.5
! Criterion for choosing ZOFF:
!    One edge of tetrahedron is parallel to z axis
!    Choose ZOFF in order to have number of dipoles along
!    this edge as close as possible to S
!    e.g., if S=odd integer, then take ZOFF=0 to place dipoles
!          at integral locations
!          if S=even integer, then take ZOFF=1/2 to place dipoles
!          at half-integral locations

! B.T.Draine, Princeton Univ. Obs., 90.10.30
! History:
! 90.11.08 (BTD): Changed DO...ENDDO to DO #...# CONTINUE
! 90.11.09 (BTD): Corrected error in location of center of mass.
! 91.01.05 (BTD): Change I4 -> I7 when printing NAT.
! 91.04.22 (BTD): Added XV,YV,ZV and A1,A2 to WRITE(IOSHP,FMT=9020)
! 95.07.13 (BTD): Eliminated unused arguments AY,AZ
! 96.01.26 (BTD): Replaced IX,IY,IZ by IXYZ
! 97.12.26 (BTD): Added DX(3) to argument list to support noncubic
!                 lattice.
!                 *** Note: additional code modifications required
!                     to support noncubic lattice!
! 98.03.06 (BTD): If called with noncubic lattice (not yet supported)
!                 call ERRMSG and halt with fatal error.
! 98.03.07 (BTD): Modify to write DX to file "target.out" for
!                 compatibility with REASHP
! 00.11.02 (BTD): Add ICOMP to argument list
!                 set ICOMP=1 at occupied sites
!                 write ICOMP to target.out
! 07.09.11 (BTD): changed IXYZ from INTEGER*2 to INTEGER
! 08.08.30 (BTD): added X0 to argument list
!                 added code to evaluate X0 for origin of TF at centroid
!                 modified format 9020
! end history

! Copyright (C) 1993,1995,1996,1997,1998,2000,2007,2008
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

! Current version of TARTET is restricted to cubic lattices

      IF(DX(1)/=1._WP.OR.DX(2)/=1._WP)THEN
         CALL ERRMSG('FATAL','TARTET',' tartet does not support noncubic lattice')
      ENDIF
      S=AX

! Set XOFF (and IMIN,IMAX):

      IMIN=-INT(S*SQRT(1._WP/24._WP))
      XOFF=0.5_WP-S*SQRT(1._WP/24._WP)-IMIN
      IMAX=IMIN+INT(S*SQRT(2._WP/3._WP)+0.5_WP)-1

! Set YOFF (and JMIN,JMAX):

      JMIN=-INT(S/SQRT(12._WP))
      YOFF=0.5_WP-S/SQRT(12._WP)-JMIN
      JMAX=JMIN+INT(S*SQRT(.75_WP)+0.5_WP)-1

! Set ZOFF (and KMIN,KMAX):
! Determine whether S is closest to even or odd integer.
!  (Temporarily let KMIN be integer which S is closest to)

      KMIN=INT(S+0.5_WP)

! If KMIN is even, then ZOFF=0.5
! If KMIN is odd, then ZOFF=0.

      ZOFF=0._WP
      IF(KMIN-2*(KMIN/2)==0)ZOFF=0.5_WP
      KMIN=-INT(0.5_WP*S+ZOFF)
      KMAX=KMIN+INT(S+0.5_WP)-1

! Determine list of occupied sites.

      NAT=0
      DO I=IMIN,IMAX
         X=REAL(I,KIND=WP)+XOFF

! YMAX=largest value of Y which can occur for this X value
! YMIN=smallest value of Y which can occur for this X value
! ZMAX0=largest value of Z which can occur for this X value

         YMAX=S*SQRT(3._WP/16._WP)-X/SQRT(2._WP)
         YMIN=-S*SQRT(3._WP/64._WP)+X/SQRT(8._WP)
         ZMAX0=3._WP*S/8._WP-X*SQRT(3._WP/8._WP)
         DO J=JMIN,JMAX
            Y=REAL(J,KIND=WP)+YOFF
            IF(Y>=YMIN.AND.Y<=YMAX)THEN

! ZMAX=largest value of Z which can occur for this (X,Y)

               FY=(Y-YMIN)/(YMAX-YMIN)
               ZMAX=(1._WP-FY)*ZMAX0
               DO K=KMIN,KMAX
                  Z=REAL(K,KIND=WP)+ZOFF
                  IF(ABS(Z)<=ZMAX)THEN

! Site is occupied:

                     NAT=NAT+1
                     IF(NAT>MXNAT)THEN
                        CALL ERRMSG('FATAL','TARTET',' NAT.GT.MXNAT ')
                     ENDIF
                     IXYZ(NAT,1)=I
                     IXYZ(NAT,2)=J
                     IXYZ(NAT,3)=K
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO

! Initialize composition:
! Set X0 so that origin of TF is at centroid

      DO K=1,3
         X0(K)=0._WP
         DO I=1,NAT
            X0(K)=X0(K)+REAL(IXYZ(I,K))
            ICOMP(I,K)=1
         ENDDO
         X0(K)=-X0(K)/REAL(NAT)
      ENDDO

!*** Specify vectors A1 and A2 which will define target orientation
!    A1 is initially along x-axis
!    A2 is initially along y-axis

      A1(1)=1._WP
      A1(2)=0._WP
      A1(3)=0._WP
      A2(1)=0._WP
      A2(2)=1._WP
      A2(3)=0._WP

      WRITE(CDESCR,FMT='(A,I7,A)')' Tetrahedron of NAT=',NAT,' dipoles'
      IF(IOSHP>=0)THEN
         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')
         WRITE(IOSHP,FMT=9020)S,NAT,A1,A2,DX,X0
         DO JX=1,NAT
            WRITE(IOSHP,FMT=9030)JX,IXYZ(JX,1),IXYZ(JX,2),IXYZ(JX,3), &
                                 ICOMP(JX,1),ICOMP(JX,2),ICOMP(JX,3)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF
      RETURN
9020  FORMAT(' >TARTET  tetrahedral grain: S=',F9.4,/,             &
         I10,' = NAT ',/,                                          &
         3F10.6,' = A_1 vector',/,                                 &
         3F10.6,' = A_2 vector',/,                                 &
         3F10.6,' = lattice spacings (d_x,d_y,d_z)/d',/,           &
         3F10.5,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ', &
               'for dipole 0 0 0',/,                               &
         '     JA  IX  IY  IZ ICOMP(x,y,z)')
9030  FORMAT (I7,3I4,3I2)
    END SUBROUTINE TARTET
