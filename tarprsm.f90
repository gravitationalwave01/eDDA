    SUBROUTINE TARPRSM(A1,A2,A,BA,CA,LA,DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ, &
                       ICOMP)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

! Arguments:

      CHARACTER :: CDESCR*67
      INTEGER :: IOSHP, MXNAT, NAT
      INTEGER*2 :: ICOMP(MXNAT,3)
      INTEGER :: IXYZ(MXNAT,3)
      REAL(WP) :: A,BA,CA,LA

! Array arguments:

      REAL(WP) :: A1(3),A2(3),DX(3),X0(3)

! Local scalars:

      CHARACTER :: CMSGNM*70
      INTEGER :: JA,JX,JY,JZ,NFAC,NLAY,NX,NY,NZ

      REAL(WP) :: B,BETA,C,COTBETA,COTGAMMA,GAMMA,L,XMAX,XMIN,Y, &
                  YMAX,Z,ZMAX,ZMIN

!***********************************************************************
! Routine to construct triangular prism from "atoms".
! triangle side lengths = a,b,c ; prism length = L
! prism axis is assumed to be in x direction
! normal to prism face of width a : (0,1,0)
! normal to prism face of width b : (0,-cos(gamma),sin(gamma))
! normal to prism face of width c : (0,-cos(beta),-sin(beta))
! angles alpha,beta,gamma are opposite sides a,b,c

! first layer of target array is at x=0 (target boundary at x=-0.5)
! layer along side a is at int(b*sin(gamma))
!      boundary is at ymax=int(b*sin(gamma))+0.5
! Input:
!        A      = a/d
!        BA     = b/a
!        CA     = c/a
!        LA     = L/a
!        DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z
!                directions, and d=(dx*dy*dz)**(1/3)=effective lattice
!                spacing
!        IOSHP  =device number for "target.out" file
!               =-1 to suppress printing of "target.out"
!        MXNAT  =dimensioning information (max number of atoms)
! Output:
!        A1(1-3)=unit vector (1,0,0) defining target axis 1 (prism axis)
!        A2(1-3)=unit vector (0,1,0) defining target axis 2 (normal to
!                rectangular faces with sides a,L
!        NAT=number of atoms in target
!        IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
!        ICOMP(1-NAT,1-3)=dielectric function identifier for dipole
!                         locations and 3 directions in space
!        CDESCR=string describing target (up to 67 chars.)
!        X0(3)=location/d in TF of dipole with IXYZ=0 0 0
!              we set the TF origin to be the centroid of the prism
! B.T.Draine, Princeton Univ. Obs., 2002.02.12

! History:
! 02.02.12 (BTD): adapted from tarprsm.f for ver5a10
! 04.02.25 (BTD): corrected typographical error found by In-Ho Lee
! 06.09.30 (BTD): remove unused variable YMIN
! 07.09.11 (BTD): changed IXYZ from INTEGER*2 to INTEGER
! 08.08.30 (BTD): added X0 to argument list
!                 add code to evaluate X0(1-3)
!                 modified format 9020
! end history

! Copyright (C) 2002,2004,2006,2007,2008 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!-----------------------------------------------------------------------

! Current version of TARPRSM is restricted to cubic lattices

      IF(DX(1)/=1._WP.OR.DX(2)/=1._WP)THEN
         CALL ERRMSG('FATAL','TARPRSM',                          &
                     ' tarprsm does not support noncubic lattice')
      ENDIF

      DO JA=1,3
         A1(JA)=0._WP
         A2(JA)=0._WP
      END DO
      A1(1)=1._WP
      A2(2)=1._WP

! A = a/d
! B = b/d
! C = c/d
! L = L/d

      B=BA*A
      C=CA*A
      L=LA*A

! sanity check:

      IF((A>=B+C).OR.(B>=A+C).OR.(C>=A+B))CALL ERRMSG('FATAL','TARPRSM',   &
         'prism sides a,b,c do not satisfy triangle inequality:check input')

      BETA=ACOS((A*A+C*C-B*B)/(2._WP*A*C))
      GAMMA=ACOS((A*A+B*B-C*C)/(2._WP*A*B))
      COTBETA=COS(BETA)/SIN(BETA)
      COTGAMMA=COS(GAMMA)/SIN(GAMMA)

! ideal prism extent:
!    x = xmin to x=xmax
!    triangle vertices at (0,ymax,zmin),(0,ymax,zmax),(0,ymin,za)
!    where xmin=-0.5
!          xmax=xmin+L
!          ymax=int[b*sin(gamma)]+0.5
!          ymin=ymax-b*sin(gamma)
!          zmax=a/2
!          zmin=-a/2
!          za=-a/2+c*cos(beta)

      XMIN=-0.5_WP
      XMAX=XMIN+L
      YMAX=INT(B*SIN(GAMMA))+0.5_WP
      ZMAX=0.5_WP*A
      ZMIN=-ZMAX

! Now determine limits for testing x,y,z values
! Along axis (x), run from 0 to NX=INT(XMAX)
! Along y direction, run from 0 to NY=INT(YMAX)
! Along z direction, run from -NZ to NZ=INT(ZMAX)

      NX=INT(XMAX)
      NY=INT(YMAX)
      NZ=INT(ZMAX)

      NAT=0
      DO JZ=-NZ,NZ
         Z=REAL(JZ,KIND=WP)
         DO JY=0,NY
            Y=REAL(JY,KIND=WP)
            IF((Z<=ZMAX+(Y-YMAX)*COTGAMMA).AND.(Z>=ZMIN-(Y-YMAX)*COTBETA))THEN
               DO JX=0,NX
                  NAT=NAT + 1
                  IXYZ(NAT,1)=JX
                  IXYZ(NAT,2)=JY
                  IXYZ(NAT,3)=JZ
               ENDDO
            ENDIF
         ENDDO
      ENDDO

      NLAY=NX+1
      NFAC=NAT/NLAY

      DO JX=1,3
         X0(JX)=0._WP
         DO JA=1,NAT
            X0(JX)=X0(JX)+REAL(IXYZ(JA,JX))
            ICOMP(JA,JX)=1
         ENDDO
         X0(JX)=-X0(JX)/REAL(NAT)
      ENDDO

! Note: string CDESCR will be printed by subr. TARGET

      WRITE(CDESCR,FMT='(A,I7,A)')' Hexagonal prism of NAT=',NAT,' dipoles'

! Here print any additional target info which is desired by using
! subroutine WRIMSG

      WRITE(CMSGNM,FMT='(A,3F7.4,A)')' A1 = (',A1,') = prism axis direction'
      CALL WRIMSG('TARPRSM',CMSGNM)
      WRITE(CMSGNM,FMT='(A,3F7.4,A)')' A2 = (',A2,') = normal to face of width a'
      CALL WRIMSG('TARPRSM',CMSGNM)
      WRITE(CMSGNM,FMT='(A,I7,A,I4,A,I4)')' NAT=',NAT, ' NFAC=',NFAC, &
                                          ' NLAY=',NLAY
      CALL WRIMSG('TARPRSM',CMSGNM)
      WRITE(CDESCR,FMT='(A,F7.4,A,F8.4,A,F7.4,A,F8.4)')' tri.prism: a/d=', &
                                          A,' b/a=',BA,' c/a=',CA,' L/a=',LA

! Now store atom coordinates in file

      IF(IOSHP>0)THEN
         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')

!*** 3 lines of general description allowed:

         WRITE(IOSHP,9020)A,BA,CA,LA,NAT,A1,A2,DX,X0

!***

         DO JX=1,NAT
            WRITE(IOSHP,FMT=9030)JX,IXYZ(JX,1),IXYZ(JX,2),IXYZ(JX,3), &
                                 ICOMP(JX,1),ICOMP(JX,2),ICOMP(JX,3)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF
      RETURN
9020  FORMAT(' Triangular Prism: a/d=',F7.4,' b/a=',F7.4,          &
             ' c/a=',F7.4,' L/a=',F7.4,/,                          &
         I10,' = NAT ',/,                                          &
         3F10.6,' = A_1 vector',/,                                 &
         3F10.6,' = A_2 vector',/,                                 &
         3F10.6,' = lattice spacings (d_x,d_y,d_z)/d',/,           &
         3F10.5,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ', &
               'for dipole 0 0 0',/,                               &
         '     JA  IX  IY  IZ ICOMP(x,y,z)')
9030  FORMAT(I7,3I4,3I2)
    END SUBROUTINE TARPRSM
