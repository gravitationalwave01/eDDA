    SUBROUTINE TARLYRSLAB(A1,A2,XV,YV,ZV,F1,F2,F3,F4,DX,X0,CDESCR,IOSHP, &
                          MXNAT,NAT,IXYZ,ICOMP)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE
! Arguments:
      REAL(WP) :: F1,F2,F3,F4,XV,YV,ZV
      INTEGER :: IOSHP,MXNAT,NAT
      CHARACTER :: CDESCR*67
      INTEGER*2 ::    &
         ICOMP(MXNAT,3)
      INTEGER ::     &
         IXYZ(MXNAT,3)
      REAL(WP) :: &
         A1(3),   &
         A2(3),   &
         DX(3),   &
         X0(3)

! Local variables:

      INTEGER :: IC,JA,JD,JX,JY,JZ,NCOMP,NX,NX1,NX2,NX3,NY,NZ
      CHARACTER :: CMSGNM*70

! External subroutines:

      EXTERNAL ERRMSG

! Intrinsic functions:

      INTEGER :: NINT
      INTRINSIC NINT

!***********************************************************************
! Routine to construct rectangular prism from "atoms"
! with layered structure
! slab has dimension XV*d in x direction
!                    YV*d   in y direction
!                    ZV*d   in z direction
!
! current version allows for up to 4 layers
! slab is layered in x direction
! top surface has x=0
! x0(1-3) = point at center of top surface

!      0       > x > -f1*a is composition 1
!   -f1*a      > x > -(f2+f1)*a is composition 2
! -(f2+f1)*a   > x > -(1-f4)*a is composition 3
!   -f4*a      > x >  -a   is composition 4

! for 1 layers, let f2=f3=f4=0
!     2             f3=f4=0
!     3             f4=0

! It is required that f1+f2+f3+f4=1.

! Input:
!        XV=(x-length)/d    (d=lattice spacing)
!        YV=(y-length)/d    (d=lattice spacing)
!        ZV=(z-length)/d    (d=lattice spacing)
!        F1,F2,F3,F4 = fraction of thickness contributed by compositions
!                      1,2,3,4
!        DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z
!                directions, and d=(dx*dy*dz)**(1/3)=effective lattice
!                spacing
!        IOSHP=device number for "target.out" file
!             =-1 to suppress printing of "target.out"
!        MXNAT=dimensioning information (max number of atoms)
! Output:
!        A1(1-3)=unit vector (1,0,0) defining target axis 1 in Target Fr
!        A2(1-3)=unit vector (0,1,0) defining target axis 2 in Target Fr
!        NAT=number of atoms in target
!        IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
!        ICOMP(1-NAT,1-3)=composition
!        X0(3)=location/d in TF of lattice site with IXYZ=(0,0,0)

! B.T.Draine, Princeton Univ. Obs.
! History:
! 07.02.23 (BTD): Created using TARSLB as starting point
! 07.09.11 (BTD): Changed IXYZ from INTEGER*2 to INTEGER
! 08.03.16 (BTD): Correct handling of NY>1, NZ>1
! 08.08.29 (BTD): Modified format 9020
! end history

! Copyright (C) 2007 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
! Sanity check:
      IF(ABS(F1+F2+F3+F4-1._WP)>1.E-6_WP)THEN
         WRITE(0,*)'F1,F2,F3,F4=',F1,F2,F3,F4
         CALL ERRMSG('FATAL','TARLYRPBC',' F1+F2+F3+F4.ne.1')
      ENDIF
      IF(F1<0.OR.F2<0.OR.F3<0.OR.F4<0._WP)THEN
         CALL ERRMSG('FATAL','TARLYRPBC',' one or more F are negative')
      ENDIF
      IF(F4==0._WP)NCOMP=3
      IF(F3==0._WP)NCOMP=2
      IF(F2==0._WP)NCOMP=1

      NX=NINT(XV/DX(1))
      NY=NINT(YV/DX(2))
      NZ=NINT(ZV/DX(3))

      WRITE(CMSGNM,FMT='(A,3I4)')' Layered slab; NX,NY,NZ=',NX,NY,NZ
      WRITE(CDESCR,FMT='(A,3I4)')' Layered slab; NX,NY,NZ=',NX,NY,NZ

      NX1=NINT(F1*NX)
      IF(NCOMP>=2)NX2=NINT((F1+F2)*NX)
      IF(NCOMP>=3)NX3=NINT((F1+F2+F3)*NX)

! Specify target axes A1 and A2
! Convention: A1 will be (1,0,0) in target frame
!             A2 will be (0,1,0) in target frame

      DO JA=1,3
         A1(JA)=0._WP
         A2(JA)=0._WP
      ENDDO

      A1(1)=1._WP
      A2(2)=1._WP

      NAT=NX*NY*NZ
      IF(NAT>MXNAT)THEN
         CALL ERRMSG('FATAL','TARSLB',' NAT.GT.MXNAT')
      ENDIF

! Now populate lattice:
! Top of lattice (in TF) consists of composition 1
! as we descend in x (in TF) we pass through compositions 2,3,4
! Top layer will have JX=-1
! Lowest layer will have JX=-NX

! Specify X0(3) = (x,y,z)/d in TF corresponding to IXYZ=(0,0,0)
! Top surface of slab (0.5d above topmost dipole layer) is assumed to have x=0
! Dipoles with JX=-1 are assumed to be at x=-dx/2
! Therefore X0(1)=0.5*DX(1)
! Set y=0 and z=0 to be at middle of slab
! JY runs from 1 to NY -> X0(2) = -(1+NY)/2
! JZ runs from 1 to NZ -> X0(3) = -(1+NZ)/2

      X0(1)=0.5*DX(1)
      X0(2)=-0.5*REAL(1+NY)
      X0(3)=-0.5*REAL(1+NZ)

      JA=0
      DO JX=-1,-NX,-1
         IF(JX.LE.NX1)IC=1
         IF(NCOMP>=2.AND.JX.GT.NX1.AND.JX.LE.NX2)IC=2
         IF(NCOMP>=3.AND.JX.GT.NX2.AND.JX.LE.NX3)IC=3
         IF(NCOMP>=4.AND.JX.GT.NX3)IC=4
         DO JY=1,NY
            DO JZ=1,NZ
               JA=JA+1
               IXYZ(JA,1)=JX
               IXYZ(JA,2)=JY
               IXYZ(JA,3)=JZ
               DO JD=1,3
                  ICOMP(JA,JD)=IC
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      WRITE(CMSGNM,FMT='(A,I7,A)') '  Total slab thickness =',NAT,' dipoles'
      IF(IOSHP>=0)THEN
         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')
         WRITE(IOSHP,FMT=9020)XV,YV,ZV,NAT,A1,A2,DX,X0
         DO JA=1,NAT
            WRITE(IOSHP,FMT=9030)JA,IXYZ(JA,1),IXYZ(JA,2),IXYZ(JA,3), &
                                 ICOMP(JA,1),ICOMP(JA,2),ICOMP(JA,3)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF

      RETURN
9020  FORMAT (' >TARREC   rectangular prism; AX,AY,AZ=',3F8.4,/,   &
         I10,' = NAT ',/,                                          &
         3F10.6,' = A_1 vector',/,                                 &
         3F10.6,' = A_2 vector',/,                                 &
         3F10.6,' = lattice spacings (d_x,d_y,d_z)/d',/,           &
         3F10.5,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ', &
               'for dipole 0 0 0',/,                               &
         '     JA  IX  IY  IZ ICOMP(x,y,z)')

9030  FORMAT (I7,3I4,3I2)
    END SUBROUTINE TARLYRSLAB
