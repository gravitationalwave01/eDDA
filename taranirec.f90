    SUBROUTINE TARANIREC(A1,A2,XV,YV,ZV,DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

! Arguments:

      REAL(WP) :: XV,YV,ZV
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

      INTEGER :: JA,JX,JY,JZ,NX,NY,NZ
      CHARACTER :: CMSGNM*70

! External subroutines:

      EXTERNAL ERRMSG

! Intrinsic functions:
! intrinsic function: do not need to explicitly type:      INTEGER IFIX
!***********************************************************************
! Subroutine TARANIREC
! Purpose: to construct rectangular target with anisotropic dielectric
!          tensor.
!          Dielectric tensor is assumed to be diagonalized in target
!          frame, with
!          material 1 for x-direction
!          material 2 for y-direction
!          material 3 for z-direction
!          Note: Must specify NCOMP=3 in ddscat.par, and provide
!                file names for three materials
! Input:
!        XV=(x-length)/d    (d=lattice spacing)
!        YV=(y-length)/d
!        ZV=(z-length)/d
!        DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z
!                directions, and d=(dx*dy*dz)**(1/3)=effective lattice
!                spacing
!        IOSHP=device number for "target.out" file
!             =-1 to suppress printing of "target.out"
!        MXNAT=dimensioning information (max number of atoms)
! Output:
!        A1(1-3)=unit vector (1,0,0) defining target axis 1 in Target Frame
!        A2(1-3)=unit vector (0,1,0) defining target axis 2 in Target Frame
!        NAT=number of atoms in target
!        IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
!        X0(1-3)=(location/d) in TF of dipole with IXYZ=(0,0,0).  This
!                will henceforth be treated as the origin of physical
!                coordinates in the TF.
!                origin is assumed to be located at midpoint of upper surface
!                normal to the A1 vector.
!                Upper surface is presumed to be 0.5*d above uppermost
!                layer of dipoles.

! B.T.Draine, Princeton Univ. Obs.
! History:
! 91.01.05 (BTD): Change I4 -> I7 when printing NAT.
! 91.04.22 (BTD): Added XV,YV,ZV and A1,A2 to WRITE(IOSHP,FMT=9020)
! 93.03.12 (BTD): Specified CDESCR*67, and now use it.
! 95.07.12 (BTD): Changed definition of axes A1,A2.  They used to be
!                 longest,next longest dimention.  They are now
!                 fixed to always be (1,0,0) and (0,1,0).
! 96.01.26 (BTD): Replaced IX,IY,IZ by IXYZ
! 97.12.26 (BTD): Added DX(1-3) to argument list to support nonuniform
!                 lattice spacing.
! 98.02.11 (BTD): Modify to generate target on noncubic lattice.
! 00.11.02 (BTD): Add ICOMP to argument list, set ICOMP=1,
!                 write ICOMP to target.out
! 04.04.29 (BTD): Created from TARREC with only minor modifications
! 07.09.11 (BTD): Changed IXYZ from INTEGER*2 to INTEGER
! 08.01.13 (BTD): Cosmetic changes
!                 Added X0 to argument list
!                 Added code to define X0
! 08.08.29 (BTD): Modified format 9020
! 12.05.14 (BTD): further modified format 9020 to provide enough space
!                 for x0(1-3) to be written out.
! end history

! Copyright (C) 1993,1995,1996,1997,1998,2000,2004,2007,2008,2010
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
!***********************************************************************
      NX=INT(XV/DX(1)+0.5_WP)
      NY=INT(YV/DX(2)+0.5_WP)
      NZ=INT(ZV/DX(3)+0.5_WP)
      WRITE(CMSGNM,FMT='(A,3I4)')' Anisotropic rectangular prism; NX,NY,NZ=', &
                                 NX,NY,NZ
      WRITE(CDESCR,FMT='(A,3I4)')' Anisotropic rectangular prism; NX,NY,NZ=', &
                                 NX,NY,NZ

! Specify target axes A1 and A2
! Convention: A1 will be (1,0,0) in target frame
!             A2 will be (0,1,0) in target frame

      DO JA=1,3
         A1(JA)=0._WP
         A2(JA)=0._WP
      ENDDO
      A1(1)=1._WP
      A2(2)=1._WP

! Now populate lattice:
      NAT=0
      DO JZ=1,NZ
         DO JY=1,NY
            DO JX=1,NX
               NAT=NAT+1
               IF(NAT>MXNAT) THEN
                  CALL ERRMSG('FATAL','TARREC',' NAT.GT.MXNAT')
               ENDIF
               IXYZ(NAT,1)=JX
               IXYZ(NAT,2)=JY
               IXYZ(NAT,3)=JZ
            ENDDO
         ENDDO
      ENDDO
      IF(NAT>MXNAT)THEN
         CALL ERRMSG('FATAL','TARREC',' NAT.GT.MXNAT ')
      ENDIF

! set TF origin of coordinates to be located at midpoint of
! upper surface normal to A1.  Upper surface is presumed to be located
! 0.5*d above uppermost dipole layer

      X0(1)=-REAL(NX)-0.5
      X0(2)=-REAL(1+NY)/2._WP
      X0(3)=-REAL(1+NZ)/2._WP

! homogeneous, anisotropic composition:

      DO JX=1,3
         DO JA=1,NAT
            ICOMP(JA,JX)=JX
         ENDDO
      ENDDO

      WRITE(CMSGNM,FMT='(A, I7,A)') &
            '  Anisotropic rectangular prism of NAT=',NAT,' dipoles'
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
9020  FORMAT (' >TARREC anisotropic rectangular prism; AX,AY,AZ=',3F9.4,/, &
         I9,' = NAT',/,                                                    &
         3F10.6,' = A_1 vector',/,                                         &
         3F10.6,' = A_2 vector',/,                                         &
         3F10.6,' = lattice spacings (d_x,d_y,d_z)/d',/,                   &
         3F10.4,' = lattice offset x0(1-3)=(x_TF,y_TF,z_TF)/d ',           &
               'for dipole 0 0 0',/,                                       &
         '     JA   IX   IY   IZ ICOMP(x,y,z)')
9030  FORMAT (I7,3I5,3I2)
    END SUBROUTINE TARANIREC
