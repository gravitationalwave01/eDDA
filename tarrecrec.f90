    SUBROUTINE TARRECREC(A1,A2,XV1,YV1,ZV1,XV2,YV2,ZV2,DX,X0, &
                         CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

! Arguments:

      REAL(WP) :: XV1,YV1,ZV1,XV2,YV2,ZV2
      INTEGER :: IOSHP,MXNAT,NAT
      CHARACTER :: CDESCR*67
      INTEGER*2 :: ICOMP(MXNAT,3)
      INTEGER ::     &
         IXYZ(MXNAT,3)
      REAL(WP) :: &
         A1(3),   &
         A2(3),   &
         DX(3),   &
         X0(3)

! Local variables:

      INTEGER :: JA,JX,JX1MAX,JX1MIN,JX2MAX,JX2MIN,JY,JY1MAX,JY1MIN,JY2MAX, &
         JY2MIN,JZ,JZ1MAX,JZ1MIN,JZ2MAX,JZ2MIN,NX1,NX2,NY1,NY2,NZ1,NZ2
      CHARACTER :: CMSGNM*70

! External subroutines:

      EXTERNAL ERRMSG

! Intrinsic functions:

!***********************************************************************
! Routine to generate target consisting of rectangular block of
!        composition 1 resting on top of rectangular block of composition 2
!        TF origin = center of top surface of block 1
!                    i.e. 0.5d above the top layer of dipoles
! Input:
!        XV1=(x-length)/d  of block 1    (d=lattice spacing)
!        YV1=(y-length)/d           1
!        ZV1=(z-length)/d           1
!        XV2=(x-length)/d  of block 2
!        YV2=(y-length)/d  of block 2
!        ZV2=(z-length)/d  of block 2
!
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
!        X0(3)=(x,y,z)_TF corresponding to lattice site IXYZ=(0,0,0)
!
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
! 07.09.11 (BTD): changed IXYZ from INTEGER*2 to INTEGER
! 07.10.27 (BTD): created starting from TARREC as template
! 08.08.30 (BTD): modified format 9020
! 12.04.02 (BTD): v2
!                 corrected typo noted by Georges Levi 
!                 (Universite Paris-Diderot)
! 12.06.13 (BTD): v3
!                 corrected two typos noted by Georges Levi
!                 (Universite Paris-Diderot)
! end history

! Copyright (C) 1993,1995,1996,1997,1998,2000,2007,2008,2012
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
      NX1=NINT(XV1/DX(1))
      NY1=NINT(YV1/DX(2))
      NZ1=NINT(ZV1/DX(3))
      NX2=NINT((XV1+XV2)/DX(1))
      NX2=NX2-NX1
      NY2=NINT(YV2/DX(2))
      NZ2=NINT(ZV2/DX(3))

      WRITE(CMSGNM,FMT='(A,3I4,A,3I4)')             &
            'Two rect.slabs NX1,NY1,NZ1=', &
            NX1,NY1,NZ1,' NX2,NYZ,NZ2=',NX2,NY2,NZ2
      WRITE(CDESCR,FMT='(A,3I4,A,3I4)')             &
            'Two rect.slabs NX1,NY1,NZ1=', &
            NX1,NY1,NZ1,' NX2,NYZ,NZ2=',NX2,NY2,NZ2

! Specify target axes A1 and A2
! Convention: A1 will be (1,0,0) in target frame
!             A2 will be (0,1,0) in target frame

      DO JA=1,3
         A1(JA)=0._WP
         A2(JA)=0._WP
      END DO
      A1(1)=1._WP
      A2(2)=1._WP

! Bottom of upper box and top of lower box will be in x_TF=0 plane
!                 JX=0 corresponds to x_TF=dx/2
!     upper box:  JX runs from 0 to NX1-1 for material 1
!                 example: NX1=1, JX runs from 0 to 1
!     lower box:  JX runs from -NX2 to -1 for material 2
!                 example: NX2=2, JX runs from -2 to -1

! If NY1 is even, JY=0 corresponds to y_TF = -dy/2
!     upper box:  JY runs from 1-(NY1/2) to (NY1/2)
!                 example: NY1=2, JY runs from 0 to +1
!     lower box:  JY runs from 1-int(NY2/2) to NY2-int(NY2/2)
!                 example1: NY2=2, JY runs from 0 to 1
!                 example2: NY2=3, JY runs from 0 to 2

!    NY1 is  odd, JY=0 corresponds to y_TF=0
!     upper box:  JY runs from -(NY1-1)/2 to (NY1-1)/2
!                 example: NY1=3, JY runs from -1 to +1
!     lower box:  JY runs from -int[(NY2-1)/2] to NY2-1-int[(NY2-1)/2]

! If NZ1 is even, JZ=0 corresponds to z_TF = -dz/2
!     upper box:  JZ runs from 1-(NZ1/2) to (NZ1/2)
!     lower box:  JZ runs from 1-int(NZ2/2) to NZ2-int(NZ2/2)

!    NZ1 is odd,  JZ=0 corresponds to z_TF=0
!     upper box:  JZ runs from -(NZ1-1)/2 to (NZ1-1)/2
!     lower box:  JZ runs from -int[(NZ2-1)/2] to NZ2-1-int[(NZ2-1)/2]

      X0(1)=0.5_WP
      JX1MIN=0
      JX1MAX=NX1-1
      JX2MIN=-NX2
      JX2MAX=-1
      IF(MOD(NY1,2)==0)THEN ! NY1 is even
         X0(2)=-0.5_WP
         JY1MIN=1-NY1/2
         JY1MAX=NY1/2
         JY2MIN=1-INT(NY2/2)
         JY2MAX=NY2-INT(NY2/2)
      ELSE                  ! NY1 is odd
         X0(2)=0.0_WP
         JY1MIN=-(NY1-1)/2
         JY1MAX=(NY1-1)/2
         JY2MIN=-INT((NY2-1)/2)
         JY2MAX=NY2-1-INT((NY2-1)/2)
      ENDIF

      IF(MOD(NZ1,2)==0)THEN ! NZ1 is even
! 12.04.02 correction
!         X0(3)=0.5_WP
         X0(3)=-0.5_WP
         JZ1MIN=1-NZ1/2
         JZ1MAX=NZ1/2
         JZ2MIN=1-INT(NZ2/2)
         JZ2MAX=NZ2-INT(NZ2/2)
      ELSE                  ! NZ2 is odd
         X0(3)=0.0_WP
         JZ1MIN=-(NZ1-1)/2
         JZ1MAX=(NZ1-1)/2
         JZ2MIN=-INT(NZ2-1)/2
         JZ2MAX=NZ2-1-INT((NZ2-1)/2)
      ENDIF

! Now populate upper rectangular slab:

      NAT=0
      DO JZ=JZ1MIN,JZ1MAX
         DO JY=JY1MIN,JY1MAX
            DO JX=JX1MIN,JX1MAX
               NAT=NAT+1
               IF(NAT>MXNAT)THEN
                  CALL ERRMSG('FATAL','TARREC',' NAT.GT.MXNAT')
               ENDIF
               IXYZ(NAT,1)=JX
               IXYZ(NAT,2)=JY
               IXYZ(NAT,3)=JZ

! homogeneous, isotropic composition:

               ICOMP(NAT,1)=1
               ICOMP(NAT,2)=1
               ICOMP(NAT,3)=1
            ENDDO
         ENDDO
      ENDDO

! Lower rectangular slab

      DO JZ=JZ2MIN,JZ2MAX
         DO JY=JY2MIN,JY2MAX
            DO JX=JX2MIN,JX2MAX
               NAT=NAT+1
               IF(NAT>MXNAT)THEN
                  CALL ERRMSG('FATAL','TARREC',' NAT.GT.MXNAT')
               ENDIF
               IXYZ(NAT,1)=JX
               IXYZ(NAT,2)=JY
               IXYZ(NAT,3)=JZ

! homogeneous, isotropic composition:

               ICOMP(NAT,1)=2
               ICOMP(NAT,2)=2
               ICOMP(NAT,3)=2
            ENDDO
         ENDDO
      ENDDO
      WRITE(CMSGNM,FMT='(A,I7,A)')'Two rectangular solids with NAT=',NAT, &
                                  ' dipoles'
      IF(IOSHP>=0)THEN
         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')
         WRITE(IOSHP,FMT=9020)NX1,NY1,NZ1,NX2,NY2,NZ2,NAT,A1,A2,DX,X0
         DO JA=1,NAT
            WRITE(IOSHP,FMT=9030)JA,IXYZ(JA,1),IXYZ(JA,2),IXYZ(JA,3), &
                                 ICOMP(JA,1),ICOMP(JA,2),ICOMP(JA,3)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF
      RETURN
9020  FORMAT(' >TARRECREC: two rectangular solids:',                 &
         ' NX1,NY1,NZ1=',3I4,' NX2,NY2,NZ2=',3I4,/,                  &
         I10,' = NAT=',/,                                            &
         3F10.6,' = A_1 vector',/,                                   &
         3F10.6,' = A_2 vector',/,                                   &
         3F10.5,' = lattice spacings (d_x,d_y,d_z)/d',/,             &
         3F10.5,' = x0(1-3)= (x_TF,y_TF,z_TF)/d for dipole 0 0 0',/, &
         '     JA  IX  IY  IZ ICOMP(x,y,z)')
9030  FORMAT(I7,3I4,3I2)
    END SUBROUTINE TARRECREC
