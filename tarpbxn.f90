    SUBROUTINE TARPBXN(A1,A2,XD,YD,XS2,XS3,YS,ZS,DX,X0,CDESCR,IOSHP, &
                       MXNAT,NAT,IXYZ,ICOMP)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

! Arguments:

      REAL(WP) :: XD,XS2,XS3,YD,YS,ZS
      INTEGER :: IOSHP,MXNAT,NAT
      CHARACTER :: CDESCR*67
      INTEGER*2 :: &
         ICOMP(MXNAT,3)
      INTEGER :: &
         IXYZ(MXNAT,3)
      REAL(WP) :: &
         A1(3),   &
         A2(3),   &
         DX(3),   &
         X0(3)

! Local variables:

      INTEGER :: JA,JX,JXMAX,JXMAX1,JXMAX2,JXMIN,JXMIN1,JXMIN2,        &
                 JY,JYMAX,JYMIN,JZ,JZMAX,JZMIN,NATD,NATS1,NATS2,NX,NY,NZ
      REAL(WP) :: R2,Y,Y2,Z,Z2
      CHARACTER :: CMSGNM*70

! External subroutines:

      EXTERNAL ERRMSG

! Intrinsic functions:

!***********************************************************************
! subroutine TARPBXN
! Purpose: to construct "pillbox" target: disk on top of a bilayer slab
!          disk has composition 1
!          slab has compositions 2 and 3: 2 above, 3 below
!          disk-slab interface is in the y-z plane (x=0) (in TF)
!          disk symmetry axis is in x-direction (in TF)
! Input:
!        XD=(disk thickness)/d
!        YD=(disk diameter)/d
!        XS2=(x thickness of composition 2 of slab)/d   (d=lattice spacing)
!        XS3=(x thickness of composition 3 of slab)/d   
!        YS=(y length of slab)/d
!        ZS=(z length of slab)/d
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
!        X0(1-3)=offset vector, such that
!                x/d=IXYZ(J,1)+X0(1)
!                y/d=IXYZ(J,2)+X0(2)
!                z/d=IXYZ(J,3)+X0(3)
!                where (x,y,z)=location of dipole in Target Frame
!                with (0,0,0)=point where axis of disk
!                             intersects upper surface of slab
!                             upper surface of slab is assumed to
!                             located 0.5*d above dipoles defining slab

! Target consists of a rectangular slab of composition 1
!        extent XS x YS x ZS
!        center of upper surface of slab is at (0,0,0)
!        upper surface of slab is in the x=0 plane
!        lower surface of slab is in the x=-XS plane
!        slab runs from -XS   to 0
!                       -YS/2 to +YS/2
!                       -ZS/2 to +ZS/2
!        plus a disk of composition 2, diameter YD, thickness XD
!        with symmetry axis in x direction
!        centroid of disk located at (XD/2,0,0)

!        XS2,XS3,YS,ZS should ideally be integers
!        if YS is even, lattice will be located at y = +/-0.5, +/-1.5, .
!                 odd                              y = 0, +/-1, +/-2, ..
!        if ZS is even, lattice will be located at z = +/-0.5, +/-1.5, .
!                 odd                              z = 0, +/-1, +/-2, ..

! B.T.Draine, Princeton Univ. Obs.
! History:
! 06.12.08 (BTD) first written for collaboration with Ward Johnson and
!                Zhandos Utegulov (NIST, Boulder Colorado).
! 07.01.17 (BTD) modified to center pillbox on slab
! 07.01.18 (BTD) dimensions not quite right: fix
! 07.06.19 (BTD) major modification
!                * add X0(1-3) to argument list
!                * different lattice offset for YS even or odd
!                * different lattice offset for ZS even or odd
!                * revamp procedure for choosing dipoles
! 07.09.11 (BTD) changed IXYZ from INTEGER*2 to INTEGER
! 07.10.26 (BTD) major modification
!                * now allow one or two layers in slab
!                * changed argument list
! 07.10.30 (BTD) renamed XS2->XS3,XS1->XS2
!                reordered argument list
! 08.01.13 (BTD) cosmetic changes
! end history

! Copyright (C) 2006,2007,2008
!               B.T. Draine and P.J. Flatau

! This code is covered by the GNU General Public License.
!***********************************************************************
!***********************************************************************
      IF(XS2.GT.0._WP.AND.XS3.GT.0._WP)THEN
         WRITE(CMSGNM,FMT='(A)')' Two-layer slab+disk'
         WRITE(CDESCR,FMT='(A)')' Two-layer slab+disk'
      ELSE
         WRITE(CMSGNM,FMT='(A)')' Slab+disk'
         WRITE(CDESCR,FMT='(A)')' Slab+disk'
      ENDIF

! Specify target axes A1 and A2
! Convention: A1 will be (1,0,0) in Target Frame
!             A2 will be (0,1,0) in Target Frame

      DO JA=1,3
         A1(JA)=0._WP
         A2(JA)=0._WP
      ENDDO
      A1(1)=1._WP
      A2(2)=1._WP

! Target reference point =intersection of disk axis with upper surface
! of slab: (x,y,z)=(0,0,0) in Target Frame

! Set lattice offset in x direction:
! Dipole with IX=0 is located at x/d=0.5

      X0(1)=0.5_WP

! Determine lattice offset in y and z directions
!   if NY is even, dipole with IY=0 is located at y/d=0.5
!            odd                                      0
!   if NZ is even, dipole with IZ=0 is located at z/d=0.5
!            odd                                      0

      NY=NINT(YS)
      IF(MOD(NY,2)==0)THEN
         X0(2)=0.5_WP
      ELSE
         X0(2)=0._WP
      ENDIF

      NZ=NINT(ZS)
      IF(MOD(NZ,2)==0)THEN
         X0(3)=0.5_WP
      ELSE
         X0(3)=0._WP
      ENDIF

! Now populate lattice:

      NAT=0

! First do the pillbox (material 1)
! disk axis is in x direction
! XD=thickness
! YD=diameter

      NX=NINT(XD)

      JXMIN=0
      JXMAX=NX-1

      JYMIN=NINT(-0.5_WP*YD-X0(2))
      JYMAX=NINT(0.5_WP*YD-X0(2))
      JZMIN=NINT(-0.5_WP*YD-X0(3))
      JZMAX=NINT(0.5_WP*YD-X0(3))

      R2=(0.5_WP*YD)**2

!*** diagnostic
!      write(0,*)'disk: jx runs from ',jxmin,' to ',jxmax
!      write(0,*)'     diameter yd=',yd
!      write(0,*)'jymin,jymax=',jymin,jymax
!      write(0,*)'jzmin,jzmax=',jzmin,jzmax
!***
      DO JY=JYMIN,JYMAX
         Y=JY+X0(2)
         Y2=Y**2
         DO JZ=JZMIN,JZMAX
            Z=JZ+X0(3)
            Z2=Z**2
            IF(Y2+Z2<R2)THEN
!*** diagnostic
!               write(0,*)'disk: jy,jz=',jy,jz
!***
               DO JX=JXMIN,JXMAX
                  NAT=NAT+1
                  IF(NAT>MXNAT)THEN
                     CALL ERRMSG('FATAL','TARPBX',' NAT.GT.MXNAT')
                  ENDIF
                  IXYZ(NAT,1)=JX
                  IXYZ(NAT,2)=JY
                  IXYZ(NAT,3)=JZ
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      NATD=NAT
      DO JX=1,3
         DO JA=1,NATD
            ICOMP(JA,JX)=1
         ENDDO
      ENDDO

! completed definition of disk
! begin definition of slab
! layer 1: extends from     -xs1   to  0         [composition 2]
!          dipoles from    JXMIN1  to -1
! layer 2: extends from -(xs1+xs2) to -xs1       [composition 3]
!          dipoles from   JXMIN2   to (JXMIN1-1)

      JYMIN=NINT(-0.5_WP*YS-X0(2))
      JYMAX=NINT(0.5_WP*YS-X0(2))
      JZMIN=NINT(-0.5_WP*ZS-X0(3))
      JZMAX=NINT(0.5_WP*ZS-X0(3))

      JXMIN2=-NINT(XS2+XS3)
      JXMAX2=-NINT(XS2)-1
      JXMIN1=-NINT(XS2)
      JXMAX1=-1

! layer 1 = composition 2:
!*** diagnostic
!      write(0,*)'slab layer 1 runs from ',jxmin1,' to ',jxmax1
!***
      IF(JXMIN1.LE.JXMAX1)THEN
         DO JX=JXMIN1,JXMAX1
            DO JY=JYMIN,JYMAX
               IF(ABS(REAL(JY)+X0(2))<0.5_WP*YS)THEN
                  DO JZ=JZMIN,JZMAX
                     IF(ABS(REAL(JZ)+X0(3))<0.5_WP*ZS)THEN
                        NAT=NAT+1
                        IF(NAT>MXNAT)THEN
                           CALL ERRMSG('FATAL','TARPBX',' NAT.GT.MXNAT ')
                        ENDIF
                        IXYZ(NAT,1)=JX
                        IXYZ(NAT,2)=JY
                        IXYZ(NAT,3)=JZ
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         DO JX=1,3
            DO JA=NATD+1,NAT
               ICOMP(JA,JX)=2
            ENDDO
         ENDDO
      ENDIF
      NATS1=NAT-NATD

      IF(JXMIN2.LE.JXMAX2)THEN
!*** diagnnostic
!         write(0,*)'slab lower layer runs from ',jxmin2,' to ',jxmax2
!***
         DO JX=JXMIN2,JXMAX2
            DO JY=JYMIN,JYMAX
               IF(ABS(REAL(JY)+X0(2))<0.5_WP*YS)THEN
                  DO JZ=JZMIN,JZMAX
                     IF(ABS(REAL(JZ)+X0(3))<0.5_WP*ZS)THEN
                        NAT=NAT+1
                        IF(NAT>MXNAT)THEN
                           CALL ERRMSG('FATAL','TARPBX',' NAT.GT.MXNAT ')
                        ENDIF
                        IXYZ(NAT,1)=JX
                        IXYZ(NAT,2)=JY
                        IXYZ(NAT,3)=JZ
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         DO JX=1,3
            DO JA=NATD+NATS1+1,NAT
               ICOMP(JA,JX)=3
            ENDDO
         ENDDO
      ENDIF
      NATS2=NAT-NATS1

      IF(NATS1*NATS2>0)THEN
         WRITE(CMSGNM,FMT='(A,I7,A)')                              &
                       ' Bilayer slab+disk with NAT=',NAT,' dipoles'
      ELSE
         WRITE(CMSGNM,FMT='(A,I7,A)')' Slab+disk with NAT=',NAT,' dipoles'
      ENDIF
      IF(IOSHP>=0)THEN
         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')
         WRITE(IOSHP,FMT=9020)XS2,XS3,YS,ZS,NAT,A1,A2,DX,X0
         DO JA=1,NAT
            WRITE(IOSHP,FMT=9030)JA,IXYZ(JA,1),IXYZ(JA,2),IXYZ(JA,3), &
                                 ICOMP(JA,1),ICOMP(JA,2),ICOMP(JA,3)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF
!*** diagnostic/sanity check
!      jxmin=ixyz(1,1)
!      jxmax=ixyz(1,1)
!      do ja=2,nat
!         if(ixyz(ja,1)<jxmin)jxmin=ixyz(ja,1)
!         if(ixyz(ja,1)>jxmax)jxmax=ixyz(ja,1)
!      enddo
!      write(0,*)'ckpt 1 in tarpbxn: jxmin= ',jxmin,' jxmax= ',jxmax
!***
      RETURN
9020  FORMAT(' >TARPBXN : slab+disk; XS2,XS3,YS,ZS=',4F8.4,/,      &
         I10,' = NAT ',/,                                          &
         3F10.6,' = A_1 vector',/,                                 &
         3F10.6,' = A_2 vector',/,                                 &
         3F10.6,' = lattice spacings (d_x,d_y,d_z)/d',/,           &
         3F10.5,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ', &
               'for dipole 0 0 0',/,                               &
         '     JA  IX  IY  IZ ICOMP(x,y,z)')
9030  FORMAT(I7,3I4,3I2)
    END SUBROUTINE TARPBXN
