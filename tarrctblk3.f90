    SUBROUTINE TARRCTBLK3(A1,A2,XS1,YS1,ZS1,XS2,YS2,ZS2,XS3,YS3,ZS3,DX, &
                          X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

! Arguments:

      REAL(WP) :: XS1,XS2,XS3,YS1,YS2,YS3,ZS1,ZS2,ZS3
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

      INTEGER :: JA,JD,JX,JXMAX1,JXMAX2,JXMAX3,JXMIN1,JXMIN2,JXMIN3,  &
                 JY,JYMAX1,JYMAX2,JYMAX3,JYMIN1,JYMIN2,JYMIN3,        &
                 JZ,JZMAX1,JZMAX2,JZMAX3,JZMIN1,JZMIN2,JZMIN3,        &
                 NATS1,NATS2,NATS3,NX1,NX2,NX3,NY1,NY2,NY3,NZ1,NZ2,NZ3
      CHARACTER :: CMSGNM*70
      REAL(WP) :: DELTA,TERM

! External subroutines:

      EXTERNAL ERRMSG

! Intrinsic functions:

!***********************************************************************
! subroutine TARRCTBLK3
! Purpose: to construct target consisting of 3 stacked rectangular blocks:
! boundary between blocks 1 and 2 lies in x=0 plane

!      block 1, of composition 1, extends from x/d= 0 to XS1
!                                              y/d=-YS1/2 to +YS1/2
!                                              z/d=-ZS1/2 to +ZS1/2
!      block 2, of composition 2, extends from x/d=-XS2 to     0
!                                              y/d=-YS2/2 to +YS2/2
!                                              z/d=-ZS2/2 to +ZS2/2
!      block 3, of composition 3, extends from x/d=-(XS2+XS3) to -XS2
!                                              y/d=-YS3/2 to +YS3/2
!                                              z/d=-ZS3/2 to +ZS3/2

! Input:
!        XS1=(x thickness of composition 1 of slab)/d   (d=lattice spacing)
!        YS1=(y extent of composition 1 of slab)/d   (d=lattice spacing)
!        ZS1=(z extent of composition 1 of slab)/d   (d=lattice spacing)
!        XS2=(x thickness of composition 2 of slab)/d
!        YS2=(y extent of composition 2 of slab)/d
!        ZS2=(z extent of composition 2 of slab)/d
!        XS3=(x thickness of composition 3 of slab)/d
!        YS3=(y extent of composition 3 of slab)/d
!        ZS3=(z extent of composition 3 of slab)/d
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
!        ICOMP(1-NAT,1-3)=composition number for locations 1-NAT 
!                         and directions 1-3
!        X0(1-3)=offset vector, such that
!                x/d=(IXYZ(J,1)+X0(1))*DX(1)
!                y/d=(IXYZ(J,2)+X0(2))*DX(2)
!                z/d=(IXYZ(J,3)+X0(3))*DX(3)
!                where (x,y,z)=location of dipole in Target Frame
!                with (0,0,0)=point where axis of disk
!                             intersects upper surface of slab

! Target consists of three stacked rectangular slabs with x-axis running
!        through center of each slab
!        slab 1 is of extent XS1 x YS1 x ZS1 (in units of d)
!        slab 2 is of extent XS2 x YS2 x ZS2 (in units of d)
!        slab 3 is of extent XS3 x YS3 x ZS3 (in units of d)

!        center of upper surface of slab 1 is at (XS1*d*DX(1), 0 , 0)
!        center of upper surface of slab 2 is at (0,0,0)
!        center of lower surface of slab 2 = upper surface of slab 3
!                                          is at (-XS2*d*DX(2) , 0 , 0 )
!        boundary between compositions 1 and 2 is in the x=0 plane
!        boundary between compositions 2 and 3 is in the x=-XS2*d*DX(1) plane
!        lower surface of slab 3 is in the x=-(XS2+XS3)*d*DX(1) plane

!        if max(YS1,YS2,YS3)/DX(2)
!           is even: lattice will be located at y/(d*DX(2)) = +/-0.5, +/-1.5, .
!               odd:                            y/(d*DX(2)) = 0, +/-1, +/-2, ..
!        if max(ZS1,ZS2,ZS3)/DX(3) 
!           is even: lattice will be located at z/(d*DX(3)) = +/-0.5, +/-1.5, .
!               odd:                            z/(d*DX(3)) = 0, +/-1, +/-2, ..

! B.T.Draine, Princeton Univ. Obs.
! History:
! 08.01.31 (BTD) created
! 08.03.23 (BTD) fixed errors
! 08.03.30 (BTD) replace 1.E-8 -> 1.E-5 so that it will work with
!                single as well as double precision arithmetic
! 08.03.30 (BTD) replace 1.E-5 by DELTA=0.499999 to obtain proper
!                behavior as YS1,YS2,YS3 vary, etc.
! 08.08.30 (BTD) modified format 9020
! 10.01.24 (BTD) ver7.0.8
!                ver7.1.0
!                reordered argument list
!                XS1,XS2,XS3,YS1,YS2,YS3,ZS1,ZS2,ZS3 ->
!                XS1,YS1,ZS1,XS2,YS2,ZS2,XS3,YS3,ZS3
!                for consistency with UserGuide
!                problem reported by
!                Bala Krishna Juluri, Penn State Univ.
! end history

! Copyright (C) 2008,2010
!               B.T. Draine and P.J. Flatau

! This code is covered by the GNU General Public License.
!***********************************************************************
!***********************************************************************
      IF(XS1*XS2*XS3.GT.0._WP)THEN
         WRITE(CMSGNM,FMT='(A)')' Three stacked rectangular slabs'
         WRITE(CDESCR,FMT='(A)')' Three stacked rectangular slabs'
      ELSEIF(XS1==0..OR.XS2==0..OR.XS3==0.)THEN
         WRITE(CMSGNM,FMT='(A)')' Two stacked rectangular slabs'
         WRITE(CDESCR,FMT='(A)')' Two stacked rectangular slabs'
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

! Target reference point =intersection of x_TF axis with plane separating
! slabs 1 and 2

! Set lattice offset in x direction:
! Dipole with IX=0 is located at x/(d*DX(1))=0.5

      X0(1)=0.5_WP

! Determine lattice offset in y and z directions
!   if NY is even, dipole with IY=0 is located at y/d=0.5
!            odd                                      0
!   if NZ is even, dipole with IZ=0 is located at z/d=0.5
!            odd                                      0

      TERM=MAX(YS1,YS2)
      TERM=MAX(TERM,YS3)
      NY1=NINT(TERM/DX(2))
      IF(MOD(NY1,2)==0)THEN
         X0(2)=0.5_WP
      ELSE
         X0(2)=0._WP
      ENDIF

      TERM=MAX(ZS1,ZS2)
      TERM=MAX(TERM,ZS3)
      NZ1=NINT(TERM/DX(3))
      IF(MOD(NZ1,2)==0)THEN
         X0(3)=0.5_WP
      ELSE
         X0(3)=0._WP
      ENDIF

! Now populate lattice:

      NAT=0
      DELTA=0.49999_WP

! First do the top block (material 1)

      JXMIN1=0
      JXMAX1=NINT(XS1/DX(1))-1
      NX1=JXMAX1-JXMIN1+1

! number of dipole layers 
!           NX1 = jxmax1-jxmin1+1 = nint(xs1/dx(1))

      JYMIN1=NINT(-0.5_WP*YS1/DX(2)-X0(2)+DELTA)
      JYMAX1=NINT(0.5_WP*YS1/DX(2)-X0(2)-DELTA)
      JZMIN1=NINT(-0.5_WP*ZS1/DX(3)-X0(3)+DELTA)
      JZMAX1=NINT(0.5_WP*ZS1/DX(3)-X0(3)-DELTA)
!*** diagnostic
!      write(0,*)'tarrctblk3 ckpt 1'
!      write(0,*)'ys1=',ys1,' x0(2)=',x0(2)
!      write(0,*)'jymin1=',jymin1,' jymax1=',jymax1
!      write(0,*)'zs1=',zs1,' x0(3)=',x0(3)
!      write(0,*)'jzmin1=',jzmin1,' jzmax1=',jzmax1
!***
      NY1=JYMAX1-JYMIN1+1
      NZ1=JZMAX1-JZMIN1+1

! number of dipoles in y direction 
!     NY1 = jymax1-jymin1+1 
!         = nint(0.5*ys1/dx(2)-x0(2)+delta)-nint(-0.5*ys1/dx(2)-x0(2)-delta)+1
!         = nint(0.5*ys1/dx(2)-x0(2)+delta)+nint(0.5*ys1/dx(2)+x0(2)+delta)+1
!
! number of dipoles in Z direction 
!     NZ1 = jzmax1-jzmin1+1 
!         = nint(0.5*zs1/dx(3)-x0(3)+delta)-nint(-0.5*zs1/dx(3)-x0(3)-delta)+1
!         = nint(0.5*zs1/dx(3)-x0(3)+delta)+nint(0.5*zs1/dx(3)+x0(3)+delta)+1

      IF(JXMIN1.LE.JXMAX1)THEN
         DO JY=JYMIN1,JYMAX1
            DO JZ=JZMIN1,JZMAX1
               DO JX=JXMIN1,JXMAX1
                  NAT=NAT+1
                  IF(NAT>MXNAT)THEN
                     CALL ERRMSG('FATAL','TARRCTBLK3',' NAT.GT.MXNAT')
                  ENDIF
                  IXYZ(NAT,1)=JX
                  IXYZ(NAT,2)=JY
                  IXYZ(NAT,3)=JZ
                  DO JD=1,3
                     ICOMP(NAT,JD)=1
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         NATS1=NAT
      ENDIF
      NATS1=NAT

!*** diagnostic
      write(0,*)'block 1 runs from JX= ',JXMIN1,' to ',JXMAX1
!***
! sanity check

      IF(NATS1/=NX1*NY1*NZ1)THEN
         CALL ERRMSG('FATAL','TARRCTBLK3',' NAT1.ne.NX1*NY1*NZ1')
      ENDIF

! now do block 2 (material 2)
! layer 1: extends from     0   to  xs1         [composition 2]
!          dipoles from     1   to NAT1
! layer 2: extends from -(xs1+xs2) to -xs1       [composition 3]
!          dipoles from   JXMIN2   to (JXMIN1-1)

      JXMIN2=-NINT(XS2/DX(1))
      JXMAX2=-1
      NX2=JXMAX2-JXMIN2+1
      JYMIN2=NINT(-0.5_WP*YS2/DX(2)-X0(2)+DELTA)
      JYMAX2=NINT(0.5_WP*YS2/DX(2)-X0(2)-DELTA)
      JZMIN2=NINT(-0.5_WP*ZS2/DX(3)-X0(3)+DELTA)
      JZMAX2=NINT(0.5_WP*ZS2/DX(3)-X0(3)-DELTA)
      NY2=JYMAX2-JYMIN2+1
      NZ2=JZMAX2-JZMIN2+1

! middle layer = composition 2:

      IF(JXMIN2.LT.JXMIN1)THEN
         DO JX=JXMIN2,JXMAX2
            DO JY=JYMIN2,JYMAX2
               DO JZ=JZMIN2,JZMAX2
                  NAT=NAT+1
                  IF(NAT>MXNAT)THEN
                     CALL ERRMSG('FATAL','TARSLBLIN',' NAT.GT.MXNAT ')
                  ENDIF
                  IXYZ(NAT,1)=JX
                  IXYZ(NAT,2)=JY
                  IXYZ(NAT,3)=JZ
                  DO JD=1,3
                     ICOMP(NAT,JD)=2
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      NATS2=NAT-NATS1

!*** diagnostic
      write(0,*)'block 2 runs from JX= ',JXMIN2,' to ',JXMAX2
!***
! sanity check

      IF(NATS2/=NX2*NY2*NZ2)THEN
         CALL ERRMSG('FATAL','TARRCTBLK3',' NATS2.ne.NX2*NY2*NZ2')
      ENDIF

! lower slab layer = composition 3

      JXMIN3=-NINT((XS2+XS3)/DX(1))
      JXMAX3=-NINT(XS2/DX(1))-1
      NX3=JXMAX3-JXMIN3+1

      JYMIN3=NINT(-0.5_WP*YS3/DX(2)-X0(2)+DELTA)
      JYMAX3=NINT(0.5_WP*YS3/DX(2)-X0(2)-DELTA)
      JZMIN3=NINT(-0.5_WP*ZS3/DX(3)-X0(3)+DELTA)
      JZMAX3=NINT(0.5_WP*ZS3/DX(3)-X0(3)-DELTA)
      NY3=JYMAX3-JYMIN3+1
      NZ3=JZMAX3-JZMIN3+1

      IF(JXMIN3.LT.JXMIN2)THEN
         DO JX=JXMIN3,JXMAX3
            DO JY=JYMIN3,JYMAX3
               DO JZ=JZMIN3,JZMAX3
                  NAT=NAT+1
                  IF(NAT>MXNAT)THEN
                     CALL ERRMSG('FATAL','TARSLBLIN',' NAT.GT.MXNAT ')
                  ENDIF
                  IXYZ(NAT,1)=JX
                  IXYZ(NAT,2)=JY
                  IXYZ(NAT,3)=JZ
                  DO JD=1,3
                     ICOMP(NAT,JD)=3
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
!***
         write(0,*)'layer 3 runs from JX= ',JXMIN3,' to ',JXMAX3
         write(0,*)'xs2=',xs2,' xs3=',xs3
!***
      ENDIF
      NATS3=NAT-NATS1-NATS2

! sanity check:

      IF(NATS3/=NX3*NY3*NZ3)THEN
         CALL ERRMSG('FATAL','TARRCTBLK3',' NATS3.ne.NX3*NY3*NZ3')
      ENDIF

      IF(NATS1*NATS2*NATS3>0)THEN
         WRITE(CMSGNM,FMT='(A,I7,A)')                          &
            ' Trilayer block structure with NAT=',NAT,' dipoles'
      ELSE
         WRITE(CMSGNM,FMT='(A,I7,A)')                         &
            ' Bilayer block structure with NAT=',NAT,' dipoles'
      ENDIF
      IF(IOSHP>=0)THEN
         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')
         WRITE(IOSHP,FMT=9020)XS1,YS1,ZS1,XS2,YS2,ZS2,XS3,YS3,ZS3, &
            NAT,NATS1,NATS2,NATS3,A1,A2,DX,X0
         DO JA=1,NAT
            WRITE(IOSHP,FMT=9030)JA,IXYZ(JA,1),IXYZ(JA,2),IXYZ(JA,3), &
                                 ICOMP(JA,1),ICOMP(JA,2),ICOMP(JA,3)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF
      RETURN
9020  FORMAT(' >TARRCTBLK3: XS1,YS1,ZS1,XS2,YS2,ZS2,XS3,YS3,ZS3=',9F8.4,/, &
         I10,3I7,' = NAT,NATS1,NATS2,NATS3',/,                             &
         3F10.6,' = A_1 vector',/,                                         &
         3F10.6,' = A_2 vector',/,                                         &
         3F10.6,' = lattice spacings (d_x,d_y,d_z)/d',/,                   &
         3F10.5,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ',         &
               'for dipole 0 0 0',/,                                       &
         '     JA  IX  IY  IZ ICOMP(x,y,z)')
9030  FORMAT(I7,3I4,3I2)
    END SUBROUTINE TARRCTBLK3
