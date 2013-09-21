    SUBROUTINE TAR2SP(A1,A2,A_1,B_1,A_2,B_2,PHI,PRINAX,DX,X0,CDESCR,IOSHP, &
                      MXNAT,NAT,IXYZ,ICOMP)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

! Arguments:

      CHARACTER :: CDESCR*67
      INTEGER :: IOSHP, MXNAT, NAT
      INTEGER*2 :: ICOMP(MXNAT,3)
      INTEGER :: IXYZ(MXNAT,3)
      REAL(WP) :: A_1,A_2,B_1,B_2,PHI,PRINAX
      REAL(WP) :: A1(3),A2(3),DX(3),X0(3)

! Local variables:

      CHARACTER :: CMSGNM*70
      INTEGER :: JA,JX,JY,JZ,LMX1,LMX2,LMY1,LMY2,LMZ1,LMZ2,NAT1,NAT2
      REAL(WP) :: AX2,AY2,AZ2,COSPHI,SINPHI,U,V,R,RX2,RYZ2,RZ2, &
                  X,XC1,XC2,Y,YC1,YC2,Z,ZC1,ZC2
      REAL(WP) :: EIGVAL(3)

!***********************************************************************
! Routine to construct two touching spheroids from "atoms"
!    First spheroid has symmetry axis in y direction.
!    Second spheroid is displaced in x direction, with symmetry
!    axis in yz plane, in direction ey*cos(phi)+ez*sin(phi).
!    Separation between centroids=(B_1+B_2)/2

! Input:
!        A_1=length/d of 1st spheroid along symm.axis  (d=lattice spacing)
!        A_2=length/d of 2nd spheroid along symm.axis
!        B_1=diameter/d of 1st spheroid
!        B_2=diameter/d of 2nd spheroid
!        PHI=angle (deg) specifying relative orientation of spheroids
!        PRINAX=0. to set A1=(1,0,0),A2=(0,1,0) in Target Framee
!              =1. to set A1,A2=principal axes of largest and second
!                               largest moment of inertia
!        DX(1-3)=(dx/d,dy/d,dz/d) where dx,dy,dz=lattice spacing in x,y,z
!            directions, and d=(dx*dy*dz)**(1/3)=effective lattice spacing
!        IOSHP=device number for "target.out" file
!             =-1 to suppress printing of "target.out"
!        MXNAT=dimensioning information (max number of atoms)
! Output:
!        A1(1-3)=(1,0,0)=unit vector defining target axis 1 in Target Frame
!        A2(1-3)=(0,1,0)=unit vector defining target axis 2 in Target Frame
!        X0(1-3)=location in TF of dipole with IXYZ=0 0 0
!        NAT=number of atoms in target
!        IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
!        ICOMP(1-NAT,1-3)=1 for sites within first spheroid
!                        =2 for sites within second spheroid
!        CDESCR=description of target (up to 67 characters)

! B.T.Draine, Princeton Univ. Obs.
! History:
! 96.01.25 (BTD): First written.
! 96.01.26 (BTD): Replaced IX,IY,IZ by IXYZ
! 96.01.29 (BTD): Modified to add PRINAX to argument list, and
!                 call to PRINAXIS to determine principal axes if
!                 PRINAX=1.
! 97.12.26 (BTD): Added DX(3) to argument list to support nonuniform
!                 lattice spacing.
!                 *** Note: additional code modifications required
!                     to support noncubic lattice!
! 98.03.06 (BTD): If called with noncubic lattice (not yet supported)
!                 call ERRMSG and halt with fatal error.
! 98.03.07 (BTD): Modified to write DX to file "target.out" for
!                 compatibility with REASHP.
! 99.06.30 (BTD): Revised procedure for constructing target.
!                 With new prescription,
!                 should obtain a symmetric structure when called
!                 with A_1,B_1 = A_2,B_2 and PHI=0 or 90
!                 Note that we have not yet made modifications
!                 necessary to support noncubic lattice.
! 00.11.02 (BTD): Modified to write ICOMP to target.out
! 03.11.06 (BTD): Modified to use new PRINAXIS with eigenvalues in
!                 argument list.
! 04.03.19 (BTD): Revert to previous version of PRINAXIS
! 04.05.23 (BTD): Modify argument list of PRINAXIS to conform to
!                 change in PRINAXIS (eigenvalues added to arg list)
! 07.09.11 (BTD): Change IXYZ from INTEGER*2 to INTEGER
! 08.08.08 (BTD): Add X0(3) to argument list
!                 Add code to determine X0
! 08.08.29 (BTD): modified format 9020
! end history

! Copyright (C) 1996,1997,1998,1999,2000,2003,2004,2007,2008
!                B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

! Current version of TAR2SP is restricted to cubic lattices

      IF(DX(1)/=1._WP.OR.DX(2)/=1._WP)THEN
         CALL ERRMSG('FATAL','TAR2SP',                              &
                     ' tar2sp does not yet support noncubic lattice')
      ENDIF

! Routine to construct target array representing two spheroids in
! contact.
! line connecting spheroid centers is parallel to x
! symmetry axis of first spheroid is parallel to y
! symmetry axis of second spheroid is in y-z plane, at angle PHI to y

! First spheroid:
! Array sites contained within spheroidal surface defined by
!    ((X-XC1)/B_1)**2+((Y-YC1)/A_1)**2+((Z-ZC1)/B_1)**2=0.25
! Ideal volume V=(pi/6)*A_1*B_1*B_1

! Second spheroid:
!    ((X-XC2)/B_2)**2+((U-UC1)/A_1)**2+((V-VC1)/B_1)**2=0.25
! where
!    U  =  Y*cos(PHI)  +  Z*sin(PHI)
!    V  = -Y*sin(PHI)  +  Z*cos(PHI)
!   UC1 = YC1*cos(PHI) + ZC1*sin(PHI)
!   VC1 =-YC1*sin(PHI) + ZC1*cos(PHI)

! Dipoles are located at sites
! (x,y,z)=(I,J,K)*d, I,J,K=integers

! Criterion for choosing XC1:
! Point where spheroids contact each other should be midway between
! lattice points.  Thus take XC1+0.5*B_1 = 0.5
! or XC1 = 0.5(1.-B_1)

! Criterion for choosing YC1 and ZC1:
! If XC1 is close to an integer, take YC1=ZC1=0
! If XC1 is close to a half-integer, take YC1=ZC1=0.5

! Note that other criteria are obviously possible.  For example,
! could have chosen YC1 and ZC1 to maximize number of lattice
! sites falling within spheroidal surfaces.

      XC1=0.5_WP*(1._WP-B_1)
      IF(XC1-INT(XC1)<0.25_WP.OR.XC1-INT(XC1)>0.75_WP)THEN
         YC1=0._WP
         ZC1=0._WP
      ELSE
         YC1=0.5_WP
         ZC1=0.5_WP
      ENDIF

      NAT=0

! Determine list of occupied sites in first spheroid

      LMX1=INT(XC1-.5_WP*B_1)
      LMX2=INT(XC1+.5_WP*B_1)
      LMY1=INT(YC1-.5_WP*A_1)
      LMY2=INT(YC1+.5_WP*A_1)
      LMZ1=INT(ZC1-.5_WP*B_1)
      LMZ2=INT(ZC1+.5_WP*B_1)
      AX2=B_1*B_1
      AY2=A_1*A_1
      AZ2=B_1*B_1
      DO JZ=LMZ1,LMZ2
         Z=REAL(JZ,KIND=WP)-ZC1
         RZ2=Z*Z/AZ2
         IF(RZ2<0.25_WP)THEN
            DO JY=LMY1,LMY2
               Y=REAL(JY,KIND=WP)-YC1
               RYZ2=RZ2+Y*Y/AY2
               IF(RYZ2<0.25_WP)THEN
                  DO JX=LMX1,LMX2
                     X=REAL(JX,KIND=WP)-XC1
                     R=RYZ2+X*X/AX2
                     IF(R<0.25_WP)THEN
                        NAT=NAT+1
                        IXYZ(NAT,1)=JX
                        IXYZ(NAT,2)=JY
                        IXYZ(NAT,3)=JZ
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      NAT1=NAT

! Determine occupied sites in second spheroid

      XC2=XC1+.5_WP*(B_1+B_2)
      YC2=YC1
      ZC2=ZC1
      LMX1=INT(XC2-.5_WP*B_2)
      LMX2=INT(XC2+.5_WP*B_2)

! This is not time-consuming, so no need to optimize choices of
! LMY1,LMY2,LMZ1,LMZ2

      LMY1=INT(YC1-.5_WP*MAX(A_2,B_2))
      LMY2=INT(YC1+.5_WP*MAX(A_2,B_2))
      LMZ1=INT(ZC1-.5_WP*MAX(A_2,B_2))
      LMZ2=INT(ZC1+.5_WP*MAX(A_2,B_2))
      SINPHI=SIN(3.1415927_WP*PHI/180._WP)
      COSPHI=COS(3.1415927_WP*PHI/180._WP)

! transform (y,z) -> (u,v) with u=y*cos(phi)-z*sin(phi)
!                               v=y*sin(phi)+z*cos(phi)

      DO JX=LMX1,LMX2
         X=REAL(JX,KIND=WP)-XC2
         RX2=(X/B_2)**2
         IF(RX2<0.25_WP)THEN
            DO JY=LMY1,LMY2
               Y=REAL(JY,KIND=WP)-YC2
               DO JZ=LMZ1,LMZ2
                  Z=REAL(JZ,KIND=WP)-ZC2
                  U=(COSPHI*Y-SINPHI*Z)/A_2
                  V=(SINPHI*Y+COSPHI*Z)/B_2
                  R=RX2+U**2+V**2
                  IF(R<0.25_WP)THEN
! Site is occupied:
                     NAT=NAT+1
                     IXYZ(NAT,1)=JX
                     IXYZ(NAT,2)=JY
                     IXYZ(NAT,3)=JZ
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      NAT2=NAT-NAT1
      IF(NAT>MXNAT)THEN
        CALL ERRMSG('FATAL','TAR2SP',' NAT.GT.MXNAT ')
      ENDIF

      DO JX=1,3
         DO JA=1,NAT1
            ICOMP(JA,JX)=1
         ENDDO
         DO JA=NAT1+1,NAT
            ICOMP(JA,JX)=2
         ENDDO
      ENDDO

! Determine X0(1-3)=TF coordinates/d of dipole with IXYZ=0 0 0

      X0(1)=-0.5_WP
      X0(2)=-YC1
      X0(3)=-ZC1

! Specify target axes A1 and A2
! If PRINAX=0. then
!    A1=(1,0,0) in target frame
!    A2=(0,1,0) in target frame
! If PRINAX=1. then
!    A1,A2 are principal axes of largest,second largest moment of
!    inertia

      IF(PRINAX<=0._WP)THEN
         DO JX=1,3
            A1(JX)=0._WP
            A2(JX)=0._WP
         ENDDO
         A1(1)=1._WP
         A2(2)=1._WP
      ELSE
         CALL PRINAXIS(MXNAT,NAT,ICOMP,IXYZ,DX,A1,A2,EIGVAL)
      ENDIF

!***********************************************************************
! Write target description into string CDESCR

      WRITE(CDESCR,FMT='(A,I7,A,I7,A,I7,A)')' Spheroids with',NAT1,' +', &
                                            NAT2,'=',NAT,' dipoles'

!***********************************************************************

      CMSGNM=CDESCR
      CALL WRIMSG('TAR2SP',CMSGNM)
      IF(IOSHP>=0)THEN
         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')
         WRITE(IOSHP,FMT=9020)A_1,B_1,A_2,B_2,PHI,NAT,NAT1,NAT2,A1,A2,DX,X0
         DO JX=1,NAT
            WRITE(IOSHP,FMT=9030)JX,IXYZ(JX,1),IXYZ(JX,2),IXYZ(JX,3), &
                                 ICOMP(JX,1),ICOMP(JX,2),ICOMP(JX,3)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF

      RETURN
9020  FORMAT(' >TAR2SP: touching spheroids, A_1,B_2,A_2,B_2,PHI=',5F8.4,/, &
         3I10,' = NAT, NAT1, NAT2',/,                                       &
         3F10.6,' = A_1 vector',/,                                          &
         3F10.6,' = A_2 vector',/,                                          &
         3F10.6,' = lattice spacings (d_x,d_y,d_z)/d',/,                    &
         3F10.5,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ',          &
               'for dipole 0 0 0',/,                                       &
         '     JA  IX  IY  IZ ICOMP(x,y,z)')
9030  FORMAT(I7,3I4,3I2)
    END SUBROUTINE TAR2SP
