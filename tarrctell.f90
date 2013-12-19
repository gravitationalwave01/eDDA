    SUBROUTINE TARRCTELL(A1,A2,AX,AY,AZ,BX,BY,BZ,DX,X0,CDESCR,IOSHP, &
                         MXNAT,NAT,IXYZ,ICOMP)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

!** Arguments:

      CHARACTER :: CDESCR*67
      INTEGER :: IOSHP,MXNAT,NAT
      INTEGER*2 ::    &
         ICOMP(MXNAT,3)
      INTEGER ::     &
         IXYZ(MXNAT,3)
      REAL(WP) :: AX,AY,AZ,BX,BY,BZ
      REAL(WP) :: & 
         A1(3),   &
         A2(3),   &
         DX(3),   &
         X0(3)

!** Local variables:

      INTEGER :: JX,JY,JZ,LMX1,LMX2,LMY1,LMY2,LMZ1,LMZ2,NIN
      REAL(WP) :: AX2,AY2,AZ2,BX2,BY2,BZ2,R,RR,RYZ2,RZ2,X,XOFF,Y,YOFF,Z,ZOFF
!***********************************************************************
! Subroutine TARRCTELL
! Routine to construct target consisting of ellipsoid of dimensions
! AX*d, AY*d, AZ*d
! embedded in a concentric rectangular volume of dimensions
! BX*d, BY*d, BZ*d
! Input:
!        AX=(x-length)/d of ellipsoid  (d=lattice spacing) of material 1
!        AY=(y-length)/d "   "     "
!        AZ=(z-length)/d "   "     "
!        BX=(x-length)/d of larger rectangular volume of material 2
!        BY=(y-length)/d "   "     "
!        BZ=(z-length)/d "   "     "
!        DX(1-3)=lattice spacing (dx/d,dy/d,dz/d), d=(dx*dy*dz)**(1/3)
!        IOSHP=device number for "target.out" file
!             =-1 to suppress printing of "target.out"
!        MXNAT=dimensioning information (max number of atoms)
! Output:
!        A1(1-3)=(1,0,0)=unit vector defining target axis 1 in Target Frame
!        A2(1-3)=(0,1,0)=unit vector defining target axis 2 in Target Frame
!        NAT=number of atoms in target
!        IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
!        ICOMP(1-NAT,1-3)=1 for sites within inner ellipsoid,
!                        =2 for sites between inner and outer ellipsoids
!        CDESCR=description of target (up to 67 characters)
!        X0(1-3)=(location/d) in Target Frame corresponding to dipole with
!                IXYZ=(0,0,0).  This will be treated at the origin of physical
!                coordinates in the TF.
!                Here set to be centroid of the cube and enclosed ellipsoid 

! B.T.Draine, Princeton Univ. Obs.
! History:
! 11.10.18 (BTD): Adapted from TARCEL
! end history

! Copyright (C) 2011 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
!*** diagnostic
!      write(0,*)'tarrctell ckpt 1'
!      write(0,*)'  ax,ay,az=',ax,ay,az
!      write(0,*)'  bx,by,bz=',bx,by,bz
!***
! Current version of TARRCTELL is restricted to cubic lattices

      IF(DX(1)/=1._WP.OR.DX(2)/=1._WP)THEN
         CALL ERRMSG('FATAL','TARRCTELL',                          &
                     ' tarrctell does not support noncubic lattice')
      ENDIF

! Inner ellipsoid is defined by
!      (X/AX)**2+(Y/AY)**2+(Z/AZ)**2 = 0.25
! Outer rectangular volume has sides BX,BY,BZ
! Material within ellipsoid is of composition 1
! Material between ellipsoids and outer surface is of composition 2

! Dipoles are located at sites
! (x,y,z)=(I+XOFF,J+YOFF,Z+KOFF), I,J,K=integers
!                                 XOFF,YOFF,ZOFF=constants

!***********************************************************************
!*** diagnostic
!      write(0,*)'tarrctell ckpt 2'
!***
! Check that input parameters have "inner" ellipsoid smaller than
! "outer" rectangular volume:
      IF(AX>BX)CALL ERRMSG('FATAL','TARCEL',' AX > BX ')
      IF(AY>BY)CALL ERRMSG('FATAL','TARCEL',' AY > BY ')
      IF(AZ>BZ)CALL ERRMSG('FATAL','TARCEL',' AZ > BZ ')

      X0(1)=0.
      X0(2)=0.
      X0(3)=0.
      JX=NINT(BX)
      IF(2*(JX/2)==JX)THEN
         XOFF=0.5_WP
      ELSE
         XOFF=0._WP
      ENDIF
! Similar criterion for YOFF:
      JY=NINT(BY)
      IF(2*(JY/2)==JY)THEN
         YOFF=0.5_WP
      ELSE
         YOFF=0._WP
      ENDIF
! Similar criterion for ZOFF:
      JZ=NINT(BZ)
      IF(2*(JZ/2)==JZ)THEN
         ZOFF=0.5_WP
      ELSE
         ZOFF=0._WP
      ENDIF
!*** diagnostic
!      write(0,*)'tarrctell ckpt 3'
!***
! JX even: run from -JX/2 to JX/2-1      e.g. -2 to +1   if JX=4, XOFF=0.5
!    odd :          -(JX/2) to (JX/2)    e.g. -2 to +2   if JX=5
      LMX1=-(JX/2)
      LMX2=(JX/2)-NINT(2*XOFF)
      LMY1=-(JY/2)
      LMY2=(JY/2)-NINT(2*YOFF)
      LMZ1=-(JZ/2)
      LMZ2=(JZ/2)-NINT(2*ZOFF)

      AX2=AX*AX
      AY2=AY*AY
      AZ2=AZ*AZ

      NAT=0
      NIN=0

! Specify target axes A1 and A2
!   A1=(1,0,0) in target frame
!   A2=(0,1,0) in target frame

      DO JX=1, 3
         A1(JX)=0._WP
         A2(JX)=0._WP
      ENDDO
      A1(1)=1._WP
      A2(2)=1._WP

! Determine list of occupied sites.
!*** diagnostic
!      write(0,*)'tarrctell ckpt 4'
!      write(0,*)'  lmz1,lmz2=',lmz1,lmz2
!      write(0,*)'  lmy1,lmy2=',lmy1,lmy2
!      write(0,*)'  lmx1,lmx2=',lmx1,lmx2
!***
      DO JZ=LMZ1,LMZ2
!*** diagnostic
!         write(0,*)'tarrctell ckpt 5, jz=',jz
!***
         Z=REAL(JZ,KIND=WP)+ZOFF
         RZ2=Z*Z/AZ2
         DO JY=LMY1,LMY2
            Y=REAL(JY,KIND=WP)+YOFF
            RYZ2=RZ2+Y*Y/AY2
            DO JX=LMX1,LMX2
               NAT=NAT+1
               X=REAL(JX,KIND=WP)+XOFF
               R=RYZ2+X*X/AX2
!*** diagnostic
!               write(0,*)'tarrctell ckpt 8, jx,jy,jz=',jx,jy,jz
!               write(0,*)'   x=',x,' R=',R
!               write(0,*)'   nat=',nat
!***
               IF(R<0.25_WP)THEN
! Site is within ellipsoid
                  ICOMP(NAT,1)=1
                  ICOMP(NAT,2)=1
                  ICOMP(NAT,3)=1
                  NIN=NIN+1
               ELSE
! in volume outside ellipsoid

!*** diagnostic
!                  write(0,*)'tarrctell ckpt 8.2'
!***
                  ICOMP(NAT,1)=2
                  ICOMP(NAT,2)=2
                  ICOMP(NAT,3)=2
               ENDIF
!*** diagnostic
!               write(0,*)'tarrctell ckpt 9'
!***
               IXYZ(NAT,1)=JX
               IXYZ(NAT,2)=JY
               IXYZ(NAT,3)=JZ
            ENDDO
         ENDDO
      ENDDO
!*** diagnostic
!      write(0,*)'tarrctell ckpt 10'
!***
! locate centroid = origin of coordinates in TF

      DO JX=1,3
         X0(JX)=0._WP
      ENDDO
      DO JX=1,3
         DO JY=1,NAT
            X0(JX)=X0(JX)+REAL(IXYZ(JY,JX))
         ENDDO
         X0(JX)=X0(JX)/REAL(NAT)
      ENDDO
      
      IF(NAT>MXNAT)THEN
         CALL ERRMSG('FATAL','TARELL',' NAT.GT.MXNAT ')
      ENDIF

!***********************************************************************
! Write target description into string CDESCR

      IF(AX==AY.AND.AX==AZ)THEN
         WRITE(CDESCR,FMT='(A,I7,A)')' spherical with',NIN, &
                                     ' dipoles'
      ELSE
         WRITE(CDESCR,FMT='(A,I7,A)')' Ellipsoid with',NIN, &
                                     ' dipoles'
      ENDIF

!***********************************************************************

      IF(IOSHP>=0)THEN
         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')
         WRITE(IOSHP,FMT=9020)AX,AY,AZ,BX,BY,BZ,NAT,NIN,A1,A2,DX,X0
         DO JX=1,NAT
            WRITE(IOSHP,FMT=9030)JX,IXYZ(JX,1),IXYZ(JX,2),IXYZ(JX,3), &
                                 ICOMP(JX,1),ICOMP(JX,2),ICOMP(JX,3)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF
      RETURN
9020  FORMAT(' >TARRCTELL: ellipsoid; AX,AY,AZ=',3F8.4,            &
         ' in rect. vol. BX,BY,BZ=',3F8.4,/,                       &
         2I10,' = NAT, NIN,',/,                                    &
         3F10.6,' = A_1 vector',/,                                 &
         3F10.6,' = A_2 vector',/,                                 &
         3F10.6,' = lattice spacings (d_x,d_y,d_z)/d',/,           &
         3F10.5,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ', &
               'for dipole 0 0 0',/,                               &
         '     JA  IX  IY  IZ ICOMP(x,y,z)')
9030  FORMAT(I7,3I4,3I2)
    END SUBROUTINE TARRCTELL
