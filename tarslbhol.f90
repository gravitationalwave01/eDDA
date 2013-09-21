    SUBROUTINE TARSLABHOLE(A1,A2,A,BA,CA,RA,DX,X0,CDESCR,IOSHP, &
                           MXNAT,NAT,IXYZ,ICOMP)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE
! Arguments:
      REAL(WP) :: A,BA,CA,RA
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

      CHARACTER :: CMSGNM*70
      INTEGER :: IC,JA,JD,JX,JY,JZ,NCOMP,NX,NX1,NX2,NX3,NY,NZ
      REAL(WP) :: R2,RAD2,RADX2

! External subroutines:

      EXTERNAL ERRMSG

! Intrinsic functions:

      INTEGER :: NINT
      INTRINSIC NINT

!***********************************************************************
! Routine to construct rectangular slab with dimensions A*d, B*d, C*d
! in x_tf,y_tf,z_tf directions
! with circular hole of radius R*d
! hole axis is in x_tf direction, through center of slab
!
! x0(1-3) = point at center of hole, at top surface of slab
!
!
! Input:
!        A=(x-length)/d     (d=lattice spacing)
!        BA=B/A  
!        CA=C/A
!        RA=R/A
!
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
!        IXYZ(1-NAT,1-3)=location indices for atoms of target
!        ICOMP(1-NAT,1-3)=composition
!        X0(3)=location/d in TF of lattice site with IXYZ=(0,0,0)
!              X(1-3)=[IXYZ(J,1-3)+X0(1-3)]*d

! B.T.Draine, Princeton Univ. Obs.
! History:
! 08.04.24 (BTD): First written
! 10.02.01 (BTD): Added to ddscat 7.1.0
! 10.02.02 (BTD): Further mods
! 10.02.06 (BTD): added X0 to target.out output
! end history

! Copyright (C) 2008,2010 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
! Sanity check:

      IF(2._WP*RA>SQRT(BA**2+CA**2)) &
      CALL ERRMSG('FATAL','TARSLABHOLE','hole diam > block diag')

      CDESCR='SLAB_HOLE: rect. block with cylindrical hole'

      NX=NINT(A/DX(1))
      NY=NINT(BA*A/DX(2))
      NZ=NINT(CA*A/DX(3))

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
! Lowest layer will have JX=-NX

! Specify X0(3) = (x,y,z)/d in TF corresponding to IXYZ=(0,0,0)
! Top surface of slab (0.5d above topmost dipole layer) is assumed to have x=0
! Dipoles with JX=-1 are assumed to be at x=-dx/2 (top layer of target).
! Therefore X0(1)=0.5*DX(1)
! Set y=0 and z=0 to be at middle of slab
! JY runs from 1 to NY -> X0(2) = -(1+NY)/2
! JZ runs from 1 to NZ -> X0(3) = -(1+NZ)/2

      X0(1)=0.5*DX(1)
      X0(2)=-0.5*REAL(1+NY)
      X0(3)=-0.5*REAL(1+NZ)

      JA=0
      R2=(RA*A)**2
      NAT=0

! diagnostic
      write(0,*)'tarslbhol ckpt 1,  a=',a
      write(0,*)'                  ba=',ba
      write(0,*)'                  ca=',ca
      write(0,*)'                  ra=',ra
      write(0,*)'                  dx=',dx
      write(0,*)'                  nx,ny,nz=',nx,ny,nz
      write(0,*)'                  r2=',r2
      write(0,*)'                  x0=',x0
!
      DO JY=1,NY
         RADX2=((JY+X0(2))/DX(2))**2
         DO JZ=1,NZ
            RAD2=RADX2+((JZ+X0(3))/DX(3))**2
            IF(RAD2>R2)THEN
! sweep downward starting from upper surface
               DO JX=-1,-NX,-1
! diagnostic
!                     write(0,*)'tarslbhol ckpt 2, jx=',jx
!
                  NAT=NAT+1
                  IF(NAT>MXNAT)THEN
                     CALL ERRMSG('FATAL','TARHOLESLAB',' NAT > MXNAT')
                  ENDIF
                  IXYZ(NAT,1)=JX
                  IXYZ(NAT,2)=JY
                  IXYZ(NAT,3)=JZ
                  DO JD=1,3
                     ICOMP(NAT,JD)=1
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

      WRITE(CMSGNM,FMT='(A,I7,A)')' total slab x-thickness =',NX,' dipoles'
      CALL WRIMSG('TARSLABHOLE',CMSGNM)
      WRITE(CMSGNM,FMT='(A,I7,A)')' slab y-width =',NY,' dipoles' 
      CALL WRIMSG('TARSLABHOLE',CMSGNM)
      WRITE(CMSGNM,FMT='(A,I7,A)')' slab z-width =',NZ,' dipoles'
      CALL WRIMSG('TARSLABHOLE',CMSGNM)
      IF(IOSHP>=0)THEN
         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')
         WRITE(IOSHP,FMT=9020)A,BA,CA,RA,NAT,A1,A2,X0,DX
         DO JA=1,NAT
            WRITE(IOSHP,FMT=9030)JA,IXYZ(JA,1),IXYZ(JA,2),IXYZ(JA,3), &
                                 ICOMP(JA,1),ICOMP(JA,2),ICOMP(JA,3)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF

      RETURN
9020  FORMAT(                                                              &
        ' >TARREC rectangular block with circular hole; A/d,B/d/C/d,R/d=', &
           4F8.4,/,                                                        &
        I7,' = NAT',/,                                                    &
        3F9.4,' = A_1 vector',/,                                           &
        3F9.4,' = A_2 vector',/,                                           &
        3F9.6,' = lattice spacings (d_x,d_y,d_z)/d',/,                     &
        3F10.5,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ',          &
           'for dipole 0 0 0',/,                                           &
        '     JA  IX  IY  IZ ICOMP(x,y,z)')
9030  FORMAT (I7,3I5,3I2)
    END SUBROUTINE TARSLABHOLE
