    SUBROUTINE TARBLOCKS(A1,A2,DX,NBLOCKS,BLKSIZ,XYZB,X0,IPRINAX,IOSHP, &
        CDESCR,MXNAT,NAT,IXYZ,ICOMP)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

! Arguments:

      INTEGER :: BLKSIZ,IPRINAX,NBLOCKS,IOSHP,MXNAT,NAT
      INTEGER*2 :: ICOMP(MXNAT,3)
      INTEGER :: IXYZ(MXNAT,3)
      INTEGER :: XYZB(3,100)
      REAL(WP) :: & 
         A1(3),   &
         A2(3),   &
         DX(3),   &
         X0(3)
      CHARACTER :: CDESCR*67

! Local variables:

      INTEGER :: J,JB,JD,JJX,JJY,JJZ,JOFF,JOFFX,JOFFY,JOFFZ,JX,JY,JZ
      REAL(WP) :: EIGVAL(3)

!***********************************************************************

! Subroutine to construct a target out of "cubic" blocks.
! Blocks are only cubic for case of a cubic lattice DX=(1.,1.,1.)
! In case of a noncubic lattice, "blocks" are rectangular, with an equal
! number (=BLKSIZ) of dipole spacings in each direction.

! Given:
! DX(1-3)     =(dx/d,dy/d,dz/d) where dx,dy,dz = x,y,z lattice spacing,
!              and d = (dx*dy*dz)**(1/3) = effective lattice spacing
! NBLOCKS     =number of cubic blocks in target
! BLKSIZ      =(integer) ratio of block width/dipole spacing
!              (there will be BLKSIZ**3 dipoles per block)
! XYZB(1-3,J) =x,y,z coordinates of block J, in units of block width
! IPRINAX     =0. to return A1=(1,0,0), A2=(0,1,0)
!             =1. to return A1,A2=principal axes with largest, 2nd large
!              moment of inertia
! IOSHP       =device number for output file "target.out"
!             =-1 to suppress generation of file "target.out"
! MXNAT       =limit on largest allowed value of NAT=number of dipoles

! Returns:
! NAT         =number of dipoles in target
! IXYZ(J,1-3) =x,y,z location of dipole J on lattice
! ICOMP(J,1-3)=composition for dipole J; x,y,z directions
! A1(1-3)     =principal axis A1 in target frame (normalized)
! A2(1-3)     =principal axis A2 in target frame (normalized)
! X0(1-3)     =location in TF of dipole with IXYZ=0 0 0
!              TF origin is taken to be at centroid of target
! CDESCR      =string describing target

! where A1 = principal axis with largest moment of inertia
!       A2 = principal axis with second-largest moment of inertia

! B.T. Draine, Princeton Univ. Observatory, 95.12.11
! History:
! 96.01.02 (BTD) modified to compute eigenvalues and eigenvectors of
!                moment of inertia tensor, and to define vectors A1 and
!                A2 to be along eigenvectors with largest and second
!                largest eigenvalues, respectively.
! 96.01.26 (BTD) Replaced IX,IY,IZ by IXYZ
!                Removed principal axis calculation to separate
!                subroutine PRINAXIS
! 96.01.29 (BTD) Added variable IPRINAX to allow control over definition
!                of axes A1,A2
! 96.02.23 (BTD) Remove computation of principal axes and eigenvalues
!                as this is now carried out by routine PRINAXIS.
! 97.12.26 (BTD) Add DX(3) to argument list to support nonuniform
!                lattice spacing.
!                Add DX to argument list of PRINAXIS.
! 98.03.07 (BTD) Add DX to file target.out for compatibility with
!                new version of REASHP
! 03.11.06 (BTD) Modified to use new version of PRINAXIS with eigenvalue
!                in argument list
! 04.03.19 (BTD) Revert to previous version of PRINAXIS
! 04.05.23 (BTD) Added EIGVAL to argument list of PRINAXIS to conform
!                to change in PRINAXIS
! 07.09.11 (BTD) Changed IXYZ from INTEGER*2 to INTEGER
! 08.08.08 (BTD) Add X0(3) to argument list
!                Add code to define X0
! 08.08.29 (BTD) Modified format 9020
! End history.
! Copyright (C) 1996,1997,1998,2003,2004,2007,2008 
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

! For the moment, let us assume that blocks do not overlap.

      JOFF=1-((BLKSIZ+1)/2)
      JD=0
      DO JB=1,NBLOCKS
         JOFFX=JOFF+BLKSIZ*XYZB(1,JB)
         JOFFY=JOFF+BLKSIZ*XYZB(2,JB)
         JOFFZ=JOFF+BLKSIZ*XYZB(3,JB)
         DO JX=1,BLKSIZ
            JJX=JX+JOFFX
            DO JY=1,BLKSIZ
               JJY=JY+JOFFY
               DO JZ=1,BLKSIZ
                  JJZ=JZ+JOFFZ
                  JD=JD+1
                  IXYZ(JD,1)=JJX
                  IXYZ(JD,2)=JJY
                  IXYZ(JD,3)=JJZ

! Homogeneous target:

                  ICOMP(JD,1)=1
                  ICOMP(JD,2)=1
                  ICOMP(JD,3)=1
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      NAT=JD

      DO J=1,3
         X0(J)=0._WP
         DO JD=1,NAT
            X0(J)=X0(J)+REAL(IXYZ(JD,J),KIND=WP)
         ENDDO
! Set value of X0:
         X0(J)=-X0(J)/REAL(NAT,KIND=WP)
      ENDDO

      IF(IPRINAX/=0)THEN

! Call PRINAXIS to compute principal axes

         CALL PRINAXIS(MXNAT,NAT,ICOMP,IXYZ,DX,A1,A2,EIGVAL)
      ELSE

! Do not compute principal axes: simply set
! a1=(1,0,0) and a2=(0,1,0) in target frame.

         A1(1)=1._WP
         A1(2)=0._WP
         A1(3)=0._WP
         A2(1)=0._WP
         A2(2)=1._WP
         A2(3)=0._WP
      ENDIF

!***********************************************************************
! Write target description into string CDESCR
! Note: this will replace whatever CDESCR was read in from the
! 'blocks.par' file. If you want to retain that description,
! then comment out the following statement:
!      WRITE(CDESCR,FMT='(A,I7,A,I3,A)')'Target containing',
!     &      NAT,' dipoles arranged in ',NBLOCKS,' cubic blocks'
!***********************************************************************

      IF(IOSHP>=0)THEN
         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')
         WRITE(IOSHP,FMT=9020)NAT,A1,A2,DX,X0
         DO JX=1,NAT
            WRITE(IOSHP,FMT=9030)JX,IXYZ(JX,1),IXYZ(JX,2),IXYZ(JX,3), &
                                 ICOMP(JX,1),ICOMP(JX,2),ICOMP(JX,3)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF
      RETURN
9020  FORMAT(' >TARBLOCKS: cubic building blocks',/,              &
         I9,' = NAT',/,                                           &
         3F9.4,' = A_1 vector',/,                                 &
         3F9.4,' = A_2 vector',/,                                 &
         3F9.6,' = lattice spacings (d_x,d_y,d_z)/d',/,           &
         3F9.6,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ', &
               'for dipole 0 0 0',/,                              &
         '     JA  IX  IY  IZ ICOMP(x,y,z)')
9030  FORMAT(I7,3I4,3I2)
    END SUBROUTINE TARBLOCKS
