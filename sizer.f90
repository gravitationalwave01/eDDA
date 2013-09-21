    SUBROUTINE SIZER(MXNAT,NAT,IXYZ,LXYZ)
      IMPLICIT NONE
! Arguments:
      INTEGER :: MXNAT,NAT
      INTEGER :: &
         LXYZ(3)
      INTEGER :: &
         IXYZ(MXNAT,3)
! Local variables:
      CHARACTER :: CMSGNM*70
      INTEGER :: JA,JD,JX
      INTEGER ::   &
         LXMIN(3), &
         LXMAX(3)
!***********************************************************************

! Purpose: to determine full extent of target in x,y,z directions
! in order to determine what "computational" volume will be.
! B.T.Draine, Princeton University Observatory, 96.01.26
! History:
! 96.11.21 (BTD) modified to replace WRITE(0,... with CALL WRIMSG
! 07.09.11 (BTD) changed IXYZ from INTEGER*2 to INTEGER
! end history
! Copyright (C) 1996,2007 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
      DO JD=1,3
        LXMIN(JD)=IXYZ(1,JD)
        LXMAX(JD)=LXMIN(JD)
      ENDDO
      DO JD=1,3
        DO JA=2,NAT
          JX=IXYZ(JA,JD)
          IF(JX<LXMIN(JD))LXMIN(JD)=JX
          IF(JX>LXMAX(JD))LXMAX(JD)=JX
        ENDDO
      ENDDO
      DO JD=1,3
        LXYZ(JD)=LXMAX(JD)-LXMIN(JD)+1
      ENDDO
      JA=LXYZ(1)*LXYZ(2)*LXYZ(3)
      WRITE(CMSGNM,6000) LXYZ,JA
      CALL WRIMSG('SIZER ',CMSGNM)
6000  FORMAT('min. comp. vol.: L_x*L_y*L_z=',I5,' *',I5,' *',I5,' =',I8)
      RETURN
    END SUBROUTINE SIZER
