    SUBROUTINE NULLER(CXVEC,IOCC,MXNAT,MXN3,NAT)
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE
! Arguments:
      INTEGER :: MXNAT,MXN3,NAT
      INTEGER*2 :: &
         IOCC(MXNAT)
      COMPLEX(WP) :: &
         CXVEC(NAT,3)
! Local variables:
      INTEGER :: J, L
!***********************************************************************
! This routine takes as input a complex vector with NAT3 elements
! and sets to zero those elements corresponding to "vacuum" sites.
! This is accomplished by multiplying by vector IOCC(J), whose
! elements are either 1 or 0 depending on whether site is occupied
! or "vacuum".
! Note: CXVEC should be given dimension CXVEC(NAT,3) here
!       so that first 3*NAT elements of CXVEC are employed.

! B.T.Draine, Princeton Univ. Obs., 90/11/1
! History:
! 90/12/15 (BTD): Corrected error in dimensioning of CXVEC.

! Copyright (C) 1993, B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
      DO L=1, 3
         DO J=1,NAT
            IF(IOCC(J)==0)THEN
               CXVEC(J,L)=IOCC(J)*CXVEC(J,L)
            ENDIF
         ENDDO
      ENDDO
      RETURN
    END SUBROUTINE NULLER
