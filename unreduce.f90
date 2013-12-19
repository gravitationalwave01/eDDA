    SUBROUTINE UNREDUCE(CXV,IOCC,MXN3,MXNAT,NAT,NAT0)
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE
! ------------------------ unreduce_v3 ---------------------------------
! Arguments

      INTEGER :: MXN3,MXNAT,NAT,NAT0
      INTEGER*2 :: IOCC(MXNAT)
      COMPLEX(WP) :: CXV(MXN3)

! Local variables:

      COMPLEX(WP) :: CXZERO
      INTEGER :: J,JOCC,M

!***********************************************************************
! subroutine UNREDUCE
! given
!   CXV(1:MXN3)   = complex vector in "reduced" form
!   IOCC(1:MXNAT) = 0/1 if site is vacant/occupied
!   MXN3          = dimensioning information
!   MXNAT         = dimensioning information
!   NAT           = number of lattice sites
!   NAT0          = number of occupied lattice sites
! returns
!   CXV(1:MXN3)   = complex vector in "natural" form

! history:
! 11.08.16 (BTD) created for use in version 7.2.1
! 12.08.09 (BTD) added code to ensure cxv is initialized to zero
!                if not defined by input CXV(1:MXN3)
! 13.04.22 (BTD) modified to correct problem noted by C. Vasyl
! end history
! Copyright (C) 2011,2012,2013
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
      CXZERO=(0._WP,0._WP)

!*** diagnostic
!      write(0,*)'unreduce ckpt 1'
!***
      IF(NAT>NAT0)THEN
         DO J=3*NAT0+1,3*NAT
            CXV(J)=CXZERO
         ENDDO
         DO M=2,0,-1
            JOCC=NAT0
            DO J=NAT,1,-1
               IF(IOCC(J)>=1)THEN
                  CXV(J+M*NAT)=CXV(JOCC+M*NAT0)
                  JOCC=JOCC-1
               ENDIF
            ENDDO
         ENDDO
         DO J=1,NAT
            IF(IOCC(J)==0)THEN
               DO M=0,2
                  CXV(J+M*NAT)=CXZERO
               ENDDO
            ENDIF
         ENDDO
      ENDIF


!*** diagnostic
!      write(0,*)'unreduce ckpt 9, jocc=',jocc,' (should be zero...)'
!***
      RETURN
    END SUBROUTINE UNREDUCE
