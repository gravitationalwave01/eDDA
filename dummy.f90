! These routines were created by Art Lazanoff in order to
! satisfy dependencies when passing dummy routines as argument
! to various subroutines.
! Added to DDSCAT 2008.05.12
      SUBROUTINE DUMMY(D1,D2,I1)
      USE DDPRECISION,ONLY : WP
      REAL(WP)::D1,D2
      INTEGER I1
      END SUBROUTINE DUMMY
      SUBROUTINE CDUMMY(D1,D2,I1)
      USE DDPRECISION,ONLY : WP
      REAL(WP)::D1,D2
      INTEGER I1
      END SUBROUTINE CDUMMY
      SUBROUTINE CDUMMY_1(D1,D2,I1)
      USE DDPRECISION,ONLY : WP
      REAL(WP)::D1,D2
      INTEGER I1
      END SUBROUTINE CDUMMY_1
