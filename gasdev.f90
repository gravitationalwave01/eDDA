    FUNCTION GASDEV(IDUM)
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE

! arguments:

      INTEGER :: IDUM
      REAL(WP) :: GASDEV

! local variables:

      INTEGER :: ISET
      REAL(WP) :: FAC,GSET,ONE,RSQ,TWO,V1,V2,ZERO

! external functions

      REAL(WP) :: RAN3

      SAVE ISET,GSET
      DATA ISET/0/

!=======================================================================
! function GASDEV
!
! given:

!    IDUM = seed for random number generator
!            if IDUM < 0, reinitialize RAN1 with seed IDUM
!            if IDUM .ge. 0, continue sequence from RAN1
! returns:

!    GASDEV(IDUM) = gaussian random deviate
!                   with mean <gasdev>= 0 
!                   and variance <gasdev**2>=1

! This subroutine follows closely the structure of the routine gasdev
! described in Numerical Recipes by Press, W.H., Flannery, B.P.,
! Teukolsky, S.A., and Vetterling, W.T. (Cambridge Univ. Press).

! NR version uses RAN1, but it was found that under RH7.1/pgf77 
! gasdev does not have <gasdev> = 0 when we use RAN1.
! Therefore have changed to use RAN3

! Modifications by B.T. Draine, Princeton University Observatory,
! 2003.09.10
! 2007.08.06 (BTD) converted to f90

!=======================================================================

      ZERO=0._WP
      ONE=1._WP
      TWO=2._WP

      IF(ISET==0)THEN

! when ISET = 0, we do not have an extra deviate on hand
!                so generate two new deviates

 0100   CONTINUE

! generate random variates V1 and V2 on interval -1,1

        V1=TWO*RAN3(IDUM)-ONE
        V2=TWO*RAN3(IDUM)-ONE
        RSQ=V1**2+V2**2
        IF(RSQ>=ONE.OR.RSQ==ZERO)GOTO 0100
        FAC=SQRT(-TWO*LOG(RSQ)/RSQ)
        GSET=V1*FAC
        GASDEV=V2*FAC
        ISET=1

      ELSE

! when ISET = 1, we have an extra deviate on hand
!                so use it, and change ISET back to 0

         GASDEV=GSET
         ISET=0
      ENDIF
      RETURN
      END

