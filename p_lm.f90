    FUNCTION P_LM(L,M,X)
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE

! arguments

      INTEGER :: L,M
      REAL(WP) :: P_LM,X

! local variables

      CHARACTER :: CMSGNM*70
      INTEGER :: I,LL
      REAL(WP) :: FACT,ONE,PLL,PMM,PMMP1,SINTH,TWO

!-----------------------------------------------------------------------
! FUNCTION P_LM  = Associated Legendre polynomial P_{lm}(x=cos(theta))
! given:
!    L,M
!    X
! returns:
!    P_LM = P_lm(x)
!             associated Legendre polynomial
! this routine follows closely the structure of function PLGNDR 
! described in section 6.6 of "Numerical Recipes", by Press, W.H.,
! Flannery, B.P., Teukolsky, S.A., and Vetterling, W.T. (1986, Cambridge
! Univ. Press).
!
! Modifications by B.T. Draine, Princeton Univ. Observatory
! history
! 02.11.12 (BTD) first written
! 04.03.31 (BTD) changed WRITE(0 to CALL WRIMSG(
! 07.08.06 (BTD) converted to f90
! end history
!-----------------------------------------------------------------------
      ONE=1._WP
      TWO=2._WP
      PMM=ONE
      IF(M>0)THEN
         SINTH=SQRT(ONE-X*X)
         FACT=ONE
         DO I=1,M
            PMM=PMM*FACT*SINTH
            FACT=FACT+TWO
         ENDDO
      ELSEIF(M<0.OR.M>=L)THEN
         WRITE(CMSGNM,'(A)')' Fatal error in P_LM: M=',M,' < 0'
         CALL WRIMSG(' P_LM ',CMSGNM)
         STOP
      ENDIF

      IF(L==M)THEN
         P_LM=PMM
      ELSE
         PMMP1=X*REAL(2*M+1)*PMM
         IF(L==M+1)THEN
            P_LM=PMMP1
         ELSE
            DO LL=M+2,L
               PLL=(X*REAL(2*LL-1)*PMMP1-REAL(LL+M-1)*PMM)/REAL(LL-M)
               PMM=PMMP1
               PMMP1=PLL
            ENDDO
            P_LM=PLL
         ENDIF
      ENDIF
      RETURN
      END
