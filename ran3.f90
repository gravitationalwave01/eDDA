    FUNCTION RAN3(IDUM)
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE

! Arguments:

      INTEGER :: IDUM
      REAL(WP) :: RAN3

! Local variables:

      INTEGER :: MBIG,MSEED,MZ
      REAL(WP) :: FAC
      PARAMETER(MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
      INTEGER :: I,IFF,II,INEXT,INEXTP,K
      INTEGER :: MJ,MK,MA(55)
      SAVE IFF,INEXT,INEXTP,MA

!***********************************************************************
! Function RAN3
! Given:
!     IDUM = integer seed (always used on first call, but
!            disregarded on subsequent calls unless IDUM < 0)
! Returns:
!     RAN3 = uniform random deviate between 0.0 and 1.0.
!
! Usage:
!     On first call, will be initialized with whatever seed
!     IDUM is provided.
!     On subsequent calls, result does not depend on IDUM 
!     unless IDUM<0, in which case IDUM is used to reinitialize
!     the sequence.
!
! This routine is based on Knuth's recommendation for a portable
! random number generator.  Present Fortran 90 implementation follows
! closely the f77 FUNCTION RAN3 described in section 7.1 of the book 
! "Numerical Recipes" by W.H. Press, B.P. Flannery, S.A. Teukolsky, and 
! W.T. Vetterling (1986; Cambridge University Press).
!
! B.T. Draine, Princeton University Observatory
!
!***********************************************************************

      DATA IFF/0/

! Initialization:

      IF(IDUM<0.OR.IFF==0)THEN
         IFF=1
         MJ=MSEED-IABS(IDUM)
         MJ=MOD(MJ,MBIG)
         MA(55)=MJ
         MK=1
         DO I=1,54
            II=MOD(21*I,55)
            MA(II)=MK
            MK=MJ-MK
            IF(MK.LT.MZ)MK=MK+MBIG
            MJ=MA(II)
         ENDDO
         DO K=1,4
            DO I=1,55
               MA(I)=MA(I)-MA(1+MOD(I+30,55))
               IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
            ENDDO
         ENDDO
         INEXT=0
         INEXTP=31
         IDUM=1
      ENDIF

! Under normal use (second and subsequent calls), this is the
! first executed statement:

      INEXT=INEXT+1
      IF(INEXT==56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP==56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ<MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAN3=MJ*FAC

      RETURN
    END FUNCTION RAN3
