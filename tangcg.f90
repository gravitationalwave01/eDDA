      SUBROUTINE TANGCG(TOL,MAXIT,XI,XR,B,MATVEC,WRK,LDA,  &
                        NDIM,NLAR,TOLE,NLOOP,NCOMPTE)
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE
! arguments:
      INTEGER :: LDA,MAXIT,NCOMPTE,NDIM,NLAR,NLOOP
      REAL(WP) :: TOL,TOLE
      COMPLEX(WP) :: B(LDA),WRK(LDA,NLAR),XI(LDA),XR(LDA) 
! Interface code to Tang et. al. conjugate gradient
! tol - tolerance to be achived
! maxit - maximum number of iterations
! xi - on input initial guess, on output vector x 
! xr - scratch array
! b - right hand side; i.e. Ax=b
! matvec - external subroutine calculating matrix times vector multiplication Ax
! wrk - scratch array wrk(lda, 12). NOTE that this has to be set properly
! lda - maximum leading dimension
! ndim - actuall dimension of the linear problem
! tole - achived tolerance
! nloop - number of iterations
! ncompte - number of  Ax multiplications (this is where calculations are expensive)
! History
! (PJF) = Piotr Jacek Flatau
! (PCC) = P. C. Chaumet 
! (AR)  = A. Rahmani
! April 2, 2010 Original code by (PCC). (AR)
! May 6, 2010  (PJF)  converted to Fortran90, changed interface,
!                     single/double precision kinds
!-------------------------------------------------------------------------- 
! correspondence with DDSCAT
! xi - initial guess of cxpol on input; cxpol  on output
! b  - cxe (right hand side, not changed)
! xr - A xi (use matvec)  work vector
! lda -  nat3
! ndim - nat3 ?
! nlar  work array dimension >= 12
! wrk - cxsc NOTE!!! THAT MXCXSC has to be set =12 in main DDSCAT(large!)
! maxit - itermx ? mxiter?
! nloop - itern
! tol -   tolr (input)
! tole -  achieved relative error
! ipar(12) - convergence;  ipar(12)=0 converged

!     .. Parameters ..

      INTEGER :: MULTIPLICATIONS

!     .. Local Scalars ..

      COMPLEX(WP) :: ALPHA,BETA,DZETA,ETA,ICOMP,R0RN
      REAL(WP) :: NORM,RESIDU
      INTEGER :: I,J,NOU,STATUS,STEPERR
      INTEGER :: IPAR

!     .. External Subroutines ..
      EXTERNAL GPBICG,MATVEC

!     .. Intrinsic Functions ..
      INTRINSIC ABS,REAL,SQRT
!     ..
!     initialization
      NLOOP=0
      NOU=0
      NCOMPTE=0
10    CONTINUE
      CALL GPBICG(XI,XR,B,LDA,NDIM,NLAR,NOU,WRK,NLOOP,MAXIT,TOL,NORM,ALPHA, &
                  BETA,ETA,DZETA,R0RN,STATUS,STEPERR,RESIDU)
      IF(STATUS.LT.0)THEN
         WRITE(*,FMT=*)'STOP NSTAT',STATUS,STEPERR
         STOP
      ENDIF
      NCOMPTE=NCOMPTE+1
      CALL MATVEC(XI,XR,IPAR)
      IF(STATUS.NE.1) GO TO 10

! finished iterations

      IF(STEPERR.EQ.0)THEN
        WRITE(*,FMT=*)'NLOOP HAS REACHED MAXIT',NLOOP,MAXIT
      ENDIF

!     Compute the relative error

      TOLE=0._WP
      DO I=1,NDIM
         TOLE=TOLE+ABS(XR(I)-B(I))**2._WP
      ENDDO
      TOLE=SQRT(TOLE)/NORM
    RETURN
    END

