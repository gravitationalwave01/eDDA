    SUBROUTINE PIMCBICG(X,B,WRK,IPAR,SPAR,MATVEC,TMATVEC,PRECONL,PRECONR, &
        PCSUM,PSCNRM,PROGRESS)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

!           PIM -- The Parallel Iterative Methods package
!           ---------------------------------------------

!                      Rudnei Dias da Cunha
!     National Supercomputing Centre and Mathematics Institute
!         Universidade Federal do Rio Grande do Sul, Brasil

!                          Tim Hopkins
!     Computing Laboratory, University of Kent at Canterbury, U.K.

! ----------------------------------------------------------------------

!     .. Parameters ..
      REAL (WP) :: ONE
      PARAMETER (ONE=1.0E0_WP)
      COMPLEX (WP) :: CZERO
      PARAMETER (CZERO=(0.0E0_WP,0.0E0_WP))
      COMPLEX (WP) :: CONE
      PARAMETER (CONE=(1.0E0_WP,0.0E0_WP))
      INTEGER :: IPARSIZ
      PARAMETER (IPARSIZ=13)
      INTEGER :: SPARSIZ
!Art      PARAMETER (SPARSIZ=2)
      PARAMETER (SPARSIZ=6)
!     ..
!     .. Array Arguments ..
      COMPLEX (WP) :: B(*), WRK(*), X(*)
      REAL (WP) :: SPAR(SPARSIZ)
      INTEGER :: IPAR(IPARSIZ)
!     ..
!     .. Function Arguments ..
      REAL (WP) :: PSCNRM
      EXTERNAL PSCNRM
!     ..
!     .. Subroutine Arguments ..
      EXTERNAL MATVEC, PCSUM, PRECONL, PRECONR, PROGRESS, TMATVEC
!     ..
!     .. Local Scalars ..
      COMPLEX (WP) :: ALPHA, BETA, DELTA, RHO, RHO0, XI
      REAL (WP) :: EPSILON, EXITNORM, RHSSTOP
      INTEGER :: BASISDIM, BLKSZ, CNVRTX, IP, IPTILDE, IR, IRTILDE, IS, ITNO, &
        IW, IXOLD, IZ, LDA, LOCLEN, MAXIT, N, NPROCS, PRECONTYPE, PROCID, &
        STATUS, STEPERR, STOPTYPE
!     ..
!     .. Local Arrays ..
      COMPLEX (WP) :: DOTS(2)
!     ..
!     .. External Functions ..
      COMPLEX (WP) :: CDOTC
      REAL (WP) :: SCSETRHSSTOP
      EXTERNAL CDOTC, SCSETRHSSTOP
!     ..
!     .. External Subroutines ..
      EXTERNAL CAXPY, CCOPY, PIMSGETPAR, STOPCRIT
!     ..
      CALL PIMSGETPAR(IPAR,SPAR,LDA,N,BLKSZ,LOCLEN,BASISDIM,NPROCS,PROCID, &
        PRECONTYPE,STOPTYPE,MAXIT,ITNO,STATUS,STEPERR,EPSILON,EXITNORM)

!  Check consistency of preconditioning and stop types
      IF (((PRECONTYPE==0) .OR. (PRECONTYPE==2)) .AND. (STOPTYPE==6)) THEN
        ITNO = 0
        STATUS = -4
        STEPERR = 0
        GO TO 9999

      END IF

!  Does not need conversion Y=Q2X for residual
      CNVRTX = 0

!  Set indices for mapping local vectors into wrk
      IR = 1
      IRTILDE = IR + LOCLEN
      IPTILDE = IRTILDE + LOCLEN
      IP = IPTILDE + LOCLEN
      IW = IP + LOCLEN
      IZ = IW + LOCLEN
      IS = IZ + LOCLEN
      IXOLD = IS + LOCLEN

!  Set rhs of stopping criteria
      RHSSTOP = SCSETRHSSTOP(B,WRK(IR),EPSILON,IPAR,PRECONL,PSCNRM)

!  1. r=Q1(b-AQ2x)
      IF (STOPTYPE/=6) THEN
        IF (PRECONTYPE==0) THEN
!     r=b-Ax
          CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==1) THEN
!     r=Q1(b-Ax)
          CALL CCOPY(LOCLEN,B,1,WRK(IZ),1)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IZ),1)
          CALL PRECONL(WRK(IZ),WRK(IR),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
!     r=b-AQ2x
          CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL PRECONR(X,WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(b-AQ2x)
          CALL CCOPY(LOCLEN,B,1,WRK(IP),1)
          CALL PRECONR(X,WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IP),1)
          CALL PRECONL(WRK(IP),WRK(IR),IPAR)
        END IF

      ELSE
!     r has been set to Qb in the call to ssetrhsstop
        IF (PRECONTYPE==1) THEN
!     r=Q1(b-Ax)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL PRECONL(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(b-AQ2x)
          CALL PRECONR(X,WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IW),IPAR)
          CALL PRECONL(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)
        END IF

      END IF

!  2. rtilde=ptilde=p=r
      CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IRTILDE),1)
      CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IPTILDE),1)
      CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IP),1)

!  3. rho=dot(rtilde,r)
      DOTS(1) = CDOTC(LOCLEN,WRK(IRTILDE),1,WRK(IR),1)
      CALL PCSUM(1,DOTS)
      RHO = DOTS(1)

!  4. w=Q1AQ2p
      IF (PRECONTYPE==0) THEN
        CALL MATVEC(WRK(IP),WRK(IW),IPAR)

      ELSE IF (PRECONTYPE==1) THEN
        CALL MATVEC(WRK(IP),WRK(IZ),IPAR)
        CALL PRECONL(WRK(IZ),WRK(IW),IPAR)

      ELSE IF (PRECONTYPE==2) THEN
        CALL PRECONR(WRK(IP),WRK(IZ),IPAR)
        CALL MATVEC(WRK(IZ),WRK(IW),IPAR)

      ELSE IF (PRECONTYPE==3) THEN
        CALL PRECONR(WRK(IP),WRK(IW),IPAR)
        CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
        CALL PRECONL(WRK(IZ),WRK(IW),IPAR)
      END IF

!  5. xi=dot(ptilde,w)
      DOTS(1) = CDOTC(LOCLEN,WRK(IPTILDE),1,WRK(IW),1)
      CALL PCSUM(1,DOTS)
      XI = DOTS(1)

!  Loop
      STATUS = 0
      EXITNORM = -ONE
      STEPERR = -1
      DO ITNO = 1, MAXIT

!  6. alpha=rho/xi
        IF (XI==CZERO) THEN
          STATUS = -3
          STEPERR = 6
          GO TO 9999

        END IF

        ALPHA = RHO/XI

!  7. x=x+alpha*p
        CALL CCOPY(LOCLEN,X,1,WRK(IXOLD),1)
        CALL CAXPY(LOCLEN,ALPHA,WRK(IP),1,X,1)

!  8. r=r-alpha*w
        CALL CAXPY(LOCLEN,-ALPHA,WRK(IW),1,WRK(IR),1)

!  9. check stopping criterion
        CALL STOPCRIT(B,WRK(IR),WRK(IZ),X,WRK(IXOLD),WRK(IW),RHSSTOP,CNVRTX, &
          EXITNORM,STATUS,IPAR,MATVEC,TMATVEC,PRECONR,PCSUM,PSCNRM)

!  Call monitoring routine
        CALL PROGRESS(LOCLEN,ITNO,EXITNORM,X,WRK(IR),WRK(IZ))

        IF (STATUS==0) THEN
          GO TO 9999
        END IF
! 10. rtilde=rtilde-alpha*Q1A^{T}Q2ptilde
        IF (PRECONTYPE==0) THEN
          CALL TMATVEC(WRK(IPTILDE),WRK(IW),IPAR)

        ELSE IF (PRECONTYPE==1) THEN
          CALL TMATVEC(WRK(IPTILDE),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IW),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
          CALL PRECONR(WRK(IPTILDE),WRK(IZ),IPAR)
          CALL TMATVEC(WRK(IZ),WRK(IW),IPAR)

        ELSE IF (PRECONTYPE==3) THEN
          CALL PRECONR(WRK(IPTILDE),WRK(IW),IPAR)
          CALL TMATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IW),IPAR)
        END IF

        CALL CAXPY(LOCLEN,-ALPHA,WRK(IW),1,WRK(IRTILDE),1)

! 11. s=Q1AQ2r
        IF (PRECONTYPE==0) THEN
          CALL MATVEC(WRK(IR),WRK(IS),IPAR)

        ELSE IF (PRECONTYPE==1) THEN
          CALL MATVEC(WRK(IR),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IS),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
          CALL PRECONR(WRK(IR),WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IS),IPAR)

        ELSE IF (PRECONTYPE==3) THEN
          CALL PRECONR(WRK(IR),WRK(IS),IPAR)
          CALL MATVEC(WRK(IS),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IS),IPAR)
        END IF

! 12. rho=dot(rtilde,r)
        RHO0 = RHO
        DOTS(1) = CDOTC(LOCLEN,WRK(IRTILDE),1,WRK(IR),1)

! 13. delta=dot(rtilde,s)
        DOTS(2) = CDOTC(LOCLEN,WRK(IRTILDE),1,WRK(IS),1)

!  Accumulate simultaneously partial values
        CALL PCSUM(2,DOTS)
        RHO = DOTS(1)
        DELTA = DOTS(2)

! 14. beta=rho/rho0
        IF (RHO0==CZERO) THEN
          STATUS = -3
          STEPERR = 14
          GO TO 9999

        END IF

        BETA = RHO/RHO0

! 15. p=r+beta*p
        CALL CCOPY(LOCLEN,WRK(IP),1,WRK(IZ),1)
        CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IP),1)
        CALL CAXPY(LOCLEN,BETA,WRK(IZ),1,WRK(IP),1)

! 16. ptilde=rtilde+beta*ptilde
        CALL CCOPY(LOCLEN,WRK(IPTILDE),1,WRK(IZ),1)
        CALL CCOPY(LOCLEN,WRK(IRTILDE),1,WRK(IPTILDE),1)
        CALL CAXPY(LOCLEN,BETA,WRK(IZ),1,WRK(IPTILDE),1)

! 17. w=s+beta*w
        CALL CCOPY(LOCLEN,WRK(IW),1,WRK(IZ),1)
        CALL CCOPY(LOCLEN,WRK(IS),1,WRK(IW),1)
        CALL CAXPY(LOCLEN,BETA,WRK(IZ),1,WRK(IW),1)

! 18. xi=delta-beta^2*xi
        XI = DELTA - BETA**2*XI

      END DO

      IF (ITNO>MAXIT) THEN
        STATUS = -1
        ITNO = MAXIT
      END IF

9999  CONTINUE

      IF ((PRECONTYPE==2) .OR. (PRECONTYPE==3)) THEN
        CALL CCOPY(LOCLEN,X,1,WRK(IZ),1)
        CALL PRECONR(WRK(IZ),X,IPAR)
      END IF

!  Set output parameters
      IPAR(11) = ITNO
      IPAR(12) = STATUS
      IPAR(13) = STEPERR
      SPAR(2) = EXITNORM

      RETURN

    END SUBROUTINE PIMCBICG

    SUBROUTINE PIMCBICGSTAB(X,B,WRK,IPAR,SPAR,MATVEC,PRECONL,PRECONR,PCSUM, &
        PSCNRM,PROGRESS)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

!           PIM -- The Parallel Iterative Methods package
!           ---------------------------------------------

!                      Rudnei Dias da Cunha
!     National Supercomputing Centre and Mathematics Institute
!         Universidade Federal do Rio Grande do Sul, Brasil

!                          Tim Hopkins
!     Computing Laboratory, University of Kent at Canterbury, U.K.

! ----------------------------------------------------------------------

!     .. Parameters ..
      REAL (WP) :: ONE
      PARAMETER (ONE=1.0E0_WP)
      COMPLEX (WP) :: CZERO
      PARAMETER (CZERO=(0.0E0_WP,0.0E0_WP))
      COMPLEX (WP) :: CONE
      PARAMETER (CONE=(1.0E0_WP,0.0E0_WP))
      INTEGER :: IPARSIZ
      PARAMETER (IPARSIZ=13)
      INTEGER :: SPARSIZ
!Art      PARAMETER (SPARSIZ=2)
      PARAMETER (SPARSIZ=6)
!     ..
!     .. Array Arguments ..
      COMPLEX (WP) :: B(*), WRK(*), X(*)
      REAL (WP) :: SPAR(SPARSIZ)
      INTEGER :: IPAR(IPARSIZ)
!     ..
!     .. Function Arguments ..
      REAL (WP) :: PSCNRM
      EXTERNAL PSCNRM
!     ..
!     .. Subroutine Arguments ..
      EXTERNAL MATVEC, PCSUM, PRECONL, PRECONR, PROGRESS
!     ..
!     .. Local Scalars ..
      COMPLEX (WP) :: ALPHA, BETA, KAPPA, OMEGA, RHO, RHO0, XI
      REAL (WP) :: EPSILON, EXITNORM, MACHEPS, RHSSTOP
      INTEGER :: BASISDIM, BLKSZ, CNVRTX, IP, IQ, IR, IRTILDE, IS, IT, ITNO, &
        IV, IW, IXOLD, IZ, LDA, LOCLEN, MAXIT, N, NPROCS, PRECONTYPE, PROCID, &
        STATUS, STEPERR, STOPTYPE
!     ..
!     .. Local Arrays ..
      COMPLEX (WP) :: DOTS(2)
!     ..
!     .. External Functions ..
      COMPLEX (WP) :: CDOTC
      REAL (WP) :: SCSETRHSSTOP
      EXTERNAL CDOTC, SCSETRHSSTOP
!     ..
!     .. External Subroutines ..
      EXTERNAL CAXPY, CCOPY, CINIT, PIMSGETPAR, SMACHCONS, STOPCRIT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS
!     ..
!*** diagnostic
!      write(0,*)'entered pibcbigcstab'
!      write(0,*)'about to call smachcons'
!****

      CALL SMACHCONS('M',MACHEPS)

!*** diagnostic
!      write(0,*)'returned from smachcons with MACHEPS=',macheps
!      write(0,*)'about to call pimsgetpar with'
!      write(0,*)'ipar=',ipar
!      write(0,*)'spar=',spar
!      write(0,*)'N=',n
!***

      CALL PIMSGETPAR(IPAR,SPAR,LDA,N,BLKSZ,LOCLEN,BASISDIM,NPROCS,PROCID, &
        PRECONTYPE,STOPTYPE,MAXIT,ITNO,STATUS,STEPERR,EPSILON,EXITNORM)

!*** diagnostic
!      write(0,*)'returned from pisgetpar'
!      write(0,*)'ipar=',ipar
!      write(0,*)'spar=',spar
!      write(0,*)'precontype=',precontype
!      write(0,*)'stoptype=',stoptype
!***

!Check consistency of preconditioning and stop types
      IF (((PRECONTYPE==0) .OR. (PRECONTYPE==2)) .AND. (STOPTYPE==6)) THEN
        ITNO = 0
        STATUS = -4
        STEPERR = 0
        GO TO 9999

      END IF

!  Does not need conversion Y=Q2X for residual
      CNVRTX = 0

!  Set indices for mapping local vectors into wrk
      IR = 1
      IRTILDE = IR + LOCLEN
      IP = IRTILDE + LOCLEN
      IQ = IP + LOCLEN
      IS = IQ + LOCLEN
      IT = IS + LOCLEN
      IV = IT + LOCLEN
      IW = IV + LOCLEN
      IZ = IW + LOCLEN
      IXOLD = IZ + LOCLEN

!  Set rhs of stopping criteria
      RHSSTOP = SCSETRHSSTOP(B,WRK(IR),EPSILON,IPAR,PRECONL,PSCNRM)
!*** diagnostic
!      write(0,*)'pimcbigstab ckpt delta, rhsstop=',rhsstop
!***

!  1. r=Q1(b-AQ2x)
      IF (STOPTYPE/=6) THEN
        IF (PRECONTYPE==0) THEN
!     r=b-Ax
          CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==1) THEN
!     r=Q1(b-Ax)
!*** diagnostic
!             write(0,*)'pimcbigstab ckpt delta.2, about to call ccopy'
!             write(0,*)'loclen=',loclen,' IZ=',IZ
!             write(0,*)'b(1)=',b(1)
!             write(0,*)'b(17000)=',b(17000)
!***
          CALL CCOPY(LOCLEN,B,1,WRK(IZ),1)
!*** diagnostic
!              write(0,*)'about to call matvec with x(17000)=',x(17000)
!***
          CALL MATVEC(X,WRK(IW),IPAR)
!*** diagnostic
!              write(0,*)'ckpt delta.3, returned from matvec with'
!              write(0,*)'iw=',iw
!              write(0,*)'wrk(iw+1000)=',wrk(iw+1000)
!              write(0,*)'wrk(iw+2000)=',wrk(iw+2000)
!              write(0,*)'wrk(iw+3000)=',wrk(iw+3000)
!              write(0,*)'wrk(iw+4000)=',wrk(iw+4000)
!              write(0,*)'wrk(iw+5000)=',wrk(iw+5000)
!              write(0,*)'wrk(iw+6000)=',wrk(iw+6000)
!              write(0,*)'wrk(iw+7000)=',wrk(iw+7000)
!              write(0,*)'wrk(iw+8000)=',wrk(iw+8000)
!              write(0,*)'wrk(iw+9000)=',wrk(iw+9000)
!              write(0,*)'wrk(iw+10000)=',wrk(iw+10000)
!              write(0,*)'wrk(iw+17000)=',wrk(iw+17000)
!              write(0,*)'x(1000)=',x(1000)
!              write(0,*)'x(2000)=',x(2000)
!              write(0,*)'x(3000)=',x(3000)
!              write(0,*)'x(4000)=',x(4000)
!              write(0,*)'x(5000)=',x(5000)
!              write(0,*)'x(6000)=',x(6000)
!              write(0,*)'x(7000)=',x(7000)
!              write(0,*)'x(8000)=',x(8000)
!              write(0,*)'x(9000)=',x(9000)
!***
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IZ),1)
!*** diagnostic
!          write(0,*)'ckpt n1: about to call PRECONL with IZ=',IZ,' IR=',IR
!          write(0,*)'wrk(iw+1)=',wrk(iw+1)
!          write(0,*)'wrk(iz+1)=',wrk(iz+1)
!***
          CALL PRECONL(WRK(IZ),WRK(IR),IPAR)
!***
!          write(0,*)'returned from PRECONL'
!***

        ELSE IF (PRECONTYPE==2) THEN
!     r=b-AQ2x
          CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL PRECONR(X,WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(b-AQ2x)
          CALL CCOPY(LOCLEN,B,1,WRK(IP),1)
          CALL PRECONR(X,WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IP),1)
          CALL PRECONL(WRK(IP),WRK(IR),IPAR)
        END IF

      ELSE
!     r has been set to Qb in the call to dsetrhsstop
        IF (PRECONTYPE==1) THEN
!     r=Q1(b-Ax)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL PRECONL(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(b-AQ2x)
          CALL PRECONR(X,WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IW),IPAR)
          CALL PRECONL(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)
        END IF

      END IF

!  2. rtilde=r
      CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IRTILDE),1)

!  3. p=v=0
      CALL CINIT(LOCLEN,CZERO,WRK(IP),1)
      CALL CINIT(LOCLEN,CZERO,WRK(IV),1)

!  4. rho=alpha=omega=1
      RHO = CONE
      ALPHA = CONE
      OMEGA = CONE

!  Loop
      STATUS = 0
      EXITNORM = -ONE
      STEPERR = -1
      DO ITNO = 1, MAXIT

!  5. rho=dot(rtilde,r)
        RHO0 = RHO
        DOTS(1) = CDOTC(LOCLEN,WRK(IRTILDE),1,WRK(IR),1)
!*** diagnostic
!          write(0,*)'itno=',itno,' dots(1)=',dots(1)
!***
        CALL PCSUM(1,DOTS)
        RHO = DOTS(1)

!  6. beta=rho*alpha/(rho0*omega)
        KAPPA = RHO0*OMEGA
!*** diagnostic
!          write(0,*)'rho0=',rho0
!          write(0,*)'omega=',omega
!          write(0,*)'kappa=',kappa
!***
        IF (KAPPA==CZERO) THEN
          STATUS = -3
          STEPERR = 6
          GO TO 9999

        END IF

        BETA = RHO*ALPHA/KAPPA

!*** diagnostic
!          write(0,*)'ckpt m100, alpha=',alpha,' beta=',beta
!***
!  7. p=r+beta*(p-omega*v)
        CALL CAXPY(LOCLEN,-OMEGA,WRK(IV),1,WRK(IP),1)
        CALL CCOPY(LOCLEN,WRK(IP),1,WRK(IW),1)
        CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IP),1)
        CALL CAXPY(LOCLEN,BETA,WRK(IW),1,WRK(IP),1)
!*** diagnostic
!          write(0,*)'ckpt m101, precontype=',precontype
!***
!  8. v=Q1AQ2p
        IF (PRECONTYPE==0) THEN
          CALL MATVEC(WRK(IP),WRK(IV),IPAR)

        ELSE IF (PRECONTYPE==1) THEN
!*** diagnostic
!             write(0,*)'ckpt m201, about to call matvec'
!             write(0,*)'with ipar=',ipar
!             write(0,*)'ip=',ip
!             write(0,*)'wrk(ip+1)=',wrk(ip+1)
!             write(0,*)'iz=',iz
!             write(0,*)'wrk(iz+1)=',wrk(iz+1)
! ************ wrk(iz+1) is ok at this point
!              problem must be in next call to matvec:
!***
          CALL MATVEC(WRK(IP),WRK(IZ),IPAR)
!*** diagnostic
!              write(0,*)'returned from matvec with'
!              write(0,*)'ipar=',ipar
!              write(0,*)'ip=',ip,' iz=',iz
!              write(0,*)'wrk(iw+1000)=',wrk(iw+1000)
!              write(0,*)'wrk(iw+2000)=',wrk(iw+2000)
!              write(0,*)'wrk(iw+3000)=',wrk(iw+3000)
!              write(0,*)'wrk(iw+4000)=',wrk(iw+4000)
!              write(0,*)'wrk(iw+5000)=',wrk(iw+5000)
!              write(0,*)'wrk(iw+6000)=',wrk(iw+6000)
!              write(0,*)'wrk(iw+7000)=',wrk(iw+7000)
!              write(0,*)'wrk(iw+8000)=',wrk(iw+8000)
!              write(0,*)'wrk(iw+9000)=',wrk(iw+9000)
!              write(0,*)'wrk(iw+10000)=',wrk(iw+10000)
!              write(0,*)'wrk(iw+17000)=',wrk(iw+17000)
!              write(0,*)'wrk(589825)=',wrk(589825)
!              write(0,*)'wkr(iz+1)=',wrk(iz+1)
!***
!*** diagnostic
!              write(0,*)'ckpt n2: about to call PRECONL with iz=',iz,' iv=',iv
! ************** problem occurs before this point: wrk(iz+1) is NaN
!                where iz=786433
!              write(0,*)'wrk(iz+1)=',wrk(iz+1)
!              write(0,*)'wkr(iv+1)=',wrk(iv+1)
!              write(0,*)'call PRECONL'
!***
          CALL PRECONL(WRK(IZ),WRK(IV),IPAR)

! PRECONL has messed with wrk(589825)...
!*** diagnostic
!              write(0,*)'ckpt m210'
!              write(0,*)'wrk(589825)=',wrk(589825)
!***
        ELSE IF (PRECONTYPE==2) THEN
          CALL PRECONR(WRK(IP),WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IV),IPAR)

        ELSE IF (PRECONTYPE==3) THEN
          CALL PRECONR(WRK(IP),WRK(IV),IPAR)
          CALL MATVEC(WRK(IV),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IV),IPAR)
        END IF

!  9. xi=dot(rtilde,v)
!*** diagnostic
!          write(0,*)'ckpt m230, loclen=',loclen
!          write(0,*)'irtilde=',irtilde
!          write(0,*)'wrk(irtilde)=',wrk(irtilde)
!          write(0,*)'iv=',iv
!          write(0,*)'wrk(iv)=',wrk(iv)
!***
        DOTS(1) = CDOTC(LOCLEN,WRK(IRTILDE),1,WRK(IV),1)
!*** diagnostic
!          write(0,*)'ckpt m250, dots(1)=',dots(1)
!***
        CALL PCSUM(1,DOTS)
        XI = DOTS(1)
!*** diagnostic
!          write(0,*)'ckpt m270, xi=',xi
!***
! 10. alpha=rho/xi
        IF (XI==CZERO) THEN
          STATUS = -3
          STEPERR = 10
          GO TO 9999

        END IF

        ALPHA = RHO/XI

! 11. s=r-alpha*v
        CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IS),1)
        CALL CAXPY(LOCLEN,-ALPHA,WRK(IV),1,WRK(IS),1)

! 12. if ||s||<breaktol then soft-breakdown has occurred
        KAPPA = PSCNRM(LOCLEN,WRK(IS))
        IF (ABS(KAPPA)<MACHEPS) THEN
          STATUS = -2
          STEPERR = 12
          CALL CAXPY(LOCLEN,ALPHA,WRK(IP),1,X,1)
          GO TO 9999

        END IF

! 13. t=Q1AQ2s
        IF (PRECONTYPE==0) THEN
          CALL MATVEC(WRK(IS),WRK(IT),IPAR)

        ELSE IF (PRECONTYPE==1) THEN
          CALL MATVEC(WRK(IS),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IT),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
          CALL PRECONR(WRK(IS),WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IT),IPAR)

        ELSE IF (PRECONTYPE==3) THEN
          CALL PRECONR(WRK(IS),WRK(IT),IPAR)
          CALL MATVEC(WRK(IT),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IT),IPAR)
        END IF

! 14. omega=dot(t,s)/dot(t,t)
        DOTS(1) = CDOTC(LOCLEN,WRK(IT),1,WRK(IT),1)
        DOTS(2) = CDOTC(LOCLEN,WRK(IT),1,WRK(IS),1)

!  Accumulate simultaneously partial values
        CALL PCSUM(2,DOTS)

        IF (DOTS(1)==CZERO) THEN
          STATUS = -3
          STEPERR = 14
          GO TO 9999

        END IF

        OMEGA = DOTS(2)/DOTS(1)

! 15. x=x+alpha*p+omega*s
        CALL CCOPY(LOCLEN,X,1,WRK(IXOLD),1)
        CALL CAXPY(LOCLEN,ALPHA,WRK(IP),1,X,1)
        CALL CAXPY(LOCLEN,OMEGA,WRK(IS),1,X,1)

! 16. r=s-omega*t
        CALL CCOPY(LOCLEN,WRK(IS),1,WRK(IR),1)
        CALL CAXPY(LOCLEN,-OMEGA,WRK(IT),1,WRK(IR),1)

! 17. check stopping criterion
        CALL STOPCRIT(B,WRK(IR),WRK(IZ),X,WRK(IXOLD),WRK(IW),RHSSTOP,CNVRTX, &
          EXITNORM,STATUS,IPAR,MATVEC,MATVEC,PRECONR,PCSUM,PSCNRM)

!  Call monitoring routine
        CALL PROGRESS(LOCLEN,ITNO,EXITNORM,X,WRK(IR),WRK(IZ))

        IF (STATUS==0) THEN
          GO TO 9999
        END IF
      END DO

      IF (ITNO>MAXIT) THEN
        STATUS = -1
        ITNO = MAXIT
      END IF

9999  CONTINUE

      IF ((PRECONTYPE==2) .OR. (PRECONTYPE==3)) THEN
        CALL CCOPY(LOCLEN,X,1,WRK(IZ),1)
        CALL PRECONR(WRK(IZ),X,IPAR)
      END IF

!  Set output parameters
      IPAR(11) = ITNO
      IPAR(12) = STATUS
      IPAR(13) = STEPERR
      SPAR(2) = EXITNORM

      RETURN

    END SUBROUTINE PIMCBICGSTAB

    SUBROUTINE PIMCCG(X,B,WRK,IPAR,SPAR,MATVEC,PRECONL,PRECONR,PCSUM,PSCNRM, &
        PROGRESS)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

!           PIM -- The Parallel Iterative Methods package
!           ---------------------------------------------

!                      Rudnei Dias da Cunha
!     National Supercomputing Centre and Mathematics Institute
!         Universidade Federal do Rio Grande do Sul, Brasil

!                          Tim Hopkins
!     Computing Laboratory, University of Kent at Canterbury, U.K.

! ----------------------------------------------------------------------

!     .. Parameters ..
      REAL (WP) :: ONE
      PARAMETER (ONE=1.0E0_WP)
      COMPLEX (WP) :: CZERO
      PARAMETER (CZERO=(0.0E0_WP,0.0E0_WP))
      COMPLEX (WP) :: CONE
      PARAMETER (CONE=(1.0E0_WP,0.0E0_WP))
      INTEGER :: IPARSIZ
      PARAMETER (IPARSIZ=13)
      INTEGER :: SPARSIZ
!Art      PARAMETER (SPARSIZ=2)
      PARAMETER (SPARSIZ=6)
!     ..
!     .. Array Arguments ..
      COMPLEX (WP) :: B(*), WRK(*), X(*)
      REAL (WP) :: SPAR(SPARSIZ)
      INTEGER :: IPAR(IPARSIZ)
!     ..
!     .. Function Arguments ..
      REAL (WP) :: PSCNRM
      EXTERNAL PSCNRM
!     ..
!     .. Subroutine Arguments ..
      EXTERNAL MATVEC, PCSUM, PRECONL, PRECONR, PROGRESS
!     ..
!     .. Local Scalars ..
      COMPLEX (WP) :: ALPHA, BETA, DELTA, RDOTR, RDOTR0, XI
      REAL (WP) :: EPSILON, EXITNORM, RHSSTOP
      INTEGER :: BASISDIM, BLKSZ, CNVRTX, IP, IR, IS, ITNO, IW, IXOLD, IZ, &
        LDA, LOCLEN, MAXIT, N, NPROCS, PRECONTYPE, PROCID, STATUS, STEPERR, &
        STOPTYPE
!     ..
!     .. Local Arrays ..
      COMPLEX (WP) :: DOTS(2)
!     ..
!     .. External Functions ..
      COMPLEX (WP) :: CDOTC
      REAL (WP) :: SCSETRHSSTOP
      EXTERNAL CDOTC, SCSETRHSSTOP
!     ..
!     .. External Subroutines ..
      EXTERNAL CAXPY, CCOPY, PIMSGETPAR, STOPCRIT
!     ..

      CALL PIMSGETPAR(IPAR,SPAR,LDA,N,BLKSZ,LOCLEN,BASISDIM,NPROCS,PROCID, &
        PRECONTYPE,STOPTYPE,MAXIT,ITNO,STATUS,STEPERR,EPSILON,EXITNORM)

!  Check consistency of preconditioning and stop types
      IF (((PRECONTYPE==0) .OR. (PRECONTYPE==2)) .AND. (STOPTYPE==6)) THEN
        ITNO = 0
        STATUS = -4
        STEPERR = 0
        GO TO 9999

      END IF

!  Does not need conversion Y=Q2X for residual
      CNVRTX = 0

!  Set indices for mapping local vectors into wrk
      IR = 1
      IP = IR + LOCLEN
      IW = IP + LOCLEN
      IZ = IW + LOCLEN
      IS = IZ + LOCLEN
      IXOLD = IS + LOCLEN

!  Set rhs of stopping criteria
      RHSSTOP = SCSETRHSSTOP(B,WRK(IR),EPSILON,IPAR,PRECONL,PSCNRM)

!  1. r=Q1(b-AQ2x)
      IF (STOPTYPE/=6) THEN
        IF (PRECONTYPE==0) THEN
!     r=b-Ax
          CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==1) THEN
!     r=Q1(b-Ax)
          CALL CCOPY(LOCLEN,B,1,WRK(IZ),1)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IZ),1)
          CALL PRECONL(WRK(IZ),WRK(IR),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
!     r=b-AQ2x
          CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL PRECONR(X,WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(b-AQ2x)
          CALL CCOPY(LOCLEN,B,1,WRK(IP),1)
          CALL PRECONR(X,WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IP),1)
          CALL PRECONL(WRK(IP),WRK(IR),IPAR)
        END IF

      ELSE
!     r has been set to Qb in the call to dsetrhsstop
        IF (PRECONTYPE==1) THEN
!     r=Q1(b-Ax)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL PRECONL(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(b-AQ2x)
          CALL PRECONR(X,WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IW),IPAR)
          CALL PRECONL(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)
        END IF

      END IF

!  2. p=r
      CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IP),1)

!  3. rdotr=dot(r,r)
      DOTS(1) = CDOTC(LOCLEN,WRK(IR),1,WRK(IR),1)
      CALL PCSUM(1,DOTS)
      RDOTR = DOTS(1)

!  4. w=Q1AQ2p
      IF (PRECONTYPE==0) THEN
        CALL MATVEC(WRK(IP),WRK(IW),IPAR)

      ELSE IF (PRECONTYPE==1) THEN
        CALL MATVEC(WRK(IP),WRK(IZ),IPAR)
        CALL PRECONL(WRK(IZ),WRK(IW),IPAR)

      ELSE IF (PRECONTYPE==2) THEN
        CALL PRECONR(WRK(IP),WRK(IZ),IPAR)
        CALL MATVEC(WRK(IZ),WRK(IW),IPAR)

      ELSE IF (PRECONTYPE==3) THEN
        CALL PRECONR(WRK(IP),WRK(IW),IPAR)
        CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
        CALL PRECONL(WRK(IZ),WRK(IW),IPAR)
      END IF

!  5. xi=dot(p,w)
      DOTS(1) = CDOTC(LOCLEN,WRK(IP),1,WRK(IW),1)
      CALL PCSUM(1,DOTS)
      XI = DOTS(1)

!  Loop
      STATUS = 0
      EXITNORM = -ONE
      STEPERR = -1
      DO ITNO = 1, MAXIT

!  6. alpha=rdotr/xi
        IF (XI==CZERO) THEN
          STATUS = -3
          STEPERR = 6
          GO TO 9999

        END IF

        ALPHA = RDOTR/XI

!  7. x=x+alpha*p
        CALL CCOPY(LOCLEN,X,1,WRK(IXOLD),1)
        CALL CAXPY(LOCLEN,ALPHA,WRK(IP),1,X,1)

!  8. r=r-alpha*w
        CALL CAXPY(LOCLEN,-ALPHA,WRK(IW),1,WRK(IR),1)

!  9. check stopping criterion
        CALL STOPCRIT(B,WRK(IR),WRK(IZ),X,WRK(IXOLD),WRK(IS),RHSSTOP,CNVRTX, &
          EXITNORM,STATUS,IPAR,MATVEC,MATVEC,PRECONR,PCSUM,PSCNRM)

!  Call monitoring routine
        CALL PROGRESS(LOCLEN,ITNO,EXITNORM,X,WRK(IR),WRK(IZ))

        IF (STATUS==0) THEN
          GO TO 9999
        END IF

! 10. s=Q1AQ2r
        IF (PRECONTYPE==0) THEN
          CALL MATVEC(WRK(IR),WRK(IS),IPAR)

        ELSE IF (PRECONTYPE==1) THEN
          CALL MATVEC(WRK(IR),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IS),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
          CALL PRECONR(WRK(IR),WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IS),IPAR)

        ELSE IF (PRECONTYPE==3) THEN
          CALL PRECONR(WRK(IR),WRK(IS),IPAR)
          CALL MATVEC(WRK(IS),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IS),IPAR)
        END IF

! 11. rdotr=dot(r,r)
        RDOTR0 = RDOTR
        DOTS(1) = CDOTC(LOCLEN,WRK(IR),1,WRK(IR),1)

! 12. delta=dot(r,s)
        DOTS(2) = CDOTC(LOCLEN,WRK(IR),1,WRK(IS),1)

!  Accumulate simultaneously partial values
        CALL PCSUM(2,DOTS)
        RDOTR = DOTS(1)
        DELTA = DOTS(2)

! 13. beta=rdotr/rdotr0
        IF (RDOTR0==CZERO) THEN
          STATUS = -3
          STEPERR = 13
          GO TO 9999

        END IF

        BETA = RDOTR/RDOTR0

! 14. p=r+beta*p
        CALL CCOPY(LOCLEN,WRK(IP),1,WRK(IZ),1)
        CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IP),1)
        CALL CAXPY(LOCLEN,BETA,WRK(IZ),1,WRK(IP),1)

! 15. w=s+beta*w
        CALL CCOPY(LOCLEN,WRK(IW),1,WRK(IZ),1)
        CALL CCOPY(LOCLEN,WRK(IS),1,WRK(IW),1)
        CALL CAXPY(LOCLEN,BETA,WRK(IZ),1,WRK(IW),1)

! 16. xi=delta-beta^2*xi
        XI = DELTA - BETA**2*XI

      END DO

      IF (ITNO>MAXIT) THEN
        STATUS = -1
        ITNO = MAXIT
      END IF

9999  CONTINUE

      IF ((PRECONTYPE==2) .OR. (PRECONTYPE==3)) THEN
        CALL CCOPY(LOCLEN,X,1,WRK(IZ),1)
        CALL PRECONR(WRK(IZ),X,IPAR)
      END IF

!  Set output parameters
      IPAR(11) = ITNO
      IPAR(12) = STATUS
      IPAR(13) = STEPERR
      SPAR(2) = EXITNORM

      RETURN

    END SUBROUTINE PIMCCG
    SUBROUTINE PIMCCGNE(X,B,WRK,IPAR,SPAR,MATVEC,TMATVEC,PRECONL,PRECONR, &
        PCSUM,PSCNRM,PROGRESS)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

!           PIM -- The Parallel Iterative Methods package
!           ---------------------------------------------

!                      Rudnei Dias da Cunha
!     National Supercomputing Centre and Mathematics Institute
!         Universidade Federal do Rio Grande do Sul, Brasil

!                          Tim Hopkins
!     Computing Laboratory, University of Kent at Canterbury, U.K.

! ----------------------------------------------------------------------

!     .. Parameters ..
      REAL (WP) :: ONE
      PARAMETER (ONE=1.0E0_WP)
      COMPLEX (WP) :: CZERO
      PARAMETER (CZERO=(0.0E0_WP,0.0E0_WP))
      COMPLEX (WP) :: CONE
      PARAMETER (CONE=(1.0E0_WP,0.0E0_WP))
      INTEGER :: IPARSIZ
      PARAMETER (IPARSIZ=13)
      INTEGER :: SPARSIZ
!Art      PARAMETER (SPARSIZ=2)
      PARAMETER (SPARSIZ=6)
!     ..
!     .. Array Arguments ..
      COMPLEX (WP) :: B(*), WRK(*), X(*)
      REAL (WP) :: SPAR(SPARSIZ)
      INTEGER :: IPAR(IPARSIZ)
!     ..
!     .. Function Arguments ..
      REAL (WP) :: PSCNRM
      EXTERNAL PSCNRM
!     ..
!     .. Subroutine Arguments ..
      EXTERNAL MATVEC, PCSUM, PRECONL, PRECONR, PROGRESS, TMATVEC
!     ..
!     .. Local Scalars ..
      COMPLEX (WP) :: ALPHA, BETA, DELTA, RDOTR, RDOTR0, XI
      REAL (WP) :: EPSILON, EXITNORM, RHSSTOP
      INTEGER :: BASISDIM, BLKSZ, CNVRTX, IP, IR, IS, ITNO, IW, IXOLD, IZ, &
        LDA, LOCLEN, MAXIT, N, NPROCS, PRECONTYPE, PROCID, STATUS, STEPERR, &
        STOPTYPE
!     ..
!     .. Local Arrays ..
      COMPLEX (WP) :: DOTS(2)
!     ..
!     .. External Functions ..
      COMPLEX (WP) :: CDOTC
      REAL (WP) :: SCSETRHSSTOP
      EXTERNAL CDOTC, SCSETRHSSTOP
!     ..
!     .. External Subroutines ..
      EXTERNAL CAXPY, CCOPY, PIMSGETPAR, STOPCRIT
!     ..

      CALL PIMSGETPAR(IPAR,SPAR,LDA,N,BLKSZ,LOCLEN,BASISDIM,NPROCS,PROCID, &
        PRECONTYPE,STOPTYPE,MAXIT,ITNO,STATUS,STEPERR,EPSILON,EXITNORM)

!  Check consistency of preconditioning and stop types
      IF (((PRECONTYPE==0) .OR. (PRECONTYPE==2)) .AND. (STOPTYPE==6)) THEN
        ITNO = 0
        STATUS = -4
        STEPERR = 0
        GO TO 9999

      END IF

!  Needs conversion Y=Q2X for residual
      CNVRTX = 1

!  Set indices for mapping local vectors into wrk
      IR = 1
      IP = IR + LOCLEN
      IW = IP + LOCLEN
      IZ = IW + LOCLEN
      IS = IZ + LOCLEN
      IXOLD = IS + LOCLEN

!  Set rhs of stopping criteria
      RHSSTOP = SCSETRHSSTOP(B,WRK(IR),EPSILON,IPAR,PRECONL,PSCNRM)

!  1. r=Q1(b-AA^{T}Q2x)
      IF (STOPTYPE/=6) THEN
        IF (PRECONTYPE==0) THEN
!     r=b-Ax
          CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL TMATVEC(X,WRK(IP),IPAR)
          CALL MATVEC(WRK(IP),WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==1) THEN
!     r=Q1(b-AA^{T})
          CALL CCOPY(LOCLEN,B,1,WRK(IZ),1)
          CALL TMATVEC(X,WRK(IP),IPAR)
          CALL MATVEC(WRK(IP),WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IZ),1)
          CALL PRECONL(WRK(IZ),WRK(IR),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
!     r=b-AA^{T}Q2x
          CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL PRECONR(X,WRK(IP),IPAR)
          CALL TMATVEC(WRK(IP),WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(b-AA^{T}Q2x)
          CALL CCOPY(LOCLEN,B,1,WRK(IP),1)
          CALL PRECONR(X,WRK(IW),IPAR)
          CALL TMATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IP),1)
          CALL PRECONL(WRK(IP),WRK(IR),IPAR)
        END IF

      ELSE
!     r has been set to Qb in the call to dsetrhsstop
        IF (PRECONTYPE==1) THEN
!     r=Q1(b-AA^{T}x)
          CALL TMATVEC(X,WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IP),IPAR)
          CALL PRECONL(WRK(IP),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(b-AA^{T}Q2x)
          CALL PRECONR(X,WRK(IZ),IPAR)
          CALL TMATVEC(WRK(IZ),WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IR),1)
        END IF

      END IF

!  2. p=r
      CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IP),1)

!  3. rdotr=dot(r,r)
      DOTS(1) = CDOTC(LOCLEN,WRK(IR),1,WRK(IR),1)
      CALL PCSUM(1,DOTS)
      RDOTR = DOTS(1)

!  4. w=Q1AA^{T}Q2p
      IF (PRECONTYPE==0) THEN
        CALL TMATVEC(WRK(IP),WRK(IZ),IPAR)
        CALL MATVEC(WRK(IZ),WRK(IW),IPAR)

      ELSE IF (PRECONTYPE==1) THEN
        CALL TMATVEC(WRK(IP),WRK(IW),IPAR)
        CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
        CALL PRECONL(WRK(IZ),WRK(IW),IPAR)

      ELSE IF (PRECONTYPE==2) THEN
        CALL PRECONR(WRK(IP),WRK(IW),IPAR)
        CALL TMATVEC(WRK(IW),WRK(IZ),IPAR)
        CALL MATVEC(WRK(IZ),WRK(IW),IPAR)

      ELSE IF (PRECONTYPE==3) THEN
        CALL PRECONR(WRK(IP),WRK(IZ),IPAR)
        CALL TMATVEC(WRK(IZ),WRK(IW),IPAR)
        CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
        CALL PRECONL(WRK(IZ),WRK(IW),IPAR)
      END IF

!  5. xi=dot(p,w)
      DOTS(1) = CDOTC(LOCLEN,WRK(IP),1,WRK(IW),1)
      CALL PCSUM(1,DOTS)
      XI = DOTS(1)

!  Loop
      STATUS = 0
      EXITNORM = -ONE
      STEPERR = -1
      DO ITNO = 1, MAXIT

!  6. alpha=rdotr/xi
        IF (XI==CZERO) THEN
          STATUS = -3
          STEPERR = 6
          GO TO 9999

        END IF

        ALPHA = RDOTR/XI

!  7. x=x+alpha*p
        CALL CCOPY(LOCLEN,X,1,WRK(IXOLD),1)
        CALL CAXPY(LOCLEN,ALPHA,WRK(IP),1,X,1)

!  8. r=r-alpha*w
        CALL CAXPY(LOCLEN,-ALPHA,WRK(IW),1,WRK(IR),1)

!  9. check stopping criterion
        CALL STOPCRIT(B,WRK(IR),WRK(IZ),X,WRK(IXOLD),WRK(IS),RHSSTOP,CNVRTX, &
          EXITNORM,STATUS,IPAR,MATVEC,TMATVEC,PRECONR,PCSUM,PSCNRM)

!  Call monitoring routine
        CALL PROGRESS(LOCLEN,ITNO,EXITNORM,X,WRK(IR),WRK(IZ))

        IF (STATUS==0) THEN
          GO TO 9999
        END IF

! 10. s=Q1AA^{T}Q2r
        IF (PRECONTYPE==0) THEN
          CALL TMATVEC(WRK(IR),WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IS),IPAR)

        ELSE IF (PRECONTYPE==1) THEN
          CALL TMATVEC(WRK(IR),WRK(IS),IPAR)
          CALL MATVEC(WRK(IS),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IS),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
          CALL PRECONR(WRK(IR),WRK(IS),IPAR)
          CALL TMATVEC(WRK(IS),WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IS),IPAR)

        ELSE IF (PRECONTYPE==3) THEN
          CALL PRECONR(WRK(IR),WRK(IZ),IPAR)
          CALL TMATVEC(WRK(IZ),WRK(IS),IPAR)
          CALL MATVEC(WRK(IS),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IS),IPAR)
        END IF

! 11. rdotr=dot(r,r)
        RDOTR0 = RDOTR
        DOTS(1) = CDOTC(LOCLEN,WRK(IR),1,WRK(IR),1)

! 12. delta=dot(r,s)
        DOTS(2) = CDOTC(LOCLEN,WRK(IR),1,WRK(IS),1)

!  Accumulate simultaneously partial values
        CALL PCSUM(2,DOTS)
        RDOTR = DOTS(1)
        DELTA = DOTS(2)

! 13. beta=rdotr/rdotr0
        IF (RDOTR0==CZERO) THEN
          STATUS = -3
          STEPERR = 10
          GO TO 9999

        END IF

        BETA = RDOTR/RDOTR0

! 14. p=r+beta*p
        CALL CCOPY(LOCLEN,WRK(IP),1,WRK(IZ),1)
        CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IP),1)
        CALL CAXPY(LOCLEN,BETA,WRK(IZ),1,WRK(IP),1)

! 15. w=s+beta*w
        CALL CCOPY(LOCLEN,WRK(IW),1,WRK(IZ),1)
        CALL CCOPY(LOCLEN,WRK(IS),1,WRK(IW),1)
        CALL CAXPY(LOCLEN,BETA,WRK(IZ),1,WRK(IW),1)

! 16. xi=delta-beta^2*xi
        XI = DELTA - BETA**2*XI

      END DO

      IF (ITNO>MAXIT) THEN
        STATUS = -1
        ITNO = MAXIT
      END IF

9999  CONTINUE

      IF ((PRECONTYPE==2) .OR. (PRECONTYPE==3)) THEN
        CALL CCOPY(LOCLEN,X,1,WRK(IZ),1)
        CALL PRECONR(WRK(IZ),WRK(IXOLD),IPAR)
        CALL TMATVEC(WRK(IXOLD),X,IPAR)
      ELSE
        CALL CCOPY(LOCLEN,X,1,WRK(IZ),1)
        CALL TMATVEC(WRK(IZ),X,IPAR)
      END IF

!  Set output parameters
      IPAR(11) = ITNO
      IPAR(12) = STATUS
      IPAR(13) = STEPERR
      SPAR(2) = EXITNORM

      RETURN

    END SUBROUTINE PIMCCGNE
    SUBROUTINE PIMCCGNR(X,B,WRK,IPAR,SPAR,MATVEC,TMATVEC,PRECONL,PRECONR, &
        PCSUM,PSCNRM,PROGRESS)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

!           PIM -- The Parallel Iterative Methods package
!           ---------------------------------------------

!                      Rudnei Dias da Cunha
!     National Supercomputing Centre and Mathematics Institute
!         Universidade Federal do Rio Grande do Sul, Brasil

!                          Tim Hopkins
!     Computing Laboratory, University of Kent at Canterbury, U.K.

! ----------------------------------------------------------------------

!     .. Parameters ..
      REAL (WP) :: ONE
      PARAMETER (ONE=1.0E0_WP)
      COMPLEX (WP) :: CZERO
      PARAMETER (CZERO=(0.0E0_WP,0.0E0_WP))
      COMPLEX (WP) :: CONE
      PARAMETER (CONE=(1.0E0_WP,0.0E0_WP))
      INTEGER :: IPARSIZ
      PARAMETER (IPARSIZ=13)
      INTEGER :: SPARSIZ
!Art      PARAMETER (SPARSIZ=2)
      PARAMETER (SPARSIZ=6)
!     ..
!     .. Array Arguments ..
      COMPLEX (WP) :: B(*), WRK(*), X(*)
      REAL (WP) :: SPAR(SPARSIZ)
      INTEGER :: IPAR(IPARSIZ)
!     ..
!     .. Function Arguments ..
      REAL (WP) :: PSCNRM
      EXTERNAL PSCNRM
!     ..
!     .. Subroutine Arguments ..
      EXTERNAL MATVEC, PCSUM, PRECONL, PRECONR, PROGRESS, TMATVEC
!     ..
!     .. Local Scalars ..
      COMPLEX (WP) :: ALPHA, BETA, DELTA, RDOTR, RDOTR0, XI
      REAL (WP) :: EPSILON, EXITNORM, RHSSTOP
      INTEGER :: BASISDIM, BLKSZ, CNVRTX, IP, IR, IS, ITNO, IW, IXOLD, IZ, &
        LDA, LOCLEN, MAXIT, N, NPROCS, PRECONTYPE, PROCID, STATUS, STEPERR, &
        STOPTYPE
!     ..
!     .. Local Arrays ..
      COMPLEX (WP) :: DOTS(2)
!     ..
!     .. External Functions ..
      COMPLEX (WP) :: CDOTC
      REAL (WP) :: SCSETRHSSTOP
      EXTERNAL CDOTC, SCSETRHSSTOP
!     ..
!     .. External Subroutines ..
      EXTERNAL CAXPY, CCOPY, PIMSGETPAR, STOPCRIT
!     ..

      CALL PIMSGETPAR(IPAR,SPAR,LDA,N,BLKSZ,LOCLEN,BASISDIM,NPROCS,PROCID, &
        PRECONTYPE,STOPTYPE,MAXIT,ITNO,STATUS,STEPERR,EPSILON,EXITNORM)

!  Check consistency of preconditioning and stop types
      IF (((PRECONTYPE==0) .OR. (PRECONTYPE==2)) .AND. (STOPTYPE==6)) THEN
        ITNO = 0
        STATUS = -4
        STEPERR = 0
        GO TO 9999

      END IF

!  Does not need conversion Y=Q2X for residual
      CNVRTX = 0

!  Set indices for mapping local vectors into wrk
      IR = 1
      IP = IR + LOCLEN
      IW = IP + LOCLEN
      IZ = IW + LOCLEN
      IS = IZ + LOCLEN
      IXOLD = IS + LOCLEN

!  Set rhs of stopping criteria
!  On (P)CGNR, the rhs vector is actually A^{T}b
      CALL TMATVEC(B,WRK(IZ),IPAR)
      RHSSTOP = SCSETRHSSTOP(WRK(IZ),WRK(IR),EPSILON,IPAR,PRECONL,PSCNRM)

!  1. r=Q1(A^{T}b-A^{T}AQ2x)
      IF (STOPTYPE/=6) THEN
        IF (PRECONTYPE==0) THEN
!     r=b-Ax
          CALL MATVEC(X,WRK(IP),IPAR)
          CALL TMATVEC(WRK(IP),WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IZ),1)
          CALL CCOPY(LOCLEN,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==1) THEN
!     r=Q1(A^{T}b-A^{T}Ax)
          CALL MATVEC(X,WRK(IP),IPAR)
          CALL TMATVEC(WRK(IP),WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IZ),1)
          CALL PRECONL(WRK(IZ),WRK(IR),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
!     r=A^{T}b-A^{T}AQ2x
          CALL PRECONR(X,WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IP),IPAR)
          CALL TMATVEC(WRK(IP),WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IZ),1)
          CALL CCOPY(LOCLEN,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(A^{T}b-A^{T}AQ2x)
          CALL PRECONR(X,WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IP),IPAR)
          CALL TMATVEC(WRK(IP),WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IZ),1)
          CALL PRECONL(WRK(IZ),WRK(IR),IPAR)
        END IF

      ELSE
!     r has been set to Q(A^{T}b) in the call to dsetrhsstop
        IF (PRECONTYPE==1) THEN
!     r=Q1(A^{T}b-A^{T}Ax)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL TMATVEC(WRK(IW),WRK(IP),IPAR)
          CALL PRECONL(WRK(IP),WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(A^{T}b-A^{T}AQ2x)
          CALL PRECONR(X,WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IP),IPAR)
          CALL TMATVEC(WRK(IP),WRK(IW),IPAR)
          CALL PRECONL(WRK(IW),WRK(IP),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IP),1,WRK(IR),1)
        END IF

      END IF

!  2. p=r
      CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IP),1)

!  3. rdotr=dot(r,r)
      DOTS(1) = CDOTC(LOCLEN,WRK(IR),1,WRK(IR),1)
      CALL PCSUM(1,DOTS)
      RDOTR = DOTS(1)

!  4. w=Q1A^{T}AQ2p
      IF (PRECONTYPE==0) THEN
        CALL MATVEC(WRK(IP),WRK(IZ),IPAR)
        CALL TMATVEC(WRK(IZ),WRK(IW),IPAR)

      ELSE IF (PRECONTYPE==1) THEN
        CALL MATVEC(WRK(IP),WRK(IW),IPAR)
        CALL TMATVEC(WRK(IW),WRK(IZ),IPAR)
        CALL PRECONL(WRK(IZ),WRK(IW),IPAR)

      ELSE IF (PRECONTYPE==2) THEN
        CALL PRECONR(WRK(IP),WRK(IW),IPAR)
        CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
        CALL TMATVEC(WRK(IZ),WRK(IW),IPAR)

      ELSE IF (PRECONTYPE==3) THEN
        CALL PRECONR(WRK(IP),WRK(IZ),IPAR)
        CALL MATVEC(WRK(IZ),WRK(IW),IPAR)
        CALL TMATVEC(WRK(IW),WRK(IZ),IPAR)
        CALL PRECONL(WRK(IZ),WRK(IW),IPAR)
      END IF

!  5. xi=dot(p,w)
      DOTS(1) = CDOTC(LOCLEN,WRK(IP),1,WRK(IW),1)
      CALL PCSUM(1,DOTS)
      XI = DOTS(1)

!  Loop
      STATUS = 0
      EXITNORM = -ONE
      STEPERR = -1
      DO ITNO = 1, MAXIT

!  6. alpha=rdotr/xi
        IF (XI==CZERO) THEN
          STATUS = -3
          STEPERR = 6
          GO TO 9999

        END IF

        ALPHA = RDOTR/XI

!  7. x=x+alpha*p
        CALL CCOPY(LOCLEN,X,1,WRK(IXOLD),1)
        CALL CAXPY(LOCLEN,ALPHA,WRK(IP),1,X,1)

!  8. r=r-alpha*w
        CALL CAXPY(LOCLEN,-ALPHA,WRK(IW),1,WRK(IR),1)

!  9. check stopping criterion
        CALL STOPCRIT(B,WRK(IR),WRK(IZ),X,WRK(IXOLD),WRK(IS),RHSSTOP,CNVRTX, &
          EXITNORM,STATUS,IPAR,MATVEC,TMATVEC,PRECONR,PCSUM,PSCNRM)

!  Call monitoring routine
        CALL PROGRESS(LOCLEN,ITNO,EXITNORM,X,WRK(IR),WRK(IZ))

        IF (STATUS==0) THEN
          GO TO 9999
        END IF

! 10. s=Q1A^{T}AQ2r
        IF (PRECONTYPE==0) THEN
          CALL MATVEC(WRK(IR),WRK(IZ),IPAR)
          CALL TMATVEC(WRK(IZ),WRK(IS),IPAR)

        ELSE IF (PRECONTYPE==1) THEN
          CALL MATVEC(WRK(IR),WRK(IS),IPAR)
          CALL TMATVEC(WRK(IS),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IS),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
          CALL PRECONR(WRK(IR),WRK(IS),IPAR)
          CALL MATVEC(WRK(IS),WRK(IZ),IPAR)
          CALL TMATVEC(WRK(IZ),WRK(IS),IPAR)

        ELSE IF (PRECONTYPE==3) THEN
          CALL PRECONR(WRK(IR),WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IS),IPAR)
          CALL TMATVEC(WRK(IS),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IS),IPAR)
        END IF

! 11. rdotr=dot(r,r)
        RDOTR0 = RDOTR
        DOTS(1) = CDOTC(LOCLEN,WRK(IR),1,WRK(IR),1)

! 12. delta=dot(r,s)
        DOTS(2) = CDOTC(LOCLEN,WRK(IR),1,WRK(IS),1)

!  Accumulate simultaneously partial values
        CALL PCSUM(2,DOTS)
        RDOTR = DOTS(1)
        DELTA = DOTS(2)

! 13. beta=rdotr/rdotr0
        IF (RDOTR0==CZERO) THEN
          STATUS = -3
          STEPERR = 13
          GO TO 9999

        END IF

        BETA = RDOTR/RDOTR0

! 14. p=r+beta*p
        CALL CCOPY(LOCLEN,WRK(IP),1,WRK(IZ),1)
        CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IP),1)
        CALL CAXPY(LOCLEN,BETA,WRK(IZ),1,WRK(IP),1)

! 15. w=s+beta*w
        CALL CCOPY(LOCLEN,WRK(IW),1,WRK(IZ),1)
        CALL CCOPY(LOCLEN,WRK(IS),1,WRK(IW),1)
        CALL CAXPY(LOCLEN,BETA,WRK(IZ),1,WRK(IW),1)

! 16. xi=delta-beta^2*xi
        XI = DELTA - BETA**2*XI

      END DO

      IF (ITNO>MAXIT) THEN
        STATUS = -1
        ITNO = MAXIT
      END IF

9999  CONTINUE

      IF ((PRECONTYPE==2) .OR. (PRECONTYPE==3)) THEN
        CALL CCOPY(LOCLEN,X,1,WRK(IZ),1)
        CALL PRECONR(WRK(IZ),X,IPAR)
      END IF

!  Set output parameters
      IPAR(11) = ITNO
      IPAR(12) = STATUS
      IPAR(13) = STEPERR
      SPAR(2) = EXITNORM

      RETURN

    END SUBROUTINE PIMCCGNR
    SUBROUTINE PIMCCGS(X,B,WRK,IPAR,SPAR,MATVEC,PRECONL,PRECONR,PCSUM,PSCNRM, &
        PROGRESS)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

!           PIM -- The Parallel Iterative Methods package
!           ---------------------------------------------

!                      Rudnei Dias da Cunha
!     National Supercomputing Centre and Mathematics Institute
!         Universidade Federal do Rio Grande do Sul, Brasil

!                          Tim Hopkins
!     Computing Laboratory, University of Kent at Canterbury, U.K.

! ----------------------------------------------------------------------

!     .. Parameters ..
      REAL (WP) :: ONE
      PARAMETER (ONE=1.0E0_WP)
      COMPLEX (WP) :: CZERO
      PARAMETER (CZERO=(0.0E0_WP,0.0E0_WP))
      COMPLEX (WP) :: CONE
      PARAMETER (CONE=(1.0E0_WP,0.0E0_WP))
      INTEGER :: IPARSIZ
      PARAMETER (IPARSIZ=13)
      INTEGER :: SPARSIZ
!Art      PARAMETER (SPARSIZ=2)
      PARAMETER (SPARSIZ=6)
!     ..
!     .. Array Arguments ..
      COMPLEX (WP) :: B(*), WRK(*), X(*)
      REAL (WP) :: SPAR(SPARSIZ)
      INTEGER :: IPAR(IPARSIZ)
!     ..
!     .. Function Arguments ..
      REAL (WP) :: PSCNRM
      EXTERNAL PSCNRM
!     ..
!     .. Subroutine Arguments ..
      EXTERNAL MATVEC, PCSUM, PRECONL, PRECONR, PROGRESS
!     ..
!     .. Local Scalars ..
      COMPLEX (WP) :: ALPHA, BETA, RHO, RHO0, XI
      REAL (WP) :: EPSILON, EXITNORM, RHSSTOP
      INTEGER :: BASISDIM, BLKSZ, CNVRTX, IP, IQ, IR, IRTILDE, IS, IT, ITNO, &
        IU, IW, IXOLD, IZ, LDA, LOCLEN, MAXIT, N, NPROCS, PRECONTYPE, PROCID, &
        STATUS, STEPERR, STOPTYPE
!     ..
!     .. Local Arrays ..
      COMPLEX (WP) :: DOTS(1)
!     ..
!     .. External Functions ..
      COMPLEX (WP) :: CDOTC
      REAL (WP) :: SCSETRHSSTOP
      EXTERNAL CDOTC, SCSETRHSSTOP
!     ..
!     .. External Subroutines ..
      EXTERNAL CAXPY, CCOPY, PIMSGETPAR, STOPCRIT
!     ..
      CALL PIMSGETPAR(IPAR,SPAR,LDA,N,BLKSZ,LOCLEN,BASISDIM,NPROCS,PROCID, &
        PRECONTYPE,STOPTYPE,MAXIT,ITNO,STATUS,STEPERR,EPSILON,EXITNORM)

!  Check consistency of preconditioning and stop types
      IF (((PRECONTYPE==0) .OR. (PRECONTYPE==2)) .AND. (STOPTYPE==6)) THEN
        ITNO = 0
        STATUS = -4
        STEPERR = 0
        GO TO 9999

      END IF

!  Does not need conversion Y=Q2X for residual
      CNVRTX = 0

!  Set indices for mapping local vectors into wrk
      IR = 1
      IRTILDE = IR + LOCLEN
      IP = IRTILDE + LOCLEN
      IQ = IP + LOCLEN
      IS = IQ + LOCLEN
      IT = IS + LOCLEN
      IU = IT + LOCLEN
      IW = IU + LOCLEN
      IZ = IW + LOCLEN
      IXOLD = IZ + LOCLEN

!  Set rhs of stopping criteria
      RHSSTOP = SCSETRHSSTOP(B,WRK(IR),EPSILON,IPAR,PRECONL,PSCNRM)

!  1. r=Q1(b-AQ2x)
      IF (STOPTYPE/=6) THEN
        IF (PRECONTYPE==0) THEN
!     r=b-Ax
          CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==1) THEN
!     r=Q1(b-Ax)
          CALL CCOPY(LOCLEN,B,1,WRK(IZ),1)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IZ),1)
          CALL PRECONL(WRK(IZ),WRK(IR),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
!     r=b-AQ2x
          CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL PRECONR(X,WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(b-AQ2x)
          CALL CCOPY(LOCLEN,B,1,WRK(IP),1)
          CALL PRECONR(X,WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IP),1)
          CALL PRECONL(WRK(IP),WRK(IR),IPAR)
        END IF

      ELSE
!     r has been set to Qb in the call to ssetrhsstop
        IF (PRECONTYPE==1) THEN
!     r=Q1(b-Ax)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL PRECONL(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(b-AQ2x)
          CALL PRECONR(X,WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IW),IPAR)
          CALL PRECONL(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)
        END IF

      END IF

!  2. p=s=rtilde=r
      CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IRTILDE),1)
      CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IP),1)
      CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IS),1)

!  3. rho=dot(rtilde,r)
      DOTS(1) = CDOTC(LOCLEN,WRK(IRTILDE),1,WRK(IR),1)
      CALL PCSUM(1,DOTS)
      RHO = DOTS(1)

!  Loop
      STATUS = 0
      EXITNORM = -ONE
      STEPERR = -1
      DO ITNO = 1, MAXIT

!  4. w=Q1AQ2p
        IF (PRECONTYPE==0) THEN
          CALL MATVEC(WRK(IP),WRK(IW),IPAR)

        ELSE IF (PRECONTYPE==1) THEN
          CALL MATVEC(WRK(IP),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IW),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
          CALL PRECONR(WRK(IP),WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IW),IPAR)

        ELSE IF (PRECONTYPE==3) THEN
          CALL PRECONR(WRK(IP),WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IW),IPAR)
        END IF

!  5. xi=dot(rtilde,w)
        DOTS(1) = CDOTC(LOCLEN,WRK(IRTILDE),1,WRK(IW),1)
        CALL PCSUM(1,DOTS)
        XI = DOTS(1)

!  6. alpha=rho/xi
        IF (XI==CZERO) THEN
          STATUS = -3
          STEPERR = 6
          GO TO 9999

        END IF

        ALPHA = RHO/XI

!  7. t=s-alpha*w
        CALL CCOPY(LOCLEN,WRK(IS),1,WRK(IT),1)
        CALL CAXPY(LOCLEN,-ALPHA,WRK(IW),1,WRK(IT),1)

!  8. w=s+t
        CALL CCOPY(LOCLEN,WRK(IS),1,WRK(IW),1)
        CALL CAXPY(LOCLEN,CONE,WRK(IT),1,WRK(IW),1)

!  9. x=x+alpha*w
        CALL CCOPY(LOCLEN,X,1,WRK(IXOLD),1)
        CALL CAXPY(LOCLEN,ALPHA,WRK(IW),1,X,1)

! 10. r=r-alpha*Q1AQ2w
        IF (PRECONTYPE==0) THEN
          CALL MATVEC(WRK(IW),WRK(IU),IPAR)

        ELSE IF (PRECONTYPE==1) THEN
          CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IU),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
          CALL PRECONR(WRK(IW),WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IU),IPAR)

        ELSE IF (PRECONTYPE==3) THEN
          CALL PRECONR(WRK(IW),WRK(IU),IPAR)
          CALL MATVEC(WRK(IU),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IU),IPAR)
        END IF

        CALL CAXPY(LOCLEN,-ALPHA,WRK(IU),1,WRK(IR),1)

! 11. check stopping criterion
        CALL STOPCRIT(B,WRK(IR),WRK(IZ),X,WRK(IXOLD),WRK(IU),RHSSTOP,CNVRTX, &
          EXITNORM,STATUS,IPAR,MATVEC,MATVEC,PRECONR,PCSUM,PSCNRM)

!  Call monitoring routine
        CALL PROGRESS(LOCLEN,ITNO,EXITNORM,X,WRK(IR),WRK(IZ))

        IF (STATUS==0) THEN
          GO TO 9999
        END IF
! 12. rho=dot(rtilde0,r)
        RHO0 = RHO
        DOTS(1) = CDOTC(LOCLEN,WRK(IRTILDE),1,WRK(IR),1)
        CALL PCSUM(1,DOTS)
        RHO = DOTS(1)

! 13. beta=rho/rho0
        IF (RHO0==CZERO) THEN
          STATUS = -3
          STEPERR = 13
          GO TO 9999

        END IF

        BETA = RHO/RHO0

! 14. s=r+beta*t
        CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IS),1)
        CALL CAXPY(LOCLEN,BETA,WRK(IT),1,WRK(IS),1)

! 15. w=t+beta*p
        CALL CCOPY(LOCLEN,WRK(IT),1,WRK(IW),1)
        CALL CAXPY(LOCLEN,BETA,WRK(IP),1,WRK(IW),1)

! 16. p=s+beta*w
        CALL CCOPY(LOCLEN,WRK(IS),1,WRK(IP),1)
        CALL CAXPY(LOCLEN,BETA,WRK(IW),1,WRK(IP),1)

      END DO

      IF (ITNO>MAXIT) THEN
        STATUS = -1
        ITNO = MAXIT
      END IF

9999  CONTINUE

      IF ((PRECONTYPE==2) .OR. (PRECONTYPE==3)) THEN
        CALL CCOPY(LOCLEN,X,1,WRK(IZ),1)
        CALL PRECONR(WRK(IZ),X,IPAR)
      END IF

!  Set output parameters
      IPAR(11) = ITNO
      IPAR(12) = STATUS
      IPAR(13) = STEPERR
      SPAR(2) = EXITNORM

      RETURN

    END SUBROUTINE PIMCCGS
    SUBROUTINE PIMCCHEBYSHEV(X,B,WRK,IPAR,SPAR,MATVEC,PRECONL,PRECONR,PCSUM, &
        PSCNRM,PROGRESS)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

!           PIM -- The Parallel Iterative Methods package
!           ---------------------------------------------

!                      Rudnei Dias da Cunha
!     National Supercomputing Centre and Mathematics Institute
!         Universidade Federal do Rio Grande do Sul, Brasil

!                          Tim Hopkins
!     Computing Laboratory, University of Kent at Canterbury, U.K.

! ----------------------------------------------------------------------

!     .. Parameters ..
      REAL (WP) :: ZERO
      PARAMETER (ZERO=0.0E0_WP)
      COMPLEX (WP) :: ZONE
      PARAMETER (ZONE=(1.0E0_WP,0.0E0_WP))
      REAL (WP) :: ONE
      PARAMETER (ONE=1.0E0_WP)
      REAL (WP) :: TWO
      PARAMETER (TWO=2.0_WP)
      INTEGER :: IPARSIZ
      PARAMETER (IPARSIZ=13)
      INTEGER :: SPARSIZ
      PARAMETER (SPARSIZ=6)
!     ..
!     .. Array Arguments ..
      COMPLEX (WP) :: B(*), WRK(*), X(*)
      REAL (WP) :: SPAR(SPARSIZ)
      INTEGER :: IPAR(IPARSIZ)
!     ..
!     .. Function Arguments ..
      REAL (WP) :: PSCNRM
      EXTERNAL PSCNRM
!     ..
!     .. Subroutine Arguments ..
      EXTERNAL MATVEC, PCSUM, PRECONL, PRECONR, PROGRESS
!     ..
!     .. Local Scalars ..
      REAL (WP) :: AXISISQ, AXISRSQ, D, DELTA, EPSILON, EXITNORM, GAMMA, &
        LENGTHI, LENGTHR, RHO, RHSSTOP, SIGMA, SIGMASQ
      INTEGER :: BASISDIM, BLKSZ, CNVRTX, IK, IR, ITNO, IW, IXOLD, IZ, LDA, &
        LOCLEN, MAXIT, N, NPROCS, PRECONTYPE, PROCID, STATUS, STEPERR, &
        STOPTYPE
!     ..
!     .. External Functions ..
      REAL (WP) :: SCSETRHSSTOP
      EXTERNAL SCSETRHSSTOP
!     ..
!     .. External Subroutines ..
      EXTERNAL CAXPY, CCOPY, CSCAL, CSWAP, PIMSGETPAR, STOPCRIT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX, MAX
!     ..
      CALL PIMSGETPAR(IPAR,SPAR,LDA,N,BLKSZ,LOCLEN,BASISDIM,NPROCS,PROCID, &
        PRECONTYPE,STOPTYPE,MAXIT,ITNO,STATUS,STEPERR,EPSILON,EXITNORM)

!  Check consistency of stop types
      IF ((STOPTYPE/=1) .AND. (STOPTYPE/=2) .AND. (STOPTYPE/=7)) THEN
        ITNO = 0
        STATUS = -6
        STEPERR = 0
        GO TO 9999

      END IF

!  Does not need conversion Y=Q2X for residual
      CNVRTX = 0

!  Set indices for mapping local vectors into wrk
      IW = 1
      IK = IW + LOCLEN
      IZ = IK + LOCLEN
      IR = IZ + LOCLEN
      IXOLD = IR + LOCLEN

!  Set rhs of stopping criteria
      RHSSTOP = SCSETRHSSTOP(B,WRK(IR),EPSILON,IPAR,PRECONL,PSCNRM)

!  1. Set parameters for iteration
      IF ((SPAR(3)==ZERO) .AND. (SPAR(4)==ZERO) .AND. (SPAR( &
          5)==ZERO) .AND. (SPAR(6)==ZERO)) THEN
        STATUS = -7
        STEPERR = 1
        GO TO 9999

      ELSE IF (SPAR(5)==SPAR(6)) THEN
!     Eigenvalues are contained in the interval [SPAR(3),SPAR(4)] on
!     the real axis:
!         sigma=(dpar(4)-dpar(3))/(2-dpar(4)-dpar(3))
!         gamma=2/(2-dpar(4)-dpar(3))
        SIGMA = (SPAR(4)-SPAR(3))/(TWO-SPAR(4)-SPAR(3))
        SIGMASQ = SIGMA*SIGMA
        GAMMA = TWO/(TWO-SPAR(4)-SPAR(3))

      ELSE IF (SPAR(3)==SPAR(4)) THEN
!     Eigenvalues are contained in the interval [SPAR(5),SPAR(6)] on
!     the imaginary axis:
!         sigma^2=-max(dpar(5),dpar(6))
!         gamma=1
        SIGMASQ = -MAX(SPAR(5),SPAR(6))
        GAMMA = ONE

      ELSE
!     Eigenvalues are complex and contained in the box
!     SPAR(3)<= Real(e) <= SPAR(4) and SPAR(5)<= Imag(e) <= SPAR(6).
!     Compute the minimum bounding ellipse that circumscribes the box;
!     this is defined by its axes a=sqrt(2)*(dpar(4)-dpar(3))/2 (along
!     the real axis) and b=sqrt(2)*(dpar(6)-dpar(5))/2 (along the
!     imaginary axis). The center of the ellipse is d.
!         sigma^2=(a^2+b^2)/(1-d)^2
!         gamma=1/(1-d)
        LENGTHR = (SPAR(4)-SPAR(3))/TWO
        LENGTHI = (SPAR(6)-SPAR(5))/TWO
        AXISRSQ = LENGTHR*LENGTHR*TWO
        AXISISQ = LENGTHI*LENGTHI*TWO
        D = (SPAR(6)+SPAR(5))/TWO
        SIGMASQ = (AXISRSQ-AXISISQ)/(ONE-D)**2
        GAMMA = ONE/(ONE-D)

      END IF

!  2. k=gamma*Q1b
      IF (PRECONTYPE==0) THEN
        CALL CCOPY(LOCLEN,B,1,WRK(IK),1)
        CALL CSCAL(LOCLEN,CMPLX(GAMMA,KIND=WP),WRK(IK),1)

      ELSE IF ((PRECONTYPE==1) .OR. (PRECONTYPE==3)) THEN
        CALL PRECONL(B,WRK(IK),IPAR)
        CALL CSCAL(LOCLEN,CMPLX(GAMMA,KIND=WP),WRK(IK),1)
      END IF

!    xold=x
      CALL CCOPY(LOCLEN,X,1,WRK(IXOLD),1)

!  Loop
      STATUS = 0
      EXITNORM = -ONE
      STEPERR = -1
      DO ITNO = 1, MAXIT

!  3. rho
        IF (ITNO==1) THEN
          RHO = ONE

        ELSE IF (ITNO==2) THEN
          RHO = ONE/(ONE-SIGMASQ/TWO)

        ELSE
          RHO = ONE/(ONE-RHO*SIGMASQ/4.0_WP)
        END IF

!  4. w=(I-Q1AQ2)x
        IF (PRECONTYPE==0) THEN
          CALL MATVEC(X,WRK(IZ),IPAR)

        ELSE IF (PRECONTYPE==1) THEN
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL PRECONL(WRK(IW),WRK(IZ),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
          CALL PRECONR(X,WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IZ),IPAR)

        ELSE IF (PRECONTYPE==3) THEN
          CALL PRECONR(X,WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IW),IPAR)
          CALL PRECONL(WRK(IW),WRK(IZ),IPAR)
        END IF

        CALL CCOPY(LOCLEN,X,1,WRK(IW),1)
        CALL CAXPY(LOCLEN,-ZONE,WRK(IZ),1,WRK(IW),1)

!  5. x=rho*(gamma*((I-Q1A)x+Q1b)+(1-gamma)*x)+(1-rho)*xold
        DELTA = RHO*GAMMA
        CALL CSCAL(LOCLEN,CMPLX(ONE-RHO,KIND=WP),WRK(IXOLD),1)
        CALL CAXPY(LOCLEN,CMPLX(RHO,KIND=WP),WRK(IK),1,WRK(IXOLD),1)
        CALL CAXPY(LOCLEN,CMPLX(RHO-DELTA,KIND=WP),X,1,WRK(IXOLD),1)
        CALL CAXPY(LOCLEN,CMPLX(DELTA,KIND=WP),WRK(IW),1,WRK(IXOLD),1)
        CALL CSWAP(LOCLEN,WRK(IXOLD),1,X,1)

!  6. check stopping criterion
        CALL STOPCRIT(B,WRK(IZ),WRK(IR),X,WRK(IXOLD),WRK(IW),RHSSTOP,CNVRTX, &
          EXITNORM,STATUS,IPAR,MATVEC,MATVEC,PRECONR,PCSUM,PSCNRM)

!  Call monitoring routine
        CALL PROGRESS(LOCLEN,ITNO,EXITNORM,X,WRK(IR),WRK(IR))

        IF (STATUS==0) THEN
          GO TO 9999
        END IF

      END DO

      IF (ITNO>MAXIT) THEN
        STATUS = -1
        ITNO = MAXIT
      END IF

9999  CONTINUE

      IF ((PRECONTYPE==2) .OR. (PRECONTYPE==3)) THEN
        CALL CCOPY(LOCLEN,X,1,WRK(IZ),1)
        CALL PRECONR(WRK(IZ),X,IPAR)
      END IF

!  Set output parameters
      IPAR(11) = ITNO
      IPAR(12) = STATUS
      IPAR(13) = STEPERR
      SPAR(2) = EXITNORM

      RETURN

    END SUBROUTINE PIMCCHEBYSHEV
    SUBROUTINE PIMCQMR(X,B,WRK,IPAR,SPAR,MATVEC,TMATVEC,PRECONL,PRECONR,PCSUM, &
        PSCNRM,PROGRESS)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

!           PIM -- The Parallel Iterative Methods package
!           ---------------------------------------------

!                      Rudnei Dias da Cunha
!     National Supercomputing Centre and Mathematics Institute
!         Universidade Federal do Rio Grande do Sul, Brasil

!                          Tim Hopkins
!     Computing Laboratory, University of Kent at Canterbury, U.K.

! ----------------------------------------------------------------------

!     .. Parameters ..
      REAL (WP) :: ONE
      PARAMETER (ONE=1.0E0_WP)
      COMPLEX (WP) :: CZERO
      PARAMETER (CZERO=(0.0E0_WP,0.0E0_WP))
      COMPLEX (WP) :: CONE
      PARAMETER (CONE=(1.0E0_WP,0.0E0_WP))
      INTEGER :: IPARSIZ
      PARAMETER (IPARSIZ=13)
      INTEGER :: SPARSIZ
!Art      PARAMETER (SPARSIZ=2)
      PARAMETER (SPARSIZ=6)
!     ..
!     .. Array Arguments ..
      COMPLEX (WP) :: B(*), WRK(*), X(*)
      REAL (WP) :: SPAR(SPARSIZ)
      INTEGER :: IPAR(IPARSIZ)
!     ..
!     .. Function Arguments ..
      REAL (WP) :: PSCNRM
      EXTERNAL PSCNRM
!     ..
!     .. Subroutine Arguments ..
      EXTERNAL MATVEC, PCSUM, PRECONL, PRECONR, PROGRESS, TMATVEC
!     ..
!     .. Local Scalars ..
      COMPLEX (WP) :: BETA, C, C0, DE, DELTA, DI, EPS, ETA, KSI, OMEGA, RHO, &
        RHO0, THETA, THETA0
      REAL (WP) :: EPSILON, EXITNORM, OVERFLOW, RHSSTOP
      INTEGER :: BASISDIM, BLKSZ, CNVRTX, ID, IP, IQ, IR, ITNO, IV, IVTILDE, &
        IW, IWTILDE, IXOLD, IY, IZ, LDA, LOCLEN, MAXIT, N, NPROCS, PRECONTYPE, &
        PROCID, STATUS, STEPERR, STOPTYPE
!     ..
!     .. Local Arrays ..
      COMPLEX (WP) :: DOTS(2)
!     ..
!     .. External Functions ..
      COMPLEX (WP) :: CDOTC
      REAL (WP) :: SCSETRHSSTOP
      EXTERNAL CDOTC, SCSETRHSSTOP
!     ..
!     .. External Subroutines ..
      EXTERNAL CAXPY, CCOPY, CINIT, CSCAL, PIMSGETPAR, SMACHCONS, STOPCRIT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS, SQRT
!     ..
      CALL SMACHCONS('O',OVERFLOW)

      CALL PIMSGETPAR(IPAR,SPAR,LDA,N,BLKSZ,LOCLEN,BASISDIM,NPROCS,PROCID, &
        PRECONTYPE,STOPTYPE,MAXIT,ITNO,STATUS,STEPERR,EPSILON,EXITNORM)

!  Check consistency of preconditioning and stop types
      IF (((PRECONTYPE==0) .OR. (PRECONTYPE==2)) .AND. (STOPTYPE==6)) THEN
        ITNO = 0
        STATUS = -4
        STEPERR = 0
        GO TO 9999

      END IF

!  Does not need conversion Y=Q2X for residual
      CNVRTX = 1

!  Set indices for mapping local vectors into wrk
      IR = 1
      IV = IR + LOCLEN
      IW = IV + LOCLEN
      IP = IW + LOCLEN
      IQ = IP + LOCLEN
      ID = IQ + LOCLEN
      IVTILDE = ID + LOCLEN
      IWTILDE = IVTILDE + LOCLEN
      IXOLD = IWTILDE + LOCLEN
      IZ = IXOLD + LOCLEN
      IY = IZ + LOCLEN

!  Set RHS of stopping criteria
      RHSSTOP = SCSETRHSSTOP(B,WRK(IR),EPSILON,IPAR,PRECONL,PSCNRM)

!  1. r=Q1(b-AQ2x)
      IF (STOPTYPE/=6) THEN
        IF (PRECONTYPE==0) THEN
!     r=b-Ax
          CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==1) THEN
!     r=Q1(b-Ax)
          CALL CCOPY(LOCLEN,B,1,WRK(IZ),1)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IZ),1)
          CALL PRECONL(WRK(IZ),WRK(IR),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
!     r=b-AQ2x
          CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL PRECONR(X,WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(b-AQ2x)
          CALL CCOPY(LOCLEN,B,1,WRK(IP),1)
          CALL PRECONR(X,WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IP),1)
          CALL PRECONL(WRK(IP),WRK(IR),IPAR)
        END IF

      ELSE
!     r has been set to Q1b in the call to dsetrhsstop
        IF (PRECONTYPE==1) THEN
!     r=Q1(b-Ax)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL PRECONL(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(b-AQ2x)
          CALL PRECONR(X,WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IW),IPAR)
          CALL PRECONL(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)
        END IF

      END IF

!  2. rho=||r||_{2}
      RHO = PSCNRM(LOCLEN,WRK(IR))

!  3. v=r/rho
      IF (RHO==CZERO) THEN
        STATUS = -3
        STEPERR = 3
        GO TO 9999

      END IF

      CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IV),1)
      CALL CSCAL(LOCLEN,CONE/RHO,WRK(IV),1)

!  4. w=-r/rho
      CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IW),1)
      CALL CSCAL(LOCLEN,-CONE/RHO,WRK(IW),1)

!  5. p=q=d=0
      CALL CINIT(LOCLEN,CZERO,WRK(IP),1)
      CALL CINIT(LOCLEN,CZERO,WRK(IQ),1)
      CALL CINIT(LOCLEN,CZERO,WRK(ID),1)

!  6. c0=1, eps=1, ksi=1, theta0=0, eta=-1, omega=1
      C0 = CONE
      EPS = CONE
      KSI = CONE
      THETA0 = CZERO
      ETA = -CONE
      OMEGA = CONE

!  Loop
      STATUS = 0
      EXITNORM = -ONE
      STEPERR = -1
      DO ITNO = 1, MAXIT

!  7. delta=w^{T}v
        DOTS(1) = CDOTC(LOCLEN,WRK(IW),1,WRK(IV),1)
        CALL PCSUM(1,DOTS)
        DELTA = DOTS(1)

!  8. if eps=0 then breakdown
        IF (EPS==CZERO) THEN
          STATUS = -3
          STEPERR = 8
          GO TO 9999

        END IF

!  9. if delta=0 then breakdown
        IF (DELTA==CZERO) THEN
          STATUS = -3
          STEPERR = 9
          GO TO 9999

        END IF

! 10. p=v-(ksi*delta/eps)*p
        DE = DELTA/EPS
        CALL CCOPY(LOCLEN,WRK(IP),1,WRK(IZ),1)
        CALL CCOPY(LOCLEN,WRK(IV),1,WRK(IP),1)
        CALL CAXPY(LOCLEN,-KSI*DE,WRK(IZ),1,WRK(IP),1)

! 11. q=w-(rho*delta/eps)*q
        CALL CCOPY(LOCLEN,WRK(IQ),1,WRK(IZ),1)
        CALL CCOPY(LOCLEN,WRK(IW),1,WRK(IQ),1)
        CALL CAXPY(LOCLEN,-RHO*DE,WRK(IZ),1,WRK(IQ),1)

! 12. vtilde=Q1AQ2p
        IF (PRECONTYPE==0) THEN
          CALL MATVEC(WRK(IP),WRK(IVTILDE),IPAR)

        ELSE IF (PRECONTYPE==1) THEN
          CALL MATVEC(WRK(IP),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IVTILDE),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
          CALL PRECONR(WRK(IP),WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IVTILDE),IPAR)

        ELSE IF (PRECONTYPE==3) THEN
          CALL PRECONR(WRK(IP),WRK(IVTILDE),IPAR)
          CALL MATVEC(WRK(IVTILDE),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IVTILDE),IPAR)
        END IF

! 13. eps=q^{T}(Q1AQ2p)=q^{T}vtilde
        DOTS(1) = CDOTC(LOCLEN,WRK(IQ),1,WRK(IVTILDE),1)
        CALL PCSUM(1,DOTS)
        EPS = DOTS(1)

! 14. beta=eps/delta
        BETA = EPS/DELTA

! 15. vtilde=vtilde-beta*v
        CALL CAXPY(LOCLEN,-BETA,WRK(IV),1,WRK(IVTILDE),1)

! 16. wtilde=Q1A^{T}Q2q-beta*w
        IF (PRECONTYPE==0) THEN
          CALL TMATVEC(WRK(IQ),WRK(IWTILDE),IPAR)

        ELSE IF (PRECONTYPE==1) THEN
          CALL TMATVEC(WRK(IQ),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IWTILDE),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
          CALL PRECONR(WRK(IQ),WRK(IZ),IPAR)
          CALL TMATVEC(WRK(IZ),WRK(IWTILDE),IPAR)

        ELSE IF (PRECONTYPE==3) THEN
          CALL PRECONR(WRK(IQ),WRK(IWTILDE),IPAR)
          CALL TMATVEC(WRK(IWTILDE),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IWTILDE),IPAR)
        END IF

        CALL CAXPY(LOCLEN,-BETA,WRK(IW),1,WRK(IWTILDE),1)

! 17. rho=norm(vtilde)
        RHO0 = RHO
        DOTS(1) = CDOTC(LOCLEN,WRK(IVTILDE),1,WRK(IVTILDE),1)

! 18. ksi=norm(wtilde)
        DOTS(2) = CDOTC(LOCLEN,WRK(IWTILDE),1,WRK(IWTILDE),1)

!  Accumulate simultaneously partial values
        CALL PCSUM(2,DOTS)
        IF (ABS(DOTS(1))>=OVERFLOW) THEN
          STATUS = -3
          STEPERR = 17
          GO TO 9999

        END IF
        RHO = SQRT(DOTS(1))

        IF (ABS(DOTS(2))>=OVERFLOW) THEN
          STATUS = -3
          STEPERR = 18
          GO TO 9999

        END IF
        KSI = SQRT(DOTS(2))

! 19. theta=(omega*rho)/(omega*c0*abs(beta))
        DI = OMEGA*C0*ABS(BETA)
        IF (DI==CZERO) THEN
          STATUS = -3
          STEPERR = 19
          GO TO 9999

        END IF
        THETA = (OMEGA*RHO)/DI

! 20. c=1/sqrt(1+theta^2)
        C = CONE/SQRT(CONE+THETA**2)

! 21. eta=-eta*rho0*c^2/(beta*c0^2)
        DI = BETA*C0**2
        IF (DI==CZERO) THEN
          STATUS = -3
          STEPERR = 21
          GO TO 9999

        END IF
        ETA = -ETA*RHO0*C**2/DI

! 22. d=p*eta+(theta0*c)^2*d
        CALL CSCAL(LOCLEN,(THETA0*C)**2,WRK(ID),1)
        CALL CAXPY(LOCLEN,ETA,WRK(IP),1,WRK(ID),1)

! 23. x=x+d
        CALL CCOPY(LOCLEN,X,1,WRK(IXOLD),1)
        CALL CAXPY(LOCLEN,CONE,WRK(ID),1,X,1)

! 24. r=b-Ax
        IF (PRECONTYPE==0) THEN
!     r=b-Ax
          CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL MATVEC(X,WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==1) THEN
!     r=Q1(b-Ax)
          CALL CCOPY(LOCLEN,B,1,WRK(IZ),1)
          CALL MATVEC(X,WRK(IY),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IY),1,WRK(IZ),1)
          CALL PRECONL(WRK(IZ),WRK(IR),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
!     r=b-AQ2x
          CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL PRECONR(X,WRK(IY),IPAR)
          CALL MATVEC(WRK(IY),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(b-AQ2x)
          CALL CCOPY(LOCLEN,B,1,WRK(IY),1)
          CALL PRECONR(X,WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IR),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IR),1,WRK(IY),1)
          CALL PRECONL(WRK(IY),WRK(IR),IPAR)
        END IF

! 25. check stopping criterion
        CALL STOPCRIT(B,WRK(IR),WRK(IZ),X,WRK(IXOLD),WRK(IW),RHSSTOP,CNVRTX, &
          EXITNORM,STATUS,IPAR,MATVEC,TMATVEC,PRECONR,PCSUM,PSCNRM)

!  Call monitoring routine
        CALL PROGRESS(LOCLEN,ITNO,EXITNORM,X,WRK(IR),WRK(IZ))

        IF (STATUS==0) THEN
          GO TO 9999
        END IF

! 26. if rho=0 then breakdown
        IF (RHO==CZERO) THEN
          STATUS = -3
          STEPERR = 26
          GO TO 9999

        END IF

! 27. if ksi=0 then breakdown
        IF (KSI==CZERO) THEN
          STATUS = -3
          STEPERR = 27
          GO TO 9999

        END IF

! 28. v=vtilde/rho
        CALL CCOPY(LOCLEN,WRK(IVTILDE),1,WRK(IV),1)
        CALL CSCAL(LOCLEN,CONE/RHO,WRK(IV),1)

! 29. w=wtilde/rho
        CALL CCOPY(LOCLEN,WRK(IWTILDE),1,WRK(IW),1)
        CALL CSCAL(LOCLEN,CONE/RHO,WRK(IW),1)

! c0=c,theta0=theta
        C0 = C
        THETA0 = THETA

      END DO

      IF (ITNO>MAXIT) THEN
        STATUS = -1
        ITNO = MAXIT
      END IF

9999  CONTINUE

      IF ((PRECONTYPE==2) .OR. (PRECONTYPE==3)) THEN
        CALL CCOPY(LOCLEN,X,1,WRK(IZ),1)
        CALL PRECONR(WRK(IZ),X,IPAR)
      END IF

!  Set output parameters
      IPAR(11) = ITNO
      IPAR(12) = STATUS
      IPAR(13) = STEPERR
      SPAR(2) = EXITNORM

      RETURN

    END SUBROUTINE PIMCQMR

    SUBROUTINE PIMCRBICGSTAB(X,B,WRK,IPAR,SPAR,MATVEC,PRECONL,PRECONR,PCSUM, &
        PSCNRM2,PROGRESS)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

!           PIM -- The Parallel Iterative Methods package
!           ---------------------------------------------

!                      Rudnei Dias da Cunha
!     National Supercomputing Centre and Mathematics Institute
!         Universidade Federal do Rio Grande do Sul, Brasil

!                          Tim Hopkins
!     Computing Laboratory, University of Kent at Canterbury, U.K.

! ----------------------------------------------------------------------

!     .. Parameters ..
      REAL (WP) :: ZERO
      PARAMETER (ZERO=0.0E0_WP)
      REAL (WP) :: ONE
      PARAMETER (ONE=1.0E0_WP)

! 04.05.21 (BTD) corrected REAL CZERO,CONE to COMPLEX CZERO,CONE

      COMPLEX (WP) :: CZERO
      PARAMETER (CZERO=(0.0E0_WP,0.0E0_WP))
      COMPLEX (WP) :: CONE
      PARAMETER (CONE=(1.0E0_WP,0.0E0_WP))
      INTEGER :: IPARSIZ
      PARAMETER (IPARSIZ=13)
      INTEGER :: SPARSIZ
!Art      PARAMETER (SPARSIZ=2)
      PARAMETER (SPARSIZ=6)
      INTEGER :: IBDIM
      PARAMETER (IBDIM=10)
!     ..
!     .. Array Arguments ..
      COMPLEX (WP) :: B(*), WRK(*), X(*)
      REAL (WP) :: SPAR(SPARSIZ)
      INTEGER :: IPAR(IPARSIZ)
!     ..
!     .. Function Arguments ..
      REAL (WP) :: PSCNRM2
      EXTERNAL PSCNRM2
!     ..
!     .. Subroutine Arguments ..
      EXTERNAL MATVEC, PCSUM, PRECONL, PRECONR, PROGRESS
!     ..
!     .. Local Scalars ..
      REAL (WP) :: ALPHA, BETA, EPSILON, EXITNORM, KSI, OMEGA, RHO0, RHO1, &
        RHSSTOP, S
      INTEGER :: BASISDIM, BLKSZ, CNVRTX, I, I0, I1, I2, I3, I4, IR, IRTILDE, &
        ITNO, IU, IW, IXOLD, IZ, J, LDA, LOCLEN, MAXIT, N, NPROCS, PRECONTYPE, &
        PROCID, STATUS, STEPERR, STOPTYPE

!     ..
!     .. Local Arrays ..
      COMPLEX (WP) :: DOTS(IBDIM)
      REAL (WP) :: GAMMA(IBDIM), GAMMA1(IBDIM), GAMMA2(IBDIM), SIGMA(IBDIM), &
        TAU(IBDIM,IBDIM)
!     ..
!     .. External Functions ..
      COMPLEX (WP) :: CDOTC
      REAL (WP) :: SCSETRHSSTOP
      EXTERNAL CDOTC, SCSETRHSSTOP
!     ..
!     .. External Subroutines ..
      EXTERNAL CAXPY, CCOPY, CINIT, PIMSGETPAR, STOPCRIT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX
!     ..
      CALL PIMSGETPAR(IPAR,SPAR,LDA,N,BLKSZ,LOCLEN,BASISDIM,NPROCS,PROCID, &
        PRECONTYPE,STOPTYPE,MAXIT,ITNO,STATUS,STEPERR,EPSILON,EXITNORM)

!  Check consistency of preconditioning and stop types
      IF (((PRECONTYPE==0) .OR. (PRECONTYPE==2)) .AND. (STOPTYPE==6)) THEN
        ITNO = 0
        STATUS = -4
        STEPERR = 0
        GO TO 9999

      END IF

!  Does not need conversion Y=Q2X for residual
      CNVRTX = 0

!  Set indices for mapping local vectors into wrk
      IRTILDE = 1
      IW = IRTILDE + LOCLEN
      IZ = IW + LOCLEN
      IXOLD = IZ + LOCLEN
      IR = IXOLD + LOCLEN
      IU = IR + (BASISDIM+1)*LOCLEN

!  Set rhs of stopping criteria
      RHSSTOP = SCSETRHSSTOP(B,WRK(IR),EPSILON,IPAR,PRECONL,PSCNRM2)

!  1. r=Q1(b-AQ2x)
      IF (PRECONTYPE==0) THEN
!     r=b-Ax
        CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
        CALL MATVEC(X,WRK(IW),IPAR)
        CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IR),1)

      ELSE IF (PRECONTYPE==1) THEN
!     r=Q1(b-Ax)
        CALL CCOPY(LOCLEN,B,1,WRK(IZ),1)
        CALL MATVEC(X,WRK(IW),IPAR)
        CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IZ),1)
        CALL PRECONL(WRK(IZ),WRK(IR),IPAR)

      ELSE IF (PRECONTYPE==2) THEN
!     r=b-AQ2x
        CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
        CALL PRECONR(X,WRK(IW),IPAR)
        CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
        CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

      ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(b-AQ2x)
        CALL CCOPY(LOCLEN,B,1,WRK(IW),1)
        CALL PRECONR(X,WRK(IR),IPAR)
        CALL MATVEC(WRK(IR),WRK(IZ),IPAR)
        CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IW),1)
        CALL PRECONL(WRK(IW),WRK(IR),IPAR)
      END IF

!  2. rtilde=r
      CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IRTILDE),1)

!  3. u0=0
      CALL CINIT(LOCLEN,CZERO,WRK(IU),1)

!  4. rho0=1, alpha=0, omega=1
      RHO0 = ONE
      ALPHA = ZERO
      OMEGA = ONE

!  Loop
      STATUS = 0
      STEPERR = -1
      EXITNORM = -ONE
      DO ITNO = 1, MAXIT

!  5. rho0=-omega*rho0
        RHO0 = -OMEGA*RHO0

!  BiCG loop
        I1 = 0
        I2 = LOCLEN
        DO J = 0, BASISDIM - 1

!  6. rho1=r(j)^{T}rtilde
          DOTS(1) = CDOTC(LOCLEN,WRK(IR+I1),1,WRK(IRTILDE),1)
          CALL PCSUM(1,DOTS)
          RHO1 = DOTS(1)

!  7. beta=alpha*rho1/rho0
          IF (RHO0==ZERO) THEN
            STATUS = -3
            STEPERR = 7
            GO TO 9999

          END IF

          BETA = ALPHA*RHO1/RHO0

!  8. rho0=rho1
          RHO0 = RHO1

!  9. u(i)=r(i)-beta*u(i), i=0:j
          I3 = 0
          DO I = 0, J
            CALL CCOPY(LOCLEN,WRK(IU+I3),1,WRK(IZ),1)
            CALL CCOPY(LOCLEN,WRK(IR+I3),1,WRK(IU+I3),1)
            CALL CAXPY(LOCLEN,CMPLX(-BETA,KIND=WP),WRK(IZ),1,WRK(IU+I3),1)
            I3 = I3 + LOCLEN
          END DO

! 10. u(j+1)=Q1AQ2u(j)
          IF (PRECONTYPE==0) THEN
            CALL MATVEC(WRK(IU+I1),WRK(IU+I2),IPAR)

          ELSE IF (PRECONTYPE==1) THEN
            CALL MATVEC(WRK(IU+I1),WRK(IW),IPAR)
            CALL PRECONL(WRK(IW),WRK(IU+I2),IPAR)

          ELSE IF (PRECONTYPE==2) THEN
            CALL PRECONR(WRK(IU+I1),WRK(IW),IPAR)
            CALL MATVEC(WRK(IW),WRK(IU+I2),IPAR)

          ELSE IF (PRECONTYPE==3) THEN
            CALL PRECONR(WRK(IU+I1),WRK(IZ),IPAR)
            CALL MATVEC(WRK(IZ),WRK(IW),IPAR)
            CALL PRECONL(WRK(IW),WRK(IU+I2),IPAR)
          END IF

! 11. ksi=u(j+1)^{T}rtilde
          DOTS(1) = CDOTC(LOCLEN,WRK(IU+I2),1,WRK(IRTILDE),1)
          CALL PCSUM(1,DOTS)
          KSI = DOTS(1)

! 12. alpha=rho0/ksi
          IF (KSI==ZERO) THEN
            STATUS = -3
            STEPERR = 12
            GO TO 9999

          END IF

          ALPHA = RHO0/KSI

! 13. r(i)=r(i)-alpha*u(i+1), i=0:j
          I3 = 0
          I4 = LOCLEN
          DO I = 0, J
            CALL CAXPY(LOCLEN,CMPLX(-ALPHA,KIND=WP),WRK(IU+I4),1,WRK(IR+I3),1)
            I3 = I3 + LOCLEN
            I4 = I4 + LOCLEN
          END DO

! 14. r(j+1)=Q1AQ2r(j)
          IF (PRECONTYPE==0) THEN
            CALL MATVEC(WRK(IR+I1),WRK(IR+I2),IPAR)

          ELSE IF (PRECONTYPE==1) THEN
            CALL MATVEC(WRK(IR+I1),WRK(IW),IPAR)
            CALL PRECONL(WRK(IW),WRK(IR+I2),IPAR)

          ELSE IF (PRECONTYPE==2) THEN
            CALL PRECONR(WRK(IR+I1),WRK(IW),IPAR)
            CALL MATVEC(WRK(IW),WRK(IR+I2),IPAR)

          ELSE IF (PRECONTYPE==3) THEN
            CALL PRECONR(WRK(IR+I1),WRK(IZ),IPAR)
            CALL MATVEC(WRK(IZ),WRK(IW),IPAR)
            CALL PRECONL(WRK(IW),WRK(IR+I2),IPAR)
          END IF

! 15. x0=x0+alpha*u0
          CALL CCOPY(LOCLEN,X,1,WRK(IXOLD),1)
          CALL CAXPY(LOCLEN,CMPLX(ALPHA,KIND=WP),WRK(IU),1,X,1)
          I1 = I1 + LOCLEN
          I2 = I2 + LOCLEN
        END DO

! 16. check stopping criterion
        CALL STOPCRIT(B,WRK(IR),WRK(IZ),X,WRK(IXOLD),WRK(IW),RHSSTOP,CNVRTX, &
          EXITNORM,STATUS,IPAR,MATVEC,MATVEC,PRECONR,PCSUM,PSCNRM2)

!  Call monitoring routine
        CALL PROGRESS(LOCLEN,ITNO,EXITNORM,X,WRK(IR),WRK(IZ))

        IF (STATUS==0) THEN
          GO TO 9999
        END IF

!  MR loop

! 17. sigma(1)=r(1)^{T}r(1), gamma'(1)=r(0)^{T}r(1)/sigma(1)
        DOTS(1) = CDOTC(LOCLEN,WRK(IR+LOCLEN),1,WRK(IR+LOCLEN),1)
        DOTS(2) = CDOTC(LOCLEN,WRK(IR),1,WRK(IR+LOCLEN),1)
        CALL PCSUM(2,DOTS)
        SIGMA(1) = DOTS(1)

        IF (SIGMA(1)==ZERO) THEN
          STATUS = -3
          STEPERR = 17
          GO TO 9999

        END IF

        GAMMA1(1) = DOTS(2)/SIGMA(1)

        I0 = LOCLEN + LOCLEN
        DO J = 2, BASISDIM

! 18. tau(i,j)=r(j)^{T}r(i)/sigma(i), r(j)=r(j)-tau(i,j)r(i)
          I1 = LOCLEN
          DO I = 1, J - 1
            DOTS(I) = CDOTC(LOCLEN,WRK(IR+I0),1,WRK(IR+I1),1)
            I1 = I1 + LOCLEN
          END DO
          CALL PCSUM(J-1,DOTS)
          I1 = LOCLEN
          DO I = 1, J - 1
            TAU(I,J) = DOTS(I)/SIGMA(I)
            CALL CAXPY(LOCLEN,CMPLX(-TAU(I,J),KIND=WP),WRK(IR+I1),1,WRK(IR+I0) &
              ,1)
          END DO

! 19. sigma(j)=r(j)^{T}r(j), gamma'(j)=r(0)^{T}r(j)/sigma(j)
          DOTS(1) = CDOTC(LOCLEN,WRK(IR+I0),1,WRK(IR+I0),1)
          DOTS(2) = CDOTC(LOCLEN,WRK(IR),1,WRK(IR+I0),1)
          CALL PCSUM(2,DOTS)
          SIGMA(J) = DOTS(1)

          IF (SIGMA(J)==ZERO) THEN
            STATUS = -3
            STEPERR = 19
            GO TO 9999

          END IF

          GAMMA1(J) = DOTS(2)/SIGMA(J)
          I0 = I0 + LOCLEN
        END DO

! 20. gamma_{l}=omega=gamma'_{l}
!     gamma_{j}=gamma'_{j}-\sum_{i=j+1}^{l}{tau_{j,i}gamma_{i}}
        GAMMA(BASISDIM) = GAMMA1(BASISDIM)
        OMEGA = GAMMA(BASISDIM)
        DO J = BASISDIM - 1, 1, -1
          S = ZERO
          DO I = J + 1, BASISDIM
            S = S + TAU(J,I)*GAMMA(I)
          END DO
          GAMMA(J) = GAMMA1(J) - S
        END DO

! 21. gamma''=gamma_{j+1}+\sum_{i=j+1}^{l-1}{tau_{j,i}gamma_{i+1}}
        DO J = 1, BASISDIM - 1
          S = ZERO
          DO I = J + 1, BASISDIM - 1
            S = S + TAU(J,I)*GAMMA(I+1)
          END DO
          GAMMA2(J) = GAMMA(J+1) + S
        END DO

!  Update

! 22. x(0)=x(0)+gamma(1)r(0)
        CALL CAXPY(LOCLEN,CMPLX(GAMMA(1),KIND=WP),WRK(IR),1,X,1)

! 23. r(0)=r(0)-gamma'(l)r(l)
        CALL CAXPY(LOCLEN,CMPLX(-GAMMA1(BASISDIM),KIND=WP), &
          WRK(IR+BASISDIM*LOCLEN),1,WRK(IR),1)

! 24. u(0)=u(0)-gamma(l)u(l)
        CALL CAXPY(LOCLEN,CMPLX(-GAMMA(BASISDIM),KIND=WP), &
          WRK(IU+BASISDIM*LOCLEN),1,WRK(IU),1)

        I0 = LOCLEN
        DO J = 1, BASISDIM - 1

! 25. u(0)=u(0)-gamma(j)u(j), j=1:l-1
          CALL CAXPY(LOCLEN,CMPLX(-GAMMA(J),KIND=WP),WRK(IU+I0),1,WRK(IU),1)

! 26. x(0)=x(0)+gamma''(j)r(j), j=1:l-1
          CALL CAXPY(LOCLEN,CMPLX(GAMMA2(J),KIND=WP),WRK(IR+I0),1,X,1)

! 27. r(0)=r(0)-gamma'(j)r(j), j=1:l-1
          CALL CAXPY(LOCLEN,CMPLX(-GAMMA1(J),KIND=WP),WRK(IR+I0),1,WRK(IR),1)
          I0 = I0 + LOCLEN
        END DO

      END DO

      IF (ITNO>MAXIT) THEN
        STATUS = -1
        ITNO = MAXIT
      END IF

9999  CONTINUE

      IF ((PRECONTYPE==2) .OR. (PRECONTYPE==3)) THEN
        CALL CCOPY(LOCLEN,X,1,WRK(IZ),1)
        CALL PRECONR(WRK(IZ),X,IPAR)
      END IF

!  Set output parameters
      IPAR(11) = ITNO
      IPAR(12) = STATUS
      IPAR(13) = STEPERR
      SPAR(2) = EXITNORM

      RETURN

    END SUBROUTINE PIMCRBICGSTAB
    SUBROUTINE PIMCRGCR(X,B,WRK,IPAR,SPAR,MATVEC,PRECONL,PRECONR,PCSUM,PSCNRM, &
        PROGRESS)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

!           PIM -- The Parallel Iterative Methods package
!           ---------------------------------------------

!                      Rudnei Dias da Cunha
!     National Supercomputing Centre and Mathematics Institute
!         Universidade Federal do Rio Grande do Sul, Brasil

!                          Tim Hopkins
!     Computing Laboratory, University of Kent at Canterbury, U.K.

! ----------------------------------------------------------------------

!     .. Parameters ..
      REAL (WP) :: ONE
      PARAMETER (ONE=1.0E0_WP)
      COMPLEX (WP) :: CZERO
      PARAMETER (CZERO=(0.0E0_WP,0.0E0_WP))
      COMPLEX (WP) :: CONE
      PARAMETER (CONE=(1.0E0_WP,0.0E0_WP))
      INTEGER :: IPARSIZ
      PARAMETER (IPARSIZ=13)
      INTEGER :: SPARSIZ
!Art      PARAMETER (SPARSIZ=2)
      PARAMETER (SPARSIZ=6)
!     ..
!     .. Array Arguments ..
      COMPLEX (WP) :: B(*), WRK(*), X(*)
      REAL (WP) :: SPAR(SPARSIZ)
      INTEGER :: IPAR(IPARSIZ)
!     ..
!     .. Function Arguments ..
      REAL (WP) :: PSCNRM
      EXTERNAL PSCNRM
!     ..
!     .. Subroutine Arguments ..
      EXTERNAL MATVEC, PCSUM, PRECONL, PRECONR, PROGRESS
!     ..
!     .. Local Scalars ..
      COMPLEX (WP) :: ALPHA, BETA, XI
      REAL (WP) :: EPSILON, EXITNORM, RHSSTOP
      INTEGER :: BASISDIM, BLKSZ, CNVRTX, I, IDOTS, IP, IQ, IR, ITNO, IW, &
        IXOLD, IZ, IZETA, J, J0, K, K1, LDA, LOCLEN, MAXIT, N, NPROCS, &
        PRECONTYPE, PROCID, STATUS, STEPERR, STOPTYPE
!     ..
!     .. Local Arrays ..
      COMPLEX (WP) :: DOTS(2)
!     ..
!     .. External Functions ..
      COMPLEX (WP) :: CDOTC
      REAL (WP) :: SCSETRHSSTOP
      EXTERNAL CDOTC, SCSETRHSSTOP
!     ..
!     .. External Subroutines ..
      EXTERNAL CAXPY, CCOPY, CINIT, PIMSGETPAR, STOPCRIT
!     ..

      CALL PIMSGETPAR(IPAR,SPAR,LDA,N,BLKSZ,LOCLEN,BASISDIM,NPROCS,PROCID, &
        PRECONTYPE,STOPTYPE,MAXIT,ITNO,STATUS,STEPERR,EPSILON,EXITNORM)

!  Check consistency of preconditioning and stop types
      IF (((PRECONTYPE==0) .OR. (PRECONTYPE==2)) .AND. (STOPTYPE==6)) THEN
        ITNO = 0
        STATUS = -4
        STEPERR = 0
        GO TO 9999

      END IF

!  Does not need conversion Y=Q2X for residual
      CNVRTX = 0

!  Set indices for mapping local vectors into wrk
      IR = 1
      IP = IR + LOCLEN
      IW = IP + BASISDIM*LOCLEN
      IZETA = IW + BASISDIM*LOCLEN
      IZ = IZETA + BASISDIM
      IQ = IZ + LOCLEN
      IXOLD = IQ + LOCLEN
      IDOTS = IXOLD + BASISDIM

!  Set rhs of stopping criteria
      RHSSTOP = SCSETRHSSTOP(B,WRK(IR),EPSILON,IPAR,PRECONL,PSCNRM)

!  1. r=Q1(b-AQ2x)
      IF (STOPTYPE/=6) THEN
        IF (PRECONTYPE==0) THEN
!     r=b-Ax
          CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==1) THEN
!     r=Q1(b-Ax)
          CALL CCOPY(LOCLEN,B,1,WRK(IZ),1)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IZ),1)
          CALL PRECONL(WRK(IZ),WRK(IR),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
!     r=b-AQ2x
          CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL PRECONR(X,WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(b-AQ2x)
          CALL CCOPY(LOCLEN,B,1,WRK(IP),1)
          CALL PRECONR(X,WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IP),1)
          CALL PRECONL(WRK(IP),WRK(IR),IPAR)
        END IF

      ELSE
!     r has been set to Qb in the call to dsetrhsstop
        IF (PRECONTYPE==1) THEN
!     r=Q1(b-Ax)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL PRECONL(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(b-AQ2x)
          CALL PRECONR(X,WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IW),IPAR)
          CALL PRECONL(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)
        END IF

      END IF

!  Loop
      STATUS = 0
      EXITNORM = -ONE
      STEPERR = -1
      DO ITNO = 1, MAXIT

!  2. p(1)=r
        CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IP),1)

        K = 0
        DO J = 1, BASISDIM
          J0 = J - 1

!  3. w(j)=Q1AQ2p(j)
          IF (PRECONTYPE==0) THEN
            CALL MATVEC(WRK(IP+K),WRK(IW+K),IPAR)

          ELSE IF (PRECONTYPE==1) THEN
            CALL MATVEC(WRK(IP+K),WRK(IZ),IPAR)
            CALL PRECONL(WRK(IZ),WRK(IW+K),IPAR)

          ELSE IF (PRECONTYPE==2) THEN
            CALL PRECONR(WRK(IP+K),WRK(IZ),IPAR)
            CALL MATVEC(WRK(IZ),WRK(IW+K),IPAR)

          ELSE IF (PRECONTYPE==3) THEN
            CALL PRECONR(WRK(IP+K),WRK(IW+K),IPAR)
            CALL MATVEC(WRK(IW+K),WRK(IZ),IPAR)
            CALL PRECONL(WRK(IZ),WRK(IW+K),IPAR)
          END IF

!  4. zeta(j)=dot(w(j),w(j))
          DOTS(1) = CDOTC(LOCLEN,WRK(IW+K),1,WRK(IW+K),1)

!  5. xi=dot(r,w(j))
          DOTS(2) = CDOTC(LOCLEN,WRK(IR),1,WRK(IW+K),1)

!  Accumulate simultaneously partial values
          CALL PCSUM(2,DOTS)
          WRK(IZETA+J0) = DOTS(1)
          XI = DOTS(2)

!  5. alpha=xi/zeta(j)
          IF (WRK(IZETA+J0)==CZERO) THEN
            STATUS = -3
            STEPERR = 5
            GO TO 9999

          END IF

          ALPHA = XI/WRK(IZETA+J0)

!  6. x=x+alpha*p(j)
          CALL CCOPY(LOCLEN,X,1,WRK(IXOLD),1)
          CALL CAXPY(LOCLEN,ALPHA,WRK(IP+K),1,X,1)

!  7. r=r-alpha*w(j)
          CALL CAXPY(LOCLEN,-ALPHA,WRK(IW+K),1,WRK(IR),1)

!  8. check stopping criterion
          CALL STOPCRIT(B,WRK(IR),WRK(IZ),X,WRK(IXOLD),WRK(IQ),RHSSTOP,CNVRTX, &
            EXITNORM,STATUS,IPAR,MATVEC,MATVEC,PRECONR,PCSUM,PSCNRM)

!  Call monitoring routine
          CALL PROGRESS(LOCLEN,ITNO,EXITNORM,X,WRK(IR),WRK(IZ))

          IF (STATUS==0) THEN
            GO TO 9999
          END IF

!  9. q=Q1AQ2r
          IF (PRECONTYPE==0) THEN
            CALL MATVEC(WRK(IR),WRK(IQ),IPAR)

          ELSE IF (PRECONTYPE==1) THEN
            CALL MATVEC(WRK(IR),WRK(IZ),IPAR)
            CALL PRECONL(WRK(IZ),WRK(IQ),IPAR)

          ELSE IF (PRECONTYPE==2) THEN
            CALL PRECONR(WRK(IR),WRK(IZ),IPAR)
            CALL MATVEC(WRK(IZ),WRK(IQ),IPAR)

          ELSE IF (PRECONTYPE==3) THEN
            CALL PRECONR(WRK(IR),WRK(IQ),IPAR)
            CALL MATVEC(WRK(IQ),WRK(IZ),IPAR)
            CALL PRECONL(WRK(IZ),WRK(IQ),IPAR)
          END IF

! 10. p(j+1)=r-sum_{i=1}^{j}{dot(q,w(i))/zeta(i)*p(i)}
          IF (J<BASISDIM) THEN

!  Compute partial inner-products
            K1 = 0
            DO I = 0, J - 1
              WRK(IDOTS+I) = CDOTC(LOCLEN,WRK(IQ),1,WRK(IW+K1),1)
              K1 = K1 + LOCLEN
            END DO

!  Accumulate simultaneously partial values
            CALL PCSUM(J,WRK(IDOTS))

!  Compute summation
            CALL CINIT(LOCLEN,CZERO,WRK(IZ),1)
            K1 = 0
            DO I = 0, J - 1
              BETA = WRK(IDOTS+I)/WRK(IZETA+I)
              CALL CAXPY(LOCLEN,BETA,WRK(IP+K1),1,WRK(IZ),1)
              K1 = K1 + LOCLEN
            END DO

!  Compute p(j+1)
            K = K + LOCLEN
            CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IP+K),1)
            CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IP+K),1)
          END IF

        END DO

      END DO

      IF (ITNO>MAXIT) THEN
        STATUS = -1
        ITNO = MAXIT
      END IF

9999  CONTINUE

      IF ((PRECONTYPE==2) .OR. (PRECONTYPE==3)) THEN
        CALL CCOPY(LOCLEN,X,1,WRK(IZ),1)
        CALL PRECONR(WRK(IZ),X,IPAR)
      END IF

!  Set output parameters
      IPAR(11) = ITNO
      IPAR(12) = STATUS
      IPAR(13) = STEPERR
      SPAR(2) = EXITNORM

      RETURN

    END SUBROUTINE PIMCRGCR
    SUBROUTINE PIMCTFQMR(X,B,WRK,IPAR,SPAR,MATVEC,PRECONL,PRECONR,PCSUM, &
        PSCNRM2,PROGRESS)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

!           PIM -- The Parallel Iterative Methods package
!           ---------------------------------------------

!                      Rudnei Dias da Cunha
!     National Supercomputing Centre and Mathematics Institute
!         Universidade Federal do Rio Grande do Sul, Brasil

!                          Tim Hopkins
!     Computing Laboratory, University of Kent at Canterbury, U.K.

! ----------------------------------------------------------------------

!     .. Parameters ..
      REAL (WP) :: ONE
      PARAMETER (ONE=1.0E0_WP)
      COMPLEX (WP) :: CZERO
      PARAMETER (CZERO=(0.0E0_WP,0.0E0_WP))
      COMPLEX (WP) :: CONE
      PARAMETER (CONE=(1.0E0_WP,0.0E0_WP))
      INTEGER :: IPARSIZ
      PARAMETER (IPARSIZ=13)
      INTEGER :: SPARSIZ
!Art      PARAMETER (SPARSIZ=2)
      PARAMETER (SPARSIZ=6)
!     ..
!     .. Array Arguments ..
      COMPLEX (WP) :: B(*), WRK(*), X(*)
      REAL (WP) :: SPAR(SPARSIZ)
      INTEGER :: IPAR(IPARSIZ)
!     ..
!     .. Function Arguments ..
      REAL (WP) :: PSCNRM2
      EXTERNAL PSCNRM2
!     ..
!     .. Subroutine Arguments ..
      EXTERNAL MATVEC, PCSUM, PRECONL, PRECONR, PROGRESS
!     ..
!     .. Local Scalars ..
      COMPLEX (WP) :: ALPHA, BETA, C, ETA, ETA0, KAPPA, RHO, RHO0, SIGMA, TAU, &
        TAU0, THETA, THETA0
      REAL (WP) :: CHANGETEST, EPSILON, EXITNORM, RHSSTOP
      INTEGER :: BASISDIM, BLKSZ, CNVRTX, ID, IM, IM0, IP, IR, IRTILDE, ITNO, &
        IV, IW, IXOLD, IY, IY0, IZ, LDA, LOCLEN, MAXIT, N, NPROCS, PRECONTYPE, &
        PROCID, STATUS, STEPERR, STOPTYPE
!     ..
!     .. Local Arrays ..
      COMPLEX (WP) :: DOTS(1)
!     ..
!     .. External Functions ..
      COMPLEX (WP) :: CDOTC
      REAL (WP) :: SCSETRHSSTOP
      EXTERNAL CDOTC, SCSETRHSSTOP
!     ..
!     .. External Subroutines ..
      EXTERNAL CAXPY, CCOPY, CINIT, PIMSGETPAR, STOPCRIT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS, REAL, SQRT
!     ..

      CALL PIMSGETPAR(IPAR,SPAR,LDA,N,BLKSZ,LOCLEN,BASISDIM,NPROCS,PROCID, &
        PRECONTYPE,STOPTYPE,MAXIT,ITNO,STATUS,STEPERR,EPSILON,EXITNORM)

!  Check consistency of preconditioning and stop types
      IF (((PRECONTYPE==0) .OR. (PRECONTYPE==2)) .AND. (STOPTYPE==6)) THEN
        ITNO = 0
        STATUS = -4
        STEPERR = 0
        GO TO 9999

      END IF

!  Needs conversion Y=Q2X for residual
      CNVRTX = 1

!  Set indices for mapping local vectors into wrk
      IR = 1
      IRTILDE = IR + LOCLEN
      IY = IRTILDE + LOCLEN
      IY0 = IY + LOCLEN
      IW = IY0 + LOCLEN
      IV = IW + LOCLEN
      ID = IV + LOCLEN
      IZ = ID + LOCLEN
      IP = IZ + LOCLEN
      IXOLD = IP + LOCLEN

!  Set rhs of stopping criteria
      RHSSTOP = SCSETRHSSTOP(B,WRK(IR),EPSILON,IPAR,PRECONL,PSCNRM2)

!  1. r=Q1(b-AQ2x)
      IF (STOPTYPE/=6) THEN
        IF (PRECONTYPE==0) THEN
!     r=b-Ax
          CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==1) THEN
!     r=Q1(b-Ax)
          CALL CCOPY(LOCLEN,B,1,WRK(IZ),1)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IW),1,WRK(IZ),1)
          CALL PRECONL(WRK(IZ),WRK(IR),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
!     r=b-AQ2x
          CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
          CALL PRECONR(X,WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(b-AQ2x)
          CALL CCOPY(LOCLEN,B,1,WRK(IP),1)
          CALL PRECONR(X,WRK(IW),IPAR)
          CALL MATVEC(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IP),1)
          CALL PRECONL(WRK(IP),WRK(IR),IPAR)
        END IF

      ELSE
!     r has been set to Qb in the call to dsetrhsstop
        IF (PRECONTYPE==1) THEN
!     r=Q1(b-Ax)
          CALL MATVEC(X,WRK(IW),IPAR)
          CALL PRECONL(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

        ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(b-AQ2x)
          CALL PRECONR(X,WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IW),IPAR)
          CALL PRECONL(WRK(IW),WRK(IZ),IPAR)
          CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)
        END IF

      END IF

!  2. w=y=r
      CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IW),1)
      CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IY),1)

!  3. v=Q1AQ2y
      IF (PRECONTYPE==0) THEN
        CALL MATVEC(WRK(IY),WRK(IV),IPAR)

      ELSE IF (PRECONTYPE==1) THEN
        CALL MATVEC(WRK(IY),WRK(IZ),IPAR)
        CALL PRECONL(WRK(IZ),WRK(IV),IPAR)

      ELSE IF (PRECONTYPE==2) THEN
        CALL PRECONR(WRK(IY),WRK(IZ),IPAR)
        CALL MATVEC(WRK(IZ),WRK(IV),IPAR)

      ELSE IF (PRECONTYPE==3) THEN
        CALL PRECONR(WRK(IY),WRK(IV),IPAR)
        CALL MATVEC(WRK(IV),WRK(IZ),IPAR)
        CALL PRECONL(WRK(IZ),WRK(IV),IPAR)
      END IF

!  4. d=0
      CALL CINIT(LOCLEN,CZERO,WRK(ID),1)

!  5. tau=||r||2
      TAU = PSCNRM2(LOCLEN,WRK(IR))

!  6. theta=eta=0
      THETA = CZERO
      ETA = CZERO

!  7. rtilde=r
      CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IRTILDE),1)

!  8. rho=dot(rtilde,r)
      DOTS(1) = CDOTC(LOCLEN,WRK(IRTILDE),1,WRK(IR),1)
      CALL PCSUM(1,DOTS)
      RHO = DOTS(1)

!  Loop
      STATUS = 0
      EXITNORM = -ONE
      STEPERR = -1
      CHANGETEST = -1
      IM0 = 1
      DO ITNO = 1, MAXIT
!  9. sigma=dot(rtilde,v)
        DOTS(1) = CDOTC(LOCLEN,WRK(IRTILDE),1,WRK(IV),1)
        CALL PCSUM(1,DOTS)
        SIGMA = DOTS(1)

! 10. alpha=rho/sigma
        IF (SIGMA==CZERO) THEN
          STATUS = -3
          STEPERR = 10
          GO TO 9999

        END IF

        ALPHA = RHO/SIGMA

! 11. y=y0-alpha*v
        CALL CCOPY(LOCLEN,WRK(IY),1,WRK(IY0),1)
        CALL CAXPY(LOCLEN,-ALPHA,WRK(IV),1,WRK(IY),1)

        DO IM = IM0, IM0 + 1

! 12. w=w-alpha*Q1AQ2y0
          IF (PRECONTYPE==0) THEN
            CALL MATVEC(WRK(IY0),WRK(IP),IPAR)

          ELSE IF (PRECONTYPE==1) THEN
            CALL MATVEC(WRK(IY0),WRK(IZ),IPAR)
            CALL PRECONL(WRK(IZ),WRK(IP),IPAR)

          ELSE IF (PRECONTYPE==2) THEN
            CALL PRECONR(WRK(IY0),WRK(IZ),IPAR)
            CALL MATVEC(WRK(IZ),WRK(IP),IPAR)

          ELSE IF (PRECONTYPE==3) THEN
            CALL PRECONR(WRK(IY0),WRK(IP),IPAR)
            CALL MATVEC(WRK(IP),WRK(IZ),IPAR)
            CALL PRECONL(WRK(IZ),WRK(IP),IPAR)
          END IF

          CALL CAXPY(LOCLEN,-ALPHA,WRK(IP),1,WRK(IW),1)

! 13. theta=||w||_2/tau0
          THETA0 = THETA
          TAU0 = TAU
          IF (TAU0==CZERO) THEN
            STATUS = -3
            STEPERR = 13
            GO TO 9999

          END IF

          THETA = PSCNRM2(LOCLEN,WRK(IW))/TAU0

! 14. c=1/sqrt(1+theta^2)
          C = CONE/SQRT(CONE+THETA*THETA)

! 15. tau=tau0*theta*c
          TAU = TAU0*THETA*C

! 16. eta=(c^2)*alpha
          ETA0 = ETA
          ETA = C*C*ALPHA

! 17. d=y0+((theta0^2)*eta0/alpha)*d
          IF (ALPHA==CZERO) THEN
            STATUS = -3
            STEPERR = 17
            GO TO 9999

          END IF

          CALL CCOPY(LOCLEN,WRK(ID),1,WRK(IP),1)
          CALL CCOPY(LOCLEN,WRK(IY0),1,WRK(ID),1)
          CALL CAXPY(LOCLEN,THETA0*THETA0*ETA0/ALPHA,WRK(IP),1,WRK(ID),1)

! 18. x=x+eta*d
          CALL CCOPY(LOCLEN,X,1,WRK(IXOLD),1)
          CALL CAXPY(LOCLEN,ETA,WRK(ID),1,X,1)

! 19. kappa=tau*sqrt(m+1)
          KAPPA = SQRT(REAL(IM+1,KIND=WP))*TAU

! 20. check stopping criterion
          IF (ABS(KAPPA)<EPSILON) THEN

!     r=Q1(b-AQ2x)
            IF (PRECONTYPE==0) THEN
!     r=b-Ax
              CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
              CALL MATVEC(X,WRK(IP),IPAR)
              CALL CAXPY(LOCLEN,-CONE,WRK(IP),1,WRK(IR),1)
              CALL CCOPY(LOCLEN,WRK(IR),1,WRK(IZ),1)

            ELSE IF (PRECONTYPE==1) THEN
!     r=Q1(b-Ax)
              CALL CCOPY(LOCLEN,B,1,WRK(IZ),1)
              CALL MATVEC(X,WRK(IP),IPAR)
              CALL CAXPY(LOCLEN,-CONE,WRK(IP),1,WRK(IZ),1)
              CALL PRECONL(WRK(IZ),WRK(IR),IPAR)

            ELSE IF (PRECONTYPE==2) THEN
!     r=b-AQ2x
              CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
              CALL PRECONR(X,WRK(IP),IPAR)
              CALL MATVEC(WRK(IP),WRK(IZ),IPAR)
              CALL CAXPY(LOCLEN,-CONE,WRK(IZ),1,WRK(IR),1)

            ELSE IF (PRECONTYPE==3) THEN
!     r=Q1(b-AQ2x)
              CALL CCOPY(LOCLEN,B,1,WRK(IR),1)
              CALL PRECONR(X,WRK(IZ),IPAR)
              CALL MATVEC(WRK(IZ),WRK(IP),IPAR)
              CALL CAXPY(LOCLEN,-CONE,WRK(IP),1,WRK(IR),1)
              CALL PRECONL(WRK(IR),WRK(IZ),IPAR)
              CALL CCOPY(LOCLEN,WRK(IZ),1,WRK(IR),1)
            END IF

            CALL STOPCRIT(B,WRK(IR),WRK(IZ),X,WRK(IXOLD),WRK(IP),RHSSTOP, &
              CNVRTX,EXITNORM,STATUS,IPAR,MATVEC,MATVEC,PRECONR,PCSUM,PSCNRM2)

!  Call monitoring routine
            CALL PROGRESS(LOCLEN,ITNO,EXITNORM,X,WRK(IR),WRK(IZ))

            IF (STATUS==0) THEN
              GO TO 9999
            END IF

          ELSE
!  Call monitoring routine
! 04.05.21 (BTD) correction: replace KAPPA by ABS(KAPPA)

            CALL PROGRESS(LOCLEN,ITNO,ABS(KAPPA),X,WRK(IR),WRK(IZ))

          END IF

!  y0=y
          CALL CCOPY(LOCLEN,WRK(IY),1,WRK(IY0),1)

        END DO

! 21. rho=dot(rtilde,w)
        RHO0 = RHO
        DOTS(1) = CDOTC(LOCLEN,WRK(IRTILDE),1,WRK(IW),1)
        CALL PCSUM(1,DOTS)
        RHO = DOTS(1)

! 22. beta=rho/rho0
        IF (RHO0==CZERO) THEN
          STATUS = -3
          STEPERR = 22
          GO TO 9999

        END IF

        BETA = RHO/RHO0

! 23. y=w+beta*y0
        CALL CCOPY(LOCLEN,WRK(IW),1,WRK(IY),1)
        CALL CAXPY(LOCLEN,BETA,WRK(IY0),1,WRK(IY),1)

! 24. v=Q1AQ2y+beta*(Q1AQ2y0+beta*v)

        IF (PRECONTYPE==0) THEN
          CALL MATVEC(WRK(IY0),WRK(IP),IPAR)

        ELSE IF (PRECONTYPE==1) THEN
          CALL MATVEC(WRK(IY0),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IP),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
          CALL PRECONR(WRK(IY0),WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IP),IPAR)

        ELSE IF (PRECONTYPE==3) THEN
          CALL PRECONR(WRK(IY0),WRK(IP),IPAR)
          CALL MATVEC(WRK(IP),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IP),IPAR)
        END IF

        CALL CAXPY(LOCLEN,BETA,WRK(IV),1,WRK(IP),1)

        IF (PRECONTYPE==0) THEN
          CALL MATVEC(WRK(IY),WRK(IV),IPAR)

        ELSE IF (PRECONTYPE==1) THEN
          CALL MATVEC(WRK(IY),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IV),IPAR)

        ELSE IF (PRECONTYPE==2) THEN
          CALL PRECONR(WRK(IY),WRK(IZ),IPAR)
          CALL MATVEC(WRK(IZ),WRK(IV),IPAR)

        ELSE IF (PRECONTYPE==3) THEN
          CALL PRECONR(WRK(IY),WRK(IV),IPAR)
          CALL MATVEC(WRK(IV),WRK(IZ),IPAR)
          CALL PRECONL(WRK(IZ),WRK(IV),IPAR)
        END IF

        CALL CAXPY(LOCLEN,BETA,WRK(IP),1,WRK(IV),1)

        IM0 = IM0 + 2

      END DO

      IF (ITNO>MAXIT) THEN
        STATUS = -1
        ITNO = MAXIT
        IF (EXITNORM==(-ONE)) THEN
          EXITNORM = KAPPA
        END IF

      END IF

9999  CONTINUE

      IF ((PRECONTYPE==2) .OR. (PRECONTYPE==3)) THEN
        CALL CCOPY(LOCLEN,X,1,WRK(IZ),1)
        CALL PRECONR(WRK(IZ),X,IPAR)
      END IF

!  Set output parameters
      IPAR(11) = ITNO
      IPAR(12) = STATUS
      IPAR(13) = STEPERR
      SPAR(2) = EXITNORM

      RETURN

    END SUBROUTINE PIMCTFQMR
