! correspondence with DDSCAT
! xi - initial guess of cxpol on input; cxpol  on output
! b  - cxe (right hand side, not changed)
! xr - A xi (use matvec)  work vector
! lda -  nat3
! ndim - nat3 ?
! nlar  work array dimension >= 12
! wrk - cxsc NOTE!!! THAT MXCXSC has to be set =10 in main DDSCAT
! maxit - itermx ? mxiter?
! nloop - itern
! tol -   tolr (input)
! tole -  achieved relative error
! ipar(12) - convergence;  ipar(12)=0 converged
! Example call from  getfml.f90 (DDSCAT)
!          nlar=10
!          call pimqmrcg(nat3,matvec,cxe,nat3,nlar,cxpol,cxscr1, cxsc, mxiter,&
!                 itern,tol,tolr,multiplications)

! IMPORTANT NOTES
!
! Change DDSCAT.f90 
!         MXCXSC=10*MXN3
! 
! (PJF) NOTE: nlar needs to be handled more gracefully i.e. test should be done if NLAR is
! sufficiently large. Nlar has to be at least 10 because I use 1 for vector xs and 9
! vectors for conjugate gradient. This is how it was implemented by PCC and AR
! 
! Bruce - we need to discuss the nlar issue for all CG routines. 
! it needs to be set in DDSCAT and transfered to getfml so one can add some tests
! this is simple

!(PJF) in pimqmrcg I assumed that A matrix is symmetric. I.e. the code  is not
! general for non-symmetric cases. I hope that A is indeed symmetric.

      subroutine pimqmrcg(ndim,matvec,                           &
                 b,lda,nlar,xi,xr,wrk,maxit,itno,tol,tole,ncompte)

!
! Interface to QMR solver
! ndim - dimension
! matvec - external Ax routine
! b - right hand side
! lda - leading dimension
! nlar - number of needed scratch vectors (9+1)
! xi - on input initial guess, on output the result
! xr - scratch vector
! wrk - work array
! maxit - maximum number of iterations
! itno - number of iterations
! tol - needed relative error
! tole - achived  relative error
! ncompte - cumber of Ax multiplications
!
! History
! (PJF) = Piotr J. Flatau
! (PCC) = P. C. Chaumet
! (AR) = A. Rahmani
! February 4, 2010  P. C. Chaumet and A. Rahmani
! May 6, 2010 (PJF) converted to Fortran90, introduce  DDPRECISION to
!                    handle single/double precision easily, 
!                    introduced pointer/target to split work array
           
!     license: GNU GPL
     USE DDPRECISION,ONLY: WP
     IMPLICIT NONE

!     .. Parameters ..
      INTEGER :: LDA,NLAR
! PARAMETER NLAR has to be 9+1=10 to accomodate vector xs

!     .. Local Scalars ..
      CHARACTER :: CMSGNM*70
      COMPLEX(WP) :: EPSILON,GAMMA,ICOMP,KAPPA,KSI,LAMBDA,MU,RHO,TAU,THETA
      REAL(WP) :: NORM,TOL,TOLE
      INTEGER :: I,ITLAST,ITNO,J,MAXIT,NCOMPTE,NDIM,NOU,NT,STATUS,STEPERR
      INTEGER :: IPAR
!     ..
!     .. Local Arrays ..
      COMPLEX(WP) :: & 
         B(LDA),     &
         DOTS(4),    &
         XI(LDA),    &
         XR(LDA)     
      COMPLEX(WP), TARGET:: WRK(LDA,NLAR)
      COMPLEX(WP), POINTER:: XS(:), WRK2(:,:)
!     ..
!     .. External Subroutines ..
      EXTERNAL pimzqmr, matvec
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC real
!     ..
! FLATAU split wrk array to 2 arrays which are needed by pimzqmr
      xs => wrk(1:lda,1)
      wrk2 => wrk(1:lda,2:nlar)
!
      nou = 0
      ncompte = 0
      ITLAST=1
10    CALL pimzqmr(xs,xi,xr,b,wrk2,norm,lda,ndim,nlar,lambda,kappa,theta,gamma, &
        ksi,rho,epsilon,mu,tau,dots,nou,nt,itno,maxit,tole,tol,status,steperr)

      IF(ITNO.GT.ITLAST)THEN
         WRITE(CMSGNM,FMT='(A,I8,A,1P,E10.3)') &
              'IT=',ITNO-2,' f.err=',TOLE
         CALL WRIMSG('QMRCCG',CMSGNM)
         ITLAST=ITNO
      ENDIF
!      print*, itno, tole

      IF (status.LT.0) THEN
        WRITE (*,fmt=*) 'stop nstat', status, steperr
        STOP
      END IF
      ncompte = ncompte + 1

      IF (nt.EQ.1) THEN
! original code
!        DO i = 1, ndim
!          xr(i) = 0._wp
!          DO j = 1, ndim
!            xr(i) = xr(i) + mat(i,j)*xi(j)
!          END DO
!        END DO
        call matvec(xi,xr,ipar)
      ELSE IF (nt.EQ.2) THEN
! Flatau. NOTE this only work for symmetric problems
! I am not 100% sure if DDSCAT is always symmetric

!original code (is this right, mat is not conjugated here?)
!        DO i = 1, ndim
!          xr(i) = 0._wp
!          DO j = 1, ndim
!            xr(i) = xr(i) + (mat(j,i))*xi(j)
!          END DO
!        END DO
        call matvec(xi,xr,ipar)
      END IF
      IF (status.NE.1) GO TO 10

      IF(STEPERR.EQ.0)THEN
!        WRITE (*,fmt=*) 'ITNO has reached MAXIT', itno, maxit
         WRITE(CMSGNM,FMT='(A,I6,A,I6)')'IT=',ITNO,' has reached MAXIT=',MAXIT
         CALL ERRMSG('FATAL','pimqmrcg',CMSGNM)
      ENDIF
! FLATAU after all is done assign xs (solution) to xi 
     xi(1:lda)=xs(1:lda)
    return
    END
!****************************************
    SUBROUTINE pimzqmr(xs,xi,xr,b,wrk,norm,lda,ndim,nlar,lambda,kappa,theta, &
        gamma,ksi,rho,epsilon,mu,tau,dots,nou,nt,itno,maxit,tole,tol,status, &
        steperr)
     USE DDPRECISION,ONLY: WP
!****************************************************************
!     Iterative solver QMR
!****************************************************************
!     Authors: P. C. Chaumet and A. Rahmani
!     Date: 04/02/2010

!     Purpose: iterative solver for linear system Ax=b. There is no
!     condition on the matrix. Notice that the products A x and At x are
!     provided by the user.

!     Reference: if you use this routine in your research, please
!     reference, as appropriate: P. C. Chaumet and A. Rahmani, Efficient
!     discrete dipole approximation for magneto-dielectric scatterers
!     Opt. Lett. 34, 917 (2009). R. D. Da Cunha and T. Hopkins,
!     Appl. Numer. Math. 19, 33 (1995).
! History
! (PJF) = Piotr J. Flatau
! (PCC) = P. C. Chaumet
! (AR) = A. Rahmani
! Originally written by  R. D. Da Cunha and T. Hopkins
! February 4, 2010 Modified by P. C. Chaumet and A. Rahmani
! May 6, 2010  (PJF) converted to Fortran90, introduce  DDPRECISION to
!                    handle sing/double precision easily, introduced pointer/target
!                    to split work array

      IMPLICIT NONE

!     .. Array Arguments ..
      COMPLEX (wp) :: b(lda), dots(4), wrk(lda,nlar), xi(lda), xr(lda), &
        xs(lda)
!     ..
!     .. Local Scalars ..
      COMPLEX (wp) :: absgamma2, abstau02, den, epsilon0, gamma0, kappa0, &
        ksi0, lambda0, mu0, rho0, tau0, tmp1
      INTEGER :: i
!     ..
!     .. Local Arrays ..
!     ..

!     .. Scalar Arguments ..
      COMPLEX (wp) :: epsilon, gamma, kappa, ksi, lambda, mu, rho, tau, theta
      REAL (wp) :: norm, tol, tole
      INTEGER :: itno, lda, maxit, ndim, nlar, nou, nt, status, steperr
!     ..
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC abs, conjg, sqrt
!     ..
      IF (nou.EQ.0) GO TO 10
      IF (nou.EQ.1) GO TO 20
      IF (nou.EQ.2) GO TO 30
      IF (nou.EQ.3) GO TO 40
      IF (nou.EQ.4) GO TO 50
      IF (nou.EQ.5) GO TO 100

!     1. lambda=1, kappa=-1, theta=-1
10    lambda = (1._wp,0._wp)
      kappa = -(1._wp,0._wp)
      theta = -(1._wp,0._wp)
      norm = 0._wp
      DO i = 1, ndim
        norm = norm + b(i)*conjg(b(i))
      END DO
      norm = sqrt(norm)

!     Loop
      status = 0
      steperr = -1
      itno = 0

!     2. wtilde=vtilde=r=b-Ax
!     r=b-Ax

!     A*x=wrk(3)
      nou = 1
      nt = 1
      DO i = 1, ndim
        xs(i) = xi(i)
      END DO

!     compute A*xi

      RETURN

20    DO i = 1, ndim
        wrk(i,3) = xr(i)
        wrk(i,1) = b(i) - wrk(i,3)
        wrk(i,7) = wrk(i,1)
        wrk(i,8) = wrk(i,1)
      END DO

!     3. p=q=d=s=0
      DO i = 1, ndim
        wrk(i,2) = 0._wp
        wrk(i,4) = 0._wp
        wrk(i,5) = 0._wp
        wrk(i,6) = 0._wp
      END DO

!     4. gamma=||vtilde||_{2}, ksi=||wtilde||_{2},
!     rho=wtilde^{T}vtilde, epsilon=(A^{T}wtilde)^{T}vtilde, mu=0

      dots(1) = 0._wp
      dots(2) = 0._wp
      dots(3) = 0._wp

      DO i = 1, ndim
        dots(1) = dots(1) + wrk(i,7)*conjg(wrk(i,7))
        dots(2) = dots(2) + wrk(i,8)*conjg(wrk(i,8))
        dots(3) = dots(3) + wrk(i,7)*wrk(i,8)
      END DO
!     Compute A^{T}wtilde
!     CALL TMATVEC(WRK(IWTILDE),WRK(IATWTILDE),IPAR)
      DO i = 1, ndim
        xi(i) = wrk(i,8)
      END DO

      nt = 2
      nou = 2
      RETURN


30    dots(4) = 0._wp
      DO i = 1, ndim
        wrk(i,9) = xr(i)
        dots(4) = dots(4) + wrk(i,7)*wrk(i,9)
      END DO

!     Accumulate simultaneously partial inner-products
!     CALL PDZSUM(4,DOTS)

      gamma = sqrt(dots(1))
      ksi = sqrt(dots(2))
      rho = dots(3)
      epsilon = dots(4)
      mu = 0._wp
!     5. tau=epsilon/rho
      IF (rho.EQ.0._wp) THEN
        itno = 0
        status = -3
        steperr = 5
        GO TO 200
      END IF

      tau = epsilon/rho
100   itno = itno + 1

!     6. p=1/gamma*vtilde-mu*p

      IF (gamma.EQ.0._wp) THEN
        status = -3
        steperr = 6
        GO TO 200

      END IF
      DO i = 1, ndim
        wrk(i,2) = wrk(i,7)/gamma - mu*wrk(i,2)
      END DO

!     7. q=1/ksi*A^{T}wtilde-(gamma*mu)/ksi*q
      IF (ksi.EQ.0._wp) THEN
        status = -3
        steperr = 7
        GO TO 200
      END IF

      DO i = 1, ndim
        wrk(i,4) = (wrk(i,9)-gamma*mu*wrk(i,4))/ksi
      END DO

!     8. vtilde=Ap-tau/gamma*vtilde
      DO i = 1, ndim
        xi(i) = wrk(i,2)
      END DO

      nt = 1
      nou = 3
      RETURN

40    DO i = 1, ndim
        wrk(i,3) = xr(i)
        wrk(i,7) = wrk(i,3) - tau/gamma*wrk(i,7)
      END DO

!     9. wtilde=q-tau/ksi*wtilde

      DO i = 1, ndim
        wrk(i,8) = wrk(i,4) - tau/ksi*wrk(i,8)
      END DO

!     11. gamma=||vtilde||_{2}, ksi=||wtilde||_{2},
!     rho=wtilde^{T}vtilde, epsilon=(A^{T}wtilde)^{T}vtilde

      dots(1) = 0._wp
      dots(2) = 0._wp
      dots(3) = 0._wp

      DO i = 1, ndim
        dots(1) = dots(1) + wrk(i,7)*conjg(wrk(i,7))
        dots(2) = dots(2) + wrk(i,8)*conjg(wrk(i,8))
        dots(3) = dots(3) + wrk(i,7)*wrk(i,8)
      END DO

!     Compute A^{T}wtilde
      DO i = 1, ndim
        xi(i) = wrk(i,8)
      END DO

      nt = 2
      nou = 4
      RETURN
50    dots(4) = 0._wp
      DO i = 1, ndim
        wrk(i,9) = xr(i)
        dots(4) = dots(4) + wrk(i,7)*wrk(i,9)
      END DO

!     Accumulate simultaneously partial inner-products
!     CALL PDZSUM(4,DOTS)

      gamma0 = gamma
      gamma = sqrt(dots(1))
      ksi0 = ksi
      ksi = sqrt(dots(2))
      rho0 = rho
      rho = dots(3)
      epsilon0 = epsilon
      epsilon = dots(4)
!     12. mu=(gamma0*ksi0*rho)/(gamma*tau*rho0)

      den = gamma*tau*rho0
      IF (den.EQ.0._wp) THEN
        status = -3
        steperr = 12
        GO TO 200

      END IF
      mu0 = mu
      mu = (gamma0*ksi0*rho)/den

!     13. tau=epsilon/rho-gamma*mu
      IF (rho.EQ.0._wp) THEN
        status = -3
        steperr = 13
        GO TO 200

      END IF
      tau0 = tau
      tau = epsilon/rho - gamma*mu
!     14. theta=(|tau0|^2*(1-lambda))/(lambda*|tau|^2+|gamma|^2)

      abstau02 = abs(tau0)**2._wp
      absgamma2 = abs(gamma)**2._wp
      den = lambda*abstau02 + absgamma2
      IF (den.EQ.0._wp) THEN
        status = -3
        steperr = 14
        GO TO 200

      END IF
      theta = (abstau02*((1._wp,0._wp)-lambda))/den

!     15. kappa=(-gamma0*CONJG(tau0)*kappa0)/(gamma0*|tau|^2+|gamma|^2)
      kappa0 = kappa
      kappa = -(gamma0*conjg(tau0)*kappa0)/den

!     16. lambda=(lambda0*|tau0|^2)/(gamma0*|tau|^2+|gamma|^2)
      lambda0 = lambda
      lambda = lambda0*abstau02/den

!     17. d=theta*d+kappa*p

      DO i = 1, ndim
        wrk(i,5) = theta*wrk(i,5) + kappa*wrk(i,2)
      END DO

!     18. s=theta*s+kappa*A*p

      DO i = 1, ndim
        wrk(i,6) = theta*wrk(i,6) + kappa*wrk(i,3)
      END DO

!     19. x=x+d

      DO i = 1, ndim
        xs(i) = xs(i) + wrk(i,5)
      END DO

!     20. r=r-s

      DO i = 1, ndim
        wrk(i,1) = wrk(i,1) - wrk(i,6)
      END DO

!     criterion to stop
      tmp1 = 0._wp
      DO i = 1, ndim
        tmp1 = tmp1 + wrk(i,1)*conjg(wrk(i,1))
      END DO
      tole = sqrt(abs(tmp1))/norm

      IF (tole.LE.tol) THEN
        nou = 5
        status = 1
        RETURN
      END IF

      IF (itno.GT.maxit) THEN
        status = 1
        steperr = 0
        RETURN
      END IF
      GO TO 100

200   RETURN

    END SUBROUTINE pimzqmr
