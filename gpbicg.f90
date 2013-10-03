    SUBROUTINE gpbicg(xi,xr,b,lda,ndim,nlar,nou,wrk,itno,maxit,tol,norm,alpha, &
        beta,eta,dzeta,r0rn,status,steperr, residu)
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE

!****************************************************************
!     Iterative solver GPBICG
!****************************************************************
!     Authors: P. C. Chaumet and A. Rahmani
!     Date: 04/02/2010
!     Purpose: iterative solver for linear system Ax=b. There is no
!     condition on the matrix. Notice that the product A x is provided
!     by the user.
!     Reference: if you use this routine in your research, please
!     reference, as appropriate: P. C. Chaumet and A. Rahmani, Efficient
!     discrete dipole approximation for magneto-dielectric scatterers
!     Opt. Lett. 34, 917 (2009). J. Tang, Y. Shen, Y. Zheng, and D. Qiu,
!     Coastal Eng. 51, 143 (2004).
!     license: GNU GPL
!     We cannot guarantee correctness of the programs...
!     Main program and example. The main program if provided for testing
!     purposes and should be commented out to use only the routine GPBICG
!     We want to solve Ax=b
!     mat:  matrix A
!     lda: Input: leading dimension array of the matrix
!     ndim: Input: dimension  of the matrix: ndim.le.lda
!     nlar: Input: size of the work vector: nlar.ge.12
!     MAXIT: Input: Maximum of iteration for the iterative solver
!     NLOOP: Output: number of iteration for the iterative solver Should
!     be initialize to 0.
!     ncompte: number of Ax products.
!     STATUS: Output: if STATUS.lt.0 a problem occured in GPBICG
!     STATUS=1 the requested tolerance has been reached or the maximum n
!     iterations has been reached
!     STEPERR: Output: if STEPERR.lt.0: indicates where the problem
!     occurs in GPBICG. STEPERR=0 the maximum number of iterations has
!     been reached.  STEPERR=-1 routine completed witho

!     tol: Input: tolerance requested by the user. At the end we have:
!     r=||Ax-b||/||b||<tol
!     b: Input: right-hand side of Ax=b
!     norm: Output: norm of b
!     xi: Input: initial guess; output:solution of the linear equation
!     xr: Input: xr = A xi
!     NOU: Local integer used by GPBICG. Should be initialized to 0.
!     ALPHA,BETA,ETA,DZETA,R0RN: Local complex needs for GPBICG
!     WRK: local array used by for GPBICG
! History:
! (PJF) = Piotr J. Flatau
! May 6, 2010 - added wrmsg, converted to Fortran 90, added single/double precision kinds
!               corrected "40 elseif" statement such that it is on separate lines 

!     .. Local Scalars ..

      COMPLEX (wp) :: ctmp, ctmp1, ctmp2, ctmp3, ctmp4, ctmp5
      REAL (wp) :: residu
      INTEGER :: ii
      character*70 CMSGNM
!     ..
!     .. Scalar Arguments ..
      COMPLEX (wp) :: alpha, beta, dzeta, eta, r0rn
      REAL (wp) :: norm, tol
      INTEGER :: itno, lda, maxit, ndim, nlar, nou, status, steperr
!     ..
!     .. Array Arguments ..
      COMPLEX (wp) :: b(lda), wrk(lda,nlar), xi(lda), xr(lda)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC abs, conjg, sqrt
!     ..
      IF (nou.EQ.1) GO TO 10
      IF (nou.EQ.2) GO TO 20
      IF (nou.EQ.3) GO TO 30
      IF (nou.EQ.4) GO TO 40

      status = 0
      steperr = -1

!     Column index of the various variables stored in WRK array.
!     1:r0
!     2:p
!     3:r
!     4:y
!     5:t
!     6:Ap
!     7:At
!     8:u
!     9:w
!     10:z
!     11:x
!     12:t old

!     Compute norm and Ax0 (x0 initial guess)
      norm = 0._wp
      DO ii = 1, ndim
        wrk(ii,11) = xi(ii)
        norm = norm + abs(b(ii))**2._wp
      END DO
      norm = sqrt(norm)

      nou = 1
      RETURN

!     initialize r0=b-Ax0,rOt=rO,p0=v0=d0
10    DO ii = 1, ndim
        wrk(ii,1) = b(ii) - xr(ii)
        wrk(ii,2) = 0._wp
        wrk(ii,3) = wrk(ii,1)
        wrk(ii,5) = 0._wp
        wrk(ii,9) = 0._wp
        wrk(ii,8) = 0._wp
        wrk(ii,10) = 0._wp
      END DO
!     Compute the initial residue
      ctmp = 0._wp
      DO ii = 1, ndim
        ctmp = ctmp + wrk(ii,1)*conjg(wrk(ii,1))
      END DO

      r0rn = ctmp

!     initialize rho,alpha,w=1,tau=norm,theta=0,eta=0
      beta = 0._wp

!     begin the iteration sequence
      itno = -1

! btd add wrimsg call
      residu=1._wp
      WRITE(CMSGNM,FMT='(A,I8,A,1P,E10.3)') &
                 'IT=',ITNO+1,' f.err=',RESIDU
      CALL WRIMSG('GPBICG',CMSGNM)

!      WRITE (*,fmt=*) 'Residu Initial', abs(sqrt(ctmp))/norm

100   itno = itno + 1

!     compute p=r+beta*(p-u)
      DO ii = 1, ndim
        wrk(ii,2) = wrk(ii,3) + beta*(wrk(ii,2)-wrk(ii,8))
        xi(ii) = wrk(ii,2)
      END DO
      nou = 2
      RETURN

!     compute Ap
20    DO ii = 1, ndim
        wrk(ii,6) = xr(ii)
      END DO

!     compute alpha=r0r/r0Ap
      ctmp = 0._wp
      DO ii = 1, ndim
        ctmp = ctmp + conjg(wrk(ii,1))*wrk(ii,6)
      END DO

      IF (ctmp.EQ.0._wp) THEN
        status = -1
        steperr = 1
        RETURN
      END IF

      alpha = r0rn/ctmp

!     compute y=t-r-alpha*w+alpha*Ap et de t=r-alpha AP
      DO ii = 1, ndim
        wrk(ii,4) = wrk(ii,5) - wrk(ii,3) - alpha*wrk(ii,9) + alpha*wrk(ii,6)
        wrk(ii,12) = wrk(ii,5)
        wrk(ii,5) = wrk(ii,3) - alpha*wrk(ii,6)
        xi(ii) = wrk(ii,5)
      END DO
      nou = 3
      RETURN

!     compute At
30    DO ii = 1, ndim
        wrk(ii,7) = xr(ii)
      END DO

!     compute dzeta and eta

      IF (itno.EQ.0) THEN

        eta = 0._wp
        dzeta = 0._wp
        ctmp = 0._wp
        DO ii = 1, ndim
          dzeta = dzeta + conjg(wrk(ii,7))*wrk(ii,5)
          ctmp = ctmp + conjg(wrk(ii,7))*wrk(ii,7)
        END DO

        IF (ctmp.EQ.0._wp) THEN
          status = -1
          steperr = 2
          RETURN
        END IF
        dzeta = dzeta/ctmp

      ELSE

        ctmp1 = 0._wp
        ctmp2 = 0._wp
        ctmp3 = 0._wp
        ctmp4 = 0._wp
        ctmp5 = 0._wp

        DO ii = 1, ndim
          ctmp1 = ctmp1 + conjg(wrk(ii,7))*wrk(ii,7)
          ctmp2 = ctmp2 + conjg(wrk(ii,4))*wrk(ii,4)
          ctmp3 = ctmp3 + conjg(wrk(ii,7))*wrk(ii,4)
          ctmp4 = ctmp4 + conjg(wrk(ii,7))*wrk(ii,5)
          ctmp5 = ctmp5 + conjg(wrk(ii,4))*wrk(ii,5)
        END DO

        ctmp = ctmp1*ctmp2 - ctmp3*conjg(ctmp3)

        IF (ctmp.EQ.0._wp) THEN
          status = -1
          steperr = 3
          RETURN
        END IF

        dzeta = (ctmp2*ctmp4-ctmp5*ctmp3)/ctmp
        eta = (ctmp1*ctmp5-conjg(ctmp3)*ctmp4)/ctmp

      END IF

!     compute u
      DO ii = 1, ndim
        wrk(ii,8) = dzeta*wrk(ii,6) + eta*(wrk(ii,12)-wrk(ii,3)+beta*wrk(ii,8) &
          )
      END DO

!     compute z
      DO ii = 1, ndim
        wrk(ii,10) = dzeta*wrk(ii,3) + eta*wrk(ii,10) - alpha*wrk(ii,8)
      END DO

!     compute x and r
      residu = 0._wp
      DO ii = 1, ndim
        wrk(ii,11) = wrk(ii,11) + alpha*wrk(ii,2) + wrk(ii,10)
        wrk(ii,3) = wrk(ii,5) - eta*wrk(ii,4) - dzeta*wrk(ii,7)
        residu = residu + abs(wrk(ii,3))**2._wp
      END DO
      residu = sqrt(residu)/norm


!flatau add wrimsg call
!btd modify
           WRITE(CMSGNM,FMT='(A,I8,A,1P,E10.3)') &
                 'IT=',itno+1,' f.err=',residu
           CALL WRIMSG('GPBICG',CMSGNM)

      IF (residu.LE.tol) THEN
        status = 1
        DO ii = 1, ndim
          xi(ii) = wrk(ii,11)
        END DO
        nou = 4
        RETURN
      END IF
40    CONTINUE
!     compute beta
      ctmp = 0._wp
      DO ii = 1, ndim
        ctmp = ctmp + conjg(wrk(ii,1))*wrk(ii,3)
      END DO

      IF (r0rn.EQ.0._wp) THEN
        status = -1
        steperr = 4
        RETURN
      END IF

      beta = alpha*ctmp/dzeta/r0rn
      r0rn = ctmp

!     compute w
      DO ii = 1, ndim
        wrk(ii,9) = wrk(ii,7) + beta*wrk(ii,6)
      END DO

      IF (itno.LE.maxit) GO TO 100

      status = 1
      steperr = 0
      DO ii = 1, ndim
        xi(ii) = wrk(ii,11)
      END DO

    END SUBROUTINE gpbicg
