!notes for DDSCAT integration
!changed petr to petr90

! 
!                              CCGPAK 2.0
!  Conjugate gradient package for solving complex matrix equations (Fortran90)
!                           Piotr J. Flatau
!                        last change January 8, 2013
!
! If you use this library in publication please reference:
! Flatau, P. J., 2012, CCGPAC 2.0 - Conjugate gradient package for solving complex matrix equations, google code
! http://code.google.com/p/conjugate-gradient-lib/
!
! Copyright (C) 2012 P.J. Flatau
! This code is covered by the GNU General Public License.
!  ----------------------------------------------------------------
! library of CG implementations
!  CGSQR     - Complex conjugate gradient-squared algorithm of Sonneveld,
!  CGSTAB    - Complex conjugate gradient squared stabilized
!  PETR      - Petravic,M and Kuo-Petravic,G. 1979, JCP, 32,263
!  CORS      - The BiCOR and CORS Iterative Algorithms 
!



    SUBROUTINE PETR90ver2(x,b,n,matvec,cxsc, mxcxsc, cmatvec, cg)
      USE ddprecision, ONLY : wp
      use cgmodule
      implicit none
      type (cgstruct) cg

! solves Ax=b

! Input:
! B(n) - RHS (doesn't change on output)
! MATVEC -- name of procedure for computing y=Ax
! CMATVEC -- name of procedure for computing y= conj(A') x

! Output:
! X(n)

! Reference:
!    Petravic,M and Kuo-Petravic,G. 1979, JCP, 32,263
! History:
!    Late 80's  (Coded by Draine)
!    Changes to include "cprod" (PJF+BTD)
!    95.06.01 (PJF) re-written to conform to PIM
!    95.08.11 (BTD) minor editing (comments, etc.)
!    08.05.12 (BTD) change CDUMMY -> CDUMMY(1) following suggestion
!                   by Art Lazanoff, NASA Ames
!    2012 July 30 (PJF) converted to Fortran90 syntax, removed all PIM
! end history
!=======================================================================

!     .. Parameters ..
!     .. Scalar Arguments ..
      REAL (wp) ::epsilon_err
      INTEGER :: iter, maxit,  n, alloc_error, dealloc_error,ioerr
!     .. Array Arguments ..
      COMPLEX (wp) :: x(n), b(n)
!     .. Local Scalars ..
      REAL (wp) ::  rnorm,bnorm, gigi, gi1gi1, qiqi, betai1, alphai
      integer maxiter, itno
      integer mxcxsc
      complex(wp),target:: cxsc(n,7)
!     .. Local Arrays ..

     COMPLEX (wp), pointer :: qi(:), gi(:), pi(:), cr(:), axi(:), r(:), ace(:)

!     .. External Functions ..
!     .. External Subroutines ..
      EXTERNAL matvec, cmatvec
!     .. Intrinsic Functions ..
      INTRINSIC sqrt
      maxiter=cg%maxit
      epsilon_err=cg%epsilon_err
      ioerr=cg%ioerr

      if(n*7 .gt. mxcxsc) stop 'stop in cg routine petr90 n*7> mxcxsc not enough memory'
      qi  => cxsc(:,1)
      gi  => cxsc(:,2)
      pi  => cxsc(:,3)
      cr  => cxsc(:,4)
      axi => cxsc(:,5)
      r   => cxsc(:,6)
      ace => cxsc(:,7)


!      allocate(qi(n), gi(n), pi(n), cr(n), axi(n), r(n), ace(n),STAT=alloc_error)
!        IF (alloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)') 'Allocation Error Detected in conjugate gradient petr'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' petr '
!        END IF
! Compute conjg(A')*B
      CALL CMATVEC(B,ACE,n)
      BNORM=dot_product(b(1:n),b(1:n))
      GI(1:n)=ACE(1:n)
      PI(1:n)=GI(1:n)
! Compute |QI>=A|PI>
      CALL MATVEC(PI,QI,n)
      GIGI=dot_product(GI(1:n),GI(1:n))
      QIQI=dot_product(qi(1:n),qi(1:n))
      ALPHAI=GIGI/QIQI

      X(1:n)=X(1:n)+ALPHAI*PI(1:n)
! compute |AX1>:
      CALL MATVEC(X,AXI,n)
      DO ITNO=1, maxiter
! Transfer <GI|GI> -> <GI-1|GI-1>:
         GI1GI1=GIGI
! compute |GI>=AC|E>-AC|AXI>:
         CALL CMATVEC(AXI,GI, n)
! Compute GIGI=<GI|GI>:
         GI(1:n)=ace(1:n)-gi(1:n)
         GIGI=dot_product(gi(1:N),gi(1:n))
! Compute BETAI-1=<GI|GI>/<GI-1|GI-1>:
         BETAI1=GIGI/GI1GI1
! Compute |PI>=|GI>+BETAI-1|PI-1>:
         PI(1:n)=GI(1:n)+BETAI1*PI(1:n)
! Compute |QI>=A|PI>:
         CALL MATVEC(PI, QI, n)
! Compute <QI|QI>:
         QIQI=dot_product(QI(1:n),QI(1:n))
! Compute ALPHAI=<GI|GI>/<QI|QI>:
         ALPHAI=GIGI/QIQI
! Compute |XI+1>=|XI>+ALPHAI*|PI>:
            X(1:n)=X(1:n)+ALPHAI*PI(1:n)
! Except every 10TH iteration, compute |AXI+1>=|AXI>+ALPHAI*|QI>
! Draine: warning this part perhaps not checked
         IF(ITNO/=10*(ITNO/10))THEN
         AXI(1:n)=AXI(1:n)+ALPHAI*QI(1:n)
         ELSE
            CALL MATVEC(X,AXI,n)
         ENDIF
! Compute residual vector |RI>=A|XI>-|B>
         R(1:n)=AXI(1:n)-B(1:n)
         rnorm=dot_product(r(1:n),r(1:n))
         if(cg%print.eq.'print')  WRITE (*,fmt=*) 'sqrt(rnorm/bnorm)= ', itno, sqrt(rnorm/bnorm)
         IF ((sqrt(rnorm/bnorm)).LT.epsilon_err) GO TO 50
      ENDDO
 50 continue
      cg%itno=itno
!      deallocate(qi, gi, pi, cr, axi, r, ace, STAT=dealloc_error)
!        IF (dealloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)') 'Deallcation Error Detected in conjugate gradient petr'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' petr '
!        END IF
      RETURN
    END SUBROUTINE PETR90ver2


    SUBROUTINE cgsqr90ver2(cx,cy,n,cxsc, mxcxsc, matvec, cg)
      USE ddprecision, ONLY : wp
      use cgmodule
      implicit none
      type (cgstruct) cg
!    CGSQR: Complex conjugate gradient-squared algorithm of Sonneveld,
!           for unsymmetric cases.  
!      H. A. van der Vorst,
!      "BI-CGSTAB: A FAST AND SMOOTHLY CONVERGING VARIANT OF BI-CG FOR
!      THE SOLUTION OF NONSYMMETRIC LINEAR SYSTEMS", SIAM J. Sci. Stat.
!      Comput., V.13, #2, pp. 631-644, 1992.
!  History: 92/09/03 
!     .. Parameters ..
!     .. Scalar Arguments ..
      REAL (wp) :: epsilon_err
      INTEGER :: iter, maxit, n, alloc_error, dealloc_error,ioerr
!     .. Array Arguments ..
      COMPLEX (wp) ::  cx(n), cy(n)
!     .. Local Scalars ..
      COMPLEX (wp) :: alpha, beta, rho, rho0
      REAL (wp) :: rnorm, ynorm
      integer mxcxsc
      complex(wp),target:: cxsc(n,9)
 !     .. Local Arrays ..
      COMPLEX (wp),pointer:: aw(:), cax(:), p(:), q(:), r(:), r0(:), u(:), v(:), w(:)
!     .. External Functions ..
!     .. External Subroutines ..
      EXTERNAL matvec
!     .. Intrinsic Functions ..
      INTRINSIC cmplx, sqrt
      maxit=cg%maxit
      epsilon_err=cg%epsilon_err
      ioerr=cg%ioerr

      if(n*9 .gt. mxcxsc) stop 'stop in cg routine cgsqr90ver90 n*9 > mxcxsc not enough memory'
      aw  => cxsc(:,1)
      cax  => cxsc(:,2)
      p  => cxsc(:,3)
      q  => cxsc(:,4)
      r => cxsc(:,5)
      r0   => cxsc(:,6)
      u => cxsc(:,7)
      v => cxsc(:,8)
      w => cxsc(:,9)

!      allocate (aw(n), cax(n), p(n), q(n), r(n), r0(n), u(n), v(n), w(n),STAT=alloc_error)
!        IF (alloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)') 'Allocation Error Detected in conjugate gradient cgsqr'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' cgsqr '
!        END IF

      ynorm=dot_product(cy(1:n),cy(1:n))
!      CALL prod(n,'u',cx,'u',cax)
      call matvec(cx,cax,n)
      r(1:n) = cy(1:n) - cax(1:n)
      r0(1:n) = r(1:n)
      p(1:n) = cmplx(0.0_wp,0.0_wp,kind=wp)
      q(1:n) = cmplx(0.0_wp,0.0_wp,kind=wp)
      rho0 = cmplx(1.0_wp,1.0_wp,kind=wp)
      DO iter = 1, maxit
        rho = sum(conjg(r0(1:n))*r(1:n))
        beta = rho/rho0
        u(1:n) = r(1:n) + beta*q(1:n)
        p(1:n) = u(1:n) + beta*(q(1:n)+beta*p(1:n))
!        CALL prod(n,'u',p,'u',v)
        call matvec(p,v,n)
        alpha=rho/sum(conjg(r0(1:n))*v(1:n))
        q(1:n) = u(1:n) - alpha*v(1:n)
        w(1:n) = u(1:n) + q(1:n)
        cx(1:n) = cx(1:n) + alpha*w(1:n)
!        CALL prod(n,'u',w,'u',aw)
        call matvec(w,aw,n)
        r(1:n) = r(1:n) - alpha*aw(1:n)
        rnorm=dot_product(r(1:n),r(1:n))
        if(cg%print.eq.'print') WRITE (*,fmt=*) 'sqrt(||r||/||y||) =  ', iter, sqrt(rnorm/ynorm)
        cg%itno=iter
        IF (sqrt(rnorm/ynorm).LT.epsilon_err) GO TO 60
        rho0 = rho
      END DO
60    CONTINUE
!      deallocate (aw, cax, p, q, r, r0, u, v, w, STAT=dealloc_error)
!        IF (dealloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)') 'Deallcation Error Detected in conjugate gradient cgsqr'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' cgsqr '
!        END IF

      RETURN
    END SUBROUTINE cgsqr90ver2

    SUBROUTINE cgstab90ver2(cx,cy,n,cxsc, mxcxsc, matvec, cg)
      USE ddprecision, ONLY : wp
      use cgmodule
      implicit none
      type (cgstruct) cg
!  CGSTAB:  Complex conjugate gradient squared stabilized
!           encoded from HAV's algorithm Bi-CGSTAB.
! HISTORY:  92/09/08 (TLS)
! TLS= T. L. Schneider, Colorado State University, Dept. Atmos. Sci.
!      For Collins, Colo.
!     .. Parameters ..
!     .. Scalar Arguments ..
      REAL (wp) ::epsilon_err
      INTEGER :: iter, maxit, n, alloc_error, dealloc_error,ioerr
!     .. Array Arguments ..
      COMPLEX (wp) :: cx(n), cy(n)
!     .. Local Scalars ..
      COMPLEX (wp) :: alpha, beta, omega, rho, rho0
      REAL (wp) :: rnorm, ynorm
      integer mxcxsc
      complex(wp),target:: cxsc(n,7)
!     .. Local Arrays ..
      COMPLEX (wp), pointer,dimension(:)  :: cax, p, r, r0, s, t, v
!     .. External Functions ..
!     .. External Subroutines ..
      EXTERNAL matvec
!     .. Intrinsic Functions ..
      INTRINSIC cmplx, sqrt
      maxit=cg%maxit
      epsilon_err=cg%epsilon_err
      ioerr=cg%ioerr

      if(n*7 .gt. mxcxsc) stop 'stop in cg routine cgstab90ver n*9> mxcxsc not enough memory'
      cax  => cxsc(:,1)
      p  => cxsc(:,2)
      r  => cxsc(:,3)
      r0  => cxsc(:,4)
      s => cxsc(:,5)
      t   => cxsc(:,6)
      v => cxsc(:,7)

!      allocate(cax(n), p(n), r(n), r0(n), s(n), t(n), v(n),STAT=alloc_error)
!        IF (alloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)') 'Allocation Error Detected in conjugate gradient cgstab'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' cgstab '
!        END IF

      ynorm=dot_product(cy(1:n),cy(1:n))
!      CALL prod(n,'u',cx,'u',cax)
      call matvec(cx,cax,n)
      r(1:n) = cy(1:n) - cax(1:n)
      r0(1:n) = r(1:n)
      p(1:n) = cmplx(0.0_wp,0.0_wp,kind=wp)
      v(1:n) = cmplx(0.0_wp,0.0_wp,kind=wp)
      rho0 = cmplx(1.0_wp,1.0_wp,kind=wp)
      alpha = cmplx(1.0_wp,1.0_wp,kind=wp)
      omega = cmplx(1.0_wp,1.0_wp,kind=wp)
      DO iter = 1, maxit
        rho = sum(conjg(r0(1:n))*r(1:n))
        beta = (rho/rho0)*(alpha/omega)
        p(1:n) = r(1:n) + beta*(p(1:n)-omega*v(1:n))
!        CALL prod(n,'u',p,'u',v)
        call matvec(p,v,n)
        alpha = rho/sum(conjg(r0(1:n))*v(1:n))
        s(1:n) = r(1:n) - alpha*v(1:n)
!        CALL prod(n,'u',s,'u',t)
        call matvec(s,t,n)
        omega = sum(conjg(t(1:n))*s(1:n))/dot_product(t(1:n),t(1:n))
        cx(1:n) = cx(1:n) + alpha*p(1:n) + omega*s(1:n)
        r(1:n) = s(1:n) - omega*t(1:n)
        rnorm=dot_product(r(1:n),r(1:n))
        if(cg%print.eq.'print')  WRITE (*,fmt=*) 'sqrt(||r||/||y||) =  ', iter, sqrt(rnorm/ynorm)
        cg%itno=iter
        IF (sqrt(rnorm/ynorm).LT.epsilon_err) GO TO 60
        rho0 = rho
      END DO
60    CONTINUE
!      deallocate(cax, p, r, r0, s, t, v, STAT=dealloc_error)
!        IF (dealloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)') 'Deallcation Error Detected in conjugate gradient cgstab'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' cgstab '
!        END IF
      RETURN
    END SUBROUTINE cgstab90ver2

        SUBROUTINE cors90ver2(x,b,n,cxsc, mxcxsc, matvec, cg)
!
!	CORS routine for solving a linear system Ax=b	  *	
! input parameters:
!
! - n is the number (INTEGER) of rows/columns in A
! - A is the COMPLEX(wp) array storing the coefficient matrix
!   A is defined as A(lda,lda)
! - b is the COMPLEX(wp) vector storing the right-hand side of the
!   linear system
! - x is the COMPLEX(wp) vector storing the initial guess in input.
!   On output, it stores the approximate solution.

! matvec     (input) EXTERNAL name of matrix vector subroutine
!            to deliver y:=A*x by CALL matvec(n,x,y)
! 

! B. Carpentieri, Y.-F. Jing and T.-Z. Huang, 2011, The BiCOR and CORS Iterative Algorithms 
! for Solving Nonsymmetric Linear Systems, SIAM J. Sci. Comput. 33, 
! Special Section: 2010 Copper Mountain Conference, 3020–3036. (17 pages) 

        USE DDPRECISION,ONLY : WP
        use cgmodule
        implicit none      
        type (cgstruct) cg
!	Variable definition
	INTEGER n, alloc_error, dealloc_error,ioerr
!        COMPLEX(wp):: aa(*)
!	INTEGER ja(*),ia(*)
        COMPLEX(wp)::  b(n),x(n)
        INTEGER flag

      integer mxcxsc
      complex(wp),target:: cxsc(n,21)

        INTEGER i,lda,nnz,j,nmv        
        INTEGER iter,max_it,printing,prectype
        REAL (wp) :: tol,error
        COMPLEX(wp) :: rho,rho1,beta,alpha,dpr,dp1,dp2,omega        

        INTEGER ios        
        COMPLEX(wp), pointer, dimension(:) ::  p,q_hat,q,r_hat,r,y                
        COMPLEX(wp), pointer, dimension(:) ::  p_hat,s_hat,s                        
        COMPLEX(wp), pointer, dimension(:) ::  r_s,t,tmp,tmp1
	COMPLEX(wp), pointer, dimension(:) ::  f,z,q1,e,ze,d,zh,h                         
        COMPLEX(wp) zero, one
        external matvec
        INTRINSIC cmplx
        PARAMETER(ZERO=cmplx(0.0_wp,0.0_wp),ONE=cmplx(1.0_wp,0.0_wp))

        real, dimension(2) :: tarray
        REAL(wp) bnrm2, dnorm2

! allocate temporary vectors
      ioerr=cg%ioerr

      if(n*21 .gt. mxcxsc) stop 'stop in cg routine cors90ver2 n*21> mxcxsc not enough memory'
      p  => cxsc(:,1)
      q_hat  => cxsc(:,2)
      q  => cxsc(:,3)
      r_hat  => cxsc(:,4)
      r => cxsc(:,5)
      y   => cxsc(:,6)
      p_hat  => cxsc(:,7)
      s_hat  => cxsc(:,8)
      s  => cxsc(:,9)
      r_s  => cxsc(:,10)
      t => cxsc(:,11)
      tmp   => cxsc(:,12)
      tmp1 => cxsc(:,13)
      f  => cxsc(:,14)
      z  => cxsc(:,15)
      q1  => cxsc(:,16)
      e  => cxsc(:,17)
      ze => cxsc(:,18)
      d   => cxsc(:,19)
      zh => cxsc(:,20)
      h => cxsc(:,21)

!        allocate(p(n),q_hat(n),q(n),r_hat(n),r(n),y(n), p_hat(n),s_hat(n),s(n), &                        
!             r_s(n),t(n),tmp(n),tmp1(n), f(n),z(n),q1(n),e(n),ze(n),d(n),zh(n),h(n),STAT=alloc_error)
!        IF (alloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)') 'Allocation Error Detected in conjugate gradient cors'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' cors '
!        END IF

!       Initialization
!	iter = nr of iteration
        iter     = 0          
!	flag = information on success/failure of iterative solution                        
        flag     = 0
!	nmv = number of M-V products
	nmv      = 0
!	max_it = maximum number of iterations allowed
        max_it   = cg%maxit
!	tol = tolerance for iterative solution
	tol      = cg%epsilon_err
        rho      = ONE 
        bnrm2=sqrt(dot_product(b(1:n),b(1:n)))

       
        IF  (bnrm2.eq.ZERO) THEN
            bnrm2 = 1.0        
	    write(*,*) 'norm of right-hand side set to 1.0'       
        ENDIF
!	write(*,*) 'problem size = ',n
!	Initial residual
!        call camux (n, x, y, aa, ja, ia) 
!        call matvec(n, x, y) 
!        call prod(n,'u',x,'u',y)
        call matvec(x,y,n)

	nmv=nmv+1
	r(1:n) = b(1:n) - y(1:n)
        dnorm2 = sqrt(dot_product(r(1:n),r(1:n)))
        error = dnorm2/ bnrm2 

        IF ( error.le.tol ) go to 50
!	Initial shadow residual: default r_s = A*r
!        call matvec (n, r, r_s) 
!        call prod(n,'u',r,'u',r_s)
        call matvec(r,r_s,n)
 	nmv=nmv+1            
        dnorm2= sqrt(dot_product(r_s(1:n),r_s(1:n)))
        IF (dnorm2.eq.(0.0_wp)) r_s(1:n)=r(1:n) 

!        call matvec (n, r, r_s) 
!        call prod(n,'u',r,'u',r_s)
        call matvec(r,r_s,n)
 	nmv=nmv+1            
        dnorm2= sqrt(dot_product(r_s(1:n),r_s(1:n)))
        IF (dnorm2.eq.(0.0_wp)) r_s(1:n)=r(1:n)

!       Start iterations
        DO iter = 1,max_it
            z(1:n)=r(1:n)
!            call matvec (n, z, r_hat) 
!             call prod(n,'u',z,'u',r_hat)
            call matvec(z,r_hat,n)
	    nmv=nmv+1               
            rho1 = rho            
            rho  = dot_product(r_s(1:n),r_hat(1:n))     

            IF (rho.eq.ZERO) EXIT 
            IF (iter.eq.1) then 
                e(1:n)=r(1:n)
		ze(1:n)=e(1:n)
	        d(1:n) = r_hat(1:n)
	        q(1:n) = r_hat(1:n)
            ELSE
                beta=(rho/rho1)
	        e(1:n) = r(1:n) + beta*h(1:n)
	        ze(1:n) = z(1:n) + beta*zh(1:n)
		d(1:n) = r_hat(1:n) + beta*f(1:n)
	        q(1:n) = d(1:n) + beta*(f(1:n)+beta*q(1:n))
            ENDIF
!	    Precondition the direction q to obtain the new direction q1:
!	    i.e. solve  M q1 = q.
!
!	    By default no preconditioner is used, so that:	


		    q1(1:n)=q(1:n)

!           Matrix-vector product: q_hat = A*q1            
!            call camux (n, q1, q_hat, aa, ja, ia) 
!            call matvec(n, q1, q_hat) 
!            call prod(n,'u',q1,'u',q_hat)
            call matvec(q1, q_hat,n)

	    nmv=nmv+1                        
            dpr = dot_product(r_s(1:n),q_hat(1:n))            
            alpha=rho/dpr
	    h(1:n) = e(1:n) - alpha*q(1:n)
	    zh(1:n) = ze(1:n) - alpha*q1(1:n)
	    f(1:n) = d(1:n) - alpha*q_hat(1:n)
	    x(1:n) = x(1:n) + alpha*(2*ze(1:n)-alpha*q1(1:n))
	    r(1:n) = r(1:n) - alpha*(2*d(1:n)-alpha*q_hat(1:n))
!           Compute the residual reduction
            dnorm2= sqrt(dot_product(r(1:n),r(1:n)))
            error = dnorm2/ bnrm2                     
!	    Print at screen iteration number, residual reduction, elapsed CPU time	
             if(cg%print.eq.'print') WRITE(*,*) ' iter cors ', iter,error
!	    Check for convergence
            IF ( error.le.tol ) EXIT
        ENDDO
 
!	Set the flag for the iterative solution process
        IF (error.le.tol) THEN                  
!           we have converged        
            flag =  0
!	    print *,'Total nr of M-V products = ',nmv 	
        ELSEIF ( rho.eq.ZERO) THEN
	    print *,'Breakdown in the algorithm'
            flag = -3
        ELSE                                              
!           no convergence
            print*, ' no convergence '
            flag = -1
        ENDIF	
        cg%status=flag
        cg%itno=iter
 50     continue
!        deallocate(p,q_hat,q,r_hat,r,y,p_hat,s_hat,s, r_s,t,tmp,tmp1,f,z,q1,e,ze,d,zh,h, STAT=dealloc_error)
!        IF (dealloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)') 'Deallcation Error Detected in conjugate gradient cors'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' cors '
!        END IF
      return
      END subroutine cors90ver2

