    SUBROUTINE PRECOND
      RETURN
    END SUBROUTINE PRECOND

    SUBROUTINE ZBCG2(PRINT_RESID,L,N,X,NONZERO_X,RHS,MATVEC,PRECOND,TOLER, &
                     MXMATVEC,WORK,INFO)
      USE DDPRECISION,ONLY : WP
      USE DDCOMMON_9,ONLY: ITERMX,ITERN

! subroutine zbcg2 (print_resid,l,n,x,nonzero_x,rhs,matvec,precond,toler, &
!                   mxmatvec,work,info)

! Improved "vanilla" BiCGStab(2) iterative method

! Copyright (c) 2001 by M.A.Botchev (http://www.math.utwente.nl/~botchev/),
!                       University of Twente
! Permission to copy all or part of this work is granted,
! provided that the copies are not made or distributed
! for resale, and that the copyright notice and this
! notice are retained.

! This is the "vanilla" version of BiCGstab(\ell) as described
! in PhD thesis of D.R.Fokkema, Chapter 3 (also available as
! Preprint 976, Dept. of Mathematics, Utrecht University, URL
! http://www.math.uu.nl/publications/).  It includes two enhancements 
! to BiCGstab(\ell) proposed by G.Sleijpen and H.van der Vorst in
! 1) G.Sleijpen and H.van der Vorst "Maintaining convergence 
!    properties of BiCGstab methods in finite precision arithmetic",
!    Numerical Algorithms, 10, 1995, pp.203-223
! 2) G.Sleijpen and H.van der Vorst "Reliable updated residuals in
!    hybrid BiCG methods", Computing, 56, 1996, pp.141-163

! {{ This code based on original work of D.R.Fokkema:

! subroutine zbistbl v1.1 1998    
! Copyright (c) 1995-1998 by D.R. Fokkema.
! Permission to copy all or part of this work is granted, 
! provided that the copies are not made or distributed 
! for resale, and that the copyright notice and this 
! notice are retained.

! }}

! Your bug reports, comments, etc. are welcome: 
! m.a.botchev@math.utwente.nl

! ------------------------------
! Description of the parameters:
! ------------------------------

! print_resid (input) LOGICAL. If print_resid=.true. the number of 
!            matrix-vector multiplications done so far and residual norm will
!            be printed to the standard output each iteration

! l          (input) INTEGER the dimension \ell of BiCGstab(\ell)
!            in this simple version it is required that l <= 2
!            l=2 is often useful for systems with nonsymmetric matrices

! n          (input) INTEGER size of the linear system to solve 

! x          (input/output) COMPLEX*16 array dimension n
!            initial guess on input, solution on output

! rhs        (input) COMPLEX*16 array dimension n
!            the right-hand side (r.h.s.) vector

! matvec     (input) EXTERNAL name of matrix vector subroutine
!            to deliver y:=A*x by CALL matvec(n,x,y)

! nonzero_x  (input) LOGICAL tells
!            BiCGstab(\ell) if the initial guess x is zero or not. 
!            If nonzero_x is .FALSE., initial residual equals r.h.s. vector
!            and one MATVEC call is avoided

! toler      (input/output) DOUBLE PRECISION tolerance: the iterations are 
!            stopped as soon as || residual ||/|| initial residual|| <= toler,
!            the norm is Euclidean.  On output, if info>=0, the value of 
!            toler is set to the actually achieved residual reduction

! mxmatvec   (input/output) INTEGER.  On input: maximum number of matrix 
!            vector multiplications allowed to be done.  On output: 
!            if info>=0, mxmatvec is set to the actual number of matrix 
!            vector multiplications done

! work       (workspace) COMPLEX*16 array of dimension (n,2*l+5)

! info       (output) INTEGER.  info = 0 in case of succesful computations
!            and 
!            info = -m (<0) - means paramater number m has an illegal value
!            info = 1 - means no convergence achieved (stopping criterion
!            is not fulfilled)
!            info = 2 - means breakdown of the algorithm (taking a larger
!            value of parameter l usually helps)

! WARNING: If the iterations are ended normally (info=0 or info=1),
! the true residual norm is computed and returned as an output value 
! of the parameter toler.  The true residual norm can be slightly larger
! than the projected residual norm used by the algorithm to stop the
! iterations.  It may thus happen that on output info=0 but the value
! of toler is (slightly) larger than tolerance prescribed on input.
! ----------------------------------------------------------
! history:
! 08.03.14 (BTD) added declaration of dummy array IPAR(*)
!                changed
!                CALL MATVEC(X,WORK(1:N,R),N) -> 
!                CALL MATVEC(X,WORK(1:N,R),IPAR)
! 08.05.12 (BTD) Following suggestion by Art Lazanoff, NASA Ames:
!                CALL MATVEC(X,WORK(1:N,R),IPAR)
!                -> CALL MATVEC(X,WORK(1,R),IPAR) ->
!                RNRM0=DNORM2_BCG(N,WORK(1:N,R))
!                -> RNRM0=DNORM2_BCG(N,WORK(1,R))
!                RHO1=ZDOT_BCG(N,WORK(1:N,RR),WORK(1:N,R+K-1))
!                -> RHO1=ZDOT_BCG(N,WORK(1,RR),WORK(1,R+K-1))
!                CALL MATVEC(WORK(1:N,U+K-1),WORK(1:N,U+K),IPAR)
!                -> CALL MATVEC(WORK(1,U+K-1),WORK(1,U+K),IPAR)
!                SIGMA=ZDOT_BCG(N,WORK(1:N,RR),WORK(1:N,U+K))
!                -> SIGMA=ZDOT_BCG(N,WORK(1,RR),WORK(1,U+K))
!                CALL MATVEC(WORK(1:N,R+K-1),WORK(1:N,R+K),IPAR)
!                -> CALL MATVEC(WORK(1,R+K-1),WORK(1,R+K),IPAR)
!                MATRIX_Z(I,J)=
!                            CONJG(ZDOT_BCG(N,WORK(1:N,R+J-1),WORK(1:N,R+I-1)))
!                -> MATRIX_Z(I,J)=
!                            CONJG(ZDOT_BCG(N,WORK(1,R+J-1),WORK(1,R+I-1)))
!                CALL MATVEC(X,WORK(1:N,R),IPAR)
!                -> CALL MATVEC(X,WORK(1,R),IPAR)
!                CALL MATVEC(X,WORK(1:N,R),IPAR)
!                -> CALL MATVEC(X,WORK(1,R),IPAR)
!                RNRM=DNORM2_BCG(N,WORK(1:N,R))
!                -> RNRM=DNORM2_BCG(N,WORK(1,R))
! 08.07.16 (BTD) Added DDCOMMON_9 to communicate ITERMX and ITERN
! end history
      IMPLICIT NONE

! Parameters:

      LOGICAL,INTENT(IN) :: PRINT_RESID,NONZERO_X
      INTEGER,INTENT(IN) :: L,N
      INTEGER,INTENT(INOUT) :: MXMATVEC
      INTEGER,INTENT(OUT) :: INFO
      COMPLEX(WP),INTENT(INOUT) :: X(N)
      COMPLEX(WP),INTENT(IN) :: RHS(N)
      REAL(WP),INTENT(INOUT) :: TOLER
      COMPLEX(WP),INTENT(OUT) :: WORK(N,3+2*(L+1))
      EXTERNAL MATVEC, PRECOND

! Local variables:

      CHARACTER :: CMSGNM*70
      COMPLEX(WP) :: MATRIX_Z(L+1,L+1),Y0(L+1),YL(L+1),ZY0(L+1),ZYL(L+1)
      LOGICAL :: RCMP,XPDT
      INTEGER :: I,J,K,NMATVEC
      COMPLEX(WP) :: ALPHA,BETA,OMEGA,RHO0,RHO1,SIGMA
      COMPLEX(WP) :: VARRHO,HATGAMMA
      REAL(WP) :: RNRM0,RNRM
      REAL(WP) :: MXNRMX,MXNRMR
      COMPLEX(WP) :: KAPPA0,KAPPAL

! Define dummy array IPAR

      INTEGER IPAR(13)

! Aliases for the parts of the work array:

      INTEGER :: RR,R,U,XP,BP

! Constants:

      REAL(WP),PARAMETER :: DELTA=1D-2
      COMPLEX(WP),PARAMETER :: ZZERO=(0D0,0D0),ZONE=(1D0,0D0)

! Functions:

      REAL(WP) :: DNORM2_BCG
      COMPLEX(WP) :: ZDOT_BCG

!---------------------------------------------------------------------------
!*** diagnostic
!      write(0,*)'zbcg2wp ckpt 1'
!***
      INFO=0

      IF(L<1 .OR. L>2)INFO=-2
      IF(N<1)INFO=-3
      IF(TOLER<=0E0_WP)INFO=-9
      IF(MXMATVEC<0)INFO=-10

      RR=1
      R=RR+1
      U=R+(L+1)
      XP=U+(L+1)
      BP=XP+1

      IF(INFO/=0)RETURN

! Initialize first residual

      IF(NONZERO_X)THEN
!*** diagnostic
!         write(0,*)'zbcg2wp ckpt 2'
!****
!        CALL MATVEC(X,WORK(1:N,R),N)
!         CALL MATVEC(X,WORK(1:N,R),IPAR)
         CALL MATVEC(X,WORK(1,R),IPAR)
!*** diagnostic
!         write(0,*)'zbcg2wp ckpt 3'
!****

         WORK(1:N,R)=RHS-WORK(1:N,R)
         NMATVEC=1
      ELSE
         WORK(1:N,R)=RHS
         NMATVEC=0
      ENDIF

!call precond (n,work(1:n,r))

! Initialize iteration loop

      WORK(1:N,RR)=WORK(1:N,R)
      WORK(1:N,BP)=WORK(1:N,R)
      WORK(1:N,XP)=X
      X=ZZERO

!      RNRM0=DNORM2_BCG(N,WORK(1:N,R))
      RNRM0=DNORM2_BCG(N,WORK(1,R))
      RNRM=RNRM0

      MXNRMX=RNRM0
      MXNRMR=RNRM0
      RCMP=.FALSE.
      XPDT=.FALSE.

      ALPHA=ZZERO
      OMEGA=ZONE
      SIGMA=ZONE
      RHO0=ZONE

!BTD 080716:
      ITERN=0
!-----------

! diagnostic
!      write(0,*)'zbcg2wp ckpt 4: itermx=',itermx
!----------

! Iterate

      DO WHILE (RNRM>TOLER*RNRM0 .AND. NMATVEC<MXMATVEC)

!BTD 080716:
         ITERN=ITERN+1
         IF(ITERN>ITERMX)CALL ERRMSG('FATAL','ZBCG2WP',' ITERN>ITERMX')
!-----------
         
! =====================
! The BiCG part ---
! =====================

         RHO0=-OMEGA*RHO0
         DO K=1,L
!            RHO1=ZDOT_BCG(N,WORK(1:N,RR),WORK(1:N,R+K-1))
            RHO1=ZDOT_BCG(N,WORK(1,RR),WORK(1,R+K-1))
            IF(RHO0==ZZERO)THEN
               INFO=2
               TOLER=RNRM/RNRM0
               MXMATVEC=NMATVEC
               RETURN
            ENDIF
            BETA=ALPHA*(RHO1/RHO0)
            RHO0=RHO1
            DO J=0,K-1
               WORK(1:N,U+J)=WORK(1:N,R+J)-BETA*WORK(1:N,U+J)
            ENDDO
!*** diagnostic
!            write(0,*)'zbcq2wp ckpt 5'
!***
!          CALL MATVEC(WORK(1:N,U+K-1),WORK(1:N,U+K),N)
!           CALL MATVEC(WORK(1:N,U+K-1),WORK(1:N,U+K),IPAR)
           CALL MATVEC(WORK(1,U+K-1),WORK(1,U+K),IPAR)
!*** diagnostic
!           write(0,*)'zbcq2wp ckpt 6'
!***
!          call precond (n, work(1:n,u+k))
           NMATVEC=NMATVEC+1

!           SIGMA=ZDOT_BCG(N,WORK(1:N,RR),WORK(1:N,U+K))
           SIGMA=ZDOT_BCG(N,WORK(1,RR),WORK(1,U+K))
           IF(SIGMA==ZZERO)THEN
              INFO=2
              TOLER=RNRM/RNRM0
              MXMATVEC=NMATVEC
              RETURN
           ENDIF
           ALPHA=RHO1/SIGMA
           X(1:N)=ALPHA*WORK(1:N,U)+X(1:N)
           DO J=0,K-1
              WORK(1:N,R+J)=-ALPHA*WORK(1:N,U+J+1)+WORK(1:N,R+J)
           ENDDO
!           CALL MATVEC(WORK(1:N,R+K-1),WORK(1:N,R+K),IPAR)
           CALL MATVEC(WORK(1,R+K-1),WORK(1,R+K),IPAR)
!      call precond (n, work(1:n,r+k))
           NMATVEC=NMATVEC+1
           RNRM=DNORM2_BCG(N,WORK(1,R))
           MXNRMX=MAX(MXNRMX,RNRM)
           MXNRMR=MAX(MXNRMR,RNRM)
        ENDDO

! ==================================
! The convex polynomial part ---
! ==================================

!  --- Z = R'R

        DO I=1,L+1
           DO J=1,I
!              MATRIX_Z(I,J)=CONJG(ZDOT_BCG(N,WORK(1:N,R+J-1),WORK(1:N,R+I-1)))
              MATRIX_Z(I,J)=CONJG(ZDOT_BCG(N,WORK(1,R+J-1),WORK(1,R+I-1)))
           ENDDO
        ENDDO

!  lower triangular part of Z is computed; compute the rest knowing that Z^H=Z
        DO J=2,L+1
           MATRIX_Z(1:J-1,J)=CONJG(MATRIX_Z(J,1:J-1))
        ENDDO

!  small vectors y0 and yl

        Y0(1)=-ZONE
        Y0(2)=(MATRIX_Z(2,1)/MATRIX_Z(2,2)) ! works only for l=2
        Y0(L+1)=ZZERO

        YL(1)=ZZERO
        YL(2)=(MATRIX_Z(2,3)/MATRIX_Z(2,2)) ! works only for l=2
        YL(L+1)=-ZONE

!  --- Convex combination

! compute Z*y0 and Z*yl
        ZY0=ZZERO
        ZYL=ZZERO
        DO J=1,L+1
           ZY0=ZY0+MATRIX_Z(:,J)*Y0(J)
           ZYL=ZYL+MATRIX_Z(:,J)*YL(J)
        ENDDO

        KAPPA0=SQRT(ABS(ZDOT_BCG(L+1,Y0,ZY0)))
        KAPPAL=SQRT(ABS(ZDOT_BCG(L+1,YL,ZYL)))

        VARRHO=ZDOT_BCG(L+1,YL,ZY0)/(KAPPA0*KAPPAL)

        HATGAMMA=VARRHO/ABS(VARRHO)*MAX(ABS(VARRHO),7E-1_WP)

        Y0=Y0-(HATGAMMA*KAPPA0/KAPPAL)*YL


!  --- Update

        OMEGA=Y0(L+1)

        DO J=1,L
           WORK(1:N,U)=WORK(1:N,U)- Y0(J+1)*WORK(1:N,U+J)
           X(1:N)=X(1:N)+Y0(J+1)*WORK(1:N,R+J-1)
           WORK(1:N,R)=WORK(1:N,R)- Y0(J+1)*WORK(1:N,R+J)
        ENDDO

! y0 has changed; compute Z*y0 once more

        ZY0=ZZERO
        DO J=1,L+1
           ZY0=ZY0+MATRIX_Z(:,J)*Y0(J)
        ENDDO

        RNRM=SQRT(ABS(ZDOT_BCG(L+1,Y0,ZY0)))

! ================================
! The reliable update part ---
! ================================

        MXNRMX=MAX(MXNRMX,RNRM)
        MXNRMR=MAX(MXNRMR,RNRM)
        XPDT=(RNRM<DELTA*RNRM0 .AND. RNRM0<MXNRMX)
        RCMP=((RNRM<DELTA*MXNRMR .AND. RNRM0<MXNRMR).OR. XPDT)
        IF(RCMP)THEN
!           CALL MATVEC(X,WORK(1:N,R),N)
!           CALL MATVEC(X,WORK(1:N,R),IPAR)
           CALL MATVEC(X,WORK(1,R),IPAR)
!   call precond (n, work(1:n,r))
           NMATVEC=NMATVEC+1
           WORK(1:N,R)=WORK(1:N,BP)- WORK(1:N,R)
           MXNRMR=RNRM
           IF(XPDT)THEN

              WORK(1:N,XP)=X(1:N)+WORK(1:N,XP)
              X=ZZERO
              WORK(1:N,BP)=WORK(1:N,R)

              MXNRMX=RNRM
           ENDIF
        ENDIF

!        IF(print_resid)PRINT *,nmatvec,' ',rnrm
        
        IF(PRINT_RESID)THEN
           WRITE(CMSGNM,FMT='(A,I8,A,1P,E10.3)') &
                 'IT=',NMATVEC/4,' f.err=',RNRM/RNRM0
           CALL WRIMSG('ZBCG2 ',CMSGNM)
         ENDIF

      ENDDO

! =========================
! End of iterations ---
! =========================

      X(1:N)=X(1:N)+WORK(1:N,XP)

      IF(RNRM>TOLER*RNRM0)INFO=1

! compute the true residual:

! --------------------- One matvec can be saved by commenting out this:
!      CALL MATVEC(X,WORK(1:N,R),N)
!      CALL MATVEC(X,WORK(1:N,R),IPAR)
      CALL MATVEC(X,WORK(1,R),IPAR)

      WORK(1:N,R)=RHS(1:N)- WORK(1:N,R)
!call precond (n,work(1:n,r))
!      RNRM=DNORM2_BCG(N,WORK(1:N,R))
      RNRM=DNORM2_BCG(N,WORK(1,R))
      NMATVEC=NMATVEC+1
! --------------------- One matvec can be saved by commenting out this^

      TOLER=RNRM/RNRM0
      MXMATVEC=NMATVEC

    END SUBROUTINE ZBCG2


    FUNCTION ZDOT_BCG(N,ZX,ZY)
      USE DDPRECISION,ONLY : WP

! complex inner product function

      IMPLICIT NONE
      INTEGER,INTENT (IN) :: N
      COMPLEX(WP),INTENT (IN) :: ZX(N),ZY(N)
      COMPLEX(WP) ZDOT_BCG
      ZDOT_BCG=SUM(CONJG(ZX)*ZY)
    END FUNCTION ZDOT_BCG


    FUNCTION DNORM2_BCG(N,ZX)
      USE DDPRECISION,ONLY : WP

! L2 norm function

      IMPLICIT NONE
      INTEGER,INTENT (IN) :: N
      COMPLEX(WP),INTENT (IN) :: ZX(N)
      COMPLEX(WP),EXTERNAL :: ZDOT_BCG
      REAL(WP) :: DNORM2_BCG
      DNORM2_BCG=SQRT(ABS(ZDOT_BCG(N,ZX,ZX)))

    END FUNCTION DNORM2_BCG
