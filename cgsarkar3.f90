!integration with DDSCAT
!changed subroutine names by adding 90 at the end

    SUBROUTINE gacg90ver2(cx,cy,n,cxsc,mxcxsc,matvec, cmatvec, cg)
      USE ddprecision, ONLY : wp
      use cgmodule
      implicit none
      type (cgstruct) cg

! Purpose: GACG --- generalized augumented conjugate gradient method
!          for unsymmetric case

! History: 1988     (TKS); original version
!          92/06/09 (PJF); re-written using "BLAS"-type style
!          92/09/11 (TLS); Removed unnecessary zeroing.
!     .. Parameters ..
!     .. Scalar Arguments ..
      REAL (wp) ::epsilon_err
      INTEGER :: iter, maxit, n, alloc_error, dealloc_error,ioerr
!     .. Array Arguments ..
      COMPLEX (wp) :: cx(n), cy(n)
      integer mxcxsc
      complex(wp),target:: cxsc(n,8)
!     .. Local Scalars ..
      REAL (wp) :: ak, ay, bk, bw2, q2, sk, sk2

!     .. Local Arrays ..
      COMPLEX (wp), pointer :: cbq(:), cbw(:), cq(:), cw(:)
!     .. External Functions ..

!     .. External Subroutines ..
      EXTERNAL matvec, cmatvec
!     .. Intrinsic Functions ..
      INTRINSIC conjg, sqrt
      maxit=cg%maxit
      epsilon_err=cg%epsilon_err
      ioerr=cg%ioerr

      if(n*8 .gt. mxcxsc)                                                     &
         stop 'stop in sarkar routine gacg90ver2 n*8> mxcxsc not enough memory'
      cbq  => cxsc(:,1)
      cbw  => cxsc(:,3)
      cq   => cxsc(:,5)
      cw   => cxsc(:,7)

!      allocate(cbq(2*n), cbw(2*n), cq(2*n), cw(2*n),STAT=alloc_error)
!        IF (alloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)')                                 &
!             'Allocation Error Detected in conjugate gradient gacg'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' gacg '
!        END IF

      iter = 0
      ay=2._wp*dot_product(cy(1:n),cy(1:n))
!      CALL prod(n,'u',cx,'u',cq)
      call matvec(cx,cq,n)
!      CALL prod(n,'c',cx,'u',cq(n+1))
      call cmatvec(cx,cq(n+1),n)
!
      cq(1:n) = cy(1:n) - cq(1:n)
      cq(n+1:2*n) = conjg(cy(1:n)) - cq(n+1:2*n)
      cw(1:n) = cq(1:n)
      cw(n+1:2*n) = cq(n+1:2*n)
!      CALL prod(n,'u',cq(n+1),'u',cbq)
      call matvec(cq(n+1),cbq,n)
      cbq(n+1:2*n)=cmplx(0.0_wp,0.0_wp,kind=wp)
!      CALL prod(n,'c',cq,'u',cbq(n+1))
      call matvec(cq,cbq(n+1),n)
!
      sk = 2.0_wp*sum(cbq(1:n)*conjg(cq(1:n)))
      cbw(1:2*n) = cbq(1:2*n)
      bw2=dot_product(cbw(1:2*n),cbw(1:2*n))
30    CONTINUE
      ak = sk/bw2
      cx(1:n)=cx(1:n)+ak*cw(n+1:2*n)
      cq(1:2*n)=cq(1:2*n)-ak*cbw(1:2*n)
      q2=dot_product(cq(1:2*n),cq(1:2*n))
      if(cg%print.eq.'print')WRITE (*,fmt=*) &
         'sqrt(q2/ay)= ', iter, sqrt(q2/ay)
      cg%itno=iter
      IF ((sqrt(q2/ay)).LT.epsilon_err) GO TO 50
!      CALL prod(n,'u',cq(n+1),'u',cbq)
      call matvec(cq(n+1),cbq,n)
!      CALL prod(n,'c',cq,'u',cbq(n+1))
      call cmatvec(cq,cbq(n+1),n)
!
      sk2 = 2.0_wp*sum(cbq(1:n)*conjg(cq(1:n)))
      bk = sk2/sk
! Warning there are differences in results (very small) when the code is run in 
! double precision two times. First time we use section 1, which is "fortran 77"
! construct with do loops. Second time I run section 2 which is "fortran 90" 
! type construct. I do not know what is the problem. This problems is not 
! occuring if the code is run in single precision. One would think that there
! should not be any difference.
! PJF
!section 1
!      DO i = 1, 2*n
!        cw(i) = cq(i) + bk*cw(i)
!        cbw(i) = cbq(i) + bk*cbw(i)
!      END DO
! section 2
       cw(1:2*n) = cq(1:2*n) + bk*cw(1:2*n)
       cbw(1:2*n) = cbq(1:2*n) + bk*cbw(1:2*n)
       bw2=dot_product(cbw(1:2*n),cbw(1:2*n))
      sk = sk2
      iter = iter + 1
      IF (iter.GT.maxit) GO TO 50
      GO TO 30
50    CONTINUE
!      deallocate(cbq, cbw, cq, cw,STAT=dealloc_error)
!        IF (dealloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)')                                  &
!             'Deallcation Error Detected in conjugate gradient gacg'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' gacg '
!        END IF
      RETURN
    END SUBROUTINE gacg90ver2





    SUBROUTINE gbicg90ver2(cx,cy,n,cxsc,mxcxsc,matvec, cmatvec, cg)
      USE ddprecision, ONLY : wp
      use cgmodule
      implicit none
      type (cgstruct) cg
! Purpose: GBICG --- generalized bi-conjugate gradient
! History: 1988     (TKS); original version
!          92/06/09 (PJF); re-written using "BLAS"-type style
!          92/09/11 (TLS); Removed unnecessary zeroing.
!     .. Parameters ..
!     .. Scalar Arguments ..
      REAL (wp) ::epsilon_err
      INTEGER :: iter, maxit,  n, alloc_error, dealloc_error,ioerr
!     .. Array Arguments ..
      COMPLEX (wp) ::  cx(n), cy(n)
!     .. Local Scalars ..
      COMPLEX (wp) :: cak, cbk, csk, csk2
      REAL (wp) :: ay, ek
      INTEGER :: i

      integer mxcxsc
      complex(wp),target:: cxsc(n,6)

!     .. Local Arrays ..
      COMPLEX (wp),pointer :: cap(:), caw(:), cp(:), cq(:), cr(:), cw(:)
!     .. External Functions ..
!     .. External Subroutines ..
      EXTERNAL matvec, cmatvec
!     .. Intrinsic Functions ..
      INTRINSIC conjg, sqrt
      maxit=cg%maxit
      epsilon_err=cg%epsilon_err
      ioerr=cg%ioerr
      if(n*6 .gt. mxcxsc)                                                      &
         stop 'stop in sarkar routine gbicg90ver2 n*6> mxcxsc not enough memory'
      cap  => cxsc(:,1)
      caw  => cxsc(:,2)
      cp   => cxsc(:,3)
      cq   => cxsc(:,4)
      cr   => cxsc(:,5)
      cw   => cxsc(:,6)
!      allocate(cap(n), caw(n), cp(n), cq(n), cr(n), cw(n),STAT=alloc_error)
!        IF (alloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)')
!             'Allocation Error Detected in conjugate gradient gbicg'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' gbicg '
!        END IF
      iter = 0._wp
      ay=dot_product(cy(1:n),cy(1:n))
!      CALL prod(n,'u',cx,'u',cr)
      call matvec(cx,cr,n)
      csk = (0._wp,0._wp)
      cr(1:n) = cy(1:n) - cr(1:n)
      cq(1:n) = conjg(cr(1:n))
      cw(1:n) = cq(1:n)
      cp(1:n) = cr(1:n)
      csk=  sum(cr(1:n)*cr(1:n))
20    CONTINUE
!      CALL prod(n,'c',cw,'u',caw)
      call cmatvec(cw,caw,n)
!      CALL prod(n,'u',cp,'u',cap)
      call matvec(cp,cap,n)
      cak = sum(cap(1:n)*conjg(cw(1:n)))
      cak = csk/cak
      cx(1:n)=cx(1:n)+cak*cp(1:n)
      cr(1:n)=cr(1:n)-cak*cap(1:n)
      cq(1:n)=cq(1:n)-conjg(cak)*caw(1:n)
      csk2 = sum(cr(1:n)*conjg(cq(1:n)))
      ek=dot_product(cr(1:n),cr(1:n))
      if(cg%print.eq.'print') WRITE (*,fmt=*) 'sqrt(ek/ay)= ', iter, sqrt(ek/ay)
      cg%itno=iter
      IF ((sqrt(ek/ay)).LT.epsilon_err) GO TO 40
      cbk = csk2/csk
      cp(1:n) = cr(1:n) + cbk*cp(1:n)
      cw(1:n) = cq(1:n) + conjg(cbk)*cw(1:n)
      csk = csk2
      iter = iter + 1
      IF (iter.GT.maxit) GO TO 40
      GO TO 20
40    CONTINUE
!      deallocate(cap, caw, cp, cq, cr, cw, STAT=dealloc_error)
!        IF (dealloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)')
!             'Deallcation Error Detected in conjugate gradient gbicg'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' gbicg '
!        END IF
      RETURN
    END SUBROUTINE gbicg90ver2

    SUBROUTINE gmcg90ver2(cx,cy,n,cxsc,mxcxsc,matvec, cmatvec, cg)
      USE ddprecision, ONLY : wp
      use cgmodule
      implicit none
      type (cgstruct) cg
! Purpose: GMCG --- modified conjugate method for unsymmetric complex
!          matrices

! History: 1988     (TKS); original version
!          92/06/09 (PJF); re-written using "BLAS"-type style
!          92/09/11 (TLS); Removed unnecessary zeroing.
!     .. Parameters ..
!     .. Scalar Arguments ..
      REAL (wp) ::epsilon_err
      INTEGER :: iter, maxit,  n, alloc_error, dealloc_error,ioerr
!     .. Array Arguments ..
      COMPLEX (wp) ::  cx(n), cy(n)
!     .. Local Scalars ..
      REAL (wp) :: ak, ay, bk, bk2, q2, w2
      integer mxcxsc
      complex(wp),target:: cxsc(n,8)
!     .. Local Arrays ..
      COMPLEX (wp), pointer :: cbw(:), cq(:), cw(:), cxh(:), cz(:)
!     .. External Functions ..
!     .. External Subroutines ..
      EXTERNAL matvec, cmatvec
!     .. Intrinsic Functions ..
      INTRINSIC conjg, sqrt
      maxit=cg%maxit
      epsilon_err=cg%epsilon_err
      ioerr=cg%ioerr

      if(n*8 .gt. mxcxsc)                                                     &
         stop 'stop in sarkar routine gmcg90ver2 n*8> mxcxsc not enough memory'
      cbw  => cxsc(:,1)
      cq  => cxsc(:,3)
      cw   => cxsc(:,5)
      cxh   => cxsc(:,7)
      cz   => cxsc(:,8)

!      allocate(cbw(2*n), cq(2*n), cw(2*n), cxh(n), cz(n),STAT=alloc_error)
!        IF (alloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)')                                 &
!             'Allocation Error Detected in conjugate gradient gmcg'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' gmcg '
!        END IF
      iter = 0
      cz(1:n) = cx(1:n)
      cxh(1:n) = cx(1:n)
      ay=2._wp*dot_product(cy(1:n),cy(1:n))
!      CALL prod(n,'u',cx,'u',cq)
      call matvec(cx,cq,n)
!      CALL prod(n,'c',cz,'u',cq(n+1))
      call cmatvec(cz,cq(n+1),n)
!
      cq(1:n) = cy(1:n) - cq(1:n)
      cq(n+1:2*n) = conjg(cy(1:n)) - cq(n+1:2*n)
      cw(1:n) = cq(1:n)
      cw(n+1:2*n) = cq(n+1:2*n)
      w2=dot_product(cw(1:2*n),cw(1:2*n))
30    CONTINUE
!      CALL prod(n,'u',cw(n+1),'u',cbw)
      call matvec(cw(n+1),cbw,n)
!      CALL prod(n,'c',cw,'u',cbw(n+1))
      call cmatvec(cw,cbw(n+1),n)
!
      ak = 2.0_wp*sum(cbw(1:n)*conjg(cw(1:n)))
      IF (ak.EQ.0._wp) THEN
        PRINT *, ' gmcg.  ak = 0 '
        STOP
      END IF
      ak = w2/ak
      cx(1:n)=cx(1:n)+ak*cw(n+1:2*n)
      cz(1:n)=cz(1:n)+ak*cw(1:n)
      cq(1:2*n)=cq(1:2*n)-ak*cbw(1:2*n)
      q2=dot_product(cq(1:2*n),cq(1:2*n))
      bk = q2/w2
      bk2 = 1._wp + bk
      cxh(1:n) = (cx(1:n)+bk*cxh(1:n))/bk2
      cw(1:2*n) = (cq(1:2*n)+bk*cw(1:2*n))/bk2
      w2=dot_product(cw(1:2*n),cw(1:2*n))
      if(cg%print.eq.'print')                             &
         WRITE (*,fmt=*) 'sqrt(w2/ay)= ', iter, sqrt(w2/ay)
      cg%itno=iter
      IF ((sqrt(w2/ay)).LT.epsilon_err) GO TO 60
      iter = iter + 1
      IF (iter.GT.maxit) GO TO 60
      GO TO 30
60    CONTINUE
        cx(1:n) = cxh(1:n)
!      deallocate(cbw, cq, cw, cxh, cz, STAT=dealloc_error)
!        IF (dealloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)') 'Deallcation Error Detected in conjugate gradient gmcg'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' gmcg '
!        END IF
      RETURN
    END SUBROUTINE gmcg90ver2





    SUBROUTINE rcg90ver2(cx,cy,n,cxsc,mxcxsc,matvec, cmatvec, cg)
      USE ddprecision, ONLY : wp
      use cgmodule
      implicit none
      type (cgstruct) cg
! Purpose: RCG --- residual minimized; the residuals are minimized
!          at each iteration, no scaling is introduced.

! History: 1988     (TKS); original version
!          92/06/09 (PJF); re-written using "BLAS"-type style
!          92/09/11 (TLS); Removed unnecessary zeroing.
! Warning there are differences in results (very small) when the code is run in
! double precision two times. First time we use section 1, which is "fortran 77"
! construct with do loops. Second time I run section 2 which is "fortran 90"
! type construct. I do not know what is the problem. This problems is not 
! occuring if the code is run in single precision. One would think that there 
! should not be any difference.
! PJF
!     .. Parameters ..
!     .. Scalar Arguments ..
      REAL (wp) ::epsilon_err
      INTEGER :: iter, maxit,  n, alloc_error, dealloc_error,ioerr
!     .. Array Arguments ..
      COMPLEX (wp) :: cx(n), cy(n)
!     .. Local Scalars ..
      REAL (wp) :: ak, ay, ek, qk, sk, sk2
      integer mxcxsc
      complex(wp),target:: cxsc(n,3)
!     .. Local Arrays ..
      COMPLEX (wp), pointer :: cp(:), cprod(:), cr(:)
!     .. External Functions ..
!     .. External Subroutines ..
      EXTERNAL matvec, cmatvec
!     .. Intrinsic Functions ..
      INTRINSIC sqrt
      maxit=cg%maxit
      epsilon_err=cg%epsilon_err
      ioerr=cg%ioerr
      if(n*3 .gt. mxcxsc)                                                    &
         stop 'stop in sarkar routine rcg90ver2 n*3> mxcxsc not enough memory'
      cp  => cxsc(:,1)
      cprod  => cxsc(:,2)
      cr   => cxsc(:,3)
!      allocate(cp(n), cprod(n), cr(n),STAT=alloc_error)
!        IF (alloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)')                                &
!             'Allocation Error Detected in conjugate gradient rcg'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' rcg '
!        END IF
      iter = 0
      ay=dot_product(cy(1:n),cy(1:n))
!      CALL prod(n,'u',cx,'u',cprod)
      call matvec(cx,cprod,n)
      cr(1:n) = cy(1:n) - cprod(1:n)
!      CALL prod(n,'c',cr,'u',cp)
      call cmatvec(cr,cp,n)
      sk=dot_product(cp(1:n),cp(1:n))
20    CONTINUE
!      CALL prod(n,'u',cp,'u',cprod)
      call matvec(cp,cprod,n)
      ak=dot_product(cprod(1:n),cprod(1:n))
      ak = sk/ak
      cx(1:n)=cx(1:n)+ak*cp(1:n)
      cr(1:n)=cr(1:n)-ak*cprod(1:n)
      ek=dot_product(cr(1:n),cr(1:n))
      if(cg%print.eq.'print')WRITE (*,fmt=*) 'sqrt(ek/ay)= ', iter, sqrt(ek/ay)
      cg%itno=iter
      IF ((sqrt(ek/ay)).LT.epsilon_err) GO TO 40
!      CALL prod(n,'c',cr,'u',cprod)
      call cmatvec(cr,cprod,n)
!
      sk2=dot_product(cprod(1:n),cprod(1:n))
      qk = sk2/sk
      cp(1:n) = cprod(1:n) + qk*cp(1:n)
      iter = iter + 1
      sk = sk2
      IF (iter.GT.maxit) GO TO 40
      GO TO 20
40    CONTINUE
!      deallocate(cp, cprod, cr, STAT=dealloc_error)
!        IF (dealloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)')                                &
!            'Deallcation Error Detected in conjugate gradient rcg'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop '  rcg '
!        END IF
      RETURN
    END SUBROUTINE rcg90ver2


    SUBROUTINE sacg90ver2(cx,cy,n,cxsc,mxcxsc,matvec, cg)
      USE ddprecision, ONLY : wp
      use cgmodule
      implicit none
      type (cgstruct) cg
! Purpose: SACG --- augumented conjugate gradient method for symmetric
!          case.
! History: 1988     (TKS); original version
!          92/06/09 (PJF); re-written using "BLAS"-type style
!          92/09/11 (TLS); Removed unnecessary zeroing.
! Warning there are differences in results (very small) when the code is run in
! double precision two times. First time we use section 1, which is "fortran 77"
! construct with do loops. Second time I run section 2 which is "fortran 90"
! type construct. I do not know what is the problem. This problems is not 
! occuring if the code is run in single precision. One would think that there
! should not be any difference.
! PJF
!     .. Parameters ..
!     .. Scalar Arguments ..
      REAL (wp) ::epsilon_err
      INTEGER :: iter, maxit,  n, alloc_error, dealloc_error,ioerr
!     .. Array Arguments ..
      COMPLEX (wp) ::  cx(n), cy(n)
      integer mxcxsc
      complex(wp),target:: cxsc(n,4)
!     .. Local Scalars ..
      REAL (wp) :: ak, ap2, ay, bk, r2, sk, sk2
!     .. Local Arrays ..
      COMPLEX (wp), pointer :: cap(:), car(:), cp(:), cr(:)
!     .. External Functions ..
!     .. External Subroutines ..
      EXTERNAL matvec
!     .. Intrinsic Functions ..
      INTRINSIC sqrt
      maxit=cg%maxit
      epsilon_err=cg%epsilon_err
      ioerr=cg%ioerr
      if(n*4 .gt. mxcxsc)                                                     &
         stop 'stop in sarkar routine sacg90ver2 n*4> mxcxsc not enough memory'
      cap  => cxsc(:,1)
      car  => cxsc(:,2)
      cp   => cxsc(:,3)
      cr   => cxsc(:,4)

!      allocate(cap(n), car(n), cp(n), cr(n),STAT=alloc_error)
!        IF (alloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)')                                 &
!             'Allocation Error Detected in conjugate gradient sacg'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' sacg '
!        END IF

      iter = 0
      ay=dot_product(cy(1:n),cy(1:n))
!      CALL prod(n,'u',cx,'u',cr)
      call matvec(cx,cr,n)
      cr(1:n) = cy(1:n) - cr(1:n)
      cp(1:n) = cr(1:n)
!      CALL prod(n,'u',cr,'c',car)
      call matvec(conjg(cr),car,n)
!
      sk = sum(car(1:n)*conjg(cr(1:n)))
      cap(1:n) = car(1:n)
      ap2=dot_product(cap(1:n),cap(1:n))
30    CONTINUE
      ak = sk/ap2
      cx(1:n)=cx(1:n)+ak*conjg(cp(1:n))
      cr(1:n)=cr(1:n)-ak*cap(1:n)
      r2=dot_product(cr(1:n),cr(1:n))
      car(1:n)=cmplx(0.0_wp,0.0_wp,kind=wp)
      if(cg%print.eq.'print')WRITE (*,fmt=*) 'sqrt(r2/ay)= ', iter, sqrt(r2/ay)
        cg%itno=iter
      IF ((sqrt(r2/ay)).LT.epsilon_err) GO TO 50
!      CALL prod(n,'u',cr,'c',car)
      call matvec(conjg(cr),car,n)
!
      sk2 = sum(car(1:n)*conjg(cr(1:n)))
      bk = sk2/sk
      cp(1:n) = cr(1:n) + bk*cp(1:n)
      cap(1:n) = car(1:n) + bk*cap(1:n)
      ap2=dot_product(cap(1:n),cap(1:n))
      sk = sk2
      iter = iter + 1
      IF (iter.GT.maxit) GO TO 50
      GO TO 30
50    CONTINUE
!      deallocate(cap, car, cp, cr, STAT=dealloc_error)
!        IF (dealloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)')                                  &
!             'Deallcation Error Detected in conjugate gradient sacg'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' sacg '
!        END IF
      RETURN
    END SUBROUTINE sacg90ver2





    SUBROUTINE sbicg90ver2(cx,cy,n,cxsc,mxcxsc,matvec, cg)
      USE ddprecision, ONLY : wp
      use cgmodule
      implicit none
      type (cgstruct) cg
! Purpose: SBICG --- simplified biconjugate gradient method specialized
!          to symmetric matrices

! History: 1988     (TKS); original version
!          92/06/09 (PJF); re-written using "BLAS"-type style
!          92/09/11 (TLS); Removed unnecessary zeroing.
!          13.01.10 (BTD); added CMSGNM variable and call to WRIMSG
!                          to report progress in DDSCAT standard way...
!     .. Parameters ..
!     .. Scalar Arguments ..
      CHARACTER :: CMSGNM*70
      REAL (wp) :: epsilon_err, time
      INTEGER :: iter, maxit, n, alloc_error, dealloc_error,ioerr
!     .. Array Arguments ..
      COMPLEX (wp) ::  cx(n), cy(n)
!     .. Local Scalars ..
      COMPLEX (wp) :: cak, cbk, csk, csk2
      REAL (wp) :: ay, ek
      integer mxcxsc
      complex(wp),target:: cxsc(n,3)
!     .. Local Arrays ..
      COMPLEX (wp), pointer :: cap(:), cp(:), cr(:)
!     .. External Functions ..

!     .. External Subroutines ..
      EXTERNAL  matvec
!     .. Intrinsic Functions ..
      INTRINSIC sqrt
      maxit=cg%maxit
      epsilon_err=cg%epsilon_err
      ioerr=cg%ioerr
      if(n*3 .gt. mxcxsc)                                                     &
         stop 'stop in sarkar routine sbicgver2 n*3 > mxcxsc not enough memory'
      cap  => cxsc(:,1)
      cp   => cxsc(:,2)
      cr   => cxsc(:,3)
!      allocate(cap(n), cp(n), cr(n),STAT=alloc_error)
!        IF (alloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)')                                  &
!             'Allocation Error Detected in conjugate gradient sbicg'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' sbicg '
!        END IF
      iter = 0
      ay=dot_product(cy(1:n),cy(1:n))
!      CALL prod(n,'u',cx,'u',cr)
      call matvec(cx,cr,n)
      cr(1:n) = cy(1:n) - cr(1:n)
      cp(1:n) = cr(1:n)
      csk = sum(cr(1:n)*cr(1:n))
20    CONTINUE
!      CALL prod(n,'u',cp,'u',cap)
      call matvec(cp,cap,n)
      cak = sum(cap(1:n)*cp(1:n))
      cak = csk/cak
      cx(1:n)=cx(1:n)+cak*cp(1:n)
      cr(1:n)=cr(1:n)-cak*cap(1:n)
      csk2 = sum(cr(1:n)*cr(1:n))
      ek=dot_product(cr(1:n),cr(1:n))
      IF(CG%PRINT.EQ.'print')THEN
         WRITE(CMSGNM,FMT='(A,I8,A,1P,E10.3)') &
               'IT=',ITER,' f.err=',SQRT(EK/AY)  
         CALL WRIMSG('SBICGM',CMSGNM)
!         WRITE (*,fmt=*) 'sqrt(ek/ay)= ', iter, sqrt(ek/ay)
      ENDIF
      cg%itno=iter
      IF ((sqrt(ek/ay)).LT.epsilon_err) GO TO 40
      cbk = csk2/csk
      cp(1:n) = cr(1:n) + cbk*cp(1:n)
      csk = csk2
      iter = iter + 1
      IF (iter.GT.maxit) GO TO 40
      GO TO 20
40    CONTINUE
!      deallocate(cap, cp, cr, STAT=dealloc_error)
!        IF (dealloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)') 'Deallcation Error Detected in conjugate gradient sbicg'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' sbicg '
!        END IF
      RETURN
    END SUBROUTINE sbicg90ver2


    SUBROUTINE smcg90ver2(cx,cy,n,cxsc,mxcxsc,matvec, cg)
      USE ddprecision, ONLY : wp
      use cgmodule
      implicit none
      type (cgstruct) cg
! Purpose: SMCG --- modified conjugate gradient method for symmetric
!          case.

! History: 1988     (TKS); original version
!          92/06/09 (PJF); re-written using "BLAS"-type style
!          92/09/11 (TLS); Removed unnecessary zeroing.
!     .. Parameters ..
!     .. Scalar Arguments ..
      REAL (wp) ::epsilon_err
      INTEGER :: iter, maxit, n, alloc_error, dealloc_error,ioerr
!     .. Array Arguments ..
      COMPLEX (wp) ::  cx(n), cy(n)
!     .. Local Scalars ..
      REAL (wp) :: ak, ay, bk, bk2, p2, r2
      integer mxcxsc
      complex(wp),target:: cxsc(n,4)
!     .. Local Arrays ..
      COMPLEX (wp), pointer :: cap(:), cp(:), cr(:), cxh(:)
!     .. External Functions ..
!     .. External Subroutines ..
      EXTERNAL matvec
!     .. Intrinsic Functions ..
      INTRINSIC sqrt
      maxit=cg%maxit
      epsilon_err=cg%epsilon_err
      ioerr=cg%ioerr
      if(n*4 .gt. mxcxsc) stop 'stop in sarkar routine smcg90ver2 n*4> mxcxsc not enough memory'
      cap  => cxsc(:,1)
      cp   => cxsc(:,2)
      cr   => cxsc(:,3)
      cxh  => cxsc(:,4)
!      allocate(cap(n), cp(n), cr(n), cxh(n),STAT=alloc_error)
!        IF (alloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)') 'Allocation Error Detected in conjugate gradient smcg'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' smcg '
!        END IF
      iter = 0
      cxh(1:n) = cx(1:n)
      ay=dot_product(cy(1:n),cy(1:n))
!      CALL prod(n,'u',cx,'u',cr)
      call matvec(cx,cr,n)
      cr(1:n) = cy(1:n) - cr(1:n)
      cp(1:n) = cr(1:n)
      p2=dot_product(cp(1:n),cp(1:n))
30    CONTINUE
!      CALL prod(n,'u',cp,'c',cap)
      call matvec(conjg(cp),cap,n)
!
      ak = sum(cap(1:n)*conjg(cp(1:n)))
      ak = p2/ak
      cx(1:n)=cx(1:n)+ak*conjg(cp(1:n))
      cr(1:n)=cr(1:n)-ak*cap(1:n)
      r2=dot_product(cr(1:n),cr(1:n))
      bk = r2/p2
      bk2 = 1._wp + bk
      cxh(1:n) = (cx(1:n)+bk*cxh(1:n))/bk2
      cp(1:n) = (cr(1:n)+bk*cp(1:n))/bk2
      p2=dot_product(cp(1:n),cp(1:n))
      if(cg%print.eq.'print')  WRITE (*,fmt=*) 'sqrt(p2/ay)= ', iter, sqrt(p2/ay)
      cg%itno=iter
      IF ((sqrt(p2/ay)).LT.epsilon_err) GO TO 60
      iter = iter + 1
      IF (iter.GT.maxit) GO TO 60
      GO TO 30
60    CONTINUE
      cx(1:n) = cxh(1:n)
!      deallocate(cap, cp, cr, cxh, STAT=dealloc_error)
!        IF (dealloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)') 'Deallcation Error Detected in conjugate gradient smcg'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' smcg '
!        END IF
      RETURN
    END SUBROUTINE smcg90ver2





    SUBROUTINE srcg90ver2(cx,cy,n,cxsc,mxcxsc,matvec, cg)
      USE ddprecision, ONLY : wp
      use cgmodule
      implicit none
      type (cgstruct) cg
! Purpose: SRCG --- search directions scaled at each iteration
!          and the residuals are minimized
! History: 1988     (TKS); original version
!          92/06/09 (PJF); re-written using "BLAS"-type style
!          92/09/11 (TLS); Removed unnecessary zeroing.
!     .. Parameters ..
!     .. Scalar Arguments ..
      REAL (wp) ::epsilon_err
      INTEGER :: iter, maxit, n, alloc_error, dealloc_error,ioerr
!     .. Array Arguments ..
      COMPLEX (wp) ::  cx(n), cy(n)
      integer mxcxsc
      complex(wp),target:: cxsc(n,3)
!     .. Local Scalars ..
      REAL (wp) :: ak, ay, ek, sk
!     .. Local Arrays ..
      COMPLEX (wp), pointer :: cp(:), cprod(:), cr(:)
!     .. External Functions ..
!     .. External Subroutines ..
      EXTERNAL matvec, cmatvec
!     .. Intrinsic Functions ..
      INTRINSIC sqrt
      maxit=cg%maxit
      epsilon_err=cg%epsilon_err
      ioerr=cg%ioerr
      if(n*3 .gt. mxcxsc) stop 'stop in sarkar routine srcg90ver2 n*3> mxcxsc not enough memory'
      cp  => cxsc(:,1)
      cprod  => cxsc(:,2)
      cr   => cxsc(:,3)
!      allocate(cp(n), cprod(n), cr(n),STAT=alloc_error)
!        IF (alloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)') 'Allocation Error Detected in conjugate gradient srcg'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' srcg '
!        END IF
      iter = 0
      ay=dot_product(cy(1:n),cy(1:n))
!      CALL prod(n,'u',cx,'u',cprod)
      call matvec(cx,cprod,n)
      cr(1:n) = cy(1:n) - cprod(1:n)
!      CALL prod(n,'c',cr,'u',cprod)
      call cmatvec(cr,cprod,n)
!
      sk=dot_product(cprod(1:n),cprod(1:n))
      cp(1:n) = cprod(1:n)/sk
30    CONTINUE
!      CALL prod(n,'u',cp,'u',cprod)
      call matvec(cp,cprod,n)
      ak=dot_product(cprod(1:n),cprod(1:n))
      ak = 1._wp/ak
      cx(1:n)=cx(1:n)+ak*cp(1:n)
      cr(1:n)=cr(1:n)-ak*cprod(1:n)
      ek=dot_product(cr(1:n),cr(1:n))
      if(cg%print.eq.'print')  WRITE (*,fmt=*) 'sqrt(ek/ay)= ', iter, sqrt(ek/ay)
        cg%itno=iter
      IF ((sqrt(ek/ay)).LT.epsilon_err) GO TO 40
!      CALL prod(n,'c',cr,'u',cprod)
      call cmatvec(cr,cprod,n)
!
      sk=dot_product(cprod(1:n),cprod(1:n))
      cp(1:n)=cp(1:n)+(1._wp/sk)*cprod(1:n)
      iter = iter + 1
      IF (iter.GT.maxit) GO TO 40
      GO TO 30
40    CONTINUE
!      deallocate(cp, cprod, cr, STAT=dealloc_error)
!        IF (dealloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)') 'Deallcation Error Detected in conjugate gradient srcg'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' srcg '
!        END IF
      RETURN
    END SUBROUTINE srcg90ver2





    SUBROUTINE xcg90ver2(cx,cy,n,cxsc,mxcxsc,matvec, cmatvec, cg)
      USE ddprecision, ONLY : wp
      use cgmodule
      implicit none
      type (cgstruct) cg

! Purpose: XCG --- the error between the true solution and the
!          approximate solution is minimized at each iteration

! History: 1988     (TKS); original version
!          92/06/09 (PJF); re-written using "BLAS"-type style
!          92/09/11 (TLS); Removed unnecessary zeroing.

! Warning there are differences in results (very small) when the code is run in double precision
! two times. First time we use section 1, which is "fortran 77" construct with do loops. Second time I run
! section 2 which is "fortran 90" type construct. I do not know what is the problem. This problems is not 
! occuring if the code is run in single precision. One would think that there should not be any difference.
! PJF

!     .. Parameters ..
!     .. Scalar Arguments ..
      REAL (wp) ::epsilon_err
      INTEGER :: iter, maxit,  n, alloc_error, dealloc_error,ioerr
!     .. Array Arguments ..
      COMPLEX (wp) :: cx(n), cy(n)
!     .. Local Scalars ..
      REAL (wp) :: ak, ay, qk, sk, sk2
      integer mxcxsc
      complex(wp),target:: cxsc(n,3)

!     .. Local Arrays ..
      COMPLEX (wp), pointer :: cp(:), cprod(:), cr(:)

!     .. External Functions ..
!     .. External Subroutines ..
      EXTERNAL matvec, cmatvec
!     .. Intrinsic Functions ..
      INTRINSIC sqrt
      maxit=cg%maxit
      epsilon_err=cg%epsilon_err
      ioerr=cg%ioerr
      if(n*3 .gt. mxcxsc) stop 'stop in sarkar routine xcg90ver2 n*3> mxcxsc not enough memory'
      cp  => cxsc(:,1)
      cprod  => cxsc(:,2)
      cr   => cxsc(:,3)
!      allocate(cp(n), cprod(n), cr(n),STAT=alloc_error)
!        IF (alloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)') 'Allocation Error Detected in conjugate gradient xcg'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' xcg '
!        END IF
      iter = 0
      ay=dot_product(cy(1:n),cy(1:n))
!      CALL prod(n,'u',cx,'u',cprod)
      call matvec(cx,cprod,n)
      cr(1:n) = cy(1:n) - cprod(1:n)
      sk=dot_product(cr(1:n),cr(1:n))
!      CALL prod(n,'c',cr,'u',cp)
      call cmatvec(cr,cp,n)
!
20    CONTINUE
!      CALL prod(n,'u',cp,'u',cprod)
      call matvec(cp,cprod,n)
      ak=dot_product(cp(1:n),cp(1:n))
      ak = sk/ak
      cx(1:n)=cx(1:n)+ak*cp(1:n)
      cr(1:n)=cr(1:n)-ak*cprod(1:n)
      sk2=dot_product(cr(1:n),cr(1:n))
      if(cg%print.eq.'print')  WRITE (*,fmt=*) 'sqrt(sk2/ay)= ', iter, sqrt(sk2/ay)
        cg%itno=iter
      IF ((sqrt(sk2/ay)).LT.epsilon_err) GO TO 40
!      CALL prod(n,'c',cr,'u',cprod)
      call cmatvec(cr,cprod,n)
      qk = sk2/sk
      cp(1:n) = cprod(1:n) + qk*cp(1:n)
      iter = iter + 1
      sk = sk2
      IF (iter.GT.maxit) GO TO 40
      GO TO 20
40    CONTINUE
!      deallocate(cp, cprod, cr, STAT=dealloc_error)
!        IF (dealloc_error/=0) THEN
!          WRITE (ioerr,fmt='(a)') 'Deallcation Error Detected in conjugate gradient xcg'
!          WRITE (ioerr,fmt='(''Aborting'')')
!          stop ' xcg '
!        END IF
      RETURN
    END SUBROUTINE xcg90ver2
