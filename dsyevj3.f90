! ----------------------------------------------------------------------
! Numerical diagonalization of 3x3 matrcies
! Copyright (C) 2006  Joachim Kopp
! ----------------------------------------------------------------------
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.

! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.

! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-13
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
    SUBROUTINE DSYEVJ3(A,Q,W)
      USE DDPRECISION, ONLY : WP
! ----------------------------------------------------------------------
! Calculates the eigenvalues and normalized eigenvectors of a symmetric
! matrix A using the Jacobi algorithm.
! The upper triangular part of A is destroyed during the calculation,
! the diagonal elements are read but not destroyed, and the lower
! triangular elements are not referenced at all.
! ----------------------------------------------------------------------
! Parameters:
!   A: The symmetric input matrix
!   Q: Storage buffer for eigenvectors
!   W: Storage buffer for eigenvalues
! ----------------------------------------------------------------------
!     .. Arguments ..
      REAL(WP) :: a(3,3)
      REAL(WP) :: q(3,3)
      REAL(WP) :: w(3)

!     .. Parameters ..
      INTEGER :: n
      PARAMETER (n=3)

!     .. Local Variables ..
      REAL(WP) :: sd, so
      REAL(WP) :: s, c, t
      REAL(WP) :: g, h, z, theta
      REAL(WP) :: thresh
      INTEGER :: i, x, y, r

!     Initialize Q to the identitity matrix
!     --- This loop can be omitted if only the eigenvalues are desired -
      DO 10 x = 1, n
        q(x,x) = 1.0E0_wp
        DO 11 y = 1, x - 1
          q(x,y) = 0.0E0_wp
          q(y,x) = 0.0E0_wp
11      CONTINUE
10    CONTINUE

!     Initialize W to diag(A)
      DO 20 x = 1, n
        w(x) = a(x,x)
20    CONTINUE

!     Calculate SQR(tr(A))
      sd = 0.0E0_wp
      DO 30 x = 1, n
        sd = sd + abs(w(x))
30    CONTINUE
      sd = sd**2

!     Main iteration loop
      DO 40 i = 1, 50
!       Test for convergence
        so = 0.0E0_wp
        DO 50 x = 1, n
          DO 51 y = x + 1, n
            so = so + abs(a(x,y))
51        CONTINUE
50      CONTINUE
        IF (so==0.0E0_wp) THEN
!*** diagnostic
           write(0,*)'>DSYEVJ3: converged.'
!***
          RETURN
        END IF

        IF (i<4) THEN
          thresh = 0.2E0_wp*so/n**2
        ELSE
          thresh = 0.0E0_wp
        END IF

!       Do sweep
        DO 60 x = 1, n
          DO 61 y = x + 1, n
            g = 100.0E0_wp*(abs(a(x,y)))
            IF (i>4 .AND. abs(w(x))+g==abs(w(x)) .AND. abs(w(y))+g==abs(w(y))) &
                THEN
              a(x,y) = 0.0E0_wp
            ELSE IF (abs(a(x,y))>thresh) THEN
!             Calculate Jacobi transformation
              h = w(y) - w(x)
              IF (abs(h)+g==abs(h)) THEN
                t = a(x,y)/h
              ELSE
                theta = 0.5E0_wp*h/a(x,y)
                IF (theta<0.0E0_wp) THEN
                  t = -1.0E0_wp/(sqrt(1.0E0_wp+theta**2)-theta)
                ELSE
                  t = 1.0E0_wp/(sqrt(1.0E0_wp+theta**2)+theta)
                END IF
              END IF

              c = 1.0E0_wp/sqrt(1.0E0_wp+t**2)
              s = t*c
              z = t*a(x,y)

!             Apply Jacobi transformation
              a(x,y) = 0.0E0_wp
              w(x) = w(x) - z
              w(y) = w(y) + z
              DO 70 r = 1, x - 1
                t = a(r,x)
                a(r,x) = c*t - s*a(r,y)
                a(r,y) = s*t + c*a(r,y)
70            CONTINUE
              DO 80 r = x + 1, y - 1
                t = a(x,r)
                a(x,r) = c*t - s*a(r,y)
                a(r,y) = s*t + c*a(r,y)
80            CONTINUE
              DO 90 r = y + 1, n
                t = a(x,r)
                a(x,r) = c*t - s*a(y,r)
                a(y,r) = s*t + c*a(y,r)
90            CONTINUE

!             Update eigenvectors
!             --- This loop can be omitted if only the eigenvalues are desired
              DO 100 r = 1, n
                t = q(r,x)
                q(r,x) = c*t - s*q(r,y)
                q(r,y) = s*t + c*q(r,y)
100           CONTINUE
            END IF
61        CONTINUE
60      CONTINUE
40    CONTINUE

      WRITE(0,*)'>DSYEVJ3: No convergence.'
      RETURN

    END SUBROUTINE dsyevj3
