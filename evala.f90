    SUBROUTINE EVALA(CXADIA,CXAOFF,CXALPH,CXALOF,MXN3,NAT)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE
!*** Arguments:
      INTEGER :: MXN3, NAT
      COMPLEX (WP) :: CXADIA(MXN3), CXALOF(MXN3), CXALPH(MXN3), CXAOFF(MXN3)

!*** Local variables:
      INTEGER :: J1, J2, J3
      COMPLEX (WP) :: CX1, CX2, CX3, DENOM

!*********************************************************************
! Subroutine EVALA
! Purpose: to evaluate 3x3 ***diagonal*** blocks of matrix A for use in
! discrete-dipole scattering calculation.

! Given:   CXALPH(J,1-3)=(alpha_11,alpha_22,alpha_33) where alpha=
!                         complex polarizability tensor
!                         of dipole J (J=1-NAT, first 3*NAT elements
!                         of CXALPH used)
!          CXALOF(J,1-3)=(alpha_23,alpha_31,alpha_12) for dipole J
!                        (J=1-NAT, first 3*NAT elements of CXALOF used)
!          NAT=number of dipoles
!          MXN3,MXNAT=dimensioning information

! Returns: CXADIAG(1-3*NAT)=diagonal elements of matrix A
!          CXAOFF(J,1-3)=(a_23,a_31,a_12) where a_ij=3x3 matrix
!                        equal to inverse of 3x3 matrix alpha_ij
!                        for dipole J (first 3*NAT elements of CXAOFF
!                        are used)

! Note that it is assumed that vectors E and P will have data
! in order E_1x,E_2x,...,E_Nx,E_1y,E_2y,...,E_Ny,E_1z,E_2z,...,E_Nz

! B. T. Draine, Princeton Univ. Obs., 87.01.07
! History:
! 88.05.05 (BTD): Modifications...
! 90.11.06 (BTD): Modified to pass array dimensions
! 97.12.25 (BTD): Added CXAOFF,CXALOF to argument list
!                 Modified to allow nondiagonal polarizabilities
!                 alpha_ij
! 98.01.01 (BTD): Correct inconsistency in assumed data ordering.

! End history
! Copyright (C) 1993,1997,1998 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
!***********************************************************************
! Given symmetric matrix m with six independent elements a1, a2, a3, b1,
! b2, b3 (see below), we seek inverse matrix M with six independent
! elements A1, A2, A3, B1, B2, B3 such that

!      (a1 b3 b2)   (A1 B3 B2)   (1 0 0)
!      (b3 a2 b1) x (B3 A2 B1) = (0 1 0)
!      (b2 b1 a3)   (B2 B1 A3)   (0 0 1)

! solution:
!                     a2*a3 - b1*b1
! A1 = ----------------------------------------------
!      a1*a2*a3+2*b1*b2*b3-a1*b1**2-a2*b2**2-a3*b3**2

!                     a3*a1 - b2*b2
! A2 = ----------------------------------------------
!      a1*a2*a3+2*b1*b2*b3-a1*b1**2-a2*b2**2-a3*b3**2

!                     a1*a2 - b3*b3
! A3 = ----------------------------------------------
!      a1*a2*a3+2*b1*b2*b3-a1*b1**2-a2*b2**2-a3*b3**2

!                     b2*b3 - a1*b1
! B1 = ----------------------------------------------
!      a1*a2*a3+2*b1*b2*b3-a1*b1**2-a2*b2**2-a3*b3**2

!                     b3*b1 - a2*b2
! B2 = ----------------------------------------------
!      a1*a2*a3+2*b1*b2*b3-a1*b1**2-a2*b2**2-a3*b3**2

!                     b1*b2 - a3*b3
! B3 = ----------------------------------------------
!      a1*a2*a3+2*b1*b2*b3-a1*b1**2-a2*b2**2-a3*b3**2

!***********************************************************************

!*** diagnostic
!      write(0,*)'evala ckpt 0'
!***

      DO J1=1,NAT
         J2=J1+NAT
         J3=J2+NAT
         CX1=CXALOF(J1)**2
         CX2=CXALOF(J2)**2
         CX3=CXALOF(J3)**2
         DENOM=CXALPH(J1)*CXALPH(J2)*CXALPH(J3)+          &
               2._WP*CXALOF(J1)*CXALOF(J2)*CXALOF(J3)-    &
               CXALPH(J1)*CX1-CXALPH(J2)*CX2-CXALPH(J3)*CX3
         CXADIA(J1)=(CXALPH(J2)*CXALPH(J3)-CX1)/DENOM
         CXADIA(J2)=(CXALPH(J3)*CXALPH(J1)-CX2)/DENOM
         CXADIA(J3)=(CXALPH(J1)*CXALPH(J2)-CX3)/DENOM
         CXAOFF(J1)=(CXALOF(J2)*CXALOF(J3)-CXALPH(J1)*CXALOF(J1))/DENOM
         CXAOFF(J2)=(CXALOF(J3)*CXALOF(J1)-CXALPH(J2)*CXALOF(J2))/DENOM
         CXAOFF(J3)=(CXALOF(J1)*CXALOF(J2)-CXALPH(J3)*CXALOF(J3))/DENOM
!*** diagnostic
!         if(j1==133.or.j1==134)then
!            write(0,*)'evala ckpt 1:'
!            write(0,*)'       j1=',j1
!            write(0,*)'       j2=',j2
!            write(0,*)'       j3=',j3
!            write(0,*)'       cx1=',cx1
!            write(0,*)'       cx2=',cx2
!            write(0,*)'       cx3=',cx3
!            write(0,*)'     denom=',denom
!            write(0,*)'cxadia(j1)=',cxadia(j1)
!            write(0,*)'cxadia(j2)=',cxadia(j2)
!            write(0,*)'cxadia(j3)=',cxadia(j3)
!            write(0,*)'cxaoff(j1)=',cxaoff(j1)
!            write(0,*)'cxaoff(j2)=',cxaoff(j2)
!            write(0,*)'cxaoff(j3)=',cxaoff(j3)
!         endif
!***
      ENDDO
!*** diagnostic
!      write(0,*)'evala ckpt 99 : about to return'
!***
      RETURN

    END SUBROUTINE EVALA
