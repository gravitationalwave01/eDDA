    SUBROUTINE PRINAXIS(MXNAT,NAT,ICOMP,IXYZ,DX,A1,A2,EIGVAL)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

! Arguments:

      INTEGER :: MXNAT,NAT
      INTEGER*2 ::   &
        ICOMP(MXNAT,3)
      INTEGER ::    &
        IXYZ(MXNAT,3)
      REAL (WP) :: A1(3),A2(3),DX(3),EIGVAL(3)

! Local variables:

      INTEGER :: LWORK
      PARAMETER (LWORK=12)
      INTEGER :: J,J1,J2,J3,JA,JD,K
      REAL (WP) :: DENOM,DSUM,PI,ZNORM,ZNORM2
      REAL (WP) :: &
        A(3,3),    &
        A3(3),     &
        CM(3),     &
        Q(3,3),    &
        W(3)
      CHARACTER :: CMSGNM*70
!***********************************************************************
! Subroutine to determine principal axes of moment of inertia tensor
! of target.
! Target is assumed to be homogeneous.
! Given:

! MXNAT       =limit on largest allowed value of NAT=number of dipoles
! NAT         =number of dipoles in target
! IXYZ(J,1-3) =x,y,z location of dipole Jon lattice
! DX(1-3)     =(dx/d,dy/d,dz/d) where (dx,dy,dz)=lattice spacing in
!              x,y,z directions,
!              and d = (dx*dy*dz)**(1/3) = effective lattice spacing
! ICOMP(J,1-3)=composition for dipole J; x,y,z directions

! Returns:

! A1(1-3)     =principal axis A1 in target frame (normalized)
! A2(1-3)     =principal axis A2 in target frame (normalized)
! EIGV(1-3)   =(eigenvalues of moment of inertia tensor)/(0.4*M*aeff^2)
!              definition of aeff = (3*V/4*pi)^{1/3} , V = solid volume
!              EIGV are in decreasing order:
!              EIGV(1) is largest, EIGV(3) is smallest

!              A1 = principal axis with largest moment of inertia
!              A2 = principal axis with second-largest moment of inertia

! B.T. Draine, Princeton Univ. Observatory, 96.01.26
! History:
! 96.01.26 (BTD) extracted from routine TARBLOCKS
!                IX,IY,IZ replaced by IXYZ
! 96.02.23 (BTD) corrected typo in 2 format statements
!                add code computing axis a3 and printing out
!                eigenvalue.
! 96.11.01 (PJF) modified to use LAPACK code rather than EIGV
!                and associated routines from NSWC package.
! 96.11.21 (BTD) modified to remove WRITE(0 and replace with
!                calls to WRIMSG
! 97.12.26 (BTD) added DX(3) to argument list, and modified to
!                compute principal axes with unequal lattice
!                spacing.
! 98.04.27 (BTD) deleted unused variables I,ZI
! 98.12.15 (BTD) Appears to be a problem with LAPACK routine SGEEV
!                when compiled using g77 under Linux.
!                Add trap to catch this problem if it occurs,
!                print out warning, and terminate execution.
! 00.06.13 (BTD) corrected error in computation of diagonal elements
!                of moment of inertia tensor (NAT/6 -> NAT/18)
! 04.05.23 (BTD) added normalized eigenvalues EIGV to argument list
! 07.08.03 (BTD) modified to replace LAPACK routine SGEEV with
!                routine DSYEVJ3 from Joachim Kopp
!                arXiv.org:physics/0610206
!                use f90 version of DSYEVJ3
! 07.09.11 (BTD) changed IXYZ from INTEGER*2 to INTEGER
! End history.

! Copyright (C) 1993,1994,1995,1996,1997,1998,2000,2004,2007
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
      PI=4._WP*ATAN(1._WP)

! Compute moment of inertia tensor.

      DO JD=1,3
        CM(JD)=0._WP
      ENDDO
      DO JA=1,NAT
        DO JD=1,3
          CM(JD)=CM(JD)+IXYZ(JA,JD)*DX(JD)
        ENDDO
      ENDDO
      DO JD=1,3
        CM(JD)=CM(JD)/NAT
      ENDDO

! Compute moment of inertia tensor

      DO J=1,3
        DO K=1,3
          A(J,K)=0._WP
        ENDDO
      ENDDO
      DSUM=0._WP
      DO K=1,3
        DO JA=1,NAT
          DSUM=DSUM+((IXYZ(JA,K)-CM(K))*DX(K))**2
        ENDDO
      ENDDO
      DO J=1,3
        A(J,J)=REAL(DSUM+(NAT/18._WP)*(DX(1)**2+DX(2)**2+DX(3)**2),KIND=WP)
      ENDDO
      DO J=1,3
        DO K=1,3
          DSUM=0._WP
          DO JA=1,NAT
            DSUM=DSUM+(IXYZ(JA,J)-CM(J))*(IXYZ(JA,K)-CM(K))*DX(J)*DX(K)
          ENDDO
          A(J,K)=A(J,K)-REAL(DSUM,KIND=WP)
        ENDDO
      ENDDO

! Normalize moment of inertia tensor to units of 0.4*mass*a_eff**2

      DENOM=0.4_WP*NAT*(.75_WP*NAT/PI)**(2._WP/3._WP)
      DO K=1,3
        DO J=1,3
          A(J,K)=A(J,K)/DENOM
        ENDDO
      ENDDO

! Compute eigenvalues and eigenvectors of A(J,K)

      CALL DSYEVJ3(A,Q,W)

! on return from DSYEVJ3:
!    W(i) = eigenvalues i=1,2,3
!    Q(i,j=1,3) = eigenvector corresponding to eigenvalue W(i) 

! Normalize the eigenvectors:

      DO J=1,3
        ZNORM2=0._WP
        DO K=1,3
          ZNORM2=ZNORM2+Q(K,J)**2
        ENDDO

! Make first component of each eigenvector positive

        ZNORM=SQRT(ZNORM2)
        IF(Q(1,J)<0._WP)THEN
          ZNORM=-ZNORM
        ELSEIF(Q(1,J)==0._WP)THEN
          IF(Q(2,J)<0._WP)ZNORM=-ZNORM
        ENDIF
        DO K=1,3
          Q(K,J)=Q(K,J)/ZNORM
        ENDDO
      ENDDO

! Order the eigenvectors by eigenvalue, largest first:

      J1=1
      IF(W(2)>W(1))J1=2
      IF(W(3)>W(J1))J1=3
      IF(J1==1)THEN
        J2=2
        IF(W(3)>W(2))J2=3
      ELSEIF(J1==2)THEN
        J2=1
        IF(W(3)>W(1))J2=3
      ELSEIF(J1==3)THEN
        J2=1
        IF(W(2)>W(1))J2=2
      ENDIF
      DO K=1,3
        A1(K)=Q(K,J1)
        A2(K)=Q(K,J2)
      ENDDO

! Set J3=3 if J1+J2=1+2=3
!        2         =1+3=4
!        1         =2+3=5

      J3=6-J1-J2

! Having fixed eigenvectors a_1 and a_2, compute a_3 = a_1 x a_2

      A3(1)=A1(2)*A2(3)-A1(3)*A2(2)
      A3(2)=A1(3)*A2(1)-A1(1)*A2(3)
      A3(3)=A1(1)*A2(2)-A1(2)*A2(1)

      EIGVAL(1)=W(J1)
      EIGVAL(2)=W(J2)
      EIGVAL(3)=W(J3)

!*** enable following lines to print out eigenvalues

      WRITE (CMSGNM,FMT='(A,0P,F8.5)')'I_1/.4Ma_eff^2 =',W(J1)
      CALL WRIMSG('PRINAXIS',CMSGNM)
      WRITE (CMSGNM,FMT='(A,0P,F8.5)')'I_2/.4Ma_eff^2 =',W(J2)
      CALL WRIMSG('PRINAXIS',CMSGNM)
      WRITE (CMSGNM,FMT='(A,0P,F8.5)')'I_3/.4Ma_eff^2 =',W(J3)
      CALL WRIMSG('PRINAXIS',CMSGNM)
      WRITE (CMSGNM,FMT='(A,0P,3F9.5,A)')'axis a_1 = (',A1,')'
      CALL WRIMSG('PRINAXIS',CMSGNM)
      WRITE (CMSGNM,FMT='(A,0P,3F9.5,A)')'     a_2 = (',A2,')'
      CALL WRIMSG('PRINAXIS',CMSGNM)
      WRITE (CMSGNM,FMT='(A,0P,3F9.5,A)')'     a_3 = (',A3,')'
      CALL WRIMSG('PRINAXIS',CMSGNM)
!***
      RETURN
    END SUBROUTINE PRINAXIS
