    SUBROUTINE SCAVEC(MXSCA,NSCAT,THETAN,PHIN,CXE01,ENSC,EM1,EM2)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE
!***********************************************************************
! Given:
!      MXSCA=dimension information for ENSC,PHIN,THETAN
!      NSCAT=number of scattering directions
!      THETAN(1-NSCAT)=scattering angles theta
!      PHIN(1-NSCAT)=scattering angles phi
!      CXE01(1-3)=complex polarization vector 1 (phi=0 direction
!                 is defined by x,y plane (Lab Frame), where incident
!                 radiation propagates along the x axis.
! Returns:
!      ENSC(1-3,1-NSCAT)=scattering vectors in Lab Frame
!      EM1(1-3,1-NSCAT)=scattered pol vectors parallel to scat. plane
!                       in Lab Frame
!      EM2(1-3,1-NSCAT)=scattered pol vectors perp. to scat. plane
!                       in Lab Frame
! It is assumed that incident propagation vector is in x-direction
! in Lab Frame

! History:
! 96.11.06 (BTD): Changed definition of scattering angle phi
!                 Previously, phi was measured from plane containing
!                 incident k vector (i.e., Lab x-axis) and Re(CXE01)
!                 Henceforth, phi is measured from Lab x,y plane.
! 10.01.30 (BTD): cosmetic changes
! end history
! Copyright (C) 1993,1996,2010 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

! Arguments:
      INTEGER :: MXSCA,NSCAT
      COMPLEX(WP) :: & 
         CXE01(3)
      REAL(WP) ::       &
         EM1(3,MXSCA),  &
         EM2(3,MXSCA),  &
         ENSC(3,MXSCA), &
         PHIN(MXSCA),   &
         THETAN(MXSCA)

! Local variables:
      INTEGER :: I
      REAL(WP) :: COSPHI,COSTHE,SINPHI,SINTHE

!***********************************************************************
      DO I=1,NSCAT
        COSPHI=COS(PHIN(I))
        SINPHI=SIN(PHIN(I))
        COSTHE=COS(THETAN(I))
        SINTHE=SIN(THETAN(I))
        ENSC(1,I)=COSTHE
        ENSC(2,I)=SINTHE*COSPHI
        ENSC(3,I)=SINTHE*SINPHI
        EM1(1,I)=-SINTHE
        EM1(2,I)=COSTHE*COSPHI
        EM1(3,I)=COSTHE*SINPHI
        EM2(1,I)=0._WP
        EM2(2,I)=-SINPHI
        EM2(3,I)=COSPHI
      ENDDO
      RETURN
    END SUBROUTINE SCAVEC
