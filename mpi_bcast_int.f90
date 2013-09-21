    SUBROUTINE MPI_BCAST_INT(INTVAR,ICOUNT,IROOT,IERR)
      IMPLICIT NONE
!-----------------------------------------------------------------------

      INCLUDE 'mpif.h'

!-----------------------------------------------------------------------
! arguments
      INTEGER :: ICOUNT,IERR,IROOT
      INTEGER :: INTVAR(*)
!-----------------------------------------------------------------------
! Subroutine MPI_BCAST_INT
! Purpose: to serve as a "jacket" for passing integer variables to
!          MPI_BCAST so that MPI_BCAST will be called with only a single 
!          variable type in a given fortran routine.
! B.T. Draine, Princeton University Observatory, 2004.04.09
! history
! 04.04.09 (BTD) first written
! 08.01.17 (BTD) f90 version
! end history
! Copyright (C) 2004,2008 B.T.Draine and P.J. Flatau
! This code is covered by the GNU General Public License
!-----------------------------------------------------------------------

      CALL MPI_BCAST(INTVAR,ICOUNT,MPI_INTEGER, &
                     IROOT,MPI_COMM_WORLD,IERR)
      RETURN
    END SUBROUTINE MPI_BCAST_INT
