    SUBROUTINE MPI_BCAST_CHAR(CHARVAR,ICOUNT,IROOT,IERR)
      IMPLICIT NONE
!-----------------------------------------------------------------------

      INCLUDE 'mpif.h'

!-----------------------------------------------------------------------
! arguments
      CHARACTER :: CHARVAR*(*)
     INTEGER :: ICOUNT,IERR,IROOT
!-----------------------------------------------------------------------
! Subroutine MPI_BCAST_CHAR
! Purpose: to serve as a "jacket" for passing character variables to
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

      CALL MPI_BCAST(CHARVAR,ICOUNT,MPI_CHARACTER,IROOT,MPI_COMM_WORLD,IERR)
      RETURN
      END
