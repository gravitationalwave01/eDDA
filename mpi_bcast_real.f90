      SUBROUTINE MPI_BCAST_REAL(SINGLE,REALVAR,ICOUNT,IROOT,IERR)
      IMPLICIT NONE
!-----------------------------------------------------------------------

      INCLUDE 'mpif.h'

!-----------------------------------------------------------------------
! arguments
      LOGICAL SINGLE
      INTEGER ICOUNT,IERR,IROOT
      REAL REALVAR(*)
!-----------------------------------------------------------------------
! Subroutine MPI_BCAST_REAL
! Purpose: to serve as a "jacket" for passing real variables to
!          MPI_BCAST so that MPI_BCAST will be called with only a single 
!          variable type in a given fortran routine.

! if SINGLE = .TRUE. : real variables are single precision
!             .FALSE. : real variables are double precision

! B.T. Draine, Princeton University Observatory, 2004.04.09
! history
! 04.04.09 (BTD) first written
! 08.01.17 (BTD) f90 version
! end history
! Copyright (C) 2004,2008 B.T.Draine and P.J. Flatau
! This code is covered by the GNU General Public License
!-----------------------------------------------------------------------
      IF(SINGLE)THEN
         CALL MPI_BCAST(REALVAR,ICOUNT,MPI_REAL, &
                        IROOT,MPI_COMM_WORLD,IERR)
      ELSE
         CALL MPI_BCAST(REALVAR,ICOUNT,MPI_DOUBLE_PRECISION, &
                        IROOT,MPI_COMM_WORLD,IERR)
      ENDIF
      RETURN
    END SUBROUTINE MPI_BCAST_REAL
