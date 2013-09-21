    SUBROUTINE MPI_BCAST_INT2(INT2VAR,ICOUNT,IROOT,IERR)
      IMPLICIT NONE
!-----------------------------------------------------------------------

      INCLUDE 'mpif.h'

!-----------------------------------------------------------------------
! arguments
      INTEGER*2 :: INT2VAR(*)
      INTEGER :: ICOUNT,IERR,IROOT
!-----------------------------------------------------------------------
! Subroutine MPI_BCAST_INT2
! Purpose: to serve as a "jacket" for passing integer*2 variables
!          to MPI_BCAST so that MPI_BCAST will be called with only a 
!          single variable type in a given fortran routine.
! B.T. Draine, Princeton University Observatory, 2004.04.09
! history
! 04.04.09 (BTD) first written
! 08.01.17 (BTD) f90 version
! end history
!-----------------------------------------------------------------------
!*** diagnostic
!      write(0,*)'mpi_bcast_int2 ckpt 1: ICOUNT=',icount
!***
      CALL MPI_BCAST(INT2VAR,ICOUNT,MPI_INTEGER2, &
                     IROOT,MPI_COMM_WORLD,IERR)
!*** diagnostic
!      write(0,*)'mpi_bcast_int2 ckpt 2'
!***
      RETURN
    END SUBROUTINE MPI_BCAST_INT2
