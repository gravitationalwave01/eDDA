!!    MODULE DDPRECISION
! Define precision of variable
!if you want to run in single-precision set  WP=KIND(0.E0)
!if you want to run in double-precision set  WP=KIND(0.D0)
!!      INTEGER,PARAMETER :: WP=KIND(0.E0)
!       INTEGER,PARAMETER :: WP=KIND(0.D0)
!!    END MODULE DDPRECISION

    module arraymodule
      use ddprecision
      implicit none
      COMPLEX(wp),allocatable::  cmat(:,:)
    end module arraymodule


    module cgmodule
    use ddprecision
    INTEGER, PARAMETER :: continue_iter = -99, converged = 0, &
        diff_iter_eps = 7, hard_breakdown = -3, left = 1, none = 0, &
        no_convergence = -1, no_eig_estimate = -7, precon_stop_error = -4, &
        right = 2, r_eps = 1, r_eps_b = 2, soft_breakdown = -2, &
        sqrt_rz_eps_b = 3, stop_error = -5, stop_error_pimcheb = -6, &
        success = 0, symmetric = 3, z_eps = 4, z_eps_b = 5, z_eps_qb = 6
    type cgstruct
     character(5) :: print
     integer:: ioerr              ! output unit 
        integer n                 ! n  - System size
        integer blksz             ! blksz - 'Size of block of data (cyclic mode) ='
        integer loclen            ! system size (typically equal to n)
        integer basisdim          ! basisdim - 'Dimension of orthogonal basis ='
        integer nprocs            ! nprocs - 'Number of processors used ='
        integer procid            ! procid - 'Processor identification ='
        integer precontype
! none=0        'No preconditioning'
! left=1        'Left preconditioning'
! right=2       'Right preconditioning'
! symmetric=3   'Symmetric preconditioning
        integer stoptype
! r_eps=1             'stopping criterion: ||r||<epsilon'
! r_eps_b=2           'Stopping criterion: ||r||<epsilon||b||'
! sqrt_rz_eps_b=3     'Stopping criterion : sqrt(r''z)<epsilon ||b||'
! z_eps=4             'Stopping criterion : ||z||<epsilon'
! z_eps_b=5           'Stopping criterion : ||z||<epsilon ||b||'
! z_eps_qb=6          'Stopping criterion : ||z||<epsilon ||qb||'
! diff_iter_eps=7     'Stopping criterion : ||x(k)-x(k-1)||<epsilon'
        integer maxit            ! maxit - maximum number of iterations
        integer itno             ! itno - actual number of iterations
        integer status
! converged           = 0   'Exit status: converged to solution'
! no_convergence      =-1   'Exit status: no convergence achieved'
! soft_breakdown      =-2   'Exit status: soft breakdown solution may have been found'
! hard_breakdown      =-3   'Exit status: hard breakdown no solution found'
! precon_stop_error   =-4   'Conflict between preconditioner and selecte stopping criterion'
! stop_error          =-5   'Error in stopping criterion, r''z<0'
        integer  steperr          ! steperr - 'Breakdown occurred at step' (step number))
        real(wp) epsilon_err      ! espilon_err - 'Value of epsilon ='
        real(wp) exitnorm         ! exitnorm - 'Value of left hand side of selected stopping criterion'
        real(wp) mineig_real     ! mineig_real - 'Minimum eigenvalue of Q1*A*Q2''
        real(wp) maxeig_real
        real(wp) mineig_imag
        real(wp) maxeig_imag
    end type cgstruct
    end module cgmodule
