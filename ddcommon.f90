    MODULE DDCOMMON_1
      USE DDPRECISION,ONLY: WP
      REAL(WP) :: &
         AKR(3),  &
         DX(3)
    END MODULE DDCOMMON_1
!------------------------------------------------------------------------------
    MODULE DDCOMMON_2
      USE DDPRECISION,ONLY: WP
      COMPLEX(WP),ALLOCATABLE :: CXADIA(:,:)
    END MODULE DDCOMMON_2
!------------------------------------------------------------------------------
    MODULE DDCOMMON_3
      USE DDPRECISION,ONLY: WP
      COMPLEX(WP),ALLOCATABLE :: CXZC(:,:,:,:)
    END MODULE DDCOMMON_3
!------------------------------------------------------------------------------
    MODULE DDCOMMON_4
      USE DDPRECISION,ONLY: WP
      COMPLEX(WP),ALLOCATABLE :: CXZW(:,:,:,:)
    END MODULE DDCOMMON_4
!------------------------------------------------------------------------------
    MODULE DDCOMMON_5
      INTEGER*2,ALLOCATABLE :: IOCC(:)
    END MODULE DDCOMMON_5
!------------------------------------------------------------------------------
    MODULE DDCOMMON_6
      USE DDPRECISION,ONLY: WP
      INTEGER :: MXNATF,MXNXF,MXNYF,MXNZF,NAT,NAT3,NAT0,NX,NY,NZ,MXN3F, &
                 IDVOUT,IPBC
      REAL(WP) :: GAMMA,PYD,PZD
    END MODULE DDCOMMON_6
!------------------------------------------------------------------------------
    MODULE DDCOMMON_7
      USE DDPRECISION,ONLY: WP
      COMPLEX(WP),ALLOCATABLE :: CXAOFF(:,:)
    END MODULE DDCOMMON_7
!------------------------------------------------------------------------------
    MODULE DDCOMMON_8
      CHARACTER(6) :: CMDFFT
    END MODULE DDCOMMON_8
!------------------------------------------------------------------------------
    MODULE DDCOMMON_9
      USE DDPRECISION,ONLY: WP
      INTEGER :: IDVOUT2,ITERMX,ITERN
      REAL(WP) :: ERRSCAL
    END MODULE DDCOMMON_9
!------------------------------------------------------------------------------
    MODULE DDCOMMON_10
      INTEGER :: MYID
    END MODULE DDCOMMON_10
!------------------------------------------------------------------------------
