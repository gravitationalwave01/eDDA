     MODULE READNF_ECOM
       USE DDPRECISION,ONLY: WP
       INTEGER*2,ALLOCATABLE :: &
          ICOMP(:,:,:,:)        !
       COMPLEX(WP),ALLOCATABLE :: &
          CXADIA(:,:,:,:),        &
          CXEINC(:,:,:,:),        &
          CXEPS(:),               &
          CXESCA(:,:,:,:),        &
          CXPOL(:,:,:,:)          !
     END MODULE READNF_ECOM      
