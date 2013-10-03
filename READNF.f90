    SUBROUTINE READNF(CFLENAME,IDVOUT,CSTAMP,VERSNUM,NRFLDB,AEFF,DPHYS, &
                      NX,NY,NZ,NAT0,X0,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,   &
                      NAMBIENT,WAVE,AK_TF,CXE0_TF,CXB0_TF,NCOMP)        !

!------------------------- subroutine readnf v2 -------------------------
      USE DDPRECISION,ONLY: WP

! modules used to transfer allocatable variables:

      USE READNF_ECOM,ONLY: CXADIA,CXEINC,CXEPS,CXESCA,CXPOL,ICOMP
      USE READNF_BCOM,ONLY: CXBINC,CXBSCA
      IMPLICIT NONE

! arguments

      CHARACTER*60 :: CFLENAME
      CHARACTER*26 :: CSTAMP

      INTEGER :: IDVOUT,NAT0,NCOMP,NRFLDB,NX,NY,NZ,VERSNUM

      REAL(WP) ::                                             &
         AEFF,DPHYS,NAMBIENT,WAVE,XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN

      REAL(WP) ::  &
         AK_TF(3), &
         X0(3)     !

      COMPLEX(WP) :: &
         CXB0_TF(3), &
         CXE0_TF(3)  !

! local variables

      LOGICAL INIT

      INTEGER IC,ILINE,IOBIN,IX1,IY1,IZ1,J,J1,JX,JY,JZ, &
         K,NAT3,RWORD,NRWORD,NRWORD_NF,NXY,NXYZ         !

      REAL(WP) ::                               &
         EINC2,PHI0,PHIYZ,PHIZ,PI,SUMERR2,TINY, &
         W1,W2,W3,W4,W5,W6,W7,W8,WX,WY,WZ,      &
         XA,XB,YA,YB,ZA,ZB,ZETA                 !

      COMPLEX(WP) :: CXERR,CXFAC,CXI,CXPHAS
!==========================================================================
! subroutine READNF
! purpose: to read file with stored polarization and EM field
!          produced by DDSCAT v7.3.0
!          to support postprocessing
! given:
!    CFLENAME = name of file with stored polarization and EM field
!    IDVOUT   = unit number for output (e.g., 7)

! returns
!
! via arguments:
!    CSTAMP     = CHARACTER*26 string with DDSCAT version used to create file
!                 CFLENAME
!                 e.g., 'DDSCAT 7.3.0 [13.03.18]'
!    VERSNUM    = integer string with version number
!                 e.g., 730 for version 7.3.0
!    NRFLDB     = 0 if scattered magnetic field was not stored in file CFLENAME
!               = 1 if scattered magnetic field was stored in file CFLENAME
!    AEFF       = effective radius of target (physical units)
!    DPHYS      = d = dipole spacing (physical units)
!    NX,NY,NZ   = dimensions of computational volume
!    NAT0       = number of occupied sites (dipoles)
!    X0(3)      = x/d,y/d,z/d for site with indices (I,J,K)=(0,0,0)
!    XMIN,XMAX  = x_min,x_max (physical units) for computational volume (in TF)
!    YMIN,YMAX  = y_min,y_max (physical units) for computational volume (in TF)
!    ZMIN,ZMAX  = z_min,z_max (physical units) for computational volume (in TF)
!    NAMBIENT   = refractive index of ambient medium (real)
!    WAVE       = wavelength **in vacuo** (physical units)
!    AK_TF(3)   = (k_x,k_y,k_z)*d in the Target Frame, where 
!                 (k_x,k_y,k_z) = 2*pi*NAMBIENT/WAVE = wavevector in the medium
!    CXE0_TF(3) = complex (E_x,E_y,E_z)_TF at (x,y,z)_TF=0 and t=0
!                 for the incident plane wave
!    CXB0_TF(3) = complex (B_x,B_y,B_z)_TF at (x,y,z)_TF=0 and t=0
!                 for the incident plane wave
!    NCOMP      = number of compositions
!
! via module READNF_ECOM:
!    CXADIA(J,3)= complex diagonal elements of the "A matrix"
!               = d^3/diagonal elements of complex polarizability tensor
!                 at lattice sites J=1-NX*NY*NZ
!    CXEINC(J,3)= complex "macroscopic" (E_x,E_y,E_z)_TF of incident wave at
!                 lattice  site J
!    CXEPS(IC)  = complex dielectric function for composition IC=1-NCOMP
!    CXESCA(J,3)= complex "macroscopic" (E_x,E_y,E_z)_TF of scattered wave at
!                 lattice site J
!    CXPOL(J,3) = complex polarization (P_x,P_y,P_z)_TF of dipole J
!    ICOMP(J,3) = integer*2 composition identifier for lattice site J
!               = 0 for ambient medium
!
! via module READNF_BCOM (only if NRFLDB=1):
!    CXBINC(J,3)= complex (B_x,B_y,B_z)_TF of incident wave at lattice site J
!    CXBSCA(J,3)= complex (B_x,B_y,B_z)_TF of scattered wave at lattice site J
!
! NB: macroscopic and microscopic E fields are related by
!                   3
!    E_macro = ----------- * E_micro
!              (epsilon+2)
!
! The dipoles respond to E_micro : P = alpha*E_micro
!
!===============================================================================
! B.T. Draine, Princeton University, 2013.03.20
! history
! 13.03.20 (BTD) v1 written to support postprocessing DDSCAT output
!                by program POSTPROCESS
! end history
!===============================================================================
      DATA INIT/.FALSE./,CXI/(0.,1._WP)/
      SAVE CXI,INIT,NRWORD,PI
!=================================================================
! determine word length
      IF(INIT)THEN

! elements of READNF_ECOM:

         DEALLOCATE(CXADIA)
         DEALLOCATE(CXEINC)
         DEALLOCATE(CXEPS)
         DEALLOCATE(CXESCA)
         DEALLOCATE(CXPOL)
         DEALLOCATE(ICOMP)

! elements of READNF_BCOM

         IF(NRFLDB==1)THEN
            DEALLOCATE(CXBINC)
            DEALLOCATE(CXBSCA)
         ENDIF
      ENDIF
      IF(.NOT.INIT)THEN
         IF(WP==KIND(0.E0))THEN
            NRWORD=4
         ELSEIF(WP==KIND(0.D0))THEN
            NRWORD=8
         ELSE
            WRITE(0,*)'Fatal error determining word length in READNF'
            STOP
         ENDIF
         WRITE(IDVOUT,FMT='(A,I2)')'>READNF word length=',NRWORD
         PI=4._WP*ATAN(1._WP)
         INIT=.TRUE.
      ENDIF

!              >>>>> Important Note! <<<<<
! The structure of the READ statements below *must* conform to the
! structure of the corresponding WRITE statements in nearfield.f90 
! Any changes must be made in both modules.

      IOBIN=17
      OPEN(UNIT=IOBIN,FILE=CFLENAME,ACCESS='STREAM')
!*** diagnostic
!      write(0,*)'readnf ckpt 3'
!***
      READ(IOBIN)CSTAMP,VERSNUM
      WRITE(IDVOUT,FMT='(2A)')'>READNF data from ',CSTAMP
!*** diagnostic
!      write(0,*)'readnf ckpt 4, cstamp=',cstamp
!      write(0,*)'              versnum=',versnum
!***
      IF(VERSNUM.EQ.730)THEN
!*** diagnostic
!         write(0,*)'readnf ckpt 5, about to read file'
!***
         READ(IOBIN)NRWORD_NF,NRFLDB,NXYZ,NAT0,NAT3,NCOMP,NX,NY,NZ,X0,AEFF, &
                    NAMBIENT,WAVE,AK_TF,CXE0_TF,CXB0_TF
!*** diagnostic
!         write(0,*)'readnf ckpt 6, nrword_nf=',nrword_nf
!***
      ELSE
         WRITE(0,FMT='(3A,I4)')'file=',CFLENAME,                 &
                                ' was written by version=',VERSNUM
         WRITE(0,FMT='(2A)')'file is incompatible with present version of ', &
                             'subroutine SUBREADNF'                           !
         STOP
      ENDIF
      IF(NRWORD_NF.NE.NRWORD)THEN
         WRITE(0,*)'READNF fatal error:'
         WRITE(0,*)'  word length=',NRWORD_NF,' in file',CFLENAME
         WRITE(0,*)'  word length=',NRWORD,' in subroutine READNF'
         STOP
      ENDIF
!*** diagnostic
!      write(0,*)'readnf ckpt 7, begin allocation'
!***
      ALLOCATE(CXEPS(1:NCOMP))
      ALLOCATE(ICOMP(1:NX,1:NY,1:NZ,1:3))
      ALLOCATE(CXPOL(1:NX,1:NY,1:NZ,1:3))
      ALLOCATE(CXESCA(1:NX,1:NY,1:NZ,1:3))
      ALLOCATE(CXEINC(1:NX,1:NY,1:NZ,1:3))
      ALLOCATE(CXADIA(1:NX,1:NY,1:NZ,1:3))
!*** diagnostic
!      write(0,*)'readnf ckpt 8, end allocation, begin reading arrays'
!***
      READ(IOBIN)CXEPS
!*** diagnostic
!      write(0,*)'readnf ckpt 9, have read cxeps'
!***
      READ(IOBIN)ICOMP
      READ(IOBIN)CXPOL
      READ(IOBIN)CXESCA
      READ(IOBIN)CXADIA
      IF(NRFLDB==1)THEN
         ALLOCATE(CXBINC(1:NX,1:NY,1:NZ,1:3))
         ALLOCATE(CXBSCA(1:NX,1:NY,1:NZ,1:3))
         READ(IOBIN)CXBSCA
      ENDIF
      CLOSE(IOBIN)
!*** diagnostic
!      write(0,*)'readnf ckpt 10'
!***

! compute phase for JX=JY=JZ=0

      PHI0=AK_TF(1)*X0(1)+AK_TF(2)*X0(2)+AK_TF(3)*X0(3)
      DO JZ=1,NZ
         PHIZ=PHI0+JZ*AK_TF(3)
         DO JY=1,NY
            PHIYZ=PHIZ+JY*AK_TF(2)
            DO JX=1,NX
!*** diagnostic
!               write(0,*)'readnf ckpt 12,j=',j
!***
               CXPHAS=EXP(CXI*(JX*AK_TF(1)+PHIYZ))
!*** diagnostic
!               write(0,*)'readnf ckpt 13,j=',j
!***
               DO K=1,3

! compute E_macro contributed by incident wave

                  IC=ICOMP(JX,JY,JZ,K)
                  CXFAC=1._WP
                  IF(IC>0)CXFAC=3._WP/(CXEPS(IC)+2._WP)
                  CXEINC(JX,JY,JZ,K)=CXFAC*CXE0_TF(K)*CXPHAS
!*** diagnostic
!                  write(0,*)'readnf ckpt 14.2,k=',k
!***
               ENDDO
               IF(NRFLDB==1)THEN
                  DO K=1,3
                     CXBINC(JX,JY,JZ,K)=CXB0_TF(K)*CXPHAS
                  ENDDO
               ENDIF
!*** diagnostic
!                  write(0,*)'readnf ckpt 14.3,k=',k
!***
            ENDDO
         ENDDO
      ENDDO
!*** diagnostic
!      write(0,*)'readnf ckpt 15'
!      write(0,*)'aeff=',aeff
!      write(0,*)'nat0=',nat0
!***
      DPHYS=AEFF*(4._WP*PI/(3._WP*NAT0))**(1._WP/3._WP)
      NXY=NX*NY
      XMIN=(X0(1)+1.-0.5001)*DPHYS
      XMAX=(X0(1)+NX+0.5001)*DPHYS
      YMIN=(X0(2)+1.-0.5001)*DPHYS
      YMAX=(X0(2)+NY+0.5001)*DPHYS
      ZMIN=(X0(3)+1.-0.5001)*DPHYS
      ZMAX=(X0(3)+NZ+0.5001)*DPHYS

!*** diagnostic
!      write(0,*)'readnf ckpt 16'
!      write(0,*)'dphys=',dphys
!      write(0,fmt='(A,1P6E11.3)')'xmin,xmax,ymin,ymax,zmin,zmax=', &
!                                 xmin,xmax,ymin,ymax,zmin,zmax
!***
      WRITE(IDVOUT,FMT='(A,1PE12.4,A)')'>READNF lambda=',WAVE,' physical units'
      WRITE(IDVOUT,FMT='(A,1PE12.4,A)')'>READNF  aeff =',AEFF,' physical units'
      WRITE(IDVOUT,FMT='(A,I8,A)')'>READNF  NAT0 = ',NAT0,' target dipoles'
      WRITE(IDVOUT,FMT='(A,1PE12.4,A)')'>READNF    d  =',DPHYS,' physical units'
      WRITE(IDVOUT,FMT='(A,1P2E12.4,A)')'>READNF target xmin,xmax=', &
                                        XMIN,XMAX,' physical units'  !
      WRITE(IDVOUT,FMT='(A,1P2E12.4,A)')'>READNF target ymin,ymax=', &
                                        YMIN,YMAX,' physical units'  !
      WRITE(IDVOUT,FMT='(A,1P2E12.4,A)')'>READNF target zmin,zmax=', &
                                        ZMIN,ZMAX,' physical units'  !
                                        
! check solution
! this check neglects the off-diagonal elements of A
! this will not be a valid assumption for anisotropic materials
! that do not have optical axes aligned with TF xyz axes

      EINC2=0.
      DO K=1,3
         EINC2=EINC2+CXE0_TF(K)*CONJG(CXE0_TF(K))
      ENDDO
      SUMERR2=0.
      J1=0
      DO JZ=1,NZ
         DO JY=1,NY
            DO JX=1,NX
               IF(ICOMP(JX,JY,JZ,1)>0)THEN
                  DO K=1,3

! CXEINC and CXESCA are macroscopic E fields.
! DDA equation is for microscopic E field
! E_micro = E_macro*(epsilon+2)/3
                     IC=ICOMP(JX,JY,JZ,K)
                     CXERR=CXPOL(JX,JY,JZ,K)*CXADIA(JX,JY,JZ,K)-    &
                           (CXEINC(JX,JY,JZ,K)+CXESCA(JX,JY,JZ,K))* &
                           (CXEPS(IC)+2._WP)/3._WP                  !
                     SUMERR2=SUMERR2+CXERR*CONJG(CXERR)
                  ENDDO
! count number of occupied sites as a sanity check 
                  J1=J1+1
               ENDIF
            ENDDO
         ENDDO
      ENDDO

!*** sanity check
      if(j1.ne.nat0)then
         write(0,*)'readnf sanity failure: inconsistent j1=',j1, &
                   ' and nat0=',nat0                             !
         stop
      endif
!***
      SUMERR2=SUMERR2/(REAL(NAT0)*EINC2)

! write information to unit IDVOUT (usually a log file)

      WRITE(IDVOUT,FMT='(A,1PE11.4,A)')                                &
         '>READNF',SUMERR2,' = normalized error |P/alpha-E|^2/|E_inc|^2'
      WRITE(IDVOUT,FMT='(A,1PE11.4,A)')                           &
         '>READNF',AEFF,' = AEFF (vol. equiv. radius, phys. units)'
      WRITE(IDVOUT,FMT='(A,I11,A)')                                    &
         '>READNF',NAT0,' = NAT0 (number of physical dipoles in target)'
      WRITE(IDVOUT,FMT='(A,1PE11.4,A)')                              &
         '>READNF',DPHYS,' = d = interdipole separation (phys. units)'
      WRITE(IDVOUT,FMT='(A,1PE11.4,A)')                      &
         '>READNF',WAVE,' = wavelength in vacuo (phys. units)'
      WRITE(IDVOUT,FMT='(A,1PE11.4,A)')                                        &
         '>READNF',WAVE/NAMBIENT,' = wavelength in ambient medium (phys. units)'

      RETURN
      END
