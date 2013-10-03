      SUBROUTINE NEARFIELD(CSTAMP,VERSNUM,NRWORD,CMDFFT,CFLE1,CFLE2, &
                           CFLEB1,CFLEB2,MYID,AKD,AK_TF,IDVOUT,IOCC, &
                           ICOMP,MXCOMP,MXPBC,NAT0,NCOMP,NX,NY,NZ,   &
                           IPBC,IORTH,NRFLDB,GAMMA,NAMBIENT,PYD,PZD, &
                           DX,X0,AEFF,WAVE,CXADIA,CXEPS,CXPOL1,      &
                           CXPOL2,CXE01_TF,CXE02_TF,CXBSCA1,         &
                           CXBSCA2,CXESCA1,CXESCA2,CXZC,CXZW)        !

!-------------------------------- v8 -------------------------------
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE

! Arguments:

      CHARACTER :: CMDFFT*6,CFLE1*15,CFLE2*15,CFLEB1*16,CFLEB2*16,CSTAMP*26

      INTEGER :: IDVOUT,IORTH,IPBC,JPBC,MXCOMP,MXPBC,MYID, &
                 NAT0,NCOMP,NRFLDB,NRWORD,NX,NY,NZ,VERSNUM

      INTEGER :: MXITER
      INTEGER :: & 
         ITNUM(2)

      INTEGER*2 ::          &
         ICOMP(3*NX*NY*NZ), &
         IOCC(NX*NY*NZ)

      INTEGER ::     &
         IXYZ0(NAT0,3)

      REAL(WP) :: AEFF,AK3,AKD,ETASCA,GAMMA,NAMBIENT,PIA2,PYD,PZD,TOL,TOLR,WAVE

      REAL(WP) ::  &
         AK_TF(3), &
         DX(3),    &
         X0(3)

      COMPLEX(WP) :: CXI,CXFAC
      COMPLEX(WP) ::                &
         CXADIA(1:3*NX*NY*NZ),      &
	 CXBSCA1(1:3*NX*NY*NZ),     &
         CXBSCA2(1:3*NX*NY*NZ),     &
         CXE01_TF(1:3),             &
         CXE02_TF(1:3),             &
         CXEPS(1:MXCOMP),           &
         CXESCA1(1:3*NX*NY*NZ),     &
         CXESCA2(1:3*NX*NY*NZ),     &
         CXPOL1(1:3*NX*NY*NZ),      &
         CXPOL2(1:3*NX*NY*NZ),      &
         CXZC(NX+1+MXPBC*(NX-1),    &
              NY+1+MXPBC*(NY-1),    &
              NZ+1+MXPBC*(NZ-1),6), &
         CXZW(NX*NY*NZ,24)

!-----------------------------------------------------------------------

! Local variables

      LOGICAL :: NONZERO_X

      INTEGER :: IDVOUT2,IOBIN,J1,J2,J,JJ,JX,JY,JZ,K,NAT3,NXY,NXYZ

      REAL(WP) :: DTIME,PHASYZ,PHASZ

      REAL(WP) :: &
         X(1:3)

      COMPLEX(WP) ::  &
         CXB01_TF(3), &
         CXB02_TF(3)

      COMPLEX(WP),ALLOCATABLE :: &
         CXZG(:,:,:,:)

      CHARACTER :: CMSGNM*70

!***********************************************************************
! Subroutine NEARFIELD
! given
!   NRWORD = length (bytes) of real word
!          = 4 for single precision, 8 for double precision
!   CMDFFT = FFT method
!   CFLB1  = name of output file for nearfield B, for polarization state 1
!   CFLB2  = name of output file for nearfield B, for polarization state 2
!   CFLE1  = name of output file for nearfield E, for polarization state 1
!   CFLE2  = name of output file for nearfield E, for polarization state 2
!   MYID
!   AKD        = k*d
!   AK_TF(1:3)   = (k_x,k_y,k_z)*d in target frame
!   IDVOUT     = unit to use for output
!   IOCC(1:NAT)=0/1 if site is vacant/occupied
!   ICOMP(1:3*NAT) = ICOMP(NAT,3)
!              = composition identifier for sites 1-NAT, directions x,y,z
!   MXCOMP     = dimensioning info for array CXEPS
!   MXPBC      = dimensioning info
!   NAT0       = number of occupied sites
!   NCOMP      = number of distinct compositions
!   NX,NY,NZ   = extent of lattice in x,y,z directions
!   IPBC       = 0 to do isolated target
!                1 to use periodic boundary conditions
!                  in either y or z direction, or both y and z
!   IORTH      = 1 if only doing a single polarization
!                2 if doing two polarizations for each orientation
!   NRFLDB     = 0 to omit calculation of B
!                1 to use BSELF to compute B
!   GAMMA      = convergence parameter used in PBC calculations
!   NAMBIENT   = (real) refractive index of ambient medium
!   PYD,PZD    = PBC period/d in y and z directions
!   DX(1:3)    = lattice spacing in x,y,z directions/d
!                d^3 = DX(1)*DX(2)*DX(3)
!   X0(1:3)    = location/d in TF for (I,J,K)=(0,0,0)
!   AEFF       = effective radius (phyical units) 
!                of target or target unit cell
!                aeff = (3*Volume/4*pi)^{1/3}
!   WAVE       = wavelength in vacuo (physical units)
!   CXADIA(1:3*NAT)=diagonal elements of polarizability tensor
!                for dipoles at locations J=1-NAT
!   CXEPS(1:MXCOMP)=complex dielectric function for compositions IC=1-MXCOMP
!   CXPOL1(1:3*NAT)=polarization (in TF) at locations J=1-NAT
!                produced by incident E polarization CXE01_TF
!                propagating in direction AK_TF
!   CXPOL2(1:3*NAT)=polarization (in TF) at locations J=1-NAT
!                produced by incident E polarization CXE02_TF
!                propagating in directin AK_TF [used only if IORTH=2]
!   CXE01_TF(1:3)= incident polarization vector in TF
!   CXE02_TF(1:3)= incident orthogonal polarization vector in TF
!   CXZC       = work space needed by ESELF
!   CXZW       = work space needed by ESELF

! returns

!   CXESCA1(1:3*NAT)=macroscopic electric field at lattice sites 1-NAT
!                    generated by the polarization CXPOL1
!   CXESCA2(1:3*NAT)=macroscopic electric field at lattice sites 1-NAT
!                    generated by the polarization CXPOL2 [only if IORTH=2]
!
!   CXBSCA1(1:3*NAT)=macroscopic=microscopic magnetic field at 
!                    lattice sites 1-NAT
!                    generated by the polarization CXPOL1
!                    [only if NRFLB=1]
!   CXBSCA2(1:3*NAT)=macroscopic=microscopic magnetic field at
!                    lattice sites 1-NAT
!                    generated by the polarization CXPOL2 
!                    [only if NRFLB=1 and IORTH=2]
! and writes to file
! if NRFLDB=0
!      CFLE1  : P, E_inc, E_sca for incident polarization 1
!      CFLE2                    for incident polarization 2 [only if IORTH=2]
! or
! if NRFLDB=1
!      CFLEB1 : P, E_inc, E_sca, B_inc, B_sca for incident polarization 1 
!      CFLEB2 :                               for incident polarization 2
!                                             [only if IORTH=2]


! history: created to perform nearfield calculations for version 7.2.1
! 11.08.16 (BTD) v7.2.a
!                * NEARFIELD cast into final form
!                * diagnostic write statements disabled
! 11.08.17 (BTD) * fixed bug
! 11.08.18 (BTD) * add UNREDUCE(CXADIA...)
! 11.08.30 (BTD) v7.2.b
!                * add NAMBIENT to argument list
!                * add NAMBIENT to WRITE(IDVOUT)....
!                * change from FORM='UNFORMATTED' to ACCESS='STREAM'
! 11.11.20 (BTD) v7.2.0
!                * added missing NAMBIENT to write list when IORTH=2
! 12.02.11 (BTD) v2
!                * add NRWORD to argument list
!                * reorder data in stored output
!                  (more essential first, least essential last)
!                * add AKD and CXE01_TF (or CXE02_TF) to stored output
! 12.07.06 (BTD) v2 edited comments
! 12.07.10 (BTD) v3 for DDSCAT 7.3
!                changed notation CXE01R -> CXE01_TF, CXE02R -> CXE02_TF
!                incorporated changes required to support use of
!                subroutine BSELF written by Ian Wong
!                * NRFLDB added to arg list
!                * CFLB1,CFLB2 added to arg list
!                * CXBSCA1,CXBSCA2 added to arg list
! 12.07.12 (BTD) v4
!                * additional changes
! 12.08.02 (IYW) modified by Ian Y. Wong
!                added DIPINT to arg list
! 12.08.09 (IYW) v4 
!                * added WRITE(IOBIN)ICOMP to output statements
!                  when writing out B fields
! 12.08.11 (BTD) v4 
!                * removed DIPINT from argument list
!                  IDIPINT is now passed to ESELF through module DDCOMMON_0
! 12.12.23 (BTD) v5
!                * added MXCOMP and CXEPS to argument list to allow
!                  calculation of macroscopic E field within target
! 12.12.25 (BTD) rewrite to consolidate E and B into single output file
! 13.01.03 (BTD) v6
!                * modified to *not* deallocate CXZG after first
!                  call to BSELF if IORTH=2
!                  This way BSELF can skip recomputation of Green function
!                  coefficients when called for second polarization
! 13.03.18 (BTD) v7
!                * add CSTAMP and VERSNUM to argument list
!                * write CSTAMP and VERSNUM to nearfield binary output
! 13.03.25 (BTD) v8
!                * remove CXEINC and CXBINC -- do not need to store these
!                  because they are easily generated by DDPOSTPROCESS
! end history
! Copyright (C) 2011,2012,2013 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

      CALL TIMEIT('NEARFIELD',DTIME)

      IOBIN=66
      CXI=(0._WP,1._WP)
      IDVOUT2=IDVOUT
      NXY=NX*NY
      NXYZ=NXY*NZ
      NAT3=3*NXYZ

      CALL UNREDUCE(CXADIA,IOCC,3*NXYZ,NXYZ,NXYZ,NAT0)

!*** diagnostic
!	write(0,*)'nearfield_v6 ckpt 1'
!***

      CALL UNREDUCE(CXPOL1,IOCC,3*NXYZ,NXYZ,NXYZ,NAT0)

!*** diagnostic
!	write(0,*)'nearfield_v6 ckpt 2: about to call eself'
!***

      CALL ESELF(CMDFFT,CXPOL1,NX,NY,NZ,IPBC,GAMMA,PYD,PZD,AK_TF,AKD,DX, &
                 CXZC,CXZW,CXESCA1)

!*** diagnostic
!	write(0,*)'nearfield_v6 ckpt 3'
!***

      IF(NRFLDB==1)THEN

! allocate CXZG = array that will contain Green function for computing B_sca
!                 from P

         ALLOCATE(CXZG(NX+1+IPBC*(NX-1),NY+1+IPBC*(NY-1),NZ+1+IPBC*(NZ-1),3))

!*** diagnostic
!         write(0,*)'nearfield_v6 ckpt 3.5: about to call bself'
!***
         CALL BSELF(CMDFFT,CXPOL1,NX,NY,NZ,IPBC,GAMMA,PYD,PZD,AK_TF,AKD,DX, &
                    CXZG,CXZW,CXBSCA1)

!*** diagnostic
!	 write(0,*)'nearfield_v6 ckpt 4'
!***

! If memory use is an issue, could deallocate CXZG and recalculate for IORTH=2
! In present version, deallocate only if IORTH=1

         IF(IORTH==1)DEALLOCATE(CXZG)

      ENDIF

!*** diagnostic
!	write(0,*)'nearfield_v6 ckpt 5'
!***

! evaluate incident wave for incident polarization state 1

! for adopted convention E02 = khat x E01
! we have                B01 = khat x E01 = E02

      DO K=1,3
         CXB01_TF(K)=CXE02_TF(K)
      ENDDO

! have calculated *microscopic* E field CXESCA1
! now convert to macroscopic E field
! only changes are within the target material

      DO JZ=1,NZ
         DO JY=1,NY
            DO JX=1,NX
               J=JX+(JY-1)*NX+(JZ-1)*NXY
               IF(IOCC(J)>0)THEN
                  DO K=1,3
                     JJ=J+(K-1)*NXYZ
                     CXFAC=3._WP/(2._WP+CXEPS(ICOMP(JJ)))
!*** diagnostic
!                     write(0,fmt='(a,3i3,i2,a,2f6.2)')                &
!                        'nearfield_v6 ckpt 6,jx,jy,jz,k=',jx,jy,jz,k, &
!                        ' cxeps=',cxeps(icomp(jj))
!***
                     CXESCA1(JJ)=CXFAC*CXESCA1(JJ)
                  ENDDO   ! enddo k=1,3
               ENDIF   ! endif(iocc...)
            ENDDO   ! enddo jx=1,nx
         ENDDO   ! enddo jy=1,ny
      ENDDO   ! enddo jz=1,nz

!              >>>>> Important Note! <<<<<
! The structure of the WRITE statements below *must* conform to the
! structure of the corresponding READ statements in readnf.f90 
! Any changes must be made in both modules.

      IF(NRFLDB==0)THEN
         OPEN(UNIT=IOBIN,FILE=CFLE1,ACCESS='STREAM')
      ELSE
         OPEN(UNIT=IOBIN,FILE=CFLEB1,ACCESS='STREAM')
      ENDIF

!*** diagnostic
!      write(0,*)'nearfield_v6 ckpt 7'
!      write(0,*)'    CSTAMP=',CSTAMP
!      write(0,*)'   VERSNUM=',VERSNUM
!***

      WRITE(IOBIN)CSTAMP,VERSNUM
      WRITE(IOBIN)NRWORD,NRFLDB,NXYZ,NAT0,NAT3,NCOMP,NX,NY,NZ,X0,AEFF, &
                  NAMBIENT,WAVE,AK_TF,CXE01_TF,CXB01_TF

!*** diagnostic
!      write(0,*)'nearfield_v6 ckpt 8'
!      write(0,*)'   nxyz=',nxyz
!      write(0,*)'   nat0=',nat0
!      write(0,*)'   nat3=',nat3
!      write(0,*)' cxe01r=',cxe01r
!***
      WRITE(IOBIN)CXEPS(1:NCOMP)
      WRITE(IOBIN)ICOMP(1:3*NXYZ)
      WRITE(IOBIN)CXPOL1(1:3*NXYZ)
      WRITE(IOBIN)CXESCA1(1:3*NXYZ)
      WRITE(IOBIN)CXADIA(1:3*NXYZ)

!*** diagnostic
!      write(0,*)'nearfield_v6 ckpt 10'
!***
      IF(NRFLDB==1)THEN
         WRITE(IOBIN)CXBSCA1(1:3*NXYZ)
      ENDIF
      CLOSE(IOBIN)

      IF(IORTH==2)THEN

         CALL UNREDUCE(CXPOL2,IOCC,3*NXYZ,NXYZ,NXYZ,NAT0)

!*** diagnostic
!         write(0,*)'nearfield_v6 ckpt 10.5: about to call eself'
!***
         CALL ESELF(CMDFFT,CXPOL2,NX,NY,NZ,IPBC,GAMMA,PYD,PZD,AK_TF,AKD,DX, &
                    CXZC,CXZW,CXESCA2)
!*** diagnostic
!         write(0,*)'nearfield_v6 ckpt 10.6'
!***
         IF(NRFLDB==1)THEN

! in present version, do not need to reallocate CXZG
! because deallocation above (in version v5) has been suppressed in v6

!            ALLOCATE(CXZG(NX+1+IPBC*(NX-1),NY+1+IPBC*(NY-1), &
!                          NZ+1+IPBC*(NZ-1),3))

!*** diagnostic
!            write(0,*)'nearfield_v6 ckpt 10.7: about to call bself'
!***
            CALL BSELF(CMDFFT,CXPOL2,NX,NY,NZ,IPBC,GAMMA,PYD,PZD,AK_TF,AKD,DX, &
                       CXZG,CXZW,CXBSCA2)
!*** diagnostic
!            write(0,*)'nearfield_v6 ckpt 10.8'
!***
            DEALLOCATE(CXZG)
         ENDIF

! evaluate incident wave

! for adopted convention E02 = khat x E01
! we have                B02 = khat x E02 = khat x (khat x E01) = -E01
         DO K=1,3
            CXB02_TF(K)=-CXE01_TF(K)
         ENDDO

! have calculated *microscopic* E field CXESCA1
! now convert to macroscopic E field
! only changes are within the target material

         DO JZ=1,NZ
            DO JY=1,NY
               DO JX=1,NX
                  J=JX+(JY-1)*NX+(JZ-1)*NXY
                  IF(IOCC(J)>0)THEN
                     DO K=1,3
                        JJ=J+(K-1)*NXYZ
                        CXFAC=3._WP/(2._WP+CXEPS(ICOMP(JJ)))
                        CXESCA2(JJ)=CXFAC*CXESCA2(JJ)
                     ENDDO
                  ENDIF   ! endif(iocc...)
               ENDDO   ! enddo jx=1,nx
            ENDDO   ! enddo jy=1,ny
         ENDDO   ! enddo jz=1,nz

         IF(NRFLDB==0)THEN
            OPEN(UNIT=IOBIN,FILE=CFLE2,ACCESS='STREAM')
         ELSE
            OPEN(UNIT=IOBIN,FILE=CFLEB2,ACCESS='STREAM')
         ENDIF
         WRITE(IOBIN)CSTAMP,VERSNUM
         WRITE(IOBIN)NRWORD,NRFLDB,NXYZ,NAT0,NAT3,NCOMP,NX,NY,NZ,X0,AEFF, &
                     NAMBIENT,WAVE,AK_TF,CXE02_TF,CXB02_TF
         WRITE(IOBIN)CXEPS(1:NCOMP)
         WRITE(IOBIN)ICOMP(1:3*NXYZ)
         WRITE(IOBIN)CXPOL2(1:3*NXYZ)
         WRITE(IOBIN)CXESCA2(1:3*NXYZ)
         WRITE(IOBIN)CXADIA(1:3*NXYZ)
         IF(NRFLDB==1)THEN
            WRITE(IOBIN)CXBSCA2(1:3*NXYZ)
         ENDIF
         CLOSE(IOBIN)

!*** diagnostic sanity check...
!      write(0,*)'nearfield_v6 ckpt 11: check CXESCA2 for NaN...'
!      j2=0
!      do j1=1,nat3
!         if(.not.(abs(cxesca2(j1))>=0.d0))then
!            write(0,*)'j1=',j1,' cxesca2(j1)=',cxesca2(j1)
!            j2=j2+1
!         endif
!      enddo
!      write(0,*)'... CXESCA2 checked for NaN'
!      write(0,*)'    j2=',j2,' instances'
!***

      ENDIF
      CALL TIMEIT('NEARFIELD',DTIME)
      RETURN
      END SUBROUTINE NEARFIELD

