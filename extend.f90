    SUBROUTINE EXTEND(CMETHD,ICOMP,ICOMP2,IDVERR,IDVOUT,IOCC,IX,IY,IZ,BETADF, &
                      PHIDF,THETADF,SCRRS2,X0,MXNX,MXNY,MXNZ,MXNAT,MXNAT0,    &
                      MXN03,NAT,NAT0,NAT3,NX,NY,NZ,IPNF)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

! Arguments:

      CHARACTER(6) :: CMETHD
      INTEGER :: IDVERR,IDVOUT,MXNAT,MXNAT0,MXN03,MXNX,MXNY,MXNZ,NAT,NAT0, &
         NAT3,NX,NY,NZ
      INTEGER*2 ::        &
         ICOMP(3*MXNAT),  &
         ICOMP2(3*MXNAT), &
         IOCC(MXNAT)
      INTEGER ::    &
         IPNF(6),   &
         IX(MXNAT0), &
         IY(MXNAT0), &
         IZ(MXNAT0)
      REAL(WP) ::        &
         BETADF(MXNAT),  &
         PHIDF(MXNAT),   &
         SCRRS2(MXNAT),  &
         THETADF(MXNAT), &
         X0(3)

! Local parameters:

      INTEGER :: MX235
      PARAMETER(MX235=137)

! Local variables:

      INTEGER :: IXA,IXMAX,IXMIN,IYA,IYMAX,IYMIN,IZA,IZMAX,IZMIN,JA,JA2,NXY
      INTEGER ::    &
         NF235(MX235)
      SAVE NF235

!***********************************************************************
! Given:
!    NAT0 = number of dipoles in real target
!    IX(1-NAT0) = x-index of dipoles in desired target
!    IY(1-NAT0) = y-index of dipoles in desired target
!    IZ(1-NAT0) = z-index of dipoles in desired target
!    X0(1-3)    = x,y,z offset/d from target origin for dipole with IX=0
!                 IY=0, IZ=0
!                 i.e. x_TF = (IX + X0(1))*d
!
!    ICOMP      = composition info. for all locations in physical target
!                 on input, composition is stored according to
!                 dimensioning ICOMP(MXNAT,3):
!                 elements 1 - NAT0                 = comps 1x,2x,3x,..,NAT0x,
!                 elements NAT0+1 - MXNAT           = 0
!                 elements MXNAT+1 - MXNAT+NAT0     = comps 1y,2y,3y,..,NAT0y
!                 elements MXNAT+NAT0+1 - 2*MXNAT   = 0
!                 elements 2*MXNAT+1 - 2*MXNAT+NAT0 = comps 1z,2z,3z,..,NAT0z
!                 elements 2*MXNAT+NAT0+1 - 3*MXNAT = 0
!
!    BETADF(1-NAT0) = constituent material orientation angle BETA
!                     (radians) for Dielectric Frame (DF) relative to
!                     Target Frame (TF)
!    PHIDF(1-NAT0)  = constituent material orientation angle PHI (radian
!                     for DF relative to TF
!    THETADF(1-NAT0)= constituent material orientation angle THETA
!                     (radians) for DF relative to TF
!    SCRRS2         = scratch array
!    MXNX,MXNY,MXNZ = maximum allowed target size in x,y,z directions
!    MXNAT,MXNAT0,MXN03 = dimensioning information
!    IPNF(1)        = padding in -x directions to create near-field zone
!    IPNF(2)        = padding in +x directions to create near-field zone
!    IPNF(3)        = padding in -Y directions to create near-field zone
!    IPNF(4)        = padding in +Y directions to create near-field zone
!    IPNF(5)        = padding in -Z directions to create near-field zone
!    IPNF(6)        = padding in +Z directions to create near-field zone

! Returns:
!    NAT = number of dipoles in extended target = NX*NY*NZ
!    NAT3 = 3*NAT
!    IX(1-NAT0) = x-index of dipoles in reordered physical target
!    IY(1-NAT0) = y-index of dipoles in reordered physical target
!    IZ(1-NAT0) = z-index of dipoles in reordered physical target
!    ICOMP(  1     -  NAT) = x-composition ident. for sites (ICOMP=0 if vacant)
!    ICOMP(NAT+1   - 2*NAT)= y-composition ident. for sites (ICOMP=0 if vacant)
!    ICOMP(2*NAT+1 - 3*NAT)= z-composition ident. for sites (ICOMP=0 if vacant)
!                     = composition inf. for dipoles (ICOMP=0 if vacant)
!                       [with storage locations according to
!                        dimensioning ICOMP(NAT,3) ]
!    IOCC(J=1-NAT3) = 0 or 1 depending on whether site (JX,JY,JZ) is
!                            unoccupied or occupied
!                     J=NX*NY(JZ-1)+NX*(JY-1)+JX
!    BETADF(1-NAT0)= material orientation angle BETA (radians) for
!                    reordered physical target
!    PHIDF(1-NAT0)=material orientation angle PHI (radians) for
!                  reordered physical target
!    THETADF(1-NAT0)=material orientation angle THETA (radians) for
!                    reordered target
!             NB: BETADF(JA)=PHIDF(JA)=THETADF(JA)=0 at unoccupied sites
!    NX = range of IX values (1 to NX)
!    NY = range of IY values (1 to NY)
!    NZ = range of IZ values (1 to NZ)

! Note: it is necessary to reorder IX,IY,IZ because DDSCAT calls routine
!    REDUCE to take list of polarizations for "extended" target and
!    collapse it to list of polarizations for "physical" target.
!    Similarly, REDUCE is also used to collapse vector of incident E
!    and vector of polarizabilities.

! B.T.Draine, Princeton Univ. Obs., 90.10.30
! History:
! 90.11. 8 (BTD): Changed DO...ENDDO to DO #...# CONTINUE
! 90.11.30 (BTD): Cosmetic changes.
! 90.11.30 (BTD): Modify to choose NX,NY,NZ from NF235.
! 90.12.05 (BTD): Change argument list; correct bug.
! 90.12.07 (BTD): Previous version was stupidly written. Rewritten.
! 90.12.13 (BTD): Corrected data list for NF235. Modified to allow
!                 CMETHD='CONVEX' option.
! 92.06.01 (BTD+PJF): Corrected major flaw in EXTEND: previous version
!                 neglected to reorder elements of IX,IY,IZ to
!                 correspond to possible reordering of occupied sites
!                 when extended target is constructed.  This error was
!                 previously overlooked because ordering was fortuitousl
!                 unchanged by extend for targets generated by TARREC,
!                 TARELL, TARCYL, TARHEX.  Note that differential
!                 scattering cross sections computed previously using
!                 DDSCAT ver4a and using TARTET are not correct.
! 92.09.09 (BTD): extended values of NF235 up to 1024, and added
!                 warning if calling with IXMAX,IYMAX,or IZMAX > 1024
! 93.01.20 (BTD): added MXNX,MXNY,MXNZ to argument list and
!                 added tests to ensure that target is consistent with
!                 limits MXNX,MXNY,MXNZ
! 93.03.12 (BTD): extended NF235 up to 4096
! 94.06.20 (PJF): added code for NEWTMP 3D FFT
! 96.10.18 (BTD): changed NEWTMP to GPFAFT (Generalized Prime Factor
!                 Algorithm for Fourier Transform).
! 98.08.09 (BTD): Modified to permit execution to continue with warning
!                 if one or more of NX,NY,NZ exceeds MXNX,MXNY,MXNZ
!                 but NX*NY*NZ < MXNX*MXNY*MXNZ
!                 since we believe that memory management should still
!                 be ok in this case
!                 Also eliminate superfluous warning/error message
!                 about NX,NY,NZ
! 00.06.23 (BTD): Modified to support option FFTWFJ
! 04.09.15 (BTD): Modified to treat material orientation angles
!                 BETADF,PHIDF,THETADF
!                 Add scratch array SCRRS1 to use in doing this.
! 06.12.24 (BTD): Increased size of array NF235 to MX235=137
!                 Added 1 to array NF235
! 07.06.19 (BTD): Added dipole offset X0(3) to argument list
!                 Modified to shift X0 when changing IX,IY,IZ
!                 [Question: do we really need to shift IX,IY,IZ?]
! 07.09.11 (BTD): Corrected shift of X0(3)
!                 Changed IX,IY,IZ from INTEGER*2 to INTEGER
! 11.08.16 (BTD) v7.2.1
!                * add IPNF to argument list, and use this to
!                  pad space around target for nearfield calculation
! end history
! Copyright (C) 1993,1996,1998,2000,2004,2006,2007,2011
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
! Elements of NF235 are numbers.le.4096  which can be factored by 2,3,5:

      DATA NF235/1,2,3,4,5,6,8,9,10,12,15,16,18,20,24,25,27,      &
         30,32,36,40,45,48,50,54,60,64,72,75,80,81,90,96,100,     &
         108,120,125,128,135,144,150,160,162,180,192,200,216,225, &
         240,243,250,256,270,288,300,320,324,360,375,384,400,405, &
         432,450,480,486,500,512,540,576,600,625,640,648,675,720, &
         729,750,768,800,810,864,900,960,972,1000,1024,1080,1125, &
         1152,1200,1215,1250,1280,1296,1350,1440,1458,1500,1536,  &
         1600,1620,1728,1800,1875,1920,1944,2000,2025,2048,2160,  &
         2187,2250,2304,2400,2430,2500,2560,2592,2700,2880,2916,  &
         3000,3072,3125,3200,3240,3375,3456,3600,3645,3750,3840,  &
         3888,4000,4050,4096/

!*** First determine max., min values of input IX,IY,IZ for occupied sites:

!*** diagnostic
!      write(0,*)'extend ckpt 0'
!      write(0,*)'   entered extend with'
!      write(0,*)'   mxnx,mxny,mxnz=',mxnx,mxny,mxnz
!      write(0,*)'      mxnat=',mxnat
!      write(0,*)'     mxnat0=',mxnat0
!      write(0,*)'   nat0,nat=',nat0,nat
!      write(0,*)'       nat3=',nat3
!      write(0,*)'   nx,ny,nz=',nx,ny,nz
!      write(0,*)'  ipnf(1-6)=',ipnf
!      write(0,*)'      ix(1)=',ix(1)
!      write(0,*)'      iy(1)=',iy(1)
!      write(0,*)'      iz(1)=',iz(1)
!      write(0,*)'   j icomp(j), j=1 - nat3=',nat3, &
!                ' [may not yet be in natural order]'
!      do ja=1,3*nat
!         write(0,fmt='(i5,i2)')ja,icomp(ja)
!      enddo
!***

      IXMIN=IX(1)
      IXMAX=IX(1)
      IYMIN=IY(1)
      IYMAX=IY(1)
      IZMIN=IZ(1)
      IZMAX=IZ(1)
      DO JA=2,NAT0
         IF(IX(JA)<IXMIN)IXMIN=IX(JA)
         IF(IX(JA)>IXMAX)IXMAX=IX(JA)
         IF(IY(JA)<IYMIN)IYMIN=IY(JA)
         IF(IY(JA)>IYMAX)IYMAX=IY(JA)
         IF(IZ(JA)<IZMIN)IZMIN=IZ(JA)
         IF(IZ(JA)>IZMAX)IZMAX=IZ(JA)
      ENDDO

!*** diagnostic
!      write(0,*)'extend ckpt 2'
!***

!*** Now shift so that 
!   IX for physical sites runs from IPNF(1) to IPNF(1)+IXMAX-IXMIN+1
!   IY for physical sites runs from IPNF(3) to IPNF(3)+IYMAX-IYMIN+1
!   IZ for physical sites runs from IPNF(5) to IPNF(5)+IZMAX-IZMIN+1

      DO JA=1,NAT0
         IX(JA)=IX(JA)+1-IXMIN+IPNF(1)
         IY(JA)=IY(JA)+1-IYMIN+IPNF(3)
         IZ(JA)=IZ(JA)+1-IZMIN+IPNF(5)
      ENDDO
      IXMAX=IXMAX-IXMIN+1
      IYMAX=IYMAX-IYMIN+1
      IZMAX=IZMAX-IZMIN+1

! we have shifted IX,IY,IZ
! therefore need to change X0(1-3)
! where X0(1-3)*d = physical location in TF of dipole with IX=IY=IZ=0

      X0(1)=X0(1)+REAL(IXMIN-1-IPNF(1),KIND=WP)
      X0(2)=X0(2)+REAL(IYMIN-1-IPNF(3),KIND=WP)
      X0(3)=X0(3)+REAL(IZMIN-1-IPNF(5),KIND=WP)

!*** diagnostic
!      write(0,*)'extend ckpt 3'
!      write(0,*)'   ixmax,iymax,izmax=',ixmax,iymax,izmax
!      write(0,*)'   ipnf(1-6)=',ipnf
!      write(0,fmt='(a,f7.3,a,f7.3,a,f7.3)')'   new x0(1-3)=', &
!                  x0(1),',',x0(2),',',x0(3)
!***

      WRITE(IDVOUT,6100)IXMAX,IYMAX,IZMAX

!*** Now determine NX,NY,NZ

      IF(CMETHD=='GPFAFT')THEN

! Using GPFAFT : Temperton's new GPFA routine
! which requireS that NX,NY,NZ be factorable by 2,3,5
! We require that NX,NY,NZ be selected from list NF235
! Check that list is long enough:

         IF(IXMAX>NF235(MX235))THEN
            WRITE(IDVERR,*)'EXTEND> Error: IXMAX=',IXMAX, &
                           '.gt.NF235(MX235)=',NF235(MX235)
            CALL ERRMSG('FATAL','EXTEND',' need to add to NF235 in EXTEND')
         ELSEIF(IYMAX>NF235(MX235))THEN
            WRITE(IDVERR,*)'EXTEND> Error: IYMAX=',IYMAX, &
                           '.gt.NF235(MX235)=',NF235(MX235)
            CALL ERRMSG('FATAL','EXTEND',' need to add to NF235 in EXTEND')
         ELSEIF(IZMAX>NF235(MX235))THEN
           WRITE(IDVERR,*)'EXTEND> Error: IZMAX=',IZMAX, &
                          '.gt.NF235(MX235)=',NF235(MX235)
           CALL ERRMSG('FATAL','EXTEND',' need to add to NF235 in EXTEND')
         ENDIF
         NX=0
         NY=0
         NZ=0
         DO JA=1,MX235
            IF(IXMAX+IPNF(1)+IPNF(2)<=NF235(JA).AND.NX==0)NX=NF235(JA)
            IF(IYMAX+IPNF(3)+IPNF(4)<=NF235(JA).AND.NY==0)NY=NF235(JA)
            IF(IZMAX+IPNF(5)+IPNF(6)<=NF235(JA).AND.NZ==0)NZ=NF235(JA)
         ENDDO

!*** diagnostic
!         write(0,*)'extend ckpt 4'
!         write(0,*)'   nx,ny,nz=',nx,ny,nz
!***
      ELSE

! Using a routine with no restriction on NX,NY,NZ:
! e.g.,
! -- FFTWFJ : FFTW routine from Frigo and Johnson

         NX=IXMAX+IPNF(1)+IPNF(2)
         NY=IYMAX+IPNF(3)+IPNF(4)
         NZ=IZMAX+IPNF(5)+IPNF(6)
      ENDIF
      IF(NX>MXNX.OR.NY>MXNY.OR.NZ>MXNZ)THEN
         IF(NX*NY*NZ>MXNX*MXNY*MXNZ)THEN
            WRITE(IDVOUT,6300)NX,NY,NZ,MXNX,MXNY,MXNZ
            CALL ERRMSG('FATAL','EXTEND', &
                 ' need to change MXNX,MXNY,MXNZ in DDSCAT and recompile')
         ELSE
            WRITE(IDVOUT,6301)NX,NY,NZ,MXNX,MXNY,MXNZ
         ENDIF
      ENDIF
      NAT=NX*NY*NZ
      NAT3=3*NAT

!*** Now extend target: generate vectors of dipole properties
!    for all sites in extended target
!        ICOMP(jx,jy,jz,1-3)=ICOMP(ja,1-3)=ICOMP(ja3)
!        IOCC(jx,jy,jz)=IOCC(ja)
!    with indices jx,jy,jz running from 1-NX,1-NY,1-NZ
!              or ja          "     "   1-(NX*NY*NZ)
!              or ja3         "     "   1-(3*NX*NY*NZ)


! First initialize ICOMP, IOCC, and SCRRS1 to zero:

!*** diagnostic
!      write(0,*)'extend ckpt 5'
!      write(0,*)'    NAT=',nat,' NAT0=',nat0
!***

      DO JA2=1,NAT
         IOCC(JA2)=0
         ICOMP2(JA2)=0
         ICOMP2(JA2+NAT)=0
         ICOMP2(JA2+2*NAT)=0
         SCRRS2(JA2)=0._WP
      ENDDO

!*** diagnostic
!      write(0,*)'extend ckpt 6, nat0=',nat0
!      write(0,*)'       ix(1)=',ix(1)
!      write(0,*)'       iy(1)=',iy(1)
!      write(0,*)'       iz(1)=',iz(1)
!***
! Now reset values at occupied sites

      NXY=NX*NY
      DO JA=1,NAT0
!*** diagnostic
!         write(0,*)'extend ckpt 6.1 ja=',ja,' ix,iy,iz=',ix(ja),iy(ja),iz(ja)
!***
         JA2=IX(JA)+NX*(IY(JA)-1)+NXY*(IZ(JA)-1)
!*** diagnostic
!         write(0,*)'extend ckpt 6.2 ja2=',ja2
!***
         IOCC(JA2)=1
         ICOMP2(JA2)=ICOMP(JA)
         ICOMP2(JA2+NAT)=ICOMP(JA+MXNAT)
         ICOMP2(JA2+2*NAT)=ICOMP(JA+2*MXNAT)
         SCRRS2(JA2)=BETADF(JA)
      ENDDO

!*** diagnostic
!      write(0,*)'extend ckpt 7'
!***

! Now overwrite original ICOMP

      DO JA=1,NAT3
         ICOMP(JA)=ICOMP2(JA)
      ENDDO

! Overwrite original BETADF

      DO JA=1,NAT
         BETADF(JA)=SCRRS2(JA)
      ENDDO

!*** diagnostic
!      write(0,*)'extend ckpt 8'
!***

! SCRRS2 should still be zero at unoccupied sites, so no need
! to rezero this array.

! Overwrite original PHIDF

      DO JA=1,NAT0
         JA2=IX(JA)+NX*(IY(JA)-1)+NXY*(IZ(JA)-1)
         SCRRS2(JA2)=PHIDF(JA)
      ENDDO
      DO JA=1,NAT
         PHIDF(JA)=SCRRS2(JA)
      ENDDO

! Overwrite original THETADF

      DO JA=1,NAT0
         JA2=IX(JA)+NX*(IY(JA)-1)+NXY*(IZ(JA)-1)
         SCRRS2(JA2)=THETADF(JA)
      ENDDO
      DO JA=1,NAT
         THETADF(JA)=SCRRS2(JA)
      ENDDO

!*** diagnostic
!      write(0,*)'extend ckpt 9'
!***

! Now reorder IX,IY,IZ to correspond to new ordering of occupied sites

      JA=0
      DO JA2=1,NAT
         IF(IOCC(JA2)>=1)THEN
            JA=JA+1
            IXA=1+MOD(JA2-1,NX)
            IYA=1+MOD(JA2-IXA,NXY)/NX
            IZA=1+(JA2-IXA-NX*(IYA-1))/NXY
            IX(JA)=IXA
            IY(JA)=IYA
            IZ(JA)=IZA
         ENDIF
      ENDDO

!*** diagnostic
!      write(0,*)'extend ckpt 10'
!***

      RETURN
6100  FORMAT(' >EXTEND target extent in x,y,z directions =',I5,',',I5,',', &
        I5,' (in Target Frame)')
6300  FORMAT(' >EXTEND: Fatal error after extending target!:',/,            &
        '          Extended target extent in x,y,z directions =',I5,',',I5, &
        ',',I5,/,'          exceeds dimension limits MXNX,MXNY,MXNZ=',I5,   &
        ',',I5,',',I5)
6301  FORMAT(' >EXTEND: Warning after extending target!:',/,                &
        '          Extended target extent in x,y,z directions =',I5,',',I5, &
        ',',I5,/,'          exceeds dimension limits MXNX,MXNY,MXNZ=',I5,   &
        ',',I5,',',I5,/,'          execution allowed to continue since ',   &
        'MXNX*MXNY*MXNZ > NX*NY*NZ')
    END SUBROUTINE EXTEND
