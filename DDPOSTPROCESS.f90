      PROGRAM DDPOSTPROCESS

!---------------------------DDPOSTPROCESS v2 -----------------------------
! purpose: 
! to use subroutine READNF to read data from near-field files written 
! by DDSCAT and then
! postprocess as desired for visualization, etc.

! allocatable arrays are passed through modules READNF_ECOM and READNF_BCOM
! other information is passed through READNF argument list

! To the user: if desired, add additional processing code
! at bottom of program.  See comments there.

      USE DDPRECISION,ONLY: WP
      USE READNF_ECOM,ONLY: CXADIA,CXEINC,CXEPS,CXESCA,CXPOL,ICOMP
      USE READNF_BCOM,ONLY: CXBINC,CXBSCA

!VTR (this is defined in vtr.f90)
      USE VTR
      IMPLICIT NONE

      CHARACTER :: CFLENAME*60,CFLPAR*60,CFLPAR_DEFAULT*60,COMMAND(10)*60
      CHARACTER :: CSTAMP*26

!VTR output file name

      CHARACTER*60 :: CFLVTR

      INTEGER ::                                                        &
         IDVOUT,ILINE,IOBIN,IVTR,IX1,IX2,IY1,IY2,IZ1,IZ2,               &
         JA,JX,JY,JZ,K,                                                 &
         NAB,NAT0,NAT3,NCOMP,NRFLDB,NRWORD,NRWORD_NF,NX,NXY,NXYZ,NY,NZ, &
         VERSNUM                                                        !

      REAL(WP) ::                                           &
         AEFF,DPHYS,E2,EINC2,                               &
         NAMBIENT,PI,SNORM,SUMERR2,TINY,                    &
         W1,W2,W3,W4,W5,W6,W7,W8,WAVE,WX,WY,WZ,             &
         XA,XB,XMAX,XMIN,YA,YB,YMAX,YMIN,ZA,ZB,ZETA,ZMAX,ZMIN

      REAL(WP) ::    &
         AK_TF(1:3), &
         SVEC(1:3),  &
         S_INC(1:3), &
         X0(1:3),    &
         XTF(1:3)    !

      REAL(WP),ALLOCATABLE :: &
         S(:,:,:,:)

      COMPLEX(WP) :: CXBX,CXBY,CXBZ,CXERR,CXEX,CXEY,CXEZ
      COMPLEX(WP) ::   &
         CXB(1:3),     &
         CXB0_TF(1:3), &
         CXB_INC(1:3), &
         CXB_SCA(1:3), &
         CXE(1:3),     &
         CXE_INC(1:3), &
         CXE_SCA(1:3), &
         CXE0_TF(1:3), &
         CXP(1:3)      !
 
!VTR related arrays
      TYPE(VTR_FILE_HANDLE) :: fd
!VTR note that VTR graphics fields have to be in double precision
      REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: VTRX,VTRY,VTRZ        ! mesh
      REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:) :: VTR8              ! data
      REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:) :: VTRVX,VTRVY,VTRVZ ! data
!VTR supplementary file of problem geometry
      REAL(KIND=8),DIMENSION(1):: XVECT,YVECT,ZVECT
      REAL(KIND=8),DIMENSION(1,1,1):: U,V,W
      COMPLEX(WP):: CXBB(3),CXEE(3)
      INTEGER::IX,IY,IZ

!=======================================================================
! Program DDPOSTPROCESS v2
! Purpose: to read binary output files created by nearfield calculation

! given (from parameter file ddpostprocess.par)
!    CFLENAME = name of file with precomputed P, E, and possibly B
!    IVTR     = 0 to skip creation of VTR file
!             = 1 to create VTR file with |E|
!             = 2 to create VTR file with |E|^2
!             = 3 to create VTR file with time-average <Re(E)xRe(B)>/(4*pi)
!    ILINE    = 1 to evaluate E field at points along 1 or more lines
!    XA,YA,ZA = (x,y,z)_TF (physical units) for starting point of line
!    XB,YB,ZB = (x,y,z)_TF (physical units) for endpoint of line
!    NAB      = number of points along line, including A and B

! extracts from file (by using subroutine READNF):
!    AEFF             = effective radius of target (phys. units)
!    NAMBIENT         = (real) refractive index of ambient medium
!    WAVE             = wavelength in vacuo of incident wave (phys. units)
!    DPHYS            = interdipole separation (phys. units)
!    NAT0             = number of dipoles in physical target
!    NCOMP            = number of distinct compositions present
!    NX,NY,NZ         = dimensions/d of computational volume
!                       (computational volume has NXYZ=NX*NY*NZ points)
!    X0(1-3)          = (x/d,y/d,z/d) in Target frame for index I,J,K=0,0,0
!                       thus (x,y,z)=[ x0(1,2,3) + (I,J,K) ]*d
!    AK_TF(1-3)         = (k_x,k_y,k_z)*d in the Target Frame
!    CXE0_TF(1-3)       = E_inc (complex) in the Target Frame
!                       at (x_TF,y_TF,z_TF)=(0,0,0)
!    ICOMP(1-3*NXYZ)  = composition identifier for all points and
!                       directions
!    CXEINC(1-NXYZ,3) = complex incident macroscopic E field at all points
!    CXESCA(1-NXYZ,3) = complex radiated macroscopic E field at all points
!    CXPOL(1-NXYZ,3)  = complex polarization/d^3 at all points
!    CXADIA(1-NXYZ,3) = diagonal element of polarizability/d^3 at all pts

! if the stored file contained magnetic field information (NRFLDB=1)
! then also return
!
!    CXBINC(1-NXYZ,3) = complex incident B field at all points
!    CXBSCA(1-NXYZ,3) = complex radiated B field at all points
!    ICOMP(1-3*NXYZ)  = composition identifier at all points
!                       (= 0 for vacuum)
!
! using extracted information, ddpostprocess uses simple interpolation to 
! evaluate, for each of NAB points on line:
!
!    CXE_INC(1-3)= (complex) incident E field at location (x_TF,y_TF,z_TF)
!    CXE_SCA(1-3)= (complex) radiated E field at location (x_TF,y_TF,z_TF)
!    CXP(1-3)    = (complex) polarization/d^3 at location (x_TF,y_TF,z_TF)
! 
! and, if NRFLDB=1:
!
!    CXB_INC(1-3) = (complex) incident B field at location (x_TF,y_TF,z_TF)
!    CXB_SCA(1-3) = (complex) radiated B field at location (x_TF,y_TF,z_TF)
!
! current version writes out textfile ddpostprocess_E.out with one line per
! point along track:

!    x_TF,y_TF,z_TF,CXE(1-3)           [if NRFLDB=0: only E field is available]
!    x_TF,y_TF,z_TF,CXE(1-3),CXB(1-3)  [if NRFLDB=1: both E and B are available]

! where CXE=CXE_INC+CXE_SCA = total macroscopic E field at point
!       CXB=CXB_INC+CXB_SCA = total B field at point
!
! If IVTR > 0, then generate VTR output files
!
!*** NB: Existing code simply writes out E and B at points along a 
!        straight line as a simple example with limited output.
!
!        Users who wish to write out other information for purposes
!        of display or analysis should go to the end of the existing
!        program and add additional code to write out whatever is desired
!        (e.g., you may want E or B at points on a 2-D plane, or 3-D volume).
!          
!---------------------------------------------------------------------- 
! DDPOSTPROCESS is adapted from program originally first named READE, 
! then renamed READNF.  READE was first written 2011.08.30
! B.T. Draine, Princeton University Observatory
! history
! 11.08.29 (BTD) completed working version of subroutine READE
! 11.08.31 (BTD) * added NAMBIENT to argument list
!                * added NAMBIENT to READ statement
!                * changed FORM='UNFORMATTED' to ACCESS='STREAM'
! 12.02.11 (BTD) merge subroutine readE and program callreadE
!                into single program readnf to simplify array allocation
! 12.02.12 (BTD) testing, debugging, and improvements
! 12.02.13 (PJF) added VTR capabilities
! 12.02.13 (BTD) create options IVTR=1 and 2 for |E| or |E^2|
!                with new output file name = Esquared
! 12.02.14 (PJF) added Real(CXE0_TF) to VTR output
! 12.12.25 (BTD) v3 and DDSCAT 7.3.0
!                * modify to read modified structure of data file
!                * data file includes integer NRFLDB=0 or 1 that
!                  indicates if B_inc and B_sca are in data file
!                * data file includes composition information and
!                  complex dielectric function(s) for possible use
!                  in analyzing results
!                * if NRFLDB=0, output E along a line as before
!                * if NRFLDB=1, output both E and B along line
! 13.01.01 (BTD) Modified comments to be in accord with latest
!                modifications to READNF
! 13.03.20 (BTD) Program DDPOSTPROCESS.f90 created to
!                use subroutine READNF to read data files
!                and allocate memory as needed.
!                Arrays are now allocated with 3 spatial indices
!                e.g., CXPOL(1:NX,1:NY,1:NZ,1:3)
! 13.03.21 (BTD,PJF) bug fix
! 13.03.23 (BTD) Add calculation of time-averaged Poynting vector S
! 13.03.24 (BTD) Add output of Poynting vector to VTR
! 13.03.25 (BTD) v2 -- no significant change as of yet
! 13.05.05 (BTD) * corrected typo reported by Choliy Vasyl
! 13.05.19 (BTD) * calculate S_INC always, not just when NRFLDB=1
!                * modified formatting of output to ddpostprocess.out
! end history
! Copyright (C) 2011,2012,2013
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!=======================================================================
      DATA CFLPAR_DEFAULT/'ddpostprocess.par'/

      CFLPAR=CFLPAR_DEFAULT
      IF(IARGC().EQ.1)THEN
         CALL GETARG(1,COMMAND(1))
         CFLPAR=COMMAND(1)
      ENDIF
      IDVOUT=0
      WRITE(IDVOUT,FMT='(A)')'>DDPOSTPROCESS: using parameter file='
      WRITE(IDVOUT,FMT='(A,A)')'       ',CFLPAR

      OPEN(UNIT=3,FILE=CFLPAR)
!*** diagnostic
!      write(0,*)'ddpostprocess_v2 ckpt 1'
!***

      READ(3,*)CFLENAME
      WRITE(IDVOUT,FMT='(A,A)')'>DDPOSTPROCESS: input data from ',CFLENAME

      READ(3,*)CFLVTR
      WRITE(IDVOUT,FMT='(A,A)')'>DDPOSTPROCESS: VTK output name= ',CFLVTR

      READ(3,*)IVTR

! if IVTR > 0, then create VTR files for subsequent visualization

      READ(3,*)ILINE


      PI=4._WP*ATAN(1._WP)
      WRITE(IDVOUT,FMT='(A,A)')'>DDPOSTPROCESS: now read file = ',CFLENAME

!==============================================================================

      CALL READNF(CFLENAME,IDVOUT,CSTAMP,VERSNUM,NRFLDB,AEFF,DPHYS, &
                  NX,NY,NZ,NAT0,X0,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,   &
                  NAMBIENT,WAVE,AK_TF,CXE0_TF,CXB0_TF,NCOMP)        !

!==============================================================================
! determine E_inc, E_sca, and P at points along defined track

! if NRFLDB=1, also calculate B and time-averaged S along track, where
! S=Poynting vector normalized by Poynting vector of incident plane wave
! Let E=Re[Ec*e^(-iwt)]=E1cos(wt)+E2sin(wt)  where Ec=E1+iE2
!     B=Re[Bc*e^(-iwt)]=B1cos(wt)+B2sin(wt)  where Bc=B1+iB2
!     S=(1/4pi)*ExB
!    <S>=(1/4pi)*(1/2)(E1xB1+E2xB2)
!       =(1/8pi)*( E1xB1 + E2xB2 )
!       =(1/8pi)*Re( Ec x conjg(Bc) )
! we will omit the (1/8pi) because we will always normalize by incident S

! calculate S_inc=8*pi*|<S>| for incident wave

!*** diagnostic
!      write(0,*)'ddpostprocess_v2 ckpt 2, NRFLDB=',NRFLDB
!      write(0,fmt='(a,1p6e10.3)')'  cxe0_tf=',cxe0_tf
!      write(0,fmt='(a,1p6e10.3)')'  cxb0_tf=',cxb0_tf
!***
      S_INC(1)=REAL(CXE0_TF(2)*CONJG(CXB0_TF(3))- &
                    CXE0_TF(3)*CONJG(CXB0_TF(2)))
      S_INC(2)=REAL(CXE0_TF(3)*CONJG(CXB0_TF(1))- &
                    CXE0_TF(1)*CONJG(CXB0_TF(3)))
      S_INC(3)=REAL(CXE0_TF(1)*CONJG(CXB0_TF(2))- &
                    CXE0_TF(2)*CONJG(CXB0_TF(1)))
      SNORM=SQRT(S_INC(1)**2+S_INC(2)**2+S_INC(3)**2)
      WRITE(IDVOUT,FMT='(A,1P3E11.3)')                        &
         '>DDPOSTPROCESS: 8*pi*|<S>| for incident wave =',SNORM
!*** 
      IF(NRFLDB==1)THEN
!*** diagnostic
!         write(0,*)'ddpostprocess_v2 ckpt 3'
!         write(0,fmt='(a,1p6e10.3)')'  cxe0_tf=',cxe0_tf
!         write(0,fmt='(a,1p6e10.3)')'  cxb0_tf=',cxb0_tf
!***

! calculate normalized Poynting vector at all points
         ALLOCATE(S(1:NX,1:NY,1:NZ,1:3))
         DO JZ=1,NZ
            DO JY=1,NY
               DO JX=1,NX
                  CXEX=CXEINC(JX,JY,JZ,1)+CXESCA(JX,JY,JZ,1)
                  CXEY=CXEINC(JX,JY,JZ,2)+CXESCA(JX,JY,JZ,2)
                  CXEZ=CXEINC(JX,JY,JZ,3)+CXESCA(JX,JY,JZ,3)
                  CXBX=CXBINC(JX,JY,JZ,1)+CXBSCA(JX,JY,JZ,1)
                  CXBY=CXBINC(JX,JY,JZ,2)+CXBSCA(JX,JY,JZ,2)
                  CXBZ=CXBINC(JX,JY,JZ,3)+CXBSCA(JX,JY,JZ,3)
                  S(JX,JY,JZ,1)=REAL(CXEY*CONJG(CXBZ)-CXEZ*CONJG(CXBY))/SNORM
                  S(JX,JY,JZ,2)=REAL(CXEZ*CONJG(CXBX)-CXEX*CONJG(CXBZ))/SNORM
                  S(JX,JY,JZ,3)=REAL(CXEX*CONJG(CXBY)-CXEY*CONJG(CXBX))/SNORM
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      WRITE(IDVOUT,FMT='(A,I3,A)')'>DDPOSTPROCESS: NCOMP=',NCOMP,' compositions'

      IF(ILINE>0)THEN

         OPEN(UNIT=7,FILE='ddpostprocess.out')

 3000    READ(3,*,END=4000)XA,YA,ZA,XB,YB,ZB,NAB

! check that track is entirely within computational volume

         IF(XA>=XMIN.AND.XB<=XMAX.AND. &
            YA>=YMIN.AND.YB<=YMAX.AND. &
            ZA>=ZMIN.AND.ZB<=ZMAX)THEN
            WRITE(7,FMT='(1PE10.4,A)')                       &
               AEFF,' = a_eff (radius of equal-volume sphere)'
            WRITE(7,FMT='(1PE10.4,A)')DPHYS,' = d = dipole spacing'
            WRITE(7,FMT='(I10,A)')                            &
               NAT0,' = N = number of dipoles in target or TUC'
            WRITE(7,FMT='(3I8,A)')                                   &
               NX,NY,NZ,' = NX,NY,NZ = extent of computational volume'
            WRITE(7,FMT='(1PE10.4,A)')WAVE,' = wavelength in vacuo'
            WRITE(7,FMT='(F10.5,A)')                          &
               NAMBIENT,' = refractive index of ambient medium'
            WRITE(7,FMT='(3F10.5,A)')AK_TF,' = k_{inc,TF} * d'
            WRITE(7,FMT='(5(2F10.5,A,/),2F10.5,A)')                     &
              CXE0_TF(1),                                               &
              ' = Re(E_inc,x) Im(E_inc,x) at x_TF=0,y_TF=0,z_TF=0,t=0', &
              CXE0_TF(2),' = Re(E_inc,y) Im(E_inc,y) "',                &
              CXE0_TF(3),' = Re(E_inc,z) Im(E_inc,z) "',                &
              CXB0_TF(1),' = Re(B_inc,x) Im(B_inc,x) "',                &
              CXB0_TF(2),' = Re(B_inc,y) Im(B_inc,y) "',                &
              CXB0_TF(3),' = Re(B_inc,z) Im(B_inc,z) "'                 !
!*** diagnostic
!            write(0,*)'ddpostprocess_v2 ckpt 4: S_INC=',S_INC
!***
            WRITE(7,FMT='(3F10.5,A,A,/,A,A)')S_INC,                          &
               ' = 2*(4pi/c)*<S_inc> where <S_inc>=time-averaged incident ', &
               'Poynting vector','[Poynting vector S = (Sx,Sy,Sz) = ',       &
               ' (c/4pi)*Re(E)xRe(B) ]'                                      !
            WRITE(7,FMT='(F10.5,A)')SNORM, &
               ' = 2*(4pi/c)*|<S_inc>|'    !
            IF(NRFLDB==0)THEN
               WRITE(7,FMT='(3X,A,7X,A,7X,A,6X,A,3X,A,3X,A,3X,A,3X,A,3X,A)') &
                  'x_TF','y_TF','z_TF','Re(E_x)','Im(E_x)','Re(E_y)',        &
                  'Im(E_y)','Re(E_z)','Im(E_z)'                              !
            ELSE
               WRITE(7,FMT='(33X,A,2X,A,2X,A)')                                &
                '----------------------- E field ---------------------------', &
                '----------------------- B field --------------------------',  &
                '-normalized Poynting vector-'                                 !
               WRITE(7,FMT='(3X,A,7X,A,7X,A,6X,A,11(3X,A),2X,A)')     &
                  'x_TF','y_TF','z_TF','Re(E_x)','Im(E_x)','Re(E_y)', &
                  'Im(E_y)','Re(E_z)','Im(E_z)','Re(B_x)','Im(B_x)',  &
                  'Re(B_y)','Im(B_y)','Re(B_z)','Im(B_z)',            &
                  '-- <(Sx,Sy,Sz)>/|<S_inc>| --'                      !
            ENDIF

            DO JA=0,NAB
               ZETA=REAL(JA)/REAL(NAB)
               XTF(1)=XA+(XB-XA)*ZETA
               XTF(2)=YA+(YB-YA)*ZETA
               XTF(3)=ZA+(ZB-ZA)*ZETA
          
! XTF is assumed to be in physical units
! XTF/DPHYS is in dipole units
! I + X0 = XTF/DPHYS
! I = XTF/DPHYS - X0

               IX1=INT(XTF(1)/DPHYS-X0(1))
               IY1=INT(XTF(2)/DPHYS-X0(2))
               IZ1=INT(XTF(3)/DPHYS-X0(3))
               IF(IX1.LT.1)IX1=1
               IF(IX1.GE.NX)IX1=NX-1
               IF(IY1.LT.1)IY1=1
               IF(IY1.GE.NY)IY1=NY-1
               IF(IZ1.LT.1)IZ1=1
               IF(IZ1.GE.NZ)IZ1=NZ-1
               IX2=IX1+1
               IY2=IY1+1
               IZ2=IZ1+1
               WX=XTF(1)/DPHYS-X0(1)-IX1
               WY=XTF(2)/DPHYS-X0(2)-IY1
               WZ=XTF(3)/DPHYS-X0(3)-IZ1
! ------------ determine weights ---------- xyz zyx
               W1=(1.-WX)*(1.-WY)*(1.-WZ) ! 000 000
               W2=WX*(1.-WY)*(1.-WZ)      ! 100 001
               W3=(1.-WX)*WY*(1.-WZ)      ! 010 010
               W4=WX*WY*(1.-WZ)           ! 110 011
               W5=(1.-WX)*(1.-WY)*WZ      ! 001 100
               W6=WX*(1.-WY)*WZ           ! 101 101
               W7=(1.-WX)*WY*WZ           ! 011 110
               W8=WX*WY*WZ                ! 111 111
! -------- evaluate weighted averages -----------------------------------
               DO K=1,3
                  CXE_INC(K)=W1*CXEINC(IX1,IY1,IZ1,K)+ &
                             W2*CXEINC(IX2,IY1,IZ1,K)+ &
                             W3*CXEINC(IX1,IY2,IZ1,K)+ &
                             W4*CXEINC(IX2,IY2,IZ1,K)+ &
                             W5*CXEINC(IX1,IY1,IZ2,K)+ &
                             W6*CXEINC(IX2,IY1,IZ2,K)+ &
                             W7*CXEINC(IX1,IY2,IZ2,K)+ &
                             W8*CXEINC(IX2,IY2,IZ2,K)  !

                  CXE_SCA(K)=W1*CXESCA(IX1,IY1,IZ1,K)+ &
                             W2*CXESCA(IX2,IY1,IZ1,K)+ &
                             W3*CXESCA(IX1,IY2,IZ1,K)+ &
                             W4*CXESCA(IX2,IY2,IZ1,K)+ &
                             W5*CXESCA(IX1,IY1,IZ2,K)+ &
                             W6*CXESCA(IX2,IY1,IZ2,K)+ &
                             W7*CXESCA(IX1,IY2,IZ2,K)+ &
                             W8*CXESCA(IX2,IY2,IZ2,K)  !

                  CXP(K)=W1*CXPOL(IX1,IY1,IZ1,K)+W2*CXPOL(IX2,IY1,IZ1,K)+ &
                         W3*CXPOL(IX1,IY2,IZ1,K)+W4*CXPOL(IX2,IY2,IZ1,K)+ &
                         W5*CXPOL(IX1,IY1,IZ2,K)+W6*CXPOL(IX2,IY1,IZ2,K)+ &
                         W7*CXPOL(IX1,IY2,IZ2,K)+W8*CXPOL(IX2,IY2,IZ2,K)

                  CXE(K)=CXE_INC(K)+CXE_SCA(K)
               ENDDO

               IF(NRFLDB==0)THEN

! write out total macroscopic E field at point on track
  
                  WRITE(7,FMT='(1PE10.3,1P2E11.3,0P6F10.5)')XTF,CXE

               ELSEIF(NRFLDB==1)THEN
                  DO K=1,3
                     CXB_INC(K)=W1*CXBINC(IX1,IY1,IZ1,K)+ &
                                W2*CXBINC(IX2,IY1,IZ1,K)+ &
                                W3*CXBINC(IX1,IY2,IZ1,K)+ &
                                W4*CXBINC(IX2,IY2,IZ1,K)+ &
                                W5*CXBINC(IX1,IY1,IZ2,K)+ &
                                W6*CXBINC(IX2,IY1,IZ2,K)+ &
                                W7*CXBINC(IX1,IY2,IZ2,K)+ &
                                W8*CXBINC(IX2,IY2,IZ2,K)  !
                     CXB_SCA(K)=W1*CXBSCA(IX1,IY1,IZ1,K)+ &
                                W2*CXBSCA(IX2,IY1,IZ1,K)+ &
                                W3*CXBSCA(IX1,IY2,IZ1,K)+ &
                                W4*CXBSCA(IX2,IY2,IZ1,K)+ &
                                W5*CXBSCA(IX1,IY1,IZ2,K)+ &
                                W6*CXBSCA(IX2,IY1,IZ2,K)+ &
                                W7*CXBSCA(IX1,IY2,IZ2,K)+ &
                                W8*CXBSCA(IX2,IY2,IZ2,K)  !
                     CXB(K)=CXB_INC(K)+CXB_SCA(K)

! calculate time-averaged Poynting vector at each point, normalized by
! magnitude of time-averaged Poynting vector of incident plane wave

                     SVEC(K)=W1*S(IX1,IY1,IZ1,K)+ &
                             W2*S(IX2,IY1,IZ1,K)+ &
                             W3*S(IX1,IY2,IZ1,K)+ &
                             W4*S(IX2,IY2,IZ1,K)+ &
                             W5*S(IX1,IY1,IZ2,K)+ &
                             W6*S(IX2,IY1,IZ2,K)+ &
                             W7*S(IX1,IY2,IZ2,K)+ &
                             W8*S(IX2,IY2,IZ2,K)  !
                  ENDDO

                  WRITE(7,FMT='(1PE10.3,1P2E11.3,0P15F10.5)')XTF,CXE,CXB,SVEC
               ELSE
                  WRITE(0,*)'ddpostprocess fatal error: NRFLDB=',NRFLDB
                  STOP
               ENDIF   ! endif(nrfldb=0 or 1)

            ENDDO
         ELSE
            WRITE(0,FMT='(A,A)')'ddpostprocess fatal error: ',       &
               'requested track extends beyond computational volume'
            WRITE(0,FMT='(A,1P3E11.3)')'XA,YA,ZA=',XA,YA,ZA
            WRITE(0,FMT='(A,1P3E11.3)')'XB,YB,ZB=',XB,YB,ZB
            WRITE(0,FMT='(A,1P3E11.3)')'XMIN,YMIN,ZMIN=',XMIN,YMIN,ZMIN
            WRITE(0,FMT='(A,1P3E11.3)')'XMAX,YMAX,ZMAX=',XMAX,YMAX,ZMAX
            STOP
         ENDIF
         GOTO 3000
 4000    CLOSE(3)
         CLOSE(7)
         WRITE(IDVOUT,FMT='(A)')                             &
            '>DDPOSTPROCESS: completed computing E along tracks'
     ENDIF   ! endif(iline>0)
 
! have completed:
! 1. reading data from file
! 2. evaluating E_inc, E_sca, and P at points along track
! 3. writing track data to ascii output file

!*************************************************************************

! *** Additional code below here to write out additional files for display
!     or analysis.

      IF(IVTR<=0)STOP

! if IVTR > 0: now write VTK file for graphics
! define mesh for graphics

!VTK supplementary arrays 
! mesh x,y,z assuming that we are on rectangular grid)
! variables are dimensioned nx,ny,nz

      ALLOCATE(VTRX(NX),VTRY(NY),VTRZ(NZ))
      IF(IVTR<3)THEN
         ALLOCATE(VTR8(NX,NY,NZ))
      ELSEIF(IVTR==3)THEN
         IF(NRFLDB==1)THEN
            ALLOCATE(VTRVX(NX,NY,NZ))
            ALLOCATE(VTRVY(NX,NY,NZ))
            ALLOCATE(VTRVZ(NZ,NY,NZ))
         ELSE
            WRITE(0,FMT='(A,A)')'Fatal Error: cannot execute option IVTR=3', &
               ' because nearfield B was not calculated (NRFLDB=0)'
            STOP
         ENDIF
      ELSE
         WRITE(0,FMT='(A,I10)')'Fatal Error: unknown option IVTR=',IVTR
      ENDIF
         

!write VTK file for graphics
!define mesh for graphics

      CALL VTR_OPEN_FILE(PREFIX=CFLVTR, FD=FD)

! calculate x,y,z in dimensional units

      DO JX=1,NX
         VTRX(JX)=(JX+X0(1))*DPHYS
      ENDDO
      DO JY=1,NY
         VTRY(JY)=(JY+X0(2))*DPHYS
      ENDDO
      DO JZ=1,NZ
         VTRZ(JZ)=(JZ+X0(3))*DPHYS
      ENDDO

      CALL VTR_WRITE_MESH(FD=fd, X=VTRX, Y=VTRY, Z=VTRZ)

! electric field intensity or energy density

      DO JZ=1,NZ
         DO JY=1,NY
            DO JX=1,NX
               CXEE(1:3)=CXEINC(JX,JY,JZ,1:3)+CXESCA(JX,JY,JZ,1:3)
         
!intrinsic function dot_product(cx,cy)=sum(conjg(cx),cy)
!         vtr8(ix,iy,iz)=sqrt(sum(abs(cxee)**2))

               IF(IVTR==1)THEN
                  VTR8(JX,JY,JZ)=SQRT(DOT_PRODUCT(CXEE(1:3),CXEE(1:3)))
               ELSEIF(IVTR==2)THEN
                  VTR8(JX,JY,JZ)=DOT_PRODUCT(CXEE(1:3),CXEE(1:3))
               ELSEIF(IVTR==3)THEN
                  CXBB(1:3)=CXBINC(JX,JY,JZ,1:3)+CXBSCA(JX,JY,JZ,1:3)
                  VTRVX(JX,JY,JZ)=S(JX,JY,JZ,1)
                  VTRVY(JX,JY,JZ)=S(JX,JY,JZ,2)
                  VTRVZ(JX,JY,JZ)=S(JX,JY,JZ,3)
               ELSE
                  WRITE(IDVOUT,FMT='(A,I6)')                             &
                     '>READNF Fatal error: unrecognized option IVTR=',IVTR
                  STOP
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF(IVTR==1)THEN
         CALL VTR_WRITE_VAR(FD=fd,NAME="Intensity",FIELD=VTR8)
      ELSEIF(IVTR==2)THEN
         CALL VTR_WRITE_VAR(FD=fd,NAME="Esquared",FIELD=VTR8)
      ELSEIF(IVTR==3)THEN
         CALL VTR_WRITE_VAR(FD=fd,NAME="Poynting",VX=VTRVX,VY=VTRVY,VZ=VTRVZ)
      ENDIF

! composition 
! This can be used to display nicely inhomogeneous objects
! ICOMP=0 outside the target

      DO JZ=1,NZ
         DO JY=1,NY
            DO JX=1,NX
! let us output only one component of icomp
               VTR8(JX,JY,JZ)=ICOMP(JX,JY,JZ,1)
            ENDDO
         ENDDO
      ENDDO

      CALL VTR_WRITE_VAR(FD=fd, NAME="Composition", FIELD=VTR8)
      CALL VTR_CLOSE_FILE(FD=fd)
      CALL VTR_COLLECT_FILE(fd)  ! Produces "output.vtr" and "output.pvd" files

!We are done with 3d objects (large files)
!Now VTR output for smaller (3)-vectors defining geometry
!we probably should also output CXE0 (but it is complex, what to do?)
!copy to kind=8 precision for graphics and output

        XVECT(1)=0.
        YVECT(1)=0.
        ZVECT(1)=0.
        CALL VTR_OPEN_FILE(PREFIX='akr',FD=fd)
        CALL VTR_WRITE_MESH(FD=fd,X=XVECT,Y=YVECT,Z=ZVECT)
        U(1,1,1)=AK_TF(1)
        V(1,1,1)=AK_TF(2)
        W(1,1,1)=AK_TF(3)
        CALL VTR_WRITE_VAR(FD=fd,NAME='akr',VX=U,VY=V,VZ=W)
        U(1,1,1)=REAL(CXE0_TF(1))
        V(1,1,1)=REAL(CXE0_TF(2))
        W(1,1,1)=REAL(CXE0_TF(3))
        CALL VTR_WRITE_VAR(FD=fd,NAME='E0',VX=U,VY=V,VZ=W)
        CALL VTR_CLOSE_FILE(FD=fd)
        CALL VTR_COLLECT_FILE(fd)

      WRITE(IDVOUT,FMT='(A)')'>DDPOSTPROCESS: completed writing VTK file'

!========================== Add code below here as desired ===================
! Following arrays are available for use in additional calculations or output:

! CXEINC(1:NX,1:NY,1:NZ,1:3)= E_macro due to incident wave
!                           = [3/(CXEPS+2)]*(E_micro of incident wave)
! CXESCA(1:NX,1:NY,1:NZ,1:3)= E_macro radiated by dipoles
!                           = [3/(CXEPS+2)]*(E_micro radiated by dipoles)
!                             total E = E_inc + E_sca
! CXPOL(1:NX,1:NY,1:NZ,1:3) = P/d^3 at lattice sites
!                           = 0 at unoccupied lattice sites
! CXEPS(1:NCOMP)            = complex dielectric function for
!                             compositions 1-NCOMP
! CXADIA(1:NX,1:NY,1:NZ,1:3)= diagonal elements of A matrix
!                           = d^3/polarizability    at occupied sites
!                           = 1 at unoccupied sites (where ICOMP(jx,jy,jz,k)=0)
! ICOMP(1:NX,1:NY,1:NZ,1:3) = composition identifier at each site
!                           = 0 at unoccupied sites
!                             where refractive index = NAMBIENT
! if original calculation was done with NRFLDB=1, then following
! arrays are also available:
! CXBINC(1:NX,1:NY,1:NZ,1:3)= B due to incident wave
! CXBSCA(1:NX,1:NY,1:NZ,1:3)= B radiated by dipoles
!                             total B = B_inc + B_sca
! S(1:NX,1:NY,1:NZ,1:3)     = time-averaged Poynting vector, normalized
!                             by Poynting flux in incident wave
 
      STOP
    END PROGRAM DDPOSTPROCESS
