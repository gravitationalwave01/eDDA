    SUBROUTINE TARGSPHER(AEFF,PRINAX,BETA,RLMAX,S2,RSEED,   &
                         A1,A2,DX,X0,CFLSHP,CDESCR,         &
                         IOSHP,MXNAT,NAT,IXYZ,ICOMP)
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE

! Parameters:
! note: this must agree with corresponding parameter in FUNCTION RSPHARM
!       below

      INTEGER NMAX
      PARAMETER(NMAX=100)

! Arguments:

      CHARACTER(67) :: CDESCR
      CHARACTER(*) :: CFLSHP
      INTEGER :: IOSHP,MXNAT,NAT
      INTEGER*2 :: ICOMP(MXNAT,3)
      INTEGER :: IXYZ(MXNAT,3)
      REAL(WP) :: AEFF,BETA,PRINAX,RLMAX,RSEED,S2,S2LMAX
      REAL(WP) :: A1(3),A2(3),DX(3),X0(3)

! Local variables

      INTEGER NRES
      PARAMETER(NRES=20)
      CHARACTER CMSGNM*70
      INTEGER ::                              &
         ISEED,JA,JPHI,JTH,JX,JY,JZ,          &
         L,LMAX,LMX1,LMX2,LMY1,LMY2,LMZ1,LMZ2,M
      REAL(WP) ::                                       &
         AR,AI,COSPH,COSTH,F1,FOURPI,ONE,PHI,PHIARG,PI, &
         R,RAD,ROOTTWO,RYZ,RYZ2,                        &
         SIGLM,SINPH,SINTH,SUM,                         &
         TERM,THETA,TWOPI,X,XMAX,XMIN,XX,               &
         Y,YMAX,YMIN,YY,Z,ZERO,ZMAX,ZMIN,ZZ
      REAL(WP) ::             &
         ALPHA(1:3),          &
         PHIS(0:NRES*NMAX,3), &
         RADS(0:NRES*NMAX,3), &
         THETS(0:NRES*NMAX,3)
      COMPLEX(WP) :: CXI
      COMPLEX(WP) :: ALM(1:NMAX,0:NMAX)

! External functions

      REAL(WP) :: GASDEV,RGSPHER
      EXTERNAL GASDEV,RGSPHER

!***********************************************************************
! Routine to construct irregular grain using spherical harmonics
! to define surface.
!
! Input:
!        AEFF  =effective radius/d (d=lattice spacing)
!        PRINAX=0 to use (1,0,0) and (0,1,0) in TF as a1,a2
!              =1 to use principal axes of moment of inertia tensor
!                        for a1 and a2
!        BETA  = power law index
!                if BETA = 0 : read ALM values from file CFLSHP
!                if BETA > 1 : generate ALM values here using
!                              random number generator, with
!
!                <|a_LM|^2> propto [1/(L+1)]*L**(-BETA) for L>1

!        RLMAX = maximum L value to use in Y_LM expansion
!                LMAX=NINT(RLMAX)
!        RSEED = to choose seed for random number generator
!              = 1, 2, 3, 4, ...
!                ISEED=-NINT(RSEED) is passed to gaussian deviate program
!        S2 = expectation value of <s^2>, where < > denotes average
!             over grain surface

!        CFLSHP=name of file containing spherical harmonic amplitudes
!               (used only if BETA = 0)
!        IOSHP=device number for "target.out" file
!             =-1 to suppress printing of "target.out"
!        MXNAT=dimensioning information (max number of atoms)
! Output:
!        A1(1-3)=unit vector (1,0,0) defining target axis 1
!        A2(1-3)=unit vector (0,1,0) defining target axis 2
!        NAT=number of atoms in target
!        IXYZ(1-NAT,1-3)=(x-x0(1))/d,(y-x0(2))/d,(z-x0(3))/d
!                        for atoms of target
!        X0(1-3)=(location/d) in Target Frame corresponding to dipole with
!                IXYZ=(0,0,0).  This will be treated at the origin of physical
!                coordinates in the TF.
!                Here origin is set to be centroid of the target.
!        CDESCR=string describing target (up to 67 chars.)
!

!=======================================================================
!                   Gaussian Sphere Geometry

! The target is defined in polar coordinates by a radius function

!                                Lmax  l
!    R(theta,phi).= const * exp[ sum  sum  a_lm Y_lm(theta,phi) ]
!                                l=1  m=-l
!
!    We assume that the a_lm are gaussian random variables,
!    with
!                   4*pi   
!    <|a_{1M}|^2> = ---- * f1* <s^2>
!                    3     
!
!                    4*pi     (1-f1)*<s^2>     1
!    <|a_{LM}|^2> = ------- * ------------ * ------
!                   (2*L+1)   zeta(beta)-1   L^beta

! After generating the a_lm values according to above prescription for
! L=1,...,NMAX, we calculate the actual <s^2> averaged over angles.
! Because the a_lm have been chosen stochastically, this will not be
! exactly equal to the input S2.
! We then rescale all coefficients a_lm by a fixed factor so that
! the resulting a_lm have angle-averaged <s^2> = S2 (summing over
! L values up to NMAX).
! For reference, we also calculate <s^2> with the sum restricted to
! L values up to LMAX (as is used in the target generation).

! For s to be real it is also required that

!    a_{L,-M} = (-1)^m conjg(a_{L,M})    ---> Im[a_{L0}] = 0

! Because of the above requirement, we do not calculate or store
! the a_lm for m < 0, since all that we require to calculate the
! target shape are the combination

!    a_{L,m}Y_{L,m} + a_{L,-m}Y_{L,-m}

!                  (2L+1)*(L-m)!
!          = sqrt[ ------------- ] * P_Lm * 
!                    4*pi*(L+m)!

!              [ a_Lm e^{i*m*phi} + (-1)^m a_{L,-m}*e^{-i*m*phi) ]
!
! where P_Lm(theta,phi) is the associated Legendre function.
! With above requirement on a_{L,-M} it follows that
!
!                  (2L+1)*(L-m)!
!          = sqrt[ ------------- ] * P_Lm * 2 Re(a_Lm e^{i*m*phi})
!                    4*pi*(L+m)!


! B.T.Draine, Princeton Univ. Obs.  Original version: 02.11.12
!
! History:
! 02.11.12 (BTD): Start to write this routine
! 02.11.13 (BTD): further modification
! 03.05.21 (BTD): bug fix (set NAT=0)
!                 bug fix to suppress cross section info when IOSHP < 1
! 03.11.06 (BTD): modify to use random number generator to choose
!                 coefficients a_lm
! 03.11.08 (BTD): Using TARPHARM as template, write routine to
!                 generate gaussian random spheres using
!                 r prop exp(s), with s = sum a_lm Y_lm
! 04.03.19 (BTD): revert to previous version of PRINAXIS
! 04.04.01 (BTD): add DX to argument list
! 04.04.05 (BTD): Replace WRITE(0, with CALL WRIMSG
! 04.05.23 (BTD): Obtain eigenvalues ALPHA from PRINAXIS
!                 modified argument list in call to PRINAXIS to
!                 conform to modification in PRINAXIS
! 07.08.06 (BTD): Converted to f90
! 07.09.11 (BTD): Changed IXYZ from INTEGER*2 to INTEGER
! 08.05.12 (BTD): Added X0 to argument list, and add code to calculate X0
! 08.08.29 (BTD): Modified format 9020
! end history
!
! Copyright (C) 2002,2003,2004,2007,2008 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!-----------------------------------------------------------------------

      ZERO=0._WP
      ONE=1._WP
      CXI=(0._WP,1._WP)
      PI=4._WP*ATAN(ONE)
      TWOPI=2._WP*PI
      FOURPI=4._WP*PI
      ROOTTWO=SQRT(2._WP)

! set Lmax

      LMAX=NINT(RLMAX)

      IF(BETA==ZERO)THEN
         WRITE(CMSGNM,'(A,A)')'going to read alm from file = ',CFLSHP
         CALL WRIMSG('TARGSPHER',CMSGNM)
         OPEN(UNIT=IOSHP,FILE=CFLSHP)
         READ(IOSHP,*)S2LMAX,S2
         READ(IOSHP,*)
         READ(IOSHP,*)
         READ(IOSHP,*)
         READ(IOSHP,*)LMAX
         WRITE(CMSGNM,'(A,I3)')'Lmax=',LMAX
         CALL WRIMSG('TARGSPHER',CMSGNM)
         READ(IOSHP,*)

! read ALM values for nonnegative m
! LMAX = 1 : 2    (m=0 and 1)
! LMAX = 2 : 2+3 = 5
! LMAX = 3 : 2+3+4 = 9
! general: (LMAX^2+3*LMAX)/2

         JTH=(LMAX**2+3*LMAX)/2
         DO JA=1,JTH
            READ(IOSHP,*)L,M,AR,AI
            ALM(L,M)=AR+CXI*AI
         ENDDO
         WRITE(CMSGNM,'(A,I4,A,I3,I3)')                  &
               'Read ',JTH,' a_lm values through L,M=',L,M
         CALL WRIMSG('TARGSPHER',CMSGNM)
         CLOSE(IOSHP)

! calculate S2

         SUM=ZERO
         DO L=1,LMAX
            SUM=SUM+ABS(ALM(L,0))**2
            DO M=1,L
                SUM=SUM+2._WP*ABS(ALM(L,M))**2
            ENDDO
         ENDDO
         WRITE(CMSGNM,FMT='(A,1PE10.3)')'input S2LMAX=',s2lmax
         CALL WRIMSG('TARGSPHER',CMSGNM)
         S2LMAX=SUM/FOURPI
         WRITE(CMSGNM,FMT='(A,1PE10.3)')'S2LMAX=',S2LMAX
         CALL WRIMSG('TARGSPHER',CMSGNM)

      ELSEIF(BETA>ONE)THEN

! set value of F1 = fraction of <s^2> contributed by L=1
! note that to leading order the L=1 terms produce merely a
! translation with no deformation, so we normally suppress the
! L=1 terms by setting F1=0.0

         F1=0.0_WP

! Assign a_lm values using gaussian random deviates
! set seed

         ISEED=-NINT(RSEED)

! sanity check:

         IF(LMAX>NMAX)THEN
            WRITE(CMSGNM,FMT='(A,I4,A,I3)')                      &
               'Fatal problem: called tarspharm with LMAX=',LMAX, &
               ' > NMAX=',NMAX
            CALL WRIMSG('TARGSPHER',CMSGNM)
            STOP
         ENDIF

! Calculate random numbers and store in appropriate elements of ALM
! Calculate them in appropriate order so that changes in NMAX will 
! leave lower L values unaffected.  
! Specifying RSEED, F1, and BETA therefore serves to define a shape, 
! with adjustment to LMAX allowing one to see important to the shape
! of high L components.

! initialize GASDEV

         AR=GASDEV(ISEED)

         DO L=1,NMAX
            ALM(L,0)=GASDEV(ISEED)
            DO M=1,L
               AR=GASDEV(ISEED)/ROOTTWO
               AI=GASDEV(ISEED)/ROOTTWO
               ALM(L,M)=AR+CXI*AI
            ENDDO
         ENDDO

! ALM(L,M) now contains random complex numbers a_lm , with
!    <a_lm> = 0 
!    <|a_lm|^2> = <[Re(a_lm)]^2> + <[Im(a_lm)]^2> = 1
!    Im[a_lm(L,0)] = 0
!    <[Re(a_lm)]^2> = <[Im(a_lm)]^2> = 1/2 for M > 0

! Now multiply the A(L,M) random deviates by
! appropriate factors to obtain intended <|a_lm|^2> values

! First treal L=1
! SIGLM = <|a_lm|^2>^{1/2} 

         SIGLM=SQRT(F1*S2*FOURPI/3.)

         ALM(1,0)=SIGLM*ALM(1,0)
         ALM(1,1)=SIGLM*ALM(1,1)

! set coefficients for L > 1:
! estimate SUM = [Riemann zeta function - 1]

         SUM=ZERO
         DO L=NMAX,2,-1
            SUM=SUM+ONE/REAL(L)**BETA
         ENDDO
         TERM=(ONE-F1)*S2*FOURPI/SUM

! Note: we calculate coefficients up to NMAX
! We then use values only up to LMAX to generate shape, but we
! want to have correct amount of power in high order coefficients,
! even if they are later suppressed.

         DO L=2,NMAX
            SIGLM=SQRT(TERM/(REAL(2*L+1)*REAL(L)**BETA))
            DO M=0,L
               ALM(L,M)=SIGLM*ALM(L,M)
            ENDDO
         ENDDO

! compute actual <s^2> for this realization

         SUM=ZERO
         DO L=1,NMAX
            SUM=SUM+ABS(ALM(L,0))**2
            DO M=1,L
               SUM=SUM+2._WP*ABS(ALM(L,M))**2
            ENDDO
         ENDDO
         SUM=SUM/FOURPI

! At this point, we rescale the a_lm values so that we will have
! <s^2> = S2 exactly including all terms up to NMAX

         TERM=SQRT(S2/SUM)
         DO L=1,NMAX
            DO M=0,L
               ALM(L,M)=TERM*ALM(L,M)
            ENDDO
         ENDDO

! Now calculate S2LMAX = <s^2> contributed by terms up to LMAX

         SUM=ZERO
         DO L=1,LMAX
            SUM=SUM+ABS(ALM(L,0))**2
            DO M=1,L
               SUM=SUM+2._WP*ABS(ALM(L,M))**2
            ENDDO
         ENDDO
         S2LMAX=SUM/FOURPI
      ELSE
         CALL ERRMSG('FATAL','TARGSPHER',' Invalid BETA')
      ENDIF

!=======================================================================
! ALM values are now defined.  
! Proceed to generate target array, using a_lm values up to Lmax

      DO JA=1,3
         A1(JA)=0.
         A2(JA)=0.
      ENDDO
      A1(1)=ONE
      A2(2)=ONE

! if IOSHP > 0 , calculate three cross-sectional slices

! polar slice with phi=0

      IF(IOSHP>0)THEN
         PHI=ZERO
         DO JTH=0,NRES*LMAX
            THETA=TWOPI*REAL(JTH)/(REAL(NRES*LMAX))
            COSTH=COS(THETA)
            THETS(JTH,1)=THETA
            PHIS(JTH,1)=PHI
            PHIARG=PHI
            IF(THETA>PI)PHIARG=PHI+PI
            RADS(JTH,1)=RGSPHER(COSTH,PHIARG,LMAX,NMAX,S2,ALM)
         ENDDO

! polar slice with phi=pi/2

         PHI=PI/2._WP
         DO JTH=0,NRES*LMAX
            THETA=TWOPI*REAL(JTH)/(REAL(NRES*LMAX))
            COSTH=COS(THETA)
            THETS(JTH,2)=THETA
            PHIS(JTH,2)=PHI
            PHIARG=PHI
            IF(THETA>PI)PHIARG=PHI+PI
            RADS(JTH,2)=RGSPHER(COSTH,PHIARG,LMAX,NMAX,S2,ALM)
         ENDDO

! equatorial slice

         THETA=PI/2._WP
         COSTH=ZERO
         DO JPHI=0,NRES*LMAX
            PHI=TWOPI*REAL(JPHI)/(REAL(NRES*LMAX))
            THETS(JPHI,3)=THETA
            COSTH=COS(THETA)
            PHIS(JPHI,3)=PHI
            RADS(JPHI,3)=RGSPHER(COSTH,PHI,LMAX,NMAX,S2,ALM)
         ENDDO
         OPEN(UNIT=IOSHP,FILE='spharm.out')
         WRITE(IOSHP,7600)LMAX
         DO JTH=0,NRES*LMAX
            WRITE(IOSHP,7700)(THETS(JTH,L),PHIS(JTH,L),RADS(JTH,L),L=1,3)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF

! sample surface to determine total extent of target in x,y,z directions

      XMIN=ZERO
      XMAX=ZERO
      YMIN=ZERO
      YMAX=ZERO
      ZMIN=ZERO
      ZMAX=ZERO

      DO JTH=1,NRES*LMAX
         THETA=PI*(REAL(JTH)-0.5_WP)/REAL(NRES*LMAX)
         COSTH=COS(THETA)
         SINTH=SIN(THETA)

         DO JPHI=1,NRES*LMAX
            PHI=TWOPI*REAL(JPHI)/REAL(NRES*LMAX)
            RAD=AEFF*RGSPHER(COSTH,PHI,LMAX,NMAX,S2,ALM)
            XX=COSTH*RAD
            YY=SINTH*COS(PHI)*RAD
            ZZ=SINTH*SIN(PHI)*RAD
            IF(XX<XMIN)XMIN=XX
            IF(XX>XMAX)XMAX=XX
            IF(YY<YMIN)YMIN=YY
            IF(YY>YMAX)YMAX=YY
            IF(ZZ<ZMIN)ZMIN=ZZ
            IF(ZZ>ZMAX)ZMAX=ZZ
         ENDDO
      ENDDO

      LMX1=-INT(-XMIN)-1
      LMY1=-INT(-YMIN)-1
      LMZ1=-INT(-ZMIN)-1
      LMX2=INT(XMAX)  +1
      LMY2=INT(YMAX)  +1
      LMZ2=INT(ZMAX)  +1

! construct list of occupied sites

      NAT=0
      DO JZ=LMZ1,LMZ2
         Z=REAL(JZ)
         DO JY=LMY1,LMY2
            Y=REAL(JY)
            RYZ2=Y**2+Z**2
            RYZ=SQRT(RYZ2)
            IF(RYZ2>ZERO)THEN
               COSPH=Y/RYZ
               SINPH=Z/RYZ
            ELSE
               COSPH=ONE
               SINPH=ZERO
            ENDIF
            DO JX=LMX1,LMX2
               X=REAL(JX)
               R=SQRT(RYZ2+X**2)
               IF(R>0.)THEN
                  COSTH=X/R
               ELSE
                  COSTH=ONE
               ENDIF
               IF(COSPH>=0.)THEN
                  PHI=ASIN(SINPH)
               ELSEIF(COSPH<0)THEN
                  PHI=PI-ASIN(SINPH)
               ENDIF

               RAD=AEFF*RGSPHER(COSTH,PHI,LMAX,NMAX,S2,ALM)
               IF(RAD>=R)THEN

! Site is occupied:

                  NAT=NAT+1
                  IF(NAT>MXNAT)THEN
!***
!                     write(0,*)'nat=',nat
!                     write(0,*)'mxnat=',mxnat
!***
                     CALL ERRMSG('FATAL','TARGSPHER',' NAT > MXNAT ')
                  ENDIF
                  IXYZ(NAT,1)=JX
                  IXYZ(NAT,2)=JY
                  IXYZ(NAT,3)=JZ
               ENDIF
            ENDDO
         ENDDO
      ENDDO

! Homogeneous target:

      DO JA=1,NAT
         DO JX=1,3
            ICOMP(JA,JX)=1
         ENDDO
      ENDDO

! Specify target axes A1 and A2
! If PRINAX=0, then
!     A1=(1,0,0) in target frame
!     A2=(0,1,0) in target frame
! If PRINAX=1., then
!     A1,A2 are principal axes of largest, second largest moment
!     of inertia

      IF(PRINAX<=ZERO)THEN
         DO JX=1,3
            A1(JX)=ZERO
            A2(JX)=ZERO
         ENDDO
         A1(1)=ONE
         A2(2)=ONE
      ELSE

         CALL PRINAXIS(MXNAT,NAT,ICOMP,IXYZ,DX,A1,A2,ALPHA)

      ENDIF

! Set X0 so that origin is at the centroid

      DO JY=1,3
         X0(JY)=0._WP
         DO JX=1,NAT
            X0(JY)=X0(JY)+REAL(IXYZ(JX,JY))
         ENDDO
         X0(JY)=-X0(JY)/REAL(NAT)
      ENDDO

!-----------------------------------------------------------------------
! Write target description into string CDESCR
! Note: string CDESCR will be printed by subr. TARGET

      WRITE(CDESCR,FMT='(A,I7,A)')                         &
            ' gaussian sphere target of NAT=',NAT,' dipoles'

!-----------------------------------------------------------------------

! Here print any additional information which is desired by using
! subroutine WRIMSG


      WRITE(CMSGNM,FMT='(A,I7,A)')' Gaussian sphere: NAT=',NAT,' dipoles'

      CALL WRIMSG('TARGSPHER',CMSGNM)

      IF(IOSHP>0)THEN

! if this is a new shape, write the ALM value out to a file

         IF(BETA>0)THEN

            OPEN(UNIT=IOSHP,FILE='targspher.alm')
            WRITE(IOSHP,9100)S2LMAX,S2,F1,BETA,LMAX,NINT(RSEED),    &
                             NAT,AEFF,(0.75*NAT/PI)**(1./3.),ALPHA, &
                             A1,A2,LMAX
            DO L=1,LMAX
               DO M=0,L
                  WRITE(IOSHP,9110)L,M,ALM(L,M)
               ENDDO
            ENDDO
            CLOSE(IOSHP)

         ENDIF

         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')

!*** 3 lines of general description allowed:

         WRITE(IOSHP,9020)S2LMAX,S2,F1,BETA,LMAX,NINT(RSEED),           &
               NAT,AEFF,(0.75_WP*NAT/PI)**(1._WP/3._WP),ALPHA,A1,A2,DX,X0

!***

         DO JX=1,NAT
            WRITE(IOSHP,FMT=9030)JX,IXYZ(JX,1),IXYZ(JX,2),IXYZ(JX,3)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF
      RETURN
 7600 FORMAT(                                                          &
        'three cross sectional slices for spherical harmonic grain',/, &
        I3,' = LMAX (highest order for Y_lm)',/,                       &
        '------ slice 1 ------ ------ slice 2 ------ ',                &
        '------ slice 3 ------',/,                                     &
        ' theta  phi   r/a_eff  theta  phi   r/a_eff  ',               &
        'theta  phi   r/a_eff')
 7700 FORMAT(F6.4,F7.4,F8.4,2(2F7.4,F8.4))
 9020 FORMAT('> TARGSPHER:',F7.5,F7.4,F6.3,F5.2,I3,I4,            &
         ' = <s^2>_Lmax, <s^2>, f_1, beta, Lmax, ISEED',          &
         ' for gaussian sphere',/,                                &
         I10,5F8.4,' = NAT, AEFF(in) a_eff/d, alpha_1-3',/,        &
         3F10.6,' = A_1 vector',/,                                 &
         3F10.6,' = A_2 vector',/,                                 &
         3F10.6,' = lattice spacings (d_x,d_y,d_z)/d',/,           &
         3F10.5,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ', &
               'for dipole 0 0 0',/,                              &
         '     JA  IX  IY  IZ ICOMP(x,y,z)')
 9030 FORMAT(I7,3I4,' 1 1 1')
 9100 FORMAT(                                                          &
        F7.5,F7.4,F6.3,F5.2,I3,I4,' = <s^2>_lmax, <s^2>, f_1, beta, ', &
        'Lmax, ISEED for gaussian sphere',/,                           &
        I7,2F8.4,3F8.4,' = NAT, a_eff/d, alpha_1-3',/,                 &
        3F9.4,' = A_1 vector',/,                                       &
        3F9.4,' = A_2 vector',/,                                       &
        I4,' = Lmax',/,                                                &
        ' L   M   Re(a_LM)  Im(a_LM)')
 9110 FORMAT(I3,I4,2F10.6)
    END SUBROUTINE TARGSPHER

    FUNCTION RGSPHER(COSTH,PHI,LMAX,NMAX,S2,ALM)
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE

! Arguments:

      INTEGER :: LMAX,NMAX
      REAL(WP) :: COSTH,PHI,RGSPHER,S2
      COMPLEX(WP) :: ALM(1:NMAX,0:NMAX)

! Local variables:

      INTEGER :: L,M
      REAL(WP) :: FAC,FOURPI,SIGN,SUM,Y_L0
      COMPLEX(WP) :: CXI,CXTERM

! External functions

      REAL(WP) :: P_LM
      EXTERNAL P_LM

!-----------------------------------------------------------------------
! FUNCTION RGSPHER
! Given:
!    COSTH = cos(theta)
!    PHI   = phi (radians)
!    LMAX = maximum order of spherical harmonics
!    S2   = <s^2>
!    ALM  = complex array of spherical harmonic amplitudes a_lm
!           L=1,...,LMAX , M=-0,..,L
!           Note that we do not require the ALM for m<0, since these
!           are obtained from the values for m>0 using
!           a_{l,-m} = (-1)^m conjg(a_{l,m})

! Returns
!    RGSPHER = distance from "center" to surface for target with
!              a_eff=1.
!                    1             Lmax  l
!            = ------------ * exp[ sum  sum a_lm Y_lm ]
!              sqrt(1+<s^2>)       l=1  m=-l
!
! Requires:
! External function P_LM(L,M,X) to return associated Legendre
!                               polynomial
!                               at the moment this function resides
!                               in file tarspharm.f
!
! Note: we assume that a_{l,-m} = (-1)^m conjg(a_{l,m})
!       so that sum over spherical harmonics gives real function
!
! History
! 02.11.12 (BTD) first written
! 03.11.08 (BTD) gaussian sphere version first created
! end history
!
! Copyright (C) 2002,2003 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!-----------------------------------------------------------------------
      FOURPI=16._WP*ATAN(1._WP)
      CXI=(0._WP,1._WP)

      SUM=0._WP
     
      DO L=1,LMAX
         FAC=SQRT(REAL(2*L+1)/FOURPI)
         Y_L0=FAC*P_LM(L,0,COSTH)

         SUM=SUM+REAL(ALM(1,0))*Y_L0
         SIGN=1._WP
         DO M=1,L
            FAC=FAC/SQRT(REAL((L+M)*(L+1-M)))
            SIGN=-SIGN
            CXTERM=ALM(L,M)*EXP(CXI*REAL(M)*PHI)
            SUM=SUM+FAC*P_LM(L,M,COSTH)*2._WP*REAL(CXTERM)
         ENDDO
      ENDDO

      RGSPHER=EXP(SUM)/SQRT(1._WP+S2)

!*** diagnostic
!      write(17,*)'cos(theta),phi=',costh,phi,' rgspher=',rgspher
!***
      RETURN
      END
