    SUBROUTINE DIELEC(WAVE,IDVOUT,CFLEPS,CXEPS,MXCOMP,MXWAVT,NCOMP,E1A,E2A,WVA)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

! Arguments:

      INTEGER :: MXCOMP,MXWAVT
      CHARACTER(60) :: &
         CFLEPS(MXCOMP)
      COMPLEX(WP) :: &
         CXEPS(MXCOMP)
      INTEGER :: IDVOUT,NCOMP
      REAL(WP) :: WAVE
      REAL(WP) ::     &
         E1A(MXWAVT), &
         E2A(MXWAVT), &
         WVA(MXWAVT)

! Local Variables:

      LOGICAL :: INIT
      CHARACTER(70) :: CMSGNM
      CHARACTER(80) :: CDESCR
      COMPLEX(WP) :: CXI
      INTEGER :: I,J,JJ,NWAVT,IWV,IREN,IIMN,IEPS1,IEPS2,ICOL
      REAL(WP) :: DUM1,DUM2,E1,E2,TEMP
      REAL(WP) :: & 
         XIN(10)
      EXTERNAL ERRMSG
      SAVE CXI,IEPS1,INIT,NWAVT
!**********************************************************************
! Given:
!       WAVE = wavelength (micron)
!       IDVOUT=output unit number
!       CFLEPS(1-NCOMP) = names of files containing dielectric data
!       MXCOMP,MXWAVT = dimensioning information
!       NCOMP = number of components
!       NWAV = number of wavelengths
!       E1A, E2A, WVA = scratch arrays
! Returns:
!       CXEPS(1-NCOMP) = dielectric constant for materials
!                        1-NCOMP
! NOTE:
!       It is assumed that file(s) containing
!       the table(s) will have following format:
! line1 = description (CHARACTER*80)
! line2 = IWV, IREN, IIMN, IEPS1, IEPS2
!         where IWV = column in which wavelengths are tabulated
!               IREN = col. in which Re(N) is tabulated (0 if not)
!               IIMN = col. in which Im(N) is tabulated (0 if not)
!               IEPS1= col. in which Re(EPS) is tabulated (0 if not)
!               IEPS2= col. in which Im(EPS) is tabulated (0 if not)
! line3 = header line (will not be read)
! line4... = data, according to columns specified on line2
! wavelengths must be monotonic, either increasing or decreasing

! B.T.Draine, Princeton Univ. Obs.
! History:
! 90.12.04 (BTD): Rewritten to allow H2OICE and H2OLIQ options.
! 90.12.21 (BTD): Correct handling of 'TABLES' option when more than
!                 one wavelength is considered.
! 91.05.02 (BTD): Added IDVOUT to argument list for subr. INTERP
! 91.09.12 (BTD): Moved declaration of MXCOMP and MWAVT ahead of
!                 declaration of CFLEPS and CXEPS
! 93.12.06 (BTD): Added IEPS1,INIT, and NWAVT to SAVE statement
!                 (without this did not run properly on SGI Indigo)
! 96.12.16 (BTD): Corrected temperature to T=250K in output statement
!                 for H2OICE option
! 98.12.21 (BTD): changed dimension of CFLPAR from CHARACTER*40 to
!                 CHARACTER*60 to allow longer file names
!                 (also changed in reapar.f and DDSCAT.f)
! 04.10.14 (BTD): added check to ensure that user does not specify
!                 refractive index = 1 (this would cause division
!                 by zero elsewhere in code)
! 07.07.30 (PJF): Converted to f90
! 07.07.31 (BTD): Removed calls to REFWAT and REFICE
! 07.10.28 (BTD): Eliminated CDIEL -- reading from tables is standard
!                 Added output line just before table read as clue in
!                 case of failure during table read.
! 09.10.19 (BTD): ver7.0.8
!                 Added some more output to report begin/end reading file
! end history
! Copyright (C) 1993,1996,1998,2004,2007,2009 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
      DATA CXI/(0._WP,1._WP)/,INIT/.TRUE./

      IF(INIT.OR.NCOMP>1)THEN
         DO J=1,NCOMP
            WRITE(CMSGNM,FMT='(A)')'about to read file ='
            CALL WRIMSG('DIELEC',CMSGNM)
            CALL WRIMSG('DIELEC',CFLEPS(J))
            OPEN(UNIT=3,FILE=CFLEPS(J),STATUS='OLD')

! Read header line:

            READ(3,9000)CDESCR
            CALL WRIMSG('DIELEC',CDESCR)
!            WRITE(IDVOUT,9000)CDESCR

! Read line specifying columns for wavelength,Re(n),Im(n),Re(eps),Im(eps
!            write(0,*)'dielec ckpt 1, j=',j
            READ(3,*)IWV,IREN,IIMN,IEPS1,IEPS2
!            write(0,*)'dielec ckpt 1.1, iwv,iren,iimn,ieps1,ieps2=', &
!                      iwv,iren,iimn,ieps1,ieps2
            ICOL=IWV
            IF(IREN>ICOL)ICOL=IREN
            IF(IIMN>ICOL)ICOL=IIMN
            IF(IEPS1>ICOL)ICOL=IEPS1
            IF(IEPS2>ICOL)ICOL=IEPS2

! Skip header line:

            READ(3,*)
            DO I=1,MXWAVT
               READ(3,*,END=600)(XIN(JJ),JJ=1,ICOL)
               WVA(I)=XIN(IWV)
               IF(IEPS1>0)THEN
                  E1A(I)=XIN(IEPS1)
                  E2A(I)=XIN(IEPS2)
               ELSE
                  E1A(I)=XIN(IREN)
                  E2A(I)=XIN(IIMN)
               ENDIF
               NWAVT=I
            ENDDO

! Check whether there is unread data remaining in file

            READ(3,*,END=600)(XIN(JJ),JJ=1,ICOL)

! If this point is reached, apparently unread data remains, so
! issue warning:

            CALL ERRMSG('WARNING','DIELEC', &
                        'parameter MXWAVT not large enough to read full dielec file')
600         CLOSE(3)

! 091019 BTD added
!            write(0,*)' dielec ckpt 2, j=',j,' cfleps(j)=',cfleps(j)
! ------
            WRITE(CMSGNM,FMT='(A,A)')' completed reading file ='
            CALL WRIMSG('DIELEC',CMSGNM)
            CALL WRIMSG('DIELEC',CFLEPS(J))

!*** Have completed reading in table for this composition.
!    Now interpolate

            CALL INTERP(WAVE,E1,E2,WVA,E1A,E2A,IDVOUT,MXWAVT,NWAVT)
            IF(IEPS1>0)THEN
               CXEPS(J)=E1+CXI*E2
            ELSE
               CXEPS(J)=(E1**2-E2**2)+2._WP*CXI*E1*E2
            ENDIF

! check that user has not specified refractive index = 1
! issue fatal warning if this occurs

            IF(CXEPS(J)==1._WP)CALL ERRMSG('FATAL','DIELEC', &
               'Refractive index = 1 is not allowed: should leave unoccupied')

         ENDDO
         INIT=.FALSE.
      ELSE

!*** Perform this only if NCOMP=1 and previously initialized:

         CALL INTERP(WAVE,E1,E2,WVA,E1A,E2A,IDVOUT,MXWAVT,NWAVT)
         IF(IEPS1>0)THEN
            CXEPS(1)=E1+CXI*E2
         ELSE
            CXEPS(1)=(E1**2-E2**2)+2._WP*CXI*E1*E2
         ENDIF
      ENDIF
      RETURN
9000  FORMAT (A80)
    END SUBROUTINE DIELEC
