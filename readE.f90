    SUBROUTINE READE(CFLENAME,XTF,AEFF,NAMBIENT,WAVE,DPHYS,NAT0,NX,NY,NZ, &
                     CXE_INC,CXE_SCA,CXP)
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE

! arguments

      CHARACTER*60 :: CFLENAME
      INTEGER NAT0,NX,NY,NZ
      REAL(WP) :: AEFF,DPHYS,NAMBIENT,WAVE
      REAL(WP) ::  &
         XTF(3)
      COMPLEX(WP) :: &
         CXE_INC(3), &
         CXE_SCA(3), &
         CXP(3)

! local variables

      LOGICAL :: INIT
      CHARACTER*60 :: CFLENAME_OLD

      INTEGER :: IDVOUT,IX1,IY1,IZ1,                               &
         J,J1,J1X,J1Y,J1Z,J2,J2X,J2Y,J2Z,J3X,J3Y,J3Z,J4X,J4Y,J4Z,  &
         J5X,J5Y,J5Z,J6X,J6Y,J6Z,J7X,J7Y,J7Z,J8X,J8Y,J8Z,JX,JY,JZ, &
         NAT3,NXY,NXYZ

      INTEGER*2,        &
         ALLOCATABLE :: &
         ICOMP(:)

      REAL(WP) :: E2,EINC2,PI,SUMERR2,TINY,W1,W2,W3,W4,W5,W6,W7,W8,WX,WY,WZ, &
         XMAX,XMIN,YMAX,YMIN,ZMAX,ZMIN
      REAL(WP) :: X0(3)

      COMPLEX(WP) :: CXERR
      COMPLEX(WP),      &
         ALLOCATABLE :: &
         CXADIA(:),     &
         CXPOL(:),      &
         CXEINC(:),     &
         CXESCA(:)

!=======================================================================
! subroutine READE
! given
!    CFLENAME = name of file with precomputed E and P
!    XTF(3)   = (x_TF,y_TF,z_TF) = coordinates in Target Frame (phys. units)
! returns
!    AEFF        = effective radius of target (phys. units)
!    NAMBIENT    = (real) refractive index of ambient medium
!    WAVE        = wavelength in vacuo of incident wave (phys. units)
!    DPHYS       = interdipole separation (phys. units)
!    NAT0        = number of dipoles in physical target
!    NX,NY,NZ    = dimensions/d of computational volume
!    CXE_INC(1-3)= (complex) incident E field at location (x_TF,y_TF,z_TF)
!    CXE_SCA(1-3)= (complex) scattered E field at location (x_TF,y_TF,z_TF)
!    CXP(1-3)    = (complex) polarization/d^3 at location (x_TF,y_TF,z_TF)
!
! B.T. Draine, Princeton Univ. Observatory, 2011.08.30
! history
! 11.08.29 (BTD) completed working version
! 11.08.31 (BTD) * added NAMBIENT to argument list
!                * added NAMBIENT to READ statement
!                * changed FORM='UNFORMATTED' to ACCESS='STREAM'
! end history
! Copyright (C) 2011
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!=======================================================================
      DATA CFLENAME_OLD/'null'/
      SAVE NXYZ
      SAVE CFLENAME_OLD,NXY,PI,X0,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
      SAVE CXEINC,CXESCA,CXPOL
      IDVOUT=0
!*** diagnostic
!      write(idvout,*)'readE ckpt 0, x_TF(1-3)=',xtf
!***
      IF(CFLENAME/=CFLENAME_OLD)THEN
         PI=4._WP*ATAN(1._WP)
         WRITE(IDVOUT,FMT='(A,A)')'>READE read file = ',CFLENAME
         CFLENAME_OLD=CFLENAME
         OPEN(UNIT=17,FILE=CFLENAME,ACCESS='STREAM')
         READ(17)NXYZ,NAT0,NAT3,NX,NY,NZ,X0,AEFF,NAMBIENT,WAVE
         DPHYS=AEFF*(4._WP*PI/(3._WP*NAT0))**(1._WP/3._WP)
         NXY=NX*NY
         XMIN=(X0(1)+1.-0.501)*DPHYS
         XMAX=(X0(1)+NX+0.501)*DPHYS
         YMIN=(X0(2)+1.-0.501)*DPHYS
         YMAX=(X0(2)+NY+0.501)*DPHYS
         ZMIN=(X0(3)+1.-0.501)*DPHYS
         ZMAX=(X0(3)+NZ+0.501)*DPHYS

!*** diagnostic
!         write(idvout,*)'readE ckpt 1'
!         write(idvout,*)' nat0=',nat0
!         write(idvout,*)' nxyz=',nxyz
!         write(idvout,*)' nat3=',nat3
!         write(idvout,*)'   nx=',nx
!         write(idvout,*)'   ny=',ny
!         write(idvout,*)'   nz=',nz
!         write(idvout,*)'   x0=',x0
!         write(idvout,*)' aeff=',aeff
!         write(idvout,*)' wave=',wave
!         write(idvout,*)'nambient=',nambient
!***
         ALLOCATE(CXPOL(1:3*NXYZ))
         ALLOCATE(CXEINC(1:3*NXYZ))
         ALLOCATE(CXESCA(1:3*NXYZ))

!*** diagnostic
!         write(idvout,*)'readE ckpt 2'
!***

         READ(17)CXPOL,CXEINC,CXESCA

!*** diagnostic
!         write(idvout,*)'readE ckpt 3'
!***

         ALLOCATE(CXADIA(1:3*NXYZ))
         ALLOCATE(ICOMP(1:3*NXYZ))

!*** diagnostic
!         write(idvout,*)'readE ckpt 4'
!***

         READ(17)CXADIA,ICOMP

!*** diagnostic
!         write(idvout,*)'readE ckpt 5'
!***

         CLOSE(17)

!*** diagnostic
!         write(idvout,*)'readE ckpt 6: check CXPOL for NaN...'
!         j2=0
!         do j1=1,nat3
!            if(.not.(abs(cxpol(j1))>=0.d0))then
!               write(idvout,*)'j1=',j1,' cxpol(j1)=',cxpol(j1)
!               j2=j2+1
!            endif
!         enddo
!         write(idvout,*)'... CXPOL checked for NaN'
!         write(idvout,*)'    j2=',j2,' instances'
!***
!*** diagnostic
!         write(idvout,*)'readE ckpt 7: check CXEINC for NaN...'
!         j2=0
!         do j1=1,nat3
!            if(.not.(abs(cxeinc(j1))>=0.d0))then
!               write(idvout,*)'j1=',j1,' cxeinc(j1)=',cxeinc(j1)
!               j2=j2+1
!            endif
!         enddo
!         write(idvout,*)'... CXEINC checked for NaN'
!         write(idvout,*)'    j2=',j2,' instances'
!***
!*** diagnostic
!         write(idvout,*)'readE ckpt 8: check CXESCA for NaN...'
!         j2=0
!         do j1=1,nat3
!            if(.not.(abs(cxesca(j1))>=0.d0))then
!               write(idvout,*)'j1=',j1,' cxesca(j1)=',cxesca(j1)
!               j2=j2+1
!            endif
!         enddo
!         write(idvout,*)'... CXESCA checked for NaN'
!         write(idvout,*)'    j2=',j2,' instances'
!***
! check solution
         EINC2=0.
         DO J=0,2
            EINC2=EINC2+CXEINC(1+J*NXYZ)*CONJG(CXEINC(1+J*NXYZ))
         ENDDO
!*** diagnostic
!         write(0,*)'readE ckpt 9 : EINC2=',EINC2
!***
         SUMERR2=0.
         J1=0
         DO J=1,NAT3
            IF(ICOMP(J)>0)THEN
               CXERR=CXPOL(J)*CXADIA(J)-(CXEINC(J)+CXESCA(J))
               SUMERR2=SUMERR2+CXERR*CONJG(CXERR)
               J1=J1+1
!*** diagnostic
!               IF(J1<=3)THEN
!                  WRITE(IDVOUT,*)'readE ckpt 10'
!                  WRITE(IDVOUT,*)'------------------------ j1 =',j1, &
!                            ' ---------------'
!                  WRITE(IDVOUT,FMT='(1P6E11.3,A)')CXEINC(J),CXEINC(J+NXYZ), &
!                                             CXEINC(J+2*NXYZ),' = E_inc'
!                  WRITE(IDVOUT,FMT='(1P6E11.3,A)')CXESCA(J),CXESCA(J+NXYZ), &
!                                             CXESCA(J+2*NXYZ),' = E_sca'
!                  WRITE(IDVOUT,FMT='(1P6E11.3,A)')(CXEINC(J)+CXESCA(J)),      &
!                                         (CXEINC(J+NXYZ)+CXESCA(J+NXYZ)),     &
!                                         (CXEINC(J+2*NXYZ)+CXESCA(J+2*NXYZ)), &
!                                          ' = E'
!                  WRITE(IDVOUT,FMT='(1P6E11.3,A)')(CXADIA(J)*CXPOL(J)),       &
!                                          (CXADIA(J+NXYZ)*CXPOL(J+NXYZ)),     &
!                                          (CXADIA(J+2*NXYZ)*CXPOL(J+2*NXYZ)), &
!                                          ' = A*P'
!                  WRITE(IDVOUT,FMT='(1P6E11.3,A)')CXPOL(J),CXPOL(J+NXYZ), &
!                                             CXPOL(J+2*NXYZ),' = P'
!                  WRITE(IDVOUT,FMT='(1P6E11.3,A)')CXADIA(J),CXADIA(J+NXYZ), &
!                                             CXADIA(J+2*NXYZ),' = A'
!               ENDIF
!***
            ENDIF
         ENDDO
         SUMERR2=SUMERR2/(REAL(NAT0)*EINC2)
         WRITE(IDVOUT,FMT='(A,1PE11.4,A)')'>READE', &
               SUMERR2,' = normalized error |P/alpha-E|^2/|E_inc|^2'
         WRITE(IDVOUT,FMT='(A,1PE11.4,A)')'>READE',AEFF,     &
               ' = AEFF (vol. equiv. radius, phys. units)'
         WRITE(IDVOUT,FMT='(A,I11,A)')'>READE',NAT0,          &
               ' = NAT0 (number of physical dipoles in target)'
         WRITE(IDVOUT,FMT='(A,1PE11.4,A)')'>READE',              &
               DPHYS,' = d = interdipole separation (phys. units)'
         WRITE(IDVOUT,FMT='(A,1PE11.4,A)')'>READE',      &
               WAVE,' = wavelength in vacuo (phys. units)'
         WRITE(IDVOUT,FMT='(A,1PE11.4,A)')'>READE',                        &
               WAVE/NAMBIENT,' = wavelength in ambient medium (phys. units)'
! write out the lattice site locations
!         write(9,fmt='(1pe11.3,a)')dphys,' = dphys'
!         do j=1,nx
!            w1=(x0(1)+j)*dphys
!            write(9,fmt='(i3,1pe11.3)')j,w1
!         enddo
!***

      ENDIF

!==============================================================================
! check that XTF(1-3) is within computational volume

!*** diagnostic
!      write(idvout,*)'readE ckpt 11'
!***
      
      IF(XTF(1).LT.XMIN)THEN
         WRITE(IDVOUT,'(A,F8.4,A,F8.4,A)')'Fatal error in readE: x_TF=', &
                                          XTF(1),' < x_min=',XMIN,' STOP'
         STOP
      ENDIF
      IF(XTF(1).GT.XMAX)THEN
         WRITE(IDVOUT,'(A,F8.4,A,F8.4,A)')'Fatal error in readE: x_TF=', &
                                          XTF(1),' > x_max=',XMAX,' STOP'
         STOP
      ENDIF
      IF(XTF(2).LT.YMIN)THEN
         WRITE(IDVOUT,'(A,F8.4,A,F8.4,A)')'Fatal error in readE: y_TF=', &
                                          XTF(2),' < y_min=',YMIN,' STOP'
         STOP
      ENDIF
      IF(XTF(2).GT.YMAX)THEN
         WRITE(IDVOUT,'(A,F8.4,A,F8.4,A)')'Fatal error in readE: y_TF=', &
                                          XTF(2),' > y_max=',YMAX,' STOP'
         STOP
      ENDIF
      IF(XTF(3).LT.ZMIN)THEN
         WRITE(IDVOUT,'(A,F8.4,A,F8.4,A)')'Fatal error in readE: z_TF=', &
                                          XTF(3),' < z_min=',ZMIN,' STOP'
         STOP
      ENDIF
      IF(XTF(3).GT.ZMAX)THEN
         WRITE(IDVOUT,'(A,F8.4,A,F8.4,A)')'Fatal error in readE: z_TF=', &
                                          XTF(3),' > z_max=',ZMAX,' STOP'
         STOP
      ENDIF

!*** diagnostic
!      write(idvout,*)'readE ckpt 12'
!***

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

!*** diagnostic
!      write(idvout,*)'readE ckpt 13'
!      write(idvout,*)'  xmin,xmax=',xmin,xmax
!      write(idvout,*)'  xtf(1-3)=',xtf
!      write(idvout,*)'  ix1,iy1,iz1=',ix1,iy1,iz1
!***
!                                      zyx
      J1X=IX1+NX*(IY1-1)+NXY*(IZ1-1) ! 000
      J2X=J1X+1                      ! 001
      J3X=J1X+NX                     ! 010
      J4X=J3X+1                      ! 011
      J5X=J1X+NXY                    ! 100
      J6X=J5X+1                      ! 101
      J7X=J5X+NX                     ! 110
      J8X=J7X+1                      ! 111
      J1Y=J1X+NXYZ
      J2Y=J2X+NXYZ
      J3Y=J3X+NXYZ
      J4Y=J4X+NXYZ
      J5Y=J5X+NXYZ
      J6Y=J6X+NXYZ
      J7Y=J7X+NXYZ
      J8Y=J8X+NXYZ
      J1Z=J1Y+NXYZ
      J2Z=J2Y+NXYZ
      J3Z=J3Y+NXYZ
      J4Z=J4Y+NXYZ
      J5Z=J5Y+NXYZ
      J6Z=J6Y+NXYZ
      J7Z=J7Y+NXYZ
      J8Z=J8Y+NXYZ
      WX=XTF(1)/DPHYS-X0(1)-IX1
      WY=XTF(2)/DPHYS-X0(2)-IY1
      WZ=XTF(3)/DPHYS-X0(3)-IZ1
      W1=(1.-WX)*(1.-WY)*(1.-WZ) ! 000
      W2=WX*(1.-WY)*(1.-WZ)      ! 001
      W3=(1.-WX)*WY*(1.-WZ)      ! 010
      W4=WX*WY*(1.-WZ)           ! 011
      W5=(1.-WX)*(1.-WY)*WZ      ! 100
      W6=WX*(1.-WY)*WZ           ! 101
      W7=(1.-WX)*WY*WZ           ! 110
      W8=WX*WY*WZ                ! 111

!*** diagnostic
!      write(idvout,*)'readE ckpt 14'
!      write(idvout,*)'   wx,wy,wz=',wx,wy,wz
!      write(idvout,fmt='(a,1i8)')' 3*nxyz=',(3*nxyz)
!      write(idvout,fmt='(a,8i8)')'j1x-j8x=',j1x,j2x,j3x,j4x,j5x,j6x,j7x,j8x
!      write(idvout,fmt='(a,8i8)')'j1y-j8y=',j1y,j2y,j3y,j4y,j5y,j6y,j7y,j8y
!      write(idvout,fmt='(a,8i8)')'j1z-j8z=',j1z,j2z,j3z,j4z,j5z,j6z,j7z,j8z
!      write(idvout,fmt='(a,8f8.4)')'  w1-w8=',w1,w2,w3,w4,w5,w6,w7,w8
!      write(idvout,*)' sum(w1-w8)=',(w1+w2+w3+w4+w5+w6+w7+w8)
!      write(idvout,*)'   nxyz=',nxyz
!
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j1x=',j1x,' cxeinc(j1x)=',cxeinc(j1x)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j2x=',j2x,' cxeinc(j2x)=',cxeinc(j2x)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j3x=',j3x,' cxeinc(j3x)=',cxeinc(j3x)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j4x=',j4x,' cxeinc(j4x)=',cxeinc(j4x)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j5x=',j5x,' cxeinc(j5x)=',cxeinc(j5x)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j6x=',j6x,' cxeinc(j6x)=',cxeinc(j6x)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j7x=',j7x,' cxeinc(j7x)=',cxeinc(j7x)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j8x=',j8x,' cxeinc(j8x)=',cxeinc(j8x)
!      write(idvout,*)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j1y=',j1y,' cxeinc(j1y)=',cxeinc(j1y)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j2y=',j2y,' cxeinc(j2y)=',cxeinc(j2y)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j3y=',j3y,' cxeinc(j3y)=',cxeinc(j3y)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j4y=',j4y,' cxeinc(j4y)=',cxeinc(j4y)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j5y=',j5y,' cxeinc(j5y)=',cxeinc(j5y)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j6y=',j6y,' cxeinc(j6y)=',cxeinc(j6y)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j7y=',j7y,' cxeinc(j7y)=',cxeinc(j7y)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j8y=',j8y,' cxeinc(j8y)=',cxeinc(j8y)
!      write(idvout,*)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j1z=',j1z,' cxeinc(j1z)=',cxeinc(j1z)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j2z=',j2z,' cxeinc(j2z)=',cxeinc(j2z)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j3z=',j3z,' cxeinc(j3z)=',cxeinc(j3z)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j4z=',j4z,' cxeinc(j4z)=',cxeinc(j4z)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j5z=',j5z,' cxeinc(j5z)=',cxeinc(j5z)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j6z=',j6z,' cxeinc(j6z)=',cxeinc(j6z)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j7z=',j7z,' cxeinc(j7z)=',cxeinc(j7z)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j8z=',j8z,' cxeinc(j8z)=',cxeinc(j8z)
!      write(idvout,*)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j1x=',j1x,' cxesca(j1x)=',cxesca(j1x)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j2x=',j2x,' cxesca(j2x)=',cxesca(j2x)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j3x=',j3x,' cxesca(j3x)=',cxesca(j3x)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j4x=',j4x,' cxesca(j4x)=',cxesca(j4x)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j5x=',j5x,' cxesca(j5x)=',cxesca(j5x)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j6x=',j6x,' cxesca(j6x)=',cxesca(j6x)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j7x=',j7x,' cxesca(j7x)=',cxesca(j7x)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j8x=',j8x,' cxesca(j8x)=',cxesca(j8x)
!      write(idvout,*)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j1y=',j1y,' cxesca(j1y)=',cxesca(j1y)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j2y=',j2y,' cxesca(j2y)=',cxesca(j2y)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j3y=',j3y,' cxesca(j3y)=',cxesca(j3y)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j4y=',j4y,' cxesca(j4y)=',cxesca(j4y)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j5y=',j5y,' cxesca(j5y)=',cxesca(j5y)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j6y=',j6y,' cxesca(j6y)=',cxesca(j6y)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j7y=',j7y,' cxesca(j7y)=',cxesca(j7y)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j8y=',j8y,' cxesca(j8y)=',cxesca(j8y)
!      write(idvout,*)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j1z=',j1z,' cxesca(j1z)=',cxesca(j1z)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j2z=',j2z,' cxesca(j2z)=',cxesca(j2z)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j3z=',j3z,' cxesca(j3z)=',cxesca(j3z)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j4z=',j4z,' cxesca(j4z)=',cxesca(j4z)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j5z=',j5z,' cxesca(j5z)=',cxesca(j5z)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j6z=',j6z,' cxesca(j6z)=',cxesca(j6z)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j7z=',j7z,' cxesca(j7z)=',cxesca(j7z)
!      write(idvout,fmt='(a,i7,a,2f8.4)')' j8z=',j8z,' cxesca(j8z)=',cxesca(j8z)
!***

      CXE_INC(1)=W1*CXEINC(J1X)+W2*CXEINC(J2X)+W3*CXEINC(J3X)+W4*CXEINC(J4X)+ &
                 W5*CXEINC(J5X)+W6*CXEINC(J6X)+W7*CXEINC(J7X)+W8*CXEINC(J8X)
      CXE_INC(2)=W1*CXEINC(J1Y)+W2*CXEINC(J2Y)+W3*CXEINC(J3Y)+W4*CXEINC(J4Y)+ &
                 W5*CXEINC(J5Y)+W6*CXEINC(J6Y)+W7*CXEINC(J7Y)+W8*CXEINC(J8Y)
      CXE_INC(3)=W1*CXEINC(J1Z)+W2*CXEINC(J2Z)+W3*CXEINC(J3Z)+W4*CXEINC(J4Z)+ &
                 W5*CXEINC(J5Z)+W6*CXEINC(J6Z)+W7*CXEINC(J7Z)+W8*CXEINC(J8Z)

!*** diagnostic
!      write(idvout,*)'readE ckpt 15'
!***

      CXE_SCA(1)=W1*CXESCA(J1X)+W2*CXESCA(J2X)+W3*CXESCA(J3X)+W4*CXESCA(J4X)+ &
                 W5*CXESCA(J5X)+W6*CXESCA(J6X)+W7*CXESCA(J7X)+W8*CXESCA(J8X)
      CXE_SCA(2)=W1*CXESCA(J1Y)+W2*CXESCA(J2Y)+W3*CXESCA(J3Y)+W4*CXESCA(J4Y)+ &
                 W5*CXESCA(J5Y)+W6*CXESCA(J6Y)+W7*CXESCA(J7Y)+W8*CXESCA(J8Y)
      CXE_SCA(3)=W1*CXESCA(J1Z)+W2*CXESCA(J2Z)+W3*CXESCA(J3Z)+W4*CXESCA(J4Z)+ &
                 W5*CXESCA(J5Z)+W6*CXESCA(J6Z)+W7*CXESCA(J7Z)+W8*CXESCA(J8Z)

!*** diagnostic
!      write(idvout,*)'readE ckpt 16'
!***

      CXP(1)=W1*CXPOL(J1X)+W2*CXPOL(J2X)+W3*CXPOL(J3X)+W4*CXPOL(J4X)+ &
             W5*CXPOL(J5X)+W6*CXPOL(J6X)+W7*CXPOL(J7X)+W8*CXPOL(J8X)
      CXP(2)=W1*CXPOL(J1Y)+W2*CXPOL(J2Y)+W3*CXPOL(J3Y)+W4*CXPOL(J4Y)+ &
             W5*CXPOL(J5Y)+W6*CXPOL(J6Y)+W7*CXPOL(J7Y)+W8*CXPOL(J8Y)
      CXP(3)=W1*CXPOL(J1Z)+W2*CXPOL(J2Z)+W3*CXPOL(J3Z)+W4*CXPOL(J4Z)+ &
             W5*CXPOL(J5Z)+W6*CXPOL(J6Z)+W7*CXPOL(J7Z)+W8*CXPOL(J8Z)

!*** diagnostic
!      write(idvout,*)'readE ckpt 17'
!      write(idvout,*)' cxe_inc(1)=',cxe_inc(1)
!      write(idvout,*)' cxe_inc(2)=',cxe_inc(2)
!      write(idvout,*)' cxe_inc(3)=',cxe_inc(3)
!      w1=cxe_inc(1)*conjg(cxe_inc(1))+cxe_inc(2)*conjg(cxe_inc(2))+ &
!         cxe_inc(3)*conjg(cxe_inc(3))
!      write(idvout,*)' |cxe_inc|=',w1
!      write(idvout,*)' cxe_sca(1)=',cxe_sca(1)
!      write(idvout,*)' cxe_sca(2)=',cxe_sca(2)
!      write(idvout,*)' cxe_sca(3)=',cxe_sca(3)
!      write(idvout,*)' cxp(1)=',cxp(1)
!      write(idvout,*)' cxp(2)=',cxp(2)
!      write(idvout,*)' cxp(3)=',cxp(3)
!***

!*** diagnostic
!      e2=(cxe_inc(1)+cxe_sca(1))*conjg(cxe_inc(1)+cxe_sca(1))+ &
!         (cxe_inc(2)+cxe_sca(2))*conjg(cxe_inc(2)+cxe_sca(2))+ &
!         (cxe_inc(3)+cxe_sca(3))*conjg(cxe_inc(3)+cxe_sca(3))
!      write(10,fmt='(f7.4,3i3,8f7.3,f7.4)')xtf(1),ix1,iy1,iz1,   &
!                                        w1,w2,w3,w4,w5,w6,w7,w8,e2
!***      
      RETURN

   END SUBROUTINE READE

