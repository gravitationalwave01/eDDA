    SUBROUTINE PBCSCAVEC(MXSCA,JPBC,NSCAT,PYDDX,PZDDX,A1,A2,THETA,BETA,    &
                         XLR,YLR,ZLR,AK1,EN0R,CXE01R,ORDERM,ORDERN,THETAN, &
                         PHIN,AKSR,EM1R,EM2R,ENSC,EM1,EM2)
! --------------------------------- v4 ------------------------------------
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE
      INTEGER :: JPBC,MXSCA,NSCAT
      REAL(WP) :: AK1,BETA,THETA,PYDDX,PZDDX
      REAL(WP) ::      &
        A1(3),         &
        A2(3),         &
        AKSR(3,MXSCA), &
        EM1(3,MXSCA),  &
        EM2(3,MXSCA),  &
        EM1R(3,MXSCA), &
        EM2R(3,MXSCA), &
        EN0R(3),       &
        ENSC(3,MXSCA), &
        ORDERM(MXSCA), &
        ORDERN(MXSCA), &
        PHIN(MXSCA),   &
        THETAN(MXSCA), &
        XLR(3),        &
        YLR(3),        &
        ZLR(3)
      COMPLEX(WP) :: &
        CXE01R(3)

! local variables:

      INTEGER :: J,JR,JT,K
      REAL(WP) :: AK2,AKPERP,AKPERP0,AKPERP2,                &
                  COSALPHAI,COSPHI,COSTHETA,DENOM,PI,SINALPHAI,SINPHI,SUM
      REAL(WP) :: &
        EM01R(3)

!-----------------------------------------------------------------------
! subroutine PBCSCAVEC
! given:

!     JPBC = 1 for target periodic in y_TF only (PZD=0)
!            2                        z_TF only (PYD=0)
!            3                        both y_TF and z_TF

!     MXSCA = dimensioning information
!     PYDDX = periodicity in y direction/d
!     PZDDX = periodicity in z direction/d
!     A1(3) = target axis a1 in TF (TF = Target Frame)    [**NOT USED**]
!     A2(3) = target axis a2 in TF                        [**NOT USED**]
!     THETA = rotation angle (rad) of a1 relative to x_LF
!     BETA  = rotation angle (rad) of target around axis a1
!     XLR(3)= unit vector xlab in TF = direction of incident radiation
!     YLR(3)= unit vector ylab in TF
!     ZLR(3)= unit vector zlab in TF

!     AK1   = k*d in vacuo
!     NSCAT = number of scattering directions
!     EN0R(1-3) = unit vector in incident wave direction in TF
!     CXE01R(1-3) = (complex) incident pol.vector 1 in TF
!     ORDERM(1-NSCAT)=diffraction order in y_TF direction if JPBC=1 or 3
!                                       in z_TF direction if JPBC=2
!     ORDERN(1-NSCAT)=rotation angle (rad) around scattering cone if
!                     JPBC=1 or 2
!                    =diffraction order in z_TF direction if JPBC=3

! returns:

!     AKSR(3,1-NSCAT) = scattering vectors in TF
!     EM1R(3,1-NSCAT) = scattered polarization vector parallel to
!                       scattering plane in TF
!     EM2R(3,1-NSCAT) = scattered polarization vector perpendicular to
!                       scattering plane in TF
!     ENSC(3,1-NSCAT) = scattering unit vectors in Lab Frame
!     EM1(3,1-NSCAT)  = scattered polarization vector parallel to
!                       scattering plane in LF
!     EM2(3,1-NSCAT)  = scattered polarization vector perpendicular to
!                       scattering plane in LF
! if JPBC=1 or 2:
!     THETAN(1-NSCAT) = alpha_s (rad) in TF
!     PHIN(1-NSCAT)   = zeta_s (rad) in TF

! where the scattering directions are determined as follows:

! if JPBC=1:
!    input ORDERM(J) = scattering order for y (integer)
!          ORDERN(J) = azimuthal angle zeta (radians) around y_TF 
!                      where ORDERM=0 and zeta=0 for forward scattering
!    computes THETAN(J) = angle between k_s and y_TF
!             PHIN(J)   = azimuthal angle zeta (radians)

! if JPBC=2:
!    input ORDERM(J) = scattering order for z (integer)
!          ORDERN(J) = azimuthal angle zeta (radians) around z_TF 
!                      where ORDERM=0 and zeta=0 for forward scattering
!    computes THETAN(J) = angle between k_s and z_TF
!             PHIN(J)   = azimuthal angle zeta (radians)

! if JPBC=3:
!    input ORDERM(J) = scattering order for y (integer)
!          ORDERN(J) = scattering order for z (integer)
!                      where ORDERM=0 and ORDERM=0 for forward scattering
!
!     NSCAT is even: each scattering order (M,N) appears twice:
!                    first for transmission, later for reflection
!                    (this is established in REAPAR)
!       J = 1         -> NSCAT/2 correspond to transmission
!       J = NSCAT/2+1 -> NSCAT   correspond to reflection

!     computes THETAN(J) = angle (rad) between k_s and x_TF
!              PHIN(J)   = azimuthal angle (around x_TF)
!                          k_sx=k_s*cos(thetan)
!                          k_sy=k_s*sin(thetan)*cos(phin)
!                          k_sz=k_s*sin(thetan)*sin(phin)
!
!     in special case of zeroth-order transmission, the angle
!     THETAN(J)=0, and the angle PHIN(J) is taken to be the
!     same as for (0,0) reflection.
!     In the special case of (0,0) normal incidence, where k_sy=k_sz=0,
!     we set PHIN=0 

! B.T. Draine, Princeton University Observatory, 2006.10.08
! history
! 06.10.08 (BTD) First written
! 06.10.13 (BTD) fixed bug when JPBC=2
!                added fatal error warning
! 06.10.13 (BTD) added new arguments ORDERM,ORDERN
!                to preserve scattering order information
!                as orientation is changed
! 06.10.14 (BTD) added code to guard against numerical problem
!                evaluating angle PHIN when incident radiation is
!                normal (or near-normal) to target axis (for JPBC=1 or
!                2) or target plane (for JPBC=3)
!                Now treat near-normal as exactly normal when
!                calculating angle PHIN
! 06.12.28 (BTD) corrected error in computation of unit scattered pol.
!                vector in special cases
!                added arrays XLR,YLR,ZLR to argument list
!                use XLR,YLR,ZLR in defining unit pol vectors in special
!                case of scattering angle theta = 0 or pi (in which
!                case scattering plane is not uniquely defined).
! 07.01.20 (BTD) corrected bug in assignment of transmitted/reflected
!                unit vectors when JPBC=3 (2-d periodic target)
! 07.06.22 (BTD) modified so that azimuthal angle = 0 for forward
!                scattering when JPBC=1 or 2
!                introduce new variable AKP1(1-3) and AKP2(1-3)
!                revise calculation of THETAN(J) and PHIN(J) for the
!                JPBC=1, 2, and 3
! 07.06.28 (BTD) corrected error in computation of AKSR when JPBC=1 or 2
!                and incident radiation is not perpendicular to
!                periodicity direction.
! 07.07.01 (BTD) corrected error in computation of EM1R
! 07.07.08 (BTD) added EM1,EM2 to argument list
!                added code to compute EM1,EM2
! 07.10.03 (BTD) corrected error in computation of angle PHIN(J) when 
!                JPBC=3
! 07.10.04 (BTD) corrected error in evalation of AKPERP when JPBC=3
! 07.10.06 (BTD) modified handling of JPBC=3:
!                * use (0,0) reflection to establish PHIN(J) for (0,0)
!                  transmission
!                * in case of (0,0) and normal incidence,
!                  set PHIN(J)=0
! 07.10.07 (BTD) * added THETA,BETA to argument list
!                * rewrote treatment of JPBC=1 and 2
! 08.01.02 (BTD) v2 Revised code to calculate scattered pol vectors
!                   when JPBC=3:
! 08.01.03 (BTD) v2 Another major revision for JPBC=3 case
! 08.01.10 (BTD) v3 Revised treatment for JPBC=3 to conform to
!                   Bohren & Huffman convention for ehat_{spar}
!                   and ehat_{sperp}
! 08.04.09 (BTD) v4 Revised creation of scattering vectors for JPBC=1 or 2
! 08.06.30 (BTD) v4 Fixed bug in determining PHIN for special case
!                   THETA=0.
! end history
!-----------------------------------------------------------------------
!*** diagnostic
!      write(0,*)'pbcscavec ckpt 0, en0r(1-3)=',en0r
!***

      PI=4._WP*ATAN(1._WP)

! for use in this module, define
! EM01R = unit polarization vector for incident polarization state 1 in TF
! AK2 = |kd|^2

      DO K=1,3
         EM01R(K)=REAL(CXE01R(K))
      ENDDO
      SUM=SQRT(EM01R(1)**2+EM01R(2)**2+EM01R(3)**2)
      DO K=1,3
         EM01R(K)=EM01R(K)/SUM
      ENDDO
      AK2=AK1*AK1

      IF(JPBC==1)THEN

! JPBC=1 : target periodic in y_TF direction
!       -> input ORDERM = diffraction order in y_TF direction
!                ORDERN = azimuthal angle around y_TF axis
!                         with ORDERN=0 for forward scattering
!*** diagnostic
!         write(0,*)'pbcscavec ckpt 1'
!***

         AKPERP0=AK1*SQRT(EN0R(1)**2+EN0R(3)**2)
         DO J=1,NSCAT
            AKSR(2,J)=AK1*EN0R(2)+ORDERM(J)*2._WP*PI/PYDDX
            IF(ABS(AKSR(2,J))<=AK1)THEN
               PHIN(J)=ORDERN(J)
               THETAN(J)=ACOS(AKSR(2,J)/AK1)
               AKPERP=SQRT(AK2-AKSR(2,J)**2)
               AKSR(3,J)=(AKPERP/AKPERP0)*AK1*                     &
                         (EN0R(3)*COS(PHIN(J))-EN0R(1)*SIN(PHIN(J)))
               AKSR(1,J)=(AKPERP/AKPERP0)*AK1*                     &
                         (EN0R(1)*COS(PHIN(J))+EN0R(3)*SIN(PHIN(J)))
            ELSE
               CALL ERRMSG('FATAL','PBCSCAVEC','invalid diffraction order')
            ENDIF
         ENDDO

      ELSEIF(JPBC==2)THEN

! JPBC=2 : target periodic in z direction
!       -> input ORDERM = diffraction order in z_TF direction
!                ORDERN = azimuthal angle around z_TF axis
!                         with ORDERN=0 for forward scattering

!*** diagnostic
!         write(0,*)'pbcscavec ckpt 2'
!***

         AKPERP0=AK1*SQRT(EN0R(1)**2+EN0R(2)**2)
         DO J=1,NSCAT
            AKSR(3,J)=AK1*EN0R(3)+ORDERM(J)*2._WP*PI/PZDDX
            IF(ABS(AKSR(3,J))<=AK1)THEN
               PHIN(J)=ORDERN(J)
               THETAN(J)=ACOS(AKSR(3,J)/AK1)
               AKPERP=SQRT(AK2-AKSR(3,J)**2)
               AKSR(1,J)=(AKPERP/AKPERP0)*AK1*                     &
                         (EN0R(1)*COS(PHIN(J))-EN0R(2)*SIN(PHIN(J)))
               AKSR(2,J)=(AKPERP/AKPERP0)*AK1*                     &
                         (EN0R(2)*COS(PHIN(J))+EN0R(1)*SIN(PHIN(J)))
            ELSE
               CALL ERRMSG('FATAL','PBCSCAVEC','invalid diffraction order')
            ENDIF
         ENDDO

!*** diagnostic
!         write(0,*)'pbcscavec ckpt 2.5'
!***

      ELSEIF(JPBC==3)THEN

! JPBC=3 : infinite and periodic in y and z directions
!       -> input ORDERM = diffraction order in y_TF direction
!                ORDERN = diffraction order in z_TF direction
!          NSCAT is assumed to be an even number (this is set by
!                subroutine REAPAR)
!          every (M,N) case appears twice: 
!             first for transmission (directions J = 1 -> NSCAT/2)
!             then for reflection (directions J = NSCAT/2+1 -> NSCAT)
!          EN0R(1-3) = components of unit vector for incident direction
!                      in TF
!          calculate
!             AKSR(1-3,J) = scattered vector in TF for J=1,NSCAT
!             THETAN(J)   = angle between k_s and x_TF
!             PHIN(J)     = azimuthal angle in TF
!             where  k_sx = k_s*cos(thetan)
!                    k_sy = k_s*sin(thetan)*cos(phin)
!                    k_sz = k_s*sin(thetan)*sin(phin)
!             Note: for special case k_sy=k_sz=0, we set phin=0

!*** diagnostic
!         write(0,*)'pbcscavec ckpt 3'
!***

         DO JT=1,NSCAT/2
            JR=NSCAT/2+JT
            AKSR(2,JT)=AK1*EN0R(2)+ORDERM(JT)*2._WP*PI/PYDDX
            AKSR(3,JT)=AK1*EN0R(3)+ORDERN(JT)*2._WP*PI/PZDDX
            AKSR(2,JR)=AKSR(2,JT)
            AKSR(3,JR)=AKSR(3,JT)
            AKPERP2=AKSR(2,JR)**2+AKSR(3,JR)**2
!*** diagnostic
!            write(0,*)'pbcscavec ckpt 3.1, jt=',jt
!            write(0,*)'                    ak1=',ak1
!            write(0,*)'                en0r(2)=',en0r(2)
!            write(0,*)'                en0r(3)=',en0r(3)
!            write(0,*)'             orderm(jt)=',orderm(jt)
!            write(0,*)'             ordern(jt)=',ordern(jt)
!            write(0,*)'                     pi=',pi
!            write(0,*)'                  pyddx=',pyddx
!            write(0,*)'                  pzddx=',pzddx
!            write(0,*)'                    aksr(2,jt)=',aksr(2,jt)
!            write(0,*)'                    aksr(3,jt)=',aksr(3,jt)
!***
            IF(AKPERP2<=AK2)THEN
               AKSR(1,JT)=SIGN(SQRT(AK2-AKPERP2),EN0R(1))
               AKSR(1,JR)=-AKSR(1,JT)
!*** diagnostic
!               write(0,*)'pbcscavec ckpt 3.2'
!
               THETAN(JT)=ACOS(AKSR(1,JT)/AK1)
!*** diagnostic
!               write(0,*)'pbcscavec ckpt 3.3'
!               write(0,*)'          aksr(1,jr)=',aksr(1,jr)
!               write(0,*)'                 ak1=',ak1
!***
               THETAN(JR)=ACOS(AKSR(1,JR)/AK1)
!*** diagnostic
!               write(0,*)'pbcscavec ckpt 3.31, thetan(jr)=',thetan(jr)
!***
               AKPERP=SQRT(AKPERP2)
               IF(AKPERP.GT.0._WP)THEN
                  COSPHI=AKSR(2,JT)/AKPERP
                  SINPHI=AKSR(3,JT)/AKPERP
                  IF(SINPHI<0._WP)THEN
                     PHIN(JT)=2._WP*PI-ACOS(COSPHI)
                  ELSE
                     PHIN(JT)=ACOS(COSPHI)
                  ENDIF
               ELSE
                  PHIN(JT)=0._WP
               ENDIF
!*** diagnostic
!               write(0,*)'pbcscavec ckpt 3.34'
!***
               PHIN(JR)=PHIN(JT)
            ELSE
!***
!               write(0,*)'pbcscavec ckpt 3.37, orderm(j)=',orderm(j)
!               write(0,*)'                     ordern(j)=',ordern(j)
!***
               CALL ERRMSG('FATAL','PBCSCAVEC','invalid diffraction order')
            ENDIF
!*** diagnostic
!            write(0,*)'pbcscavec ckpt 3.4'
!***
         ENDDO ! end loop over JT

!*** diagnostic
!         write(0,*)'pbcscavec ckpt 3.5'
!***
      ELSE
         CALL ERRMSG('FATAL','PBCSCAVEC',' invalid value of JPBC')
      ENDIF

!=======================================================================

! scattering directions AKSR are now defined.
! For each transmission direction that is not either parallel or 
! antiparallel to the incident beam there is a scattering plane defined
! by vectors k_0 and k_s.
! Construct unit vector for scattered polarization
! parallel to plane (EM1R) and perpendicular to plane (EM2R)
! We adopt convention of Bohren & Huffman 1984 sec 3.2 for
! unit vectors parallel and perpendicular to the scattering plane:
!   em1 = ehat_{para,s} is in direction of increasing scattering angle theta
!   em2 = ehat_{perp,s} = ehat_{para,s} cross khat_s

! For special case of |k_s cross k_0| = 0  
!   (i.e., k_s dot k_0 = +/- |k_0|^2 , or k_s = +/- k_0: 
!   forward scattering or 180deg backscattering):
!   Let k_sperp = component of k_s perpendicular to target normal
!   If |k_sperp| = 0, then incident and scattered rays are normal
!                     to target.  For this special case, let
!                     
!                     ehat_spar =  ehat_01
!                     ehat_sperp = ehat_spar cross khat_s
!                               
!   If |k_s cross k_0} > 0, then
!                                k_s cross k_0 
!                     e_sperp = ---------------
!                               |k_s cross k_0|
!
! For all cases, we take:
!                     e_iperp = e_sperp
!                     e_iparr = (k_0 cross e_iperp) / |k_0|
!                     e_sparr = (k_s cross e_sperp) / |k_s|

! Here we evaluate em1r(1-3,J) = (x,y,z) components of e_sparr in TF
!                  em2r(1-3,J) = (x,y,z) components of e_sperp in TF

!*** diagnostic
!      write(0,*)'pbcscavec ckpt 4'
!***
      DO J=1,NSCAT

! determine which case we have

         COSTHETA=0._WP
         DO K=1,3
            COSTHETA=COSTHETA+AKSR(K,J)*EN0R(K)
         ENDDO
         COSTHETA=COSTHETA/AK1

         IF(ABS(COSTHETA).LE.0.9999_WP)THEN

! scattering angle is neither zero nor pi
! scattering plane is therefore defined by k_0 and k_s
! obtain EM2R(1-3,J) = components of em2 = (ks cross k0)/|ks cross k0|= in TF

!*** diagnostic
!            write(0,*)'pbcscavec ckpt 5'
!***

            EM2R(1,J)=(AKSR(2,J)*EN0R(3)-AKSR(3,J)*EN0R(2))
            EM2R(2,J)=(AKSR(3,J)*EN0R(1)-AKSR(1,J)*EN0R(3))
            EM2R(3,J)=(AKSR(1,J)*EN0R(2)-AKSR(2,J)*EN0R(1))
            DENOM=SQRT(EM2R(1,J)**2+EM2R(2,J)**2+EM2R(3,J)**2)
            DO K=1,3
               EM2R(K,J)=EM2R(K,J)/DENOM
            ENDDO

!*** diagnostic
!            if(j<=20)then
!            write(0,9970)j,sum,aksr(1,j),aksr(2,j),aksr(3,j),en0r,     &
!                         ak1,em1r(1,j),em1r(2,j),em1r(3,j),            &
!                         sqrt(em1r(1,j)**2+em1r(2,j)**2+em1r(3,j)**2), &
!                         em2r(1,j),em2r(2,j),em2r(3,j),                &
!                         sqrt(em2r(1,j)**2+em2r(2,j)**2+em2r(3,j)**2)
!            endif
! 9970 format('in pbcscavec ckpt alpha: J=',I2,'  sum=',F10.6,/, &
!             '   aksr=',3F10.6,/,       &
!             '   en0r=',3F10.6,/,       &
!             '    ak1=',f10.6,/,        &
!             '   em1r=',3f10.6,/,       &
!             ' |em1r|=',f10.6,/,        &
!             '   em2r=',3f10.6,/,       &
!             ' |em2r|=',f10.6)
!***

         ELSE

! Scattering angle is either zero or pi

!*** diagnostic
!            write(0,*)'pbcscavec ckpt 6'
!***

            IF(JPBC==1.OR.JPBC==2)THEN

! JPBC=1 or 2 scattering angle = 0  -> ORDERM=0 and PHIN=0
!             scattering angle = pi -> ORDERM=0 and PHIN=pi and alpha=pi/2
! Check to see whether scattering angle theta = 0 or pi

!*** diagnostic
!               write(0,*)'pbcscavec ckpt 7'
!***

               IF(COSTHETA>0)THEN

! forward scattering: em2 = [chat - nhat (nhat dot chat)]/sin(alpha)

!*** diagnostic
!                  write(0,*)'pbcscavec ckpt 8'
!***
                  IF(JPBC==1)THEN
                     SINALPHAI=SQRT(1._WP-EN0R(2)**2)
                     COSALPHAI=EN0R(2)
                     EM2R(1,J)=EN0R(1)*COSALPHAI/SINALPHAI
                     EM2R(2,J)=(EN0R(2)*COSALPHAI-1._WP)/SINALPHAI
                     EM2R(3,J)=EN0R(3)*COSALPHAI/SINALPHAI
!*** diagnostic
!                     if(j<=20)then
!                        write(0,*)'pbcscavec ckpt 2, j=',j
!                        write(0,*)'sinalphai=',SINALPHAI
!                        write(0,*)'em2r(1-3,j)=',em2r(1,j),em2r(2,j),em2r(3,j)
!                     endif
!*** end diagnostic
                  ELSE
                     SINALPHAI=SQRT(1._WP-EN0R(3)**2)
                     COSALPHAI=EN0R(3)
                     EM2R(1,J)=EN0R(1)*COSALPHAI/SINALPHAI
                     EM2R(2,J)=EN0R(2)*COSALPHAI/SINALPHAI
                     EM2R(3,J)=(EN0R(3)-1._WP)/SINALPHAI
                  ENDIF

               ELSE

! backscattering: em2 = chat

!*** diagnostic
!                  write(0,*)'pbcscavec ckpt 9'
!***
                  IF(JPBC==1)THEN
                     EM2R(1,J)=0._WP
                     EM2R(2,J)=-1._WP
                     EM2R(3,J)=0._WP
                  ELSE
                     EM2R(1,J)=0._WP
                     EM2R(2,J)=0._WP
                     EM2R(3,J)=-1._WP
                  ENDIF

               ENDIF

            ELSEIF(JPBC==3)THEN

!*** diagnostic
!               write(0,*)'pbcscavec ckpt 10'
!***
! Check to see if k_0 is normal to target plane.

               IF((EN0R(2)**2+EN0R(3)**2)<0.0001_WP)THEN

! incident and scattered ray are both normal to the target plane:

                  DO K=1,3
                     EM1R(K,J)=REAL(CXE01R(K))
                  ENDDO
                  DENOM=SQRT(EM1R(1,J)**2+EM1R(2,J)**2+EM1R(3,J)**2)
                  DO K=1,3
                     EM1R(K,J)=EM1R(K,J)/DENOM
                  ENDDO
                  EM2R(1,J)=(EM1R(2,J)*AKSR(3,J)-EM1R(3,J)*AKSR(2,J))/AK1
                  EM2R(2,J)=(EM1R(3,J)*AKSR(1,J)-EM1R(1,J)*AKSR(3,J))/AK1
                  EM2R(3,J)=(EM1R(1,J)*AKSR(2,J)-EM1R(2,J)*AKSR(1,J))/AK1

               ELSE

! incident ray is not normal to the target plane.  Therefore we can
! use the reflected ray to construct a scattering plane, and use this
! to define the perp polarization direction.
! Note that AKSR for the reflected ray differs from AKSR for the
! transmitted ray only by change in sign of AKSR(1,J)

                  EM2R(1,J)=AKSR(2,J)*EN0R(3)-AKSR(3,J)*EN0R(2)
                  EM2R(2,J)=AKSR(3,J)*EN0R(1)+AKSR(1,J)*EN0R(3)
                  EM2R(3,J)=-AKSR(1,J)*EN0R(2)-AKSR(2,J)*EN0R(1)
                  DENOM=SQRT(EM2R(1,J)**2+EM2R(2,J)**2+EM2R(3,J)**2)
                  DO K=1,3
                     EM2R(K,J)=EM2R(K,J)/DENOM
                  ENDDO
               ENDIF

            ENDIF


         ENDIF ! end IF(ABS(SUM).LE.0.9999_WP) ...
      
! Now use EM2R to obtain EM1R from e_par = khat cross e_perp

         EM1R(1,J)=(AKSR(2,J)*EM2R(3,J)-AKSR(3,J)*EM2R(2,J))/AK1
         EM1R(2,J)=(AKSR(3,J)*EM2R(1,J)-AKSR(1,J)*EM2R(3,J))/AK1

!*** diagnostic
!         write(0,9811)j,aksr(3,j),em2r(1,j),aksr(1,j),em2r(3,j), &
!                         ak1,em1r(2,j)
! 9811 format('in pbcscavec ckpt 11: j=',i2,/,                          &
!             'aksr(3,j)=',f10.6,' em2r(1,j)=',f10.6,' aksr(1,j)=',f10.6, &
!             ' em2r(3,j)=',f10.6,/,'ak1=',f10.6,' em1r(2,j)=',f10.6)
!*** end diagnostic

         EM1R(3,J)=(AKSR(1,J)*EM2R(2,J)-AKSR(2,J)*EM2R(1,J))/AK1

!*** diagnostic
!         write(0,9980)j,sum,aksr(1,j),aksr(2,j),aksr(3,j),en0r,     &
!                         ak1,em1r(1,j),em1r(2,j),em1r(3,j),            &
!                         sqrt(em1r(1,j)**2+em1r(2,j)**2+em1r(3,j)**2), &
!                         em2r(1,j),em2r(2,j),em2r(3,j),                &
!                         sqrt(em2r(1,j)**2+em2r(2,j)**2+em2r(3,j)**2)
! 9980 format('in pbcscavec ckpt 12: J=',I2,'  sum=',F10.6,/, &
!             '   aksr=',3F10.6,/,       &
!             '   en0r=',3F10.6,/,       &
!             '    ak1=',f10.6,/,        &
!             '   em1r=',3f10.6,/,       &
!             ' |em1r|=',f10.6,/,        &
!             '   em2r=',3f10.6,/,       &
!             ' |em2r|=',f10.6)
!***

! now evaluate scattering unit vector in Lab Frame

         ENSC(1,J)=AKSR(1,J)*XLR(1)+AKSR(2,J)*XLR(2)+AKSR(3,J)*XLR(3)
         ENSC(2,J)=AKSR(1,J)*YLR(1)+AKSR(2,J)*YLR(2)+AKSR(3,J)*YLR(3)
         ENSC(3,J)=AKSR(1,J)*ZLR(1)+AKSR(2,J)*ZLR(2)+AKSR(3,J)*ZLR(3)
         DO K=1,3
            ENSC(K,J)=ENSC(K,J)/AK1
         ENDDO

! evaluate scattering polarization vectors in Lab Frame

         EM1(1,J)=EM1R(1,J)*XLR(1)+EM1R(2,J)*XLR(2)+EM1R(3,J)*XLR(3)
         EM1(2,J)=EM1R(1,J)*YLR(1)+EM1R(2,J)*YLR(2)+EM1R(3,J)*YLR(3)
         EM1(3,J)=EM1R(1,J)*ZLR(1)+EM1R(2,J)*ZLR(2)+EM1R(3,J)*ZLR(3)

         EM2(1,J)=EM2R(1,J)*XLR(1)+EM2R(2,J)*XLR(2)+EM2R(3,J)*XLR(3)
         EM2(2,J)=EM2R(1,J)*YLR(1)+EM2R(2,J)*YLR(2)+EM2R(3,J)*YLR(3)
         EM2(3,J)=EM2R(1,J)*ZLR(1)+EM2R(2,J)*ZLR(2)+EM2R(3,J)*ZLR(3)

      ENDDO

!*** diagnostic
!         write(0,*)'pbcscavec ckpt 99'
!         write(0,*)'aksr(1-3,1)=',aksr(1,1),aksr(2,1),aksr(3,1)
!         write(0,*)'aksr(1-3,2)=',aksr(1,2),aksr(2,2),aksr(3,2)
!         write(0,*)'aksr(1-3,3)=',aksr(1,3),aksr(2,3),aksr(3,3)
!         write(0,*)'aksr(1-3,4)=',aksr(1,4),aksr(2,4),aksr(3,4)
!         write(0,*)'aksr(1-3,9)=',aksr(1,9),aksr(2,9),aksr(3,9)
!         write(0,*)'aksr(1-3,10)=',aksr(1,10),aksr(2,10),aksr(3,10)
!         write(0,*)'aksr(1-3,11)=',aksr(1,11),aksr(2,11),aksr(3,11)
!
!         write(0,*)'em1r(1-3,1)=',em1r(1,1),em1r(2,1),em1r(3,1)
!         write(0,*)'em1r(1-3,2)=',em1r(1,2),em1r(2,2),em1r(3,2)
!         write(0,*)'em1r(1-3,3)=',em1r(1,3),em1r(2,3),em1r(3,3)
!         write(0,*)'em1r(1-3,4)=',em1r(1,4),em1r(2,4),em1r(3,4)
!         write(0,*)'em1r(1-3,9)=',em1r(1,9),em1r(2,9),em1r(3,9)
!         write(0,*)'em1r(1-3,10)=',em1r(1,10),em1r(2,10),em1r(3,11)
!         write(0,*)'em1r(1-3,11)=',em1r(1,11),em1r(2,10),em1r(3,11)
!
!         write(0,*)'em2r(1-3,1)=',em2r(1,1),em2r(2,1),em2r(3,1)
!         write(0,*)'em2r(1-3,2)=',em2r(1,2),em2r(2,2),em2r(3,2)
!         write(0,*)'em2r(1-3,3)=',em2r(1,3),em2r(2,3),em2r(3,3)
!         write(0,*)'em2r(1-3,4)=',em2r(1,4),em2r(2,4),em2r(3,4)
!         write(0,*)'em2r(1-3,9)=',em2r(1,9),em2r(2,9),em2r(3,9)
!         write(0,*)'em2r(1-3,10)=',em2r(1,10),em2r(2,10),em2r(3,10)
!         write(0,*)'em2r(1-3,11)=',em2r(1,11),em2r(2,11),em2r(3,11)
!
!         write(0,*)'ensc(1-3,1)=',ensc(1,1),ensc(2,1),ensc(3,1)
!         write(0,*)'ensc(1-3,2)=',ensc(1,2),ensc(2,2),ensc(3,2)
!         write(0,*)'ensc(1-3,3)=',ensc(1,3),ensc(2,3),ensc(3,3)
!         write(0,*)'ensc(1-3,4)=',ensc(1,4),ensc(2,4),ensc(3,4)
!         write(0,*)'ensc(1-3,9)=',ensc(1,9),ensc(2,9),ensc(3,9)
!         write(0,*)'ensc(1-3,10)=',ensc(1,10),ensc(2,10),ensc(3,10)
!         write(0,*)'ensc(1-3,11)=',ensc(1,11),ensc(2,11),ensc(3,11)
!
!         write(0,*)'em1(1-3,1)=',em1(1,1),em1(2,1),em1(3,1)
!         write(0,*)'em1(1-3,2)=',em1(1,2),em1(2,2),em1(3,2)
!         write(0,*)'em1(1-3,3)=',em1(1,3),em1(2,3),em1(3,3)
!         write(0,*)'em1(1-3,4)=',em1(1,4),em1(2,4),em1(3,4)
!         write(0,*)'em1(1-3,9)=',em1(1,9),em1(2,9),em1(3,9)
!         write(0,*)'em1(1-3,10)=',em1(1,10),em1(2,10),em1(3,10)
!         write(0,*)'em1(1-3,11)=',em1(1,11),em1(2,11),em1(3,11)
!
!         write(0,*)'em2(1-3,1)=',em2(1,1),em2(2,1),em2(3,1)
!         write(0,*)'em2(1-3,2)=',em2(1,2),em2(2,2),em2(3,2)
!         write(0,*)'em2(1-3,3)=',em2(1,3),em2(2,3),em2(3,3)
!         write(0,*)'em2(1-3,4)=',em2(1,4),em2(2,4),em2(3,4)
!         write(0,*)'em2(1-3,9)=',em2(1,9),em2(2,9),em2(3,9)
!         write(0,*)'em2(1-3,10)=',em2(1,10),em2(2,10),em2(3,10)
!         write(0,*)'em2(1-3,11)=',em2(1,11),em2(2,11),em2(3,11)
!*** end diagnostic

      RETURN
    END SUBROUTINE PBCSCAVEC
