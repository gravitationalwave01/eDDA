    SUBROUTINE GETMUELLER(IBETA,IORTH,IPHI,ITHETA,JPBC,MXBETA,MXSCA,MXPHI,   &
        MXTHET,NSCAT,ORDERM,ORDERN,CMDTRQ,AK1,AKSR,ENSC,ENSCR,PYDDX,PZDDX,   &
        PHIN,SM,SMORI,S1111,S2121,CX1121,CXE01,CXE02,CXF11,CXF12,CXF21,      &
        CXF22,CXS1,CXS2,CXS3,CXS4,QABS,QABSUM,QBKSCA,QBKSUM,QEXSUM,QEXT,     &
        QPHA,QPHSUM,QSCAG,QSCAG2,QSCAT,QSCG2SUM,QSCGSUM,QSCSUM,QTRQAB,       &
        QTRQABSUM,QTRQSC,QTRQSCSUM,WGTA,WGTB,EM1,EM2,EM1R,EM2R)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

! Arguments:

      INTEGER :: IBETA,IORTH,IPHI,ITHETA,JPBC,MXBETA,MXSCA,MXPHI,MXTHET,NSCAT
      CHARACTER :: CMDTRQ*6
      REAL(WP) :: AK1,PYDDX,PZDDX,WG
      REAL(WP) ::            &
         AKSR(3,MXSCA),      &
         EM1(3,MXSCA),       &
         EM1R(3,MXSCA),      &
         EM2(3,MXSCA),       &
         EM2R(3,MXSCA),      &
         ENSC(3,MXSCA),      &
         ENSCR(3,MXSCA),     &
         ORDERM(MXSCA),      &
         ORDERN(MXSCA),      &
         PHIN(MXSCA),        &
         QABS(2),            &
         QABSUM(2),          &
         QBKSCA(2),          &
         QBKSUM(2),          &
         QEXSUM(2),          &
         QEXT(2),            &
         QPHA(2),            &
         QPHSUM(2),          &
         QSCAG(3,2),         &
         QSCAG2(2),          &
         QSCAT(2),           &
         QSCG2SUM(2),        &
         QSCGSUM(3,2),       &
         QSCSUM(2),          &
         QTRQAB(3,2),        &
         QTRQABSUM(3,2),     &
         QTRQSC(3,2),        &
         QTRQSCSUM(3,2),     &
         S1111(MXSCA),       &
         S2121(MXSCA),       &
         SM(4,4,MXSCA),      &
         SMORI(4,4,MXSCA),   &
         WGTA(MXTHET,MXPHI), &
         WGTB(MXBETA)
      COMPLEX(WP) ::    &
         CX1121(MXSCA), &
         CXE01(3),      &
         CXE02(3),      &
         CXF11(MXSCA),  &
         CXF12(MXSCA),  &
         CXF21(MXSCA),  &
         CXF22(MXSCA),  &
         CXS1(MXSCA),   &
         CXS2(MXSCA),   &
         CXS3(MXSCA),   &
         CXS4(MXSCA)

! Local variables:

      INTEGER :: J,K,JO,ND,ND2
      REAL(WP) :: COSPHI,FAC,FAC0,PI,SINPHI
      COMPLEX(WP) :: CXA,CXB,CXC,CXD,CXI,CXTRM1,CXTRM2

!***********************************************************************
! Given:
!   IBETA =index specifying rotation of target around axis a1
!   IORTH =number of incident polarization states being calculated for
!   IPHI  =index specifying phi value for orientation of target axis a1
!   ITHETA=index specifying theta value for orientation of target axis
!          a1
!   JPBC  = 0 if PBC not used
!           1 if PBC used in y direction only
!           2 if PBC used in z direction only
!           3 if PBC used in both y and z directions
!   MXBETA=dimensioning information
!   MXSCA =dimensioning information
!   MXPHI =dimensioning information
!   MXTHET=dimensioning information
!   NSCAT =number of scattering directions for evaluation of Mueller
!          matrix
!   CMDTRQ='DOTORQ' or 'NOTORQ' depending on whether or not torque
!          calculation is to be carried out
!   ENSC(1-3,1-NSCAT)=normalized scattering direction vector in Lab Frame
!                     for scattering directions 1-NSCAT
!   ENSCR(1-3,1-NSCAT)=normalized scattering direction vector in Target
!                      Frame for scattering directions 1-NSCA
!   EM1(1-3,1-NSCAT)=unit scattering polarization vector 1 in Lab Frame
!   EM2(1-3,1-NSCAT)=unit scattering polarization vector 2 in Lab Frame
!   EM1R(1-3,1-NSCAT)=unit scattering polarization vector 1 in Target Frame
!   EM2R(1-3,1-NSCAT)=unit scattering polarization vector 2 in Target Frame
!   PHIN(1-NSCAT)=values of PHI for scattering directions 1-NSCAT in
!                 Lab Frame
!                 [only used if JPBC=0: scattering by finite target]
!   AKSR(1-3,1-NSCAT)=(kx,ky,kz)*d for NSCAT scattering directions
!                    in Target Frame
!   ORDERM(1-NSCAT)=scattering order M for case of 1-d or 2-d target
!   ORDERN(1-NSCAT)=scattering order N for case of 2-d target
!   S1111(1-NSCAT)=weighted sum of |f_11|^2 over previous orientations
!   S2121(1-NSCAT)=weighted sum of |f_21|^2 over previous orientations
!   CXE01(1-3)    =incident pol state 1 in Lab Frame
!   CXE02(1-3)    =                   2
!   CXF11(1-NSCAT)=f_11 for current orientation
!   CXF12(1-NSCAT)=f_12 for current orientation
!   CXF21(1-NSCAT)=f_21 for current orientation
!   CXF22(1-NSCAT)=f_22 for current orientation
!   QABS(2)       =Q_abs for two incident polarizations
!   QABSUM(2)     =weighted sum of Qabs over previous orientations
!   QBKSCA(2)     =Q_bk for two incident polarizations
!   QBKSUM(2)     =weighted sum of Q_bk over previous orientations
!   QEXSUM(2)     =weighted sum of Q_ext over previous orientations
!   QEXT(2)       =Q_ext for two incident polarizations
!   QPHA(2)       =Q_pha for two incident polarizations
!   QPHSUM(2)     =weighted sum of Q_ph over previous orientations
!   QSCAG(3,2)    =Q_pr vector for two incident polarizations
!   QSCAT(2)      =Q_sca for two incident polarizations
!   QSCGSUM(3,2)  =weighted sum of Q_pr vector over prev. orients.
!   QSCSUM(2)     =weighted sum of Q_sca over previous orientations
!   QTRQAB(3,2)   =absorptive part of Q_gamma for two inc. pols.
!   QTRQABSUM(3,2)=weighted sum of Q_gamma over previous orientations
!   QTRQSC(3,2)   =scattering part of Q_gamma for two inc. pols.
!   QTRQSCSUM(3,2)=weighted sum of Q_gamma over previous orientations
!   WGTA(ITHETA,IPHI)=weight factor for theta,phi orientation
!   WGTB(IBETA)      =weight factor for beta orientation
!   SMORI(1-NSCAT,4,4)=weighted sum of Mueller matrix over previous
!                      orientations
! If IORTH=1, returns:
!   S1111(1-NSCAT)=updated weighted sum of |f_11|^2 over NSCAT
!                  scattering directions
!   S2121(1-NSCAT)=updated weighted sum of |f_21|^2
!   CX1121(1-NSCAT)=updated weighted summ of f_11*conjg(f_21)

! If IORTH=2, returns:
!   SM(1-NSCAT,4,4)=Mueller matrix for current orientation, and NSCAT
!                   scattering directions
!   CXS1(1-NSCAT)=amplitude scattering matrix element S_1 for current
!                 orientation
!   CXS2(1-NSCAT)=amplitude scattering matrix element S_2 for current
!                 orientation
!   CXS3(1-NSCAT)=amplitude scattering matrix element S_3 for current
!                 orientation
!   CXS4(1-NSCAT)=amplitude scattering matrix element S_4 for current
!                 orientations
!   QABSUM(2)     =updated weighted sum of Q_abs
!   QBKSUM(2)     =updated weighted sum of Q_bk
!   QEXSUM(2)     =updated weighted sum of Q_ext
!   QSCGSUM(3,2)  =updated weighted sum of Q_pr vector
!   QSCSUM(2)     =updated weighted sum of Q_sca
!   QTRQABSUM(3,2)=updated weighted sum of absorptive contribution to
!                  Q_gamma
!   QTRQSCSUM(3,2)=updated weighted sum of scattering contribution to
!                  Q_gamma

!   SMORI(1-NSCAT,4,4)=updated weighted sum of Mueller matrix elements

! History:
! 96.11.11 (BTD) ready for use
! 97.07.24 (BTD) corrected error in computation of amplitude scattering
!                matrix elements (CXS1,CXS2,CXS3,CXS4) and Mueller matri
!                elements SM.  Until now code mistakenly used PHI
!                defining target orientation rather than PHIN
!                defining scattering direction.  Correction of error
!                required replacing PHI by PHIN in argument list,
!                thereby requiring corresponding change in DDSCAT.
!                All calculated Mueller matrix elements until now were
!                therefore incorrect except for special cases where
!                PHI and PHIN were coincidentally the same (e.g.
!                target orientation PHI=0 and scattering plane PHIN=0).
!                This error was pointed out by Henriette Lemke (Inst.
!                for Atmospheric Physics, GKSS - Research Center,
!                Geesthacht.
! 97.07.24 (BTD) also corrected sign mistake in evaluation of S_14
!                and typo in evaluation of S_31.  Note that since S_31
!                had been evaluated incorrectly, PPOL (evaluated by
!                subroutine WRITESCA) was then numerically incorrect.
! 98.11.01 (BTD) Timo Nousiainen (tpnousia@lumi.meteo.helsinki.fi)
!                reported a problem with S_12 computed by DDSCAT for
!                a pseudosphere target when the scattering plane has
!                e.g. phi=90.  Problem is real and must be tracked down
!                and corrected.
! 98.11.16 (BTD) Modified code calculating complex amplitude scattering
!                matrix elements S1,S2,S3,S4 from f_ij.
! 00.06.13 (BTD) make data conversions explicit
! 03.10.23 (BTD) add QSCAG2 and QSCAG2SUM to argument list
!                add code to update QSCAG2SUM
! 06.10.07 (BTD) extend to handle 1-d and 2-d periodic targets
!                * add JPBC to arg list
!                * add PYDDX and PZDDX to arg list
! 06.12.24 (BTD) ver 7.0.1
!                * redefine normalization of Mueller matrix for 1d targe
! 06.12.25 (BTD) * for 2d targets, add contribution of incident wave to
!                  forward scattering if M=N=0
! 07.06.22 (BTD) ver 7.0.2
!                correct calculation of S_{ij} for targets with 1-d
!                periodicity
! 07.07.03 (BTD) introduce new treatment for evaluation of CXS1,CXS2,
!                CXS3,CXS4 when JPBC > 0 (previous treatment was not
!                correct).
!                add ENSC to argument list to support this new approach
! 07.07.08 (BTD) added EM1,EM2 to argument list
!                modified code to determine S matrix in case of
!                forward scattering (use EM1,EM2)
! 07.09.13 (BTD) Cosmetic changes
! 07.10.08 (BTD) Rewritten
! 07.10.09 (BTD) Added EM1R,EM2R to argument list
! 07.10.24 (BTD) Added ENSCR to argument list
! 07.12.29 (BTD) Experimentally determined correct phase relationship
!                between incident wave and scattered wave in case
!                of forward scattering from target that is periodic
!                in 2 dimensions.
! 08.05.29 (BTD) ver7.0.6
!                * changed declaration of arguments
!                  SM(MXSCA,4,4)    -> SM(4,4,MXSCA)
!                  SMORI(MXSCA,4,4) -> SMORI(4,4,MXSCA)
!                * made corresponding reordering of indices in lines
!                  evaluating SM and SMORI
! end history

! Copyright (C) 1996,1997,1998,2000,2003,2006,2007
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!*********  Sum scattering properties over orientations ****************

! Correspondence to notation of Bohren & Huffman (1983):
! Our f_jk correspond to notation of Draine (1988).
! For the special case where incident polarization modes 1 and 2 are
! parallel and perpendicular to the scattering plane, we have the
! relationship between the f_jk and the elements S_j of the amplitude
! scattering matrix as defined by Bohren & Huffman (1983) is as follows:
!     S_1 = -i*f_22
!     S_2 = -i*f_11
!     S_3 =  i*f_12
!     S_4 =  i*f_21
! In this case (incident pol. states parallel,perpendicular to the
! scattering plane), the elements of the 4x4 Mueller matrix (see Bohren
! Huffman 1983) are related to our f_jk as follows:
!     S_11 = [|f_11|^2 + |f_22|^2 + |f_12|^2 + |f_21|^2]/2
!     S_12 = [|f_11|^2 + |f_21|^2 - |f_22|^2 - |f_12|^2]/2
!     S_13 = -Re[f_11*conjg(f_12) + f_22*conjg(f_21)]
!     S_14 =  Im[f_22*conjg(f_21) - f_11*conjg(f_12)]
!     S_21 = [|f_11|^2 + |f_12|^2 - |f_22|^2 - |f_21|^2]/2
!     S_22 = [|f_11|^2 + |f_22|^2 - |f_21|^2 - |f_12|^2]/2
!     S_23 =  Re[f_22*conjg(f_21) - f_11*conjg(f_12)]
!     S_24 = -Im[f_11*conjg(f_12) + f_22*conjg(f_21)]
!     S_31 = -Re[f_11*conjg(f_21) + f_22*conjg(f_22)]
!     S_32 =  Re[f_22*conjg(f_12) - f_22*conjg(f_21)]
!     S_33 =  Re[f_22*conjg(f_11) + f_12*conjg(f_21)]
!     S_34 =  Im[f_11*conjg(f_22) + f_21*conjg(f_22)]
!     S_41 = -Im[f_21*conjg(f_11) + f_22*conjg(f_12)]
!     S_42 =  Im[f_22*conjg(f_12) - f_22*conjg(f_11)]
!     S_43 =  Im[f_22*conjg(f_11) - f_12*conjg(f_21)]
!     S_44 =  Re[f_22*conjg(f_11) - f_12*conjg(f_21)]

! Notation internal to this program:
! ij.ne.kl: CX_ijkl = \f_ij* \times f_kl  (summed over orientations)
!            S_ijij = \f_ij* \times f_ij  (summed over orientations)
! 4 pure real elements: S1111,S1212,S2121,S2222
! 8 complex elements:   CX1112,CX1121,CX1122,CX1221,CX1222,CX2122
! 4 redundant elts.:    CX1211 = Conjg(CX1112)
!                       CX2111 = Conjg(CX1121)
!                       CX2112 = Conjg(CX1221)
!                       CX2211 = Conjg(CX1122)

! In this routine we will carry out fully general calculation of the
! 4x4 Mueller matrix elements, valid for arbitrary incident polarization
! states j=1,2 for which we have the scattering matrix f_ij

! First, we compute the amplitude scattering matrix elements S_1,S_2,
! S_3,S_4 as defined by Bohren & Huffman:
!=======================================================================
      CXI=(0._WP,1._WP)
      PI=4._WP*ATAN(1._WP)

! WG=weighting factor for this orientation

      WG=WGTA(ITHETA,IPHI)*WGTB(IBETA)
      IF(IORTH==2)THEN

! CXE01,CXE02 are unit incident polarization vectors in Lab Frame.
! We assume CXE01 and CXE02 are normalized (REAPAR makes sure of this).
! CXA,CXB,CXC,CXD are complex coefficients corresponding to dot product
! of complex polarization vectors CXE01,CXE02 with y and z unit vectors
! in Lab Frame.

         CXA=CONJG(CXE01(2))
         CXB=CONJG(CXE01(3))
         CXC=CONJG(CXE02(2))
         CXD=CONJG(CXE02(3))

!*** diagnostic
!         write(0,*)'getmueller ckp a'
!         write(0,*)'cxa=',cxa
!         write(0,*)'cxb=',cxb
!         write(0,*)'cxc=',cxc
!         write(0,*)'cxd=',cxd
!***

! Compute complex scattering amplitudes S_1,S_2,S_3,S_4 for this
! particular target orientation and NSCAT scattering directions:

         IF(JPBC==0)THEN
            DO ND=1,NSCAT
               SINPHI=SIN(PHIN(ND))
               COSPHI=COS(PHIN(ND))
               CXS1(ND)=CXI*(CXF21(ND)*(CXA*SINPHI-CXB*COSPHI)+ &
                             CXF22(ND)*(CXC*SINPHI-CXD*COSPHI))
               CXS2(ND)=-CXI*(CXF11(ND)*(CXB*SINPHI+CXA*COSPHI)+ &
                              CXF12(ND)*(CXD*SINPHI+CXC*COSPHI))
               CXS3(ND)=-CXI*(CXF11(ND)*(CXA*SINPHI-CXB*COSPHI)+ &
                              CXF12(ND)*(CXC*SINPHI-CXD*COSPHI))
               CXS4(ND)=CXI*(CXF21(ND)*(CXB*SINPHI+CXA*COSPHI)+ &
                             CXF22(ND)*(CXD*SINPHI+CXC*COSPHI))
            ENDDO

         ELSEIF(JPBC==1.OR.JPBC==2)THEN

! JPBC = 1 or 2:
! ENSC(1-3,ND) = components of scattering unit vector in Lab Frame
! CXA,CXB = y,z components of incident polarization state 1 in Lab Frame
! CXC,CXD = y,z                                           2 in Lab Frame

            DO ND=1,NSCAT
               FAC=SQRT(ENSC(2,ND)**2+ENSC(3,ND)**2)
               IF(FAC<=1.E-3_WP)THEN

! if FAC = 0, then scattering is in either forward (theta=0) or
!                  backward (theta=pi) directions

                  FAC=SQRT(EM1(2,ND)**2+EM2(2,ND)**2)
                  FAC=FAC*SQRT(REAL(CXA)**2+REAL(CXB)**2)
                  COSPHI=REAL(CXA*EM1(2,ND)+CXB*EM1(3,ND))/FAC
                  SINPHI=-REAL(CXA*EM2(2,ND)+CXB*EM2(3,ND))/FAC

               ELSE

! if scattering angle is neither 0 nor pi:

                  COSPHI=REAL(CXA*ENSC(2,ND)+CXB*ENSC(3,ND))/FAC
                  SINPHI=REAL(CXC*ENSC(2,ND)+CXD*ENSC(3,ND))/FAC

               ENDIF
               CXS1(ND)=CXI*(CXF21(ND)*(CXA*SINPHI-CXB*COSPHI)+ &
                             CXF22(ND)*(CXC*SINPHI-CXD*COSPHI))
               CXS2(ND)=-CXI*(CXF11(ND)*(CXB*SINPHI+CXA*COSPHI)+ &
                              CXF12(ND)*(CXD*SINPHI+CXC*COSPHI))
               CXS3(ND)=-CXI*(CXF11(ND)*(CXA*SINPHI-CXB*COSPHI)+ &
                              CXF12(ND)*(CXC*SINPHI-CXD*COSPHI))
               CXS4(ND)=CXI*(CXF21(ND)*(CXB*SINPHI+CXA*COSPHI)+ &
                             CXF22(ND)*(CXD*SINPHI+CXC*COSPHI))
            ENDDO

         ELSEIF(JPBC==3)THEN

! JPBC = 3
!     Target is periodic in y and z directions in Lab Frame
!     Allowed scattering directions are identified by scattering
!     order (M,N).  
!     For each allowed (M,N) there is one forward (transmitted) 
!     direction and one backward (reflected) direction.
!     Following convention set in subroutine PBCSCAVEC for JPBC=3, 
!     it is assumed that NSCAT is even, with 
!     the first NSCAT/2 directions corresponding to transmission, and
!     the remaining NSCAT/2 directions corresponding to reflection

! CXA,CXB = y,z components of incident polarization state 1 in Lab Frame
! CXC,CXD = y,z                                           2 in Lab Frame
! ENSC(1-3,ND) = components of scattering unit vector in Lab Frame
! EM1(1-3,ND) = components of scattering polarization 1 in Lab Frame
! EM2(1-3,ND) = components of scattering polarization 2 in Lab Frame

            DO ND=1,NSCAT

               FAC=SQRT(ENSC(2,ND)**2+ENSC(3,ND)**2)

!*** diagnostic
!               write(0,7013)nd,ensc(1,nd),ensc(2,nd),ensc(3,nd),fac
! 7013 format('getmueller checkpoint baker: nd=',i2,/, &
!             'ensc(1-3,nd)=',3f8.5,' fac=',f10.7)
!*** end diagnostic

               IF(FAC<1.E-3_WP)THEN

! either zero-deg forward scattering, or 180 deg backscattering:

                  FAC=SQRT(EM1(2,ND)**2+EM2(2,ND)**2)
                  FAC=FAC*SQRT(REAL(CXA)**2+REAL(CXB)**2)
                  COSPHI=REAL(CXA*EM1(2,ND)+CXB*EM1(3,ND))/FAC
! 080111 (BTD) change
!                  SINPHI=-REAL(CXA*EM2(2,ND)+CXB*EM2(3,ND))/FAC
                  SINPHI=REAL(CXA*EM2(2,ND)+CXB*EM2(3,ND))/FAC

!*** diagnostic
!                  write(0,*)'in getmueller, ckpt ashanti, IBETA=',IBETA
!                  write(0,*)'nd=',nd,' --- zero deg forward scattering ---'
!                  write(0,7700)ensc(1,nd),ensc(2,nd),ensc(3,nd),    &
!                               em1(1,nd),em1(2,nd),em1(3,nd),       &
!                               em2(1,nd),em2(2,nd),em2(3,nd),       &
!                               enscr(1,nd),enscr(2,nd),enscr(3,nd), &
!                               em1r(1,nd),em1r(2,nd),em1r(3,nd),    &
!                               em2r(1,nd),em2r(2,nd),em2r(3,nd)
!                  write(0,7701)cxf11(nd),cxf21(nd),cxf12(nd),cxf22(nd), &
!                               cxa,cxb,cxc,cxd
!                  write(0,7702)cosphi,sinphi
! 7700 format('ensc(1-3,nd) =',3f9.5,' [Lab Frame]',/,     &
!             'em1(1-3,nd)  =',3f9.5,' [Lab Frame]',/,     &
!             'em2(1-3,nd)  =',3f9.5,' [Lab Frame]',/,     &
!             'enscr(1-3,nd)=',3f9.5,' [Target Frame]',/,  &
!             'em1r(1-3,nd) =',3f9.5,' [Target Frame]',/,  &
!             'em2r(1-3,nd) =',3f9.5,' [Target Frame]')
! 7701 format('cxf11(nd)=',2f10.6,/,'cxf21(nd)=',2f10.6,/,   &
!             'cxf12(nd)=',2f10.6,/,'cxf22(nd)=',2f10.6,/,   &
!             'cxa=',2f10.6,/,'cxb=',2f10.6,/,'cxc=',2f10.6,/,'cxd=',2f10.6)
! 7702 format('cosphi=',f10.6,' sinphi=',f10.6)
!*** end diagnostic

               ELSE

! scattering plane is well-defined:

                  COSPHI=REAL(CXA*ENSC(2,ND)+CXB*ENSC(3,ND))/FAC
                  SINPHI=REAL(CXC*ENSC(2,ND)+CXD*ENSC(3,ND))/FAC
!*** diagnostic
!                  write(0,*)'in getmueller, ckpt zulu'
!                  write(0,*)'nd=',nd,' --- not forward scattering ---'
!                  write(0,7700)ensc(1,nd),ensc(2,nd),ensc(3,nd),    &
!                               em1(1,nd),em1(2,nd),em1(3,nd),       &
!                               em2(1,nd),em2(2,nd),em2(3,nd),       &
!                               enscr(1,nd),enscr(2,nd),enscr(3,nd), &
!                               em1r(1,nd),em1r(2,nd),em1r(3,nd),    &
!                               em2r(1,nd),em2r(2,nd),em2r(3,nd)
!                  write(0,7701)cxf11(nd),cxf21(nd),cxf12(nd),cxf22(nd), &
!                               cxa,cxb,cxc,cxd
!                  write(0,7702)cosphi,sinphi
!*** end diagnostic

               ENDIF

               CXS1(ND)=CXI*(CXF21(ND)*(CXA*SINPHI-CXB*COSPHI)+ &
                             CXF22(ND)*(CXC*SINPHI-CXD*COSPHI))
               CXS2(ND)=-CXI*(CXF11(ND)*(CXB*SINPHI+CXA*COSPHI)+ &
                              CXF12(ND)*(CXD*SINPHI+CXC*COSPHI))
               CXS3(ND)=-CXI*(CXF11(ND)*(CXA*SINPHI-CXB*COSPHI)+ &
                              CXF12(ND)*(CXC*SINPHI-CXD*COSPHI))
               CXS4(ND)=CXI*(CXF21(ND)*(CXB*SINPHI+CXA*COSPHI)+ &
                             CXF22(ND)*(CXD*SINPHI+CXC*COSPHI))
            ENDDO

!*** diagnostic
!            if(nd.lt.20)then
!               write(0,*)'ensc(1-3,nd)=',ensc(1,nd),ensc(2,nd),ensc(3,nd)
!               write(0,*)'fac=',fac
!               write(0,*)'sinphi=',sinphi
!               write(0,*)'cosphi=',cosphi
!               write(0,*)'cxf11(nd)=',cxf11(nd)
!               write(0,*)'cxf21(nd)=',cxf21(nd)
!               write(0,*)'cxf12(nd)=',cxf12(nd)
!               write(0,*)'cxf22(nd)=',cxf22(nd)
!               write(0,9700)nd,cxs1(nd),cxs2(nd),cxs3(nd),cxs4(nd)
! 9700          format(i2,' s1,s2,s3,s4=',1p,  &
!                      2e11.3,1x,2e11.3,1x,2e11.3,1x,2e11.3)
!               fac=cxs1(nd)*conjg(cxs1(nd))
!               write(0,*)'|S1|^2=',fac
!               fac=cxs2(nd)*conjg(cxs2(nd))
!               write(0,*)'|S2|^2=',fac
!               fac=cxs3(nd)*conjg(cxs3(nd))
!               write(0,*)'|S3|^2=',fac
!               fac=cxs4(nd)*conjg(cxs4(nd))
!               write(0,*)'|S4|^2=',fac
!            endif
!***
         ENDIF

! if JPBC=1,2,3 (target periodic in y, z, or y and z directions):
! Check for special case: JPBC=3 (2-d target) and forward scattering wit
! ORDERM=ORDERN=0.  For this case we need to add incident wave to
! S_1 and S_2.
! Note that S_1 and S_2 are defined such that for M=N=0 we change
! iS_1 -> iS_1 +1  and iS_2 -> iS_2 +1
!  S_1 ->  S_1 -i       S_2 ->  S_2 -i
! note that CXS1 and CXS2 have yet to be multiplied by factor
! 2*pi/(ak1*aksr(1,nd)*pyddx*pzddx)

         IF(JPBC==3)THEN
!*** diagnostic
!            write(0,9089)ibeta
!            do nd=1,nscat
!               write(0,9090)nd,em1(1,nd),em1(2,nd),em1(3,nd),        &
!                            em2(1,nd),em2(2,nd),em2(3,nd),           &
!                            em1r(1,nd),em1r(2,nd),em1r(3,nd),        &
!                            em2r(1,nd),em2r(2,nd),em2r(3,nd),        &
!                            cxf11(nd),cxf21(nd),cxf12(nd),cxf22(nd), &
!                            cosphi,sinphi,cxs1(nd),cxs2(nd),         &
!                            cxs3(nd),cxs4(nd),cxa,cxb,cxc,cxd
!            enddo
! 9089 format('------------ in getmueller, ibeta=',i3,' ------------------')
! 9090 format(' nd=',i2,/,                            &
!             'em1=',3f10.6,' in lab frame',/,        &
!             'em2=',3f10.6,' in lab frame',/,        &
!             'em1r=',3f10.6,' in TF',/,              &
!             'em2r=',3f10.6,' in TF',/,              &
!             'cxf11=',2F12.8,/,'cxf21=',2F12.8,/,    &
!             'cxf12=',2F12.8,/,'cxf22=',2F12.8,/,    &
!             'cosphi=',f10.6,' sinphi=',f10.6,/,     &
!             'cxs1=',2f12.8,/,'cxs2=',2f12.8,/,      &
!             'cxs3=',2f12.8,/,'cxs4=',2f12.8,/,      &
!             'cxa=',2f13.8,/,'cxb=',2f13.8,/,        &
!             'cxc=',2f13.8,/,'cxd=',2f13.8)
!*** end diagnostic

            DO ND=1,NSCAT/2
               IF(NINT(ORDERM(ND))==0.AND.NINT(ORDERN(ND))==0)THEN
                  FAC=AK1*ABS(AKSR(1,ND))*PYDDX*PZDDX/(2._WP*PI)
                  CXTRM1= FAC*(COSPHI*(CXD*COSPHI-CXC*SINPHI)  &
                               +SINPHI*(-CXB*COSPHI+CXA*SINPHI))
                  CXTRM2= FAC*(COSPHI*(-CXA*COSPHI-CXB*SINPHI) &
                               -SINPHI*(CXC*COSPHI+CXD*SINPHI))

!*** begin diagnostic
!
!                  write(0,9050)nd,cxs1(nd),cxtrm1,(cxs1(nd)+cxtrm1)
!                  write(0,9051)cxs2(nd),cxtrm2,(cxs2(nd)+cxtrm2)
! 9050 format('nd --- radiated cxs --- ',   &
!                 '----- direct ------- ',  &
!                 '------ sum ---------',/, &
!              i2,1x,2f10.6,1x,2f10.6,1x,2f10.6)
! 9051 format(3x,2f10.6,1x,2f10.6,1x,2f10.6)
!*** end diagnostic

                  CXS1(ND)=CXS1(ND)+CXTRM1
                  CXS2(ND)=CXS2(ND)+CXTRM2
               ENDIF
!*** diagnostic
!                  write(0,9091)nd,aksr(1,nd),  &
!                  cxa,cxb,cxc,cxd,cosphi,sinphi,cxs1(nd),cxs2(nd), &
!                  cxs3(nd),cxs4(nd)
!                  nd2=nd+nscat/2
!                  write(0,9091)nd2,aksr(1,nd2),  &
!                  cxa,cxb,cxc,cxd,cosphi,sinphi,cxs1(nd2),cxs2(nd2), &
!                  cxs3(nd2),cxs4(nd2)
! 9091 format('in getmueller: nd=',i2,            &
!             ' aksr(1,nd)=',f8.5,/,                           &
!             'cxa=',2f8.5,/,                                  &
!             'cxb=',2f8.5,/,                                  &
!             'cxc=',2f8.5,/,                                  &
!             'cxd=',2f8.5,/,                                  &
!             'cosphi=',f8.5,' sinphi=',f8.5,/,                &
!             'cxs1=',2f12.8,/,'cxs2=',2f12.8,/,               &
!             'cxs3=',2f12.8,/,'cxs4=',2f12.8)
!*** end diagnostic
            ENDDO
         ENDIF

! Calculation of scattering amplitudes CXS1, CXS2, CXS3, CXS4
! is complete.
! Now compute Mueller matrix elements for this particular target
! orientation:

         DO ND=1,NSCAT
            SM(1,1,ND)=0.5_WP*REAL(CXS1(ND)*CONJG(CXS1(ND))+ &
                                   CXS2(ND)*CONJG(CXS2(ND))+ &
                                   CXS3(ND)*CONJG(CXS3(ND))+ &
                                   CXS4(ND)*CONJG(CXS4(ND)))
            SM(1,2,ND)=0.5_WP*REAL(CXS2(ND)*CONJG(CXS2(ND))- &
                                   CXS1(ND)*CONJG(CXS1(ND))+ &
                                   CXS4(ND)*CONJG(CXS4(ND))- &
                                   CXS3(ND)*CONJG(CXS3(ND)))
            SM(1,3,ND)=REAL(CXS2(ND)*CONJG(CXS3(ND))+CXS1(ND)*CONJG(CXS4(ND)))
            SM(1,4,ND)=AIMAG(CXS2(ND)*CONJG(CXS3(ND))- &
                             CXS1(ND)*CONJG(CXS4(ND)))
            SM(2,1,ND)=0.5_WP*REAL(CXS2(ND)*CONJG(CXS2(ND))- &
                                   CXS1(ND)*CONJG(CXS1(ND))- &
                                   CXS4(ND)*CONJG(CXS4(ND))+ &
                                   CXS3(ND)*CONJG(CXS3(ND)))
            SM(2,2,ND)=0.5_WP*REAL(CXS2(ND)*CONJG(CXS2(ND))+ &
                                   CXS1(ND)*CONJG(CXS1(ND))- &
                                   CXS4(ND)*CONJG(CXS4(ND))- &
                                   CXS3(ND)*CONJG(CXS3(ND)))
            SM(2,3,ND)=REAL(CXS2(ND)*CONJG(CXS3(ND))-CXS1(ND)*CONJG(CXS4(ND)))
            SM(2,4,ND)=AIMAG(CXS2(ND)*CONJG(CXS3(ND))+ &
                             CXS1(ND)*CONJG(CXS4(ND)))
            SM(3,1,ND)=REAL(CXS2(ND)*CONJG(CXS4(ND))+CXS1(ND)*CONJG(CXS3(ND)))
            SM(3,2,ND)=REAL(CXS2(ND)*CONJG(CXS4(ND))-CXS1(ND)*CONJG(CXS3(ND)))
            SM(3,3,ND)=REAL(CXS1(ND)*CONJG(CXS2(ND))+CXS3(ND)*CONJG(CXS4(ND)))
            SM(3,4,ND)=AIMAG(CXS2(ND)*CONJG(CXS1(ND))+ &
                             CXS4(ND)*CONJG(CXS3(ND)))
            SM(4,1,ND)=AIMAG(CXS4(ND)*CONJG(CXS2(ND))+ &
                             CXS1(ND)*CONJG(CXS3(ND)))
            SM(4,2,ND)=AIMAG(CXS4(ND)*CONJG(CXS2(ND))- &
                             CXS1(ND)*CONJG(CXS3(ND)))
            SM(4,3,ND)=AIMAG(CXS1(ND)*CONJG(CXS2(ND))- &
                             CXS3(ND)*CONJG(CXS4(ND)))
            SM(4,4,ND)=REAL(CXS1(ND)*CONJG(CXS2(ND))-CXS3(ND)*CONJG(CXS4(ND)))
         ENDDO

         IF(JPBC>0)THEN

            IF(JPBC==1)THEN

! prepare to calculate S^{(1d)} for target with 1-d periodicity in y

               FAC0=2._WP*PI/(AK1*PYDDX**2)

            ELSEIF(JPBC==2)THEN

! prepare to calculate S^{(1d)} for target with 1-d periodicity in z

               FAC0=2._WP*PI/(AK1*PZDDX**2)

            ELSEIF(JPBC==3)THEN

! Prepare to calculate S^{(2d)} for target with periodicity in y and z

               FAC0=(2._WP*PI/(AK1*PYDDX*PZDDX))**2
            ENDIF

! k*sin(alpha_s)=component of k_s perpendicular to target

            DO ND=1,NSCAT
               IF(JPBC==1)THEN
                  FAC=FAC0/SQRT(AKSR(1,ND)**2+AKSR(3,ND)**2)
               ELSEIF(JPBC==2)THEN
                  FAC=FAC0/SQRT(AKSR(1,ND)**2+AKSR(2,ND)**2)
               ELSEIF(JPBC==3)THEN
                  FAC=FAC0/AKSR(1,ND)**2
               ELSE
                  CALL ERRMSG('FATAL','GETMUELLER',' invalid JPBC')
               ENDIF

               DO K=1,4
                  DO J=1,4
                     SM(J,K,ND)=FAC*SM(J,K,ND)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

! Augment Mueller matrix for orientational average

         DO ND=1,NSCAT
            DO K=1,4
               DO J=1,4
                  SMORI(J,K,ND)=SMORI(J,K,ND)+WG*SM(J,K,ND)
               ENDDO
            ENDDO
         ENDDO

!***********************************************************************

      ELSE

! When IORTH=1, cannot compute full Mueller matrix.
! Compute three scattering properties:

         DO ND=1,NSCAT
            S1111(ND)=S1111(ND)+WG*REAL(CONJG(CXF11(ND))*CXF11(ND))
            S2121(ND)=S2121(ND)+WG*REAL(CONJG(CXF21(ND))*CXF21(ND))
            CX1121(ND)=CX1121(ND)+WG*REAL(CONJG(CXF11(ND))*CXF21(ND))
         ENDDO
      ENDIF

!**** Now augment sums of QABS,QBKSCA,QEXT,QSCA,G*QSCA over incident
!     directions and polarizations.

      DO JO=1,IORTH
         QEXSUM(JO)=QEXSUM(JO)+QEXT(JO)*WG
         QABSUM(JO)=QABSUM(JO)+QABS(JO)*WG
         QBKSUM(JO)=QBKSUM(JO)+QBKSCA(JO)*WG
         QPHSUM(JO)=QPHSUM(JO)+QPHA(JO)*WG
         QSCSUM(JO)=QSCSUM(JO)+QSCAT(JO)*WG
         QSCG2SUM(JO)=QSCG2SUM(JO)+QSCAG2(JO)*WG
         DO J=1,3
            QSCGSUM(J,JO)=QSCGSUM(J,JO)+QSCAG(J,JO)*WG
         ENDDO
         IF(CMDTRQ=='DOTORQ')THEN
            DO J=1,3
               QTRQABSUM(J,JO)=QTRQABSUM(J,JO)+QTRQAB(J,JO)*WG
               QTRQSCSUM(J,JO)=QTRQSCSUM(J,JO)+QTRQSC(J,JO)*WG
            ENDDO
         ENDIF
      ENDDO
!*** diagnostic
!      write(0,*)'returning from getmueller'
!***
      RETURN
    END SUBROUTINE GETMUELLER
