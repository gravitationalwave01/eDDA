    SUBROUTINE GETMUELLER(IBETA,IORTH,IPHI,ITHETA,JPBC,MXBETA,MXSCA,MXPHI,   &
         MXTHET,NSCAT,ORDERM,ORDERN,CMDTRQ,AK1,AKS_TF,ENSC_LF,ENSC_TF,PYDDX, &
         PZDDX,PHIN,SM,SMORI,S1111,S2121,CX1121,CXE01_LF,CXE02_LF,CXE01_TF,  &
         CXE02_TF,CXF11,CXF12,CXF21,CXF22,CXS1,CXS2,CXS3,CXS4,QABS,QABSUM,   &
         QBKSCA,QBKSUM,QEXSUM,QEXT,QPHA,QPHSUM,QSCAG,QSCAG2,QSCAT,QSCG2SUM,  &
         QSCGSUM,QSCSUM,QTRQAB,QTRQABSUM,QTRQSC,QTRQSCSUM,WGTA,WGTB,EM1_LF,  &
         EM2_LF,EM1_TF,EM2_TF)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

!                       getmueller v5
! Arguments:

      INTEGER :: IBETA,IORTH,IPHI,ITHETA,JPBC,MXBETA,MXSCA,MXPHI,MXTHET,NSCAT
      CHARACTER :: CMDTRQ*6
      REAL(WP) :: AK1,PYDDX,PZDDX,WG
      REAL(WP) ::            &
         AKS_TF(3,MXSCA),    &
         EM1_LF(3,MXSCA),    &
         EM1_TF(3,MXSCA),    &
         EM2_LF(3,MXSCA),    &
         EM2_TF(3,MXSCA),    &
         ENSC_LF(3,MXSCA),   &
         ENSC_TF(3,MXSCA),   &
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
         CXE01_LF(3),   &
         CXE02_LF(3),   &
         CXE01_TF(3),   &
         CXE02_TF(3),   &
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
      COMPLEX(WP) :: CXA,CXAA,CXB,CXBB,CXC,CXCC,CXD,CXDD, &
                     CXFAC0,CXFAC1,CXFAC2,CXI

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
!   ENSC_LF(1-3,1-NSCAT)=normalized scattering direction vector in Lab Frame
!                     for scattering directions 1-NSCAT
!   ENSC_TF(1-3,1-NSCAT)=normalized scattering direction vector in Target
!                      Frame for scattering directions 1-NSCA
!   EM1_LF(1-3,1-NSCAT)=unit scattering polarization vector 1 in Lab Frame
!   EM2_LF(1-3,1-NSCAT)=unit scattering polarization vector 2 in Lab Frame
!   EM1_TF(1-3,1-NSCAT)=unit scattering polarization vector 1 in Target Frame
!   EM2_TF(1-3,1-NSCAT)=unit scattering polarization vector 2 in Target Frame
!   PHIN(1-NSCAT)=values of PHI for scattering directions 1-NSCAT in
!                 Lab Frame
!                 [only used if JPBC=0: scattering by finite target]
!   AKS_TF(1-3,1-NSCAT)=(kx,ky,kz)*d for NSCAT scattering directions
!                    in Target Frame
!   ORDERM(1-NSCAT)=scattering order M for case of 1-d or 2-d target
!   ORDERN(1-NSCAT)=scattering order N for case of 2-d target
!   S1111(1-NSCAT)=weighted sum of |f_11|^2 over previous orientations
!   S2121(1-NSCAT)=weighted sum of |f_21|^2 over previous orientations
!   CXE01_LF(1-3)    =incident pol state 1 in Lab Frame
!   CXE02_LF(1-3)    =                   2
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
!                * redefine normalization of Mueller matrix for 1d target
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
! 11.11.04 (BTD) ddscat ver7.1.1
!                v3
!                * change notation: 
!                  AKSR  -> AKS_TF
!                  CXE01 -> CXE01_LF
!                  CXE02 -> CXE02_LF
!                  ENSC  -> ENSC_LF
!                  ENSCR -> ENSC_TF
!                  EM1   -> EM1_LF
!                  EM2   -> EM2_LF
!                  EM1R  -> EM1_TF
!                  EM2R  -> EM2_TF
!                * for JPBC=3, change ENSC_LF -> ENSC_TF in calculation
!                  of SINPHI and COSPHI for transformation of incident
!                  polarization states CXE01_LF, CXE02_LF 
!                  to parallel or perp to scattering plane
! 11.11.15 (BTD) * add CXE01_TF,CXE02_TF to argument list here and in
!                  DDSCAT
! 11.11.16 (BTD) v4
!                * finally got the transformation from f_ij to S_i
!                  correct!
!                * cleaned up code
! 12.04.25 (BTD) v5
!                * corrected error in computation of S_i and S_ij
!                  for JPBC=0
! end history

! Copyright (C) 1996,1997,1998,2000,2003,2006,2007,2008,2011,2012
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

! Notation:
!        khat_0     = unit vector along incident direction
!        ehat_01    = CXE01 = complex incident polarization vector 1
!        ehat_02    = CXE02 = complex incident polarization vector 2
!        ehat_0perp = incident pol vector perp to scattering plane
!                   = scattered pol vector perp to scattering plane
!                   = em2
!        ehat_0para = incident pol vector para to scattering plane
!                   = khat_0 cross ehat_0perp   [conventional def.]
!                   = xhat_LF cross ehat_0perp
!                   = xhat_LF cross [yhat_LF (yhat_LF dot ehat_0perp)
!                                    + zhat_LF (zhat_LF dot ehat_0perp)]
!                   = - yhat_LF (ehat_0perp dot zhat_LF)
!                     + zhat_LF (ehat_0perp dot yhat_LF)
!                   = - yhat_LF em2_LF(3) + zhat_LF em2_LF(2)
!        em2        = scattered pol vector perp to scattering plane
!        cxa        = conjg(ehat_01) dot yhat_LF
!        cxb        = conjg(ehat_01) dot zhat_LF
!        cxc        = conjg(ehat_02) dot yhat_LF
!        cxd        = conjg(ehat_02) dot zhat_LF

! WG=weighting factor for this orientation

      WG=WGTA(ITHETA,IPHI)*WGTB(IBETA)
      IF(IORTH==2)THEN

! Compute complex scattering amplitudes S_1,S_2,S_3,S_4 for this
! particular target orientation and NSCAT scattering directions:

         IF(JPBC==0)THEN
            CXA=CONJG(CXE01_LF(2))
            CXB=CONJG(CXE01_LF(3))
            CXC=CONJG(CXE02_LF(2))
            CXD=CONJG(CXE02_LF(3))
            DO ND=1,NSCAT

! 2012.04.25 (BTD) replace ----------------------------------------------
!               CXA=-CXE01_LF(2)*EM2_LF(3,ND)+CXE01_LF(3)*EM2_LF(2,ND)
!               CXB=CXE01_LF(2)*EM2_LF(2,ND)+CXE01_LF(3)*EM2_LF(3,ND)
!               CXC=-CXE02_LF(2)*EM2_LF(3,ND)+CXE02_LF(3)*EM2_LF(2,ND)
!               CXD=CXE02_LF(2)*EM2_LF(2,ND)+CXE02_LF(3)*EM2_LF(3,ND)
!               CXS1(ND)=CXI*(CXF11(ND)*CONJG(CXA)+CXF12(ND)*CONJG(CXC))
!               CXS2(ND)=CXI*(CXF21(ND)*CONJG(CXB)+CXF22(ND)*CONJG(CXD))
!               CXS3(ND)=CXI*(CXF11(ND)*CONJG(CXB)+CXF12(ND)*CONJG(CXD))
!               CXS4(ND)=CXI*(CXF21(ND)*CONJG(CXA)+CXF22(ND)*CONJG(CXC))

               COSPHI=EM2_LF(3,ND)
               SINPHI=-EM2_LF(2,ND)
               CXAA=CXA*COSPHI+CXB*SINPHI
               CXBB=CXB*COSPHI-CXA*SINPHI
               CXCC=CXC*COSPHI+CXD*SINPHI
               CXDD=CXD*COSPHI-CXC*SINPHI
               CXS1(ND)=-CXI*(CXF21(ND)*CXBB+CXF22(ND)*CXDD)
               CXS2(ND)=-CXI*(CXF11(ND)*CXAA+CXF12(ND)*CXCC)
               CXS3(ND)=CXI*(CXF11(ND)*CXBB+CXF12(ND)*CXDD)
               CXS4(ND)=CXI*(CXF21(ND)*CXAA+CXF22(ND)*CXCC)
!-----------------------------------------------------------------------
!*** diagnostic
!               if(nd.eq.10)then
!                  write(0,*)'getmueller_v5 ckpt 2, nd=',nd
!                  write(0,fmt='(a,i2,a,3f9.5)')'nd=',nd,' ensc_lf(1-3,nd)=', &
!                     ensc_lf(1,nd),ensc_lf(2,nd),ensc_lf(3,nd)
!                  write(0,fmt='(a,1pe10.3,1pe11.3,a)')' s1 = (',cxs1(nd),')'
!                  write(0,fmt='(a,1pe10.3,1pe11.3,a)')' s2 = (',cxs2(nd),')'
!                  write(0,fmt='(a,1pe10.3,1pe11.3,a)')' s3 = (',cxs3(nd),')'
!                  write(0,fmt='(a,1pe10.3,1pe11.3,a)')' s4 = (',cxs4(nd),')'
!               endif
!***
            ENDDO

         ELSEIF(JPBC==1.OR.JPBC==2)THEN

! JPBC = 1 or 2:
! ENSC_LF(1-3,ND) = components of scattering unit vector in Lab Frame
! CXA,CXB = y,z components of incident polarization state 1 in Lab Frame
! CXC,CXD = y,z                                           2 in Lab Frame

            CXA=CONJG(CXE01_LF(2))
            CXB=CONJG(CXE01_LF(3))
            CXC=CONJG(CXE02_LF(2))
            CXD=CONJG(CXE02_LF(3))
            DO ND=1,NSCAT

! 2012.04.25 (BTD) replace ----------------------------------------------
!               CXA=-CXE01_LF(2)*EM2_LF(3,ND)+CXE01_LF(3)*EM2_LF(2,ND)
!               CXB=CXE01_LF(2)*EM2_LF(2,ND)+CXE01_LF(3)*EM2_LF(3,ND)
!               CXC=-CXE02_LF(2)*EM2_LF(3,ND)+CXE02_LF(3)*EM2_LF(2,ND)
!               CXD=CXE02_LF(2)*EM2_LF(2,ND)+CXE02_LF(3)*EM2_LF(3,ND)
!               CXS1(ND)=CXI*(CXF11(ND)*CONJG(CXA)+CXF12(ND)*CONJG(CXC))
!               CXS2(ND)=CXI*(CXF21(ND)*CONJG(CXB)+CXF22(ND)*CONJG(CXD))
!               CXS3(ND)=CXI*(CXF11(ND)*CONJG(CXB)+CXF12(ND)*CONJG(CXD))
!               CXS4(ND)=CXI*(CXF21(ND)*CONJG(CXA)+CXF22(ND)*CONJG(CXC))

               COSPHI=EM2_LF(3,ND)
               SINPHI=-EM2_LF(2,ND)
               CXAA=CXA*COSPHI+CXB*SINPHI
               CXBB=CXB*COSPHI-CXA*SINPHI
               CXCC=CXC*COSPHI+CXD*SINPHI
               CXDD=CXD*COSPHI-CXC*SINPHI
               CXS1(ND)=-CXI*(CXF21(ND)*CXBB+CXF22(ND)*CXDD)
               CXS2(ND)=-CXI*(CXF11(ND)*CXAA+CXF12(ND)*CXCC)
               CXS3(ND)=CXI*(CXF11(ND)*CXBB+CXF12(ND)*CXDD)
               CXS4(ND)=CXI*(CXF21(ND)*CXAA+CXF22(ND)*CXCC)
!-----------------------------------------------------------------------
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

! ENSC_LF(1-3,ND) = components of scattering unit vector in Lab Frame
! EM1_LF(1-3,ND) = components of scattering polarization 1 in Lab Frame
! EM2_LF(1-3,ND) = components of scattering polarization 2 in Lab Frame

            CXA=CONJG(CXE01_LF(2))
            CXB=CONJG(CXE01_LF(3))
            CXC=CONJG(CXE02_LF(2))
            CXD=CONJG(CXE02_LF(3))
            DO ND=1,NSCAT
! 2012.04.25 (BTD) replace ----------------------------------------------
!               CXA=-CXE01_LF(2)*EM2_LF(3,ND)+CXE01_LF(3)*EM2_LF(2,ND)
!               CXB=CXE01_LF(2)*EM2_LF(2,ND)+CXE01_LF(3)*EM2_LF(3,ND)
!               CXC=-CXE02_LF(2)*EM2_LF(3,ND)+CXE02_LF(3)*EM2_LF(2,ND)
!               CXD=CXE02_LF(2)*EM2_LF(2,ND)+CXE02_LF(3)*EM2_LF(3,ND)

!*** diagnostic
!               write(0,*)'getmueller_v5 ckpt 2.9, nd=',nd
!               write(0,fmt='(a,2f9.5,a,2f9.5)') &
!                   'cxe01_lf(2)=',cxe01_lf(2),' cxe01_lf(3)=',cxe01_lf(3)
!               write(0,fmt='(a,f9.5,a,f9.5)') &
!                   'em1_lf(2)=',em1_lf(2,nd),' em1_lf(3)=',em1_lf(3,nd)
!               write(0,fmt='(a,f9.5,a,f9.5)') &
!                   'em2_lf(2)=',em2_lf(2,nd),' em2_lf(3)=',em2_lf(3,nd)
!               write(0,fmt='(a,2f9.5)')'cxfac1=',cxfac1
!               write(0,fmt='(a,2f9.5)')'cxfac2=',cxfac2
!***
!*** diagnostic
!               write(0,*)'getmueller_v5 ckpt 3, IBETA=',IBETA
!               write(0,*)'nd=',nd,' --- zero deg forward scattering ---'
!               write(0,7700)ensc_lf(1,nd),ensc_lf(2,nd),ensc_lf(3,nd), &
!                            em1_lf(1,nd),em1_lf(2,nd),em1_lf(3,nd),    &
!                            em2_lf(1,nd),em2_lf(2,nd),em2_lf(3,nd),    &
!                            ensc_tf(1,nd),ensc_tf(2,nd),ensc_tf(3,nd), &
!                            em1_tf(1,nd),em1_tf(2,nd),em1_tf(3,nd),    &
!                            em2_tf(1,nd),em2_tf(2,nd),em2_tf(3,nd),    &
!                            cxa,cxb,cxc,cxd,cxe01_lf,cxe02_lf
!               write(0,7701)cxf11(nd),cxf21(nd),cxf12(nd),cxf22(nd)
!               write(0,7702)cxfac1,cxfac2
! 7700 format('ensc_lf(1-3,nd)=',3f9.5,' [Lab Frame]',/,    &
!             'em1_lf(1-3,nd) =',3f9.5,' [Lab Frame]',/,    &
!             'em2_lf(1-3,nd) =',3f9.5,' [Lab Frame]',/,    &
!             'ensc_tf(1-3,nd)=',3f9.5,' [Target Frame]',/, &
!             'em1_tf(1-3,nd) =',3f9.5,' [Target Frame]',/, &
!             'em2_tf(1-3,nd) =',3f9.5,' [Target Frame]',/, &
!             '  cxa   = (',f10.6,',',f10.6,')',/,   &
!             '  cxb   = (',f10.6,',',f10.6,')',/,   &
!             '  cxc   = (',f10.6,',',f10.6,')',/,   &
!             '  cxd   = (',f10.6,',',f10.6,')',/,   &
!             'cxe01_LF(1-3,nd)=(',f9.5,',',f9.5,')(', &
!             f9.5,',',f9.5,')(',f9.5,',',f9.5')',/,   &
!             'cxe02_LF(1-3,nd)=(',f9.5,',',f9.5,')(', &
!             f9.5,',',f9.5,')(',f9.5,',',f9.5')')
! 7701 format('cxf11(nd)=',2f10.6,/,'cxf21(nd)=',2f10.6,/, &
!             'cxf12(nd)=',2f10.6,/,'cxf22(nd)=',2f10.6)
! 7702 format('cxfac1=',2f10.6,/,'cxfac2=',2f10.6)
!*** end diagnostic

!               CXS1(ND)=CXI*(CXF11(ND)*CONJG(CXA)+CXF12(ND)*CONJG(CXC))
!               CXS2(ND)=CXI*(CXF21(ND)*CONJG(CXB)+CXF22(ND)*CONJG(CXD))
!               CXS3(ND)=CXI*(CXF11(ND)*CONJG(CXB)+CXF12(ND)*CONJG(CXD))
!               CXS4(ND)=CXI*(CXF21(ND)*CONJG(CXA)+CXF22(ND)*CONJG(CXC))

               COSPHI=EM2_LF(3,ND)
               SINPHI=-EM2_LF(2,ND)
               CXAA=CXA*COSPHI+CXB*SINPHI
               CXBB=CXB*COSPHI-CXA*SINPHI
               CXCC=CXC*COSPHI+CXD*SINPHI
               CXDD=CXD*COSPHI-CXC*SINPHI
               CXS1(ND)=-CXI*(CXF21(ND)*CXBB+CXF22(ND)*CXDD)
               CXS2(ND)=-CXI*(CXF11(ND)*CXAA+CXF12(ND)*CXCC)
               CXS3(ND)=CXI*(CXF11(ND)*CXBB+CXF12(ND)*CXDD)
               CXS4(ND)=CXI*(CXF21(ND)*CXAA+CXF22(ND)*CXCC)
!-------------------------------------------------------------------------
!*** diagnostic
!               if(nd.lt.20)then
!                  write(0,*)'getmueller_v5 ckpt 5'
!                  write(0,fmt='(a,i2,a,3f9.5)')'nd=',nd,' ensc_lf(1-3,nd)=', &
!                     ensc_lf(1,nd),ensc_lf(2,nd),ensc_lf(3,nd)
!                  write(0,fmt='(a,f10.6)')'fac=',fac
!                  write(0,fmt='(a,2f10.6)')'cxf11(nd)=',cxf11(nd)
!                  write(0,fmt='(a,2f10.6)')'cxf21(nd)=',cxf21(nd)
!                  write(0,fmt='(a,2f10.6)')'cxf12(nd)=',cxf12(nd)
!                  write(0,fmt='(a,2f10.6)')'cxf22(nd)=',cxf22(nd)
!                  write(0,fmt='(a,1pe10.3,1pe11.3,a)')' s1 = (',cxs1(nd),')'
!                  write(0,fmt='(a,1pe10.3,1pe11.3,a)')' s2 = (',cxs2(nd),')'
!                  write(0,fmt='(a,1pe10.3,1pe11.3,a)')' s3 = (',cxs3(nd),')'
!                  write(0,fmt='(a,1pe10.3,1pe11.3,a)')' s4 = (',cxs4(nd),')'
!                  fac=cxs1(nd)*conjg(cxs1(nd))
!                  write(0,'(a,f12.8)')'|S1|^2=',fac
!                  fac=cxs2(nd)*conjg(cxs2(nd))
!                  write(0,'(a,f12.8)')'|S2|^2=',fac
!                  fac=cxs3(nd)*conjg(cxs3(nd))
!                  write(0,'(a,f12.8)')'|S3|^2=',fac
!                  fac=cxs4(nd)*conjg(cxs4(nd))
!                  write(0,'(a,f12.8)')'|S4|^2=',fac
!               endif
!***
            ENDDO   ! loop over ND

         ENDIF   ! JPBC==3

! if JPBC=1,2,3 (target periodic in y, z, or y and z directions):
! Check for special case: JPBC=3 (2-d target) and forward scattering wit
! ORDERM=ORDERN=0.  For this case we need to add incident wave to
! S_1 and S_2.
! Note that S_1 and S_2 are defined such that for M=N=0 we change
! iS_1 -> iS_1 +1  and iS_2 -> iS_2 +1
!  S_1 ->  S_1 -i       S_2 ->  S_2 -i
! note that CXS1 and CXS2 have yet to be multiplied by factor
! 2*pi/(ak1*aks_tf(1,nd)*pyddx*pzddx)

         IF(JPBC==3)THEN
!*** diagnostic
!            write(0,9089)ibeta,cxe01_lf,cxe02_lf,cxe01_tf,cxe02_tf
!            do nd=1,nscat
!--------
!               write(0,9090)nd,em1_lf(1,nd),em1_lf(2,nd),em1_lf(3,nd), &
!                            em2_lf(1,nd),em2_lf(2,nd),em2_lf(3,nd),    &
!                            em1_tf(1,nd),em1_tf(2,nd),em1_tf(3,nd),    &
!                            em2_tf(1,nd),em2_tf(2,nd),em2_tf(3,nd),    &
!                            cxfac1,cxfac2,                             &
!                            cxf11(nd),cxf21(nd),cxf12(nd),cxf22(nd),   &
!                            cxs1(nd),cxs2(nd),           &
!                            cxs3(nd),cxs4(nd)
!            enddo
! 9089 format('---------- getmueller_v5 ckpt 6, ibeta=',i3,' -------_-------',/, & 
!             'cxe01_LF =',3('(',f10.6,',',f10.6,') '),/,                      &
!             'cxe02_LF =',3('(',f10.6,',',f10.6,') '),/,                      &
!             'cxe01_TF =',3('(',f10.6,',',f10.6,') '),/,                      &
!             'cxe02_TF =',3('(',f10.6,',',f10.6,') '))
! 9090 format(' nd=',i2,/,                            &
!             'em1_lf=',3f10.6,' in lab frame',/,     &
!             'em2_lf=',3f10.6,' in lab frame',/,     &
!             'em1_tf=',3f10.6,' in TF',/,            &
!             'em2_tf=',3f10.6,' in TF',/,            &
!             'cxfac1=a=(',f10.6,',',f10.6,')',/,     &
!             'cxfac2=b=(',f10.6,',',f10.6,')',/,     &
!             'cxf11=',2F12.8,/,'cxf21=',2F12.8,/,    &
!             'cxf12=',2F12.8,/,'cxf22=',2F12.8,/,    &
!             'cxs1=',2f12.8,/,'cxs2=',2f12.8,/,      &
!             'cxs3=',2f12.8,/,'cxs4=',2f12.8)
!*** end diagnostic

            DO ND=1,NSCAT/2
               IF(NINT(ORDERM(ND))==0.AND.NINT(ORDERN(ND))==0)THEN

! zero-degree forward-scattering
! need to add incident wave to radiated wave
! 2012.04.26 (BTD) change
!                  CXFAC0=AKS_TF(1,ND)*AK1*PYDDX*PZDDX/(2._WP*PI)
                  CXFAC0=ABS(AKS_TF(1,ND))*AK1*PYDDX*PZDDX/(2._WP*PI)
!*** diagnostic
!                  write(0,fmt='(a,2f11.7)') &
!                     'getmueller_v5 ckpt 7:   s1=',cxs1(nd)
!                  write(0,fmt='(a,2f11.7)') &
!                     '                     s2=',cxs2(nd)
!                  write(0,fmt='(a,2f11.7)') &
!                     '                     s3=',cxs3(nd)
!                  write(0,fmt='(a,2f11.7)') &
!                     '                     s4=',cxs4(nd)
!                  write(0,fmt='(a,2f11.7)') &
!                     '               cxfac0  =',cxfac0
!***
! 2012.04.25 (BTD): for reasons not understood, incident wave
!                   added for S2 has different sign from S1
!                   ... does this indicate a sign mistake in
!                       calculation of either S1 or S2?
!                   ... is it possible that there is an error in eq. 68 of
!                       Draine & Flatau 2008?. Something like this could
!                       arise from a 180deg shift in directions between
!                       the incident and scattered "parallel" basis vectors,
!                       or between the incident and scattered "perpendicular"
!                       basis vectors (should check for this in output files!)
!                   This appears to indicate that for JPBC=3 the
!                   sign of either CXS1 or CXS2 is incorrect
!                   The Mueller matrix elements S11, S22, S12, and S21 appear
!                   to be correct, but they would not be affected by
!                   a sign error in either CXS1 or CXS2
!                   Elements S13,S14,S23,S24,S31,S32,S33,S34,S41,S42,S43,S44
!                   would be sensitive to a sign error, and it would be
!                   a good idea to find a way to test them.

                  CXS1(ND)=CXS1(ND)-CXFAC0
                  CXS2(ND)=CXS2(ND)+CXFAC0
!---
!*** diagnostic
!                  write(0,fmt='(a,2f11.7)') &
!                     'getmueller_v5 ckpt 7.1: s1=',cxs1(nd)
!                  write(0,fmt='(a,2f11.7)') &
!                     '                     s2=',cxs2(nd)
!*** end diagnostic
               ENDIF
            ENDDO   ! loop over ND
         ENDIF   ! JPBC==3

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
                  FAC=FAC0/SQRT(AKS_TF(1,ND)**2+AKS_TF(3,ND)**2)
               ELSEIF(JPBC==2)THEN
                  FAC=FAC0/SQRT(AKS_TF(1,ND)**2+AKS_TF(2,ND)**2)
               ELSEIF(JPBC==3)THEN
                  FAC=FAC0/AKS_TF(1,ND)**2
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

!*** diagnostic
!         write(0,fmt='(a,7x,a,7x,a,7x,a,7x,a,7x,a,7x,a,2(/,20x,6f11.7))')     &
!            'getmueller_v5 ckpt 8','S_11','S_12','S_21','S_22','S_31','S_41', &
!            SM(1,1,1),SM(1,2,1),SM(2,1,1),SM(2,2,1),SM(3,1,1),SM(4,1,1),      &
!            SM(1,1,2),SM(1,2,2),SM(2,1,2),SM(2,2,2),SM(3,1,2),SM(4,1,2)
!***
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
