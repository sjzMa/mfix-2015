!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_MU_g                                               C
!  Purpose: Calculate the effective viscosity for a turbulent flow,    C
!           which is the sum of molecular and eddy viscosities         C
!                                                                      C
!  Author: W. Sams/M. Syamlal                         Date: 18-JUL-94  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: MFIX 2.0 mods (previous name CALC_MU_gt)                   C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: allow SI unit                                              C
!  Author: S. Dartevelle                              Date: 01-Jul-02  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: compute turbulent eddy viscosity                           C
!  Author: S. Benyahia                                Date: May-13-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!     Perry, R. H., and Chilton, C. H., Chemical Engineers' Handbook,  C
!        5th Edition, McGraw-Hill Inc., 1973, pp. 248, eqn. 3-133.     C
!     Arnold, J. H., Vapor viscosities and the Sutherland equation,    C
!        Journal of Chemical Physics, 1 (2), 1933, pp. 170-176.        C
!     Sutherland, W., The Viscosity of Gases and Molecular Force,      C
!        Phil. Mag. 5:507-531, 1893.                                   C
!                                                                      C
!     Cao, J. and Ahmadi, G., 1995, Gas-particle two-phase turbulent   C
!        flow in a vertical duct. Int. J. Multiphase Flow, vol. 21,    C
!        No. 6, pp. 1203-1228.                                         C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_MU_G()

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE physprop
      USE geometry
      USE fldvar
      USE visc_g
      USE visc_s
      USE indices
      USE constant
      USE toleranc
      USE compar
      USE drag
      USE run          !S. Dartevelle
      USE turb
      USE sendrecv
      USE mms
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      DOUBLE PRECISION, PARAMETER :: F2O3 = 2.D0/3.D0
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Cell indices
      INTEGER :: I, J, K, IM, JM, KM
      INTEGER :: IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
      INTEGER :: IMJPK, IMJMK, IMJKP, IMJKM, IPJKM, IPJMK
      INTEGER :: IJMKP, IJMKM, IJPKM
! Solids phase index
      INTEGER :: M
! Strain rate tensor components for mth solids phase
      DOUBLE PRECISION :: D_g(3,3)
! Second invariant of the deviator of D_g
      DOUBLE PRECISION :: I2_devD_g
! Constant in turbulent viscosity formulation
      DOUBLE PRECISION :: C_MU
! particle relaxation time
      DOUBLE PRECISION :: Tau_12_st
! U_g at the north face of the THETA cell-(i, j+1/2, k)
      DOUBLE PRECISION :: U_g_N
! U_g at the south face of the THETA cell-(i, j-1/2, k)
      DOUBLE PRECISION :: U_g_S
! U_g at the top face of the THETA cell-(i, j, k+1/2)
      DOUBLE PRECISION :: U_g_T
! U_g at the bottom face of the THETA cell-(i, j, k-1/2)
      DOUBLE PRECISION :: U_g_B
! U_g at the center of the THETA cell-(i, j, k)
! Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION :: U_g_C
! V_g at the east face of the THETA cell-(i+1/2, j, k)
      DOUBLE PRECISION :: V_g_E
! V_g at the west face of the THETA cell-(i-1/2, j, k)
      DOUBLE PRECISION :: V_g_W
! V_g at the top face of the THETA cell-(i, j, k+1/2)
      DOUBLE PRECISION :: V_g_T
! V_g at the bottom face of the THETA cell-(i, j, k-1/2)
      DOUBLE PRECISION :: V_g_B
! W_g at the east face of the THETA cell-(i+1/2, j, k)
      DOUBLE PRECISION :: W_g_E
! W_g at the west face of the THETA cell-(1-1/2, j, k)
      DOUBLE PRECISION :: W_g_W
! W_g at the north face of the THETA cell-(i, j+1/2, k)
      DOUBLE PRECISION :: W_g_N
! W_g at the south face of the THETA cell-(i, j-1/2, k)
      DOUBLE PRECISION :: W_g_S
! W_g at the center of the THETA cell-(i, j, k).
! Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION :: W_g_C
!-----------------------------------------------

! JFD: Calling SET_EP_FACTORS here to make sure EPG_IFAC is defined
!      when calc_mu_g is called for the firt time.
!      There may be a better way to do this...
      CALL SET_EP_FACTORS

! solids phase index used throughout routine...
! may be inappropriate for multiple solids phases
      M = 1 ! for solids phase


!!$omp parallel do private(ijk) schedule(dynamic,chunk_size)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN


! Gas viscosity   (in Poise or Pa.s)
! Calculating gas viscosity using Sutherland's formula with
! Sutherland's constant (C) given by Vogel's equation C = 1.47*Tb.
! For air  C = 110 (Tb=74.82)
!         mu = 1.71*10-4 poise at T = 273K

            IF (MU_G0 == UNDEFINED) MU_G(IJK) = to_SI*1.7D-4 * &
               (T_G(IJK)/273.0D0)**1.5D0 * (383.D0/(T_G(IJK)+110.D0))

            MU_GT(IJK) = MU_G(IJK)
            LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)
! if ishii then multiply by void fraction otherwise multiply by 1
            EPMU_GT(IJK)= EPG_IFAC(IJK)*MU_GT(IJK)
            EPLAMBDA_GT(IJK) = EPG_IFAC(IJK)*LAMBDA_GT(IJK)

! K_epsilon model
! ---------------------------------------------------------------->>>
            IF (K_Epsilon) THEN

               C_MU = 9D-02

! I'm not very confident about this correction in Peirano paper, but it's made
! available here, uncomment to use it. sof@fluent.com --> 02/01/05
!               IF(SIMONIN .AND. F_GS(IJK,1) > SMALL_NUMBER) THEN
!                  Tau_12_st = Ep_s(IJK,M)*RO_S(IJK,M)/F_GS(IJK,1)
!                  X_21 = Ep_s(IJK,M)*RO_S(IJK,M)/(EP_g(IJK)*RO_g(IJK))
! new definition of C_mu (equation A.12, Peirano et al. (2002) Powder tech. 122,69-82)
!                  IF( K_12(IJK)/(2.0D0*K_Turb_G(IJK)) < ONE) &
!                     C_MU = C_MU/(ONE+ 0.314D0*X_21*Tau_12_st / Tau_1(IJK) * &
!                            (ONE - K_12(IJK)/(2.0D0*K_Turb_G(IJK))) )
!               ENDIF
! Correction in Ahmadi paper (Cao and Ahmadi)
               IF(AHMADI .AND. F_GS(IJK,1) > SMALL_NUMBER) THEN
                  Tau_12_st = Ep_s(IJK,M)*RO_S(IJK,M)/F_GS(IJK,1)
                  C_MU = C_MU/(ONE+ Tau_12_st/Tau_1(IJK) * &
                         (EP_s(IJK,M)/(ONE-EP_star_array(IJK)))**3)
               ENDIF

! Definition of the turbulent viscosity
               MU_GT(IJK) = MU_G(IJK) + RO_G(IJK)*C_MU*&
                  K_Turb_G(IJK)**2 / (E_Turb_G(IJK) + SMALL_NUMBER)
               MU_GT(IJK) = MIN(MU_GMAX, MU_GT(IJK))
               LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)
! if ishii then multiply by void fraction otherwise multiply by 1
               EPMU_GT(IJK) = EPG_IFAC(IJK)*MU_GT(IJK)
               EPLAMBDA_GT(IJK) = EPG_IFAC(IJK)*LAMBDA_GT(IJK)
            ENDIF
! ----------------------------------------------------------------<<<

         ELSE
            MU_G(IJK)  = ZERO
            MU_GT(IJK) = ZERO
            LAMBDA_GT(IJK) = ZERO
            EPMU_GT(IJK) = ZERO
            EPLAMBDA_GT(IJK) = ZERO
         ENDIF   ! end if (fluid_at(ijk))

      ENDDO   ! end do (ijk=ijkstart3,ijkend3)


! MMS: Force constant gas viscosity at all cells including ghost cells.
      IF (USE_MMS) THEN
         DO IJK = ijkstart3, ijkend3
            MU_G(IJK) = MU_G0
            MU_GT(IJK) = MU_G(IJK)
            LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)
! if ishii then multiply by void fraction otherwise multiply by 1
            EPMU_GT(IJK) = EPG_IFAC(IJK)*MU_GT(IJK)
            EPLAMBDA_GT(IJK) = EPG_IFAC(IJK)*LAMBDA_GT(IJK)
         ENDDO
      END IF ! end if (USE_MMS)



! L_scale0 model
! ---------------------------------------------------------------->>>
!!$omp parallel do &
!!$omp$ schedule(dynamic,chunk_size) &
!!$omp$ private(IJK, I,J,K,IM,JM,KM, &
!!$omp& IMJK,IPJK,IJMK,IJPK,IJKM,IJKP,IMJPK,IMJMK,IMJKP, &
!!$omp& IMJKM,IPJKM,IPJMK,IJMKP,IJMKM,IJPKM, &
!!$omp& U_G_N,U_G_S,U_G_T,U_G_B,V_G_E,V_G_W,V_G_T,V_G_B, &
!!$omp$ W_G_N,W_G_S,W_G_E,W_G_W,  U_G_C,W_G_C, D_G,I2_DEVD_G )
      DO IJK = ijkstart3, ijkend3
         IF ( FLUID_AT(IJK) .AND. L_SCALE(IJK)/=ZERO) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IM = IM1(I)
            JM = JM1(J)
            KM = KM1(K)

            IMJK = IM_OF(IJK)
            IPJK = IP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKM = KM_OF(IJK)
            IJKP = KP_OF(IJK)
            IMJPK = IM_OF(IJPK)
            IMJMK = IM_OF(IJMK)
            IMJKP = IM_OF(IJKP)
            IMJKM = IM_OF(IJKM)
            IPJKM = IP_OF(IJKM)
            IPJMK = IP_OF(IJMK)
            IJMKP = JM_OF(IJKP)
            IJMKM = JM_OF(IJKM)
            IJPKM = JP_OF(IJKM)

! Find fluid velocity values at faces of the cell
            U_G_N = AVG_Y(AVG_X_E(U_G(IMJK),U_G(IJK),I),&
                          AVG_X_E(U_G(IMJPK),U_G(IJPK),I),J)        !i, j+1/2, k
            U_G_S = AVG_Y(AVG_X_E(U_G(IMJMK),U_G(IJMK),I),&
                          AVG_X_E(U_G(IMJK),U_G(IJK),I),JM)         !i, j-1/2, k
            U_G_T = AVG_Z(AVG_X_E(U_G(IMJK),U_G(IJK),I),&
                          AVG_X_E(U_G(IMJKP),U_G(IJKP),I),K)        !i, j, k+1/2
            U_G_B = AVG_Z(AVG_X_E(U_G(IMJKM),U_G(IJKM),I),&
                          AVG_X_E(U_G(IMJK),U_G(IJK),I),KM)         !i, j, k-1/2
            V_G_E = AVG_X(AVG_Y_N(V_G(IJMK),V_G(IJK)),&
                          AVG_Y_N(V_G(IPJMK),V_G(IPJK)),I)          !i+1/2, j, k
            V_G_W = AVG_X(AVG_Y_N(V_G(IMJMK),V_G(IMJK)),&
                          AVG_Y_N(V_G(IJMK),V_G(IJK)),IM)           !i-1/2, j, k
            V_G_T = AVG_Z(AVG_Y_N(V_G(IJMK),V_G(IJK)),&
                          AVG_Y_N(V_G(IJMKP),V_G(IJKP)),K)          !i, j, k+1/2
            V_G_B = AVG_Z(AVG_Y_N(V_G(IJMKM),V_G(IJKM)),&
                          AVG_Y_N(V_G(IJMK),V_G(IJK)),KM)           !i, j, k-1/2
            W_G_N = AVG_Y(AVG_Z_T(W_G(IJKM),W_G(IJK)),&
                          AVG_Z_T(W_G(IJPKM),W_G(IJPK)),J)          !i, j+1/2, k
            W_G_S = AVG_Y(AVG_Z_T(W_G(IJMKM),W_G(IJMK)),&
                          AVG_Z_T(W_G(IJKM),W_G(IJK)),JM)           !i, j-1/2, k
            W_G_E = AVG_X(AVG_Z_T(W_G(IJKM),W_G(IJK)),&
                          AVG_Z_T(W_G(IPJKM),W_G(IPJK)),I)          !i+1/2, j, k
            W_G_W = AVG_X(AVG_Z_T(W_G(IMJKM),W_G(IMJK)),&
                          AVG_Z_T(W_G(IJKM),W_G(IJK)),IM)           !i-1/2, j, k

            IF (CYLINDRICAL) THEN
               U_G_C = AVG_X_E(U_G(IMJK),U_G(IJK),I)  !i, j, k
               W_G_C = AVG_Z_T(W_G(IJKM),W_G(IJK))    !i, j, k
            ELSE
               U_G_C = ZERO
               W_G_C = ZERO
            ENDIF

! Find components of fluid phase strain rate tensor, D_g, at center of
! the cell - (i,j,k)
            D_G(1,1) = (U_G(IJK)-U_G(IMJK))*ODX(I)
            D_G(1,2) = HALF*((U_G_N - U_G_S)*ODY(J)+(V_G_E-V_G_W)*&
                       ODX(I))
            D_G(1,3) = HALF*((W_G_E - W_G_W)*ODX(I)+(U_G_T-U_G_B)*&
                       (OX(I)*ODZ(K))-W_G_C*OX(I))
            D_G(2,1) = D_G(1,2)
            D_G(2,2) = (V_G(IJK)-V_G(IJMK))*ODY(J)
            D_G(2,3) = HALF*((V_G_T-V_G_B)*(OX(I)*ODZ(K))+&
                       (W_G_N-W_G_S)*ODY(J))
            D_G(3,1) = D_G(1,3)
            D_G(3,2) = D_G(2,3)
            D_G(3,3) = (W_G(IJK)-W_G(IJKM))*(OX(I)*ODZ(K)) + U_G_C*OX(I)

! Calculate the second invariant of the deviator of D_g
            I2_DEVD_G = ((D_G(1,1)-D_G(2,2))**2+(D_G(2,2)-D_G(3,3))**2+&
                         (D_G(3,3)-D_G(1,1))**2)/6.D0 + &
                        D_G(1,2)**2 + D_G(2,3)**2 + D_G(3,1)**2

            MU_GT(IJK) = MIN(MU_GMAX, MU_G(IJK)+2.0*&
               L_SCALE(IJK)*L_SCALE(IJK)*RO_G(IJK)*SQRT(I2_DEVD_G))
            LAMBDA_GT(IJK) = -F2O3*MU_GT(IJK)
! if ishii then multiply by void fraction otherwise multiply by 1
            EPMU_GT(IJK) = EPG_IFAC(IJK)*MU_GT(IJK)
            EPLAMBDA_GT(IJK) = EPG_IFAC(IJK)*LAMBDA_GT(IJK)
         ENDIF ! end if (fluid_at(ijk) and l_scale(ijk)/=0))
      ENDDO    ! end loop (ijk=ijkstart3,ijkend3)
! end calculations for L_scale0 model
! ----------------------------------------------------------------<<<

      RETURN
      END SUBROUTINE CALC_MU_G

