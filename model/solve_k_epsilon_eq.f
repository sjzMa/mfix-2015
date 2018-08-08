!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOLVE_K_Epsilon_EQ(IER)                                C
!  Purpose: Solve K & Epsilon equations for a turbulent flow           C
!                                                                      C
!                                                                      C
!  Author: S. Benyahia                                Date: MAY-13-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SOLVE_K_Epsilon_EQ(IER)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      !USE matrix
      USE bc
      USE compar
      USE constant
      USE cutcell
      USE drag
      USE energy
      USE fldvar
      USE fun_avg
      USE functions
      USE geometry
      USE indices
      USE leqsol
      USE mflux
      USE output
      USE param
      USE param1
      USE pgcor
      USE physprop
      USE pscor
      USE residual
      USE run
      USE rxns
      USE toleranc
      USE turb
      USE ur_facs
      USE usr
      Use ambm
      Use tmp_array, S_p => Array1, S_c => Array2, EPs => Array3, VxGama => Array4

      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!
!                      phase index
      INTEGER          m, I, J, K, LC
!
      DOUBLE PRECISION apo

!                      temporary variables in residual computation
      DOUBLE PRECISION res1, mres1, num_res, den_res
      INTEGER          ires1
!
!                      Indices
      INTEGER          IJK
!
!                      linear equation solver method and iterations
      INTEGER          LEQM, LEQI

!                      A default zero flux will be defined for both K & Epsilon at walls
      DOUBLE PRECISION BC_hw_K_Turb_G (DIMENSION_BC),  BC_hw_E_Turb_G (DIMENSION_BC)
      DOUBLE PRECISION BC_K_Turb_GW (DIMENSION_BC),   BC_E_Turb_GW (DIMENSION_BC)
      DOUBLE PRECISION BC_C_K_Turb_G (DIMENSION_BC),  BC_C_E_Turb_G (DIMENSION_BC)
!
!                      small value of K or E, 1 cm2/s2 = 1e-4 m2/s2 = 1e-4 m2/s3
      DOUBLE PRECISION smallTheta
!
      character(LEN=8) :: Vname
!-----------------------------------------------

      IF( .NOT. K_Epsilon) RETURN

      call lock_ambm
      call lock_tmp_array

!
      smallTheta = (to_SI)**4 * ZERO_EP_S

      RESID(RESID_ke,0) = ZERO
      NUM_RESID(RESID_ke,0) = ZERO
      DEN_RESID(RESID_ke,0) = ZERO
      MAX_RESID(RESID_ke,0) = ZERO
      IJK_RESID(RESID_ke,0) = 0

! Setting default zero flux for K & Epsilon since we use wall functions.
! If an expert user want to use Low Re K-Epilon model and needs to set
! the turbulence quatities to zero at walls, then set the hw's to UNDEFINE will
! do it. All the variables below can be changed in the same way as in the
! MFIX data file in the boundary conditions section.
!
      DO LC = 1, DIMENSION_BC
        BC_hw_K_Turb_G (LC) = ZERO
        BC_hw_E_Turb_G (LC) = ZERO
        BC_K_Turb_GW (LC) = ZERO
        BC_E_Turb_GW (LC) = ZERO
        BC_C_K_Turb_G (LC) = ZERO
        BC_C_E_Turb_G (LC) = ZERO
      ENDDO
! End of setting default zero flux for K & Epsilon wall boundary conditions
!
! Equations solved for gas phase, thus M = 0
      M = 0
          CALL INIT_AB_M (A_M, B_M, IJKMAX2, M)

! Solve first fot the K_Turb_G Equation

          DO IJK = IJKSTART3, IJKEND3
!
             I = I_OF(IJK)
             J = J_OF(IJK)
             K = K_OF(IJK)
               IF (FLUID_AT(IJK)) THEN

                  APO = ROP_G(IJK)*VOL(IJK)*ODT
                  S_P(IJK) = APO + (  ZMAX(SUM_R_G(IJK)) &
                                    + K_Turb_G_p(IJK) )*VOL(IJK)
                  S_C(IJK) =   APO*K_Turb_GO(IJK) &
                            + K_Turb_G(IJK)*ZMAX((-SUM_R_G(IJK)))*VOL(IJK) &
                            + K_Turb_G_c(IJK) *VOL(IJK)
               ELSE
!
                  S_P(IJK) = ZERO
                  S_C(IJK) = ZERO
!
               ENDIF

            END DO

            IF(.NOT.ADDED_MASS) THEN
               CALL CONV_DIF_PHI (K_Turb_G, DIF_K_Turb_G, DISCRETIZE(9), &
                               U_G, V_G, W_G, Flux_gE, Flux_gN, Flux_gT, M, A_M, B_M)
            ELSE
               CALL CONV_DIF_PHI (K_Turb_G, DIF_K_Turb_G, DISCRETIZE(9), &
                               U_G, V_G, W_G, Flux_gSE, Flux_gSN, Flux_gST, M, A_M, B_M)
            ENDIF
!
!
            CALL BC_PHI (K_Turb_G, BC_K_Turb_G, BC_K_Turb_GW, BC_HW_K_Turb_G, &
                         BC_C_K_Turb_G, M, A_M, B_M)
!
!
            CALL SOURCE_PHI (S_P, S_C, EP_G, K_Turb_G, M, A_M, B_M)
!
            CALL CALC_RESID_S (K_Turb_G, A_M, B_M, M, num_res, den_res, res1, &
                               mres1, ires1, ZERO)

            RESID(RESID_ke,0) = RESID(RESID_ke,0)+res1
            NUM_RESID(RESID_ke,0) = NUM_RESID(RESID_ke,0)+num_res
            DEN_RESID(RESID_ke,0) = DEN_RESID(RESID_ke,0)+den_res
            if(mres1 .gt. MAX_RESID(RESID_ke,0))then
              MAX_RESID(RESID_ke,0) = mres1
              IJK_RESID(RESID_ke,0) = ires1
            endif
!
            CALL UNDER_RELAX_S (K_Turb_G, A_M, B_M, M, UR_FAC(9))
!
!          call check_ab_m(a_m, b_m, m, .false., ier)
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)
!          write(*,*) res1, mres1, &
!           ires1
!
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, DO_K, &
!          ier)
!
            CALL ADJUST_LEQ (RESID(RESID_ke,0), LEQ_IT(9), LEQ_METHOD(9), &
               LEQI, LEQM)
!
            write(Vname, '(A,I2)')'K_Turb_G'
            CALL SOLVE_LIN_EQ (Vname, 9, K_Turb_G, A_M, B_M, M, LEQI, LEQM, &
                             LEQ_SWEEP(9), LEQ_TOL(9),  LEQ_PC(9), IER)
!          call out_array(K_Turb_G, Vname)
!
! remove small negative K values generated by linear solver
! same as adjust_theta.f
!
           DO IJK = IJKSTART3, IJKEND3
            IF (FLUID_AT(IJK)) THEN
             IF(K_Turb_G(IJK) < smallTheta) K_Turb_G(IJK) = smallTheta
            ENDIF
           END DO

!!!!!!!!!!!!!
! Now, we'll solve for the E_Turb_G (dissipation) Equation.
!!!!!!!!!!!!!
!
! Initiate (again) the Am Bm matrix. This has to be done for every scalar equation.
          CALL INIT_AB_M (A_M, B_M, IJKMAX2, M)

          DO IJK = IJKSTART3, IJKEND3
!
             I = I_OF(IJK)
             J = J_OF(IJK)
             K = K_OF(IJK)
               IF (FLUID_AT(IJK)) THEN

                  APO = ROP_G(IJK)*VOL(IJK)*ODT
                  S_P(IJK) = APO + (  ZMAX(SUM_R_G(IJK)) &
                                    + E_Turb_G_p(IJK)  )*VOL(IJK)
                  S_C(IJK) =   APO*E_Turb_GO(IJK) &
                            + E_Turb_G(IJK)*ZMAX((-SUM_R_G(IJK)))*VOL(IJK) &
                            + E_Turb_G_c(IJK) *VOL(IJK)
               ELSE
!
                  S_P(IJK) = ZERO
                  S_C(IJK) = ZERO
!
               ENDIF

            END DO

            IF(.NOT.ADDED_MASS) THEN
               CALL CONV_DIF_PHI (E_Turb_G, DIF_E_Turb_G, DISCRETIZE(9), &
                               U_G, V_G, W_G, Flux_gE, Flux_gN, Flux_gT, M, A_M, B_M)
            ELSE
               CALL CONV_DIF_PHI (E_Turb_G, DIF_E_Turb_G, DISCRETIZE(9), &
                               U_G, V_G, W_G, Flux_gSE, Flux_gSN, Flux_gST, M, A_M, B_M)
            ENDIF
!
!
            CALL BC_PHI (E_Turb_G, BC_E_Turb_G, BC_E_Turb_GW, BC_HW_E_Turb_G, &
                         BC_C_E_Turb_G, M, A_M, B_M)
!
!
            CALL SOURCE_PHI (S_P, S_C, EP_G, E_Turb_G, M, A_M, B_M)
!
! When implementing the wall functions, The Epsilon (dissipation) value at the fluid cell
! near the walls needs to be set.

            DO IJK = IJKSTART3, IJKEND3
             I = I_OF(IJK)
             J = J_OF(IJK)
             K = K_OF(IJK)
!
               IF (FLUID_AT(IJK)) THEN
!
                  IF(WALL_AT(JP_OF(IJK)).OR.WALL_AT(JM_OF(IJK))) THEN
                     A_M(IJK,1,M) = ZERO
                     A_M(IJK,-1,M) = ZERO
                     A_M(IJK,2,M) = ZERO
                     A_M(IJK,-2,M) = ZERO
                     A_M(IJK,3,M) = ZERO
                     A_M(IJK,-3,M) = ZERO
                     A_M(IJK,0,M) = -ONE
                     B_M(IJK,M) =-((0.09D+0)**0.75*K_Turb_G(IJK)**1.5)/DY(J) &
                                 *2.0D+0/0.42D+0

                  ELSE IF(WALL_AT(KP_OF(IJK)).OR.WALL_AT(KM_OF(IJK))) THEN
                     A_M(IJK,1,M) = ZERO
                     A_M(IJK,-1,M) = ZERO
                     A_M(IJK,2,M) = ZERO
                     A_M(IJK,-2,M) = ZERO
                     A_M(IJK,3,M) = ZERO
                     A_M(IJK,-3,M) = ZERO
                     A_M(IJK,0,M) = -ONE
                     B_M(IJK,M) =-((0.09D+0)**0.75*K_Turb_G(IJK)**1.5)* &
                                  (ODZ(K)*OX(I)*2.0D+0)/0.42D+0
                  ENDIF  !for identifying wall cells in J or K direction

                  IF(CYLINDRICAL) THEN
                     IF (WALL_AT(IP_OF(IJK)))  THEN
                        A_M(IJK,1,M) = ZERO
                        A_M(IJK,-1,M) = ZERO
                        A_M(IJK,2,M) = ZERO
                        A_M(IJK,-2,M) = ZERO
                        A_M(IJK,3,M) = ZERO
                        A_M(IJK,-3,M) = ZERO
                        A_M(IJK,0,M) = -ONE
                        B_M(IJK,M) =-((0.09D+0)**0.75*K_Turb_G(IJK)**1.5)/DX(I) &
                                    *2.0D+0/0.42D+0

                     ENDIF! for wall cells in I direction

                  ELSE IF (WALL_AT(IP_OF(IJK)).OR.WALL_AT(IM_OF(IJK))) THEN
                     A_M(IJK,1,M) = ZERO
                     A_M(IJK,-1,M) = ZERO
                     A_M(IJK,2,M) = ZERO
                     A_M(IJK,-2,M) = ZERO
                     A_M(IJK,3,M) = ZERO
                     A_M(IJK,-3,M) = ZERO
                     A_M(IJK,0,M) = -ONE
                     B_M(IJK,M) =-((0.09D+0)**0.75*K_Turb_G(IJK)**1.5)/DX(I) &
                                   *2.0D+0/0.42D+0

                  ENDIF ! for cylindrical

                  IF(CUT_CELL_AT(IJK)) THEN
                     A_M(IJK,1,M) = ZERO
                     A_M(IJK,-1,M) = ZERO
                     A_M(IJK,2,M) = ZERO
                     A_M(IJK,-2,M) = ZERO
                     A_M(IJK,3,M) = ZERO
                     A_M(IJK,-3,M) = ZERO
                     A_M(IJK,0,M) = -ONE

                     B_M(IJK,M) =-((0.09D+0)**0.75*K_Turb_G(IJK)**1.5)/(0.42D+0*DELH_Scalar(IJK))
                  ENDIF

!
               ENDIF  !for fluid at ijk
            ENDDO

            CALL CALC_RESID_S (E_Turb_G, A_M, B_M, M, num_res, den_res, res1, &
                               mres1, ires1, ZERO)

            RESID(RESID_ke,0) = RESID(RESID_ke,0)+res1
            NUM_RESID(RESID_ke,0) = NUM_RESID(RESID_ke,0)+num_res
            DEN_RESID(RESID_ke,0) = DEN_RESID(RESID_ke,0)+den_res
            if(mres1 .gt. MAX_RESID(RESID_ke,0))then
              MAX_RESID(RESID_ke,0) = mres1
              IJK_RESID(RESID_ke,0) = ires1
            endif
!
            CALL UNDER_RELAX_S (E_Turb_G, A_M, B_M, M, UR_FAC(9))
!
!          call check_ab_m(a_m, b_m, m, .false., ier)
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)
!          write(*,*) res1, mres1, &
!           ires1
!
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, DO_K, &
!          ier)
!
            CALL ADJUST_LEQ (RESID(RESID_ke,0), LEQ_IT(9), LEQ_METHOD(9), &
               LEQI, LEQM)
!
            write(Vname, '(A,I2)')'E_Turb_G'
            CALL SOLVE_LIN_EQ (Vname, 9, E_Turb_G, A_M, B_M, M, LEQI, LEQM, &
                             LEQ_SWEEP(9), LEQ_TOL(9), LEQ_PC(9), IER)
!          call out_array(E_Turb_G, Vname)
!
! remove small negative Epsilon values generated by linear solver
! same as adjust_theta.f
!
           DO IJK = IJKSTART3, IJKEND3
            IF (FLUID_AT(IJK)) THEN
             IF(E_Turb_G(IJK) < smallTheta) E_Turb_G(IJK) = smallTheta
            ENDIF
           END DO

      call unlock_ambm
      call unlock_tmp_array

      RETURN
      END SUBROUTINE SOLVE_K_Epsilon_EQ


!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
