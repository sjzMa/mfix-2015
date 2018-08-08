!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_GRANULAR_ENERGY                                   C
!  Purpose: Solve granular energy equations in matrix equation         C
!     form Ax=b. The center coefficient (ap) and source vector (b)     C
!     are negative.  The off-diagonal coefficients are positive.       C
!                                                                      C
!  Purpose: Solve granular energy equations                            C
!                                                                      C
!                                                                      C
!  Author: K. Agrawal                                 Date:            C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOLVE_GRANULAR_ENERGY(IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE toleranc
      USE run
      USE physprop
      USE visc_s
      USE geometry
      USE fldvar
      USE constant
      USE output
      USE indices
      USE drag
      USE residual
      USE ur_facs
      USE pgcor
      USE pscor
      USE leqsol
      USE bc
      USE energy
      USE rxns
      Use ambm
      Use tmp_array, S_p => Array1, S_c => Array2, EPs => Array3
      USE compar
      USE mflux
      USE mpi_utility
      USE mms
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Error index
      INTEGER :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! phase index
      INTEGER :: M, L
! Cp * Flux
      DOUBLE PRECISION CpxFlux_E(DIMENSION_3), CpxFlux_N(DIMENSION_3), CpxFlux_T(DIMENSION_3)
! previous time step term
      DOUBLE PRECISION :: apo
! source terms which appear appear in the center coefficient (lhs) and
! the rhs vector
      DOUBLE PRECISION :: sourcelhs, sourcerhs
! Indices
      INTEGER :: IJK
! linear equation solver method and iterations
      INTEGER :: LEQM, LEQI
! particle mass and diameter
      DOUBLE PRECISION :: M_PM,D_PM
! total solids volume fraction
      DOUBLE PRECISION :: TOT_EPS(DIMENSION_3)
! net production of all solids phase
      DOUBLE PRECISION :: TOT_SUM_RS(DIMENSION_3)
! old value of total solids number density
      DOUBLE PRECISION :: TOT_NO(DIMENSION_3)
! small value of theta_m, 1 cm2/s2 = 1e-4 m2/s2
      DOUBLE PRECISION :: smallTheta
! temporary variable for volume x interphase transfer
! coefficient between solids phases
      DOUBLE PRECISION :: VxTC_ss(DIMENSION_3,DIMENSION_LM)

! Arrays for storing errors:
! 140 - Linear Eq diverged
! 141 - Unphysical values
! 14x - Unclassified
      INTEGER :: Err_l(0:numPEs-1)  ! local
      INTEGER :: Err_g(0:numPEs-1)  ! global

! temporary use of global arrays:
! array1 (locally s_p)
! source vector: coefficient of dependent variable
! becomes part of a_m matrix; must be positive
!      DOUBLE PRECISION :: S_P(DIMENSION_3)
! array2 (locally s_c)
! source vector: constant part becomes part of b_m vector
!      DOUBLE PRECISION :: S_C(DIMENSION_3)
! array3 (locally eps)
!      DOUBLE PRECISION :: eps(DIMENSION_3)

! Septadiagonal matrix A_m, vector b_m
!      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------

      call lock_ambm       ! locks arrays a_m and b_m
      call lock_tmp_array  ! locks array1,array2,array3
                           ! (locally s_p, s_c, eps)

! Initialize error flags.
      Err_l = 0

      smallTheta = (to_SI)**4 * ZERO_EP_S

      DO M = 0, MMAX
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M)
      ENDDO

      SELECT CASE (KT_TYPE_ENUM)
        CASE(LUN_1984, AHMADI_1995, SIMONIN_1996, GD_1999, GTSH_2012)
! ---------------------------------------------------------------->>>
          DO M = 1, SMAX

            DO IJK = ijkstart3, ijkend3

               IF(.NOT.ADDED_MASS .OR. M /= M_AM) THEN
                  CpxFlux_E(IJK) = 1.5D0 * Flux_sE(IJK,M)
                  CpxFlux_N(IJK) = 1.5D0 * Flux_sN(IJK,M)
                  CpxFlux_T(IJK) = 1.5D0 * Flux_sT(IJK,M)
               ELSE
                  CpxFlux_E(IJK) = 1.5D0 * Flux_sSE(IJK)
                  CpxFlux_N(IJK) = 1.5D0 * Flux_sSN(IJK)
                  CpxFlux_T(IJK) = 1.5D0 * Flux_sST(IJK)
               ENDIF

               IF (FLUID_AT(IJK)) THEN
! calculate the source terms to be used in the a matrix and b vector
                  IF (KT_TYPE_ENUM == GD_1999 .OR. &
                      KT_TYPE_ENUM == GTSH_2012) THEN
                     CALL SOURCE_GRANULAR_ENERGY_GD(SOURCELHS, &
                        SOURCERHS, IJK, M)
                  ELSE
                     CALL SOURCE_GRANULAR_ENERGY(SOURCELHS, &
                        SOURCERHS, IJK, M)
                  ENDIF
                  APO = 1.5D0*ROP_SO(IJK,M)*VOL(IJK)*ODT
                  S_P(IJK) = APO + SOURCELHS + 1.5D0* &
                     ZMAX(SUM_R_S(IJK,M)) * VOL(IJK)
                  S_C(IJK) = APO*THETA_MO(IJK,M) + SOURCERHS + 1.5d0* &
                      THETA_M(IJK,M)*ZMAX((-SUM_R_S(IJK,M))) * VOL(IJK)
                  EPS(IJK) = EP_S(IJK,M)
! MMS Source term.
                  IF(USE_MMS) S_C(IJK) = S_C(IJK) + &
                     MMS_THETA_M_SRC(IJK)*VOL(IJK)
               ELSE
                  EPS(IJK) = ZERO
                  S_P(IJK) = ZERO
                  S_C(IJK) = ZERO
                  IF(USE_MMS) EPS(IJK) = EP_S(IJK,M)
               ENDIF
            ENDDO   ! end do loop (ijk=ijkstart3,ijkend3)

! calculate the convection-diffusion terms
            CALL CONV_DIF_PHI (THETA_M(1,M), KTH_S(1,M), DISCRETIZE(8), &
               U_S(1,M), V_S(1,M), W_S(1,M), &
               CpxFlux_E, CpxFlux_N, CpxFlux_T, M, A_M, B_M)

! calculate standard bc
            CALL BC_PHI (THETA_M(1,M), BC_THETA_M(1,M), BC_THETAW_M(1,M), &
               BC_HW_THETA_M(1,M), BC_C_THETA_M(1,M), M, A_M, B_M)

! override bc settings if Johnson-Jackson bcs are specified
            CALL BC_THETA (M, A_M, B_M)

! set the source terms in a and b matrix form
            CALL SOURCE_PHI (S_P, S_C, EPS, THETA_M(1,M), M, A_M, B_M)

! Adjusting the values of theta_m to zero when Ep_g < EP_star
! (Shaeffer, 1987). This is done here instead of calc_mu_s to
! avoid convergence problems. (sof)
            IF (SCHAEFFER) THEN
               DO IJK = ijkstart3, ijkend3
                  IF (FLUID_AT(IJK) .AND. &
                      EP_g(IJK) .LT. EP_g_blend_start(ijk)) THEN

                     A_M(IJK,1,M) = ZERO
                     A_M(IJK,-1,M) = ZERO
                     A_M(IJK,2,M) = ZERO
                     A_M(IJK,-2,M) = ZERO
                     A_M(IJK,3,M) = ZERO
                     A_M(IJK,-3,M) = ZERO
                     A_M(IJK,0,M) = -ONE
                     B_M(IJK,M) = -smallTheta
                  ENDIF
               ENDDO
            ENDIF

            CALL CALC_RESID_S (THETA_M(1,M), A_M, B_M, M, &
               NUM_RESID(RESID_TH,M), &
               DEN_RESID(RESID_TH,M), RESID(RESID_TH,M), &
               MAX_RESID(RESID_TH,M), IJK_RESID(RESID_TH,M), ZERO)

            CALL UNDER_RELAX_S (THETA_M(1,M), A_M, B_M, M, UR_FAC(8))

!            call check_ab_m(a_m, b_m, m, .true., ier)
!            write(*,*) resid(resid_th, m), max_resid(resid_th, m),&
!               ijk_resid(resid_th, m)
!            call write_ab_m(a_m, b_m, ijkmax2, m, ier)

!            call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), &
!               1, DO_K, ier)

            CALL ADJUST_LEQ (RESID(RESID_TH,M), LEQ_IT(8), LEQ_METHOD(8),&
               LEQI, LEQM)

            CALL SOLVE_LIN_EQ ('Theta_m', 8, THETA_M(1,M), A_M, B_M, M,&
               LEQI, LEQM, LEQ_SWEEP(8), LEQ_TOL(8),  LEQ_PC(8), IER)
!            call out_array(Theta_m(1,m), 'Theta_m')

! Check for linear solver divergence.
            IF(ier == -2) Err_l(myPE) = 140

! Remove very small negative values of theta caused by leq solvers
            CALL ADJUST_THETA (M, IER)
! large negative granular temp -> divergence
            IF (IER /= 0) Err_l(myPE) = 141

          ENDDO   ! end do loop (m=1,smax)
! end case(lun_1984, simonin_1996, ahmadi_1995, gd_1999, gtsh_2012)
! ----------------------------------------------------------------<<<

        CASE(IA_2005)
! ---------------------------------------------------------------->>>
          DO M = 1, SMAX
            DO IJK = ijkstart3, ijkend3

! Skip walls where some values are undefined.
               IF(WALL_AT(IJK)) cycle

               D_PM = D_P(IJK,M)
               M_PM = (PI/6.d0)*(D_PM**3)*RO_S(IJK,M)

! In Iddir & Arastoopour (2005) the granular temperature includes
! mass of the particle in the definition.
               IF(.NOT.ADDED_MASS .OR. M /= M_AM) THEN
                  CpxFlux_E(IJK) = (1.5D0/M_PM) * Flux_sE(IJK,M)
                  CpxFlux_N(IJK) = (1.5D0/M_PM) * Flux_sN(IJK,M)
                  CpxFlux_T(IJK) = (1.5D0/M_PM) * Flux_sT(IJK,M)
               ELSE ! in case added mass is used.
                  CpxFlux_E(IJK) = (1.5D0/M_PM) * Flux_sSE(IJK)
                  CpxFlux_N(IJK) = (1.5D0/M_PM) * Flux_sSN(IJK)
                  CpxFlux_T(IJK) = (1.5D0/M_PM) * Flux_sST(IJK)
               ENDIF

               IF (FLUID_AT(IJK)) THEN
! calculate the source terms to be used in the a matrix and b vector
                  CALL SOURCE_GRANULAR_ENERGY_IA(SOURCELHS, &
                     SOURCERHS, IJK, M)
                  APO = (1.5D0/M_PM)*ROP_SO(IJK,M)*VOL(IJK)*ODT
                  S_P(IJK) = APO + SOURCELHS + (1.5d0/M_PM)*&
                     ZMAX(SUM_R_S(IJK,M)) * VOL(IJK)
                  S_C(IJK) = APO*THETA_MO(IJK,M) + SOURCERHS + &
                     (1.50d0/M_PM)*THETA_M(IJK,M)*&
                     ZMAX((-SUM_R_S(IJK,M))) * VOL(IJK)
                  EPS(IJK) = EP_S(IJK,M)
               ELSE
                  EPS(IJK) = ZERO
                  S_P(IJK) = ZERO
                  S_C(IJK) = ZERO
               ENDIF
            ENDDO    ! end do loop (ijk=ijkstart3,ijkend3)

! calculate the convection-diffusion terms
            CALL CONV_DIF_PHI (THETA_M(1,M), KTH_S(1,M), DISCRETIZE(8),&
               U_S(1,M), V_S(1,M), W_S(1,M), &
               CpxFlux_E, CpxFlux_N, CpxFlux_T, M, A_M, B_M)

! calculate standard bc
            CALL BC_PHI (THETA_M(1,M), BC_THETA_M(1,M), BC_THETAW_M(1,M),&
               BC_HW_THETA_M(1,M),BC_C_THETA_M(1,M), M, A_M, B_M)

! override bc settings if Johnson-Jackson bcs are specified
            CALL BC_THETA (M, A_M, B_M)

! set the source terms in a and b matrix form
            CALL SOURCE_PHI (S_P, S_C, EPS, THETA_M(1,M), M, A_M, B_M)
          ENDDO   ! end do loop (m = 1, smax)

! use partial elimination on collisional dissipation term: SUM(Nip)
!     SUM( ED_s_ip* (Theta_p-Theta_i))
          IF (SMAX > 1) THEN
            CALL CALC_VTC_SS (VXTC_SS)
            CALL PARTIAL_ELIM_IA (THETA_M, VXTC_SS, A_M, B_M)
          ENDIF

! Adjusting the values of theta_m to zero when Ep_g < EP_star
! (Shaeffer, 1987). This is done here instead of calc_mu_s to
! avoid convergence problems. (sof)
          IF (SCHAEFFER) THEN
            DO M = 1, SMAX
               DO IJK = ijkstart3, ijkend3
                  IF (FLUID_AT(IJK) .AND. &
                      EP_g(IJK) .LT. EP_g_blend_start(ijk)) THEN

                     D_PM = D_P(IJK,M)
                     M_PM = (PI/6.d0)*(D_PM**3)*RO_S(IJK,M)
                     A_M(IJK,1,M) = ZERO
                     A_M(IJK,-1,M) = ZERO
                     A_M(IJK,2,M) = ZERO
                     A_M(IJK,-2,M) = ZERO
                     A_M(IJK,3,M) = ZERO
                     A_M(IJK,-3,M) = ZERO
                     A_M(IJK,0,M) = -ONE
! In Iddir & Arastoopour (2005) the granular temperature includes
! mass of the particle in the definition. Systems with small
! particles can give rise to smaller temperatures than the standard
! zero_ep_s
                     B_M(IJK,M) = -smallTheta*M_PM
                  ENDIF
               ENDDO
            ENDDO ! for M
          ENDIF  ! end if (schaeffer)

          DO M = 1, SMAX
            CALL CALC_RESID_S (THETA_M(1,M), A_M, B_M, M, &
               NUM_RESID(RESID_TH,M), DEN_RESID(RESID_TH,M), RESID(RESID_TH,M),&
               MAX_RESID(RESID_TH,M), IJK_RESID(RESID_TH,M), ZERO)

            CALL UNDER_RELAX_S (THETA_M(1,M), A_M, B_M, M, UR_FAC(8))

            CALL ADJUST_LEQ (RESID(RESID_TH,M), LEQ_IT(8), &
               LEQ_METHOD(8), LEQI, LEQM)

            CALL SOLVE_LIN_EQ ('Theta_m', 8, THETA_M(1,M), A_M, B_M, M,&
               LEQI, LEQM, LEQ_SWEEP(8), LEQ_TOL(8),  LEQ_PC(8), IER)
!            call out_array(Theta_m(1,m), 'Theta_m')

! Check for linear solver divergence.
            IF(ier == -2) Err_l(myPE) = 140

! Remove very small negative values of theta caused by leq solvers
            CALL ADJUST_THETA (M, IER)
! large negative granular temp -> divergence
            IF (IER /= 0) Err_l(myPE) = 141

          ENDDO   ! end do loop (m=1,smax)
! end case(ia_2005)
! ----------------------------------------------------------------<<<


        CASE (GHD_2007)
! ---------------------------------------------------------------->>>
! do not loop over solids phases. solve for the mixture phase
          M = MMAX

! initialize
          TOT_NO(:) = ZERO
          TOT_SUM_RS(:) = ZERO
          TOT_EPS(:) = ZERO

          CALL CALC_NFLUX ()

          DO IJK = ijkstart3, ijkend3

! total number density is used for GHD theory
            CpxFlux_E(IJK) = 1.5D0 * Flux_nE(IJK)
            CpxFlux_N(IJK) = 1.5D0 * Flux_nN(IJK)
            CpxFlux_T(IJK) = 1.5D0 * Flux_nT(IJK)

            IF (FLUID_AT(IJK)) THEN
               DO L = 1,SMAX
                  M_PM = (PI/6.d0)*(D_P(IJK,L)**3)*RO_S(IJK,L)
                  TOT_SUM_RS(IJK) = TOT_SUM_RS(IJK) + SUM_R_S(IJK,L)/M_PM
                  TOT_NO(IJK) = TOT_NO(IJK) + ROP_SO(IJK,L)/M_PM
                  TOT_EPS(IJK) = TOT_EPS(IJK) + EP_S(IJK,L)
               ENDDO

! calculate the source terms to be used in the a matrix and b vector
               CALL SOURCE_GHD_GRANULAR_ENERGY (SOURCELHS, &
                  SOURCERHS, IJK)
               APO = 1.5D0 * TOT_NO(IJK)*VOL(IJK)*ODT
               S_P(IJK) = APO + SOURCELHS + 1.5d0 *&
                  ZMAX(TOT_SUM_RS(IJK)) * VOL(IJK)
               S_C(IJK) = APO*THETA_MO(IJK,M) + SOURCERHS + 1.50d0 *&
                  THETA_M(IJK,M)*ZMAX(-TOT_SUM_RS(IJK)) * VOL(IJK)
               EPS(IJK) = TOT_EPS(IJK)
            ELSE
               EPS(IJK) = ZERO
               S_P(IJK) = ZERO
               S_C(IJK) = ZERO
            ENDIF
          ENDDO   ! end do loop (ijk=ijkstart3,ijkend3)

! calculate the convection-diffusion terms
          CALL CONV_DIF_PHI (THETA_M(1,M), KTH_S(1,M), DISCRETIZE(8),&
            U_S(1,M), V_S(1,M), W_S(1,M), &
            CpxFlux_E, CpxFlux_N, CpxFlux_T, M, A_M, B_M)

! calculate standard bc
          CALL BC_PHI (THETA_M(1,M), BC_THETA_M(1,M), BC_THETAW_M(1,M), &
            BC_HW_THETA_M(1,M), BC_C_THETA_M(1,M), M, A_M, B_M)

! override bc settings if Johnson-Jackson bcs are specified
          CALL BC_THETA (M, A_M, B_M)

! set the source terms in a and b matrix form
          CALL SOURCE_PHI (S_P, S_C, EPS, THETA_M(1,M), M, A_M, B_M)

          CALL CALC_RESID_S (THETA_M(1,M), A_M, B_M, M, &
            NUM_RESID(RESID_TH,M), DEN_RESID(RESID_TH,M), RESID(RESID_TH,M),&
            MAX_RESID(RESID_TH,M), IJK_RESID(RESID_TH,M), ZERO)

          CALL UNDER_RELAX_S (THETA_M(1,M), A_M, B_M, M, UR_FAC(8))

          CALL ADJUST_LEQ (RESID(RESID_TH,M), LEQ_IT(8), LEQ_METHOD(8),&
            LEQI, LEQM)

          CALL SOLVE_LIN_EQ ('Theta_m', 8, THETA_M(1,M), A_M, B_M, M, &
            LEQI, LEQM, LEQ_SWEEP(8), LEQ_TOL(8),  LEQ_PC(8), IER)
!         call out_array(Theta_m(1,m), 'Theta_m')

! Check for linear solver divergence.
          IF(ier == -2) Err_l(myPE) = 140

! Remove very small negative values of theta caused by leq solvers
          CALL ADJUST_THETA (M, IER)

! large negative granular temp -> divergence
          IF (IER /= 0) Err_l(myPE) = 141

! end case(ghd_2007)
! ----------------------------------------------------------------<<<

        CASE DEFAULT
! should never hit this
          WRITE (*, '(A)') 'ADJUST_THETA'
          WRITE (*, '(A,A)') 'Unknown KT_TYPE: ', KT_TYPE
          call mfix_exit(myPE)
      END SELECT   ! end selection of kt_type_enum

      call unlock_ambm
      call unlock_tmp_array


      CALL global_all_sum(Err_l, Err_g)
      IER = maxval(Err_g)


      RETURN
      END SUBROUTINE SOLVE_GRANULAR_ENERGY


