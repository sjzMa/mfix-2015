!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_VEL_STAR                                          C
!  Purpose: Solve starred velocity components                          C
!                                                                      C
!  Author: M. Syamlal                                 Date: 25-APR-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: allow multiparticle in D_E, D_N and D_T claculations       C
!           and account for the averaged Solid-Solid drag              C
!  Author: S. Dartevelle, LANL                        Date: 28-FEb-04  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 3                                                  C
!  Purpose: To flag the solids calculations: Continuum or DES          C
!           and to change the gas arrays incorporating the drag        C
!           when doing DES                                             C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!                                                                      C
!  Revision Number: 4                                                  C
!  Purpose: Incorporation of QMOM for the solution of the particle     C
!           kinetic equation                                           C
!  Author: Alberto Passalacqua - Fox Research Group   Date: 02-Dec-09  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOLVE_VEL_STAR(IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE toleranc
      USE run
      USE physprop
      USE geometry
      USE fldvar
      USE ghdtheory
      USE output
      USE indices
      USE drag
      USE residual
      USE ur_facs
      USE pgcor
      USE pscor
      USE leqsol
      Use ambm
      Use tmp_array1,  VxF_gs => Arraym1
      Use tmp_array,  VxF_ss => ArrayLM
      USE compar
      USE discretelement
      USE qmom_kinetic_equation
      use ps

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! solids phase index
      INTEGER :: M
! fluid cell index
      INTEGER :: IJK
! temporary velocity arrays
      DOUBLE PRECISION, DIMENSION(:), allocatable :: U_gtmp,  V_gtmp, W_gtmp
      DOUBLE PRECISION, DIMENSION(:,:), allocatable :: U_stmp, V_stmp, W_stmp
! linear equation solver method and iterations
      INTEGER :: LEQM, LEQI

      LOGICAL :: DO_SOLIDS

! temporary use of global arrays:
! arraym1 (locally vxf_gs)
! the volume x average gas-solids drag at momentum cell centers
!      DOUBLE PRECISION :: VXF_GS(DIMENSION_3, DIMENSION_M)
! arraylm (locally vxf_ss)
! the volume x average solids-solids drag at momentum cell centers
!      DOUBLE PRECISION :: VXF_SS(DIMENSION_3, DIMENSION_LM)
! Septadiagonal matrix A_m, vector b_m
!      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)

!-----------------------------------------------

      allocate(U_gtmp(DIMENSION_3))
      allocate(V_gtmp(DIMENSION_3))
      allocate(W_gtmp(DIMENSION_3))
      allocate(U_stmp(DIMENSION_3,DIMENSION_M))
      allocate(V_stmp(DIMENSION_3,DIMENSION_M))
      allocate(W_stmp(DIMENSION_3,DIMENSION_M))

      call lock_ambm        ! locks arrys a_m and b_m
      call lock_tmp_array1  ! locks arraym1 (locally vxf_gs)
      call lock_tmp_array2  ! locks arraylm (locally vxf_ss)

! Store the velocities so that the order of solving the momentum
! equations does not matter
      DO IJK = ijkstart3, ijkend3
         U_gtmp(IJK) = U_g(IJK)
         V_gtmp(IJK) = V_g(IJK)
         W_gtmp(IJK) = W_g(IJK)
      ENDDO
      DO M = 1, MMAX
         IF(TRIM(KT_TYPE) /= 'GHD' .OR. &
            (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN
            DO IJK = ijkstart3, ijkend3
               U_stmp(IJK, M) = U_s(IJK, M)
               V_stmp(IJK, M) = V_s(IJK, M)
               W_stmp(IJK, M) = W_s(IJK, M)
            ENDDO
         ENDIF
      ENDDO

      DO_SOLIDS = .NOT.(DISCRETE_ELEMENT .OR. QMOMK) .OR. &
         DES_CONTINUUM_HYBRID


! Calculate U_m_star and residuals
! ---------------------------------------------------------------->>>
      DO M = 0, MMAX
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M)
         IF (M >= 1) VXF_GS(:,M) = ZERO
      ENDDO

! calculate the convection-diffusion terms for the gas and solids phase
! u-momentum equations
      CALL CONV_DIF_U_G (A_M, B_M)
      IF(DO_SOLIDS) CALL CONV_DIF_U_S (A_M, B_M)

! calculate the source terms for the gas and solids phase u-momentum
! equations
      CALL SOURCE_U_G (A_M, B_M)
      IF(POINT_SOURCE) CALL POINT_SOURCE_U_G (A_M, B_M)
      IF(DO_SOLIDS) THEN
         CALL SOURCE_U_S (A_M, B_M)
         IF(POINT_SOURCE) CALL POINT_SOURCE_U_S (A_M, B_M)
      ENDIF

! evaluate local variable vxf_gs and vxf_ss.  both terms are sent to the
! subroutine calc_d (pressure correction equation coefficients).  the
! former is also used in the subroutine partial_elim_u while the latter
! is effectively re-evaluated within said subroutine
      CALL VF_GS_X (VXF_GS)
      IF(DO_SOLIDS .AND. (TRIM(KT_TYPE) /= 'GHD')) THEN
         IF (MMAX > 0) CALL VF_SS_X (VXF_SS)
      ENDIF

! calculate coefficients for the pressure correction equation
      IF(TRIM(KT_TYPE) == 'GHD') THEN
         CALL CALC_D_GHD_E (A_M, VXF_GS, D_E)
      ELSE
         CALL CALC_D_E (A_M, VXF_GS, VXF_SS, D_E, IER)
      ENDIF

      IF(DO_SOLIDS) THEN
! calculate coefficients for a solids volume correction equation
         IF (MMAX > 0) CALL CALC_E_E (A_M, MCP, E_E)

! calculate modifications to the A matrix center coefficient and B
! source vector for partial elimination
         IF(TRIM(KT_TYPE) == 'GHD') THEN
            IF (MMAX > 0) CALL PARTIAL_ELIM_GHD_U (U_G,U_S,VXF_GS,A_M,B_M)
         ELSE
            IF (MMAX > 0) CALL PARTIAL_ELIM_U (U_G,U_S,VXF_GS,A_M,B_M)
         ENDIF
      ENDIF

! handle special case where center coefficient is zero
      CALL ADJUST_A_U_G (A_M, B_M)
      IF(DO_SOLIDS) CALL ADJUST_A_U_S (A_M, B_M)

! calculate modifications to the A matrix center coefficient and B
! source vector for treating DEM drag terms
      IF(DES_CONTINUUM_COUPLED) THEN
         CALL GAS_DRAG_U(A_M, B_M, IER)
         IF (DES_CONTINUUM_HYBRID) &
            CALL SOLID_DRAG_U(A_M, B_M)
      ENDIF

      IF(QMOMK .AND. QMOMK_COUPLED) THEN
         CALL QMOMK_GAS_DRAG(A_M, B_M, IER, 1, 0, 0)
      ENDIF

      IF (MOMENTUM_X_EQ(0)) THEN
         CALL CALC_RESID_U (U_G, V_G, W_G, A_M, B_M, 0, &
            NUM_RESID(RESID_U,0), DEN_RESID(RESID_U,0), &
            RESID(RESID_U,0), MAX_RESID(RESID_U,0), &
            IJK_RESID(RESID_U,0))
         CALL UNDER_RELAX_U (U_G, A_M, B_M, 0, UR_FAC(3))
!         call check_ab_m(a_m, b_m, 0, .false., ier)
!         call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!         write(*,*) &
!            resid(resid_u, 0), max_resid(resid_u, 0), &
!            ijk_resid(resid_u, 0)
      ENDIF

      IF(DO_SOLIDS) THEN
         DO M = 1, MMAX
            IF(TRIM(KT_TYPE) /= 'GHD' .OR. &
              (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN
               IF (MOMENTUM_X_EQ(M)) THEN
                  CALL CALC_RESID_U (U_S(1,M), V_S(1,M), W_S(1,M), A_M,&
                     B_M, M, NUM_RESID(RESID_U,M), &
                     DEN_RESID(RESID_U,M), RESID(RESID_U,M), &
                     MAX_RESID(RESID_U,M), IJK_RESID(RESID_U,M))
                  CALL UNDER_RELAX_U (U_S(1,M), A_M, B_M, M, &
                     UR_FAC(3))
!                  call check_ab_m(a_m, b_m, m, .false., ier)
!                  write(*,*) &
!                     resid(resid_u, m), max_resid(resid_u, m), &
!                     ijk_resid(resid_u, m)
!                  call write_ab_m(a_m, b_m, ijkmax2, m, ier)
               ENDIF   ! end if (momentum_x_eq(m))
            ENDIF ! end if check for GHD Theory
         ENDDO   ! end do (m=1,mmax)
      ENDIF

      IF (MOMENTUM_X_EQ(0)) THEN
!         call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1,&
!            DO_K,ier)
         CALL ADJUST_LEQ (RESID(RESID_U,0), LEQ_IT(3), LEQ_METHOD(3), &
            LEQI, LEQM)
         CALL SOLVE_LIN_EQ ('U_g', 3, U_Gtmp, A_M, B_M, 0, LEQI, LEQM, &
            LEQ_SWEEP(3), LEQ_TOL(3),  LEQ_PC(3), IER)
!         call out_array(u_g, 'u_g')
      ENDIF

      IF(DO_SOLIDS) THEN
         DO M = 1, MMAX
            IF(TRIM(KT_TYPE) /= 'GHD' .OR. &
              (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN
               IF (MOMENTUM_X_EQ(M)) THEN
!                  call test_lin_eq(ijkmax2, ijmax2, imax2, &
!                     a_m(1, -3, M), 1, DO_K,ier)
                  CALL ADJUST_LEQ (RESID(RESID_U,M), LEQ_IT(3),&
                     LEQ_METHOD(3), LEQI, LEQM)
                  CALL SOLVE_LIN_EQ ('U_s', 3, U_Stmp(1,M), A_M, &
                     B_M, M, LEQI, LEQM, LEQ_SWEEP(3), LEQ_TOL(3),&
                     LEQ_PC(3), IER)
!                  call out_array(u_s(1,m), 'u_s')
               ENDIF   ! end if (momentum_x_eq(m))
            ENDIF ! end if check for GHD Theory
         ENDDO   ! end do (m=1,mmax)
      ENDIF
! End U_m_star and residuals
! ----------------------------------------------------------------<<<

! Calculate V_m_star and residuals
! ---------------------------------------------------------------->>>
      DO M = 0, MMAX
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M)
         IF (M >= 1) VXF_GS(:,M) = ZERO
      ENDDO

      CALL CONV_DIF_V_G (A_M, B_M, IER)
!      call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
      IF(DO_SOLIDS) CALL CONV_DIF_V_S (A_M, B_M, IER)

      CALL SOURCE_V_G (A_M, B_M)
      IF(POINT_SOURCE) CALL POINT_SOURCE_V_G (A_M, B_M)

!      call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
      IF(DO_SOLIDS) THEN
         CALL SOURCE_V_S (A_M, B_M)
         IF(POINT_SOURCE) CALL POINT_SOURCE_V_S (A_M, B_M)
      END IF

      CALL VF_GS_Y (VXF_GS)
      IF(DO_SOLIDS .AND. (TRIM(KT_TYPE) /= 'GHD')) THEN
         IF (MMAX > 0) CALL VF_SS_Y (VXF_SS)
      ENDIF

! calculate coefficients for the pressure correction equation
      IF(TRIM(KT_TYPE) == 'GHD') THEN
         CALL CALC_D_GHD_N (A_M, VXF_GS, D_N)
      ELSE
         CALL CALC_D_N (A_M, VXF_GS, VXF_SS, D_N, IER)
      ENDIF

      IF(DO_SOLIDS) THEN
! calculate coefficients for a solids pressure correction equation
         IF (MMAX > 0) CALL CALC_E_N (A_M, MCP, E_N)

! calculate modifications to the A matrix center coefficient and B
! source vector for partial elimination
         IF(TRIM(KT_TYPE) == 'GHD') THEN
           IF (MMAX > 0) CALL PARTIAL_ELIM_GHD_V (V_G,V_S,VXF_GS,A_M,B_M)
         ELSE
           IF (MMAX > 0) CALL PARTIAL_ELIM_V (V_G,V_S,VXF_GS,A_M,B_M)
         ENDIF
      ENDIF

!      call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
      CALL ADJUST_A_V_G (A_M, B_M)
!      call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
      IF(.NOT.(DISCRETE_ELEMENT .OR. QMOMK) .OR. &
         DES_CONTINUUM_HYBRID) THEN
         CALL ADJUST_A_V_S (A_M, B_M)
      ENDIF

      IF(DES_CONTINUUM_COUPLED) THEN
         CALL GAS_DRAG_V(A_M, B_M, IER)
         IF (DES_CONTINUUM_HYBRID) &
            CALL SOLID_DRAG_V(A_M, B_M)
      ENDIF

      IF(QMOMK .AND. QMOMK_COUPLED) THEN
         CALL QMOMK_GAS_DRAG(A_M, B_M, IER, 0, 1, 0)
      ENDIF


      IF (MOMENTUM_Y_EQ(0)) THEN
         CALL CALC_RESID_V (U_G, V_G, W_G, A_M, B_M, 0, &
            NUM_RESID(RESID_V,0), DEN_RESID(RESID_V,0), &
            RESID(RESID_V,0), MAX_RESID(RESID_V,0), &
            IJK_RESID(RESID_V,0))
         CALL UNDER_RELAX_V (V_G, A_M, B_M, 0, UR_FAC(4))
!         call check_ab_m(a_m, b_m, 0, .false., ier)
!         call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!         write(*,*) &
!            resid(resid_v, 0), max_resid(resid_v, 0), &
!            ijk_resid(resid_v, 0)
      ENDIF

      IF(DO_SOLIDS) THEN
         DO M = 1, MMAX
            IF(TRIM(KT_TYPE) /= 'GHD' .OR. &
              (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN
               IF (MOMENTUM_Y_EQ(M)) THEN
                  CALL CALC_RESID_V (U_S(1,M), V_S(1,M), W_S(1,M), A_M,&
                     B_M, M, NUM_RESID(RESID_V,M), &
                     DEN_RESID(RESID_V,M),RESID(RESID_V,M), &
                     MAX_RESID(RESID_V,M), IJK_RESID(RESID_V,M))
                  CALL UNDER_RELAX_V (V_S(1,M),A_M,B_M,M,UR_FAC(4))
!                  call check_ab_m(a_m, b_m, m, .false., ier)
!                  write(*,*) &
!                     resid(resid_v, m), max_resid(resid_v, m),
!                     ijk_resid(resid_v, m)
!                  call write_ab_m(a_m, b_m, ijkmax2, m, ier)
               ENDIF   ! end if (momentum_y_eq(m))
            ENDIF ! end if check for GHD Theory
         ENDDO   ! end do (m=1,mmax)
      ENDIF

      IF (MOMENTUM_Y_EQ(0)) THEN
!         call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), &
!            1, DO_K, ier)
         CALL ADJUST_LEQ (RESID(RESID_V,0), LEQ_IT(4), LEQ_METHOD(4), &
            LEQI, LEQM)
         CALL SOLVE_LIN_EQ ('V_g', 4, V_Gtmp, A_M, B_M, 0, LEQI, LEQM, &
            LEQ_SWEEP(4), LEQ_TOL(4),  LEQ_PC(4), IER)
!         call out_array(v_g, 'v_g')
      ENDIF

      IF(DO_SOLIDS) THEN
         DO M = 1, MMAX
            IF(TRIM(KT_TYPE) /= 'GHD' .OR. &
              (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN
               IF (MOMENTUM_Y_EQ(M)) THEN
!                  call test_lin_eq(ijkmax2, ijmax2, imax2, &
!                     a_m(1, -3, M), 1, DO_K, ier)
                  CALL ADJUST_LEQ (RESID(RESID_V,M), LEQ_IT(4), &
                     LEQ_METHOD(4), LEQI, LEQM)
                  CALL SOLVE_LIN_EQ ('V_s', 4, V_Stmp(1,M), A_M, &
                     B_M, M, LEQI, LEQM, LEQ_SWEEP(4), LEQ_TOL(4), &
                     LEQ_PC(4), IER)
!                  call out_array(v_s(1,m), 'v_s')
               ENDIF   ! end if (momentum_y_eq(m))
            ENDIF ! end if check for GHD Theory
         ENDDO   ! end do (m=1,mmax)
      ENDIF
! End V_m_star and residuals
! ----------------------------------------------------------------<<<


! Calculate W_m_star and residuals
! ---------------------------------------------------------------->>>
      IF (DO_K)THEN
         DO M = 0, MMAX
            CALL INIT_AB_M (A_M, B_M, IJKMAX2, M)
            IF (M >= 1) VXF_GS(:,M) = ZERO
         ENDDO

         CALL CONV_DIF_W_G (A_M, B_M)
         IF(DO_SOLIDS) CALL CONV_DIF_W_S (A_M, B_M)

         CALL SOURCE_W_G (A_M, B_M)
         IF(POINT_SOURCE) CALL POINT_SOURCE_W_G (A_M, B_M)
!         call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
         IF(DO_SOLIDS) THEN
            CALL SOURCE_W_S (A_M, B_M)
            IF(POINT_SOURCE) CALL POINT_SOURCE_W_S (A_M, B_M)
         ENDIF
!        call write_ab_m(a_m, b_m, ijkmax2, 0, ier)

         CALL VF_GS_Z (VXF_GS)
         IF(DO_SOLIDS .AND. (TRIM(KT_TYPE) /= 'GHD')) THEN
            IF (MMAX > 0) CALL VF_SS_Z (VXF_SS)
         ENDIF

! calculate coefficients for the pressure correction equation
         IF(TRIM(KT_TYPE) == 'GHD') THEN
            CALL CALC_D_GHD_T (A_M, VXF_GS, D_T)
         ELSE
            CALL CALC_D_T (A_M, VXF_GS, VXF_SS, D_T, IER)
         ENDIF

         IF(DO_SOLIDS) THEN
! calculate coefficients for a solids pressure correction equation
            IF (MMAX > 0) CALL CALC_E_T (A_M, MCP, E_T)

! calculate modifications to the A matrix center coefficient and B
! source vector for partial elimination
            IF(TRIM(KT_TYPE) == 'GHD') THEN
               IF (MMAX > 0) CALL PARTIAL_ELIM_GHD_W (W_G, W_S, VXF_GS, A_M, B_M)
            ELSE
               IF (MMAX > 0) CALL PARTIAL_ELIM_W (W_G, W_S, VXF_GS, A_M, B_M)
            ENDIF
         ENDIF

         CALL ADJUST_A_W_G (A_M, B_M)
         IF(.NOT.(DISCRETE_ELEMENT .OR. QMOMK) .OR. &
            DES_CONTINUUM_HYBRID) THEN
            CALL ADJUST_A_W_S (A_M, B_M)
         ENDIF

         IF(DES_CONTINUUM_COUPLED) THEN
            CALL GAS_DRAG_W(A_M, B_M, IER)
            IF (DISCRETE_ELEMENT .AND. DES_CONTINUUM_HYBRID) &
               CALL SOLID_DRAG_W(A_M, B_M)
         ENDIF

         IF(QMOMK .AND. QMOMK_COUPLED) THEN
            CALL QMOMK_GAS_DRAG(A_M, B_M, IER, 0, 0, 1)
         ENDIF

         IF (MOMENTUM_Z_EQ(0)) THEN
            CALL CALC_RESID_W (U_G, V_G, W_G, A_M, B_M, 0, &
               NUM_RESID(RESID_W,0), DEN_RESID(RESID_W,0), &
               RESID(RESID_W,0), MAX_RESID(RESID_W,0), &
               IJK_RESID(RESID_W,0))
            CALL UNDER_RELAX_W (W_G, A_M, B_M, 0, UR_FAC(5))
!            call check_ab_m(a_m, b_m, 0, .false., ier)
!            write(*,*) &
!               resid(resid_w, 0), max_resid(resid_w, 0), &
!               ijk_resid(resid_w, 0)
         ENDIF

         IF(DO_SOLIDS) THEN
            DO M = 1, MMAX
               IF(TRIM(KT_TYPE) /= 'GHD' .OR. &
                 (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN
                  IF (MOMENTUM_Z_EQ(M)) THEN
                     CALL CALC_RESID_W (U_S(1,M), V_S(1,M), W_S(1,M),&
                        A_M, B_M, M, NUM_RESID(RESID_W,M), &
                        DEN_RESID(RESID_W,M), RESID(RESID_W,M), &
                        MAX_RESID(RESID_W,M), IJK_RESID(RESID_W,M))
                     CALL UNDER_RELAX_W (W_S(1,M), A_M, B_M, M, &
                        UR_FAC(5))
!                     call check_ab_m(a_m, b_m, m, .false., ier)
!                     write(*,*) &
!                        resid(resid_w, m), max_resid(resid_w, m), &
!                        ijk_resid(resid_w, m)
!                     call write_ab_m(a_m, b_m, ijkmax2, m, ier)
                  ENDIF   ! end if (momentum_z_eq(m))
               ENDIF ! end if check for GHD Theory
            ENDDO   ! end do (m=1,mmax)
         ENDIF

         IF (MOMENTUM_Z_EQ(0)) THEN
!            call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), &
!               1, DO_K, ier)
            CALL ADJUST_LEQ (RESID(RESID_W,0), LEQ_IT(5), &
               LEQ_METHOD(5), LEQI, LEQM)
            CALL SOLVE_LIN_EQ ('W_g', 5, W_Gtmp, A_M, B_M, 0, LEQI,&
               LEQM, LEQ_SWEEP(5), LEQ_TOL(5), LEQ_PC(5), IER)
!            call out_array(w_g, 'w_g')
         ENDIF

         IF(DO_SOLIDS) THEN
            DO M = 1, MMAX
               IF(TRIM(KT_TYPE) /= 'GHD' .OR. &
                 (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN
                  IF (MOMENTUM_Z_EQ(M)) THEN
!                     call test_lin_eq(ijkmax2, ijmax2, imax2, &
!                        a_m(1, -3, M), 1, DO_K, ier)
                     CALL ADJUST_LEQ (RESID(RESID_W,M), LEQ_IT(5), &
                        LEQ_METHOD(5), LEQI, LEQM)
                     CALL SOLVE_LIN_EQ ('W_s', 5, W_Stmp(1,M), &
                        A_M, B_M, M, LEQI, LEQM, LEQ_SWEEP(5), &
                        LEQ_TOL(5), LEQ_PC(5), IER)
!                     call out_array(w_s(1,m), 'w_s')
                  ENDIF   ! end if (momentum_z_eq(m))
               ENDIF ! end if check for GHD Theory
            ENDDO   ! end do (m=1,mmax)
         ENDIF
      ENDIF   ! end if (do_k)
! End W_m_star and residuals
! ----------------------------------------------------------------<<<

! Now update all velocity components
      DO IJK = ijkstart3, ijkend3
         U_g(IJK) = U_gtmp(IJK)
         V_g(IJK) = V_gtmp(IJK)
         W_g(IJK) = W_gtmp(IJK)
      ENDDO
      DO M = 1, MMAX
         IF(TRIM(KT_TYPE) /= 'GHD' .OR. &
           (TRIM(KT_TYPE) == 'GHD' .AND. M==MMAX)) THEN
            DO IJK = ijkstart3, ijkend3
               U_s(IJK, M) = U_stmp(IJK, M)
               V_s(IJK, M) = V_stmp(IJK, M)
               W_s(IJK, M) = W_stmp(IJK, M)
            ENDDO
         ENDIF
      ENDDO

! modification for GHD theory to compute species velocity: Ui = Joi/(mi ni) + U.
      IF(TRIM(KT_TYPE) == 'GHD') THEN
         CALL calc_external_forces()
         CALL GHDMassFlux() ! to compute solid species mass flux
         CALL UpdateSpeciesVelocities() ! located at end of ghdMassFlux.f file
      ENDIF

      call unlock_ambm
      call unlock_tmp_array1
      call unlock_tmp_array2

      deallocate(U_gtmp)
      deallocate(V_gtmp)
      deallocate(W_gtmp)
      deallocate(U_stmp)
      deallocate(V_stmp)
      deallocate(W_stmp)

      RETURN
      END SUBROUTINE SOLVE_VEL_STAR
