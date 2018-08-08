
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_Pp_g
!  Purpose: Solve fluid pressure correction equation                   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOLVE_PP_G(NORMG, RESG, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE physprop
      USE geometry
      USE pgcor
      USE residual
      USE leqsol
      USE run
      Use ambm
      Use tmp_array1, B_mMAX => ARRAYm1
      use ps

      IMPLICIT NONE
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! Parameter to make tolerance for residual scaled with max value
! compatible with residual scaled with first iteration residual.
! Increase it to tighten convergence.
      DOUBLE PRECISION, PARAMETER :: DEN = 1.0D1   !5.0D2
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Normalization factor for gas pressure correction residual.
! At start of the iterate loop normg will either be 1 (i.e. not
! normalized) or a user defined value given by norm_g.  If norm_g
! was set to zero then the normalization is based on dominate
! term in the equation
      DOUBLE PRECISION, INTENT(IN) :: NORMg
! gas pressure correction residual
      DOUBLE PRECISION, INTENT(OUT) :: RESg
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! phase index
      INTEGER :: M
! Normalization factor for gas pressure correction residual
      DOUBLE PRECISION :: NORMGloc
! linear equation solver method and iterations
      INTEGER :: LEQM, LEQI

! temporary use of global arrays:
! arraym1 (locally b_mmax)
! vector B_M based on dominate term in correction equation
!      DOUBLE PRECISION :: B_MMAX(DIMENSION_3, DIMENSION_M)
! Septadiagonal matrix A_m, vector B_m
!      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------

      call lock_ambm
      call lock_tmp_array1

! initializing
      CALL ZERO_ARRAY (PP_G)
      DO M = 0, MMAX
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M)
      ENDDO

! If gas momentum equations in x and y directions are not solved return
      IF (.NOT.(MOMENTUM_X_EQ(0) .OR. MOMENTUM_Y_EQ(0)) .AND.&
          RO_G0 .NE. UNDEFINED) THEN
        call unlock_ambm
        call unlock_tmp_array1
        RETURN
      ENDIF

! Forming the sparse matrix equation.
      CALL CONV_PP_G (A_M, B_M)
      CALL SOURCE_PP_G (A_M, B_M, B_MMAX)
      IF(POINT_SOURCE) CALL POINT_SOURCE_PP_G (B_M, B_MMAX)

!      call check_ab_m(a_m, b_m, 0, .false., ier)
!      call write_ab_m(a_m, b_m, ijkmax2, 0, ier)


! Find average residual, maximum residual and location
      NORMGloc = NORMG
      IF(NORMG == ZERO) THEN
! calculating the residual based on dominate term in correction equation
! and use this to form normalization factor
        CALL CALC_RESID_PP (B_MMAX, ONE, NUM_RESID(RESID_P,0), &
         DEN_RESID(RESID_P,0), RESID(RESID_P,0), MAX_RESID(RESID_P,0), &
         IJK_RESID(RESID_P,0))
         NORMGloc = RESID(RESID_P,0)/DEN
      ENDIF
      CALL CALC_RESID_PP (B_M, NORMGloc, NUM_RESID(RESID_P,0),  &
         DEN_RESID(RESID_P,0), RESID(RESID_P,0), MAX_RESID(RESID_P,0), &
         IJK_RESID(RESID_P,0))
      RESG = RESID(RESID_P,0)
!      write(*,*) resid(resid_p, 0), max_resid(resid_p, 0), &
!         ijk_resid(resid_p, 0)


! Solve P_g_prime equation
       LEQI = LEQ_IT(1)
       LEQM = LEQ_METHOD(1)
!      CALL ADJUST_LEQ(RESID(RESID_P,0),LEQ_IT(1),LEQ_METHOD(1),LEQI,LEQM,IER)

!     call check_symmetry(A_m, 0, IER)
!     call test_lin_eq(A_M, LEQ_IT(1),LEQ_METHOD(1), LEQ_SWEEP(1), LEQ_TOL(1), LEQ_PC(1),0,IER)
      CALL SOLVE_LIN_EQ ('Pp_g', 1, PP_G, A_M, B_M, 0, LEQI, LEQM, &
                         LEQ_SWEEP(1), LEQ_TOL(1), LEQ_PC(1), IER)

!      call out_array(Pp_g, 'Pp_g')

      call unlock_tmp_array1
      call unlock_ambm

      RETURN
      END SUBROUTINE SOLVE_PP_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_Pp_g                                       C
!  Purpose: Adds point sources to the Pressure correction equation.    C
!                                                                      C
!  Notes: The off-diagonal coefficients are positive. The center       C
!         coefficient and the source vector are negative. See          C
!         conv_Pp_g                                                    C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE POINT_SOURCE_PP_G(B_M, B_mmax)

      use compar
      use constant
      use geometry
      use indices
      use param1, only: small_number
      use physprop
      use ps
      use run
      use functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! maximum term in b_m expression
      DOUBLE PRECISION, INTENT(INOUT) :: B_mmax(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------

! Indices
      INTEGER :: IJK, I, J, K
      INTEGER :: PSV

! terms of bm expression
      DOUBLE PRECISION pSource

!-----------------------------------------------
      PS_LP: do PSV = 1, DIMENSION_PS

         if(.NOT.PS_DEFINED(PSV)) cycle PS_LP
         if(PS_MASSFLOW_G(PSV) < small_number) cycle PS_LP

         do k = PS_K_B(PSV), PS_K_T(PSV)
         do j = PS_J_S(PSV), PS_J_N(PSV)
         do i = PS_I_W(PSV), PS_I_E(PSV)

            if(.NOT.IS_ON_myPE_plus2layers(I,J,K)) cycle
            IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells

            ijk = funijk(i,j,k)
            if(fluid_at(ijk)) then
               pSource = PS_MASSFLOW_G(PSV) * (VOL(IJK)/PS_VOLUME(PSV))

               B_M(IJK,0) = B_M(IJK,0) - pSource
               B_MMAX(IJK,0) = max(abs(B_MMAX(IJK,0)), abs(B_M(IJK,0)))
            endif

         enddo
         enddo
         enddo

      enddo PS_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_PP_G
