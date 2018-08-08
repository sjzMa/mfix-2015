!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_CONTINUITY                                        C
!  Purpose: Solve for solids bulk density (i.e., material density      C
!           multiplied by volume fraction).                            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 2-JUL-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOLVE_CONTINUITY(M,IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE physprop
      USE geometry
      USE fldvar
      USE indices
      USE residual
      USE cont
      USE leqsol
      Use ambm
      USE ur_facs
      use run
      use ps

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! phase index
      INTEGER, INTENT(IN) :: M
! error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! solids volume fraction residual
      DOUBLE PRECISION :: RESs
! linear equation solver method and iterations
      INTEGER :: LEQM, LEQI

! temporary use of global arrays:
! Septadiagonal matrix A_m, vector B_m
!      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------

      call lock_ambm


      IF (M==0) THEN
! solve gas phase continuity equation. note that this branch will never
! be entered given that the calling subroutine (iterate) only calls
! solve_continuity when M>1

! initializing
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, 0)

! forming the matrix equation
         CALL CONV_ROP_G (A_M, B_M)
         CALL SOURCE_ROP_G (A_M, B_M)

! calculating the residual
         CALL CALC_RESID_C (ROP_G, A_M, B_M, 0, NUM_RESID(RESID_RO,0), &
            DEN_RESID(RESID_RO,0), RESID(RESID_RO,0), MAX_RESID(&
            RESID_RO,0), IJK_RESID(RESID_RO,0))

!         call check_ab_m(a_m, b_m, 0, .true., ier)
!         call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!         write(*,*) resid(resid_ro, 0), max_resid(resid_ro, 0), &
!            ijk_resid(resid_ro, 0)
!         call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0),&
!            1, DO_K, ier)

! solving gas continuity
         CALL ADJUST_LEQ (RESID(RESID_RO,0), LEQ_IT(2), LEQ_METHOD(2),&
            LEQI, LEQM)
         CALL SOLVE_LIN_EQ ('ROP_g', 2, ROP_G, A_M, B_M, 0, LEQI, &
            LEQM, LEQ_SWEEP(2), LEQ_TOL(2), LEQ_PC(2), IER)
         CALL ADJUST_ROP (ROP_G)
!        call out_array(ROP_g, 'rop_g')

      ELSE
! solve solids phase M continuity equation.

! initializing
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M)

! forming the matrix equation
         CALL CONV_ROP_S (A_M, B_M, M)
         CALL SOURCE_ROP_S (A_M, B_M, M)
         IF(POINT_SOURCE) CALL POINT_SOURCE_ROP_S (B_M, M)

         CALL CALC_RESID_C (ROP_S(1,M), A_M, B_M, M, &
            NUM_RESID(RESID_RO,M), DEN_RESID(RESID_RO,M), &
            RESID(RESID_RO,M), MAX_RESID(RESID_RO,M), &
            IJK_RESID(RESID_RO,M))
         RESS = RESID(RESID_RO,M)

!          call check_ab_m(a_m, b_m, m, .true., ier)
!          write(*,*)'solve_cont=',resid(resid_ro, m),max_resid(resid_ro, m),&
!        ijk_resid(resid_ro, m),m
!          call write_ab_m(a_m, b_m, ijkmax2, m, ier)
!          call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, M), 1,
!     &      DO_K, ier)
!
         CALL ADJUST_LEQ (RESID(RESID_RO,M), LEQ_IT(2), LEQ_METHOD(2),&
            LEQI, LEQM)
         CALL SOLVE_LIN_EQ ('ROP_s', 2, ROP_S(1,M), A_M, B_M, M, LEQI,&
            LEQM,LEQ_SWEEP(2), LEQ_TOL(2), LEQ_PC(2), IER)
         CALL ADJUST_ROP (ROP_S(1,M))
!          call out_array(rop_s(1,m), 'rop_s')
      ENDIF

      call unlock_ambm


      RETURN
      END SUBROUTINE SOLVE_CONTINUITY
