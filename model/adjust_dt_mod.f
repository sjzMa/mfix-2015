      MODULE ADJUST_DT

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: ADJUST_DT(IER, NIT)                                    !
!  Author: M. Syamlal                                 Date: FEB-10-97  !
!                                                                      !
!  Purpose: Automatically adjust time step.                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      LOGICAL FUNCTION ADJUSTDT (IER, NIT)

! Global Variables:
!---------------------------------------------------------------------//
! User defined type of run: new or restart
      use run, only: RUN_TYPE
! User defined aximum number of iterations
      use leqsol, only: MAX_NIT
! User defined: min, max DT and adjustment factor
      use run, only: DT_MIN, DT_MAX, DT_FAC
! Flag: Use stored DT value for advancing TIME
      use run, only: USE_DT_PREV
! Flag: 2nd order time implementation
      use run, only: CN_ON
! Flag: Continue to run at DT_MIN
      use run, only: PERSISTENT_MODE
! Flag: Running in interactive mode
      use run, only: INTERACTIVE_MODE
! Flag: Interupt the code in interactive mode
      use run, only: INTERUPT
! The current number of time steps (value at restart).
      use run, only: NSTEP, NSTEPRST
! Current DT (1/DT) and direction of last change (+/-)
      use run, only: DT, oDT, DT_DIR

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: ZERO, ONE, UNDEFINED


! Module proceedures:
!---------------------------------------------------------------------//
! Routine to break successive time step reductions.
      use interactive, only: CHECK_TIMESTEP_FAIL_RATE
      use error_manager

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Integer flag: 0=Good, 100=initialize, otherwise bad.
      INTEGER, INTENT(INOUT) :: IER
! Number of iterations for current time step
      INTEGER, INTENT(IN) :: NIT


! Local Variables:
!---------------------------------------------------------------------//
! Number of steps in between DT adjustments.
      INTEGER, PARAMETER :: STEPS_MIN = 5
! Number of time steps since last DT adjustment
      INTEGER, SAVE :: STEPS_TOT=0
! number of iterations since last DT adjustment
      INTEGER, SAVE :: NIT_TOT=0
! Iterations per second for last dt
      DOUBLE PRECISION, SAVE :: NIToS=0.0
! Current number of iterations per second
      DOUBLE PRECISION :: NITOS_NEW
! Flag to half/double the current time step
      LOGICAL :: CN_ADJUST_DT
!......................................................................!

! Initialize the function result.
      ADJUSTDT = .FALSE.
      USE_DT_PREV = .FALSE.

! Steady-state simulation.
      IF (DT==UNDEFINED .OR. DT<ZERO) RETURN

! Local flag for adjusting the time step for CN implementation.
      CN_ADJUST_DT = CN_ON .AND. ((RUN_TYPE=='NEW' .AND. NSTEP>1) .OR. &
         (RUN_TYPE/='NEW' .AND. NSTEP >= (NSTEPRST+1)))

! Stop sequential DT reductions by reseting the time-step and returning
! control waiting for new user input.
      IF(INTERACTIVE_MODE) &
         CALL CHECK_TIMESTEP_FAIL_RATE(IER)

! Iterate successfully converged.
!---------------------------------------------------------------------//
      IF(IER == 0) THEN

! Set back the timestep to original size which was halved previously for
! 2nd order accurate time implementation.
         IF(CN_ADJUST_DT) DT = 2.0D0*DT

! Calculate a new DT every STEPS_MIN time steps.
         IF(STEPS_TOT >= STEPS_MIN) THEN
            NITOS_NEW = DBLE(NIT_TOT)/(STEPS_TOT*DT)
            IF (NITOS_NEW > NITOS) DT_DIR = DT_DIR*(-1)
            STEPS_TOT = 0
            NITOS = NITOS_NEW
            NIT_TOT = 0
            IF (DT_DIR > 0) THEN
               IF(NIT < MAX_NIT) DT = MIN(DT_MAX,DT/DT_FAC)
            ELSE
               DT = DT*DT_FAC
               IF(PERSISTENT_MODE) DT = max(DT, DT_MIN)
            ENDIF

! DT was modified. Use the stored DT should be used to update TIME.
            USE_DT_PREV = .TRUE.

! Write the convergence stats to the screen/log file.
            WRITE(ERR_MSG,"('DT=',g11.4,3x,'NIT/s=',A)")  &
               DT, trim(iVal(nint(NITOS)))
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., &
               FOOTER=.FALSE., LOG=.FALSE.)

         ELSE
            STEPS_TOT = STEPS_TOT + 1
            NIT_TOT = NIT_TOT + NIT
         ENDIF
! No need to iterate again
         ADJUSTDT = .FALSE.
! Cut the timestep into half for 2nd order accurate time implementation.
         IF(CN_ADJUST_DT) DT = 0.5d0*DT

! Iterate failed to converge.
!---------------------------------------------------------------------//
      ELSE

! Clear the error flag.
         IER = 0

! Reset the timestep to original size which was halved previously for
! 2nd order accurate time implementation.
         IF(CN_ADJUST_DT) DT = 2.0d0*DT

! Reset counters.
         STEPS_TOT = 0
         NITOS = 0.
         NIT_TOT = 0

! Reduce the step size.
         DT = DT*DT_FAC

! The step size has decreased to the minimum.
         IF (DT_FAC >= ONE) THEN

            IF(PERSISTENT_MODE) THEN
               ADJUSTDT = .FALSE.
            ELSE
               WRITE(ERR_MSG,"(3X,A)") &
                  'DT_FAC >= 1. Recovery not possible!'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE., &
                  HEADER=.FALSE., FOOTER=.FALSE.)
            ENDIF

         ELSEIF (DT > DT_MIN) THEN

            IF(.NOT.INTERUPT) THEN
               WRITE(ERR_MSG,"(3X,'Recovered: Dt=',G12.5,' :-)')") DT
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            ENDIF

            CALL RESET_NEW

! Iterate again with new dt
            ADJUSTDT = .TRUE.

! Cut the timestep for 2nd order accurate time implementation.
            IF(CN_ADJUST_DT) DT = 0.5d0*DT

! Set the return flag stop iterating.
         ELSE

! Prevent DT from dropping below DT_MIN.
            IF(PERSISTENT_MODE) DT = max(DT, DT_MIN)
            ADJUSTDT = .FALSE.
         ENDIF

      ENDIF

! Update ONE/DT variable.
      ODT = ONE/DT

! Break out of the iterate loop if INTERUPT signal is given
      IF(INTERUPT) ADJUSTDT = .FALSE.

      RETURN
      END FUNCTION ADJUSTDT

      END MODULE ADJUST_DT
