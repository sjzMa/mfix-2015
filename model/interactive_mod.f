      MODULE INTERACTIVE

      use param1, only: UNDEFINED_I, UNDEFINED, UNDEFINED_C
      use error_manager

      IMPLICIT NONE

      PUBLIC :: INTERACT
      PUBLIC :: CHECK_INTERACT_ITER
      PUBLIC :: CHECK_TIMESTEP_FAIL_RATE
      PUBLIC :: INIT_INTERACTIVE_MODE

      PRIVATE

! Actions supported in interactive mode.
!----------------------------------------------------------------------!
      INTEGER :: ACTION
! Do nothing (release interactive control)
      INTEGER, PARAMETER :: NULL_ENUM=0
! Pause the simulation to await further commands.
      INTEGER, PARAMETER :: WAIT_ENUM=1
! Exit the simulation cleanly.
      INTEGER, PARAMETER :: EXIT_ENUM=2
! Exit the simulation without save current state.
      INTEGER, PARAMETER :: ABORT_ENUM=3

! Take N time-steps then wait for next command.
      INTEGER, PARAMETER :: NSTEPS_ENUM=10
! Run to a user-specified time then wait for next command.
      INTEGER, PARAMETER :: WAIT_AT_ENUM=11
! Create a RES file backup.
      INTEGER, PARAMETER :: BCKUP_RES_ENUM=20
! Wite a RES file
      INTEGER, PARAMETER :: WRITE_RES_ENUM=21
! Read a RES file
      INTEGER, PARAMETER :: READ_RES_ENUM=22
! Refresh the simulation by processing the mfix.dat file.
      INTEGER, PARAMETER :: REINITIALIZE_ENUM=100


! Generic variables for interactive communication.
!---------------------------------------------------------------------//
      INTEGER, PARAMETER :: ACT_DIMN=50
! Integer input for action items
      INTEGER :: ACTION_INT(ACT_DIMN)
! Float input for action items
      DOUBLE PRECISION :: ACTION_DP(ACT_DIMN)
! String input for action items
      CHARACTER(LEN=256), DIMENSION(ACT_DIMN) :: ACTION_STR


! Internal data storage variables.
!---------------------------------------------------------------------//
      INTEGER :: USER_DEFINED_STEPS = UNDEFINED_I
      INTEGER :: USER_DEFINED_NITS = UNDEFINED_I
      INTEGER :: STEP_CNT

      DOUBLE PRECISION :: USER_WAIT_AT_TIME = UNDEFINED

      INTEGER :: HISTORY_POS
      DOUBLE PRECISION :: FAILRATE, MAX_FAILRATE
      INTEGER, ALLOCATABLE :: TIMESTEP_HISTORY(:)

      CONTAINS
!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: INTERACT                                                !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE INTERACT(pEXIT_SIGNAL, pABORT_SIGNAL)

      use compar, only: myPE, PE_IO

      use mpi_utility
      use error_manager

      LOGICAL :: INTERACTING

      LOGICAL, INTENT(INOUT) :: pEXIT_SIGNAL
      LOGICAL, INTENT(OUT) :: pABORT_SIGNAL

      pABORT_SIGNAL = .FALSE.

      write(*,*) 'In the INTERACT ROUTINE'

      CALL CHECK_INTERACT_SIGNAL(INTERACTING)

      DO WHILE(INTERACTING)

         CALL GET_INTERACTIVE_DATA

         SELECT CASE(ACTION)
         CASE(NULL_ENUM)
            INTERACTING = .FALSE.

         CASE(WAIT_ENUM)
            INTERACTING = .TRUE.

         CASE(EXIT_ENUM)
            pEXIT_SIGNAL = .TRUE.
            INTERACTING = .FALSE.

         CASE(ABORT_ENUM)
            pABORT_SIGNAL = .TRUE.
            INTERACTING = .FALSE.

         CASE(NSTEPS_ENUM)
            CALL CHECK_INTERACTIVE_NSTEPS(INTERACTING)

         CASE(WAIT_AT_ENUM)
            CALL CHECK_INTERACTIVE_WAIT_AT(INTERACTING)

         CASE(BCKUP_RES_ENUM)
            CALL CHECK_RES_ACTION
            INTERACTING = .TRUE.

         CASE(WRITE_RES_ENUM)
            CALL CHECK_RES_ACTION
            INTERACTING = .TRUE.

         CASE(READ_RES_ENUM)
            CALL CHECK_RES_ACTION
            INTERACTING = .TRUE.

         CASE(REINITIALIZE_ENUM)
            CALL REINITIALIZE
            INTERACTING = .TRUE.

         CASE DEFAULT
            CALL UNKNOWN_INTERACTIVE_CMD
            INTERACTING = .TRUE.

         END SELECT

         IF(INTERACTING) CALL INTERACTIVE_WAIT
      ENDDO

      RETURN
      END SUBROUTINE INTERACT


!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: INTERACT                                                !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_INTERACT_SIGNAL(pINTERACT)

      use compar, only: myPE, PE_IO
      use param1, only: UNDEFINED_I, UNDEFINED
      use run, only: INTERUPT
      use run, only: TIME, DT

      use mpi_utility

      LOGICAL, INTENT(OUT) :: pINTERACT

! Root checks if there is there is information from the user.
      IF(myPE == PE_IO) THEN

         IF(INTERUPT) THEN
            pINTERACT = .TRUE.
            INTERUPT = .FALSE.
            CALL INTERACTIVE_WAIT
         ELSE
            INQUIRE(file="interact.dat", exist=pINTERACT)
         ENDIF

         IF(.NOT.pINTERACT .AND. USER_DEFINED_STEPS /= UNDEFINED_I) THEN
            STEP_CNT = STEP_CNT + 1
            IF(STEP_CNT >= USER_DEFINED_STEPS) THEN
               STEP_CNT = UNDEFINED_I
               USER_DEFINED_STEPS = UNDEFINED_I
               pINTERACT = .TRUE.
               CALL INTERACTIVE_WAIT
            ENDIF
         ENDIF

         IF(.NOT.pINTERACT .AND. USER_WAIT_AT_TIME /= UNDEFINED) THEN
            IF(TIME + 0.1d0*DT >= USER_WAIT_AT_TIME) THEN
               USER_WAIT_AT_TIME = UNDEFINED
               pINTERACT = .TRUE.
               CALL INTERACTIVE_WAIT
            ENDIF
         ENDIF


      ENDIF

      CALL BCAST(pINTERACT, PE_IO)

      RETURN
      END SUBROUTINE CHECK_INTERACT_SIGNAL


!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: CHECK_INTERACTIVE_NSTEPS                                !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_INTERACT_ITER(iMUSTIT)

      use run, only: INTERACTIVE_NITS

      INTEGER, INTENT(OUT) :: iMUSTIT
      LOGICAL :: INTERACTING

      INTERACTIVE_NITS = INTERACTIVE_NITS - 1

      IF(INTERACTIVE_NITS > 0) THEN
         iMUSTIT = 1
         RETURN
      ENDIF

      iMUSTIT = 0
      INTERACTIVE_NITS = UNDEFINED_I

      CALL INTERACTIVE_WAIT

      INTERACTING = .TRUE.
      DO WHILE(INTERACTING)

         CALL GET_INTERACTIVE_DATA

         SELECT CASE(ACTION)
         CASE(NULL_ENUM)
            INTERACTING = .FALSE.

         CASE(WAIT_ENUM)
            INTERACTING = .TRUE.

         CASE(NSTEPS_ENUM)
            CALL CHECK_INTERACTIVE_NSTEPS(INTERACTING)
            IF(INTERACTIVE_NITS /= UNDEFINED_I) iMUSTIT = 1

         CASE(WAIT_AT_ENUM)
            CALL CHECK_INTERACTIVE_WAIT_AT(INTERACTING)

         CASE DEFAULT
            CALL UNAVAILABLE_INTERACTIVE_CMD
            INTERACTING = .TRUE.

         END SELECT

         IF(INTERACTING) CALL INTERACTIVE_WAIT
      ENDDO


      RETURN
      END SUBROUTINE CHECK_INTERACT_ITER

!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: CHECK_INTERACTIVE_WAIT_AT                               !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_INTERACTIVE_WAIT_AT(pINTERACT)

      use compar, only: myPE, PE_IO

      LOGICAL, INTENT(OUT) :: pINTERACT

      IF(ACTION_DP(1) == UNDEFINED_I) THEN
         IF(myPE == PE_IO) THEN
            WRITE(*,*) ' '
            WRITE(*,*) ' ACTION: WAIT_AT is not fully specified.'
            WRITE(*,*) ' I will wait until you fix it...'
            WRITE(*,*) ' '
         ENDIF
         pINTERACT = .TRUE.
      ELSE
         USER_WAIT_AT_TIME = ACTION_DP(1)
         pINTERACT = .FALSE.
      ENDIF

      RETURN
      END SUBROUTINE CHECK_INTERACTIVE_WAIT_AT



!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: UNKNOWN_INTERACTIVE_CMD                                 !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE UNKNOWN_INTERACTIVE_CMD

      use compar, only: myPE, PE_IO

      IF(myPE == PE_IO) THEN
         WRITE(*,*) ' '
         WRITE(*,*) ' '
         WRITE(*,*) ' Unknown interactive command.'
         WRITE(*,*) ' Enter a different command.'
         WRITE(*,*) ' '
      ENDIF

      RETURN
      END SUBROUTINE UNKNOWN_INTERACTIVE_CMD

!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: UNKNOWN_AVAILABLE_CMD                                   !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE UNAVAILABLE_INTERACTIVE_CMD

      use compar, only: myPE, PE_IO

      IF(myPE == PE_IO) THEN
         WRITE(*,*) ' '
         WRITE(*,*) ' '
         WRITE(*,*) ' Unavailable interactive command.'
         WRITE(*,*) ' Enter a different command.'
         WRITE(*,*) ' '
      ENDIF

      RETURN
      END SUBROUTINE UNAVAILABLE_INTERACTIVE_CMD

!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: CHECK_INTERACTIVE_NSTEPS                                !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_INTERACTIVE_NSTEPS(pINTERACT)

      use run, only: INTERACTIVE_NITS
      use compar, only: myPE, PE_IO

      LOGICAL, INTENT(OUT) :: pINTERACT

      IF(ACTION_INT(1) == UNDEFINED_I) THEN
         IF(myPE == PE_IO) THEN
            WRITE(*,*) ' '
            WRITE(*,*) ' ACTION: NSTEPS is not fully specified.'
            WRITE(*,*) ' I will wait until you fix it...'
            WRITE(*,*) ' '
         ENDIF
         pINTERACT = .TRUE.
      ELSE
         USER_DEFINED_STEPS = ACTION_INT(1)
         STEP_CNT = 0
         pINTERACT = .FALSE.
      ENDIF

      IF(USER_DEFINED_STEPS > 0) RETURN

      IF(ACTION_INT(2) == UNDEFINED_I) THEN
         IF(myPE == PE_IO)THEN
            WRITE(*,*) ' '
            WRITE(*,*) ' ACTION: NSTEPS=0 is not fully specified.'
            WRITE(*,*) ' I will wait until you fix it...'
            WRITE(*,*) ' '
         ENDIF
         pINTERACT=.TRUE.
      ELSE
         INTERACTIVE_NITS = ACTION_INT(2)
      ENDIF


      RETURN
      END SUBROUTINE CHECK_INTERACTIVE_NSTEPS

!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: CHECK_RES_ACTION                                        !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_RES_ACTION

      use run, only: RUN_NAME, TIME

      use compar, only: myPE, PE_IO

      LOGICAL :: EXISTS

      CHARACTER(LEN=128) :: tFNAME
      CHARACTER(LEN=256) :: CMD

      DOUBLE PRECISION :: lTIME

      IF(ACTION_STR(1) == UNDEFINED_C) THEN
         IF(ACTION == BCKUP_RES_ENUM .OR. &
            ACTION == READ_RES_ENUM) THEN
            IF(myPE == PE_IO) THEN
               WRITE(*,*) ' '
               WRITE(*,*) ' ACTION: BCKUP_RES is not fully specified.'
               WRITE(*,*) ' I will wait until you fix it...'
               WRITE(*,*) ' '
            ENDIF
            RETURN
         ENDIF
      ELSEIF(myPE == PE_IO) THEN
         tFNAME=''; WRITE(tFNAME,"(A,'.RES')") trim(RUN_NAME)
         WRITE(*,*) 'Create RES backup: ',trim(ACTION_STR(1))
         CMD=''; WRITE(CMD,"('cp ',2(1x,A))") &
            trim(tFNAME), trim(ACTION_STR(1))

         write(*,*) trim(CMD)
         CALL SYSTEM(trim(CMD))
      ENDIF

      IF(ACTION == WRITE_RES_ENUM) THEN
         IF(myPE==PE_IO) WRITE(*,*) 'Writing RES file.'
         CALL WRITE_RES1

      ELSEIF(ACTION == READ_RES_ENUM) THEN
         IF(ACTION_STR(2) == UNDEFINED_C) THEN
            IF(myPE == PE_IO) THEN
               WRITE(*,*) ' '
               WRITE(*,*) ' ACTION: READ_RES is not fully specified.'
               WRITE(*,*) ' I will wait until you fix it...'
               WRITE(*,*) ' '
            ENDIF
            RETURN
         ELSE
            INQUIRE(FILE=trim(ACTION_STR(2)), EXIST=EXISTS)
            IF(.NOT.EXISTS)THEN
               IF(myPE == PE_IO) THEN
                  WRITE(*,*) ' '
                  WRITE(*,*) ' ACTION: READ_RES is missing.'
                  WRITE(*,*) ' I will wait until you fix it...'
                  WRITE(*,*) ' '
                  RETURN
               ENDIF
            ENDIF
         ENDIF
         IF(myPE==PE_IO) WRITE(*,*) 'Reading RES file.'
         tFNAME=''; WRITE(tFNAME,"(A,'.RES')") trim(RUN_NAME)
         CMD=''; WRITE(CMD,"('cp ',2(1x,A))") &
            trim(ACTION_STR(2)), trim(tFNAME)
         write(*,*) trim(CMD)
         CALL SYSTEM(trim(CMD))

         lTIME = TIME
         CALL READ_RES1
         TIME = lTIME
      ENDIF

      RETURN
      END SUBROUTINE CHECK_RES_ACTION

!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: GET_INTERACTIVE_DATA                                    !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE GET_INTERACTIVE_DATA

      CALL INTERACTIVE_DATA_INIT
      CALL INTERACTIVE_DATA_READ
      CALL INTERACTIVE_DATA_FINL

      RETURN
      END SUBROUTINE GET_INTERACTIVE_DATA


!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: INTERACTIVE_DATA_INIT                                   !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE INTERACTIVE_DATA_INIT

      use error_manager

      ACTION = NULL_ENUM

      ACTION_INT = UNDEFINED_I
      ACTION_DP = UNDEFINED
      ACTION_STR = UNDEFINED_C

      RETURN
      END SUBROUTINE INTERACTIVE_DATA_INIT


!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: INTERACTIVE_DATA_READ                                   !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE INTERACTIVE_DATA_READ

      use funits, only: UNIT_DAT
      use error_manager

      INTEGER :: IOS
      CHARACTER(LEN=512) :: LINE_STRING, STRING

      NAMELIST / INTERACT_INPUT / ACTION, &
         ACTION_INT, ACTION_DP, ACTION_STR

      OPEN(UNIT=UNIT_DAT, FILE='interact.dat', STATUS='OLD', IOSTAT=IOS)

      READ_LP: DO
         READ(UNIT_DAT,"(A)",IOSTAT=IOS) LINE_STRING
         IF(IOS < 0) EXIT READ_LP
! Setup the namelist format as a string.
         STRING=''; STRING = '&INTERACT_INPUT '//&
            trim(adjustl(LINE_STRING))//'/'
! Process the input line.
         READ(STRING, NML=INTERACT_INPUT, IOSTAT=IOS)
      ENDDO READ_LP

      CLOSE(UNIT_DAT)

      RETURN
      END SUBROUTINE INTERACTIVE_DATA_READ


!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: INIT_INTERACTIVE_DATA                                   !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE INTERACTIVE_DATA_FINL

      use compar, only: myPE, PE_IO
      use mpi_utility

      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
      IF(myPE == PE_IO) CALL SYSTEM('rm interact.dat')


      RETURN
      END SUBROUTINE INTERACTIVE_DATA_FINL

!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: INTERACTIVE_WAIT                                        !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE INTERACTIVE_WAIT

      use compar, only: myPE, PE_IO
      use mpi_utility

      use run, only: INTERUPT

      LOGICAL :: FOUND

      IF(myPE == PE_IO) THEN
         DO; INQUIRE(file="interact.dat", exist=FOUND)
            IF(FOUND) EXIT
            IF(.NOT.INTERUPT) EXIT
            CALL SLEEP(1)
            write(*,*) 'I am waiting...'
         ENDDO
      ENDIF

      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)

      RETURN
      END SUBROUTINE INTERACTIVE_WAIT


!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: UPDATE_FAILURE_RATE                                     !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_TIMESTEP_FAIL_RATE(IER)

      use run, only: INTERUPT
      use run, only: TIMESTEP_FAIL_RATE

      INTEGER, INTENT(IN) :: IER
      INTEGER :: WINDOW, NEW, OLD

      WINDOW = TIMESTEP_FAIL_RATE(2)

      NEW = merge(0,1, IER == 0)
      OLD = TIMESTEP_HISTORY(HISTORY_POS)

      TIMESTEP_HISTORY(HISTORY_POS) = NEW
      HISTORY_POS = MOD(HISTORY_POS, WINDOW) + 1

      FAILRATE = FAILRATE + DBLE(NEW-OLD)/DBLE(WINDOW)

      IF(FAILRATE >= MAX_FAILRATE) THEN
         INTERUPT = .TRUE.
         FAILRATE = 0.0
         WRITE(ERR_MSG, 1000) trim(iVal(TIMESTEP_FAIL_RATE(1))), &
            trim(iVal(TIMESTEP_FAIL_RATE(2)))
         CALL FLUSH_ERR_MSG(HEADER=.FALSE.)
      ENDIF

 1000 FORMAT(2/,70('*'),/'DT reduced ',A,' times over the last ',A,1x, &
         'steps. The current time step was',/'reset. User action is ', &
         'required.')

      RETURN
      END SUBROUTINE CHECK_TIMESTEP_FAIL_RATE


!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: INIT_INTERACTIVE_MODE                                   !
!  Author: J.Musser                                   Date: APR-15-02  !
!                                                                      !
!  Purpose: Interact with MFIX via text files at runtime.              !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE INIT_INTERACTIVE_MODE

      use run, only: TIMESTEP_FAIL_RATE

      ALLOCATE(TIMESTEP_HISTORY(TIMESTEP_FAIL_RATE(2)))
      TIMESTEP_HISTORY = 0

      HISTORY_POS = 1
      FAILRATE = 0.0
      MAX_FAILRATE = dble(TIMESTEP_FAIL_RATE(1)) / &
         dble(TIMESTEP_FAIL_RATE(2))

      RETURN
      END SUBROUTINE INIT_INTERACTIVE_MODE



      END MODULE INTERACTIVE
