!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ERROR_ROUTINE (CALL_ROUTINE,MESSAGE,ACTION_CODE,       C
!                              MESSAGE_CODE)                           C
!  Purpose:  Assist in printing error messages during input processing C
!                                                                      C
!  Author: P.NICOLETTI                                Date: 25-NOV-91  C
!  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 29-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:  Added myPE identifier in outputs and also replace STOP    C
!            with mfix_exit() to abort all processors                  C
!                                                                      C
!  Author:   Aeolus Res. Inc.                         Date: 04-SEP-99  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: ABORT_CONT                                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE ERROR_ROUTINE(CALL_ROUTINE,MESSAGE,ACTION_CODE,MESSAGE_CODE)

      USE funits
      USE compar
      USE mpi_utility

      IMPLICIT NONE

! Action to be taken after error messgase is given.
! > 0 :: Continute
! > 1 :: Abort
      INTEGER, intent(in) :: ACTION_CODE
! Message code
! - 1 :: Write header, message, and footer.
! - 2 :: Write only the header and message.
! - 3 :: Write only the message and footer.
      INTEGER, intent(in) ::  MESSAGE_CODE

! Name of routine calling ERROR_ROUTINE. Used in constructing
! the error message header.
      CHARACTER, intent(in) :: CALL_ROUTINE*(*)

! Message to be written.
      CHARACTER, intent(in) :: MESSAGE*(*)

! Local Variables
! String used to format integers to characters.
      CHARACTER(len=16) :: int2char

! Write the header information and passed messaged.
      IF (MESSAGE_CODE /= 3 .AND. DMP_LOG) THEN
         write(*,1000); write(UNIT_LOG,1000)
         IF(numPEs > 1) THEN
            int2char=''; write(int2char,*) myPE
            write(*,1001) trim(adjustl(int2Char)), CALL_ROUTINE
            write(UNIT_LOG,1001) trim(adjustl(int2Char)), CALL_ROUTINE
         ELSE
            write(*,1002) CALL_ROUTINE
            write(UNIT_LOG,1002) CALL_ROUTINE
         ENDIF

         WRITE (*, 1005) MESSAGE
         WRITE (UNIT_LOG, 1005) MESSAGE
      ENDIF

! WRITE OUT TRAILER INFO, UNLESS MESSAGE_CODE = 2
      IF (MESSAGE_CODE /= 2 .AND. DMP_LOG) THEN
         IF(ACTION_CODE == 0) THEN
            WRITE (*, 1100)
            WRITE (UNIT_LOG, 1100)
         ELSE
            WRITE (*, 1101)
            WRITE (UNIT_LOG, 1101)
         ENDIF
      ENDIF

      IF (ACTION_CODE == 1) CALL MFIX_EXIT(myPE)

      RETURN

 1000 FORMAT(2/,1X,70('*'))
 1001 FORMAT(1X,'(PE ',A,'): From : ',A)
 1002 FORMAT(1X,'From : ',A)
 1005 FORMAT(1X,'Message : ',A)

 1100 FORMAT(1X,70('*'),2/)
 1101 FORMAT(/1X,'Aborting execution.',/1X,70('*'),2/)

      END SUBROUTINE ERROR_ROUTINE
