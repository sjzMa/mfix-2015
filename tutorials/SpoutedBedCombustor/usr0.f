!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR0                                                   C
!  Purpose: This routine is called before the time loop starts and is  C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.  This        C
!           can be used for setting constants and checking errors in   C
!           data.  This routine is not called from an IJK loop, hence  C
!           all indices are undefined.                                 C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
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
      SUBROUTINE USR0
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
      USE constant
      USE funits
      USE param
      USE param1
      USE physprop
      USE toleranc
      USE usr

      IMPLICIT NONE
      INCLUDE 'usrnlst.inc'
!-----------------------------------------------
!
!  Include files defining common blocks here
!
!
!  Define local variables here
!
      DOUBLE PRECISION SUM
!
!  Include files defining statement functions here
!
!
!  Insert user-defined code here
!

!     allocate array declared in usrnlst.inc
      Allocate(  N_Sh (DIMENSION_3, DIMENSION_M) )

        IF(PAFC .EQ. UNDEFINED ) &
          CALL ERROR_ROUTINE ('USR0', 'PAFC not specified', 1, 1)
!
        IF(PAA .NE. UNDEFINED)THEN
          SUM = PAFC + PAA
          IF( .NOT.COMPARE(ONE,SUM) )THEN
            WRITE(UNIT_LOG,'(A,F10.5/A)') &
              ' *** PAFC + PAA = ',SUM,' It should be equal to 1.0'
            CALL EXIT
          ENDIF
        ELSE
          PAA = 1.0 - PAFC
        ENDIF
!
!  Function of the ash-layer void fraction
      f_EP_A = (0.25 + 0.75 * ( 1.0 - PAA )) ** 2.5
      RETURN
      END SUBROUTINE USR0
