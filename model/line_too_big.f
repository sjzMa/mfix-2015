!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: LINE_TOO_BIG (LINE,LINE_LEN,MAXCOL)                    C
!  Purpose: return an error condition if input data is located past    C
!           column MAXCOL in the data input file                       C
!                                                                      C
!  Author: P.Nicoletti                                Date: 25-NOV-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: LINE_TOO_BIG                                    C
!                                                                      C
!  Local variables: L                                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      INTEGER FUNCTION LINE_TOO_BIG (LINE, LINE_LEN, MAXCOL)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                   input data line
      CHARACTER(LEN=*) :: LINE
!
!                   length of input data line
      INTEGER       LINE_LEN
!
!                   maximum column that non-blank charcaters are
!                   are in the input data line
      INTEGER       MAXCOL
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!               loop index
      INTEGER :: L
!-----------------------------------------------
!
!     The function LINE_TOO_BIG returns a value greater than 0 to
!     indicate an error condition (data passed column MAXCOL in LINE)
!
!
      DO L = MAXCOL + 1, LINE_LEN
         IF (LINE(L:L) /= ' ') THEN
            LINE_TOO_BIG = L
            RETURN
         ENDIF
      END DO
      LINE_TOO_BIG = 0
      RETURN
      END FUNCTION LINE_TOO_BIG
