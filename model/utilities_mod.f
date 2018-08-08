MODULE utilities

  IMPLICIT NONE

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  function: mfix_isnan                                                !
!  Purpose: check whether argument is NAN                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      LOGICAL FUNCTION mfix_isnan(x)

!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      double precision x
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      CHARACTER(LEN=80) :: notnumber
!-----------------------------------------------

      mfix_isnan = .False.
      WRITE(notnumber,*) x
! To check for NaN's in x, see if x (a real number) contain a letter "N"
! "n" or symbol "?", in which case it is a NaN (Not a Number)

      IF(INDEX(notnumber,'?') > 0 .OR.     &
         INDEX(notnumber,'n') > 0 .OR.     &
         INDEX(notnumber,'N') > 0 ) THEN
        mfix_isnan = .TRUE.
         RETURN
      ENDIF

      RETURN
    END FUNCTION mfix_isnan

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  function: MAX_VEL_INLET                                             C
!  Purpose: Find maximum velocity at inlets.                           C
!                                                                      C
!  Author: S. Benyahia                                Date: 26-AUG-05  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      DOUBLE PRECISION FUNCTION MAX_VEL_INLET()

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE bc
      USE fldvar
      USE geometry
      USE physprop
      USE indices
      USE constant
      USE run
      USE compar
      USE discretelement
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: L, I, J, K, IJK, IJK2, M
!-----------------------------------------------

! initializing
      MAX_VEL_INLET = ZERO

      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN
            IF (BC_TYPE(L) == 'MASS_INFLOW' .OR. BC_TYPE(L) == 'P_INFLOW') THEN

               DO K = BC_K_B(L), BC_K_T(L)
                  DO J = BC_J_S(L), BC_J_N(L)
                     DO I = BC_I_W(L), BC_I_E(L)
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)

                        SELECT CASE (BC_PLANE(L))
                        CASE ('S')
                           IJK2 = JM_OF(IJK)
                           IF( ABS(V_G(IJK2)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(V_G(IJK2))
                        CASE ('N')
                           IF( ABS(V_G(IJK)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(V_G(IJK))
                        CASE ('W')
                           IJK2 = IM_OF(IJK)
                           IF( ABS(U_G(IJK2)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(U_G(IJK2))
                        CASE ('E')
                           IF( ABS(U_G(IJK)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(U_G(IJK))
                        CASE ('B')
                           IJK2 = KM_OF(IJK)
                           IF( ABS(W_G(IJK2)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(W_G(IJK2))
                        CASE ('T')
                           IF( ABS(W_G(IJK)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(W_G(IJK))
                        END SELECT

                       IF (.NOT.DES_CONTINUUM_COUPLED .OR. DES_CONTINUUM_HYBRID) THEN
                          SELECT CASE (BC_PLANE(L))
                           CASE ('S')
                              IJK2 = JM_OF(IJK)
                              DO M = 1, MMAX
                                 IF( ABS(V_s(IJK2, M)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(V_s(IJK2, M))
                              ENDDO
                           CASE ('N')
                              DO M = 1, MMAX
                                IF( ABS(V_s(IJK, M)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(V_s(IJK, M))
                              ENDDO
                           CASE ('W')
                              IJK2 = IM_OF(IJK)
                              DO M = 1, MMAX
                                IF( ABS(U_s(IJK2, M)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(U_s(IJK2, M))
                              ENDDO
                           CASE ('E')
                              DO M = 1, MMAX
                                IF( ABS(U_s(IJK, M)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(U_s(IJK, M))
                              ENDDO
                           CASE ('B')
                              IJK2 = KM_OF(IJK)
                              DO M = 1, MMAX
                                IF( ABS(W_s(IJK2, M)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(W_s(IJK2, M))
                              ENDDO
                           CASE ('T')
                              DO M = 1, MMAX
                                IF( ABS(W_s(IJK, M)) > MAX_VEL_INLET ) MAX_VEL_INLET = ABS(W_s(IJK, M))
                              ENDDO
                           END SELECT
                        ENDIF   ! end if (.not.des_continuum_coupled .or. des_continuum_hybrid)

                     ENDDO
                  ENDDO
               ENDDO

           ENDIF
         ENDIF
      ENDDO

      RETURN
      END FUNCTION MAX_VEL_INLET


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function: CHECK_VEL_BOUND()                                         C
!  Purpose: Check velocities upper bound to be less than speed of      C
!           sound                                                      C
!                                                                      C
!  Author: S. Benyahia                                Date: 25-AUG-05  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      LOGICAL FUNCTION CHECK_VEL_BOUND ()

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE bc
      USE geometry
      USE physprop
      USE indices
      USE run
      USE toleranc
      USE compar
      USE mpi_utility
      USE discretelement
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: M
! Indices
      INTEGER :: IJK
      LOGICAL :: ALL_IS_ERROR
!-----------------------------------------------

!!$omp   parallel do private(IJK)
! initializing
      CHECK_VEL_BOUND = .FALSE.
      ALL_IS_ERROR    = .FALSE.

LOOP_FLUID : DO IJK = IJKSTART3, IJKEND3

         IF (FLUID_AT(IJK)) THEN
            IF(ABS(U_G(IJK)) > MAX_INLET_VEL .OR. &
               ABS(V_G(IJK)) > MAX_INLET_VEL .OR. &
               ABS(W_G(IJK)) > MAX_INLET_VEL) THEN
               CHECK_VEL_BOUND = .TRUE.
               WRITE(*,1000) MAX_INLET_VEL, I_OF(IJK), J_OF(IJK), K_OF(IJK), &
                             EP_g(IJK), U_G(IJK), V_G(IJK), W_G(IJK)
               EXIT LOOP_FLUID
            ENDIF

            IF (.NOT.DES_CONTINUUM_COUPLED .OR. DES_CONTINUUM_HYBRID) THEN
               DO M = 1, MMAX
                 IF(ABS(U_S(IJK,M)) > MAX_INLET_VEL .OR. &
                    ABS(V_S(IJK,M)) > MAX_INLET_VEL .OR. &
                    ABS(W_S(IJK,M)) > MAX_INLET_VEL) THEN
                   CHECK_VEL_BOUND = .TRUE.
                   WRITE(*,1010) MAX_INLET_VEL, I_OF(IJK), J_OF(IJK), K_OF(IJK), M, &
                                 EP_s(IJK, M), U_S(IJK,M), V_S(IJK,M), W_S(IJK,M)
                   EXIT LOOP_FLUID
                 ENDIF
               ENDDO
            ENDIF   ! end if(.not.des_continuum_coupled or des_continuum_hybrid)
         ENDIF

      ENDDO LOOP_FLUID

      CALL GLOBAL_ALL_OR(CHECK_VEL_BOUND, ALL_IS_ERROR)
      IF(ALL_IS_ERROR) CHECK_VEL_BOUND = .TRUE.

      RETURN
 1000 FORMAT(1X,'Message from: CHECK_VEL_BOUND',/&
            'WARNING: velocity higher than maximum allowed velocity: ', &
            G12.5, '(to change this adjust the scale factor MAX_INLET_VEL_FAC)'/&
            'in this cell: ','I = ',I4,2X,' J = ',I4,2X,' K = ',I4, /&
            '  ','Epg = ', G12.5, 'Ug = ', G12.5, 'Vg = ', G12.5, 'Wg = ', G12.5)
 1010 FORMAT(1X,'Message from: CHECK_VEL_BOUND',/&
            'WARNING: velocity higher than maximum allowed velocity: ', &
            G12.5,/&
            'in this cell: ','I = ',I4,2X,' J = ',I4,2X,' K = ',I4,' M = ',I4, /&
            '  ','Eps = ', G12.5,'Us = ', G12.5, 'Vs = ', G12.5, 'Ws = ', G12.5)

      END FUNCTION CHECK_VEL_BOUND



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function name: SEEK_COMMENT (LINE_MAXCOL)                           C
!  Purpose: determine if (and where) a comment character appears       C
!           in a data input line                                       C
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
!  Variables modified: SEEK_COMMENT                                    C
!                                                                      C
!  Local variables: DIM_COMMENT, COMMENT_CHAR, L, COMMENT, L2          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      INTEGER FUNCTION SEEK_COMMENT (LINE, MAXCOL)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                   input data line
      CHARACTER(len=*) LINE
!
!                   maximum column of input data line to search
      INTEGER       MAXCOL
!
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
!                   the number of designated comment characters
      INTEGER, PARAMETER :: DIM_COMMENT = 2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                   loop indicies
      INTEGER :: L, L2
!
!                   the comment characters
      CHARACTER, DIMENSION(DIM_COMMENT) :: COMMENT_CHAR
!-----------------------------------------------
!
!     The function SEEK_COMMENT returns the index to where a comment
!     character was found in the input data line.  Equals MAXCOL + 1
!     if no-comment characters in the line
!
!
      DATA COMMENT_CHAR/'#', '!'/
!
      DO L = 1, MAXCOL
         DO L2 = 1, DIM_COMMENT
            IF (LINE(L:L) == COMMENT_CHAR(L2)) THEN
               SEEK_COMMENT = L
               RETURN
            ENDIF
         END DO
      END DO
      SEEK_COMMENT = MAXCOL + 1
!
      RETURN
      END FUNCTION SEEK_COMMENT

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function name: SEEK_END (LINE, MAXCOL)                              C
!  Purpose: determine where trailing blanks begin in a line            C
!                                                                      C
!  Author: P.Nicoletti, M. Syamlal                    Date: 7-AUG-92   C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: SEEK_END                                        C
!                                                                      C
!  Local variables: L                                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      INTEGER FUNCTION SEEK_END (LINE, MAXCOL)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                   maximum column of input data line to search
      INTEGER MAXCOL
!
!                   input data line
      CHARACTER LINE*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: L
!-----------------------------------------------
!
!     The function SEEK_END returns the index to where the last
!     character was found in the input data line.  Equals MAXCOL
!     if no trailing blank characters in the line
!
!
      SEEK_END = 0
      DO L = 1, MAXCOL
         IF (LINE(L:L) /= ' ') SEEK_END = L
      END DO
      RETURN
      END FUNCTION SEEK_END

!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function name: LINE_TOO_BIG (LINE,LINE_LEN,MAXCOL)                  C
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

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Function: BLANK_LINE                                                !
! Author: P. Nicoletti                                Date: 25-NOV-91  !
!                                                                      !
! Purpose: Return .TRUE. if a line contains no input or only spaces.   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      LOGICAL FUNCTION BLANK_LINE (line)

      IMPLICIT NONE

      CHARACTER :: LINE*(*)

      INTEGER :: L

      BLANK_LINE = .FALSE.
      DO L=1, len(line)
         IF(line(L:L)/=' ' .and. line(L:L)/='    ')RETURN
      ENDDO

      BLANK_LINE = .TRUE.
      RETURN
      END FUNCTION BLANK_LINE

END MODULE utilities
