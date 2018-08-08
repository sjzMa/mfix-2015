!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_CORNER_CELLS(IER)                                  C
!  Purpose: Identify wall cells with more than one fulid cell as       C
!           a neighbor.  No heat mass, momentum, or energy transfer    C          is allowed to such cells to avoid ambiguity.
!                                                                      C
!  Author: M. Syamlal                                 Date: 08-JUL-98  C
!  Reviewer:                                          Date:            C
!                                                                      C
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
      SUBROUTINE GET_CORNER_CELLS()
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE compar
      USE corner
      USE functions
      USE funits
      USE geometry
      USE indices
      USE machine, only: start_log, end_log
      USE matrix
      USE param
      USE param1
      USE physprop
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Loop index
      INTEGER          L

!
!                      indices
      INTEGER          IJK, IMJK, IJMK, IJKM, IPJK, IJPK, &
                      IJKP
!
!                      number of faces adjacent to a fluid cell
      INTEGER          NUM
!
!                      fluid face location, whether not a corner
      LOGICAL          dir(-3:3), NotCorner
!
!-----------------------------------------------

      NCORN = 0
!
      DO IJK = ijkstart3, ijkend3
         IF (WALL_AT(IJK).AND..NOT.CYCLIC_AT(IJK)) THEN
!
!----------------------------------------------------------------
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKP = KP_OF(IJK)
!----------------------------------------------------------------
            NUM = 0
!
            IF (FLUID_AT(IMJK)) THEN
               NUM = NUM + 1
               DIR(W) = .TRUE.
            ELSE
               DIR(W) = .FALSE.
            ENDIF
!
            IF (FLUID_AT(IPJK)) THEN
               NUM = NUM + 1
               DIR(E) = .TRUE.
            ELSE
               DIR(E) = .FALSE.
            ENDIF
!
            IF (FLUID_AT(IJMK)) THEN
               NUM = NUM + 1
               DIR(S) = .TRUE.
            ELSE
               DIR(S) = .FALSE.
            ENDIF
!
            IF (FLUID_AT(IJPK)) THEN
               NUM = NUM + 1
               DIR(N) = .TRUE.
            ELSE
               DIR(N) = .FALSE.
            ENDIF
!
            IF (FLUID_AT(IJKM)) THEN
               NUM = NUM + 1
               DIR(B) = .TRUE.
            ELSE
               DIR(B) = .FALSE.
            ENDIF
!
            IF (FLUID_AT(IJKP)) THEN
               NUM = NUM + 1
               DIR(T) = .TRUE.
            ELSE
               DIR(T) = .FALSE.
            ENDIF
!
            IF (NUM > 1) THEN
!
!
               NOTCORNER = .TRUE.
!
!           check for single cell thick internal walls
               IF (DIR(W) .AND. DIR(E) .OR. DIR(S) .AND. DIR(N) .OR. DIR(T)&
                   .AND. DIR(B)) THEN
!
                  IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
               ENDIF
!
!           check for corner cells
!
               IF (DIR(E)) THEN
!
                  IF (DIR(N)) THEN
                     IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                  ENDIF
!
                  IF (DIR(S)) THEN
                     IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                  ENDIF
!
                  IF (DO_K) THEN
                     IF (DIR(T)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
!
                     IF (DIR(B)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
                  ENDIF
!
               ENDIF
!
               IF (DIR(W)) THEN
!
                  IF (DIR(N)) THEN
                     IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                  ENDIF
!
                  IF (DIR(S)) THEN
                     IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                  ENDIF
!
                  IF (DO_K) THEN
                     IF (DIR(T)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
!
                     IF (DIR(B)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
                  ENDIF
!
               ENDIF
!
               IF (DIR(N)) THEN
!
!
                  IF (DO_K) THEN
                     IF (DIR(T)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
!
                     IF (DIR(B)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
                  ENDIF
!
               ENDIF
!
               IF (DIR(S)) THEN
!
!
                  IF (DO_K) THEN
                     IF (DIR(T)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
!
                     IF (DIR(B)) THEN
                        IF (NOTCORNER) CALL ADDCORN (NOTCORNER, NCORN)
                     ENDIF
                  ENDIF
!
               ENDIF
!
               IF (.NOT.NOTCORNER) THEN
                  IJK_CORN(NCORN) = IJK
!
                  AXZ(IJK) = ZERO
                  AXZ(IJMK) = ZERO
                  AYZ(IJK) = ZERO
                  AYZ(IMJK) = ZERO
                  AXY(IJK) = ZERO
                  AXY(IJKM) = ZERO
!
                  AXZ_U(IM_OF(IJMK)) = ZERO
                  AXZ_U(IJMK) = ZERO
                  AXZ_U(IMJK) = ZERO
                  AXZ_U(IJK) = ZERO
!
                  AXY_U(IM_OF(IJKM)) = ZERO
                  AXY_U(IJKM) = ZERO
                  AXY_U(IMJK) = ZERO
                  AXY_U(IJK) = ZERO
!
                  AYZ_V(JM_OF(IMJK)) = ZERO
                  AYZ_V(IMJK) = ZERO
                  AYZ_V(IJMK) = ZERO
                  AYZ_V(IJK) = ZERO
!
                  AXY_V(JM_OF(IJKM)) = ZERO
                  AXY_V(IJKM) = ZERO
                  AXY_V(IJMK) = ZERO
                  AXY_V(IJK) = ZERO
!
                  AYZ_W(KM_OF(IMJK)) = ZERO
                  AYZ_W(IMJK) = ZERO
                  AYZ_W(IJKM) = ZERO
                  AYZ_W(IJK) = ZERO
!
                  AXZ_W(KM_OF(IJMK)) = ZERO
                  AXZ_W(IJMK) = ZERO
                  AXZ_W(IJKM) = ZERO
                  AXZ_W(IJK) = ZERO
               ENDIF
!
            ENDIF
!
         ENDIF
      END DO
      IF (NCORN > 0) THEN
            CALL START_LOG
            IF(DMP_LOG)WRITE (UNIT_LOG, 1000)
!
         DO L = 1, NCORN
            IJK = IJK_CORN(L)
            IF(DMP_LOG)WRITE (UNIT_LOG, 1100) IJK, I_OF(IJK), J_OF(IJK), K_OF(IJK)
         END DO
         IF(DMP_LOG)WRITE (UNIT_LOG, 1300)
         CALL END_LOG
      ENDIF
!
      RETURN
!
 1000 FORMAT(/1X,70('*')//' From: Get_Corner_Cells',/&
         ' Warning: The following wall-cells are adjacent to two or',/,&
         ' more fluid-cells.  Mass, momentum, and energy transfer ',/,&
         ' to these wall-cells have been set to zero.',/,&
         '     IJK     I     J     K')
 1100 FORMAT(3X,I6,2X,I4,2X,I4,2X,I4)
!
 1300 FORMAT(/1X,70('*')/)
      END SUBROUTINE GET_CORNER_CELLS
!
      SUBROUTINE ADDCORN(NOTCORNER, NCORN)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER NCORN
      LOGICAL NOTCORNER
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(LEN=80) :: LINE
!-----------------------------------------------
!
!                      error message
!
      NCORN = NCORN + 1
      IF (NCORN > MAX_NCORN) THEN
         WRITE (LINE, '(A)') 'Error: Increase MAX_NCORN in param1.inc.'
         CALL WRITE_ERROR ('AddCorn', LINE, 1)
         CALL MFIX_EXIT(myPE)
      ENDIF
!
      NOTCORNER = .FALSE.
!
      RETURN
      END SUBROUTINE ADDCORN

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3

