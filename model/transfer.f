!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: TRANSFER(B_m, Var, M, IER)                             C
!  Purpose: Transfer the result of linear equation solver to variable  C
!           array.                                                     C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 28-MAY-96  C
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
      SUBROUTINE TRANSFER(B_M, VAR, M)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      phase index
      INTEGER          M
!
!                      index
      INTEGER          IJK
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!
!                      Variable
      DOUBLE PRECISION Var(DIMENSION_3)
!
!-----------------------------------------------
!
      IJK = 1
      IF (IJKMAX2 > 0) THEN
         VAR(:IJKMAX2) = B_M(:IJKMAX2,M)
         IJK = IJKMAX2 + 1
      ENDIF
      RETURN
      END SUBROUTINE TRANSFER
