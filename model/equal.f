!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: EQUAL                                                   C
!  Purpose: Loop on the number of solids phases to set a variable      C
!           equal to the value or negative value of another variable   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 29-JAN-92  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE EQUAL(ARRAY1, IJK1, SIGN0, ARRAY2, IJK2)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE indices
      USE physprop, only : MMAX
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! First array
      DOUBLE PRECISION, INTENT(OUT) :: ARRAY1 (DIMENSION_3, *)
! Second array
      DOUBLE PRECISION, INTENT(IN) :: ARRAY2 (DIMENSION_3, *)
! IJK index for the first array
      INTEGER, INTENT(IN) :: IJK1
! IJK index for the second array
      INTEGER, INTENT(IN) :: IJK2
! Sign to be used when setting ARRAY1.  Legal values
! are + or - 1.0.
      DOUBLE PRECISION, INTENT(IN) :: SIGN0
!-----------------------------------------------
! Local variables
!-----------------------------------------------
!-----------------------------------------------

      IF (MMAX > 0) THEN
         ARRAY1(IJK1,:MMAX) = SIGN0*ARRAY2(IJK2,:MMAX)
      ENDIF

      RETURN
      END SUBROUTINE EQUAL
