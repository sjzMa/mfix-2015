!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: ADJUST_ROP                                              C
!  Purpose: Remove small negative values of density.                   C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 7-AUG-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE ADJUST_ROP(ROP)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE indices
      USE compar
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! density
      DOUBLE PRECISION, INTENT(INOUT) :: ROP(DIMENSION_3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK
!-----------------------------------------------

      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) ROP(IJK) = DMAX1(ZERO,ROP(IJK))
      ENDDO

      RETURN
      END SUBROUTINE ADJUST_ROP


