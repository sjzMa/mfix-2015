!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Zero_array(Array, IER)                                 C
!  Purpose: Zero out an array                                          C
!                                                                      C
!  Author: M. Syamlal                                 Date: 14-AUG-96  C
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
      SUBROUTINE ZERO_ARRAY(ARRAY)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!//11-02-99/Modified by Sreekanth to Zero-out the entire array passed in.
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Array
!// 504 1120 Modified the inherent dimensioning to absolute
!      double precision, intent(inout), dimension(:) :: ARRAY
      double precision, intent(inout), dimension(DIMENSION_3) :: ARRAY
!
!-----------------------------------------------
!
      ARRAY(:) = ZERO

      RETURN
      END SUBROUTINE ZERO_ARRAY
