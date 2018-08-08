!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_DISTANCE (DIST_MIN,D,DIM_MAX2,DIST_SC,DIST_VEC)   C
!  Purpose: Calculate the scalar and vector distance for the passed    C
!           direction.                                                 C
!                                                                      C
!  Author: P. Nicoletti                               Date: 07-JUL-92  C
!  Reviewer:                                                           C
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
!  Local variables: I                                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_DISTANCE (DIST_MIN, D,DIM_MAX2,DIST_SC,DIST_VEC)
!
      IMPLICIT NONE
!
      REAL               DIST_SC(*) , DIST_VEC(*)
      DOUBLE PRECISION   DIST_MIN, D(*)
      INTEGER            DIM_MAX2 , I
!
      DIST_SC(1)  = DIST_MIN - 0.5 * D(1)
      DIST_VEC(1) = DIST_MIN
!
      DO I = 2,DIM_MAX2
         DIST_SC(I)  = DIST_VEC(I-1) + 0.5 * D(I)
         DIST_VEC(I) = DIST_VEC(I-1) + D(I)
      END DO
!
      RETURN
      END
