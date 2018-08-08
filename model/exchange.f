!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: EXCHANGE                                                C
!  Purpose: Calls routines to drive calculations of the interphase     C
!           mass, momentum, and energy exchange coefficients/terms     C
!           if directed to do so by the corresponding flags            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 25-APR-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE EXCHANGE(IER)

! Global Variables:
!-----------------------------------------------------------------------
! Flags for calculating drag coefficient.
      use coeff, only: DRAGCOEF
! Flags for calculating heat transfer coefficient.
      use coeff, only: HEAT_TR

      implicit none

! Dummy arguments:
!-----------------------------------------------------------------------
      INTEGER, INTENT(INOUT) :: IER ! Error index


! Calculate drag coefficients
      CALL CALC_DRAG (DRAGCOEF, IER)

! Calculate interphase heat transfer coefficients
      CALL CALC_GAMA (HEAT_TR)


      return
      END SUBROUTINE EXCHANGE
