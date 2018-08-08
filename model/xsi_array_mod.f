!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Use convection weighting factors for higher order          C
!  discretization.                                                     C
!                                                                      C
!  IMPORTANT:  For using these arrays in a subroutine                  C
!   -lock the module in the beginning of the subroutine                C
!    call lock_xsi_array                                               C
!   -and unlock the module at the end of the subroutine                C
!    call unlock_xsi_array                                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE xsi_array
      IMPLICIT NONE

! discretization factors xsi of dimension_3
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: XSI_e, XSI_n, &
                                                     XSI_t

      LOGICAL :: xsi_array_locked = .false.


      CONTAINS


      SUBROUTINE lock_xsi_array
      use compar, only: mype
      IMPLICIT NONE
      IF(xsi_array_locked) THEN
         Write(*,*)'Error:  Multiple use of xsi_array ',&
               '(xsi_array_mod.f)'
         CALL MFIX_EXIT(myPE)
      Else
         xsi_array_locked = .true.
      Endif
      END SUBROUTINE lock_xsi_array


      SUBROUTINE unlock_xsi_array
      IMPLICIT NONE
      xsi_array_locked = .false.
      END SUBROUTINE unlock_xsi_array

      END MODULE xsi_array
