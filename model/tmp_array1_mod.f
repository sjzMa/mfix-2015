!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: tmp_array1                                                  !
!  Purpose:                                                            !
!     IMPORTANT:  For using these arrays in a subroutine               !
!     -lock the module in the beginning of the subroutine              !
!      call lock_tmp_array1                                            !
!     -and unlock the module at the end of the subroutine              !
!      call unlock_tmp_array1                                          !
! Contains the following subroutines:                                  !
!      lock_tmp_array1, unlock_tmp_array1                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      MODULE tmp_array1

! temporary storage of dimension_3, dimension_m
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: &
         ARRAYm1

      LOGICAL :: tmp_array1_locked = .false.

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE lock_tmp_array1
      USE compar, only: myPE
      USE funits, only: dmp_log
      IMPLICIT NONE
      IF(tmp_array1_locked) THEN
         IF(DMP_LOG) WRITE(*,*) &
            'Error:  Multiple use of tmp_array1 (tmp_array1_mod.f)'
         CALL MFIX_EXIT(myPE)
      ELSE
         tmp_array1_locked = .true.
      ENDIF
      END SUBROUTINE lock_tmp_array1

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE unlock_tmp_array1
      IMPLICIT NONE
      tmp_array1_locked = .false.
      END SUBROUTINE unlock_tmp_array1

      END MODULE tmp_array1
