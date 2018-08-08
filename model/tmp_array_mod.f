!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: tmp_array                                                   !
!  Purpose:                                                            !
!     IMPORTANT:  For using these arrays in a subroutine               !
!     -lock the module in the beginning of the subroutine              !
!      call lock_tmp_array                                             !
!     -and unlock the module at the end of the subroutine              !
!      call unlock_tmp_array                                           !
! Contains the following subroutines:                                  !
!      lock_tmp_array, unlock_tmp_array                                !
!      lock_tmp_array2, unlock_tmp_array2                              !
!      lock_tmp4_array, unlock_tmp4_array                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      MODULE tmp_array
      IMPLICIT NONE

! temporary storage of dimension_3
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: &
         ARRAY1, ARRAY2, ARRAY3, ARRAY4

! temporary storage array of dimension_3, dimension_lm
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: &
         ArrayLM

! temporary storage for 4th order scheme of dimension_4
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TMP4

! temporary storage of dimension of dimension_3
      INTEGER, DIMENSION(:), ALLOCATABLE :: ARRAY1I

      character(LEN=3), dimension(:), pointer :: ARRAY1C

      LOGICAL :: tmp_array_locked = .false.
      LOGICAL :: tmp_array2_locked = .FALSE.
      LOGICAL :: tmp4_array_locked = .false.

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE lock_tmp_array
      USE compar, only: myPE
      USE funits, only: dmp_log
      IMPLICIT NONE
      IF(tmp_array_locked) THEN
         IF (DMP_LOG) WRITE(*,*) &
            'Error:  Multiple use of tmp_array (tmp_array_mod.f)'
         CALL MFIX_EXIT(myPE)
      ELSE
         tmp_array_locked = .true.
      ENDIF
      END SUBROUTINE lock_tmp_array

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE unlock_tmp_array
      IMPLICIT NONE
      tmp_array_locked = .false.
      END SUBROUTINE unlock_tmp_array

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE lock_tmp_array2
      USE compar, only: myPE
      USE funits, only: dmp_log
      IMPLICIT NONE
      IF(tmp_array2_locked) THEN
         IF(DMP_LOG) WRITE(*,*) &
            'Error:  Multiple use of tmp_array2 (tmp_array_mod.f)'
         CALL MFIX_EXIT(myPE)
      ELSE
         tmp_array2_locked = .true.
      ENDIF
      END SUBROUTINE lock_tmp_array2

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE unlock_tmp_array2
      IMPLICIT NONE
      tmp_array2_locked = .false.
      END SUBROUTINE unlock_tmp_array2

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE lock_tmp4_array
      USE compar, only: myPE
      USE funits, only: dmp_log
      IMPLICIT NONE
      IF(tmp4_array_locked) THEN
         IF(DMP_LOG) WRITE(*,*) &
            'Error:  Multiple use of tmp_array4 (tmp_array_mod.f)'
         CALL MFIX_EXIT(myPE)
      ELSE
         tmp4_array_locked = .true.
      ENDIF
      END SUBROUTINE lock_tmp4_array

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE unlock_tmp4_array
      IMPLICIT NONE
      tmp4_array_locked = .false.
      END SUBROUTINE unlock_tmp4_array


      END MODULE tmp_array



