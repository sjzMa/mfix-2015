!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR0                                                   !
!  Author: J.Musser                                   Date: dd-mmm-yy  !
!  Purpose: This routine is called before the time loop starts and is  !
!           user-definable.  The user may insert code in this routine  !
!           or call appropriate user defined subroutines.  This        !
!           can be used for setting constants and checking errors in   !
!           data.  This routine is not called from an IJK loop, hence  !
!           all indices are undefined.                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE USR0

      use usr
      use discretelement
      use constant

      IMPLICIT NONE

      if(particles /= 2) then
         write(*,"(3x, 'invalid setup for test case')")
         call mfix_exit(0)
      endif

      gY1 = des_pos_new(2,1)
      gY2 = des_pos_new(2,2)
      gX1 = des_vel_new(2,1)
      gX1 = des_vel_new(2,2)

! Body forces. (Gravity)
      F1b = -GRAVITY
      F2b = -GRAVITY

      return
      END SUBROUTINE USR0
