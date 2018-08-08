!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_CONSTANTS                                           C
!  Purpose: Set various constants                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_CONSTANTS

! Modules
!---------------------------------------------------------------------//
      USE param1, only: zero, one, undefined
      USE constant, only: gas_const
      USE constant, only: gravity, gravity_x, gravity_y, gravity_z
      USE constant, only: to_SI
      USE constant, only: k_scale
      USE run, only: LAM_HYS, UNITS
      USE error_manager, only: err_msg, init_err_msg, finl_err_msg
      USE error_manager, only: flush_err_msg

      IMPLICIT NONE
!---------------------------------------------------------------------//

! Note that the cell flags are not set when this routine is called.

! Dimensionless constants
      K_SCALE = .08D0   ! this actually isn't used anywhere...

! Enter the value of all constants in various units (CGS or SI)
      IF (UNITS == 'SI') THEN
         IF (GRAVITY == UNDEFINED) GRAVITY = 9.80665D0 ! m/s2
         GAS_CONST = 8314.56D0                !Pa.m3/kmol.K, or kg m2/s2 kmol K (Perry and Green, 1984)
         to_SI = 0.1D0                        !to convert dyne/cm2 to Pa. e.g. calc_mu_g.f
         IF (LAM_HYS == UNDEFINED) LAM_HYS = 0.000001d0    ! m
      ELSEIF (UNITS == 'CGS') THEN
         IF (GRAVITY == UNDEFINED) GRAVITY = 980.665D0 !cm/s2
         GAS_CONST = 8314.56D4                !g.cm2/s2.mol.K
         to_SI = ONE                          !does not do anything in CGS. e.g.: calc_mu_g.f
         IF (LAM_HYS == UNDEFINED) LAM_HYS = 0.0001d0    ! cm
      ELSE
! Initialize the error manager.
         CALL INIT_ERR_MSG("SET_CONSTANTS")
         WRITE(ERR_MSG,1005) UNITS
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         CALL FINL_ERR_MSG
 1005 FORMAT('Error 1005: Unknown UNITS = ',A,/ &
         'Please correct the mfix.dat file.')
      ENDIF

! If all components of the gravitational acceleration vector are
! undefined (zero), then use the default value for the negative
! y-direction. This ensures backward compatibility with the old
! (legacy) GRAVITY keyword. At this point GRAVITY is defined,
! either from mfix.dat or by default above
      IF(GRAVITY_X==ZERO.AND.GRAVITY_Y==ZERO.AND.GRAVITY_Z==ZERO) THEN
         GRAVITY_X = ZERO
         GRAVITY_Y = - GRAVITY
         GRAVITY_Z = ZERO
      ENDIF

      RETURN
      END SUBROUTINE SET_CONSTANTS
