!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS3_DES                                               !
!                                                                      !
!  Purpose: This routine is called after the discrete phase time loop  !
!  and is user-definable. The user may insert code in this routine or  !
!  call appropriate user defined subroutines.                          !
!                                                                      !
!  This routien is not called from a loop, hence all indicies are      !
!  undefined.                                                          !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR3_DES

      use rxns, only: R_gp
      use rxns, only: RoX_gc
      use rxns, only: R_PHASE
      use rxns, only: SUM_R_G
      use energy, only: HOR_g

      use param1, only: ZERO

      IMPLICIT NONE

      R_gp(:,:) = ZERO
      RoX_gc(:,:) = ZERO
      R_PHASE(:,:) = ZERO
      HOR_G(:) = ZERO
      SUM_R_G(:) = ZERO

      RETURN
      END SUBROUTINE USR3_DES
