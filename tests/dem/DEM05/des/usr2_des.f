!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS2_DES                                               !
!                                                                      !
!  Purpose: This routine is called within the discrete phase time loop !
!  after the source terms are applied and the time step updated. The   !
!  The user may insert code in this routine or call user defined       !
!  subroutines.                                                        !
!                                                                      !
!  This routien is called from the time loop, but no indicies (fluid   !
!  cell or particle) are defined.                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR2_DES

      Use discretelement, only: DES_VEL_NEW
      Use discretelement, only: DES_POS_NEW
      Use discretelement, only: OMEGA_NEW

      IMPLICIT NONE

      INTEGER :: LL

! Move particles 63-93 below particles 32-62 to fake a wall.
      DO LL=63, 93
         DES_VEL_NEW(:,LL) = 0.0d0
         DES_POS_NEW(1,LL) = DES_POS_NEW(1,LL-31)
         DES_POS_NEW(3,LL) = DES_POS_NEW(3,LL-31)
         DES_POS_NEW(2,LL) = 0.0475d0
         OMEGA_NEW(:,LL) = 0.0d0
      ENDDO

      RETURN
      END SUBROUTINE USR2_DES
