!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: PHYSICAL_PROP_ROg                                       !
!  Purpose: User hook for calculating the gas phase density.           !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PHYSICAL_PROP_ROg

! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! Function to identify wall cells
      use functions, only: WALL_AT

      use error_manager

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! Loop indicies
      INTEGER :: IJK   ! Computational cell
!......................................................................!


! The following error message is used to make sure that if a user
! defined gas density is invoked, that this routine has been modified.

!- REMOVE THE FOLLOWING ---------------------------------------------->>

      CALL INIT_ERR_MSG('USR_PHYSICAL_PROP_ROg')
      WRITE(ERR_MSG,9999)
      CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 9999 FORMAT('ERROR 9999: The user-defined drag routine was invoked ', &
         'but this',/'generic error message exits. Either choose a ',  &
         'different drag law',/'or correct mfix/model/usr_drag.f')

!- END REMOVE --------------------------------------------------------<<


      IJK_LP: DO IJK = IJKSTART3, IJKEND3
         IF(WALL_AT(IJK)) cycle IJK_LP

! Calculate the fluid density and bulk density
!         RO_G(IJK) = 
!         ROP_G(IJK) = 

      ENDDO IJK_LP

      RETURN
      END SUBROUTINE USR_PHYSICAL_PROP_ROg


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: PHYSICAL_PROP_ROs                                       !
!  Purpose: User hook for calculating solids phase density.            !
!                                                                      !
!  Author: J. Musser                                  Date: 28-JUN-13  !
!  Reviewer:                                          Date:            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PHYSICAL_PROP_ROs

! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! Function to identify wall cells
      use functions, only: WALL_AT
! Number of solids phases
      use physprop, only: MMAX

      use error_manager


      implicit none


! Local Variables:
!---------------------------------------------------------------------//
! Loop indicies
      INTEGER :: M, IJK
!......................................................................!

! The following error message is used to make sure that if a user
! defined gas density is invoked, that this routine has been modified.

!- REMOVE THE FOLLOWING ---------------------------------------------->>

      CALL INIT_ERR_MSG('USR_PHYSICAL_PROP_ROs')
      WRITE(ERR_MSG,9999)
      CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 9999 FORMAT('ERROR 9999: The user-defined drag routine was invoked ', &
         'but this',/'generic error message exits. Either choose a ',  &
         'different drag law',/'or correct mfix/model/usr_drag.f')

!- END REMOVE --------------------------------------------------------<<

      M_LP: DO M=1, MMAX
         IJK_LP: DO IJK = IJKSTART3, IJKEND3
            IF(WALL_AT(IJK)) cycle IJK_LP

! Calculate the solids density.
!            RO_S(IJK,M) = 

         ENDDO IJK_LP
      ENDDO M_LP

      RETURN
      END SUBROUTINE USR_PHYSICAL_PROP_ROs


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: PHYSICAL_PROP_CPg                                       !
!  Purpose: User hook for calculating the gas phase constant pressure  !
!  specific heat.                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PHYSICAL_PROP_CPg

! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! Function to identify wall cells
      use functions, only: WALL_AT

      use error_manager

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! Loop indicies
      INTEGER :: IJK
!......................................................................!

! The following error message is used to make sure that if a user
! defined gas density is invoked, that this routine has been modified.

!- REMOVE THE FOLLOWING ---------------------------------------------->>

      CALL INIT_ERR_MSG('USR_PHYSICAL_PROP_CPg')
      WRITE(ERR_MSG,9999)
      CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 9999 FORMAT('ERROR 9999: The user-defined drag routine was invoked ', &
         'but this',/'generic error message exits. Either choose a ',  &
         'different drag law',/'or correct mfix/model/usr_drag.f')

!- END REMOVE --------------------------------------------------------<<

      IJK_LP: DO IJK = IJKSTART3, IJKEND3
         IF(WALL_AT(IJK)) cycle IJK_LP

! Calculate the fluid density and bulk density
!         C_PG(IJK) = 

      ENDDO IJK_LP

      RETURN
      END SUBROUTINE USR_PHYSICAL_PROP_CPg


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: PHYSICAL_PROP_CPs                                       !
!  Purpose: User hook for calculating solids phase constant pressure   !
!  specific heat.                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_PHYSICAL_PROP_CPs

! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! Function to identify wall cells
      use functions, only: WALL_AT
! Number of solids phases
      use physprop, only: MMAX

      use error_manager

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! Loop indicies
      INTEGER :: M, IJK
!......................................................................!

! The following error message is used to make sure that if a user
! defined gas density is invoked, that this routine has been modified.

!- REMOVE THE FOLLOWING ---------------------------------------------->>

      CALL INIT_ERR_MSG('USR_PHYSICAL_PROP_CPs')
      WRITE(ERR_MSG,9999)
      CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 9999 FORMAT('ERROR 9999: The user-defined drag routine was invoked ', &
         'but this',/'generic error message exits. Either choose a ',  &
         'different drag law',/'or correct mfix/model/usr_drag.f')

!- END REMOVE --------------------------------------------------------<<

      M_LP: DO M=1, MMAX
         IJK_LP: DO IJK = IJKSTART3, IJKEND3
            IF(WALL_AT(IJK)) cycle IJK_LP

! Calculate the solids phase specific heat.
!            C_PS(IJK,M) = 

         ENDDO IJK_LP
      ENDDO M_LP

      RETURN
      END SUBROUTINE USR_PHYSICAL_PROP_CPs
