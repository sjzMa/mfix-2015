!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_RADIATION                                          !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 25-Jun-10  !
!                                                                      !
!  Commen:                                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_RADIATION(I, iM, iIJK, FOCUS)

      USE constant
      USE des_thermo
      USE discretelement
      USE fldvar
      USE param1
      USE physprop, only: SMAX
      USE toleranc

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
! Global index of particle
      INTEGER, INTENT(IN) :: I
! Solid phase of particle I
      INTEGER, INTENT(IN) :: iM
! Fluid cell index containing particle I
      INTEGER, INTENT(IN) :: iIJK
! Logical used for debugging
      LOGICAL, INTENT(IN) :: FOCUS

! Local variables
!---------------------------------------------------------------------//
! Surface area of particle
      DOUBLE PRECISION :: A_S
! Radiative heat transfer
      DOUBLE PRECISION :: Qrd
! Environment temperature
      DOUBLE PRECISION :: Tenv
! Particle Emmisivity
      DOUBLE PRECISION :: lEm

! Set the environment temperature.
      IF(COMPARE(EP_g(iIJK),ONE)) THEN
         Tenv = T_g(iIJK)
      ELSE
         Tenv = EP_g(iIJK)*T_g(iIJK) + (ONE-EP_g(iIJK))*avgDES_T_s(iIJK)
      ENDIF

! Set the particle emmisivity. Phase shift needed for TFM/DEM hybrid.
      lEM= DES_Em(iM + SMAX)

! Calculate the surface area of the particle
      A_S = 4.0d0 * Pi * DES_RADIUS(I) * DES_RADIUS(I)
! Calculate the heat source.
      Qrd = SB_CONST * A_s * lEm * (Tenv**4 - (DES_T_s_NEW(I))**4)
! Update the thermal source term.
      Q_Source(I) = Q_Source(I) + Qrd

      IF(FOCUS)THEN
         WRITE(*,"(//5X,A)")'From: DES_RADIATION -'
         WRITE(*,"(8X,A,D13.6)")'Tp: ',DES_T_s_NEW(I)
         WRITE(*,"(8X,A,D13.6)")'Tenv: ',Tenv
         WRITE(*,"(8X,A,D13.6)")'Qrd: ',Qrd
         WRITE(*,"(5X,25('-')/)")
      ENDIF

      RETURN

      END SUBROUTINE DES_RADIATION
