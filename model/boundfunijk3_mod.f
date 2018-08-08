  MODULE BOUNDFUNIJK3

  CONTAINS

      INTEGER FUNCTION BOUND_FUNIJK3(pLI, pLJ, pLK)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE physprop
      USE geometry
      USE compar
      USE fldvar
      USE indices
      USE function3
      IMPLICIT NONE
!
!                      Dummy indices for I, J, K
      INTEGER          pLI, pLJ, pLK

!-----------------------------------------------

      BOUND_FUNIJK3  = FUNIJK3 ( MIN( IEND4, MAX (ISTART4, pLI) ),&
                               MIN( JEND4, MAX (JSTART4, pLJ) ),&
                               MIN( KEND4, MAX (KSTART4, pLK) ) )

      END FUNCTION BOUND_FUNIJK3

  END MODULE BOUNDFUNIJK3
