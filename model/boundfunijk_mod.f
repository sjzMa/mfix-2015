  MODULE BOUNDFUNIJK

  CONTAINS

      INTEGER FUNCTION BOUND_FUNIJK(pLI, pLJ, pLK)

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
      USE functions
      IMPLICIT NONE
!
!                      Dummy indices for I, J, K
      INTEGER          pLI, pLJ, pLK

!-----------------------------------------------

      BOUND_FUNIJK  = FUNIJK ( MIN( IEND3, MAX (ISTART3, pLI) ),&
                               MIN( JEND3, MAX (JSTART3, pLJ) ),&
                               MIN( KEND3, MAX (KSTART3, pLK) ) )

      END FUNCTION BOUND_FUNIJK

  END MODULE BOUNDFUNIJK
