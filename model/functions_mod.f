MODULE functions

  USE compar
  USE discretelement
  USE geometry
  USE indices
  USE param1

! Functions for generating IJK indices for indicated basis:
!---------------------------------------------------------------------//
! INTEGER :: FUNIJK, FUNIJK_0
! INTEGER :: FUNIJK_PROC
! INTEGER :: FUNIJK_GL
! INTEGER :: FUNIJK_IO


! Logical functions to determine whether index is on my PE's domain or
! indicated subset:
!---------------------------------------------------------------------//
! LOGICAL :: IS_ON_myPE_owns
! LOGICAL :: IS_ON_myPE_plus1layer
! LOGICAL :: IS_ON_myPE_plus2layers
! LOGICAL :: IS_ON_myPE_wobnd


! Functions for calculating indicated directional shift in given IJK
! index. This will return the ijk index of the computational cell
! corresponding to the indicated shift when that computational cell
! is NOT a wall cell. If the computational cell is a wall cell then
! this will return its own ijk index. For example, east_of will return
! IPJK when IPJK is a fluid or flow at cell. However, if IPJK is a
! wall cell east_of will return IJK.
!---------------------------------------------------------------------//
! INTEGER :: EAST_OF, EAST_OF_0
! INTEGER :: WEST_OF, WEST_OF_0
! INTEGER :: SOUTH_OF, SOUTH_OF_0
! INTEGER :: NORTH_OF, NORTH_OF_0
! INTEGER :: BOTTOM_OF, BOTTOM_OF_0
! INTEGER :: TOP_OF, TOP_OF_0


! Functions for calculating indicated directional shift in given IJK
! index. This will generally return the ijk index of the computational
! cell corresponding to the indicated shift regardless of the wall
! status of that computational cell. It may not return corner cells
! unless the ijk cell itself is a corner cell.
!---------------------------------------------------------------------//
! INTEGER :: IP_OF, IP_OF_0
! INTEGER :: IM_OF, IM_OF_0
! INTEGER :: JP_OF, JP_OF_0
! INTEGER :: JM_OF, JM_OF_0
! INTEGER :: KP_OF, KP_OF_0
! INTEGER :: KM_OF, KM_OF_0


! Logical functions to identify indicated condition:
!---------------------------------------------------------------------//
! logical function to identify various fluid/flow cells
! LOGICAL :: FLUID_AT
! LOGICAL :: P_FLOW_AT
! LOGICAL :: P_OUTFLOW_AT
! LOGICAL :: MASS_OUTFLOW_AT
! LOGICAL :: OUTFLOW_AT
! LOGICAL :: FLOW_AT
! LOGICAL :: FLUIDorP_FLOW_AT

! logical function to identify various wall cells
! LOGICAL :: WALL_AT
! LOGICAL :: NS_WALL_AT, FS_WALL_AT, PS_WALL_AT
! LOGICAL :: DEFAULT_WALL_AT
! LOGICAL :: WALL_ICBC_FLAG

! Logical function to identify a cyclic cell and different
! cyclic flow boundaries
! LOGICAL :: CYCLIC_AT
! LOGICAL :: CYCLIC_AT_E, CYCLIC_AT_N, CYCLIC_AT_T

! logical function to identify different flow at boundaries
! LOGICAL :: FLOW_AT_E, FLOW_AT_N, FLOW_AT_T
! LOGICAL :: PFLOW_AT_E, PFLOW_AT_N, PFLOW_AT_T
! LOGICAL :: MFLOW_AT_E, MFLOW_AT_N, MFLOW_AT_T

! Logical functions to identify different impermeable/semipermeable
! surface boundaries (specific type of internal surface)
! LOGICAL :: IP_AT_E, IP_AT_N, IP_AT_T
! LOGICAL :: SIP_AT_E, SIP_AT_N, SIP_AT_T
! LOGICAL :: SP_AT_E, SP_AT_N, SP_AT_T

! Logical functions concerning general internal surfaces
! LOGICAL :: IS_AT_E, IS_AT_N, IS_AT_T
! LOGICAL :: NO_IS_AT_E, NO_IS_AT_N, NO_IS_AT_T
! Integer function to return internal surface ID
! INTEGER :: IS_ID_AT_E, IS_ID_AT_N, IS_ID_AT_T


! Additional functions
!---------------------------------------------------------------------//
! DOUBLE PRECISION :: ZMAX
! INTEGER FUNCTION :: FUNLM


CONTAINS

  INCLUDE 'functions.inc'

END MODULE functions
