      MODULE stl


      Use param
      Use param1

!     Maximum number of facets that can be read
      INTEGER, PARAMETER          :: DIM_STL = 10000000   !10 Million
!     Number of facets
      INTEGER                     :: N_FACETS
      INTEGER, PARAMETER :: WEST_FACEID = 9000000
      INTEGER, PARAMETER :: EAST_FACEID = 9000001
      INTEGER, PARAMETER :: SOUTH_FACEID = 9000002
      INTEGER, PARAMETER :: NORTH_FACEID = 9000003
      INTEGER, PARAMETER :: BOTTOM_FACEID = 9000004
      INTEGER, PARAMETER :: TOP_FACEID = 9000005

!     Number of facets for des. This could be a diiferent number from
! N_FACETS if the outer boundary is triangulated here
      INTEGER                     :: N_FACETS_DES
!     Vertex Coordinates X ,Y and Z
      DOUBLE PRECISION, DIMENSION(3,3,DIM_STL) :: VERTEX
!     Face normal vector (normalized)
      DOUBLE PRECISION, DIMENSION(3,DIM_STL) :: NORM_FACE
!     TRANSLATION COMPONENTS
      DOUBLE PRECISION :: TX_STL,TY_STL,TZ_STL
      DOUBLE PRECISION :: TX_MSH,TY_MSH,TZ_MSH
!     SCALING FACTOR
      DOUBLE PRECISION :: SCALE_STL
      DOUBLE PRECISION :: SCALE_MSH
!     RANGE OF STL FILE:
      DOUBLE PRECISION :: XMIN_STL,XMAX_STL
      DOUBLE PRECISION :: YMIN_STL,YMAX_STL
      DOUBLE PRECISION :: ZMIN_STL,ZMAX_STL
      DOUBLE PRECISION :: XMIN_MSH,XMAX_MSH
      DOUBLE PRECISION :: YMIN_MSH,YMAX_MSH
      DOUBLE PRECISION :: ZMIN_MSH,ZMAX_MSH
!     VALUE OF F_STL OUTSIDE OF STL BOUNDING BOX
      DOUBLE PRECISION :: OUT_STL_VALUE
      DOUBLE PRECISION :: OUT_MSH_VALUE
!     SMALLEST ANGLE FOR DETECTION OF SMALL TRIANGLES
      DOUBLE PRECISION :: STL_SMALL_ANGLE
      DOUBLE PRECISION :: MSH_SMALL_ANGLE
!     DIRECTION OF RAY TRACED TO DETERMINE WHETHER A POINT IS INSIDE/OUTSIDE
      CHARACTER(LEN=3) :: RAY_DIR
!     Tolerance for polygone edge detection
      DOUBLE PRECISION :: TOL_STL
      DOUBLE PRECISION :: TOL_MSH
      DOUBLE PRECISION :: TOL_STL_DP
!     Boundary condition ID
      INTEGER :: STL_BC_ID

      INTEGER, DIMENSION(DIM_STL) :: BC_ID_STL_FACE

!     Maximum number of facets per cell. The arrays below are used
! to define cut-cells under the CG modules
      INTEGER          :: DIM_FACETS_PER_CELL
      INTEGER, DIMENSION (:), ALLOCATABLE ::  N_FACET_AT
      INTEGER, DIMENSION (:,:), ALLOCATABLE ::  LIST_FACET_AT


!RG: Since Lagrangian requires facets that do no intersect at any edge of a cell,
!a separate facet list is maintained for Lagrangian modules, identfied by _DES
!appended to the key word
!     Maximum number of facets per cell
! Dynamic variable. for each ijk computational fluid cell store the
! total number of facets and the id's of the facets in that cell
      INTEGER :: MAX_FACETS_PER_CELL_DES
! in order to facilitate the parallel processing the PIC is defined
! as single array IJK
      TYPE FACETS_TO_CELL
         INTEGER :: COUNT_FACETS
         INTEGER, DIMENSION(:), ALLOCATABLE ::  FACET_LIST
      END TYPE FACETS_TO_CELL

      TYPE (FACETS_TO_CELL), DIMENSION (:), ALLOCATABLE ::  LIST_FACET_AT_DES
      CHARACTER(LEN=3) :: CAD_PROPAGATE_ORDER

      Logical, DIMENSION (:), ALLOCATABLE ::  NO_NEIGHBORING_FACET_DES

      END MODULE stl


