!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: des_stl_functions                                      !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose: This module containd routines for geometric interaction    !
!  required for STL files.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE STL_PREPROC_DES

      use des_stl_functions, only: TestTriangleAABB

      use desgrid
      use error_manager

      IMPLICIT NONE

! Use this module only to define functions and subroutines.
      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: des_stl_preprocessing                                   !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_STL_PREPROCESSING

      USE cutcell, only: use_stl

      USE discretelement, only: STL_FACET_TYPE
      USE discretelement, only: FACET_TYPE_NORMAL
      USE discretelement, only: FACET_TYPE_MI
      USE discretelement, only: FACET_TYPE_PO
      USE discretelement, only: COUNT_FACET_TYPE_NORMAL
      USE discretelement, only: COUNT_FACET_TYPE_MI
      USE discretelement, only: COUNT_FACET_TYPE_PO

      USE desgrid

      USE param
      USE stl

      implicit none

      INTEGER :: IJK, COUNT, NF

! Pre-procssing for the des in order to assign facets to grid cells.
      WRITE(ERR_MSG,"('Pre-Processing geometry for DES.')")
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      CALL ALLOCATE_DES_STL_ARRAYS

      N_FACETS_DES = 0

! Set N_facets_des to add any more facets needed by dem and not to.
! contaminate the Eulerian-Eulerian CG stuff
      IF(USE_STL) N_FACETS_DES = N_FACETS

! Set all the facets generated so far as normal facet types
      STL_FACET_TYPE(1:N_FACETS) = FACET_TYPE_NORMAL


      CALL BIN_FACETS_TO_GRID_DES

      DO IJK = 1, DG_IJKSIZE2
         COUNT = LIST_FACET_AT_DES(IJK)%COUNT_FACETS
         IF(COUNT.eq.0) DEALLOCATE(LIST_FACET_AT_DES(IJK)%FACET_LIST)
      ENDDO

      COUNT_FACET_TYPE_NORMAL = 0
      COUNT_FACET_TYPE_PO     = 0
      COUNT_FACET_TYPE_MI     = 0
      DO NF = 1, N_FACETS_DES
         IF(STL_FACET_TYPE(NF) == FACET_TYPE_NORMAL) &
            COUNT_FACET_TYPE_NORMAL = COUNT_FACET_TYPE_NORMAL + 1

         IF(STL_FACET_TYPE(NF).EQ.FACET_TYPE_MI) &
            COUNT_FACET_TYPE_MI = COUNT_FACET_TYPE_MI + 1

         IF(STL_FACET_TYPE(NF).EQ.FACET_TYPE_PO) &
            COUNT_FACET_TYPE_PO = COUNT_FACET_TYPE_PO + 1
      ENDDO

!      CALL DEBUG_WRITE_GRID_FACEINFO
!      CALL DEBUG_WRITE_STL_FROM_GRID_FACET(WRITE_FACETS_EACH_CELL=.FALSE.)

!      CALL DEBUG_WRITE_ALL_READIN_FACETS

      RETURN
      END SUBROUTINE DES_STL_PREPROCESSING

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ALLOCATE_DES_STL_ARRAYS                                 !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ALLOCATE_DES_STL_ARRAYS

      USE desgrid
      USE discretelement, only: STL_Facet_type
      USE discretelement, only: CELLNEIGHBOR_FACET
      USE discretelement, only: CELLNEIGHBOR_FACET_MAX
      USE discretelement, only: CELLNEIGHBOR_FACET_NUM
      USE discretelement, only: CELLNEIGHBOR_FACET_NUM
      USE stl, only: LIST_FACET_AT_DES
      USE stl, only: NO_NEIGHBORING_FACET_DES
      USE param
      USE stl, only: DIM_STL, MAX_FACETS_PER_CELL_DES

      IMPLICIT NONE

      integer :: ii,ijk

      ALLOCATE(LIST_FACET_AT_DES(DG_IJKSIZE2))

      Allocate(  CELLNEIGHBOR_FACET (DG_IJKSIZE2) )
      Allocate(  CELLNEIGHBOR_FACET_MAX (DG_IJKSIZE2) )
      Allocate(  CELLNEIGHBOR_FACET_NUM (DG_IJKSIZE2) )
      DO II = 1, DG_IJKSIZE2
         CELLNEIGHBOR_FACET_MAX(II) = 4
         allocate(CELLNEIGHBOR_FACET(II)%P(CELLNEIGHBOR_FACET_MAX(II)))
         allocate(CELLNEIGHBOR_FACET(II)%EXTENTDIR(CELLNEIGHBOR_FACET_MAX(II)))
         allocate(CELLNEIGHBOR_FACET(II)%EXTENTMIN(CELLNEIGHBOR_FACET_MAX(II)))
         allocate(CELLNEIGHBOR_FACET(II)%EXTENTMAX(CELLNEIGHBOR_FACET_MAX(II)))
         CELLNEIGHBOR_FACET_NUM(II) = 0
      ENDDO

      DO IJK = 1, DG_IJKSIZE2
         LIST_FACET_AT_DES(IJK)%COUNT_FACETS = 0
         ALLOCATE(LIST_FACET_AT_DES(IJK)%FACET_LIST(MAX_FACETS_PER_CELL_DES))
      ENDDO

      ALLOCATE(NO_NEIGHBORING_FACET_DES(DG_IJKSIZE2))
      NO_NEIGHBORING_FACET_DES = .false.

      Allocate(stl_facet_type(DIM_STL))

      END SUBROUTINE ALLOCATE_DES_STL_ARRAYS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: bin_facets_to_grid_des                                  !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine BIN_FACETS_TO_GRID_DES

      USE desgrid
      USE discretelement, only: dimn
      USE functions
      USE geometry
      USE indices
      USE mpi_utility
      USE param1, only: zero
      Use stl, only: tol_stl, n_facets_des
      Use stl, only: vertex, NO_NEIGHBORING_FACET_DES, LIST_FACET_AT_DES

      IMPLICIT NONE

! DES Grid cell index.
      INTEGER :: IJK, IJK2
! Loop counters:
      INTEGER :: I, I1, I2, II  ! X-axis
      INTEGER :: J, J1, J2, JJ  ! Y-axis
      INTEGER :: K, K1, K2, KK  ! Z-axis
      INTEGER :: N              ! STLs
! Generic accumulator
      INTEGER :: COUNT_FAC

! Maximum and minimum extents of the indexed STL
      DOUBLE PRECISION:: X1,Y1,Z1
      DOUBLE PRECISION:: X2,Y2,Z2



      DO N = 1,N_FACETS_DES

         X1 = minval(VERTEX(1:3,1,N))
         X2 = maxval(VERTEX(1:3,1,N))
         Y1 = minval(VERTEX(1:3,2,N))
         Y2 = maxval(VERTEX(1:3,2,N))
         Z1 = minval(VERTEX(1:3,3,N))
         Z2 = maxval(VERTEX(1:3,3,N))

         I1 = DG_IEND2
         I2 = DG_ISTART2
         IF(X2>=ZERO .AND. X1<=XLENGTH+TOL_STL) THEN
            I1 = iofpos(X1)
            I2 = iofpos(1/dg_dxinv+X2)
         ENDIF

         J1 = DG_JEND2
         J2 = DG_JSTART2
         IF(Y2>=ZERO .AND. Y1<=YLENGTH+TOL_STL) THEN
            J1 = jofpos(Y1)
            J2 = jofpos(1/dg_dyinv+Y2)
         ENDIF

         K1 = DG_KEND2
         K2 = DG_KSTART2
         IF(DO_K) THEN
            IF(Z2>=ZERO .AND. Z1<=ZLENGTH+TOL_STL) THEN
               K1 = kofpos(Z1)
               K2 = kofpos(1/dg_dzinv+Z2)
            ENDIF
         ENDIF


         DO K=K1,K2
         DO J=J1,J2
         DO I=I1,I2
            IF(dg_is_ON_myPE_plus1layers(I,J,K)) THEN
               IJK = DG_FUNIJK(I,J,K)
               CALL ADD_FACET_FOR_DES(I,J,K,IJK,N)
            ENDIF
         ENDDO
         ENDDO
         ENDDO

      ENDDO


! JM: I don't think this is needed.
      DO K = DG_KSTART1, DG_KEND1
      DO J = DG_JSTART1, DG_JEND1
      DO I = DG_ISTART1, DG_IEND1

         IJK = DG_FUNIJK(I, J, K)

         I1 =  MAX( I - 1, DG_ISTART2)
         I2 =  MIN( I + 1, DG_IEND2)

         J1 =  MAX( J - 1, DG_JSTART2)
         J2 =  MIN (J + 1, DG_JEND2)

         K1 =  MAX( K - 1, DG_KSTART2)
         K2 =  MIN (K + 1, DG_KEND2)

         NO_NEIGHBORING_FACET_DES(IJK)  = .FALSE.

         COUNT_FAC = 0
         DO KK = K1, k2
         DO JJ = J1, j2
         DO II = I1, i2
            IJK2 = DG_FUNIJK(II, JJ, KK)
            COUNT_FAC = COUNT_FAC + LIST_FACET_AT_DES(IJK2)%COUNT_FACETS
         ENDDO
         ENDDO
         ENDDO

         IF(COUNT_FAC == 0) NO_NEIGHBORING_FACET_DES(IJK) = .TRUE.

      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE BIN_FACETS_TO_GRID_DES

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ADD_FACET_FOR_DES                                       !
!  Author: Rahul Garg                                  Date: 24-Oct-13 !
!                                                                      !
!  Purpose: Add facets to DES grid cells.                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ADD_FACET_FOR_DES(I,J,K,IJK,N)

      USE param1, only: zero, one
      use geometry, only: DO_K

      use desgrid, only: dg_dxinv
      use desgrid, only: dg_dyinv
      use desgrid, only: dg_dzinv

      Use discretelement, only: MAX_RADIUS

      Use stl
      use error_manager

      IMPLICIT NONE

! DES grid index and facet index
      INTEGER, INTENT(IN) :: I,J,K,IJK, N

! Center of DES grid cell and half size. Note that a buffer is added to
! the half size to make the cell appear a little larger. This ensures
! that paricles near the edge 'see' STLs that are nearby but do not
! directly intersect the DES grid cell contain the particle center.
      DOUBLE PRECISION :: CENTER(3), HALFSIZE(3)
! Flag: STL intersects the DES grid cell
      LOGICAL :: OVERLAP
! DES grid cell dimensions
      DOUBLE PRECISION :: lDX, lDY, lDZ
! Buffer to ensure all particle-STL collisions are captured.
      DOUBLE PRECISION :: BUFFER
! Legacy variable - should be removed
      INTEGER ::  CURRENT_COUNT


      BUFFER = 1.1d0*MAX_RADIUS

      lDX = ONE/DG_DXINV
      lDY = ONE/DG_DYINV
      lDZ = ONE/DG_DZINV

      CENTER(1) = (dble(I-2)+HALF)*lDX
      HALFSIZE(1) = HALF*lDX + BUFFER

      CENTER(2) = (dble(J-2)+HALF)*lDY
      HALFSIZE(2) = HALF*lDY + BUFFER

      IF(DO_K)THEN
         CENTER(3) = (dble(K-2)+HALF)*lDZ
         HALFSIZE(3) = HALF*lDZ + BUFFER
      ELSE
         CENTER(3) = HALF*lDZ
         HALFSIZE(3) = HALF*lDZ
      ENDIF

      CALL TRI_BOX_OVERLAP(CENTER, HALFSIZE, VERTEX(:,:,N), OVERLAP)

      IF(OVERLAP) THEN

         CALL ADD_FACET(IJK, N)
         CURRENT_COUNT = LIST_FACET_AT_DES(IJK)%COUNT_FACETS

! These variables are not really needed. They are currently used
! by the PIC model and should be moved to the new forms used by
! the DEM wall collision routine.
         IF(CURRENT_COUNT .LT. MAX_FACETS_PER_CELL_DES) THEN
            LIST_FACET_AT_DES(IJK)%COUNT_FACETS = CURRENT_COUNT+1
            LIST_FACET_AT_DES(IJK)%FACET_LIST(CURRENT_COUNT+1) = N
         ELSE

         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE ADD_FACET_FOR_DES


!......................................................................!
! Subroutine TRI_BOX_OVERLAP                                           !
! Author: J.Musser                                   Date: 10-22-2015  !
!                                                                      !
! Purpose: Determine if a box (DES grid cell) intersects the triangle  !
!    (SLT). Note that the DES grid size is slightly increased to       !
!    capture STLs near the boarder of the cell. Otherwise, some        !
!    collisions could ve over looked.                                  !
!                                                                      !
! Author: Tomas Akenine-Moller                   Accessed: 10-22-2015  !
! REF: http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/    !
!         code/tribox2.txt                                             !
!......................................................................!
      SUBROUTINE TRI_BOX_OVERLAP(pCENTER, pHALFSIZE, pVERTS, pOVERLAP)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: pCENTER(3), pHALFSIZE(3)
      DOUBLE PRECISION, INTENT(IN) :: pVERTS(3,3)
      LOGICAL, INTENT(OUT) :: pOVERLAP

      DOUBLE PRECISION :: v0(3), v1(3), v2(3)
      DOUBLE PRECISION :: fex, fey, fez
      DOUBLE PRECISION :: normal(3), e0(3), e1(3), e2(3)

      pOVERLAP = .FALSE.

      v0 = pVERTS(1,:) - pCENTER
      v1 = pVERTS(2,:) - pCENTER
      v2 = pVERTS(3,:) - pCENTER

      e0 = v1-v0
      e1 = v2-v1
      e2 = v0-v2

      fex = abs(e0(1))
      fey = abs(e0(2))
      fez = abs(e0(3))

      if(ATEST_X01(e0(3),e0(2),fez,fey)) return
      if(ATEST_Y02(e0(3),e0(1),fez,fex)) return
      if(ATEST_Z12(e0(2),e0(1),fey,fex)) return

      fex = abs(e1(1))
      fey = abs(e1(2))
      fez = abs(e1(3))

      if(ATEST_X01(e1(3),e1(2),fez,fey)) return
      if(ATEST_Y02(e1(3),e1(1),fez,fex)) return
      if(ATEST_Z0 (e1(2),e1(1),fey,fex)) return

      fex = abs(e2(1))
      fey = abs(e2(2))
      fez = abs(e2(3))

      if(ATEST_X2 (e2(3),e2(2),fez,fey)) return
      if(ATEST_Y1 (e2(3),e2(1),fez,fex)) return
      if(ATEST_Z12(e2(2),e2(1),fey,fex)) return

      if(findMin(v0(1),v1(1),v2(1)) > phalfsize(1) .OR. &
         findMax(v0(1),v1(1),v2(1)) <-phalfsize(1)) return

      if(findMin(v0(2),v1(2),v2(2)) > phalfsize(2) .OR. &
         findMax(v0(2),v1(2),v2(2)) <-phalfsize(2)) return

      if(findMin(v0(3),v1(3),v2(3)) > phalfsize(3) .OR. &
         findMax(v0(3),v1(3),v2(3)) <-phalfsize(3)) return


      normal(1) = e0(2)*e1(3)-e0(3)*e1(2)
      normal(2) = e0(3)*e1(1)-e0(1)*e1(3)
      normal(3) = e0(1)*e1(2)-e0(2)*e1(1)

      if(.NOT.planeBoxOverlap(normal,v0,phalfsize)) return

      pOVERLAP = .TRUE.

      RETURN

      CONTAINS

!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!                                                                      !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION planeBoxOverlap(norm, vert, maxbox)

      double precision :: norm(3), vert(3), maxbox(3)

      integer :: lc
      double precision :: vmin(3), vmax(3), v

      do lc=1,3
         v=vert(lc)
         if(norm(lc) > 0.0d0) then
            vmin(lc) = -maxbox(lc) - v
            vmax(lc) =  maxbox(lc) - v
         else
            vmin(lc) = maxbox(lc) - v
            vmax(lc) =-maxbox(lc) - v
         endif
      enddo

      if(dot_product(norm,vmin) > 0.0d0) RETURN
      planeBoxOverlap=(dot_product(norm,vmax) >= 0.0d0)

      RETURN
      END FUNCTION planeBoxOverlap

!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!                                                                      !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      DOUBLE PRECISION FUNCTION findMin(x0,x1,x2)

      double precision :: x0,x1,x2

      findMin = x0

      if(x1<findMin) findMin=x1
      if(x2<findMin) findMin=x2

      RETURN
      END FUNCTION findMin

!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!                                                                      !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      DOUBLE PRECISION FUNCTION findMax(x0,x1,x2)

      double precision :: x0,x1,x2

      findMax = x0

      if(x1>findMax) findMax=x1
      if(x2>findMax) findMax=x2

      RETURN
      END FUNCTION findMax

!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!                                                                      !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION ATEST_X01(a,b,fa,fb)

      double precision :: a, b, fa, fb
      double precision :: lMin, lMax, p0, p2, rad

      p0 = a*v0(2) - b*v0(3)
      p2 = a*v2(2) - b*v2(3)

      if(p0<p2) then; lMIN=p0; lMAX=p2
      else; lMIN=p2; lMAX=p0; endif

      rad=fa*phalfsize(2) + fb*phalfsize(3)
      ATEST_X01=(lmin>rad .OR. lmax<-rad)

      END FUNCTION ATEST_X01

!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!                                                                      !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION ATEST_X2(a,b,fa,fb)

      double precision :: a, b, fa, fb
      double precision :: lMin, lMax, p0, p1, rad

      p0 = a*v0(2) - b*v0(3)
      p1 = a*v1(2) - b*v1(3)

      if(p0<p1) then; lMIN=p0; lMAX=p1
      else; lMIN=p1; lMAX=p0; endif

      rad=fa*phalfsize(2) + fb*phalfsize(3)
      ATEST_X2=(lmin>rad .OR. lmax<-rad)

      END FUNCTION ATEST_X2

!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!                                                                      !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION ATEST_Y02(a,b,fa,fb)

      double precision :: a, b, fa, fb
      double precision :: lMin, lMax, p0, p2, rad

      p0 = -a*v0(1) + b*v0(3)
      p2 = -a*v2(1) + b*v2(3)

      if(p0<p2) then; lMIN=p0; lMAX=p2
      else; lMIN=p2; lMAX=p0; endif

      rad=fa*phalfsize(1) + fb*phalfsize(3)
      ATEST_Y02=(lmin>rad .OR. lmax<-rad)

      END FUNCTION ATEST_Y02

!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!                                                                      !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION ATEST_Y1(a,b,fa,fb)

      double precision :: a, b, fa, fb
      double precision :: lMin, lMax, p0, p1, rad

      p0 = -a*v0(1) + b*v0(3)
      p1 = -a*v1(1) + b*v1(3)

      if(p0<p1) then; lMIN=p0; lMAX=p1
      else; lMIN=p1; lMAX=p0; endif

      rad=fa*phalfsize(1) + fb*phalfsize(3)
      ATEST_Y1=(lmin>rad .OR. lmax<-rad)

      END FUNCTION ATEST_Y1

!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!                                                                      !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION ATEST_Z12(a,b,fa,fb)

      double precision :: a, b, fa, fb
      double precision :: lMin, lMax, p1, p2, rad

      p1 = a*v1(1) - b*v1(2)
      p2 = a*v2(1) - b*v2(2)

      if(p2<p1) then; lMIN=p2; lMAX=p1
      else; lMIN=p1; lMAX=p2; endif

      rad=fa*phalfsize(1) + fb*phalfsize(2)
      ATEST_Z12=(lmin>rad .OR. lmax<-rad)

      END FUNCTION ATEST_Z12

!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!                                                                      !
!                                                                      !
!``````````````````````````````````````````````````````````````````````!
      LOGICAL FUNCTION ATEST_Z0(a,b,fa,fb)

      double precision :: a, b, fa, fb
      double precision :: lMin, lMax, p0, p1, rad

      p0 = a*v0(1) - b*v0(2)
      p1 = a*v1(1) - b*v1(2)

      if(p0<p1) then; lMIN=p0; lMAX=p1
      else; lMIN=p1; lMAX=p0; endif

      rad=fa*phalfsize(1) + fb*phalfsize(2)
      ATEST_Z0=(lmin>rad .OR. lmax<-rad)

      END FUNCTION ATEST_Z0

      END SUBROUTINE TRI_BOX_OVERLAP


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DEBUG_WRITE_GRID_FACEINFO                               !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DEBUG_WRITE_GRID_FACEINFO
      USE param
      USE param1
      USE run
      USE geometry
      USE indices
      USE compar
      USE stl
      USE functions

      IMPLICIT NONE

      INTEGER :: CELL_ID, I, J, K, COUNT, COUNT_FACETS

      CHARACTER(LEN=100) :: FN

      IF(NODESI*NODESJ*NODESK == 1) THEN
         WRITE(FN,'("FACETS_DG_GRID.DAT")')
      ELSE
         WRITE(FN,'("FACETS_DG_GRID_",I5.5,".DAT")') myPE
      ENDIF

      OPEN(1001, file = TRIM(FN), form ='formatted',CONVERT='BIG_ENDIAN')

      DO K=DG_KSTART2, DG_KEND2
      DO J=DG_JSTART2, DG_JEND2
      DO I=DG_ISTART2, DG_IEND2

         CELL_ID = DG_FUNIJK(I,J,K)
         COUNT_FACETS =  LIST_FACET_AT_DES(CELL_ID)%COUNT_FACETS
         IF(COUNT_FACETS.eq.0) cycle

         WRITE(1001,2000) CELL_ID, I, J, K,  COUNT_FACETS

 2000 FORMAT(50('*'),/2X,'CELL IJK, I, J, K =        = ', i20, 2x,     &
         4(2x,i10),/2X,'TOTAL FACETS',18(' '),'= ',3(2x, i10))

         DO COUNT = 1, COUNT_FACETS
            WRITE(1001, '(2x, i20)')                                   &
               LIST_FACET_AT_DES(CELL_ID)%FACET_LIST(COUNT)
         ENDDO

      ENDDO
      ENDDO
      ENDDO

      CLOSE(1001, STATUS = "keep")

      RETURN
      END SUBROUTINE DEBUG_WRITE_GRID_FACEINFO



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DEBUG_write_all_readin_facets                           !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine  DEBUG_write_all_readin_facets

      USE stl
      USE compar

      IMPLICIT NONE

      INTEGER :: N
!      CHARACTER(LEN=100) :: FN


      WRITE(ERR_MSG, 2000) N_FACETS, N_FACETS_DES
      CALL FLUSH_ERR_MSG(HEADER=.FALSE.,FOOTER=.FALSE.)

 2000 FORMAT('Saving STL geometry on DES grid: FACETS_TO_DG.stl',      &
         11x,'Facet Count:',I10,/2x,'DES Grid Facet Count:',I10)

      OPEN(UNIT=444, FILE='FACETS_TO_DG.stl',CONVERT='BIG_ENDIAN')
      write(444,*)'solid vcg'

      DO N = 1, N_FACETS_DES
         write(444,*) '   facet normal ', NORM_FACE(1:3,N)
         write(444,*) '      outer loop'
         write(444,*) '         vertex ', VERTEX(1,1:3,N)
         write(444,*) '         vertex ', VERTEX(2,1:3,N)
         write(444,*) '         vertex ', VERTEX(3,1:3,N)
         write(444,*) '      endloop'
         write(444,*) '   endfacet'
      ENDDO
      write(444,*)'endsolid vcg'
      close(444)



      END Subroutine DEBUG_write_all_readin_facets


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DEBUG_write_stl_from_grid_facet                         !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DEBUG_WRITE_STL_FROM_GRID_FACET(WRITE_FACETS_EACH_CELL)

      use run
      USE stl

      USE functions
      USE geometry
      USE indices
      USE compar
      USE discretelement, only: STL_FACET_TYPE, FACET_TYPE_NORMAL, &
           FACET_TYPE_MI, FACET_TYPE_PO, &
           COUNT_FACET_TYPE_PO

      IMPLICIT NONE

      LOGICAL, INTENT(IN),optional  :: WRITE_FACETS_EACH_CELL
      INTEGER ::  CELL_ID, N, I, J, K, COUNT, COUNT_FACETS, W_UNIT
      CHARACTER(LEN=200) :: FN, FN_PO
      LOGICAL :: write_each_cell
      LOGICAL, DIMENSION(:), allocatable :: FACET_WRITTEN

      ALLOCATE (FACET_WRITTEN(DIM_STL))

      write_each_cell = .false.
      if(present(WRITE_FACETS_EACH_CELL)) then
         write_each_cell = WRITE_FACETS_EACH_CELL
      endif

      FACET_WRITTEN = .false.

      IF(NODESI*NODESJ*NODESK == 1) THEN
         WRITE(FN,'("GEOM_DG.stl")')
         WRITE(FN_PO,'("GEOM_DG_PO.stl")')
      ELSE
         WRITE(FN,'("GEOM_DG_",I5.5,".stl")') MYPE
         WRITE(FN_PO,'("GEOM_DG_PO_",I5.5,".stl")')  MYPE
      ENDIF

      IF(COUNT_FACET_TYPE_PO.GE.1) THEN
         OPEN(UNIT=443, FILE=TRIM(FN_PO),CONVERT='BIG_ENDIAN')
         write(443,*)'solid vcg'
      endif

      OPEN(UNIT=444, FILE=trim(FN),CONVERT='BIG_ENDIAN')
      write(444,*)'solid vcg'

      DO K=DG_KSTART2, DG_KEND2
      DO J=DG_JSTART2, DG_JEND2
      DO I=DG_ISTART2, DG_IEND2

         CELL_ID = DG_FUNIJK(I,J,K)
         COUNT_FACETS =  LIST_FACET_AT_DES(CELL_ID)%COUNT_FACETS

         IF(COUNT_FACETS.EQ.0) CYCLE

         IF(WRITE_EACH_CELL) THEN

            WRITE(FN,'("GEOM_DG",3("_",I3.3),"_",I8.8,".stl")')        &
               I, J, K, CELL_ID


            WRITE(FN_PO,'("GEOM_DG_PO",3("_",I3.3),"_",I8.8,".stl")') &
               I,J,K,CELL_ID

            OPEN(UNIT=446, FILE=FN,CONVERT='BIG_ENDIAN')
            WRITE(446,*)'solid vcg'

            OPEN(UNIT=445, FILE=FN_PO,CONVERT='BIG_ENDIAN')
            WRITE(445,*)'solid vcg'
         ENDIF

         DO COUNT = 1, COUNT_FACETS
            N = LIST_FACET_AT_DES(CELL_ID)%FACET_LIST(COUNT)

            IF(WRITE_EACH_CELL) THEN
               IF(STL_FACET_TYPE(N).EQ.FACET_TYPE_NORMAL) THEN
                  W_UNIT  = 446
               ELSE
                  W_UNIT  = 445
               ENDIF

               write(w_unit,*) '   facet normal ', NORM_FACE(:,N)
               write(w_unit,*) '      outer loop'
               write(w_unit,*) '         vertex ', VERTEX(1,1:3,N)
               write(w_unit,*) '         vertex ', VERTEX(2,1:3,N)
               write(w_unit,*) '         vertex ', VERTEX(3,1:3,N)
               write(w_unit,*) '      endloop'
               write(w_unit,*) '   endfacet'

            ENDIF

            IF (FACET_WRITTEN(N)) CYCLE

            IF(STL_FACET_TYPE(N).EQ.FACET_TYPE_NORMAL) THEN
               W_UNIT  = 444
            ELSE
               W_UNIT  = 443
            ENDIF

            write(w_unit,*) '   facet normal ', NORM_FACE(:,N)
            write(w_unit,*) '      outer loop'
            write(w_unit,*) '         vertex ', VERTEX(1,1:3,N)
            write(w_unit,*) '         vertex ', VERTEX(2,1:3,N)
            write(w_unit,*) '         vertex ', VERTEX(3,1:3,N)
            write(w_unit,*) '      endloop'
            write(w_unit,*) '   endfacet'
            facet_written(N) = .true.
         ENDDO

         if(write_each_cell) then
            write(445,*)'endsolid vcg'
            close(445)

            write(446,*)'endsolid vcg'
            close(446)
         ENDIF

         ENDDO
      ENDDO
      ENDDO
      write(444,*)'endsolid vcg'
      close(444)

      if(count_facet_type_po.gt.0) then
         write(443,*)'endsolid vcg'
         close(443)
      ENDIF

      DEALLOCATE (FACET_WRITTEN)

!      IF(MyPE == PE_IO) THEN
!         WRITE(*,*) ' The file geometry_from_grid_facets.stl was sucessfully written.'
!         WRITE(*,*) ' This is the based on facets stored grid wise for DES modules'
!         WRITE(*,*) ' and is provided for convenience and debugging (it is not used).'
!      ENDIF
      END Subroutine  DEBUG_write_stl_from_grid_facet



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ADD_FACET                                               !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ADD_FACET(CELL_ID, FACET_ID)

      Use discretelement
      USE stl

      implicit none

      INTEGER, INTENT(IN) :: cell_id, facet_id

      INTEGER, DIMENSION(:), ALLOCATABLE :: int_tmp
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: real_tmp

      INTEGER :: lSIZE2, ii
      DOUBLE PRECISION :: smallest_extent, min_temp, max_temp

      IF(STL_FACET_TYPE(facet_id) /= FACET_TYPE_NORMAL) RETURN

      DO II = 1, CELLNEIGHBOR_FACET_NUM(CELL_ID)
         IF(FACET_ID .EQ. CELLNEIGHBOR_FACET(CELL_ID)%P(II)) RETURN
      ENDDO

      CELLNEIGHBOR_FACET_NUM(CELL_ID) = &
         CELLNEIGHBOR_FACET_NUM(CELL_ID) + 1

      NO_NEIGHBORING_FACET_DES(CELL_ID)  = .FALSE.

      IF(cellneighbor_facet_num(cell_id) > &
         cellneighbor_facet_max(cell_id)) THEN

         cellneighbor_facet_max(cell_id) = &
         2*cellneighbor_facet_max(cell_id)

         lSIZE2 = size(cellneighbor_facet(cell_id)%p)
         allocate(int_tmp(cellneighbor_facet_max(cell_id)))
         int_tmp(1:lSIZE2) = cellneighbor_facet(cell_id)%p(1:lSIZE2)
         call move_alloc(int_tmp,cellneighbor_facet(cell_id)%p)

         lSIZE2 = size(cellneighbor_facet(cell_id)%extentdir)
         allocate(int_tmp(cellneighbor_facet_max(cell_id)))
         int_tmp(1:lSIZE2) = &
            cellneighbor_facet(cell_id)%extentdir(1:lSIZE2)
         call move_alloc(int_tmp,cellneighbor_facet(cell_id)%extentdir)

         lSIZE2 = size(cellneighbor_facet(cell_id)%extentmin)
         allocate(real_tmp(cellneighbor_facet_max(cell_id)))
         real_tmp(1:lSIZE2) = &
            cellneighbor_facet(cell_id)%extentmin(1:lSIZE2)
         call move_alloc(real_tmp,cellneighbor_facet(cell_id)%extentmin)

         lSIZE2 = size(cellneighbor_facet(cell_id)%extentmax)
         allocate(real_tmp(cellneighbor_facet_max(cell_id)))
         real_tmp(1:lSIZE2) = &
            cellneighbor_facet(cell_id)%extentmax(1:lSIZE2)
         call move_alloc(real_tmp,cellneighbor_facet(cell_id)%extentmax)

      ENDIF

      CELLNEIGHBOR_FACET(CELL_ID)%&
         P(CELLNEIGHBOR_FACET_NUM(CELL_ID)) = FACET_ID
      SMALLEST_EXTENT = HUGE(0.0)

      DO II=1,3
         MIN_TEMP = MINVAL(VERTEX(:,II,FACET_ID))
         MAX_TEMP = MAXVAL(VERTEX(:,II,FACET_ID))

         IF(ABS(MAX_TEMP - MIN_TEMP) < SMALLEST_EXTENT ) THEN
             CELLNEIGHBOR_FACET(CELL_ID)%&
                EXTENTDIR(CELLNEIGHBOR_FACET_NUM(CELL_ID)) = II
             CELLNEIGHBOR_FACET(CELL_ID)%&
                EXTENTMIN(CELLNEIGHBOR_FACET_NUM(CELL_ID)) = MIN_TEMP
             CELLNEIGHBOR_FACET(CELL_ID)%&
                EXTENTMAX(CELLNEIGHBOR_FACET_NUM(CELL_ID)) = MAX_TEMP
             SMALLEST_EXTENT = ABS(MAX_TEMP - MIN_TEMP)
         ENDIF
      ENDDO

      END SUBROUTINE ADD_FACET

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_IF_PARTICLE_OVERLAPS_STL                          !
!  Authors: Rahul Garg                               Date: 21-Mar-2014 !
!                                                                      !
!  Purpose: This subroutine is special written to check if a particle  !
!          overlaps any of the STL faces. The routine exits on         !
!          detecting an overlap. It is called after initial            !
!          generation of lattice configuration to remove out of domain !
!          particles                                                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_IF_PARTICLE_OVERLAPS_STL(POS, fI, fJ, fK, REMOVE)

      USE run
      USE param1
      USE discretelement, only: dimn, xe, yn, zt
      USE geometry
      USE constant
      USE cutcell
      USE indices
      USE stl
      USE compar
      USE des_stl_functions
      use desgrid
      USE functions

      Implicit none

      DOUBLE PRECISION, INTENT(IN) :: POS(DIMN)
      INTEGER, INTENT(IN) :: fI, fJ, fK
      LOGICAL, INTENT(OUT) :: REMOVE

! Integers mapping the fluid cell corners to DES Grid indices.
      INTEGER :: I1, I2, J1, J2, K1, K2

      INTEGER I, J, K, IJK, NF, LC

      DOUBLE PRECISION :: LINE_T
      DOUBLE PRECISION :: RADSQ, DIST(3)

      REMOVE = .TRUE.

      I1 = IofPOS(XE(fI-1))
      I2 = IofPOS(XE( fI ))

      J1 = JofPOS(YN(fJ-1))
      J2 = JofPOS(YN( fJ ))

      K1 = KofPOS(ZT(fK-1))
      K2 = KofPOS(ZT( fK ))

      RADSQ = (1.05d0*MAX_RADIUS)**2

      DO K = K1, K2
      DO J = J1, J2
      DO I = I1, I2

         IF(.NOT.DG_is_ON_myPE_plus1layers(I,J,K)) CYCLE

         IJK = DG_FUNIJK(I,J,K)

! The point is on the non-fluid side of the plane if t>0
         DO LC = 1, LIST_FACET_AT_DES(IJK)%COUNT_FACETS
            NF = LIST_FACET_AT_DES(IJK)%FACET_LIST(LC)

            CALL INTERSECTLNPLANE(POS, NORM_FACE(:,NF), &
               VERTEX(1,:,NF), NORM_FACE(:,NF), LINE_T)

! Orthogonal projection puts the point outside of the domain or less than
! one particle radius to the facet.
            DIST = LINE_T*NORM_FACE(:,NF)
            IF(LINE_T > ZERO .OR. dot_product(DIST,DIST)<=RADSQ) RETURN

         ENDDO
      ENDDO
      ENDDO
      ENDDO

      REMOVE = .FALSE.

      RETURN
      END SUBROUTINE CHECK_IF_PARTICLE_OVERLAPS_STL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE write_this_facet_and_part(FID, PID)
      USE run
      USE param1
      USE discretelement
      USE geometry
      USE compar
      USE constant
      USE cutcell
      USE funits
      USE indices
      USE physprop
      USE parallel
      USE stl
      USE des_stl_functions
      Implicit none
      !facet id and particle id
      Integer, intent(in) :: fid, pid
      Integer :: stl_unit, vtp_unit , k
      CHARACTER(LEN=100) :: stl_fname, vtp_fname
      real :: temp_array(3)

      stl_unit = 1001
      vtp_unit = 1002

      WRITE(vtp_fname,'(A,"_OFFENDING_PARTICLE",".vtp")') TRIM(RUN_NAME)
      WRITE(stl_fname,'(A,"_STL_FACE",".stl")') TRIM(RUN_NAME)

      open(vtp_unit, file = vtp_fname, form='formatted',convert='big_endian')
      open(stl_unit, file = stl_fname, form='formatted',convert='big_endian')

      write(vtp_unit,"(a)") '<?xml version="1.0"?>'
      write(vtp_unit,"(a,es24.16,a)") '<!-- time =',s_time,'s -->'
      write(vtp_unit,"(a,a)") '<VTKFile type="PolyData"',&
           ' version="0.1" byte_order="LittleEndian">'
      write(vtp_unit,"(3x,a)") '<PolyData>'
      write(vtp_unit,"(6x,a,i10.10,a,a)")&
           '<Piece NumberOfPoints="',1,'" NumberOfVerts="0" ',&
           'NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">'
      write(vtp_unit,"(9x,a)")&
           '<PointData Scalars="Diameter" Vectors="Velocity">'
      write(vtp_unit,"(12x,a)")&
           '<DataArray type="Float32" Name="Diameter" format="ascii">'
      write (vtp_unit,"(15x,es13.6)") (2*des_radius(pid))
      write(vtp_unit,"(12x,a)") '</DataArray>'

      temp_array = zero
      temp_array(:) = des_vel_new(:,pid)
      write(vtp_unit,"(12x,a,a)") '<DataArray type="Float32" ',&
           'Name="Velocity" NumberOfComponents="3" format="ascii">'
      write (vtp_unit,"(15x,3(es13.6,3x))")&
           ((temp_array(k)),k=1,3)
      write(vtp_unit,"(12x,a,/9x,a)") '</DataArray>','</PointData>'
      ! skip cell data
      write(vtp_unit,"(9x,a)") '<CellData></CellData>'

      temp_array = zero
      temp_array(1:dimn) = des_pos_new(1:dimn, pid)
      write(vtp_unit,"(9x,a)") '<Points>'
      write(vtp_unit,"(12x,a,a)") '<DataArray type="Float32" ',&
           'Name="Position" NumberOfComponents="3" format="ascii">'
      write (vtp_unit,"(15x,3(es13.6,3x))")&
           ((temp_array(k)),k=1,3)
      write(vtp_unit,"(12x,a,/9x,a)")'</DataArray>','</Points>'
      ! Write tags for data not included (vtp format style)
      write(vtp_unit,"(9x,a,/9x,a,/9x,a,/9x,a)")'<Verts></Verts>',&
           '<Lines></Lines>','<Strips></Strips>','<Polys></Polys>'
      write(vtp_unit,"(6x,a,/3x,a,/a)")&
           '</Piece>','</PolyData>','</VTKFile>'

      !Now write the facet info

      write(stl_unit,*)'solid vcg'

      write(stl_unit,*) '   facet normal ', NORM_FACE(1:3,FID)
      write(stl_unit,*) '      outer loop'
      write(stl_unit,*) '         vertex ', VERTEX(1,1:3,FID)
      write(stl_unit,*) '         vertex ', VERTEX(2,1:3,FID)
      write(stl_unit,*) '         vertex ', VERTEX(3,1:3,FID)
      write(stl_unit,*) '      endloop'
      write(stl_unit,*) '   endfacet'

      write(stl_unit,*)'endsolid vcg'

      close(vtp_unit, status = 'keep')
      close(stl_unit, status = 'keep')

      end SUBROUTINE write_this_facet_and_part








!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE write_stls_this_dg(dg)

      Use usr

! STL Vertices
      use stl, only: VERTEX
! STL Facet normals
      use stl, only: NORM_FACE

      use discretelement, only: CELLNEIGHBOR_FACET
      use discretelement, only: CELLNEIGHBOR_FACET_NUM

      IMPLICIT NONE
!-----------------------------------------------
      integer, intent(in) :: dg
      integer :: lc, nf

      logical :: EXISTS
      character(len=6) :: IDX

      write(idx,"(I6.6)") dg
      open(unit=555, file='dg_'//idx//'.stl', status='UNKNOWN')

      write(555,*) 'solid vcg'
      DO lc=1, cellneighbor_facet_num(dg)

         NF = cellneighbor_facet(dg)%p(lc)

         write(555,*) '   facet normal ', NORM_FACE(:,NF)
         write(555,*) '      outer loop'
         write(555,*) '         vertex ', VERTEX(1,1:3,NF)
         write(555,*) '         vertex ', VERTEX(2,1:3,NF)
         write(555,*) '         vertex ', VERTEX(3,1:3,NF)
         write(555,*) '      endloop'
         write(555,*) '   endfacet'
      enddo
      close(555)

      RETURN
      END SUBROUTINE write_stls_this_dg


!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE write_this_stl(lc)


! STL Vertices
      use stl, only: VERTEX
! STL Facet normals
      use stl, only: NORM_FACE
      use compar, only: myPE

      IMPLICIT NONE
!-----------------------------------------------
      integer, intent(in) :: lc

      logical :: EXISTS
      character(len=4) :: IDX
      character(len=4) :: IPE


      write(idx,"(I4.4)") LC
      write(ipe,"(I4.4)") myPE
      open(unit=555, file='geo_'//idx//'_'//IPE//'.stl',&
         status='UNKNOWN')
      write(555,*) 'solid vcg'
      write(555,*) '   facet normal ', NORM_FACE(:,LC)
      write(555,*) '      outer loop'
      write(555,*) '         vertex ', VERTEX(1,1:3,LC)
      write(555,*) '         vertex ', VERTEX(2,1:3,LC)
      write(555,*) '         vertex ', VERTEX(3,1:3,LC)
      write(555,*) '      endloop'
      write(555,*) '   endfacet'
      close(555)


      RETURN
      END SUBROUTINE write_this_stl

      END MODULE STL_PREPROC_DES


