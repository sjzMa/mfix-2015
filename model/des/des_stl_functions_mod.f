!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: des_stl_functions                                      !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose: This module containd routines for geometric interaction    !
!  required for STL files.                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE des_stl_functions

      IMPLICIT NONE

! Use this module only to define functions and subroutines.
      CONTAINS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: TestTriangleAABB                                        !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE TestTriangleAABB(vertices, tri_norm, box_origin, &
         box_extents, sa_exist, sa, i,j,k)

      USE stl

      Implicit none

      !Separating axis test algorithm from Realtimecollision book by
      !Christer Ericson
      double precision, intent(in), dimension(3) :: tri_norm
      double precision, intent(in), dimension(3,3) :: vertices
      double precision, intent(in), dimension(3) :: box_origin, box_extents
      logical, intent(out) :: sa_exist
!      double precision, intent(out) :: proj_box, proj_tri(3)
      Integer, intent(out) :: sa
      Integer, intent(in) :: i,j,k

      double precision :: p0, p1, p2, r, e0, e1, e2

      double precision :: c_orig(3), c(3)
      double precision, dimension(3) :: vert0, vert1, vert2

      !box face normals
      double precision, dimension(3) :: u0, u1, u2

      double precision, dimension(3) :: v0, v1, v2
      double precision, dimension(3) :: f0, f1, f2

      !13 possible Separating axes
      double precision, dimension(9, 3) :: sep_axis
      !triangle plane n and d
      double precision :: tol_tri_aabb_proj

      double precision :: axis_magsq, s
      integer :: count
      sa_exist = .false.
      sa = 0
      tol_tri_aabb_proj = 1e-06

      vert0 = vertices(1,:)
      vert1 = vertices(2,:)
      vert2 = vertices(3,:)

      !Initially set e0, e1, and e2 to box extents

      !Compute box center and extents (if not already given in that format)
      c_orig(:) = box_origin(:) + box_extents(:)* 0.5d0

      !Translate triangle as conceptually moving AABB's center to origin

      v0 = vert0 - c_orig
      v1 = vert1 - c_orig
      v2 = vert2 - c_orig

      !Scale everything by box_extents
      v0(1) = v0(1)/box_extents(1)
      v1(1) = v1(1)/box_extents(1)
      v2(1) = v2(1)/box_extents(1)

      v0(2) = v0(2)/box_extents(2)
      v1(2) = v1(2)/box_extents(2)
      v2(2) = v2(2)/box_extents(2)

      v0(3) = v0(3)/box_extents(3)
      v1(3) = v1(3)/box_extents(3)
      v2(3) = v2(3)/box_extents(3)

      !set the box origin to (0,0,0)
      c = 0.d0
      !e0, e1, and e2 are half box extents which equals 0.5d0 in the scaled co-ordinates
      e0 = 0.5d0
      e1 = 0.5d0
      e2 = 0.5d0
      u0 = (/ 1.d0, 0.d0, 0.d0/)
      u1 = (/ 0.d0, 1.d0, 0.d0/)
      u2 = (/ 0.d0, 0.d0, 1.d0/)

!!$  write(*,'(A20,5(2x,g17.8))') 'box center orgi = ',c_orig(:)
!!$  write(*,'(A20,5(2x,g17.8))') 'box center new = ',c(:)
!!$  write(*,'(A20,5(2x,g17.8))') 'box extents = ',e0,e1,e2
!!$  write(*,'(A20,5(2x,g17.8))') 'tri vert1 = ',v0(:)
!!$  write(*,'(A20,5(2x,g17.8))') 'tri vert2 = ',v1(:)
!!$  write(*,'(A20,5(2x,g17.8))') 'tri vert3 = ',v2(:)

      !Compute edge vectors for triangle
      f0 = v1 - v0
      f1 = v2 - v1
      f2 = v0 - v2
      !Category 3 (edge-edge cross products)...
      !total of 9 (3 from triangle cross 3 independent from box)
      sep_axis(1, 1:3) = (/0.d0, -f0(3), f0(2)/)
      sep_axis(2, 1:3) = (/0.d0, -f1(3), f1(2)/)
      sep_axis(3, 1:3) = (/0.d0, -f2(3), f2(2)/)
      sep_axis(4, 1:3) = (/f0(3),  0.d0, -f0(1)/)
      sep_axis(5, 1:3) = (/f1(3),  0.d0, -f1(1)/)
      sep_axis(6, 1:3) = (/f2(3),  0.d0, -f2(1)/)
      sep_axis(7, 1:3) = (/-f0(2), f0(1), 0.d0/)
      sep_axis(8, 1:3) = (/-f1(2), f1(1), 0.d0/)
      sep_axis(9, 1:3) = (/-f2(2), f2(1), 0.d0/)

      do count = 1, 9
         axis_magsq = (sep_axis(count,1)**2+sep_axis(count,2)**2 &
         +sep_axis(count,3)**2)
         if(axis_magsq < TOL_STL*TOL_STL) THEN
            !nearly a zero vector. This will happen if the cross-products
            !from category 3 are formed from co-planar vectors
            !WRITE(*,'(/5x,A20, i2,/5x, A, 3(2x,g17.8))') 'Axis number:', count, &
            !     'Ignored due to zero vec:', sep_axis(count,:)
            CYCLE
         endif

         r = e0*abs(dot_product(u0(:), sep_axis(count,:)))             &
            + e1*abs(dot_product(u1(:), sep_axis(count,:)))            &
            + e2*abs(dot_product(u2(:), sep_axis(count,:)))

         p0 = dot_product(v0(:), sep_axis(count, :))
         p1 = dot_product(v1(:), sep_axis(count, :))
         p2 = dot_product(v2(:), sep_axis(count, :))

        IF(.FALSE.) THEN
            WRITE(*,'(A6,5(2X,G17.8))')'C:',C_ORIG
            WRITE(*,'(A6,5(2X,G17.8))')'VORIG:',VERT0
            WRITE(*,'(A6,5(2X,G17.8))')'V0:',V0
            WRITE(*,'(A6,5(2X,G17.8))')'V1:',V1
            WRITE(*,'(A6,5(2X,G17.8))')'EXTENTS:',E0,E1,E2
            WRITE(*,'(A6,5(2X,G17.8))')'SA:', COUNT, SEP_AXIS(COUNT,:)
            WRITE(*,'(A6,5(2X,G17.8))')'DATA:',  R, P0, P1, P2
            WRITE(*,'(A6,(2X,G17.8),"<",2(2X,G17.8), ">", 2X,G17.8)')  &
               'DATA2:',  MAX(P0, P1, P2)*1E8, -R*1E8, &
               MIN(P0, P1, P2)*1E8,R*1E8

            WRITE(*,'(A6,5(2X,L2))') 'DATA2:', MAX(P0, P1, P2)< -R,&
               MIN(P0, P1, P2) > R

            WRITE(*,'(A6,(2X,G17.8),"<",2(2X,G17.8), ">", 2X,G17.8)')&
               'DATA2:', (MAX(P0, P1, P2)+TOL_TRI_AABB_PROJ)*1E8,&
               -R*1E8, MIN(P0, P1, P2)*1E8, (R+TOL_TRI_AABB_PROJ)*1E8
            WRITE(*,'(A6,5(2X,L2))') 'DATA2:', MAX(P0, P1, P2)+&
               TOL_TRI_AABB_PROJ< -R, MIN(P0, P1, P2) > R+TOL_TRI_AABB_PROJ
         ENDIF

         IF (MAX(P0, P1, P2) +TOL_TRI_AABB_PROJ < -R .OR. &
            MIN(P0,P1,P2) > R+TOL_TRI_AABB_PROJ) THEN
            SA = COUNT
            SA_EXIST = .TRUE.
            RETURN
         ENDIF
      ENDDO

! Test the three axes corresponding to the face normals of AABB b (category 1).
! Could have been through dot products like above, but explotiing the
! the simple definition of box face normals to optimize the code.

! Exit if...

! ... [-e0, e0] and [min(v0.x,v1.x,v2.x), max(v0.x,v1.x,v2.x)] do not overlap
      IF (MAX(V0(1), V1(1), V2(1)) + TOL_TRI_AABB_PROJ< -E0 .OR.       &
         MIN(V0(1), V1(1), V2(1)) > E0+TOL_TRI_AABB_PROJ) THEN
         SA = 10
         SA_EXIST = .TRUE.
         RETURN
      ENDIF

!... [-e1, e1] and [min(v0.y,v1.y,v2.y), max(v0.y,v1.y,v2.y)] do not overlap
      IF (MAX(V0(2), V1(2), V2(2)) +TOL_TRI_AABB_PROJ< -E1 .OR.        &
         MIN(V0(2), V1(2), V2(2)) > E1+TOL_TRI_AABB_PROJ) THEN
         SA = 11
         SA_EXIST = .TRUE.
         RETURN
      ENDIF

! ... [-e2, e2] and [min(v0.z,v1.z,v2.z), max(v0.z,v1.z,v2.z)] do not overlap
      IF (MAX(V0(3), V1(3), V2(3)) +TOL_TRI_AABB_PROJ< -E2 .OR.        &
         MIN(V0(3), V1(3), V2(3)) > E2+TOL_TRI_AABB_PROJ) THEN
         SA = 12
         SA_EXIST = .TRUE.
         RETURN
      ENDIF

      R = E0*ABS(DOT_PRODUCT(U0(:), TRI_NORM(:)))                      &
         + E1*ABS(DOT_PRODUCT(U1(:), TRI_NORM(:)))                     &
         + E2*ABS(DOT_PRODUCT(U2(:), TRI_NORM(:)))

! Compute distance of box center from plane
! s = Dot_product(tri_norm(:), c(:) - v0(:))
! c = zero in the shfted co-ordinates
      s = Dot_product(tri_norm(:), v0(:))

!INTERSECTION OCCURS WHEN DISTANCE S FALLS WITHIN [-R,+R] INTERVAL
      IF(ABS(S) > R+ TOL_TRI_AABB_PROJ) THEN
         !PROJECTIONS DO NOT INTERSECT, SO THE SEPARATING AXIS EXISTS
         SA_EXIST = .TRUE.
         SA = 13
      ENDIF
      END SUBROUTINE TestTriangleAABB





!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: TestTriangleAABB                                        !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine ClosestPtPointTriangle(pointp, points, closest_point)
      USE param1, only: zero, one
      USE discretelement, only: dimn

      IMPLICIT NONE
      !point a, pointb, and pointc are the three nodes of the triangle
      !point p is the sphere center
      double precision, intent(in), dimension(3,3) :: points
      double precision, intent(in), dimension(dimn) :: pointp
      double precision, intent(out), dimension(dimn) ::  closest_point
      !Local variables
      double precision, dimension(3) :: pointa, pointb, pointc
      double precision, dimension(dimn) :: ab, ac, ap, bp,cp
      double precision :: d1, d2, d3, d4, vc, v, d5, d6, vb, w, va, denom

      pointa = points(1,:)
      pointb = points(2,:)
      pointc = points(3,:)

      ab = pointb - pointa
      ac = pointc - pointa
      ap = pointp - pointa
      d1 = DOT_PRODUCT(ab, ap)
      d2 = DOT_PRODUCT(ac, ap)

      IF(d1 <= Zero .AND. d2 <= zero) then
         closest_point = pointa ! barycentric coordinates (1,0,0)
         return
      end if

      !Check if P in vertex region outside B
      bp = pointp - pointb
      d3 = DOT_PRODUCT(ab, bp);
      d4 = DOT_PRODUCT(ac, bp);
      if (d3 >= zero .and. d4 <= d3) then
         closest_point = pointb !barycentric coordinates (0,1,0)
         return
      endif

      ! Check if P in edge region of AB, if so return projection of P onto AB
      vc = d1*d4 - d3*d2;
      if (vc <= zero .and. d1 >= zero .and. d3 <= zero) then
         v = d1 / (d1 - d3);
         closest_point =  pointa + v * ab; ! barycentric coordinates (1-v,v,0)
         return
      end if

      !Check if P in vertex region outside C
      cp = pointp - pointc
      d5 = DOT_PRODUCT(ab, cp)
      d6 = DOT_PRODUCT(ac, cp)
      if (d6 >= zero .and. d5 <= d6) then
         closest_point  = pointc ! barycentric coordinates (0,0,1)
         return
      endif

      !Check if P in edge region of AC, if so return projection of P onto AC
      vb = d5*d2 - d1*d6

      if (vb <= zero .and. d2 >= zero .and. d6 <= zero) then
         w = d2 / (d2 - d6)
         closest_point = pointa + w * ac ! barycentric coordinates (1-w,0,w)
         return
      end if

      !Check if P in edge region of BC, if so return projection of P onto BC
      va = d3*d6 - d5*d4
      if (va <= zero .and.(d4 - d3) >= zero .and. (d5 - d6) >= zero) then
         w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
         closest_point = pointb + w * (pointc - pointb) !barycentric coordinates (0,1-w,w)
         return
      end if


      !P inside face region. Compute Q through its barycentric coordinates (u,v,w)
      denom = one / (va + vb + vc)
      v = vb * denom
      w = vc * denom
      closest_point = pointa + ab * v + ac * w; ! = u*a + v*b + w*c, u = va * denom = 1.0f - v - w
      return
      end Subroutine ClosestPtPointTriangle

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: TestTriangleAABB                                        !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECKPTONTRIANGLE(POINTP, POINTS, ON_TRIAN)

! Global parameters:
!---------------------------------------------------------------------//
! DEM array spatial array dimension.
      use discretelement, only: DIMN
! Tolerance surround STL
      use stl, only: tol_stl
! Common values in working percision
      use param1, only : ONE


      IMPLICIT NONE


! Dummy arguments
!---------------------------------------------------------------------//
! the three nodes of the triangle
      DOUBLE PRECISION, INTENT(IN) :: POINTS(3,3)
! the test point
      DOUBLE PRECISION, INTENT(IN) :: POINTP(DIMN)
! Flag that the point lies on the triangle
      LOGICAL, INTENT(OUT)::  ON_TRIAN

! Local variables
!---------------------------------------------------------------------//
! triangle edges
      DOUBLE PRECISION, DIMENSION(DIMN) :: V0, V1, V2
      DOUBLE PRECISION :: D00, D01, D11, D20, D21, DENOM
! barcycentric coordinates
      DOUBLE PRECISION :: V, W


      V0 = POINTS(2,:) - POINTS(1,:)
      V1 = POINTS(3,:) - POINTS(1,:)
      V2 = POINTP - POINTS(1,:)

      D00 = dot_product(V0, V0)
      D01 = dot_product(V0, V1)
      D11 = dot_product(V1, V1)
      D20 = dot_product(V2, V0)
      D21 = dot_product(V2, V1)

      DENOM = D00*D11 - D01*D01

      V = (D11*D20 - D01*D21)/DENOM
      W = (D00*D21 - D01*D20)/DENOM

      ON_TRIAN = ((V>=-TOL_STL) .AND. (W>=-TOL_STL) .AND.              &
         ((V+W)<=ONE+TOL_STL))

      RETURN
      END SUBROUTINE CHECKPTONTRIANGLE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: intersectLnPlane                                        !
!  Author: Rahul Garg                                 Date: 24-Oct-13  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      Subroutine intersectLnPlane(ref_line, dir_line, ref_plane,       &
         norm_plane, line_param)

      USE discretelement, only: dimn
      USE param1, only: zero
      IMPLICIT NONE
      !reference point and direction of the line
      double precision, intent(in), dimension(dimn) :: ref_line,  dir_line
      !reference point and normal of the plane
      double precision, intent(in), dimension(dimn) :: ref_plane, norm_plane

      !line is parameterized as p = p_ref + t * dir_line, t is line_param
      double precision, intent(out) :: line_param

      !local vars
      double precision :: denom

      denom = DOT_PRODUCT(dir_line, norm_plane)
      if(denom*denom.gt.zero) then
         line_param = DOT_PRODUCT(ref_plane(:) - ref_line(:), norm_plane(:))
         line_param = line_param/denom
      endif
      return
      end subroutine intersectLnPlane

      subroutine set_facet_type_normal(facet_type)
        USE discretelement, only: Facet_type_normal
        Implicit none
        integer, intent(out) :: facet_type
        facet_type = FACET_type_NORMAL
      end subroutine set_facet_type_normal

      subroutine set_facet_type_po(facet_type)
        USE discretelement, only: Facet_type_po
        Implicit none
        integer, intent(out) :: facet_type
        facet_type = FACET_type_PO
      end subroutine set_facet_type_po

      subroutine set_facet_type_mi(facet_type)
        USE discretelement, only: Facet_type_mi
        Implicit none
        integer, intent(out) :: facet_type
        facet_type = FACET_type_MI
      end subroutine set_facet_type_mi


      END MODULE DES_STL_FUNCTIONS


