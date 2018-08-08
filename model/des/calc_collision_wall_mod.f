!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_COLLISION_WALL                                    C
!  Author: Rahul Garg                               Date: 1-Dec-2013   C
!                                                                      C
!  Purpose: subroutines for particle-wall collisions when cutcell is   C
!           used. Also contains rehack of routines for cfslide and     C
!           cffctow which might be different from the stand alone      C
!           routines. Eventually all the DEM routines will be          C
!           consolidated.                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      MODULE CALC_COLLISION_WALL

      PRIVATE
      PUBLIC :: CALC_DEM_FORCE_WITH_WALL_STL

      CONTAINS

!VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV!
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_DEM_FORCE_WITH_WALL_STL

      USE run
      USE param1
      USE desgrid
      USE discretelement
      USE geometry
      USE compar
      USE constant
      USE indices
      USE stl
      USE des_stl_functions
      USE functions
      Implicit none

      INTEGER :: LL
      INTEGER IJK, NF
      DOUBLE PRECISION OVERLAP_N, SQRT_OVERLAP

      DOUBLE PRECISION V_REL_TRANS_NORM, DISTSQ, RADSQ, CLOSEST_PT(DIMN)
! local normal and tangential forces
      DOUBLE PRECISION NORMAL(DIMN), VREL_T(DIMN), DIST(DIMN), DISTMOD
      DOUBLE PRECISION, DIMENSION(DIMN) :: FTAN, FNORM, OVERLAP_T

      LOGICAL :: DES_LOC_DEBUG
      INTEGER :: CELL_ID, cell_count
      INTEGER :: PHASELL

      DOUBLE PRECISION :: TANGENT(DIMN)
      DOUBLE PRECISION :: FTMD, FNMD
! local values used spring constants and damping coefficients
      DOUBLE PRECISION ETAN_DES_W, ETAT_DES_W, KN_DES_W, KT_DES_W

      double precision :: MAG_OVERLAP_T

      double precision :: line_t
! flag to tell if the orthogonal projection of sphere center to
! extended plane detects an overlap

      DOUBLE PRECISION :: MAX_DISTSQ, DISTAPART, FORCE_COH, R_LM
      INTEGER :: MAX_NF, axis
      DOUBLE PRECISION, DIMENSION(3) :: PARTICLE_MIN, PARTICLE_MAX

      DES_LOC_DEBUG = .false. ;      DEBUG_DES = .false.
      FOCUS_PARTICLE = -1


!$omp parallel default(none) private(LL,ijk,MAG_OVERLAP_T,             &
!$omp    cell_id,radsq,particle_max,particle_min,tangent,              &
!$omp    axis,nf,closest_pt,dist,r_lm,distapart,force_coh,distsq,      &
!$omp    line_t,max_distsq,max_nf,normal,distmod,overlap_n,VREL_T,     &
!$omp    v_rel_trans_norm,phaseLL,sqrt_overlap,kn_des_w,kt_des_w,      &
!$omp    etan_des_w,etat_des_w,fnorm,overlap_t,ftan,ftmd,fnmd)         &
!$omp shared(max_pip,focus_particle,debug_des,no_neighboring_facet_des,&
!$omp    pijk,dg_pijk,list_facet_at_des,i_of,j_of,k_of,des_pos_new,    &
!$omp    des_radius,cellneighbor_facet_num,cellneighbor_facet,vertex,  &
!$omp    hert_kwn,hert_kwt,kn_w,kt_w,des_coll_model_enum,mew_w,tow,    &
!$omp    des_etan_wall,des_etat_wall,dtsolid,fc,norm_face,             &
!$omp    wall_collision_facet_id,wall_collision_PFT,use_cohesion,      &
!$omp    van_der_waals,wall_hamaker_constant,wall_vdw_outer_cutoff,    &
!$omp    wall_vdw_inner_cutoff,asperities,surface_energy)
!$omp do
      DO LL = 1, MAX_PIP

         IF(LL.EQ.FOCUS_PARTICLE) DEBUG_DES = .TRUE.

! skipping non-existent particles or ghost particles
! make sure the particle is not classified as a new 'entering' particle
! or is already marked as a potential exiting particle
         IF(.NOT.IS_NORMAL(LL)) CYCLE

         CELL_ID = DG_PIJK(LL)

! If no neighboring facet in the surrounding 27 cells, then exit
!        IF (NO_NEIGHBORING_FACET_DES(DG_PIJK(LL))) THEN
         IF(CELLNEIGHBOR_FACET_NUM(CELL_ID) < 1) THEN
            WALL_COLLISION_FACET_ID(:,LL) = -1
            WALL_COLLISION_PFT(:,:,LL) = 0.0d0
            CYCLE
         ENDIF

! Check particle LL for wall contacts
         RADSQ = DES_RADIUS(LL)*DES_RADIUS(LL)

         particle_max(:) = des_pos_new(:, LL) + des_radius(LL)
         particle_min(:) = des_pos_new(:, LL) - des_radius(LL)

         DO CELL_COUNT = 1, cellneighbor_facet_num(cell_id)

            axis = cellneighbor_facet(cell_id)%extentdir(cell_count)

            NF = cellneighbor_facet(cell_id)%p(cell_count)

! Compute particle-particle VDW cohesive short-range forces
            IF(USE_COHESION .AND. VAN_DER_WAALS) THEN

               CALL ClosestPtPointTriangle(DES_POS_NEW(:,LL),          &
                  VERTEX(:,:,NF), CLOSEST_PT(:))

               DIST(:) = CLOSEST_PT(:) - DES_POS_NEW(:,LL)
               DISTSQ = DOT_PRODUCT(DIST, DIST)
               R_LM = 2*DES_RADIUS(LL)

               IF(DISTSQ < (R_LM+WALL_VDW_OUTER_CUTOFF)**2) THEN
                  IF(DISTSQ > (WALL_VDW_INNER_CUTOFF+R_LM)**2) THEN
                     DistApart = (SQRT(DISTSQ)-R_LM)
                     FORCE_COH = WALL_HAMAKER_CONSTANT*DES_RADIUS(LL) /&
                        (12d0*DistApart**2)*(Asperities/(Asperities +  &
                        DES_RADIUS(LL)) + ONE/(ONE+Asperities/         &
                        DistApart)**2)
                  ELSE
                     FORCE_COH = 2d0*PI*SURFACE_ENERGY*DES_RADIUS(LL)* &
                        (Asperities/(Asperities+DES_RADIUS(LL)) + ONE/ &
                        (ONE+Asperities/WALL_VDW_INNER_CUTOFF)**2 )
                  ENDIF
                  FC(:,LL) = FC(:,LL) + DIST(:)*FORCE_COH/SQRT(DISTSQ)
               ENDIF
            ENDIF

            if (cellneighbor_facet(cell_id)%extentmin(cell_count) >    &
               particle_max(axis)) then
               call remove_collision(LL, nf, wall_collision_facet_id)
               cycle
            endif

            if (cellneighbor_facet(cell_id)%extentmax(cell_count) <    &
               particle_min(axis)) then
               call remove_collision(LL, nf, wall_collision_facet_id)
               cycle
            endif

! Checking all the facets is time consuming due to the expensive
! separating axis test. Remove this facet from contention based on
! a simple orthogonal projection test.

! Parametrize a line as p = p_0 + t normal and intersect with the
! triangular plane. If t>0, then point is on the non-fluid side of
! the plane. If the plane normal is assumed to point toward the fluid.

! -undefined, because non zero values will imply the sphere center
! is on the non-fluid side of the plane. Since the testing
! is with extended plane, this could very well happen even
! when the particle is well inside the domain (assuming the plane
! normal points toward the fluid). See the pic below. So check
! only when line_t is negative

!                 \   Solid  /
!                  \  Side  /
!                   \      /
!                    \    /
! Wall 1, fluid side  \  /  Wall 2, fluid side
!                      \/
!                        o particle
!
! line_t will be positive for wall 1 (incorrectly indicating center
! is outside the domain) and line_t will be negative for wall 2.
!
! Therefore, only stick with this test when line_t is negative and let
! the separating axis test take care of the other cases.

! Since this is for checking static config, line's direction is the
! same as plane's normal. For moving particles, the line's normal will
! be along the point joining new and old positions.

            line_t = DOT_PRODUCT(VERTEX(1,:,NF) - des_pos_new(:,LL),&
               NORM_FACE(:,NF))

! k - rad >= tol_orth, where k = -line_t, then orthogonal
! projection is false. Substituting for k
! => line_t + rad <= -tol_orth
! choosing tol_orth = 0.01% of des_radius = 0.0001*des_radius

! Orthogonal projection will detect false positives even
! when the particle does not overlap the triangle.
! However, if the orthogonal projection shows no overlap, then
! that is a big fat negative and overlaps are not possible.
            if((line_t.le.-1.0001d0*des_radius(LL))) then  ! no overlap
               call remove_collision(LL,nf,wall_collision_facet_id)
               CYCLE
            ENDIF

            CALL ClosestPtPointTriangle(DES_POS_NEW(:,LL),             &
               VERTEX(:,:,NF), CLOSEST_PT(:))

            DIST(:) = CLOSEST_PT(:) - DES_POS_NEW(:,LL)
            DISTSQ = DOT_PRODUCT(DIST, DIST)

            IF(DISTSQ .GE. RADSQ) THEN !No overlap exists
               call remove_collision(LL,nf,wall_collision_facet_id)
               CYCLE
            ENDIF

            MAX_DISTSQ = DISTSQ
            MAX_NF = NF

! Assign the collision normal based on the facet with the
! largest overlap.
            NORMAL(:) = DIST(:)/sqrt(DISTSQ)

! Facet's normal is correct normal only when the intersection is with
! the face. When the intersection is with edge or vertex, then the
! normal is based on closest pt and sphere center. The definition above
! of the normal is generic enough to account for differences between
! vertex, edge, and facet.

! Calculate the particle/wall overlap.
            DISTMOD = SQRT(MAX_DISTSQ)
            OVERLAP_N = DES_RADIUS(LL) - DISTMOD

! Calculate the translational relative velocity
            CALL CFRELVEL_WALL(LL, V_REL_TRANS_NORM,VREL_T,           &
               NORMAL, DISTMOD)

! Calculate the spring model parameters.
            phaseLL = PIJK(LL,5)

! Hertz vs linear spring-dashpot contact model
            IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN) THEN
               sqrt_overlap = SQRT(OVERLAP_N)
               KN_DES_W = hert_kwn(phaseLL)*sqrt_overlap
               KT_DES_W = hert_kwt(phaseLL)*sqrt_overlap
               sqrt_overlap = SQRT(sqrt_overlap)
               ETAN_DES_W = DES_ETAN_WALL(phaseLL)*sqrt_overlap
               ETAT_DES_W = DES_ETAT_WALL(phaseLL)*sqrt_overlap
            ELSE
               KN_DES_W = KN_W
               KT_DES_W = KT_W
               ETAN_DES_W = DES_ETAN_WALL(phaseLL)
               ETAT_DES_W = DES_ETAT_WALL(phaseLL)
            ENDIF

! Calculate the normal contact force
            FNORM(:) = -(KN_DES_W * OVERLAP_N * NORMAL(:) + &
               ETAN_DES_W * V_REL_TRANS_NORM * NORMAL(:))

! Calculate the tangential displacement.
            OVERLAP_T(:) = DTSOLID*VREL_T(:) + GET_COLLISION(LL,       &
               NF, WALL_COLLISION_FACET_ID, WALL_COLLISION_PFT)
            MAG_OVERLAP_T = sqrt(DOT_PRODUCT(OVERLAP_T, OVERLAP_T))

! Check for Coulombs friction law and limit the maximum value of the
! tangential force on a particle in contact with a wall.
            IF(MAG_OVERLAP_T > 0.0) THEN
! Tangential froce from spring.
               FTMD = KT_DES_W*MAG_OVERLAP_T
! Max force before the on set of frictional slip.
               FNMD = MEW_W*sqrt(DOT_PRODUCT(FNORM,FNORM))
! Direction of tangential force.
               TANGENT = OVERLAP_T/MAG_OVERLAP_T
               IF(FTMD < FNMD) THEN
                  FTAN = -FTMD * TANGENT
               ELSE
                  FTAN = -FNMD * TANGENT
                  OVERLAP_T = (FNMD/KT_DES_W) * TANGENT
               ENDIF
            ELSE
               FTAN = 0.0
            ENDIF
! Add in the tangential dashpot damping force
            FTAN = FTAN - ETAT_DES_W*VREL_T(:)

! Save the tangential displacement.
            CALL UPDATE_COLLISION(OVERLAP_T, LL, NF,                   &
               WALL_COLLISION_FACET_ID, WALL_COLLISION_PFT)

! Add the collision force to the total forces acting on the particle.
            FC(:,LL) = FC(:,LL) + FNORM(:) + FTAN(:)

! Add the torque force to the toal torque acting on the particle.
            TOW(:,LL) = TOW(:,LL) + DISTMOD*DES_CROSSPRDCT(NORMAL,FTAN)

         ENDDO

      ENDDO
!$omp end do
!$omp end parallel

      RETURN

       contains

!......................................................................!
!  Function: GET_COLLISION                                             !
!                                                                      !
!  Purpose: Return the integrated (t0->t) tangential displacement.     !
!......................................................................!
      FUNCTION GET_COLLISION(LLL,FACET_ID,WALL_COLLISION_FACET_ID,     &
          WALL_COLLISION_PFT)
      use stl_preproc_des
      use error_manager

      IMPLICIT NONE

      DOUBLE PRECISION :: GET_COLLISION(DIMN)
      INTEGER, INTENT(IN) :: LLL,FACET_ID
      INTEGER, INTENT(INOUT) :: WALL_COLLISION_FACET_ID(:,:)
      DOUBLE PRECISION, INTENT(INOUT) :: WALL_COLLISION_PFT(:,:,:)
      INTEGER :: CC, FREE_INDEX, LC, dgIJK

      free_index = -1

      do cc = 1, COLLISION_ARRAY_MAX
         if (facet_id == wall_collision_facet_id(cc,LLL)) then
            get_collision(:) = wall_collision_PFT(:,cc,LLL)
            return
         else if (-1 == wall_collision_facet_id(cc,LLL)) then
            free_index = cc
         endif
      enddo

! Overwrite old data. This is needed because a particle moving from
! one dg cell to another may no longer 'see' an STL before it moved
! out of contact range. Therefore, the 'remove_collision' function
! does not get called to cleanup the stale data.
      if(-1 == free_index) then
         dgIJK=DG_PIJK(LLL)
         cc_lp: do cc=1, COLLISION_ARRAY_MAX
            do lc=1, cellneighbor_facet_num(dgIJK)
               if(wall_collision_facet_id(cc,LLL) == &
                  cellneighbor_facet(dgIJK)%p(LC))  cycle cc_lp
            enddo
            free_index = cc
            exit cc_lp
         enddo cc_lp
      endif


      if(-1 == free_index) then

         do cc = 1, COLLISION_ARRAY_MAX
            call write_this_stl(wall_collision_facet_id(cc,LLL))
         enddo
         call write_stls_this_dg(dg_pijk(LLL))

         CALL INIT_ERR_MSG("CALC_COLLISION_WALL_MOD: GET_COLLISION")
         WRITE(ERR_MSG, 1100) LLL, CC
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

      else
         wall_collision_facet_id(free_index,LLL) = facet_id
         wall_collision_PFT(:,free_index,LLL) = ZERO
         get_collision(:) = wall_collision_PFT(:,free_index,LLL)
         return
      endif

 1100 FORMAT('Error: COLLISION_ARRAY_MAX too small.'/'Particle: ',I9,/ &
           'Facet:    ',I9)

      END FUNCTION GET_COLLISION


!......................................................................!
!  Function: UPDATE_COLLISION                                          !
!                                                                      !
!  Purpose: Update the integrated (t0->t) tangential displacement.     !
!......................................................................!
      SUBROUTINE UPDATE_COLLISION(PFT, LLL, FACET_ID,                  &
         WALL_COLLISION_FACET_ID, WALL_COLLISION_PFT)

      use error_manager
      implicit none

      DOUBLE PRECISION, INTENT(IN) :: PFT(DIMN)
      INTEGER, INTENT(IN) :: LLL,FACET_ID
      INTEGER, INTENT(IN) :: WALL_COLLISION_FACET_ID(:,:)
      DOUBLE PRECISION, INTENT(INOUT) :: WALL_COLLISION_PFT(:,:,:)
      INTEGER :: CC

      do cc = 1, COLLISION_ARRAY_MAX
         if (facet_id == wall_collision_facet_id(cc,LLL)) then
            wall_collision_PFT(:,cc,LLL) = PFT(:)
            return
         endif
      enddo

      CALL INIT_ERR_MSG("CALC_COLLISION_WALL_MOD: UPDATE_COLLISION")
      WRITE(ERR_MSG, 1100)
      CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1100 FORMAT('Error: COLLISION_ARRAY_MAX too small. ')

      END SUBROUTINE UPDATE_COLLISION

!......................................................................!
!  Function: REMOVE_COLLISION                                          !
!                                                                      !
!  Purpose: Clear the integrated (t0->t) tangential displacement once  !
!  the collision is over (contact ended).                              !
!......................................................................!
      SUBROUTINE REMOVE_COLLISION(LLL,FACET_ID,WALL_COLLISION_FACET_ID)

      use error_manager

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: LLL,FACET_ID
      INTEGER, INTENT(INOUT) :: WALL_COLLISION_FACET_ID(:,:)
      INTEGER :: CC

      DO CC = 1, COLLISION_ARRAY_MAX
         IF (FACET_ID == WALL_COLLISION_FACET_ID(CC,LLL)) THEN
            WALL_COLLISION_FACET_ID(CC,LLL) = -1
            RETURN
         ENDIF
      ENDDO

      END SUBROUTINE REMOVE_COLLISION

      END SUBROUTINE CALC_DEM_FORCE_WITH_WALL_STL


!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: CFRELVEL_WALL                                           !
!                                                                      !
!  Purpose: Calculate the normal and tangential components of the      !
!  relative velocity between a particle and wall contact.              !
!                                                                      !
!  Comments: Only the magnitude of the normal component is returned    !
!  whereas the full tangential vector is returned.                     !
!----------------------------------------------------------------------!
      SUBROUTINE CFRELVEL_WALL(LL, VRN, VRT, NORM, DIST)

! Particle translational velocity
      use discretelement, only: DES_VEL_NEW
! Particle rotational velocity
      use discretelement, only: OMEGA_NEW
! Spatial array size (parameter)
      use discretelement, only: DIMN
! Function for calculating the cross prodcut
      use discretelement, only: DES_CROSSPRDCT

      IMPLICIT NONE

! Dummy arguments:
!---------------------------------------------------------------------//
! Particle index.
      INTEGER, INTENT(IN) :: LL
! Magnitude of the total relative translational velocity.
      DOUBLE PRECISION, INTENT(OUT):: VRN
! Total relative translational velocity (vector).
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(OUT):: VRT
! Unit normal from particle center to closest point on stl (wall)
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: NORM
! Distance between particle center and stl (wall).
      DOUBLE PRECISION, INTENT(IN) :: DIST

! Local variables
!---------------------------------------------------------------------//
! Additional relative translational motion due to rotation
      DOUBLE PRECISION, DIMENSION(DIMN) :: V_ROT
! Total relative velocity at contact point
      DOUBLE PRECISION, DIMENSION(DIMN) :: VRELTRANS

! Total relative velocity + rotational contribution
      V_ROT = DIST*OMEGA_NEW(:,LL)
      VRELTRANS(:) =  DES_VEL_NEW(:,LL) + DES_CROSSPRDCT(V_ROT, NORM)

! magnitude of normal component of relative velocity (scalar)
      VRN = DOT_PRODUCT(VRELTRANS,NORM)

! total relative translational slip velocity at the contact point
! Equation (8) in Tsuji et al. 1992
      VRT(:) =  VRELTRANS(:) - VRN*NORM(:)

      RETURN
      END SUBROUTINE CFRELVEL_WALL

      end module CALC_COLLISION_WALL

