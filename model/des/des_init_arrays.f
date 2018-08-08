!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                         !
!  Subrourtine: DES_INIT_ARRAYS                                           !
!  Author: Jay Boyalakuntla                              Date: 12-Jun-04  !
!                                                                         !
!  Purpose: Initialize arrays at the start of the simulation. Note that   !
!  arrays based on the number of particles (MAX_PIP) should be added to   !
!  the DES_INIT_PARTICLE_ARRAYS as they need to be reinitialized after    !
!  particle arrays are grown (see PARTICLE_GROW).                         !
!                                                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_INIT_ARRAYS

      USE param
      USE param1
      USE discretelement
      USE indices
      USE geometry
      USE compar
      USE physprop
      USE des_bc
      USE run
      use desgrid
      use desmpi
      USE des_thermo
      USE des_rxns

      IMPLICIT NONE

!-----------------------------------------------

      PINC(:) = ZERO

      DES_U_s(:,:) = ZERO
      DES_V_s(:,:) = ZERO
      DES_W_s(:,:) = ZERO
      DES_ROP_S(:,:) = ZERO
      DES_ROP_SO(:,:) = ZERO

      P_FORCE(:,:) = ZERO

      IF(allocated(DRAG_AM)) DRAG_AM = ZERO
      IF(allocated(DRAG_BM)) DRAG_BM = ZERO

      F_GDS = ZERO
      VXF_GDS = ZERO

      IF (DES_CONTINUUM_HYBRID) THEN
         F_SDS = ZERO
         VXF_SDS = ZERO
         SDRAG_AM = ZERO
         SDRAG_BM = ZERO
      ENDIF

      GRAV(:) = ZERO

      IF(ENERGY_EQ)THEN
         avgDES_T_s(:) = ZERO
         DES_ENERGY_SOURCE(:) = ZERO
      ENDIF

      CALL DES_INIT_PARTICLE_ARRAYS(1,MAX_PIP)

      RETURN
      END SUBROUTINE DES_INIT_ARRAYS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                         !
!  Subrourtine: DES_INIT_PARTICLE_ARRAYS                                  !
!  Author: Jay Boyalakuntla                              Date: 12-Jun-04  !
!                                                                         !
!  Purpose: Initialize particle arrays. The upper and lower bounds are    !
!  passed so that after resizing particle arrays (see GROW_PARTICLE) the  !
!  new portions of the arrays can be initialized.                         !
!                                                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_INIT_PARTICLE_ARRAYS(LB,UB)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use des_rxns
      use des_thermo
      use desgrid
      use desmpi
      use discretelement
      use functions
      use mfix_pic, only: AVGSOLVEL_P, EPG_P
      use mfix_pic, only: MPPIC, DES_STAT_WT, PS_GRAD
      use particle_filter, only: FILTER_CELL, FILTER_WEIGHT
      use particle_filter, only: FILTER_SIZE
      use run, only: ANY_SPECIES_EQ
      use run, only: ENERGY_EQ

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: LB, UB
      INTEGER :: II

      IGLOBAL_ID(LB:UB) = 0

! Physical properties:
      DES_RADIUS(LB:UB) = ZERO
      RO_Sol(LB:UB) = ZERO
      PVOL(LB:UB) = ZERO
      PMASS(LB:UB) = ZERO
      OMOI(LB:UB) = ZERO

! Particle position, velocity, etc
      DES_POS_NEW(:,LB:UB) = ZERO
      DES_VEL_NEW(:,LB:UB) = ZERO
      OMEGA_NEW(:,LB:UB) = ZERO
      IF(PARTICLE_ORIENTATION) THEN
         ORIENTATION(1,:) = INIT_ORIENTATION(1)
         ORIENTATION(2,:) = INIT_ORIENTATION(2)
         ORIENTATION(3,:) = INIT_ORIENTATION(3)
      ENDIF

! Particle state flag
      DO II = LB, UB
         call set_nonexistent(II)
      ENDDO
      NEIGHBOR_INDEX(:) = 0

! DES grid bin information
      DG_PIJK(LB:UB) = -1
      DG_PIJKPRV(LB:UB) = -1
      IGHOST_UPDATED(LB:UB) = .false.

! Fluid cell bin information
      PIJK(LB:UB,:) = 0

! Translation and rotational forces
      FC(:,LB:UB) = ZERO
      TOW(:,LB:UB) = ZERO

! Collision data
      WALL_COLLISION_FACET_ID(:,LB:UB) = -1
      WALL_COLLISION_PFT(:,:,LB:UB) = ZERO

! Initializing user defined array
      IF(DES_USR_VAR_SIZE > 0) &
         DES_USR_VAR(:,LB:UB) = ZERO

! Particle center drag coefficient and explicit drag force
      F_GP(LB:UB) = ZERO
      DRAG_FC(:,LB:UB) = ZERO


! Interpolation variables.
      IF(FILTER_SIZE > 0)THEN
         FILTER_CELL(:,LB:UB) = -1
         FILTER_WEIGHT(:,LB:UB) = ZERO
      ENDIF

! MPPIC variables
      IF(MPPIC) THEN
         DES_STAT_WT(LB:UB) = ZERO
         PS_GRAD(LB:UB,:) = ZERO
         AVGSOLVEL_P(:,LB:UB) = ZERO
         EPG_P(LB:UB) = ZERO
      ENDIF

! Higher order time integration variables.
      IF (DO_OLD) THEN
         DES_POS_OLD(:,LB:UB) = ZERO
         DES_VEL_OLD(:,LB:UB) = ZERO
         DES_ACC_OLD(:,LB:UB) = ZERO
         OMEGA_OLD(:,LB:UB) = ZERO
         ROT_ACC_OLD(:,LB:UB) = ZERO
      ENDIF

! Energy equation variables.
      IF(ENERGY_EQ)THEN
         DES_T_s_OLD(LB:UB) = ZERO
         DES_T_s_NEW(LB:UB) = ZERO
         DES_C_PS(LB:UB) = ZERO
         DES_X_s(LB:UB,:) = ZERO
         Q_Source(LB:UB) = ZERO
         IF (INTG_ADAMS_BASHFORTH) &
            Q_Source0(LB:UB) = ZERO
      ENDIF

! Chemical reaction variables.
      IF(ANY_SPECIES_EQ)THEN
         DES_R_sp(LB:UB,:) = ZERO
         DES_R_sc(LB:UB,:) = ZERO
         IF (INTG_ADAMS_BASHFORTH) THEN
            dMdt_OLD(LB:UB) = ZERO
            dXdt_OLD(LB:UB,:) = ZERO
         ENDIF
         Qint(LB:UB) = ZERO
      ENDIF


! Cohesion VDW forces
      IF(USE_COHESION) THEN
         PostCohesive (LB:UB) = ZERO
      ENDIF


      RETURN
      END SUBROUTINE DES_INIT_PARTICLE_ARRAYS



