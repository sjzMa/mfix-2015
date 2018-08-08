      MODULE GENERATE_PARTICLES

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: GENERATE_PARTICLE_CONFIG                                C
!                                                                      C
!  Purpose: Generate particle configuration based on maximum particle  C
!           radius and filling from top to bottom within specified     C
!           bounds                                                     C
!                                                                      C
!                                                                      C
!  Authors: Rahul Garg                                Date: 19-Mar-14  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GENERATE_PARTICLE_CONFIG

      use mfix_pic, only: MPPIC
      use discretelement, only: PIP, PARTICLES
! Flag indicating that the IC region is defined.
      use ic, only: IC_DEFINED
! Parameter for detecting unspecified values, zero, and one
      use param1, only: ONE
! Maximum number of initial conditions
      use param, only: DIMENSION_IC
! IC Region gas volume fraction.
      use ic, only: IC_EP_G

      use mpi_utility

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

      INTEGER :: ICV

! Initialize the error manager.
      CALL INIT_ERR_MSG("Generate_Particle_Config")

      DO ICV = 1, DIMENSION_IC

         IF(.NOT.IC_DEFINED(ICV)) CYCLE
         IF(IC_EP_G(ICV) == ONE) CYCLE

         IF(MPPIC) THEN
            CALL GENERATE_PARTICLE_CONFIG_MPPIC(ICV)
         ELSE
            CALL GENERATE_PARTICLE_CONFIG_DEM(ICV)
         ENDIF

      ENDDO

      CALL GLOBAL_SUM(PIP,PARTICLES)

      WRITE(ERR_MSG, 1004) PARTICLES
 1004 FORMAT(/,'Total number of particles in the system: ',I15)

      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE GENERATE_PARTICLE_CONFIG



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: GENERATE_PARTICLE_CONFIG                                !
!  Authors: Rahul Garg                               Date: 21-Mar-2014 !
!                                                                      !
!  Purpose: Generate particle configuration for DEM solids for each IC !
!           region. Now using the particle linked lists for initial    !
!           build                                                      !
!           This routine will ultimately supersede the older rouine    !
!           that has not been deleted yet                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GENERATE_PARTICLE_CONFIG_DEM(ICV)


! Global Variables:
! particle radius and density
      use discretelement, only: DES_RADIUS, RO_Sol
! particle position new and old
      use discretelement, only: DES_POS_NEW, DES_POS_OLD
! particle velocity new and old
      use discretelement, only: DES_VEL_NEW, DES_VEL_OLD
! Simulation dimension (2D/3D)
      use discretelement, only: DIMN
! Number of particles in the system (current)
      use discretelement, only: PIP
! Number of DEM solids phases.
      use discretelement, only: DES_MMAX
! Flag to use _OLD variables
      use discretelement, only: DO_OLD
! Angular velocity
      use discretelement, only: OMEGA_OLD, OMEGA_NEW, PIJK
! DEM solid phase diameters and densities.
      use physprop, only: D_p0, RO_s0, SMAX
! IC Region solids volume fraction.
      use ic, only: IC_EP_S

! Constant: 3.14159...
      use constant, only: PI
! min and max physical co-ordinates of IC regions in each direction
      use ic, only: IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t
! initally specified velocity field and granular temperature
      use ic, only: IC_U_s, IC_V_s, IC_W_s, IC_Theta_M
! Flag to extend the lattice distribution in a given IC to available area
      use ic, only: IC_DES_FIT_TO_REGION
! Parameter for detecting unspecified values, zero, and one
      use param1, only: UNDEFINED, UNDEFINED_I, ZERO, ONE, Half
! Parameter for small and large numbers
      use param1, only: SMALL_NUMBER, LARGE_NUMBER

! to access random number generator subroutines
      use randomno
      use mpi_utility
      use functions, only: SET_NORMAL, FUNIJK, FLUID_AT

      use error_manager

! direction wise spans of the domain and grid spacing in each direction
      use geometry, only: xlength, ylength, zlength

      use cutcell, only : CARTESIAN_GRID, CUT_CELL_AT
      use STL_PREPROC_DES, only: CHECK_IF_PARTICLE_OVERLAPS_STL
      use run, only: solids_model
      use des_allocate, only: PARTICLE_GROW

      use discretelement, only: MAX_RADIUS

      use discretelement, only: XE, YN, ZT

      use param, only: DIM_M, DIMENSION_I, DIMENSION_J, DIMENSION_K
      use functions, only: IS_ON_MYPE_WOBND

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ICV

! Local variables
!---------------------------------------------------------------------//
! Starting positions in the axial directions
      DOUBLE PRECISION :: xINIT, yINIT, zINIT
! Fractor used to scale particle diameter
      DOUBLE PRECISION :: lFAC
! Particle position and velocity
      DOUBLE PRECISION :: POS(3), VEL(3)
! Number of particles in the lattice
      INTEGER :: SEED_X, SEED_Y, SEED_Z
! Loop indices phase/fluid cell
      INTEGER :: M, MM, I, J, K, IJK, LB, UB
! Loop indicies for seeding
      INTEGER :: II, JJ, KK
! Start and end bound for IC region.
      DOUBLE PRECISION :: IC_START(3), IC_END(3)
! Volume and lengths of the IC Region
      DOUBLE PRECISION :: DOM_VOL, DOML(3)
! Flag to skip the particle
      LOGICAL :: SKIP
! Diameter adjusted for space padding
      DOUBLE PRECISION :: ADJ_DIA
! Number of particles calculated from volume fracton
      INTEGER :: rPARTS(DIM_M), tPARTS
! Spacing between particles.
      DOUBLE PRECISION :: lDEL, lDX, lDY, lDZ
! Flag that the setup failed to fit the particles to the IC region
      LOGICAL :: FIT_FAILED
! Number of seeded particles
      INTEGER :: pCOUNT(DIM_M), tCOUNT

      DOUBLE PRECISION :: SOLIDS_DATA(0:DIM_M)

      LOGICAL :: VEL_FLUCT
      DOUBLE PRECISION :: VEL_SIG
      DOUBLE PRECISION, ALLOCATABLE :: randVEL(:,:)

!......................................................................!

      CALL INIT_ERR_MSG("GENERATE_PARTICLE_CONFIG_DEM")

      WRITE(ERR_MSG,"(2/,'Generating initial particle configuration:')")
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      SOLIDS_DATA = ZERO
      CALL GET_IC_VOLUME(ICV, SOLIDS_DATA(0))

! setting particle seed spacing grid to be slightly greater than
! the maximum particle diameter. seed at ~particle radii
      lFAC = 1.05D0

! Setup local arrays with IC region bounds.
      IC_START(1)=IC_X_W(ICV);   IC_END(1)=IC_X_E(ICV)
      IC_START(2)=IC_Y_S(ICV);   IC_END(2)=IC_Y_N(ICV)
      IC_START(3)=IC_Z_B(ICV);   IC_END(3)=IC_Z_T(ICV)

      DOML = IC_END-IC_START
      IF(NO_K) DOML(3)=DZ(1)

! Volume of the IC region
      DOM_VOL = DOML(1)*DOML(2)*DOML(3)

      rPARTS=0
      VEL_FLUCT = .FALSE.
      DO M=1,SMAX+DES_MMAX
         IF(SOLIDS_MODEL(M) == 'DEM') THEN
! Number of particles for phase M
            rPARTS(M) = &
               floor((6.0d0*IC_EP_S(ICV,M)*DOM_VOL)/(PI*(D_P0(M)**3)))
! Set flags for random velocity fluctuations
            VEL_FLUCT = VEL_FLUCT .OR. (IC_Theta_M(ICV,M) /= 0.0d0)
         ENDIF
      ENDDO

! Total number of particles in this IC region.
      tPARTS = sum(rPARTS)
      IF(tPARTS == 0) RETURN

      ADJ_DIA = 2.0d0*MAX_RADIUS*lFAC

! Calculate velocity fluctuations if specified.
      IF(VEL_FLUCT) THEN
         allocate(randVEL(tPARTS,3))
         LB=1
         DO M=1,SMAX+DES_MMAX
            IF(SOLIDS_MODEL(M) == 'DEM' .AND. rPARTS(M) > 0.0d0) THEN
               UB=LB+int(rPARTS(M))-1
               IF(IC_Theta_M(ICV,1) /= 0.0d0)THEN
                  VEL_SIG = sqrt(IC_Theta_M(ICV,1))
                  CALL NOR_RNO(randVEL(LB:UB,1), IC_U_s(ICV,1), VEL_SIG)
                  CALL NOR_RNO(randVEL(LB:UB,2), IC_V_s(ICV,2), VEL_SIG)
                  IF(DO_K) CALL NOR_RNO(randVEL(LB:UB,3), &
                     IC_W_s(ICV,3), VEL_SIG)
               ELSE
                  randVEL(LB:UB,1) = ZERO
                  randVEL(LB:UB,2) = ZERO
                  IF(DO_K) randVEL(LB:UB,3) = ZERO
               ENDIF
               LB=UB+1
            ENDIF
         ENDDO
      ENDIF


! Attempt to seed particle throughout the IC region
      FIT_FAILED=.FALSE.
      IF(IC_DES_FIT_TO_REGION(ICV)) THEN
         IF(NO_K) THEN
            lDEL = (DOML(1)-ADJ_DIA)*(DOML(2)-ADJ_DIA)
            lDEL = (lDEL/dble(tPARTS))**(1.0/2.0)
            SEED_X = max(1,ceiling((DOML(1)-ADJ_DIA)/lDEL))
            SEED_Y = max(1,ceiling((DOML(2)-ADJ_DIA)/lDEL))
            SEED_Z = 1
         ELSE
            lDEL = (DOML(1)-ADJ_DIA)*(DOML(2)-ADJ_DIA)*(DOML(3)-ADJ_DIA)
            lDEL = (lDEL/dble(tPARTS))**(1.0/3.0)
            SEED_X = max(1,ceiling((DOML(1)-ADJ_DIA)/lDEL))
            SEED_Y = max(1,ceiling((DOML(2)-ADJ_DIA)/lDEL))
            SEED_Z = max(1,ceiling((DOML(3)-ADJ_DIA)/lDEL))
         ENDIF
         FIT_FAILED=(dble(SEED_X*SEED_Y*SEED_Z) < tPARTS)
      ENDIF

! Generic filling
      IF(.NOT.IC_DES_FIT_TO_REGION(ICV) .OR. FIT_FAILED) THEN
         SEED_X = max(1,floor((IC_END(1)-IC_START(1)-ADJ_DIA)/ADJ_DIA))
         SEED_Y = max(1,floor((IC_END(2)-IC_START(2)-ADJ_DIA)/ADJ_DIA))
         SEED_Z = max(1,floor((IC_END(3)-IC_START(3)-ADJ_DIA)/ADJ_DIA))
      ENDIF

      lDX = DOML(1)/dble(SEED_X)
      lDY = DOML(2)/dble(SEED_Y)
      IF(DO_K) THEN
         lDZ = DOML(3)/dble(SEED_Z)
      ELSE
         lDZ = 0.0d0
      ENDIF

      xINIT = IC_START(1)+HALF*lDX
      yINIT = IC_START(2)+HALF*lDY
      zINIT = IC_START(3)+HALF*lDZ


      M=1
      pCOUNT = 0
      tCOUNT = 0
      JJ_LP: DO  JJ=1, SEED_Y
      KK_LP: DO  KK=1, SEED_Z
      II_LP: DO  II=1, SEED_X

! Exit if all particles were seeded.
         IF(tCOUNT > int(tPARTS)) THEN
            EXIT JJ_LP
! Find the next phase that needs to be seeded
         ELSEIF(pCOUNT(M) > int(rPARTS(M))) THEN
            MM_LP: DO MM=M+1,SMAX+DES_MMAX
               IF(rPARTS(MM) > 0.0) THEN
                  M=MM
                  EXIT MM_LP
               ENDIF
            ENDDO MM_LP
            IF(MM > SMAX+DES_MMAX) THEN
               EXIT JJ_LP
            ELSEIF(IC_Theta_M(ICV,MM) /= 0.0d0) THEN
               VEL_SIG = sqrt(IC_Theta_M(ICV,MM))
               CALL NOR_RNO(randVEL(:,1), 0.0d0, VEL_SIG)
               CALL NOR_RNO(randVEL(:,2), 0.0d0, VEL_SIG)
               IF(DO_K) CALL NOR_RNO(randVEL(:,3), 0.0d0, VEL_SIG)
            ENDIF
         ENDIF

         pCOUNT(M) = pCOUNT(M) + 1
         tCOUNT = tCOUNT + 1

! Position the particle
         POS(1) = xINIT + (II-1)*lDX
         POS(2) = YINIT + (JJ-1)*lDY
         POS(3) = ZINIT + (KK-1)*lDZ

! Bin the parcel to the fuild grid.
         K=1
         IF(DO_K) CALL PIC_SEARCH(K, POS(3), ZT, DIMENSION_K, KMIN2, KMAX2)
         CALL PIC_SEARCH(J, POS(2), YN, DIMENSION_J, JMIN2, JMAX2)
         CALL PIC_SEARCH(I, POS(1), XE, DIMENSION_I, IMIN2, IMAX2)

! Skip cells that are not part of the local fuild domain.
         IF(.NOT.IS_ON_MYPE_WOBND(I,J,K)) CYCLE
         IF(DEAD_CELL_AT(I,J,K)) CYCLE

         IJK = FUNIJK(I,J,K)
         IF(.NOT.FLUID_AT(IJK)) CYCLE

         IF(CARTESIAN_GRID) THEN
            CALL CHECK_IF_PARTICLE_OVERLAPS_STL(POS, I, J, K, SKIP)
            IF(SKIP) CYCLE
         ENDIF

         PIP = PIP + 1
         CALL PARTICLE_GROW(PIP)

         CALL SET_NORMAL(PIP)

         IF(VEL_FLUCT) THEN
            VEL(1) = randVEL(tCOUNT,1)
            VEL(2) = randVEL(tCOUNT,2)
            VEL(3) = randVEL(tCOUNT,3)
         ELSE
            VEL(1) = IC_U_s(ICV,M)
            VEL(2) = IC_V_s(ICV,M)
            VEL(3) = IC_W_s(ICV,M)
         ENDIF
         IF(NO_K) VEL(3) = 0.0d0


         DES_POS_NEW(:,PIP) = POS(:)
         DES_VEL_NEW(:,PIP) = VEL(:)
         OMEGA_NEW(:,PIP) = 0.0d0

         DES_RADIUS(PIP) = D_P0(M)*HALF
         RO_SOL(PIP) =  RO_S0(M)

         PIJK(PIP,1) = I
         PIJK(PIP,2) = J
         PIJK(PIP,3) = K
         PIJK(PIP,4) = IJK
         PIJK(PIP,5) = M

         IF(DO_OLD) THEN
            DES_VEL_OLD(:,PIP) = DES_VEL_NEW(:,PIP)
            DES_POS_OLD(:,PIP) = DES_POS_NEW(:,PIP)
            OMEGA_OLD(:,PIP) = ZERO
         ENDIF

         SOLIDS_DATA(M) = SOLIDS_DATA(M) + 1.0

      ENDDO II_LP
      ENDDO KK_LP
      ENDDO JJ_LP

! Collect the data
      CALL GLOBAL_ALL_SUM(SOLIDS_DATA)

! Verify that the IC region volume is not zero.
      IF(SOLIDS_DATA(0) <= 0.0d0) THEN
         WRITE(ERR_MSG,1000) ICV, SOLIDS_DATA(0)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
 
1000 FORMAT('Error 1000: Invalid IC region volume: IC=',I3,' VOL=',&
         ES15.4,/'Please correct the mfix.dat file.')

      WRITE(ERR_MSG,2000) ICV
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      DO M=1, SMAX+DES_MMAX
         IF(SOLIDS_DATA(M) < SMALL_NUMBER) CYCLE
         WRITE(ERR_MSG,2010) M, int(SOLIDS_DATA(M)), IC_EP_S(ICV,M),   &
            (dble(SOLIDS_DATA(M))*(Pi/6.0d0)*D_P0(M)**3)/SOLIDS_DATA(0)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDDO

      IF(allocated(randVEL)) deallocate(randVEL)

      CALL FINL_ERR_MSG

      RETURN

 2000 FORMAT(/2x,'|',43('-'),'|',/2x,'| IC Region: ',I3,28x,'|',/2x,   &
         '|',43('-'),'|',/2x,'| Phase | Number of |    EPs    |    EP',&
         's    |',/2x,'|   ID  | Particles | Specified |   Actual  |', &
         /2x,'|-------|',3(11('-'),'|'))

 2010 FORMAT(2x,'|  ',I3,'  |',1x,I9,1x,'|',2(1x,ES9.2,1x,'|'),/2x,    &
         '|-------|',3(11('-'),'|'))

      END SUBROUTINE GENERATE_PARTICLE_CONFIG_DEM




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GENERATE_PARTICLE_CONFIG_MMPPIC                         !
!  Author: Rahul Garg                                 Date: 3-May-2011 !
!                                                                      !
!  Purpose: Generates particle position distribution for MPPIC.        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GENERATE_PARTICLE_CONFIG_MPPIC(ICV)

! Number of DES solids phases.
      use discretelement, only: DES_MMAX
! Flag indicating that the IC region is defined.
      use ic, only: IC_DEFINED
! IC Region bulk density (RO_s * EP_s)
      use ic, only: IC_ROP_s
! IC Region gas volume fraction.
      use ic, only: IC_EP_G
! MPPIC specific IC region specification.
      use ic, only: IC_PIC_CONST_NPC, IC_PIC_CONST_STATWT

      use param1, only: UNDEFINED, UNDEFINED_I
      use param1, only: ZERO, ONE, HALF
! Maximum number of IC regions and solids phases
      use param, only: DIMENSION_IC
      use param, only: DIM_M

      use mpi_utility, only: GLOBAL_ALL_SUM

      use error_manager


      IMPLICIT NONE


      INTEGER, INTENT(IN) :: ICV

! Local variables
!---------------------------------------------------------------------//
! Generic loop counters
      INTEGER :: M
! Actual volume of IC region
      DOUBLE PRECISION :: IC_VOL
! Solids data in IC Region by phase:
      DOUBLE PRECISION :: SOLIDS_DATA(0:4*DIM_M)
!......................................................................!

      CALL INIT_ERR_MSG("GENERATE_PARTICLE_CONFIG_MPPIC")

      WRITE(ERR_MSG,"(2/,'Generating initial parcel configuration:')")
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)


      SOLIDS_DATA = ZERO
      CALL GET_IC_VOLUME(ICV, SOLIDS_DATA(0))

      WRITE(ERR_MSG,2000) ICV
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

! Set up the individual solids phases.
      DO M=1, DES_MMAX
         IF(IC_ROP_S(ICV,M) == ZERO) CYCLE
! Seed parcels with constant stastical weight.
         IF(IC_PIC_CONST_STATWT(ICV,M) /= ZERO) THEN
            CALL GPC_MPPIC_CONST_STATWT(ICV, M, SOLIDS_DATA(0), &
               SOLIDS_DATA((4*M-3):(4*M)))
! Seed parcels with a constant number per cell
         ELSEIF(IC_PIC_CONST_NPC(ICV,M) /= 0) THEN
            CALL GPC_MPPIC_CONST_NPC(ICV, M, SOLIDS_DATA(0), &
               SOLIDS_DATA((4*M-3):(4*M)))
         ENDIF
      ENDDO

! Collect the data
      CALL GLOBAL_ALL_SUM(SOLIDS_DATA)

! Verify that the IC region volume is not zero.
      IF(SOLIDS_DATA(0) <= 0.0d0) THEN
         WRITE(ERR_MSG,1000) ICV, SOLIDS_DATA(0)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1000 FORMAT('Error 1000: Invalid IC region volume: IC=',I3,' VOL=',&
         ES15.4,/'Please correct the mfix.dat file.')

! Report solids information for the IC region.
      DO M=1, DES_MMAX
         WRITE(ERR_MSG,2010) M, int(SOLIDS_DATA(4*M-3)),&
            int(SOLIDS_DATA(4*M-3)*SOLIDS_DATA(4*M-2)), &
            SOLIDS_DATA(4*M-2), SOLIDS_DATA(4*M-1),     &
            SOLIDS_DATA(4*M)/SOLIDS_DATA(0)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDDO

      CALL FINL_ERR_MSG

      RETURN

 2000 FORMAT(/2x,'|',67('-'),'|',/2x,'| IC Region: ',I3,52x,'|',/2x,   &
         '|',67('-'),'|',/2x,'| Phase | Num. Comp | Num. Real ',       &
         '| Stastical |    EPs    |    EPs    |',/2x,'|   ID  |  ',    &
         'Parcels  | Particles |   Weight  | Specified |   Actual  |', &
         /2x,'|-------|',5(11('-'),'|'))

 2010 FORMAT(2x,'|  ',I3,'  |',2(1x,I9,1x,'|'),3(1x,ES9.2,1x,'|'),/2x,&
         '|-------|',5(11('-'),'|'))

      END SUBROUTINE GENERATE_PARTICLE_CONFIG_MPPIC


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GPC_MPPIC_CONST_STATWT                                  !
!  Author: J.Musser                                 Date:  26-Aug-2015 !
!                                                                      !
!  Purpose: generates particle position distribution for MPPIC         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GPC_MPPIC_CONST_STATWT(ICV, M, IC_VOL, sDATA)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use cutcell, only : CARTESIAN_GRID

! IC Region Bounds
      use ic, only: IC_X_w, IC_X_e
      use ic, only: IC_Y_s, IC_Y_n
      use ic, only: IC_Z_b, IC_Z_t
! IC Region bulk density (RO_s * EP_s)
      use ic, only: IC_ROP_s
! IC Region velocity field and granular temperature
      use ic, only: IC_U_s, IC_V_s, IC_W_s, IC_Theta_M
! MPPIC specific IC region specification.
      use ic, only: IC_PIC_CONST_STATWT
! Flag for Cartesian cut-cell 
      use cutcell, only : CARTESIAN_GRID
! DES solid phase diameters and densities.
      use discretelement, only: DES_D_p0, DES_RO_s
! Total number of active parcles
      use discretelement, only: PIP
! Parcel position and velocity
      use discretelement, only: DES_POS_NEW, DES_VEL_NEW
! Parcel radius and density
      use discretelement, only: DES_RADIUS, RO_SOL
! Stastical weight of each parcel
      use mfix_pic, only: DES_STAT_WT
! Paricle fluid I/J/K/IJK/Phase information
      use discretelement, only: PIJK
! The east/north/top location of grid partition
      use discretelement, only: XE, YN, ZT
! Fluid grid cell dimensions and mesh size
      USE geometry, only: IMIN2, IMAX2
      USE geometry, only: JMIN2, JMAX2
      USE geometry, only: KMIN2, KMAX2
      USE geometry, only: DO_K, NO_K, DZ

! Global parameters
!----------------------------------------------------------------------//
      use param1, only: ZERO, HALF, ONE
! Maximum cell partitions in the axial directions
      use param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K
! Maximum number of IC regions and solids phases
      use param, only: DIMENSION_IC
! Constant: 3.14159...
      use constant, only: PI


! Module procedures
!----------------------------------------------------------------------//
      use functions, only: FUNIJK, FLUID_AT, DEAD_CELL_AT
      use stl_preproc_des, only: CHECK_IF_PARTICLE_OVERLAPS_STL
      use cutcell, only: cut_cell_at
      use randomno, only: UNI_RNO, NOR_RNO
      use functions, only: SET_NORMAL
      use functions, only: IS_ON_MYPE_WOBND
      use des_allocate, only: PARTICLE_GROW

      use error_manager


      IMPLICIT NONE


! Dummy Arguments
!----------------------------------------------------------------------//
! Index of IC region and solids phase
      INTEGER, INTENT(IN) :: ICV, M
! Specific volume of IC region (accounts for blocked cells)
      DOUBLE PRECISION, INTENT(IN) :: IC_VOL
! Data about solids in the IC region.
      DOUBLE PRECISION, INTENT(OUT) :: sDATA(4)

! Local variables
!----------------------------------------------------------------------//
! Solids phase M volume fraction in IC region
      DOUBLE PRECISION :: EP_SM
! Start/End bounds of IC region
      DOUBLE PRECISION :: IC_START(3), IC_END(3)
! IC region lenghts and volume
      DOUBLE PRECISION :: DOML(3), DOM_VOL
! Number of real particles calculated from volume fracton
      DOUBLE PRECISION ::  rPARTS
! Number of computational parcels
      INTEGER :: cPARTS
! Calculated statistical weight of parcels
      DOUBLE PRECISION :: STAT_WT
! Volume of a parcel and total solids volume
      DOUBLE PRECISION :: sVOL, sVOL_TOT
! Parcel position and velocity
      DOUBLE PRECISION :: POS(3), VEL(3)
! Arrays for assigning random position and velocities
      DOUBLE PRECISION, ALLOCATABLE :: randPOS(:,:), randVEL(:,:)
! Standard deivation of initial velocities
      DOUBLE PRECISION :: VEL_SIG
! Flag to skip parcel
      LOGICAL :: SKIP
! Counter for seeded parcles.
      INTEGER :: SEEDED
! Generic fluid cell indices and loop counters
      INTEGER :: I, J, K, IJK, LC
!......................................................................!

      CALL INIT_ERR_MSG("GENERATE_PARTICLE_CONFIG_MPPIC")

! Setup local arrays with IC region bounds.
      IC_START(1)=IC_X_W(ICV);   IC_END(1)=IC_X_E(ICV)
      IC_START(2)=IC_Y_S(ICV);   IC_END(2)=IC_Y_N(ICV)
      IC_START(3)=IC_Z_B(ICV);   IC_END(3)=IC_Z_T(ICV)

      DOML = IC_END-IC_START
      IF(NO_K) DOML(3)=DZ(1)

! Volume of the IC region
      DOM_VOL = DOML(1)*DOML(2)*DOML(3)
! Solids volume fraction in IC region
      EP_SM = IC_ROP_S(ICV,M)/DES_RO_s(M)

! Total number of real and computational particles
      rPARTS = 6.d0*EP_SM*DOM_VOL/(PI*(Des_D_P0(M)**3.d0))
      cPARTS = max(1,int(rPARTS/real(IC_PIC_CONST_STATWT(ICV,M))))

! The 'actual' statistical weight
      STAT_WT = rPARTS/real(cPARTS)

! Obtain random variations for parcel position 
      ALLOCATE(randPOS(cPARTS,3))

      DO LC = 1, merge(2,3,NO_K)
         CALL UNI_RNO(RANDPOS(1:cPARTS,LC))
      ENDDO
      IF(NO_K) randPOS(:,3) = 0.0d0

! Obtain initial velocities with random variations
      ALLOCATE(randVEL(cPARTS,3))

      randVEL(:,1) = IC_U_s(ICV,M)
      randVEL(:,2) = IC_V_s(ICV,M)
      randVEL(:,3) = merge(IC_W_s(ICV,M), 0.0d0, DO_K)

      VEL_SIG = sqrt(IC_Theta_M(ICV,M))
      IF(VEL_SIG /= ZERO) THEN
         CALL NOR_RNO(randVEL(:,1), IC_U_S(ICV,M), VEL_SIG)
         CALL NOR_RNO(randVEL(:,2), IC_V_S(ICV,M), VEL_SIG)
         IF(DO_K) CALL NOR_RNO(randVEL(:,3), IC_W_S(ICV,M), VEL_SIG)
      ENDIF

! Volume occupied by one parcel
      sVOL = (Pi/6.0d0)*(DES_D_P0(M)**3.d0)*STAT_WT
      sVOL_TOT = IC_VOL*EP_SM

      SEEDED = 0
      DO LC = 1, cPARTS

! Generate the initial parcel position.
         POS = IC_START + DOML*RANDPOS(LC,:)
         VEL = 0.0d0

! Bin the parcel to the fuild grid.
         K=1
         IF(DO_K) CALL PIC_SEARCH(K, POS(3), ZT, DIMENSION_K, KMIN2, KMAX2)
         CALL PIC_SEARCH(J, POS(2), YN, DIMENSION_J, JMIN2, JMAX2)
         CALL PIC_SEARCH(I, POS(1), XE, DIMENSION_I, IMIN2, IMAX2)

! Skip cells that are not part of the local fuild domain.
         IF(.NOT.IS_ON_MYPE_WOBND(I,J,K)) CYCLE
         IF(DEAD_CELL_AT(I,J,K)) CYCLE

         IJK = FUNIJK(I, J, K)
         IF(.NOT.FLUID_AT(IJK)) CYCLE

         IF(CARTESIAN_GRID) THEN
            CALL CHECK_IF_PARCEL_OVERLAPS_STL(POS, SKIP)
            IF(SKIP) CYCLE
         ENDIF

         PIP = PIP + 1
         CALL PARTICLE_GROW(PIP)

         DES_POS_NEW(:,PIP) = POS(:)
         DES_VEL_NEW(:,PIP) = VEL(:)

         DES_RADIUS(PIP) = DES_D_P0(M)*HALF
         RO_SOL(PIP) =  DES_RO_S(M)

         PIJK(PIP,1) = I
         PIJK(PIP,2) = J
         PIJK(PIP,3) = K
         PIJK(PIP,4) = IJK
         PIJK(PIP,5) = M

         DES_STAT_WT(PIP) = STAT_WT

         CALL SET_NORMAL(PIP)

         SEEDED = SEEDED + 1

         IF(sVOL_TOT <= sVOL*dble(SEEDED)) EXIT
      ENDDO

      sDATA(1) = dble(SEEDED)
      sDATA(2) = STAT_WT
      sDATA(3) = EP_SM
      sDATA(4) = sVOL*dble(SEEDED)

      IF(allocated(randPOS)) deallocate(randPOS)
      IF(allocated(randPOS)) deallocate(randVEL)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE GPC_MPPIC_CONST_STATWT



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GENERATE_PARTICLE_CONFIG_MMPPIC                         !
!  Author: Rahul Garg                                 Date: 3-May-2011 !
!                                                                      !
!  Purpose: generates particle position distribution for MPPIC         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GPC_MPPIC_CONST_NPC(ICV, M, IC_VOL, sDATA)

! Global variables
!---------------------------------------------------------------------//
      use cutcell, only : CARTESIAN_GRID
      use discretelement, only: XE, YN, ZT

! Number of DES solids phases.
      use discretelement, only: DES_MMAX

! DEM solid phase diameters and densities.
      use discretelement, only: DES_D_p0, DES_RO_s

! Implied total number of physical particles
      use mfix_pic, only: rnp_pic
! Total number of computational particles/parcels
      use mfix_pic, only: cnp_pic

      use mfix_pic, only: des_stat_wt

! Flag indicating that the IC region is defined.
      use ic, only: IC_DEFINED
! IC Region bulk density (RO_s * EP_s)
      use ic, only: IC_EP_s

      use ic, only: IC_I_w, IC_I_e, IC_J_s, IC_J_n, IC_K_b, IC_K_t

! initally specified velocity field and granular temperature
      use ic, only: IC_U_s, IC_V_s, IC_W_s, IC_Theta_M
      use ic, only: IC_PIC_CONST_NPC

! Cut_cell identifier array
      use cutcell, only: cut_cell_at

! Maximum number of IC regions and solids phases
      use param, only: DIMENSION_IC
      use param, only: DIM_M
      use param1, only: UNDEFINED, UNDEFINED_I, ZERO, ONE, HALF

! Constant: 3.14159...
      use constant, only: PI

      use mpi_utility

      use randomno
      use error_manager
      use functions
!      use STL_PREPROC_DES, only: CHECK_IF_PARTICLE_OVERLAPS_STL
      use run, only: solids_model
      use des_allocate, only: PARTICLE_GROW

      IMPLICIT NONE

! Dummy Arguments
!----------------------------------------------------------------------//
! Index of IC region and solids phase
      INTEGER, INTENT(IN) :: ICV, M
! Specific volume of IC region (accounts for blocked cells)
      DOUBLE PRECISION, INTENT(IN) :: IC_VOL
! Data about solids in the IC region.
      DOUBLE PRECISION, INTENT(OUT) :: sDATA(4)

! Local variables
!----------------------------------------------------------------------//
      DOUBLE PRECISION :: EP_SM


! Number of real and comp. particles in a cell.
      DOUBLE PRECISION ::  rPARTS
      INTEGER :: cPARTS
      DOUBLE PRECISION :: DOML(3), IC_START(3)
! Parcel position with random
      DOUBLE PRECISION :: POS(3)
! Average velocity and standard deivation
      DOUBLE PRECISION :: IC_VEL(3), VEL_SIG
! Flag to not keep parcel.
      LOGICAL :: SKIP
! Arrasy for assigning random position and velocities
      DOUBLE PRECISION, ALLOCATABLE :: randPOS(:,:), randVEL(:,:)
! Statistical weights
      DOUBLE PRECISION :: STAT_WT, SUM_STAT_WT
! Volume of a parcel and total solids volume
      DOUBLE PRECISION :: sVOL, sVOL_TOT
! Counter for seeded parcles.
      INTEGER :: SEEDED
! Generic loop indices and loop counters
      INTEGER :: I, J, K, IJK, LC, LC_MAX
!......................................................................!

      CALL INIT_ERR_MSG("GPC_MPPIC_CONST_NPC")

      cPARTS = IC_PIC_CONST_NPC(ICV,M)

      allocate(randPOS(cPARTS,3))
      allocate(randVEL(cPARTS,3))

      IC_VEL(1) = IC_U_s(ICV,M)
      IC_VEL(2) = IC_V_s(ICV,M)
      IC_VEL(3) = merge(IC_W_s(ICV,M),0.0d0,DO_K)

      VEL_SIG = sqrt(IC_Theta_M(ICV,M))

! Volume occupied by one particle
      sVOL = (Pi/6.0d0)*(DES_D_P0(M)**3.d0)

      SEEDED = 0
      sVOL_TOT = 0.0d0
      SUM_STAT_WT = 0.0d0

      DO K = IC_K_B(ICV), IC_K_T(ICV)
      DO J = IC_J_S(ICV), IC_J_N(ICV)
      DO I = IC_I_W(ICV), IC_I_E(ICV)

         IF(.not.IS_ON_myPE_wobnd(I,J,K)) cycle
         IF(DEAD_CELL_AT(I,J,K)) cycle

         IJK = FUNIJK(I,J,K)
         IF(.not.FLUID_AT(IJK)) cycle

         rPARTS = IC_EP_s(ICV,M)*VOL(IJK)/sVOL

         LC_MAX = cPARTS
         IF(CUT_CELL_AT(IJK)) LC_MAX = &
            max(1, int(VOL(IJK)*dble(cPARTS)/(DX(I)*DY(J)*DZ(K))))

         DO LC=1, merge(2,3,NO_K)
            CALL UNI_RNO(RANDPOS(1:LC_MAX,LC))
            IF(VEL_SIG > ZERO) THEN
               CALL NOR_RNO(randVEL(1:LC_MAX,LC), IC_VEL(LC), VEL_SIG)
            ELSE
               randVEL(1:LC_MAX,LC) = IC_VEL(LC)
            ENDIF
         ENDDO
         IF(NO_K) THEN
            randPOS(1:LC_MAX,3) = 0.0d0
            randVEL(1:LC_MAX,3) = 0.0d0
         ENDIF

         IC_START(1) = XE(I-1)
         IC_START(2) = YN(J-1)
         IC_START(3) = ZERO;  IF(DO_K) IC_START(3) = ZT(K-1)

         DOML(1) = DX(I)
         DOML(2) = DY(J)
         DOML(3) = ZERO;  IF(DO_K) DOML(3) = DZ(K)

         DO LC=1,LC_MAX

            POS(:) = IC_START + DOML*randPOS(LC,:)

            IF(CARTESIAN_GRID) THEN
               CALL CHECK_IF_PARCEL_OVERLAPS_STL(POS, SKIP)
               IF(SKIP) CYCLE
            ENDIF

            PIP = PIP + 1
            CALL PARTICLE_GROW(PIP)
 
            DES_POS_NEW(:,PIP) = POS(:)
            DES_VEL_NEW(:,PIP) = randVEL(LC,:)
 
            DES_RADIUS(PIP) = DES_D_P0(M)*HALF
            RO_SOL(PIP) =  DES_RO_S(M)
 
            PIJK(PIP,1) = I
            PIJK(PIP,2) = J
            PIJK(PIP,3) = K
            PIJK(PIP,4) = IJK
            PIJK(PIP,5) = M
 
            STAT_WT = rPARTS/dble(LC_MAX)
            SUM_STAT_WT = SUM_STAT_WT + STAT_WT

            DES_STAT_WT(PIP) = STAT_WT
            sVOL_TOT = sVOL_TOT + sVOL*STAT_WT
 
            CALL SET_NORMAL(PIP)
 
            SEEDED = SEEDED + 1

         ENDDO

      ENDDO
      ENDDO
      ENDDO

      IF(allocated(randPOS)) deallocate(randPOS)
      IF(allocated(randPOS)) deallocate(randVEL)

      sDATA(1) = dble(SEEDED)
      sDATA(2) = SUM_STAT_WT/dble(SEEDED)
      sDATA(3) = IC_EP_S(ICV,M)
      sDATA(4) = sVOL_TOT

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE GPC_MPPIC_CONST_NPC

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GET_IC_VOLUME                                           !
!  Author: J.Musser                                 Date: 26-Aug-2015  !
!                                                                      !
!  Purpose: Calculate the actual volume of the IC region.              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GET_IC_VOLUME(ICV, IC_VOL)

! IC region index bounds
      use ic, only: IC_I_w, IC_I_e
      use ic, only: IC_J_s, IC_J_n
      use ic, only: IC_K_b, IC_K_t

! Volume of computational cells.
      use geometry, only: VOL

      use functions

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! Index of IC region
      INTEGER, INTENT(IN) :: ICV
! Total calculated volume of IC region
      DOUBLE PRECISION, INTENT(OUT) :: IC_VOL

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: I, J, K, IJK
!......................................................................!


      IC_VOL = 0.0d0
      DO K = IC_K_B(ICV), IC_K_T(ICV)
      DO J = IC_J_S(ICV), IC_J_N(ICV)
      DO I = IC_I_W(ICV), IC_I_E(ICV)

         IF(.NOT.IS_ON_MYPE_WOBND(I,J,K)) CYCLE
         IF(DEAD_CELL_AT(I,J,K)) CYCLE

         IJK = FUNIJK(I,J,K)
         IF(FLUID_AT(IJK)) IC_VOL = IC_VOL + VOL(IJK)

      ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE GET_IC_VOLUME

      END MODULE GENERATE_PARTICLES
