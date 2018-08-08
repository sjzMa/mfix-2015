!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNEWVALUES                                            C
!
!  Purpose: DES - Calculate the new values of particle velocity,
!           position, angular velocity etc
!
!                                                                      C
!  Comments: Implements Eqns 1, 2, 3, 4 & 5  from the following paper:
!    Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical
!    simulation of plug glow of cohesionless particles in a
!    horizontal pipe", Powder technology, 71, 239-250, 1992
!
!  pradeep : changes for parallel processing
!          1. periodic boundaries might lie in different proc. so adjust
!             particle position for periodic removed
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFNEWVALUES

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE constant
      USE des_bc
      USE discretelement
      USE fldvar
      USE mfix_pic
      USE mpi_utility
      USE parallel
      USE param
      USE param1
      USE physprop
      use geometry, only: DO_K, NO_K
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER :: L
      DOUBLE PRECISION :: DD(3), NEIGHBOR_SEARCH_DIST
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.
      DOUBLE PRECISION :: OMEGA_MAG,OMEGA_UNIT(3),ROT_ANGLE
!-----------------------------------------------

      IF(MPPIC) THEN
         IF(MPPIC_SOLID_STRESS_SNIDER) THEN
            CALL CFNEWVALUES_MPPIC_SNIDER
         ELSE
! call the coloring function like approach
            CALL CFNEWVALUES_MPPIC
         ENDIF
         RETURN
      ENDIF

! Adams-Bashforth defaults to Euler for the first time step.
      IF(FIRST_PASS .AND. INTG_ADAMS_BASHFORTH) THEN
         DO L =1, MAX_PIP
            IF(IS_NONEXISTENT(L)) CYCLE                       ! Only real particles
            IF(IS_ENTERING(L).or.IS_ENTERING_GHOST(L)) CYCLE  ! Only non-entering
            IF(IS_GHOST(L)) CYCLE                             ! Skip ghost particles
            DES_ACC_OLD(:,L) = FC(:,L)/PMASS(L) + GRAV(:)
            ROT_ACC_OLD(:,L) = TOW(:,L)
         ENDDO
      ENDIF

!$omp parallel do if(max_pip .ge. 10000) default(none)                    &
!$omp shared(MAX_PIP,INTG_EULER,INTG_ADAMS_BASHFORTH,fc,tow,              &
!$omp       omega_new,omega_old,pmass,grav,des_vel_new,des_pos_new,       &
!$omp       des_vel_old,des_pos_old,dtsolid,omoi,des_acc_old,rot_acc_old, &
!$omp       ppos,neighbor_search_rad_ratio,des_radius,DO_OLD, iGlobal_ID, &
!$omp       particle_orientation, orientation) &
!$omp private(l,dd,neighbor_search_dist,rot_angle,omega_mag,omega_unit)   &
!$omp reduction(.or.:do_nsearch) schedule (auto)

      DO L = 1, MAX_PIP
! only process particles that exist
         IF(IS_NONEXISTENT(L)) CYCLE
! skip ghost particles
         IF(IS_GHOST(L).or.IS_ENTERING_GHOST(L).or.IS_EXITING_GHOST(L)) CYCLE

! If a particle is classified as new, then forces are ignored.
! Classification from new to existing is performed in routine
! des_check_new_particle.f
         IF(.NOT.IS_ENTERING(L) .AND. .NOT.IS_ENTERING_GHOST(L))THEN
            FC(:,L) = FC(:,L)/PMASS(L) + GRAV(:)
         ELSE
            FC(:,L) = ZERO
            TOW(:,L) = ZERO
         ENDIF

! Advance particle position, velocity
        IF (INTG_EULER) THEN
! first-order method
            DES_VEL_NEW(:,L) = DES_VEL_NEW(:,L) + FC(:,L)*DTSOLID
            DD(:) = DES_VEL_NEW(:,L)*DTSOLID
            DES_POS_NEW(:,L) = DES_POS_NEW(:,L) + DD(:)
            OMEGA_NEW(:,L)   = OMEGA_NEW(:,L) + TOW(:,L)*OMOI(L)*DTSOLID
         ELSEIF (INTG_ADAMS_BASHFORTH) THEN

! Second-order Adams-Bashforth/Trapezoidal scheme
            DES_VEL_NEW(:,L) = DES_VEL_OLD(:,L) + 0.5d0*&
               ( 3.d0*FC(:,L)-DES_ACC_OLD(:,L) )*DTSOLID
            OMEGA_NEW(:,L)   =  OMEGA_OLD(:,L) + 0.5d0*&
               ( 3.d0*TOW(:,L)*OMOI(L)-ROT_ACC_OLD(:,L) )*DTSOLID
            DD(:) = 0.5d0*( DES_VEL_OLD(:,L)+DES_VEL_NEW(:,L) )*DTSOLID
            DES_POS_NEW(:,L) = DES_POS_OLD(:,L) + DD(:)
            DES_ACC_OLD(:,L) = FC(:,L)
            ROT_ACC_OLD(:,L) = TOW(:,L)*OMOI(L)
         ENDIF

! Update particle orientation - Always first order
! When omega is non-zero, compute the rotation angle, and apply the
! Rodrigues' rotation formula

         IF(PARTICLE_ORIENTATION) THEN
            OMEGA_MAG = OMEGA_NEW(1,L)**2 +OMEGA_NEW(2,L)**2 + OMEGA_NEW(3,L)**2

            IF(OMEGA_MAG>ZERO) THEN
               OMEGA_MAG=DSQRT(OMEGA_MAG)
               OMEGA_UNIT(:) = OMEGA_NEW(:,L)/OMEGA_MAG
               ROT_ANGLE = OMEGA_MAG * DTSOLID

               ORIENTATION(:,L) = ORIENTATION(:,L)*DCOS(ROT_ANGLE) &
                                 + DES_CROSSPRDCT(OMEGA_UNIT,ORIENTATION(:,L))*DSIN(ROT_ANGLE) &
                                 + OMEGA_UNIT(:)*DOT_PRODUCT(OMEGA_UNIT,ORIENTATION(:,L))*(ONE-DCOS(ROT_ANGLE))
            ENDIF
         ENDIF

! Check if the particle has moved a distance greater than or equal to
! its radius during one solids time step. if so, call stop
         IF(dot_product(DD,DD).GE.DES_RADIUS(L)**2) THEN
            WRITE(*,1002) iGlobal_ID(L), sqrt(dot_product(DD,DD)), DES_RADIUS(L)
            IF (DO_OLD) WRITE(*,'(5X,A,3(ES17.9))') &
               'old particle pos = ', DES_POS_OLD(:,L)
            WRITE(*,'(5X,A,3(ES17.9))') &
               'new particle pos = ', DES_POS_NEW(:,L)
            WRITE(*,'(5X,A,3(ES17.9))')&
               'new particle vel = ', DES_VEL_NEW(:,L)
            WRITE(*,1003)
            STOP 1
         ENDIF

! Check if the particle has moved a distance greater than or equal to
! its radius since the last time a neighbor search was called. if so,
! make sure that neighbor is called in des_time_march
         IF(.NOT.DO_NSEARCH) THEN
            DD(:) = DES_POS_NEW(:,L) - PPOS(:,L)
            NEIGHBOR_SEARCH_DIST = NEIGHBOR_SEARCH_RAD_RATIO*&
               DES_RADIUS(L)
            IF(dot_product(DD,DD).GE.NEIGHBOR_SEARCH_DIST**2) DO_NSEARCH = .TRUE.
         ENDIF

! Reset total contact force and torque
         FC(:,L) = ZERO
         TOW(:,L) = ZERO

      ENDDO
!$omp end parallel do

      FIRST_PASS = .FALSE.


 1002 FORMAT(/1X,70('*')//&
         ' From: CFNEWVALUES -',/&
         ' Message: Particle ',I10, ' moved a distance ', ES17.9, &
         ' during a',/10X, 'single solids time step, which is ',&
         ' greater than',/10X,'its radius: ', ES17.9)
 1003 FORMAT(1X,70('*')/)

      RETURN

      contains

        include 'functions.inc'

      END SUBROUTINE CFNEWVALUES
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNEWVALUES_MPPIC_SNIDER                               C
!
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFNEWVALUES_MPPIC_SNIDER

      USE param
      USE param1
      USE parallel
      USE physprop
      USE constant
      USE fldvar
      USE discretelement
      USE des_bc
      USE mpi_utility
      USE fldvar
      USE cutcell
      USE mfix_pic
      USE randomno
      USE geometry, only: DO_K, NO_K
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER L, M, IDIM
      INTEGER I, J, K, IJK, IJK_OLD

      DOUBLE PRECISION DD(3), DIST, &
                       DP_BAR, COEFF_EN, MEANVEL(3), D_GRIDUNITS(3)

      DOUBLE PRECISION DELUP(3), UPRIMETAU(3), UPRIMETAU_INT(3), MEAN_FREE_PATH, PS_FORCE(3)
! index to track accounted for particles
      INTEGER PC

! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG

! maximum distance particles can move in MPPIC
      DOUBLE PRECISION UPRIMEMOD

! dt's in each direction  based on cfl_pic for the mppic case

      DOUBLE PRECISION DTPIC_TMPX, DTPIC_TMPY , DTPIC_TMPZ, THREEINTOSQRT2, RAD_EFF

      LOGICAL :: DELETE_PART, INSIDE_DOMAIN
      INTEGER :: PIP_DEL_COUNT

      INTEGER :: LPIP_DEL_COUNT_ALL(0:numPEs-1), PIJK_OLD(5)

      double precision  sig_u, mean_u
      double precision, allocatable, dimension(:,:) ::  rand_vel
!-----------------------------------------------

      PC = 1
      DTPIC_CFL = LARGE_NUMBER

      if(NO_K) THREEINTOSQRT2 = 2.D0*SQRT(2.D0)
      if(DO_K) THREEINTOSQRT2 = 3.D0*SQRT(2.D0)
      THREEINTOSQRT2 = 3.D0*SQRT(2.D0)
      DES_VEL_MAX(:) = ZERO
      PIP_DEL_COUNT = 0

      !EPG_MIN2 = MINVAL(EP_G(:))
      !epg_min_loc = MINLOC(EP_G)
      !IJK = epg_min_loc(1)

      allocate(rand_vel(3, MAX_PIP))
      do idim = 1, merge(2,3,NO_K)
         mean_u = zero
         sig_u = 1.d0
         CALL NOR_RNO(RAND_VEL(IDIM, 1:MAX_PIP), MEAN_U, SIG_U)
      enddo

      DO L = 1, MAX_PIP
         DELETE_PART = .false.
! pradeep skip ghost particles
         if(pc.gt.pip) exit
         if(is_nonexistent(l)) cycle
         pc = pc+1
         if(is_ghost(l) .or. is_entering_ghost(l) .or. is_exiting_ghost(l)) cycle

         DES_LOC_DEBUG = .FALSE.

         !if(mppic) then
         !   IJK = PIJK(L, 4)

         !   COUNT_BC = DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC
         !   IF(COUNT_BC.GE.1) CALL MPPIC_ADD_FRIC_FORCE(L)
         !ENDIF
! If a particle is classified as new, then forces are ignored.
! Classification from new to existing is performed in routine
! des_check_new_particle.f
         IF(.NOT.IS_ENTERING(L) .AND. .NOT.IS_ENTERING_GHOST(L))THEN
            FC(:,L) = FC(:,L)/PMASS(L) + GRAV(:)
         ELSE
            FC(:,L) = ZERO
         ENDIF

         !DP_BAR is the D_p in the Snider paper (JCP, vol 170, 523-549, 2001)
         !By comparing MFIX equations and equations in those papers,
         !D_p = Beta/(EP_S*RHOP)
         !F_gp in drag_fgs.f  = Beta*PVOL/EP_S
         !Therefore, D_p = F_gp/(PVOL*RHOP) = F_gp/PMASS
         DP_BAR = F_gp(L)/(PMASS(L))
         IF(.NOT.DES_CONTINUUM_COUPLED) DP_BAR = ZERO

         if(.not.MPPIC_PDRAG_IMPLICIT) DP_BAR = ZERO

         M = PIJK(L,5)
         IJK = PIJK(L,4)

         IJK_OLD = IJK

         PIJK_OLD(:) = PIJK(L,:)

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         COEFF_EN =  MPPIC_COEFF_EN1
         UPRIMETAU(:) = ZERO

         DES_VEL_NEW(:,L) = (DES_VEL_NEW(:,L) + &
         & FC(:,L)*DTSOLID)/(1.d0+DP_BAR*DTSOLID)

         do idim = 1, merge(2,3,NO_K)
            SIG_U = 0.05D0
            rand_vel(idim, L)  = sig_u*DES_VEL_NEW(IDIM,L)*rand_vel(idim, L)
         enddo

         !MEANVEL(1) = DES_U_S(IJK_OLD,M)
         !MEANVEL(2) = DES_V_S(IJK_OLD,M)
         !IF(DO_K) MEANVEL(3) = DES_W_S(IJK_OLD,M)

         PS_FORCE(:) = PS_GRAD(L, :)
         DELUP(:) = -( DTSOLID*PS_FORCE(:))/((1.d0+DP_BAR*DTSOLID))
         DELUP(:) = DELUP(:)/ROP_S(IJK_OLD,M)

         MEANVEL(:) = AVGSOLVEL_P(:,L)

         DO IDIM = 1, merge(2,3,NO_K)
            IF(PS_FORCE(IDIM).LE.ZERO) THEN
               UPRIMETAU(IDIM) = MIN(DELUP(IDIM), (1+COEFF_EN)*(MEANVEL(IDIM)-DES_VEL_NEW(IDIM,L)))
               UPRIMETAU_INT(IDIM) = UPRIMETAU(IDIM)
               UPRIMETAU(IDIM) = MAX(UPRIMETAU(IDIM), ZERO)
            ELSE
               UPRIMETAU(IDIM) = MAX(DELUP(IDIM), (1+COEFF_EN)*(MEANVEL(IDIM)-DES_VEL_NEW(IDIM,L)))
               UPRIMETAU_INT(IDIM) = UPRIMETAU(IDIM)
               UPRIMETAU(IDIM) = MIN(UPRIMETAU(IDIM), ZERO)
            END IF

         ENDDO

         DES_VEL_NEW(:,L) = DES_VEL_NEW(:,L) + UPRIMETAU(:)
         DD(:) = DES_VEL_NEW(:,L)*DTSOLID

         UPRIMEMOD = SQRT(DOT_PRODUCT(DES_VEL_NEW(1:,L), DES_VEL_NEW(1:,L)))

         RAD_EFF = DES_RADIUS(L)
               !RAD_EFF = (DES_STAT_WT(L)**(1.d0/3.d0))*DES_RADIUS(L)
         MEAN_FREE_PATH  = MAX(1.d0/(1.d0-EP_STAR), 1.d0/(1.D0-EP_G(IJK_OLD)))
         MEAN_FREE_PATH = MEAN_FREE_PATH*RAD_EFF/THREEINTOSQRT2

         !IF(UPRIMEMOD*DTSOLID.GT.MEAN_FREE_PATH) then
         !   DES_VEL_NEW(:,L) = (DES_VEL_NEW(:,L)/UPRIMEMOD)*MEAN_FREE_PATH/DTSOLID
         !ENDIF

          DIST = SQRT(DOT_PRODUCT(DD,DD))

         CALL PIC_FIND_NEW_CELL(L)

         IJK = PIJK(L,4)

         INSIDE_DOMAIN = .true.

         INSIDE_DOMAIN = FLUID_AT(IJK)

         IF(EP_G(IJK).LT.EP_STAR.and.INSIDE_DOMAIN.and.IJK.NE.IJK_OLD) then
            DD(:) = rand_vel(:, L)*dtsolid
            DES_VEL_NEW(:,L) = 0.8d0*DES_VEL_NEW(:,L)
         ENDIF



         PIJK(L,:) = PIJK_OLD(:)
         DIST = SQRT(DOT_PRODUCT(DD,DD))

          IF(DIST.GT.MEAN_FREE_PATH) THEN
             DD(:) =  DES_VEL_NEW(:,L)*DTSOLID*MEAN_FREE_PATH/DIST
          ENDIF

         DES_POS_NEW(:,L) = DES_POS_NEW(:,L) + DD(:)

         D_GRIDUNITS(1) = ABS(DD(1))/DX(PIJK(L,1))
         D_GRIDUNITS(2) = ABS(DD(2))/DY(PIJK(L,2))
         D_GRIDUNITS(3) = 1.d0
         IF(DO_K) D_GRIDUNITS(3) = ABS(DD(3))/DZ(PIJK(L,3))

         DIST = SQRT(DOT_PRODUCT(DD,DD))
         DTPIC_TMPX = (CFL_PIC*DX(PIJK(L,1)))/(ABS(DES_VEL_NEW(1,L))+SMALL_NUMBER)
         DTPIC_TMPY = (CFL_PIC*DY(PIJK(L,2)))/(ABS(DES_VEL_NEW(2,L))+SMALL_NUMBER)
         DTPIC_TMPZ = LARGE_NUMBER
         IF(DO_K) DTPIC_TMPZ = (CFL_PIC*DZ(PIJK(L,3)))/(ABS(DES_VEL_NEW(3,L))+SMALL_NUMBER)


! Check if the particle has moved a distance greater than or equal to grid spacing
! if so, then exit

         DO IDIM = 1, merge(2,3,NO_K)
            IF(D_GRIDUNITS(IDIM).GT.ONE) THEN
               IF(DMP_LOG.OR.myPE.eq.pe_IO) THEN

                  WRITE(UNIT_LOG, 2001) L, D_GRIDUNITS(:), DES_VEL_NEW(:,L)
                  WRITE(*, 2001) L, D_GRIDUNITS(:), DES_VEL_NEW(:,L)
                  DELETE_PART = .true.

               ENDIF
               !CALL mfix_exit(myPE)
            END IF

         END DO
         IF(.not.DELETE_PART) DTPIC_CFL = MIN(DTPIC_CFL, DTPIC_TMPX, DTPIC_TMPY, DTPIC_TMPZ)
         FC(:,L) = ZERO

         IF(DELETE_PART) THEN
            CALL SET_NONEXISTENT(l)
            PIP_DEL_COUNT = PIP_DEL_COUNT + 1
         ENDIF
         IF (DES_LOC_DEBUG) WRITE(*,1001)
      ENDDO


      IF(MPPIC) THEN
         CALL global_all_max(DTPIC_CFL)
         PIP = PIP - PIP_DEL_COUNT

         LPIP_DEL_COUNT_ALL(:) = 0
         LPIP_DEL_COUNT_ALL(mype) = PIP_DEL_COUNT
         CALL GLOBAL_ALL_SUM(LPIP_DEL_COUNT_ALL)
         IF((DMP_LOG).AND.SUM(LPIP_DEL_COUNT_ALL(:)).GT.0) THEN
            IF(PRINT_DES_SCREEN) THEN
               WRITE(*,'(/,2x,A,2x,i10,/,A)') &
                    'TOTAL NUMBER OF PARTS STEPPING MORE THAN ONE GRID SPACE = ', SUM(LPIP_DEL_COUNT_ALL(:)), &
                    'THIS SHOULD NOT HAPPEN FREQUENTLY: MONITOR THIS MESSAGE'
            ENDIF

            WRITE(UNIT_LOG,'(/,2x,A,2x,i10,/,A)') &
                 'TOTAL NUMBER OF PARTS  STEPPING MORE THAN ONE GRID SPACEC= ', SUM(LPIP_DEL_COUNT_ALL(:)), &
                 'THIS SHOULD NOT HAPPEN FREQUENTLY: MONITOR THIS MESSAGE'
            !DO IPROC = 0, NUMPES-1
            !   WRITE(UNIT_LOG, '(/,A,i4,2x,A,2x,i5)') &
            !   'PARTICLES OUTSIDE DOMAIN (PIC)  ON PROC:', IPROC,' EQUAL TO', LPIP_DEL_COUNT_ALL(IPROC)
            !ENDDO

         ENDIF

         DTPIC_MAX = MIN(DTPIC_CFL, DTPIC_TAUP)

         IF(DTSOLID.GT.DTPIC_MAX) THEN
            !IF(DMP_LOG) WRITE(UNIT_LOG, 2001) DTSOLID
            IF(myPE.eq.pe_IO) WRITE(*, 2004) DTSOLID
         ELSEIF(DTSOLID.LT.DTPIC_MAX) THEN

            !IF(DMP_LOG) WRITE(UNIT_LOG, 2002) DTSOLID
            IF(myPE.eq.pe_IO) WRITE(*, 2002) DTPIC_MAX
         ELSE
            !WRITE(*,'(A40,2x,g17.8)') 'DT
            !IF(mype.eq.pe_IO) WRITE(*,2003) DTSOLID
         END IF


      ENDIF
 1001 FORMAT(3X,'<---------- END CFNEWVALUES ----------')

2001  FORMAT(/1X,70('*'),//,10X,  &
           & 'MOVEMENT UNDESIRED IN CFNEWVALUES: PARTICLE', i5, /,10X, &
           & ' MOVED MORE THAN A GRID SPACING IN ONE TIME STEP', /,10X, &
           & 'MOVEMENT IN GRID UNITS = ', 3(2x, g17.8),/,10X,  &
           & 'TERMINAL ERROR: NOT STOPPING, BUT DELETING THE PARTICLE', &
           & /1X,70('*'), /10X, &
           & 'DES_VEL_NEW = ',  3(2x, g17.8))

 2004 FORMAT(/10x, &
      & 'DTSOLID SHUD BE REDUCED TO', g17.8)

 2002 FORMAT(/10x, &
      & 'DTSOLID CAN BE INCREASED TO', g17.8)

      RETURN
      END SUBROUTINE CFNEWVALUES_MPPIC_SNIDER


      SUBROUTINE CFNEWVALUES_MPPIC

      USE param
      USE param1
      USE parallel
      USE physprop
      USE constant
      USE fldvar
      USE discretelement
      USE des_bc
      USE mpi_utility
      USE mfix_pic
      USE randomno
      USE cutcell
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER L, M, IDIM
      INTEGER I, J, K, IJK, IJK_OLD

      DOUBLE PRECISION DD(3), DIST, &
                       DP_BAR, MEANVEL(3), D_GRIDUNITS(3)

      DOUBLE PRECISION MEAN_FREE_PATH
! index to track accounted for particles
      INTEGER PC

! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG

! maximum distance particles can move in MPPIC
      DOUBLE PRECISION UPRIMEMOD

! dt's in each direction  based on cfl_pic for the mppic case

      DOUBLE PRECISION DTPIC_TMPX, DTPIC_TMPY , DTPIC_TMPZ, THREEINTOSQRT2, RAD_EFF, POS_Z
      DOUBLE PRECISION :: VELP_INT(3)

      LOGICAL :: DELETE_PART, INSIDE_DOMAIN
      INTEGER :: PIP_DEL_COUNT

      INTEGER :: LPIP_DEL_COUNT_ALL(0:numPEs-1), PIJK_OLD(5)

      double precision  sig_u, mean_u, DTPIC_MIN_X,  DTPIC_MIN_Y,  DTPIC_MIN_Z
      double precision, allocatable, dimension(:,:) ::  rand_vel

      double precision :: norm1, norm2, norm3
      Logical :: OUTER_STABILITY_COND, DES_FIXED_BED
!-----------------------------------------------

      OUTER_STABILITY_COND = .false.
      DES_FIXED_BED = .false.
      PC = 1
      DTPIC_CFL = LARGE_NUMBER
      DTPIC_MIN_X =  LARGE_NUMBER
      DTPIC_MIN_Y =  LARGE_NUMBER
      DTPIC_MIN_Z =  LARGE_NUMBER

      if(NO_K) THREEINTOSQRT2 = 2.D0*SQRT(2.D0)
      if(DO_K) THREEINTOSQRT2 = 3.D0*SQRT(2.D0)
      THREEINTOSQRT2 = 3.D0*SQRT(2.D0)
      DES_VEL_MAX(:) = ZERO
      PIP_DEL_COUNT = 0

      !EPG_MIN2 = MINVAL(EP_G(:))
      !epg_min_loc = MINLOC(EP_G)
      !IJK = epg_min_loc(1)
      allocate(rand_vel(3, MAX_PIP))
      do idim = 1, merge(2,3,NO_K)
         mean_u = zero
         sig_u = 1.d0
         CALL NOR_RNO(RAND_VEL(IDIM, 1:MAX_PIP), MEAN_U, SIG_U)
      enddo
      !WRITE(*, '(A,2x,3(g17.8))') 'FC MIN AND MAX IN Z = ', MINVAL(FC(:,3)), MAXVAL(FC(:,3))
      !WRITE(UNIT_LOG, '(A,2x,3(g17.8))') 'FC MIN AND MAX IN Z = ', MINVAL(FC(:,3)), MAXVAL(FC(:,3))

      DO L = 1, MAX_PIP
         DELETE_PART = .false.
! pradeep skip ghost particles
         if(pc.gt.pip) exit
         if(is_nonexistent(l)) cycle
         pc = pc+1
         if(is_ghost(l) .or. is_entering_ghost(l) .or. is_exiting_ghost(l)) cycle

         DES_LOC_DEBUG = .FALSE.

         IF(.NOT.IS_ENTERING(L))THEN
            FC(:,L) = FC(:,L)/PMASS(L) + GRAV(:)
         ELSE
            FC(:,L) = ZERO
         ENDIF

         IF(DES_FIXED_BED) FC(:,L) = ZERO

         !DP_BAR is the D_p in the Snider paper (JCP, vol 170, 523-549, 2001)
         !By comparing the MFIX and equations in those papers,
         !D_p = Beta/(EP_S*RHOP)
         !F_gp in drag_fgs.f  = Beta*PVOL/EP_S
         !Therefore, D_p = F_gp/(PVOL*RHOP) = F_gp/PMASS
         DP_BAR = F_gp(L)/(PMASS(L))
         IF(.NOT.DES_CONTINUUM_COUPLED) DP_BAR = ZERO

         if(.not.MPPIC_PDRAG_IMPLICIT) DP_BAR = ZERO
         M = PIJK(L,5)
         IJK = PIJK(L,4)
         IJK_OLD = IJK
         PIJK_OLD(:) = PIJK(L,:)

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         DES_VEL_NEW(:,L) = (DES_VEL_NEW(:,L) + &
         & FC(:,L)*DTSOLID)/(1.d0+DP_BAR*DTSOLID)

         IF(DES_FIXED_BED) DES_VEL_NEW(:,L) = ZERO

         !MPPIC_VPTAU(L,:) = DES_VEL_NEW(:,L)

         VELP_INT(:) = DES_VEL_NEW(:,L)

         MEANVEL(1) = DES_U_S(IJK_OLD,M)
         MEANVEL(2) = DES_V_S(IJK_OLD,M)
         IF(DO_K) MEANVEL(3) = DES_W_S(IJK_OLD,M)

         RAD_EFF = DES_RADIUS(L)
         !RAD_EFF = (DES_STAT_WT(L)**(1.d0/3.d0))*DES_RADIUS(L)
         IF(.not.DES_ONEWAY_COUPLED) then
            MEAN_FREE_PATH  = MAX(1.d0/(1.d0-EP_STAR), 1.d0/(1.D0-EP_G(IJK_OLD)))
            MEAN_FREE_PATH  = MEAN_FREE_PATH*RAD_EFF/THREEINTOSQRT2
         endif

         DO IDIM = 1, merge(2,3,NO_K)
            !SIG_U = 0.05D0*MEANVEL(IDIM)
            !DES_VEL_NEW(IDIM,L) = DES_VEL_NEW(IDIM,L) + SIG_U*RAND_VEL(IDIM, L )
            !PART_TAUP = RO_Sol(L)*((2.d0*DES_RADIUS(L))**2.d0)/(18.d0* MU_G(IJK))
            SIG_U = 0.005D0      !*MEANVEL(IDIM)
            RAND_VEL(IDIM, L)  = SIG_U*RAND_VEL(IDIM, L)*DES_VEL_NEW(IDIM,L)
            IF(DES_FIXED_BED) RAND_VEL(IDIM,L) = ZERO
            !rand_vel(idim, L)  = sig_u*rand_vel(idim, L)/part_taup
            !rand_vel(idim, L)  = sig_u* mean_free_path*rand_vel(idim, L)/part_taup
            DES_VEL_NEW(idim,L) = DES_VEL_NEW(idim,L) + rand_vel(idim, L)
         enddo

         IF(.not.DES_ONEWAY_COUPLED.and.(.not.des_fixed_bed)) CALL MPPIC_APPLY_PS_GRAD_PART(L)

         UPRIMEMOD = SQRT(DOT_PRODUCT(DES_VEL_NEW(1:,L), DES_VEL_NEW(1:,L)))

         !IF(UPRIMEMOD*DTSOLID.GT.MEAN_FREE_PATH) then
         !   DES_VEL_NEW(:,L) = (DES_VEL_NEW(:,L)/UPRIMEMOD)*MEAN_FREE_PATH/DTSOLID
         !ENDIF

         IF(DES_FIXED_BED) THEN
            DD(:) = ZERO
         ELSE
            DD(:) = DES_VEL_NEW(:,L)*DTSOLID !+ rand_vel(:, L)*dtsolid
         ENDIF

         CALL PIC_FIND_NEW_CELL(L)

         IJK = PIJK(L,4)

         !IF((EP_G(IJK).LT.0.35.and.fluid_at(ijk)).or.(ep_g(ijk_old).lt.0.35)) then !.and.(ijk.ne.ijk_old)) then
         !IF((EP_G(IJK).LT.EP_STAR.and.fluid_at(ijk)).and.(ijk.ne.ijk_old)) then
         INSIDE_DOMAIN = .true.
         INSIDE_DOMAIN = FLUID_AT(IJK)

         IF(CUT_CELL_AT(IJK)) THEN
            POS_Z = zero
            IF(DO_K) POS_Z = DES_POS_NEW(3,L)
            CALL GET_DEL_H_DES(IJK,'SCALAR', &
            & DES_POS_NEW(1,L),  DES_POS_NEW(2,L), &
            & POS_Z, &
            & DIST, NORM1, NORM2, NORM3, .true.)

            IF(DIST.LE.ZERO) INSIDE_DOMAIN = .false.
         ENDIF

         !IF((EP_G(IJK).LT.EP_STAR.and.INSIDE_DOMAIN)) then
         !IF(EP_G(IJK).LT.EP_STAR.and.INSIDE_DOMAIN) then
         IF(1.d0 - EP_G(IJK).GT. 1.3d0*(1.d0 - EP_STAR).and.INSIDE_DOMAIN.and.IJK.NE.IJK_OLD.and.OUTER_STABILITY_COND) then
            DD(:) = ZERO
            DES_VEL_NEW(:,L) = 0.8d0*DES_VEL_NEW(:,L)
         ENDIF

         PIJK(L,:) = PIJK_OLD(:)

         DES_POS_NEW(:,L) = DES_POS_NEW(:,L) + DD(:)

         DIST = SQRT(DOT_PRODUCT(DD,DD))

         D_GRIDUNITS(1) = ABS(DD(1))/DX(PIJK(L,1))
         D_GRIDUNITS(2) = ABS(DD(2))/DY(PIJK(L,2))
         D_GRIDUNITS(3) = 1.d0
         IF(DO_K) D_GRIDUNITS(3) = ABS(DD(3))/DZ(PIJK(L,3))

         DIST = SQRT(DOT_PRODUCT(DD,DD))
         DTPIC_TMPX = (CFL_PIC*DX(PIJK(L,1)))/(ABS(DES_VEL_NEW(1,L))+SMALL_NUMBER)
         DTPIC_TMPY = (CFL_PIC*DY(PIJK(L,2)))/(ABS(DES_VEL_NEW(2,L))+SMALL_NUMBER)
         DTPIC_TMPZ = LARGE_NUMBER
         IF(DO_K) DTPIC_TMPZ = (CFL_PIC*DZ(PIJK(L,3)))/(ABS(DES_VEL_NEW(3,L))+SMALL_NUMBER)


! Check if the particle has moved a distance greater than or equal to grid spacing
! if so, then exit

         DO IDIM = 1, merge(2,3,NO_K)
            IF(D_GRIDUNITS(IDIM).GT.ONE) THEN
               IF(DMP_LOG.OR.myPE.eq.pe_IO) THEN

                  WRITE(UNIT_LOG, 2001) L, D_GRIDUNITS(:), DES_VEL_NEW(:,L)
                  WRITE(UNIT_LOG, '(A,2x,3(g17.8))') 'rand_vel = ', rand_vel(:,L)
                  IF (DO_OLD) WRITE(UNIT_LOG, '(A,2x,3(g17.8))') 'des_pos_old = ', des_pos_old(:,l)
                  WRITE(UNIT_LOG, '(A,2x,3(g17.8))') 'des_pos_new = ', des_pos_new(:,L)
                  WRITE(UNIT_LOG, '(A,2x,3(g17.8))') 'FC          = ', FC(:,L)

                  WRITE(*, 2001) L, D_GRIDUNITS(:), DES_VEL_NEW(:,L)

                  WRITE(*, '(A,2x,3(g17.8))') 'rand_vel = ', rand_vel(:,L)
                  IF (DO_OLD) WRITE(*, '(A,2x,3(g17.8))') 'des_pos_old = ', des_pos_old(:,l)
                  WRITE(*, '(A,2x,3(g17.8))') 'des_pos_new = ', des_pos_new(:,L)
                  WRITE(*, '(A,2x,3(g17.8))') 'FC          = ', FC(:,L)
                  read(*,*)
                  DELETE_PART = .true.

               ENDIF
               !CALL mfix_exit(myPE)
            END IF

         END DO

         IF(.not.DELETE_PART) then
            DTPIC_CFL = MIN(DTPIC_CFL, DTPIC_TMPX, DTPIC_TMPY, DTPIC_TMPZ)
            DTPIC_MIN_X = MIN(DTPIC_MIN_X, DTPIC_TMPX)
            DTPIC_MIN_Y = MIN(DTPIC_MIN_Y, DTPIC_TMPY)
            DTPIC_MIN_Z = MIN(DTPIC_MIN_Z, DTPIC_TMPZ)
         ENDIF
         FC(:,L) = ZERO

         IF(DELETE_PART) THEN
            CALL SET_NONEXISTENT(L)
            PIP_DEL_COUNT = PIP_DEL_COUNT + 1
         ENDIF
         IF (DES_LOC_DEBUG) WRITE(*,1001)
      ENDDO


      DEALLOCATE(RAND_VEL)

      CALL global_all_max(DTPIC_CFL)
      PIP = PIP - PIP_DEL_COUNT

      LPIP_DEL_COUNT_ALL(:) = 0
      LPIP_DEL_COUNT_ALL(mype) = PIP_DEL_COUNT
      CALL GLOBAL_ALL_SUM(LPIP_DEL_COUNT_ALL)
      IF((DMP_LOG).AND.SUM(LPIP_DEL_COUNT_ALL(:)).GT.0) THEN
         IF(PRINT_DES_SCREEN) THEN
            WRITE(*,'(/,2x,A,2x,i10,/,A)') &
                 'TOTAL NUMBER OF PARTS STEPPING MORE THAN ONE GRID SPACE = ', SUM(LPIP_DEL_COUNT_ALL(:)), &
                 'THIS SHOULD NOT HAPPEN FREQUENTLY: MONITOR THIS MESSAGE'
         ENDIF

         WRITE(UNIT_LOG,'(/,2x,A,2x,i10,/,A)') &
              'TOTAL NUMBER OF PARTS  STEPPING MORE THAN ONE GRID SPACEC= ', SUM(LPIP_DEL_COUNT_ALL(:)), &
              'THIS SHOULD NOT HAPPEN FREQUENTLY: MONITOR THIS MESSAGE'
         !DO IPROC = 0, NUMPES-1
         !   WRITE(UNIT_LOG, '(/,A,i4,2x,A,2x,i5)') &
         !     'PARTICLES OUTSIDE DOMAIN (PIC)  ON PROC:', IPROC,' EQUAL TO', LPIP_DEL_COUNT_ALL(IPROC)
         !ENDDO

      ENDIF

      DTPIC_MAX = MIN(DTPIC_CFL, DTPIC_TAUP)

      RETURN
      IF(DTSOLID.GT.DTPIC_MAX) THEN
         !IF(DMP_LOG) WRITE(UNIT_LOG, 2001) DTSOLID
         IF(myPE.eq.pe_IO) WRITE(*, 2004) DTSOLID
      ELSEIF(DTSOLID.LT.DTPIC_MAX) THEN

         !IF(DMP_LOG) WRITE(UNIT_LOG, 2002) DTSOLID
         IF(myPE.eq.pe_IO) WRITE(*, 2002) DTPIC_MAX
      ELSE
        !WRITE(*,'(A40,2x,g17.8)') 'DT
        !IF(mype.eq.pe_IO) WRITE(*,2003) DTSOLID
      END IF

      WRITE(UNIT_LOG, '(10x,A,2x,3(g17.8))') 'DTPIC MINS IN EACH DIRECTION = ', DTPIC_MIN_X, DTPIC_MIN_Y, DTPIC_MIN_Z
      WRITE(*, '(10x,A,2x,3(g17.8))') 'DTPIC MINS IN EACH DIRECTION = ', DTPIC_MIN_X, DTPIC_MIN_Y, DTPIC_MIN_Z

 1001 FORMAT(3X,'<---------- END CFNEWVALUES ----------')

2001  FORMAT(/1X,70('*'),//,10X,  &
           & 'MOVEMENT UNDESIRED IN CFNEWVALUES: PARTICLE', i5, /,10X, &
           & ' MOVED MORE THAN A GRID SPACING IN ONE TIME STEP', /,10X, &
           & 'MOVEMENT IN GRID UNITS = ', 3(2x, g17.8),/,10X,  &
           & 'TERMINAL ERROR: NOT STOPPING, BUT DELETING THE PARTICLE', &
           & /1X,70('*'), /10X, &
           & 'DES_VEL_NEW = ',  3(2x, g17.8))

 2004 FORMAT(/10x, &
      & 'DTSOLID SHUD BE REDUCED TO', g17.8)

 2002 FORMAT(/10x, &
      & 'DTSOLID CAN BE INCREASED TO', g17.8)

      RETURN
      END SUBROUTINE CFNEWVALUES_MPPIC


!------------------------------------------------------------------------
! subroutine       : des_dbgpic
! Author           : Pradeep G.
! Purpose          : For debugging the pic values
! Parameters       : pstart - start indices of the particle
!                    pend - end indices of the particle
!                    pfreq - optional frequency (when the local count matches the
!                    frequency the filw will be written)
!                    if not send then it prints the file
!------------------------------------------------------------------------

      subroutine des_dbgpic (pstart,pend,pfreq)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use discretelement
      use fldvar
      use functions
      implicit none
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      INTEGER, INTENT(IN) :: pstart,pend
      INTEGER, INTENT(IN), OPTIONAL :: pfreq
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer lp,lijk
      integer, save :: lfcount = 0 ,lfreq =0
      character(255) :: filename
!-----------------------------------------------
      if (present(pfreq)) then
         lfreq = lfreq+1
         if (lfreq .ne. pfreq) return
         lfreq =0
      end if
      lfcount = lfcount + 1
      write(filename,'("debug",I3.3)') lfcount
      open (unit=100,file=filename)
      do lp = pstart,pend
         if (is_normal(lp) .or. is_entering(lp) .or. is_exiting(lp)) then
            lijk = pijk(lp,4)
            write(100,*)"positon =",lijk,pijk(lp,1),pijk(lp,2), &
               pijk(lp,3),ep_g(lijk),DES_U_s(lijk,1)
            write(100,*)"forces =", FC(2,lp),tow(1,lp)
         end if
      end do
      close (100)

      RETURN
      END SUBROUTINE des_dbgpic

!------------------------------------------------------------------------
! subroutine       : des_dbgtecplot
! Author           : Pradeep G.
! Purpose          : prints the tecplot file for particle location
! Parameters       : pstart - start indices of the particle
!                    pend - end indices of the particle
!                    pfreq - optional frequency (when the local count matches the
!                    frequency the filw will be written)
!                    if not send then it prints the file
!------------------------------------------------------------------------

      subroutine des_dbgtecplot (pstart,pend,pfreq)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement
      USE fldvar
      USE functions
      implicit none
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      INTEGER, INTENT(IN) :: pstart,pend
      INTEGER, INTENT(IN), OPTIONAL :: pfreq
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer lp,lijk
      integer, save :: lfcount = 0 ,lfreq =0
      character(255) :: filename
!-----------------------------------------------

      if (present(pfreq)) then
         lfreq = lfreq+1
         if (lfreq .ne. pfreq) return
         lfreq =0
      end if
      lfcount = lfcount + 1
      write(filename,'("new_tec",I3.3,".dat")') lfcount
      open (unit=100,file=filename)
      write(100,'(9(A,3X),A)') &
         'VARIABLES = ', '"ijk"', '"x"', '"y"', '"vx"', '"vy"', &
         '"ep_g"', '"FCX"' ,'"FCY"', '"TOW"'
      write(100,'(A,F14.7,A)') 'zone T = "' , s_time , '"'
      do lp = pstart,pend
         if (.not.is_nonexistent(lp)) then
            lijk = pijk(lp,4)
            write(100,*)lijk,des_pos_new(1,lp),des_pos_new(2,lp), &
               des_vel_new(1,lp),des_vel_new(2,lp),ep_g(lijk),&
               fc(1,lp),fc(2,lp),tow(lp,1)
         endif
      enddo
      close (100)
      RETURN
      END SUBROUTINE DES_DBGTECPLOT
