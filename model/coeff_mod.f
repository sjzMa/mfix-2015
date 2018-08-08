!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  MODULE: COEFF                                                       !
!  Purpose: Contains logic flags that tells the code whether to        !
!           perform the indicated type of calculation when the         !
!           value is true                                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE coeff

! Flags used by PHYSICAL_PROP :: (0:DIMENSION_M)
!```````````````````````````````````````````````````````````````````````
      LOGICAL, ALLOCATABLE :: DENSITY(:)  ! Density
      LOGICAL, ALLOCATABLE :: SP_HEAT(:)  ! Specific heat
      LOGICAL, ALLOCATABLE :: PSIZE(:)    ! Particle diameter


! Flags used by TRANSPORT_PROP :: (0:DIMENSION_M)
!```````````````````````````````````````````````````````````````````````
      LOGICAL, ALLOCATABLE :: VISC(:)      ! Viscosity
      LOGICAL, ALLOCATABLE :: COND(:)      ! Conductivity
      LOGICAL, ALLOCATABLE :: DIFF(:)      ! Diffusivity
      LOGICAL, ALLOCATABLE :: GRAN_DISS(:) ! Granular energy dissipation


! Flags used by EXCHANGE :: (0:DIMENSION_M)x(0:DIMENSION_M)
!```````````````````````````````````````````````````````````````````````
      LOGICAL, ALLOCATABLE :: DRAGCOEF(:,:) ! Drag coefficient
      LOGICAL, ALLOCATABLE :: HEAT_TR(:,:)  ! Heat transfer coeff


      contains

!**********************************************************************!
!  SUBROUTINE: INIT_COEFF                                              !
!                                                                      !
!  Purpose: Initialize logical flags.                                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE INIT_COEFF(IER)

! Global Variables:
!-----------------------------------------------------------------------
      use param, only: DIMENSION_M
      use param1, only: UNDEFINED
! Kinetic theory model.
      USE run, only: kt_type_enum
      USE run, only: gd_1999
      USE run, only: gtsh_2012
      USE run, only: ia_2005
! Run-time flag for invoking DQMOM
      use run, only: CALL_DQMOM
! Real number of solids phases (GHD theory)
      use physprop, only: SMAX
! Run-time flag for to solve energy equations.
      use run, only: ENERGY_EQ
! Run-time flag for to solve species equations.
      use run, only: SPECIES_EQ
! Flag to recalculate gas viscosity.
      use visc_g, only: RECALC_VISC_G
! Run-time flag for invoking discrete element model
      use discretelement, only: DISCRETE_ELEMENT
! Run-time flag for gas/DEM coupling
      use discretelement, only: DES_CONTINUUM_COUPLED
! Run-time flag for invoking TFM/DEM hybrid model
      use discretelement, only: DES_CONTINUUM_HYBRID
! Run-time flag invoking QMOM theory
      use qmom_kinetic_equation, only: QMOMK
! Specified constant gas phase density (incompressible)
      use physprop, only: RO_G0
! Specified constant specific heat.
      use physprop, only: C_PG0, C_PS0
! Specified constant thermal conductivity.
      use physprop, only: K_G0, K_S0
! Specified number of solids phases.
      use physprop, only: MMAX
! Specified constant viscosity.
      use physprop, only: MU_g0
! Variable solids density flag.
      use run, only: SOLVE_ROs
! UDF flags for physical properties
      use run, only: USR_ROg, USR_ROs, USR_CPg, USR_CPs
! MMS flag
      use mms, only: USE_MMS
      implicit none

! Dummy Arguments:
!-----------------------------------------------------------------------
! Error flag.
      INTEGER, intent(inout) :: IER


! Local Variables.
!-----------------------------------------------------------------------
! Invoke debug routine:
      LOGICAL, parameter :: dbg_coeffs = .FALSE.
! Loop counter for solids phases
      INTEGER :: M

! Allocate and initialize:
!```````````````````````````````````````````````````````````````````````
      IF(.NOT.allocated(DENSITY)) allocate( DENSITY(0:DIMENSION_M))
      IF(.NOT.allocated(SP_HEAT)) allocate( SP_HEAT(0:DIMENSION_M))
      IF(.NOT.allocated(PSIZE)) allocate( PSIZE(0:DIMENSION_M))

      DENSITY = .FALSE.
      SP_HEAT = .FALSE.
      PSIZE   = .FALSE.

      IF(.NOT.allocated(VISC)) allocate( VISC(0:DIMENSION_M))
      IF(.NOT.allocated(COND)) allocate( COND(0:DIMENSION_M))
      IF(.NOT.allocated(DIFF)) allocate( DIFF(0:DIMENSION_M))
      IF(.NOT.allocated(GRAN_DISS)) allocate( GRAN_DISS(0:DIMENSION_M))

      VISC = .FALSE.
      COND = .FALSE.
      DIFF = .FALSE.
      GRAN_DISS = .FALSE.

      IF(.NOT.allocated(DRAGCOEF)) &
         allocate( DRAGCOEF(0:DIMENSION_M,0:DIMENSION_M))
      IF(.NOT.allocated(HEAT_TR)) &
         allocate( HEAT_TR(0:DIMENSION_M,0:DIMENSION_M))

      DRAGCOEF = .FALSE.
      HEAT_TR = .FALSE.

! Coefficients for gas phase parameters.
!```````````````````````````````````````````````````````````````````````
! Compressible flow.
      if(RO_G0 == UNDEFINED .OR. USR_ROg) DENSITY(0) = .TRUE.
! Viscosity is recalculated iteration-to-iteration if:
! 1) the energy equations are solved
! 2) a turbulace length scale is defined (L_SCALE0 /= ZERO)
! 3) K-Epsilon model is used.
      VISC(0) = RECALC_VISC_G
! Specific heat and thermal conductivity.
      if(ENERGY_EQ) then
         if(C_PG0 == UNDEFINED) SP_HEAT(0) = .TRUE.
         if(K_G0  == UNDEFINED) COND(0) = .TRUE.
      endif
      if(USR_CPg) SP_HEAT(0) = .TRUE.
! Species diffusivity.
      if(SPECIES_EQ(0)) DIFF(0) = .TRUE.


! Interphase transfer terms.
!```````````````````````````````````````````````````````````````````````
       if(.NOT.QMOMK .AND. .NOT.USE_MMS) DRAGCOEF(0:MMAX,0:MMAX)=.TRUE.

! Coefficients for solids phase parameters.
!```````````````````````````````````````````````````````````````````````
      IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID) THEN
! Interphase heat transfer coefficient (GAMA)
         if(ENERGY_EQ .AND. .NOT.USE_MMS) HEAT_TR(0:MMAX,0:MMAX)=.TRUE.

! Variable solids density.
         if(any(SOLVE_ROs)) DENSITY(1:MMAX) = .TRUE.
         if(USR_ROs) DENSITY(1:MMAX) = .TRUE.

! Solids viscosity.
!         DO M = 1, MMAX
!            IF (MU_s0(M) == UNDEFINED) THEN
!               VISC(M) = .TRUE.
!            ENDIF
!         ENDDO
! Calc_mu_s must be invoked every iteration even if constant viscosity
! (mu_s0 /= undefined) to incorporate ishii form of governing equations
! wherein the viscosity is multiplied by the phase volume fraction
         VISC(1:MMAX) = .TRUE.

! Specific heat and thermal conductivity.
         do M=1,MMAX
            if(ENERGY_EQ) THEN
               if(C_PS0(M) == UNDEFINED) SP_HEAT(M) = .TRUE.
               if(K_S0(M)  == UNDEFINED) COND(M) = .TRUE.
            endif
            if(USR_CPS) SP_HEAT(M) = .TRUE.
         enddo

! Species diffusivity. There is no reason to invoke this routine as the
! diffusion coefficient for solids is always zero.
         DIFF(1:MMAX) = .FALSE.

! Particle-Particle Energy Dissipation
         IF (KT_TYPE_ENUM == IA_2005 .OR. &
             KT_TYPE_ENUM == GD_1999 .OR. &
             KT_TYPE_ENUM == GTSH_2012) THEN
            GRAN_DISS(:MMAX) = .TRUE.
         ENDIF

! Particle diameter.
         if(Call_DQMOM) PSIZE(1:SMAX)=.TRUE.

      ENDIF   ! end if (.not.discrete_element .or des_continuum_hybrid)

      if(dbg_coeffs) CALL DEBUG_COEFF

! Invoke calc_coeff.
      IF(.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_COUPLED) THEN
! If gas viscosity is undefined and the flag for calculating gas
! viscosity is turned off: Turn it on and make the call to calc_coeff.
! Once viscosity values have been calculated (i.e., an initial value
! is calculated), turn the flag off again so it isn't recalculated.
         IF(MU_g0 == UNDEFINED .AND. .NOT.VISC(0)) THEN
            VISC(0) = .TRUE.; CALL CALC_COEFF(IER, 2)
            VISC(0) = .FALSE.
         ELSE
            CALL CALC_COEFF(IER, 2)
         ENDIF
      ENDIF

      END SUBROUTINE INIT_COEFF

!**********************************************************************!
!  SUBROUTINE: DEBUG_COEFF                                             !
!                                                                      !
!  Purpose: Dump the coefficient arrays for debugging.                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DEBUG_COEFF

      use compar
      use physprop, only: MMAX

      implicit none

      INTEGER :: M, MM

      if(myPE /= PE_IO) return

      write(*,"(/3x,'From DEBUG_COEFF:')")

      write(*,"(/3x,'Gas phase coefficients:')")
      write(*,"( 5x,'Density (RO_g):',1x,1L1)") DENSITY(0)
      write(*,"( 5x,'Specific heat (C_pg):',1x,1L1)") SP_HEAT(0)
      write(*,"( 5x,'Viscosity: (MU_g)',1x,1L1)") VISC(0)
      write(*,"( 5x,'Thermal conductivity (K_g):',1x,1L1)") COND(0)
      write(*,"( 5x,'Species diffusivity: (DIF_G)',1x,1L1)") DIFF(0)


      DO M=1, MMAX
         write(*,"(/3x,'Solids ',I1,' phase coefficients:')") M
         write(*,"( 5x,'Density: (RO_s)',1x,1L1)") DENSITY(M)
         write(*,"( 5x,'Specific heat (C_ps):',1x,1L1)") SP_HEAT(M)
         write(*,"( 5x,'Viscosity (MU_s):',1x,1L1)") VISC(M)
         write(*,"( 5x,'Thermal conductivity (K_s):',1x,1L1)") COND(M)
         write(*,"( 5x,'Species diffusivity (DIF_s):',1x,1L1)") DIFF(M)
         write(*,"( 5x,'Gran. Dissipation (D_p):',1x,1L1)") GRAN_DISS(M)
         write(*,"( 5x,'Diameter (D_p):',1x,1L1)") PSIZE(M)
      ENDDO


      write(*,"(/3x,'Interphase drag:')")
      write(*,"( 5x,'ref')",ADVANCE="NO")
      DO M=0, MMAX
         write(*,"(2x,I3)",ADVANCE="NO")M
      ENDDO
      write(*,"('')")

      DO M=0, MMAX
         write(*,"( 5x,I3)",ADVANCE="NO") M
         DO MM=0, MMAX
            write(*,"(2x,L3)",ADVANCE="NO")DRAGCOEF(M, MM)
         ENDDO
         write(*,"('')")
      ENDDO

      write(*,"(/3x,'Interphase heat transfer:')")
      write(*,"( 5x,'ref')",ADVANCE="NO")
      DO M=0, MMAX
         write(*,"(2x,I3)",ADVANCE="NO")M
      ENDDO
      write(*,"('')")
      DO M=0, MMAX
         write(*,"( 5x,I3)",ADVANCE="NO") M
         DO MM=0, MMAX
            write(*,"(2x,L3)",ADVANCE="NO")HEAT_TR(M, MM)
         ENDDO
         write(*,"('')")
      ENDDO

      write(*,"(/3x,'DEBUG_COEFF - Exit',3/)")

      END SUBROUTINE DEBUG_COEFF

      END MODULE coeff
