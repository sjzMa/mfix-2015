!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_THERMO                                             !
!                                                                      !
!  Purpose: Common elements for MFIX-DEM heat transfer.                !
!  condition.                                                          !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE DES_THERMO

      USE param

! Heat transfer correlation specified in mfix.dat
! Default [RANZ_1952]
      CHARACTER(LEN=24) :: DES_CONV_CORR

      INTEGER :: DES_CONV_CORR_ENUM
      INTEGER, PARAMETER :: RANZ_1952 = 0

! Particle properties
!--------------------------
! Particle temperature at current time step (S_TIME)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_T_s_OLD !(PARTICLES)
! Particle temperature at previous time step (S_TIME - DT_SOLID)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_T_s_NEW !(PARTICLES)

! DES specific heat of particles by particle
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_C_ps

! Emissivity of particles
      DOUBLE PRECISION DES_Em(DIM_M)
! Stefan-Boltzmann Constant
      DOUBLE PRECISION SB_CONST
! Bulk solids temperature for radiative heat transfer
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: avgDES_T_s


! Rates of heat transfer
!--------------------------
! Internal heat generation resulting from a chemical reaction
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Qint !(PARTICLES)
! Heat transfer source to the particle.
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Q_Source
! Previous time step's rate of heat transfer. Used for Adams-Bashforth
! time integration scheme
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Q_Source0

! Thermodynamic Neighborhood
!--------------------------
! Fluid Lens Proportion Constant used to calculate the radius of the
! fluid lens that surrounds the particle for particle-fluid-particle
! conduction.  Default [ 0.2 ]
      DOUBLE PRECISION FLPC

! Mininum separation distance between the surface of two contacting
! particles. This value is used to remove the singluarity that the
! particle-fluid-particle conduciton model develops at the contact
! interface. [4.0x10^(-10) meters]
      DOUBLE PRECISION DES_MIN_COND_DIST


! Run time flags for calculating the various modes of heat transfer
      LOGICAL :: CALC_CONV_DES
      LOGICAL :: CALC_COND_DES(DIM_M)
      LOGICAL :: CALC_RADT_DES(DIM_M)

! Fluid/Particle coupling
!---------------------------------------------------------------------//
! Source term for TFM gas phase energy equation.
      DOUBLE PRECISION, ALLOCATABLE :: DES_ENERGY_SOURCE(:)


      END MODULE DES_THERMO
