!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP                                           C
!  Purpose: Calculate the indicated physical properties that vary      C
!           with time if directed to do so by the corresponding flag   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUL-92  C
!  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Mods for MFIX 2.0 (old name CALC_PHYSPROP)                 C
!  Author: M. Syamlal                                 Date: 23-APR-96  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: allow SI                                                   C
!  Author: S. Dartevelle                              Date: 01-Jul-02  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!    Perry, R.H., and C.H. Chilton, Chemical Engineer's Handbook, 5th  C
!      edition, McGraw-Hill Kogakusha, Tokyo, 1973.                    C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PHYSICAL_PROP(IER, LEVEL)

      use compar
      use funits
      use geometry
      use indices
      use mpi_utility
      use param, only: dimension_m
      use param1
      use physprop

      use coeff, only: DENSITY  ! Density
      use coeff, only: SP_HEAT  ! Specific heat
      use coeff, only: PSIZE    ! Particle diameter

      implicit none

! Dummy arguments
!-----------------------------------------------------------------------
! Global error Flag.
      INTEGER, intent(inout) :: IER
      INTEGER, intent(in) :: LEVEL

! Local variables
!-----------------------------------------------------------------------
! Arrays for storing errors:
! 100 - Negative gas phase density
! 101 - Negative solids phase density
! 10x - Unclassified
! 900 - Invalid temperature in calc_CpoR
      INTEGER :: Err_l(0:numPEs-1)  ! local
      INTEGER :: Err_g(0:numPEs-1)  ! global

! Initialize error flags.
      Err_l = 0

! Calculate density only. This is invoked several times within iterate,
! making it the most frequently called.
      if(LEVEL == 0) then
         if(DENSITY(0)) CALL PHYSICAL_PROP_ROg
         if(any(DENSITY(1:DIMENSION_M))) CALL PHYSICAL_PROP_ROs

! Calculate everything except density. This is called at the start of
! each iteration.
      elseif(LEVEL == 1) then
         if(SP_HEAT(0)) CALL PHYSICAL_PROP_CPg
!         if(any(DENSITY(1:DIMENSION_M))) CALL PHYSICAL_PROP_ROs
         if(any(SP_HEAT(1:DIMENSION_M))) CALL PHYSICAL_PROP_CPs
         if(any(PSIZE(1:DIMENSION_M))) CALL PHYSICAL_PROP_Dp


! Calculate everything. This is invoked via calc_coeff_all as part of
! the initialization (before starting the time march) and at the start
! of each step step thereafter.
      elseif(LEVEL == 2) then
         if(DENSITY(0)) CALL PHYSICAL_PROP_ROg
         if(SP_HEAT(0)) CALL PHYSICAL_PROP_CPg
         if(any(DENSITY(1:DIMENSION_M))) CALL PHYSICAL_PROP_ROs
         if(any(SP_HEAT(1:DIMENSION_M))) CALL PHYSICAL_PROP_CPs
         if(any(PSIZE(1:DIMENSION_M))) CALL PHYSICAL_PROP_Dp
      endif


! In case of negative density force exit from the physical property
! calculation routine and reduce the time step
      CALL global_all_sum(Err_l, Err_g)
      IER = maxval(Err_g)
      if(IER == 0) return


! Error handeling. - Local.
!-----------------------------------------------------------------------
! An invalid temperature was found by calc_CpoR. This is a fatal run-
! time error and forces a call to MFIX_EXIT.
      IF(IER == 901 .OR. IER == 902) then
         if(myPE == PE_IO) then
            write(*,2000) IER
            write(UNIT_LOG,2000) IER
         endif
         CALL MFIX_EXIT(myPE)
      ENDIF

      return

 2000 FORMAT(/1X,70('*')/' From: PHYSICAL_PROP',/' Fatal Error 2000:', &
         ' calc_CpoR reporetd an invalid temperature: 0x0', I3/,       &
         'See Cp.log for details. Calling MFIX_EXIT.',/1X,70('*')/)

      contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP_ROg                                       C
!  Purpose: Calculate the gas phase density.                           C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PHYSICAL_PROP_ROg


! Global variables:
!-----------------------------------------------------------------------
! Gas phase species mass fractions.
      use fldvar, only: X_g
! Gas phase temperature.
      use fldvar, only: T_g
! Gas phase density (compressible).
      use fldvar, only: RO_g
! Gas phase pressure.
      use fldvar, only: P_g
! Gas phase volume fraction.
      use fldvar, only: EP_g
! Gas phase material density.
      use fldvar, only: ROP_g
! Maximum value for molecular weight (divided by one)
      use toleranc, only: OMW_MAX
! Run time flag for generating negative gas density log files
      use run, only: REPORT_NEG_DENSITY
! Equation of State - GAS
      use eos, only: EOSG
! Flag for user defined function
      use run, only: USR_ROg

      use functions

      implicit none

! Local Variables:
!-----------------------------------------------------------------------
! Average molecular weight
      DOUBLE PRECISION :: MW
! Loop indicies
      INTEGER :: IJK   ! Computational cell
! Flag to write log header
      LOGICAL :: wHeader
!......................................................................!

! User-defined function
      IF(USR_ROg) THEN
         CALL USR_PHYSICAL_PROP_ROg
         RETURN
      ENDIF

! Initialize:
      wHeader = .TRUE.

! Average molecular weight: Xg1/Mw1 + Xg2/Mw2 + Xg3/Mw3 + ....
      IF(.NOT.database_read) call read_database0(IER)

      IJK_LP: DO IJK = IJKSTART3, IJKEND3
         IF(WALL_AT(IJK)) cycle IJK_LP
         IF (MW_AVG == UNDEFINED) THEN
! Calculating the average molecular weight of the fluid.
            MW = SUM(X_G(IJK,:NMAX(0))/MW_G(:NMAX(0)))
            MW = ONE/MAX(MW,OMW_MAX)
            MW_MIX_G(IJK) = MW
! Calculate the fluid density and bulk density
            RO_G(IJK) = EOSG(MW,P_G(IJK),T_G(IJK))
            ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
         ELSE
            RO_G(IJK) = EOSG(MW_AVG,P_G(IJK),T_G(IJK))
            ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
         ENDIF

         IF(RO_G(IJK) < ZERO) THEN
            Err_l(myPE) = 100
            IF(REPORT_NEG_DENSITY)CALL ROgErr_LOG(IJK, wHeader)
         ENDIF
      ENDDO IJK_LP


      RETURN
      END SUBROUTINE PHYSICAL_PROP_ROg


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP_ROs                                       C
!  Purpose: Calculate solids phase (variable) density.                 C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PHYSICAL_PROP_ROs

! Global variables:
!-----------------------------------------------------------------------
! Solid phase species mass fractions.
      use fldvar, only: X_s
! Solid phase density (variable).
      use fldvar, only: ROP_s, RO_s
! Baseline/Unreaced solids density
      use physprop, only: BASE_ROs
! Initial mass fraction of inert species
      use physprop, only: X_S0
! Index of inert solids phase species.
      use physprop, only: INERT_SPECIES
! Inert solids phase species mass fraction in dilute region.
      use physprop, only: DIL_INERT_X_VSD
! Factor to define dilute region where DIL_INERT_X_VSD is used
      use physprop, only: DIL_FACTOR_VSD
! Flag for variable solids density.
      use run, only: SOLVE_ROs
! Run time flag for generating negative density log files.
      use run, only: REPORT_NEG_DENSITY
! Minimum solids volume fraction
      use toleranc, only: DIL_EP_s
! Equation of State - Solid
      use eos, only: EOSS
! Flag for user defined function
      use run, only: USR_ROs

      use functions

      implicit none

! Local Variables:
!-----------------------------------------------------------------------
! Loop indicies
      INTEGER :: IJK, M
! Index of inert species
      INTEGER :: IIS
! Flag to write log header
      LOGICAL :: wHeader
! Minimum bulk density
      DOUBLE PRECISION :: minROPs
!......................................................................!

! User defined function
      IF(USR_ROs) THEN
         CALL USR_PHYSICAL_PROP_ROs
         RETURN
      ENDIF

      M_LP: DO M=1, MMAX
         IF(.NOT.SOLVE_ROs(M)) cycle M_LP
! Initialize header flag.
         wHeader = .TRUE.
! Set the index of the inert species
         IIS = INERT_SPECIES(M)
! Calculate the minimum solids denisty.
         minROPs = BASE_ROs(M)*(DIL_FACTOR_VSD*DIL_EP_s)

! Calculate the solids denisty over all cells.
         IJK_LP: DO IJK = IJKSTART3, IJKEND3
            IF(WALL_AT(IJK)) cycle IJK_LP
            IF(ROP_s(IJK,M) > minROPs) THEN
               RO_S(IJK,M) = EOSS(BASE_ROs(M), X_s0(M,IIS),            &
                  X_s(IJK,M,IIS))
            ELSE
!               RO_s(IJK,M) = BASE_ROs(M)
               RO_S(IJK,M) = EOSS(BASE_ROs(M), X_s0(M,IIS),            &
                  DIL_INERT_X_VSD(M))
            ENDIF

! Report errors.
            IF(RO_S(IJK,M) <= ZERO) THEN
               Err_l(myPE) = 101
               IF(REPORT_NEG_DENSITY) CALL ROsErr_LOG(IJK, M, wHeader)
            ENDIF
         ENDDO IJK_LP
      ENDDO M_LP


      RETURN
      END SUBROUTINE PHYSICAL_PROP_ROs


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP_CPg                                       C
!  Purpose: Calculate the gas phase constant pressure specific heat.   C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
! Notes:                                                               C
!  > Unit conversion: 1 cal = 4.183925 J                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PHYSICAL_PROP_CPg

! Global Variables:
!-----------------------------------------------------------------------
! Universal gas constant in cal/mol.K
      use constant, only: RGAS => GAS_CONST_cal
! Maximum value for molecular weight (divided by one)
      use toleranc, only: OMW_MAX
! Gas phase species mass fractions.
      use fldvar, only: X_g
! Gas phase temperature.
      use fldvar, only: T_g
! Units: CGS/SI
      use run, only: UNITS
! Function to calculate Cp over gas constant R
      use read_thermochemical, only: calc_CpoR
! Flag for user defined function
      use run, only: USR_CPg

      use functions

      implicit none

! Local Variables:
!-----------------------------------------------------------------------
! Species specific heat.
      DOUBLE PRECISION :: lCp
! Error flag returned from calc_CpoR
      INTEGER :: lCP_Err
      INTEGER :: gCP_Err
! Loop indicies
      INTEGER :: IJK, N
!......................................................................!

! User defined function
      IF(USR_CPg) THEN
         CALL USR_PHYSICAL_PROP_CPg
         RETURN
      ENDIF

!-----------------------------------------------------------------------

! Ensure that the database was read. This *should* have been caught by
! check_data_05 but this call remains to prevent an accident.
      IF(.NOT.database_read) CALL read_database0(IER)

      gCP_Err = 0
      lCP_Err = 0

      IJK_LP: DO IJK = IJKSTART3, IJKEND3
         IF(WALL_AT(IJK)) CYCLE IJK_LP
! Calculating an average specific heat for the fluid.
         C_PG(IJK) = ZERO
         DO N = 1, NMAX(0)
            lCp = calc_CpoR(T_G(IJK), 0, N, lCP_Err)
            C_PG(IJK) = C_PG(IJK) + X_g(IJK,N) * lCp * RGAS / MW_g(N)
            gCP_Err = max(gCP_Err, lCP_Err)
         ENDDO
      ENDDO IJK_LP

! The database calculation always returns cal/g.K thus the following
! conversion is needed if using SI units.
! Only convert useful part of the array. When the array is re-indexed,
! elements past IJKEND3 are UNDEFINED and would keep growing if 
! multiplied, leading to overflow.
      IF(UNITS == 'SI') THEN
         DO IJK = IJKSTART3, IJKEND3
           C_PG(IJK) = 4.183925d3 * C_PG(IJK)
         ENDDO
      ENDIF

! Increment the error to 900+ to invoke fatal exit.
      IF(gCP_Err /= 0) Err_l(myPE) = 800 + gCP_Err

      RETURN
      END SUBROUTINE PHYSICAL_PROP_CPg




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP_CPs                                       C
!  Purpose: Calculate solids phase constant pressure specific heat.    C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PHYSICAL_PROP_CPs

! Universal gas constant in cal/mol.K
      use constant, only: RGAS => GAS_CONST_cal
! Units: CGS/SI
      use run, only: UNITS
! Max value for molecular weight inverse
      use toleranc, only: OMW_MAX
! Solids temperature
      use fldvar, only: T_s
! Solids species mass fractions
      use fldvar, only: X_s
! Function to calculate Cp over gas constant R
      use read_thermochemical, only: calc_CpoR
! Flag for user defined function
      use run, only: USR_CPs

      use functions

      implicit none

! Local variables:
!---------------------------------------------------------------------//
! Local value for Cp
      DOUBLE PRECISION :: lCp
! Loop indicies
      INTEGER :: IJK, M, N
! Local error flag indicating that the Cp is out of range.
      INTEGER :: lCP_Err
!......................................................................!

! User defined function
      IF(USR_CPs) THEN
         CALL USR_PHYSICAL_PROP_CPs
         RETURN
      ENDIF

! Ensure that the database was read. This *should* have been caught by
! check_data_05 but this call remains to prevent an accident.
      IF(.NOT.database_read) CALL read_database0(IER)

      lCP_Err = 0

      M_LP: DO M=1, MMAX
         IJK_LP: DO IJK = IJKSTART3, IJKEND3
            IF(WALL_AT(IJK)) CYCLE IJK_LP
! Calculating an average specific heat for the fluid.
            C_PS(IJK, M) = ZERO

            DO N = 1, NMAX(M)
               lCp = calc_CpoR(T_S(IJK,M), M, N, lCP_Err)
               C_PS(IJK,M) = C_PS(IJK,M) + X_s(IJK,M,N) * &
                  (lCp * RGAS / MW_s(M,N))
            ENDDO

         ENDDO IJK_LP
      ENDDO M_LP

! The database calculation always returns cal/g.K thus the following
! conversion is needed if using SI units.
! Only convert useful part of the array. When the array is re-indexed,
! elements past IJKEND3 are UNDEFINED and would keep growing if 
! multiplied, leading to overflow.

      IF(UNITS == 'SI') THEN
         DO M=1, MMAX
            DO IJK = IJKSTART3, IJKEND3
              C_PS(IJK,M) = 4.183925d3 * C_PS(IJK,M)
            ENDDO
         ENDDO
      ENDIF
      END SUBROUTINE PHYSICAL_PROP_CPs

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP_CPs                                       C
!  Purpose: Calculate solids phase constant pressure specific heat.    C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PHYSICAL_PROP_Dp

      use run, only: CALL_DQMOM

      use fldvar, only: scalar
      use scalars, only: phase4scalar

      use fldvar, only: D_p, EP_S
      use functions

      implicit none

! Loop indicies
      INTEGER :: IJK   ! Computational cell
      INTEGER :: M     ! Solids phase

! Map from true index to map.
      INTEGER :: lM

      IF(.NOT.CALL_DQMOM) return

      M_LP: DO M=1, MMAX

         lM = phase4scalar(M) ! Map from scalar eq to solids phase

         IF(.NOT.PSIZE(M)) CYCLE M_LP
         IJK_LP: DO IJK = IJKSTART3, IJKEND3
            IF(WALL_AT(IJK)) CYCLE IJK_LP

            IF(EP_s(IJK,lM) > small_number) D_p(IJK,M)= Scalar(IJK,lM)

         ENDDO IJK_LP
      ENDDO M_LP

      return
      END SUBROUTINE PHYSICAL_PROP_Dp



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: NegROg_LOG                                              C
!  Purpose: Record information about the location and conditions that  C
!           resulted in a negative gas phase density.                  C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE ROgErr_LOG(IJK, tHeader)

! Simulation time
      use run, only: TIME
! Gas phase temperature.
      use fldvar, only: T_g
! Gas phase density (compressible).
      use fldvar, only: RO_g
! Gas phase pressure.
      use fldvar, only: P_g
      use cutcell

      INTEGER, intent(in) :: IJK
      LOGICAL, intent(inout) :: tHeader

      LOGICAL :: lExists
      CHARACTER(LEN=255) :: lFile
      INTEGER, parameter :: lUnit = 4868
      LOGICAL, save :: fHeader = .TRUE.


      lFile = '';
      if(numPEs > 1) then
         write(lFile,"('ROgErr_',I4.4,'.log')") myPE
      else
         write(lFile,"('ROgErr.log')")
      endif
      inquire(file=trim(lFile),exist=lExists)
      if(lExists) then
         open(lUnit,file=trim(adjustl(lFile)),                         &
            status='old', position='append', convert='big_endian')
      else
         open(lUnit,file=trim(adjustl(lFile)), status='new', convert='big_endian')
      endif

      if(fHeader) then
         write(lUnit,1000)
         fHeader = .FALSE.
      endif

      if(tHeader) then
         write(lUnit,"(/2x,'Simulation time: ',g12.5)") TIME
         tHeader = .FALSE.
      endif

      write(lUnit,1001) IJK, I_OF(IJK), J_OF(IJK), K_OF(IJK)
      write(lUnit,"(6x,A,1X,g12.5)",ADVANCE='NO') 'RO_g:', RO_g(IJK)
      write(lUnit,"(2x,A,1X,g12.5)",ADVANCE='NO') 'P_g:', P_g(IJK)
      write(lUnit,"(2x,A,1X,g12.5)") 'T_g:', T_g(IJK)
      if(CARTESIAN_GRID) then
         write(lUnit,"(6x,A,1X,L1)",ADVANCE='NO') 'Cut Cell:', CUT_CELL_AT(IJK)
         write(lUnit,"(6x,A,1X,L1)") 'Small Cell:', SMALL_CELL_AT(IJK)
         write(lUnit,"(6x,'Coordinates (E/N/T): ',1X,3(2x, g17.8))") &
            xg_e(I_OF(IJK)), yg_n(J_of(ijk)), zg_t(k_of(ijk))
      endif

      close(lUnit)

      RETURN

 1000 FORMAT(2X,'One or more cells have reported a negative gas',      &
         ' density (RO_g(IJK)). If',/2x,'this is a persistent issue,', &
         ' lower UR_FAC(1) in mfix.dat.')

 1001 FORMAT(/4X,'IJK: ',I8,7X,'I: ',I4,'  J: ',I4,'  K: ',I4)

      END SUBROUTINE ROgErr_LOG


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: NegROs_LOG                                              !
!  Purpose: Record information about the location and conditions that  !
!           resulted in a negative solids phase density.               !
!                                                                      !
!  Author: J. Musser                                  Date: 09-Oct-13  !
!  Reviewer:                                          Date:            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ROsErr_LOG(IJK, M, tHeader)

! Global variables:
!-----------------------------------------------------------------------
! Simulation time
      use run, only: TIME
! Solid phase species mass fractions.
      use fldvar, only: X_s
! Solid phase density (variable).
      use fldvar, only: RO_s
! Baseline/Unreaced solids density
      use physprop, only: BASE_ROs
! Initial mass fraction of inert species
      use physprop, only: X_S0
! Index of inert solids phase species.
      use physprop, only: INERT_SPECIES

! Full Access to cutcell data.
      use cutcell


! Local Variables:
!-----------------------------------------------------------------------
! Fluid cell and solids phase indices
      INTEGER, intent(in) :: IJK, M
! Flag to output header
      LOGICAL, intent(inout) :: tHeader
! Local aliase for inert species index
      INTEGER :: N
! Local file values.
      LOGICAL :: lExists
      CHARACTER(LEN=255) :: lFile
      INTEGER, parameter :: lUnit = 4868
      LOGICAL, save :: fHeader = .TRUE.


      lFile = '';
      if(numPEs > 1) then
         write(lFile,"('ROsErr_',I4.4,'.log')") myPE
      else
         write(lFile,"('ROsErr.log')")
      endif
      inquire(file=trim(lFile),exist=lExists)
      if(lExists) then
         open(lUnit,file=trim(adjustl(lFile)),                         &
            status='old', position='append',convert='big_endian')
      else
         open(lUnit,file=trim(adjustl(lFile)), status='new',convert='big_endian')
      endif

      if(fHeader) then
         write(lUnit,1000)
         fHeader = .FALSE.
      endif

      if(tHeader) then
         write(lUnit,"(/2x,'Simulation time: ',g12.5,'  Phase: ',I2)")&
             TIME, M
         tHeader = .FALSE.
      endif

      N = INERT_SPECIES(M)
      write(lUnit,1001) IJK, I_OF(IJK), J_OF(IJK), K_OF(IJK)
      write(lUnit,"(6x,A,1X,g12.5)",advance='no') 'RO_s:', RO_s(IJK,M)
      write(lUnit,"(2x,A,1X,g12.5)",advance='no') 'Base:', BASE_ROs(M)
      write(lUnit,"(2x,A,1X,g12.5)",advance='no') 'X_s0:', X_s0(M,N)
      write(lUnit,"(2x,A,1X,g12.5)") 'X_s:', X_s(IJK,M,N)

      if(CARTESIAN_GRID) then
         write(lUnit,"(6x,A,1X,L1)",ADVANCE='NO') 'Cut Cell:', CUT_CELL_AT(IJK)
         write(lUnit,"(6x,A,1X,L1)") 'Small Cell:', SMALL_CELL_AT(IJK)
         write(lUnit,"(6x,'Coordinates (E/N/T): ',1X,3(2x, g17.8))") &
            xg_e(I_OF(IJK)), yg_n(J_of(ijk)), zg_t(k_of(ijk))
      endif

      close(lUnit)

      RETURN

 1000 FORMAT(2X,'One or more cells have reported a negative gas',      &
         ' density (RO_g(IJK)). If',/2x,'this is a persistent issue,', &
         ' lower UR_FAC(1) in mfix.dat.')

 1001 FORMAT( 4X,'IJK: ',I8,7X,'I: ',I4,'  J: ',I4,'  K: ',I4)

      END SUBROUTINE ROsErr_LOG

      END SUBROUTINE PHYSICAL_PROP
