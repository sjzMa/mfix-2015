!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine name: PHYSICAL_PROP                                      !
!                                                                      !
!  Purpose: Calculate physical properties that vary with time.         !
!                                                                      !
!  Author: J.Musser                                   Date: 09-May-11  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_PHYSICAL_PROP(NP, FOCUS)

      Use des_rxns
      Use des_thermo
      Use discretelement
      Use funits
      Use param
      Use param1
      Use physprop
      Use run

! Calculate the specific heat from polynomical data obtained from the
! thermodynamic databases.
      use read_thermochemical, only: DES_calc_CpoR

! Universal gas constant in cal/mol.K
      use constant, only: RGAS => GAS_CONST_cal

      IMPLICIT NONE

! Passed Variables
!-----------------------------------------------------------------------
! Index of particle
      INTEGER, INTENT(IN) :: NP
! Logical indicating to write additional data about the particle for
! debugging purposes.
      LOGICAL, INTENT(IN) :: FOCUS

! Local Variables
!-----------------------------------------------------------------------
! Dummy indices
      INTEGER M, N
! error indicator
      INTEGER IER

      DOUBLE PRECISION :: lCPoR

! Get the solids phase of the particle
      M = PIJK(NP,5) + SMAX

! Specific heat
!-----------------------------------------------------------------------
! This only needs calculated when solving the energy equations.
      IF(ENERGY_EQ) THEN
! If a constant value specific heat has not been defined in the mfix.dat
! file, calculate the temperature dependent value based on data from the
! thermodynamic databases.
         IF (C_PS0(M) == UNDEFINED) THEN
! Read the thermodynamic database if it has not already been done.
            DES_C_PS(NP) = ZERO
! Calculate the specific heat based on the species composition of the
! particle and the data from the thermodynamic databases.
            DO N = 1, NMAX_s(M)
               lCPoR = DES_calc_CpoR(DES_T_s_NEW(NP), M, N, IER)
               DES_C_PS(NP) = DES_C_PS(NP) +                          &
                  lCpoR * DES_X_s(NP,N) * RGAS / MW_s(M,N)
            ENDDO
! Convert to SI units if needed.
            IF (UNITS == 'SI') DES_C_PS(NP) = 4183.925D0*DES_C_PS(NP)
         ELSE
! If a constant value specific heat has been assigned to the particle
! in the mfix.dat file, use this value.
            DES_C_PS(NP) = C_PS0(M)
         ENDIF
      ENDIF

!      WRITE(*,"(3X,A,I2,A,F10.6)")'DES_C_PS(',NP,'): ',DES_C_PS(NP)

      RETURN

      END SUBROUTINE DES_PHYSICAL_PROP
