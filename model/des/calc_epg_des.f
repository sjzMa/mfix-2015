!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: CALC_EPG_DES                                            !
!  Author: R.Garg                                     Date: ??-???-??  !
!                                                                      !
!  Purpose: Calculate the gas phase volume fraction (and in turn the   !
!  gas phase bulk density) from the sum of the solids volume fractions.!
!                                                                      !
!  NOTE: This routine uses a global communication to notify all ranks  !
!  of potential errors. Therefore all ranks can call MFIX_EXIT and     !
!  prevent dead-lock. Communications may be reduced by passing the     !
!  flag back to the caller and combining with other error checks.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_EPG_DES

! Global Variables:
!---------------------------------------------------------------------//
! Flag: Coupled gas-solids DEM simulation
      use discretelement, only: DES_CONTINUUM_HYBRID
! Flag: Discrete and continuum solids co-exist
      use discretelement, only: DES_CONTINUUM_COUPLED
! Number of discrete solids phases
      use discretelement, only: DES_MMAX
! Global ID of particles
      use discretelement, only: iGLOBAL_ID
! Particle positions
      use discretelement, only: DES_POS_NEW
! Number of continuum solids phases
      use physprop, only: SMAX
! Discrete particle material and bulk densities
      use discretelement, only: DES_RO_S, DES_ROP_s
! Number of particles in indexed fluid cell
      use discretelement, only: PINC
! List of particles in each cell.
      use discretelement, only: PIC
! Gas phae volume fraction, density, and build density
      use fldvar, only: EP_G, RO_G, ROP_G
! Bulk density of continuum solids phases
      use fldvar, only: EP_S
! Volume of scalar grid cell.
      use geometry, only: VOL
! Flag: Status of indexed cell
      use cutcell, only: CUT_CELL_AT
! Flag: Indexed cell contains fluid
      USE functions, only: FLUID_AT
! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! Flag: Fluid exists at indexed cell
      use functions, only: FLUID_AT
! The I, J, and K values that comprise an IJK
      use indices, only: I_OF, J_OF, K_OF
! Rank ID of current process
      use compar, only: myPE
! Global communication function to sum to all ranks.
      use mpi_utility, only: GLOBAL_ALL_SUM

! Global Parameters:
!---------------------------------------------------------------------//
      USE param1, only: ZERO, ONE

      use error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop indices
      INTEGER :: IJK, M, LC
! Total solids volume fraction
      DOUBLE PRECISION SUM_EPS
! Integer Error Flag
      INTEGER :: IER
!......................................................................!

! Initialize error flag.
      IER = 0

! Calculate gas volume fraction from solids volume fraction:
!---------------------------------------------------------------------//
!$omp parallel do if(ijkend3 .ge. 2000) default(none) reduction(+:IER) &
!$omp shared(IJKSTART3, IJKEND3, DES_CONTINUUM_COUPLED, DES_MMAX, SMAX,&
!$omp    EP_G, RO_G, ROP_G, DES_ROP_S, DES_RO_S, DES_CONTINUUM_HYBRID) &
!$omp private(IJK, SUM_EPs, M)
      DO IJK = IJKSTART3, IJKEND3
! Skip wall cells.
         IF(.NOT.FLUID_AT(IJK)) CYCLE
! Initialize EP_g and the accumulator.
         EP_G(IJK) = ONE
         SUM_EPS = ZERO
! Sum the DES solids volume fraction.
         DO M = 1, DES_MMAX
            SUM_EPS = SUM_EPS + DES_ROP_S(IJK,M)/DES_RO_S(M)
         ENDDO
! Add in TFM solids contributions.
         IF(DES_CONTINUUM_HYBRID) THEN
            DO M = 1,SMAX
               SUM_EPS = SUM_EPS + EP_S(IJK,M)
            ENDDO
         ENDIF
! Calculate the gas phase volume fraction and bulk density.
         EP_G(IJK) = ONE - SUM_EPS
         ROP_G(IJK) = RO_G(IJK) * EP_G(IJK)
! Flag an error if gas volume fraction is unphysical.
         IF(DES_CONTINUUM_COUPLED) THEN
            IF(EP_G(IJK) <= ZERO .OR. EP_G(IJK) > ONE) IER = IER + 1
         ENDIF
      ENDDO
!omp end parallel do


      CALL GLOBAL_ALL_SUM(IER)
      IF(IER == 0) RETURN


! Report any errors. Volume fraction errors are fatal.
!---------------------------------------------------------------------//
      CALL INIT_ERR_MSG("CALC_EPG_DES")
      CALL OPEN_PE_LOG(IER)

      WRITE(ERR_MSG, 1100)
      CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

 1100 FORMAT('Error 1100: Unphysical gas phase volume fraction ',      &
         'calculated. A .vtp',/'file will be written and the code ',   &
         'will exit. Fluid cell details:')

         DO IJK=IJKSTART3, IJKEND3
            IF(.NOT.FLUID_AT(IJK)) CYCLE
            IF(EP_G(IJK) > ZERO .AND. EP_G(IJK) <= ONE) CYCLE

            WRITE(ERR_MSG,1101) trim(iVal(IJK)), trim(iVal(I_OF(IJK))),&
               trim(iVal(J_OF(IJK))), trim(iVal(K_OF(IJK))),EP_G(IJK), &
               CUT_CELL_AT(IJK), trim(iVal(PINC(IJK))), VOL(IJK)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

            WRITE(ERR_MSG,1102)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            DO LC=1,PINC(IJK)
               M=PIC(IJK)%P(LC)
               WRITE(ERR_MSG,1103) iGlobal_ID(M), trim(iVal(           &
                  DES_POS_NEW(1,M))), trim(iVal(DES_POS_NEW(2,M))),    &
                  trim(iVal(DES_POS_NEW(3,M)))
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            ENDDO
         ENDDO

 1101 FORMAT(/3x,'Fluid Cell IJK: ',A,6x,'I/J/K: (',A,',',A,',',A,')',/&
         T6,'EP_G = ',g11.4,T30,'CUT_CELL_AT = ',L1,/T6,'PINC: ',A,T30,&
         'VOL = ',g11.4)

 1102 FORMAT(/T6,'Global ID',T30,'Position')

 1103 FORMAT(T6,I9,3x,'(',A,', ',A,', ',A,')')

      WRITE(ERR_MSG, 1104)
      CALL FLUSH_ERR_MSG(HEADER=.FALSE.)
 1104 FORMAT('This is a fatal error. A particle output file (vtp) ',   &
         'will be written',/'to aid debugging.')

      CALL WRITE_DES_DATA
      CALL MFIX_EXIT(myPE)

      END SUBROUTINE CALC_EPG_DES
