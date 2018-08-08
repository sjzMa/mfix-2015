!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CALC_DRAG                                               !
!  Author: M. Syamlal                                 Date: 29-JAN-92  !
!                                                                      !
!  Purpose: Calculate the gas solids and solids-solids drag terms if   !
!           directed to do so by the corresponding flags               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE CALC_DRAG(DRAGD, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE drag
      USE compar
      USE qmom_kinetic_equation

      use discretelement, only: DES_EXPLICITLY_COUPLED
      use discretelement, only: DES_CONTINUUM_COUPLED
      use discretelement, only: DES_CONTINUUM_HYBRID

      IMPLICIT NONE

!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Error index
      INTEGER :: IER
! Flag for exchange functions
      LOGICAL :: DRAGD(0:DIMENSION_M, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Solids phase indices
      INTEGER :: M, L
!-----------------------------------------------

! Alberto Passalacqua:  QMOMB
      IF (QMOMK) RETURN


! calculate drag between continuum phases (gas-solid & solids-solids)
      IF (.NOT.DES_CONTINUUM_COUPLED .OR. DES_CONTINUUM_HYBRID) THEN
         DO M = 1, SMAX
            IF (DRAGD(0,M) .AND. RO_G0/=ZERO) CALL DRAG_GS(M, IER)
         ENDDO

         IF (GRANULAR_ENERGY) THEN
            SELECT CASE (KT_TYPE_ENUM)
               CASE(LUN_1984, SIMONIN_1996, AHMADI_1995)
                  DO M = 1, SMAX
                     DO L = 1, M - 1
                        IF (DRAGD(L,M)) CALL DRAG_SS (L, M, IER)
                     ENDDO
                  ENDDO
               CASE(IA_2005)
                  DO M= 1, SMAX
                     DO L = 1, SMAX
                        IF (DRAGD(L,M)) CALL DRAG_SS_IA (L,M,IER)
                     ENDDO
                  ENDDO
               CASE (GD_1999, GTSH_2012, GHD_2007)
! strictly speaking gd and gtsh are monodisperse theories and so
! do not have solids-solids drag. ghd theory is a polydisperse
! theory but is self contained and does not invoke this routine

               CASE DEFAULT
! should never hit this
                  WRITE (*, '(A)') 'CALC_DRAG'
                  WRITE (*, '(A,A)') 'Unknown KT_TYPE: ', KT_TYPE
               call mfix_exit(myPE)
            END SELECT
         ELSE
            DO M = 1,SMAX
               DO L = 1, M - 1
                  IF (DRAGD(L,M)) CALL DRAG_SS (L, M, IER)
               ENDDO
            ENDDO
         ENDIF
      ENDIF

! calculate drag between continuum phases and discrete particles
! (gas-particle & solids-particle)
      IF (DES_CONTINUUM_COUPLED .AND. .NOT.DES_EXPLICITLY_COUPLED) &
         CALL CALC_DRAG_DES_2FLUID


      RETURN
      END SUBROUTINE CALC_DRAG
