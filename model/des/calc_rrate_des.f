!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine name: CALC_RRATE_DES                                     !
!                                                                      !
!  Author: J.Musser                                   Date: 16-May-11  !
!                                                                      !
!  Purpose: This routine manages gas-solid reactions for the continuum !
!  phase. Reactions are calculated from the DEM and the continuum      !
!  variables are integrated over the solids time step. These values    !
!  are then applied to the continuum phase by distributing them over   !
!  the full fluid time step. This ensures that mass is conserved       !
!  between the TFM/DEM and observes the variable time step within the  !
!  TFM.                                                                !
!                                                                      !
!  Although mass is conserved, the two models (TFM/DEM) are out of     !
!  sync. The general approach is as follows.                           !
!                                                                      !
!  Time: T0                                                            !
!  TFM-A) Hydrodynamcics (gas phase) [stationary solids]               !
!  TFM-B) Homogeneous (gas phase) reactions                            !
!  TFM-C)  < Do Nothing >                                              !
!  TFM-C) Apply gas/DEM reactions from T0 to TFM  <<  [out of sync]    !
!  TFM-D) Updated TFM variables: Time (T0 --> T1)                      !
!  DEM-E) Particle/Particle/Wall collisions                            !
!  DEM-F) Hydrodynamics [constant gas field]                           !
!  DEM-G) Calculate heterogeneous (gas/DEM) reactions                  !
!  DEM-H) Update DEM variables: Time (T0 --> T0.1)                     !
!  * Repeat DEM-[E-H] until T0.x=T1 (T0 -> T0.1 -> T0.2 -> ... -> T1)  !
!                                                                      !
!  TFM/DEM are at the same physical time (T1).                         !
!                                                                      !
!  Time: T1                                                            !
!  TFM-A) Hydrodynamcics (gas phase) [stationary solids]               !
!  TFM-B) Homogeneous (gas phase) reactions                            !
!  TFM-C) Apply gas/DEM reactions from T0 to TFM  <<  [out of sync]    !
!    ...                                                               !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_RRATE_DES(CLEAR)

      USE compar
      USE des_rxns
      USE discretelement
      USE energy
      USE fldvar
      USE geometry
      USE interpolation
      USE param1
      USE physprop
      USE run
      USE rxns
      USE usr

      IMPLICIT NONE

! Passed Variables
!---------------------------------------------------------------------//
      LOGICAL, INTENT(IN) :: CLEAR

! Local variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION DEM_to_TFM
      INTEGER IJK


! Initialize global storage arrays to zero
!---------------------------------------------------------------------//
      IF(CLEAR) THEN
         SUM_R_G(:) = ZERO
         HOR_G(:) = ZERO
         R_GP(:,:) = ZERO
         ROX_GC(:,:) = ZERO
         R_PHASE(:,:) = ZERO
      ENDIF


      DO IJK = IJKSTART3, IJKEND3
         DEM_to_TFM = DT * VOL(IJK)
         R_gp(IJK,:) = R_gp(IJK,:) + &
            DES_R_gp(IJK,:)/DEM_to_TFM
         WHERE(X_g(IJK,:) > SMALL_NUMBER)
            RoX_gc(IJK,:) = RoX_gc(IJK,:) + &
               DES_R_gc(IJK,:) / (DEM_to_TFM * X_g(IJK,:))
         ELSEWHERE
            RoX_gc(IJK,:) = RoX_gc(IJK,:) + 1.0d-9/DEM_to_TFM
         ENDWHERE
         R_PHASE(IJK,:) = R_PHASE(IJK,:) + &
            DES_R_PHASE(IJK,:) / DEM_to_TFM
         SUM_R_g(IJK) = SUM_R_g(IJK) + &
            DES_SUM_R_g(IJK) / DEM_to_TFM
         HOR_g(IJK) = HOR_g(IJK) + &
            DES_HOR_g(IJK) / DEM_to_TFM
      ENDDO

      RETURN


      END SUBROUTINE CALC_RRATE_DES



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine name: CALC_RRATE_DES                                     !
!                                                                      !
!  Purpose: This routine manages gas-solid reactions for the continuum !
!  phase.                                                              !
!                                                                      !
!  Author: J.Musser                                   Date: 16-May-11  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ZERO_RRATE_DES

      USE des_rxns

      IMPLICIT NONE

      DES_R_gp(:,:) = ZERO
      DES_R_gc(:,:) = ZERO
      DES_R_PHASE(:,:) = ZERO
      DES_HOR_G(:) = ZERO
      DES_SUM_R_g(:) = ZERO

      RETURN
      END SUBROUTINE ZERO_RRATE_DES
