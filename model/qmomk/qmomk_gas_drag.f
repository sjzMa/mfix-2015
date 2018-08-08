!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: QMOMK_GAS_DRAG(A_M, B_M, IER, UV, VV, WV)              C
!  Purpose: QMOMK - Accounting for the equal and opposite drag force   C
!           on gas due to particles by introducing the drag            C
!           as a source term. Face centered                            C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Alberto Passalacqua                      Date:            C
!            Re-used to develop QMOMK_GAS_DRAG from the                C
!            original DES_GAS_DRAG                                     C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE QMOMK_GAS_DRAG(A_M, B_M, IER, UV, VV, WV)
!-----------------------------------------------
!     M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE matrix
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE visc_g
      USE rxns
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE is
      USE tau_g
      USE bc
      USE compar
      USE sendrecv
      USE discretelement
      USE qmom_kinetic_equation
      USE drag
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!
!     Error index
      INTEGER          IER
!
!     Indices
      INTEGER          IJK
!
!     Phase index
      INTEGER          M, UV, VV, WV
!     Averaging Factor
      DOUBLE PRECISION :: AVG_FACTOR
!     Grid indices
      INTEGER          I,J,K,IN
!
      DOUBLE PRECISION USFCM, VSFCM, WSFCM
      DOUBLE PRECISION A_M(DIMENSION_3, -3:3, 0:DIMENSION_M)
      DOUBLE PRECISION B_M(DIMENSION_3, 0:DIMENSION_M)
      DOUBLE PRECISION tmp_A, tmp_B

!-----------------------------------------------
      AVG_FACTOR = 0.25D0*(DIMN-2) + 0.5D0*(3-DIMN)
      IF(UV.EQ.1) THEN
         DO M = 1, MMAX
            DO IJK = IJKSTART3, IJKEND3
               IF(FLUID_AT(IJK)) THEN
                  I = I_OF(IJK)
                  J = J_OF(IJK)
                  K = K_OF(IJK)

                  DO IN = 1, QMOMK_NN
                    USFCM = AVG_X(QMOMK_U1(IN,IJK,M),QMOMK_U1(IN,EAST_OF(IJK),M),I_OF(IJK))
                    tmp_A = -AVG_X(QMOMK_F_GS(IN,IJK,M),QMOMK_F_GS(IN,EAST_OF(IJK),M),I)*VOL_U(IJK)
                    tmp_B = tmp_A*USFCM

                    A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A
                    B_M(IJK,0) = B_M(IJK,0) + tmp_B
                  END DO
               END IF
            END DO
         END DO


      ELSE IF(VV.EQ.1) THEN

         DO M = 1, MMAX
            DO IJK = IJKSTART3, IJKEND3
               IF(FLUID_AT(IJK)) THEN
                  I = I_OF(IJK)
                  J = J_OF(IJK)
                  K = K_OF(IJK)

                  DO IN = 1, QMOMK_NN
                     VSFCM = AVG_Y(QMOMK_V1(IN,IJK,M), QMOMK_V1(IN,NORTH_OF(IJK),M),J_OF(IJK))
                     tmp_A = -AVG_Y(QMOMK_F_GS(IN,IJK,M),QMOMK_F_GS(IN,NORTH_OF(IJK),M),J)*VOL_V(IJK)
                     tmp_B = tmp_A*VSFCM

                     A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A
                     B_M(IJK,0) = B_M(IJK,0) + tmp_B
                  END DO
               END IF

            END DO
         END DO
90       FORMAT(3(1x,i2),10(2x,E12.5))
         close(900)

      ELSE IF(WV.EQ.1) THEN
         DO M = 1, MMAX
            DO IJK = IJKSTART3, IJKEND3
               IF(FLUID_AT(IJK)) THEN
                  I = I_OF(IJK)
                  J = J_OF(IJK)
                  K = K_OF(IJK)

                  DO IN = 1, QMOMK_NN
                     WSFCM = AVG_Z(QMOMK_W1(IN,IJK,M),QMOMK_W1(IN,TOP_OF(IJK),M),K_OF(IJK))
                     tmp_A = -AVG_Z(QMOMK_F_GS(IN,IJK,M),QMOMK_F_GS(IN,TOP_OF(IJK),M),K)*VOL_W(IJK)
                     tmp_B = tmp_A*WSFCM

                     A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A
                     B_M(IJK,0) = B_M(IJK,0) + tmp_B
                  END DO
               END IF
            END DO
         END DO

      END IF

      RETURN
      END SUBROUTINE QMOMK_GAS_DRAG
