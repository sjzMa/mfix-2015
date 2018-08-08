!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_GAMA                                               !
!  Author: M. Syamlal                                 Date: 15-JUL-92  !
!                                                                      !
!  Purpose: Calculate gas-solids heat transfer coefficients            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_GAMA(HEAT_TR)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE physprop
      USE geometry
      USE fldvar
      USE energy
      USE rxns
      USE indices
      USE compar
      USE sendrecv
      USE discretelement
      USE functions

      IMPLICIT NONE

!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
!                      Flag for exchange functions
      LOGICAL          HEAT_TR(0:DIMENSION_M, 0:DIMENSION_M)

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
!                      Indices
      INTEGER          I, IJK

!                      Solids phase
      INTEGER          M

!                      Cell center value of U_g, Vg, and Wg [Recall
!                      U_g(IJK) refers to U_g at I+1/2, J, K
!                      V_g(IJK) refers to V_g at I, J+1/2, K and
!                      W_g [W_G(IJK) refers to W_G at I, J, K+1/2]
      DOUBLE PRECISION UGC, VGC, WGC

!                      Cell center value of U_sm, V_sm, W_sm
      DOUBLE PRECISION USCM, VSCM, WSCM

!                      Gas-solids relative velocity
      DOUBLE PRECISION VREL

!                      function of Prandtl number, Pr^(1/3)
      DOUBLE PRECISION Pr1o3

!                      Reynolds number, Re
      DOUBLE PRECISION Re

!                      EP_g^2
      DOUBLE PRECISION EP_g2

!                      a factor
      DOUBLE PRECISION FAC

      INTEGER          LM

!-----------------------------------------------

      IF (DISCRETE_ELEMENT) RETURN

! this needs to be generalized
      DO M = 1, SMAX
         IF (HEAT_TR(0,M)) THEN

            DO IJK = ijkstart3, ijkend3

               IF (FLUIDorP_FLOW_AT(IJK)) THEN
                  I = I_OF(IJK)
                  EP_G2 = EP_G(IJK)*EP_G(IJK)

! Calculate Prandtl number to the 1/3 power
                  if(K_G(IJK) > ZERO) then
                    PR1O3 = (C_PG(IJK)*MU_G(IJK)/K_G(IJK))**(1.D0/3.D0)
                  else
                    PR1O3 = LARGE_NUMBER
                  endif

! Calculate velocity components at the cell center for gas phase
                  UGC = HALF*(U_G(IJK)+U_G(IM_OF(IJK)))
                  VGC = HALF*(V_G(IJK)+V_G(JM_OF(IJK)))
                  WGC = HALF*(W_G(IJK)+W_G(KM_OF(IJK)))

! Calculate velocity components at the cell center for solids phase m
                  USCM = HALF*(U_S(IJK,M)+U_S(IM_OF(IJK),M))
                  VSCM = HALF*(V_S(IJK,M)+V_S(JM_OF(IJK),M))
                  WSCM = HALF*(W_S(IJK,M)+W_S(KM_OF(IJK),M))

! Calculate the magnitude of gas-solids relative velocity
                  VREL=SQRT((UGC-USCM)**2+(VGC-VSCM)**2+(WGC-WSCM)**2)

                  if(MU_G(IJK) > ZERO)then
                     RE = EP_G(IJK)*D_P(IJK,M)*VREL*RO_G(IJK)/MU_G(IJK)
                  else
                     RE = LARGE_NUMBER
                  endif

! Calculate gas-solids heat transfer coefficient (Gunn 1978)
                  GAMA_GS(IJK,M) = ((7.D0 - 10.D0*EP_G(IJK)+5.D0*EP_G2)*(ONE+0.7D0*RE**&
                     0.2D0*PR1O3)+(1.33D0 - 2.4D0*EP_G(IJK)+1.2D0*EP_G2)*RE**0.7D0*PR1O3)*(&
                     K_G(IJK)/D_P(IJK,M))*(6.D0*EP_S(IJK,M)/D_P(IJK,M))

! Correct the heat transfer coefficient for transpiration
! Bird, Stewart, and Lightfoot (1960, p.663)
                  IF (GAMA_GS(IJK,M) > SMALL_NUMBER) THEN
                     LM = 1
                     FAC = R_PHASE(IJK,LM)*C_PG(IJK)/GAMA_GS(IJK,M)
                     IF (ABS(FAC) < 0.1D0) THEN
                        GAMA_GS(IJK,M)=GAMA_GS(IJK,M)/(ONE+FAC/2.D0+FAC*FAC/6.D0)
                     ELSE
                        IF (R_PHASE(IJK,LM) > ZERO) THEN
                           GAMA_GS(IJK,M) = R_PHASE(IJK,LM)*C_PG(IJK)*EXP((-FAC&
                              ))/(ONE - EXP((-FAC)))
                        ELSE
                           GAMA_GS(IJK,M) = R_PHASE(IJK,LM)*C_PG(IJK)/(EXP(FAC)&
                               - ONE)
                        ENDIF
                     ENDIF
                  ENDIF   ! end if gama_gs(ijk,M) > 0

               ENDIF   ! end if (fluidorp_flow_at(ijk)
            ENDDO    ! end do loop (ijk=ijkstart3,ijkend3)
         ENDIF    ! end if (heat_tr(0,m))
      ENDDO    ! end do loop (m=1,smax)

      call send_recv(GAMA_GS,2)

      RETURN
      END SUBROUTINE CALC_GAMA
