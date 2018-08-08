!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                           C
!  Module name: CALC_K_s(M, IER)                                            C
!  Purpose: Calculate the effective conductivity of solids phases           C
!                                                                           C
!  Author:M. Syamlal                                  Date: 24-APR-96       C
!  Reviewer:                                          Date: dd-mmm-yy       C
!                                                                           C
!  Revision Number: 01                                                      C
!  Purpose: (1) allow to use Bauer & Schlunder's (1978) model in CGS or SI  C
!           (2) If fluid_at(IJK) condition for the Bauer & Schlunder's modelC
!  Author:  S. Dartevelle                             Date: 10-July-02      C
!  Reviewer:                                          Date: dd-mmm-yy       C
!                                                                           C
!  Literature/Document References:                                          C
!                                                                           C
!  Variables referenced:                                                    C
!  Variables modified:                                                      C
!                                                                           C
!  Local variables:                                                         C
!                                                                           C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_K_S(M)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE physprop
      USE fldvar
      USE geometry
      USE indices
      USE constant
      USE toleranc
      USE compar
      USE sendrecv
      USE run
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------

! define two
      DOUBLE PRECISION, PARAMETER          :: TWO = 2.0d0

! microscopic conductivity of ash in cal/s.cm.K
! (not modified by the gas phase)
      DOUBLE PRECISION Ks_micro
      PARAMETER (Ks_micro = 0.5258D-2)    !(2.2 J/s.m.K)

! constant in conductivity equation
      DOUBLE PRECISION PHI_k
      PARAMETER (PHI_k = 7.26D-3)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!                      Indices
      INTEGER          IJK

!                      Solids phase
      INTEGER          M

!                      Quantities in solids conductivity formula
      DOUBLE PRECISION B, R_km, BoR, L_rm

!                      Transform K_g(IJK) into the CGS if we work with SI
      DOUBLE PRECISION Kg_micro
!-----------------------------------------------

      IF (K_S0(M) /= UNDEFINED) RETURN

!!!!$omp parallel do private(IJK,B,R_km,BoR,L_rm,Kg_micro) &
!!!!$omp& schedule(dynamic,chunk_size)
      DO IJK = ijkstart3, ijkend3

! All calculations are in CGS (1 cal = 4.183925J)
         IF (FLUID_AT(IJK)) THEN
            IF (UNITS == 'SI') THEN
               Kg_micro = K_g(IJK)/418.3925D0  !convert K_g to CGS units (cal/s.cm.K)
            ELSE
               Kg_micro = K_g(IJK)     ! K_g already in CGS units (cal/s.cm.K)
            ENDIF

! Bauer & Schlunder's (1978) theory:
            IF( EP_s(IJK,M) >  DIL_EP_s) THEN
               B = 1.25D0 * ((ONE - EP_g(IJK))/EP_g(IJK))**(10.D0/9.D0)
               R_km = Ks_micro/Kg_micro
               BoR  = B/R_km
               L_rm = -(TWO/(ONE-BoR)) * &
                  ( ((R_km-ONE)/(ONE-BoR)**2)*BoR*LOG(BoR) + &
                    (B-ONE)/(ONE-BoR) + (B+ONE)/TWO )
! K_s is the macroscopic conductivity that has been modified by the presence of
! the gas phase (cal/s.cm.K)
               K_S(IJK,M) = (Phi_k*R_km + (ONE-Phi_k)*L_rm)*&
                  Kg_micro/SQRT(ONE - EP_g(IJK))
            ELSE
               K_S(IJK, M) = ZERO
            ENDIF

! An approximate average value for the solids conductivity is 2.5*K_g
!         K_S(IJK,M) = 2.5*Kg_micro            !in CGS system

         ELSE   ! else branch if(fluid_at(ijk))
            K_S(IJK,M) = ZERO
         ENDIF   ! end if/else (fluid_at(ijk))

         IF (UNITS == 'SI') K_s(IJK, M) = 418.3925D0*K_s(IJK, M)   !J/s.m.K

      ENDDO   ! end do (ijk=ijkstart3,ijkend3)

      CALL send_recv(K_S, 2)

      RETURN
      END SUBROUTINE CALC_K_S


