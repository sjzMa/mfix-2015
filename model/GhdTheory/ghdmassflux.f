!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GHDMASSFLUX                                             C
!  Purpose: Calculate the species mass flux 3-components of Joi at cellC
!           faces to compute species velocities and source terms in T  C
!           equation.                                                  C
!                                                                      C
!  Author: S. Benyahia                              Date: 13-MAR-09    C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!     C. Hrenya handnotes and Garzo, Hrenya, Dufty papers (PRE, 2007)  C
!                                                                      C
!  Variables modified:   JoiX  JoiY   JoiZ                             C
!                                                                      C
!  Local variables: all terms in mass flux                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE GHDMASSFLUX ()

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE fldvar
      USE indices
      USE visc_s
      USE ghdtheory
      USE physprop
      USE run
      USE constant
      USE toleranc
      USE drag
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Index
      INTEGER :: IJK, I, J, K
! Index
      INTEGER :: IJKE, IJKN, IJKT
! Solids phase
      INTEGER :: M, L
! number densities to compute del(Nj)
      DOUBLE PRECISION NjC, NjE, NjN, NjT

! mass of species
      DOUBLE PRECISION :: Mi, Mj
! number density of species
      DOUBLE PRECISION :: Ni
! mixture density and temperature at cell faces
      DOUBLE PRECISION :: ropsE, ropsN, ropsT, ThetaE, ThetaN, ThetaT
      DOUBLE PRECISION :: EPSA
! transport coefficient at cell faces
      DOUBLE PRECISION :: DiTE, DiTN, DiTT
      DOUBLE PRECISION :: DijE, DijN, DijT
      DOUBLE PRECISION :: DijFE, DijFN, DijFT
! Terms in the calculation of Joi-X,Y,Z
      DOUBLE PRECISION :: ordinDiffTermX, ordinDiffTermY, ordinDiffTermZ
      DOUBLE PRECISION :: massMobilityTermX, massMobilityTermY, &
                          massMobilityTermZ
      DOUBLE PRECISION :: massMobilityTermXvelupdate, &
                          massMobilityTermYvelupdate, &
                          massMobilityTermZvelupdate
      DOUBLE PRECISION :: thermalDiffTermX, thermalDiffTermY, &
                          thermalDiffTermZ
      DOUBLE PRECISION :: ropsme, ropsmn, ropsmt

      DOUBLE PRECISION :: addtermx, &
                          addtermy, addtermz
      DOUBLE PRECISION :: massMobilityTermNoDragX, &
                          massMobilityTermNoDragY, &
                          massMobilityTermNoDragZ
      DOUBLE PRECISION :: DiTE_H, DiTE_A, DiTN_H, DiTN_A, DiTT_H, DiTT_A
      DOUBLE PRECISION :: DijE_H, DijE_A, DijN_H, DijN_A, DijT_H, DijT_A
      DOUBLE PRECISION :: DijFE_H, DijFE_A, DijFN_H, DijFN_A, DijFT_H, DijFT_A
!-----------------------------------------------
!     Function subroutines
!-----------------------------------------------

      DO M = 1, SMAX
        DO 200 IJK = ijkstart3, ijkend3
          I = I_OF(IJK)
          J = J_OF(IJK)
          K = K_OF(IJK)

          IF ( FLUID_AT(IJK) ) THEN
            Mi = (PI/6.d0)*D_P(IJK,M)**3 * RO_S(IJK,M)
            Ni = ROP_s(IJK,M) / Mi

            IJKE = EAST_OF(IJK)
            IJKN = NORTH_OF(IJK)
            IJKT = TOP_OF(IJK)

! mixture density and temperature evaluated at cell faces
            ropsE = AVG_X(ROP_S(IJK,MMAX),ROP_S(IJKE,MMAX),I)
            ropsN = AVG_Y(ROP_S(IJK,MMAX),ROP_S(IJKN,MMAX),J)
            ropsT = AVG_Z(ROP_S(IJK,MMAX),ROP_S(IJKT,MMAX),K)

            ThetaE = AVG_X(THETA_M(IJK,MMAX),THETA_M(IJKE,MMAX),I)
            ThetaN = AVG_Y(THETA_M(IJK,MMAX),THETA_M(IJKN,MMAX),J)
            ThetaT = AVG_Z(THETA_M(IJK,MMAX),THETA_M(IJKT,MMAX),K)

! Thermal diffusion evaluated at cell faces (all used transport coef.
! will be evaluated this way)
            DiTE_H = AVG_X_S(DiT(IJK,M)*ROP_S(IJK,MMAX)/Theta_m(IJK,MMAX),&
               DiT(IJKE,M)*ROP_S(IJKE,MMAX)/Theta_m(IJKE,MMAX),I)
            DiTE_A = AVG_X(DiT(IJK,M)*ROP_S(IJK,MMAX)/Theta_m(IJK,MMAX),&
                DiT(IJKE,M)*ROP_S(IJKE,MMAX)/Theta_m(IJKE,MMAX),I)

            DiTN_H = AVG_Y_S(DiT(IJK,M)*ROP_S(IJK,MMAX)/Theta_m(IJK,MMAX),&
                DiT(IJKN,M)*ROP_S(IJKN,MMAX)/Theta_m(IJKN,MMAX),J)
            DiTN_A = AVG_Y(DiT(IJK,M)*ROP_S(IJK,MMAX)/Theta_m(IJK,MMAX),&
                 DiT(IJKN,M)*ROP_S(IJKN,MMAX)/Theta_m(IJKN,MMAX),J)

            DiTT_H = AVG_Z_S(DiT(IJK,M)*ROP_S(IJK,MMAX)/Theta_m(IJK,MMAX),&
                DiT(IJKT,M)*ROP_S(IJKT,MMAX)/Theta_m(IJKT,MMAX),K)
            DiTT_A = AVG_Z(DiT(IJK,M)*ROP_S(IJK,MMAX)/Theta_m(IJK,MMAX),&
                DiT(IJKT,M)*ROP_S(IJKT,MMAX)/Theta_m(IJKT,MMAX),K)

            IF(M .eq. 1)THEN
              IF(MIN(ABS(DiTE_H),ABS(DiTE_A)) .eq. ABS(DiTE_H))THEN
                DiTE = DiTE_H
                DiT_HarmE(IJK) = .TRUE.
              ELSE
                DiTE = DiTE_A
                DiT_HarmE(IJK) = .FALSE.
              ENDIF

              IF(MIN(ABS(DiTN_H),ABS(DiTN_A)) .eq. ABS(DiTN_H))THEN
                DiTN = DiTN_H
                DiT_HarmN(IJK) = .TRUE.
              ELSE
                DiTN = DiTN_A
                DiT_HarmN(IJK) = .FALSE.
              ENDIF

              IF(MIN(ABS(DiTT_H),ABS(DiTT_A)) .eq. ABS(DiTT_H))THEN
                DiTT = DiTT_H
                DiT_HarmT(IJK) = .TRUE.
              ELSE
                DiTT = DiTT_A
                DiT_HarmT(IJK) = .FALSE.
              ENDIF
            ELSE
              IF(DiT_HarmE(IJK))THEN
                DiTE = DiTE_H
              ELSE
                DiTE = DiTE_A
              ENDIF

              IF(DiT_HarmN(IJK))THEN
                DiTN = DiTN_H
              ELSE
                DiTN = DiTN_A
              ENDIF

              IF(DiT_HarmT(IJK))THEN
                DiTT = DiTT_H
              ELSE
                DiTT = DiTT_A
              ENDIF
            ENDIF   ! end if/else M=1

! initializing variables for summation over L
            ordinDiffTermX = ZERO
            ordinDiffTermY = ZERO
            ordinDiffTermZ = ZERO
            massMobilityTermX = ZERO
            massMobilityTermY = ZERO
            massMobilityTermZ = ZERO
            massMobilityTermXvelUpdate = ZERO
            massMobilityTermYvelUpdate = ZERO
            massMobilityTermZvelUpdate = ZERO
            addtermx = ZERO
            addtermy = ZERO
            addtermz = ZERO
            massMobilityTermNoDragX = ZERO
            massMobilityTermNoDragY = ZERO
            massMobilityTermNoDragZ = ZERO

            DO L = 1, SMAX
              Mj  = (PI/6.d0)*D_P(IJK,L)**3 * RO_S(IJK,L)

              NjC = ROP_s(IJK,L) / Mj
              NjE = ROP_S(IJKE,L) / Mj
              NjN = ROP_S(IJKN,L) / Mj
              NjT = ROP_S(IJKT,L) / Mj

              IF((ROP_S(IJK,MMAX)/RO_S(IJK,M) > DIL_EP_S) .and. &
                (ROP_S(IJKE,MMAX)/RO_S(IJK,M) > DIL_EP_S))THEN
                DijE_H = AVG_X_S(Dij(IJK,M,L)*Mi*Mj/ROP_S(IJK,MMAX),&
                   Dij(IJKE,M,L)*Mi*Mj/ROP_S(IJKE,MMAX),I)
                DijE_A = AVG_X(Dij(IJK,M,L)*Mi*Mj/ROP_S(IJK,MMAX),&
                   Dij(IJKE,M,L)*Mi*Mj/ROP_S(IJKE,MMAX),I)

                IF(M .eq. 1)THEN
                  IF(MIN(ABS(DijE_H),ABS(DijE_A)) .eq. ABS(DijE_H))THEN
                    DijE = DijE_H
                    Dij_HarmE(IJK,L) = .TRUE.
                  ELSE
                    DijE = DijE_A
                    Dij_HarmE(IJK,L) = .FALSE.
                  ENDIF
                ELSE
                  IF(Dij_HarmE(IJK,L))THEN
                    DijE = DijE_H
                  ELSE
                    DijE = DijE_A
                  ENDIF
                ENDIF
              ELSE
                DijE = ZERO
              ENDIF

              IF((ROP_S(IJK,MMAX)/RO_S(IJK,M) > DIL_EP_S) .and. &
                 (ROP_S(IJKN,MMAX)/RO_S(IJK,M) > DIL_EP_S))THEN
                DijN_H = AVG_Y_S(Dij(IJK,M,L)*Mi*Mj/ROP_S(IJK,MMAX),&
                  Dij(IJKN,M,L)*Mi*Mj/ROP_S(IJKN,MMAX),J)
                DijN_A = AVG_Y(Dij(IJK,M,L)*Mi*Mj/ROP_S(IJK,MMAX),&
                  Dij(IJKN,M,L)*Mi*Mj/ROP_S(IJKN,MMAX),J)
                IF(M .eq. 1)THEN
                  IF(MIN(ABS(DijN_H),ABS(DijN_A)) .eq. ABS(DijN_H))THEN
                    DijN = DijN_H
                    Dij_HarmN(IJK,L) = .TRUE.
                  ELSE
                    DijN = DijN_A
                    Dij_HarmN(IJK,L) = .FALSE.
                  ENDIF
                ELSE
                  IF(Dij_HarmN(IJK,L))THEN
                    DijN = DijN_H
                  ELSE
                    DijN = DijN_A
                  ENDIF
                ENDIF
              ELSE
                DijN = ZERO
              ENDIF

              IF((ROP_S(IJK,MMAX)/RO_S(IJK,M) > DIL_EP_S) .and. &
                 (ROP_S(IJKT,MMAX)/RO_S(IJK,M) > DIL_EP_S))THEN
                DijT_H = AVG_Z_S(Dij(IJK,M,L)*Mi*Mj/ROP_S(IJK,MMAX),&
                  Dij(IJKT,M,L)*Mi*Mj/ROP_S(IJKT,MMAX),K)
                DijT_A = AVG_Z(Dij(IJK,M,L)*Mi*Mj/ROP_S(IJK,MMAX),&
                  Dij(IJKT,M,L)*Mi*Mj/ROP_S(IJKT,MMAX),K)
                IF(M .eq. 1)THEN
                  IF(MIN(ABS(DijT_H),ABS(DijT_A)) .eq. ABS(DijT_H))THEN
                    DijT = DijT_H
                    Dij_HarmT(IJK,L) = .TRUE.
                  ELSE
                    DijT = DijT_A
                    Dij_HarmT(IJK,L) = .FALSE.
                  ENDIF
                ELSE
                  IF(Dij_HarmT(IJK,L))THEN
                    DijT = DijT_H
                  ELSE
                    DijT = DijT_A
                  ENDIF
                ENDIF
              ELSE
                DijT = ZERO
              ENDIF


              DijFE_H = AVG_X_S(DijF(IJK,M,L),DijF(IJKE,M,L),I)
              DijFE_A = AVG_X(DijF(IJK,M,L),DijF(IJKE,M,L),I)
              DijFN_H = AVG_Y_S(DijF(IJK,M,L),DijF(IJKN,M,L),J)
              DijFN_A = AVG_Y(DijF(IJK,M,L),DijF(IJKN,M,L),J)
              DijFT_H = AVG_Z_S(DijF(IJK,M,L),DijF(IJKT,M,L),K)
              DijFT_A = AVG_Z(DijF(IJK,M,L),DijF(IJKT,M,L),K)

              IF(M .eq. 1)THEN
                IF(MIN(ABS(DijFE_H),ABS(DijFE_A)) .eq. ABS(DijFE_H))THEN
                  DijFE = DijFE_H
                  DijF_HarmE(IJK,L) = .TRUE.
                ELSE
                  DijFE = DijFE_A
                  DijF_HarmE(IJK,L) = .FALSE.
                ENDIF

                IF(MIN(ABS(DijFN_H),ABS(DijFN_A)) .eq. ABS(DijFN_H))THEN
                  DijFN = DijFN_H
                  DijF_HarmN(IJK,L) = .TRUE.
                ELSE
                  DijFN = DijFN_A
                  DijF_HarmN(IJK,L) = .FALSE.
                ENDIF

                IF(MIN(ABS(DijFT_H),ABS(DijFT_A)) .eq. ABS(DijFT_H))THEN
                  DijFT = DijFT_H
                  DijF_HarmT(IJK,L) = .TRUE.
                ELSE
                  DijFT = DijFT_A
                  DijF_HarmT(IJK,L) = .FALSE.
                ENDIF
              ELSE
                IF(DijF_HarmE(IJK,L))THEN
                  DijFE = DijFE_H
                ELSE
                  DijFE = DijFE_A
                ENDIF

                IF(DijF_HarmN(IJK,L))THEN
                  DijFN = DijFN_H
                ELSE
                  DijFN = DijFN_A
                ENDIF

                IF(DijF_HarmT(IJK,L))THEN
                  DijFT = DijFT_H
                ELSE
                  DijFT = DijFT_A
                ENDIF
              ENDIF   ! end if/else (M=1)

              ordinDiffTermX = ordinDiffTermX + DijE * (NjE - NjC) * oDX_E(I)
              ordinDiffTermY = ordinDiffTermY + DijN * (NjN - NjC) * oDY_N(J)
              ordinDiffTermZ = ordinDiffTermZ + DijT * (NjT - NjC) * (oX_E(I)*oDZ_T(K))

              massMobilityTermX = massMobilityTermX + DijFE * FiX(IJK,L)
              massMobilityTermY = massMobilityTermY + DijFN * FiY(IJK,L)
              massMobilityTermZ = massMobilityTermZ + DijFT * FiZ(IJK,L)

              massMobilityTermXvelUpdate = massMobilityTermXvelUpdate + DijFE * FiXvel(IJK,L)
              massMobilityTermYvelUpdate = massMobilityTermYvelUpdate + DijFN * FiYvel(IJK,L)
              massMobilityTermZvelUpdate = massMobilityTermZvelUpdate + DijFT * FiZvel(IJK,L)

              massMobilityTermNoDragX = massMobilityTermNoDragX + DijFE * FiMinusDragX(IJK,L)
              massMobilityTermNoDragY = massMobilityTermNoDragY + DijFN * FiMinusDragY(IJK,L)
              massMobilityTermNoDragZ = massMobilityTermNoDragZ + DijFT * FiMinusDragZ(IJK,L)

            ENDDO   ! end do (l=1,smax)


            thermalDiffTermX = DiTE * ( THETA_M(IJKE,MMAX) - THETA_M(IJK,MMAX) )  * oDX_E(I)
            thermalDiffTermY = DiTN * ( THETA_M(IJKN,MMAX) - THETA_M(IJK,MMAX) )  * oDY_N(J)
            thermalDiffTermZ = DiTT * ( THETA_M(IJKT,MMAX) - THETA_M(IJK,MMAX) )  * (oX_E(I)*oDZ_T(K))

            JoiX(IJK,M) = -(ordinDiffTermX + thermalDiffTermX + massMobilityTermX)
            JoiY(IJK,M) = -(ordinDiffTermY + thermalDiffTermY + massMobilityTermY)
            JoiZ(IJK,M) = -(ordinDiffTermZ + thermalDiffTermZ + massMobilityTermZ)

            JoiMinusDragX(IJK,M) = (ordinDiffTermX + thermalDiffTermX + massMobilityTermNoDragX)
            JoiMinusDragY(IJK,M) = (ordinDiffTermY + thermalDiffTermY + massMobilityTermNoDragY)
            JoiMinusDragZ(IJK,M) = (ordinDiffTermZ + thermalDiffTermZ + massMobilityTermNoDragZ)

            ropsME=AVG_X(ROP_S(IJK,M),ROP_S(IJKE,M),I)
            ropsMN=AVG_Y(ROP_S(IJK,M),ROP_S(IJKN,M),J)
            ropsMT=AVG_Z(ROP_S(IJK,M),ROP_S(IJKT,M),K)

            DELTAU(IJK,M) = -(ordinDiffTermX+thermalDiffTermX+massMobilityTermXvelupdate)
            DELTAV(IJK,M) = -(ordinDiffTermY+thermalDiffTermY+massMobilityTermYvelupdate)
            DELTAW(IJK,M) = -(ordinDiffTermZ+thermalDiffTermz+massMobilityTermZvelUpdate)

! set fluxes to zero in case very dilute conditions
            EPSA = AVG_X(ROP_S(IJK,M),ROP_S(IJKE,M),I) / RO_S(IJK,M)
            IF(EPSA <= ZERO_EP_S) JoiX(IJK,M) = ZERO

            EPSA = AVG_Y(ROP_S(IJK,M),ROP_S(IJKN,M),J) / RO_S(IJK,M)
            IF(EPSA <= ZERO_EP_S) JoiY(IJK,M) = ZERO

            EPSA = AVG_Z(ROP_S(IJK,M),ROP_S(IJKT,M),K) / RO_S(IJK,M)
            IF(EPSA <= ZERO_EP_S) JoiZ(IJK,M) = ZERO

! set flux components to zero in case of walls
            IF (IP_AT_E(IJK))  JoiX(IJK,M) = ZERO
            IF (IP_AT_N(IJK))  JoiY(IJK,M) = ZERO
            IF (IP_AT_T(IJK))  JoiZ(IJK,M) = ZERO

          ELSE
            JoiX(IJK,M) = ZERO
            JoiY(IJK,M) = ZERO
            JoiZ(IJK,M) = ZERO
          ENDIF   ! end if/else (fluid_at(ijk))

 200  CONTINUE   ! outer IJK loop
      ENDDO   ! for m = 1,smax

      RETURN
      END SUBROUTINE GHDMASSFLUX


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: UpdateSpeciesVelocities                                 C
!  Purpose: Update solids velocities at celll faces based on the       C
!           formula Ui = U + Joi/(mi ni); also calculate averaged      C
!           velocities for dilute conditions.                          C
!                                                                      C
!  Author: S. Benyahia                              Date: 6-MAR-09     C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!     C. Hrenya handnotes and Garzo, Hrenya, Dufty papers (PRE, 2007)  C
!                                                                      C
!  Variables modified: solid species velocity components: Us, Vs, Ws   C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE UpdateSpeciesVelocities ()

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      USE fldvar
      USE indices
      USE is
      USE drag
      USE visc_s
      USE ghdtheory
      USE physprop
      USE run
      USE constant
      USE toleranc
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Index
      INTEGER :: IJK, I, J, K
! Index
      INTEGER :: IJKE, IJKN, IJKT, IJKW, IJKS, IJKB, IMJK, IJMK, &
                 IJKM
! Solids phase
      INTEGER :: S
! species density at cell faces
      integer :: kk, maxFluxS
      double precision :: epgN, rogN, mugN, Vg
      double precision :: Ur(smax), vrelSq(smax), vel, velup(smax)
      double precision :: maxFlux
      double precision :: rosN(smax), dp(smax)
      double precision :: DijN(smax,smax), JoiM(smax), &
                          DijN_H(smax,smax), DijN_A(smax,smax)
      double precision :: beta_cell(smax), beta_ij_cell(smax,smax)

      integer :: ntrial
      double precision tolx, tolf
!-----------------------------------------------
!     Function subroutines
!-----------------------------------------------

! for Newton method
      ntrial = 100
      tolx = 1d-14
      tolf = 1d-14

      DO 200 IJK = ijkstart3, ijkend3
        I = I_OF(IJK)
        J = J_OF(IJK)
        K = K_OF(IJK)

        IF ( FLUID_AT(IJK) ) THEN

          IJKE = EAST_OF(IJK)
          IJKW = WEST_OF(IJK)
          IJKN = NORTH_OF(IJK)
          IJKS = SOUTH_OF(IJK)
          IJKT = TOP_OF(IJK)
          IJKB = BOTTOM_OF(IJK)
          IMJK = IM_OF(IJK)
          IJMK = JM_OF(IJK)
          IJKM = KM_OF(IJK)

! First Compute Us
! ---------------------------------------------------------------->>>
          IF(.NOT.IP_AT_E(IJK) .OR. .NOT.SIP_AT_E(IJK)) THEN
            IF(RO_g0 == ZERO) THEN ! special case of a granular flow
              do s = 1, smax
                rosN(s) = AVG_X(ROP_S(IJK,s),ROP_S(IJKE,s),I)/ RO_S(IJK,s) ! this is ep_s
                if(rosN(s) > zero_ep_s) then
                  u_s(ijk,s) = u_s(ijk,mmax) + JoiX(IJK,s)/(rosN(s)*ro_s(IJK,s))
                else
                  u_s(ijk,s) = u_s(ijk,mmax) ! case of zero flux
                endif
              enddo
            ELSE
              Vg   = U_g(ijk) - u_s(ijk,mmax)
              epgN = AVG_X(EP_g(IJK),EP_g(IJKE),I)
              rogN = AVG_X(ROP_g(IJK),ROP_g(IJKE),I)
              mugN = AVG_X(MU_g(IJK),MU_g(IJKE),I)

              do s = 1, smax
                Ur(s) = u_g(ijk)-u_s(ijk,s)
! vrel must not include Ur, which is being solved for iteratively and must be updated.
                vrelSq(s) = (v_g(ijk)-v_s(ijk,s))**2 + (w_g(ijk)-w_s(ijk,s))**2
                rosN(s) = AVG_X(ROP_S(IJK,s),ROP_S(IJKE,s),I)
                velup(s) = 0.d0
                beta_cell(s) = beta_cell_X(IJK,s)
                dp(s)   = D_P(IJK,s)

                IF(DRAG_TYPE_ENUM .eq. HYS)THEN
                  JoiM(s) = DELTAU(IJK,s)
                ELSE
                  JoiM(s) = JoiMinusDragX(ijk,s)
                ENDIF

                do kk = 1, smax
                  DijN_H(s,kk) = AVG_X_S(DijF(IJK,s,kk),DijF(IJKE,s,kk),I)
                  DijN_A(s,kk) = AVG_X(DijF(IJK,s,kk),DijF(IJKE,s,kk),I)
                  if(DijF_HarmE(IJK,kk))THEN
                    DijN(s,kk) = DijN_H(s,kk)
                  ELSE
                    DijN(s,kk) = DijN_A(s,kk)
                  ENDIF
                  if(s .eq. kk)then
                    beta_ij_cell(s,kk)=0.d0
                  else
                    beta_ij_cell(s,kk)=beta_ij_cell_X(IJK,s,kk)
                  endif
                enddo
              enddo

              IF(DRAG_TYPE_ENUM .eq. HYS)THEN
                vel=U_S(IJK,MMAX)
                CALL VELOCITY_UPDATE(velup, smax, rosN, DijN, &
                  JoiM, beta_cell, beta_ij_cell,vel)
              ELSE
                CALL UrNEWT(ntrial, Ur, smax, ijk, tolx, tolf, &
                  epgN, rogN, mugN, Vg, vrelSq, rosN, dp, DijN, JoiM)
              ENDIF

! species velocity and flux update.
              do s = 1, smax
                IF(DRAG_TYPE_ENUM .eq. HYS)THEN
                  U_S(IJK,s)=velup(s)
                  JoiX(IJK,s) = rosN(s) * (u_s(ijk,s)-u_s(ijk,mmax))
                ELSE
                  u_s(ijk,s) =  u_g(ijk) - Ur(s)
                  JoiX(IJK,s) = rosN(s) * (u_s(ijk,s)-u_s(ijk,mmax))
                ENDIF
              enddo
            ENDIF ! for a granular case (no gas and no drag)
          ENDIF   ! end if (.not.ip_at_e or .not.sip_at_e)

          if(smax==2) then ! only for binary, how to implement for smax > 2?
            JoiX(IJK,2)=-JoiX(IJK,1)
          elseif(smax > 2) then
            maxFlux = JoiX(IJK,1)
            maxFluxS = 1
            do s = 2, smax  ! finding species with maximum flux in a cell
              if( abs(JoiX(IJK,s)) > abs(maxFlux) ) then
                maxFlux = JoiX(IJK,s)
                maxFluxS = s
              endif
            enddo
            JoiX(IJK,maxFluxS) = 0d0 ! reset max. flux to zero
            do s = 1, smax ! re-calc species with max. flux to satisfy SUM(fluxes) = 0
              if(s /= maxFluxS) JoiX(IJK,maxFluxS) = JoiX(IJK,maxFluxS) - JoiX(IJK,s)
            enddo
          endif
! End Compute Us
! ----------------------------------------------------------------<<<



! Now Compute Vs
! ---------------------------------------------------------------->>>
          IF (.NOT.IP_AT_N(IJK) .OR. .NOT.SIP_AT_N(IJK)) THEN
            IF(RO_g0 == ZERO) THEN ! special case of a granular flow
              do s = 1, smax
                rosN(s) = AVG_Y(ROP_S(IJK,s),ROP_S(IJKN,s),J)/ RO_S(IJK,s) ! this is ep_s
                if(rosN(s) > zero_ep_s) then
                  v_s(ijk,s) = v_s(ijk,mmax) + JoiY(IJK,s)/(rosN(s)*ro_s(IJK,s))
                else
                   v_s(ijk,s) = v_s(ijk,mmax) ! case of zero flux
                endif
              enddo
            ELSE
              Vg   = V_g(ijk) - v_s(ijk,mmax)
              epgN = AVG_Y(EP_g(IJK),EP_g(IJKN),J)
              rogN = AVG_Y(ROP_g(IJK),ROP_g(IJKN),J)
              mugN = AVG_Y(MU_g(IJK),MU_g(IJKN),J)

              do s = 1, smax
                Ur(s) = v_g(ijk)-v_s(ijk,s)
! vrel must not include Ur, which is being solved for iteratively and must be updated.
                vrelSq(s) = (u_g(ijk)-u_s(ijk,s))**2 + (w_g(ijk)-w_s(ijk,s))**2
                rosN(s) = AVG_Y(ROP_S(IJK,s),ROP_S(IJKN,s),J)
                velup(s) = 0.d0
                beta_cell(s) = beta_cell_Y(IJK,s)
                dp(s)   = D_P(IJK,s)

                IF(DRAG_TYPE_ENUM .eq. HYS)THEN
                  JoiM(s) = DELTAV(IJK,s)
                ELSE
                  JoiM(s) = JoiMinusDragY(ijk,s)
                ENDIF

                do kk = 1, smax

                  DijN_H(s,kk) = AVG_Y_S(DijF(IJK,s,kk),DijF(IJKN,s,kk),J)
                  DijN_A(s,kk) = AVG_Y(DijF(IJK,s,kk),DijF(IJKN,s,kk),J)

                  if(DijF_HarmN(IJK,kk))THEN
                    DijN(s,kk) = DijN_H(s,kk)
                  ELSE
                    DijN(s,kk) = DijN_A(s,kk)
                  ENDIF

                  if(s .eq. kk)then
                     beta_ij_cell(s,kk)=0.d0
                  else
                    beta_ij_cell(s,kk)=beta_ij_cell_Y(IJK,s,kk)
                  endif
                enddo
              enddo

              IF(DRAG_TYPE_ENUM .eq. HYS)THEN
                vel=V_S(IJK,MMAX)
                CALL VELOCITY_UPDATE(velup, smax, rosN, DijN, JoiM, &
                  beta_cell, beta_ij_cell,vel)
              ELSE
                CALL UrNEWT(ntrial, Ur, smax, ijk, tolx, tolf, &
                  epgN, rogN, mugN, Vg, vrelSq, rosN, dp, DijN, JoiM)
              ENDIF

! species velocity and flux update
              do s = 1, smax
                IF(DRAG_TYPE_ENUM .eq. HYS)THEN
                  V_S(IJK,s)=velup(s)
                  JoiY(IJK,s) = rosN(s) * (v_s(ijk,s)-v_s(ijk,mmax))
                ELSE
                  v_s(ijk,s) =  v_g(ijk) - Ur(s)
                  JoiY(IJK,s) = rosN(s) * (v_s(ijk,s)-v_s(ijk,mmax))
                ENDIF
              enddo
            ENDIF ! for a granular case (no gas and no drag)
          ENDIF   ! end if (.not.ip_at_n or .not.sip_at_n)

          if(smax==2) then ! only for binary, how to implement for smax > 2?
            JoiY(IJK,2)=-JoiY(IJK,1)
          elseif(smax > 2) then
            maxFlux = JoiY(IJK,1)
            maxFluxS = 1
            do s = 2, smax  ! finding species with maximum flux in a cell
              if( abs(JoiY(IJK,s)) > abs(maxFlux) ) then
                maxFlux = JoiY(IJK,s)
                maxFluxS = s
              endif
            enddo
            JoiY(IJK,maxFluxS) = 0d0 ! reset max. flux to zero
            do s = 1, smax ! re-calc species with max. flux to satisfy SUM(fluxes) = 0
              if(s /= maxFluxS) JoiY(IJK,maxFluxS) = JoiY(IJK,maxFluxS) - JoiY(IJK,s)
            enddo
          endif
! End Compute Vs
! ----------------------------------------------------------------<<<


! Finaly Compute Ws
! ---------------------------------------------------------------->>>
          IF(.NOT.NO_K .AND. (.NOT.IP_AT_T(IJK) .OR. .NOT.SIP_AT_T(IJK))) THEN
            IF(RO_g0 == ZERO) THEN ! special case of a granular flow
              do s = 1, smax
                rosN(s) = AVG_Z(ROP_S(IJK,s),ROP_S(IJKT,s),K)/ RO_S(IJK,s) ! this is ep_s
                if(rosN(s) > zero_ep_s) then
                  w_s(ijk,s) = w_s(ijk,mmax) + JoiZ(IJK,s)/(rosN(s)*ro_s(IJK,s))
                else
                  w_s(ijk,s) = w_s(ijk,mmax) ! case of zero flux
                endif
              enddo
            ELSE
              Vg   = W_g(ijk) - W_s(ijk,mmax)
              epgN = AVG_Z(EP_g(IJK),EP_g(IJKT),K)
              rogN = AVG_Z(ROP_g(IJK),ROP_g(IJKT),K)
              mugN = AVG_Z(MU_g(IJK),MU_g(IJKT),K)

              do s = 1, smax
                Ur(s) = w_g(ijk)-w_s(ijk,s)
! vrel must not include Ur, which is being solved for iteratively and must be updated.
                vrelSq(s) = (u_g(ijk)-u_s(ijk,s))**2 + (v_g(ijk)-v_s(ijk,s))**2
                rosN(s) = AVG_Z(ROP_S(IJK,s),ROP_S(IJKT,s),K)
                velup(s) = 0.d0
                beta_cell(s) = beta_cell_Z(IJK,s)
                dp(s)   = D_P(IJK,s)

                IF(DRAG_TYPE_ENUM .eq. HYS)THEN
                  JoiM(s) = DELTAW(IJK,s)
                ELSE
                  JoiM(s) = JoiMinusDragZ(ijk,s)
                ENDIF

                do kk = 1, smax
                  DijN_H(s,kk) = AVG_Z_S(DijF(IJK,s,kk),DijF(IJKT,s,kk),K)
                  DijN_A(s,kk) = AVG_Z(DijF(IJK,s,kk),DijF(IJKT,s,kk),K)
                  if(DijF_HarmT(IJK,kk))THEN
                    DijN(s,kk) = DijN_H(s,kk)
                  ELSE
                    DijN(s,kk) = DijN_A(s,kk)
                  ENDIF
                  if(s .eq. kk)then
                     beta_ij_cell(s,kk)=0.d0
                  else
                    beta_ij_cell(s,kk)=beta_ij_cell_Z(IJK,s,kk)
                  endif
                enddo
              enddo

              IF(DRAG_TYPE_ENUM .eq. HYS)THEN
                vel=W_S(IJK,MMAX)
                CALL VELOCITY_UPDATE(velup, smax, rosN, DijN, &
                  JoiM, beta_cell, beta_ij_cell,vel)
              ELSE
                CALL UrNEWT(ntrial, Ur, smax, ijk, tolx, tolf, &
                  epgN, rogN, mugN, Vg, vrelSq, rosN, dp, DijN, JoiM)
              ENDIF

! species velocity and flux update
              do s = 1, smax
                IF(DRAG_TYPE_ENUM .eq. HYS)THEN
                  W_S(IJK,s)=velup(s)
                  JoiZ(IJK,s) = rosN(s) * (w_s(ijk,s)-w_s(ijk,mmax))
                ELSE
                  w_s(ijk,s) =  w_g(ijk) - Ur(s)
                  JoiZ(IJK,s) = rosN(s) * (w_s(ijk,s)-w_s(ijk,mmax))
                ENDIF
              enddo
            ENDIF ! for a granular case (no gas and no drag)
          ENDIF   ! end if (.not.ip_at_t or .not.sip_at_t)

          if(smax==2) then ! only for binary, how to implement for smax > 2?
            JoiZ(IJK,2)=-JoiZ(IJK,1)
          elseif(smax > 2) then
            maxFlux = JoiZ(IJK,1)
            maxFluxS = 1
            do s = 2, smax  ! finding species with maximum flux in a cell
              if( abs(JoiZ(IJK,s)) > abs(maxFlux) ) then
                maxFlux = JoiZ(IJK,s)
                maxFluxS = s
              endif
            enddo
            JoiZ(IJK,maxFluxS) = 0d0 ! reset max. flux to zero
            do s = 1, smax ! re-calc species with max. flux to satisfy SUM(fluxes) = 0
              if(s /= maxFluxS) JoiZ(IJK,maxFluxS) = JoiZ(IJK,maxFluxS) - JoiZ(IJK,s)
            enddo
          endif
! End Compute Ws
! ----------------------------------------------------------------<<<

! if .not. fluid_at(ijk), do nothing (keep old values of velocities).

          ENDIF     ! Fluid_at
 200  CONTINUE     ! outer IJK loop

      RETURN
      END SUBROUTINE UpdateSpeciesVelocities


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: UrNEWT                                                  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      subroutine UrNEWT(ntrial, x, s, ijk, tolx, tolf, epgN, rogN, &
                        mugN, Vg, vrelSq, rosN, dp, DijN, JoiM)

      Implicit NONE

!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      integer, intent(in) :: s
      integer, intent(in) :: ijk
      integer, intent(in) :: ntrial
      double precision, intent(in) :: tolx, tolf

      DOUBLE PRECISION, intent(inout) :: X(s)
      double precision, intent(in) :: epgN, rogN, mugN, Vg
      double precision, intent(in) :: vrelSq(s), rosN(s), dp(s)
      double precision, intent(in) :: DijN(s,s), JoiM(s)
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      INTEGER :: NP
      PARAMETER (NP=15)  ! solves up to NP variable/equations;
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: I, K, INDX(s)
      DOUBLE PRECISION :: D, ERRF, ERRX, FJAC(s,s), FVEC(s), P(s)
!-----------------------------------------------

      DO K = 1, NTRIAL
        CALL Ur_JACOBI_EVAL(X, s, ijk, FVEC, FJAC, epgN, rogN, mugN, Vg, &
                             vrelSq, rosN, dp, DijN, JoiM)
        ERRF = 0d0
        DO I = 1, s
             ERRF = ERRF + DABS(FVEC(I))
        ENDDO

        IF(ERRF <= TOLF) RETURN
        DO I = 1, s
            P(I) = -FVEC(I) ! RHS of linear equations.
            IF(rosN(I) .eq. 0.d0)THEN
              P(I) = 0.d0
            ENDIF
        ENDDO

        CALL LUDCMP(fjac, s, NP, indx, d, 'UrNewt')  ! solve system of s linear equations using
        CALL LUBKSB(fjac, s, NP, indx, p)  ! LU decomposition

        ERRX = 0d0
        DO I = 1, s
            ERRX = ERRX + DABS(P(I))
            X(I) = X(I) + P(I)
        ENDDO
        IF(ERRX <= TOLX) RETURN
      ENDDO  ! for K trials

      RETURN
      END SUBROUTINE UrNEWT


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: Ur_JACOBI_EVAL                                          C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE Ur_JACOBI_EVAL(X, s, ijk, FVEC, FJAC, epgN, rogN, &
                                mugN, Vg, vrelSq, rosN, dp, DijN, JoiM)
      Implicit NONE

!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      integer, intent(in) :: s
      integer, intent(in) :: ijk
      double precision, intent(in) :: X(s)
! vector function and matrix jacobian
      DOUBLE PRECISION, INTENT(OUT) :: FVEC(s), FJAC(s,s)
      DOUBLE PRECISION, INTENT(IN) :: epgN, rogN, mugN, Vg
      DOUBLE PRECISION, INTENT(IN) :: vrelSq(s), rosN(s), dp(s)
      DOUBLE PRECISION, INTENT(IN) :: DijN(s,s), JoiM(s)
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      double precision :: pi
      parameter (pi=3.14159265458979323846d0)
      double precision :: one
      parameter (one=1.d0)
      double precision :: zero
      parameter (zero=0.d0)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J
! various quantities
      DOUBLE PRECISION :: RE_G, C_d, DgA, Vi, vrel
      DOUBLE PRECISION :: FgsOni(s), dFgsdVi(s), sum(s)
!-----------------------------------------------

      DO i = 1, s
        vrel = dsqrt(vrelSq(i) + x(i)**2)
        Vi = pi * dp(i)**3 / 6d0
        RE_G = dp(i)*vrel*rogN/mugN
        IF(Re_G <= 1000D0 .and. Re_G> zero)THEN
           C_d = (24.D0/Re_G) * (ONE + 0.15D0 * Re_G**0.687D0)
        ELSE
           C_d = 0.44D0
        ENDIF
        DgA = 0.75D0 * C_d * vrel * rogN / (epgN**2.65D0 * dp(i))
        FgsOni(i) = DgA * Vi  ! this is the wen-Yu drag
        IF(vrel == ZERO)THEN
           dFgsdVi(i) = ZERO
        ELSEIF(Re_G <= 1000D0)THEN
           dFgsdVi(i) = 1.8549d0*mugN*Re_G**0.687d0*Vi*dabs(x(i))/ &
                            (dp(i)**2*epgN**2.65d0*vrel**2)
        ELSE
           dFgsdVi(i) = 0.33d0*rogN*Vi*dabs(x(i))/ &
                             (dp(i)*epgN**2.65d0*vrel)
        ENDIF
      ENDDO

! Start computing values of the function FVEC(i)
      DO i = 1, s
         sum(i) = zero
         DO j = 1, s
            sum(i) = sum(i) + DijN(i,j) * FgsOni(j) * x(j)
         ENDDO
      ENDDO
      DO i = 1, s
         FVEC(i) = sum(i) - rosN(i)*x(i) + rosN(i)*Vg + JoiM(i)
      ENDDO

!
! Evaluation of functions FVEC(i) is done above, and their Jacobian is computed below.
!
      DO i = 1, s
         DO j = 1, s
            if(j == i) THEN
               FJAC(i,j) = DijN(i,i) * (FgsOni(i) + dFgsdVi(i)*x(i)) - rosN(i)
               IF(rosN(i) .eq. 0.d0)THEN
                 FJAC(i,j) = 1.d0
               ENDIF
            else
               FJAC(i,j) = DijN(i,j) * (FgsOni(j) + dFgsdVi(j)*x(j))
            endif
         ENDDO ! for j
      ENDDO ! for i
!
! End of computing jacobian (J_ij).
!
      RETURN
      END SUBROUTINE Ur_JACOBI_EVAL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: Velocity_update                                         C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE VELOCITY_UPDATE(FVEC, s, rosi, Diji, Joii, &
                                 beta_celli, beta_ij_celli, &
                                 velocity)
      USE toleranc
      USE fldvar
      Implicit NONE

!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      integer, intent (in) :: s
! various quantities
      DOUBLE PRECISION, INTENT(IN) :: rosi(s), Diji(s,s), &
                                      Joii(s), beta_celli(s), &
                                      beta_ij_celli(s,s), velocity
! vector function
      DOUBLE PRECISION, INTENT(OUT) :: FVEC(s)
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
      integer :: NP
      PARAMETER (NP=15)  ! solves up to NP variable/equations;
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: i, k, l
      INTEGER :: indx(s)
      DOUBLE PRECISION :: D
! matrix jacobian
      DOUBLE PRECISION :: FJAC(s,s)
!-----------------------------------------------

      DO i=1,s
         FVEC(i)= 0.d0
         DO k=1,s
           FJAC(i,k) = 0.d0
         ENDDO
      ENDDO

      DO i=1,s

        IF(rosi(i) .ne. 0.d0)THEN
           FVEC(i) = rosi(i)*velocity+Joii(i)
        ELSE
           FVEC(i) = 0.d0
        ENDIF

        DO l=1,s

          if(i .eq. l)then
            FJAC(i,l)=FJAC(i,l)+rosi(i)
            IF(rosi(i) .eq. 0.d0)THEN   ! This is done to avoid a singular matrix if ros(i) is ZERO
              FJAC(i,l) = 1.d0
            ENDIF
          endif

          FJAC(i,l)=FJAC(i,l)-Diji(i,l)*beta_celli(l)

          DO k=1,s

           if(l .ne. k)then
              FJAC(i,k)=FJAC(i,k)+Diji(i,l)*beta_ij_celli(l,k)
           endif

          enddo

        enddo
      enddo

      CALL LUDCMP(fjac, s, NP, indx, d, 'velocity_update')  ! solve system of s linear equations using
      CALL LUBKSB(fjac, s, NP, indx, FVEC)  ! LU decomposition

      RETURN

      END SUBROUTINE VELOCITY_UPDATE

