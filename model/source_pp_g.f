!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_Pp_g                                             C
!  Purpose: Determine source terms for Pressure correction equation.   C
!                                                                      C
!  Notes: The off-diagonal coefficients are positive. The center       C
!         coefficient and the source vector are negative. See          C
!         conv_Pp_g                                                    C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

SUBROUTINE SOURCE_PP_G(A_M, B_M, B_MMAX)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc, ONLY: SMALL_NUMBER, ONE, ZERO, UNDEFINED, IJK_P_G, DIMENSION_3, DIMENSION_M
      USE compar, ONLY: IJKSTART3, IJKEND3
      USE cutcell, ONLY: CARTESIAN_GRID, A_UPG_E, A_VPG_N, A_WPG_T
      USE eos, ONLY: DROODP_G
      USE fldvar, ONLY: U_G, V_G, W_G, U_S, V_S, W_S, ROP_G, ROP_S, ROP_GO, ROP_SO, RO_G, P_G, EP_G
      USE geometry, ONLY: VOL
      USE matrix, ONLY: E, W, N, S, T, B
      USE pgcor, ONLY: D_E, D_N, D_T
      USE physprop, ONLY: CLOSE_PACKED, MMAX, RO_G0
      USE run, ONLY: SHEAR, ODT, UNDEFINED_I
      USE rxns, ONLY: SUM_R_G, SUM_R_S
      USE ur_facs, ONLY: UR_FAC
      USE vshear, ONLY: VSH
      USE xsi_array, ONLY: LOCK_XSI_ARRAY, UNLOCK_XSI_ARRAY

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! maximum term in b_m expression
      DOUBLE PRECISION, INTENT(INOUT) :: B_mmax(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! solids phase index
      INTEGER :: M
! Indices
      INTEGER :: IJK, IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
! under relaxation factor for pressure
      DOUBLE PRECISION fac
! terms of bm expression
      DOUBLE PRECISION bma, bme, bmw, bmn, bms, bmt, bmb, bmr
! loezos: used for including shearing
      integer :: incr
! error message
      CHARACTER(LEN=80) :: LINE(1)
! temporary use of global arrays:
! xsi_array: convection weighting factors
!      DOUBLE PRECISION XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),&
!                       XSI_t(DIMENSION_3)
!-----------------------------------------------
      call lock_xsi_array
! loezos
! update to true velocity
      IF (SHEAR) THEN
!!$omp parallel do private(IJK)
         DO IJK = IJKSTART3, IJKEND3
            IF (FLUID_AT(IJK)) THEN
               V_G(IJK)=V_G(IJK)+VSH(IJK)
            ENDIF
         ENDDO
      ENDIF

! Calculate convection-diffusion fluxes through each of the faces

!$omp parallel default(none) &
!$omp          private(IJK, IMJK, IJMK, IJKM, M, bma,bme,bmw,bmn,bms,bmt,bmb,bmr,line)  &
!$omp          shared(ijkstart3,ijkend3,cartesian_grid,rop_g,rop_go,rop_s,rop_so,vol,odt,u_g,v_g,w_g,u_s,v_s,w_s,b_m, &
!$omp                 b_mmax,d_e,d_n,d_t,a_m,a_upg_e,a_vpg_n,a_wpg_t,mmax,close_packed,sum_r_s,sum_r_g,ro_g0)
!$omp do
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)

            bma = (ROP_G(IJK)-ROP_GO(IJK))*VOL(IJK)*ODT
            bme = A_M(IJK,E,0)*U_G(IJK)
            bmw = A_M(IJK,W,0)*U_G(IMJK)
            bmn = A_M(IJK,N,0)*V_G(IJK)
            bms = A_M(IJK,S,0)*V_G(IJMK)
            bmt = A_M(IJK,T,0)*W_G(IJK)
            bmb = A_M(IJK,B,0)*W_G(IJKM)
            bmr = SUM_R_G(IJK)*VOL(IJK)
            B_M(IJK,0) = -((-(bma + bme - bmw + bmn - bms + bmt - bmb ))+ bmr )
            B_MMAX(IJK,0) = max(abs(bma), abs(bme), abs(bmw), abs(bmn), abs(bms), abs(bmt), abs(bmb), abs(bmr) )

            A_M(IJK,E,0) = A_M(IJK,E,0)*D_E(IJK,0)
            A_M(IJK,W,0) = A_M(IJK,W,0)*D_E(IMJK,0)
            A_M(IJK,N,0) = A_M(IJK,N,0)*D_N(IJK,0)
            A_M(IJK,S,0) = A_M(IJK,S,0)*D_N(IJMK,0)
            A_M(IJK,T,0) = A_M(IJK,T,0)*D_T(IJK,0)
            A_M(IJK,B,0) = A_M(IJK,B,0)*D_T(IJKM,0)

            IF(CARTESIAN_GRID) THEN
               A_M(IJK,E,0) = A_M(IJK,E,0) * A_UPG_E(IJK)
               A_M(IJK,W,0) = A_M(IJK,W,0) * A_UPG_E(IMJK)
               A_M(IJK,N,0) = A_M(IJK,N,0) * A_VPG_N(IJK)
               A_M(IJK,S,0) = A_M(IJK,S,0) * A_VPG_N(IJMK)
               A_M(IJK,T,0) = A_M(IJK,T,0) * A_WPG_T(IJK)
               A_M(IJK,B,0) = A_M(IJK,B,0) * A_WPG_T(IJKM)
            ENDIF

! solids phase terms are needed for those solids phases that do not
! become close packed. these terms are used to essentially make the
! matrix equation for gas pressure correction a mixture pressure
! correction equation by adding solids phases that do not become
! close packed
            DO M = 1, MMAX
               IF (.NOT.CLOSE_PACKED(M)) THEN
                  B_M(IJK,0) = B_M(IJK,0) - &
                     ((-((ROP_S(IJK,M)-ROP_SO(IJK,M))*VOL(IJK)*ODT+&
                     A_M(IJK,E,M)*U_S(IJK,M)-A_M(IJK,W,M)*U_S(IMJK,M)+&
                     A_M(IJK,N,M)*V_S(IJK,M)-A_M(IJK,S,M)*V_S(IJMK,M)+&
                     A_M(IJK,T,M)*W_S(IJK,M)-A_M(IJK,B,M)*W_S(IJKM,M)))+&
                     SUM_R_S(IJK,M)*VOL(IJK))

                  IF(.NOT.CARTESIAN_GRID) THEN
                     A_M(IJK,E,0) = A_M(IJK,E,0) + A_M(IJK,E,M)*D_E(IJK,M)
                     A_M(IJK,W,0) = A_M(IJK,W,0) + A_M(IJK,W,M)*D_E(IMJK,M)
                     A_M(IJK,N,0) = A_M(IJK,N,0) + A_M(IJK,N,M)*D_N(IJK,M)
                     A_M(IJK,S,0) = A_M(IJK,S,0) + A_M(IJK,S,M)*D_N(IJMK,M)
                     A_M(IJK,T,0) = A_M(IJK,T,0) + A_M(IJK,T,M)*D_T(IJK,M)
                     A_M(IJK,B,0) = A_M(IJK,B,0) + A_M(IJK,B,M)*D_T(IJKM,M)
                  ELSE
                     A_M(IJK,E,0) = A_M(IJK,E,0) + A_M(IJK,E,M)*D_E(IJK,M)  * A_UPG_E(IJK)
                     A_M(IJK,W,0) = A_M(IJK,W,0) + A_M(IJK,W,M)*D_E(IMJK,M) * A_UPG_E(IMJK)
                     A_M(IJK,N,0) = A_M(IJK,N,0) + A_M(IJK,N,M)*D_N(IJK,M)  * A_VPG_N(IJK)
                     A_M(IJK,S,0) = A_M(IJK,S,0) + A_M(IJK,S,M)*D_N(IJMK,M) * A_VPG_N(IJMK)
                     A_M(IJK,T,0) = A_M(IJK,T,0) + A_M(IJK,T,M)*D_T(IJK,M)  * A_WPG_T(IJK)
                     A_M(IJK,B,0) = A_M(IJK,B,0) + A_M(IJK,B,M)*D_T(IJKM,M) * A_WPG_T(IJKM)
                  ENDIF
               ENDIF
            ENDDO

            A_M(IJK,0,0) = -(A_M(IJK,E,0)+A_M(IJK,W,0)+A_M(IJK,N,0)+A_M(IJK,S,0&
               )+A_M(IJK,T,0)+A_M(IJK,B,0))

            IF (ABS(A_M(IJK,0,0)) < SMALL_NUMBER) THEN
               IF (ABS(B_M(IJK,0)) < SMALL_NUMBER) THEN
                  A_M(IJK,0,0) = -ONE
                  B_M(IJK,0) = ZERO
               ELSEIF (RO_G0 .NE. UNDEFINED) THEN !This is an error only in incompressible flow
!!$omp             critical
                  WRITE (LINE, '(A,I6,A,I1,A,G12.5)') 'Error: At IJK = ', IJK, &
                     ' M = ', 0, ' A = 0 and b = ', B_M(IJK,0)
                  CALL WRITE_ERROR ('SOURCE_Pp_g', LINE, 1)
!!$omp             end critical
               ENDIF
            ENDIF

         ELSE   ! if/else branch .not.fluid_at(ijk)
! set the value (correction) in all wall and flow boundary cells to zero
! note the matrix coefficients and source vector should already be zero
! from the initialization of A and B but the following ensures the zero
! value.
            A_M(IJK,E,0) = ZERO
            A_M(IJK,W,0) = ZERO
            A_M(IJK,N,0) = ZERO
            A_M(IJK,S,0) = ZERO
            A_M(IJK,T,0) = ZERO
            A_M(IJK,B,0) = ZERO
            A_M(IJK,0,0) = -ONE
            B_M(IJK,0) = ZERO
         ENDIF   ! end if/else branch fluid_at(ijk)
      ENDDO    ! end do loop (ijk=ijkstart3,ijkend3)
!$omp end parallel


! loezos
      IF (SHEAR) THEN
!!$omp parallel do private(IJK)
         DO IJK = IJKSTART3, IJKEND3
            IF (FLUID_AT(IJK)) THEN
               V_G(IJK)=V_G(IJK)-VSH(IJK)
            ENDIF
         ENDDO
      ENDIF


! make correction for compressible flows
      IF (RO_G0 == UNDEFINED) THEN
         fac = UR_FAC(1)  !since p_g = p_g* + ur_fac * pp_g

! loezos
         incr=0
!         CALL CALC_XSI(DISCRETIZE(1),ROP_G,U_G,V_G,W_G,XSI_E,XSI_N,XSI_T,incr)

!!$omp    parallel do                                                     &
!!$omp&   private(IJK,I,J,K,                                       &
!!$omp&            IMJK,IJMK,IJKM,IJKE,IJKW,IJKN,IJKS,IJKT,IJKB)
         DO IJK = ijkstart3, ijkend3
            IF (FLUID_AT(IJK)) THEN

               A_M(IJK,0,0) = A_M(IJK,0,0) - &
                  fac*DROODP_G(RO_G(IJK),P_G(IJK))*&
                  EP_G(IJK)*VOL(IJK)*ODT

! Although the following is a better approximation for high speed flows because
! it considers density changes in the neighboring cells, the code runs faster
! without it for low speed flows.  The gas phase mass balance cannot be
! maintained to machine precision with the following approximation. If the
! following lines are uncommented, the calc_xsi call above should also be
! uncommented
!               IMJK = IM_OF(IJK)
!               IJMK = JM_OF(IJK)
!               IJKM = KM_OF(IJK)
!               IJKE = EAST_OF(IJK)
!               IJKW = WEST_OF(IJK)
!               IJKN = NORTH_OF(IJK)
!               IJKS = SOUTH_OF(IJK)
!               IJKT = TOP_OF(IJK)
!               IJKB = BOTTOM_OF(IJK)
!               A_M(IJK,0,0) = A_M(IJK,0,0) - fac*DROODP_G(RO_G(IJK),P_G(IJK))*EP_G(&
!                  IJK)*((ONE - XSI_E(IJK))*U_G(IJK)*AYZ(IJK)-XSI_E(IMJK)*U_G(&
!                  IMJK)*AYZ(IMJK)+(ONE-XSI_N(IJK))*V_G(IJK)*AXZ(IJK)-XSI_N(IJMK&
!                  )*V_G(IJMK)*AXZ(IJMK))

!               A_M(IJK,E,0) = A_M(IJK,E,0) - EP_G(IJKE)*fac*DROODP_G(RO_G(IJKE),P_G&
!                  (IJKE))*XSI_E(IJK)*U_G(IJK)*AYZ(IJK)
!               A_M(IJK,W,0) = A_M(IJK,W,0) + EP_G(IJKW)*fac*DROODP_G(RO_G(IJKW),P_G&
!                  (IJKW))*(ONE - XSI_E(IMJK))*U_G(IMJK)*AYZ(IMJK)
!               A_M(IJK,N,0) = A_M(IJK,N,0) - EP_G(IJKN)*fac*DROODP_G(RO_G(IJKN),P_G&
!                  (IJKN))*XSI_N(IJK)*V_G(IJK)*AXZ(IJK)
!               A_M(IJK,S,0) = A_M(IJK,S,0) + EP_G(IJKS)*fac*DROODP_G(RO_G(IJKS),P_G&
!                  (IJKS))*(ONE - XSI_N(IJMK))*V_G(IJMK)*AXZ(IJMK)
!               IF (DO_K) THEN
!                  A_M(IJK,0,0) = A_M(IJK,0,0) - fac*DROODP_G(RO_G(IJK),P_G(IJK))*&
!                     EP_G(IJK)*((ONE - XSI_T(IJK))*W_G(IJK)*AXY(IJK)-XSI_T(IJKM&
!                     )*W_G(IJKM)*AXY(IJKM))
!                  A_M(IJK,T,0) = A_M(IJK,T,0) - EP_G(IJKT)*fac*DROODP_G(RO_G(IJKT),&
!                     P_G(IJKT))*XSI_T(IJK)*W_G(IJK)*AXY(IJK)
!                  A_M(IJK,B,0) = A_M(IJK,B,0) + EP_G(IJKB)*fac*DROODP_G(RO_G(IJKB),&
!                     P_G(IJKB))*(ONE - XSI_T(IJKM))*W_G(IJKM)*AXY(IJKM)
!               ENDIF
!
            ENDIF   !end if (fluid_at(ijk))
         ENDDO    ! end do (ijk=ijkstart3,ijkend3)
      ENDIF   ! end if (ro_g0 == undefined); i.e., compressible flow


! Remove the asymmetry in matrix caused by the pressure outlet or inlet
! boundaries.  Because the P' at such boundaries is zero we may set the
! coefficient in the neighboring fluid cell to zero without affecting
! the linear equation set.
!!$omp    parallel do                                                     &
!!$omp&   private(IJK,IMJK, IPJK, IJMK, IJPK, IJKM, IJKP)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN
            IMJK = IM_OF(IJK)
            IPJK = IP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKM = KM_OF(IJK)
            IJKP = KP_OF(IJK)
! Cutting the neighbor link between fluid cell and adjacent p_flow_at cell
            if(p_flow_at(imjk)) A_m(IJK, w, 0) = ZERO
            if(p_flow_at(ipjk)) A_m(IJK, e, 0) = ZERO
            if(p_flow_at(ijmk)) A_m(IJK, s, 0) = ZERO
            if(p_flow_at(ijpk)) A_m(IJK, n, 0) = ZERO
            if(p_flow_at(ijkm)) A_m(IJK, b, 0) = ZERO
            if(p_flow_at(ijkp)) A_m(IJK, t, 0) = ZERO
         ENDIF
      ENDDO

! Specify P' to zero for incompressible flows. Check set_bc0
! for details on selection of IJK_P_g.
      IF (IJK_P_G /= UNDEFINED_I) THEN
         B_M(IJK_P_G,0) = ZERO
         A_M(IJK_P_G,:,0) = ZERO
         A_M(IJK_P_G,0,0) = -ONE
      ENDIF

      call unlock_xsi_array

      RETURN

CONTAINS

      INCLUDE 'functions.inc'

END SUBROUTINE SOURCE_PP_G
