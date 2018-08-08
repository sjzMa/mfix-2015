!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_trD_g(trD_g, IER)                                 C
!  Purpose: Calculate the trace of gas phase rate of strain tensor     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-DEC-96  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_TRD_G(TRD_G)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!     Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE geometry
      USE fldvar
      USE indices
      USE compar
      USE sendrecv
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      USE bc
      USE cutcell
      USE quadric
      USE functions
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Indices
      INTEGER          I, J, K, IJK, IMJK, IJMK, IJKM, IM
!
!                      Strain rate tensor components for mth solids phase
      DOUBLE PRECISION trD_g(DIMENSION_3)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
      DOUBLE PRECISION :: DEL_H,Nx,Ny,Nz
      DOUBLE PRECISION :: dudx,dvdy,dwdz
      DOUBLE PRECISION :: Xi,Yi,Zi,Ui,Vi,Wi,Sx,Sy,Sz
      DOUBLE PRECISION :: UW_g,VW_g,WW_g

      LOGICAL :: U_NODE_AT_E, U_NODE_AT_W
      LOGICAL :: V_NODE_AT_N, V_NODE_AT_S
      LOGICAL :: W_NODE_AT_T, W_NODE_AT_B
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
!-----------------------------------------------
!
!
!!!!$omp  parallel do private( IJK, I,J,K, IM,IMJK,IJMK,IJKM ) &
!!!!$omp& schedule(dynamic,chunk_size)
      DO IJK = ijkstart3, ijkend3
         IF (.NOT.WALL_AT(IJK)) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IM = IM1(I)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)
!=======================================================================
! JFD: START MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
            IF(.NOT.CUT_CELL_AT(IJK)) THEN

               TRD_G(IJK) = (X_E(I)*U_G(IJK)-X_E(IM)*U_G(IMJK))*OX(I)*ODX(I) + (&
                 V_G(IJK)-V_G(IJMK))*ODY(J) + (W_G(IJK)-W_G(IJKM))*(OX(I)*ODZ(K))

            ELSE  ! CUT CELL

               BCV = BC_ID(IJK)

               IF(BCV > 0 ) THEN
                  BCT = BC_TYPE(BCV)
               ELSE
                  BCT = 'NONE'
               ENDIF

               SELECT CASE (BCT)
                  CASE ('CG_NSW')
                     NOC_TRDG = .TRUE.
                     UW_g = ZERO
                     VW_g = ZERO
                     WW_g = ZERO
                  CASE ('CG_FSW')
                     NOC_TRDG = .FALSE.
                     UW_g = ZERO
                     VW_g = ZERO
                     WW_g = ZERO
                  CASE('CG_PSW')
                     IF(BC_HW_G(BCV)==UNDEFINED) THEN   ! same as NSW
                        NOC_TRDG = .TRUE.
                        UW_g = BC_UW_G(BCV)
                        VW_g = BC_VW_G(BCV)
                        WW_g = BC_WW_G(BCV)
                     ELSEIF(BC_HW_G(BCV)==ZERO) THEN   ! same as FSW
                        NOC_TRDG = .FALSE.
                        UW_g = ZERO
                        VW_g = ZERO
                        WW_g = ZERO
                     ELSE                              ! partial slip
                        NOC_TRDG = .FALSE.
                     ENDIF
                  CASE ('CG_MI')
                     TRD_G(IJK) = ZERO
                     RETURN
                  CASE ('CG_PO')
                     TRD_G(IJK) = ZERO
                     RETURN
                  CASE ('NONE')
                     TRD_G(IJK) = ZERO
                     RETURN
               END SELECT


               IF(FLOW_AT(IJK)) THEN
                  TRD_G(IJK) = ZERO
                  RETURN
               ENDIF

!              du/dx

               U_NODE_AT_E = ((.NOT.BLOCKED_U_CELL_AT(IJK)) .AND.(.NOT.WALL_U_AT(IJK)))
               U_NODE_AT_W = ((.NOT.BLOCKED_U_CELL_AT(IMJK)).AND.(.NOT.WALL_U_AT(IMJK)))

               IF(U_NODE_AT_E.AND.U_NODE_AT_W) THEN

                  Ui = HALF * (U_G(IJK) + U_G(IMJK))
                  Xi = HALF * (X_U(IJK) + X_U(IMJK))
                  Yi = HALF * (Y_U(IJK) + Y_U(IMJK))
                  Zi = HALF * (Z_U(IJK) + Z_U(IMJK))
                  Sx = X_U(IJK) - X_U(IMJK)
                  Sy = Y_U(IJK) - Y_U(IMJK)
                  Sz = Z_U(IJK) - Z_U(IMJK)

                  CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

                  IF(abs(Sx) > ZERO) THEN
                     dudx =  (U_G(IJK) - U_G(IMJK))/Sx
                     IF(NOC_TRDG) dudx = dudx - ((Ui-UW_g)/(Sx*DEL_H)*(Sy*Ny+Sz*Nz))
                  ELSE
                     dudx = ZERO
                  ENDIF

               ELSE
                  dudx = ZERO
               ENDIF

!              dv/dy

               V_NODE_AT_N = ((.NOT.BLOCKED_V_CELL_AT(IJK)) .AND.(.NOT.WALL_V_AT(IJK)))
               V_NODE_AT_S = ((.NOT.BLOCKED_V_CELL_AT(IJMK)).AND.(.NOT.WALL_V_AT(IJMK)))

               IF(V_NODE_AT_N.AND.V_NODE_AT_S) THEN

                  Vi = HALF * (V_G(IJK) + V_G(IJMK))
                  Xi = HALF * (X_V(IJK) + X_V(IJMK))
                  Yi = HALF * (Y_V(IJK) + Y_V(IJMK))
                  Zi = HALF * (Z_V(IJK) + Z_V(IJMK))
                  Sx = X_V(IJK) - X_V(IJMK)
                  Sy = Y_V(IJK) - Y_V(IJMK)
                  Sz = Z_V(IJK) - Z_V(IJMK)

                  CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

                  IF(abs(Sy) > ZERO) THEN
                     dvdy =  (V_G(IJK) - V_G(IJMK))/Sy
                     IF(NOC_TRDG) dvdy = dvdy - ((Vi-VW_g)/(Sy*DEL_H)*(Sx*Nx+Sz*Nz))
                  ELSE
                     dvdy =  ZERO
                  ENDIF

               ELSE IF (V_NODE_AT_N.AND.(.NOT.V_NODE_AT_S).AND.NOC_TRDG) THEN
                  CALL GET_DEL_H(IJK,'SCALAR',X_V(IJK),Y_V(IJK),Z_V(IJK),DEL_H,Nx,Ny,Nz)
                  dvdy = (V_g(IJK) - VW_g) / DEL_H * Ny

               ELSE IF ((.NOT.V_NODE_AT_N).AND.V_NODE_AT_S.AND.NOC_TRDG) THEN
                  CALL GET_DEL_H(IJK,'SCALAR',X_V(IJMK),Y_V(IJMK),Z_V(IJMK),DEL_H,Nx,Ny,Nz)
                  dvdy = (V_g(IJMK) - VW_g) / DEL_H * Ny


               ELSE
                  dvdy = ZERO
               ENDIF

!              dw/dz

               IF(NO_K) THEN

                  dwdz = ZERO

               ELSE

                  W_NODE_AT_T = ((.NOT.BLOCKED_W_CELL_AT(IJK)) .AND.(.NOT.WALL_W_AT(IJK)))
                  W_NODE_AT_B = ((.NOT.BLOCKED_W_CELL_AT(IJKM)).AND.(.NOT.WALL_W_AT(IJKM)))

                  IF(W_NODE_AT_T.AND.W_NODE_AT_B) THEN

                     Wi = HALF * (W_G(IJK) + W_G(IJKM))
                     Xi = HALF * (X_W(IJK) + X_W(IJKM))
                     Yi = HALF * (Y_W(IJK) + Y_W(IJKM))
                     Zi = HALF * (Z_W(IJK) + Z_W(IJKM))
                     Sx = X_W(IJK) - X_W(IJKM)
                     Sy = Y_W(IJK) - Y_W(IJKM)
                     Sz = Z_W(IJK) - Z_W(IJKM)

                     CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

                     IF(abs(Sz) > ZERO) THEN
                        dwdz =  (W_G(IJK) - W_G(IJKM))/Sz
                        IF(NOC_TRDG) dwdz = dwdz - ((Wi-WW_g)/(Sz*DEL_H)*(Sx*Nx+Sy*Ny))
                     ELSE
                        dwdz = ZERO
                     ENDIF

                  ELSE IF (W_NODE_AT_T.AND.(.NOT.W_NODE_AT_B).AND.NOC_TRDG) THEN
                     CALL GET_DEL_H(IJK,'SCALAR',X_W(IJK),Y_W(IJK),Z_W(IJK),DEL_H,Nx,Ny,Nz)
                     dwdz = (W_g(IJK) - WW_g) / DEL_H * Nz

                  ELSE IF ((.NOT.W_NODE_AT_T).AND.W_NODE_AT_B.AND.NOC_TRDG) THEN
                     CALL GET_DEL_H(IJK,'SCALAR',X_W(IJKM),Y_W(IJKM),Z_W(IJKM),DEL_H,Nx,Ny,Nz)
                     dwdz = (W_g(IJKM) - WW_g) / DEL_H * Nz

                  ELSE
                     dwdz = ZERO
                  ENDIF

               ENDIF


               TRD_G(IJK) = dudx + dvdy + dwdz


            ENDIF  ! CUT CELL

! Original term:
!            TRD_G(IJK) = (X_E(I)*U_G(IJK)-X_E(IM)*U_G(IMJK))*OX(I)*ODX(I) + (&
!               V_G(IJK)-V_G(IJMK))*ODY(J) + (W_G(IJK)-W_G(IJKM))*(OX(I)*ODZ(K))
!=======================================================================
! JFD: END MODIFICATION FOR CARTESIAN GRID IMPLEMENTATION
!=======================================================================
         ENDIF
      END DO

      RETURN
      END SUBROUTINE CALC_TRD_G

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_CALC_VEL_G_GRAD(IJK,DELV, IER)                      C
!  Purpose: Calculate velocity derivatives in scalar cut-cell          C
!           Gas phase                                                  C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 25-JAN-96  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CG_CALC_VEL_G_GRAD(IJK,DELV)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!     Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE geometry
      USE fldvar
      USE indices
      USE compar
      USE sendrecv
      USE bc
      USE cutcell
      USE quadric
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Indices
      INTEGER          I, J, K, IJK, IMJK, IJMK, IJKM
!
      DOUBLE PRECISION :: DEL_H,Nx,Ny,Nz
      DOUBLE PRECISION :: dudx,dudy,dudz
      DOUBLE PRECISION :: dvdx,dvdy,dvdz
      DOUBLE PRECISION :: dwdx,dwdy,dwdz
      DOUBLE PRECISION :: Xi,Yi,Zi,Ui,Vi,Wi,Sx,Sy,Sz
      DOUBLE PRECISION :: UW_g,VW_g,WW_g
      DOUBLE PRECISION, DIMENSION (3,3) :: DELV

!              |  du/dx    du/dy   du/dz  |
!      DELV =  |  dv/dx    dv/dy   dv/dz  |  =  dUi/dxj
!              |  dw/dx    dw/dy   dw/dz  |

      LOGICAL :: U_NODE_AT_E, U_NODE_AT_W
      LOGICAL :: V_NODE_AT_N, V_NODE_AT_S
      LOGICAL :: W_NODE_AT_T, W_NODE_AT_B
      INTEGER :: BCV
      CHARACTER(LEN=9) :: BCT
!-----------------------------------------------
!
!
!!!!$omp  parallel do private( IJK, I,J,K, IM,IMJK,IJMK,IJKM ) &
!!!!$omp& schedule(dynamic,chunk_size)

      DELV = ZERO

      IF (.NOT.CUT_CELL_AT(IJK)) RETURN

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)

      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)

      BCV = BC_ID(IJK)

      IF(BCV > 0 ) THEN
         BCT = BC_TYPE(BCV)
      ELSE
         BCT = 'NONE'
      ENDIF

      SELECT CASE (BCT)
         CASE ('CG_NSW')
            NOC_TRDG = .TRUE.
            UW_g = ZERO
            VW_g = ZERO
            WW_g = ZERO
         CASE ('CG_FSW')
            NOC_TRDG = .FALSE.
            UW_g = ZERO
            VW_g = ZERO
            WW_g = ZERO
         CASE('CG_PSW')
            IF(BC_HW_G(BCV)==UNDEFINED) THEN   ! same as NSW
               NOC_TRDG = .TRUE.
               UW_g = BC_UW_G(BCV)
               VW_g = BC_VW_G(BCV)
               WW_g = BC_WW_G(BCV)
            ELSEIF(BC_HW_G(BCV)==ZERO) THEN   ! same as FSW
               NOC_TRDG = .FALSE.
               UW_g = ZERO
               VW_g = ZERO
               WW_g = ZERO
            ELSE                              ! partial slip
               NOC_TRDG = .FALSE.
            ENDIF
         CASE ('CG_MI')
            DELV = ZERO
            RETURN
         CASE ('CG_PO')
            DELV = ZERO
            RETURN
         CASE ('NONE')
            DELV = ZERO
            RETURN
      END SELECT

      IF(FLOW_AT(IJK)) THEN
         DELV = ZERO
         RETURN
      ENDIF

!=======================================================================
!              du/dx, du/dy, du/dz
!=======================================================================

      U_NODE_AT_E = ((.NOT.BLOCKED_U_CELL_AT(IJK)) .AND.(.NOT.WALL_U_AT(IJK)))
      U_NODE_AT_W = ((.NOT.BLOCKED_U_CELL_AT(IMJK)).AND.(.NOT.WALL_U_AT(IMJK)))

      IF(U_NODE_AT_E.AND.U_NODE_AT_W) THEN

         Ui = HALF * (U_G(IJK) + U_G(IMJK))
         Xi = HALF * (X_U(IJK) + X_U(IMJK))
         Yi = HALF * (Y_U(IJK) + Y_U(IMJK))
         Zi = HALF * (Z_U(IJK) + Z_U(IMJK))
         Sx = X_U(IJK) - X_U(IMJK)
         Sy = Y_U(IJK) - Y_U(IMJK)
         Sz = Z_U(IJK) - Z_U(IMJK)

         CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

         IF(abs(Sx) > ZERO) THEN
            dudx =  (U_G(IJK) - U_G(IMJK))/Sx
            dudy =  ZERO
            dudz =  ZERO
            IF(NOC_TRDG) THEN
               dudx = dudx - ((Ui-UW_g)/(Sx*DEL_H)*(Sy*Ny+Sz*Nz))
               dudy = (Ui-UW_g) / DEL_H * Ny
               dudz = (Ui-UW_g) / DEL_H * Nz
            ENDIF
         ELSE
            dudx = ZERO
            dudy = ZERO
            dudz = ZERO
         ENDIF


      ELSE IF (U_NODE_AT_E.AND.(.NOT.U_NODE_AT_W).AND.NOC_TRDG) THEN

         CALL GET_DEL_H(IJK,'SCALAR',X_U(IJK),Y_U(IJK),Z_U(IJK),DEL_H,Nx,Ny,Nz)

         dudx = (U_g(IJK) - UW_g) / DEL_H * Nx
         dudy = (U_g(IJK) - UW_g) / DEL_H * Ny
         dudz = (U_g(IJK) - UW_g) / DEL_H * Nz
      ELSE IF ((.NOT.U_NODE_AT_E).AND.U_NODE_AT_W.AND.NOC_TRDG) THEN
         CALL GET_DEL_H(IJK,'SCALAR',X_U(IMJK),Y_U(IMJK),Z_U(IMJK),DEL_H,Nx,Ny,Nz)
         dudx = (U_g(IMJK) - UW_g) / DEL_H * Nx
         dudy = (U_g(IMJK) - UW_g) / DEL_H * Ny
         dudz = (U_g(IMJK) - UW_g) / DEL_H * Nz
      ELSE
         dudx = ZERO
         dudy = ZERO
         dudz = ZERO
      ENDIF


      DELV(1,1) = dudx
      DELV(1,2) = dudy
      DELV(1,3) = dudz
!=======================================================================
!              dv/dx, dv/dy, dv/dz
!=======================================================================

      V_NODE_AT_N = ((.NOT.BLOCKED_V_CELL_AT(IJK)) .AND.(.NOT.WALL_V_AT(IJK)))
      V_NODE_AT_S = ((.NOT.BLOCKED_V_CELL_AT(IJMK)).AND.(.NOT.WALL_V_AT(IJMK)))

      IF(V_NODE_AT_N.AND.V_NODE_AT_S) THEN

         Vi = HALF * (V_G(IJK) + V_G(IJMK))
         Xi = HALF * (X_V(IJK) + X_V(IJMK))
         Yi = HALF * (Y_V(IJK) + Y_V(IJMK))
         Zi = HALF * (Z_V(IJK) + Z_V(IJMK))
         Sx = X_V(IJK) - X_V(IJMK)
         Sy = Y_V(IJK) - Y_V(IJMK)
         Sz = Z_V(IJK) - Z_V(IJMK)

         CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)




         IF(abs(Sy) > ZERO) THEN
            dvdx =  ZERO
            dvdy =  (V_G(IJK) - V_G(IJMK))/Sy
            dvdz =  ZERO
            IF(NOC_TRDG) THEN
               dvdx = (Vi-VW_g) / DEL_H * Nx
               dvdy = dvdy - ((Vi-VW_g)/(Sy*DEL_H)*(Sx*Nx+Sz*Nz))
               dvdz = (Vi-VW_g) / DEL_H * Nz
            ENDIF
         ELSE
            dvdx =  ZERO
            dvdy =  ZERO
            dvdz =  ZERO
         ENDIF


      ELSE IF (V_NODE_AT_N.AND.(.NOT.V_NODE_AT_S).AND.NOC_TRDG) THEN
         CALL GET_DEL_H(IJK,'SCALAR',X_V(IJK),Y_V(IJK),Z_V(IJK),DEL_H,Nx,Ny,Nz)
         dvdx = (V_g(IJK) - VW_g) / DEL_H * Nx
         dvdy = (V_g(IJK) - VW_g) / DEL_H * Ny
         dvdz = (V_g(IJK) - VW_g) / DEL_H * Nz
      ELSE IF ((.NOT.V_NODE_AT_N).AND.V_NODE_AT_S.AND.NOC_TRDG) THEN
         CALL GET_DEL_H(IJK,'SCALAR',X_V(IJMK),Y_V(IJMK),Z_V(IJMK),DEL_H,Nx,Ny,Nz)
         dvdx = (V_g(IJMK) - VW_g) / DEL_H * Nx
         dvdy = (V_g(IJMK) - VW_g) / DEL_H * Ny
         dvdz = (V_g(IJMK) - VW_g) / DEL_H * Nz
      ELSE
         dvdx =  ZERO
         dvdy =  ZERO
         dvdz =  ZERO
      ENDIF


      DELV(2,1) = dvdx
      DELV(2,2) = dvdy
      DELV(2,3) = dvdz

!=======================================================================
!              dw/dx, dw/dy, dw/dz
!=======================================================================

      IF(NO_K) THEN

         dwdx = ZERO
         dwdy = ZERO
         dwdz = ZERO

      ELSE

         W_NODE_AT_T = ((.NOT.BLOCKED_W_CELL_AT(IJK)) .AND.(.NOT.WALL_W_AT(IJK)))
         W_NODE_AT_B = ((.NOT.BLOCKED_W_CELL_AT(IJKM)).AND.(.NOT.WALL_W_AT(IJKM)))

         IF(W_NODE_AT_T.AND.W_NODE_AT_B) THEN

            Wi = HALF * (W_G(IJK) + W_G(IJKM))
            Xi = HALF * (X_W(IJK) + X_W(IJKM))
            Yi = HALF * (Y_W(IJK) + Y_W(IJKM))
            Zi = HALF * (Z_W(IJK) + Z_W(IJKM))
            Sx = X_W(IJK) - X_W(IJKM)
            Sy = Y_W(IJK) - Y_W(IJKM)
            Sz = Z_W(IJK) - Z_W(IJKM)

            CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)

            IF(abs(Sz) > ZERO) THEN
               dwdx =  ZERO
               dwdy =  ZERO
               dwdz =  (W_G(IJK) - W_G(IJKM))/Sz
               IF(NOC_TRDG) THEN
                  dwdx = (Wi-WW_g) / DEL_H * Nx
                  dwdy = (Wi-WW_g) / DEL_H * Ny
                  dwdz = dwdz - ((Wi-WW_g)/(Sz*DEL_H)*(Sx*Nx+Sy*Ny))
               ENDIF
            ELSE
               dwdx = ZERO
               dwdy = ZERO
               dwdz = ZERO
            ENDIF

         ELSE IF (W_NODE_AT_T.AND.(.NOT.W_NODE_AT_B).AND.NOC_TRDG) THEN
            CALL GET_DEL_H(IJK,'SCALAR',X_W(IJK),Y_W(IJK),Z_W(IJK),DEL_H,Nx,Ny,Nz)
            dwdx = (W_g(IJK) - WW_g) / DEL_H * Nx
            dwdy = (W_g(IJK) - WW_g) / DEL_H * Ny
            dwdz = (W_g(IJK) - WW_g) / DEL_H * Nz

         ELSE IF ((.NOT.W_NODE_AT_T).AND.W_NODE_AT_B.AND.NOC_TRDG) THEN
            CALL GET_DEL_H(IJK,'SCALAR',X_W(IJKM),Y_W(IJKM),Z_W(IJKM),DEL_H,Nx,Ny,Nz)
            dwdx = (W_g(IJKM) - WW_g) / DEL_H * Nx
            dwdy = (W_g(IJKM) - WW_g) / DEL_H * Ny
            dwdz = (W_g(IJKM) - WW_g) / DEL_H * Nz

         ELSE
            dwdx = ZERO
            dwdy = ZERO
            dwdz = ZERO
         ENDIF

      ENDIF

      DELV(3,1) = dwdx
      DELV(3,2) = dwdy
      DELV(3,3) = dwdz

      RETURN
      END SUBROUTINE CG_CALC_VEL_G_GRAD

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3



