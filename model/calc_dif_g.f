!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_DIF_g(IER)
!  Purpose: Calculate the effective diffusivity of fluid phase        C
!                                                                      C
!  Author:M. Syamlal                                  Date: 13-FEB-98  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!
!  Updated with the dilute mixture approximation for calculation of    C
!  multicomponent diffusion coefficients                               C
!  Author:N. Reuge                                    Date: 11-APR-07  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_DIF_G()
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
      USE scales
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

!                      Indices
      INTEGER          IJK, N

      DOUBLE PRECISION Dab(3,3), Tg0, Pg0

!-----------------------------------------------

! Default gas diffusion coefficient
! Bird, Stewart, and Lightfoot (1960) -- CO2--N2 at 298.2 K
      Dab(1,2) = 0.165D0       !cm^2/s
      Tg0 = 298.2D0
      Pg0 = 1.01D6             !dyne

      IF(UNITS == 'SI') THEN
         Dab(1,2) = Dab(1,2)*0.0001D0   !m^2/s
         Pg0 = Pg0/10.D0                !Pa
      ENDIF


!********Gas diffusion coef for a system of 3 species***************
! Species: SiH4, H2 & N2
! Calculated using relation derived from Chapman and Enskog's theory
! of gases - Reid, Prausnitz and Poling (1987)

! Binary diffusion coefficient SiH4/H2 at 873 K
!      Dab(1,2) = 3.78      ! cm^2/s
!      Dab(2,1) = Dab(1,2)

! Binary diffusion coefficient SiH4/N2 at 873 K
!      Dab(1,3) = 1.02      ! cm^2/s
!      Dab(3,1) = Dab(1,3)

! Binary diffusion coefficient H2/N2 at 873 K
!      Dab(2,3) = 4.52      ! cm^2/s
!      Dab(3,2) = Dab(2,3)

!      Tg0 = 873.0
!      Pg0 = 1.01e6         ! dyne

!      IF(UNITS == 'SI') THEN
!         DO N1 = 1, NMAX(0)-1
!            DO N2 = N1+1, NMAX(0)
!               Dab(N1,N2) = Dab(N1,N2)*0.0001D0   !m^2/s
!               Dab(N2,N1) = Dab(N1,N2)
!            ENDDO
!         ENDDO
!         Pg0 = Pg0/10.D0                    !Pa
!      ENDIF

!*******************************************************************

      IF (DIF_G0 /= UNDEFINED) RETURN

!!!!$omp  parallel do private(ijk) &
!!!!$omp& schedule(dynamic,chunk_size)

! Default calculation of diffusivities
! Influence of gas temperature and gas pressure from Fuller relation
      DO N = 1, NMAX(0)
         DO IJK = IJKSTART3, IJKEND3
            IF (FLUID_AT(IJK)) THEN
               DIF_G(IJK,N) = ROP_G(IJK)*Dab(1,2)*(T_g(IJK)/Tg0)**1.75 * &
                              Pg0/(P_g(IJK)+P_REF)
            ELSE
               DIF_G(IJK,N) = ZERO
            ENDIF
         ENDDO
      ENDDO


!*******************************************************************
! Calculation of diffusivities using the dilute mixture approximation for
! multicomponent diffusion - Curtiss-Hirschfelder, Wilke & Blanc
! Valid if the mass fraction of the carrier species > 0.9

! Influence of gas temperature and gas pressure from Fuller relation
!      DO N = 1, NMAX(0)
!         DO IJK = IJKSTART3, IJKEND3
!            IF (FLUID_AT(IJK)) THEN
!               IF ((1.0-X_g(IJK,N)) > 1.e-8) THEN
!                  SUMJ = ZERO
!                  DO N2 = 1, NMAX(0)
!                     IF (N2 /= N) SUMJ = SUMJ+X_g(IJK,N2)/Dab(N,N2)
!                  ENDDO
!                  DIF_G(IJK,N) = ROP_G(IJK)*(1-X_g(IJK,N))/SUMJ * &
!                                 (T_g(IJK)/Tg0)**1.75*Pg0/P_g(IJK)
!               ELSE
!                  DIF_G(IJK,N) = ROP_G(IJK)*Dab(1,2)
!               ENDIF
!            ELSE
!               DIF_G(IJK,N) = ZERO
!            ENDIF
!         ENDDO
!      ENDDO
!****************************************************************************


      RETURN
      END SUBROUTINE CALC_DIF_G

