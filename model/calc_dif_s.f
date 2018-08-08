!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_DIF_s(M, IER)
!  Purpose: Calculate the effective diffusivity of solids phases       C
!                                                                      C
!  Author:M. Syamlal                                  Date: 13-EFB-98  C
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

      SUBROUTINE CALC_DIF_S(M)
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

!                      Indices
      INTEGER          IJK, N

!                      Solids phase
      INTEGER          M

      DOUBLE PRECISION Dab

!-----------------------------------------------

      Dab = ZERO       !cm^2/s
      IF(UNITS == 'SI') Dab = Dab*0.0001D0   !m^2/s

      IF (DIF_s0  /= UNDEFINED) RETURN

!!!!$omp  parallel do private(n,ijk) &
!!!!$omp& schedule(dynamic,chunk_size)

      DO N = 1, NMAX(M)
         DO IJK = IJKSTART3, IJKEND3
            IF (FLUID_AT(IJK)) THEN
               DIF_S(IJK,M,N) = ROP_S(IJK,M)*Dab
            ELSE
               DIF_S(IJK,M,N) = ZERO
            ENDIF
         ENDDO
      ENDDO

      CALL SEND_RECV(DIF_S, 2)

      RETURN
      END SUBROUTINE CALC_DIF_S


