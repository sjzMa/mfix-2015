!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_K_g(IER)                                          C
!  Purpose: Calculate the effective conductivity of fluid phase        C
!                                                                      C
!  Author:M. Syamlal                                  Date: 24-APR-96  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: allow SI unit                                              C
!  Author: S. Dartevelle                              Date: 01-Jul-02  C
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

      SUBROUTINE CALC_K_G()
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
      USE compar
      USE run
      USE sendrecv
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      Indices
      INTEGER          IJK
!-----------------------------------------------

      IF (K_G0 /= UNDEFINED) RETURN

!!!!$omp parallel do private(ijk) &
!!!!$omp& schedule(dynamic,chunk_size)

      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN
! Gas conductivity (air)
! Bird, Stewart, and Lightfoot (1960) -- Temperature dependence from formula
! 8.3-12 on p. 255 and conductivity value at 300 K from p. 263
            K_G(IJK) = 6.02D-5*SQRT(T_G(IJK)/300.D0)           ! cal/(s.cm.K)
         ELSE
            K_G(IJK) = ZERO
         ENDIF
! 1 cal = 4.183925D0 J
         IF (UNITS == 'SI') K_G(IJK) = 418.3925D0*K_G(IJK)      !J/s.m.K

      ENDDO

      CALL send_recv(K_G, 2)

      RETURN
      END SUBROUTINE CALC_K_G

