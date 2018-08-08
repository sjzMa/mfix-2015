!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_H (IJK, M, N)                                       C
!  Purpose: Calculate specific enthalpy of species N in phase M        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 27-DEC-2007C
!  Reviewer:                                         Date:   C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:None                                           C
!  Variables modified:None                                             C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      DOUBLE PRECISION FUNCTION CALC_H(refT, M, N)
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE physprop
      USE fldvar

      USE constant, only: RGAS => GAS_CONST_cal
      USE read_thermochemical, only: calc_ICpoR

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!                      cell, phase and species indices

      DOUBLE PRECISION, INTENT(IN) :: refT   ! Temperature

      INTEGER, INTENT(IN) :: M ! Phase index
      INTEGER, INTENT(IN) :: N ! Species index

      DOUBLE PRECISION ICpoR
      DOUBLE PRECISION lMW

      INTEGER :: IER

!-----------------------------------------------
!
      IER = 0

      if(M == 0)then
         lMW = MW_g(N)
      else
         lMW = MW_s(M,N)
      endif

! Integrate the specific heat from zero to refT
      ICpoR = calc_ICpoR(refT, M, N, IER)

! Evaluate the enthalpy of speices N at refT
      CALC_H = (HfrefoR(M,N)  + ICpoR) * (RGAS / lMW)

      RETURN
      END FUNCTION CALC_H
