!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Init_Ab_m(A_m, b_m, IJKMAX2, M, IER)                   C                     C
!  Purpose:Initialiize the sparse matrix coefficients and the          C
!           source vector.                                             C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 16-MAY-96  C
!  Reviewer:                                          Date:            C
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
      SUBROUTINE INIT_AB_M(A_M, B_M, IJKMAX2A, M)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE matrix
      USE parallel
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Phase index
      INTEGER          M
!
!                      cell index
      INTEGER          IJK
!
!                      Maximum dimension
      INTEGER          IJKMAX2A
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Source vector
      DOUBLE PRECISION b_m(DIMENSION_3, 0:DIMENSION_M)
!
!-----------------------------------------------
!
!      IJK = 1
      IF (IJKMAX2A > 0) THEN
       IF (USE_DOLOOP) THEN
!!!$omp    parallel do private( IJK )
         DO IJK = ijkstart3, ijkend3
           A_M(IJK,B,M) = ZERO
           A_M(IJK,S,M) = ZERO
           A_M(IJK,W,M) = ZERO
           A_M(IJK,0,M) = -ONE
           A_M(IJK,E,M) = ZERO
           A_M(IJK,N,M) = ZERO
           A_M(IJK,T,M) = ZERO

           B_M(IJK,M) = ZERO
         ENDDO
       ELSE
         A_M(:,B,M) = ZERO
         A_M(:,S,M) = ZERO
         A_M(:,W,M) = ZERO
         A_M(:,0,M) = -ONE
         A_M(:,E,M) = ZERO
         A_M(:,N,M) = ZERO
         A_M(:,T,M) = ZERO
         B_M(:,M) = ZERO

!         IJK = IJKMAX2A + 1
       ENDIF
      ENDIF
      RETURN
      END SUBROUTINE INIT_AB_M

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 120 Replaced the index for initialization, :IJKMAX2A => :
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!
