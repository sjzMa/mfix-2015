!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_CONDUCTION                                         !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 25-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!  REF: Batchelor and O'Brien, "Thermal or electrical conduction       !
!       through a granular material," Proceedings of the Royal Society !
!       of London. Series A, Mathematical and Physical Sciences,       !
!       Vol. 355, no 1682, pp. 313-333, July 1977.                     !
!                                                                      !
!  REF: Rong and Horio, "DEM simulation of char combustion in a        !
!       fluidized bed," in Second International Conference on CFD in   !
!       the materials and Process Industries, Melbourne, 1999, pp.     !
!       65-70.                                                         !
!                                                                      !
!  REF: Xavier and Davidson, "Heat transfer to surfaces immersed in    !
!       fluidised beds, particularly tub arrays," in Fluidization:     !
!       Proceedings of the Second Engineering Foundation Conference,   !
!       Cambridge, England, 2-6 April, 1978, pp. 333-338.              !
!                                                                      !
!  REF: Zhou, Yu, and Zulli, "Particle scale study of heat transfer in !
!       packed and bubbling fluidized beds," AIChE Journal, Vol. 55,   !
!       no 4, pp 868-884, 2009.                                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
MODULE DES_THERMO_COND

CONTAINS

      FUNCTION DES_CONDUCTION(I, J, CENTER_DIST, iM, iIJK)

      use constant
      use des_thermo
      use discretelement
      use funits
      use physprop

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
! Index of particle being looped over
      INTEGER, INTENT(IN) :: I,J
! Index of fluid cell containing particle I
      INTEGER, INTENT(IN) :: iIJK
! Solid phase indices for the given particles
      INTEGER, INTENT(IN) :: iM
! Distance between the centers of particle I and particle J (magnitude)
      DOUBLE PRECISION, INTENT(IN) :: CENTER_DIST
      DOUBLE PRECISION :: DES_CONDUCTION

! Local variables
!---------------------------------------------------------------------//
! Solid phase indices for the given particles
      INTEGER jM
! Radius of smaller particle
      DOUBLE PRECISION MIN_RAD
! Radius of larger particle
      DOUBLE PRECISION MAX_RAD
! Rate of particle-particle conduction
      DOUBLE PRECISION Q_pp
! Rate of particle-fluid-particle conduction
      DOUBLE PRECISION Q_pfp
! Outer radius of region delineating particle-fluid-particle conduction
      DOUBLE PRECISION RD_OUT
! Inner radius of region delineating particle-fluid-particle conduction
      DOUBLE PRECISION RD_IN
! The radius of the fluid lens containing the larger particle
      DOUBLE PRECISION LENS_RAD
! Temperature difference between two particles
      DOUBLE PRECISION DeltaTp
! Effective thermal conductivity
      DOUBLE PRECISION lK_eff
! Radius of contact area
      DOUBLE PRECISION lRadius

! Functions
!---------------------------------------------------------------------//

! Identify the solid phases of the neighbor particle
         jM = PIJK(J,5) + SMAX

! Determine the radius of the larger and smaller particle
         MIN_RAD = MIN(DES_RADIUS(I), DES_RADIUS(J))
         MAX_RAD = MAX(DES_RADIUS(I), DES_RADIUS(J))

         DeltaTp = DES_T_s_NEW(J) - DES_T_s_NEW(I)

! Calculate the particle-particle conduction
! REF: Batchelor and O'Brien, 1977 (MODIFIED)
!---------------------------------------------------------------------//
         Q_pp = 0.0
         IF(CENTER_DIST < (MAX_RAD + MIN_RAD)) THEN
! Effective thermal conductivity
            lK_eff = K_eff(K_s0(iM),K_s0(jM))
! Effective contact area's radius
            lRadius = RADIUS(MAX_RAD, MIN_RAD)
! Inter-particle heat transfer
            Q_pp = 2.0d0 * lK_eff * lRadius * DeltaTp
! Assign the inter-particle heat transfer to both particles.
         ENDIF

! Calculate the particle-fluid-particle conduction
! REF: Rong and Horio, 1999 (MODIFIED)
!---------------------------------------------------------------------//
! Calculate the radius of the fluid lens surrounding the larger particle
! Default FLPC = 0.2
         LENS_RAD = MAX_RAD * (1.0D0 + FLPC)

! Calculate the outer radial distance of the region for particle-fluid-
! particle heat conduction.
         RD_OUT = RADIUS( LENS_RAD, MIN_RAD)

! If the value returned is less than zero, then the fluid lens
! surrounding the larger particle does not intersect with the surface
! of the smaller particle. In this case, particle-fluild-particle
! conduction does not occur.
         Q_pfp = 0.0
         IF(RD_OUT .GT. ZERO)THEN
! Calculate the distance from the line connecting the particles' centers
! to the point of contact between the two particles. This value is
! zero if the particles are not touching and is the radius of the
! shared contact area otherwise.
            RD_IN = ZERO
            IF(CENTER_DIST < (MAX_RAD + MIN_RAD) ) &
               RD_IN = RADIUS(MAX_RAD, MIN_RAD)
! Calculate the rate of heat transfer between the particles through the
! fluid using adaptive Simpson's rule to manage the integral.
            Q_pfp = K_g(iIJK) * DeltaTp * &
               ADPT_SIMPSON(RD_IN,RD_OUT)

         ENDIF ! RD_OUT > ZERO

         ! Return total heat flux
         DES_CONDUCTION = Q_pp + Q_pfp

      RETURN

      CONTAINS

!......................................................................!
! Subroutine Procedure: FUNCTION RADIUS                                !
!                                                                      !
! Propose: Calculate the center line radius between two particles.     !
! Depending on what is given as input, this function calculates:       !
!   1) radius of contact area between two particles                    !
!   2) radius delineating particle-fluid-particle region.              !
!......................................................................!
      DOUBLE PRECISION FUNCTION RADIUS(R1, R2)

      IMPLICIT NONE

! Radius values
      DOUBLE PRECISION, INTENT(IN) :: R1, R2
! Distance between particle centers
      DOUBLE PRECISION BETA
! Generic value
      DOUBLE PRECISION VALUE

! Calculate
      VALUE = (R1**2 - R2**2 + CENTER_DIST**2)/(2.0d0 * R1 * CENTER_DIST)
! Check to ensure that VALUE is less than or equal to one. If VALUE is
! greater than one, the triangle inequality has been violated. Therefore
! there is no intersection between the fluid lens surrounding the larger
! particle and the surface of the smaller particle.
! Thus, there is no particle-fluid-particle heat transfer.
      IF( VALUE .GT. 1.0d0) THEN
         RADIUS = -1.0d0
      ELSE
! Calculate beta (Law of cosines)
         BETA = ACOS( VALUE )
! Calculate the radius
         RADIUS = R1 * SIN(BETA)
      ENDIF

      RETURN
      END FUNCTION RADIUS

!......................................................................!
! Subroutine Procedure: FUNCTION K_eff                                 !
!                                                                      !
! Propose: Calculate the effective thermal conductivity of two         !
! particles with conductivities K1 and K2                              !
!......................................................................!
      DOUBLE PRECISION FUNCTION K_eff(K1, K2)

      IMPLICIT NONE

! Thermal conductivities
      DOUBLE PRECISION, INTENT(IN) :: K1, K2

      K_eff = 2.0d0 * (K1*K2)/(K1 + K2)

      RETURN
      END FUNCTION K_eff

!``````````````````````````````````````````````````````````````````````!
! Function: F                                                          !
!                                                                      !
! Purpose: This function defines the region of particle-fluid-particle !
!          heat transfer. Note that this function develops a           !
!          singularity at the boundary of the contact region when two  !
!          particles are touching.  This singularity is resolved by    !
!          making the assumption that the surfaces of two particles    !
!          never directly touch (c.f. Rong and Horio, 1999). By        !
!          default, it is assumed that the particles are separated by  !
!          a minimum length of 4.0x10^(-10) meters. This value may be  !
!          modified by specifying a different value in the mfix.dat    !
!          file under the variable name DES_MIN_COND_DIST.             !
!                                                                      !
!          Note, the closer this value is to zero, the harder it will  !
!          be for the quadrature routine to converge.                  !
!......................................................................!
      DOUBLE PRECISION FUNCTION F(R)
        USE utilities
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: R

      F = (2.0d0*Pi*R)/MAX((CENTER_DIST - SQRT(MAX_RAD**2-R**2) - &
         SQRT(MIN_RAD**2-R**2)), DES_MIN_COND_DIST)

      IF( mfix_isnan(F))PRINT*,'F IS NAN'

      END FUNCTION F

!``````````````````````````````````````````````````````````````````````!
! Function: ADPT_SIMPSON                                               !
!                                                                      !
! Purpose: This function applies adaptive Simpson's Rule to solve the  !
!          definite integral of the function F.                        !
!                                                                      !
!         *NOTE: This route will return the result of the method even  !
!          if convergence has not been achieved. If this occurs, a     !
!          message will be written to the log file to notify the user. !
!......................................................................!
      DOUBLE PRECISION FUNCTION ADPT_SIMPSON(A, B)

      IMPLICIT NONE

! Bounds of integration
      DOUBLE PRECISION, INTENT(IN) :: A, B
! Maximum recursive depth
      INTEGER, PARAMETER :: Kmax = 25
! Current depth of recursion
      INTEGER :: K
! Array storage for recursion
      DOUBLE PRECISION :: V(Kmax,6)
! Step size
      DOUBLE PRECISION :: H
! Simpson's Rule evaluated on the left and right intervals
      DOUBLE PRECISION :: lS, rS
! Value of the function F evaluate at the midpoints of the left and
! right intervals
      DOUBLE PRECISION :: rF, lF
! Local error bound
      DOUBLE PRECISION :: Err_BND
! Convergence bound
      DOUBLE PRECISION :: EPS = 10.0**(-2)

! Error indicating that an error has been specified to the user.
      LOGICAL, SAVE :: ADPT_SIMPSON_ERR = .FALSE.

! Initialize variables
      V(:,:) = 0.0d0   !
      ADPT_SIMPSON = 0.0d0   ! Integral value
      H = (B-A)/2.0d0 ! Dynamic interval length
! Calculate Simpson's Rule over the interval [A,B]
      lS = (F(A) + 4.0d0*F((A+B)/2.0d0) + F(B))*(H/3.0d0)
! Initialize the storage vector for the first pass
      V(1,:) = (/ A, H, F(A), F((A+B)/2.0d0), F(B), lS/)
      K = 1
! Enter recursion calculations
      DO
! Establish the new interval length
         H = V(K,2)/2.0d0
! Evaluate the function on the left interval
         lF = F(V(K,1) + H)
! Calculate Simpson's Rule over the left interval
         lS = (V(K,3) + 4.0d0*lF + V(K,4))*(H/3.0d0)
! Evaluate the function on the right interval
         rF = F(V(K,1) + 3.0d0*H)
! Calculate Simpson's Rule over the right interval
         rS = (V(K,4) + 4.0d0*rF + V(K,5))*(H/3.0d0)
! Evaluate the error bound
         Err_BND = (30.0d0*EPS*H)/(B-A)
! Check to see if conversion on the interval has been met
         IF( ABS(lS + rS - V(K,6)) .LT. Err_BND)THEN
! Update the integral value
            ADPT_SIMPSON = ADPT_SIMPSON + lS + rS + &
               (1.0d0/15.0d0)*(lS + rS - V(K,6))
! Decrement the recursion counter
            K = K - 1
! If K=0, then integration has been successfully completed over the
! entire interval [A,B].
            IF( K == 0) RETURN
         ELSEIF( (K .GE. Kmax) .OR. &
            (H == (B-A)*(1.0d0/2.0d0)**(Kmax+3))) THEN
! Flag that the method did not converge.
            IF(.NOT.ADPT_SIMPSON_ERR)THEN
               WRITE(*,1000)
               WRITE(UNIT_LOG,1000)
               ADPT_SIMPSON_ERR = .TRUE.
            ENDIF
! Update the integral value
            ADPT_SIMPSON = ADPT_SIMPSON + lS + rS + &
               (1.0d0/15.0d0)*(lS + rS - V(K,6))
! Decrement the recursion counter
            K = K - 1
         ELSE
! Refine the subintervals through recursive splitting of the intervals
            V(K+1,:) = (/V(K,1) + 2.0d0*H, H, V(K,4), rF, V(K,5), rS/)
            V(K,:) = (/V(K,1), H, V(K,3), lF, V(K,4), lS/)
! Increment level counter
            K = K+1
         ENDIF
      ENDDO

 1000 FORMAT(/1X,70('*'),/' From: DES_COND_EQ',/, ' Message: ',        &
         'Integration of the particle-fluid-particle equation did ',   &
         'not',/' converge! No definite bound can be placed on the ',  &
         'error.',/' Future convergence messages will be suppressed!', &
         /1X,70('*'))

      END FUNCTION ADPT_SIMPSON

    END FUNCTION DES_CONDUCTION

  END MODULE DES_THERMO_COND
