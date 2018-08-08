!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_CELL_MOVEMENT                                     !
!                                                                      !
!  Purpose: Check to see if particles have moved into ghost cells.     !
!                                                                      !
!  Note: This routine needs a global communicator to identify errors.  !
!  The collection could get expensive so the call frequency of this    !
!  routine should probably be reduced.                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_CELL_MOVEMENT

! Global Variables:
!---------------------------------------------------------------------//
! Max number of particles in process
      use discretelement, only: MAX_PIP
! The I/J/K/IJK indicies of the fluid cell
      use discretelement, only: PIJK
! Run time flag indicating DEM or PIC solids.
      use run, only: DEM_SOLIDS, PIC_SOLIDS

      use functions, only: IS_NORMAL
      use mpi_utility
      use error_manager

      IMPLICIT NONE

! Local Variables:
!----------------------------------------------------------------------!
! Loop indicies:
      INTEGER :: L, I, J, K
! Integer error flag.
      INTEGER :: IER

! Initialize local variables.
      IER = 0

! Set an error flag if any errors are found. Preform a global collection
! to sync error flags. If needed, reort errors.
!.......................................................................
!!$omp parallel default(shared) private(L, I, J, K, IJK)
!!$omp do reduction(+:IER) schedule (guided,50)
      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(.NOT.IS_NORMAL(L)) CYCLE

! assigning local aliases for particle i, j, k fluid grid indices
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)

         IF(I > IEND1 .OR. I < ISTART1) IER = 1
         IF(J > JEND1 .OR. J < JSTART1) IER = 1
         IF(DO_K .AND. (K > KEND1 .OR. K < KSTART1)) IER = 1
      ENDDO
!!$omp end parallel

      CALL GLOBAL_ALL_SUM(IER)
      IF(IER == 0) RETURN

      IF(DEM_SOLIDS) CALL CHECK_CELL_MOVEMENT_DEM
      IF(PIC_SOLIDS) CALL CHECK_CELL_MOVEMENT_PIC

      RETURN
      END SUBROUTINE CHECK_CELL_MOVEMENT

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_CELL_MOVEMENT_DEM                                 !
!                                                                      !
!  Purpose: Report which DEM particles have moved into ghost cells.    !
!  This is a dead-end routine. Once called, the simulation will exit.  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_CELL_MOVEMENT_DEM

! Global Variables:
!---------------------------------------------------------------------//
! The global ID of a particle.
      use discretelement, only: iGlobal_ID
! Max number of particles in process
      use discretelement, only: MAX_PIP
! The I/J/K/IJK indicies of the fluid cell
      use discretelement, only: PIJK
! Particle positions and velocities
      use discretelement, only: DES_POS_NEW, DES_VEL_NEW

      use functions, only: IS_NORMAL
      use mpi_utility
      USE error_manager

      IMPLICIT NONE

! Local Variables:
!----------------------------------------------------------------------!
! Loop indicies:.
      INTEGER :: L, I, J, K
! Integer error flag
      INTEGER :: IER


      CALL INIT_ERR_MSG("CHECK_CELL_MOVEMENT_DEM")
      CALL OPEN_PE_LOG(IER)

      WRITE(ERR_MSG, 1100)
      CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

 1100 FORMAT('Error 1100: Particles detected in a ghost cell:',/' ')

      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(.NOT.IS_NORMAL(L)) CYCLE

! assigning local aliases for particle i, j, k fluid grid indices
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)

         IF (I.GT.IEND1 .OR. I.LT.ISTART1) THEN
            WRITE(ERR_MSG, 1101) trim(iVal(iGlobal_ID(L))),'I',        &
               trim(iVal(I)),'X',DES_POS_NEW(1,L),'X',DES_VEL_NEW(1,L)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF

         IF(J.GT.JEND1 .OR. J.LT.JSTART1) THEN
            WRITE(ERR_MSG, 1101) trim(iVal(iGlobal_id(L))),'J',        &
               trim(iVal(J)),'Y',DES_POS_NEW(2,L),'Y',DES_VEL_NEW(2,L)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF

         IF (DO_K .AND. (K.GT.KEND1 .OR. K.LT.KSTART1)) THEN
            WRITE(ERR_MSG, 1101) trim(iVal(iGlobal_ID(L))),'K',        &
               trim(iVal(K)),'Z',DES_POS_NEW(3,L),'Z',DES_VEL_NEW(3,L)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF
      ENDDO

 1101 FORMAT('Particle ',A,' moved into cell with ',A,' index ',A,/ &
         3x,A,'-Position: ',g11.4,6x,A,'-Velocity:',g11.4,/' ')

      WRITE(ERR_MSG, 1102)
      CALL FLUSH_ERR_MSG(HEADER=.FALSE.)
 1102 FORMAT('This is a fatal error. A particle output file (vtp) ',   &
         'will be written',/'to aid debugging.')


      CALL WRITE_DES_DATA
      CALL MFIX_EXIT(myPE)

      END SUBROUTINE CHECK_CELL_MOVEMENT_DEM


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_CELL_MOVEMENT_PIC                                 !
!                                                                      !
!  Purpose: Report which PIC particles have moved into ghost cells.    !
!  Unlike DEM particles, this routine either deletes particles that    !
!  are out of the domain, or reassigns the index to try and recover    !
!  the particle.                                                       !
!                                                                      !
!  Notes:                                                              !
!                                                                      !
!  PIC particles may end up in ghost cells if the cell adjacent to the !
!  ghost cell is a cut-cell the particle was not detected outside the  !
!  system because of tolerances.                                       !
!                                                                      !
!  Future Work:                                                        !
!                                                                      !
!  1. Redistribute particle's weight among other particles in the      !
!     domain to conserve mass.                                         !
!                                                                      !
!  2. Rather than deactivating the particle, reflect the particle      !
!     inside the domain using the ghost cell bc's                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_CELL_MOVEMENT_PIC

! Global Variables:
!---------------------------------------------------------------------//
! The (max) number of particles in process
      use discretelement, only: PIP, MAX_PIP
! The I/J/K/IJK indicies of the fluid cell
      use discretelement, only: PIJK
! The number of particles in cell IJK
      use discretelement, only: PINC
! Particle positions New/Previous
      use discretelement, only: DES_POS_NEW, DES_POS_OLD
! The East/North/Top face loctions of fluid cells
      use discretelement, only: XE, YN, ZT
! Particle velocities
      use discretelement, only: DES_VEL_NEW
! Flag: identifies a fluid cell as a cut cell.
      use cutcell, only: CUT_CELL_AT

      use mpi_utility
      use error_manager
      use functions

      IMPLICIT NONE

! Local Variables:
!----------------------------------------------------------------------!
! particle no.
      INTEGER :: L, I, J, K, IJK
! Integer error flag.
      INTEGER :: IER
! Position
      DOUBLE PRECISION :: lPOS
! Number of deleted particles found on the local process
      INTEGER :: lDELETED, gDELETED
! Number of recovered particles found on the local process
      INTEGER :: lRECOVERED, gRECOVERED
! Local parameter to print verbose messages about particles.
      LOGICAL, PARAMETER :: lDEBUG = .FALSE.

!.......................................................................


      CALL INIT_ERR_MSG("CHECK_CELL_MOVEMENT_PIC")
      IF(lDEBUG) CALL OPEN_PE_LOG(IER)


! Initialize the counters for deleted/recovered parcels.
      lDELETED = 0
      lRECOVERED = 0

      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(.NOT.IS_NORMAL(L)) CYCLE

! Assigning local aliases for particle position

! assigning local aliases for particle i, j, k fluid grid indices
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)

! this ijk is still an old value as it has not been updated
         IJK=PIJK(L,4)

! in MPPIC a particle can lie on the surface of the wall as only the
! centers are tracked.
         IF (I.GT.IEND1 .OR. I.LT.ISTART1) THEN

            lPOS = DES_POS_NEW(1,L)
            IF(I.EQ.IEND1+1 .AND. &
               (lPOS >= XE(IEND1-1) .AND. lPOS <= XE(IEND1)) )THEN

               lRECOVERED = lRECOVERED + 1
               PIJK(L,1) = IEND1

               IF(lDEBUG) THEN
                  WRITE(ERR_MSG,1100) trim(iVal(L)),'I',trim(iVal(I)), &
                  'X',DES_POS_OLD(1,L),'X',lPOS,'X',DES_VEL_NEW(1,L)
                  CALL FLUSH_ERR_MSG
               ENDIF
            ELSE

               lDELETED = lDELETED + 1
               CALL SET_NONEXISTENT(L)
               PINC(IJK) = PINC(IJK) - 1

               IF(lDEBUG) THEN
                  WRITE(ERR_MSG,1110) trim(iVal(L)),'I',trim(iVal(I)), &
                  'X',DES_POS_OLD(1,L),'X',lPOS,'X',DES_VEL_NEW(1,L),  &
                  trim(iVal(IJK)), CUT_CELL_AT(IJK), FLUID_AT(IJK)
                  CALL FLUSH_ERR_MSG
               ENDIF
               CYCLE
            ENDIF
         ENDIF

         IF(J.GT.JEND1 .OR. J.LT.JSTART1) THEN
            lPOS = DES_POS_NEW(2,L)
            IF(J.EQ.JEND1+1.AND.&
              (lPOS >= YN(JEND1-1) .AND. lPOS <= YN(JEND1)) ) THEN

               lRECOVERED = lRECOVERED + 1
               PIJK(L,2) = JEND1

               IF(lDEBUG) THEN
                  WRITE(ERR_MSG,1100) trim(iVal(L)),'J',trim(iVal(J)), &
                  'Y',DES_POS_OLD(2,L),'Y',lPOS,'Y',DES_VEL_NEW(2,L)
                  CALL FLUSH_ERR_MSG
               ENDIF

            ELSE

               lDELETED = lDELETED + 1
               CALL SET_NONEXISTENT(L)
               PINC(IJK) = PINC(IJK) - 1

               IF(lDEBUG) THEN
                  WRITE(ERR_MSG,1110) trim(iVal(L)),'J',trim(iVal(J)), &
                  'Y',DES_POS_OLD(2,L),'Y',lPOS,'Y',DES_VEL_NEW(2,L),  &
                  trim(iVal(IJK)), CUT_CELL_AT(IJK), FLUID_AT(IJK)
                  CALL FLUSH_ERR_MSG
               ENDIF
               CYCLE
            ENDIF
         ENDIF

         IF(DO_K .AND. (K > KEND1 .OR. K < KSTART1)) THEN
            lPOS = DES_POS_NEW(3,L)
            IF(K == KEND1+1 .AND. &
              (lPOS >= ZT(KEND1-1) .AND. lPOS <= ZT(KEND1)) ) THEN

               lRECOVERED = lRECOVERED + 1
               PIJK(L,3) = KEND1

               IF(lDEBUG) THEN
                  WRITE(ERR_MSG,1100) trim(iVal(L)),'K',trim(iVal(K)), &
                  'Z',DES_POS_OLD(3,L),'Z',lPOS,'Z',DES_VEL_NEW(3,L)
                  CALL FLUSH_ERR_MSG
               ENDIF
            ELSE

               lDELETED = lDELETED + 1
               CALL SET_NONEXISTENT(L)
               PINC(IJK) = PINC(IJK) - 1

               IF(lDEBUG) THEN
                  WRITE(ERR_MSG,1110) trim(iVal(L)),'K',trim(iVal(K)), &
                  'Z',DES_POS_OLD(3,L),'Z',lPOS,'Z',DES_VEL_NEW(3,L),  &
                  trim(iVal(IJK)), CUT_CELL_AT(IJK), FLUID_AT(IJK)
                  CALL FLUSH_ERR_MSG
               ENDIF
               CYCLE
            ENDIF
         ENDIF
      ENDDO

 1100 FORMAT('Warninge 1100: Particle ',A,' was recovered from a ',    &
         'ghost cell.',2/2x,'Moved into cell with ',A1,' index: ',A,   &
         /2x,A1,'-Position OLD:',g11.4,/2x,A1,'-Position NEW:',g11.4,  &
         /2x,A1,'-Velocity:',g11.4)

 1110 FORMAT('Warninge 1110: Particle ',A,' was deleted from a ',      &
         'ghost cell.',2/2x,'Moved into cell with ',A1,' index: ',A,   &
         /2x,A1,'-Position OLD:',g11.4,/2x,A1,'-Position NEW:',g11.4,  &
         /2X,A1,'-Velocity:',g11.4,/2x,'Fluid Cell: ',A,/2x,           &
         'Cut cell? ',L1,/2x,'Fluid at? ',L1)

! Update the number of particles
      PIP = PIP - lDELETED

! Send the numbers to the IO process.
      CALL GLOBAL_SUM(lRECOVERED, gRECOVERED, PE_IO)
      CALL GLOBAL_SUM(lDELETED, gDELETED, PE_IO)

      IF(gRECOVERED + gDELETED > 0) THEN
         WRITE(ERR_MSG,1115) trim(iVal(gDELETED + gRECOVERED)),        &
            trim(iVal(gDELETED)), trim(iVal(gRECOVERED))
         CALL FLUSH_ERR_MSG
      ENDIF

 1115 FORMAT('Warning 1115: ',A,' particles detected outside the ',    &
         'domain.',/2x,A,' particles were deleted.',/2x,A,' particles',&
         ' were recovered.')

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_CELL_MOVEMENT_PIC


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: REPORT_PIC_STATS                                        !
!                                                                      !
!  Purpose: Output stats about PIC simulation.                         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE REPORT_PIC_STATS

! Global Variables:
!---------------------------------------------------------------------//
! Flag to report minimum EP_G
      use mfix_pic, only: PIC_REPORT_MIN_EPG
! Gas phase volume fraction
      use fldvar, only: EP_G
! Location of cell faces (East, North, Top)
      use discretelement, only: XE, YN, ZT

      use param1, only: large_number
      use mpi_utility
      USE error_manager
      USE functions

      IMPLICIT NONE

! Local Variables:
!----------------------------------------------------------------------!
! Loop counters
      INTEGER I, J, K, IJK, IPROC

      INTEGER :: EPg_MIN_loc(0:numpes-1, 4), EPg_MIN_loc2(1)
      DOUBLE PRECISION :: EPg_MIN(0:numpes-1), EPg_min2

!-----------------------------------------------

      CALL INIT_ERR_MSG("REPORT_PIC_STATS")


      IF(PIC_REPORT_MIN_EPG) THEN

         EPG_MIN(:) = 0
         EPG_MIN(mype) = LARGE_NUMBER

         EPG_MIN_LOC(:,:) = 0
         EPG_MIN_LOC(mype,:) = -1

         DO K = KSTART1, KEND1
            DO J = JSTART1, JEND1
               DO I = ISTART1, IEND1
                  IJK = funijk(I,J,K)

                  IF(EP_G(IJK) < EPG_MIN(mype)) THEN
                     EPG_MIN_LOC(mype,1) = I
                     EPG_MIN_LOC(mype,2) = J
                     EPG_MIN_LOC(mype,3) = K
                     EPG_MIN_LOC(mype,4) = IJK
                     EPG_MIN(mype) = EP_G(IJK)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         call GLOBAL_ALL_SUM(EPg_MIN)
         CALL GLOBAL_ALL_SUM(EPg_MIN_loc)

         epg_min2     = MINVAL(epg_min(0:numpes-1))
         epg_min_loc2 = MINLOC(epg_min(0:numpes-1)) - 1
         !-1, since minloc goes from 1:size of the array.
         !If not corrected by -1, then the proc id will be off by 1

         iproc = epg_min_loc2(1)

         I     = epg_min_loc(iproc, 1)
         J     = epg_min_loc(iproc, 2)
         K     = epg_min_loc(iproc, 3)
         IJK   = epg_min_loc(iproc, 4)
         WRITE(ERR_MSG,1014) EPG_MIN2, Iproc, I, J, K, IJK, &
            XE(I) - 0.5*DX(I), YN(J)-0.5*DY(J), ZT(K) - 0.5*DZ(K)

 1014       FORMAT( /, &
            &      5x,'EPGMIN                    = ', 2x,g17.8,/ &
            &      5x,'EPGMIN PROC RANK          = ', 2x, I10, / &
            &      5x,'EPGMIN (I, J, K, IJK)     = ', 3(2x,i5),2x,i10,/ &
            &      5x,'XMID, YMID, ZMID FOR CELL = ', 3(2x,g17.8))

            call flush_err_msg(header = .false., footer = .false.)

      ENDIF

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE REPORT_PIC_STATS


