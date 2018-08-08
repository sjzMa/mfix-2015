!----------------------------------------------------------------------!
! Module: MPI_FUNS_DES                                                 !
! Author: Pradeep Gopalakrishnan                                       !
!                                                                      !
! Purpose: This module contains the subroutines and functions for MPI  !
! communications in DES simulations.                                   !
!----------------------------------------------------------------------!
      module mpi_funs_des

      contains

!----------------------------------------------------------------------!
! Subroutine: DES_PAR_EXCHANGE                                         !
! Author: Pradeep Gopalakrishnan                                       !
!                                                                      !
! Purpose: This subroutine controls entire exchange of particles       !
!    between processors.                                               !
!                                                                      !
! Steps:                                                               !
! 1) Bin the particles to the DES grid.                                !
! 2) Check if the send and recv buffer size is large enough            !
! 3) Pack and send active particles located in ghost cells to the      !
!    processors that own the ghost cells. The exchange occurs in       !
!    the following order to take care of particles crossing at corners !
!    (e.g., crossing directly into the northwest block):               !
!    a.) top-bottom interface                                          !
!    b.) north-south interface                                         !
!    c.) east-west interface                                           !
! 4) Bin the particles (if required)                                   !
! 5) Pack and send particles adjacent to neighboring processes. The    !
!    exchange occurs in the following order:                           !
!    a.) east-west interface                                           !
!    b.) north-south interface                                         !
!    c.) top-bottom interface                                          !
!                                                                      !
! Comments: The DO_NSEARCH flag should be set before calling           !
!   DES_PAR_EXCHANGE; When DO_NSEARCH is true, ghost particles are     !
!   updated and later used  to generate the PAIR lists.                !
!----------------------------------------------------------------------!
      subroutine des_par_exchange()

      use discretelement, only: DO_NSEARCH, DES_PERIODIC_WALLS

      use mfix_pic, only: MPPIC
      use desmpi, only: iEXCHFLAG
      use desmpi, only: dSENDBUF, dRECVBUF
      use discretelement, only: iGHOST_UPDATED
      use desmpi, only: iSPOT

      use geometry, only: NO_K

! Module procedures
!---------------------------------------------------------------------//
      use mpi_pack_des, only: desmpi_pack_parcross
      use mpi_unpack_des, only: desmpi_unpack_parcross

      use mpi_pack_des, only: desmpi_pack_ghostpar
      use mpi_unpack_des, only: desmpi_unpack_ghostpar

      use mpi_comm_des, only: desmpi_sendrecv_init
      use mpi_comm_des, only: desmpi_sendrecv_wait

      use desgrid, only: desgrid_pic
      use desmpi_wrapper, only: des_mpi_barrier
      use compar, only: numpes

      implicit none

! Local variables:
!---------------------------------------------------------------------//
! Loop counters.
      integer :: linter, lface, ii
! Number of calls since the buffer was last checked.
      integer, save :: lcheckbuf = 0

!......................................................................!

      IF (.not.((numPEs>1) .OR. DES_PERIODIC_WALLS)) RETURN

! Check that the send/recv buffer is sufficient every 100 calls to
! avoid the related global communications.
      if (mod(lcheckbuf,100) .eq. 0) then
         call desmpi_check_sendrecvbuf
         lcheckbuf = 0
      end if
      lcheckbuf = lcheckbuf + 1

! call particle crossing the boundary exchange in T-B,N-S,E-W order
      do ii=1, size(dsendbuf)
         dsendbuf(ii)%facebuf(1) = 0
         drecvbuf(ii)%facebuf(1) = 0
      end do

      ispot = 1
      do linter = merge(2,3,NO_K),1,-1
         do lface = linter*2-1,linter*2
            if(.not.iexchflag(lface))cycle
            call desmpi_pack_parcross(lface)
            call desmpi_sendrecv_init(lface)
         end do
         do lface = linter*2-1,linter*2
            if(.not.iexchflag(lface)) cycle
            call desmpi_sendrecv_wait(lface)
            call desmpi_unpack_parcross(lface)
         end do
! update pic this is required for particles crossing corner cells
         do lface = linter*2-1,linter*2
            if(dsendbuf(1+mod(lface,2))%facebuf(1).gt.0 .or. &
               drecvbuf(1+mod(lface,2))%facebuf(1).gt.0) then
               call desgrid_pic(plocate=.false.)
               exit
            end if
         end do
      end do
      call des_mpi_barrier

!      call des_dbgmpi(5)


      IF(.NOT.MPPIC) THEN
! call ghost particle exchange in E-W, N-S, T-B order

         do ii=1, size(dsendbuf)
            dsendbuf(ii)%facebuf(1) = 0
            drecvbuf(ii)%facebuf(1) = 0
         end do

         ighost_updated(:) = .false.
         ispot = 1
         do linter = 1,merge(2,3,NO_K)
            do lface = linter*2-1,linter*2
               if(.not.iexchflag(lface))cycle
               call desmpi_pack_ghostpar(lface)
               call desmpi_sendrecv_init(lface)
            end do
            do lface = linter*2-1,linter*2
               if(.not.iexchflag(lface)) cycle
               call desmpi_sendrecv_wait(lface)
               call desmpi_unpack_ghostpar(lface)
            end do

! Rebin particles to the DES grid as ghost particles may be moved.
            do lface = linter*2-1,linter*2
               if(dsendbuf(1+mod(lface,2))%facebuf(1) .gt.0.or.&
                  drecvbuf(1+mod(lface,2))%facebuf(1).gt.0) then
                  call desgrid_pic(plocate=.false.)
                  exit
               end if
            end do
         end do
         if(do_nsearch) call desmpi_cleanup
         call des_mpi_barrier
      ENDIF   ! end if(.not.mppic)

!      call des_dbgmpi(2)
!      call des_dbgmpi(3)
!      call des_dbgmpi(4)
!      call des_dbgmpi(6)
!      call des_dbgmpi(7)

      END SUBROUTINE DES_PAR_EXCHANGE


!----------------------------------------------------------------------!
! Subroutine: DESMPI_CHECK_SENDRECVBUF                                 !
! Author: Pradeep Gopalakrishnan                                       !
!                                                                      !
! Purpose: Checks if the sendrecvbuf size is large enough. If the      !
!    buffers are not sufficent, they are resized.                      !
!----------------------------------------------------------------------!
      SUBROUTINE DESMPI_CHECK_SENDRECVBUF

      use discretelement, only: DIMN, dg_pic
      use desmpi, only: iMAXBUF
      use desmpi, only: iBUFOFFSET
      use desmpi, only: dSENDBUF, dRECVBUF
      use desmpi, only: iSENDINDICES
      use desmpi, only: iGhostPacketSize

      use mpi_utility, only: global_all_max
      use error_manager
      implicit none

! Local variables:
!---------------------------------------------------------------------//
! Loop counters
      INTEGER :: lface, lindx, lijk
! Particle count in send/recv region on current face
      INTEGER :: lparcnt
! Max particle count in send/recv region over all proc faces.
      INTEGER :: lmaxcnt
! Total number of DES grid cells on lface in send/recv
      INTEGER :: ltot_ind
! Previous Buffer
      INTEGER :: pBUF
! Growth factor when resizing send/recv buffers.
      REAL :: lfactor = 1.5
      DOUBLE PRECISION, PARAMETER :: ONEMBo8 = 131072.0
!......................................................................!

      lmaxcnt = 0
      do lface = 1,6
         ltot_ind = isendindices(1,lface)
         lparcnt = 0
         do lindx = 2,ltot_ind+1
            lijk = isendindices(lindx,lface)
            lparcnt = lparcnt + dg_pic(lijk)%isize
         enddo
         if(lparcnt.gt.lmaxcnt) lmaxcnt = lparcnt
      enddo

      call global_all_max(lmaxcnt)
      if (imaxbuf .lt. lmaxcnt*iGhostPacketSize+ibufoffset) then
         pbuf = imaxbuf
         imaxbuf = lmaxcnt*iGhostPacketSize*lfactor
         do lface = 1,2*dimn
            if(allocated(dsendbuf(1+mod(lface,2))%facebuf)) then
               deallocate(dsendbuf(1+mod(lface,2))%facebuf,drecvbuf(1+mod(lface,2))%facebuf)
            endif
            allocate(dsendbuf(1+mod(lface,2))%facebuf(imaxbuf),&
               drecvbuf(1+mod(lface,2))%facebuf(imaxbuf))
         end do

         WRITE(ERR_MSG, 1000) iMAXBUF/ONEMBo8, &
            100.0d0+100.0d0*dble(iMAXBUF-pbuf)/dble(pbuf)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 1000 FORMAT(/'Resizeing DES MPI buffers: ',F7.1,' MB  (+',F5.1'%)')

      endif

      END SUBROUTINE DESMPI_CHECK_SENDRECVBUF

!----------------------------------------------------------------------!
! Subroutine: DESMPI_CLEANUP                                           !
! Author: Pradeep Gopalakrishnan                                       !
!                                                                      !
! Purpose: Cleans the ghost particle array positions.                  !
!----------------------------------------------------------------------!
      SUBROUTINE DESMPI_CLEANUP

      use discretelement, only: DIMN
      use discretelement, only: DES_POS_NEW, DES_POS_OLD
      use discretelement, only: DES_VEL_NEW, DES_VEL_OLD
      use discretelement, only: OMEGA_NEW
      use discretelement, only: PARTICLE_ORIENTATION,ORIENTATION,INIT_ORIENTATION
      use discretelement, only: FC
      use discretelement, only: DO_OLD
      use discretelement, only: PIP
      use discretelement, only: iGHOST_CNT
      use discretelement, only: DES_USR_VAR_SIZE, DES_USR_VAR
      use discretelement, only: dg_pic, pijk

      use run, only: ENERGY_EQ,ANY_SPECIES_EQ
      use des_thermo, only: DES_T_s_NEW, DES_T_s_OLD

      use des_rxns, only: DES_X_s

      use discretelement, only: iGHOST_UPDATED
      use functions, only: SET_NONEXISTENT
      use desmpi, only: iRECVINDICES
      use desmpi, only: iEXCHFLAG

      use param, only: DIMENSION_N_s
      use param1, only: ZERO

      implicit none

! Local variables:
!---------------------------------------------------------------------//
      integer ltot_ind,lface,lindx,lijk,lcurpar,lpicloc

      do lface = 1,dimn*2
         if(.not.iexchflag(lface))cycle
         ltot_ind = irecvindices(1,lface)
         do lindx = 2,ltot_ind+1
            lijk = irecvindices(lindx,lface)
            do lpicloc =1,dg_pic(lijk)%isize
               lcurpar = dg_pic(lijk)%p(lpicloc)
               if(ighost_updated(lcurpar)) cycle
               pip = pip - 1
               ighost_cnt = ighost_cnt-1
               call set_nonexistent(lcurpar)
               fc(:,lcurpar) = 0.0
               des_pos_new(:,lcurpar)=0
               pijk(lcurpar,:) = -10
               IF (DO_OLD) THEN
                  des_pos_old(:,lcurpar)=0
                  des_vel_old(:,lcurpar)=0
               ENDIF
               des_vel_new(:,lcurpar)=0
               omega_new(:,lcurpar)=0

               IF(PARTICLE_ORIENTATION) ORIENTATION(1:3,lcurpar) = INIT_ORIENTATION

               if(ENERGY_EQ) then
                  des_t_s_new(lcurpar)=0
                  des_t_s_old(lcurpar)=0
               endif

               if(ANY_SPECIES_EQ)then
                  des_x_s(lcurpar,1:dimension_n_s)= 0
               endif

               IF(DES_USR_VAR_SIZE > 0)&
                  des_usr_var(:,lcurpar)= 0

            end do
         end do
      end do
      END SUBROUTINE DESMPI_CLEANUP


      END MODULE MPI_FUNS_DES
