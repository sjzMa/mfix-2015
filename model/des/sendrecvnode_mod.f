!------------------------------------------------------------------------
! Module           : sendrecvnode
! Purpose          : Created to send and recv values at grid nodes
!                    used for interpolation variables stored at grid nodes
!
! Author           : Pradeep.G
! Comments         : Node values used for interpolation is from istart2
!                    to iend1; mapping technique used in sendrecv mod is
!                    not applicable to nodes; hence different techniques is
!                    employed
! Contains subroutines:
!    des_setnodeindices, des_exchangenode, des_dbgnodesr
!------------------------------------------------------------------------
      module sendrecvnode

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use parallel_mpi
      use mpi_utility
      use discretelement
      use compar
      use physprop
      use sendrecv
      use desmpi_wrapper
      use desgrid
      use functions
!-----------------------------------------------

      integer :: itotalneigh,itotalindx
      integer,dimension(:),allocatable ::  itoproc,iprocsumindx,istartindx
! Following variables are used to exchange grid index values when
! des_interp_on is true
      integer,dimension(:),allocatable :: isendnodes
      double precision,dimension(:),allocatable :: dsendnodebuf,drecvnodebuf
      integer,dimension(:),allocatable :: irecvreqnode,isendreqnode


      contains


!------------------------------------------------------------------------
! Subroutine       : des_setnodesindices
! Purpose          : allocates and initializes the variables related
!                    to send and recv for grid node values
!------------------------------------------------------------------------
      subroutine des_setnodeindices
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer lijkproc,liproc,ljproc,lkproc
      integer li,lj,lk
      integer li2,lj2,lk2
      integer liproc_start,liproc_end,ljproc_start,ljproc_end,lkproc_start,lkproc_end
      integer lci,lcj,lck,lproc,lcount
      integer linode,ljnode,lknode
      integer linode_start,linode_end,ljnode_start,ljnode_end,lknode_start,lknode_end
      logical lpresent
!-----------------------------------------------

! set flags for interprocessor boundaries and set the corresponding to proc
      liproc = iofproc(mype)
      ljproc = jofproc(mype)
      lkproc = kofproc(mype)

! if not periodic then limit the processor
      if(des_periodic_walls_x.and.nodesi.gt.1) then
         liproc_start=liproc-1
         liproc_end=liproc+1
      else
         liproc_start =max(liproc-1,0)
         liproc_end=min(liproc+1,nodesi-1)
      end if
      if(des_periodic_walls_y.and.nodesj.gt.1) then
         ljproc_start=ljproc-1
         ljproc_end=ljproc+1
      else
         ljproc_start =max(ljproc-1,0)
         ljproc_end=min(ljproc+1,nodesj-1)
      end if
      if(des_periodic_walls_z.and.nodesk.gt.1) then
         lkproc_start=lkproc-1
         lkproc_end=lkproc+1
      else
         lkproc_start =max(lkproc-1,0)
         lkproc_end=min(lkproc+1,nodesk-1)
      end if
      itotalneigh = (liproc_end-liproc_start+1)*(ljproc_end-ljproc_start+1)* &
                   (lkproc_end-lkproc_start+1)-1
! allocate the variables
      allocate (itoproc(itotalneigh),iprocsumindx(itotalneigh),istartindx(itotalneigh+1), &
                irecvreqnode(itotalneigh),isendreqnode(itotalneigh))

! First loop to count the total index for each processor and count the
! neighbour processor
      itotalneigh = 0
      itoproc(:)=-1
      iprocsumindx(:) =0
      do lk = lkproc_start,lkproc_end
      do lj = ljproc_start,ljproc_end
      do li = liproc_start,liproc_end
         li2 = mod(li,nodesi);if(li2.lt.0)li2=nodesi-1
         lj2 = mod(lj,nodesj);if(lj2.lt.0)lj2=nodesj-1
         lk2 = mod(lk,nodesk);if(lk2.lt.0)lk2=nodesk-1
         lijkproc = procijk(li2,lj2,lk2)
         if (lijkproc.eq.mype) cycle
! check if the processor exits in the previous list
         lpresent = .false.
         do lproc = 1,itotalneigh
            if (lijkproc .eq.itoproc(lproc)) then
               lpresent = .true.
               exit
            end if
         end do
         if(.not.lpresent) then
            itotalneigh = itotalneigh + 1
            lproc = itotalneigh
         end if
         itoproc(lproc) = lijkproc
         lci=(liproc-li);lcj=(ljproc-lj);lck=(lkproc-lk)
         linode_start = istart2; linode_end=iend1
         ljnode_start = jstart2; ljnode_end=jend1
         lknode_start = kstart2; lknode_end=kend1
         if(lci.eq.1) linode_end = istart2
         if(lci.eq.-1)  linode_start = iend1
         if(lcj.eq.1) ljnode_end = jstart2
         if(lcj.eq.-1)  ljnode_start = jend1
         if(lck.eq.1) lknode_end = kstart2
         if(lck.eq.-1)  lknode_start = kend1
         do lknode = lknode_start,lknode_end
         do linode = linode_start,linode_end
         do ljnode = ljnode_start,ljnode_end
            IF(DEAD_CELL_AT(linode,ljnode,lknode)) CYCLE
            iprocsumindx(lproc) = iprocsumindx(lproc) + 1
         end do
         end do
         end do
      end do
      end do
      end do
!assign the start index
      do lproc =1,itotalneigh+1
         istartindx(lproc)=sum(iprocsumindx(1:lproc-1))+1
      end do
      itotalindx=istartindx(itotalneigh+1)-1

! allocate the variables
      allocate(isendnodes(itotalindx),dsendnodebuf(itotalindx),drecvnodebuf(itotalindx))

! second loop to assign actual index
      iprocsumindx(:)=0
      do lk = lkproc_start,lkproc_end
      do lj = ljproc_start,ljproc_end
      do li = liproc_start,liproc_end
         li2 = mod(li,nodesi);if(li2.lt.0)li2=nodesi-1
         lj2 = mod(lj,nodesj);if(lj2.lt.0)lj2=nodesj-1
         lk2 = mod(lk,nodesk);if(lk2.lt.0)lk2=nodesk-1
         lijkproc = procijk(li2,lj2,lk2)
         if (lijkproc.eq.mype) cycle
! find the index of the neighbour
         do lproc =1,itotalneigh
            if(lijkproc.eq.itoproc(lproc)) then
               exit
            end if
         end do
         lci=(liproc-li);lcj=(ljproc-lj);lck=(lkproc-lk)
         linode_start = istart2; linode_end=iend1
         ljnode_start = jstart2; ljnode_end=jend1
         lknode_start = kstart2; lknode_end=kend1
         if(lci.eq.1) linode_end = istart2
         if(lci.eq.-1)  linode_start = iend1
         if(lcj.eq.1) ljnode_end = jstart2
         if(lcj.eq.-1)  ljnode_start = jend1
         if(lck.eq.1) lknode_end = kstart2
         if(lck.eq.-1)  lknode_start = kend1
         lcount = istartindx(lproc)+iprocsumindx(lproc)
         do lknode = lknode_start,lknode_end
         do linode = linode_start,linode_end
         do ljnode = ljnode_start,ljnode_end
            IF(DEAD_CELL_AT(linode,ljnode,lknode)) CYCLE
            isendnodes(lcount)=funijk(linode,ljnode,lknode)
            iprocsumindx(lproc)=iprocsumindx(lproc)+1
            lcount = lcount+1
         end do
         end do
         end do
      end do
      end do
      end do

!      call  des_dbgnodesr()
      end subroutine des_setnodeindices

!------------------------------------------------------------------------
! Subroutine       : des_exchangenode
! Purpose          : calls send and recv to exchange the node values and
!                    adds based on the flag
!                    to send and recv for grid node values
! Parameters       : pvar - variable that has to be exchanged
!                    padd - if true node values will be added
!                           else node values will be replaced
!------------------------------------------------------------------------
      subroutine des_exchangenode(pvar,padd)
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      double precision,dimension(:),intent(inout) ::pvar
      logical:: padd
!-----------------------------------------------
! local variables
!-----------------------------------------------
      character(len=80), parameter :: name = 'des_exchangenode'
      integer :: lindx,lcount,lcount2,lneigh,ltag,lerr
      integer :: lstart,lend,ltotal
!-----------------------------------------------

! steps pack the buffer call isend and irecv
      do lcount = 1,itotalneigh
         lneigh = itoproc(lcount)
         lstart = istartindx(lcount);lend=istartindx(lcount+1)-1
         do lcount2 = lstart,lend
            dsendnodebuf(lcount2) = pvar(isendnodes(lcount2))
         end do
         ltag = message_tag(lneigh,mype)
         ltotal = lend-lstart+1
         call des_mpi_irecv(drecvnodebuf(lstart:lend),ltotal, &
                            lneigh,ltag,irecvreqnode(lcount),lerr)
         call mpi_check( name //':mpi_irecv ', lerr )
         ltag = message_tag(mype,lneigh)
         call des_mpi_isend(dsendnodebuf(lstart:lend),ltotal, &
                            lneigh,ltag,isendreqnode(lcount),lerr)
         call mpi_check( name //':mpi_irecv ', lerr )
      end do
! call mpi wait to complete the exchange
      do lcount = 1,itotalneigh
         call des_mpi_wait(isendreqnode(lcount),lerr)
         call mpi_check( name //':mpi_wait-send', lerr )
         call des_mpi_wait(irecvreqnode(lcount),lerr)
         call mpi_check( name //':mpi_wait-recv', lerr )
      end do
! after receiving the buffer the values are either added or
! replaced based on the flag
      if (padd) then
         do lcount = 1,itotalindx
            lindx = isendnodes(lcount)
            pvar(lindx) = pvar(lindx) + drecvnodebuf(lcount)
         end do
      else
         do lcount = 1,itotalindx
            lindx = isendnodes(lcount)
            pvar(lindx) = drecvnodebuf(lcount)
         end do
      end if
      return

      contains

        integer function message_tag(lsource,ldest)
          implicit none
          integer, intent(in) :: lsource,ldest
          message_tag = lsource+numpes*ldest+200
        end function message_tag

      end subroutine des_exchangenode

!------------------------------------------------------------------------
! Subroutine       : des_dbgnodesr
! Purpose          : For debugging prints the indices
!------------------------------------------------------------------------
      subroutine des_dbgnodesr()
!-----------------------------------------------

      use functions

      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      character (255) filename
      integer ijk
      integer lcount,lcount2,lstart,lend
!-----------------------------------------------

! pradeep remove print the flags
      write(filename,'("dbg_nodesr",I4.4,".dat")') mype
      open(44,file=filename,convert='big_endian')
      do lcount = 1,itotalneigh
         lstart = istartindx(lcount);lend=istartindx(lcount+1)-1
         write(44,*) "--------------------------------------------"
         write(44,*) "to proc start end",itoproc(lcount),lstart,lend
         write(44,*) "--------------------------------------------"
         do lcount2 = lstart,lend
          ijk = isendnodes(lcount2)
          write(44,*)ijk,i_of(ijk),j_of(ijk),k_of(ijk)
         end do
         write(44,*) "--------------------------------------------"
      end do
      close (44)
      end subroutine des_dbgnodesr

      end module


