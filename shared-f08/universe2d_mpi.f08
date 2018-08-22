   module parallel_universe
   use mpi_f08
   use universe
   use patches, only: optimal_gridmap, gridmap
   implicit none
   private
   public pp2d_mpi, pr, particle, pp2d

      type(mpi_datatype)      :: particle_mpi_type

      type, extends(pp2d)     :: pp2d_mpi
         !private
         integer, public      :: size, rank
         type(gridmap)        :: map 
         integer              :: offset(2)
         integer, allocatable :: owner(:), old_owner(:)
         integer              :: c_mine(2), w_mine(2) 
         integer              :: n_mine, n_borders, mincap
      contains
         !private
         procedure            :: mod2, cl_on
         procedure, public    :: randcc, unique_rnd, fixed_rnd
         procedure, public    :: start_parallel, end_parallel
         procedure, private   :: isend_irecv
         procedure, public    :: shift_blocks
         procedure, public    :: stage_push
         procedure, public    :: push_blocks, push_borders
      end type pp2d_mpi
   
   contains

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
   
      function mod2(pos,c1,cw) result(c2)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos 
      integer,         intent(in)    :: c1(2)
      integer, intent(in), optional  :: cw(2)
      integer                        :: c2(2)
      c2 = c1
      if( present(cw) ) c2 = c2 + cw
      c2 = mod(c2,pos%nxy) 
      where(c2<0) c2 = c2+pos%nxy
      end function mod2

      function cl_on(pos,c1,cw) result(cl)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos 
      integer,         intent(in)    :: c1(2)
      integer, intent(in), optional  :: cw(2)
      integer                        :: cl, c2(2)
      c2 = c1
      if( present(cw) ) c2 = c2 + cw
      c2 = pos%mod2(c2)
      cl = c2(2)*pos%nx + c2(1)
      end function cl_on

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      function randcc(pos) result(cc)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos 
      integer                        :: cc(2), ierr
      real(pr)                       :: rnd(2)
      call random_number(rnd)
      cc = floor(pos%map%cw(:,0)*rnd)
      call mpi_bcast( cc, 2, mpi_integer, pos%size-1, mpi_comm_world, ierr)
      end function randcc

      subroutine unique_rnd(pos)
      implicit none
      class(pp2d_mpi), intent(in) :: pos
      integer                     :: n
      integer,  allocatable       :: seed(:)
      integer                     :: time
      call random_seed(size=n)
      allocate(seed(n))
      call system_clock(time)
      seed = time + pos%rank
      call random_seed(put=seed)
      end subroutine unique_rnd

      subroutine fixed_rnd(pos)
      implicit none
      class(pp2d_mpi), intent(in) :: pos
      integer                     :: n
      integer,  allocatable       :: seed(:)
      call random_seed(size=n)
      allocate(seed(n))
      seed = 666777657 + pos%rank
      call random_seed(put=seed)
      end subroutine fixed_rnd

      !//////////////////////////////////////////////////
      !////////////////////////////////////////////////// comm
      !//////////////////////////////////////////////////

      subroutine make_particle_mpi_type()
      implicit none
      type(particle)            :: pa
      type(mpi_datatype)        :: datatypes(2)
      integer(mpi_address_kind) :: locs(2)
      integer                   :: ierr, lengths(2)
      call mpi_get_address(pa%name,locs(1),ierr); lengths(1) = 1; datatypes(1) =  mpi_integer
      call mpi_get_address(pa%pos, locs(2),ierr); lengths(2) = 2; datatypes(2) =  mpi_real
      locs = locs - locs(1)
      call mpi_type_create_struct(2,lengths,locs,datatypes,particle_mpi_type,ierr)
      call mpi_type_commit(particle_mpi_type,ierr)
      end subroutine make_particle_mpi_type

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine start_parallel(pos,cpu_grid)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos
      integer, intent(in), optional  :: cpu_grid      
      integer                        :: ierr, g(2)
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,pos%rank,ierr)
      call mpi_comm_size(mpi_comm_world,pos%size,ierr)
      call make_particle_mpi_type()
      if( present(cpu_grid) ) then
         g = cpu_grid
      else
         call optimal_gridmap(pos%nxy,pos%size,g)
      end if
      pos%map = gridmap(pos%nxy,g)
      pos%offset = 0
      allocate( pos%owner(0:pos%noc-1), &
                 pos%old_owner(0:pos%noc-1) )
      pos%owner = 0 
      call pos%shift_blocks([0,0])
      end subroutine start_parallel

   
      subroutine end_parallel(pos)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos
      integer                        :: ierr
      call mpi_finalize(ierr)
      deallocate( pos%owner, pos%old_owner )
      end subroutine end_parallel

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine shift_blocks(pos,shift)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos 
      integer,           intent(in)  :: shift(2)
      integer                        :: i, cl, cx, cy, new_owner, &
                                        count_ranks(0:pos%size-1)
      pos%offset = pos%mod2( pos%offset, shift )
      count_ranks = 0
      do i = 0, pos%map%n-1
         if( i<pos%size ) then
                 new_owner = i
         else
                 new_owner = 0
         end if
         do cy = pos%map%c1(2,i), pos%map%c2(2,i)
            do cx = pos%map%c1(1,i), pos%map%c2(1,i)
               cl = pos%cl_on( [cx,cy], pos%offset )
               pos%old_owner(cl) = pos%owner(cl)
               pos%owner(cl) = new_owner 
               count_ranks(new_owner) = count_ranks(new_owner) + 1
            end do
         end do
      end do
      pos%c_mine = pos%mod2( pos%map%c1(:,pos%rank), pos%offset )
      pos%w_mine = pos%map%cw(:,pos%rank)
      pos%n_mine = count_ranks(pos%rank)
      pos%n_borders = 2*(pos%map%nxt*pos%ny+pos%map%nyt*pos%nx) + 4*pos%map%n
      end subroutine shift_blocks

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine stage_push(pos)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos 
      integer                        :: cl, maxocc, mincap, ierr
      maxocc = maxval(pos%cup(:)%occ)
      ! prepare for sendrecv
      call mpi_allreduce(maxocc,mincap,1,mpi_integer,mpi_max,mpi_comm_world,ierr)
      do cl = 0, pos%noc-1
         do
            if( pos%cup(cl)%cap >= mincap ) exit
            call pos%cup(cl)%ext()
         end do
      end do
      pos%mincap = mincap
      end subroutine stage_push


      subroutine isend_irecv(pos,cl,rk1,rk2,ns,requests)
      implicit none
      class(pp2d_mpi), intent(inout)   :: pos 
      integer,         intent(in)      :: cl, rk1, rk2
      integer,          intent(inout)  :: ns
      type(mpi_request), intent(inout) :: requests(:)
      integer                          :: ierr
      if( rk1==rk2 ) return
      if( pos%rank==rk2 ) then
         ns = ns + 1
         call mpi_irecv( pos%cup(cl)%occ, 1, mpi_integer, rk1, pos%noc+cl, mpi_comm_world, requests(ns), ierr)
         ns = ns + 1
         call mpi_irecv( pos%cup(cl)%hld(1)%name, pos%mincap, particle_mpi_type, rk1, cl, mpi_comm_world, requests(ns), ierr)
      elseif( pos%rank==rk1 ) then
         ns = ns + 1
         call mpi_isend( pos%cup(cl)%occ, 1, mpi_integer, rk2, pos%noc+cl, mpi_comm_world, requests(ns), ierr)
         ns = ns + 1
        call mpi_isend( pos%cup(cl)%hld(1)%name, pos%cup(cl)%occ, particle_mpi_type, rk2, cl, mpi_comm_world, requests(ns), ierr)
      end if
      end subroutine isend_irecv


      subroutine push_blocks(pos)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos 
      integer                        :: cl, ierr, ns
      type(mpi_request)              :: requests(4*pos%n_mine)
      call pos%stage_push()
      ns = 0
      do cl = 0, pos%noc-1
         call pos%isend_irecv( cl, pos%old_owner(cl), pos%owner(cl), ns, requests )
      end do
      call mpi_waitall(ns,requests(1:ns),mpi_statuses_ignore,ierr)
      end subroutine push_blocks

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine push_borders(pos)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos
      integer                        :: i, j, c1, c2, cx, cy, ns, ierr
      type(mpi_request)              :: requests( pos%n_borders )
      call pos%stage_push()
      ns = 0
      do i = 1, pos%map%nxt
         cx = pos%map%xticks(i)
         do cy = 0, pos%nxy(2)-1
            c1 = pos%cl_on( [cx,cy],   pos%offset ) ! right
            c2 = pos%cl_on( [cx-1,cy], pos%offset ) ! left
            call pos%isend_irecv( c2, pos%owner(c2), pos%owner(c1), ns, requests )
            call pos%isend_irecv( c1, pos%owner(c1), pos%owner(c2), ns, requests )
         end do
      end do
      do i = 1, pos%map%nyt
         cy = pos%map%yticks(i)
         do cx = 0, pos%nxy(1)-1
            c1 = pos%cl_on( [cx,cy],   pos%offset ) ! up
            c2 = pos%cl_on( [cx,cy-1], pos%offset ) ! down
            call pos%isend_irecv( c2, pos%owner(c2), pos%owner(c1), ns, requests )
            call pos%isend_irecv( c1, pos%owner(c1), pos%owner(c2), ns, requests )
         end do
      end do
      do i = 1, pos%map%nxt
         cx = pos%map%xticks(i)
         do j = 1, pos%map%nyt
            cy = pos%map%yticks(j)
            c1 = pos%cl_on( [cx,cy],     pos%offset ) ! upright
            c2 = pos%cl_on( [cx-1,cy-1], pos%offset ) ! downleft
            call pos%isend_irecv( c2, pos%owner(c2), pos%owner(c1), ns, requests )
            call pos%isend_irecv( c1, pos%owner(c1), pos%owner(c2), ns, requests )
            c1 = pos%cl_on( [cx,cy-1], pos%offset )   ! downright
            c2 = pos%cl_on( [cx-1,cy], pos%offset )   ! upleft
            call pos%isend_irecv( c2, pos%owner(c2), pos%owner(c1), ns, requests )
            call pos%isend_irecv( c1, pos%owner(c1), pos%owner(c2), ns, requests )
         end do
      end do
      call mpi_waitall(ns,requests(1:ns),mpi_statuses_ignore,ierr)
      end subroutine push_borders

   end module parallel_universe


!   program test
!   use parallel_universe
!   use mpi_f08
!   implicit none
!   type(pp2d_mpi) :: pos
!   integer :: i, j, idx, dir
!   real(pr) :: delta
!   integer                   :: s_time, e_time, c_rate
!   pos%pp2d = pp2d([128,128],1.0,2.5)
!   call pos%start_parallel()
!   call pos%fixed_rnd()
!   ! randomly move particles
!   call system_clock(s_time,c_rate)
!   do j = 1, 1000
!      call pos%stage(pos%c_mine,pos%w_mine-1)
!      do i = 0, pos%nop
!         call pos%random(-0.1,0.1,idx,dir,delta)
!         call pos%move(idx,dir,delta) 
!      end do
!      call pos%shift_blocks([20,20])
!      call pos%push_blocks()
!      call pos%push_borders()
!   end do
!   call system_clock(e_time)
!   if(pos%rank==0) then
!      call pos%cup(0)%show()
!      write(*,*)real(e_time-s_time)/c_rate
!   end if
!   call pos%end_parallel()
!   end program test
   
