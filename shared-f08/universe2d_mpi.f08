   module parallel_universe
   use mpi_f08
   use universe
   use iso_fortran_env, only: int64, real32, real64
   implicit none
   private
   public pp2d_mpi, pr, particle, pp2d
   
      type gridmap
         integer              :: qx, qy, n
         integer              :: q(2)
         integer, allocatable :: wx(:), wy(:)
         integer, allocatable :: c1(:,:), c2(:,:), cw(:,:)
         integer              :: shift(2)
         integer              :: style
      end type gridmap
   
      type, extends(pp2d)     :: pp2d_mpi
        private
         integer, public      :: size, rank
         integer, allocatable :: owner(:)
         integer              :: c_mine(2), w_mine(2)
         type(gridmap)        :: map 
      contains
        private
         procedure            :: mod2, cl_on
         procedure, public    :: randcc, unique_rnd
         procedure            :: suggest_mapping
         procedure            :: create_mapping
         procedure            :: sendrecv_pp 
         procedure, public    :: start_parallel
         procedure, public    :: end_parallel
         procedure, public    :: shift_mapping
         procedure            :: update_sides 
         procedure, public    :: update_master 
!         procedure            :: update_all 
!         procedure            :: bcast_from_master 
!         procedure            :: metro_mpi => metropolis 
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
      cc = floor(pos%nxy*rnd)
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

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine suggest_mapping(pos,best)
      implicit none
      class(pp2d_mpi), intent(in) :: pos 
      integer,         intent(out)   :: best(2)
      integer                        :: q, qx, qy
      integer                        :: nbx, nby
      integer                        :: borders, least
      best = [pos%size,pos%size]
      least = (pos%nx+pos%ny)*pos%size
      do qx = 1, pos%size
         do qy = 1, pos%size
            q = qx*qy
            if( q>= pos%size .and. q-qx<pos%size .and. &
                     q-qy<pos%size .and. &
                    qx<=pos%nx .and. qy<=pos%ny ) then
               if( qx==1 ) then
                  nbx = 0
               else
                  nbx = qx
               end if
               if( qy==1 ) then
                  nby = 0
               else
                  nby = qy
               end if
               borders = nbx*pos%ny+nby*pos%nx
               if( borders<least ) then
                  least = borders
                  best = [qx, qy]
               end if
            end if
         end do
      end do
      end subroutine suggest_mapping
   
      subroutine create_mapping(pos,q)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos 
      integer,         intent(in)    :: q(2)
      integer                        :: w(2), r(2) 
      integer                        :: nx(q(1)), ny(q(2))
      integer                        :: i,j, n, cx1, cx2, cy1, cy2
      pos%map%q = q; pos%map%qx = q(1); pos%map%qy = q(2)
      pos%map%n = q(1)*q(2)
      w = pos%nxy/q
      r = pos%nxy - q*w
      nx = w(1); nx(1:r(1)) = nx(1:r(1)) + 1
      ny = w(2); ny(1:r(2)) = ny(1:r(2)) + 1
      pos%map%wx = nx
      pos%map%wy = ny
      allocate(pos%map%c1(2,0:pos%map%n-1))
      allocate(pos%map%c2(2,0:pos%map%n-1))
      allocate(pos%map%cw(2,0:pos%map%n-1))
      n = 0
      cy1 = 0
      do j = 1, q(2)
         cy2 = cy1 + ny(j)
         cx1 = 0
         do i = 1, q(1)
            cx2 = cx1 + nx(i)
            pos%map%c1(:,n) = [cx1,cy1]
            pos%map%c2(:,n) = [cx2,cy2] - 1
            pos%map%cw(:,n) = [nx(i),ny(j)]
            n = n+1
            cx1 = cx2
         end do
         cy1 = cy2
      end do
      pos%map%shift = 0
      pos%map%style = 0
      end subroutine create_mapping

      !//////////////////////////////////////////////////
      !////////////////////////////////////////////////// comm
      !//////////////////////////////////////////////////

      subroutine sendrecv_pp(pos,i,rk1,rk2)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos 
      integer,         intent(in)    :: i, rk1, rk2
      integer                        :: occ, tag, ierr
      if( rk1 == rk2 ) return
      tag = i
      if( pos%rank == rk1 ) then
         occ = pos%cup(i)%occ
         call mpi_send( occ,                            1, mpi_integer, rk2, tag, mpi_comm_world, ierr)
         call mpi_send( pos%cup(i)%hld(1:occ)%name  , occ, mpi_integer, rk2, tag, mpi_comm_world, ierr)
         call mpi_send( pos%cup(i)%hld(1:occ)%pos(1), occ, mpi_real,    rk2, tag, mpi_comm_world, ierr)
         call mpi_send( pos%cup(i)%hld(1:occ)%pos(2), occ, mpi_real,    rk2, tag, mpi_comm_world, ierr)
      elseif( pos%rank == rk2 ) then
         call mpi_recv( occ,                            1, mpi_integer, rk1, tag, mpi_comm_world, mpi_status_ignore, ierr)
         do
            if( occ<=pos%cup(i)%cap ) exit
            call pos%cup(i)%ext()
         end do
         call mpi_recv( pos%cup(i)%hld(1:occ)%name  , occ, mpi_integer, rk1, tag, mpi_comm_world, mpi_status_ignore, ierr)
         call mpi_recv( pos%cup(i)%hld(1:occ)%pos(1), occ, mpi_real,    rk1, tag, mpi_comm_world, mpi_status_ignore, ierr)
         call mpi_recv( pos%cup(i)%hld(1:occ)%pos(2), occ, mpi_real,    rk1, tag, mpi_comm_world, mpi_status_ignore, ierr)
         pos%cup(i)%occ = occ
      end if
      end subroutine sendrecv_pp

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine shift_mapping(pos,shift)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos 
      integer,           intent(in)  :: shift(2)
      integer                        :: i, cl, c0(2), cx, cy
      integer                        :: old_owner, new_owner
      pos%map%shift = pos%mod2( pos%map%shift, shift )
      c0 = pos%map%shift
      do i = 0, pos%map%n-1
         do cy = pos%map%c1(2,i), pos%map%c2(2,i)
            do cx = pos%map%c1(1,i), pos%map%c2(1,i)
               cl = pos%cl_on(c0,[cx,cy])
               old_owner = pos%owner(cl)
               new_owner = i; if( new_owner >= pos%size ) &
               new_owner = 0 
               call pos%sendrecv_pp(cl,old_owner,new_owner)
               pos%owner(cl) = new_owner 
            end do
         end do
      end do
      pos%c_mine = pos%mod2( pos%map%c1(:,pos%rank), c0 )
      pos%w_mine = pos%map%cw(:,pos%rank)
      call pos%update_sides()
      call pos%stage(pos%c_mine,pos%w_mine-1) ! if changed this, look at push_cl
      end subroutine shift_mapping

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine start_parallel(pos,cpu_grid)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos
      integer, intent(in), optional  :: cpu_grid      
      integer                        :: ierr, g(2)
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,pos%rank,ierr)
      call mpi_comm_size(mpi_comm_world,pos%size,ierr)
      allocate(pos%owner(0:pos%noc-1))
      pos%owner = 0 
      if( present(cpu_grid) ) then
         g = cpu_grid
      else
         call pos%suggest_mapping(g)
      end if
      call pos%create_mapping(g)
      call pos%shift_mapping([0,0])
      end subroutine start_parallel

   
      subroutine end_parallel(pos)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos
      integer                        :: ierr
      call mpi_finalize(ierr)
      deallocate(pos%owner)
      end subroutine end_parallel

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine push_cl(pos, cl)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos 
      integer,         intent(in)    :: cl
      integer                        :: ccl, k, sender, reciever
      logical                        :: upd(0:pos%size-1)
      integer                        :: i, nn(3)
      upd = .false.
      sender = pos%owner(cl)
      upd(sender) = .true.
      nn = [1,2,6]
      do i = 1, 3
         k = nn(i)
         ccl = pos%cup(cl)%su(k)
         reciever = pos%owner(ccl)
         if( upd(reciever) .eqv. .false. ) then
              call pos%sendrecv_pp(cl,sender,reciever)
              upd(reciever) = .true.
         end if
      end do
      end subroutine push_cl

      subroutine update_sides(pos)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos
      integer                        :: cl
      do cl = 0, pos%noc-1
         call push_cl(pos,cl)
      end do
      end subroutine update_sides

      subroutine update_master(pos)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos
      integer                        :: cl, sender, reciever
      reciever = 0
      do cl = 0, pos%noc-1
         sender = pos%owner(cl)
         call pos%sendrecv_pp(cl,sender,reciever)
      end do
      end subroutine update_master

!      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!      subroutine update_all(pos)
!      implicit none
!      class(pp2d_mpi), intent(inout) :: pos 
!      integer                        :: i, ierr, occ, owner
!      if( pr==4 ) then
!         do i = 0, pos%noc-1
!            owner = pos%owner(i)
!            if( owner==-1 ) cycle
!            occ = pos%cup(i)%occ
!            call mpi_bcast(occ, 1, mpi_integer, owner, mpi_comm_world, ierr)
!            do
!               if( occ<=pos%cup(i)%cap ) exit
!               call pos%cup(i)%ext()
!            end do
!            call mpi_bcast( pos%cup(i)%hld(1:occ)%name  , occ, mpi_integer, owner, mpi_comm_world, ierr)
!            call mpi_bcast( pos%cup(i)%hld(1:occ)%pos(1), occ, mpi_real,    owner, mpi_comm_world, ierr)
!            call mpi_bcast( pos%cup(i)%hld(1:occ)%pos(2), occ, mpi_real,    owner, mpi_comm_world, ierr)
!            pos%cup(i)%occ = occ
!         end do
!      elseif( pr==8 ) then
!         do i = 0, pos%noc-1
!            owner = pos%owner(i)
!            if( owner==-1 ) cycle
!            occ = pos%cup(i)%occ
!            call mpi_bcast(occ, 1, mpi_integer, owner, mpi_comm_world, ierr)
!            do
!               if( occ<=pos%cup(i)%cap ) exit
!               call pos%cup(i)%ext()
!            end do
!            call mpi_bcast( pos%cup(i)%hld(1:occ)%name  , occ, mpi_integer, owner, mpi_comm_world, ierr)
!            call mpi_bcast( pos%cup(i)%hld(1:occ)%pos(1), occ, mpi_double_precision, owner, mpi_comm_world, ierr)
!            call mpi_bcast( pos%cup(i)%hld(1:occ)%pos(2), occ, mpi_double_precision, owner, mpi_comm_world, ierr)
!            pos%cup(i)%occ = occ
!         end do
!      else
!         write(*,*) "wrong kind!"
!         stop
!      end if
!      end subroutine update_all
!
!      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!      subroutine bcast_from_master(pos)
!      implicit none
!      class(pp2d_mpi), intent(inout) :: pos 
!      integer                        :: i, ierr, occ, owner
!      owner = 0
!      if( pr==4 ) then
!         do i = 0, pos%noc-1
!            occ = pos%cup(i)%occ
!            call mpi_bcast(occ, 1, mpi_integer, owner, mpi_comm_world, ierr)
!            do
!               if( occ<=pos%cup(i)%cap ) exit
!               call pos%cup(i)%ext()
!            end do
!            call mpi_bcast( pos%cup(i)%hld(1:occ)%name  , occ, mpi_integer, owner, mpi_comm_world, ierr)
!            call mpi_bcast( pos%cup(i)%hld(1:occ)%pos(1), occ, mpi_real,    owner, mpi_comm_world, ierr)
!            call mpi_bcast( pos%cup(i)%hld(1:occ)%pos(2), occ, mpi_real,    owner, mpi_comm_world, ierr)
!            pos%cup(i)%occ = occ
!         end do
!      elseif( pr==8 ) then
!         do i = 0, pos%noc-1
!            occ = pos%cup(i)%occ
!            call mpi_bcast(occ, 1, mpi_integer, owner, mpi_comm_world, ierr)
!            do
!               if( occ<=pos%cup(i)%cap ) exit
!               call pos%cup(i)%ext()
!            end do
!            call mpi_bcast( pos%cup(i)%hld(1:occ)%name  , occ, mpi_integer, owner, mpi_comm_world, ierr)
!            call mpi_bcast( pos%cup(i)%hld(1:occ)%pos(1), occ, mpi_double_precision, owner, mpi_comm_world, ierr)
!            call mpi_bcast( pos%cup(i)%hld(1:occ)%pos(2), occ, mpi_double_precision, owner, mpi_comm_world, ierr)
!            pos%cup(i)%occ = occ
!         end do
!      else
!         write(*,*) "wrong kind!"
!         stop
!      end if
!      end subroutine bcast_from_master
!
!      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!      subroutine metropolis(pos,efunc,eargs,tem,dmax,steps,nsuccess)
!      implicit none
!      class(pp2d_mpi), intent(inout) :: pos 
!      interface 
!         function efunc(x,y,args) result(en)
!         import pr 
!         implicit none
!         real(pr), intent(in) :: x(:), y(:), args(:)
!         real(pr)             :: en(size(x))
!         end function
!      end interface
!      real(pr),    intent(in)    :: eargs(:), tem, dmax
!      integer,     intent(in)    :: steps
!      integer,     intent(out)   :: nsuccess
!      type(gridmap)              :: map
!      integer                    :: g(2), shift(2)
!      call pos%suggest_mapping(g)
!      call pos%create_mapping(g,map)
!      shift = pos%randcc() 
!      call pos%do_mapping(map,shift)
!      call pos%stage(pos%c_mine, pos%w_mine-1)
!      call pos%metro(efunc,eargs,tem,dmax,steps,nsuccess)
!      call pos%update_all()
!      end subroutine metropolis

   end module parallel_universe

