   module parallel_universe
   use mpi_f08
   use universe
   use iso_fortran_env, only: int64, real32, real64
   implicit none
   private
   public pp2d_mpi, pr, particle, pp2d, gridmap
   
      type gridmap
         integer              :: qx, qy, n
         integer              :: q(2)
         integer, allocatable :: wx(:), wy(:)
         integer, allocatable :: c1(:,:), c2(:,:), cw(:,:)
      end type gridmap
   
      type, extends(pp2d)     :: pp2d_mpi
         integer              :: size
         integer              :: rank
         integer, allocatable :: owner(:)
         integer              :: c_mine(2), w_mine(2) 
      contains
         procedure, private   :: mod2, cl_on
         procedure            :: randcc, unique_rnd
         procedure            :: start_parallel
         procedure            :: end_parallel
         procedure            :: suggest_mapping
         procedure            :: create_mapping
         procedure            :: do_mapping
         procedure            :: push 
         procedure            :: update_master 
         procedure            :: update_all 
         procedure            :: bcast_from_master 
         procedure            :: metro_mpi => metropolis 
      end type pp2d_mpi
   
   contains
   
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

      function randcc(pos) result(cc)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos 
      integer                        :: cc(2), ierr
      real(pr)                       :: rnd(2)
      call random_number(rnd)
      cc = floor(pos%nxy*rnd)
      call mpi_bcast( cc, 2, mpi_integer, pos%size-1, mpi_comm_world, ierr)
      end function randcc

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine start_parallel(pos)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos
      integer                        :: ierr
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,pos%rank,ierr)
      call mpi_comm_size(mpi_comm_world,pos%size,ierr)
      allocate(pos%owner(0:pos%noc-1))
      end subroutine start_parallel

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
      
      subroutine end_parallel(pos)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos
      integer                        :: ierr
      call mpi_finalize(ierr)
      deallocate(pos%owner)
      end subroutine end_parallel
   
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
   
      subroutine create_mapping(pos,q,map)
      implicit none
      class(pp2d_mpi), intent(in) :: pos 
      integer,         intent(in)    :: q(2)
      type(gridmap),  intent(out)    :: map
      integer                        :: w(2), r(2) 
      integer                        :: nx(q(1)), ny(q(2))
      integer                        :: i,j, n, cx1, cx2, cy1, cy2
      map%q = q; map%qx = q(1); map%qy = q(2)
      map%n = q(1)*q(2)
      w = pos%nxy/q
      r = pos%nxy - q*w
      nx = w(1); nx(1:r(1)) = nx(1:r(1)) + 1
      ny = w(2); ny(1:r(2)) = ny(1:r(2)) + 1
      map%wx = nx
      map%wy = ny
      allocate(map%c1(2,0:map%n-1))
      allocate(map%c2(2,0:map%n-1))
      allocate(map%cw(2,0:map%n-1))
      n = 0
      cy1 = 0
      do j = 1, q(2)
         cy2 = cy1 + ny(j)
         cx1 = 0
         do i = 1, q(1)
            cx2 = cx1 + nx(i)
            map%c1(:,n) = [cx1,cy1]
            map%c2(:,n) = [cx2,cy2] - 1
            map%cw(:,n) = [nx(i),ny(j)]
            n = n+1
            cx1 = cx2
         end do
         cy1 = cy2
      end do
      end subroutine create_mapping
   
      subroutine do_mapping(pos,map,shift)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos 
      type(gridmap),    intent(in)   :: map
      integer, intent(in), optional  :: shift(2)
      integer                        :: c0(2), cx, cy
      integer                        :: i, cl
      if( present(shift) ) then
         c0 = shift
      else
         c0 = 0
      end if
      do i = 0, map%n-1
         do cy = map%c1(2,i), map%c2(2,i)
            do cx = map%c1(1,i), map%c2(1,i)
               cl = pos%cl_on(c0,[cx,cy])
               pos%owner(cl) = i
            end do
         end do
      end do
      where(pos%owner>=pos%size) pos%owner=-1
      pos%c_mine = pos%mod2( map%c1(:,pos%rank), c0 )
      pos%w_mine = map%cw(:,pos%rank)
      end subroutine do_mapping

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine sendrecv(pos,i,rk1,rk2)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos 
      integer,         intent(in)    :: i, rk1, rk2
      integer                        :: occ, tag, ierr
      if( rk1 == rk2 .or. rk1==-1 .or. rk2==-1 ) return
      if( pos%rank == rk1 ) then
         tag = i
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
      end subroutine sendrecv

      subroutine push_cl(pos, cl)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos 
      integer,         intent(in)    :: cl
      integer                        :: ccl, sender, reciever, k
      sender = pos%owner(cl)
      if( sender==-1 ) return
      do k = 1, 8
         ccl = pos%cup(cl)%su(k)
         reciever = pos%owner(ccl)
         if( reciever==-1 .or. reciever == sender ) then
                 cycle
         else
                 call sendrecv(pos,cl,sender,reciever)
         end if
      end do
      end subroutine push_cl

      subroutine push(pos)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos
      integer                        :: cl
      do cl = 0, pos%noc-1
         call push_cl(pos,cl)
      end do
      end subroutine push

      subroutine update_master(pos)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos
      integer                        :: cl, sender, reciever
      reciever = 0
      do cl = 0, pos%noc-1
         sender = pos%owner(cl)
         call sendrecv(pos,cl,sender,reciever)
      end do
      end subroutine update_master

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine update_all(pos)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos 
      integer                        :: i, ierr, occ, owner
      if( pr==4 ) then
         do i = 0, pos%noc-1
            owner = pos%owner(i)
            if( owner==-1 ) cycle
            occ = pos%cup(i)%occ
            call mpi_bcast(occ, 1, mpi_integer, owner, mpi_comm_world, ierr)
            do
               if( occ<=pos%cup(i)%cap ) exit
               call pos%cup(i)%ext()
            end do
            call mpi_bcast( pos%cup(i)%hld(1:occ)%name  , occ, mpi_integer, owner, mpi_comm_world, ierr)
            call mpi_bcast( pos%cup(i)%hld(1:occ)%pos(1), occ, mpi_real,    owner, mpi_comm_world, ierr)
            call mpi_bcast( pos%cup(i)%hld(1:occ)%pos(2), occ, mpi_real,    owner, mpi_comm_world, ierr)
            pos%cup(i)%occ = occ
         end do
      elseif( pr==8 ) then
         do i = 0, pos%noc-1
            owner = pos%owner(i)
            if( owner==-1 ) cycle
            occ = pos%cup(i)%occ
            call mpi_bcast(occ, 1, mpi_integer, owner, mpi_comm_world, ierr)
            do
               if( occ<=pos%cup(i)%cap ) exit
               call pos%cup(i)%ext()
            end do
            call mpi_bcast( pos%cup(i)%hld(1:occ)%name  , occ, mpi_integer, owner, mpi_comm_world, ierr)
            call mpi_bcast( pos%cup(i)%hld(1:occ)%pos(1), occ, mpi_double_precision, owner, mpi_comm_world, ierr)
            call mpi_bcast( pos%cup(i)%hld(1:occ)%pos(2), occ, mpi_double_precision, owner, mpi_comm_world, ierr)
            pos%cup(i)%occ = occ
         end do
      else
         write(*,*) "wrong kind!"
         stop
      end if
      end subroutine update_all

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine bcast_from_master(pos)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos 
      integer                        :: i, ierr, occ, owner
      owner = 0
      if( pr==4 ) then
         do i = 0, pos%noc-1
            occ = pos%cup(i)%occ
            call mpi_bcast(occ, 1, mpi_integer, owner, mpi_comm_world, ierr)
            do
               if( occ<=pos%cup(i)%cap ) exit
               call pos%cup(i)%ext()
            end do
            call mpi_bcast( pos%cup(i)%hld(1:occ)%name  , occ, mpi_integer, owner, mpi_comm_world, ierr)
            call mpi_bcast( pos%cup(i)%hld(1:occ)%pos(1), occ, mpi_real,    owner, mpi_comm_world, ierr)
            call mpi_bcast( pos%cup(i)%hld(1:occ)%pos(2), occ, mpi_real,    owner, mpi_comm_world, ierr)
            pos%cup(i)%occ = occ
         end do
      elseif( pr==8 ) then
         do i = 0, pos%noc-1
            occ = pos%cup(i)%occ
            call mpi_bcast(occ, 1, mpi_integer, owner, mpi_comm_world, ierr)
            do
               if( occ<=pos%cup(i)%cap ) exit
               call pos%cup(i)%ext()
            end do
            call mpi_bcast( pos%cup(i)%hld(1:occ)%name  , occ, mpi_integer, owner, mpi_comm_world, ierr)
            call mpi_bcast( pos%cup(i)%hld(1:occ)%pos(1), occ, mpi_double_precision, owner, mpi_comm_world, ierr)
            call mpi_bcast( pos%cup(i)%hld(1:occ)%pos(2), occ, mpi_double_precision, owner, mpi_comm_world, ierr)
            pos%cup(i)%occ = occ
         end do
      else
         write(*,*) "wrong kind!"
         stop
      end if
      end subroutine bcast_from_master

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine metropolis(pos,efunc,eargs,tem,dmax,steps,nsuccess)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos 
      interface 
         function efunc(x,y,args) result(en)
         import pr 
         implicit none
         real(pr), intent(in) :: x(:), y(:), args(:)
         real(pr)             :: en(size(x))
         end function
      end interface
      real(pr),    intent(in)    :: eargs(:), tem, dmax
      integer,     intent(in)    :: steps
      integer,     intent(out)   :: nsuccess
      type(gridmap)              :: map
      integer                    :: g(2), shift(2)
      call pos%suggest_mapping(g)
      call pos%create_mapping(g,map)
      shift = pos%randcc() 
      call pos%do_mapping(map,shift)
      call pos%stage(pos%c_mine, pos%w_mine-1)
      call pos%metro(efunc,eargs,tem,dmax,steps,nsuccess)
      call pos%update_all()
      end subroutine metropolis

   end module parallel_universe

