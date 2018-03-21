   module parallel_universe
   use mpi_f08
   use universe
   implicit none
   private
   public pp2d_mpi, pr, particle, pp2d, gridmap
   
      type gridmap
         integer              :: qx, qy, n
         integer              :: q(2), nxy(2)
         integer, allocatable :: wx(:), wy(:)
         integer, allocatable :: c1(:,:), c2(:,:), cw(:,:)
      end type gridmap
   
      type, extends(pp2d)     :: pp2d_mpi
         integer              :: size
         integer              :: rank
         integer, allocatable :: owner(:)
         integer, allocatable :: old_owner(:)
         integer              :: c_block(2), w_block(2) 
      contains
         procedure            :: start_parallel
         procedure            :: end_parallel
         procedure            :: suggest_mapping
         procedure            :: create_mapping
         procedure            :: do_mapping
         procedure            :: update 
      end type pp2d_mpi
   
   contains
   
      subroutine start_parallel(pos)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos
      integer                        :: ierr
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,pos%rank,ierr)
      call mpi_comm_size(mpi_comm_world,pos%size,ierr)
      allocate(pos%owner(0:pos%noc-1))
      allocate(pos%old_owner(0:pos%noc-1))
      end subroutine start_parallel
      
      subroutine end_parallel(pos)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos
      integer                        :: ierr
      call mpi_finalize(ierr)
      deallocate(pos%owner)
      end subroutine end_parallel
   
      subroutine suggest_mapping(pos,best)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos 
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
      class(pp2d_mpi), intent(inout) :: pos 
      integer,         intent(in)    :: q(2)
      type(gridmap),  intent(out)    :: map
      integer                        :: w(2), r(2) 
      integer                        :: nx(q(1)), ny(q(2))
      integer                        :: i,j, n, cx1, cx2, cy1, cy2
      map%q = q; map%qx = q(1); map%qy = q(2)
      map%n = q(1)*q(2)
      map%nxy = pos%nxy
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
      integer                        :: cx0, cy0, cx, cy
      integer                        :: i, k, l, cl
      if( present(shift) ) then
         cx0 = shift(1); cy0 = shift(2)
      else
         cx0 = 0; cy0 = 0
      end if
      do i = 0, map%n-1
         do k = map%c1(2,i), map%c2(2,i)
            cy = mod(k + cy0, pos%ny)
            if( cy<0 ) cy = cy + pos%ny
            do l = map%c1(1,i), map%c2(1,i)
               cx = mod(l + cx0, pos%nx)
               if( cx<0 ) cx = cx + pos%nx
               cl = cy*pos%nx + cx 
               pos%owner(cl) = i
            end do
         end do
      end do
      pos%c_block = mod( map%c1(:,pos%rank) + shift, pos%nxy )
      where( pos%c_block<0 ) pos%c_block = pos%c_block + pos%nxy
      pos%w_block = map%cw(:,pos%rank)
      end subroutine do_mapping
   
      subroutine update(pos)
      implicit none
      class(pp2d_mpi), intent(inout) :: pos 
      integer                        :: i, ierr, occ
      do i = 0, pos%noc-1
         occ = pos%cup(i)%occ
         call mpi_bcast(occ, 1, mpi_integer, pos%owner(i), mpi_comm_world, ierr)
         do
            if( occ<=pos%cup(i)%cap ) exit
            call pos%cup(i)%ext()
         end do
         call mpi_bcast(pos%cup(i)%hld(1:occ)%name, occ, mpi_integer, pos%owner(i), mpi_comm_world, ierr)
         call mpi_bcast(pos%cup(i)%hld(1:occ)%pos(1), occ, mpi_real, pos%owner(i), mpi_comm_world, ierr)
         call mpi_bcast(pos%cup(i)%hld(1:occ)%pos(2), occ, mpi_real, pos%owner(i), mpi_comm_world, ierr)
         pos%cup(i)%occ = occ
      end do
      end subroutine update 
   
   end module parallel_universe



   program parallel_metropolis
   use parallel_universe
   use iso_fortran_env, only: int64
   implicit none
   type(pp2d_mpi)      :: pos
   type(particle)      :: part
   type(gridmap)       :: map
   integer             :: i, j, g(2), shift(2), idx, dir 
   real(pr)            :: delta, step_time
   real(pr), parameter :: dmax = 0.1
   integer, parameter  :: steps = 10**6, cycles = 10
   integer(int64)      :: s_time, e_time, c_rate
   ! build system 
   pos%pp2d = pp2d([500.0_pr,500.0_pr],3.0_pr)
   call pos%relative()
   call pos%assoc()
   do i = 1, 500*500
      part%name = i
      call random_number(part%pos)
      part%pos = part%pos * pos%ll
      call pos%add(i,part)
   end do
   call pos%reserve()
   ! write initial
   if( pos%rank==0 ) then
     open(1,file="init.txt")
        call pos%write(1,0)
     close(1)
   end if
   ! mpi
   call pos%start_parallel()
   call pos%suggest_mapping(g)
   call pos%create_mapping(g,map)
   call system_clock(s_time,c_rate)
   shift = [0,0]
   do j = 1, cycles
      shift = shift+1
      call pos%do_mapping(map,shift)
      call pos%activate(pos%c_block, pos%w_block-1)
      do i = 1, steps
         call pos%random(dmax,idx,dir,delta)
         call pos%move(idx,dir,delta)
      end do
      call pos%update()
   end do
   call system_clock(e_time)
   step_time = real(10**9*dble(e_time-s_time)/(c_rate*steps*cycles))
   write(*,*) step_time, step_time/pos%size 
   ! write last
   if( pos%rank==0 ) then
     open(1,file="fin.txt")
        call pos%write(1,-1)
     close(1)
   end if
   ! end mpi
   call pos%end_parallel()
   end program
