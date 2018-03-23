
   program parallel_metropolis
   use parallel_universe
   use iso_fortran_env, only: int64
   implicit none
   type(pp2d_mpi)      :: pos
   type(gridmap)       :: map
   integer             :: i, j, g(2), shift(2), idx, dir 
   real(pr)            :: delta, step_time
   real(pr), parameter :: dmax = 0.1
   integer, parameter  :: steps = 10**6, cycles = 4
   integer(int64)      :: s_time, e_time, c_rate
   ! build system 
   pos%pp2d = pp2d([100,100],0.962_pr,3.0_pr)
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
   call pos%unique_rnd()
   do j = 1, cycles
      shift = pos%randcc() 
      call pos%do_mapping(map,shift)
      call pos%stage(pos%c_mine, pos%w_mine-1)
      do i = 1, steps
         call pos%random(dmax,idx,dir,delta)
         call pos%move(idx,dir,delta)
      end do
      call pos%update_all()
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
