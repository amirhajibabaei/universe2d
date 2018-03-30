! note: increasing the update frequency by a factor of 10
! increases the overall cost by a factor of 2.
!
   program parallel_metropolis
   use universe, only: stack
   use parallel_universe
   use seed_md, only: seed
   use iso_fortran_env, only: int64
   implicit none
   type(pp2d_mpi)      :: pos
   type(gridmap)       :: map
   type(stack)         :: env
   type(seed)          :: sd
   integer(int64)      :: s_time, e_time, c_rate, timestamp
   integer             :: g(2), shift(2), idx, dir, n,  &
                           step, steps, cyc, cycles, uscalars
   real(pr)            :: energy, virial, delta, step_time, de, &
                           rho, tem, rc, rc2, ecut, dmax, rnd
   ! build system 
   call sd%make(pos%pp2d,timestamp)
   rc   = sd%rc
   rc2  = rc*rc
   ecut = 1.0_pr/rc2**6 - 1.0_pr/rc2**3
   dmax = sd%dmax
   tem  = sd%tem
   rho  = sd%rho
   steps = 10**5
   cycles = 10**1

   ! mpi
   call pos%start_parallel()
   call pos%suggest_mapping(g)
   call pos%create_mapping(g,map)
   call pos%unique_rnd()
   if( pos%rank==0 ) call sd%open(uscalars)
   call system_clock(s_time,c_rate)
   do cyc = 1, cycles
      shift = pos%randcc() 
      call pos%do_mapping(map,shift)
      call pos%stage(pos%c_mine, pos%w_mine-1)
      do step = 1, steps
         call pos%random(dmax,idx,dir,delta)
         call pos%zoom_on(idx,env)
         n = env%n
         env%f(1:n) = efunc(env%x(1:n),env%y(1:n))
         if( dir==1 ) then
            env%g(1:n) = efunc(env%x(1:n)-delta,env%y(1:n))
         else
            env%g(1:n) = efunc(env%x(1:n),env%y(1:n)-delta)
         end if
         de = 4*sum(env%g(1:n)-env%f(1:n))
         if( de<=0.0_pr ) then
            call pos%move(idx,dir,delta)
         else
            call random_number(rnd)
            if(rnd<=exp(-de/tem)) then
               call pos%move(idx,dir,delta)
            end if
         end if
      end do
      call pos%update_all()
      timestamp = timestamp + steps*pos%size
      ! calculate
      if( pos%rank==0) then      
         call pos%stage()
         energy = 0.0_pr
         virial = 0.0_pr
         do idx = 0, pos%lnop-1
            call pos%zoom_on(idx,env)
            n = env%n
            env%f(1:n) = efunc(env%x(1:n),env%y(1:n))
            env%g(1:n) = vfunc(env%x(1:n),env%y(1:n))
            energy = energy + sum(env%f(1:n))
            virial = virial + sum(env%g(1:n))         
         end do
         energy = 2*energy/(pos%lnop)
         virial = rho*tem + virial/(pos%lx*pos%ly)
         write(uscalars,*) timestamp, energy, virial
         call sd%dump(pos,timestamp)
      end if
   end do
   call system_clock(e_time)
   if( pos%rank==0 ) close(uscalars)
   step_time = real(10**9*dble(e_time-s_time)/(c_rate*steps*cycles))
   write(*,*) step_time, step_time/pos%size 
   ! end mpi
   call pos%end_parallel()

   contains

      elemental &
      function efunc(x,y) result(en)
      implicit none
      real(pr), intent(in) :: x, y
      real(pr)             :: en, r
      r = x*x+y*y
      if( r>=rc2 ) then
         en = 0.0_pr
      else
         r = r*r*r
         en = (1.0_pr/r-1.0_pr)/r - ecut
      end if
      end function

      elemental &
      function vfunc(x,y) result(en)
      implicit none
      real(pr), intent(in) :: x, y
      real(pr)             :: en, r
      r = x*x+y*y
      if( r>=rc2 ) then
         en = 0.0_pr
      else
         r = r*r*r
         en = 6*(2.0_pr/r-1.0_pr)/r 
      end if
      end function

   end program
