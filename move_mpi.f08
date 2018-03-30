! note: 
!      1. increasing the update frequency by a factor of 10
!         increases the overall cost by a factor of 2.
!
!      2. timestamp*nop = total mc trials
!
   program parallel_metropolis
   use universe, only: stack
   use parallel_universe
   use seed_md, only: seed
   use iso_fortran_env, only: int64
   implicit none
   real(pr), parameter :: pi = 3.14159265359_pr
   type(pp2d_mpi)      :: pos
   type(gridmap)       :: map
   type(stack)         :: env
   type(seed)          :: sd
   integer(int64)      :: s_time, e_time, c_rate, timestamp
   integer             :: g(2), shift(2), idx, dir, n,  cycle_reward, &
                           step, steps, cyc, cycles, uscalars, nrush
   real(pr)            :: energy, virial, delta, step_time, de, psi(2), &
                           rho, tem, rc, rc2, ecut, dmax, rnd, a0, rn2
   ! build system 
   call sd%make(pos%pp2d,timestamp)
   rc   = sd%rc
   rc2  = rc*rc
   ecut = 1.0_pr/rc2**6 - 1.0_pr/rc2**3
   dmax = sd%dmax
   tem  = sd%tem
   rho  = sd%rho
   a0   = sqrt(2.0_pr/(rho*sqrt(3.0_pr)))
   rn2  = (1.5_pr*a0)**2

   ! mpi setup
   call pos%start_parallel()
   call pos%unique_rnd()
   if( pos%rank==0 ) call sd%open(uscalars)

   ! scheduling
   nrush = 100
   if( pos%rank==0 ) then
      steps = (nrush + 1 - pos%size)*pos%nop 
   else
      steps = (nrush + 1 )*pos%nop 
   end if
   cycle_reward = nrush*pos%size
   cycles       = 10**1

   ! run
   call pos%suggest_mapping(g)
   call pos%create_mapping(g,map)
   call system_clock(s_time,c_rate)
   do cyc = 1, cycles

      ! main
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

      ! timestamp
      call pos%update_all()
      timestamp = timestamp + cycle_reward 

      ! calculate
      if( pos%rank==0) then      
         call pos%stage()
         energy = 0.0_pr
         virial = 0.0_pr
         psi = 0.0_pr
         do idx = 0, pos%nop-1
            call pos%zoom_on(idx,env)
            n = env%n
            env%f(1:n) = efunc(env%x(1:n),env%y(1:n))
            env%g(1:n) = vfunc(env%x(1:n),env%y(1:n))
            energy = energy + sum(env%f(1:n))
            virial = virial + sum(env%g(1:n))         
            call order(env%x(1:n),env%y(1:n),env%f(1:n),env%g(1:n),env%h(1:n))
            rnd = sum(env%f(1:n))
            if(rnd>=1.0_pr) psi = psi + [sum(env%g(1:n)),sum(env%h(1:n))]/rnd 
         end do
         energy = 2*energy/(pos%lnop)
         virial = rho*tem + virial/(pos%lx*pos%ly)
         psi = psi/pos%nop
         write(uscalars,*) timestamp, energy, virial, psi
         call sd%dump(pos,timestamp)
      end if

   end do
   call system_clock(e_time)

   ! end mpi
   if( pos%rank==0 ) close(uscalars)
   step_time = real(10**9*dble(e_time-s_time)/(c_rate*cycles*cycle_reward*pos%nop))
   write(*,*) pos%size, pos%rank, step_time
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

      elemental &
      subroutine order(x,y,n,or,oi) 
      implicit none
      real(pr), intent(in)  :: x, y
      real(pr), intent(out) :: n, or, oi
      real(pr)              :: r, theta
      r = x*x+y*y
      if( r>=rn2 ) then
         n = 0.0_pr
         or = 0.0_pr
         oi = 0.0
      else
         n = 1.0_pr
         r = sqrt(r)
         theta = acos(x/r) 
         if( y<0.0_pr ) theta = 2*pi - theta
         or = cos(6*theta)
         oi = sin(6*theta)
      end if
      end subroutine

   end program
