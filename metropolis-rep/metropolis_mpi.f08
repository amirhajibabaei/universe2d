! note: 
!      increasing the update frequency by a factor of 10
!      increases the overall cost by a factor of 2.
!
!       internal pos%metro is slower by a factor of 2
!
! scheduling:
!      a sweep is N trials
!      a cycle  is ncpu*tdamp sweeps + an update
!      timestamp is in sweeps
!      a cycle_reward is ncpu*tdamp timestamps
!
   program parallel_metropolis
   use universe,          only: stack, pr
   use parallel_universe, only: pp2d_mpi
   use seed_md,           only: seed
   use iso_fortran_env,   only: int64
   implicit none
   type(pp2d_mpi)      :: pos
   type(stack)         :: env
   type(seed)          :: sd
   integer(int64)      :: s_time, e_time, c_rate, timestamp
   integer             :: shift(2), idx, dir, n,  cycle_reward, &
                          step, steps, cyc, all_cycles, uscalars, uvectors, &
                          tdamp, dump_every, cyc_dump, dump_total, ntry, nsuccess, &
                          pow, powh
   real(pr)            :: energy, virial, delta, step_time, de, psi(2), &
                          rho, rc, rc2, ecut, dmax, rnd
   character(len=20)   :: arg

   ! build system 
   call sd%make(pos%pp2d,timestamp)

   ! parameters of simulation
   pow   = sd%pow; powh = pow/2
   if( mod(pow,2)/=0 ) stop
   rc    = sd%rc; rc2   = rc*rc
   ecut  = 1.0_pr/rc2**powh
   dmax  = sd%dmax
   rho   = sd%rho
   tdamp = sd%tdamp

   ! mpi setup
   call pos%start_parallel()
   call pos%unique_rnd()
   if( pos%rank==0 ) then
      call sd%open(uscalars,uvectors)
   end if

   ! scheduling
   dump_every = 10**6
   dump_total = 40
   call get_command_argument(1,arg)
   if( len_trim(arg)>0 ) read(arg,*) dump_total
   if( pos%rank==0 ) then
      steps = (tdamp + 1 - pos%size)*pos%nop 
   else
      steps = (tdamp + 1 )*pos%nop 
   end if
   cycle_reward = tdamp*pos%size
   cyc_dump     = max(dump_every/cycle_reward,1)
   all_cycles   = dump_total*cyc_dump

   ! run
   ntry = 0
   nsuccess = 0
   call system_clock(s_time,c_rate)
   do cyc = 1, all_cycles

      ! main
      shift = pos%randcc() 
      call pos%shift_mapping(shift)
      do step = 1, steps
         call pos%random(-dmax,dmax,idx,dir,delta)
         call pos%zoom_on(idx,env)
         n = env%n
         env%f(1:n) = efunc(env%x(1:n),env%y(1:n))
         if( dir==1 ) then
            env%g(1:n) = efunc(env%x(1:n)-delta,env%y(1:n))
         else
            env%g(1:n) = efunc(env%x(1:n),env%y(1:n)-delta)
         end if
         de = sum(env%g(1:n)-env%f(1:n))
         if( de<=0.0_pr ) then
            call pos%move(idx,dir,delta) ; nsuccess = nsuccess + 1
         else
            call random_number(rnd)
            if(rnd<=exp(-de)) then
               call pos%move(idx,dir,delta) ; nsuccess = nsuccess + 1
            end if
         end if
      end do

      ! timestamp
      call pos%update_master()
      timestamp = timestamp + cycle_reward

      ! tune dmax but keep it less than sd%dmax 
      ntry = ntry + steps
      if( ntry>1000 ) then
         rnd = real(nsuccess)/ntry
         if( rnd>0.5 ) then
              dmax = 1.05*dmax
         else
              dmax = 0.95*dmax
         end if
         if( dmax>sd%dmax ) dmax = sd%dmax
         ntry = 0
         nsuccess = 0
      end if

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
            psi = psi + env%hexorder()
         end do
         energy = 0.5_pr*energy/(pos%lnop)
         virial = rho + 0.25_pr*virial/(pos%lx*pos%ly)
         psi = psi/pos%nop
         write(uscalars,*) timestamp, energy, virial, psi
         call sd%dump(pos,timestamp)
         if( mod(cyc,cyc_dump)==0 ) then
                 call pos%write(uvectors,string="new",ints=[timestamp])
         end if
      end if

      ! mc runtime notes
      if( pos%rank==0 ) then 
         call system_clock(e_time)
         step_time = real(10**9*dble(e_time-s_time)/(c_rate*cycle_reward*pos%nop))
         open(newunit=n,file="mc_notes.txt")
            write(n,*) "timestamp", timestamp, "nxy", pos%nxy
            write(n,*)
            write(n,*) "cost of step: ", pos%size*step_time, " / "
            write(n,*) "num of procs: ", pos%size
            write(n,*) "             =", step_time
            write(n,*) 
            write(n,*) "dmax         =", dmax, sd%dmax
         close(n)
         call system_clock(s_time,c_rate)
      end if

   end do

   ! end mpi
   if( pos%rank==0 ) then 
      close(uscalars)
      close(uvectors)
   end if
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
         en = 1.0_pr/r**powh - ecut
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
         en = pow/r**powh
      end if
      end function

!      subroutine fopen(fname,u)
!      implicit none
!      character(len=*), intent(in) :: fname
!      integer, intent(out)         :: u
!      logical                      :: yes
!      inquire(file="mc_scalars.txt",exist=yes)
!      if( yes ) then
!         open(newunit=u,file=fname,status="old",action="write",access="append")
!      else
!         open(newunit=u,file=fname,status="new",action="write")
!      end if
!      end subroutine

   end program
