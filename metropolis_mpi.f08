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
   use parallel_universe, only: pp2d_mpi, gridmap
   use seed_md,           only: seed
   use histogram,         only: hist1d
   use iso_fortran_env,   only: int64
   implicit none
   type(pp2d_mpi)      :: pos
   type(gridmap)       :: map
   type(stack)         :: env
   type(seed)          :: sd
   integer(int64)      :: s_time, e_time, c_rate, timestamp
   integer             :: g(2), shift(2), idx, dir, n,  cycle_reward, &
                          step, steps, cyc, all_cycles, uscalars, uvectors, &
                          tdamp, dump_every, cyc_dump, dump_total, ntry, nsuccess, &
                          u8, u16, u32
   real(pr)            :: energy, virial, delta, step_time, de, psi(2), &
                          rho, tem, rc, rc2, ecut, dmax, rnd, alpha
   character(len=20)   :: arg
   type(hist1d)        :: hst8, hst16, hst32
   ! build system 
   call sd%make(pos%pp2d,timestamp)

   ! parameters of simulation
   rc    = sd%rc
   rc2   = rc*rc
   ecut  = 1.0_pr/rc2**6 - 1.0_pr/rc2**3
   dmax  = sd%dmax
   tem   = sd%tem
   alpha = sd%alpha
   rho   = sd%rho
   tdamp = sd%tdamp
   dump_every = 10**6
   dump_total = 100
   call get_command_argument(1,arg)
   if( len_trim(arg)>0 ) read(arg,*) dump_total

   ! mpi setup
   call pos%start_parallel()
   call pos%bcast_from_master()    ! maybe redundant 
   call pos%unique_rnd()
   if( pos%rank==0 ) then
      call sd%open(uscalars,uvectors)
      call fopen("mc_dprof_8.txt",u8)
      call fopen("mc_dprof_16.txt",u16)
      call fopen("mc_dprof_32.txt",u32)
      delta = 1.0_pr/(pos%wx*pos%wy*256)
      call hst8%init(rho-0.1_pr,rho+0.1_pr,delta)
      call hst16%init(rho-0.1_pr,rho+0.1_pr,delta)
      call hst32%init(rho-0.1_pr,rho+0.1_pr,delta)
   end if

   ! scheduling
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
   call pos%suggest_mapping(g)
   call pos%create_mapping(g,map)
   call system_clock(s_time,c_rate)
   do cyc = 1, all_cycles

      ! main
      shift = pos%randcc() 
      call pos%do_mapping(map,shift)
      call pos%stage(pos%c_mine, pos%w_mine-1)
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
         de = 4*sum(env%g(1:n)-env%f(1:n))
         if( de<=0.0_pr ) then
            call pos%move(idx,dir,delta) ; nsuccess = nsuccess + 1
         else
            call random_number(rnd)
            if(rnd<=exp(-de/tem)) then
               call pos%move(idx,dir,delta) ; nsuccess = nsuccess + 1
            end if
         end if
      end do

      ! timestamp
      call pos%update_all()
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
         energy = 2*energy/(pos%lnop)
         virial = rho*tem + virial/(pos%lx*pos%ly)
         psi = psi/pos%nop
         write(uscalars,*) timestamp, energy, virial, psi
         call hst8%gather(pos%dprof(8))
         call hst16%gather(pos%dprof(16))
         call hst32%gather(pos%dprof(32))
         call sd%dump(pos,timestamp)
         if( mod(cyc,cyc_dump)==0 ) then
                 call pos%write(uvectors,string="new",ints=[timestamp])
                 call hst8%write(u8,ints=[timestamp])
                 call hst16%write(u16,ints=[timestamp])
                 call hst32%write(u32,ints=[timestamp])
                 call hst8%reset() 
                 call hst16%reset() 
                 call hst32%reset() 
         end if
      end if

      ! speed calculator
      if( pos%rank==0 ) then 
         call system_clock(e_time)
         step_time = real(10**9*dble(e_time-s_time)/(c_rate*cycle_reward*pos%nop))
         open(newunit=n,file="mc_speed_notes.txt")
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
      close(u8)
      close(u16)
      close(u32)
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
         r = r*r*r
         en = (1.0_pr/r-alpha)/r - ecut
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
         en = 6*(2.0_pr/r-alpha)/r 
      end if
      end function

      subroutine fopen(fname,u)
      implicit none
      character(len=*), intent(in) :: fname
      integer, intent(out)         :: u
      logical                      :: yes
      inquire(file="mc_scalars.txt",exist=yes)
      if( yes ) then
         open(newunit=u,file=fname,status="old",action="write",access="append")
      else
         open(newunit=u,file=fname,status="new",action="write")
      end if
      end subroutine

   end program
