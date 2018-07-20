!
!   at rho=1.0 => P(delta) = exp(-delta/0.052)
!
!
!


    program ecmc
    use iso_fortran_env, only: int64
    use universe, only: pp2d, stack, pr
    use histogram, only: hist1d
    implicit none
    type(pp2d)     :: pos
    type(stack)    :: env
    type(hist1d)   :: h_delta
    real(pr)       :: rnd, cutoff, skin, cut2, wth, rho, &
                      powhi, delta, rj(2), max_del, sum_del, &
                      dji, dji_ave, pressure
    integer        :: idx, dir, nx, pow, powh, mloc(1), n, &
                       ch, num_chains 
    integer(int64) :: s_time, e_time, c_rate
    character(len=*) , parameter :: address = "mc_restart.txt" !"../../repulsive/p_12/c_1.8/r_1.000/n_256/mc_restart.txt"
    
    rho = 1.0_pr
    nx  = 64
    pow = 12; powh = pow/2; powhi = 1.0_pr/powh
    cutoff = 1.8_pr; cut2 = cutoff**2
    skin = 0.2_pr
    wth = cutoff + skin
    !pos = pp2d([nx,nx],rho,wth,0.01)
    pos = pp2d(address,wth)
    call pos%stage()
    call h_delta%init(0.0_pr,0.6_pr,0.004_pr)


    max_del = 0.2*nx/2


    call system_clock(s_time,c_rate)
    dji_ave = 0.0_pr
    num_chains = 10000
    do ch = 1, num_chains
    
         ! select a random direction
         call random_number(rnd)
         if( rnd<0.5_pr ) then
                 dir = 1
         else
                 dir = 2
         end if
         
         ! select a random particle
         call random_number(rnd)
         idx = floor(rnd*pos%nop)

         ! event chain    
         sum_del = 0.0_pr
         dji     = 0.0_pr
         do

             ! find the earliest collision of idx
             if( dir==1 ) then
                     call pos%zoom_r(idx,env)
                     n = env%n
                     call random_number(env%g(1:n))
                     env%g(1:n) = -log(env%g(1:n))
                     env%f(1:n) = collide(env%x(1:n),env%y(1:n),env%g(1:n))
             else
                     call pos%zoom_u(idx,env)
                     n = env%n
                     call random_number(env%g(1:n))
                     env%g(1:n) = -log(env%g(1:n))
                     env%f(1:n) = collide(env%y(1:n),env%x(1:n),env%g(1:n))
             end if
             mloc  = minloc(env%f(1:n))
             delta = env%f(mloc(1));  call h_delta%gather(delta)
             rj    = [ env%x(mloc(1)), env%y(mloc(1)) ]

             ! move
             sum_del = sum_del + delta
             if( sum_del>max_del ) then
                     delta = delta - (sum_del - max_del)
                     call pos%move(idx,dir,delta)
                     exit
             else
                     dji = dji + (rj(dir) - delta)
                     call pos%move(idx,dir,delta)
             end if

             ! transfer idx 
             idx   = env%id(mloc(1))

         end do

         dji_ave = dji_ave + dji

    end do
    call system_clock(e_time)
    write(*,*) 10**9*dble(e_time-s_time)/(c_rate*10**6), sum_del/10**6
    call h_delta%write("hst.txt")
    

    dji_ave = dji_ave/num_chains
    pressure = rho*(1.0_pr + dji_ave/max_del)
    write(*,*)"pressure = ", pressure

    contains
            
            elemental &
            function collide(x,y,de) result(delta)
            implicit none
            real(pr), intent(in) :: x, y, de
            real(pr)             :: delta, e, u, v, &
                                    y2, xx, xxx
            if( x<=0.0_pr ) then ! no collision ?
                    delta = wth; return
            end if
            y2 = y*y
            if( y2>=cut2 ) then ! no collision ?
                    delta = wth; return
            end if
            u = x*x + y2
            if( u>cut2 ) then
                   xx = sqrt( cut2-y2 )
                   delta = x - xx
                   u = cut2
            else
                    xx = x
                    delta = 0.0_pr
            end if
            e = de + 1.0_pr/u**powh
            v = (1.0_pr/e)**powhi
            if( y2>v ) then ! no collision ?
                    delta = wth; return
            end if
            xxx = sqrt(v-y2)
            delta = delta + xx - xxx
            end function collide
    
    end program

