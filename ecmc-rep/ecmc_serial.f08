!
!   at rho=1.0 => P(delta) = exp(-delta/0.052)
!
!   ~ 5 metro sweeps per sec
!


    program ecmc
    use iso_fortran_env, only: int64
    use universe, only: pp2d, stack, pr, pi
    use histogram, only: hist1d, hist2d
    implicit none
    type(pp2d)     :: pos
    type(stack)    :: env
    type(hist1d)   :: h_delta, h_xji, h_rji
    type(hist2d)   :: h_gxy
    real(pr)       :: rnd, cutoff, skin, cut2, wth, rho, &
                      powhi, delta, rji(2), max_del, sum_del, &
                      xji, dji, dji_ave, pressure, psi(2), &
                      theta, mat(2,2)
    integer        :: idx, dir, pow, powh, mloc(1), n, &
                       ch, sweep, num_sweeps, uscalars, &
                       uvecs, sample
    integer(int64) :: timestamp
!    character(len=*) , parameter :: address = "../../repulsive/p_64/c_1.8/r_0.882/n_256/mc_restart.txt"


    pow = 12; powh = pow/2; powhi = 1.0_pr/powh
    cutoff = 1.8_pr; cut2 = cutoff**2
    skin = 0.2_pr
    wth = cutoff + skin
    call h_delta%init(0.0_pr,0.6_pr,0.005_pr)
    call h_xji%init(0.0_pr,wth,0.01_pr)
    call h_rji%init(0.6_pr,wth,0.01_pr)
    call h_gxy%init(-wth,wth,0.05_pr,-wth,wth,0.05)
    call fopen("ecmc_scalars.txt",uscalars) 
    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


    open(newunit=uvecs,file="mc_vectors.txt",action="read")
    do sample = 1, 40
    
        pos = pp2d(uvecs,wth); rho = pos%nop/(pos%lx*pos%ly)
        call pos%stage()
    
        max_del = 0.05*sqrt(real(pos%nop))/2
        num_sweeps = 10 
        timestamp = 0
        do sweep = 1, num_sweeps
           !~~~~~~~~~~~ a sweep
           dji_ave = 0.0_pr
           do ch = 1, pos%nop
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ an event chain
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
                    timestamp = timestamp + 1
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
                    delta = env%f(mloc(1))                    ; call h_delta%gather(delta)
                    rji   = [ env%x(mloc(1)), env%y(mloc(1)) ]; call h_rji%gather(sqrt(sum(rji**2)))
                    xji   = rji(dir) - delta                  ; call h_xji%gather(xji)
                    ! move
                    sum_del = sum_del + delta
                    if( sum_del>=max_del ) then
                            delta = delta - (sum_del - max_del)
                            call pos%move(idx,dir,delta)
                            exit
                    else
                            dji = dji + xji 
                            call pos%move(idx,dir,delta)
                    end if
                    ! transfer idx 
                    idx   = env%id(mloc(1))
                end do!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ end of an event chain
                dji_ave = dji_ave + dji
           end do !~~~~~~~~~~~~~~~~~  end of a sweep
           ! pressure
           dji_ave = dji_ave/(pos%nop)
           pressure = rho*(1.0_pr + dji_ave/max_del)
           ! orientational order
           psi = 0.0_pr
           do idx = 0, pos%nop-1
              call pos%zoom_on(idx,env)
              psi = psi + env%hexorder()
           end do
           psi = psi/pos%nop
           ! angle of orientation 
           rnd = sum(psi**2)
           theta = acos( psi(1)/sqrt(rnd) )
           if( psi(2)<0.0_pr ) theta = 2*pi-theta
           theta = theta/6
           ! gxy
           mat = reshape( [cos(-theta), sin(-theta), &
                          -sin(-theta), cos(-theta)], [2,2] )
           do idx = 0, pos%nop-1
              call pos%zoom_on(idx,env)
              call env%rotate(matrix=mat)
              n = env%n
              call h_gxy%gather(env%x(1:n),env%y(1:n))
           end do
           !
           write(uscalars,*) sample, pressure, psi, theta
           write(0,*) sample, pressure, psi, theta
        end do
    
    end do

    !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    close(uscalars)
    call h_delta%write("hist_delta.txt")
    call h_xji%write("hist_xji.txt")
    call h_rji%write("hist_rji.txt")
    call h_gxy%write("hist_gxy.txt")


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
   

            subroutine fopen(fname,u)
            implicit none
            character(len=*), intent(in) :: fname
            integer, intent(out)         :: u
            logical                      :: yes
            inquire(file=fname,exist=yes)
            if( yes ) then
               open(newunit=u,file=fname,status="old",action="write",access="append")
            else
               open(newunit=u,file=fname,status="new",action="write")
            end if
            end subroutine
      
    end program

