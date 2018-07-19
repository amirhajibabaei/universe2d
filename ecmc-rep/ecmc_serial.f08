program ecmc
use universe, only: pp2d, stack, pr
implicit none
type(pp2d)  :: pos
type(stack) :: env
real(pr)    :: rnd, cutoff, skin, cut2, wth, rho, powhi
integer     :: idx, dir, nx, pow, powh

pow = 12; powh = pow/2; powhi = 1.0_pr/powh
cutoff = 1.8_pr; cut2 = cutoff**2
skin = 0.2_pr
wth = cutoff + skin
rho = 1.0_pr
nx  = 64
pos = pp2d([nx,nx],rho,wth)
call pos%stage()


call random_number(rnd)
if( rnd<0.5_pr ) then
        dir = 1
else
        dir = 2
end if
call random_number(rnd)
idx = floor(rnd*pos%nop)
call pos%zoom_on(idx,env)


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
