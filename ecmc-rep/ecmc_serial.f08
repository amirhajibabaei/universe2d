program ecmc
use universe, only: pp2d, pr
implicit none
type(pp2d) :: pos
pos = pp2d([64,64],1.0_pr,2.0_pr)
print*, pos%nop
end program
