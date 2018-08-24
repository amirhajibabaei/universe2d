


    program ecmc
    use universe, only: pp2d, pr
    use histogram, only: hist1d
    implicit none
    integer, parameter  :: cw(2) = [1,1]*20
    real(pr), parameter :: wth = 5.0_pr
    type(pp2d)         :: pos
    type(hist1d)       :: h_rho
    integer            :: sample, k, nx, uvecs, nums
    real(pr)           :: rho, ll(2), area
    character(len=:),     &
         allocatable   :: root
    character(len=100) :: arg


    call get_command_argument(1,arg)
    root = trim(arg)
    call read_param(rho,nx)

    ! tune params with mc_restart
    pos  = pp2d(root//"mc_restart.txt",wth)
    ll   = pos%ww*cw
    area = ll(1)*ll(2)
    call pos%stage([0,0],cw)
    call h_rho%init( real(floor(0.9*pos%lnop)), real(floor(1.1*pos%lnop)), 1.0 )

    nums = num_samples()
    open(newunit=uvecs,file=root//"mc_vectors.txt",status="old",action="read")
    do sample = 1, nums 
    
        pos = pp2d(uvecs,wth)
        if( sample<10 ) cycle

        do k = 1, 1000
           call pos%stage( pos%randcup() , cw )
           call h_rho%gather( real(pos%lnop) )
        end do

        write(*,*) sample, nums, root 

    end do
    close(uvecs)

    call h_rho%write(root//"hist_rho.txt", string=" cw, ww (block = cw*ww):", &
                                     ints=cw, & 
                                     reals=pos%ww, & 
                                     rescale_x_by_fac=1.0/area)


    contains

        function num_samples() result(k)
        implicit none
        integer :: i, j, k, num
        real    :: x, y
        open( 1, file=root//'mc_vectors.txt' )
        k = 0
        do
           read(1,*,end=173)
           read(1,*,end=173)
           read(1,*,end=173)
           read(1,*,end=173) num
           read(1,*,end=173) x, y
           do i = 1, num
              read(1,*,end=173) j, x, y
           end do 
           k = k+1
        end do 
        173 continue
        close( 1 )
        end function

        subroutine read_param(rho,nx)
        implicit none
        real, intent(out)    :: rho
        integer, intent(out) :: nx
        character(len=20)    :: line(4)
        integer              :: u
        open(newunit=u,file=root//"seed.txt",status='old')
            read(u,*) line!; read(line(4),*) pow
            read(u,*) line!; read(line(4),*) rc
            read(u,*) line; read(line(4),*) rho
            read(u,*) line; read(line(4),*) nx
            read(u,*) line!; read(line(4),*) dmax
            read(u,*) line!; read(line(4),*) tdamp
        close(u)
        end subroutine

    end program

