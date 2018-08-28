


    program ecmc
    use universe, only: pp2d, pr
    implicit none
    real(pr), parameter :: wth = 5.0_pr
    type(pp2d)          :: pos
    integer             :: nx, sample, uvecs, nums, nskip, skip(0:1000)
    real(pr)            :: rho
    character(len=:),      &
         allocatable    :: root
    character(len=100)  :: arg

    call get_command_argument(1,arg)
    root = trim(arg)
    call read_param(rho,nx)

    ! tune params with mc_restart
    pos  = pp2d(root//"mc_restart.txt",wth)

    ! pre-process input file
    call num_samples(nums,nskip,skip)
    write(*,*) root, nums, nskip

    !
    open(newunit=uvecs,file=root//"mc_vectors.txt",status="old",action="read")
    do sample = 1, nums
        ! read next configuration
        call skip_lines( uvecs, skip(sample-1) )
        pos = pp2d(uvecs,wth)
        if( sample<8 ) cycle
        !******************************
        ! your code here
        !******************************
    end do
    close(uvecs)

    write(*,*) "finished processing "//root

    contains


        subroutine num_samples(nums,nskip,skip) 
        implicit none
        integer, intent(out) :: nums, nskip, skip(0:1000)
        integer              :: i, j, k, num, line
        real                 :: x, y
        open( 1, file=root//'mc_vectors.txt' )
        skip = 0; nskip = 0
        k = 0
        do
           read(1,*,end=173)
           read(1,*,end=173)
           read(1,*,end=173)
           read(1,*,end=173) num
           read(1,*,end=173) x, y
           line = 5
           do i = 1, num
              line = line + 1
              read(1,*,end=173,err=172) j, x, y
           end do 
           k = k+1
           cycle
           172 continue
           nskip = nskip + 1
           skip(k) = line - 1
           backspace(1)
        end do 
        173 continue
        close( 1 )
        nums = k
        end subroutine


        subroutine skip_lines(u,numl)
        integer u, numl, i
        do i = 1, numl
           read(u,*)
        end do
        end subroutine


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

