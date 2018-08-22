!
! Author: Amir H. T.
!
module patches
    implicit none
    private
    public  optimal_gridmap, gridmap

    
    type gridmap
        integer              :: n ! 0, 1, ...
        integer, allocatable :: c1(:,:), c2(:,:), cw(:,:)
        integer              :: nxt      , nyt
        integer, allocatable :: xticks(:), yticks(:)
    end type gridmap


    interface gridmap
        module procedure make_gridmap
    end interface


  contains
    
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    subroutine optimal_gridmap( rect, min_parts, ticks )
    ! solves division of a rectangle with rect(1)*rect(2)
    ! discrete grids into ticks(1)*ticks(2) (>= min_parts)
    ! regions so that the length of borders is minimal
    implicit none
    integer, intent(in)  :: rect(2), min_parts
    integer, intent(out) :: ticks(2)
    integer              :: q, qx, qy, nbx, nby,&
                            length, least 
    ticks = [ min_parts, min_parts ]
    least = sum(rect)*min_parts
    do qx = 1, min_parts
        do qy = 1, min_parts
            q = qx*qy
            if( q>=min_parts .and. q-qx<min_parts .and. &
                           q-qy<min_parts .and. &
                       qx<=rect(1) .and. qy<=rect(2) ) then
               if( qx==1 ) then
                  nbx = 0
               else
                  nbx = qx
               end if
               if( qy==1 ) then
                  nby = 0
               else
                  nby = qy
               end if
               length = nbx*rect(2)+nby*rect(1)
               if( length<least ) then
                  least = length
                  ticks = [qx, qy]
               end if
            end if
        end do
    end do
    end subroutine optimal_gridmap

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    function make_gridmap( rect, q ) result(map)
    implicit none
    integer, intent(in) :: rect(2), q(2)
    type(gridmap)       :: map
    integer             :: w(2), r(2) 
    integer             :: nx(q(1)), ny(q(2))
    integer             :: i,j, n, cx1, cx2, cy1, cy2
    map%n = q(1)*q(2)
    w = rect/q
    r = rect - q*w
    nx = w(1); nx(1:r(1)) = nx(1:r(1)) + 1
    ny = w(2); ny(1:r(2)) = ny(1:r(2)) + 1
    allocate(map%c1(2,0:map%n-1))
    allocate(map%c2(2,0:map%n-1))
    allocate(map%cw(2,0:map%n-1))
    n = 0
    cy1 = 0
    do j = 1, q(2)
     cy2 = cy1 + ny(j)
     cx1 = 0
     do i = 1, q(1)
        cx2 = cx1 + nx(i)
        map%c1(:,n) = [cx1,cy1]
        map%c2(:,n) = [cx2,cy2] - 1
        map%cw(:,n) = [nx(i),ny(j)]
        n = n+1
        cx1 = cx2
     end do
     cy1 = cy2
    end do
    ! x-,y-ticks
    map%nxt = q(1); if( map%nxt==1 ) map%nxt = 0
    map%nyt = q(2); if( map%nyt==1 ) map%nyt = 0
    allocate(map%xticks(q(1)))
    allocate(map%yticks(q(2)))
    cx1 = 0
    do i = 1, q(1)
       map%xticks(i) = cx1
       cx1 = cx1 + nx(i)
    end do
    cy1 = 0
    do i = 1, q(2)
       map%yticks(i) = cy1
       cy1 = cy1 + ny(i)
    end do
    end function make_gridmap

end module patches


!    program main
!    use patches
!    implicit none
!    integer :: x,y,n,ticks(2)
!    type(gridmap) :: map
!    read(*,*) x,y,n
!    call optimal_gridmap([x,y],n,ticks)
!    map = gridmap([x,y],ticks)
!    write(*,*) ticks
!    write(*,*) map%nxt, map%xticks
!    write(*,*) map%nyt, map%yticks
!    do n = 0, map%n - 1 
!        write(*,*) map%c1(:,n),map%c2(:,n),map%cw(:,n)
!    end do
!    end program
