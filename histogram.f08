!
!  Author: Amir H. T.
!  Email: a.hajibabaei.86@gmail.com
!
!  use histpgram, only: hist1d
!  ...
!  type(hist1d) :: h
!  call h%init(x1,x2,dx)
!  call h%gather(x)        ! x is scalar or 1d array
!  call h%write(unit)      ! unit optional
!
!
   module histogram
   use iso_fortran_env, only: real32, output_unit
   implicit none
   private
   public  hist1d 
   
      integer, parameter          :: pr = real32
      
      type              hist1d
            real(pr)              :: x1, x2, dx
            integer               :: bins
            integer, allocatable  :: count(:)
            integer               :: miss1, hit, miss2, total
         contains
            procedure             :: init, write 
            procedure, private    :: gather_a, gather_s
            generic               :: gather => gather_s, gather_a
      end type          hist1d
   
   
   contains
   
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
      subroutine init(h,x1,x2,dx)
      implicit none
      class(hist1d), intent(inout)   :: h
      real(pr),      intent(in)      :: x1, x2, dx
      h%x1   = min(x1,x2)
      h%x2   = max(x1,x2)
      h%dx   = abs(dx)
      h%bins = floor( (h%x2-h%x1)/h%dx )
      allocate(h%count(0:h%bins))
      h%count = 0
      h%miss1 = 0
      h%hit   = 0
      h%miss2 = 0
      h%total = 0
      end subroutine init
   
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
      subroutine gather_s(h,x)
      implicit none
      class(hist1d), intent(inout) :: h
      real(pr),      intent(in)    :: x
      integer                      :: k
      k = floor( (x-h%x1)/h%dx )
      if( k<0 ) then
         h%miss1 = h%miss1 + 1
      else if( k>h%bins ) then
         h%miss2 = h%miss2 + 1
      else
         h%hit = h%hit + 1
         h%count(k) = h%count(k) + 1
      end if
      h%total = h%total + 1
      end subroutine gather_s
   
      subroutine gather_a(h,x)
      implicit none
      class(hist1d), intent(inout) :: h
      real(pr),      intent(in)    :: x(:)
      integer                      :: i
      do i = 1, size(x)
         call gather_s(h,x(i))
      end do
      end subroutine gather_a
   
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
      subroutine write(h,uout)
      implicit none 
      class(hist1d), intent(in)     :: h
      integer, intent(in), optional :: uout
      integer                       :: i, u
      real(pr)                      :: p
      if( present(uout) ) then
              u = uout
      else
              u= output_unit
      end if
      !
      write(u,*)
      write(u,*) "# x <", h%x1, h%miss1, real(h%miss1)/h%total
      do i = 0, h%bins
         p = real(h%count(i))/h%total
         write(u,*) h%x1 + i*h%dx, h%count(i), p, p/h%dx 
      end do
      write(u,*) "# x >=", h%x1 + (h%bins+1)*h%dx, h%miss2, real(h%miss2)/h%total
      write(u,*)
      end subroutine write
   
   end module histogram

!   ! example:
!   
!   program main
!   use histogram, only: hist1d, pr
!   implicit none
!   type(hist1d) :: h
!   real(pr) :: x(10), xx
!   integer  :: i
!   call h%init(0.1_pr, 0.6_pr, 0.1_pr)
!   do i=1,10**5
!      call random_number(x)
!      call h%gather(x)
!      call random_number(xx)
!      call h%gather(xx)
!   end do
!   call h%write()
!   end program main
