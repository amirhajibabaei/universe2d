!
!  Author: Amir H. T.
!  Email: a.hajibabaei.86@gmail.com
!
!  use histpgram, only: hist1d
!  ...
!  type(hist1d) :: h
!  call h%init(x1,x2,dx)
!  call h%gather(x)        ! x is scalar or 1d array
!  call h%write(unit,[str],[ints],...)      
!
!
   module histogram
   use iso_fortran_env, only: real32, int64
   implicit none
   private
   public  hist1d , hist2d
   
      integer, parameter          :: pr = real32
      
      type              hist1d
            real(pr)              :: x1, x2, dx
            integer               :: bins
            integer, allocatable  :: count(:)
            integer               :: miss1, hit, miss2, total
         contains
            procedure             :: init, write, reset 
            procedure, private    :: gather_a, gather_s
            generic               :: gather => gather_s, gather_a
      end type          hist1d
   
      type              hist2d
            real(pr)              :: x1, x2, dx
            real(pr)              :: y1, y2, dy
            integer               :: xbins, ybins
            integer, allocatable  :: count(:,:)
            integer               :: xmiss1, ymiss1, hit, &
                                     xmiss2, ymiss2, total
         contains
            procedure             :: init   => init_2d
            procedure             :: write  => write_2d
            procedure             :: reset  => reset_2d
            procedure, private    :: gather_a_2d, gather_s_2d
            generic               :: gather => gather_s_2d, gather_a_2d
      end type          hist2d
  
   contains

      !////////////////////////////////////////////
      !////////////////////////////////////////////    hist1d
      !////////////////////////////////////////////
   
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

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine reset(h)
      implicit none
      class(hist1d), intent(inout) :: h
      h%count = 0
      h%miss1 = 0
      h%hit   = 0
      h%miss2 = 0
      h%total = 0
      end subroutine reset
  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      subroutine write(h,handle,string,ints,reals) 
      implicit none
      class(hist1d), intent(in)      :: h
      class(*),     intent(in)       :: handle
      character(len=*), intent(in), &
                         optional    :: string
      class(*), intent(in), optional :: ints(:)
      real(pr), intent(in), optional :: reals(:)
      integer                        :: i, uout
      character(len=20)              :: str1, str2, str3
      real(pr)                       :: p
      ! 
      select type (handle)
      type is (integer)
         uout = handle
      type is (character(len=*))
         open(newunit=uout,file=handle)
      end select 
      ! 2 empty lines
      write(uout,*)
      write(uout,*)
      ! header: 3 lines
      if( present(string) ) then
         write(uout,*) "# ", string
      else
         write(uout,*) "# "
      end if
      if( present(ints) ) then
         select type (ints)
         type is (integer(int64))
            write(uout,*) "# ", ints
         type is (integer)
            write(uout,*) "# ", ints
         end select
      else
         write(uout,*) "# "
      end if
      if( present(reals) ) then
         write(uout,*) "# ", reals
      else
         write(uout,*) "# "
      end if
      ! abstract: 3 lines
      write(str1,*)h%x1  
      write(str2,*)h%x1 + (h%bins+1)*h%dx 
      str1 = "# x<"//trim(adjustl(str1))
      str2 = "# x>="//trim(adjustl(str2))
      str3 = "# else"  
      write(uout,*) str1, h%miss1, real(h%miss1)/h%total
      write(uout,*) str2, h%miss2, real(h%miss2)/h%total
      write(uout,*) str3, h%hit, real(h%hit)/h%total
      write(uout,*) "# ", h%bins + 1, h%x1, h%dx
      ! body
      do i = 0, h%bins
         p = real(h%count(i))/h%total
         write(uout,*) h%x1 + i*h%dx, h%count(i), p, p/h%dx 
      end do
      !
      select type (handle)
      type is (character(len=*))
         close(uout) 
      end select 
      end subroutine write

      !////////////////////////////////////////////
      !////////////////////////////////////////////    hist2d
      !////////////////////////////////////////////
   
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
      subroutine init_2d(h,x1,x2,dx,y1,y2,dy)
      implicit none
      class(hist2d), intent(inout)   :: h
      real(pr),      intent(in)      :: x1, x2, dx
      real(pr),      intent(in)      :: y1, y2, dy
      h%x1   = min(x1,x2)                ; h%y1   = min(y1,y2)
      h%x2   = max(x1,x2)                ; h%y2   = max(y1,y2)
      h%dx   = abs(dx)                   ; h%dy   = abs(dy)
      h%xbins = floor( (h%x2-h%x1)/h%dx ); h%ybins = floor( (h%y2-h%y1)/h%dy )
      allocate(h%count(0:h%xbins,0:h%ybins))
      h%count  = 0
      h%xmiss1 = 0; h%ymiss1 = 0
      h%hit    = 0
      h%xmiss2 = 0; h%ymiss2 = 0
      h%total  = 0
      end subroutine init_2d
   
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
      subroutine gather_s_2d(h,x,y)
      implicit none
      class(hist2d), intent(inout) :: h
      real(pr),      intent(in)    :: x, y
      integer                      :: k, j
      k = floor( (x-h%x1)/h%dx )
      if( k<0 ) then
         h%xmiss1 = h%xmiss1 + 1
      else if( k>h%xbins ) then
         h%xmiss2 = h%xmiss2 + 1
      else
          j = floor( (y-h%y1)/h%dy )
          if( j<0 ) then
             h%ymiss1 = h%ymiss1 + 1
          else if( j>h%ybins ) then
             h%ymiss2 = h%ymiss2 + 1
          else
             h%hit = h%hit + 1
             h%count(k,j) = h%count(k,j) + 1
          end if
      end if
      h%total = h%total + 1
      end subroutine gather_s_2d
   
      subroutine gather_a_2d(h,x,y)
      implicit none
      class(hist2d), intent(inout) :: h
      real(pr),      intent(in)    :: x(:),y(:)
      integer                      :: i
      do i = 1, size(x)
         call gather_s_2d(h,x(i),y(i))
      end do
      end subroutine gather_a_2d

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine reset_2d(h)
      implicit none
      class(hist2d), intent(inout) :: h
      h%count  = 0
      h%xmiss1 = 0
      h%ymiss1 = 0
      h%hit    = 0
      h%xmiss2 = 0
      h%ymiss2 = 0
      h%total  = 0
      end subroutine reset_2d
  
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      subroutine write_2d(h,handle,string,ints,reals) 
      implicit none
      class(hist2d), intent(in)      :: h
      class(*),     intent(in)       :: handle
      character(len=*), intent(in), &
                         optional    :: string
      class(*), intent(in), optional :: ints(:)
      real(pr), intent(in), optional :: reals(:)
      integer                        :: i, j, uout
      character(len=20)              :: str1, str2, str3, &
                                        str4, str5
      real(pr)                       :: p
      ! 
      select type (handle)
      type is (integer)
         uout = handle
      type is (character(len=*))
         open(newunit=uout,file=handle)
      end select 
      ! 2 empty lines
      write(uout,*)
      write(uout,*)
      ! header: 3 lines
      if( present(string) ) then
         write(uout,*) "# ", string
      else
         write(uout,*) "# "
      end if
      if( present(ints) ) then
         select type (ints)
         type is (integer(int64))
            write(uout,*) "# ", ints
         type is (integer)
            write(uout,*) "# ", ints
         end select
      else
         write(uout,*) "# "
      end if
      if( present(reals) ) then
         write(uout,*) "# ", reals
      else
         write(uout,*) "# "
      end if
      ! abstract: 3 lines
      write(str1,*)h%x1  
      write(str2,*)h%x1 + (h%xbins+1)*h%dx 
      str1 = "# x<"//trim(adjustl(str1))
      str2 = "# x>="//trim(adjustl(str2))
      write(str4,*)h%y1  
      write(str5,*)h%y1 + (h%ybins+1)*h%dy 
      str4 = "# y<"//trim(adjustl(str4))
      str5 = "# y>="//trim(adjustl(str5))
      str3 = "# else"  
      write(uout,*) str1, h%xmiss1, real(h%xmiss1)/h%total
      write(uout,*) str2, h%xmiss2, real(h%xmiss2)/h%total
      write(uout,*) "# x in range, but"
      write(uout,*) str4, h%ymiss1, real(h%ymiss1)/h%total
      write(uout,*) str5, h%ymiss2, real(h%ymiss2)/h%total
      write(uout,*) str3, h%hit, real(h%hit)/h%total
      write(uout,*) "# ", h%xbins + 1, h%x1, h%dx
      write(uout,*) "# ", h%ybins + 1, h%y1, h%dy
      ! body
      do j = 0, h%ybins
         do i = 0, h%xbins
            p = real(h%count(i,j))/h%total
            write(uout,*) h%x1 + i*h%dx, h%y1 + j*h%dy, & 
                          h%count(i,j), p, p/(h%dx*h%dy)
         end do
      end do
      !
      select type (handle)
      type is (character(len=*))
         close(uout) 
      end select 
      end subroutine write_2d

   end module histogram

   ! example:
   
!   program main
!   use histogram, only: hist1d
!   use iso_fortran_env, only: real32, output_unit
!   implicit none
!   integer, parameter :: pr = real32
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
!   call h%write(output_unit)
!   end program main


!   program main
!   use histogram, only: hist2d
!   use iso_fortran_env, only: real32, output_unit
!   implicit none
!   integer, parameter :: pr = real32
!   type(hist2d) :: h
!   real(pr) :: x(10), xx
!   real(pr) :: y(10), yy
!   integer  :: i
!   call h%init(0.1_pr, 0.6_pr, 0.1_pr, 0.2_pr, 0.8_pr, 0.2_pr)
!   do i=1,10**5
!      call random_number(x)
!      call random_number(y)
!      call h%gather(x,y)
!      call random_number(xx)
!      call random_number(yy)
!      call h%gather(xx,yy)
!   end do
!   call h%write(output_unit)
!   end program main
