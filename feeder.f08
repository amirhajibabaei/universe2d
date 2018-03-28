! rational for make:
!          check for seed.txt
!                             if yes, read seed and pos
!                             if no, 
!          check for seed.abs
!                             if yes, read seed and make pos and seed.txt
!                             if no,
!          stop
!
module seed_md
use universe, only: pr, pp2d
use iso_fortran_env, only: int64
implicit none
private
public seed

   type seed
      character(len=5) :: names(5) = ["tem  ","rho  ","nx   ","rc   ","dmax:"]
      real(pr)         :: tem 
      real(pr)         :: rho
      integer          :: nx
      real(pr)         :: rc
      real(pr)         :: dmax
   contains
      procedure        :: make, str, dump
   end type seed

contains

   function str(sd) result(string)
   implicit none
   class(seed), intent(in)    :: sd
   character(len=1000)        :: st
   character(:), allocatable  :: string
   select type (sd)
   type is (seed)
      write(st,*) sd
   end select
   string = trim(st)
   end function str

   subroutine make(sd,pos,timestamp)
   implicit none
   class(seed), intent(inout) :: sd
   type(pp2d),    intent(out) :: pos
   integer,     intent(out)   :: timestamp
   real(pr)                   :: wth
   logical                    :: yes
   integer                    :: u
   character(len=20)          :: line(4)
   inquire(file="seed.txt",exist=yes)
   if( yes ) then
      open(newunit=u,file="seed.txt",status='old')
         select type (sd)
         type is (seed)
            read(u,*) sd
         end select
         read(u,*) timestamp
      close(u)
      wth = sd%rc + sd%dmax
      pos = pp2d("seed.txt",wth)
   else
      inquire(file="seed.abs",exist=yes)
      if( yes ) then
        open(newunit=u,file="seed.abs",status='old')
            read(u,*) line; read(line(4),*) sd%tem
            read(u,*) line; read(line(4),*) sd%rho
            read(u,*) line; read(line(4),*) sd%nx
            read(u,*) line; read(line(4),*) sd%rc
            read(u,*) line; read(line(4),*) sd%dmax
         close(u)
         wth = sd%rc + sd%dmax
         pos = pp2d([sd%nx,sd%nx],sd%rho,wth)
         call pos%write("seed.txt",string=sd%str(),ints=[0])
         timestamp = 0
      else
         write(*,*) "no seed found!"
         stop
      end if
   end if
   end subroutine make

   subroutine dump(sd,pos,timestamp)
   implicit none
   class(seed), intent(in) :: sd
   class(pp2d), intent(in) :: pos
   integer,     intent(in) :: timestamp
   call pos%write("seed.txt",string=sd%str(),ints=[timestamp])
   end subroutine dump

end module seed_md
