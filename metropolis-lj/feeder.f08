! rational for make:
!          check for mc_restart.txt
!                             if yes, read seed, pos, timestamp
!                             if no, 
!          check for seed.txt
!                             if yes, read seed, make pos, put timestamp = 0
!                             if no,
!          stop
!
! modification: we do not want seed to write the restart file in make
!               because of conflict with mpi, but seed%dump can be called
!               from main (if rank is 0).
!
   module seed_md
   use universe, only: pr, pp2d
   use iso_fortran_env, only: int64
   implicit none
   private
   public seed
   
      type seed
         character(len=6) :: names(7) = ["alpha ","tem   ","rho   ", &
                             "nx    ","rc    ","dmax  ","tdamp "]
         real(pr)         :: alpha 
         real(pr)         :: tem 
         real(pr)         :: rho
         integer          :: nx
         real(pr)         :: rc
         real(pr)         :: dmax
         integer          :: tdamp
      contains
         procedure         :: make, str, dump, read
         procedure, nopass :: open
      end type seed
   
   contains

      subroutine read(sd)
      implicit none
      class(seed), intent(inout)  :: sd
      integer                    :: u
      character(len=20)          :: line(4)
      open(newunit=u,file="seed.txt",status='old')
          read(u,*) line; read(line(4),*) sd%alpha
          read(u,*) line; read(line(4),*) sd%tem
          read(u,*) line; read(line(4),*) sd%rho
          read(u,*) line; read(line(4),*) sd%nx
          read(u,*) line; read(line(4),*) sd%rc
          read(u,*) line; read(line(4),*) sd%dmax
          read(u,*) line; read(line(4),*) sd%tdamp
      close(u)
      end subroutine read

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
      class(seed), intent(inout)  :: sd
      type(pp2d),    intent(out)  :: pos
      integer(int64), intent(out) :: timestamp
      real(pr)                    :: wth
      logical                     :: yes
      integer                     :: u
      character(len=20)           :: line(4)
      inquire(file="mc_restart.txt",exist=yes)
      if( yes ) then
         open(newunit=u,file="mc_restart.txt",status='old')
            select type (sd)
            type is (seed)
               read(u,*) sd
            end select
            read(u,*) timestamp
         close(u)
         wth = sd%rc + sd%dmax
         pos = pp2d("mc_restart.txt",wth)
      else
         inquire(file="seed.txt",exist=yes)
         if( yes ) then
           open(newunit=u,file="seed.txt",status='old')
               read(u,*) line; read(line(4),*) sd%alpha
               read(u,*) line; read(line(4),*) sd%tem
               read(u,*) line; read(line(4),*) sd%rho
               read(u,*) line; read(line(4),*) sd%nx
               read(u,*) line; read(line(4),*) sd%rc
               read(u,*) line; read(line(4),*) sd%dmax
               read(u,*) line; read(line(4),*) sd%tdamp
            close(u)
            wth = sd%rc + sd%dmax
            pos = pp2d([sd%nx,sd%nx],sd%rho,wth)
            timestamp = 0
         else
            write(*,*) "no seed found!"
            stop
         end if
      end if
      end subroutine make

      subroutine open(uscalars,uvectors)
      implicit none
      integer,      intent(out)  :: uscalars, uvectors
      logical                    :: yes
      inquire(file="mc_scalars.txt",exist=yes)
      if( yes ) then
         open(newunit=uscalars,file="mc_scalars.txt",status="old",action="write",access="append")
         open(newunit=uvectors,file="mc_vectors.txt",status="old",action="write",access="append")
      else
         open(newunit=uscalars,file="mc_scalars.txt",status="new",action="write")
         open(newunit=uvectors,file="mc_vectors.txt",status="new",action="write")
      end if
      end subroutine open
   
      subroutine dump(sd,pos,timestamp)
      implicit none
      class(seed),    intent(in) :: sd
      class(pp2d),    intent(in) :: pos
      integer(int64), intent(in) :: timestamp
      call pos%write("mc_restart.txt",string=sd%str(),ints=[timestamp])
      end subroutine dump

   end module seed_md
