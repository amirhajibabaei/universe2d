
   module universe
   use iso_fortran_env, only: real32, int64, output_unit
   implicit none
   private
   public      pr,  pp2d, stack, particle, pi, latt_const
  
      integer, parameter      :: pr = real32 

      real(pr), parameter     :: pi = 3.14159265359_pr

      type particle
          integer             :: name   
          real(pr)            :: pos(2) 
      contains 
          procedure           :: move
      end type particle

      interface cube
         procedure               init_cube
      end interface cube
      
      type cube
         integer              :: cap                       
         integer              :: occ                       
         integer, allocatable :: idx(:)                    
         type(particle),                 &
                  allocatable :: hld(:)                    
         integer              :: su(0:8)                   
         integer              :: cc(2)                     
         real(pr)             :: llim(2)                   
         real(pr)             :: ulim(2)                   
         type(cube), pointer  :: sur1, sur2,&              
                                 sur3, sur4                
         integer              :: tag                       
      contains
         procedure            :: add    => add_cube        
         procedure            :: pop    => pop_cube        
         procedure            :: ext    => ext_cube        
         procedure            :: shrink => shrink_cube     
         procedure            :: show   => show_cube       
         procedure            :: rel    => relative_cube   
         procedure            :: abs    => absolute_cube   
      end type cube

      interface pp2d
         procedure               init_pp2d, load_pp2d, hex_pp2d, read_lammps_pp2d
      end interface

      type pp2d
         real(pr)             :: ll(2), lx, ly             
         real(pr)             :: ww(2), wx, wy             
         integer              :: nxy(2), nx, ny            
         integer              :: noc                       
         integer              :: nop                       
         type(cube), &
          allocatable         :: cup(:)                    
         integer              :: stat                      
         integer, allocatable :: cidx(:)                   
         integer, allocatable :: kidx(:)                   
         integer              :: lnop                      
      contains
         procedure            :: assoc    => associate_pp2d
         procedure            :: add      => add_pp2d      
         procedure            :: relative => relative_pp2d 
         procedure            :: absolute => absolute_pp2d 
         procedure            :: shrink   => shrink_pp2d   
         procedure            :: write    => write_pp2d    
         procedure            :: reserve  => reserve_pp2d  
         procedure            :: stage    => active_pp2d   
         procedure            :: original => original_pp2d 
         procedure            :: random   => random_pp2d   
         procedure            :: randcup 
         procedure            :: move     => move_pp2d
         procedure            :: zoom_on  => scube_pp2d 
         procedure            :: zoom_u   => scube_up_pp2d 
         procedure            :: zoom_r   => scube_right_pp2d 
         procedure            :: look_up 
         procedure            :: metro    => metropolis_pp2d
      end type pp2d

      integer, parameter      :: su_cubes(2,8) = &
                                 reshape( [1,0,0,1,-1,0,0,-1, &    ! 7 2 6
                                          1,-1,1,1,-1,1,-1,-1], &  ! 3 0 1
                                          [2,8] )                  ! 8 4 5
      real(pr)                :: su_ww(2,8) 

      integer, parameter      :: nstk=1000
      type stack
         integer              :: n
         integer              :: id(nstk)
         real(pr)             :: x(nstk)
         real(pr)             :: y(nstk)
         real(pr)             :: r(nstk)
         real(pr)             :: f(nstk)
         real(pr)             :: g(nstk)
         real(pr)             :: h(nstk)
      contains
         procedure            :: xx => x_stack
         procedure            :: yy => y_stack
         procedure            :: hexorder
         procedure            :: show => show_stack
         procedure            :: rotate => rotate_stack
      end type stack

   contains
     
      !//////////////////////////////////////
      !//////////////////////////////////////           particle
      !//////////////////////////////////////

      subroutine move(part,dir,delta)
      implicit none
      class(particle), intent(inout) :: part
      integer,           intent(in)  :: dir
      real(pr),          intent(in)  :: delta
      part%pos(dir) = part%pos(dir) + delta
      end subroutine move
   
      !//////////////////////////////////////
      !//////////////////////////////////////           cube
      !//////////////////////////////////////

      function init_cube(cap) result(cub)
      implicit none
      integer, intent(in) :: cap
      type(cube)          :: cub
      cub%cap = cap
      cub%occ = 0
      allocate( cub%idx(cap) )
      allocate( cub%hld(cap) )
      end function init_cube
      
      subroutine ext_cube(cub)
      implicit none
      class(cube), intent(inout)  :: cub
      integer,        allocatable :: idx(:)
      type(particle), allocatable :: hld(:)
      allocate(idx,source=[cub%idx,0])
      deallocate(cub%idx)
      call move_alloc(idx,cub%idx)
      allocate(hld,source=[cub%hld,cub%hld(1)])
      deallocate(cub%hld)
      call move_alloc(hld,cub%hld)
      cub%cap = cub%cap + 1
      end subroutine ext_cube

      subroutine shrink_cube(cub)
      implicit none
      class(cube), intent(inout)  :: cub
      integer,        allocatable :: idx(:)
      type(particle), allocatable :: hld(:)
      integer                     :: occ
      occ = cub%occ
      allocate(idx,source=cub%idx(1:occ))
      deallocate(cub%idx)
      call move_alloc(idx,cub%idx)
      allocate(hld,source=cub%hld(1:occ))
      deallocate(cub%hld)
      call move_alloc(hld,cub%hld)
      cub%cap = occ
      end subroutine shrink_cube

      subroutine add_cube(cub,idx,p)
      implicit none
      class(cube), intent(inout) :: cub
      integer,        intent(in) :: idx
      type(particle), intent(in) :: p
      integer                    :: k
      k = cub%occ + 1
      if( k>cub%cap ) call cub%ext()
      cub%hld(k) = p
      cub%idx(k) = idx
      cub%occ = k
      end subroutine add_cube
      
      subroutine pop_cube(cub,k,idx)
      implicit none
      class(cube), intent(inout) :: cub
      integer,       intent(in)  :: k
      integer,      intent(out)  :: idx
      idx = cub%idx(cub%occ)
      if( k<cub%occ ) then
         cub%idx(k) = idx
         cub%hld(k) = cub%hld(cub%occ)
      end if
      cub%occ = cub%occ - 1
      end subroutine pop_cube
       
      subroutine show_cube(cub)
      implicit none
      class(cube), intent(in) :: cub
      integer                 :: i
      write(*,*)
      write(*,'(i3,a,i3)') cub%occ, "/", cub%cap
      do i=1, cub%occ
         write(*,'(i7,6x,i7,2f8.3)') cub%idx(i), &
                   cub%hld(i)%name, &
                   cub%hld(i)%pos
      end do
      write(*,*)
      end subroutine show_cube

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ rel/abs

      subroutine relative_cube(cub)
      implicit none
      class(cube), intent(inout) :: cub
      cub%hld(1:cub%occ)%pos(1) = &
      cub%hld(1:cub%occ)%pos(1) - cub%llim(1)
      cub%hld(1:cub%occ)%pos(2) = &
      cub%hld(1:cub%occ)%pos(2) - cub%llim(2)
      end subroutine relative_cube

      subroutine absolute_cube(cub)
      implicit none
      class(cube), intent(inout) :: cub
      cub%hld(1:cub%occ)%pos(1) = &
      cub%hld(1:cub%occ)%pos(1) + cub%llim(1)
      cub%hld(1:cub%occ)%pos(2) = &
      cub%hld(1:cub%occ)%pos(2) + cub%llim(2)
      end subroutine absolute_cube


      !///////////////////////////////////////////
      !///////////////////////////////////////////        pp2d
      !///////////////////////////////////////////


      function init_pp2d(ll,wth) result(pos)   
      implicit none
      type(pp2d)           :: pos
      real(pr), intent(in) :: ll(2), wth
      integer              :: cx, cy, cl, i, j, k
      pos%ll = ll; pos%lx = ll(1); pos%ly = ll(2)
      pos%nxy = floor(ll/wth); pos%nx = pos%nxy(1); pos%ny = pos%nxy(2)
      pos%ww = ll/pos%nxy; pos%wx = pos%ww(1); pos%wy = pos%ww(2)
      pos%noc = pos%nx*pos%ny 
      allocate(pos%cup(0:pos%noc-1))
      do cy = 0, pos%ny-1
         do cx = 0, pos%nx-1
            cl = cy*pos%nx + cx
            pos%cup(cl) = cube(1) 
            pos%cup(cl)%cc = (/cx,cy/)
            pos%cup(cl)%llim = (/ cx   *pos%wx, cy   *pos%wy/)
            pos%cup(cl)%ulim = (/(cx+1)*pos%wx,(cy+1)*pos%wy/)
            pos%cup(cl)%su(0) = cl
            do k = 1, 8
               i = mod(cx+su_cubes(1,k),pos%nx)
               j = mod(cy+su_cubes(2,k),pos%ny)
               if(i<0) i = i + pos%nx
               if(j<0) j = j + pos%ny
               pos%cup(cl)%su(k) =  j*pos%nx + i
            end do
         end do
      end do
      pos%nop = 0    
      do k = 1, 8 
         su_ww(:,k) = su_cubes(:,k)*pos%ww
      end do
      pos%stat = 1 
      call pos%assoc()
      end function init_pp2d

      subroutine associate_pp2d(pos)  
      implicit none
      class(pp2d), intent(inout), target :: pos
      integer cl
      do cl = 0, pos%noc-1 
         pos%cup(cl)%sur1 => pos%cup(pos%cup(cl)%su(1))
         pos%cup(cl)%sur2 => pos%cup(pos%cup(cl)%su(2))
         pos%cup(cl)%sur3 => pos%cup(pos%cup(cl)%su(3))
         pos%cup(cl)%sur4 => pos%cup(pos%cup(cl)%su(4))
      end do
      end subroutine associate_pp2d

      subroutine add_pp2d(pos,idx,p) 
      implicit none
      class(pp2d), intent(inout)    :: pos
      integer,       intent(in)     :: idx
      type(particle), intent(inout) :: p
      integer                       :: cc(2), cl
      type(particle)                :: rp
      p%pos = mod(p%pos,pos%ll)
      where( p%pos<0.0_pr ) p%pos = p%pos+pos%ll
      cc = floor(p%pos/pos%ww)
      ! *** bug *** when r = l-eps; floor(r/w)=floor(l/w)=nxy exceeds dims
      ! fixed by this line: 
      where( cc>=pos%nxy ) cc=pos%nxy-1
      ! end bug 
      cl = cc(2)*pos%nx + cc(1)
      if( pos%stat==0 ) then
         call pos%cup(cl)%add(idx,p)
      else if( pos%stat==1 ) then
         rp = p
         rp%pos = rp%pos - pos%cup(cl)%llim
         call pos%cup(cl)%add(idx,rp)
      end if
      pos%nop = pos%nop + 1
      end subroutine add_pp2d

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

      function load_pp2d(handle, wth) result(pos) 
      implicit none
      class(*),         intent(in)    :: handle
      real(pr),         intent(in)    :: wth
      type(pp2d)                      :: pos
      integer                         :: uin, i, num
      real(pr)                        :: lx, ly
      type(particle)                  :: p
      logical                         :: doclose
      
      doclose = .false.
      select type(handle)
      type is (integer)
         uin = handle
      type is (character(len=*))
         open( newunit=uin, file=handle )
         doclose = .true.
      class default
         stop
      end select

      read(uin,*)
      read(uin,*)
      read(uin,*)
      read(uin,*) num
      read(uin,*) lx, ly
      pos = pp2d([lx,ly],wth)
      do i = 0, num-1
          read(uin,*) p
          call pos%add(i,p)
      end do

      if( doclose ) close(uin)
      end function load_pp2d


      subroutine write_pp2d(pos,handle,string,ints,reals) 
      implicit none
      class(pp2d), intent(in)        :: pos
      class(*),     intent(in)       :: handle
      character(len=*), intent(in), &
                         optional    :: string
      class(*), intent(in), optional :: ints(:)
      real(pr), intent(in), optional :: reals(:)
      integer                        :: i, j, uout
      ! 
      select type (handle)
      type is (integer)
         uout = handle
      type is (character(len=*))
         open(newunit=uout,file=handle)
      end select 
      ! header
      if( present(string) ) then
         write(uout,*) string
      else
         write(uout,*)
      end if
      if( present(ints) ) then
         select type (ints)
         type is (integer(int64))
            write(uout,*) ints
         type is (integer)
            write(uout,*) ints
         end select
      else
         write(uout,*)
      end if
      if( present(reals) ) then
         write(uout,*) reals
      else
         write(uout,*)
      end if
      ! body
      write(uout,*) pos%nop
      write(uout,*) pos%ll
      if( pos%stat==0) then
         do i = 0, pos%noc-1
            do j = 1, pos%cup(i)%occ
               write(uout,*) pos%cup(i)%hld(j)
            end do
         end do
      else if( pos%stat==1 ) then
          do i = 0, pos%noc-1
            do j = 1, pos%cup(i)%occ
               write(uout,*) pos%cup(i)%hld(j)%name, &
               pos%cup(i)%llim + pos%cup(i)%hld(j)%pos
            end do
         end do
      end if
      !
      select type (handle)
      type is (character(len=*))
         close(uout) 
      end select 
      end subroutine write_pp2d

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ lattices

      function hex_pp2d(nxy,rho,wth,f_distort) result(pos) 
      implicit none
      integer,  intent(in) :: nxy(2)
      real(pr), intent(in) :: rho, wth
      real(pr), intent(in),   &
             optional      :: f_distort
      type(pp2d)           :: pos
      real(pr)             :: v1(2), v2(2), a, ax, ay, rnd(2), delta
      integer              :: i, j, idx
      type(particle)       :: part
      ax = 0.5_pr
      ay = sqrt(3.0_pr)/2
      a = sqrt(1.0_pr/(rho*ay))
      v1 = a*[1.0_pr,0.0_pr]
      v2 = a*[ax,ay]
      pos = pp2d(nxy*[v1(1),v2(2)],wth)
      if( present(f_distort) ) then
              delta = f_distort
      else
              delta = 0.0_pr
      end if
      idx = -1
      do j = 0, nxy(2)-1
         do i = 0, nxy(1)-1
            call random_number(rnd) 
            idx = idx + 1
            part%name = idx
            part%pos  = i*v1 + j*v2 + ( rnd(1)*v1 + rnd(2)*v2 )*delta
            call pos%add(idx, part)
         end do
      end do
      end function hex_pp2d
      
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ read lammps

      function read_lammps_pp2d(fname, wth,lammps) result(pos) 
      implicit none
      character(len=*), intent(in)    :: fname,lammps
      real(pr),         intent(in)    :: wth
      type(pp2d)                      :: pos
      integer                         :: uin, i, num
      real(pr)                        :: lx, ly, skip
      type(particle)                  :: p
      if(lammps /= "lammps") then
         write(*,*)"not lammps";stop
      end if
      open( newunit=uin, file=fname )
         read(uin,*)
         read(uin,*)
         read(uin,*)
         read(uin,*) num
         read(uin,*)
         read(uin,*) skip, lx
         read(uin,*) skip, ly
         read(uin,*)
         read(uin,*)
         pos = pp2d([lx,ly],wth)
         do i = 0, num-1
             read(uin,*) p
             call pos%add(i,p)
         end do
      close(uin)
      end function read_lammps_pp2d

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      
      subroutine relative_pp2d(pos)   
      implicit none
      class(pp2d), intent(inout) :: pos
      integer                    :: cl 
      if( pos%stat==0 ) then
         do cl = 0, pos%noc-1
            call pos%cup(cl)%rel()
         end do
         pos%stat = 1
      end if
      end subroutine relative_pp2d


      subroutine absolute_pp2d(pos) 
      implicit none
      class(pp2d), intent(inout) :: pos
      integer                    :: cl
      if( pos%stat==1 ) then
         do cl = 0, pos%noc-1
            call pos%cup(cl)%abs()
         end do
         pos%stat = 0
      end if
      end subroutine absolute_pp2d

      subroutine shrink_pp2d(pos) 
      implicit none
      class(pp2d), intent(inout) :: pos
      integer                    :: cl
      do cl = 0, pos%noc-1
         call pos%cup(cl)%shrink()
      end do
      end subroutine shrink_pp2d


      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


      subroutine reserve_pp2d(pos,maxnum)    
      implicit none
      class(pp2d), intent(inout)    :: pos
      integer, intent(in), optional :: maxnum
      integer                       :: ub
      if( present(maxnum) ) then
         ub = maxnum - 1
      else
         ub = pos%nop - 1 
      end if
      if( .not.allocated(pos%cidx) ) then
          allocate(pos%cidx(0:ub))   
          allocate(pos%kidx(0:ub))  
      else if( size(pos%cidx)<ub+1 ) then
          deallocate(pos%cidx)   
          deallocate(pos%kidx)  
          allocate(pos%cidx(0:ub))   
          allocate(pos%kidx(0:ub))  
      end if
      end subroutine reserve_pp2d 


      subroutine active_pp2d(pos,c,cw)
      implicit none
      class(pp2d), intent(inout) :: pos
      integer,      intent(in),&
                      optional   :: c(2), cw(2)
      integer                    :: c1(2), c2(2)
      integer                    :: cx, cy, cl
      integer                    :: k, idx, kx, ky
      if(present(c) .and. present(cw)) then
         c1 = c; c2 = c1+cw-1
      else
         c1 = [0,0]; c2 = pos%nxy-1
      end if
      ! first count lnop and reserve
      idx = 0
      do cy = c1(2), c2(2)
         ky = mod(cy,pos%ny)
         if( ky<0 ) ky = ky + pos%ny
         do cx = c1(1), c2(1)
            kx = mod(cx,pos%nx)
            if( kx<0 ) kx = kx + pos%nx
            cl = ky*pos%nx + kx
            idx = idx + pos%cup(cl)%occ
         end do
      end do
      call pos%reserve(idx)
      pos%lnop = idx
      ! assign tags then
      pos%cup%tag = 1 
      idx = -1
      do cy = c1(2), c2(2)
         ky = mod(cy,pos%ny)
         if( ky<0 ) ky = ky + pos%ny
         do cx = c1(1), c2(1)
            kx = mod(cx,pos%nx)
            if( kx<0 ) kx = kx + pos%nx
            cl = ky*pos%nx + kx
            pos%cup(cl)%tag = 0 
            do k = 1, pos%cup(cl)%occ
               idx = idx + 1
               pos%cup(cl)%idx(k) = idx
               pos%cidx(idx) = cl
               pos%kidx(idx) = k
            end do
         end do
      end do
      end subroutine active_pp2d


      subroutine original_pp2d(pos) 
      implicit none
      class(pp2d), intent(inout) :: pos
      integer                    :: cl, k, idx
      call pos%reserve()
      pos%cup%tag = 0 
      do cl = 0, pos%noc-1
         do k = 1, pos%cup(cl)%occ
            idx = pos%cup(cl)%hld(k)%name
            pos%cup(cl)%idx(k) = idx
            pos%cidx(idx) = cl
            pos%kidx(idx) = k
         end do
      end do
      pos%lnop = pos%nop
      end subroutine original_pp2d
      
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      recursive &
      subroutine random_pp2d(pos,dmin,dmax,idx,dir,delta) 
      implicit none
      class(pp2d), intent(in) :: pos
      real(pr),    intent(in) :: dmin, dmax
      integer,    intent(out) :: idx, dir
      real(pr),   intent(out) :: delta
      real(pr)                :: rand_pr
      integer                 :: ncl, cl, k
      call random_number(rand_pr)
      idx = floor(rand_pr*pos%lnop)
      cl = pos%cidx(idx)
      k = pos%kidx(idx)
      call random_number(rand_pr)
      if(rand_pr<0.5_pr) then
         dir = 1
      else
         dir = 2
      end if
      call random_number(rand_pr)
      delta = dmin + rand_pr*(dmax-dmin)
      ! check if forbiden
      rand_pr = pos%cup(cl)%hld(k)%pos(dir) + delta 
      if( rand_pr>=pos%ww(dir) ) then
         ncl = pos%cup(cl)%su(dir)
         if( pos%cup(ncl)%tag==1 ) &
         call pos%random(dmin,dmax,idx,dir,delta)
      else if( rand_pr<0.0_pr ) then
         ncl = pos%cup(cl)%su(dir+2)
         if( pos%cup(ncl)%tag==1 ) &
         call pos%random(dmin,dmax,idx,dir,delta)
      end if
      end subroutine random_pp2d


      function randcup(pos) result(cc)
      implicit none
      class(pp2d), intent(in) :: pos
      integer                 :: cc(2)
      real(pr)                :: rnd(2)
      call random_number(rnd)
      cc = floor(rnd*pos%nxy)
      end function randcup

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     

      subroutine look_up(pos,idx,xy,env) 
      implicit none
      class(pp2d), intent(in), target :: pos
      integer,   intent(in), optional :: idx
      real(pr),      intent(in)       :: xy(2)
      type(stack), intent(out)        :: env
      integer                         :: i, j, cl, cc(2)
      real(pr)                        :: r(2)
      type(cube), pointer             :: ccub
      ! find r
      if( present(idx) ) then
         cl = pos%cidx(idx)
         j = pos%kidx(idx)
         r = pos%cup(cl)%llim + pos%cup(cl)%hld(j)%pos + xy
      else
         r = xy
      end if
      ! find cl of r
      r = mod(r,pos%ll)
      where( r<0.0_pr ) r = r+pos%ll
      cc = floor(r/pos%ww)
      where( cc>=pos%nxy ) cc=pos%nxy-1
      cl = cc(2)*pos%nx + cc(1)
      ! standardize r
      ccub => pos%cup(cl)
      r = r - ccub%llim
      ! look around then
      env%n = 0
      do i = 1, ccub%occ
          env%n = env%n + 1
          env%id(env%n) = ccub%idx(i) 
          env%x(env%n) = ccub%hld(i)%pos(1)
          env%y(env%n) = ccub%hld(i)%pos(2)
      end do
      ccub => ccub%sur1 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(1)<r(1) ) then 
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) + su_ww(1,1)
                   env%y(env%n) = ccub%hld(i)%pos(2) !+ su_ww(2,1)
                end if
             end do
      ccub => ccub%sur2 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(1)<r(1) .and. &
                     ccub%hld(i)%pos(2)<r(2) ) then
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) + su_ww(1,6)
                   env%y(env%n) = ccub%hld(i)%pos(2) + su_ww(2,6)
                end if
             end do
      ccub => ccub%sur3 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(2)<r(2) ) then 
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) !+ su_ww(1,2)
                   env%y(env%n) = ccub%hld(i)%pos(2) + su_ww(2,2)
                end if
             end do
      ccub => ccub%sur3 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(1)>r(1) .and. &
                     ccub%hld(i)%pos(2)<r(2) ) then
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) + su_ww(1,7)
                   env%y(env%n) = ccub%hld(i)%pos(2) + su_ww(2,7)
                end if
             end do
      ccub => ccub%sur4 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(1)>r(1) ) then 
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) + su_ww(1,3)
                   env%y(env%n) = ccub%hld(i)%pos(2) !+ su_ww(2,3)
                end if
             end do
      ccub => ccub%sur4 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(1)>r(1) .and. &
                     ccub%hld(i)%pos(2)>r(2) ) then
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) + su_ww(1,8)
                   env%y(env%n) = ccub%hld(i)%pos(2) + su_ww(2,8)
                end if
             end do
      ccub => ccub%sur1 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(2)>r(2) ) then 
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) !+ su_ww(1,4)
                   env%y(env%n) = ccub%hld(i)%pos(2) + su_ww(2,4)
                end if
             end do
      ccub => ccub%sur1 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(1)<r(1) .and. &
                     ccub%hld(i)%pos(2)>r(2) ) then
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) + su_ww(1,5)
                   env%y(env%n) = ccub%hld(i)%pos(2) + su_ww(2,5)
                end if
             end do
      env%x(1:env%n) = env%x(1:env%n)-r(1)
      env%y(1:env%n) = env%y(1:env%n)-r(2)
      end subroutine look_up

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     

      subroutine scube_pp2d(pos,idx,env) 
      implicit none
      class(pp2d), intent(in), &
                       target   :: pos
      integer,       intent(in) :: idx
      type(stack), intent(out)  :: env
      integer                   :: i, cl, k
      real(pr)                  :: r(2)
      type(cube), pointer       :: ccub
      cl = pos%cidx(idx)
      k = pos%kidx(idx)
      ccub => pos%cup(cl)
      r = ccub%hld(k)%pos
      env%n = 0
      do i = 1, ccub%occ
         if(i/=k) then
                 env%n = env%n + 1
                 env%id(env%n) = ccub%idx(i) 
                 env%x(env%n) = ccub%hld(i)%pos(1)
                 env%y(env%n) = ccub%hld(i)%pos(2)
         end if
      end do
      ccub => ccub%sur1 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(1)<r(1) ) then 
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) + su_ww(1,1)
                   env%y(env%n) = ccub%hld(i)%pos(2) !+ su_ww(2,1)
                end if
             end do
      ccub => ccub%sur2 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(1)<r(1) .and. &
                     ccub%hld(i)%pos(2)<r(2) ) then
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) + su_ww(1,6)
                   env%y(env%n) = ccub%hld(i)%pos(2) + su_ww(2,6)
                end if
             end do
      ccub => ccub%sur3 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(2)<r(2) ) then 
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) !+ su_ww(1,2)
                   env%y(env%n) = ccub%hld(i)%pos(2) + su_ww(2,2)
                end if
             end do
      ccub => ccub%sur3 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(1)>r(1) .and. &
                     ccub%hld(i)%pos(2)<r(2) ) then
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) + su_ww(1,7)
                   env%y(env%n) = ccub%hld(i)%pos(2) + su_ww(2,7)
                end if
             end do
      ccub => ccub%sur4 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(1)>r(1) ) then 
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) + su_ww(1,3)
                   env%y(env%n) = ccub%hld(i)%pos(2) !+ su_ww(2,3)
                end if
             end do
      ccub => ccub%sur4 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(1)>r(1) .and. &
                     ccub%hld(i)%pos(2)>r(2) ) then
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) + su_ww(1,8)
                   env%y(env%n) = ccub%hld(i)%pos(2) + su_ww(2,8)
                end if
             end do
      ccub => ccub%sur1 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(2)>r(2) ) then 
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) !+ su_ww(1,4)
                   env%y(env%n) = ccub%hld(i)%pos(2) + su_ww(2,4)
                end if
             end do
      ccub => ccub%sur1 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(1)<r(1) .and. &
                     ccub%hld(i)%pos(2)>r(2) ) then
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) + su_ww(1,5)
                   env%y(env%n) = ccub%hld(i)%pos(2) + su_ww(2,5)
                end if
             end do
      env%x(1:env%n) = env%x(1:env%n)-r(1)
      env%y(1:env%n) = env%y(1:env%n)-r(2)
      end subroutine scube_pp2d

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     

      subroutine scube_up_pp2d(pos,idx,env) 
      implicit none
      class(pp2d), intent(in), &
                       target   :: pos
      integer,       intent(in) :: idx
      type(stack), intent(out)  :: env
      integer                   :: i, cl, k
      real(pr)                  :: r(2)
      type(cube), pointer       :: ccub
      cl = pos%cidx(idx)
      k = pos%kidx(idx)
      ccub => pos%cup(cl)
      r = ccub%hld(k)%pos
      env%n = 0
      do i = 1, ccub%occ
         if( i/=k .and. ccub%hld(i)%pos(2)>r(2) ) then
                 env%n = env%n + 1
                 env%id(env%n) = ccub%idx(i) 
                 env%x(env%n) = ccub%hld(i)%pos(1)
                 env%y(env%n) = ccub%hld(i)%pos(2)
         end if
      end do
      ccub => ccub%sur1 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(1)<r(1) .and. &
                     ccub%hld(i)%pos(2)>r(2) ) then 
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) + su_ww(1,1)
                   env%y(env%n) = ccub%hld(i)%pos(2) !+ su_ww(2,1)
                end if
             end do
      ccub => ccub%sur2 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(1)<r(1) .and. &
                     ccub%hld(i)%pos(2)<r(2) ) then
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) + su_ww(1,6)
                   env%y(env%n) = ccub%hld(i)%pos(2) + su_ww(2,6)
                end if
             end do
      ccub => ccub%sur3 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(2)<r(2) ) then 
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) !+ su_ww(1,2)
                   env%y(env%n) = ccub%hld(i)%pos(2) + su_ww(2,2)
                end if
             end do
      ccub => ccub%sur3 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(1)>r(1) .and. &
                     ccub%hld(i)%pos(2)<r(2) ) then
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) + su_ww(1,7)
                   env%y(env%n) = ccub%hld(i)%pos(2) + su_ww(2,7)
                end if
             end do
      ccub => ccub%sur4 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(1)>r(1) .and. &
                     ccub%hld(i)%pos(2)>r(2) ) then 
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) + su_ww(1,3)
                   env%y(env%n) = ccub%hld(i)%pos(2) !+ su_ww(2,3)
                end if
             end do
      env%x(1:env%n) = env%x(1:env%n)-r(1)
      env%y(1:env%n) = env%y(1:env%n)-r(2)
      end subroutine scube_up_pp2d

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     

      subroutine scube_right_pp2d(pos,idx,env) 
      implicit none
      class(pp2d), intent(in), &
                       target   :: pos
      integer,       intent(in) :: idx
      type(stack), intent(out)  :: env
      integer                   :: i, cl, k
      real(pr)                  :: r(2)
      type(cube), pointer       :: ccub
      cl = pos%cidx(idx)
      k = pos%kidx(idx)
      ccub => pos%cup(cl)
      r = ccub%hld(k)%pos
      env%n = 0
      do i = 1, ccub%occ
         if(i/=k .and. ccub%hld(i)%pos(1)>r(1) ) then
                 env%n = env%n + 1
                 env%id(env%n) = ccub%idx(i) 
                 env%x(env%n) = ccub%hld(i)%pos(1)
                 env%y(env%n) = ccub%hld(i)%pos(2)
         end if
      end do
      ccub => ccub%sur4 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(2)>r(2) .and. &
                     ccub%hld(i)%pos(1)>r(1) ) then 
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) !+ su_ww(1,4)
                   env%y(env%n) = ccub%hld(i)%pos(2) + su_ww(2,4)
                end if
             end do
      ccub => ccub%sur1 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(1)<r(1) .and. &
                     ccub%hld(i)%pos(2)>r(2) ) then
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) + su_ww(1,5)
                   env%y(env%n) = ccub%hld(i)%pos(2) + su_ww(2,5)
                end if
             end do

      ccub => ccub%sur2 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(1)<r(1) ) then 
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) + su_ww(1,1)
                   env%y(env%n) = ccub%hld(i)%pos(2) !+ su_ww(2,1)
                end if
             end do
      ccub => ccub%sur2 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(1)<r(1) .and. &
                     ccub%hld(i)%pos(2)<r(2) ) then
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) + su_ww(1,6)
                   env%y(env%n) = ccub%hld(i)%pos(2) + su_ww(2,6)
                end if
             end do
      ccub => ccub%sur3 
             do i = 1, ccub%occ
                if( ccub%hld(i)%pos(2)<r(2) .and. &
                     ccub%hld(i)%pos(1)>r(1) ) then 
                   env%n = env%n + 1
                   env%id(env%n) = ccub%idx(i) 
                   env%x(env%n) = ccub%hld(i)%pos(1) !+ su_ww(1,2)
                   env%y(env%n) = ccub%hld(i)%pos(2) + su_ww(2,2)
                end if
             end do
      env%x(1:env%n) = env%x(1:env%n)-r(1)
      env%y(1:env%n) = env%y(1:env%n)-r(2)
      end subroutine scube_right_pp2d

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine move_pp2d(pos,idx,dir,delta)  
      implicit none
      class(pp2d), intent(inout), &
                          target :: pos
      integer,       intent(in)  :: idx, dir
      real(pr),      intent(in)  :: delta
      integer                    :: kk, side_eff, cl, k
      type(cube), pointer        :: cub
      cl = pos%cidx(idx)
      k = pos%kidx(idx)
      cub => pos%cup(cl)
      call cub%hld(k)%move(dir,delta)
      if( cub%hld(k)%pos(dir)>=pos%ww(dir) ) then
            call cub%hld(k)%move(dir,-pos%ww(dir))
            if( dir==1 ) then
               call cub%sur1%add(idx,cub%hld(k))
               kk = cub%sur1%occ
            else
               call cub%sur2%add(idx,cub%hld(k))
               kk = cub%sur2%occ
            end if
            call cub%pop( k, side_eff )
            pos%kidx( side_eff ) = k
            pos%cidx( idx ) = cub%su(dir)
            pos%kidx( idx ) = kk 
      else  if( cub%hld(k)%pos(dir)<0.0_pr ) then
            call cub%hld(k)%move(dir,pos%ww(dir))
            if( dir==1 ) then
               call cub%sur3%add(idx,cub%hld(k))
               kk = cub%sur3%occ
            else
               call cub%sur4%add(idx,cub%hld(k))
               kk = cub%sur4%occ
            end if
            call cub%pop( k, side_eff )
            pos%kidx( side_eff ) = k
            pos%cidx( idx ) = cub%su(dir+2)
            pos%kidx( idx ) = kk 
      end if
      end subroutine move_pp2d
   
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine metropolis_pp2d(pos,efunc,eargs,tem,dmax,steps,nsuccess)
      implicit none
      class(pp2d), intent(inout) :: pos
      interface 
         function efunc(x,y,args) result(en)
         import pr 
         implicit none
         real(pr), intent(in) :: x(:), y(:), args(:)
         real(pr)             :: en(size(x))
         end function
      end interface
      real(pr),    intent(in)    :: eargs(:), tem, dmax
      integer,     intent(in)    :: steps
      integer,     intent(out)   :: nsuccess
      integer                    :: step, idx, dir, n
      real(pr)                   :: delta, de, rnd
      type(stack)                :: env
      nsuccess = 0
      do step = 1, steps
         call pos%random(-dmax,dmax,idx,dir,delta)
         call pos%zoom_on(idx,env)
         n = env%n
         env%f(1:n) = efunc(env%x(1:n),env%y(1:n),eargs)
         if( dir==1 ) then
            env%g(1:n) = efunc(env%x(1:n)-delta,env%y(1:n),eargs)
         else
            env%g(1:n) = efunc(env%x(1:n),env%y(1:n)-delta,eargs)
         end if
         de = sum(env%g(1:n)-env%f(1:n))
         if( de<=0.0_pr ) then
            call pos%move(idx,dir,delta) ; nsuccess = nsuccess + 1
         else
            call random_number(rnd)
            if(rnd<=exp(-de/tem)) then
               call pos%move(idx,dir,delta) ; nsuccess = nsuccess + 1
            end if
         end if
      end do
      end subroutine metropolis_pp2d

      !//////////////////////////////////////
      !//////////////////////////////////////           stack
      !//////////////////////////////////////

      pure &
      function x_stack(env) result(x)
      implicit none
      class(stack), intent(in) :: env
      real(pr)                 :: x(env%n)
      x = env%x(1:env%n)
      end function x_stack

      pure &
      function y_stack(env) result(y)
      implicit none
      class(stack), intent(in) :: env
      real(pr)                 :: y(env%n)
      y = env%x(1:env%n)
      end function y_stack

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      pure &
      function hexorder(env) result(psi)
      implicit none
      class(stack), intent(in) :: env
      integer                  :: i, n, l(1)
      real(pr)                 :: psi(2), r, x, y, theta, &
                                  d(max(1,env%n))
      logical                  :: mask(max(1,env%n))
      n = env%n
      psi = 0.0_pr
      if( n>=6 ) then
         d = env%x(1:n)**2 + env%y(1:n)**2
         mask = .true.
         do i = 1, 6
            l = minloc(d,mask)
            r = sqrt(d(l(1)))
            x = env%x(l(1))
            y = env%y(l(1))
            theta = acos(x/r) 
            if( y<0.0_pr ) theta = 2*pi - theta
            psi = psi + [cos(6*theta),sin(6*theta)]
            mask(l(1)) = .false.
         end do
         psi = psi/6
      end if
      end function hexorder

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine rotate_stack(env,theta,matrix) 
      implicit none
      class(stack), intent(inout)    :: env
      real(pr), intent(in), optional :: theta, matrix(2,2)
      real(pr)                       :: mat(2,2), r(2)
      integer                        :: i
      if( present(matrix) ) then
         mat = matrix
      elseif( present(theta) ) then
         mat = reshape( [cos(theta), sin(theta), -sin(theta), cos(theta)], [2,2] )
      else
         write(output_unit,*)"rotate stack"; stop
      end if
      do i = 1, env%n
         r = matmul(mat,[env%x(i),env%y(i)])
         env%x(i) = r(1)
         env%y(i) = r(2)
      end do
      end subroutine rotate_stack

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine show_stack(env) 
      implicit none
      class(stack), intent(in) :: env
      integer                  :: i
      do i = 1, env%n
         write(*,*) env%id(i), env%x(i), env%y(i)
      end do
      end subroutine show_stack

      !//////////////////////////////////////
      !//////////////////////////////////////           other
      !//////////////////////////////////////

      function latt_const(rho) result(a)
      implicit none
      real(pr),   intent(in)  :: rho
      real(pr)                :: a
      a = sqrt(2.0_pr/(rho*sqrt(3.0_pr)))
      end function latt_const

   end module universe 
