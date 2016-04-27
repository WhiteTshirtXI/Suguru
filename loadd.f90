module mod_loadflds
use decomp_2d
use decomp_2d_io
use mod_param 
use mod_common_mpi
use mod_common
implicit none
private
public loadflds
contains
!
subroutine loadflds(in,nr)
implicit none
integer :: in,nr
integer :: fh
integer(kind=MPI_OFFSET_KIND) :: filesize,disp
character(len=7) :: istepchar
real, dimension(3) :: fldinfo
real, dimension(imax,jmax,kmax) :: temp
real, dimension(imaxc,jmaxc,kmaxc) :: tempc
!
if (in.eq.0) then
  write(istepchar,'(i7.7)') nr
  call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//'fldbg'//istepchar, &
       MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
  call MPI_FILE_CLOSE(fh,error)
  call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//'fldbg'//istepchar, &
       MPI_MODE_RDONLY, MPI_INFO_NULL,fh, error)
  disp = 0_MPI_OFFSET_KIND
  call decomp_2d_read_var(fh,disp,3,temp)
  unew(1:imax,1:jmax,1:kmax) = temp
  call decomp_2d_read_var(fh,disp,3,temp)
  vnew(1:imax,1:jmax,1:kmax) = temp
  call decomp_2d_read_var(fh,disp,3,temp)
  wnew(1:imax,1:jmax,1:kmax) = temp
  call decomp_2d_read_var(fh,disp,3,temp)
  uo(1:imax,1:jmax,1:kmax) = temp
  call decomp_2d_read_var(fh,disp,3,temp)
  vo(1:imax,1:jmax,1:kmax) = temp
  call decomp_2d_read_var(fh,disp,3,temp)
  wo(1:imax,1:jmax,1:kmax) = temp
  call decomp_2d_read_var(fh,disp,3,temp)
  pnew(1:imax,1:jmax,1:kmax) = temp
  call decomp_2d_read_scalar(fh,disp,3,fldinfo)
  time = fldinfo(1)
  nr = int(fldinfo(2))
  dt = fldinfo(3)
  call MPI_FILE_CLOSE(fh,error)

  if (ifoverlap==1) then
  call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//'fldol'//istepchar, &
       MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
  call MPI_FILE_CLOSE(fh,error)
  call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//'fldol'//istepchar, &
       MPI_MODE_RDONLY, MPI_INFO_NULL,fh, error)
  disp = 0_MPI_OFFSET_KIND
  call decomp_2d_read_var(fh,disp,3,tempc)
  unewc(1:imaxc,1:jmaxc,1:kmaxc) = tempc
  call decomp_2d_read_var(fh,disp,3,tempc)
  vnewc(1:imaxc,1:jmaxc,1:kmaxc) = tempc
  call decomp_2d_read_var(fh,disp,3,tempc)
  wnewc(1:imaxc,1:jmaxc,1:kmaxc) = tempc
  call decomp_2d_read_var(fh,disp,3,tempc)
  uoc(1:imaxc,1:jmaxc,1:kmaxc) = tempc
  call decomp_2d_read_var(fh,disp,3,tempc)
  voc(1:imaxc,1:jmaxc,1:kmaxc) = tempc
  call decomp_2d_read_var(fh,disp,3,tempc)
  woc(1:imaxc,1:jmaxc,1:kmaxc) = tempc
  call decomp_2d_read_var(fh,disp,3,tempc)
  pnewc(1:imaxc,1:jmaxc,1:kmaxc) = tempc
  call decomp_2d_read_scalar(fh,disp,3,fldinfo)
  time = fldinfo(1)
  nr = int(fldinfo(2))
  dt = fldinfo(3)
  call MPI_FILE_CLOSE(fh,error)
  endif

endif
!
if (in.eq.1) then
  write(istepchar,'(i7.7)') nr
  fldinfo = (/time,1.*nr,dt/)
  call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//'fldbg'//istepchar, &
       MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
  filesize = 0_MPI_OFFSET_KIND
  call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
  disp = 0_MPI_OFFSET_KIND
  temp = unew(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
  temp = vnew(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
  temp = wnew(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
  temp = uo(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
  temp = vo(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
  temp = wo(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
  temp = pnew(1:imax,1:jmax,1:kmax)
  call decomp_2d_write_var(fh,disp,3,temp)
  call decomp_2d_write_scalar(fh,disp,3,fldinfo)
  call MPI_FILE_CLOSE(fh,error)

  if (ifoverlap==1) then
  fldinfo = (/time,1.*nr,dt/)
  call MPI_FILE_OPEN(MPI_COMM_WORLD, datadir//'fldol'//istepchar, &
       MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, error)
  filesize = 0_MPI_OFFSET_KIND
  call MPI_FILE_SET_SIZE(fh,filesize,error)  ! guarantee overwriting
  disp = 0_MPI_OFFSET_KIND
  tempc = unewc(1:imaxc,1:jmaxc,1:kmaxc)
  call decomp_2d_write_var(fh,disp,3,tempc)
  tempc = vnewc(1:imaxc,1:jmaxc,1:kmaxc)
  call decomp_2d_write_var(fh,disp,3,tempc)
  tempc = wnewc(1:imaxc,1:jmaxc,1:kmaxc)
  call decomp_2d_write_var(fh,disp,3,tempc)
  tempc = uoc(1:imaxc,1:jmaxc,1:kmaxc)
  call decomp_2d_write_var(fh,disp,3,tempc)
  tempc = voc(1:imaxc,1:jmaxc,1:kmaxc)
  call decomp_2d_write_var(fh,disp,3,tempc)
  tempc = woc(1:imaxc,1:jmaxc,1:kmaxc)
  call decomp_2d_write_var(fh,disp,3,tempc)
  tempc = pnewc(1:imaxc,1:jmaxc,1:kmaxc)
  call decomp_2d_write_var(fh,disp,3,tempc)
  call decomp_2d_write_scalar(fh,disp,3,fldinfo)
  call MPI_FILE_CLOSE(fh,error)
  endif

endif
!
return
end subroutine loadflds
!
end module mod_loadflds
!
module mod_loadpart
use mpi
use mod_param
use mod_common
use mod_common_mpi
implicit none
private
public loadpart
contains
!
subroutine loadpart(in,nr)
implicit none
type particle_restart
 real :: x,y,z,integralv, &
         volume,volume0 
 real, dimension(nl) :: xfp,yfp,zfp,    &
                        dxdt,dydt,dzdt, & 
                        fxl,fyl,fzl,    &
                        ul,vl,wl,       &
                        xfpo,yfpo,zfpo, &
                        xfpold,yfpold,zfpold, &                        
                        xfp0,yfp0,zfp0,area
end type particle_restart
type(particle_restart), dimension(np) :: glob,glob_all

integer, parameter :: skip = 6+22*nl
integer p,i,idp
integer :: fh
integer     in,nr
integer :: proccoords(1:ndims),procrank
real :: leftbound,rightbound,frontbound,backbound
real :: xdum
integer :: lenr
real :: dist,angle
real :: xp,yp
real :: ax
real :: ay
integer :: count_mstr_all
character(len=7) :: istepchar
character(len=4) rankpr
character(len=11) :: tempstr2
integer :: counter,count_slve_loc
integer(kind=MPI_OFFSET_KIND) :: filesize,disp
!
inquire (iolength=lenr) xdum
write(istepchar,'(i7.7)') nr
!
if (in.eq.0) then
  if(myid.eq.0) then
    open(20,file=datadir//'allpartdata'//istepchar,access='direct',recl=np*skip*lenr)
    read(20,rec=1) glob
    close(20)
  endif
  !
  ! rank 0 broadcasts the particle data to all the others
  !
  call MPI_BCAST(glob(1)%x,np*skip,MPI_REAL8,0,comm_cart,error)
  ! 
  ! Determine master and slave processes for each particle.
  !
  ! initialisation
  !
ap(1:pmax)%mslv = 0
ap(1:pmax)%x = 0.
ap(1:pmax)%y = 0.
ap(1:pmax)%z = 0.
ap(1:pmax)%integralv = 0.
ap(1:pmax)%volume = 0.
ap(1:pmax)%volume0 = 0.

forall (p=1:pmax)
  ap(p)%xfp(:)  = 0.
  ap(p)%yfp(:)  = 0.
  ap(p)%zfp(:)  = 0.
  ap(p)%ul(:)   = 0.
  ap(p)%vl(:)   = 0.
  ap(p)%wl(:)   = 0.
  ap(p)%fxl(:)  = 0.
  ap(p)%fyl(:)  = 0.
  ap(p)%fzl(:)  = 0.
  ap(p)%dxdt(:) = 0.
  ap(p)%dydt(:) = 0.
  ap(p)%dzdt(:) = 0.
  ap(p)%xfpo(:)  = 0.
  ap(p)%yfpo(:)  = 0.
  ap(p)%zfpo(:)  = 0.
  ap(p)%xfpold(:)  = 0.
  ap(p)%yfpold(:)  = 0.
  ap(p)%zfpold(:)  = 0.
  ap(p)%xfp0(:)  = 0.
  ap(p)%yfp0(:)  = 0.
  ap(p)%zfp0(:)  = 0.
  ap(p)%area(:)  = 0.
  ap(p)%nb(1:8)   = 0
end forall
  ax = 0.5
  ay = 0.5
  i = 0
  npmstr = 0
  do p=1,np
    if (glob(p)%x.lt.0..or.glob(p)%x.gt.lx .or. &
        glob(p)%y.lt.0..or.glob(p)%y.gt.ly .or. &
        glob(p)%z.lt.0..or.glob(p)%z.gt.lz) then
      if (myid.eq.0) then
        write(6,*) 'Fatal error in initialisation of particle positions - '
        write(6,*) 'particle outside the domain!'
        write(6,*) 'Program aborted...'
      endif
      call mpi_finalize(error)
      stop
    endif
    if (glob(p)%x.eq.lx) ax = 0.51
    if (glob(p)%x.eq.0) ax = 0.49
    if (glob(p)%y.eq.ly) ay = 0.51
    if (glob(p)%y.eq.0) ay = 0.49

    proccoords(1) = nint(glob(p)%x*dims(1)/lx - ax)
    proccoords(2) = nint(glob(p)%y*dims(2)/ly - ay)
    leftbound     = (proccoords(1)  )*lx/(1.*dims(1)) ! left  boundary of particle's master
    rightbound    = (proccoords(1)+1)*lx/(1.*dims(1)) ! right boundary of particle's master
    frontbound    = (proccoords(2)  )*ly/(1.*dims(2)) ! front boundary of particle's master
    backbound     = (proccoords(2)+1)*ly/(1.*dims(2)) ! back  boundary of particle's master
    call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
    if (myid.eq.procrank) then
      npmstr = npmstr + 1
      i = i + 1



      ap(i)%mslv = p
      ap(i)%x = glob(p)%x
      ap(i)%y = glob(p)%y
      ap(i)%z = glob(p)%z
      ap(i)%integralv = glob(p)%integralv
      ap(i)%volume    = glob(p)%volume
      ap(i)%volume0   = glob(p)%volume0
      ap(i)%xfp(1:NL)   = glob(p)%xfp(1:NL)
      ap(i)%yfp(1:NL)   = glob(p)%yfp(1:NL)
      ap(i)%zfp(1:NL)   = glob(p)%zfp(1:NL)
      ap(i)%ul(1:NL)    = glob(p)%ul(1:NL)
      ap(i)%vl(1:NL)    = glob(p)%vl(1:NL)
      ap(i)%wl(1:NL)    = glob(p)%wl(1:NL)
      ap(i)%fxl(1:NL)   = glob(p)%fxl(1:NL)
      ap(i)%fyl(1:NL)   = glob(p)%fyl(1:NL)
      ap(i)%fzl(1:NL)   = glob(p)%fzl(1:NL)
      ap(i)%dxdt(1:NL)  = glob(p)%dxdt(1:NL)
      ap(i)%dydt(1:NL)  = glob(p)%dydt(1:NL)
      ap(i)%dzdt(1:NL)  = glob(p)%dzdt(1:NL)
      ap(i)%xfpo(1:NL)  = glob(p)%xfpo(1:NL)
      ap(i)%yfpo(1:NL)  = glob(p)%yfpo(1:NL)
      ap(i)%zfpo(1:NL)  = glob(p)%zfpo(1:NL)
      ap(i)%xfpold(1:NL)= glob(p)%xfpold(1:NL)
      ap(i)%yfpold(1:NL)= glob(p)%yfpold(1:NL)
      ap(i)%zfpold(1:NL)= glob(p)%zfpold(1:NL)
      ap(i)%xfp0(1:NL)   = glob(p)%xfp0(1:NL)
      ap(i)%yfp0(1:NL)   = glob(p)%yfp0(1:NL)
      ap(i)%zfp0(1:NL)   = glob(p)%zfp0(1:NL)
      ap(i)%area(1:NL)   = glob(p)%area(1:NL)

      do l=1,NL
 
        !neighbor 1
        if ( (ap(i)%xfp(l)+offset) .gt. rightbound ) then
          if ( ((ap(i)%yfp(l)+offset) .ge. frontbound) .and. ((ap(i)%yfp(l)-offset) .le. backbound) ) then 
             ap(i)%nb(1) = 1 !neighbor 1 is slave of particle ap(i)%mslv
          endif
        endif

      !neighbor 2
        if ( (ap(i)%xfp(l)+offset) .gt. rightbound ) then
          if  ( (ap(i)%yfp(l)-offset) .le. frontbound ) then 
            ap(i)%nb(2) = 1 !neighbor 2 is slave of particle ap(i)%mslv
          endif
        endif

      !neighbor 3
        if ( (ap(i)%yfp(l)-offset) .lt. frontbound ) then
          if ( ((ap(i)%xfp(l)+offset) .ge. leftbound) .and. ((ap(i)%xfp(l)-offset) .le. rightbound )) then 
             ap(i)%nb(3) = 1 !neighbor 3 is slave of particle ap(i)%mslv
          endif
        endif

      !neighbor 4
        if ( (ap(i)%yfp(l)-offset) .lt. frontbound ) then
          if  ( (ap(i)%xfp(l)-offset) .le. leftbound ) then 
            ap(i)%nb(4) = 1 !neighbor 4 is slave of particle ap(i)%mslv
          endif
        endif

      !neighbor 5
        if ( (ap(i)%xfp(l)-offset) .lt. leftbound ) then
          if ( ((ap(i)%yfp(l)+offset) .ge. frontbound) .and. ((ap(i)%yfp(l)-offset) .le. backbound) ) then 
             ap(i)%nb(5) = 1 !neighbor 5 is slave of particle ap(i)%mslv
          endif
        endif 

      !neighbor 6
        if ( (ap(i)%yfp(l)+offset) .gt. backbound ) then
          if  ( (ap(i)%xfp(l)-offset) .le. leftbound ) then 
            ap(i)%nb(6) = 1 !neighbor 6 is slave of particle ap(i)%mslv
          endif
        endif

      !neighbor 7
        if ( (ap(i)%yfp(l)+offset) .gt. backbound ) then
          if ( ((ap(i)%xfp(l)+offset) .ge. leftbound) .and. ((ap(i)%xfp(l)-offset) .le. rightbound )) then 
             ap(i)%nb(7) = 1 !neighbor 7 is slave of particle ap(i)%mslv
          endif
        endif

      !neighbor 8
        if ( (ap(i)%yfp(l)+offset) .gt. backbound ) then
          if  ( (ap(i)%xfp(l)+offset) .ge. rightbound ) then 
            ap(i)%nb(8) = 1 !neighbor 8 is slave of particle ap(i)%mslv
          endif
        endif
              
      enddo

    else


    count_slve_loc = 0
    !neighbor 1 of particle's master
    proccoords(1) = nint( dims(1)*glob(p)%x/lx - ax ) + 1
    proccoords(2) = nint( dims(2)*glob(p)%y/ly - ay ) 
    call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
    if (myid .eq. procrank) then

      do l=1,NL

        if ( (glob(p)%xfp(l)+offset) .gt. rightbound ) then
          if ( ((glob(p)%yfp(l)+offset) .ge. frontbound) .and. ((glob(p)%yfp(l)-offset) .le. backbound) ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            ap(i)%mslv = -p     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%nb(5) = 1      !neighbor 5 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            goto 201

          endif
        endif  
      enddo

201 continue  
    endif

    !neighbor 2 of particle's master
    proccoords(1) = nint( dims(1)*glob(p)%x/lx - ax ) + 1
    proccoords(2) = nint( dims(2)*glob(p)%y/ly - ay ) - 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
    if (myid .eq. procrank) then

      do l=1,NL

        if ( (glob(p)%xfp(l)+offset) .gt. rightbound ) then
          if  ( (glob(p)%yfp(l)-offset) .le. frontbound ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            ap(i)%mslv = -p     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%nb(6) = 1      !neighbor 6 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            goto 202

          endif
        endif  
      enddo

202 continue 
    endif



    !neighbor 3 of particle's master
    proccoords(1) = nint( dims(1)*glob(p)%x/lx - ax )
    proccoords(2) = nint( dims(2)*glob(p)%y/ly - ay ) - 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
    if (myid .eq. procrank) then


      do l=1,NL

        if ( (glob(p)%yfp(l)-offset) .lt. frontbound ) then
          if ( ((glob(p)%xfp(l)+offset) .ge. leftbound) .and. ((glob(p)%xfp(l)-offset) .le. rightbound) ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            ap(i)%mslv = -p     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%nb(7) = 1      !neighbor 7 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            goto 203

          endif
        endif  
      enddo

203 continue 
    endif


    !neighbor 4 of particle's master
    proccoords(1) = nint( dims(1)*glob(p)%x/lx - ax ) - 1
    proccoords(2) = nint( dims(2)*glob(p)%y/ly - ay ) - 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
    if (myid .eq. procrank) then

      do l=1,NL

        if ( (glob(p)%yfp(l)-offset) .lt. frontbound ) then
          if  ( (glob(p)%xfp(l)-offset) .le. leftbound ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            ap(i)%mslv = -p     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%nb(8) = 1      !neighbor 8 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            goto 204

          endif
        endif  
      enddo

204 continue 
    endif


    !neighbor 5 of particle's master
    proccoords(1) = nint( dims(1)*glob(p)%x/lx - ax ) - 1
    proccoords(2) = nint( dims(2)*glob(p)%y/ly - ay )
    call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
    if (myid .eq. procrank) then

      do l=1,NL

      if ( (glob(p)%xfp(l)-offset) .lt. leftbound ) then
        if ( ((glob(p)%yfp(l)+offset) .ge. frontbound) .and. ((glob(p)%yfp(l)-offset) .le. backbound ) ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            ap(i)%mslv = -p     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%nb(1) = 1      !neighbor 1 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            goto 205

          endif
        endif  
      enddo

205 continue 
    endif


    !neighbor 6 of particle's master
    proccoords(1) = nint( dims(1)*glob(p)%x/lx - ax ) - 1
    proccoords(2) = nint( dims(2)*glob(p)%y/ly - ay ) + 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
    if (myid .eq. procrank) then

      do l=1,NL

      if ( (glob(p)%yfp(l)+offset) .gt. backbound ) then
        if  ( (glob(p)%xfp(l)-offset) .le. leftbound ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            ap(i)%mslv = -p     !myid is slave of particle abs(ap(i)%mslv)
            ap(i)%nb(2) = 1      !neighbor 2 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            goto 206

          endif
        endif  
      enddo

206 continue 
    endif


    !neighbor 7 of particle's master
    proccoords(1) = nint( dims(1)*glob(p)%x/lx - ax )
    proccoords(2) = nint( dims(2)*glob(p)%y/ly - ay ) + 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
    if (myid .eq. procrank) then

      do l=1,NL

      if ( (glob(p)%yfp(l)+offset) .gt. backbound ) then
        if ( ((glob(p)%xfp(l)+offset) .ge. leftbound) .and. ((glob(p)%xfp(l)-offset) .le. rightbound) ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            ap(i)%mslv = -p     !myid is slave of particle p=ap(i)%mslv
            ap(i)%nb(3) = 1      !neighbor 3 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            goto 207

          endif
        endif  
      enddo

207 continue 
    endif


    !neighbor 8 of particle's master
    proccoords(1) = nint( dims(1)*glob(p)%x/lx - ax ) + 1
    proccoords(2) = nint( dims(2)*glob(p)%y/ly - ay ) + 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
    if (myid .eq. procrank) then

      do l=1,NL

      if ( (glob(p)%yfp(l)+offset) .gt. backbound ) then
        if  ( (glob(p)%xfp(l)+offset) .ge. rightbound ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            ap(i)%mslv = -p     !myid is slave of particle p=ap(i)%mslv
            ap(i)%nb(4) = 1      !neighbor 4 of myid is particle's master
            count_slve_loc = count_slve_loc + 1
            goto 208

          endif
        endif  
      enddo

208 continue 
    endif


  endif
enddo


  pmax = i 

  !if(pmax.eq.0) pmax = 1 ! for now
  print*, 'Thread ', myid, ' contains ', pmax, 'paricles.'

  call MPI_ALLREDUCE(npmstr,count_mstr_all,1,MPI_INTEGER,MPI_SUM,comm_cart,error)
  print*,count_mstr_all,npmstr,np
  if (count_mstr_all.ne.np) then
    write(6,*) 'Fatal error in initialisation of particle positions!'
    write(6,*) 'Program aborted...'
    call mpi_abort(comm_cart,error,error)
    stop
  elseif(pmax.gt.npmax) then
    write(6,*) 'Size of local particle array of process ', myid, ' is too small!'
    write(6,*) 'Program aborted... (later I will simply write a warning and allocate more memory)'
    call mpi_abort(comm_cart,error,error)
    stop
  else
    write(6,*) 'The particles were successfully initialized in thread ', myid, ' !'
  endif


endif ! in.eq.0
!
if (in.eq.1) then
!!
!! write particle related data directly in parallel with MPI-IO

  do i=1,np
    counter = 0
    do p=1,pmax
      if (i.eq.ap(p)%mslv) then
        counter = counter + 1
        glob(i)%x = ap(p)%x !WRITE THIS IN A MORE COMPACT WAY!
        glob(i)%y = ap(p)%y
        glob(i)%z = ap(p)%z
        glob(i)%integralv = ap(p)%integralv
        glob(i)%volume    = ap(p)%volume
        glob(i)%volume0   = ap(p)%volume0
        glob(i)%xfp(1:NL) = ap(p)%xfp(1:NL)
        glob(i)%yfp(1:NL) = ap(p)%yfp(1:NL)
        glob(i)%zfp(1:NL) = ap(p)%zfp(1:NL)
        glob(i)%ul(1:NL)  = ap(p)%ul(1:NL)
        glob(i)%vl(1:NL)  = ap(p)%vl(1:NL)
        glob(i)%wl(1:NL)  = ap(p)%wl(1:NL)
        glob(i)%fxl(1:NL)  = ap(p)%fxl(1:NL)
        glob(i)%fyl(1:NL)  = ap(p)%fyl(1:NL)
        glob(i)%fzl(1:NL)  = ap(p)%fzl(1:NL)
        glob(i)%dxdt(1:NL) = ap(p)%dxdt(1:NL)
        glob(i)%dydt(1:NL) = ap(p)%dydt(1:NL)
        glob(i)%dzdt(1:NL) = ap(p)%dzdt(1:NL)
        glob(i)%xfpo(1:NL) = ap(p)%xfpo(1:NL)
        glob(i)%yfpo(1:NL) = ap(p)%yfpo(1:NL)
        glob(i)%zfpo(1:NL) = ap(p)%zfpo(1:NL)
        glob(i)%xfpold(1:NL) = ap(p)%xfpold(1:NL)
        glob(i)%yfpold(1:NL) = ap(p)%yfpold(1:NL)
        glob(i)%zfpold(1:NL) = ap(p)%zfpold(1:NL)
        glob(i)%xfp0(1:NL) = ap(p)%xfp0(1:NL)
        glob(i)%yfp0(1:NL) = ap(p)%yfp0(1:NL)
        glob(i)%zfp0(1:NL) = ap(p)%zfp0(1:NL)
        glob(i)%area(1:NL) = ap(p)%area(1:NL)
      endif
    enddo
    if(counter.eq.0) then
        glob(i)%x = 0. !WRITE THIS IN A MORE COMPACT WAY!
        glob(i)%y = 0.
        glob(i)%z = 0.
        glob(i)%integralv = 0.
        glob(i)%volume  = 0.
        glob(i)%volume0 = 0.
        glob(i)%xfp(1:NL) = 0.
        glob(i)%yfp(1:NL) = 0.
        glob(i)%zfp(1:NL) = 0.
        glob(i)%ul(1:NL)  = 0.
        glob(i)%vl(1:NL)  = 0.
        glob(i)%wl(1:NL)  = 0.
        glob(i)%fxl(1:NL)  = 0.
        glob(i)%fyl(1:NL)  = 0.
        glob(i)%fzl(1:NL)  = 0.
        glob(i)%dxdt(1:NL) = 0.
        glob(i)%dydt(1:NL) = 0.
        glob(i)%dzdt(1:NL) = 0.
        glob(i)%xfpo(1:NL) = 0.
        glob(i)%yfpo(1:NL) = 0.
        glob(i)%zfpo(1:NL) = 0.
        glob(i)%xfpold(1:NL) = 0.
        glob(i)%yfpold(1:NL) = 0.
        glob(i)%zfpold(1:NL) = 0.
        glob(i)%xfp0(1:NL) = 0.
        glob(i)%yfp0(1:NL) = 0.
        glob(i)%zfp0(1:NL) = 0.
        glob(i)%area(1:NL) = 0.
    endif
  enddo
  call mpi_reduce(glob(1)%x,glob_all(1)%x,np*skip,mpi_real8,mpi_sum,0,MPI_COMM_WORLD,error)
  if (myid.eq.0) then ! myid = 0 writes the data into a single file
    open(20,file=datadir//'allpartdata'//istepchar,access='direct',recl=np*skip*lenr)
    write(20,rec=1) glob_all
    close(20)
  endif
endif

return
end subroutine loadpart
!
end module mod_loadpart
