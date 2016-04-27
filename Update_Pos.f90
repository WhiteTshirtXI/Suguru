module mod_Update_Pos
use mpi
use mod_common
use mod_common_mpi
use mod_interp_spread
use mod_param
implicit none
private
public Update_Pos
contains
!
subroutine Update_Pos
implicit none

integer:: p,g,i,j,cap,itr
real::  Cont_x,Cont_y,distx,disty,dtnew,dtold,length




type neighbour2
  real :: x,y,z
  real, dimension(nl) :: xfp,yfp,zfp
end type neighbour2
type(neighbour2), dimension(1:npmax,0:8) :: anb2


type pd
  real, dimension(1:NL) :: xfp,yfp,zfp
end type pd
type(pd), dimension(1:npmax) :: p_data

integer :: tag,l,k
integer :: nrrequests
integer :: arrayrequests(1:30)
integer :: arraystatuses(MPI_STATUS_SIZE,1:30)

real :: ax,ay

real :: leftbound,rightbound,frontbound,backbound
integer :: nbrec2,nbrec,nbsend
integer, dimension(ndims) :: proccoords
integer :: procrank
integer :: counter !delete after checking things

integer :: count_mstr,count_slve,count_mstr_all,count_slve_loc
logical :: found_mstr

integer, dimension(0:8) :: pmax_nb
integer, dimension(1:npmax,0:8) :: mslv_nb,newmaster_nb
integer :: idp,idp_nb,nb
real :: boundleftnb,boundrightnb,boundfrontnb,boundbacknb
integer,parameter :: nprocs = dims(1)*dims(2)




pmax_nb(0)=pmax
!$omp workshare
mslv_nb(1:pmax,0) = ap(1:pmax)%mslv
!$omp end workshare
do nb=1,8
  nbsend = 4+nb
  nbrec  = nb
  if (nbsend .gt. 8) nbsend = nbsend - 8
  call MPI_SENDRECV(pmax_nb(0),1,MPI_INTEGER,neighbor(nbsend),1, &
                    pmax_nb(nbrec),1,MPI_INTEGER,neighbor(nbrec),1, &
                    comm_cart,status,error)
  call MPI_SENDRECV(mslv_nb(1,0),pmax,MPI_INTEGER,neighbor(nbsend),2, &
                    mslv_nb(1,nbrec),pmax_nb(nbrec),MPI_INTEGER,neighbor(nbrec),2, &
                    comm_cart,status,error)
enddo




do p=1,pmax

  if (ap(p)%mslv.gt.0) then
!************************************periodicity for old values of xfp, yfp and zfp*********************************************
  do l=1,NL
  if ((ap(p)%xfpold(l)-ap(p)%xfp(l))>(lx*.5)) then
  ap(p)%xfpold(l)=ap(p)%xfpold(l)-lx
  elseif ((ap(p)%xfpold(l)-ap(p)%xfp(l))<-(lx*.5)) then
  ap(p)%xfpold(l)=ap(p)%xfpold(l)+lx
  endif
  if ((ap(p)%yfpold(l)-ap(p)%yfp(l))>(ly*.5)) then
  ap(p)%yfpold(l)=ap(p)%yfpold(l)-ly
  elseif ((ap(p)%yfpold(l)-ap(p)%yfp(l))<-(ly*.5)) then
  ap(p)%yfpold(l)=ap(p)%yfpold(l)+ly
  endif
  enddo
!*****************************************updating position of Lagrangian points***********************************************

    do l=1,NL
    ap(p)%xfp(l)=ap(p)%xfp(l)+(ap(p)%ul(l)*dt)
    ap(p)%yfp(l)=ap(p)%yfp(l)+(ap(p)%vl(l)*dt)
    ap(p)%zfp(l)=ap(p)%zfp(l)+(ap(p)%wl(l)*dt)
    enddo


     
  
!open(1369,file="force.dat")
!write(1369,*) 'VARIABLES = "y","z","forcey","forcez"'
!write(1369,*) 'ZONE T="Zone1"',' I=',jtot,' J=',ktot,', F=POINT'
!do j = 1,nnlon2
!do i = 1,nnlat2
!write(1369,'(6E15.6)') x_cap(i,j),y_cap(i,j),z_cap(i,j),x_n_cap(i,j),y_n_cap(i,j),z_n_cap(i,j)
!enddo
!enddo
!close(1369)

!**************************computing area for each element, total area and volume of capsule*******************************
 !call surf_vol

!*********************************seeing whether capsule crosses each wall to fix it *************************************
    do l=1,NL
    if (ap(p)%zfp(l)<0.) then
    print*,"particle crossed lower wall"
    ap(p)%zfp(l)=0.
    else if (ap(p)%zfp(l)>lz) then
    print*,"particle crossed upper wall"
    ap(p)%zfp(l)=lz
    endif
    enddo
!********************************************end of solving procedure**************************************************
   endif
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!checking periodicity!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do p=1,pmax
  if (ap(p)%mslv .gt. 0) then

    ap(p)%x=0.
    ap(p)%y=0.
    ap(p)%z=0.
    ulmean=0.
    vlmean=0.
    wlmean=0.
    do l=1,NL
    ap(p)%x = ap(p)%x+ap(p)%xfp(l)
    ap(p)%y = ap(p)%y+ap(p)%yfp(l)
    ap(p)%z = ap(p)%z+ap(p)%zfp(l)
    ulmean=ulmean+ap(p)%ul(l)
    vlmean=vlmean+ap(p)%vl(l)
    wlmean=wlmean+ap(p)%wl(l)
    enddo
    ap(p)%x = ap(p)%x/NL
    ap(p)%y = ap(p)%y/NL
    ap(p)%z = ap(p)%z/NL
    ulmean=ulmean/NL
    vlmean=vlmean/NL
    wlmean=wlmean/NL

    
    if (ap(p)%x .gt. lx) then
      ap(p)%x = ap(p)%x - lx
      do l=1,NL
        ap(p)%xfp(l) = ap(p)%xfp(l) - lx 
      enddo
    endif

    if (ap(p)%x .lt. 0.) then
      ap(p)%x = ap(p)%x + lx
      do l=1,NL
        ap(p)%xfp(l) = ap(p)%xfp(l) + lx 
      enddo
    endif

    if (ap(p)%y .gt. ly) then
      ap(p)%y = ap(p)%y - ly
      do l=1,NL
        ap(p)%yfp(l) = ap(p)%yfp(l) - ly 
      enddo
    endif

    if (ap(p)%y .lt. 0.) then
      ap(p)%y = ap(p)%y + ly
      do l=1,NL
        ap(p)%yfp(l) = ap(p)%yfp(l) + ly 
      enddo
    endif

  endif
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! particle positions were updated. Now check if there are new masters
!
do p=1,pmax
  newmaster_nb(p,0) = 0
  if (ap(p)%mslv.gt.0) then
    ! myid was master of particle ap(p)%mslv at previous time step n
    ax = 0.5
    ay = 0.5
    if (ap(p)%x.eq.lx) ax = 0.51
    if (ap(p)%x.eq.0) ax = 0.49
    if (ap(p)%y.eq.ly) ay = 0.51
    if (ap(p)%y.eq.0) ay = 0.49
    proccoords(1) = nint( dims(1)*ap(p)%x/lx - ax )
    proccoords(2) = nint( dims(2)*ap(p)%y/ly - ay ) 
    call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
    if (procrank .ne. myid) then
      ! particle ap(p)%mslv has a new master at time step n+1
      do nb=1,8
        if (procrank .eq. neighbor(nb)) then
          newmaster_nb(p,0) = nb
        endif
      enddo
    endif
  endif
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MPI_integer should be changed
!
! exchange data
!
do nb=1,8
  nbsend = nb
  nbrec  = nb+4
  if (nbrec .gt. 8) nbrec = nbrec-8
  call MPI_SENDRECV(newmaster_nb(1,0),pmax,MPI_INTEGER,neighbor(nbsend),1, &
                    newmaster_nb(1,nbrec),pmax_nb(nbrec),MPI_INTEGER,neighbor(nbrec),1, &
                    comm_cart,status,error)
enddo
!
do p=1,pmax
  nrrequests = 0
  if (newmaster_nb(p,0).gt.0) then
    nbsend = newmaster_nb(p,0)
    tag = ap(p)%mslv*10+nbsend
    nrrequests = nrrequests + 1
    call MPI_ISEND(ap(p)%x,send_real,MPI_REAL8,neighbor(nbsend),tag,comm_cart,arrayrequests((nrrequests-1)*2+1),error)
    ap(p)%mslv = -ap(p)%mslv ! master is a slave now
  endif
  if(mslv_nb(p,0).lt.0) then
    idp = -ap(p)%mslv
    do nbrec = 1,8
      nbrec2 = nbrec + 4
      if(nbrec2 .gt. 8) nbrec2 = nbrec2-8
      do i=1,pmax_nb(nbrec)
        idp_nb = mslv_nb(i,nbrec)
        if(newmaster_nb(i,nbrec) .eq. nbrec2.and.idp.eq.idp_nb) then
          ap(p)%mslv = -ap(p)%mslv ! slave became a master
          nrrequests = nrrequests + 1
          tag = ap(p)%mslv*10+nbrec2
          call MPI_IRECV(ap(p)%x,send_real,MPI_REAL8,neighbor(nbrec),tag,comm_cart,arrayrequests((nrrequests-1)*2+1),error)
        endif
      enddo
    enddo
  endif
  nrrequests = nrrequests
  call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,error)
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!
! Masters are known now: ap(p)%mslv > 0.  
! Next step: determine slaves and master/slave neighbors
!
forall(p=1:npmax,k=0:8) mslv_nb(p,k) = 0
mslv_nb(1:pmax,0) = ap(1:pmax)%mslv  ! new masters, but new slaves not all determined yet!
                                     ! process might remain master, but slave might lose particle
anb2(1:pmax,0)%x = ap(1:pmax)%x ! from new master
anb2(1:pmax,0)%y = ap(1:pmax)%y ! from new master
anb2(1:pmax,0)%z = ap(1:pmax)%z ! from new master

do l=1,NL
  anb2(1:pmax,0)%xfp(l) = ap(1:pmax)%xfp(l) ! from new master
  anb2(1:pmax,0)%yfp(l) = ap(1:pmax)%yfp(l) ! from new master
  anb2(1:pmax,0)%zfp(l) = ap(1:pmax)%zfp(l) ! from new master
enddo

do nb=1,8
  nbsend = 4+nb
  if (nbsend .gt. 8) nbsend = nbsend - 8
  nbrec  = nb
  call MPI_SENDRECV(mslv_nb(1,0),pmax,MPI_INTEGER,neighbor(nbsend),1, &
                    mslv_nb(1,nbrec),pmax_nb(nbrec),MPI_INTEGER,neighbor(nbrec),1, &
                    comm_cart,status,error)
  call MPI_SENDRECV(anb2(1,0)%x,pmax*3*(1+NL),MPI_REAL8,neighbor(nbsend),2, &
                    anb2(1,nbrec)%x,pmax_nb(nbrec)*3*(1+NL),MPI_REAL8,neighbor(nbrec),2, &
                    comm_cart,status,error)
  ! send x,y,z -> 3*pmax contiguous info
  ! (see definition of type neighbor2 in the begining of the subroutine)
enddo
!
! recompute particle positions because of periodic b.c.'s in the x and y-direction
!
nb=1
do p=1,pmax_nb(nb)
  if (mslv_nb(p,nb) .gt. 0) then
    boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left  boundary of neighbor nb
    if (anb2(p,nb)%x .lt. boundleftnb) then
      anb2(p,nb)%x = anb2(p,nb)%x + lx
      do l=1,NL
        anb2(p,nb)%xfp(l) = anb2(p,nb)%xfp(l) + lx 
      enddo
    endif
  endif
enddo
nb=2
do p=1,pmax_nb(nb)
  if (mslv_nb(p,nb) .gt. 0) then
    boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left  boundary of neighbor nb
    boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back  boundary of neighbor nb
    if (anb2(p,nb)%x .lt. boundleftnb) then
      anb2(p,nb)%x = anb2(p,nb)%x + lx
      do l=1,NL
        anb2(p,nb)%xfp(l) = anb2(p,nb)%xfp(l) + lx 
      enddo
    endif
    if (anb2(p,nb)%y .gt. boundbacknb) then
      anb2(p,nb)%y = anb2(p,nb)%y - ly
      do l=1,NL
        anb2(p,nb)%yfp(l) = anb2(p,nb)%yfp(l) - ly 
      enddo
    endif
  endif
enddo
nb=3
do p=1,pmax_nb(nb)
  if (mslv_nb(p,nb) .gt. 0) then
    boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back  boundary of neighbor nb
    if (anb2(p,nb)%y .gt. boundbacknb) then
      anb2(p,nb)%y = anb2(p,nb)%y - ly
      do l=1,NL
        anb2(p,nb)%yfp(l) = anb2(p,nb)%yfp(l) - ly 
      enddo
    endif
  endif
enddo
nb=4
do p=1,pmax_nb(nb)
  if (mslv_nb(p,nb) .gt. 0) then
    boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
    boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back  boundary of neighbor nb
    if (anb2(p,nb)%x .gt. boundrightnb) then
      anb2(p,nb)%x = anb2(p,nb)%x - lx
      do l=1,NL
        anb2(p,nb)%xfp(l) = anb2(p,nb)%xfp(l) - lx 
      enddo
    endif
    if (anb2(p,nb)%y .gt. boundbacknb) then
      anb2(p,nb)%y = anb2(p,nb)%y - ly
      do l=1,NL
        anb2(p,nb)%yfp(l) = anb2(p,nb)%yfp(l) - ly 
      enddo
    endif
  endif
enddo
nb=5
do p=1,pmax_nb(nb)
  if (mslv_nb(p,nb) .gt. 0) then
    boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
    if (anb2(p,nb)%x .gt. boundrightnb) then
      anb2(p,nb)%x = anb2(p,nb)%x - lx
      do l=1,NL
        anb2(p,nb)%xfp(l) = anb2(p,nb)%xfp(l) - lx 
      enddo
    endif
  endif
enddo
nb=6
do p=1,pmax_nb(nb)
  if (mslv_nb(p,nb) .gt. 0) then
    boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
    boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
    if (anb2(p,nb)%x .gt. boundrightnb) then
      anb2(p,nb)%x = anb2(p,nb)%x - lx
      do l=1,NL
        anb2(p,nb)%xfp(l) = anb2(p,nb)%xfp(l) - lx 
      enddo
    endif

    if (anb2(p,nb)%y .lt. boundfrontnb) then
      anb2(p,nb)%y = anb2(p,nb)%y + ly
      do l=1,NL
        anb2(p,nb)%yfp(l) = anb2(p,nb)%yfp(l) + ly 
      enddo
    endif
  endif
enddo
nb=7
do p=1,pmax_nb(nb)
  if (mslv_nb(p,nb) .gt. 0) then
    boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
    if (anb2(p,nb)%y .lt. boundfrontnb) then
      anb2(p,nb)%y = anb2(p,nb)%y + ly
      do l=1,NL
        anb2(p,nb)%yfp(l) = anb2(p,nb)%yfp(l) + ly 
      enddo
    endif
  endif
enddo
nb=8
do p=1,pmax_nb(nb)
  if (mslv_nb(p,nb) .gt. 0) then
    boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left  boundary of neighbor nb
    boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
    if (anb2(p,nb)%x .lt. boundleftnb) then
      anb2(p,nb)%x = anb2(p,nb)%x + lx
      do l=1,NL
        anb2(p,nb)%xfp(l) = anb2(p,nb)%xfp(l) + lx 
      enddo
    endif

    if (anb2(p,nb)%y .lt. boundfrontnb) then
      anb2(p,nb)%y = anb2(p,nb)%y + ly
      do l=1,NL
        anb2(p,nb)%yfp(l) = anb2(p,nb)%yfp(l) + ly 
      enddo
    endif
  endif
enddo

!
! important info that should be kept:
!
sp(1:pmax)%mslv      = ap(1:pmax)%mslv
sp(1:pmax)%x         = ap(1:pmax)%x
sp(1:pmax)%y         = ap(1:pmax)%y
sp(1:pmax)%z         = ap(1:pmax)%z
sp(1:pmax)%integralv = ap(1:pmax)%integralv
sp(1:pmax)%volume    = ap(1:pmax)%volume
sp(1:pmax)%volume0   = ap(1:pmax)%volume0

forall(l=1:NL)
  sp(1:pmax)%xfp(l)  = ap(1:pmax)%xfp(l)
  sp(1:pmax)%yfp(l)  = ap(1:pmax)%yfp(l)
  sp(1:pmax)%zfp(l)  = ap(1:pmax)%zfp(l)
  sp(1:pmax)%ul(l)   = ap(1:pmax)%ul(l)
  sp(1:pmax)%vl(l)   = ap(1:pmax)%vl(l)
  sp(1:pmax)%wl(l)   = ap(1:pmax)%wl(l)
  sp(1:pmax)%fxl(l)  = ap(1:pmax)%fxl(l) 
  sp(1:pmax)%fyl(l)  = ap(1:pmax)%fyl(l)
  sp(1:pmax)%fzl(l)  = ap(1:pmax)%fzl(l) 
  sp(1:pmax)%dxdt(l) = ap(1:pmax)%dxdt(l)
  sp(1:pmax)%dydt(l) = ap(1:pmax)%dydt(l)
  sp(1:pmax)%dzdt(l) = ap(1:pmax)%dzdt(l)
  sp(1:pmax)%xfpo(l) = ap(1:pmax)%xfpo(l)
  sp(1:pmax)%yfpo(l) = ap(1:pmax)%yfpo(l)
  sp(1:pmax)%zfpo(l) = ap(1:pmax)%zfpo(l)
  sp(1:pmax)%xfpold(l) = ap(1:pmax)%xfpold(l)
  sp(1:pmax)%yfpold(l) = ap(1:pmax)%yfpold(l)
  sp(1:pmax)%zfpold(l) = ap(1:pmax)%zfpold(l)
  sp(1:pmax)%xfp0(l) = ap(1:pmax)%xfp0(l)
  sp(1:pmax)%yfp0(l) = ap(1:pmax)%yfp0(l)
  sp(1:pmax)%zfp0(l) = ap(1:pmax)%zfp0(l)
  sp(1:pmax)%xn(l) = ap(1:pmax)%xn(l)
  sp(1:pmax)%yn(l) = ap(1:pmax)%yn(l)
  sp(1:pmax)%zn(l) = ap(1:pmax)%zn(l)
  sp(1:pmax)%area(l)   = ap(1:pmax)%area(l)
end forall

!
! clear structure ap for re-ordering:
!
ap(1:pmax)%mslv = 0
ap(1:pmax)%x = 0.
ap(1:pmax)%y = 0.
ap(1:pmax)%z = 0.
ap(1:pmax)%integralv = 0.
ap(1:pmax)%volume=0.
ap(1:pmax)%volume0=0.

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
  ap(p)%xfpo(:) = 0.
  ap(p)%yfpo(:) = 0.
  ap(p)%zfpo(:) = 0.
  ap(p)%xfpold(:) = 0.
  ap(p)%yfpold(:) = 0.
  ap(p)%zfpold(:) = 0.
  ap(p)%area(:) = 0.
  ap(p)%xfp0(:) = 0.
  ap(p)%yfp0(:) = 0.
  ap(p)%zfp0(:) = 0.
  ap(p)%xn(:) = 0.
  ap(p)%yn(:) = 0.
  ap(p)%zn(:) = 0.
  ap(p)%nb(:)   = 0
end forall
!
leftbound   = (coords(1)  )*lx/(1.*dims(1)) ! left  boundary of process myid
rightbound  = (coords(1)+1)*lx/(1.*dims(1)) ! right boundary of process myid
frontbound  = (coords(2)  )*ly/(1.*dims(2)) ! front boundary of process myid
backbound   = (coords(2)+1)*ly/(1.*dims(2)) ! back  boundary of process myid
!
i = 0
count_mstr = 0
count_slve = 0
do idp = 1,np
  found_mstr = .false.
  call binsearch(idp,mslv_nb(1:pmax_nb(0),0),pmax_nb(0),found_mstr,p)
!  do p=1,pmax
!    if(sp(p)%mslv.eq.idp) then
!      found_mstr = .true.
    if(found_mstr) then
      i = i + 1
      ap(i)%mslv      = sp(p)%mslv
      ap(i)%x         = sp(p)%x
      ap(i)%y         = sp(p)%y
      ap(i)%z         = sp(p)%z
      ap(i)%integralv = sp(p)%integralv  
      ap(i)%volume    = sp(p)%volume
      ap(i)%volume0   = sp(p)%volume0
      do l=1,NL
        ap(i)%xfp(l)  = sp(p)%xfp(l)
        ap(i)%yfp(l)  = sp(p)%yfp(l)
        ap(i)%zfp(l)  = sp(p)%zfp(l)
        ap(i)%ul(l)   = sp(p)%ul(l)
        ap(i)%vl(l)   = sp(p)%vl(l)
        ap(i)%wl(l)   = sp(p)%wl(l)
        ap(i)%fxl(l)  = sp(p)%fxl(l)
        ap(i)%fyl(l)  = sp(p)%fyl(l)
        ap(i)%fzl(l)  = sp(p)%fzl(l)
        ap(i)%dxdt(l) = sp(p)%dxdt(l)
        ap(i)%dydt(l) = sp(p)%dydt(l)
        ap(i)%dzdt(l) = sp(p)%dzdt(l)
        ap(i)%xfpo(l) = sp(p)%xfpo(l)
        ap(i)%yfpo(l) = sp(p)%yfpo(l)
        ap(i)%zfpo(l) = sp(p)%zfpo(l)
        ap(i)%xfpold(l) = sp(p)%xfpold(l)
        ap(i)%yfpold(l) = sp(p)%yfpold(l)
        ap(i)%zfpold(l) = sp(p)%zfpold(l)
        ap(i)%area(l)   = sp(p)%area(l)
        ap(i)%xfp0(l) = sp(p)%xfp0(l)
        ap(i)%yfp0(l) = sp(p)%yfp0(l)
        ap(i)%zfp0(l) = sp(p)%zfp0(l)
        ap(i)%xn(l) = sp(p)%xn(l)
        ap(i)%yn(l) = sp(p)%yn(l)
        ap(i)%zn(l) = sp(p)%zn(l)
      enddo

      count_mstr = count_mstr + 1

      do l=1,NL


      !neighbor 1
        if ( (sp(i)%xfp(l)+offset) .gt. rightbound ) then
          if ( ((sp(i)%yfp(l)+offset) .ge. frontbound) .and. ((sp(i)%yfp(l)-offset) .le. backbound) ) then 
             ap(i)%nb(1) = 1 !neighbor 1 is slave of particle ap(i)%mslv
          endif
        endif

      !neighbor 2
        if ( (sp(i)%xfp(l)+offset) .gt. rightbound ) then
          if  ( (sp(i)%yfp(l)-offset) .le. frontbound ) then 
            ap(i)%nb(2) = 1 !neighbor 2 is slave of particle ap(i)%mslv
          endif
        endif

      !neighbor 3
        if ( (sp(i)%yfp(l)-offset) .lt. frontbound ) then
          if ( ((sp(i)%xfp(l)+offset) .ge. leftbound) .and. ((sp(i)%xfp(l)-offset) .le. rightbound )) then 
             ap(i)%nb(3) = 1 !neighbor 3 is slave of particle ap(i)%mslv
          endif
        endif

      !neighbor 4
        if ( (sp(i)%yfp(l)-offset) .lt. frontbound ) then
          if  ( (sp(i)%xfp(l)-offset) .le. leftbound ) then 
            ap(i)%nb(4) = 1 !neighbor 4 is slave of particle ap(i)%mslv
          endif
        endif

      !neighbor 5
        if ( (sp(i)%xfp(l)-offset) .lt. leftbound ) then
          if ( ((sp(i)%yfp(l)+offset) .ge. frontbound) .and. ((sp(i)%yfp(l)-offset) .le. backbound) ) then 
             ap(i)%nb(5) = 1 !neighbor 5 is slave of particle ap(i)%mslv
          endif
        endif 

      !neighbor 6
        if ( (sp(i)%yfp(l)+offset) .gt. backbound ) then
          if  ( (sp(i)%xfp(l)-offset) .le. leftbound ) then 
            ap(i)%nb(6) = 1 !neighbor 6 is slave of particle ap(i)%mslv
          endif
        endif

      !neighbor 7
        if ( (sp(i)%yfp(l)+offset) .gt. backbound ) then
          if ( ((sp(i)%xfp(l)+offset) .ge. leftbound) .and. ((sp(i)%xfp(l)-offset) .le. rightbound )) then 
             ap(i)%nb(7) = 1 !neighbor 7 is slave of particle ap(i)%mslv
          endif
        endif

      !neighbor 8
        if ( (sp(i)%yfp(l)+offset) .gt. backbound ) then
          if  ( (sp(i)%xfp(l)+offset) .ge. rightbound ) then 
            ap(i)%nb(8) = 1 !neighbor 8 is slave of particle ap(i)%mslv
          endif
        endif
              
      enddo
 
    else

      count_slve_loc = 0

      nb=5
      if ( (pmax_nb(nb) .ne. 0) ) then

        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then
  
          do l=1,NL
  
            if ( (anb2(k,nb)%xfp(l)+offset) .gt. leftbound ) then
              if ( ((anb2(k,nb)%yfp(l)+offset) .ge. frontbound) .and. ((anb2(k,nb)%yfp(l)-offset) .le. backbound ) ) then 
  
                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1

!                if (myid .eq. 3) write(6,*) 'hey myid3',nb!,anb2(k,nb)%xfp(l)-offset,leftbound

                goto 206  
                          
              endif
            endif

          enddo 
        endif

206 continue
      endif


      nb=6
      if ( (pmax_nb(nb) .ne. 0) ) then

        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then

          do l=1,NL

            if ( (anb2(k,nb)%yfp(l)-offset) .lt. backbound ) then
              if  ( (anb2(k,nb)%xfp(l)+offset) .ge. leftbound ) then 

                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1

!                if (myid .eq. 3) write(6,*) 'hey myid3',nb!,anb2(k,nb)%xfp(l)-offset,leftbound

                goto 207 
           
              endif
            endif

          enddo 
        endif

207 continue
      endif


      nb=7
      if ( (pmax_nb(nb) .ne. 0) ) then

        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then

          do l=1,NL

            if ( (anb2(k,nb)%yfp(l)-offset) .lt. backbound ) then
              if ( ((anb2(k,nb)%xfp(l)+offset) .ge. leftbound) .and. ((anb2(k,nb)%xfp(l)-offset) .le. rightbound) ) then 

                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1

!                if (myid .eq. 3) write(6,*) 'hey myid3',nb!,anb2(k,nb)%xfp(l)-offset,leftbound

                goto 208 
           
              endif
            endif

          enddo 
        endif

208 continue    
      endif


      nb=8
      if ( (pmax_nb(nb) .ne. 0) ) then

        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then

          do l=1,NL

            if ( (anb2(k,nb)%yfp(l)-offset) .lt. backbound ) then
              if  ( (anb2(k,nb)%xfp(l)-offset) .le. rightbound ) then 

                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1

!                if (myid .eq. 3) write(6,*) 'hey myid3',nb!,anb2(k,nb)%xfp(l)-offset,leftbound

                goto 201 
           
              endif
            endif

          enddo 
        endif

201 continue  
      endif


      nb=1
      if ( (pmax_nb(nb) .ne. 0) ) then

        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then

          do l=1,NL

            if ( (anb2(k,nb)%xfp(l)-offset) .lt. rightbound ) then
              if ( ((anb2(k,nb)%yfp(l)+offset) .ge. frontbound) .and. ((anb2(k,nb)%yfp(l)-offset) .le. backbound) ) then 

                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1

!                if (myid .eq. 3) write(6,*) 'hey myid3',nb,anb2(k,nb)%xfp(l)-offset,rightbound

                goto 202 
           
              endif
            endif

          enddo 
        endif

202 continue 
      endif


      nb=2
      if ( (pmax_nb(nb) .ne. 0) ) then

        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then

          do l=1,NL

            if ( (anb2(k,nb)%xfp(l)-offset) .lt. rightbound ) then
              if  ( (anb2(k,nb)%yfp(l)+offset) .ge. frontbound ) then 

                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1

!                if (myid .eq. 3) write(6,*) 'hey myid3',nb!,anb2(k,nb)%xfp(l)-offset,leftbound

                goto 203 
           
              endif
            endif

          enddo 
        endif

203 continue 
      endif


      nb=3
      if ( (pmax_nb(nb) .ne. 0) ) then

        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then

          do l=1,NL

            if ( (anb2(k,nb)%yfp(l)+offset) .gt. frontbound ) then
              if ( ((anb2(k,nb)%xfp(l)+offset) .ge. leftbound) .and. ((anb2(k,nb)%xfp(l)-offset) .le. rightbound) ) then 

                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1

!                if (myid .eq. 3) write(6,*) 'hey myid3',nb!,anb2(k,nb)%xfp(l)-offset,leftbound


                goto 204 
           
              endif
            endif

          enddo 
        endif

204 continue 
      endif


      nb=4
      if ( (pmax_nb(nb) .ne. 0) ) then

        call binsearch(idp,mslv_nb(1:pmax_nb(nb),nb),pmax_nb(nb),found_mstr,k)
        if(found_mstr) then

          do l=1,NL

            if ( (anb2(k,nb)%yfp(l)+offset) .gt. frontbound ) then
              if  ( (anb2(k,nb)%xfp(l)+offset) .ge. leftbound ) then 

                if(count_slve_loc.eq. 0 ) i = i+1
                ap(i)%mslv  = -mslv_nb(k,nb) ! myid is slave of particle mslv_nb(k,nb)
                ap(i)%nb(nb) = 1             ! neighbor nb of myid is particle's master
                count_slve_loc = count_slve_loc + 1

!                if (myid .eq. 3) write(6,*) 'hey myid3',nb!,anb2(k,nb)%xfp(l)-offset,leftbound


                goto 205 
            
              endif
            endif

          enddo 
        endif

205 continue 
      endif

      if(count_slve_loc.ne.0) count_slve = count_slve + 1

    endif
enddo
!
! the new value of pmax is:
!
pmax = count_mstr + count_slve
npmstr = count_mstr
!
!check if number of masters yield np
!
call MPI_ALLREDUCE(count_mstr,count_mstr_all,1,MPI_INTEGER,MPI_SUM,comm_cart,error)




do p=1,pmax
  if (ap(p)%mslv .le. 0) then
    do l=1,NL
      ap(p)%xfp(l)  = 0.
      ap(p)%yfp(l)  = 0.
      ap(p)%zfp(l)  = 0.
      ap(p)%ul(l)   = 0.
      ap(p)%vl(l)   = 0.
      ap(p)%wl(l)   = 0.
      ap(p)%fxl(l)  = 0.
      ap(p)%fyl(l)  = 0.
      ap(p)%fzl(l)  = 0.
      ap(p)%dxdt(l) = 0.
      ap(p)%dydt(l) = 0.
      ap(p)%dzdt(l) = 0.
      ap(p)%xfpo(l)  = 0.
      ap(p)%yfpo(l)  = 0.
      ap(p)%zfpo(l)  = 0.
      ap(p)%xfpold(l)  = 0.
      ap(p)%yfpold(l)  = 0.
      ap(p)%zfpold(l)  = 0.
      ap(p)%area(l)    = 0.
      ap(p)%xfp0(l)  = 0.
      ap(p)%yfp0(l)  = 0.
      ap(p)%zfp0(l)  = 0.
      ap(p)%xn(l)  = 0.
      ap(p)%yn(l)  = 0.
      ap(p)%zn(l)  = 0.
    enddo
  endif
enddo



write(6,'(A7,I5,A8,I5,A18,I5,A11,A8,I5)') 'Thread ', myid, ' masters ', count_mstr, &
' and is slave for ', count_slve, ' particles. ', ' pmax = ', pmax

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





return
end subroutine Update_Pos


subroutine binsearch(ival,array,idim,found,index)
integer, intent(in), dimension(1:) :: array
integer, intent(in) :: ival,idim
logical, intent(out) :: found
integer, intent(out) :: index
integer :: start,finish,range,mid

start = 1
finish = idim
range = finish-start
mid = (start+finish)/2
do while( abs(array(mid)) .ne. ival .and. range .gt.  0)
  if (ival .gt. abs(array(mid))) then
    start = mid + 1
  else
    finish = mid - 1
  endif
  range = finish - start
  mid = (start + finish)/2
enddo

if(array(mid).ne.ival) then
  found = .false.
  index = 0
else
  found = .true.
  index = mid
endif

return
end subroutine binsearch


end module mod_Update_Pos
