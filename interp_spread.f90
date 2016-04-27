module mod_interp_spread
use mod_param
use mod_common
use mod_common_mpi
use mod_kernel
use mod_bound
implicit none
private
public eulr2lagr,lagr2eulr,lagr2eulrind
contains

subroutine eulr2lagr
implicit none

integer :: i,j,k,l,p
integer :: ilow,ihigh,jlow,jhigh,klow,khigh
real :: coorx,coory,coorz
real :: coorxs,coorys,coorzs
real :: coorxfp,cooryfp,coorzfp
real :: kernelx,kernely,kernelz
real :: kernelxs,kernelys,kernelzs
type pneighbor
  real :: x,y,z
  real, dimension(1:nl) :: xfp,yfp,zfp,ul,vl,wl
!  real, dimension(1:nl) :: sumu 
end type pneighbor
type(pneighbor), dimension(0:8,1:npmax) :: anb
integer :: nb,nbsend,nbrecv
integer :: nrrequests
integer :: arrayrequests(1:3) !3=3*1 (master might have 3 slaves)
integer :: arraystatuses(MPI_STATUS_SIZE,1:3)
real :: boundleftnb,boundrightnb,boundfrontnb,boundbacknb
!real, dimension(1:nl,1:npmax) :: sumu
!character*3 :: partnr,rankpr
integer :: idp,tag
real :: auxu,auxv,auxw
real :: rcoeff
!
! first step: slaves need from master x,..,xfp,.. etc
!
!$omp workshare
anb(0,1:pmax)%x = ap(1:pmax)%x
anb(0,1:pmax)%y = ap(1:pmax)%y
anb(0,1:pmax)%z = ap(1:pmax)%z
anb(1:8,1:pmax)%x = 0.
anb(1:8,1:pmax)%y = 0.
anb(1:8,1:pmax)%z = 0.
forall(l=1:nl,p=1:pmax)
  anb(0,p)%xfp(l) = ap(p)%xfp(l)
  anb(0,p)%yfp(l) = ap(p)%yfp(l)
  anb(0,p)%zfp(l) = ap(p)%zfp(l)
end forall
forall(l=1:nl,i=1:8,p=1:pmax)
  anb(i,p)%xfp(l) = 0.
  anb(i,p)%yfp(l) = 0.
  anb(i,p)%zfp(l) = 0.
end forall
!$omp end workshare
!
do p=1,pmax
  nrrequests = 0
  do nb=1,8
    nbsend = nb    ! rank of process which receives data ('sends data to neighbor nbsend')
    nbrecv = nb+4  ! rank of process which sends data ('receives data from neighbor nbrecv')
    if (nbrecv .gt. 8) nbrecv = nbrecv - 8
    if (ap(p)%mslv .gt. 0) then
      ! myid is master of particle ap(p)%mslv
      if (ap(p)%nb(nbsend) .eq. 1) then
        ! neighbor(nbsend) is rank of slave for particle ap(p)%mslv
        if ( neighbor(nbsend) .eq. myid ) then
          ! process might be both master and slave of same particle due to periodic b.c.'s
          anb(nbrecv,p)%x = anb(0,p)%x
          anb(nbrecv,p)%y = anb(0,p)%y
          anb(nbrecv,p)%z = anb(0,p)%z
!$omp parallel default(shared) private(l)
          do l=1,nl
            anb(nbrecv,p)%xfp(l) = anb(0,p)%xfp(l)
            anb(nbrecv,p)%yfp(l) = anb(0,p)%yfp(l)
            anb(nbrecv,p)%zfp(l) = anb(0,p)%zfp(l)
          enddo
!$omp barrier
          if (nbrecv .eq. 1) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
!$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
            enddo
          endif
          if (nbrecv .eq. 2) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
            anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
!$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
            enddo
          endif
          if (nbrecv .eq. 3) then
            anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
!$omp do
            do l=1,nl
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
            enddo
          endif
          if (nbrecv .eq. 4) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
            anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
!$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
            enddo
          endif
          if (nbrecv .eq. 5) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
!$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
            enddo
          endif
          if (nbrecv .eq. 6) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
            anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
!$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
            enddo
          endif
          if (nbrecv .eq. 7) then
            anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
!$omp do
            do l=1,nl
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
            enddo
          endif
          if (nbrecv .eq. 8) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
            anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
!$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
            enddo
          endif
!$omp end parallel
        else
          idp = ap(p)%mslv
          tag = idp*10+nbsend
          nrrequests = nrrequests + 1
          call MPI_ISEND(anb(0,p)%x,3*(1+nl),MPI_REAL8,neighbor(nbsend), &
                         tag,comm_cart,arrayrequests((nrrequests-1)+1),error)
          ! send x,y,z,xfp(:),yfp(:),zfp(:) -> 3*(1+nl) contiguous info
          ! (see definition of type pneighbor in the begining of the subroutine)
        endif
      endif
    endif
    if (ap(p)%mslv .lt. 0) then
      ! myid is slave of particle -ap(p)%mslv
      if (ap(p)%nb(nbrecv) .eq. 1) then
        ! ap(p)%nb(nbrecv) is rank of master of particle -ap(p)%mslv
        idp = -ap(p)%mslv
        tag = idp*10+nbsend
        nrrequests = nrrequests + 1
          call MPI_IRECV(anb(nbrecv,p)%x,3*(1+nl),MPI_REAL8,neighbor(nbrecv), &
                         tag,comm_cart,arrayrequests((nrrequests-1)+1),error)
          ! recv x,y,z,xfp(:),yfp(:),zfp(:) -> 3*(1+nl) contiguous info
          ! (see definition of type pneighbor in the begining of the subroutine)
      endif
    endif
  enddo !do nb=
  call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,error)
enddo
!
! second step: recompute particle positions for slaves due to periodic b.c.'s.
! required: (part of) particle within domain bounds of slave process.
!
!$omp parallel default(shared) &
!$omp& private(p,nbrecv,boundleftnb,boundbacknb,boundrightnb,boundfrontnb)
!!$omp do schedule(dynamic)
!$omp do 
do p=1,pmax
  if (ap(p)%mslv .lt. 0) then
    ! myid is slave of particle -ap(p)%mslv
    nbrecv=1
    if (ap(p)%nb(nbrecv) .eq. 1) then
      ! ap(p)%nb(nbrecv) is rank of master of particle -ap(p)%mslv
      boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left boundary of neighbor nb
      if (anb(nbrecv,p)%x .lt. boundleftnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
        enddo
      endif
    endif
    nbrecv=2
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left boundary of neighbor nb
      boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back boundary of neighbor nb
      if (anb(nbrecv,p)%x .lt. boundleftnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
        enddo
      endif
      if (anb(nbrecv,p)%y .gt. boundbacknb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
        enddo
      endif
    endif
    nbrecv=3
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back boundary of neighbor nb
      if (anb(nbrecv,p)%y .gt. boundbacknb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
        enddo
      endif
    endif
    nbrecv=4
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
      boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back  boundary of neighbor nb
      if (anb(nbrecv,p)%x .gt. boundrightnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
        enddo
      endif
      if (anb(nbrecv,p)%y .gt. boundbacknb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
        enddo
      endif
    endif
    nbrecv=5
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
      if (anb(nbrecv,p)%x .gt. boundrightnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
        enddo
      endif
    endif
    nbrecv=6
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
      boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
      if (anb(nbrecv,p)%x .gt. boundrightnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
        enddo
      endif
      if (anb(nbrecv,p)%y .lt. boundfrontnb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
        enddo
      endif
    endif
    nbrecv=7
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
      if (anb(nbrecv,p)%y .lt. boundfrontnb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
        enddo
      endif
    endif
    nbrecv=8
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left  boundary of neighbor nb
      boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
      if (anb(nbrecv,p)%x .lt. boundleftnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
        enddo
      endif
      if (anb(nbrecv,p)%y .lt. boundfrontnb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
        enddo
      endif
    endif
  endif
enddo
!$omp end parallel
!
! third step: perform partial integration.
!
!$omp parallel default(shared) &
!$omp&private(p,l,nb)
!$omp do 
do p=1,pmax
  do l=1,nl
    do nb=0,8
      anb(nb,p)%ul(l) = 0.
      anb(nb,p)%vl(l) = 0.
      anb(nb,p)%wl(l) = 0.
    enddo
  enddo
enddo
!$omp end parallel
!
!rcoeff = 1.*ibmiter/max(1,ibmiter)
!$omp parallel default(shared) &
!$omp&private(p,l,nbrecv,coorxfp,cooryfp,coorzfp,ilow,ihigh,jlow,jhigh,klow,khigh) &
!$omp&private(coorzs,coorz,kernelzs,kernelz,coorys,coory,kernelys,kernely,coorxs,coorx) &
!$omp&private(kernelxs,kernelx,nbsend,i,j,k,nb)
!!$omp do schedule(dynamic)
!$omp do 
do p=1,pmax
  if (ap(p)%mslv .ne. 0) then
    ! myid is master or slave of particle abs(ap(p)%mslv)
    do nb=0,8
      if((nb.gt.0.and.ap(p)%mslv.lt.0.and.ap(p)%nb(nb).eq.1) .or. & !slave
         (nb.eq.0.and.ap(p)%mslv.gt.0) .or. & !pure master
         (nb.gt.0.and.ap(p)%mslv.gt.0.and.ap(p)%nb(nb).eq.1.and.neighbor(nb) .eq. myid)) then
         ! master that looks like a slave due to periodic bcs
        nbrecv = nb
        ! ap(p)%nb(nbrecv) is rank of master of particle -ap(p)%mslv
        if ( neighbor(nb) .eq. myid .and. nb.gt.0.and.ap(p)%mslv.gt.0) then
          nbrecv = nb + 4
          if (nbrecv .gt. 8) nbrecv = nbrecv-8
        endif
        do l=1,nl
          coorzfp =  anb(nbrecv,p)%zfp(l)*dzi
          if ((coorzfp>=nf.and.ifoverlap==1).or.ifoverlap==0) then
          coorxfp = (anb(nbrecv,p)%xfp(l)-boundleftmyid)*dxi
          cooryfp = (anb(nbrecv,p)%yfp(l)-boundfrontmyid)*dyi
          ilow  = nint( coorxfp - 2. )
          ihigh = nint( coorxfp + 2. )
          jlow  = nint( cooryfp - 2. )
          jhigh = nint( cooryfp + 2. )
          klow  = nint( coorzfp - 2. )
          khigh = nint( coorzfp + 2. )
          if (ilow .lt. 1) ilow = 1
          if (jlow .lt. 1) jlow = 1
          if (klow .lt. 1) klow = 1
          if (ihigh .gt. imax) ihigh = imax
          if (jhigh .gt. jmax) jhigh = jmax
          if (khigh .gt. kmax) khigh = kmax
          do k=klow,khigh
            coorzs   = (1.*k)-coorzfp
            coorz    = coorzs-0.5           !vert.    distance in grid points
            kernelzs = kernel(coorzs)
            kernelz  = kernel(coorz)
            do j=jlow,jhigh
              coorys   = (1.*j)-cooryfp
              coory    = coorys-0.5         !spanw.   distance in grid points
              kernelys = kernel(coorys)
              kernely  = kernel(coory)
              do i=ilow,ihigh
                coorxs   = (1.*i)-coorxfp
                coorx    = coorxs-0.5       !streamw. distance in grid points
                kernelxs = kernel(coorxs)
                kernelx  = kernel(coorx)



                anb(nbrecv,p)%ul(l) = anb(nbrecv,p)%ul(l) + unew(i,j,k)*kernelxs*kernely*kernelz
                anb(nbrecv,p)%vl(l) = anb(nbrecv,p)%vl(l) + vnew(i,j,k)*kernelx*kernelys*kernelz
                anb(nbrecv,p)%wl(l) = anb(nbrecv,p)%wl(l) + wnew(i,j,k)*kernelx*kernely*kernelzs

              enddo
            enddo
          enddo
         endif
        enddo ! do l=
      endif
    enddo ! do nbrecv=
  endif
enddo
!$omp end parallel
!
! fourth step: communicate data of slaves to their masters
! note : after debugging take sumu out and change 4 to 3
!
do p=1,pmax
  nrrequests = 0
  do nb=1,8
    nbsend = nb ! rank of process which sends data ('data is received from neighbor nbsend')
    idp = abs(ap(p)%mslv)
    tag = idp*10+nbsend
    nbrecv  = nb+4 ! rank of process which receives data ('data is send to neighbor nbrecv')
    if (nbrecv .gt. 8) nbrecv = nbrecv - 8
    if (ap(p)%mslv .gt. 0) then
      ! myid is master of particle ap(p)%mslv
      if (ap(p)%nb(nbsend) .eq. 1) then
        ! neighbor(nbsend) is rank of slave for particle ap(p)%mslv
        if ( neighbor(nbsend) .ne. myid ) then
          nrrequests = nrrequests + 1
          call MPI_IRECV(anb(nbsend,p)%ul(1),3*nl,MPI_REAL8,neighbor(nbsend), &
                         tag,comm_cart,arrayrequests((nrrequests-1)+1),error)
          ! recv dudt(:),dvdt(:),dwdt(:),sumu(:) -> 4*(nl) contiguous info
          ! (see definition of type pneighbor in the begining of the subroutine)
        endif
      endif
    endif
    if (ap(p)%mslv .lt. 0) then
      ! myid is slave of particle -ap(p)%mslv
      if (ap(p)%nb(nbrecv) .eq. 1) then
        ! ap(p)%nb(nbrecv) is rank of master of particle -ap(p)%mslv
        nrrequests = nrrequests + 1
          call MPI_ISEND(anb(nbrecv,p)%ul(1),3*nl,MPI_REAL8,neighbor(nbrecv), &
                         tag,comm_cart,arrayrequests((nrrequests-1)+1),error)
          ! send dudt(:),dvdt(:),dwdt(:),sumu(:) -> 4*(nl) contiguous info
          ! (see definition of type pneighbor in the begining of the subroutine)
      endif
    endif
  enddo ! do nb=
  call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,error)
enddo
!
! Sum all contributions together.
!
!$omp parallel default(shared) &
!$omp&private(p,l,nb) 
!!$omp do schedule(dynamic)
!$omp do 
do p=1,pmax
  if (ap(p)%mslv .gt. 0) then
    do l=1,nl
      ap(p)%ul(l) = 0.
      ap(p)%vl(l) = 0.
      ap(p)%wl(l) = 0.
      do nb=0,8
        ap(p)%ul(l) = ap(p)%ul(l) + anb(nb,p)%ul(l)
        ap(p)%vl(l) = ap(p)%vl(l) + anb(nb,p)%vl(l)
        ap(p)%wl(l) = ap(p)%wl(l) + anb(nb,p)%wl(l)
      enddo
    enddo
  else
    do l=1,nl
      ap(p)%ul(l) = 0.
      ap(p)%vl(l) = 0.
      ap(p)%wl(l) = 0.
    enddo
  endif
enddo
!$omp end parallel
!
! data to file
!
!write(rankpr,'(i3.3)') myid
!do p=1,pmax
!  if (ap(p)%mslv .gt. 0) then
!    write(partnr,'(i3.3)') ap(p)%mslv
!    open(23,file=datadir//'partnr'//partnr//'_ranknr'//rankpr//'_sum_over_lfp.txt')
!    do l=1,nl
!      if ( abs(sumu(l,p)-1.) .gt. 1.e-12 ) then
!        write(6,*) 'myid,l,p,sumu = ',myid,l,p,sumu(l,p),ap(p)%xfp(l),ap(p)%yfp(l),ap(p)%zfp(l)
!        write(6,*) 'ERROR!'
!        call mpi_finalize(error)
!        stop
!      endif 
!      write(23,'(I5,I5,10E16.8)') l,p,sumu(l,p),anb(0,p)%sumu(l),anb(1,p)%sumu(l),anb(2,p)%sumu(l), &
!                          anb(3,p)%sumu(l),anb(4,p)%sumu(l),anb(5,p)%sumu(l),anb(6,p)%sumu(l), &
!                          anb(7,p)%sumu(l),anb(8,p)%sumu(l)
!    enddo
!    close(23)
!  endif
!enddo
!
return
end subroutine eulr2lagr


































subroutine lagr2eulr
implicit none
integer :: i,j,k,l,p
integer :: ilow,ihigh,jlow,jhigh,klow,khigh
real :: coorx,coory,coorz
real :: coorxs,coorys,coorzs
real :: coorxfp,cooryfp,coorzfp
real :: kernelx,kernely,kernelz
real :: kernelxs,kernelys,kernelzs
real :: forcex_sc,forcey_sc,forcez_sc
type pneighbor
  real :: x,y,z
  real, dimension(1:nl) :: xfp,yfp,zfp,fxl,fyl,fzl,area
!  real, dimension(1:nl) :: sumu
end type pneighbor
type(pneighbor), dimension(0:8,1:npmax) :: anb
integer :: nb,nbsend,nbrecv
integer :: nrrequests
integer :: arrayrequests(1:3) !3=3*1 (master might have 3 slaves)
integer :: arraystatuses(MPI_STATUS_SIZE,1:3)
real :: boundleftnb,boundrightnb,boundfrontnb,boundbacknb
real :: dVlagrdVeuli
!real, dimension(1:nl,1:npmax) :: sumu
integer :: idp,tag
dVeul=dx*dy*dz
!
dVlagrdVeuli = dVlagr/dVeul
!
!$omp workshare
anb(0,1:pmax)%x = ap(1:pmax)%x
anb(0,1:pmax)%y = ap(1:pmax)%y
anb(0,1:pmax)%z = ap(1:pmax)%z
anb(1:8,1:pmax)%x = 0.
anb(1:8,1:pmax)%y = 0.
anb(1:8,1:pmax)%z = 0.
forall(l=1:nl,p=1:pmax)
  anb(0,p)%xfp(l) = ap(p)%xfp(l)
  anb(0,p)%yfp(l) = ap(p)%yfp(l)
  anb(0,p)%zfp(l) = ap(p)%zfp(l)
  anb(0,p)%fxl(l) = ap(p)%fxl(l)
  anb(0,p)%fyl(l) = ap(p)%fyl(l)
  anb(0,p)%fzl(l) = ap(p)%fzl(l)
  anb(0,p)%area(l) = ap(p)%area(l)

end forall
forall(l=1:nl,i=1:8,p=1:pmax)
  anb(i,p)%xfp(l) = 0.
  anb(i,p)%yfp(l) = 0.
  anb(i,p)%zfp(l) = 0.
  anb(i,p)%fxl(l) = 0.
  anb(i,p)%fyl(l) = 0.
  anb(i,p)%fzl(l) = 0.
  anb(i,p)%area(l) = 0.

end forall
!$omp end workshare
!
do p=1,pmax
  nrrequests = 0
  do nb=1,8
    nbsend = nb    ! rank of process which receives data ('sends data to neighbor nbsend')
    idp = abs(ap(p)%mslv)
    tag = idp*10+nbsend
    nbrecv  = nb+4  ! rank of process which sends data ('receives data from neighbor nbrecv')
    if (nbrecv .gt. 8) nbrecv = nbrecv - 8
    if (ap(p)%mslv .gt. 0) then
      ! myid is master of particle ap(p)%mslv
      if (ap(p)%nb(nbsend) .eq. 1) then
        ! neighbor(nbsend) is rank of slave for particle ap(p)%mslv
        if ( neighbor(nbsend) .eq. myid ) then
          ! process might be both master and slave of same particle due to periodic b.c.'s
          anb(nbrecv,p)%x = anb(0,p)%x
          anb(nbrecv,p)%y = anb(0,p)%y
          anb(nbrecv,p)%z = anb(0,p)%z
!$omp parallel default(shared) &
!$omp&private(l)
!$omp do
          do l=1,nl
            anb(nbrecv,p)%xfp(l) = anb(0,p)%xfp(l)
            anb(nbrecv,p)%yfp(l) = anb(0,p)%yfp(l)
            anb(nbrecv,p)%zfp(l) = anb(0,p)%zfp(l)
            anb(nbrecv,p)%fxl(l) = anb(0,p)%fxl(l)
            anb(nbrecv,p)%fyl(l) = anb(0,p)%fyl(l)
            anb(nbrecv,p)%fzl(l) = anb(0,p)%fzl(l)
            anb(nbrecv,p)%area(l) = anb(0,p)%area(l)

          enddo
          if (nbrecv .eq. 1) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
!$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
            enddo
          endif
          if (nbrecv .eq. 2) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
            anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
!$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
            enddo
          endif
          if (nbrecv .eq. 3) then
            anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
!$omp do
            do l=1,nl
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
            enddo
          endif
          if (nbrecv .eq. 4) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
            anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
!$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
            enddo
          endif
          if (nbrecv .eq. 5) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
!$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
            enddo
          endif
          if (nbrecv .eq. 6) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
            anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
!$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
            enddo
          endif
          if (nbrecv .eq. 7) then
            anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
!$omp do
            do l=1,nl
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
            enddo
          endif
          if (nbrecv .eq. 8) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
            anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
!$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
            enddo
          endif
!$omp end parallel
        else
          nrrequests = nrrequests + 1
          call MPI_ISEND(anb(0,p)%x,(3+7*nl),MPI_REAL8,neighbor(nbsend), &
                         tag,comm_cart,arrayrequests((nrrequests-1)+1),error)
          ! send x,y,z,xfp,yfp,zfp,fxl,fyl,fzl -> 3*(1+2*nl) contiguous info
          ! (see definition of type pneighbor in the begining of the subroutine)
        endif
      endif
    endif
    if (ap(p)%mslv .lt. 0) then
      ! myid is slave of particle -ap(p)%mslv
      if (ap(p)%nb(nbrecv) .eq. 1) then
        ! ap(p)%nb(nbrecv) is rank of master of particle -ap(p)%mslv
        nrrequests = nrrequests + 1
        call MPI_IRECV(anb(nbrecv,p)%x,(3+7*nl),MPI_REAL8,neighbor(nbrecv), &
                       tag,comm_cart,arrayrequests((nrrequests-1)+1),error)
        ! recv x,y,z,xfp,yfp,zfp,fxl,fyl,fzl -> 3*(1+2*nl) contiguous info
        ! (see definition of type pneighbor in the begining of the subroutine)
      endif
    endif
  enddo ! do nb=
  call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,error)
enddo
!
! second step: recompute particle positions for slaves due to periodic b.c.'s.
! required: (part of) particle within domain bounds of slave process.
!
!$omp parallel default(shared) &
!$omp&private(p,l,nbrecv,boundleftnb,boundbacknb,boundrightnb,boundfrontnb)  
!!$omp do schedule(dynamic)
!$omp do
do p=1,pmax
  if (ap(p)%mslv .lt. 0) then
    ! myid is slave of particle -ap(p)%mslv
    nbrecv=1
    if (ap(p)%nb(nbrecv) .eq. 1) then
      ! ap(p)%nb(nbrecv) is rank of master of particle -ap(p)%mslv
      boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left boundary of neighbor nb
      if (anb(nbrecv,p)%x .lt. boundleftnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
        enddo
      endif
    endif
    nbrecv=2
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left boundary of neighbor nb
      boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back boundary of neighbor nb
      if (anb(nbrecv,p)%x .lt. boundleftnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
        enddo
      endif
      if (anb(nbrecv,p)%y .gt. boundbacknb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
        enddo
      endif
    endif
    nbrecv=3
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back boundary of neighbor nb
      if (anb(nbrecv,p)%y .gt. boundbacknb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
        enddo
      endif
    endif
    nbrecv=4
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
      boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back  boundary of neighbor nb
      if (anb(nbrecv,p)%x .gt. boundrightnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
        enddo
      endif
      if (anb(nbrecv,p)%y .gt. boundbacknb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
        enddo
      endif
    endif
    nbrecv=5
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
      if (anb(nbrecv,p)%x .gt. boundrightnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
        enddo
      endif
    endif
    nbrecv=6
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
      boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
      if (anb(nbrecv,p)%x .gt. boundrightnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
        enddo
      endif
      if (anb(nbrecv,p)%y .lt. boundfrontnb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
        enddo
      endif
    endif
    nbrecv=7
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
      if (anb(nbrecv,p)%y .lt. boundfrontnb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
        enddo
      endif
    endif
    nbrecv=8
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left  boundary of neighbor nb
      boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
      if (anb(nbrecv,p)%x .lt. boundleftnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
        enddo
      endif
      if (anb(nbrecv,p)%y .lt. boundfrontnb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
        enddo
      endif
    endif
  endif
enddo
!$omp end parallel
!
! third step: perform partial integration.
!
!$omp workshare
forcex(:,:,:) = 0.
forcey(:,:,:) = 0.
forcez(:,:,:) = 0.
!$omp end workshare
!
!$omp parallel default(shared) &
!$omp&private(p,nbrecv,l,coorxfp,cooryfp,coorzfp,ilow,ihigh,jlow,jhigh,klow,khigh) &
!$omp&private(i,j,k,nbsend,nb) &
!$omp&private(coorzs,coorz,kernelzs,kernelz,coorys,coory,kernelys,kernely,coorxs,coorx,kernelxs,kernelx,forcex_sc,forcey_sc,forcez_sc)
!!$omp do schedule(dynamic)
!$omp do 
do p=1,pmax
  if (ap(p)%mslv .ne. 0) then
    ! myid is master or slave of particle abs(ap(p)%mslv)
    do nb=0,8
      if((nb.gt.0.and.ap(p)%mslv.lt.0.and.ap(p)%nb(nb).eq.1) .or. & !slave
         (nb.eq.0.and.ap(p)%mslv.gt.0) .or. & !pure master
         (nb.gt.0.and.ap(p)%mslv.gt.0.and.ap(p)%nb(nb).eq.1.and.neighbor(nb) .eq. myid)) then
        nbrecv = nb
        ! ap(p)%nb(nbrecv) is rank of master of particle -ap(p)%mslv
        if ( neighbor(nb) .eq. myid .and. nb.gt.0.and.ap(p)%mslv.gt.0) then
          nbrecv = nb + 4
          if (nbrecv .gt. 8) nbrecv = nbrecv-8
        endif
        do l=1,nl
          coorxfp = (anb(nbrecv,p)%xfp(l)-boundleftmyid)*dxi
          cooryfp = (anb(nbrecv,p)%yfp(l)-boundfrontmyid)*dyi
          coorzfp =  anb(nbrecv,p)%zfp(l)*dzi
          ilow  = nint( coorxfp - 2. )
          ihigh = nint( coorxfp + 2. )
          jlow  = nint( cooryfp - 2. )
          jhigh = nint( cooryfp + 2. )
          klow  = nint( coorzfp - 2. )
          khigh = nint( coorzfp + 2. )
          if (ilow .lt. 1) ilow = 1
          if (jlow .lt. 1) jlow = 1
          if (klow .lt. 1) klow = 1
          if (ihigh .gt. imax) ihigh = imax
          if (jhigh .gt. jmax) jhigh = jmax
          if (khigh .gt. kmax) khigh = kmax
          do k=klow,khigh
            coorzs   = (1.*k)-coorzfp
            coorz    = coorzs-0.5           !vert.    distance in grid points
            kernelzs = kernel(coorzs)
            kernelz  = kernel(coorz)
            do j=jlow,jhigh
              coorys   = (1.*j)-cooryfp
              coory    = coorys-0.5         !spanw.   distance in grid points
              kernelys = kernel(coorys)
              kernely  = kernel(coory)
              do i=ilow,ihigh
                coorxs   = (1.*i)-coorxfp
                coorx    = coorxs-0.5       !streamw. distance in grid points
                kernelxs = kernel(coorxs)
                kernelx  = kernel(coorx)


!$omp atomic
                forcex(i,j,k) = forcex(i,j,k) + anb(nbrecv,p)%fxl(l)*kernelxs*kernely*kernelz*(anb(nbrecv,p)%area(l)/dVeul)
!$omp atomic
                forcey(i,j,k) = forcey(i,j,k) + anb(nbrecv,p)%fyl(l)*kernelx*kernelys*kernelz*(anb(nbrecv,p)%area(l)/dVeul)
!$omp atomic
                forcez(i,j,k) = forcez(i,j,k) + anb(nbrecv,p)%fzl(l)*kernelx*kernely*kernelzs*(anb(nbrecv,p)%area(l)/dVeul)

              enddo
            enddo
          enddo
        enddo !do l=
      endif
    enddo !do nbrecv=
  endif
enddo
!$omp end parallel





!do k=1,kmax
!  do j=1,jmax
!     do i=1,imax

!       write(6,*) 'forcey',forcey(i,j,k)

!              enddo
!            enddo
!          enddo


!
! communicate data in x direction (periodic b.c.'s incorporated)
!
call updthalos(forcex,1)
call updthalos(forcey,1)
call updthalos(forcez,1)
!
! communicate data in y direction (periodic b.c.'s incorporated)
!
call updthalos(forcex,2)
call updthalos(forcey,2)
call updthalos(forcez,2)

! call updthalos(dudt,1)
! call updthalos(dvdt,1)
! call updthalos(dwdt,1)
! call updthalos(dudt,2)
! call updthalos(dvdt,2)
! call updthalos(dwdt,2)
!
return
end subroutine lagr2eulr























































subroutine lagr2eulrind
implicit none
integer :: i,j,k,l,p
integer :: ilow,ihigh,jlow,jhigh,klow,khigh
real :: coorx,coory,coorz
real :: coorxs,coorys,coorzs
real :: coorxfp,cooryfp,coorzfp
real :: kernelx,kernely,kernelz
real :: kernelxs,kernelys,kernelzs
real :: forcex_sc,forcey_sc,forcez_sc
type pneighbor
  real :: x,y,z
  real, dimension(1:nl) :: xfp,yfp,zfp,xn,yn,zn,area
!  real, dimension(1:nl) :: sumu
end type pneighbor
type(pneighbor), dimension(0:8,1:npmax) :: anb
integer :: nb,nbsend,nbrecv
integer :: nrrequests
integer :: arrayrequests(1:3) !3=3*1 (master might have 3 slaves)
integer :: arraystatuses(MPI_STATUS_SIZE,1:3)
real :: boundleftnb,boundrightnb,boundfrontnb,boundbacknb
real :: dVlagrdVeuli
!real, dimension(1:nl,1:npmax) :: sumu
integer :: idp,tag
dVeul=dx*dy*dz
!
dVlagrdVeuli = dVlagr/dVeul
!
!$omp workshare
anb(0,1:pmax)%x = ap(1:pmax)%x
anb(0,1:pmax)%y = ap(1:pmax)%y
anb(0,1:pmax)%z = ap(1:pmax)%z
anb(1:8,1:pmax)%x = 0.
anb(1:8,1:pmax)%y = 0.
anb(1:8,1:pmax)%z = 0.
forall(l=1:nl,p=1:pmax)
  anb(0,p)%xfp(l) = ap(p)%xfp(l)
  anb(0,p)%yfp(l) = ap(p)%yfp(l)
  anb(0,p)%zfp(l) = ap(p)%zfp(l)
  anb(0,p)%xn(l) = ap(p)%xn(l)
  anb(0,p)%yn(l) = ap(p)%yn(l)
  anb(0,p)%zn(l) = ap(p)%zn(l)
  anb(0,p)%area(l) = ap(p)%area(l)
end forall
forall(l=1:nl,i=1:8,p=1:pmax)
  anb(i,p)%xfp(l) = 0.
  anb(i,p)%yfp(l) = 0.
  anb(i,p)%zfp(l) = 0.
  anb(i,p)%xn(l) = 0.
  anb(i,p)%yn(l) = 0.
  anb(i,p)%zn(l) = 0.
  anb(i,p)%area(l) = 0.
end forall
!$omp end workshare
!
do p=1,pmax
  nrrequests = 0
  do nb=1,8
    nbsend = nb    ! rank of process which receives data ('sends data to neighbor nbsend')
    idp = abs(ap(p)%mslv)
    tag = idp*10+nbsend
    nbrecv  = nb+4  ! rank of process which sends data ('receives data from neighbor nbrecv')
    if (nbrecv .gt. 8) nbrecv = nbrecv - 8
    if (ap(p)%mslv .gt. 0) then
      ! myid is master of particle ap(p)%mslv
      if (ap(p)%nb(nbsend) .eq. 1) then
        ! neighbor(nbsend) is rank of slave for particle ap(p)%mslv
        if ( neighbor(nbsend) .eq. myid ) then
          ! process might be both master and slave of same particle due to periodic b.c.'s
          anb(nbrecv,p)%x = anb(0,p)%x
          anb(nbrecv,p)%y = anb(0,p)%y
          anb(nbrecv,p)%z = anb(0,p)%z
!$omp parallel default(shared) &
!$omp&private(l)
!$omp do
          do l=1,nl
            anb(nbrecv,p)%xfp(l) = anb(0,p)%xfp(l)
            anb(nbrecv,p)%yfp(l) = anb(0,p)%yfp(l)
            anb(nbrecv,p)%zfp(l) = anb(0,p)%zfp(l)
            anb(nbrecv,p)%xn(l) = anb(0,p)%xn(l)
            anb(nbrecv,p)%yn(l) = anb(0,p)%yn(l)
            anb(nbrecv,p)%zn(l) = anb(0,p)%zn(l)
            anb(nbrecv,p)%area(l) = anb(0,p)%area(l)
          enddo
          if (nbrecv .eq. 1) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
!$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
            enddo
          endif
          if (nbrecv .eq. 2) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
            anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
!$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
            enddo
          endif
          if (nbrecv .eq. 3) then
            anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
!$omp do
            do l=1,nl
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
            enddo
          endif
          if (nbrecv .eq. 4) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
            anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
!$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
            enddo
          endif
          if (nbrecv .eq. 5) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
!$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
            enddo
          endif
          if (nbrecv .eq. 6) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
            anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
!$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
            enddo
          endif
          if (nbrecv .eq. 7) then
            anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
!$omp do
            do l=1,nl
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
            enddo
          endif
          if (nbrecv .eq. 8) then
            anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
            anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
!$omp do
            do l=1,nl
              anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
              anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
            enddo
          endif
!$omp end parallel
        else
          nrrequests = nrrequests + 1
          call MPI_ISEND(anb(0,p)%x,(3+7*nl),MPI_REAL8,neighbor(nbsend), &
                         tag,comm_cart,arrayrequests((nrrequests-1)+1),error)
          ! send x,y,z,xfp,yfp,zfp,fxl,fyl,fzl -> 3*(1+2*nl) contiguous info
          ! (see definition of type pneighbor in the begining of the subroutine)
        endif
      endif
    endif
    if (ap(p)%mslv .lt. 0) then
      ! myid is slave of particle -ap(p)%mslv
      if (ap(p)%nb(nbrecv) .eq. 1) then
        ! ap(p)%nb(nbrecv) is rank of master of particle -ap(p)%mslv
        nrrequests = nrrequests + 1
        call MPI_IRECV(anb(nbrecv,p)%x,(3+7*nl),MPI_REAL8,neighbor(nbrecv), &
                       tag,comm_cart,arrayrequests((nrrequests-1)+1),error)
        ! recv x,y,z,xfp,yfp,zfp,fxl,fyl,fzl -> 3*(1+2*nl) contiguous info
        ! (see definition of type pneighbor in the begining of the subroutine)
      endif
    endif
  enddo ! do nb=
  call MPI_WAITALL(nrrequests,arrayrequests,arraystatuses,error)
enddo
!
! second step: recompute particle positions for slaves due to periodic b.c.'s.
! required: (part of) particle within domain bounds of slave process.
!
!$omp parallel default(shared) &
!$omp&private(p,l,nbrecv,boundleftnb,boundbacknb,boundrightnb,boundfrontnb)  
!!$omp do schedule(dynamic)
!$omp do
do p=1,pmax
  if (ap(p)%mslv .lt. 0) then
    ! myid is slave of particle -ap(p)%mslv
    nbrecv=1
    if (ap(p)%nb(nbrecv) .eq. 1) then
      ! ap(p)%nb(nbrecv) is rank of master of particle -ap(p)%mslv
      boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left boundary of neighbor nb
      if (anb(nbrecv,p)%x .lt. boundleftnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
        enddo
      endif
    endif
    nbrecv=2
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left boundary of neighbor nb
      boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back boundary of neighbor nb
      if (anb(nbrecv,p)%x .lt. boundleftnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
        enddo
      endif
      if (anb(nbrecv,p)%y .gt. boundbacknb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
        enddo
      endif
    endif
    nbrecv=3
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back boundary of neighbor nb
      if (anb(nbrecv,p)%y .gt. boundbacknb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
        enddo
      endif
    endif
    nbrecv=4
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
      boundbacknb  = (coords(2))*ly/(1.*dims(2)) ! back  boundary of neighbor nb
      if (anb(nbrecv,p)%x .gt. boundrightnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
        enddo
      endif
      if (anb(nbrecv,p)%y .gt. boundbacknb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y - ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) - ly
        enddo
      endif
    endif
    nbrecv=5
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
      if (anb(nbrecv,p)%x .gt. boundrightnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
        enddo
      endif
    endif
    nbrecv=6
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundrightnb = (coords(1))*lx/(1.*dims(1)) ! right boundary of neighbor nb
      boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
      if (anb(nbrecv,p)%x .gt. boundrightnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x - lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) - lx
        enddo
      endif
      if (anb(nbrecv,p)%y .lt. boundfrontnb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
        enddo
      endif
    endif
    nbrecv=7
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
      if (anb(nbrecv,p)%y .lt. boundfrontnb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
        enddo
      endif
    endif
    nbrecv=8
    if (ap(p)%nb(nbrecv) .gt. 0) then
      boundleftnb  = (coords(1)+1)*lx/(1.*dims(1)) ! left  boundary of neighbor nb
      boundfrontnb = (coords(2)+1)*ly/(1.*dims(2)) ! front boundary of neighbor nb
      if (anb(nbrecv,p)%x .lt. boundleftnb) then
        anb(nbrecv,p)%x = anb(nbrecv,p)%x + lx
        do l=1,nl
          anb(nbrecv,p)%xfp(l) = anb(nbrecv,p)%xfp(l) + lx
        enddo
      endif
      if (anb(nbrecv,p)%y .lt. boundfrontnb) then
        anb(nbrecv,p)%y = anb(nbrecv,p)%y + ly
        do l=1,nl
          anb(nbrecv,p)%yfp(l) = anb(nbrecv,p)%yfp(l) + ly
        enddo
      endif
    endif
  endif
enddo
!$omp end parallel
!
! third step: perform partial integration.
!
!$omp workshare
gx(:,:,:) = 0.
gy(:,:,:) = 0.
gz(:,:,:) = 0.
!$omp end workshare
!
!$omp parallel default(shared) &
!$omp&private(p,nbrecv,l,coorxfp,cooryfp,coorzfp,ilow,ihigh,jlow,jhigh,klow,khigh) &
!$omp&private(i,j,k,nbsend,nb) &
!$omp&private(coorzs,coorz,kernelzs,kernelz,coorys,coory,kernelys,kernely,coorxs,coorx,kernelxs,kernelx,forcex_sc,forcey_sc,forcez_sc)
!!$omp do schedule(dynamic)
!$omp do 
do p=1,pmax
  if (ap(p)%mslv .ne. 0) then
    ! myid is master or slave of particle abs(ap(p)%mslv)
    do nb=0,8
      if((nb.gt.0.and.ap(p)%mslv.lt.0.and.ap(p)%nb(nb).eq.1) .or. & !slave
         (nb.eq.0.and.ap(p)%mslv.gt.0) .or. & !pure master
         (nb.gt.0.and.ap(p)%mslv.gt.0.and.ap(p)%nb(nb).eq.1.and.neighbor(nb) .eq. myid)) then
        nbrecv = nb
        ! ap(p)%nb(nbrecv) is rank of master of particle -ap(p)%mslv
        if ( neighbor(nb) .eq. myid .and. nb.gt.0.and.ap(p)%mslv.gt.0) then
          nbrecv = nb + 4
          if (nbrecv .gt. 8) nbrecv = nbrecv-8
        endif
        do l=1,nl
          coorxfp = (anb(nbrecv,p)%xfp(l)-boundleftmyid)*dxi
          cooryfp = (anb(nbrecv,p)%yfp(l)-boundfrontmyid)*dyi
          coorzfp =  anb(nbrecv,p)%zfp(l)*dzi
          ilow  = nint( coorxfp - 2. )
          ihigh = nint( coorxfp + 2. )
          jlow  = nint( cooryfp - 2. )
          jhigh = nint( cooryfp + 2. )
          klow  = nint( coorzfp - 2. )
          khigh = nint( coorzfp + 2. )
          if (ilow .lt. 1) ilow = 0
          if (jlow .lt. 1) jlow = 0
          if (klow .lt. 1) klow = 0
          if (ihigh .gt. imax) ihigh = imax+1
          if (jhigh .gt. jmax) jhigh = jmax+1
          if (khigh .gt. kmax) khigh = kmax+1
          do k=klow,khigh
            coorzs   = (1.*k)-coorzfp
            coorz    = coorzs-0.5           !vert.    distance in grid points
            kernelzs = kernel(coorzs)
            kernelz  = kernel(coorz)
            do j=jlow,jhigh
              coorys   = (1.*j)-cooryfp
              coory    = coorys-0.5         !spanw.   distance in grid points
              kernelys = kernel(coorys)
              kernely  = kernel(coory)
              do i=ilow,ihigh
                coorxs   = (1.*i)-coorxfp
                coorx    = coorxs-0.5       !streamw. distance in grid points
                kernelxs = kernel(coorxs)
                kernelx  = kernel(coorx)

!$omp atomic
                gx(i,j,k) = gx(i,j,k) + anb(nbrecv,p)%xn(l)*kernelx*kernely*kernelz*(anb(nbrecv,p)%area(l)/dVeul)
!$omp atomic
                gy(i,j,k) = gy(i,j,k) + anb(nbrecv,p)%yn(l)*kernelx*kernely*kernelz*(anb(nbrecv,p)%area(l)/dVeul)
!$omp atomic
                gz(i,j,k) = gz(i,j,k) + anb(nbrecv,p)%zn(l)*kernelx*kernely*kernelz*(anb(nbrecv,p)%area(l)/dVeul)
!!$omp atomic
              enddo
            enddo
          enddo
        enddo !do l=
      endif
    enddo !do nbrecv=
  endif
enddo
!$omp end parallel





!do k=1,kmax
!  do j=1,jmax
!     do i=1,imax

!       write(6,*) 'forcey',forcey(i,j,k)

!              enddo
!            enddo
!          enddo


!
! communicate data in x direction (periodic b.c.'s incorporated)
!

call updthalos(gx,1)
call updthalos(gy,1)
call updthalos(gz,1)
!
! communicate data in y direction (periodic b.c.'s incorporated)
!

call updthalos(gx,2)
call updthalos(gy,2)
call updthalos(gz,2)
! call updthalos(dudt,1)
! call updthalos(dvdt,1)
! call updthalos(dwdt,1)
! call updthalos(dudt,2)
! call updthalos(dvdt,2)
! call updthalos(dwdt,2)
!
return
end subroutine lagr2eulrind
!
end module mod_interp_spread
