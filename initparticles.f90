module mod_initparticles
use decomp_2d
use mod_param
use mod_common
use mod_common_mpi
use mod_sph
implicit none
private
public initparticles
contains
subroutine initparticles
implicit none
!integer,parameter extrapoints = nint(2*(radius+0.1)*dxi)! for visualisation purposes
real, dimension(np) :: xcglob,ycglob,zcglob,thetacglob,phicglob
integer :: i,j,k,p,pp,rk,icap,jcap
integer :: proccoords(1:ndims),procrank,cap
integer, dimension(2) :: sbuf,rbuf
real :: leftbound,rightbound,frontbound,backbound
real :: dist,distx,disty,distz,distzw,angle
real :: xp,yp
real :: ax
real :: ay
real :: rn
integer :: counter,crys
integer :: idp
character(len=5) rankpr 
character(len=7) :: tempstr
character(len=11) :: tempstr2
integer :: count_mstr,count_slve,count_mstr_all,count_slve_loc
real, dimension(NL) :: xfp_temp,yfp_temp,zfp_temp 
real :: d_arc,theta,phy
d_arc=(2.*pi)/m_capf


!**************************************************************************initialization for spherical harmonics*****************************************************
!    do l=1,NL      
!    icap=int(floor(real(((l-1)/m_cap)))+1)
!    jcap=l-((icap-1)*m_cap)
!    theta=(jcap-.5)*d_arc
!    phy=(icap-.5)*d_arc
!    crdyR(icap,jcap) = diam*.5*cos(phy)
!    crdxR(icap,jcap) = diam*.5*sin(phy)*sin(theta) 
!    crdzR(icap,jcap) = diam*.5*sin(phy)*cos(theta)
!    enddo
    
    if (begin==0) then

    do icap=0,(dims(1)*dims(2))-1
    if (myid==icap) then
    call capsule_init_cal_coeff   !(crdxR,crdyR,crdzR)
    intl=0
    igrida(1)=-2
    igrida(2)=1
    igridb(1)=-2
    igridb(2)=1
    call TRSSPH (intl,igrida,m_cap,n_cap,surfmetD,igridb,m_capf,n_capf, &
    surfmetDf,wsave,lssave,lsvmin,work,lswork,lwkmin,dwork,lsdwork,ier)
    do i=1,n_capf
    do j=1,m_capf 
    area_cap(i,j)=surfmetDf(i,j)*wgtdsinf(i,j)
    enddo
    enddo
    do p=1,np 
    do i=1,n_capf
    do j=1,m_capf
    cap=((i-1)*m_capf)+j
    ap(p)%area(cap)=area_cap(i,j)
    enddo
    enddo
    enddo
    endif
    enddo

!*********************************************************************************************************************************************************************
!
! position of spheres: global initialisation by root (myid = 0)
!
crys=0
if (myid .eq. 0) then
  !!!call random_seed( put = (/16112006/))
  write(6,*) 'part as crystal 1, random 0: ',crys
  open(23,file=datadir//"position_spheres.txt")
  if(crys.eq.1) then
!    p=0
!    ! pseudo-crystal
!    do k=1,4
!      do j=1,16
!        do i=1,8
!          p=p+1
!          thetacglob(p) = 0.
!          phicglob(p)   = 0.
!          call random_number(rn)
!          xcglob(p)     = (lx/4.)*((1.*i)-0.5) + 0.25*(rn-0.5)*radius
!          call random_number(rn)
!          ycglob(p)     = (ly/8.)*((1.*j)-0.5) + 0.25*(rn-0.5)*radius
!          call random_number(rn)
!          zcglob(p)     = (lz/2. )*((1.*k)-0.5) + (rn-0.5)*radius
!          write(6,*) 'Location of sphere #',p
!          write(6,*) 'x,y,z = ',xcglob(p),ycglob(p),zcglob(p)
!          write(23,'(I3,5E32.16)') p,thetacglob(p),phicglob(p), &
!                                  xcglob(p),ycglob(p),zcglob(p)
!        enddo
!      enddo
!    enddo
  else
    ! pseudo-random
    counter = 0
    do p=1,NP
      if(NP.eq.1) then
            xcglob(p)     = lx*.5
            ycglob(p)     = ly*.25
            zcglob(p)     = lz*.2
      else      
        111  continue
        call random_number(rn)
        xcglob(p)     = lx*rn
        call random_number(rn)
        ycglob(p)     = ly*rn
        call random_number(rn)
        zcglob(p)     = lz*rn
        do pp=1,p
          distzw = min(abs(lz-zcglob(p)),abs(zcglob(p)))
          distz = abs(zcglob(p)- zcglob(pp))
          if(distzw.lt.1.05*radius) goto 111
          do j=-1,1
            disty = abs(ycglob(p)- (ycglob(pp)+j*ly))
            do i=-1,1
              distx = abs(xcglob(p)- (xcglob(pp)+i*lx))
              if(distx.gt.2.05*radius.or. &
                 disty.gt.2.05*radius.or. &
                 distz.gt.2.05*radius.or.p.eq.pp) then
                ! good particle
              else
                dist = distx**2+disty**2.+distz**2.
                if((dist.lt.(4.2*radius**2.))) then 
                  !write(*,*)'RANDOM DEVIATION'
                  !write(*,*)dist,distw
                  counter=counter+1
                  goto 111
                endif
              endif
            enddo
          enddo
        222 continue
        enddo
      endif
      write(6,*) 'Location of sphere #',p
      write(6,*) 'x,y,z = ',xcglob(p),ycglob(p),zcglob(p)
      write(23,'(I3,3E32.16)') p,xcglob(p),ycglob(p),zcglob(p) 
                               
    enddo
    write(*,*)'RANDOM DEVIATIONS: ',counter
    p=p-1
!    p = np
  endif
  close(23)
  if ( (p.ne.np) ) then
    print*,counter,np
    write(6,*) 'Fatal error in initialisation of particle positions!'
    write(6,*) 'Program aborted...'
    call mpi_finalize(error)
    stop
  endif
  do rk=1,Nproc-1
    call MPI_SSEND(xcglob    ,np,MPI_REAL8,rk,rk+0*(Nproc-1),comm_cart,error)
    call MPI_SSEND(ycglob    ,np,MPI_REAL8,rk,rk+1*(Nproc-1),comm_cart,error)
    call MPI_SSEND(zcglob    ,np,MPI_REAL8,rk,rk+2*(Nproc-1),comm_cart,error)
  enddo
else ! if myid is not 0:
  call MPI_RECV(xcglob    ,np,MPI_REAL8,0,myid+0*(Nproc-1),comm_cart,status,error)
  call MPI_RECV(ycglob    ,np,MPI_REAL8,0,myid+1*(Nproc-1),comm_cart,status,error)
  call MPI_RECV(zcglob    ,np,MPI_REAL8,0,myid+2*(Nproc-1),comm_cart,status,error)
endif
!
! Determine master and slave processes for each particle.
!
! initialisation
!
ap(1:npmax)%x = 0.
ap(1:npmax)%y = 0.
ap(1:npmax)%z = 0.
ap(1:npmax)%mslv = 0.
forall(i=1:npmax) ap(i)%nb(1:8) = 0
!
count_mstr = 0
count_slve = 0
i = 0
ax = 0.5
ay = 0.5
!
pmax = 0
do p=1,np
  if (xcglob(p).lt.0..or.xcglob(p).gt.lx .or. &
      ycglob(p).lt.0..or.ycglob(p).gt.ly .or. &
      zcglob(p).lt.radius.or.zcglob(p).gt.lz-radius) then
    if (myid.eq.0) then
      write(6,*) 'Fatal error in initialisation of particle positions - '
      write(6,*) 'particle outside the domain!'
      write(6,*) 'Program aborted...'
    endif
    !call mpi_finalize(error)
    !stop
  endif
  if (xcglob(p).eq.lx) ax = 0.51
  if (xcglob(p).eq.0) ax = 0.49
  if (ycglob(p).eq.ly) ay = 0.51
  if (ycglob(p).eq.0) ay = 0.49
!  
  proccoords(1) = nint(xcglob(p)*dims(1)/lx - ax)
  proccoords(2) = nint(ycglob(p)*dims(2)/ly - ay)
  leftbound     = (proccoords(1)  )*lx/(1.*dims(1)) ! left  boundary of particle's master
  rightbound    = (proccoords(1)+1)*lx/(1.*dims(1)) ! right boundary of particle's master
  frontbound    = (proccoords(2)  )*ly/(1.*dims(2)) ! front boundary of particle's master
  backbound     = (proccoords(2)+1)*ly/(1.*dims(2)) ! back  boundary of particle's master
  call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
  if (myid.eq.procrank) then
    i = i + 1
    count_mstr = count_mstr + 1
    ap(i)%x = xcglob(p)
    ap(i)%y = ycglob(p)
    ap(i)%z = zcglob(p)
    ap(i)%mslv = p

      
    do l=1,NL
      
      icap=int(floor(real(((l-1)/m_capf)))+1)
      jcap=l-((icap-1)*m_capf)
      theta=(jcap-.5)*d_arc
      phy=(icap-.5)*d_arc
      ap(i)%xfp(l) = ap(i)%x + crdxRf(icap,jcap)
      ap(i)%yfp(l) = ap(i)%y + crdyRf(icap,jcap) 
      ap(i)%zfp(l) = ap(i)%z + crdzRf(icap,jcap) 
        

      ap(i)%xfpo(l) = ap(i)%xfp(l)
      ap(i)%yfpo(l) = ap(i)%yfp(l)
      ap(i)%zfpo(l) = ap(i)%zfp(l)    
      

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
    proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax ) + 1
    proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) 
    call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
    if (myid .eq. procrank) then

      do l=1,NL

      icap=int(floor(real(((l-1)/m_capf)))+1)
      jcap=l-((icap-1)*m_capf)
      theta=(jcap-.5)*d_arc
      phy=(icap-.5)*d_arc
      xfp_temp(l) = xcglob(p) + crdxRf(icap,jcap)
      yfp_temp(l) = ycglob(p) + crdyRf(icap,jcap) 
      zfp_temp(l) = zcglob(p) + crdzRf(icap,jcap)   

        if ( (xfp_temp(l)+offset) .gt. rightbound ) then
          if ( ((yfp_temp(l)+offset) .ge. frontbound) .and. ((yfp_temp(l)-offset) .le. backbound) ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            count_slve = count_slve + 1
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
    proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax ) + 1
    proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) - 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
    if (myid .eq. procrank) then

      do l=1,NL

      icap=int(floor(real(((l-1)/m_capf)))+1)
      jcap=l-((icap-1)*m_capf)
      theta=(jcap-.5)*d_arc
      phy=(icap-.5)*d_arc
      xfp_temp(l) = xcglob(p) + crdxRf(icap,jcap)
      yfp_temp(l) = ycglob(p) + crdyRf(icap,jcap) 
      zfp_temp(l) = zcglob(p) + crdzRf(icap,jcap)  

        if ( (xfp_temp(l)+offset) .gt. rightbound ) then
          if  ( (yfp_temp(l)-offset) .le. frontbound ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            count_slve = count_slve + 1
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
    proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax )
    proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) - 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
    if (myid .eq. procrank) then


      do l=1,NL

      icap=int(floor(real(((l-1)/m_capf)))+1)
      jcap=l-((icap-1)*m_capf)
      theta=(jcap-.5)*d_arc
      phy=(icap-.5)*d_arc
      xfp_temp(l) = xcglob(p) + crdxRf(icap,jcap)
      yfp_temp(l) = ycglob(p) + crdyRf(icap,jcap) 
      zfp_temp(l) = zcglob(p) + crdzRf(icap,jcap)  

        if ( (yfp_temp(l)-offset) .lt. frontbound ) then
          if ( ((xfp_temp(l)+offset) .ge. leftbound) .and. ((xfp_temp(l)-offset) .le. rightbound) ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            count_slve = count_slve + 1
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
    proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax ) - 1
    proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) - 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
    if (myid .eq. procrank) then

      do l=1,NL

      icap=int(floor(real(((l-1)/m_capf)))+1)
      jcap=l-((icap-1)*m_capf)
      theta=(jcap-.5)*d_arc
      phy=(icap-.5)*d_arc
      xfp_temp(l) = xcglob(p) + crdxRf(icap,jcap)
      yfp_temp(l) = ycglob(p) + crdyRf(icap,jcap) 
      zfp_temp(l) = zcglob(p) + crdzRf(icap,jcap)  

        if ( (yfp_temp(l)-offset) .lt. frontbound ) then
          if  ( (xfp_temp(l)-offset) .le. leftbound ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            count_slve = count_slve + 1
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
    proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax ) - 1
    proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay )
    call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
    if (myid .eq. procrank) then

      do l=1,NL

      icap=int(floor(real(((l-1)/m_capf)))+1)
      jcap=l-((icap-1)*m_capf)
      theta=(jcap-.5)*d_arc
      phy=(icap-.5)*d_arc
      xfp_temp(l) = xcglob(p) + crdxRf(icap,jcap)
      yfp_temp(l) = ycglob(p) + crdyRf(icap,jcap) 
      zfp_temp(l) = zcglob(p) + crdzRf(icap,jcap)  

      if ( (xfp_temp(l)-offset) .lt. leftbound ) then
        if ( ((yfp_temp(l)+offset) .ge. frontbound) .and. ((yfp_temp(l)-offset) .le. backbound ) ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            count_slve = count_slve + 1
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
    proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax ) - 1
    proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) + 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
    if (myid .eq. procrank) then

      do l=1,NL

      icap=int(floor(real(((l-1)/m_capf)))+1)
      jcap=l-((icap-1)*m_capf)
      theta=(jcap-.5)*d_arc
      phy=(icap-.5)*d_arc
      xfp_temp(l) = xcglob(p) + crdxRf(icap,jcap)
      yfp_temp(l) = ycglob(p) + crdyRf(icap,jcap) 
      zfp_temp(l) = zcglob(p) + crdzRf(icap,jcap)  

      if ( (yfp_temp(l)+offset) .gt. backbound ) then
        if  ( (xfp_temp(l)-offset) .le. leftbound ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            count_slve = count_slve + 1
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
    proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax )
    proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) + 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
    if (myid .eq. procrank) then

      do l=1,NL

      icap=int(floor(real(((l-1)/m_capf)))+1)
      jcap=l-((icap-1)*m_capf)
      theta=(jcap-.5)*d_arc
      phy=(icap-.5)*d_arc
      xfp_temp(l) = xcglob(p) + crdxRf(icap,jcap)
      yfp_temp(l) = ycglob(p) + crdyRf(icap,jcap) 
      zfp_temp(l) = zcglob(p) + crdzRf(icap,jcap)  

      if ( (yfp_temp(l)+offset) .gt. backbound ) then
        if ( ((xfp_temp(l)+offset) .ge. leftbound) .and. ((xfp_temp(l)-offset) .le. rightbound) ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            count_slve = count_slve + 1
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
    proccoords(1) = nint( dims(1)*xcglob(p)/lx - ax ) + 1
    proccoords(2) = nint( dims(2)*ycglob(p)/ly - ay ) + 1
    call MPI_CART_RANK(comm_cart,proccoords,procrank,error)
    if (myid .eq. procrank) then

      do l=1,NL

      icap=int(floor(real(((l-1)/m_capf)))+1)
      jcap=l-((icap-1)*m_capf)
      theta=(jcap-.5)*d_arc
      phy=(icap-.5)*d_arc
      xfp_temp(l) = xcglob(p) + crdxRf(icap,jcap)
      yfp_temp(l) = ycglob(p) + crdyRf(icap,jcap) 
      zfp_temp(l) = zcglob(p) + crdzRf(icap,jcap)  

      if ( (yfp_temp(l)+offset) .gt. backbound ) then
        if  ( (xfp_temp(l)+offset) .ge. rightbound ) then 

            if(count_slve_loc.eq. 0 ) i = i+1
            count_slve = count_slve + 1
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


! maximum number of particles in a thread is equal to the number of particles
! 'mastered' and 'enslaved' by it
!
pmax = count_mstr+count_slve
npmstr = count_mstr

do p=1,pmax
  if (ap(p)%mslv .le. 0) then
    ! myid is either slave of particle ap(p)%mslv or does not contain this particle
    do l=1,NL
      ap(p)%xfp(l)  = 0.
      ap(p)%yfp(l)  = 0.
      ap(p)%zfp(l)  = 0.
      ap(p)%xfpo(l) = 0.
      ap(p)%yfpo(l) = 0.
      ap(p)%zfpo(l) = 0. 
      ap(p)%xfpold(l) = 0.
      ap(p)%yfpold(l) = 0.
      ap(p)%zfpold(l) = 0.
      ap(p)%area(l)   = 0.
      ap(p)%xfp0(l) = 0.
      ap(p)%yfp0(l) = 0.
      ap(p)%zfp0(l) = 0.
      ap(p)%xn(l) = 0.
      ap(p)%yn(l) = 0.
      ap(p)%zn(l) = 0.   
    enddo

  endif
enddo

ap(1:npmax)%integralv = 0.
forall (p=1:npmax)
  ap(p)%dxdt(:) = 0.
  ap(p)%dydt(:) = 0.
  ap(p)%dzdt(:) = 0.
  ap(p)%fxl(:)  = 0.
  ap(p)%fyl(:)  = 0.
  ap(p)%fzl(:)  = 0.
  ap(p)%ul(:)   = 0.
  ap(p)%vl(:)   = 0.
  ap(p)%wl(:)   = 0.
end forall

forall (p=pmax+1:npmax)
  ap(p)%xfp(:)  = 0.
  ap(p)%yfp(:)  = 0.
  ap(p)%zfp(:)  = 0.
  ap(p)%xfpo(:) = 0.
  ap(p)%yfpo(:) = 0.
  ap(p)%zfpo(:) = 0. 
  ap(p)%xfpold(:) = 0.
  ap(p)%yfpold(:) = 0.
  ap(p)%zfpold(:) = 0. 
  ap(p)%area(:)   = 0.
  ap(p)%xfp0(:) = 0.
  ap(p)%yfp0(:) = 0.
  ap(p)%zfp0(:) = 0.
  ap(p)%xn(:) = 0.
  ap(p)%yn(:) = 0.
  ap(p)%zn(:) = 0.
end forall






write(6,'(A7,I5,A8,I5,A18,I5,A11,A8,I5)') 'Thread ', myid, ' masters ', count_mstr, &
' and is slave for ', count_slve, ' particles. ', ' pmax = ', pmax

!
call MPI_ALLREDUCE(count_mstr,count_mstr_all,1,MPI_INTEGER,MPI_SUM,comm_cart,error)
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





   else


    do icap=0,(dims(1)*dims(2))-1
    if (myid==icap) then
    call capsule_init_cal_coeff   !(crdxR,crdyR,crdzR)
    intl=0
    igrida(1)=-2
    igrida(2)=1
    igridb(1)=-2
    igridb(2)=1
    call TRSSPH (intl,igrida,m_cap,n_cap,surfmetD,igridb,m_capf,n_capf, &
    surfmetDf,wsave,lssave,lsvmin,work,lswork,lwkmin,dwork,lsdwork,ier)
    do i=1,n_capf
    do j=1,m_capf 
    area_cap(i,j)=surfmetDf(i,j)*wgtdsinf(i,j)
    enddo
    enddo
    do p=1,np 
    do i=1,n_capf
    do j=1,m_capf
    cap=((i-1)*m_capf)+j
    ap(p)%area(cap)=area_cap(i,j)
    enddo
    enddo
    enddo
    endif
    enddo



    endif
!
return
end subroutine initparticles
!
end module mod_initparticles
