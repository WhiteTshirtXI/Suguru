module mod_boundc
use mod_common
implicit none
private
public bounduvwc, boundpc, boundpcorrc, updthalosc, boundindch,boundviscc
contains
!
subroutine bounduvwc(uc,vc,wc)
use mod_param
use mod_common
implicit none
integer :: i,j
real, dimension(0:,0:,0:) :: uc,vc,wc
!
select case(iniu)
case('poi','zer','log')
if(isfreeslip) then
!$omp parallel default(none) shared(uc,vc,wc) private(i,j)
!$omp do
  do j=0,j1c
    do i=0,i1c
      



      uc(i,j,0)    = -uc(i,j,1)      ! no-slip
      vc(i,j,0)    = (2.*bulk_v_sup)-vc(i,j,1)      ! no-slip
      wc(i,j,0)    = 0.             ! no-penetration
      uc(i,j,k1c)   = -uc(i,j,kmaxc)   ! no-slip
      vc(i,j,k1c)   = (2.*bulk_v_sup)-vc(i,j,kmaxc)    ! free-slip
      wc(i,j,kmaxc) = 0.             ! no-penetration
      wc(i,j,k1c)   = wc(i,j,kmaxc-1)  ! dw/dz=0 at wall (not used) 
    

   enddo
  enddo
!$omp end parallel
else
!$omp parallel default(none) shared(u,v,w) private(i,j)
!$omp do
  do j=0,j1c
    do i=0,i1c
      uc(i,j,0)    = -uc(i,j,1)     ! no-slip
      vc(i,j,0)    = -vc(i,j,1)     ! no-slip
      wc(i,j,0)    = 0.            ! no-penetration
      uc(i,j,k1c)   = -uc(i,j,kmaxc) +(2.*ub(i,j))!overlap boundary condition 
      vc(i,j,k1c)   = -vc(i,j,kmaxc) +(2.*vb(i,j))!overlap boundary condition 
      wc(i,j,kmaxc) = wb(i,j)                    !overlap boundary condition

      
      uc(i,j,k1c+1)=-uc(i,j,kmaxc-1)+(2.*ub(i,j)) !interpolation for second ghost cell (used for Eulr to Lagr)
      vc(i,j,k1c+1)=-vc(i,j,kmaxc-1)+(2.*vb(i,j)) !interpolation for second ghost cell (used for Eulr to Lagr)
      wc(i,j,k1c)=-wc(i,j,kmaxc-1)+(2.*wb(i,j)) !interpolation for second ghost cell (used for Eulr to Lagr)
      wc(i,j,k1c+1)=-wc(i,j,kmaxc-2)+(2.*wb(i,j)) !interpolation for second ghost cell (used for Eulr to Lagr)


    enddo
  enddo
!$omp end parallel
endif
case('cou')
!$omp parallel default(none) shared(u,v,w) private(i,j)
!$omp do

  do j=0,j1c
    do i=0,i1c
      uc(i,j,0) = -uc(i,j,1) ! no-slip
      vc(i,j,0) = lz*shear - vc(i,j,1) ! no-slip
      wc(i,j,0) = 0. ! no-penetration
      uc(i,j,k1c) = -uc(i,j,kmaxc) ! no-slip
      vc(i,j,k1c) = -lz*shear - vc(i,j,kmaxc) ! no-slip
      wc(i,j,kmaxc) = 0. ! no-penetration
      wc(i,j,k1c) = wc(i,j,kmaxc-1) ! dw/dz=0 at wall (not used)
    enddo
  enddo

!$omp end parallel
end select
!
! communicate data in x direction (periodic b.c.'s incorporated)
!
call updthalosc(uc,1)
call updthalosc(vc,1)
call updthalosc(wc,1)
!
! communicate data in y direction (periodic b.c.'s incorporated)
!
call updthalosc(uc,2)
call updthalosc(vc,2)
call updthalosc(wc,2)
!
return
end subroutine bounduvwc
!
subroutine boundpc(pc)
use mod_param
implicit none
integer :: i,j
real, dimension(0:,0:,0:) :: pc
!
!$omp parallel default(none) shared(p) private(i,j)
!$omp do
do j=0,j1c
  do i=0,i1c
    pc(i,j,0) = pc(i,j,1)     ! Newmann (consistent with no/free-slip)
    pc(i,j,k1c) = pc(i,j,kmaxc)  !overlap boundary condition
  enddo
enddo
!$omp end parallel
!
call updthalosc(pc,1)
!
!  communicate data in y direction (periodic b.c.'s incorporated)
!
call updthalosc(pc,2)
!
return
end subroutine boundpc

subroutine boundpcorrc(pc)
use mod_param
implicit none
integer :: i,j
real, dimension(0:,0:,0:) :: pc
!
!$omp parallel default(none) shared(p) private(i,j)
!$omp do
do j=0,j1c
  do i=0,i1c
    pc(i,j,0) = pc(i,j,1)     ! Newmann (consistent with no/free-slip)
    pc(i,j,k1c) = pc(i,j,kmaxc) !overlap boundary condition
  enddo
enddo
!$omp end parallel
!
call updthalosc(pc,1)
!
!  communicate data in y direction (periodic b.c.'s incorporated)
!
call updthalosc(pc,2)
!
return
end subroutine boundpcorrc

subroutine boundviscc(viscc)
use mod_param
implicit none
integer :: i,j
real, dimension(0:,0:,0:) :: viscc
!
!$omp parallel default(none) shared(p) private(i,j)
!$omp do
do j=0,j1c
  do i=0,i1c
    viscc(i,j,0) = viscc(i,j,1)     ! Newmann (consistent with no/free-slip)
    viscc(i,j,k1c) = viscc(i,j,kmaxc) !overlap boundary condition
  enddo
enddo
!$omp end parallel
!
call updthalosc(viscc,1)
!
!  communicate data in y direction (periodic b.c.'s incorporated)
!
call updthalosc(viscc,2)
!
return
end subroutine boundviscc
subroutine boundindch(indich)
use mod_param
implicit none
integer :: i,j
real, dimension(0:,0:,0:) :: indich
!
!$omp parallel default(none) shared(p) private(i,j)
!$omp do
do j=0,j1c
  do i=0,i1c
    indich(i,j,0) = indich(i,j,1)     ! Newmann
    indich(i,j,k1c) = indich(i,j,kmaxc)  !overlap boundary condition
  enddo
enddo
!$omp end parallel
!
call updthalosc(indich,1)
!
!  communicate data in y direction (periodic b.c.'s incorporated)
!
call updthalosc(indich,2)
!
return
end subroutine boundindch
!
subroutine updthalosc(var,dir)
use mpi
use mod_param
use mod_common_mpi
implicit none
real, dimension(0:,0:,0:), intent(inout) :: var
integer, intent(in) :: dir
integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
!
!  This subroutine updates the halos that store info
!  from the neighboring computational sub-domain
!
select case(dir)
case(1) ! x direction
  call MPI_SENDRECV(var(1,0,0),1,xhaloc,left,0,   &
                    var(i1c,0,0),1,xhaloc,right,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(imaxc,0,0),1,xhaloc,right,0, &
                    var(0,0,0),1,xhaloc,left,0,     &
                    comm_cart,status,error)
!call MPI_IRECV(var(0,0,0),1,xhalo,left,1, &
!               comm_cart,requests(2),error)
!call MPI_IRECV(var(i1c,0,0),1,xhalo,right,0, &
!               comm_cart,requests(1),error)
!call MPI_ISSEND(var(imaxc,0,0),1,xhalo,right,1, &
!               comm_cart,requests(4),error)
!call MPI_ISSEND(var(1,0,0),1,xhalo,left,0, &
!               comm_cart,requests(3),error)
!call MPI_WAITALL(4, requests, statuses, error)
case(2) ! y direction
  call MPI_SENDRECV(var(0,1,0),1,yhaloc,front,0, &
                    var(0,j1c,0),1,yhaloc,back,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(0,jmaxc,0),1,yhaloc,back,0, &
                    var(0,0,0),1,yhaloc,front,0,   &
                    comm_cart,status,error)
!call MPI_IRECV(var(0,j1c,0),1,yhalo,back,0, &
!               comm_cart,requests(1),error)
!call MPI_IRECV(var(0,0,0),1,yhalo,front,1, &
!               comm_cart,requests(2),error)
!call MPI_ISSEND(var(0,1,0),1,yhalo,front,0, &
!               comm_cart,requests(3),error)
!call MPI_ISSEND(var(0,jmaxc,0),1,yhalo,back,1, &
!               comm_cart,requests(4),error)
!call MPI_WAITALL(4, requests, statuses, error)
end select
!
return
end subroutine updthalosc
!
end module mod_boundc
