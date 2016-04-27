module mod_bound
implicit none
private
public bounduvw, boundp, updthalos, boundind
contains
!
subroutine bounduvw(u,v,w)
use mod_param
use mod_common
implicit none
integer :: i,j
real, dimension(0:,0:,0:) :: u,v,w
!
select case(iniu)
case('poi','zer','log')
if(isfreeslip) then
!$omp parallel default(none) shared(u,v,w) private(i,j)
!$omp do
  do j=0,j1
    do i=0,i1
      



      u(i,j,0)    = -u(i,j,1)      ! no-slip
      v(i,j,0)    = (2.*bulk_v_sup)-v(i,j,1)      ! no-slip
      w(i,j,0)    = 0.             ! no-penetration
      u(i,j,k1)   = -u(i,j,kmax)   ! no-slip
      v(i,j,k1)   = (2.*bulk_v_sup)-v(i,j,kmax)    ! free-slip
      w(i,j,kmax) = 0.             ! no-penetration
      w(i,j,k1)   = w(i,j,kmax-1)  ! dw/dz=0 at wall (not used) 
    

   enddo
  enddo
!$omp end parallel
else
!$omp parallel default(none) shared(u,v,w) private(i,j)
!$omp do
  do j=0,j1
    do i=0,i1
      u(i,j,0)    = -u(i,j,1)     ! no-slip
      v(i,j,0)    = -v(i,j,1)     ! no-slip
      w(i,j,0)    = 0.            ! no-penetration
      u(i,j,k1)   = -u(i,j,kmax)  ! no-slip
      v(i,j,k1)   = -v(i,j,kmax)  ! no-slip
      w(i,j,kmax) = 0.            ! no-penetration
      w(i,j,k1)   = w(i,j,kmax-1) ! dw/dz=0 at wall (not used) 
    enddo
  enddo
!$omp end parallel
endif
case('cou')
!$omp parallel default(none) shared(u,v,w) private(i,j)
!$omp do

  do j=0,j1
    do i=0,i1
      u(i,j,0) = -u(i,j,1) ! no-slip
      v(i,j,0) = lz*shear - v(i,j,1) ! no-slip
      w(i,j,0) = 0. ! no-penetration
      u(i,j,k1) = -u(i,j,kmax) ! no-slip
      v(i,j,k1) = -lz*shear - v(i,j,kmax) ! no-slip
      w(i,j,kmax) = 0. ! no-penetration
      w(i,j,k1) = w(i,j,kmax-1) ! dw/dz=0 at wall (not used)
    enddo
  enddo

!  write(6,*) 'v(i,j,k1)',v(10,10,k1)
!$omp end parallel
end select
!
! communicate data in x direction (periodic b.c.'s incorporated)
!
call updthalos(u,1)
call updthalos(v,1)
call updthalos(w,1)
!
! communicate data in y direction (periodic b.c.'s incorporated)
!
call updthalos(u,2)
call updthalos(v,2)
call updthalos(w,2)
!
return
end subroutine bounduvw
!
subroutine boundp(p)
use mod_param
implicit none
integer :: i,j
real, dimension(0:,0:,0:) :: p
!
!$omp parallel default(none) shared(p) private(i,j)
!$omp do
do j=0,j1
  do i=0,i1
    p(i,j,0) = p(i,j,1)     ! Newmann (consistent with no/free-slip)
    p(i,j,k1) = p(i,j,kmax) ! Newmann (consistent with no/free-slip)
  enddo
enddo
!$omp end parallel
!
call updthalos(p,1)
!
!  communicate data in y direction (periodic b.c.'s incorporated)
!
call updthalos(p,2)
!
return
end subroutine boundp

subroutine boundind(indi)
use mod_param
implicit none
integer :: i,j
real, dimension(0:,0:,0:) :: indi
!
!$omp parallel default(none) shared(p) private(i,j)
!$omp do
do j=0,j1
  do i=0,i1
    indi(i,j,0) = indi(i,j,1)     ! Newmann
    indi(i,j,k1) = indi(i,j,kmax) ! Newmann
  enddo
enddo

!$omp end parallel
!
call updthalos(indi,1)

!
!  communicate data in y direction (periodic b.c.'s incorporated)
!
call updthalos(indi,2)

!
return
end subroutine boundind
!
subroutine updthalos(var,dir)
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
  call MPI_SENDRECV(var(1,0,0),1,xhalo,left,0,   &
                    var(i1,0,0),1,xhalo,right,0, &
                    comm_cart,status,error)

  call MPI_SENDRECV(var(imax,0,0),1,xhalo,right,0, &
                    var(0,0,0),1,xhalo,left,0,     &
                    comm_cart,status,error)
!call MPI_IRECV(var(0,0,0),1,xhalo,left,1, &
!               comm_cart,requests(2),error)
!call MPI_IRECV(var(i1,0,0),1,xhalo,right,0, &
!               comm_cart,requests(1),error)
!call MPI_ISSEND(var(imax,0,0),1,xhalo,right,1, &
!               comm_cart,requests(4),error)
!call MPI_ISSEND(var(1,0,0),1,xhalo,left,0, &
!               comm_cart,requests(3),error)
!call MPI_WAITALL(4, requests, statuses, error)
case(2) ! y direction
  call MPI_SENDRECV(var(0,1,0),1,yhalo,front,0, &
                    var(0,j1,0),1,yhalo,back,0, &
                    comm_cart,status,error)
  call MPI_SENDRECV(var(0,jmax,0),1,yhalo,back,0, &
                    var(0,0,0),1,yhalo,front,0,   &
                    comm_cart,status,error)
!call MPI_IRECV(var(0,j1,0),1,yhalo,back,0, &
!               comm_cart,requests(1),error)
!call MPI_IRECV(var(0,0,0),1,yhalo,front,1, &
!               comm_cart,requests(2),error)
!call MPI_ISSEND(var(0,1,0),1,yhalo,front,0, &
!               comm_cart,requests(3),error)
!call MPI_ISSEND(var(0,jmax,0),1,yhalo,back,1, &
!               comm_cart,requests(4),error)
!call MPI_WAITALL(4, requests, statuses, error)
end select
!
return
end subroutine updthalos
!
end module mod_bound
