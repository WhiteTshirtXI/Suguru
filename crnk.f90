module mod_crnk
use mod_mom
use mod_bound
use mod_common
use mod_common_mpi
implicit none
private
public crnk
contains
!
subroutine crnk(durk1,dvrk1,dwrk1)
! 
!
implicit none
integer i,j,k
real, intent(inout) :: durk1(0:,0:,0:)
real, intent(inout) :: dvrk1(0:,0:,0:)
real, intent(inout) :: dwrk1(0:,0:,0:)
real, dimension(0:i1,0:j1,0:k1):: cup,cue,cuw,cun,cus,cut,cub
real, dimension(0:i1,0:j1,0:k1):: cvp,cve,cvw,cvn,cvs,cvt,cvb
real, dimension(0:i1,0:j1,0:k1):: cwp,cwe,cww,cwn,cws,cwt,cwb
real :: wallshear_all,errorx1,errory1,errorz1
real :: unorm,vnorm,wnorm,unorm1,vnorm1,wnorm1
real :: factor
real :: dVeul
!
factor = dt
!
if (ar==1) then
call momxp(durk1,pnew)
call momyp(dvrk1,pnew)
call momzp(dwrk1,pnew)
!$omp parallel default(shared) &
!$omp&private(i,j,k)
!$omp do
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      dudt(i,j,k) = uo(i,j,k) + factor*durk1(i,j,k)
      dvdt(i,j,k) = vo(i,j,k) + factor*dvrk1(i,j,k)
      dwdt(i,j,k) = wo(i,j,k) + factor*dwrk1(i,j,k)
    enddo
  enddo
enddo
!$omp end parallel
call momxad(durk1,uo,vo,wo,uoo,voo,woo)
call momyad(dvrk1,uo,vo,wo,uoo,voo,woo)
call momzad(dwrk1,uo,vo,wo,uoo,voo,woo)
!$omp parallel default(shared) &
!$omp&private(i,j,k)
!$omp do
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      dudt(i,j,k) = dudt(i,j,k) + factor*durk1(i,j,k)
      dvdt(i,j,k) = dvdt(i,j,k) + factor*dvrk1(i,j,k)
      dwdt(i,j,k) = dwdt(i,j,k) + factor*dwrk1(i,j,k)
    enddo
  enddo
enddo
!$omp end parallel
call momxdif(durk1,uo,vo,wo,cup,cue,cuw,cun,cus,cut,cub)
call momydif(dvrk1,uo,vo,wo,cvp,cve,cvw,cvn,cvs,cvt,cvb)
call momzdif(dwrk1,uo,vo,wo,cwp,cwe,cww,cwn,cws,cwt,cwb)
!$omp parallel default(shared) &
!$omp&private(i,j,k)
!$omp do
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      dudt(i,j,k) = dudt(i,j,k) + factor*.5*durk1(i,j,k)
      dvdt(i,j,k) = dvdt(i,j,k) + factor*.5*dvrk1(i,j,k)
      dwdt(i,j,k) = dwdt(i,j,k) + factor*.5*dwrk1(i,j,k)
    enddo
  enddo
enddo

endif


do k=1,kmax
  do j=1,jmax
    do i=1,imax
      unew(i,j,k) = (dudt(i,j,k)-(cue(i,j,k)*unew(i+1,j,k))-(cuw(i,j,k)*unew(i-1,j,k))-(cun(i,j,k)*unew(i,j+1,k))-(cus(i,j,k)*unew(i,j-1,k))- &
      (cut(i,j,k)*unew(i,j,k+1))-(cub(i,j,k)*unew(i,j,k-1)))/cup(i,j,k)
      vnew(i,j,k) = (dvdt(i,j,k)-(cve(i,j,k)*vnew(i+1,j,k))-(cvw(i,j,k)*vnew(i-1,j,k))-(cvn(i,j,k)*vnew(i,j+1,k))-(cvs(i,j,k)*vnew(i,j-1,k))- &
      (cvt(i,j,k)*vnew(i,j,k+1))-(cvb(i,j,k)*vnew(i,j,k-1)))/cvp(i,j,k)
    enddo
  enddo
enddo


do k=1,kmax-1
  do j=1,jmax
    do i=1,imax
      wnew(i,j,k) = (dwdt(i,j,k)-(cwe(i,j,k)*wnew(i+1,j,k))-(cww(i,j,k)*wnew(i-1,j,k))-(cwn(i,j,k)*wnew(i,j+1,k))-(cws(i,j,k)*wnew(i,j-1,k))- &
      (cwt(i,j,k)*wnew(i,j,k+1))-(cwb(i,j,k)*wnew(i,j,k-1)))/cwp(i,j,k)
    enddo
  enddo
enddo
!***********************************************************SOR*********************************************************************************
!do k=1,kmax-1
!  do j=1,jmax
!    do i=1,imax
!    unew(i,j,k)=(unew(i,j,k)*1.2)-(upi(i,j,k)*.2)
!    vnew(i,j,k)=(vnew(i,j,k)*1.2)-(vpi(i,j,k)*.2)
!    wnew(i,j,k)=(wnew(i,j,k)*1.2)-(wpi(i,j,k)*.2)
!    enddo
!  enddo
!enddo
!***********************************************************************************************************************************************
call bounduvw(unew,vnew,wnew)
errorx1=0.
errory1=0.
errorz1=0.
unorm1=0.
vnorm1=0.
wnorm1=0.
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      errorx1=max(errorx1,abs(unew(i,j,k)-upi(i,j,k)))
      errory1=max(errory1,abs(vnew(i,j,k)-vpi(i,j,k)))
      errorz1=max(errorz1,abs(wnew(i,j,k)-wpi(i,j,k)))
      unorm1=unorm1+abs(unew(i,j,k))
      vnorm1=vnorm1+abs(vnew(i,j,k))
      wnorm1=wnorm1+abs(wnew(i,j,k))
    enddo
  enddo
enddo
call mpi_allreduce(errorx1,errorx,1,mpi_real8,mpi_max,comm_cart,error)
call mpi_allreduce(errory1,errory,1,mpi_real8,mpi_max,comm_cart,error)
call mpi_allreduce(errorz1,errorz,1,mpi_real8,mpi_max,comm_cart,error)
call mpi_allreduce(unorm1,unorm,1,mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(vnorm1,vnorm,1,mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(wnorm1,wnorm,1,mpi_real8,mpi_sum,comm_cart,error)
unorm=unorm/(itot*jtot*ktot)
vnorm=vnorm/(itot*jtot*ktot)
wnorm=wnorm/(itot*jtot*ktot)
errorx=errorx/unorm
errory=errory/vnorm
errorz=errorz/wnorm

!!
return
end subroutine crnk

!
end module mod_crnk
