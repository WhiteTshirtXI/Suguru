module mod_crnkc
use mod_momc
use mod_boundc
use mod_common
use mod_common_mpi
implicit none
private
public crnkc
contains
!
subroutine crnkc(durk1c,dvrk1c,dwrk1c)
! 
!
implicit none
integer i,j,k
real, intent(inout) :: durk1c(0:,0:,0:)
real, intent(inout) :: dvrk1c(0:,0:,0:)
real, intent(inout) :: dwrk1c(0:,0:,0:)
real, dimension(0:i1c,0:j1c,0:k1c):: cupc,cuec,cuwc,cunc,cusc,cutc,cubc
real, dimension(0:i1c,0:j1c,0:k1c):: cvpc,cvec,cvwc,cvnc,cvsc,cvtc,cvbc
real, dimension(0:i1c,0:j1c,0:k1c):: cwpc,cwec,cwwc,cwnc,cwsc,cwtc,cwbc
real :: wallshear_all,errorx1,errory1,errorz1
real :: unorm,vnorm,wnorm,unorm1,vnorm1,wnorm1
real :: factor
real :: dVeul
!
factor = dt
!
if (ar==1) then
call momxpc(durk1c,pnewc)
call momypc(dvrk1c,pnewc)
call momzpc(dwrk1c,pnewc)
!$omp parallel default(shared) &
!$omp&private(i,j,k)
!$omp do
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
      dudtc(i,j,k) = uoc(i,j,k) + factor*durk1c(i,j,k)
      dvdtc(i,j,k) = voc(i,j,k) + factor*dvrk1c(i,j,k)
      dwdtc(i,j,k) = woc(i,j,k) + factor*dwrk1c(i,j,k)
    enddo
  enddo
enddo
!$omp end parallel
call momxadc(durk1c,uoc,voc,woc,uooc,vooc,wooc)
call momyadc(dvrk1c,uoc,voc,woc,uooc,vooc,wooc)
call momzadc(dwrk1c,uoc,voc,woc,uooc,vooc,wooc)
!$omp parallel default(shared) &
!$omp&private(i,j,k)
!$omp do
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
      dudtc(i,j,k) = dudtc(i,j,k) + factor*durk1c(i,j,k)
      dvdtc(i,j,k) = dvdtc(i,j,k) + factor*dvrk1c(i,j,k)
      dwdtc(i,j,k) = dwdtc(i,j,k) + factor*dwrk1c(i,j,k)
    enddo
  enddo
enddo
!$omp end parallel
call momxdifc(durk1c,uoc,voc,woc,cupc,cuec,cuwc,cunc,cusc,cutc,cubc)
call momydifc(dvrk1c,uoc,voc,woc,cvpc,cvec,cvwc,cvnc,cvsc,cvtc,cvbc)
call momzdifc(dwrk1c,uoc,voc,woc,cwpc,cwec,cwwc,cwnc,cwsc,cwtc,cwbc)
!$omp parallel default(shared) &
!$omp&private(i,j,k)
!$omp do
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc

      dudtc(i,j,k) = dudtc(i,j,k) + factor*.5*durk1c(i,j,k)
      dvdtc(i,j,k) = dvdtc(i,j,k) + factor*.5*dvrk1c(i,j,k)
      dwdtc(i,j,k) = dwdtc(i,j,k) + factor*.5*dwrk1c(i,j,k)
    enddo
  enddo
enddo

endif

do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
      unewc(i,j,k) = (dudtc(i,j,k)-(cuec(i,j,k)*unewc(i+1,j,k))-(cuwc(i,j,k)*unewc(i-1,j,k))-(cunc(i,j,k)*unewc(i,j+1,k))-(cusc(i,j,k)*unewc(i,j-1,k))- &
      (cutc(i,j,k)*unewc(i,j,k+1))-(cubc(i,j,k)*unewc(i,j,k-1)))/cupc(i,j,k)
      vnewc(i,j,k) = (dvdtc(i,j,k)-(cvec(i,j,k)*vnewc(i+1,j,k))-(cvwc(i,j,k)*vnewc(i-1,j,k))-(cvnc(i,j,k)*vnewc(i,j+1,k))-(cvsc(i,j,k)*vnewc(i,j-1,k))- &
      (cvtc(i,j,k)*vnewc(i,j,k+1))-(cvbc(i,j,k)*vnewc(i,j,k-1)))/cvpc(i,j,k)
    enddo
  enddo
enddo

do k=1,kmaxc-1
  do j=1,jmaxc
    do i=1,imaxc
      wnewc(i,j,k) = (dwdtc(i,j,k)-(cwec(i,j,k)*wnewc(i+1,j,k))-(cwwc(i,j,k)*wnewc(i-1,j,k))-(cwnc(i,j,k)*wnewc(i,j+1,k))-(cwsc(i,j,k)*wnewc(i,j-1,k))- &
      (cwtc(i,j,k)*wnewc(i,j,k+1))-(cwbc(i,j,k)*wnewc(i,j,k-1)))/cwpc(i,j,k)
    enddo
  enddo
enddo
!***********************************************************SOR*********************************************************************************
!do k=1,kmaxc-1
!  do j=1,jmaxc
!    do i=1,imaxc
!    unewc(i,j,k)=(unewc(i,j,k)*1.2)-(upic(i,j,k)*.2)
!    vnewc(i,j,k)=(vnewc(i,j,k)*1.2)-(vpic(i,j,k)*.2)
!    wnewc(i,j,k)=(wnewc(i,j,k)*1.2)-(wpic(i,j,k)*.2)
!    enddo
!  enddo
!enddo
!***********************************************************************************************************************************************

call bounduvwc(unewc,vnewc,wnewc)


errorx1=0.
errory1=0.
errorz1=0.
unorm1=0.
vnorm1=0.
wnorm1=0.
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
      errorx1=max(errorx1,abs(unewc(i,j,k)-upic(i,j,k)))
      errory1=max(errory1,abs(vnewc(i,j,k)-vpic(i,j,k)))
      errorz1=max(errorz1,abs(wnewc(i,j,k)-wpic(i,j,k)))
      unorm1=unorm1+abs(unewc(i,j,k))
      vnorm1=vnorm1+abs(vnewc(i,j,k))
      wnorm1=wnorm1+abs(wnewc(i,j,k))
    enddo
  enddo
enddo
call mpi_allreduce(errorx1,errorx,1,mpi_real8,mpi_max,comm_cart,error)
call mpi_allreduce(errory1,errory,1,mpi_real8,mpi_max,comm_cart,error)
call mpi_allreduce(errorz1,errorz,1,mpi_real8,mpi_max,comm_cart,error)
call mpi_allreduce(unorm1,unorm,1,mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(vnorm1,vnorm,1,mpi_real8,mpi_sum,comm_cart,error)
call mpi_allreduce(wnorm1,wnorm,1,mpi_real8,mpi_sum,comm_cart,error)
unorm=unorm/(itotc*jtotc*ktotc)
vnorm=vnorm/(itotc*jtotc*ktotc)
wnorm=wnorm/(itotc*jtotc*ktotc)
errorx=errorx/unorm
errory=errory/vnorm
errorz=errorz/wnorm
!!
return
end subroutine crnkc

!
end module mod_crnkc
