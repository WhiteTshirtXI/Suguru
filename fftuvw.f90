module mod_fftuvw
use mod_mom
use mod_bound
use mod_common
use mod_common_mpi
use mod_initsolver
use mod_solver
implicit none
private
public fftuvw
contains
!
subroutine fftuvw(durk1,dvrk1,dwrk1)
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

call momxp(durk1,pnew)
call momyp(dvrk1,pnew)
call momzp(dwrk1,pnew)
!$omp parallel default(shared) &
!$omp&private(i,j,k)
!$omp do
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      dudt(i,j,k) =  -durk1(i,j,k)
      dvdt(i,j,k) =  -dvrk1(i,j,k)
      dwdt(i,j,k) =  -dwrk1(i,j,k)
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
      dudt(i,j,k) = dudt(i,j,k) - durk1(i,j,k)
      dvdt(i,j,k) = dvdt(i,j,k) - dvrk1(i,j,k)
      dwdt(i,j,k) = dwdt(i,j,k) - dwrk1(i,j,k)
    enddo
  enddo
enddo
!$omp end parallel
call momzad(dwrk1,uo,vo,wo,uoo,voo,woo)
!$omp parallel default(shared) &
!$omp&private(i,j,k)
!$omp do
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      dudt(i,j,k) = dudt(i,j,k) - (uo(i,j,k)/dt)
      dvdt(i,j,k) = dvdt(i,j,k) - (vo(i,j,k)/dt)
      dwdt(i,j,k) = dwdt(i,j,k) - (wo(i,j,k)/dt)
    enddo
  enddo
enddo
!$omp end parallel

dudt=dudt*Rep
dvdt=dvdt*Rep
dwdt=dwdt*Rep
!******************************************************setting coefficients for w***************************************************************
call initsolver(3)
!******************************************************SOurce term of w equation****************************************************************
!$omp parallel default(shared) &
!$omp&private(i,j,k)
!$omp do
  do j=1,jmax
    do i=1,imax
      dwdt(i,j,kmax) = 0.
    enddo
  enddo
!$omp end parallel
!********************************************************Sloving w equation*********************************************************************
wfft=dwdt
call solver1d(wfft,3)
!$omp parallel default(shared) &
!$omp&private(i,j,k)
!$omp do
do k=1,kmax
  do j=1,jmax
    do i=1,imax
    wnew(i,j,k)=wfft(i,j,k)
    enddo
  enddo
enddo
!************************************************setting coefficients for u and v***************************************************************
call initsolver(2)
!******************************no extra source term needed for u and v equations in main mesh***************************************************

!*************************************************Sloving u and v equations*********************************************************************
ufft=dudt
vfft=dvdt
call solver1d(ufft,2)
call solver1d(vfft,2)
do k=1,kmax
  do j=1,jmax
    do i=1,imax
    unew(i,j,k)=ufft(i,j,k)
    vnew(i,j,k)=vfft(i,j,k)
    enddo
  enddo
enddo
!***********************************************************************************************************************************************
!!
return
end subroutine fftuvw

!
end module mod_fftuvw
