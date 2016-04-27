module mod_fftuvwc
use mod_momc
use mod_boundc
use mod_common
use mod_initsolver
use mod_solverc
use mod_common_mpi
implicit none
private
public fftuvwc
contains
!
subroutine fftuvwc(durk1c,dvrk1c,dwrk1c)
! 
!
implicit none
integer i,j,k
real, intent(inout) :: durk1c(0:,0:,0:)
real, intent(inout) :: dvrk1c(0:,0:,0:)
real, intent(inout) :: dwrk1c(0:,0:,0:)
real :: wallshear_all,errorx1,errory1,errorz1
real :: unorm,vnorm,wnorm,unorm1,vnorm1,wnorm1
real :: factor
real :: dVeul
!
factor = dt
!

call momxpc(durk1c,pnewc)
call momypc(dvrk1c,pnewc)
call momzpc(dwrk1c,pnewc)
!$omp parallel default(shared) &
!$omp&private(i,j,k)
!$omp do
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
      dudtc(i,j,k) = -durk1c(i,j,k)
      dvdtc(i,j,k) = -dvrk1c(i,j,k)
      dwdtc(i,j,k) = -dwrk1c(i,j,k)
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
      dudtc(i,j,k) = dudtc(i,j,k) - durk1c(i,j,k)
      dvdtc(i,j,k) = dvdtc(i,j,k) - dvrk1c(i,j,k)
      dwdtc(i,j,k) = dwdtc(i,j,k) - dwrk1c(i,j,k)
    enddo
  enddo
enddo
!$omp end parallel
!$omp parallel default(shared) &
!$omp&private(i,j,k)
!$omp do
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
      dudtc(i,j,k) = dudtc(i,j,k) - (uoc(i,j,k)/dt)
      dvdtc(i,j,k) = dvdtc(i,j,k) - (voc(i,j,k)/dt)
      dwdtc(i,j,k) = dwdtc(i,j,k) - (woc(i,j,k)/dt)
    enddo
  enddo
enddo
!$omp end parallel

!******************************************************setting coefficients for w***************************************************************
call initsolver(3)
!******************************************************SOurce term of w equation****************************************************************
dwdtc=dwdtc*Rep
!$omp parallel default(shared) &
!$omp&private(i,j,k)
!$omp do
  do j=1,jmaxc
    do i=1,imaxc
      dwdtc(i,j,kmaxc) = wb(i,j)
    enddo
  enddo
!$omp end parallel
!********************************************************Sloving w equation*********************************************************************
wfftc=dwdtc
call solver1dc(wfftc,3)
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
    wnewc(i,j,k)=wfftc(i,j,k)
    enddo
  enddo
enddo
!************************************************setting coefficients for u and v***************************************************************
call initsolver(2)
!***********************************************SOurce term of u and v equations****************************************************************
dudtc=dudtc*Rep
dvdtc=dvdtc*Rep
!$omp parallel default(shared) &
!$omp&private(i,j,k)
!$omp do
do j=1,jmaxc
  do i=1,imaxc
  dudtc(i,j,kmaxc)=dudtc(i,j,kmaxc)-((2.*ub(i,j))/(dzc**2.))
  dvdtc(i,j,kmaxc)=dvdtc(i,j,kmaxc)-((2.*vb(i,j))/(dzc**2.))
  enddo
enddo
!$omp end parallel
!*************************************************Sloving u and v equations*********************************************************************
ufftc=dudtc
vfftc=dvdtc
call solver1dc(ufftc,2)
call solver1dc(vfftc,2)
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
    unewc(i,j,k)=ufftc(i,j,k)
    vnewc(i,j,k)=vfftc(i,j,k)
    enddo
  enddo
enddo
!***********************************************************************************************************************************************

!!
return
end subroutine fftuvwc

!
end module mod_fftuvwc
