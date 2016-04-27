module mod_correcc
use mod_param
use mod_common
implicit none
private
public correcc
contains
!
! Corrects the velocity so that it is divergence free
!
subroutine correcc(pc)
real, intent(in), dimension(0:,0:,0:) :: pc
real :: factor,factori,factorj,factork
integer :: i,j,k
!
factor = dt
factori = factor*dxic
factorj = factor*dyic
factork = factor*dzic
!
!$omp parallel default(none) &
!$omp&shared(factori,factorj,factork) &
!$omp&shared(p,unew,vnew,wnew,dudt,dvdt,dwdt) &
!$omp&private(i,j,k)
!$omp do
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
      unewc(i,j,k) = unewc(i,j,k) - factori * ( pc(i+1,j,k)-pc(i,j,k) )
      vnewc(i,j,k) = vnewc(i,j,k) - factorj * ( pc(i,j+1,k)-pc(i,j,k) )
    enddo
  enddo
enddo
do k=1,kmaxc-1
  do j=1,jmaxc
    do i=1,imaxc
      wnewc(i,j,k) = wnewc(i,j,k) - factork * ( pc(i,j,k+1)-pc(i,j,k) )
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine correcc
!
end module mod_correcc
