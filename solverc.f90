!******************************************************************
!********************  DIRECT POISSON SOLVER **********************
!*****                                                        *****
!***               p_xx + p_yy + p_zz  = f(x,y,z)               ***
!****                                                         *****
!******************************************************************
!   Fourier Transforms of the Poisson equation in the x and y
!   directions yield:
!
!   a^2 p' + b^2 p' + p'_zz = f'(kx,ky,z) = FFT_j[ FTT_i[ f(x,y,z) ] ] ,
!
!   where a and b are the known eigenvalues, kx and ky are the
!   wavenumbers in the x and y direction, and p'_zz is given by:
!
!   p'_zz =[ p'_{kx,ky,k+1} - 2*p'_{kx,ky,k} + p'_{kx,ky,k-1} ]/(dz*dz) .
!
!   The equation above results in a tridiagonal system in k:
!
!   a^2 p' + b^2 p' + p'_zz =
!   [p'_{kx,ky,k+1}-(2+a^2+b^2)*p'_{kx,ky,k}+p'_{kx,ky,k-1}]/(dz*dz)=f'(kx,ky,z) .
!
!   This can be solved with Gaussian elimination. The pressure p in
!   physical space is obtained from 2 inverse Fourier Transforms
!   according to: p = iFFT_j[ iFFT_i[ p' ] ].
!******************************************************************
!******************************************************************
!******************************************************************
!****   Programmers: Bendiks Jan Boersma                     ******
!****                Wim-Paul Breugem (1D parallellisation)  ******
!****                Pedro Costa (2D parallelisation)        ******
!****   email      : w.p.breugem@tudelft.nl                  ******
!****   USES       : VFFTPACK   (netlib)                     ******
!****                2DECOMP&FFT (www.2decomp.org)           ******
!******************************************************************
!******************************************************************
module mod_solverc
use mod_param
use mod_common
implicit none
private
public solver1dc
contains

subroutine solver1dc(pc,nec)
use mod_zredistributec
use mod_zredistribute
implicit none
integer,intent(in)::nec
integer, parameter :: nprocs = dims(1)*dims(2)
integer, parameter :: ksolc = kmaxc/nprocs
real, intent(inout) :: pc(0:,0:,0:)
real,dimension(itotc,jtotc,ksolc) :: p2c
real :: bbc
real :: zc,dc(imaxc,jmaxc,kmaxc)
real :: dic(itotc),djc(jtotc)
integer i,j,k
!
! Redistribute the pressure so that it is only distributed
! in the z-direction (starting from 0 to nprocs)
!

call zredistributec(pc,p2c,0)

!$omp parallel default(none) &
!$omp&shared(p2) &
!$omp&private(i,j,k,di,dj) &
!$omp&firstprivate(wi,wj)
!$omp do
do k=1,ksolc
  !  FFT  ---> I direction
  do j=1,jtotc
    call vrfftf(1,itotc,p2c(1:itotc,j,k),dic,1,wic)
  enddo
  !  FFT  ---> J direction
  do i=1,itotc
    call vrfftf(1,jtotc,p2c(i,1:jtotc,k),djc,1,wjc)
  enddo
enddo
!$omp end parallel
!
! Redistribute the pressure so that it is distributed
! on the cartesian grid
!

call zredistributec(pc,p2c,1)

!$omp parallel default(none) &
!$omp&shared(a,p,b,c,d,xyrtc) &
!$omp&private(i,j,k,zc,bbc)
k=1
!$omp do
do j=1,jmaxc
  do i=1,imaxc
    zc        = 1./(bc(1)+xyrtc(i,j))
    dc(i,j,1) = cc(1)*zc
    pc(i,j,1) = pc(i,j,1)*zc
  enddo
enddo
!$omp barrier
do k=2,kmaxc-1
!$omp do
  do j=1,jmaxc
    do i=1,imaxc
      bbc       = bc(k)+xyrtc(i,j)
      zc        = 1./(bbc-ac(k)*dc(i,j,k-1))
      dc(i,j,k) = cc(k)*zc
      pc(i,j,k) = (pc(i,j,k)-ac(k)*pc(i,j,k-1))*zc
    enddo
  enddo
!$omp barrier
enddo
k = kmaxc
!$omp do
do j=1,jmaxc
  do i=1,imaxc
    if (nec==3) then
    bbc       = bc(kmaxc)
    else
    bbc       = bc(kmaxc)+xyrtc(i,j)
    endif
    zc        = bbc-ac(kmaxc)*dc(i,j,kmaxc-1)
    if(zc.ne.0.) then
      pc(i,j,kmaxc) = (pc(i,j,kmaxc)-ac(kmaxc)*pc(i,j,kmaxc-1))/zc
    else
      pc(i,j,kmaxc) =0.
    endif
  enddo
enddo
do k=kmaxc-1,1,-1
!$omp do
  do j=1,jmaxc
    do i=1,imaxc
      pc(i,j,k) = pc(i,j,k)-dc(i,j,k)*pc(i,j,k+1)
    enddo
  enddo
!$omp barrier
enddo
!$omp end parallel
!
call zredistributec(pc,p2c,0)
!
!$omp parallel default(none) &
!$omp&shared(p2) &
!$omp&private(i,j,k,di,dj) &
!$omp&firstprivate(wi,wj)
!$omp do
do k=1,ksolc
  ! BACKWARD FFT ---> J direction
  do i=1,itotc
    call vrfftb(1,jtotc,p2c(i,1:jtotc,k),djc,1,wjc)
  enddo
  ! BACKWARD FFT ---> I direction
  do j=1,jtotc
    call vrfftb(1,itotc,p2c(1:itotc,j,k),dic,1,wic)
  enddo
enddo
!$omp end parallel

call zredistributec(pc,p2c,1)
!
return
end subroutine solver1dc
end module mod_solverc
