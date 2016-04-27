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
module mod_solver
use mod_param
use mod_common
implicit none
private
public solver1d
contains
subroutine solver1d(p,ne)
use mod_zredistribute
implicit none
integer,intent(in)::ne
integer, parameter :: nprocs = dims(1)*dims(2)
integer, parameter :: ksol = kmax/nprocs
real, intent(inout) :: p(0:,0:,0:)
real,dimension(itot,jtot,ksol) :: p2
real :: bb
real :: z,d(imax,jmax,kmax)
real :: di(itot),dj(jtot)
integer i,j,k
!
! Redistribute the pressure so that it is only distributed
! in the z-direction (starting from 0 to nprocs)
!
call zredistribute(p,p2,0)
!$omp parallel default(none) &
!$omp&shared(p2) &
!$omp&private(i,j,k,di,dj) &
!$omp&firstprivate(wi,wj)
!$omp do
do k=1,ksol
  !  FFT  ---> I direction
  do j=1,jtot
    call vrfftf(1,itot,p2(1:itot,j,k),di,1,wi)
  enddo
  !  FFT  ---> J direction
  do i=1,itot
    call vrfftf(1,jtot,p2(i,1:jtot,k),dj,1,wj)
  enddo
enddo
!$omp end parallel
!
! Redistribute the pressure so that it is distributed
! on the cartesian grid
!
call zredistribute(p,p2,1)
!$omp parallel default(none) &
!$omp&shared(a,p,b,c,d,xyrt) &
!$omp&private(i,j,k,z,bb)
k=1
!$omp do
do j=1,jmax
  do i=1,imax
    z        = 1./(b(1)+xyrt(i,j))
    d(i,j,1) = c(1)*z
    p(i,j,1) = p(i,j,1)*z
  enddo
enddo
!$omp barrier
do k=2,kmax-1
!$omp do
  do j=1,jmax
    do i=1,imax
      bb       = b(k)+xyrt(i,j)
      z        = 1./(bb-a(k)*d(i,j,k-1))
      d(i,j,k) = c(k)*z
      p(i,j,k) = (p(i,j,k)-a(k)*p(i,j,k-1))*z
    enddo
  enddo
!$omp barrier
enddo
k = kmax
!$omp do
do j=1,jmax
  do i=1,imax
    if (ne==3) then
    bb       = b(kmax)
    else
    bb       = b(kmax)+xyrt(i,j)
    endif
    z        = bb-a(kmax)*d(i,j,kmax-1)
    if(z.ne.0.) then
      p(i,j,kmax) = (p(i,j,kmax)-a(kmax)*p(i,j,kmax-1))/z
    else
      p(i,j,kmax) =0.
    endif
  enddo
enddo
do k=kmax-1,1,-1
!$omp do
  do j=1,jmax
    do i=1,imax
      p(i,j,k) = p(i,j,k)-d(i,j,k)*p(i,j,k+1)
    enddo
  enddo
!$omp barrier
enddo
!$omp end parallel


!
call zredistribute(p,p2,0)
!
!$omp parallel default(none) &
!$omp&shared(p2) &
!$omp&private(i,j,k,di,dj) &
!$omp&firstprivate(wi,wj)
!$omp do
do k=1,ksol
  ! BACKWARD FFT ---> J direction
  do i=1,itot
    call vrfftb(1,jtot,p2(i,1:jtot,k),dj,1,wj)
  enddo
  ! BACKWARD FFT ---> I direction
  do j=1,jtot
    call vrfftb(1,itot,p2(1:itot,j,k),di,1,wi)
  enddo
enddo
!$omp end parallel
call zredistribute(p,p2,1)
!
return
end subroutine solver1d
end module mod_solver
