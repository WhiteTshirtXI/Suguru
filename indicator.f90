module mod_indicator
use mpi
use mod_common
use mod_common_mpi
use mod_interp_spread
use mod_interp_spreadc
use mod_param
use mod_sph
use mod_bound
use mod_boundc
use mod_initsolver
use mod_solver
use mod_solverc
use mod_interpolation
implicit none
private
public indicator
contains
subroutine indicator
implicit none
integer:: i,j,k,p,ip,jp,kp,im,jm,km
!*****************************************************shifting capsules upward**********************************************************
!this is a numerical trick for prevension of unrealizable values for indicator function when cells go very close to the wall
!the idea is the domain for indicator function is shifted downward to capsule become at the center of it 
do p=1,pmax
if (ap(p)%mslv>0) then
do i=1,NL
ap(p)%zfp(i)=ap(p)%zfp(i)+((ktot/4)*dz)
enddo
endif
enddo

!*****************************************************calculating G at main mesh********************************************************
call lagr2eulrind

!********************************************returning capsules to their initial value**************************************************
do p=1,pmax
if (ap(p)%mslv>0) then
do i=1,NL
ap(p)%zfp(i)=ap(p)%zfp(i)-((ktot/4)*dz)
enddo
endif
enddo
!*****************************************************calculating G at fine mesh********************************************************
if (ifoverlap==1) call lagr2eulrindch
!*****************************************calculating divergence of G for both of meshes************************************************
!$omp parallel default(shared) &
!$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
!$omp&private(dgx,dgy,dgz,dg)
!$omp do
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      ip = i + 1
      jp = j + 1
      kp = k + 1
      im = i - 1
      jm = j - 1
      km = k - 1
     dgx(i,j,k)=(gx(ip,j,k)-gx(im,j,k))*.5*dxi
     dgy(i,j,k)=(gy(i,jp,k)-gy(i,jm,k))*.5*dyi
     dgz(i,j,k)=(gz(i,j,kp)-gz(i,j,km))*.5*dzi
     dg(i,j,k)=dgx(i,j,k)+dgy(i,j,k)+dgz(i,j,k)
     indi(i,j,k)=dg(i,j,k)
    enddo
  enddo
enddo
if (ifoverlap==1) then
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
      ip = i + 1
      jp = j + 1
      kp = k + 1
      im = i - 1
      jm = j - 1
      km = k - 1
     dgxc(i,j,k)=(gxc(ip,j,k)-gxc(im,j,k))*.5*dxic
     dgyc(i,j,k)=(gyc(i,jp,k)-gyc(i,jm,k))*.5*dyic
     dgzc(i,j,k)=(gzc(i,j,kp)-gzc(i,j,km))*.5*dzic
     dgc(i,j,k)=dgxc(i,j,k)+dgyc(i,j,k)+dgzc(i,j,k)
     indich(i,j,k)=dgc(i,j,k)
    enddo
  enddo
enddo
endif
!$omp end parallel
!***************************************************computing indicator function********************************************************
call initsolver(0)
call solver1d(indi,0)
call boundind(indi)
!*******************************************returning indicator function to main mesh***************************************************
do k=1,kmax
do j=1,jmax
do i=1,imax
if (k<=kmax-(ktot/4)) then
indic(i,j,k)=indi(i,j,k+(ktot/4))
else
indic(i,j,k)=0.
endif
enddo
enddo
enddo
call boundind(indic)
do k=1,kmax
  do j=1,jmax
    do i=1,imax
    if (indic(i,j,k)<.01) indic(i,j,k)=0.
    if (indic(i,j,k)>=.9) indic(i,j,k)=1.
    enddo
  enddo
enddo
!**********************calculating source term associated with overlap boundary condition for fine mesh*********************************  
if (ifoverlap==1) then
call interpolation
  do j=1,jmaxc
    do i=1,imaxc
      indich(i,j,kmaxc)=indich(i,j,kmaxc)-(2.*indb(i,j)*dzic*dzic)
    enddo
  enddo
call solver1dc(indich,0)
call boundindch(indich)
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
    if (indich(i,j,k)<.01) indich(i,j,k)=0.
    if (indich(i,j,k)>=.9) indich(i,j,k)=1.
    enddo
  enddo
enddo
endif
!*******************************************************calculating viscosity***********************************************************
!$omp parallel default(shared) &
!$omp&private(i,j,k) &
!$omp&private(visc)
!$omp do
do k=1,kmax
  do j=1,jmax
    do i=1,imax
    visc(i,j,k)=visc0*(1.+((lambda-1.)*indic(i,j,k)))
    enddo
  enddo
enddo
if (ifoverlap==1) then
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
    viscc(i,j,k)=visc0*(1.+((lambda-1.)*indich(i,j,k)))
    enddo
  enddo
enddo
endif
!$omp end parallel
call boundp(visc)
if (ifoverlap==1) call boundviscc(viscc)
!**************************************************calculating centrifugal force********************************************************
!$omp parallel default(shared) &
!$omp&private(i,j,k) &
!$omp&private(visc)
!$omp do
do k=1,kmax
  do j=1,jmax
    do i=1,imax
     forcey(i,j,k)=forcey(i,j,k)+(.5*(indic(i,j,k)+indic(i,j+1,k))*(sin(thetainc))*((ta/(Rep**2.))**.25))
     forcez(i,j,k)=forcez(i,j,k)+(.5*(indic(i,j,k)+indic(i,j,k+1))*(-cos(thetainc))*((ta/(Rep**2.))**.25))
    enddo
  enddo
enddo
if (ifoverlap==1) then
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
     forceyc(i,j,k)=forceyc(i,j,k)+(.5*(indich(i,j,k)+indich(i,j+1,k))*(sin(thetainc))*((ta/(Rep**2.))**.25))
     forcezc(i,j,k)=forcezc(i,j,k)+(.5*(indich(i,j,k)+indich(i,j,k+1))*(-cos(thetainc))*((ta/(Rep**2.))**.25))
    enddo
  enddo
enddo
endif
!$omp end parallel
!****************************************************end of solving procedure***********************************************************
return
end subroutine indicator
end module mod_indicator
