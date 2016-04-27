module mod_initsolver
use mod_common
use mod_common_mpi
use mod_param
use decomp_2d
implicit none
private
public initsolver
contains
!
subroutine initsolver(e)
integer,intent(in)::e
integer :: i,j,k,iv,jv
real :: xrt(itot), yrt(jtot)
real :: xrtc(itotc), yrtc(jtotc)
real :: bbc
!
! Generate tridiagonal matrix
!

do k=1,kmax
  a(k) = dzi*dzi
  c(k) = dzi*dzi
  b(k) = -(a(k) + c(k))
enddo
do k=1,kmaxc
  ac(k) = dzic*dzic
  cc(k) = dzic*dzic
  bc(k) = -(ac(k) + cc(k))
enddo

!
! boundary condition  in z-direction;
! consistent with no-penetration condition for prediction velocity 
! at solid walls (Channel)
!

if (e==0 .or.e==1) then


b(1) = b(1) + a(1)
b(kmax) = b(kmax) + c(kmax)
a(1) = 0.
 c(kmax) = 0.


if(e==0) then
bc(1) = bc(1) + ac(1)
bc(kmaxc) = bc(kmaxc) - cc(kmaxc)
ac(1) = 0.
 cc(kmaxc) = 0.
else if (e==1) then
bc(1) = bc(1) + ac(1)
bc(kmaxc) = bc(kmaxc) + cc(kmaxc)
ac(1) = 0.
 cc(kmaxc) = 0.
endif

else if (e==2 .or.e==3) then

if(e==2) then
b(1) = b(1) - a(1)
b(kmax) = b(kmax) - c(kmax)
a(1) = 0.
 c(kmax) = 0.
else if (e==3) then
b(kmax) = 1.
a(1) = 0.
 a(kmax)=0.
 c(kmax) = 0.
endif



if(e==2) then
bc(1) = bc(1) - ac(1)
bc(kmaxc) = bc(kmaxc) - cc(kmaxc)
ac(1) = 0.
 cc(kmaxc) = 0.
else if (e==3) then
bc(kmaxc) = 1.
ac(1) = 0.
 ac(kmaxc)=0.
 cc(kmaxc) = 0.
endif

do k=1,kmax
b(k)=b(k)-(Rep/dt)
enddo

do k=1,kmaxc
bc(k)=bc(k)-(Rep/dt)
enddo
if (e==3) then
b(kmax)=1.
bc(kmaxc)=1.
endif

endif

!
! set lookup tables.
!
call vrffti(itot,wi)
call vrffti(jtot,wj)
call vrffti(itotc,wic)
call vrffti(jtotc,wjc)
!
! generate eigenvalues ( xrt and yrt ).
!
!
!************************************************main mesh*******************************************
! x direction
!
do i=3,itot,2
  xrt(i-1) = -4.*dxi*dxi*(sin(float((i-1))*pi/(2.*itot)))**2
  xrt(i) = xrt(i-1)
enddo
xrt(1   ) = 0.
xrt(itot) = -4.*dxi*dxi
!
! y direction
!
do j=3,jtot,2
  yrt(j-1) = -4.*dyi*dyi*(sin(float((j-1))*pi/(2.*jtot)))**2
  yrt(j  ) = yrt(j-1)
enddo
yrt(1   ) = 0.
yrt(jtot) = -4.*dyi*dyi
!
do j=1,jmax
  jv = j + zstart(2) - 1
  do i=1,imax
    iv = i + zstart(1) - 1
    xyrt(i,j) = xrt(iv)+yrt(jv)
  enddo
enddo
!************************************************fine mesh*******************************************
! x direction
!
do i=3,itotc,2
  xrtc(i-1) = -4.*dxic*dxic*(sin(float((i-1))*pi/(2.*itotc)))**2
  xrtc(i) = xrtc(i-1)
enddo
xrtc(1   ) = 0.
xrtc(itotc) = -4.*dxic*dxic
!
! y direction
!
do j=3,jtotc,2
  yrtc(j-1) = -4.*dyic*dyic*(sin(float((j-1))*pi/(2.*jtotc)))**2
  yrtc(j  ) = yrtc(j-1)
enddo
yrtc(1   ) = 0.
yrtc(jtotc) = -4.*dyic*dyic
!
do j=1,jmaxc
  jv = j + zstartc(2) - 1
  do i=1,imaxc
    iv = i + zstartc(1) - 1
    xyrtc(i,j) = xrtc(iv)+yrtc(jv)
  enddo
enddo
!
return
end subroutine initsolver
!
end module mod_initsolver
