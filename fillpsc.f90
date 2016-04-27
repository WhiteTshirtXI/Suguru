module mod_fillpsc
use mod_param
use mod_common
implicit none
private
public fillpsc
contains
subroutine fillpsc(pc)
real, intent(out), dimension(0:,0:,0:) :: pc
real :: dti
real :: dtidxic,dtidyic,dtidzic
integer :: i,j,k,im,jm,km
!
dti = 1./dt
dtidxic = dti*dxic
dtidyic = dti*dyic
dtidzic = dti*dzic
!
!
!  fill the right-hand side of the Poisson equation for the correction pressure.
!
!  the discrete divergence is:
!
!  w(i,j,k)-w(i,j,k-1)   v(i,j,k)-v(i,j-1,k)   u(i,j,k)-u(i-1,j,k)
!  ------------------- + ------------------- + -------------------  = div
!          dz                    dy                    dx
!
!  note: in this subroutine p is not the correction pressure, but
!  the rhs of the Poisson equation.
!
do k=1,kmaxc
  km = k-1
  do j=1,jmaxc
    jm = j-1
    do i=1,imaxc
      im = i-1
      pc(i,j,k) = ( &
                ( wnewc(i,j,k)-wnewc(i,j,km))*dtidzic+ &
                ( vnewc(i,j,k)-vnewc(i,jm,k))*dtidyic+ &
                ( unewc(i,j,k)-unewc(im,j,k))*dtidxic )
    enddo
  enddo
enddo
!
return
end subroutine fillpsc
!
end module mod_fillpsc
