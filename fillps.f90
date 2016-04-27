module mod_fillps
use mod_param
use mod_common
implicit none
private
public fillps
contains
subroutine fillps(p)
real, intent(out), dimension(0:,0:,0:) :: p
real :: dti
real :: dtidxi,dtidyi,dtidzi
integer :: i,j,k,im,jm,km
!
dti = 1./dt
dtidxi = dti*dxi
dtidyi = dti*dyi
dtidzi = dti*dzi
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
do k=1,kmax
  km = k-1
  do j=1,jmax
    jm = j-1
    do i=1,imax
      im = i-1
      p(i,j,k) = ( &
                ( wnew(i,j,k)-wnew(i,j,km))*dtidzi+ &
                ( vnew(i,j,k)-vnew(i,jm,k))*dtidyi+ &
                ( unew(i,j,k)-unew(im,j,k))*dtidxi )
    enddo
  enddo
enddo

!
return
end subroutine fillps
!
end module mod_fillps
