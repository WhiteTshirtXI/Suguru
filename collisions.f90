module mod_collisions !VER O QUE ACONTECE SE TIVER NDT=40 E ITERMAX = 1
use mpi
use mod_common
use mod_common_mpi
implicit none
private
public collisions
contains
!
subroutine collisions
implicit none
integer:: i,j,k,p,jmaxcol
real:: ss,cutoff,fc(NL)
cutoff=rad0
do p=1,pmax
if (ap(p)%mslv>0) then
do i=1,NL
fc(i)=0.
jmaxcol=nint(sqrt(max(((cutoff**2.)-(ap(p)%zfp(i)**2.)),0.))*dyi)
do j=-jmaxcol,jmaxcol
do k=-jmaxcol,jmaxcol
ss=sqrt((ap(p)%zfp(i)**2.)+((j*dy)**2.)+((k*dx)**2.))
if (ss<rad0) then
fc(i)=fc(i)+(2.*(col/(Rep**2.)))*(exp(2.*betacol*(rad0-ss))-exp(betacol*(rad0-ss)))*abs(ap(p)%zfp(i)/ss)*(min(ap(p)%area(i),(dx**2.))/max(ap(p)%area(i),(dx**2.)))
endif
enddo
enddo
ap(p)%fzl(i)=ap(p)%fzl(i)+fc(i)
enddo
endif
enddo
end subroutine collisions

end module mod_collisions
