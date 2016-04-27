module mod_init
use mod_common
use decomp_2d
use mod_param
use mod_common_mpi
implicit none
private
public init
contains
subroutine init
implicit none
integer :: i,j,k,iglob,jglob,kglob
real :: coorz,rn1,rn2,rn3
real, dimension(kmax) :: v
real umean,vmean,wmean
real umean_all,vmean_all,wmean_all
!
select case(iniu)
!
case('zer') ! velocity field = 0 forall i,j,k
  forall(i=0:i1,j=0:j1,k=0:k1)
    unew(i,j,k) = 0.
    vnew(i,j,k) = 0.
    wnew(i,j,k) = 0.
    dudt(i,j,k) = 0.
    dvdt(i,j,k) = 0.
    dwdt(i,j,k) = 0.
    pnew(i,j,k) = 0.
  end forall
!
case('cou') ! Couette flow
  do k=0,k1
    coorz = ((k-0.5)/dzi)/lz ! normalised with channel height
    do j=0,j1
      do i=0,i1
        unew(i,j,k) = 0.
        vnew(i,j,k) = 0.5*lz*shear*(1.-2.*coorz)
        wnew(i,j,k) = 0.
        dudt(i,j,k) = 0.
        dvdt(i,j,k) = 0.
        dwdt(i,j,k) = 0.
        pnew(i,j,k) = 0.
      enddo
    enddo
  enddo
!
case('log') ! Logarithmic profile
  if(isfreeslip) then
    do k=1,kmax
      v(k) = 2.5*log( (1./visc0)*(k-0.5)/(2.*kmax) ) + 5.5
      if ((k-0.5)/kmax/visc0 .le. 11.6) v(k)=(k-0.5)/kmax/visc0
    enddo
  else
    do k=1,kmax/2
      v(k) = 2.5*log( (1./visc0)*(k-0.5)/(1.*kmax) ) + 5.5
      if ((k-0.5)/kmax/visc0 .le. 11.6) v(k)=(k-0.5)/kmax/visc0
    enddo
    do k=1,kmax/2
      i = kmax+1-k
      v(kmax+1-k) = 2.5*log( (1./visc0)* (k-0.5)/(1.*kmax) ) + 5.5
      if ((k-0.5)/kmax/visc0 .le. 11.6) v(kmax+1-k)=(k-0.5)/kmax/visc0
    enddo
  endif
  !
  ! Below, random numbers are generated such that the initial field
  ! does not depend in the number of mpi tasks.
  !
  umean=0.
  vmean=0.
  wmean=0.
  !!!!!call random_seed( put = (/16112006/) )
  do kglob=1,ktot
    do jglob=1,jtot
      do iglob=1,itot
        call random_number(rn1)
        call random_number(rn2)
        call random_number(rn3)
        i = iglob-imax*coords(1)
        j = jglob-jmax*coords(2)
        k = kglob
        if( &
            (1.le.i.and.i.le.imax) .and. &
            (1.le.j.and.j.le.jmax) &
          ) then
          unew(i,j,k)=15.*(rn1-0.5)
          umean=umean+unew(i,j,k)
          vnew(i,j,k)=v(k)*(1.+5.*(rn2-0.5))
          vmean=vmean+vnew(i,j,k)
          wnew(i,j,k)=15.*(rn3-0.5)
          wmean=wmean+wnew(i,j,k)
          dudt(i,j,k)=0.
          dvdt(i,j,k)=0.
          dwdt(i,j,k)=0.
          pnew(i,j,k)=0.
        endif
      enddo
    enddo
  enddo
!  do k=1,kmax
!    do j=1,jmax
!      do i=1,imax
!        call random_number(rn)
!        unew(i,j,k)=15.*(rn-0.5)
!        umean=umean+unew(i,j,k)
!        call random_number(rn)
!        vnew(i,j,k)=v(k)*(1.+5.*(rn-0.5))
!        vmean=vmean+vnew(i,j,k)
!        call random_number(rn)
!        wnew(i,j,k)=15.*(rn-0.5)
!        wmean=wmean+wnew(i,j,k)
!        dudt(i,j,k)=0.
!        dvdt(i,j,k)=0.
!        dwdt(i,j,k)=0.
!        pnew(i,j,k)=0.
!      enddo
!    enddo
!  enddo
  call mpi_allreduce(umean,umean_all,1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(vmean,vmean_all,1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(wmean,wmean_all,1,mpi_real8,mpi_sum,comm_cart,error)
  umean_all = umean_all/(1.*itot*jtot*kmax)
  vmean_all = vmean_all/(1.*itot*jtot*kmax)
  wmean_all = wmean_all/(1.*itot*jtot*kmax)
  do k=1,kmax
    do j=1,jmax
      do i=1,imax
        unew(i,j,k) = unew(i,j,k)-umean_all !average = 0
        vnew(i,j,k) = vnew(i,j,k)/vmean_all !average = 1
        wnew(i,j,k) = wnew(i,j,k)-wmean_all !average = 0
      enddo
    enddo
  enddo
!
  case('poi') ! Poiseuille
  do k = 0,k1
    coorz = (k-0.5)*dz/lz ! normalised with channel height
    do j = 0,j1
      do i = 0,i1
        unew(i,j,k) = 0.
        vnew(i,j,k) = 6.*bulk_v_sup*coorz*(1.-coorz)
        wnew(i,j,k) = 0.
        dudt(i,j,k) = 0.
        dvdt(i,j,k) = 0.
        dwdt(i,j,k) = 0.
        pnew(i,j,k) = 0.
      enddo
    enddo
  enddo
end select
!
return
end subroutine init
!
end module mod_init
