module mod_forcing
use mod_common
use mod_common_mpi
use mod_surf_vol
use mod_sph
implicit none
private
public complagrforces,updtintermediatevel
contains
!
subroutine complagrforces
implicit none
integer :: l,p,i,j,cap
real :: dti
real :: fxe(n_cap,m_cap),fye(n_cap,m_cap),fze(n_cap,m_cap)
real :: fxef(n_capf,m_capf),fyef(n_capf,m_capf),fzef(n_capf,m_capf)
real::  x_n_cap(n_cap,m_cap),y_n_cap(n_cap,m_cap),z_n_cap(n_cap,m_cap)
real::  x_n_capf(n_capf,m_capf),y_n_capf(n_capf,m_capf),z_n_capf(n_capf,m_capf)
!
dti = 1./dt
!
!$omp parallel default(shared) &
!$omp& private(p,l)
!$omp do
intl=0
igrida(1)=-2
igrida(2)=1
igridb(1)=-2
igridb(2)=1

  do p=1,pmax

  if (ap(p)%mslv .gt. 0) then



    do i=1,n_capf
    do j=1,m_capf
    cap=((i-1)*m_capf)+j
    x_cap(i,j)= ap(p)%xfp(cap)
    y_cap(i,j)= ap(p)%yfp(cap)
    z_cap(i,j)= ap(p)%zfp(cap)
    enddo
    enddo

!*******************************interpolating position of Lag points to the grid used for spherical harmonics*************************************

call TRSSPH (intl,igrida,m_capf,n_capf,x_cap,igridb,m_cap,n_cap, &
            x_cap_sph,wsave,lssave,lsvmin,work,lswork,lwkmin,dwork,lsdwork,ier)
intl=1

call TRSSPH (intl,igrida,m_capf,n_capf,y_cap,igridb,m_cap,n_cap, &
            y_cap_sph,wsave,lssave,lsvmin,work,lswork,lwkmin,dwork,lsdwork,ier)
call TRSSPH (intl,igrida,m_capf,n_capf,z_cap,igridb,m_cap,n_cap, &
            z_cap_sph,wsave,lssave,lsvmin,work,lswork,lwkmin,dwork,lsdwork,ier)
!***********************************************computing elastic force***************************************************************************
  call cal_spectral_force(x_cap_sph,y_cap_sph,z_cap_sph,fxe,fye,fze)

do i=1,n_cap
 do j=1,m_cap 
  x_n_cap(i,j)=unitO(1,i,j)
  y_n_cap(i,j)=unitO(2,i,j)
  z_n_cap(i,j)=unitO(3,i,j)
 enddo
enddo
!*******************************interpolating some variables to the grid used for spreading to Eulerian mesh*************************************
!**********************************************force*****************************************
intl=0
call TRSSPH (intl,igrida,m_cap,n_cap,fxe,igridb,m_capf,n_capf, &
            fxef,wsave,lssave,lsvmin,work,lswork,lwkmin,dwork,lsdwork,ier)
intl=1
call TRSSPH (intl,igrida,m_cap,n_cap,fye,igridb,m_capf,n_capf, &
            fyef,wsave,lssave,lsvmin,work,lswork,lwkmin,dwork,lsdwork,ier)
call TRSSPH (intl,igrida,m_cap,n_cap,fze,igridb,m_capf,n_capf, &
            fzef,wsave,lssave,lsvmin,work,lswork,lwkmin,dwork,lsdwork,ier)
!************************************surface normals*****************************************
intl=0
call TRSSPH (intl,igrida,m_cap,n_cap,x_n_cap,igridb,m_capf,n_capf, &
            x_n_capf,wsave,lssave,lsvmin,work,lswork,lwkmin,dwork,lsdwork,ier)
intl=1
call TRSSPH (intl,igrida,m_cap,n_cap,y_n_cap,igridb,m_capf,n_capf, &
            y_n_capf,wsave,lssave,lsvmin,work,lswork,lwkmin,dwork,lsdwork,ier)
call TRSSPH (intl,igrida,m_cap,n_cap,z_n_cap,igridb,m_capf,n_capf, &
            z_n_capf,wsave,lssave,lsvmin,work,lswork,lwkmin,dwork,lsdwork,ier)
!*******************************surface area metrics*****************************************
intl=0
call TRSSPH (intl,igrida,m_cap,n_cap,surfmetD,igridb,m_capf,n_capf, &
            surfmetDf,wsave,lssave,lsvmin,work,lswork,lwkmin,dwork,lsdwork,ier)


  area_cap(:,:)=surfmetDf(:,:)*wgtdsinf(:,:)

    do i=1,n_capf
    do j=1,m_capf
    cap=((i-1)*m_capf)+j
    ap(p)%area(cap)=area_cap(i,j)
    ap(p)%fxl(cap)=fxef(i,j)
    ap(p)%fyl(cap)=fyef(i,j)
    ap(p)%fzl(cap)=fzef(i,j)
    ap(p)%xn(cap)=x_n_capf(i,j)
    ap(p)%yn(cap)=y_n_capf(i,j)
    ap(p)%zn(cap)=z_n_capf(i,j)
    enddo
    enddo


  endif
enddo


!$omp end parallel
!
return
end subroutine complagrforces
!
!
subroutine updtintermediatevel
implicit none
integer i,j,k,p
real sumfx,sumfy,sumfz
real sumfx_all,sumfy_all,sumfz_all,bulk_v_periodic
real forcextot_all,forceytot_all,forceztot_all



!forceytot = sumfy + wallshearnew
!
! flux imposed
!
forcextot = 0.
forceytot = 0.
forceytot = 0.

bulk_v_periodic = 2.*sin(time)/3.

if(iniu.eq.'poi'.or.iniu.eq.'log') then

!me  forceytot = -1.0*(bulk_v_sup-v_bulk)/dt
!  forceytot = -0.05!-0.25*cos(time)/3.!-1.0*(bulk_v_sup-v_bulk)/dt
endif

!$omp parallel default(shared) &
!$omp& private(i,j,k) 
!$omp do
      
do k=0,k1
  do j=0,j1
    do i=0,i1


      dudt(i,j,k) = dudtold(i,j,k) + forcex(i,j,k)*dt
      !forceytot subtracted to get net zero forcing
       !-dp/dy needed to balance total drag force = -forceytot
      dvdt(i,j,k) = dvdtold(i,j,k) + forcey(i,j,k)*dt
      dwdt(i,j,k) = dwdtold(i,j,k) + forcez(i,j,k)*dt

    enddo
  enddo
enddo
!$omp end parallel
!
end subroutine updtintermediatevel
!
end module mod_forcing
