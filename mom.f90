module mod_mom
use mod_param
use mod_common
implicit none
private
public momxad,momxdif,momxp,momyad,momydif,momyp,momzad,momzdif,momzp
contains
!
subroutine momxad(dudt,uo,vo,wo,uoo,voo,woo)
implicit none
real, dimension(0:,0:,0:), intent(in) :: uo,vo,wo,uoo,voo,woo
real, dimension(0:,0:,0:), intent(out) :: dudt
integer :: im,ip,jm,jp,km,kp,i,j,k
real :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
real :: uuipo,uuimo,uvjpo,uvjmo,uwkpo,uwkmo
!
!$omp parallel default(shared) &
!$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
!$omp&private(uuip,uuim,uvjp,uvjm,uwkp,uwkm,dudxp,dudxm,dudyp,dudym,dudzp,dudzm,uuipo,uuimo,uvjpo,uvjmo,uwkpo,uwkmo)
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
      uuip  = 0.25 * ( Uo(ip,j,k)+Uo(i,j,k) )*( Uo(ip,j,k)+Uo(i,j,k)  )
      uuim  = 0.25 * ( Uo(im,j,k)+Uo(i,j,k) )*( Uo(im,j,k)+Uo(i,j,k)  )
      uvjp  = 0.25 * ( Uo(i,jp,k)+Uo(i,j,k) )*( Vo(ip,j,k)+Vo(i,j,k)  )
      uvjm  = 0.25 * ( Uo(i,jm,k)+Uo(i,j,k) )*( Vo(ip,jm,k)+Vo(i,jm,k))
      uwkp  = 0.25 * ( Uo(i,j,kp)+Uo(i,j,k) )*( Wo(ip,j,k) +Wo(i,j,k) )
      uwkm  = 0.25 * ( Uo(i,j,km)+Uo(i,j,k) )*( Wo(ip,j,km)+Wo(i,j,km))
      uuipo  = 0.25 * ( Uoo(ip,j,k)+Uoo(i,j,k) )*( Uoo(ip,j,k)+Uoo(i,j,k)  )
      uuimo  = 0.25 * ( Uoo(im,j,k)+Uoo(i,j,k) )*( Uoo(im,j,k)+Uoo(i,j,k)  )
      uvjpo  = 0.25 * ( Uoo(i,jp,k)+Uoo(i,j,k) )*( Voo(ip,j,k)+Voo(i,j,k)  )
      uvjmo  = 0.25 * ( Uoo(i,jm,k)+Uoo(i,j,k) )*( Voo(ip,jm,k)+Voo(i,jm,k))
      uwkpo  = 0.25 * ( Uoo(i,j,kp)+Uoo(i,j,k) )*( Woo(ip,j,k) +Woo(i,j,k) )
      uwkmo  = 0.25 * ( Uoo(i,j,km)+Uoo(i,j,k) )*( Woo(ip,j,km)+Woo(i,j,km))
      !
      ! Momentum balance
      !
      dudt(i,j,k) = dxi*((1.5*( -uuip + uuim ))-(.5*( -uuipo + uuimo )))  + &
                    dyi*((1.5*( -uvjp + uvjm ))-(.5*( -uvjpo + uvjmo )))  + &
                    dzi*((1.5*( -uwkp + uwkm ))-(.5*( -uwkpo + uwkmo )))  
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momxad
!
subroutine momxdif(dudt,uo,vo,wo,cup,cue,cuw,cun,cus,cut,cub)
implicit none
real, dimension(0:,0:,0:), intent(in) :: uo,vo,wo
real, dimension(0:,0:,0:), intent(out) :: dudt,cup,cue,cuw,cun,cus,cut,cub
integer :: im,ip,jm,jp,km,kp,i,j,k
real :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm
real :: dudxop,dudxmo,dudypo,dudymo,dudzpo,dudzmo
real :: viscip,viscim,viscjp,viscjm,visckp,visckm
!
!$omp parallel default(shared) &
!$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
!$omp&private(dudxp,dudxm,dudyp,dudym,dudzp,dudzm,viscip,viscim,viscjp,viscjm,visckp,visckm,dudxop)
!$omp&private(dudxmo,dudypo,dudymo,dudzpo,dudzmo)
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
      viscip=visc(ip,j,k)
      viscim=visc(i,j,k)
      viscjp=.25*(visc(i,j,k)+visc(i,jp,k)+visc(ip,j,k)+visc(ip,jp,k))
      viscjm=.25*(visc(i,j,k)+visc(i,jm,k)+visc(ip,j,k)+visc(ip,jm,k))
      visckp=.25*(visc(i,j,k)+visc(i,j,kp)+visc(ip,j,k)+visc(ip,j,kp))
      visckm=.25*(visc(i,j,k)+visc(i,j,km)+visc(ip,j,k)+visc(ip,j,km))
      dudxp = (Uo(ip,j,k)-Uo(i,j,k))*dxi
      dudxm = (Uo(i,j,k)-Uo(im,j,k))*dxi
      dudyp = (Uo(i,jp,k)-Uo(i,j,k))*dyi
      dudym = (Uo(i,j,k)-Uo(i,jm,k))*dyi
      dudzp = (Uo(i,j,kp)-Uo(i,j,k))*dzi
      dudzm = (Uo(i,j,k)-Uo(i,j,km))*dzi
      !
      ! Momentum balance
      !
      dudt(i,j,k) =  ((viscip*dudxp)-(viscim*dudxm))*dxi + &
                     ((viscjp*dudyp)-(viscjm*dudym))*dyi + &
                     ((visckp*dudzp)-(visckm*dudzm))*dzi 


    cup(i,j,k)=((dxi**2.)*.5*dt*(viscip+viscim+viscjp+viscjm+visckp+visckm))+1.
    cue(i,j,k)=-viscip*.5*dt*(dxi**2.)
    cuw(i,j,k)=-viscim*.5*dt*(dxi**2.)
    cun(i,j,k)=-viscjp*.5*dt*(dxi**2.)
    cus(i,j,k)=-viscjm*.5*dt*(dxi**2.)
    cut(i,j,k)=-visckp*.5*dt*(dxi**2.)
    cub(i,j,k)=-visckm*.5*dt*(dxi**2.)
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momxdif
subroutine momxp(dudt,p)
implicit none
real, dimension(0:,0:,0:), intent(in) :: p
real, dimension(0:,0:,0:), intent(out) :: dudt
integer :: i,j,k
integer :: ip
!
!$omp parallel default(shared) &
!$omp&private(i,j,k,ip)
!$omp do
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      ip = i + 1
      dudt(i,j,k) = (-dxi*(( p(ip,j,k)-p(i,j,k) )))+forcex(i,j,k)
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momxp
!
subroutine momyad(dvdt,uo,vo,wo,uoo,voo,woo)
implicit none
real, dimension(0:,0:,0:), intent(in) :: uo,vo,wo,uoo,voo,woo
real, dimension(0:,0:,0:), intent(out) :: dvdt
integer :: im,ip,jm,jp,km,kp,i,j,k
real :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
real :: uvipo,uvimo,vvjpo,vvjmo,wvkpo,wvkmo
!
!$omp parallel default(shared) &
!$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
!$omp&private(uvip,uvim,vvjp,vvjm,wvkp,wvkm,dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm,uvipo,uvimo,vvjpo,vvjmo,wvkpo,wvkmo)
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
      uvip  = 0.25 * (Uo(i ,j,k)+Uo(i ,jp,k))*( Vo(i,j,k)+Vo(ip,j,k) )
      uvim  = 0.25 * (Uo(im,j,k)+Uo(im,jp,k))*( Vo(i,j,k)+Vo(im,j,k) )
      vvjp  = 0.25 * (Vo(i,j,k )+Vo(i,jp,k) )*( Vo(i,j,k)+Vo(i,jp,k) )
      vvjm  = 0.25 * (Vo(i,j,k )+Vo(i,jm,k) )*( Vo(i,j,k)+Vo(i,jm,k) )
      wvkp  = 0.25 * (Wo(i,j,k )+Wo(i,jp,k) )*( Vo(i,j,kp)+Vo(i,j,k) )
      wvkm  = 0.25 * (Wo(i,j,km)+Wo(i,jp,km))*( Vo(i,j,km)+Vo(i,j,k) )
      uvipo  = 0.25 * (Uoo(i ,j,k)+Uoo(i ,jp,k))*( Voo(i,j,k)+Voo(ip,j,k) )
      uvimo  = 0.25 * (Uoo(im,j,k)+Uoo(im,jp,k))*( Voo(i,j,k)+Voo(im,j,k) )
      vvjpo  = 0.25 * (Voo(i,j,k )+Voo(i,jp,k) )*( Voo(i,j,k)+Voo(i,jp,k) )
      vvjmo  = 0.25 * (Voo(i,j,k )+Voo(i,jm,k) )*( Voo(i,j,k)+Voo(i,jm,k) )
      wvkpo  = 0.25 * (Woo(i,j,k )+Woo(i,jp,k) )*( Voo(i,j,kp)+Voo(i,j,k) )
      wvkmo  = 0.25 * (Woo(i,j,km)+Woo(i,jp,km))*( Voo(i,j,km)+Voo(i,j,k) )
      !
      ! Momentum balance
      !
      dvdt(i,j,k) = dxi*((1.5*( -uvip + uvim ))-(.5*( -uvipo + uvimo )))  + &
                    dyi*((1.5*( -vvjp + vvjm ))-(.5*( -vvjpo + vvjmo )))  + &
                    dzi*((1.5*( -wvkp + wvkm ))-(.5*( -wvkpo + wvkmo )))  
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momyad
!
subroutine momydif(dvdt,uo,vo,wo,cvp,cve,cvw,cvn,cvs,cvt,cvb)
implicit none
real, dimension(0:,0:,0:), intent(in) :: uo,vo,wo
real, dimension(0:,0:,0:), intent(out) :: dvdt,cvp,cve,cvw,cvn,cvs,cvt,cvb
integer :: im,ip,jm,jp,km,kp,i,j,k
real :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm
real :: dvdxpo,dvdxmo,dvdypo,dvdymo,dvdzpo,dvdzmo
real :: viscip,viscim,viscjp,viscjm,visckp,visckm
!
!$omp parallel default(shared) &
!$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
!$omp&private(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm,viscip,viscim,viscjp,viscjm,visckp,visckm)
!$omp&private(dvdxpo,dvdxmo,dvdypo,dvdymo,dvdzpo,dvdzmo)
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
      viscip=.25*(visc(i,j,k)+visc(ip,j,k)+visc(i,jp,k)+visc(ip,jp,k))
      viscim=.25*(visc(i,j,k)+visc(im,j,k)+visc(i,jp,k)+visc(im,jp,k))
      viscjp=visc(i,jp,k)              
      viscjm=visc(i,j,k)                
      visckp=.25*(visc(i,j,k)+visc(i,j,kp)+visc(i,jp,k)+visc(i,jp,kp))
      visckm=.25*(visc(i,j,k)+visc(i,j,km)+visc(i,jp,k)+visc(i,jp,km))
      dvdxp = (Vo(ip,j,k)-Vo(i,j,k))*dxi
      dvdxm = (Vo(i,j,k)-Vo(im,j,k))*dxi
      dvdyp = (Vo(i,jp,k)-Vo(i,j,k))*dyi
      dvdym = (Vo(i,j,k)-Vo(i,jm,k))*dyi
      dvdzp = (Vo(i,j,kp)-Vo(i,j,k))*dzi
      dvdzm = (Vo(i,j,k)-Vo(i,j,km))*dzi
      !
      ! Momentum balance
      !
      dvdt(i,j,k) =  ((viscip*dvdxp)-(viscim*dvdxm))*dxi + &
                     ((viscjp*dvdyp)-(viscjm*dvdym))*dyi + &
                     ((visckp*dvdzp)-(visckm*dvdzm))*dzi 


      cvp(i,j,k)=((dxi**2.)*.5*dt*(viscip+viscim+viscjp+viscjm+visckp+visckm))+1.
      cve(i,j,k)=-viscip*.5*dt*(dxi**2.)
      cvw(i,j,k)=-viscim*.5*dt*(dxi**2.)
      cvn(i,j,k)=-viscjp*.5*dt*(dxi**2.)
      cvs(i,j,k)=-viscjm*.5*dt*(dxi**2.)
      cvt(i,j,k)=-visckp*.5*dt*(dxi**2.)
      cvb(i,j,k)=-visckm*.5*dt*(dxi**2.)
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momydif
subroutine momyp(dvdt,p)
implicit none
real, dimension(0:,0:,0:), intent(in) :: p
real, dimension(0:,0:,0:), intent(out) :: dvdt
integer :: i,j,k
integer :: jp
!
!$omp parallel default(shared) &
!$omp&private(i,j,k,jp)
!$omp do
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      jp = j + 1
      !
      ! Momentum balance
      !
      dvdt(i,j,k) = (-dyi*(p(i,jp,k)-p(i,j,k)))+forcey(i,j,k) 
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momyp
!
subroutine momzad(dwdt,uo,vo,wo,uoo,voo,woo)
implicit none
real, dimension(0:,0:,0:), intent(in) :: uo,vo,wo,uoo,voo,woo
real, dimension(0:,0:,0:), intent(out) :: dwdt
integer :: im,ip,jm,jp,km,kp,i,j,k
real :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
real :: uwipo,uwimo,vwjpo,vwjmo,wwkpo,wwkmo
!
!$omp parallel default(shared) &
!$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
!$omp&private(uwip,uwim,vwjp,vwjm,wwkp,wwkm,dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm,uwipo,uwimo,vwjpo,vwjmo,wwkpo,wwkmo)
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
      uwip  = 0.25 * ( Wo(i,j,k)+Wo(ip,j,k))*(Uo(i ,j,k)+Uo(i ,j,kp) )
      uwim  = 0.25 * ( Wo(i,j,k)+Wo(im,j,k))*(Uo(im,j,k)+Uo(im,j,kp) )
      vwjp  = 0.25 * ( Wo(i,j,k)+Wo(i,jp,k))*(Vo(i ,j,k)+Vo(i,j ,kp) )
      vwjm  = 0.25 * ( Wo(i,j,k)+Wo(i,jm,k))*(Vo(i,jm,k)+Vo(i,jm,kp) )
      wwkp  = 0.25 * ( Wo(i,j,k)+Wo(i,j,kp))*(Wo(i ,j,k)+Wo(i,j,kp ) )
      wwkm  = 0.25 * ( Wo(i,j,k)+Wo(i,j,km))*(Wo(i ,j,k)+Wo(i,j,km ) )
      uwipo  = 0.25 * ( Woo(i,j,k)+Woo(ip,j,k))*(Uoo(i ,j,k)+Uoo(i ,j,kp) )
      uwimo  = 0.25 * ( Woo(i,j,k)+Woo(im,j,k))*(Uoo(im,j,k)+Uoo(im,j,kp) )
      vwjpo  = 0.25 * ( Woo(i,j,k)+Woo(i,jp,k))*(Voo(i ,j,k)+Voo(i,j ,kp) )
      vwjmo  = 0.25 * ( Woo(i,j,k)+Woo(i,jm,k))*(Voo(i,jm,k)+Voo(i,jm,kp) )
      wwkpo  = 0.25 * ( Woo(i,j,k)+Woo(i,j,kp))*(Woo(i ,j,k)+Woo(i,j,kp ) )
      wwkmo  = 0.25 * ( Woo(i,j,k)+Woo(i,j,km))*(Woo(i ,j,k)+Woo(i,j,km ) )
      !
      ! Momentum balance
      !
      dwdt(i,j,k) = dxi*((1.5*( -uwip + uwim ))-(.5*( -uwipo + uwimo )))  + &
                    dyi*((1.5*( -vwjp + vwjm ))-(.5*( -vwjpo + vwjmo )))  + &
                    dzi*((1.5*( -wwkp + wwkm ))-(.5*( -wwkpo + wwkmo )))  
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momzad
subroutine momzdif(dwdt,uo,vo,wo,cwp,cwe,cww,cwn,cws,cwt,cwb)
implicit none
real, dimension(0:,0:,0:), intent(in) :: uo,vo,wo
real, dimension(0:,0:,0:), intent(out) :: dwdt,cwp,cwe,cww,cwn,cws,cwt,cwb
integer :: im,ip,jm,jp,km,kp,i,j,k
real :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
real :: dwdxpo,dwdxmo,dwdypo,dwdymo,dwdzpo,dwdzmo
real :: viscip,viscim,viscjp,viscjm,visckp,visckm
!
!$omp parallel default(shared) &
!$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
!$omp&private(uwip,uwim,vwjp,vwjm,wwkp,wwkm,dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm,viscip,viscim,viscjp,viscjm,visckp,visckm)
!$omp&private(dwdxpo,dwdxmo,dwdypo,dwdymo,dwdzpo,dwdzmo)
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
      viscip=.25*(visc(i,j,k)+visc(ip,j,k)+visc(i,j,kp)+visc(ip,j,kp))
      viscim=.25*(visc(i,j,k)+visc(im,j,k)+visc(i,j,kp)+visc(im,j,kp))
      viscjp=.25*(visc(i,j,k)+visc(i,jp,k)+visc(i,j,kp)+visc(i,jp,kp))              
      viscjm=.25*(visc(i,j,k)+visc(i,jm,k)+visc(i,j,kp)+visc(i,jm,kp))                
      visckp=visc(i,j,kp)                   
      visckm=visc(i,j,k)                    
      dwdxp = (Wo(ip,j,k)-Wo(i,j,k))*dxi
      dwdxm = (Wo(i,j,k)-Wo(im,j,k))*dxi
      dwdyp = (Wo(i,jp,k)-Wo(i,j,k))*dyi
      dwdym = (Wo(i,j,k)-Wo(i,jm,k))*dyi
      dwdzp = (Wo(i,j,kp)-Wo(i,j,k))*dzi
      dwdzm = (Wo(i,j,k)-Wo(i,j,km))*dzi
      !
      ! Momentum balance
      !
      dwdt(i,j,k) =  ((viscip*dwdxp)-(viscim*dwdxm))*dxi + &
                     ((viscjp*dwdyp)-(viscjm*dwdym))*dyi + &
                     ((visckp*dwdzp)-(visckm*dwdzm))*dzi 


      cwp(i,j,k)=((dxi**2.)*.5*dt*(viscip+viscim+viscjp+viscjm+visckp+visckm))+1.
      cwe(i,j,k)=-viscip*.5*dt*(dxi**2.)
      cww(i,j,k)=-viscim*.5*dt*(dxi**2.)
      cwn(i,j,k)=-viscjp*.5*dt*(dxi**2.)
      cws(i,j,k)=-viscjm*.5*dt*(dxi**2.)
      cwt(i,j,k)=-visckp*.5*dt*(dxi**2.)
      cwb(i,j,k)=-visckm*.5*dt*(dxi**2.)      
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momzdif
!
subroutine momzp(dwdt,p)
implicit none
real, dimension(0:,0:,0:), intent(in) :: p
real, dimension(0:,0:,0:), intent(out) :: dwdt
integer :: kp
integer :: i,j,k
!
!$omp parallel default(shared) &
!$omp&private(i,j,k,kp)
!$omp do
do k=1,kmax
  do j=1,jmax
    do i=1,imax
      kp = k + 1
      !
      ! Momentum balance
      !
      dwdt(i,j,k) = (-dzi*(p(i,j,kp)-p(i,j,k)))+forcez(i,j,k)
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momzp
!
end module mod_mom
