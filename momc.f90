module mod_momc
use mod_param
use mod_common
implicit none
private
public momxadc,momxdifc,momxpc,momyadc,momydifc,momypc,momzadc,momzdifc,momzpc
contains
!
subroutine momxadc(dudtc,uoc,voc,woc,uooc,vooc,wooc)
implicit none
real, dimension(0:,0:,0:), intent(in) :: uoc,voc,woc,uooc,vooc,wooc
real, dimension(0:,0:,0:), intent(out) :: dudtc
integer :: im,ip,jm,jp,km,kp,i,j,k
real :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
real :: uuipo,uuimo,uvjpo,uvjmo,uwkpo,uwkmo
!
!$omp parallel default(shared) &
!$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
!$omp&private(uuip,uuim,uvjp,uvjm,uwkp,uwkm,dudxp,dudxm,dudyp,dudym,dudzp,dudzm,uuipo,uuimo,uvjpo,uvjmo,uwkpo,uwkmo)
!$omp do
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
      ip = i + 1
      jp = j + 1
      kp = k + 1
      im = i - 1
      jm = j - 1
      km = k - 1
      uuip  = 0.25 * ( uoc(ip,j,k)+uoc(i,j,k) )*( uoc(ip,j,k)+uoc(i,j,k)  )
      uuim  = 0.25 * ( uoc(im,j,k)+uoc(i,j,k) )*( uoc(im,j,k)+uoc(i,j,k)  )
      uvjp  = 0.25 * ( uoc(i,jp,k)+uoc(i,j,k) )*( voc(ip,j,k)+voc(i,j,k)  )
      uvjm  = 0.25 * ( uoc(i,jm,k)+uoc(i,j,k) )*( voc(ip,jm,k)+voc(i,jm,k))
      uwkp  = 0.25 * ( uoc(i,j,kp)+uoc(i,j,k) )*( woc(ip,j,k) +woc(i,j,k) )
      uwkm  = 0.25 * ( uoc(i,j,km)+uoc(i,j,k) )*( woc(ip,j,km)+woc(i,j,km))
      uuipo  = 0.25 * ( uooc(ip,j,k)+uooc(i,j,k) )*( uooc(ip,j,k)+uooc(i,j,k)  )
      uuimo  = 0.25 * ( uooc(im,j,k)+uooc(i,j,k) )*( uooc(im,j,k)+uooc(i,j,k)  )
      uvjpo  = 0.25 * ( uooc(i,jp,k)+uooc(i,j,k) )*( vooc(ip,j,k)+vooc(i,j,k)  )
      uvjmo  = 0.25 * ( uooc(i,jm,k)+uooc(i,j,k) )*( vooc(ip,jm,k)+vooc(i,jm,k))
      uwkpo  = 0.25 * ( uooc(i,j,kp)+uooc(i,j,k) )*( wooc(ip,j,k) +wooc(i,j,k) )
      uwkmo  = 0.25 * ( uooc(i,j,km)+uooc(i,j,k) )*( wooc(ip,j,km)+wooc(i,j,km))
      !
      ! Momentum balance
      !
      dudtc(i,j,k) = dxic*((1.5*( -uuip + uuim ))-(.5*( -uuipo + uuimo )))  + &
                    dyic*((1.5*( -uvjp + uvjm ))-(.5*( -uvjpo + uvjmo )))  + &
                    dzic*((1.5*( -uwkp + uwkm ))-(.5*( -uwkpo + uwkmo )))  
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momxadc
!
subroutine momxdifc(dudtc,uoc,voc,woc,cupc,cuec,cuwc,cunc,cusc,cutc,cubc)
implicit none
real, dimension(0:,0:,0:), intent(in) :: uoc,voc,woc
real, dimension(0:,0:,0:), intent(out) :: dudtc,cupc,cuec,cuwc,cunc,cusc,cutc,cubc
integer :: im,ip,jm,jp,km,kp,i,j,k
real :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm
real :: dudxop,dudxmo,dudypo,dudymo,dudzpo,dudzmo
real :: visccip,visccim,visccjp,visccjm,viscckp,viscckm
!
!$omp parallel default(shared) &
!$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
!$omp&private(dudxp,dudxm,dudyp,dudym,dudzp,dudzm,visccip,visccim,visccjp,visccjm,viscckp,viscckm,dudxop)
!$omp&private(dudxmo,dudypo,dudymo,dudzpo,dudzmo)
!$omp do
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
      ip = i + 1
      jp = j + 1
      kp = k + 1
      im = i - 1
      jm = j - 1
      km = k - 1
      visccip=viscc(ip,j,k)
      visccim=viscc(i,j,k)
      visccjp=.25*(viscc(i,j,k)+viscc(i,jp,k)+viscc(ip,j,k)+viscc(ip,jp,k))
      visccjm=.25*(viscc(i,j,k)+viscc(i,jm,k)+viscc(ip,j,k)+viscc(ip,jm,k))
      viscckp=.25*(viscc(i,j,k)+viscc(i,j,kp)+viscc(ip,j,k)+viscc(ip,j,kp))
      viscckm=.25*(viscc(i,j,k)+viscc(i,j,km)+viscc(ip,j,k)+viscc(ip,j,km))
      dudxp = (uoc(ip,j,k)-uoc(i,j,k))*dxic
      dudxm = (uoc(i,j,k)-uoc(im,j,k))*dxic
      dudyp = (uoc(i,jp,k)-uoc(i,j,k))*dyic
      dudym = (uoc(i,j,k)-uoc(i,jm,k))*dyic
      dudzp = (uoc(i,j,kp)-uoc(i,j,k))*dzic
      dudzm = (uoc(i,j,k)-uoc(i,j,km))*dzic
      !
      ! Momentum balance
      !
      dudtc(i,j,k) =  ((visccip*dudxp)-(visccim*dudxm))*dxic + &
                     ((visccjp*dudyp)-(visccjm*dudym))*dyic + &
                     ((viscckp*dudzp)-(viscckm*dudzm))*dzic 


    cupc(i,j,k)=((dxic**2.)*.5*dt*(visccip+visccim+visccjp+visccjm+viscckp+viscckm))+1.
    cuec(i,j,k)=-visccip*.5*dt*(dxic**2.)
    cuwc(i,j,k)=-visccim*.5*dt*(dxic**2.)
    cunc(i,j,k)=-visccjp*.5*dt*(dxic**2.)
    cusc(i,j,k)=-visccjm*.5*dt*(dxic**2.)
    cutc(i,j,k)=-viscckp*.5*dt*(dxic**2.)
    cubc(i,j,k)=-viscckm*.5*dt*(dxic**2.)
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momxdifc
subroutine momxpc(dudtc,p)
implicit none
real, dimension(0:,0:,0:), intent(in) :: p
real, dimension(0:,0:,0:), intent(out) :: dudtc
integer :: i,j,k
integer :: ip
!
!$omp parallel default(shared) &
!$omp&private(i,j,k,ip)
!$omp do
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
      ip = i + 1
      dudtc(i,j,k) = (-dxic*(( p(ip,j,k)-p(i,j,k) )))+forcexc(i,j,k)
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momxpc
!
subroutine momyadc(dvdtc,uoc,voc,woc,uooc,vooc,wooc)
implicit none
real, dimension(0:,0:,0:), intent(in) :: uoc,voc,woc,uooc,vooc,wooc
real, dimension(0:,0:,0:), intent(out) :: dvdtc
integer :: im,ip,jm,jp,km,kp,i,j,k
real :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
real :: uvipo,uvimo,vvjpo,vvjmo,wvkpo,wvkmo
!
!$omp parallel default(shared) &
!$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
!$omp&private(uvip,uvim,vvjp,vvjm,wvkp,wvkm,dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm,uvipo,uvimo,vvjpo,vvjmo,wvkpo,wvkmo)
!$omp do
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
      ip = i + 1
      jp = j + 1
      kp = k + 1
      im = i - 1
      jm = j - 1
      km = k - 1
      uvip  = 0.25 * (uoc(i ,j,k)+uoc(i ,jp,k))*( voc(i,j,k)+voc(ip,j,k) )
      uvim  = 0.25 * (uoc(im,j,k)+uoc(im,jp,k))*( voc(i,j,k)+voc(im,j,k) )
      vvjp  = 0.25 * (voc(i,j,k )+voc(i,jp,k) )*( voc(i,j,k)+voc(i,jp,k) )
      vvjm  = 0.25 * (voc(i,j,k )+voc(i,jm,k) )*( voc(i,j,k)+voc(i,jm,k) )
      wvkp  = 0.25 * (woc(i,j,k )+woc(i,jp,k) )*( voc(i,j,kp)+voc(i,j,k) )
      wvkm  = 0.25 * (woc(i,j,km)+woc(i,jp,km))*( voc(i,j,km)+voc(i,j,k) )
      uvipo  = 0.25 * (uooc(i ,j,k)+uooc(i ,jp,k))*( vooc(i,j,k)+vooc(ip,j,k) )
      uvimo  = 0.25 * (uooc(im,j,k)+uooc(im,jp,k))*( vooc(i,j,k)+vooc(im,j,k) )
      vvjpo  = 0.25 * (vooc(i,j,k )+vooc(i,jp,k) )*( vooc(i,j,k)+vooc(i,jp,k) )
      vvjmo  = 0.25 * (vooc(i,j,k )+vooc(i,jm,k) )*( vooc(i,j,k)+vooc(i,jm,k) )
      wvkpo  = 0.25 * (wooc(i,j,k )+wooc(i,jp,k) )*( vooc(i,j,kp)+vooc(i,j,k) )
      wvkmo  = 0.25 * (wooc(i,j,km)+wooc(i,jp,km))*( vooc(i,j,km)+vooc(i,j,k) )
      !
      ! Momentum balance
      !
      dvdtc(i,j,k) = dxic*((1.5*( -uvip + uvim ))-(.5*( -uvipo + uvimo )))  + &
                    dyic*((1.5*( -vvjp + vvjm ))-(.5*( -vvjpo + vvjmo )))  + &
                    dzic*((1.5*( -wvkp + wvkm ))-(.5*( -wvkpo + wvkmo )))  
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momyadc
!
subroutine momydifc(dvdtc,uoc,voc,woc,cvpc,cvec,cvwc,cvnc,cvsc,cvtc,cvbc)
implicit none
real, dimension(0:,0:,0:), intent(in) :: uoc,voc,woc
real, dimension(0:,0:,0:), intent(out) :: dvdtc,cvpc,cvec,cvwc,cvnc,cvsc,cvtc,cvbc
integer :: im,ip,jm,jp,km,kp,i,j,k
real :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm
real :: dvdxpo,dvdxmo,dvdypo,dvdymo,dvdzpo,dvdzmo
real :: visccip,visccim,visccjp,visccjm,viscckp,viscckm
!
!$omp parallel default(shared) &
!$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
!$omp&private(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm,visccip,visccim,visccjp,visccjm,viscckp,viscckm)
!$omp&private(dvdxpo,dvdxmo,dvdypo,dvdymo,dvdzpo,dvdzmo)
!$omp do
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
      ip = i + 1
      jp = j + 1
      kp = k + 1
      im = i - 1
      jm = j - 1
      km = k - 1
      visccip=.25*(viscc(i,j,k)+viscc(ip,j,k)+viscc(i,jp,k)+viscc(ip,jp,k))
      visccim=.25*(viscc(i,j,k)+viscc(im,j,k)+viscc(i,jp,k)+viscc(im,jp,k))
      visccjp=viscc(i,jp,k)              
      visccjm=viscc(i,j,k)                
      viscckp=.25*(viscc(i,j,k)+viscc(i,j,kp)+viscc(i,jp,k)+viscc(i,jp,kp))
      viscckm=.25*(viscc(i,j,k)+viscc(i,j,km)+viscc(i,jp,k)+viscc(i,jp,km))
      dvdxp = (voc(ip,j,k)-voc(i,j,k))*dxic
      dvdxm = (voc(i,j,k)-voc(im,j,k))*dxic
      dvdyp = (voc(i,jp,k)-voc(i,j,k))*dyic
      dvdym = (voc(i,j,k)-voc(i,jm,k))*dyic
      dvdzp = (voc(i,j,kp)-voc(i,j,k))*dzic
      dvdzm = (voc(i,j,k)-voc(i,j,km))*dzic
      !
      ! Momentum balance
      !
      dvdtc(i,j,k) =  ((visccip*dvdxp)-(visccim*dvdxm))*dxic + &
                     ((visccjp*dvdyp)-(visccjm*dvdym))*dyic + &
                     ((viscckp*dvdzp)-(viscckm*dvdzm))*dzic 


      cvpc(i,j,k)=((dxic**2.)*.5*dt*(visccip+visccim+visccjp+visccjm+viscckp+viscckm))+1.
      cvec(i,j,k)=-visccip*.5*dt*(dxic**2.)
      cvwc(i,j,k)=-visccim*.5*dt*(dxic**2.)
      cvnc(i,j,k)=-visccjp*.5*dt*(dxic**2.)
      cvsc(i,j,k)=-visccjm*.5*dt*(dxic**2.)
      cvtc(i,j,k)=-viscckp*.5*dt*(dxic**2.)
      cvbc(i,j,k)=-viscckm*.5*dt*(dxic**2.)
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momydifc
subroutine momypc(dvdtc,p)
implicit none
real, dimension(0:,0:,0:), intent(in) :: p
real, dimension(0:,0:,0:), intent(out) :: dvdtc
integer :: i,j,k
integer :: jp
!
!$omp parallel default(shared) &
!$omp&private(i,j,k,jp)
!$omp do
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
      jp = j + 1
      !
      ! Momentum balance
      !
      dvdtc(i,j,k) = (-dyic*(p(i,jp,k)-p(i,j,k)))+forceyc(i,j,k) 
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momypc
!
subroutine momzadc(dwdtc,uoc,voc,woc,uooc,vooc,wooc)
implicit none
real, dimension(0:,0:,0:), intent(in) :: uoc,voc,woc,uooc,vooc,wooc
real, dimension(0:,0:,0:), intent(out) :: dwdtc
integer :: im,ip,jm,jp,km,kp,i,j,k
real :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
real :: uwipo,uwimo,vwjpo,vwjmo,wwkpo,wwkmo
!
!$omp parallel default(shared) &
!$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
!$omp&private(uwip,uwim,vwjp,vwjm,wwkp,wwkm,dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm,uwipo,uwimo,vwjpo,vwjmo,wwkpo,wwkmo)
!$omp do
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
      ip = i + 1
      jp = j + 1
      kp = k + 1
      im = i - 1
      jm = j - 1
      km = k - 1                   
      uwip  = 0.25 * ( woc(i,j,k)+woc(ip,j,k))*(uoc(i ,j,k)+uoc(i ,j,kp) )
      uwim  = 0.25 * ( woc(i,j,k)+woc(im,j,k))*(uoc(im,j,k)+uoc(im,j,kp) )
      vwjp  = 0.25 * ( woc(i,j,k)+woc(i,jp,k))*(voc(i ,j,k)+voc(i,j ,kp) )
      vwjm  = 0.25 * ( woc(i,j,k)+woc(i,jm,k))*(voc(i,jm,k)+voc(i,jm,kp) )
      wwkp  = 0.25 * ( woc(i,j,k)+woc(i,j,kp))*(woc(i ,j,k)+woc(i,j,kp ) )
      wwkm  = 0.25 * ( woc(i,j,k)+woc(i,j,km))*(woc(i ,j,k)+woc(i,j,km ) )
      uwipo  = 0.25 * ( wooc(i,j,k)+wooc(ip,j,k))*(uooc(i ,j,k)+uooc(i ,j,kp) )
      uwimo  = 0.25 * ( wooc(i,j,k)+wooc(im,j,k))*(uooc(im,j,k)+uooc(im,j,kp) )
      vwjpo  = 0.25 * ( wooc(i,j,k)+wooc(i,jp,k))*(vooc(i ,j,k)+vooc(i,j ,kp) )
      vwjmo  = 0.25 * ( wooc(i,j,k)+wooc(i,jm,k))*(vooc(i,jm,k)+vooc(i,jm,kp) )
      wwkpo  = 0.25 * ( wooc(i,j,k)+wooc(i,j,kp))*(wooc(i ,j,k)+wooc(i,j,kp ) )
      wwkmo  = 0.25 * ( wooc(i,j,k)+wooc(i,j,km))*(wooc(i ,j,k)+wooc(i,j,km ) )
      !
      ! Momentum balance
      !
      dwdtc(i,j,k) = dxic*((1.5*( -uwip + uwim ))-(.5*( -uwipo + uwimo )))  + &
                    dyic*((1.5*( -vwjp + vwjm ))-(.5*( -vwjpo + vwjmo )))  + &
                    dzic*((1.5*( -wwkp + wwkm ))-(.5*( -wwkpo + wwkmo )))  
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momzadc
subroutine momzdifc(dwdtc,uoc,voc,woc,cwpc,cwec,cwwc,cwnc,cwsc,cwtc,cwbc)
implicit none
real, dimension(0:,0:,0:), intent(in) :: uoc,voc,woc
real, dimension(0:,0:,0:), intent(out) :: dwdtc,cwpc,cwec,cwwc,cwnc,cwsc,cwtc,cwbc
integer :: im,ip,jm,jp,km,kp,i,j,k
real :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
real :: dwdxpo,dwdxmo,dwdypo,dwdymo,dwdzpo,dwdzmo
real :: visccip,visccim,visccjp,visccjm,viscckp,viscckm
!
!$omp parallel default(shared) &
!$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
!$omp&private(uwip,uwim,vwjp,vwjm,wwkp,wwkm,dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm,visccip,visccim,visccjp,visccjm,viscckp,viscckm)
!$omp&private(dwdxpo,dwdxmo,dwdypo,dwdymo,dwdzpo,dwdzmo)
!$omp do
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
      ip = i + 1
      jp = j + 1
      kp = k + 1
      im = i - 1
      jm = j - 1
      km = k - 1
      visccip=.25*(viscc(i,j,k)+viscc(ip,j,k)+viscc(i,j,kp)+viscc(ip,j,kp))
      visccim=.25*(viscc(i,j,k)+viscc(im,j,k)+viscc(i,j,kp)+viscc(im,j,kp))
      visccjp=.25*(viscc(i,j,k)+viscc(i,jp,k)+viscc(i,j,kp)+viscc(i,jp,kp))              
      visccjm=.25*(viscc(i,j,k)+viscc(i,jm,k)+viscc(i,j,kp)+viscc(i,jm,kp))                
      viscckp=viscc(i,j,kp)                   
      viscckm=viscc(i,j,k)                    
      dwdxp = (woc(ip,j,k)-woc(i,j,k))*dxic
      dwdxm = (woc(i,j,k)-woc(im,j,k))*dxic
      dwdyp = (woc(i,jp,k)-woc(i,j,k))*dyic
      dwdym = (woc(i,j,k)-woc(i,jm,k))*dyic
      dwdzp = (woc(i,j,kp)-woc(i,j,k))*dzic
      dwdzm = (woc(i,j,k)-woc(i,j,km))*dzic
      !
      ! Momentum balance
      !
      dwdtc(i,j,k) =  ((visccip*dwdxp)-(visccim*dwdxm))*dxic + &
                     ((visccjp*dwdyp)-(visccjm*dwdym))*dyic + &
                     ((viscckp*dwdzp)-(viscckm*dwdzm))*dzic 


      cwpc(i,j,k)=((dxic**2.)*.5*dt*(visccip+visccim+visccjp+visccjm+viscckp+viscckm))+1.
      cwec(i,j,k)=-visccip*.5*dt*(dxic**2.)
      cwwc(i,j,k)=-visccim*.5*dt*(dxic**2.)
      cwnc(i,j,k)=-visccjp*.5*dt*(dxic**2.)
      cwsc(i,j,k)=-visccjm*.5*dt*(dxic**2.)
      cwtc(i,j,k)=-viscckp*.5*dt*(dxic**2.)
      cwbc(i,j,k)=-viscckm*.5*dt*(dxic**2.)      
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momzdifc
!
subroutine momzpc(dwdtc,p)
implicit none
real, dimension(0:,0:,0:), intent(in) :: p
real, dimension(0:,0:,0:), intent(out) :: dwdtc
integer :: kp
integer :: i,j,k
!
!$omp parallel default(shared) &
!$omp&private(i,j,k,kp)
!$omp do
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
      kp = k + 1
      !
      ! Momentum balance
      !
      dwdtc(i,j,k) = (-dzic*(p(i,j,kp)-p(i,j,k)))+forcezc(i,j,k)
    enddo
  enddo
enddo
!$omp end parallel
!
return
end subroutine momzpc
!
end module mod_momc
