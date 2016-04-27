module mod_common
use mod_param
!
COMMON /MATLAW/HF_GsCELL
!********************************************main mesh variables**********************************************************************
real ,dimension(0:i1,0:j1,0:k1) :: unew,vnew,wnew,pnew,uo,vo,wo,uoo,voo,woo,upi,vpi,wpi,dudt,dvdt,dwdt, &
                                   dudtold,dvdtold,dwdtold,ufft,vfft,wfft, &
                                   forcex,forcey,forcez, &
                                   gx,gy,gz,dgx,dgy,dgz,dg,indi,indic,visc
!*************************************************overlapped mesh variables***********************************************************
real ,dimension(0:i1c,0:j1c,0:k1c+1) :: unewc,vnewc,wnewc,uoc,voc,woc,uooc,vooc,wooc,upic,vpic,wpic
                                   
real ,dimension(0:i1c,0:j1c,0:k1c) :: pnewc,forcexc,forceyc,forcezc,dudtc,dvdtc,dwdtc, &
                                      dudtoldc,dvdtoldc,dwdtoldc,ufftc,vfftc,wfftc, &
                                      gxc,gyc,gzc,dgxc,dgyc,dgzc,dgc,indich,indicch,viscc 
!****************************************boundary values of required variables********************************************************
real ,dimension(0:i1c,0:j1c) :: ub,vb,wb,pb,pcb,indb
!*************************************************************************************************************************************
real(mytype) :: time,dt
!******************************************* arrays required for initializing FFT solver**********************************************
real(mytype) ::  wi(itot+15), wj(jtot+15)
real(mytype) ::  wic(itotc+15), wjc(jtotc+15)
! **********************************summation of eigenvalues for poisson equation(used for FFT solver)********************************
real, dimension(imax,jmax) :: xyrt
real, dimension(imaxc,jmaxc) :: xyrtc
!**********************************************coefficients of poisson equation in z direction****************************************
real, dimension(kmax) :: a,b,c
real, dimension(kmaxc) :: ac,bc,cc
!*******************************************two dimensional representation of force points******************************************** 
real, dimension(n_capf,m_capf)::x_cap,y_cap,z_cap
!Note that these two dimensional representations are used only for a number of subroutines but they cannot be used to exchange data
!between cores. The main structure for particles is described below
!********************************two dimensional representation of spherical harmonics points*****************************************
real, dimension(n_cap,m_cap)::x_cap_sph,y_cap_sph,z_cap_sph
!*********************************************total force exerted on flow in each direction*******************************************
real :: forcextot,forceytot,forcefztot
!***************************************************well known variables!*************************************************************
real :: u_bulk,v_bulk,w_bulk,v_bulkc,errorx,errory,errorz
integer :: istep,l,ar
real:: exy1,exy2
real:: area_cap(n_capf,m_capf)!surface area of each element
real:: ulmean,vlmean,wlmean !mean velocity of the force points
!*******************************************time step number for restarting from saved file*******************************************
integer:: begin=0
!*********************************************variables for spherical harmonics part**************************************************
integer::intl,igrida(2),igridb(2),work(lswork),lsvmin,lwkmin,ier
real:: wgtd(n_cap,m_cap),wgtdsin(n_cap,m_cap),wgtdsinf(n_capf,m_capf),surfmetD(n_cap,m_cap),surfmetDf(n_capf,m_capf),thetaN(n_cap2)
real:: phiN(m_cap2),thetaO(n_cap),phiO(m_cap),thetaOf(n_capf),phiOf(m_capf),wsave(lssave),dwork(lsdwork)
real:: conAR(2,2,n_cap2,m_cap2),conARd1(2,2,2,n_cap2,m_cap2),conARd2(2,2,2,2,n_cap2,m_cap2)
real:: BR(2,2,n_cap2,m_cap2),BRd1(2,2,2,n_cap2,m_cap2),BRd2(2,2,2,2,n_cap2,m_cap2)
real:: prstsmax(n_cap2,m_cap2),prstsmin(n_cap2,m_cap2),prstsmax2(n_cap,m_cap),prstsmin2(n_cap,m_cap)
real:: prstsmaxpk(2),prstsminpk(2),Eshr(n_cap2,m_cap2),Edil(n_cap2,m_cap2),Eben(n_cap2,m_cap2)
real:: E_shr(n_cap,m_cap),E_dil(n_cap,m_cap),E_ben(n_cap,m_cap),wtsO(n_cap),wtsOf(n_capf)
real:: E_shrpk(n_cap,m_cap),E_dilpk(n_cap,m_cap),E_benpk(n_cap,m_cap)
real:: crdxD(n_cap,m_cap),crdyD(n_cap,m_cap),crdzD(n_cap,m_cap),crdxR(n_cap,m_cap),crdyR(n_cap,m_cap),crdzR(n_cap,m_cap)
real:: crdxRf(n_capf,m_capf),crdyRf(n_capf,m_capf),crdzRf(n_capf,m_capf)
real:: fx(n_cap,m_cap),fy(n_cap,m_cap),fz(n_cap,m_cap),unitN(3,n_cap2,m_cap2),unitO(3,n_cap,m_cap),unitOf(3,n_capf,m_capf)
real:: evlxlat(nqlat),wgtlat(nqlat),evlxlon(nqlon),wgtlon(nqlon)
real:: evlxlonS(nqlonS),wgtlonS(nqlonS)
real:: fxL(n_cap,m_cap),fyL(n_cap,m_cap),fzL(n_cap,m_cap),acrdR(n_cap,n_cap,3),bcrdR(n_cap,n_cap,3)
real:: acrdD(n_cap,n_cap,3),bcrdD(n_cap,n_cap,3)
real:: NAED(nn3),nnvd(nn3,3),aed(nn3),nvd(ne3,3),nvr(ne3,3),in3(ne3,6),uat(ne3,4),xiT(20),etaT(20),wT(20),ie3(nn3,7)
real:: x3_eq(nn3),y3_eq(nn3),z3_eq(nn3),ccen(3),hf_gscell
!*************************************************************************************************************************************
!***********************************************main structure of particle variables**************************************************

type particle
 real :: x,y,z,integralv,volume,volume0 
 real, dimension(nl) :: xfp,yfp,zfp,    & ! force points
                        dxdt,dydt,dzdt, & ! velocity of force points (not used for inertialess membrane)
                        fxl,fyl,fzl,    & ! Lagrangian force
                        ul,vl,wl,       & ! velocity of fluid at Lagrangian points
                        xfpo,yfpo,zfpo, & ! force points 1 time step before
                        xfpold,yfpold,zfpold, & ! force points 2 time steps before
                        xfp0,yfp0,zfp0,area, &  ! initial valuesof force points and area
                        xn,yn,zn                ! surface normals

! total ammount of reals to be communicated: 6+25*nl
 integer :: mslv
 integer, dimension(8) :: nb
end type particle
type(particle), dimension(npmax) :: ap ! 'a particle' array
type(particle), dimension(npmax) :: sp !send particle array (re-ordering of masters) ! used to exchange particle data between cores     
!as an example ap(5)%fyl(100) means Lagrangian force in 100thpoint of particle number 5 in each core (it is 0 if the core doesn't have particle number 5!)
!volume of Eulerian and Lagrangian cell
real :: dVlagr,dVeul
!
integer :: pmax,npmstr!pmax is sum of number of particles which a specific processor is their master + 
!number of particles which that processor is their slave. In simpler way, it is all particles which a processor has! 
!npmstr is number of particles which a specific processor is their master 
end module mod_common
!***********************************************************some variables about position and neigbours of each processor**********************************************
module mod_common_mpi
use mpi
use decomp_2d
implicit none
integer :: myid,xhalo,yhalo,restartp,rankcw,xhaloc,yhaloc
integer :: comm_cart
!
logical periods(3),reorder
integer error,request,status(MPI_STATUS_SIZE)
integer right,rightfront,front,leftfront,left,leftback,back,rightback
integer, dimension(0:8) :: neighbor
integer, dimension(1:2) :: coords,coordsc,zstartc
real :: boundleftmyid,boundfrontmyid

end module mod_common_mpi
