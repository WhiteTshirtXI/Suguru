!                                                                                                                                                                                        
!                                                                                                                                                                                        
! IIIIIIIIII                         tttt                            RRRRRRRRRRRRRRRRR   PPPPPPPPPPPPPPPPP                                               tttt             SSSSSSSSSSSSSSS 
! I::::::::I                      ttt:::t                            R::::::::::::::::R  P::::::::::::::::P                                           ttt:::t           SS:::::::::::::::S
! I::::::::I                      t:::::t                            R::::::RRRRRR:::::R P::::::PPPPPP:::::P                                          t:::::t          S:::::SSSSSS::::::S
! II::::::II                      t:::::t                            RR:::::R     R:::::RPP:::::P     P:::::P                                         t:::::t          S:::::S     SSSSSSS
!   I::::Innnn  nnnnnnnn    ttttttt:::::ttttttt        eeeeeeeeeeee    R::::R     R:::::R  P::::P     P:::::Paaaaaaaaaaaaa  rrrrr   rrrrrrrrr   ttttttt:::::ttttttt    S:::::S            
!   I::::In:::nn::::::::nn  t:::::::::::::::::t      ee::::::::::::ee  R::::R     R:::::R  P::::P     P:::::Pa::::::::::::a r::::rrr:::::::::r  t:::::::::::::::::t    S:::::S            
!   I::::In::::::::::::::nn t:::::::::::::::::t     e::::::eeeee:::::eeR::::RRRRRR:::::R   P::::PPPPPP:::::P aaaaaaaaa:::::ar:::::::::::::::::r t:::::::::::::::::t     S::::SSSS         
!   I::::Inn:::::::::::::::ntttttt:::::::tttttt    e::::::e     e:::::eR:::::::::::::RR    P:::::::::::::PP           a::::arr::::::rrrrr::::::rtttttt:::::::tttttt      SS::::::SSSSS    
!   I::::I  n:::::nnnn:::::n      t:::::t          e:::::::eeeee::::::eR::::RRRRRR:::::R   P::::PPPPPPPPP      aaaaaaa:::::a r:::::r     r:::::r      t:::::t              SSS::::::::SS  
!   I::::I  n::::n    n::::n      t:::::t          e:::::::::::::::::e R::::R     R:::::R  P::::P            aa::::::::::::a r:::::r     rrrrrrr      t:::::t                 SSSSSS::::S 
!   I::::I  n::::n    n::::n      t:::::t          e::::::eeeeeeeeeee  R::::R     R:::::R  P::::P           a::::aaaa::::::a r:::::r                  t:::::t                      S:::::S
!   I::::I  n::::n    n::::n      t:::::t    tttttte:::::::e           R::::R     R:::::R  P::::P          a::::a    a:::::a r:::::r                  t:::::t    tttttt            S:::::S
! II::::::IIn::::n    n::::n      t::::::tttt:::::te::::::::e        RR:::::R     R:::::RPP::::::PP        a::::a    a:::::a r:::::r                  t::::::tttt:::::tSSSSSSS     S:::::S
! I::::::::In::::n    n::::n      tt::::::::::::::t e::::::::eeeeeeeeR::::::R     R:::::RP::::::::P        a:::::aaaa::::::a r:::::r                  tt::::::::::::::tS::::::SSSSSS:::::S
! I::::::::In::::n    n::::n        tt:::::::::::tt  ee:::::::::::::eR::::::R     R:::::RP::::::::P         a::::::::::aa:::ar:::::r                    tt:::::::::::ttS:::::::::::::::SS 
! IIIIIIIIIInnnnnn    nnnnnn          ttttttttttt      eeeeeeeeeeeeeeRRRRRRRR     RRRRRRRPPPPPPPPPP          aaaaaaaaaa  aaaarrrrrrr                      ttttttttttt   SSSSSSSSSSSSSSS 
! 
! InteRPTartS - Interface-Resolved Particle Simulations
! 
! Code for channel (full or half) transport
!
! Contributors: Arash Alizad Banaei, Pedro Costa & Wim-Paul Breugem
!  p.simoes.costa@gmail.com / w.p.breugem@tudelft.nl
!
! modifications for deformable particles, flow solver and adaptive mesh: Arash Alizad Banaei arash2@mech.kth.se 
!
! Last modified: April 2016
!
!**********************************************************
program interpt
use mod_param
use mod_common
use mod_common_mpi
use mod_initmpi
use mod_init
use mod_bound
use mod_boundc
use mod_chkdiv
use mod_chkdivc
use mod_chkdt
use mod_loadflds
use mod_loadpart
use mod_crnk
use mod_crnkc
use mod_fftuvw
use mod_fftuvwc
use mod_initsolver
use mod_fillps
use mod_fillpsc
use mod_solver
use mod_solverc
use mod_correc
use mod_correcc
use mod_collisions
use mod_initparticles
use mod_kernel
use mod_interp_spread
use mod_interp_spreadc
use mod_forcing
use mod_Update_Pos
use mod_surf_vol
use mod_sph
use mod_indicator
use mod_interpolation
use mod_output
use mod_outputc
implicit none
integer :: iopt,nstep,nsub,i,j,pp
real :: dtmax
real, dimension(0:i1,0:j1,0:k1) :: durkold,dvrkold,dwrkold
real, dimension(0:i1c,0:j1c,0:k1c) :: durkoldc,dvrkoldc,dwrkoldc
real, target, dimension(0:i1,0:j1,0:k1) :: p
real, target, dimension(0:i1c,0:j1c,0:k1c) :: pc
!real, pointer, dimension(:,:,:) :: pp ! uncomment iff solver2d is used
real :: maxerror,avererror,maxerrorsum,avererrorsum
integer :: ibmiter
integer :: k
real :: norm,norm_all,sumk
real :: v_bulk_all,v_bulk_allc
!character(len=4) :: number
real :: comp_time,comp_time_avrg,comp_time_min,comp_time_max ! for mpi timer
 character(10)::filename
 character(50)::fileplace

 character*3 number
 character*7 filenumber
!
!***************************************************************************************************************************************************************************************************
!*******************************************************************************************decomposition of domains********************************************************************************
!***************************************************************************************************************************************************************************************************
iopt = 0
if(iopt.eq.1) then
  call MPI_INIT(error)
  call decomp_2d_init(itot, jtot, ktot, 0,0)
  call MPI_FINALIZE(error)
  stop
else
  call initmpi
endif
!**************************************************************************************number of time steps*****************************************************************************************
nstep = 5000000  
!***************************************************************************************************************************************************************************************************
!**********************************************************************************************initialize*******************************************************************************************
!***************************************************************************************************************************************************************************************************
!begin is the number of time step to start from. it can be altered in the module "common"
if(begin.eq.0) then
!********************************to initialize flow**********************************
  call init
!********************************to initialize particles*****************************
!************************************************************************************
  call initparticles
!initial position of particle center (for 1 particle) an be specified in the module "initparticles"

  time = 0.
!  call chkdt(dtmax)
!*******************************size of time steps***********************************
  dt = .0001
!putting current values of force points into old ones for first iteration
!(pp) represents number of particle and the second index is number of force points 
  do pp=1,pmax
  ap(pp)%xfp0(:) = ap(pp)%xfp(:)
  ap(pp)%yfp0(:) = ap(pp)%yfp(:)
  ap(pp)%zfp0(:) = ap(pp)%zfp(:)
  enddo
else
  !$omp workshare
  unew(:,:,:) = 0.
  vnew(:,:,:) = 0.
  wnew(:,:,:) = 0.
  pnew(:,:,:) = 0.
  unewc(:,:,:) = 0.
  vnewc(:,:,:) = 0.
  wnewc(:,:,:) = 0.
  uoc(:,:,:) = 0.
  voc(:,:,:) = 0.
  woc(:,:,:) = 0.
  uooc(:,:,:) = 0.
  vooc(:,:,:) = 0.
  wooc(:,:,:) = 0.
  pnewc(:,:,:) = 0.
  !$omp end workshare
!********************loading flow and particles in the case of starting from saved file**********************************
  call loadflds(0,begin)
  call loadpart(0,begin)
!*********************only initializes spherical harmonics if start from saved file *************************************
  call initparticles

  if (myid .eq. 0) write(6,*) 'nr steps at beginning simulation = ',begin
endif
!
!***********************initializes coefficients for the equations which will be solved using FFT************************ 
call initsolver(1)
!pp => p(1:imax,1:jmax,1:kmax) ! uncomment iff solver2d is used
!*******************************************imposing boundary conditions*************************************************
call bounduvw(unew,vnew,wnew)
call boundp(pnew)
call chkdiv!checking RHS of continuity equation


!$omp workshare
durkold(:,:,:) = 0.
dvrkold(:,:,:) = 0.
dwrkold(:,:,:) = 0.
!$omp end workshare
!*********************************************************************************end of initialize************************************************************************************************

!***********************************************************************************************************************************************************************************************
!******************************************************************************************main loop********************************************************************************************
!***********************************************************************************************************************************************************************************************
do istep = begin+1,nstep
  time = time + dt

  if (myid.eq.0) write(6,*) 'time = ', time, 'istep = ', istep
 
!**************************putting current values of force points into old ones *****************************        
    do pp=1,pmax
    ap(pp)%xfpold(:) = ap(pp)%xfpo(:)
    ap(pp)%yfpold(:) = ap(pp)%yfpo(:)
    ap(pp)%zfpold(:) = ap(pp)%zfpo(:)
    ap(pp)%xfpo(:) = ap(pp)%xfp(:)
    ap(pp)%yfpo(:) = ap(pp)%yfp(:)
    ap(pp)%zfpo(:) = ap(pp)%zfp(:)
    enddo

!**************************putting current values of flow variables into old ones ***************************
    if (istep==begin+1) then
    uoo=unew
    voo=vnew
    woo=wnew
    if (ifoverlap==1) then
    uooc=unewc
    vooc=vnewc
    wooc=wnewc
    endif
    endif


    uoo=uo
    voo=vo
    woo=wo
    if (ifoverlap==1) then
    uooc=uoc
    vooc=voc
    wooc=woc
    endif

    uo=unew
    vo=vnew
    wo=wnew
    if (ifoverlap==1) then
    uoc=unewc
    voc=vnewc
    woc=wnewc
    endif

    if (ifoverlap==1) then
    if (istep==begin+1) then
    call interpolation
    call bounduvwc(unewc,vnewc,wnewc)
    endif
    endif
!**************************************interpolating velocities from Eulerian cells to lagrangian points****************************** 
    call eulr2lagr

    if (ifoverlap==1) call eulr2lagrc
!********************************************************computing elastic forces*****************************************************
    call complagrforces
!********************************************adding collision forces to Lagrangian forces********************************************* 
    call collisions
!*********************************************spreading elastic force back to Eulerian grid******************************************* 
    call lagr2eulr

    if (ifoverlap==1) call lagr2eulrc 
!*******************************************************computing indicator function**************************************************  
    call indicator
!********************************updating position of Lagrangian points and renewing parallelization for particles******************** 
    call Update_Pos
!*********************************************computing deformation parameter for particles******************************************* 
    call surf_vol



!***********************************************************************************************************************************************************************************************
!******************************************************************************************flow solver******************************************************************************************
!***********************************************************************************************************************************************************************************************
 if (fsolver=="c") then
!****************************************************solving for flow for main mesh***************************************************
    errorx=1.
    errory=1.
    errorz=1.
    ar=0
!***************************************************Crank Nicolson step for main mesh*************************************************
    comp_time = MPI_WTIME()!flow timer
    do while(errorx>=.0001.or.errory>=.0001.or.errorz>=.0001)
    ar=ar+1
    upi=unew
    vpi=vnew
    wpi=wnew


   !$omp workshare
    v_bulk =sum(vnew(1:imax,1:jmax,1:kmax))
   !$omp end workshare
    call mpi_allreduce(v_bulk,v_bulk_all,1,mpi_real8,mpi_sum,comm_cart,error)
    v_bulk=v_bulk_all/(itot*jtot*kmax)

    if(myid.eq.0) then

      write(6,*) 'v_bulk',v_bulk
      write(6,*) 'residuals for velocities in main mesh',errorx,errory,errorz      

!      open(29,file=datadir//'bulk_velocity.txt',position='append')
!      write(29,'(4E16.8)') 1.*istep,dt,bulk_v_sup*sin(time),v_bulk
!      close(29)

     endif

    
      if(myid.eq.0) write(6,*) 'GS step for main mesh',ar
      call crnk(durkold,dvrkold,dwrkold)
    
    
!********forcing to keep the vulk velocity at constant value (only use this for Poiseuille flow)********
!    call updtintermediatevel
!*******************************************************************************************************
  
   !
    call bounduvw(unew,vnew,wnew)

    enddo
!***********************************************end of Crank Nicolson step for main mesh********************************************

!**********************************************computing pressure correction for main mesh******************************************
    call initsolver(1)
!***************************************computing source term of pressure correction equation***************************************
    call fillps(p)
!    call solver2d(pp)
!******************************************solving pressure correction equation using FFT*******************************************
    call solver1d(p,1)
!*****************************************************boundary condition for pressure***********************************************
    call boundp(p)
!*******************************************************correctiong velocities for main mesh****************************************
    call correc(p)
    call bounduvw(unew,vnew,wnew)
!$omp workshare
!******************************************************correcting pressure for main mesh********************************************
    pnew(:,:,:) = pnew(:,:,:) + p(:,:,:)
!$omp end workshare
    call boundp(pnew)
   
!*******************************interpolating velocities from main mesh to upper boundary of  overlapped mesh*********************** 
    if (ifoverlap==1) then
    call interpolation
    call bounduvwc(unewc,vnewc,wnewc)
!*************************************************solving for flow in overlapped  mesh**********************************************
!***********************************************************************************************************************************
!***********************************************the same procedure for overlapped mesh!*********************************************
    errorx=1.
    errory=1.
    errorz=1.
    ar=0
    do while(errorx>=.0001.or.errory>=.0001.or.errorz>=.0001)
    ar=ar+1
    upic=unewc
    vpic=vnewc
    wpic=wnewc


   !$omp workshare
    v_bulkc =sum(vnewc(1:imaxc,1:jmaxc,1:kmaxc))
   !$omp end workshare
    call mpi_allreduce(v_bulkc,v_bulk_allc,1,mpi_real8,mpi_sum,comm_cart,error)
    v_bulkc=v_bulk_allc/(itotc*jtotc*kmaxc)

    if(myid.eq.0) then

      write(6,*) 'v_bulk for fine mesh',v_bulkc
      write(6,*) 'residuals for velocities in fine mesh',errorx,errory,errorz      

     endif

 
      if(myid.eq.0) write(6,*) 'GS step for fine mesh',ar
      call crnkc(durkoldc,dvrkoldc,dwrkoldc)
    
    
!!!!!!!!!!!!!!!!!!!!!!!!Forcing!!!!!!!!!!!!!!!!!!!!!!

!    call updtintermediatevel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !
    call bounduvwc(unewc,vnewc,wnewc)

    enddo

    call initsolver(1)
    call fillpsc(pc)
!    call solver2d(pp)
    call solver1dc(pc,1)
    call boundpcorrc(pc)
    call correcc(pc)
    call bounduvwc(unewc,vnewc,wnewc)
!$omp workshare
    pnewc(:,:,:) = pnewc(:,:,:) + pc(:,:,:)
!$omp end workshare
    call boundpc(pnewc)
    endif

endif
!***************************************solving for flow using fft for all equations**********************************************
!*********************************************************************************************************************************
!the procedure is the same as before except for the momentum equations, it calls fftuvw for main mesh and fftuvwc for overlapped 
!mesh, inside of these subroutnes, coefficients and source term of momentum equations become prepared and given to the subroutines
!solver and solverc to solve for main and overlapped mesh respectively

   if (fsolver=="f") then

   !$omp workshare
    v_bulk =sum(vnew(1:imax,1:jmax,1:kmax))
   !$omp end workshare
    call mpi_allreduce(v_bulk,v_bulk_all,1,mpi_real8,mpi_sum,comm_cart,error)
    v_bulk=v_bulk_all/(itot*jtot*kmax)

    if(myid.eq.0) then

      write(6,*) 'v_bulk',v_bulk

     endif



    comp_time = MPI_WTIME()
    call fftuvw(durkold,dvrkold,dwrkold)

    call bounduvw(unew,vnew,wnew)


call interpolation
    

    call initsolver(1)
    call fillps(p)
    call solver1d(p,1)
    call boundp(p)
    call correc(p)
    call bounduvw(unew,vnew,wnew)
!$omp workshare
    pnew(:,:,:) = pnew(:,:,:) + p(:,:,:)
!$omp end workshare
    call boundp(pnew)

    if (ifoverlap==1) then  

    
    call bounduvwc(unewc,vnewc,wnewc)


   !$omp workshare
    v_bulkc =sum(vnewc(1:imaxc,1:jmaxc,1:kmaxc))
   !$omp end workshare
    call mpi_allreduce(v_bulkc,v_bulk_allc,1,mpi_real8,mpi_sum,comm_cart,error)
    v_bulkc=v_bulk_allc/(itotc*jtotc*kmaxc)

    if(myid.eq.0) then

      write(6,*) 'v_bulk for fine mesh',v_bulkc     

     endif

    call fftuvwc(durkoldc,dvrkoldc,dwrkoldc)


    call bounduvwc(unewc,vnewc,wnewc)

    call initsolver(1)
    call fillpsc(pc)
    call solver1dc(pc,1)
    call boundpcorrc(pc)
    call correcc(pc)
    call interpolation
    call bounduvwc(unewc,vnewc,wnewc)
!$omp workshare
    pnewc(:,:,:) = pnewc(:,:,:) + pc(:,:,:)
!$omp end workshare
    call boundpc(pnewc)
    endif


    endif
!***********************************************************************************************************************************************************************************************
!************************************************************************************end of flow solver*****************************************************************************************
!***********************************************************************************************************************************************************************************************
  do pp=1,pmax
  if (mod(istep,100).eq.0) then
  
    if (ap(pp)%mslv .gt. 0) then !if the core is master of particle

!*******************************************************writing Lagrangian points*****************************************************
       write(number,'(i3.3)') pp
       write(filenumber,'(i7.7)') istep
       open(42,file=datadir//'Capsule'//number//'step'//filenumber//'.txt')
       do j=1,NL
        write(42,'(3E15.7)')  ap(pp)%yfp(j),ap(pp)%zfp(j),ap(pp)%xfp(j)
       enddo
       close(42)
!**********************************writing mean values of position and velocities*****************************************************
      open(58,file=datadir//'center_velocity.txt',position='append')
      write(58,'(7E16.8)') time,ap(1)%x,ap(1)%y,ap(1)%z,ulmean,vlmean,wlmean
      close(58)
!********************************************writing deformation parameter************************************************************
!exy1 is deformation parameter using moment of inertia and exy2 is obtained from  maximum radius and minimum radius
      open(59,file=datadir//'epsilonxy.txt',position='append')
      write(59,'(3E16.8)') time,exy1,exy2
      close(59)
  
     endif 

   endif

  enddo

    
  if (mod(istep,100).eq.0) then
!*************************checking source term of pressure correction equation to ensure that solution is not diverging!*************
    call chkdiv
    if (ifoverlap==1) call chkdivc
!    call chkdt(dtmax)
    dt = .0001


    if (myid .eq. 0) then
      write(6,*) 'dtmax = ', dtmax, ' dt = ', dt
      open(29,file=datadir//'timestep.txt',position='append')
      write(29,'(3E16.8)') 1.*istep,time,dt
      close(29)
    endif

 endif
!*********************************************************writing 2D data for flow**************************************************
    if (mod(istep,iout2d).eq.0) then
    call post2d(istep)
    if (ifoverlap==1) call post2dc(istep)
    endif
!****************************************************saving flow and particle data**************************************************
  if (mod(istep,ioutfld).eq.0) then
    call loadflds(1,istep)
    call loadpart(1,istep)
    endif
!****************************************************calculating computational time*************************************************
  comp_time = MPI_WTIME()-comp_time
  call mpi_allreduce(comp_time,comp_time_avrg,1,mpi_real8,mpi_sum,comm_cart,error)
  call mpi_allreduce(comp_time,comp_time_min,1,mpi_real8,mpi_min,comm_cart,error)
  call mpi_allreduce(comp_time,comp_time_max,1,mpi_real8,mpi_max,comm_cart,error)
  if(myid.eq.0) write(6,'(A,3E16.8)') 'Avrg, min & max elapsed time = ', &
                            comp_time_avrg/product(dims),comp_time_min,comp_time_max

!***********************************************************************************************************************************************************************************************
!*****************************************************************************end of main loop**************************************************************************************************
!***********************************************************************************************************************************************************************************************
enddo
!
if(myid.eq.0) write(6,*) '*** Fim ***'
!
call decomp_2d_finalize
call MPI_FINALIZE(error)
!
stop
end program interpt
