module mod_param
use decomp_2d
implicit none
!***************************************************************************************************************************************************************************************
!**********************************************************this module sets the required parameters for the problem*********************************************************************
!***************************************************************************************************************************************************************************************

!******************************************************************************************************************************
!*********************************************parallelization parameters*******************************************************
!******************************************************************************************************************************
integer, parameter :: ndims = 2 !dimensions of parallelization for this version of the code,this number should not 
!be touched since the parallelization is only in x and y direction !
integer, dimension(ndims), parameter :: dims = (/1,8/)!the left number shows the parallelization in x direction and
!the right one is for parallelization in y direction total number of cores will be dims(1)*dims(2)
!******************************************************************************************************************************
!******************************************************main domain*************************************************************
!******************************************************************************************************************************
!in this code x , y, and z directions are spanwise, streamwise and wall normal directions respectively
integer, parameter :: itot =128 , jtot = 256, ktot = 128!itot, jtot and ktot are number of grid points in spanwise,
! streamwise and wall normal directions respectively
integer, parameter :: it1 = itot+1, jt1 = jtot+1, kt1 = ktot+1
integer, parameter :: imax = itot/dims(1), jmax = jtot/dims(2),kmax = ktot! number of grid points for each block of processors
integer, parameter :: i1 = imax+1, j1 = jmax+1, k1 = kmax+1
real, parameter :: lx = 5.,ly = 10.,lz = 5.!lx, ly and lz are channel width, length and height respectively
real, parameter :: dxi = itot/lx, dyi = jtot/ly, dzi = ktot/lz
real, parameter :: dx = 1./dxi, dy = 1./dyi, dz = 1./dzi
!WARNING:in this code dx, dy and dz have to be equal i.e. lx/itot=ly/jtot=lz/ktot
real, parameter :: diam = 2. ! lengths are scaled with the particle radius
!******************************************************************************************************************************
!*********************************************overlapped domain****************************************************************
!******************************************************************************************************************************
integer, parameter ::ifoverlap=0 !if zero, the code doesn't use overlap mesh (finer grid near bottom wall)
integer, parameter :: nf=16,nr=2,ktotc =nf*nr,itotc=itot*nr,jtotc=jtot*nr!nf is the number of cells near wall in main mesh
!which we want to make them finer, nr is the ratio of main mesh size to fine mesh size
!name of variables are similar to those for main mesh with an additional letter "c" in front!
integer, parameter :: it1c = itotc+1,jt1c = jtotc+1,kt1c = ktotc+1
integer, parameter :: imaxc = itotc/dims(1), jmaxc = jtotc/dims(2),kmaxc = ktotc
integer, parameter :: i1c = imaxc+1,j1c = jmaxc+1,k1c = kmaxc+1
real, parameter :: lzc = nf*dz
real, parameter :: dxic = itotc/lx,dyic = jtotc/ly,dzic = ktotc/lzc
real, parameter :: dxc = 1./dxic,dyc = 1./dyic,dzc = 1./dzic
!******************************************************************************************************************************
!*********************************************number of particles**************************************************************
!******************************************************************************************************************************
integer, parameter :: np = 1, npmax = 1 !np is total number of particles and npmax is maximum allowed number of particles in
!each processor
!*********************************************creating PI number***************************************************************
real, parameter :: pi = acos(-1.)
!******************************************************************************************************************************
!************************************************* centrifuge parameters*******************************************************
!******************************************************************************************************************************
real,parameter :: thetaincd= 70.,thetainc=(thetaincd*pi)/180. ! thetaincd is (90)-inclination angle of centrifuge in degrees
!thetainc is representation of thetaincd in Radians
!******************************************************************************************************************************
!**************************************************flow parameters*************************************************************
!******************************************************************************************************************************
real, parameter ::Rep=.0016,ta=1160000.,Ca=.1,we=Rep*Ca,lambda=1.!Rep: particle Reynolds number, ta: Taylor number,
!Ca: Capillary number, we: Weber number, lambda: viscosity ratio
real, parameter :: visc0=1./Rep 
real,parameter:: shear = 1.!shear rate for Couette flow
!******************************************************************************************************************************
!************************************************** collision parameters*******************************************************
!******************************************************************************************************************************
real,parameter :: col= 5.68533*(10.**-2.),betacol=15.36,rad0=2.5*dz!col is non dimensional Beta*De in the Morse potential model
!for interaction with wall, betacol is non dimensional Beta, rad0 is zero force distance
!******************************************************************************************************************************
!****************************************************particle parameters*******************************************************
!******************************************************************************************************************************
character(len=1), parameter :: ifshape = "s" !"s" means spherical capsule and "b" stands for biconcave red blood cell shape
real,parameter:: arbc=.5,crbc=.6 !parameters for red blood cell shape obtained from the paper by Angelov & Mladenov
!******************************************************************************************************************************
!For the sake of efficiency, two sets of Lagrangian points have been used. One is for spreading force and normals to
!Eulerian cells (nl) and the other is for representation of capsule shape with spherical harmonics(nls). Normally, the number
!of force points (nl) is larger than those for spherical harmonic points (nls). all calculations regarding force and normals
!will be done in spherical harmonic mesh and then a spectral interpolation will be used to transfer them to force points
!******************************************************************************************************************************  
integer, parameter :: nl =4608,n_capf=48,m_capf=96 !n_capf and m_capf is number of force points on capsule in latitudinal and
!longitudinal directions respectively so nl=n_capf*m_capf  
integer,parameter:: nls=1152,n_cap=24,m_cap=48,n_cap2=48,m_cap2=96!n_cap and m_cap is number of spherical harmonics points on
!capsule in latitudinal and longitudinal directions respectively so nls=n_cap*m_cap
!n_cap2 and m_cap2 are for the mesh used for dealiasing which is used to remove aliasing errors. They should be 2 times of 
!those for spherical harmonics points
integer,parameter:: nnlat=n_cap2,nnlon=m_cap2,nnlat2=n_cap,nnlon2=m_cap!the same as spherical harmonic points with a 
!different name which are used in the modules for spherical harmonics
!******************************************************************************************************************************
!at the moment the code is using Neo-Hookian model. To switch to Skalak model, go to makefile and change -DNONHOOK to
!-DSKALAK in 7 and 8th line
real,parameter::    gscell =1./we,CSK=50.! gscell is non dimensionalized surface shear modulus and CSK is the constant used in 
!Skalak model
integer,parameter:: ibend = 1 ! if 0 then bending is zero
real,parameter::    cbend = 0.3, ebend =cbend*gscell !non dimensionalized surface bending modulus
integer,parameter:: idealias=1!if 0 the code doesn't do dealiasing
!******************************************************************************************************************************
!the parameters you see below are for creating workspace for the subroutines for spherical harmonics please never touch them!
!******************************************************************************************************************************
integer,parameter::  lldworkO=n_cap*(n_cap+4),lldworkN=n_cap2*(n_cap2+4),lldworkOf=n_capf*(n_capf+4)
integer,parameter::  nqlat = 2,nqlatS = 4
integer,parameter::  nqlon = 2*nqlat,nqlonS = 2*nqlatS
integer,parameter::  nn3=642,ne3=1280 
integer,parameter::  la1 = min0(n_cap,(m_cap+2)/2),la2 = (n_cap+1)/2 ,lb1 = min0(n_capf,(m_capf+2)/2),lb2 = (n_capf+1)/2 
integer,parameter::  lwa = n_cap*(2*la2+3*la1-2)+3*la1*(1-la1)/2+m_cap+15
integer,parameter::  lwb = 2*n_capf*lb2+3*((lb1-2)*(n_capf+n_capf-lb1-1))/2+m_capf+15
integer,parameter::  lssave=2*(lwa+lwb)
integer,parameter::  nlat = max0(n_cap,n_capf), nlon = max0(m_cap,m_capf)
integer,parameter::  l1s = min0(nlat,(nlon+2)/2), l2s = (nlat+1)/2
integer,parameter::  lswork=nlat*(4*l1s+nlon+2*nlat+4)+3*((l1s-2)*2*(2*nlat-l1s-1))/2 
integer,parameter::  lsdwork=nlon*(nlon+4)
!******************************************************************************************************************************
!**********************************************setting flow solver schemes*****************************************************
!******************************************************************************************************************************
! in this code, the Adams bashforth method has been used for discretization of convective term. For diffusive term there is a
!choice between crank nicolson and FFT(fast fourier transform) solver. FFT is more fast and accurate but don't use it when 
!having viscosity contrast 
character(len=1),parameter ::fsolver="f" !for Crank Nicolson scheme put "c" and for FFT put"f"
!******************************************************************************************************************************                                          
!******************************************************************************************************************************
integer, parameter :: send_real = 6+25*nl!amount of data for each particle to be transfered between the cores
!******************************************************************************************************************************
!************************************************flow initialization***********************************************************
!******************************************************************************************************************************
! type of initial velocity field (see init.f90)
! iniu = 'cou' --> plane Couette flow
!      = 'poi' --> plane Poiseuille flow
!      = 'zer' --> zero velocity everywhere
!      = 'log' --> logarithmic profile + random noise
 character(len=3), parameter :: iniu = 'zer' 
!******************************************************************************************************************************
!**************************************************setting free slip***********************************************************
!******************************************************************************************************************************
! type of BC z = lz is set to free-slip if
! isfreeslip is true
logical, parameter :: isfreeslip = .false.
!******************************************************************************************************************************
real, parameter :: bulk_v_sup =0.!mean velocity used for representation of Poiseuille flow. it is =1 in case
real, parameter :: radius = 1., offset = ( sqrt(3.*(1.5**2)) )/dxi + 0.01/dxi! offset is used for parallelization of particles
!to distinguish slave of each particle
!******************************************************************************************************************************
!****************************************************writing data**************************************************************
!******************************************************************************************************************************
character(len=5), parameter :: datadir = 'data/'!directory of output files 
!******************************************************************************************************************************
! after "iout2d" time steps, the code writes a 2D plot in mid yz plane
! after "ioutfld" time steps, the code writes data files for flow and particles which is necessary to restart from saved file 
integer,parameter :: iout2d = 100 ,ioutfld = 100
end module mod_param
