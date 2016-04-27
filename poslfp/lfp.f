      program init_force_points
      use mod_param
      use mod_loadd
!
!     Program computes the positions of the Lagrangian force points wrt center of sphere.
!
      implicit none
      integer l,q,j,i,k
      real rn
      real absforce
      real thetanew,phinew
!     real forcerad(1:NL)
      real forcephi(1:NL),forcetheta(1:NL)
      real forcethetaold(1:NL),forcephiold(1:NL)
      real difvectorx,difvectory,difvectorz
      real lengthdifvector
      real ethetax,ethetay,ethetaz
      real ephix,ephiy,ephiz
      real erx,ery,erz
      real forcephimean,forcephivar
      real forcethetamean,forcethetavar
      real xfp(1:NL),yfp(1:NL),zfp(1:NL)
      integer begin
      real dt
      parameter(dt=1.e-4) !1.e-3)
      real dist1max,dist2max,dist
      integer lmin1,lmin2
      real phi(1:NL),theta(1:NL)
      integer restart
      parameter(restart=0)
      real radfp
      parameter(radfp=radius-retraction) !IMPORTANT, lfp's positioned at slightly lower radius!
      integer qmax
      parameter(qmax=50)
      real distr(1:qmax),lowbound,upbound,sumq
      real sumsquares
      integer rkiter
!
!     Initialization     
!
      write(6,*) 'NL = ',NL
      write(6,*) 'Volume occupied by 1 lfp = ', 
     $ (4./3.)*picon*( (radfp+0.5/dxi)**3. - (radfp-0.5/dxi)**3. )/(1.*NL)
      write(6,*) 'Volume of Eulerian grid cell = ', 1./dxi/dyi/dzi
      do q=1,qmax
        distr(q) = 0.
      enddo
      sumq = 0
      rn = 10.
      do l=1,NL
        call random_number(rn)
        call random_number(rn)
        call random_number(rn)
        q=0
        lowbound = -1./qmax
        upbound  = 0.
 777    q=q+1
        lowbound = lowbound + 1./qmax
        upbound  = upbound + 1./qmax
        if ( (rn .gt. lowbound) .and. (rn .lt. upbound) ) then
          distr(q) = distr(q)+1.
          sumq     = sumq + 1
        else
          go to 777
        endif
        phi(l)      = rn*2.*picon !random value between 0 and 2*pi         
        call random_number(rn)
        call random_number(rn)
        call random_number(rn)
        q=0
        lowbound = -1./qmax
        upbound  = 0.
 888    q=q+1
        lowbound = lowbound + 1./qmax
        upbound  = upbound + 1./qmax
        if ( (rn .gt. lowbound) .and. (rn .lt. upbound) ) then
          distr(q) = distr(q)+1.
          sumq     = sumq + 1
        else
          go to 888
        endif
        theta(l)    = rn*picon    !random value between 0 and pi         
        xfp(l)      = radfp*sin(theta(l))*cos(phi(l))
        yfp(l)      = radfp*sin(theta(l))*sin(phi(l))
        zfp(l)      = radfp*cos(theta(l))
        forcetheta(l)    = 0.
        forcephi(l)      = 0.
      enddo
      if (sumq .ne. 2*NL) then
        write(6,*) 'Error! Failure in binning of random number.'
        write(6,*) 'Program aborted...'
        stop
      endif
      open(18,file=datadir//'distribution.txt')
      do q=1,qmax
        write(18,'(3E16.8)') (1.*q-0.5)/(1.*qmax),distr(q),sumq/(1.*qmax)
      enddo
      close(18)
!     data to file
      open(42,file=datadir//'lagrangianforcepoints_init')
      write(42,*) 'VARIABLES = "lfpx","lfpy","lfpz","theta","phi","radfp"'
      write(42,*) 'ZONE T="Zone1"',' I=',NL,', F=POINT'
      write(42,*) ''
      do l=1,NL
        write(42,'(6E16.8)') xfp(l),yfp(l),zfp(l),theta(l),phi(l),radfp
      enddo
      close(42) 
!     data for sphere
      open(42,file=datadir//'datasphere')
      write(42,*) 'VARIABLES = "x","y","z","r"'
      write(42,*) 'ZONE T="Zone1"',' I=',imax,' J=',jmax,' K=',kmax,', F=POINT'
      write(42,*) ''
      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            write(42,'(4E16.8)') (i-imax/2)/dxi,(j-jmax/2)/dyi,(k-kmax/2)/dzi,
     $ (sqrt( ((i-imax/2)/dxi)**2. + ((j-jmax/2)/dyi)**2. + ((k-kmax/2)/dzi)**2. ))/radfp
          enddo
        enddo
      enddo
      close(42)
!
!     Compute 'Coulomb force' acting on every charged particle
!
      begin = 0
      if (restart .eq. 1) then 
        call loadd(0,begin,phi,theta)
        do l=1,NL
          xfp(l)      = radfp*sin(theta(l))*cos(phi(l))
          yfp(l)      = radfp*sin(theta(l))*sin(phi(l))
          zfp(l)      = radfp*cos(theta(l))
        enddo
      endif
      write(6,*) 'begin = ',begin

      do j=begin+1,50000!100000
        write(6,*) 'Iteration step = ',j
        rkiter = 0
 999    rkiter = rkiter+1
!$omp parallel default(shared)
!$omp& private(l,q,ethetax,ethetay,ethetaz,ephix,ephiy,ephiz)
!$omp&private(difvectorx,difvectory,difvectorz,sumsquares,absforce,lengthdifvector)
!$omp do
        do l=1,NL
!         unit vector in theta direction
          ethetax    =  cos( theta(l) )*cos( phi(l) )
          ethetay    =  cos( theta(l) )*sin( phi(l) )
          ethetaz    = -sin( theta(l) )
!         unit vector in phi direction
          ephix      = -sin( phi(l) )
          ephiy      =  cos( phi(l) )
          ephiz      =  0.
!         unit vector in radial direction
!         erx        = xfp(l)/radfp
!         ery        = yfp(l)/radfp
!         erz        = zfp(l)/radfp
          forcetheta(l) = 0.
          forcephi(l)   = 0.
!         forcerad(l)   = 0.
          do q=1,NL
            if (q .ne. l) then
              difvectorx  = xfp(l)-xfp(q) 
              difvectory  = yfp(l)-yfp(q)
              difvectorz  = zfp(l)-zfp(q)
              sumsquares  = difvectorx*difvectorx + difvectory*difvectory + difvectorz*difvectorz
              absforce   = 1./( sumsquares ) !1/(distance**2)
!             normalized difference vector
              lengthdifvector = sqrt( sumsquares )
              difvectorx = difvectorx/lengthdifvector 
              difvectory = difvectory/lengthdifvector
              difvectorz = difvectorz/lengthdifvector
!             force component in theta direction              
              forcetheta(l) = forcetheta(l) + (ethetax*difvectorx + ethetay*difvectory + ethetaz*difvectorz)*absforce
!             force component in phi direction              
              forcephi(l) = forcephi(l) + (ephix*difvectorx + ephiy*difvectory + ephiz*difvectorz)*absforce
!             force component in radial direction              
!             forcerad(l) = forcerad(l) + (erx*difvectorx + ery*difvectory + erz*difvectorz)*absforce
            endif
          enddo
          forcetheta(l) = forcetheta(l)/(1.*NL-1.) !averaged force 
          forcephi(l)   = forcephi(l)/(1.*NL-1.)   !averaged force
!         forcerad(l)   = forcerad(l)/(1.*NL-1.)   !averaged force
!         The forces are expected to scale with dxi**2. Since NL is proportional to dxi**2, 
!         they scale approximately with NL. By dividing the forces by NL, the forces become independent
!         of the resolution. This then implies that the time step is insensitive to the resolution.
        enddo
!$omp end do
!$omp end parallel
!       RK3 scheme
        if (rkiter .eq. 1) then
!$omp parallel default(shared)
!$omp& private(l)
!$omp do
          do l=1,NL
            theta(l)         = theta(l) + dt*(32./60.)*forcetheta(l)
            phi(l)           = phi(l) + dt*(32./60.)*forcephi(l)
            xfp(l)           = radfp*sin(theta(l))*cos(phi(l))
            yfp(l)           = radfp*sin(theta(l))*sin(phi(l))
            zfp(l)           = radfp*cos(theta(l))
            forcethetaold(l) = forcetheta(l)
            forcephiold(l)   = forcephi(l)
          enddo
!$omp end do
!$omp end parallel
          go to 999
        endif
        if (rkiter .eq. 2) then
!$omp parallel default(shared)
!$omp& private(l)
!$omp do
          do l=1,NL
            theta(l)         = theta(l) + dt*( (25./60.)*forcetheta(l) - (17./60.)*forcethetaold(l) )
            phi(l)           = phi(l) + dt*( (25./60.)*forcephi(l)-(17./60.)*forcephiold(l) )
            xfp(l)           = radfp*sin(theta(l))*cos(phi(l))
            yfp(l)           = radfp*sin(theta(l))*sin(phi(l))
            zfp(l)           = radfp*cos(theta(l))
            forcethetaold(l) = forcetheta(l)
            forcephiold(l)   = forcephi(l)
          enddo
!$omp end do
!$omp end parallel
          go to 999
        endif
        if (rkiter .eq. 3) then
!$omp parallel default(shared)
!$omp& private(l)
!$omp do
          do l=1,NL
            theta(l)         = theta(l) + dt*( (45./60.)*forcetheta(l) - (25./60.)*forcethetaold(l) )
            phi(l)           = phi(l) + dt*( (45./60.)*forcephi(l)-(25./60.)*forcephiold(l) )
            xfp(l)           = radfp*sin(theta(l))*cos(phi(l))
            yfp(l)           = radfp*sin(theta(l))*sin(phi(l))
            zfp(l)           = radfp*cos(theta(l))
          enddo
!$omp end do
!$omp end parallel
        endif
!
        if (mod(j,10) .eq.0) then
          forcephimean = 0.
          forcethetamean = 0.
!$omp parallel default(shared)
!$omp& private(l)
!$omp& reduction(+:forcephimean,forcethetamean,forcephivar,forcethetavar)
!$omp do
          do l=1,NL
            forcephimean = forcephimean + forcephi(l)
            forcethetamean = forcethetamean + forcetheta(l)
          enddo
!$omp end do
!$omp end parallel
          forcephimean = forcephimean/(1.*NL)
          forcethetamean = forcethetamean/(1.*NL)
          forcephivar = 0.
          forcethetavar = 0.
!$omp parallel default(shared)
!$omp& private(l)
!$omp& reduction(+:forcephivar,forcethetavar)
!$omp do
          do l=1,NL
            forcephivar = forcephivar + (forcephi(l)-forcephimean)**2.
            forcethetavar = forcethetavar + (forcetheta(l)-forcethetamean)**2.
          enddo
!$omp end do
!$omp end parallel
          forcephivar = forcephivar/(1.*NL-1.)
          forcethetavar = forcethetavar/(1.*NL-1.)
          dist1max = 9999.
          do l=2,NL
            dist = sqrt( (xfp(l)-xfp(1))**2. + (yfp(l)-yfp(1))**2. + (zfp(l)-zfp(1))**2. )
            if (dist .lt. dist1max) then
              lmin1    = l
              dist1max = dist
            endif
          enddo
          dist2max = 9999.
          do l=2,NL
            dist = sqrt( (xfp(l)-xfp(1))**2. + (yfp(l)-yfp(1))**2. + (zfp(l)-zfp(1))**2. )
            if ( (dist .lt. dist2max) .and. (l .ne. lmin1) ) then
              lmin2    = l
              dist2max = dist
            endif
          enddo
          open(22,file=datadir//'variances',position='append')
          thetanew = theta(NL)-theta(1) !position of l=NL wrt l=1
          phinew   = phi(NL)-phi(1)     !position of l=NL wrt l=1
          write(22,'(I5,7E16.8,I4,E16.8,I4,E16.8)') j,forcephimean,forcephivar,forcethetamean,forcethetavar,
     $                                              radfp*sin(thetanew)*cos(phinew),
     $                                              radfp*sin(thetanew)*sin(phinew),
     $                                              radfp*cos(thetanew),lmin1,dist1max,lmin2,dist2max
          close(22)
        endif
!
        if (mod(j,50).eq.0) then
          call random_number(rn)
          open(42,file=datadir//'lagrangianforcepoints_shell')
          write(42,*) 'VARIABLES = "lfpx","lfpy","lfpz","theta","phi","radfp"'
          write(42,*) 'ZONE T="Zone1"',' I=',NL,', F=POINT'
          write(42,*) ''
          do l=1,NL
            write(42,'(6E16.8)') xfp(l),yfp(l),zfp(l),theta(l),phi(l),radfp
          enddo
          close(42)
          call loadd(1,j,phi,theta)
        endif
      enddo

      end program init_force_points
