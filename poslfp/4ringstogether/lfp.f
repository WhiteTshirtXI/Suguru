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
!     real forcerad(1:NLtot)
      real forcephi(1:NLtot),forcetheta(1:NLtot)
      real forcephiold(1:NLtot),forcethetaold(1:NLtot)
      real difvectorx,difvectory,difvectorz
      real lengthdifvector
      real ethetax,ethetay,ethetaz
      real ephix,ephiy,ephiz
      real erx,ery,erz
      real forcephimean,forcephivar
      real forcethetamean,forcethetavar
      real xfp(1:NLtot),yfp(1:NLtot),zfp(1:NLtot)
      integer begin
      real dt
      parameter(dt=1.e-3) !1.e-5
      real dist1max,dist2max,dist
      integer lmin1,lmin2
      real phi(1:NLtot),theta(1:NLtot)
      integer npoints
      integer restart
      parameter(restart=1)
      real radfp,radfp2,radfp3,radfp4
      parameter(radfp=radius-retraction,radfp2=radfp-1./dxi,radfp3=radfp2-1./dxi,radfp4=radfp3-1./dxi)
      integer qmax
      parameter(qmax=50)
      real distr(1:qmax),lowbound,upbound,sumq
      real sumsquares
      integer rkiter
!
!     Initialization     
!
      write(6,*) 'NL,NL2,NL3,NL4,NLtot = ',NL,NL2,NL3,NL4,NLtot
      write(6,*) 'Volume occupied by 1 lfp in outer shell = ', 
     $ (4./3.)*picon*( ((radfp +0.5/dxi)**3.) - (radfp -0.5/dxi)**3. )/(1.*NL)
      write(6,*) 'Volume occupied by 1 lfp in 2nd shell = ', 
     $ (4./3.)*picon*( ((radfp2+0.5/dxi)**3.) - (radfp2-0.5/dxi)**3. )/(1.*NL2)
      write(6,*) 'Volume occupied by 1 lfp in 3rd shell = ', 
     $ (4./3.)*picon*( ((radfp3+0.5/dxi)**3.) - (radfp3-0.5/dxi)**3. )/(1.*NL3)
      write(6,*) 'Volume occupied by 1 lfp in 4rd shell = ', 
     $ (4./3.)*picon*( ((radfp4+0.5/dxi)**3.) - (radfp4-0.5/dxi)**3. )/(1.*NL4)
      write(6,*) 'Volume of Eulerian grid cell = ', 1./dxi/dyi/dzi
      do q=1,qmax
        distr(q) = 0.
      enddo
      sumq = 0
      rn = 10.
      do l=1,NLtot
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
      enddo
      do l=1,NL 
        xfp(l)      = radfp*sin(theta(l))*cos(phi(l))
        yfp(l)      = radfp*sin(theta(l))*sin(phi(l))
        zfp(l)      = radfp*cos(theta(l))
      enddo
      do q=1,NL2
        l=NL+q
        xfp(l)      = radfp2*sin(theta(l))*cos(phi(l))
        yfp(l)      = radfp2*sin(theta(l))*sin(phi(l))
        zfp(l)      = radfp2*cos(theta(l))
      enddo
      do q=1,NL3
        l=NL+NL2+q
        xfp(l)      = radfp2*sin(theta(l))*cos(phi(l))
        yfp(l)      = radfp2*sin(theta(l))*sin(phi(l))
        zfp(l)      = radfp2*cos(theta(l))
      enddo
      do q=1,NL4
        l=NL+NL2+NL3+q
        xfp(l)      = radfp3*sin(theta(l))*cos(phi(l))
        yfp(l)      = radfp3*sin(theta(l))*sin(phi(l))
        zfp(l)      = radfp3*cos(theta(l))
      enddo
      do l=1,NLtot
        forcetheta(l)    = 0.
        forcephi(l)      = 0.
      enddo
      if (sumq .ne. 2*NLtot) then
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
      write(42,*) 'ZONE T="Zone1"',' I=',NLtot,', F=POINT'
      write(42,*) ''
      do l=1,NL
        write(42,'(6E16.8)') xfp(l),yfp(l),zfp(l),theta(l),phi(l),radfp
      enddo
      do q=1,NL2
        l = q+NL
        write(42,'(6E16.8)') xfp(l),yfp(l),zfp(l),theta(l),phi(l),radfp2
      enddo
      do q=1,NL3
        l = q+NL+NL2
        write(42,'(6E16.8)') xfp(l),yfp(l),zfp(l),theta(l),phi(l),radfp3
      enddo
      do q=1,NL4
        l = q+NL+NL2+NL3
        write(42,'(6E16.8)') xfp(l),yfp(l),zfp(l),theta(l),phi(l),radfp4
      enddo
      close(42)
!     data for sphere
      npoints = nint(3.*radfp*dxi)
      open(42,file=datadir//'datasphere')
      write(42,*) 'VARIABLES = "x","y","z","r1","r2","r3","r4"'
      write(42,*) 'ZONE T="Zone1"',' I=',npoints,' J=',npoints,' K=',npoints,', F=POINT'
      write(42,*) ''
      do k=1,npoints
        do j=1,npoints
          do i=1,npoints
            write(42,'(6E16.8)') (i-npoints/2)/dxi,(j-npoints/2)/dyi,(k-npoints/2)/dzi,
     $ (sqrt( ((i-npoints/2)/dxi)**2. + ((j-npoints/2)/dyi)**2. + ((k-npoints/2)/dzi)**2. ))/radfp,
     $ (sqrt( ((i-npoints/2)/dxi)**2. + ((j-npoints/2)/dyi)**2. + ((k-npoints/2)/dzi)**2. ))/radfp2,
     $ (sqrt( ((i-npoints/2)/dxi)**2. + ((j-npoints/2)/dyi)**2. + ((k-npoints/2)/dzi)**2. ))/radfp3,
     $ (sqrt( ((i-npoints/2)/dxi)**2. + ((j-npoints/2)/dyi)**2. + ((k-npoints/2)/dzi)**2. ))/radfp4
          enddo
        enddo
      enddo
      close(42)
!
!     Compute 'Coulomb force' acting on every charged particle '
!
      begin = 0
      if (restart .eq. 1) then 
        call loadd(0,begin,phi,theta)
        do l=1,NL
          xfp(l)      = radfp*sin(theta(l))*cos(phi(l))
          yfp(l)      = radfp*sin(theta(l))*sin(phi(l))
          zfp(l)      = radfp*cos(theta(l))
        enddo
        do q=1,NL2
          l=NL+q
          xfp(l)      = radfp2*sin(theta(l))*cos(phi(l))
          yfp(l)      = radfp2*sin(theta(l))*sin(phi(l))
          zfp(l)      = radfp2*cos(theta(l))
        enddo
        do q=1,NL3
          l=NL+NL2+q
          xfp(l)      = radfp3*sin(theta(l))*cos(phi(l))
          yfp(l)      = radfp3*sin(theta(l))*sin(phi(l))
          zfp(l)      = radfp3*cos(theta(l))
        enddo
        do q=1,NL4
          l=NL+NL2+NL3+q
          xfp(l)      = radfp4*sin(theta(l))*cos(phi(l))
          yfp(l)      = radfp4*sin(theta(l))*sin(phi(l))
          zfp(l)      = radfp4*cos(theta(l))
        enddo
      endif
      write(6,*) 'begin = ',begin

      do j=begin+1,120000
        write(6,*) 'Iteration step = ',j
        rkiter = 0
 999    rkiter = rkiter+1
        do l=1,NLtot
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
          do q=1,NLtot
            if (q .ne. l) then
              difvectorx = xfp(l)-xfp(q) 
              difvectory = yfp(l)-yfp(q)
              difvectorz = zfp(l)-zfp(q)
              absforce   = 1./( difvectorx**2. + difvectory**2. + difvectorz**2. ) !1/(distance**2)
!             normalized difference vector
              lengthdifvector = sqrt( difvectorx**2. + difvectory**2. + difvectorz**2. )
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
          forcetheta(l) = forcetheta(l)/(1.*NLtot-1.) !averaged force
          forcephi(l)   = forcephi(l)/(1.*NLtot-1.)   !averaged force
!         forcerad(l)   = forcerad(l)/(1.*NLtot-1.)   !averaged force
!         The forces are expected to scale with dxi**2. Since NLtot is proportional to dxi**2,
!         they scale approximately with NLtot. By dividing the forces by NLtot, the forces become independent
!         of the resolution. This then implies that the time step is insensitive to the resolution.
        enddo
!       RK3 scheme
        if (rkiter .eq. 1) then
          do l=1,NLtot
            theta(l)    = theta(l) + dt*(32./60.)*forcetheta(l)
            phi(l)      = phi(l) + dt*(32./60.)*forcephi(l)
            forcethetaold(l) = forcetheta(l)
            forcephiold(l)   = forcephi(l)
          enddo
          do l=1,NL
            xfp(l)      = radfp*sin(theta(l))*cos(phi(l))
            yfp(l)      = radfp*sin(theta(l))*sin(phi(l))
            zfp(l)      = radfp*cos(theta(l))
          enddo
          do q=1,NL2
            l=NL+q
            xfp(l)      = radfp2*sin(theta(l))*cos(phi(l))
            yfp(l)      = radfp2*sin(theta(l))*sin(phi(l))
            zfp(l)      = radfp2*cos(theta(l))
          enddo
          do q=1,NL3
            l=NL+NL2+q
            xfp(l)      = radfp3*sin(theta(l))*cos(phi(l))
            yfp(l)      = radfp3*sin(theta(l))*sin(phi(l))
            zfp(l)      = radfp3*cos(theta(l))
          enddo
          do q=1,NL4
            l=NL+NL2+NL3+q
            xfp(l)      = radfp4*sin(theta(l))*cos(phi(l))
            yfp(l)      = radfp4*sin(theta(l))*sin(phi(l))
            zfp(l)      = radfp4*cos(theta(l))
          enddo
          go to 999
        endif
        if (rkiter .eq. 2) then
          do l=1,NLtot
            theta(l)    = theta(l) + dt*( (25./60.)*forcetheta(l) - (17./60.)*forcethetaold(l) )
            phi(l)      = phi(l) + dt*( (25./60.)*forcephi(l)-(17./60.)*forcephiold(l) )
            forcethetaold(l) = forcetheta(l)
            forcephiold(l)   = forcephi(l)
          enddo
          do l=1,NL
            xfp(l)      = radfp*sin(theta(l))*cos(phi(l))
            yfp(l)      = radfp*sin(theta(l))*sin(phi(l))
            zfp(l)      = radfp*cos(theta(l))
          enddo
          do q=1,NL2
            l=NL+q
            xfp(l)      = radfp2*sin(theta(l))*cos(phi(l))
            yfp(l)      = radfp2*sin(theta(l))*sin(phi(l))
            zfp(l)      = radfp2*cos(theta(l))
          enddo
          do q=1,NL3
            l=NL+NL2+q
            xfp(l)      = radfp3*sin(theta(l))*cos(phi(l))
            yfp(l)      = radfp3*sin(theta(l))*sin(phi(l))
            zfp(l)      = radfp3*cos(theta(l))
          enddo
          do q=1,NL4
            l=NL+NL2+NL3+q
            xfp(l)      = radfp4*sin(theta(l))*cos(phi(l))
            yfp(l)      = radfp4*sin(theta(l))*sin(phi(l))
            zfp(l)      = radfp4*cos(theta(l))
          enddo
          go to 999
        endif
        if (rkiter .eq. 3) then
          do l=1,NLtot
            theta(l)    = theta(l) + dt*( (45./60.)*forcetheta(l) - (25./60.)*forcethetaold(l) )
            phi(l)      = phi(l) + dt*( (45./60.)*forcephi(l)-(25./60.)*forcephiold(l) )
          enddo
          do l=1,NL
            xfp(l)      = radfp*sin(theta(l))*cos(phi(l))
            yfp(l)      = radfp*sin(theta(l))*sin(phi(l))
            zfp(l)      = radfp*cos(theta(l))
          enddo
          do q=1,NL2
            l=NL+q
            xfp(l)      = radfp2*sin(theta(l))*cos(phi(l))
            yfp(l)      = radfp2*sin(theta(l))*sin(phi(l))
            zfp(l)      = radfp2*cos(theta(l))
          enddo
          do q=1,NL3
            l=NL+NL2+q
            xfp(l)      = radfp3*sin(theta(l))*cos(phi(l))
            yfp(l)      = radfp3*sin(theta(l))*sin(phi(l))
            zfp(l)      = radfp3*cos(theta(l))
          enddo
          do q=1,NL4
            l=NL+NL2+NL3+q
            xfp(l)      = radfp4*sin(theta(l))*cos(phi(l))
            yfp(l)      = radfp4*sin(theta(l))*sin(phi(l))
            zfp(l)      = radfp4*cos(theta(l))
          enddo
        endif
!
        if (mod(j,10) .eq.0) then
          forcephimean = 0.
          forcethetamean = 0.
          do l=1,NLtot
            forcephimean = forcephimean + forcephi(l)
            forcethetamean = forcethetamean + forcetheta(l)
          enddo
          forcephimean = forcephimean/(1.*NLtot)
          forcethetamean = forcethetamean/(1.*NLtot)
          forcephivar = 0.
          forcethetavar = 0.
          do l=1,NLtot
            forcephivar = forcephivar + (forcephi(l)-forcephimean)**2.
            forcethetavar = forcethetavar + (forcetheta(l)-forcethetamean)**2.
          enddo
          forcephivar = forcephivar/(1.*NLtot-1.)
          forcethetavar = forcethetavar/(1.*NLtot-1.)
          dist1max = 9999.
          do l=1,NLtot-1
            dist = sqrt( (xfp(l)-xfp(NLtot))**2. + (yfp(l)-yfp(NLtot))**2. + (zfp(l)-zfp(NLtot))**2. )
            if (dist .lt. dist1max) then
              lmin1    = l
              dist1max = dist
            endif
          enddo
          dist2max = 9999.
          do l=1,NLtot-1
            dist = sqrt( (xfp(l)-xfp(NLtot))**2. + (yfp(l)-yfp(NLtot))**2. + (zfp(l)-zfp(NLtot))**2. )
            if ( (dist .lt. dist2max) .and. (l .ne. lmin1) ) then
              lmin2    = l
              dist2max = dist
            endif
          enddo
          open(22,file=datadir//'variances',position='append')
          thetanew = theta(NLtot)-theta(1) !position of l=NLtot wrt l=1
          phinew   = phi(NLtot)-phi(1)     !position of l=NLtot wrt l=1
          write(22,'(I5,7E16.8,I4,E16.8,I4,E16.8)') j,forcephimean,forcephivar,forcethetamean,forcethetavar,
     $                                              radfp4*sin(thetanew)*cos(phinew),
     $                                              radfp4*sin(thetanew)*sin(phinew),
     $                                              radfp4*cos(thetanew),lmin1,dist1max,lmin2,dist2max
          close(22)
        endif
!
        if (mod(j,50).eq.0) then
          call loadd(1,j,phi,theta)
          open(42,file=datadir//'lagrangianforcepoints')
          write(42,*) 'VARIABLES = "lfpx","lfpy","lfpz","theta","phi","radfp"'
          write(42,*) 'ZONE T="Zone1"',' I=',NLtot,', F=POINT'
          write(42,*) ''
          do l=1,NL
            write(42,'(6E16.8)') xfp(l),yfp(l),zfp(l),theta(l),phi(l),radfp
          enddo
          do q=1,NL2
            l = q+NL
            write(42,'(6E16.8)') xfp(l),yfp(l),zfp(l),theta(l),phi(l),radfp2
          enddo
          do q=1,NL3
            l = q+NL+NL2
            write(42,'(6E16.8)') xfp(l),yfp(l),zfp(l),theta(l),phi(l),radfp3
          enddo
          do q=1,NL4
            l = q+NL+NL2+NL3
            write(42,'(6E16.8)') xfp(l),yfp(l),zfp(l),theta(l),phi(l),radfp4
          enddo
          close(42)
        endif
      enddo

      end program init_force_points
