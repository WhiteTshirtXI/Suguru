     module mod_sph2
     use mod_common
     use mod_param
     use mod_math
     public 
     contains



      subroutine cal_sph_3dfun_point(x,&   
                                    theta,phi,&
                                    a,      b,&
                                    nlat,nlon)
      !calculate vat1 and vat2 for a particular pair of value (theta,phi)
      !x is \vec{x} \cdot \vec{n},used for volume calculation.
      real x(3)
      real theta,phi
      real a(nlat,nlat,3),b(nlat,nlat,3)
      real cost,sint
      real cosmphi,sinmphi
      real pbar_array (nlat+1)
      real tmp1(3)
      real mphi
      x    = 0
      lmax = nlat
      cost = cos(theta)
      sint = sin(theta)
      do m = 0,nlat-1
      call cal_legendre_array(pbar_array,lmax,m,cost)

      if(m.eq.0) then
      pbar_array = pbar_array*0.5
      end if

      mphi = m*phi
      cosmphi = cos(mphi)
      sinmphi = sin(mphi)
         
      do n = m,nlat-1
      ind  = n - m + 1
      tmp1 = (a(m+1,n+1,:)*cosmphi-b(m+1,n+1,:)*sinmphi)
      x    = x+pbar_array(ind)*tmp1
      end do
      end do

      return
      end subroutine cal_sph_3dfun_point

      
      subroutine cal_sph_3dfun_vector(x,&
                                     theta,phi,&
                                     a,      b,&
                                     pbar,&
                                     nx,nlat,nlon)
      !calculate vat1 and vat2 for a particular pair of value (theta,phi)
      !x is \vec{x} \cdot \vec{n},used for volume calculation.
      real x(nx,3)
      real theta(nx),phi(nx)
      real a(nlat,nlat,3),b(nlat,nlat,3)
      real cosmphi(nx),sinmphi(nx)
      real pbar(nx,nlat+1,0:nlat-1)
      real tmp1(nx,3)
      real mphi(nx)
      x    = 0
      lmax = nlat

      tmp1 = 0
      do m = 0,nlat-1
      mphi = m*phi
      cosmphi = cos(mphi)
      sinmphi = sin(mphi)
         
      do n = m,nlat-1
      ind  = n - m + 1
      do k=1,3
         tmp1(:,k) = a(m+1,n+1,k)*cosmphi-b(m+1,n+1,k)*sinmphi
         x(:,k) = x(:,k) + pbar(:,ind,m)*tmp1(:,k)
      end do
      end do
      end do
      
      return
      end subroutine cal_sph_3dfun_vector


      subroutine cal_sph_3dfun_withderi1_vector(x,dxdthe,dxdphi,&   
                                     theta,phi,&
                                     a,      b,&
                                     pbar,pbarderi1,&
                                     nx,nlat,nlon)
      !calculate vat1 and vat2 for a particular pair of value (theta,phi)
      !x is \vec{x} \cdot \vec{n},used for volume calculation.
      real x(nx,3),dxdthe(nx,3),dxdphi(nx,3)
      real theta(nx),phi(nx)
      real a(nlat,nlat,3),b(nlat,nlat,3)
      real cost(nx),sint(nx)
      real cosmphi(nx),sinmphi(nx)
      real pbar(nx,nlat+1,0:nlat-1)
      real pbarderi1(nx,nlat+1,0:nlat-1)
      real tmp1(nx,3)
      real mphi(nx)
      x    = 0
      dxdthe = 0
      dxdphi = 0
      lmax = nlat
      sint = sin(theta)

      do m = 0,nlat-1
      mphi = m*phi
      cosmphi = cos(mphi)
      sinmphi = sin(mphi)
         
      do n = m,nlat-1
      ind  = n - m + 1
      do k=1,3
         tmp1(:,k) = a(m+1,n+1,k)*cosmphi-b(m+1,n+1,k)*sinmphi
         x(:,k) = x(:,k) + pbar(:,ind,m)*tmp1(:,k)
         dxdthe(:,k) = dxdthe(:,k) - sint*pbarderi1(:,ind,m)*tmp1(:,k)
         tmp1(:,k) = m*(-a(m+1,n+1,k)*sinmphi-b(m+1,n+1,k)*cosmphi)
         dxdphi(:,k) = dxdphi(:,k) + pbar(:,ind,m)*tmp1(:,k)
      end do
      end do
      end do
      
      return
      end subroutine cal_sph_3dfun_withderi1_vector




      subroutine cal_sph_3dfun_mesh2d(fun3d,the2d,phi2d,&
                               a,b,nthe,nphi,nlat,nlon)
      real fun3d(3,nthe,nphi)
      real the2d(nthe,nphi),phi2d(nthe,nphi)
      real a(nlat,nlat,3),b(nlat,nlat,3)
      do j=1,nphi
      do i=1,nthe
      call cal_sph_3dfun_point(fun3d(:,i,j),&   
                     the2d(i,j),phi2d(i,j),&
                                 a,      b,&
                                 nlat,nlon)
      end do
      end do
      end subroutine cal_sph_3dfun_mesh2d


      subroutine create_2dmesh_sphere(the1d,phi1d,the2d,phi2d,nthe,nphi)
      !!!!!!!nphi=2*(nthe-1)
      common /myconstants/ pimy
      real the1d(nthe),phi1d(nphi)
      real the2d(nthe,nphi),phi2d(nthe,nphi)
      dphi = 2*pimy/(nphi-1)
      dthe =   pimy/(nthe-1)
      forall(j=1:nphi) phi1d(j)=(j-1)*dphi
      forall(j=1:nthe) the1d(j)=(j-1)*dthe
      forall(j=1:nphi) the2d(:,j)=the1d
      forall(j=1:nthe) phi2d(j,:)=phi1d
      end subroutine create_2dmesh_sphere


      end module mod_sph2

