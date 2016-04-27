     module mod_sph3
     use mod_common
     use mod_param
     use mod_math
     public 
     contains


      subroutine cal_co_tang_basis2(vat1,vat2,&
     deFvat1,deFvat2,&
     deSvat,deTvat,&
     a,b,dtheta,dphi,&
     nlat2,nlon2,&
     nlat,nlon,&
     ibend)
!calculate covariant tangent basis vector 
!      include 'CELL3DGEO'
      integer,save:: icald
      real a   (nlat2,nlat2,3)
      real b   (nlat2,nlat2,3)
      real, allocatable, save,dimension(:,:,:) :: pbar,pbar_deri1,&
     pbar_deri2, &
     pbar_deri3,pbar_deri4
      real, allocatable, save,dimension(:,:) :: cosmphi,sinmphi
      real, allocatable, save,dimension(:) :: cost,sint
      real vat1(3,nlat,nlon)
      real vat2(3,nlat,nlon)
      real deFvat1(3,2,nlat,nlon)! F: first derivative
      real deFvat2(3,2,nlat,nlon) 
      real deSvat(3,4,nlat,nlon) ! S: second derivative
      real deTvat(3,5,nlat,nlon) ! T: third derivative
      integer lmax
      real mphi
      real dtheta(nlat),dphi(nlon)
      real parthe0(3),parthe1(3),parthe2(3),parthe3(3),parthe4(3)
      data icald /0/

      vat1 = 0
      vat2 = 0
      deFvat1 = 0
      deFvat2 = 0
      deSvat = 0
      deTvat = 0
      lmax = nlat2

      if(icald.eq.0) then
         allocate(pbar(nlat2+1,nlat2,nlat))
         allocate(pbar_deri1(nlat2+1,nlat2,nlat))
         allocate(pbar_deri2(nlat2+1,nlat2,nlat))
         allocate(pbar_deri3(nlat2+1,nlat2,nlat))
         allocate(pbar_deri4(nlat2+1,nlat2,nlat))
         allocate(cosmphi(nlat2,nlon),sinmphi(nlat2,nlon))
         allocate(cost(nlat),sint(nlat))
         cost = cos(dtheta(1:nlat))
         sint = sin(dtheta(1:nlat))
      if(ibend.ne.1) then
         do i = 1,nlat
            do m=0,nlat2-1
!!!!!!!!!!!!!!NO BENDING MOMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call cal_legendre_deri2_array(pbar(1,m+1,i), &
     pbar_deri1(1,m+1,i),pbar_deri2(1,m+1,i),lmax,m,cost(i))

            pbar_deri2(:,m+1,i) = pbar_deri2(:,m+1,i)*(sint(i)**2) - &
                                cost(i)*pbar_deri1(:,m+1,i)
            pbar_deri1(:,m+1,i) = pbar_deri1(:,m+1,i)*(-sint(i))
!!!!!!!!!!!!!!NO BENDING MOMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end do
         end do
         pbar(:,1,:)       = pbar(:,1,:)*0.5
         pbar_deri1(:,1,:) = pbar_deri1(:,1,:)*0.5  
         pbar_deri2(:,1,:) = pbar_deri2(:,1,:)*0.5               
      else
         do i = 1,nlat
            do m=0,nlat2-1
!!!!!!!!!!!!!!WITH BENDING MOMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call cal_legendre_deri4_array(pbar(1,m+1,i),&
     pbar_deri1(1,m+1,i),pbar_deri2(1,m+1,i),pbar_deri3(1,m+1,i),&
     pbar_deri4(1,m+1,i),lmax,m,cost(i))
            pbar_deri4(:,m+1,i) = pbar_deri4(:,m+1,i)*sint(i)**4-&
     6*pbar_deri3(:,m+1,i)*sint(i)**2*cost(i)+pbar_deri2(:,m+1,i)*&
     (-4*sint(i)**2+3*cost(i)**2)+pbar_deri1(:,m+1,i)*cost(i)
            pbar_deri3(:,m+1,i) = -pbar_deri3(:,m+1,i)*sint(i)**3+&
     3*sint(i)*cost(i)*pbar_deri2(:,m+1,i)+sint(i)*pbar_deri1(:,m+1,i)
            pbar_deri2(:,m+1,i) = pbar_deri2(:,m+1,i)*(sint(i)**2)-&
     cost(i)*pbar_deri1(:,m+1,i)
            pbar_deri1(:,m+1,i) = pbar_deri1(:,m+1,i)*(-sint(i))
!!!!!!!!!!!!!!WITH BENDING MOMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end do
         end do
         pbar(:,1,:)       = pbar(:,1,:)*0.5
         pbar_deri1(:,1,:) = pbar_deri1(:,1,:)*0.5  
         pbar_deri2(:,1,:) = pbar_deri2(:,1,:)*0.5               
         pbar_deri3(:,1,:) = pbar_deri3(:,1,:)*0.5  
         pbar_deri4(:,1,:) = pbar_deri4(:,1,:)*0.5               
      end if
      
      do j=1,nlon
         do m=0,nlat2-1
            mphi=m*dphi(j)
            cosmphi(m+1,j)=cos(mphi)
            sinmphi(m+1,j)=sin(mphi)            
         end do
      end do
      icald = icald + 1
      end if





      if(ibend.ne.1) then
      do j = 1,nlon
         do i = 1,nlat
            do m=0,nlat2-1
            do n=m,nlat2-1
            ind = n - m + 1
            parthe0 = (a(m+1,n+1,:)*cosmphi(m+1,j)-b(m+1,n+1,:)*&
     sinmphi(m+1,j))
            parthe1 = m*(a(m+1,n+1,:)*(-sinmphi(m+1,j))-b(m+1,n+1,:)*&
     cosmphi(m+1,j))
            parthe2 = -m**2*parthe0
            vat1(:,i,j)=vat1(:,i,j)+pbar_deri1(ind,m+1,i)*parthe0
            vat2(:,i,j)=vat2(:,i,j)+pbar(ind,m+1,i)*parthe1
            deFvat1(:,1,i,j)=deFvat1(:,1,i,j)+&
     pbar_deri2(ind,m+1,i)*parthe0
            deFvat1(:,2,i,j)=deFvat1(:,2,i,j)+&
     pbar_deri1(ind,m+1,i)*parthe1
            deFvat2(:,2,i,j)=deFvat2(:,2,i,j)+pbar(ind,m+1,i)*parthe2
            !!!!!!!!!!!get vat1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end do
            end do
         end do
      end do
      deFvat2(:,1,:,:) = deFvat1(:,2,:,:)

!!!!!!!!!!!!!!NO BENDING MOMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else
!!!!!!!!!!!!!!WITH BENDING MOMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j = 1,nlon
         do i = 1,nlat
            do m=0,nlat2-1
            do n=m,nlat2-1
            ind = n - m + 1
            !!!!!!!!!!!get vat1=\frac{\partial x}{\partial theta}!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!get vat2=\frac{\partial x}{\partial phi}!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !!!!!!!!!!!get deFvat1(:,1,i,j)=\frac{\partial vat1}{\partial theta} !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!get deFvat1(:,2,i,j)=\frac{\partial vat1}{\partial phi} !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
            !!!!!!!!!!!get deFvat2(:,2,i,j)=\frac{\partial vat2}{\partial phi} !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            parthe0 = (a(m+1,n+1,:)*cosmphi(m+1,j)-b(m+1,n+1,:)*&
     sinmphi(m+1,j))
            parthe1 = m*(a(m+1,n+1,:)*(-sinmphi(m+1,j))-b(m+1,n+1,:)*&
     cosmphi(m+1,j))
            parthe2 = -m**2*parthe0
            parthe3 = m**3*(a(m+1,n+1,:)*sinmphi(m+1,j)+b(m+1,n+1,:)*&
     cosmphi(m+1,j))
            parthe4 = m**4*parthe0
            vat1(:,i,j)=vat1(:,i,j)+pbar_deri1(ind,m+1,i)*parthe0
            vat2(:,i,j)=vat2(:,i,j)+pbar(ind,m+1,i)*parthe1
            deFvat1(:,1,i,j)=deFvat1(:,1,i,j)+&
     pbar_deri2(ind,m+1,i)*parthe0
            deFvat1(:,2,i,j)=deFvat1(:,2,i,j)+&
     pbar_deri1(ind,m+1,i)*parthe1
            deFvat2(:,2,i,j)=deFvat2(:,2,i,j)+pbar(ind,m+1,i)*parthe2
            !!!!!!!!!!!get vat1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           deSvat(:,1,i,j)=deSvat(:,1,i,j)+pbar_deri3(ind,m+1,i)*parthe0 
           deSvat(:,2,i,j)=deSvat(:,2,i,j)+pbar_deri2(ind,m+1,i)*parthe1
           deSvat(:,3,i,j)=deSvat(:,3,i,j)+pbar_deri1(ind,m+1,i)*parthe2 
           deSvat(:,4,i,j)=deSvat(:,4,i,j)+pbar(ind,m+1,i)*parthe3
           deTvat(:,1,i,j)=deTvat(:,1,i,j)+pbar_deri4(ind,m+1,i)*parthe0
           deTvat(:,2,i,j)=deTvat(:,2,i,j)+pbar_deri3(ind,m+1,i)*parthe1            
           deTvat(:,3,i,j)=deTvat(:,3,i,j)+pbar_deri2(ind,m+1,i)*parthe2
           deTvat(:,4,i,j)=deTvat(:,4,i,j)+pbar_deri1(ind,m+1,i)*parthe3                   
           deTvat(:,5,i,j)=deTvat(:,5,i,j)+pbar(ind,m+1,i)*parthe4
            end do
            end do
         end do
      end do
      deFvat2(:,1,:,:) = deFvat1(:,2,:,:)
      end if
      return
      end subroutine cal_co_tang_basis2


!*******************************************************calculates surface metric***********************************************************
      subroutine cal_surface_metric_sph__OPT(surfmet,&
     a,b,dtheta,dphi,&
     nlat,nlon,&
     cvat1x,cvat1y,cvat1z,&
     cvat2x,cvat2y,cvat2z)

!calculate covariant tangent basis vector 
!      include 'CELL3DGEO'
      real a   (nlat,nlat,3)
      real b   (nlat,nlat,3)
!      real pbar (nlat+1)
!      real pbar_deri1(nlat+1)
      real surfmet(nlat,nlon)
      real vat1(3,nlat,nlon)
      real vat2(3,nlat,nlon)
      integer lmax
      real mphi
      real dtheta(nlat),dphi(nlon)
      real parthe0(3),parthe1(3)
      real cross(3)
      real,save,allocatable,dimension(:,:) :: cosmphi,sinmphi
      real,save,allocatable,dimension(:) :: cost,sint
      real,save,allocatable,dimension(:,:,:) :: pbar,pbar_deri1
      real,dimension(:,:),optional :: cvat1x,cvat1y,cvat1z,&
                                     cvat2x,cvat2y,cvat2z
      integer, save :: icald
      data icald /0/
      vat1 = 0
      vat2 = 0
      lmax = nlat

      if(icald.eq.0) then
         allocate(cosmphi(nlat,nlon),sinmphi(nlat,nlon))
         allocate(cost(nlat),sint(nlat))
         allocate(pbar(1+nlat,nlat,nlat),pbar_deri1(1+nlat,nlat,nlat))
         do j=1,nlon
            do m=0,nlat-1
               mphi=m*dphi(j)
               cosmphi(m+1,j)=cos(mphi)
               sinmphi(m+1,j)=sin(mphi)            
            end do
         end do
         cost=cos(dtheta(1:nlat))
         sint=sin(dtheta(1:nlat))

         do i=1,nlat
            do m=0,nlat-1
            call cal_legendre_deri1_array(pbar(1,m+1,i),&
     pbar_deri1(1,m+1,i),lmax,m,cost(i))
            pbar_deri1(:,m+1,i) = pbar_deri1(:,m+1,i)*(-sint(i))
            end do
         end do
         pbar(:,1,:)       = pbar(:,1,:)*0.5
         pbar_deri1(:,1,:) = pbar_deri1(:,1,:)*0.5  
         icald = icald + 1
      end if


      do j = 1,nlon
         do i = 1,nlat
            do m=0,nlat-1
            do n=m,nlat-1
            ind = n - m + 1
            parthe0 = (a(m+1,n+1,:)*cosmphi(m+1,j)-b(m+1,n+1,:)*&
     sinmphi(m+1,j))
            parthe1 = m*(a(m+1,n+1,:)*(-sinmphi(m+1,j))-b(m+1,n+1,:)*&
     cosmphi(m+1,j))
            vat1(:,i,j)=vat1(:,i,j)+pbar_deri1(ind,m+1,i)*parthe0
            vat2(:,i,j)=vat2(:,i,j)+pbar(ind,m+1,i)*parthe1
            end do
            end do
         end do
      end do

      do j=1,nlon
         do i=1,nlat
!**********************************************************************subroutine vcross is in the math module********************************************************************************
            call vcross(cross(1),cross(2),cross(3),&
     vat1(1,i,j),vat1(2,i,j),vat1(3,i,j),&
     vat2(1,i,j),vat2(2,i,j),vat2(3,i,j),&
     1)
!**********************************************************************function rnorm2 is in the math module**********************************************************************************
         surfmet(i,j)=rnorm2(cross(1),cross(2),cross(3))
         end do
      end do

      if(present(cvat1x)) then
         cvat1x=vat1(1,:,:);cvat1y=vat1(2,:,:);cvat1z=vat1(3,:,:)
         cvat2x=vat2(1,:,:);cvat2y=vat2(2,:,:);cvat2z=vat2(3,:,:)
      end if
      return
      end subroutine cal_surface_metric_sph__OPT



      subroutine cal_co_vat1(vat1,vat2,&
                            x        ,  & 
                            theta,phi,&
                            a,      b,&
                            nlat,nlon)
      !calculate vat1 and vat2 for a particular pair of value (theta,phi)
      !x is \vec{x} \cdot \vec{n},used for volume calculation.
      real vat1(3),vat2(3)
      real x(3)
      real theta,phi
      real a(nlat,nlat,3),b(nlat,nlat,3)
      real cost,sint
      real cosmphi,sinmphi
!      real pbar (nlat+1),pbar_deri1 (nlat+1)
      real tmp1(3),tmp2(3)
      real mphi
      real,allocatable,save,dimension(:,:) :: pbar,pbar_deri1
      integer,save :: icald
      data icald /0/
      
      x    = 0
      vat1 = 0
      vat2 = 0
      lmax = nlat
      cost = cos(theta)
      sint = sin(theta)

      if(icald.eq.0) then
         allocate(pbar(1+nlat,nlat),pbar_deri1(1+nlat,nlat))
         icald = icald + 1
      end if

!      write(*,*) cost
      do m=0,nlat-1
         call cal_legendre_deri1_array(pbar(1,m+1),pbar_deri1(1,m+1),&
             lmax,m,cost)
         pbar_deri1(:,m+1)  = pbar_deri1(:,m+1)*(-sint)
         if(m.eq.0) then
            pbar(:,m+1)       = pbar(:,m+1)*0.5
            pbar_deri1(:,m+1) = pbar_deri1(:,m+1)*0.5
         end if
         
         mphi = m*phi
         cosmphi = cos(mphi)
         sinmphi = sin(mphi)
         do n=m,nlat-1
            ind = n - m + 1
            tmp1 = (a(m+1,n+1,:)*cosmphi-b(m+1,n+1,:)*sinmphi)
            tmp2 = (a(m+1,n+1,:)*(-sinmphi)-b(m+1,n+1,:)*cosmphi)
            vat1 = vat1+pbar_deri1(ind,m+1)*tmp1
            x    = x+pbar(ind,m+1)*tmp1
            vat2 = vat2+pbar(ind,m+1)*m*tmp2
         end do
      end do

      return
      end  subroutine cal_co_vat1


      subroutine cal_cell_area_vol_new(area,vol,a,b,wgtdsin,&
      x,y,z,the,phi,nlat,nlon,nhm)
!!!!!!!!when calculating the initial volume, we need unit normal, and no dealiasing is used for the unit normal.
      real area,vol
      real a(nlat,nlat,3),b(nlat,nlat,3)
      real x(nlat,nlon),y(nlat,nlon),z(nlat,nlon)
      real wgtdsin(nlat,nlon),surfmet(nlat,nlon),wgt(nlat,nlon)
      real the(nlat),phi(nlon)
      real vat1(nhm,3),vat2(nhm,3),van3(nhm,3),xdotn(nlat,nlon),&
     van3norm(nhm)
      call   cal_surface_metric_sph__OPT(surfmet,&
                                    a,b,the,phi,&
                                      nlat,nlon)
      wgt = surfmet*wgtdsin
      area = sum(wgt)

      call cal_co_vat1_array(vat1,vat2,&
          the, phi,&
          a, b,&
          nlat,nlon,nhm)
      call vcross(van3(1,1),van3(1,2),van3(1,3),&
          vat1(1,1),vat1(1,2),vat1(1,3),&
          vat2(1,1),vat2(1,2),vat2(1,3),&
          nhm)
      van3norm = sqrt(van3(:,1)**2+van3(:,2)**2+van3(:,3)**2)

      do i=1,3
         van3(:,i) = van3(:,i)/van3norm
      end do

      call vdot3(xdotn(1,1),x(1,1),y(1,1),z(1,1),&
          van3(1,1),van3(1,2),van3(1,3),&
          nhm) 
      vol = sum(wgt*xdotn)/3.0

      print*,"arash test new",area,vol
      end subroutine cal_cell_area_vol_new





      subroutine cal_co_vat1_array(vat1,vat2,&
                                  theta,phi,&
                                  a,      b,&
                                  nlat,nlon,nhm)
      !calculate vat1 and vat2 for a particular pair of value (theta,phi)
      !x is \vec{x} \cdot \vec{n},used for volume calculation.
      real vat1(nhm,3),vat2(nhm,3)
      real theta(nlat),phi(nlon)
      real a(nlat,nlat,3),b(nlat,nlat,3)
      real tmp1(3),tmp2(3)
      real mphi
      real,allocatable,save,dimension(:,:,:) :: pbar,pbar_deri1
      real, allocatable, save,dimension(:,:) :: cosmphi,sinmphi
      real, allocatable, save,dimension(:) :: cost,sint
      integer,save :: icald
      data icald /0/
      
      vat1 = 0
      vat2 = 0
      lmax = nlat
      cost = cos(theta)
      sint = sin(theta)

      if(icald.eq.0) then
         allocate(pbar(1+nlat,nlat,nlat),pbar_deri1(1+nlat,nlat,nlat))
         allocate(cosmphi(nlat,nlon),sinmphi(nlat,nlon))
         allocate(cost(nlat),sint(nlat))
         cost = cos(theta)
         sint = sin(theta)
         do i = 1,nlat
            do m=0,nlat-1
!!!!!!!!!!!!!!NO BENDING MOMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call cal_legendre_deri1_array(pbar(1,m+1,i),&
     pbar_deri1(1,m+1,i),lmax,m,cost(i))
            pbar_deri1(:,m+1,i) = pbar_deri1(:,m+1,i)*(-sint(i))
!!!!!!!!!!!!!!NO BENDING MOMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end do
         end do
         pbar(:,1,:)       = pbar(:,1,:)*0.5
         pbar_deri1(:,1,:) = pbar_deri1(:,1,:)*0.5  

         do j=1,nlon
            do m=0,nlat-1
               mphi=m*phi(j)
               cosmphi(m+1,j)=cos(mphi)
               sinmphi(m+1,j)=sin(mphi)            
            end do
         end do
         icald = icald + 1
      end if

      do ij=1,nhm
         call ind1to2(i,j,ij,nlat,nlon)
         do m=0,nlat-1
            do n=m,nlat-1
               ind = n - m + 1
               tmp1 = (a(m+1,n+1,:)*cosmphi(m+1,j)-&
                   b(m+1,n+1,:)*sinmphi(m+1,j))
               tmp2 = (a(m+1,n+1,:)*(-sinmphi(m+1,j))-&
                   b(m+1,n+1,:)*cosmphi(m+1,j))
               vat1(ij,:) = vat1(ij,:)+pbar_deri1(ind,m+1,i)*tmp1
               vat2(ij,:) = vat2(ij,:)+pbar(ind,m+1,i)*m*tmp2
            end do
         end do
      end do
      return
      end subroutine cal_co_vat1_array





      end module mod_sph3

