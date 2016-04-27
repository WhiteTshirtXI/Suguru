

!     a program for testing all scalar analysis and synthesis subroutines




     module mod_sph1
     use mod_common
     use mod_param
     use mod_math
     use mod_sph3
     public 
     contains

      subroutine testsphharm
      integer l,m
      real :: coeff,x
      l = 10
      m = 3
      x = 0.4
      call cal_legendre(coeff,l,m,x)
      write(*,*) "legendre coeff is:", coeff
      return
      end subroutine testsphharm



      subroutine cal_3dcrd_spectral_coeff(a,b,x,y,z,nlat,nlon)
      real :: a(nlat,nlat,3),b(nlat,nlat,3)
      real :: x(nlat,nlon),y(nlat,nlon),z(nlat,nlon)
      real :: stmp(nlat,nlon,3)
      stmp(:,:,1)=x
      stmp(:,:,2)=y
      stmp(:,:,3)=z
      call cal_spectral_coeff1(a,b,stmp,nlat,nlon,3)
      end subroutine cal_3dcrd_spectral_coeff

      subroutine cal_spectral_coeff1(a,b,s,nlat,nlon,nt)
      integer nlat,nlon,nt
!      real :: dwork(lldwork)
!      real :: work(lleng),wsave(llsav)
      real :: a(nlat,nlat,nt),b(nlat,nlat,nt),&
     s(nlat,nlon,nt)
      real :: dwts(nlat)
      real :: dtheta(nlat)
      real,allocatable :: dwork(:),work(:),wsave(:)
      lwork  =  5*nlat*nlat*nlon
      lsave  =  5*nlat*nlat*nlon
      ldwork =  nlat*(nlat+4)
      isym = 0
      a = 0
      b = 0
      allocate(dwork(lwork),work(lwork),wsave(lsave))
      dwork = 0
      work  = 0
      wsave = 0
!     
!     set equally spaced colatitude and longitude increments
!

!     compute nlat gaussian points in thetag
!

      call gaqd(nlat,dtheta,dwts,dwork,ldwork,ier)
       wsave = 0.0

!      call name("**gc")

      call shagci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)

!      call name("shai")
!      call iout(ierror,"ierr")

      call shagc(nlat,nlon,isym,nt,s,nlat,nlon,a,b,nlat,nlat,wsave,&
     +lsave,work,lwork,ierror)

      deallocate(dwork,work,wsave)
      end subroutine cal_spectral_coeff1


      subroutine cal_surface_metric_sph_Itp2(cvat1x,cvat1y,cvat1z,&
     cvat2x,cvat2y,cvat2z,&
     a,b,theoe,phioe,&
     nlat,nlon,nltoe,nlnoe)
!calculate covariant tangent basis vector 

      
      real :: a   (nlat,nlat,3)
      real :: b   (nlat,nlat,3)
      real :: vat1(3,nltoe,nlnoe)
      real :: vat2(3,nltoe,nlnoe)
      integer,save :: lmax,icald,nmodes
      real :: mphi
      real :: theoe(nltoe),phioe(nlnoe)
      real :: parthe0(3),parthe1(3)
      real :: cross(3)
      real :: cvat1x(nltoe,nlnoe),cvat1y(nltoe,nlnoe),&
     cvat1z(nltoe,nlnoe),cvat2x(nltoe,nlnoe),cvat2y(nltoe,nlnoe),&
     cvat2z(nltoe,nlnoe)
      real ,allocatable,save :: pbare(:,:,:), pbare_deri1(:,:,:),&
     cosmpe(:,:), sinmpe(:,:)
      data icald /0/
      real :: pbar (nlat+1)
      real :: pbar_deri1(nlat+1)
      lmax = nlat

      vat1 = 0; vat2 = 0
      do j = 1,nlnoe
	  phi = phioe(j)
         do i = 1,nltoe !!! when i = 1 or nltoe, sin(theta)=0, you have problem to calculate pbar_deri1,
	    theta = theoe(i)
	    cost  = cos(theta)
	    sint  = sin(theta)
            do m=0,nlat-1
!!!!!!!!!!!!!!NO BENDING MOMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(i.eq.1.or.i.eq.nltoe) then
               call cal_legendre_array(pbar,lmax,m,cost)
               call cal_deri1_array(pbar_deri1,lmax,m,cost)
            else
            call cal_legendre_deri1_array(pbar,pbar_deri1,&
           lmax,m,cost)
            pbar_deri1 = pbar_deri1*(-sint)
            end if

            if(m.eq.0) then
            pbar       = pbar*0.5
            pbar_deri1 = pbar_deri1*0.5  
            end if
!!!!!!!!!!!!!!NO BENDING MOMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            mphi = m*phi
            cosmphi = cos(mphi)
            sinmphi = sin(mphi)
               
            do n=m,nlat-1
            ind = n - m + 1
            parthe0 = (a(m+1,n+1,:)*cosmphi-b(m+1,n+1,:)*sinmphi)
            parthe1 = m*(a(m+1,n+1,:)*(-sinmphi)-b(m+1,n+1,:)*cosmphi)
            vat1(:,i,j)=vat1(:,i,j)+pbar_deri1(ind)*parthe0
            vat2(:,i,j)=vat2(:,i,j)+pbar(ind)*parthe1
            end do
            end do
         end do
      end do    
      cvat1x = vat1(1,:,:);cvat1y = vat1(2,:,:);cvat1z = vat1(3,:,:)
      cvat2x = vat2(1,:,:);cvat2y = vat2(2,:,:);cvat2z = vat2(3,:,:)
      return
      end subroutine cal_surface_metric_sph_Itp2






!***********************************************************computes covariant vectors*************************************************
      subroutine cal_co_all_basis(cvat1,cvat2,cvan3,&
     surfmet,&
     cdevat1,cdevat2,&
     cdevan3,&
     cdeSvat,cdeTvat,&
     a,b,dtheta,dphi,&
     nlat2,nlon2,&
     nlat,nlon,&
     ibend)
      real :: a   (nlat2,nlat2,3)
      real :: b   (nlat2,nlat2,3)
      real :: dtheta(nlat),dphi(nlon)
      real :: cvat1(3,nlat,nlon)
      real :: cvat2(3,nlat,nlon)
      real :: cvan3(3,nlat,nlon)
      real :: cdevat1(3,2,nlat,nlon)
      real :: cdevat2(3,2,nlat,nlon)
      real :: cdevan3(3,2,nlat,nlon)
      real :: cdeSvat(3,4,nlat,nlon) ! S: second derivative
      real :: cdeTvat(3,5,nlat,nlon) ! T: third derivative
      real :: surfmet(nlat,nlon)
!*******************************this is located in spherharm3.f cvat1 and cvat2 are coordinates of tangent vectors*************************
      call cal_co_tang_basis2(cvat1,cvat2,&
     cdevat1,cdevat2,&
     cdeSvat,cdeTvat,&
     a,b,dtheta,dphi,&
     nlat2,nlon2,&
     nlat,nlon,&
     ibend)
!************************************************cvan3 is normal and cdevan3 is its derivative*********************************************
      call cal_unit_normal_deri(cvan3,cdevan3,&
     surfmet,&
     cvat1,cvat2,&
     cdevat1,cdevat2,&
     nlat,nlon)

      return
      end subroutine cal_co_all_basis

!***********************************************************computes contravariant vectors*************************************************
      subroutine cal_contra_all_basis(convat1,convat2,convan3,&
     conderi1,conderi2,&     
     cvat1,cvat2,cvan3,&
     cderi1,cderi2,cderi3,& 
     nlat,nlon)
      use mod_capsulesolid
      real :: cvat1(3,nlat,nlon),convat1(3,nlat,nlon)
      real :: cvat2(3,nlat,nlon),convat2(3,nlat,nlon)
      real :: cvan3(3,nlat,nlon),convan3(3,nlat,nlon)
      real :: cderi1(3,2,nlat,nlon)
      real :: cderi2(3,2,nlat,nlon)
      real :: cderi3(3,2,nlat,nlon)
      real :: conderi1(3,2,nlat,nlon)
      real :: conderi2(3,2,nlat,nlon)
      real :: g,gde(2)
      real :: a1(3),a2(3),n3(3)
      real :: a1de(3,2),a2de(3,2),n3de(3,2)
      real :: a2_n3(3)     ! cvat2 \cross cvan3
      real :: n3_a1(3)     ! n \cross cvat1
      real :: n3de_a1(3,2) ! \frac{\partial cvan3}{\partial \alpha} \cross cvat1
      real :: n3_a1de(3,2) ! cvan3 \cross \frac{\partial cvat1}{\partial \alpha} 
      real :: a2de_n3(3,2) ! \frac{\partial cvat2}{\partial \alpha} \cross cvan3
      real :: a2_n3de(3,2) ! cvat2 \cross \frac{\partial cvan3}{\partial \alpha}
      integer al
      real :: stmp1,stmp2
      real :: vtmp(3)

      do j=1,nlon
      do i=1,nlat
      a1   = cvat1(:,i,j)
      a2   = cvat2(:,i,j)
      n3   = cvan3(:,i,j)
      a1de = cderi1(:,:,i,j)
      a2de = cderi2(:,:,i,j)
      n3de = cderi3(:,:,i,j)

      call co2contra(a1(1),         a1(2),         a1(3),&
                    a2(1),         a2(2),         a2(3),&
                    n3(1),         n3(2),         n3(3),&
                    convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),&
                    convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),&
                    convan3(1,i,j),convan3(2,i,j),convan3(3,i,j))

      call vcross(a2_n3(1),a2_n3(2),a2_n3(3),&
                 a2(1),   a2(2),   a2(3),&
                 n3(1),   n3(2),   n3(3),&
     1)
      call vcross(n3_a1(1),n3_a1(2),n3_a1(3),&
                 n3(1),   n3(2),   n3(3),&
                 a1(1),   a1(2),   a1(3),&
     1)

      do al = 1,2
      call vcross(a2de_n3(1,al),a2de_n3(2,al),a2de_n3(3,al),&
                 a2de(1,al),   a2de(2,al),   a2de(3,al),&
                 n3(1),        n3(2),        n3(3),&
     1)
      call vcross(a2_n3de(1,al),a2_n3de(2,al),a2_n3de(3,al),&
                 a2(1),        a2(2),        a2(3),&
                 n3de(1,al),   n3de(2,al),   n3de(3,al),&
     1)

      call vcross(n3de_a1(1,al),n3de_a1(2,al),n3de_a1(3,al),&
                 n3de(1,al),   n3de(2,al),   n3de(3,al),&
                 a1(1),        a1(2),        a1(3),&
     1)

      call vcross(n3_a1de(1,al),n3_a1de(2,al),n3_a1de(3,al),&
                 n3(1),        n3(2),        n3(3),&
                 a1de(1,al),   a1de(2,al),   a1de(3,al),&
     1)
      end do

      call vdot3(g,a1(1),a1(2),a1(3),&
     a2_n3(1),a2_n3(2),a2_n3(3),&
     1)

      do al=1,2
      vtmp = a2de_n3(:,al) + a2_n3de(:,al)
      call vdot3(stmp1,a1de(1,al),a1de(2,al),a1de(3,al),&
      a2_n3(1),a2_n3(2),a2_n3(3),&
     1)
      call vdot3(stmp2,a1(1),a1(2),a1(3),&
      vtmp(1),vtmp(2),vtmp(3),&
     1)
      gde(al) = stmp1 + stmp2
      end do
!      write(*,*) gde(1)-sph_rad**2*cos(dtheta(i)),gde(2)

      do al=1,2
      conderi1(:,al,i,j)=1.0/g*(-gde(al)*convat1(:,i,j)+ &
      a2de_n3(:,al)+a2_n3de(:,al))
      conderi2(:,al,i,j)=1.0/g*(-gde(al)*convat2(:,i,j)+& 
      n3de_a1(:,al)+n3_a1de(:,al))
      end do

      end do
      end do

      return
      end subroutine cal_contra_all_basis

      subroutine cal_contra_second_derivative(&
     condeS1,condeS2,&     
     cvat1,cvat2,cvan3,&
     cderi1,cderi2,cderi3,& 
     cdeSvat,&
     cderiS3,&
     nlat,nlon,&
     ibend_dummy)
      real :: cvat1(3,nlat,nlon)
      real :: cvat2(3,nlat,nlon)
      real :: cvan3(3,nlat,nlon)
      real :: cderi1(3,2,nlat,nlon)
      real :: cderi2(3,2,nlat,nlon)
      real :: cderi3(3,2,nlat,nlon)
      real :: cdeSvat(3,4,nlat,nlon)
      real :: cderiS3(3,3,nlat,nlon)
      real :: condeS1(3,3,nlat,nlon)
      real :: condeS2(3,3,nlat,nlon)

      real :: g,gde(2),gdeS(3),ginv1,ginv2,ginv3
      real :: a1(3),a2(3),n3(3)
      real :: a1de(3,2),a2de(3,2),n3de(3,2),adeS(3,4)
      real :: a2_n3(3)     ! cvat2 \cross cvan3
      real :: n3_a1(3)     ! n \cross cvat1
      real :: n3de_a1(3,2) ! \frac{\partial cvan3}{\partial \alpha} \cross cvat1
      real :: n3_a1de(3,2) ! cvan3 \cross \frac{\partial cvat1}{\partial \alpha} 
      real :: a2de_n3(3,2) ! \frac{\partial cvat2}{\partial \alpha} \cross cvan3
      real :: a2_n3de(3,2) ! cvat2 \cross \frac{\partial cvan3}{\partial \alpha}

      real :: a2_n3deS(3,3)! cvat2 \cross \frac{\partial^2 cvan3}{\partial \alpha \partial \beta} (:,1) to alpha=beta=1, (:,2) to alpha=beta=2, (:,3) to alpha=1(2),beta=2(1)
      real :: a2de_n3de(3,2,2) !(:,beta,alpha) = \frac{\partial cvat2}{\partial \alpha} \cross \frac{\partial cvan3}{\partial \beta}
      real :: a2deS_n3(3,3) !a2deS is \frac{\partial^2 cvat2}{\partial \alpha \partial \beta}
      real :: n3_a1deS(3,3)
      real :: n3de_a1de(3,2,2) !(:,beta,alpha) = \frac{\partial \vec{n}}{\partial \alpha} \cross \frac{\partial cvat1}{\partial \beta}
      real :: n3deS_a1(3,3) 
      integer al,al2
      real :: stmp1,stmp2
      real :: vtmp(3),vtmp1(3),vtmp2(3),vtmp3(3)

      if(ibend_dummy.eq.1) then
      do j=1,nlon
      do i=1,nlat
      a1   = cvat1(:,i,j)
      a2   = cvat2(:,i,j)
      n3   = cvan3(:,i,j)
      a1de = cderi1(:,:,i,j)
      a2de = cderi2(:,:,i,j)
      n3de = cderi3(:,:,i,j)

      adeS = cdeSvat(:,:,i,j)



      call vcross(a2_n3(1),a2_n3(2),a2_n3(3),&
                 a2(1),   a2(2),   a2(3),&
                 n3(1),   n3(2),   n3(3),&
     1)
      call vcross(n3_a1(1),n3_a1(2),n3_a1(3),&
                 n3(1),   n3(2),   n3(3),&
                 a1(1),   a1(2),   a1(3),&
     1)

      do al = 1,2
      call vcross(a2de_n3(1,al),a2de_n3(2,al),a2de_n3(3,al),&
                 a2de(1,al),   a2de(2,al),   a2de(3,al),&
                 n3(1),        n3(2),        n3(3),&
     1)
      call vcross(a2_n3de(1,al),a2_n3de(2,al),a2_n3de(3,al),&
                 a2(1),        a2(2),        a2(3),&
                 n3de(1,al),   n3de(2,al),   n3de(3,al),&
     1)

      call vcross(n3de_a1(1,al),n3de_a1(2,al),n3de_a1(3,al),&
                 n3de(1,al),   n3de(2,al),   n3de(3,al),&
                 a1(1),        a1(2),        a1(3),&
     1)

      call vcross(n3_a1de(1,al),n3_a1de(2,al),n3_a1de(3,al),&
                 n3(1),        n3(2),        n3(3),&
                 a1de(1,al),   a1de(2,al),   a1de(3,al),&
     1)
      end do

      call vdot3(g,a1(1),a1(2),a1(3),&
     a2_n3(1),a2_n3(2),a2_n3(3),&
     1)

      do al=1,2
      vtmp = a2de_n3(:,al) + a2_n3de(:,al)
      call vdot3(stmp1,a1de(1,al),a1de(2,al),a1de(3,al),&
     a2_n3(1),a2_n3(2),a2_n3(3),&
     1)
      call vdot3(stmp2,a1(1),a1(2),a1(3),&
     vtmp(1),vtmp(2),vtmp(3),&
     1)
      gde(al) = stmp1 + stmp2
      end do

      do al=1,3
      call vcross(a2_n3deS(1,al),a2_n3deS(2,al),a2_n3deS(3,al),&
                 a2(1),         a2(2),         a2(3),&
     cderiS3(1,al,i,j), cderiS3(2,al,i,j), cderiS3(3,al,i,j),&
     1)      
      end do

      do al=1,2
      do al2=1,2
      call vcross(a2de_n3de(1,al2,al),a2de_n3de(2,al2,al),&
     a2de_n3de(3,al2,al),&
     a2de(1,al),a2de(2,al),a2de(3,al),&
     n3de(1,al2),n3de(2,al2),n3de(3,al2),1)
      end do
      end do

      
      call vcross(a2deS_n3(1,1),a2deS_n3(2,1),a2deS_n3(3,1),&
     cdeSvat(1,2,i,j),cdeSvat(2,2,i,j),cdeSvat(3,2,i,j),&
     n3(1),n3(2),n3(3),1)
      call vcross(a2deS_n3(1,2),a2deS_n3(2,2),a2deS_n3(3,2),&
     cdeSvat(1,4,i,j),cdeSvat(2,4,i,j),cdeSvat(3,4,i,j),&
     n3(1),n3(2),n3(3),1)
      call vcross(a2deS_n3(1,3),a2deS_n3(2,3),a2deS_n3(3,3),&
     cdeSvat(1,3,i,j),cdeSvat(2,3,i,j),cdeSvat(3,3,i,j),&
     n3(1),n3(2),n3(3),1)
      

      call vcross(n3_a1deS(1,1),n3_a1deS(2,1),n3_a1deS(3,1),&
     n3(1),n3(2),n3(3),&
     cdeSvat(1,1,i,j),cdeSvat(2,1,i,j),cdeSvat(3,1,i,j),1)
      call vcross(n3_a1deS(1,2),n3_a1deS(2,2),n3_a1deS(3,2),&
     n3(1),n3(2),n3(3),&
     cdeSvat(1,3,i,j),cdeSvat(2,3,i,j),cdeSvat(3,3,i,j),1)
      call vcross(n3_a1deS(1,3),n3_a1deS(2,3),n3_a1deS(3,3),&
     n3(1),n3(2),n3(3),&
     cdeSvat(1,2,i,j),cdeSvat(2,2,i,j),cdeSvat(3,2,i,j),1)


      do al=1,2
      do al2=1,2
      call vcross(n3de_a1de(1,al2,al),n3de_a1de(2,al2,al),&
     n3de_a1de(3,al2,al),&
     n3de(1,al), n3de(2,al), n3de(3,al),&
     a1de(1,al2),a1de(2,al2),a1de(3,al2),1)
      end do
      end do

      
      do al=1,3
      call vcross(n3deS_a1(1,al),n3deS_a1(2,al),n3deS_a1(3,al),&
     cderiS3(1,al,i,j), cderiS3(2,al,i,j), cderiS3(3,al,i,j),&
                 a1(1),         a1(2),         a1(3),&
     1)      
      end do

      
      
      call vdot3(tmp1,adeS(1,1),adeS(2,1),adeS(3,1),&
     a2_n3(1),a2_n3(2),a2_n3(3),1)
      vtmp1 = a2de_n3(:,1)+a2_n3de(:,1)
      call vdot3(tmp2,a1de(1,1),a1de(2,1),a1de(3,1),&
     vtmp1(1),vtmp1(2),vtmp1(3),1)
      vtmp2 = a2de_n3(:,1)+a2_n3de(:,1)
      call vdot3(tmp3,a1de(1,1),a1de(2,1),a1de(3,1),&
     vtmp2(1),vtmp2(2),vtmp2(3),1)
      vtmp3 = a2deS_n3(:,1)+a2de_n3de(:,1,1)+a2de_n3de(:,1,1)+&
     a2_n3deS(:,1)
      call vdot3(tmp4,a1(1),a1(2),a1(3),&
     vtmp3(1),vtmp3(2),vtmp3(3),1)
      gdeS(1) = tmp1 + tmp2 + tmp3 + tmp4

      call vdot3(tmp1,adeS(1,3),adeS(2,3),adeS(3,3),&
     a2_n3(1),a2_n3(2),a2_n3(3),1)
      vtmp1 = a2de_n3(:,2)+a2_n3de(:,2)
      call vdot3(tmp2,a1de(1,2),a1de(2,2),a1de(3,2),&
     vtmp1(1),vtmp1(2),vtmp1(3),1)
      vtmp2 = a2de_n3(:,2)+a2_n3de(:,2)
      call vdot3(tmp3,a1de(1,2),a1de(2,2),a1de(3,2),&
     vtmp2(1),vtmp2(2),vtmp2(3),1)
      vtmp3 = a2deS_n3(:,2)+a2de_n3de(:,2,2)+a2de_n3de(:,2,2)+&
     a2_n3deS(:,2)
      call vdot3(tmp4,a1(1),a1(2),a1(3),&
     vtmp3(1),vtmp3(2),vtmp3(3),1)
      gdeS(2) = tmp1 + tmp2 + tmp3 + tmp4

      call vdot3(tmp1,adeS(1,2),adeS(2,2),adeS(3,2),&
     a2_n3(1),a2_n3(2),a2_n3(3),1)
      vtmp1 = a2de_n3(:,2)+a2_n3de(:,2)
      call vdot3(tmp2,a1de(1,1),a1de(2,1),a1de(3,1),&
     vtmp1(1),vtmp1(2),vtmp1(3),1)
      vtmp2 = a2de_n3(:,1)+a2_n3de(:,1)
      call vdot3(tmp3,a1de(1,2),a1de(2,2),a1de(3,2),&
     vtmp2(1),vtmp2(2),vtmp2(3),1)
      vtmp3 = a2deS_n3(:,3)+a2de_n3de(:,2,1)+a2de_n3de(:,1,2)+&
     a2_n3deS(:,3)
      call vdot3(tmp4,a1(1),a1(2),a1(3),&
     vtmp3(1),vtmp3(2),vtmp3(3),1)
      gdeS(3) = tmp1 + tmp2 + tmp3 + tmp4

      ginv1=1./g;ginv2=ginv1*ginv1
      ginv3=ginv1*ginv2
      !!\frac{\partial^2 cona1}{\partial alpha \partial beta}
      condeS1(:,1,i,j)=2*ginv3*gde(1)*gde(1)*a2_n3&
     -ginv2*gdeS(1)*a2_n3&
     -ginv2*gde(1)*(a2de_n3(:,1)+a2_n3de(:,1))&
     -ginv2*gde(1)*(a2de_n3(:,1)+a2_n3de(:,1))&
     +ginv1*(a2deS_n3(:,1)+a2de_n3de(:,1,1)+a2de_n3de(:,1,1)+&
     a2_n3deS(:,1))
      condeS1(:,2,i,j)=2*ginv3*gde(2)*gde(2)*a2_n3&
     -ginv2*gdeS(2)*a2_n3&
     -ginv2*gde(2)*(a2de_n3(:,2)+a2_n3de(:,2))&
     -ginv2*gde(2)*(a2de_n3(:,2)+a2_n3de(:,2))&
     +ginv1*(a2deS_n3(:,2)+a2de_n3de(:,2,2)+a2de_n3de(:,2,2)+&
     a2_n3deS(:,2))
      condeS1(:,3,i,j)=2*ginv3*gde(2)*gde(1)*a2_n3&
     -ginv2*gdeS(3)*a2_n3&
     -ginv2*gde(1)*(a2de_n3(:,2)+a2_n3de(:,2))&
     -ginv2*gde(2)*(a2de_n3(:,1)+a2_n3de(:,1))&
     +ginv1*(a2deS_n3(:,3)+a2de_n3de(:,2,1)+a2de_n3de(:,1,2)+&
     a2_n3deS(:,3))
      !!\frac{\partial^2 cona2}{\partial alpha \partial beta}
      condeS2(:,1,i,j)=2*ginv3*gde(1)*gde(1)*n3_a1&
     -ginv2*gdeS(1)*n3_a1&
     -ginv2*gde(1)*(n3de_a1(:,1)+n3_a1de(:,1))&
     -ginv2*gde(1)*(n3de_a1(:,1)+n3_a1de(:,1))&
     +ginv1*(n3deS_a1(:,1)+n3de_a1de(:,1,1)+n3de_a1de(:,1,1)+&
     n3_a1deS(:,1))
      condeS2(:,2,i,j)=2*ginv3*gde(2)*gde(2)*n3_a1&
     -ginv2*gdeS(2)*n3_a1&
     -ginv2*gde(2)*(n3de_a1(:,2)+n3_a1de(:,2))&
     -ginv2*gde(2)*(n3de_a1(:,2)+n3_a1de(:,2))&
     +ginv1*(n3deS_a1(:,2)+n3de_a1de(:,2,2)+n3de_a1de(:,2,2)+&
     n3_a1deS(:,2))
      condeS2(:,3,i,j)=2*ginv3*gde(2)*gde(1)*n3_a1&
     -ginv2*gdeS(3)*n3_a1&
     -ginv2*gde(1)*(n3de_a1(:,2)+n3_a1de(:,2))&
     -ginv2*gde(2)*(n3de_a1(:,1)+n3_a1de(:,1))&
     +ginv1*(n3deS_a1(:,3)+n3de_a1de(:,2,1)+n3de_a1de(:,1,2)+&
     n3_a1deS(:,3))

      end do
      end do
      end if
      return
      end subroutine cal_contra_second_derivative





      subroutine cal_unit_normal(vat1,vat2,van3,nlat,nlon)
      integer nlat,nlon,ntotal
      real :: vat1(3,nlat,nlon)
      real :: vat2(3,nlat,nlon)
      real :: van3(3,nlat,nlon)
      ntotal = nlat*nlon
      call vcross(van3(1,:,:),van3(2,:,:),van3(3,:,:),&
                 vat1(1,:,:),vat1(2,:,:),vat1(3,:,:),&
                 vat2(1,:,:),vat2(2,:,:),vat2(3,:,:),ntotal)
      do j=1,nlon
         do i=1,nlat
            van3(:,i,j)=van3(:,i,j)/&
           rnorm2(van3(1,i,j),van3(2,i,j),van3(3,i,j))
         end do
      end do
   
      return
      end subroutine cal_unit_normal

      subroutine cal_unit_normal1(van3,vat1,vat2)
      real :: vat1(3)
      real :: vat2(3)
      real :: van3(3)
      call vcross(van3(1),van3(2),van3(3),&
                 vat1(1),vat1(2),vat1(3),&
                 vat2(1),vat2(2),vat2(3),1)
      van3=van3/&
     rnorm2(van3(1),van3(2),van3(3))
   
      return
      end subroutine cal_unit_normal1


      subroutine cal_unit_normal_deri(van3,devan3,&
     surfmet,&
     vat1,vat2,&
     devat1,devat2,&
     nlat,nlon)
      integer nlat,nlon,ntotal
      real :: vat1(3,nlat,nlon)
      real :: vat2(3,nlat,nlon)
      real :: van3(3,nlat,nlon)
      real :: devat1(3,2,nlat,nlon)
      real :: devat2(3,2,nlat,nlon)
      real :: devan3(3,2,nlat,nlon)
      real :: surfmet(nlat,nlon)

      do j=1,nlon
      do i=1,nlat
         call cal_unit_normal_deri1(van3(:,i,j),&
                               devan3(:,:,i,j),&
                                  surfmet(i,j),&
                                   vat1(:,i,j),&
                                   vat2(:,i,j),&
                               devat1(:,:,i,j),&
                               devat2(:,:,i,j))
      end do
      end do
   
      return
      end subroutine cal_unit_normal_deri


      subroutine cal_unit_normal_deri1(n,n_deri,surfmet,a1,a2,&
     ade1,ade2)
      !\vec{n} is the unit out normal vector, defined as 
      !\vec{c} = \vec{a1} \crossproduct \vec{a2} \vec{n} = \frac{\vec{c}}{|\vec{c}|}
      !\vec{ade1_al} = \frac{\partial \vec{a1}}{\partial x} 
      !\vec{ade2_be} = \frac{\partial \vec{a2}}{\partial x}
      !\n_deri(:,1) = \frac{\partial \vec{n}}{ \partial \alpha}
      !\n_deri(:,2) = \frac{\partial \vec{n}}{ \partial \beta} 
      real :: n_deri(3,2),n(3)
      real :: surfmet
      real :: a1(3),a2(3)
      real :: ade1(3,2),ade2(3,2)
      real :: c(3)
      real :: cross1(3),cross2(3)
      real :: tmp
      integer al
      call vcross(c (1), c(2), c(3),&
                 a1(1),a1(2),a1(3),&
                 a2(1),a2(2),a2(3),&
                 1)
      surfmet = rnorm2(c(1),c(2),c(3))
      n = c/surfmet

      do al = 1,2
      !calculate \frac{\partial \vec{n}}{ \partial \alpha}
      call vcross(cross1(1),  cross1(2),   cross1(3),&
                ade1(1,al), ade1(2,al),  ade1(3,al),&
                     a2(1),      a2(2),      a2(3),&
                                                 1)
      call vcross(cross2(1),  cross2(2),   cross2(3),&
                     a1(1),      a1(2),       a1(3),&
                ade2(1,al), ade2(2,al),  ade2(3,al),&
                                                  1)
      tmp = deri_norm_cross(a1,a2,ade1(:,al),ade2(:,al))
      n_deri(:,al) = 1./surfmet*(cross1+cross2-tmp*n)
!      write(*,*) cross1(1),cross2(1),tmp
      !calculate \frac{\partial \vec{n}}{ \partial \alpha}
      end do

      return
      end subroutine cal_unit_normal_deri1





      subroutine filter_modes_unitout(xout,yout,zout,&
                              xin ,yin ,zin ,&
                              theout,phiout,&
                              nltout,nlnout,&
                              nltin, nlnin)

!      include 'SIZE'
      real :: xout(nltout,nlnout)
      real :: yout(nltout,nlnout)
      real :: zout(nltout,nlnout)
      real :: theout(nltout)
      real :: phiout(nlnout)
      real :: stmp(nltin,nlnin,3)
      real :: xin(nltin,nlnin)
      real :: yin(nltin,nlnin)
      real :: zin(nltin,nlnin)
      real :: a(nltin,nltin,3)
      real :: b(nltin,nltin,3)
      integer,save :: lmax,icald,nmodes
      real ,allocatable,save :: pbar(:,:,:), cosmp(:,:), sinmp(:,:)
      real :: pbarthe0(3)
      real :: mphi
      data icald  /0/
      if(icald.eq.0) then
         nmodes = nltin
         if(nltin.gt.nltout) nmodes=nltout
         lmax = nltin

         allocate(pbar(lmax+1,0:nmodes-1,nltout))
         allocate(cosmp(0:nmodes-1,nlnout),sinmp(0:nmodes-1,nlnout))
         do i=1,nltout
            cost = cos(theout(i))
            do m=0,nmodes-1
               call cal_legendre_array(pbar(:,m,i),lmax,m,cost)
            end do
         end do
         pbar(:,0,:)=0.5*pbar(:,0,:)
         do j=1,nlnout
            phi = phiout(j)
            do m=0,nmodes-1
               cosmp(m,j)=cos(m*phi)
               sinmp(m,j)=sin(m*phi)
            end do
         end do

         icald = 1
      end if

      xout = 0; yout = 0; zout = 0
      !!!!!!!!!!remove higher number modes!!!!!!!!! 
      call cal_3dcrd_spectral_coeff(a,b,xin,yin,zin,nltin,nlnin)

      do m = 0,nmodes-1
      do n = m,nmodes-1
         ind = n-m+1
         do j = 1,nlnout
         pbarthe0 = (a(m+1,n+1,:)*cosmp(m,j)-b(m+1,n+1,:)*sinmp(m,j))  
         xout(:,j)=xout(:,j)+pbar(ind,m,:)*pbarthe0(1)
         yout(:,j)=yout(:,j)+pbar(ind,m,:)*pbarthe0(2)
         zout(:,j)=zout(:,j)+pbar(ind,m,:)*pbarthe0(3)
         end do
      end do
      end do
      return
      end subroutine filter_modes_unitout

      subroutine filter_modes4__OPT(xout,yout,zout,&
                          prstsmaxout,prstsminout,&
                          eshrout,edilout,ebenout,&
                          pmaxoutpk,pminoutpk,&
                          esoutpk,edoutpk,eboutpk,&
                                   xin ,yin ,zin ,&
                            prstsmaxin,prstsminin,&
                             eshrin,edilin,ebenin,&
                                    theout,phiout,&
                                    nltout,nlnout,&
                                     nltin, nlnin,&
                              ifprstress,ifenergy,&
                                      xoe,yoe,zoe,&
                                      theoe,phioe)

!      include 'SIZE'
!      use mod_filter_modes
      real :: xout(nltout,nlnout),yout(nltout,nlnout),zout(nltout,nlnout)
      real :: prstsmaxout(nltout,nlnout),prstsminout(nltout,nlnout)
      real :: pmaxoutpk(2),pminoutpk(2)
      real :: esoutpk(2),edoutpk(2),eboutpk(2)
      real :: prstsmaxin(nltin,nlnin),prstsminin(nltin,nlnin)
      real :: eshrin(nltin,nlnin),edilin(nltin,nlnin),ebenin(nltin,nlnin)
      real :: eshrout(nltout,nlnout),edilout(nltout,nlnout),&
      ebenout(nltout,nlnout)
      real ,optional,dimension(:,:) :: xoe,yoe,zoe
      real,optional,dimension(:) :: theoe,phioe
      real :: theout(nltout),phiout(nlnout)
      real :: thepk(2),phipk(2)
      real :: pbarpk(2,nltout+1,0:nltout-1)
      real :: stmp(nltin,nlnin,3)
      real :: xin(nltin,nlnin),yin(nltin,nlnin),zin(nltin,nlnin)
      real :: a(nltin,nltin,3),b(nltin,nltin,3)
      real :: a2(nltin,nltin,3),b2(nltin,nltin,3)
      real :: a3(nltin,nltin,3),b3(nltin,nltin,3)
      integer,save :: lmax,icald,nmodes,ifop,nltoe,nlnoe
      real ,allocatable,save :: pbar(:,:,:), cosmp(:,:), sinmp(:,:)
      real ,allocatable,save :: pbare(:,:,:), cosmpe(:,:), sinmpe(:,:)
      real :: pbarthe0(3)
      real :: mphi
      data icald  /0/
      data ifop   /0/
      real :: xtmp(2,3)
      if(icald.eq.0) then
         nmodes = nltin
         if(nltin.gt.nltout) nmodes=nltout
         lmax = nltin
         if(present(xoe).and.present(yoe).and.present(zoe).and.&
           present(theoe).and.present(phioe)) then
            nltoe = size(xoe,1);nlnoe = size(xoe,2)
            ifop = 1
         allocate(pbare(lmax+1,0:nmodes-1,nltoe))
         allocate(cosmpe(0:nmodes-1,nlnoe),sinmpe(0:nmodes-1,nlnoe))
         end if
         allocate(pbar(lmax+1,0:nmodes-1,nltout))
         allocate(cosmp(0:nmodes-1,nlnout),sinmp(0:nmodes-1,nlnout))
         do i=1,nltout
            cost = cos(theout(i))
            do m=0,nmodes-1
               call cal_legendre_array(pbar(:,m,i),lmax,m,cost)
            end do
         end do
         pbar(:,0,:)=0.5*pbar(:,0,:)
         do j=1,nlnout
            phi = phiout(j)
            do m=0,nmodes-1
               cosmp(m,j)=cos(m*phi)
               sinmp(m,j)=sin(m*phi)
            end do
         end do
         if(ifop.eq.1) then
         do i=1,nltoe
            cost = cos(theoe(i))
            do m=0,nmodes-1
               call cal_legendre_array(pbare(:,m,i),lmax,m,cost)
            end do
         end do
         pbare(:,0,:)=0.5*pbare(:,0,:)
         do j=1,nlnoe
            phi = phioe(j)
            do m=0,nmodes-1
               cosmpe(m,j)=cos(m*phi)
               sinmpe(m,j)=sin(m*phi)
            end do
         end do
         end if
         icald = 1
      end if

      xout = 0; yout = 0; zout = 0

      if(ifprstress.eq.1) then
          prstsmaxout=0;prstsminout=0
      end if
      if(ifenergy.eq.1) then
         eshrout=0;edilout=0;ebenout=0
      end if

      if(ifop.eq.1) then
         xoe = 0; yoe = 0; zoe = 0
      end if


      !!!!!!!!!!remove higher number modes!!!!!!!!! 
      if(ifprstress.eq.1.and.ifenergy.eq.1) then
      call cal_3dcrd_spectral_coeff(a,b,xin,yin,zin,nltin,nlnin)
      call cal_3dcrd_spectral_coeff(a2,b2,prstsmaxin,prstsminin,yin,&
                                   nltin,nlnin)
      call cal_3dcrd_spectral_coeff(a3,b3,eshrin,edilin,ebenin,&
                                   nltin,nlnin)


      do m = 0,nmodes-1
      do n = m,nmodes-1
         ind = n-m+1
         do j = 1,nlnout
         pbarthe0 = (a(m+1,n+1,:)*cosmp(m,j)-b(m+1,n+1,:)*sinmp(m,j))  
         xout(:,j)=xout(:,j)+pbar(ind,m,:)*pbarthe0(1)
         yout(:,j)=yout(:,j)+pbar(ind,m,:)*pbarthe0(2)
         zout(:,j)=zout(:,j)+pbar(ind,m,:)*pbarthe0(3)
         pbarthe0(1:2) = (a2(m+1,n+1,1:2)*cosmp(m,j)-&
                         b2(m+1,n+1,1:2)*sinmp(m,j))  
         prstsmaxout(:,j)=prstsmaxout(:,j)+pbar(ind,m,:)*pbarthe0(1)
         prstsminout(:,j)=prstsminout(:,j)+pbar(ind,m,:)*pbarthe0(2)
         pbarthe0 = (a3(m+1,n+1,:)*cosmp(m,j)-b3(m+1,n+1,:)*sinmp(m,j))  
         eshrout(:,j)=eshrout(:,j)+pbar(ind,m,:)*pbarthe0(1)
         edilout(:,j)=edilout(:,j)+pbar(ind,m,:)*pbarthe0(2)
         ebenout(:,j)=ebenout(:,j)+pbar(ind,m,:)*pbarthe0(3)
         end do
      end do
      end do



      pmaxoutpk(1) = sum(prstsmaxout(1,:))/nlnout
      pmaxoutpk(2) = sum(prstsmaxout(nltout,:))/nlnout 
      pminoutpk(1) = sum(prstsminout(1,:))/nlnout
      pminoutpk(2) = sum(prstsminout(nltout,:))/nlnout 
      esoutpk(1)   = sum(eshrout(1,:))/nlnout
      esoutpk(2)   = sum(eshrout(nltout,:))/nlnout
      edoutpk(1)   = sum(edilout(1,:))/nlnout
      edoutpk(2)   = sum(edilout(nltout,:))/nlnout
      eboutpk(1)   = sum(ebenout(1,:))/nlnout
      eboutpk(2)   = sum(ebenout(nltout,:))/nlnout

      elseif(ifprstress.eq.1.and.ifenergy.eq.0) then
      call cal_3dcrd_spectral_coeff(a,b,xin,yin,zin,nltin,nlnin)
      call cal_3dcrd_spectral_coeff(a2,b2,prstsmaxin,prstsminin,yin,&
                                   nltin,nlnin)

      do m = 0,nmodes-1
      do n = m,nmodes-1
         ind = n-m+1
         do j = 1,nlnout
         pbarthe0 = (a(m+1,n+1,:)*cosmp(m,j)-b(m+1,n+1,:)*sinmp(m,j))  
         xout(:,j)=xout(:,j)+pbar(ind,m,:)*pbarthe0(1)
         yout(:,j)=yout(:,j)+pbar(ind,m,:)*pbarthe0(2)
         zout(:,j)=zout(:,j)+pbar(ind,m,:)*pbarthe0(3)
         pbarthe0(1:2) = (a2(m+1,n+1,1:2)*cosmp(m,j)-&
                         b2(m+1,n+1,1:2)*sinmp(m,j))  
         prstsmaxout(:,j)=prstsmaxout(:,j)+pbar(ind,m,:)*pbarthe0(1)
         prstsminout(:,j)=prstsminout(:,j)+pbar(ind,m,:)*pbarthe0(2)
         end do
      end do
      end do

      pmaxoutpk(1) = sum(prstsmaxout(1,:))/nlnout
      pmaxoutpk(2) = sum(prstsmaxout(nltout,:))/nlnout 
      pminoutpk(1) = sum(prstsminout(1,:))/nlnout
      pminoutpk(2) = sum(prstsminout(nltout,:))/nlnout 
      elseif(ifprstress.eq.0.and.ifenergy.eq.1) then
      call cal_3dcrd_spectral_coeff(a,b,xin,yin,zin,nltin,nlnin)
      call cal_3dcrd_spectral_coeff(a3,b3,eshrin,edilin,ebenin,&
                                   nltin,nlnin)

      do m = 0,nmodes-1
      do n = m,nmodes-1
         ind = n-m+1
         do j = 1,nlnout
         pbarthe0 = (a(m+1,n+1,:)*cosmp(m,j)-b(m+1,n+1,:)*sinmp(m,j))  
         xout(:,j)=xout(:,j)+pbar(ind,m,:)*pbarthe0(1)
         yout(:,j)=yout(:,j)+pbar(ind,m,:)*pbarthe0(2)
         zout(:,j)=zout(:,j)+pbar(ind,m,:)*pbarthe0(3)
         pbarthe0 = (a3(m+1,n+1,:)*cosmp(m,j)-b3(m+1,n+1,:)*sinmp(m,j))  
         eshrout(:,j)=eshrout(:,j)+pbar(ind,m,:)*pbarthe0(1)
         edilout(:,j)=edilout(:,j)+pbar(ind,m,:)*pbarthe0(2)
         ebenout(:,j)=ebenout(:,j)+pbar(ind,m,:)*pbarthe0(3)
         end do
      end do
      end do
      esoutpk(1)   = sum(eshrout(1,:))/nlnout
      esoutpk(2)   = sum(eshrout(nltout,:))/nlnout
      edoutpk(1)   = sum(edilout(1,:))/nlnout
      edoutpk(2)   = sum(edilout(nltout,:))/nlnout
      eboutpk(1)   = sum(ebenout(1,:))/nlnout
      eboutpk(2)   = sum(ebenout(nltout,:))/nlnout
      else
      call cal_3dcrd_spectral_coeff(a,b,xin,yin,zin,nltin,nlnin)
      do m = 0,nmodes-1
      do n = m,nmodes-1
         ind = n-m+1
         do j = 1,nlnout
         pbarthe0 = (a(m+1,n+1,:)*cosmp(m,j)-b(m+1,n+1,:)*sinmp(m,j))  
         xout(:,j)=xout(:,j)+pbar(ind,m,:)*pbarthe0(1)
         yout(:,j)=yout(:,j)+pbar(ind,m,:)*pbarthe0(2)
         zout(:,j)=zout(:,j)+pbar(ind,m,:)*pbarthe0(3)
         end do
      end do
      end do
      end if

      if(ifop.eq.1) then
      do m = 0,nmodes-1
      do n = m,nmodes-1
         ind = n-m+1
         do j = 1,nlnoe
         pbarthe0 = (a(m+1,n+1,:)*cosmpe(m,j)-b(m+1,n+1,:)*sinmpe(m,j))  
         xoe(:,j)=xoe(:,j)+pbare(ind,m,:)*pbarthe0(1)
         yoe(:,j)=yoe(:,j)+pbare(ind,m,:)*pbarthe0(2)
         zoe(:,j)=zoe(:,j)+pbare(ind,m,:)*pbarthe0(3)
         end do
      end do
      end do
      end if



      return
      end subroutine filter_modes4__OPT



 
      subroutine iout(ivar,nam)
      real :: nam
      write(6,10) nam , ivar
   10 format(1h a4, 3h = ,i8)
      return
      end subroutine iout
!
      subroutine vout(var,nam)
      real :: nam
      write(6,10) nam , var
   10 format(1h a4,3h = ,e12.5)
      return
      end subroutine vout
!
      subroutine name(nam)
      real :: nam
      write(6,100) nam
  100 format(1h a8)
      return
      end subroutine name
!
      subroutine vecout(vec,nam,len)
      dimension vec(len)
      real :: nam
      write(6,109) nam, (vec(l),l=1,len)
  109 format(1h a4,/(1h 8e11.4))
      return
      end subroutine vecout

      function deri_norm_cross(a,b,aderi,bderi)
      !a,b are vectors, varies with x. the function 
      !calculates the derivative of the norm of crossproduct of a and b
      !with respect to x.
      !aderi, bderi are derivatives of a and be with respect to x.
      real :: deri_norm_cross
      real :: a(3),b(3),aderi(3),bderi(3)
      real :: c(3)
      call vcross(c(1),c(2),c(3),&
                 a(1),a(2),a(3),&
                 b(1),b(2),b(3),&
                 1)
      deri_norm_cross=1.0/sqrt(c(1)**2+c(2)**2+c(3)**2)*&
     (c(1)*(aderi(2)*b(3)+a(2)*bderi(3)-aderi(3)*b(2)-a(3)*bderi(2))+&
      c(2)*(aderi(3)*b(1)+a(3)*bderi(1)-aderi(1)*b(3)-a(1)*bderi(3))+&
      c(3)*(aderi(1)*b(2)+a(1)*bderi(2)-aderi(2)*b(1)-a(2)*bderi(1)))
      return
      end function deri_norm_cross




      
      subroutine  surface_metric_xdotn(surfm,xdotn,&
                                        theta,phi,&
                                        a,      b,&
                                        nlat,nlon)

      real :: surfm,xdotn
      real :: vat1(3),vat2(3),van3(3),x(3)
      real :: theta,phi
      real :: a(nlat,nlat,3),b(nlat,nlat,3)
      call       cal_co_vat1(vat1,vat2,&
                            x        ,&
                            theta,phi,&
                            a,      b,&
                            nlat,nlon)      
      call vcross(van3(1),van3(2),van3(3),&
                 vat1(1),vat1(2),vat1(3),&
                 vat2(1),vat2(2),vat2(3),&
     1) 
      surfm = rnorm2(van3(1),van3(2),van3(3))

      van3 = van3/surfm !the unit normal vector
      call vdot3(xdotn,x(1),x(2),x(3), &
      van3(1),van3(2),van3(3),&
      1)
      !xdotn = \vec{x} \cdot \vec{n}
      return
      end subroutine  surface_metric_xdotn



      subroutine  surface_metric(surfm,x,&
                              theta,phi,&
                              a,      b,&
                              nlat,nlon)

      real :: surfm
      real :: vat1(3),vat2(3),van3(3),x(3)
      real :: theta,phi
      real :: a(nlat,nlat,3),b(nlat,nlat,3)
      call       cal_co_vat1(vat1,vat2,&
                            x        ,&
                            theta,phi,&
                            a,      b,&
                            nlat,nlon)      
      call vcross(van3(1),van3(2),van3(3),&
                 vat1(1),vat1(2),vat1(3),&
                 vat2(1),vat2(2),vat2(3),&
     1) 
      surfm = rnorm2(van3(1),van3(2),van3(3))
      return
      end subroutine  surface_metric



      subroutine cal_weight_vol1_old(weight   ,&
                             int_xdotn,&
                             phi      ,&
                             a,      b,&
                             indt     ,&
                             thetag   ,&
                             nlat,nlon,&
                             xi,    wt,&
                             nq       )
! weight is the nodal weight, like the surface area of the patch surrounding the node.
! int_xdotn is the integration of &\vec{x}\cdot\vec{n}& on the patch surrounding the node.
! which will be used to calculate the volumn.
      real :: weight,int_xdotn
      real :: theta,phi
      real :: a(nlat,nlat,3),b(nlat,nlat,3)
      real :: thetag(nlat)
      real :: dphi
      real :: the1,the2,phi1,phi2
      real :: xi(nq)
      real :: wt(nq)
      real :: tmpthe,tmpph
      real :: vat1(3), vat2(3)
      real :: tmpwgt,tmp_xdotn
      real :: surfm,xdotn !xdotn = \vec{x} \cdot \vec{n}
      real :: coef1,coef2
      weight = 0
      int_xdotn=0
      !theta = thetag(indt)
      !integration is performed in the range theta:[the1, the2]
      !                                      phi  :[phi1, phi2] 

      dphi = 2*pi/nlon
      phi1 = phi - 0.5*dphi
      phi2 = phi + 0.5*dphi
      if(indt.eq.1) then
      the1 = 0
      the2 = 0.5*(thetag(indt) + thetag(indt+1))
      else if (indt.eq.nlat) then
      the1 = 0.5*(thetag(indt) + thetag(indt-1))
      the2 = pi
      else
      the1 = 0.5*(thetag(indt) + thetag(indt-1))
      the2 = 0.5*(thetag(indt) + thetag(indt+1))
      end if


      do j = 1,nq !in the theta direction
      tmpthe = 0.5*(xi(j)*(the2-the1)+(the1+the2))
      tmpwgt = 0
      tmp_xdotn = 0
      do i = 1,nq !in the phi   direction
      tmpphi = 0.5*(xi(i)*(phi2-phi1)+(phi1+phi2))
!      write(*,'(5f12.8)') tmpthe,xi(j),the1,the2,cos(tmpthe)
      call surface_metric_xdotn(surfm,  xdotn,&
                               tmpthe,tmpphi,&
                               a,          b,&
                               nlat,    nlon)

      coef1 = wt(i)*surfm
      tmpwgt = tmpwgt + coef1
      tmp_xdotn = tmp_xdotn + coef1*xdotn
      end do
      coef2 = 0.5*(phi2-phi1)*wt(j)
      weight = weight + coef2*tmpwgt
      int_xdotn = int_xdotn + coef2*tmp_xdotn
      end do
      weight = weight*0.5*(the2-the1)
      int_xdotn = int_xdotn*0.5*(the2-the1)
!      write(*,*) weight
      return
      end subroutine cal_weight_vol1_old      



      subroutine cal_weight_vol1(weight,&
                             int_xdotn,&
                             phi      ,&
                             a,      b,&
                             indt     ,&
                             thetag   ,&
                             nlat,nlon,&
                        xila,wtla,nqla,&
                        xilo,wtlo,nqlo)
! weight is the nodal weight, like the surface area of the patch surrounding the node.
! int_xdotn is the integration of &\vec{x}\cdot\vec{n}& on the patch surrounding the node.
! which will be used to calculate the volumn.
      real :: weight,int_xdotn
      real :: theta,phi
      real :: a(nlat,nlat,3),b(nlat,nlat,3)
      real :: thetag(nlat)
      real :: the1,the2,phi1,phi2
      real :: xila(nqla),wtla(nqla)
      real :: xilo(nqlo),wtlo(nqlo)
      real :: tmpthe,tmpph
      real :: tmpwgt,tmp_xdotn
      real :: surfm,xdotn !xdotn = \vec{x} \cdot \vec{n}
      real :: coef1,coef2
      weight = 0
      int_xdotn=0
      !theta = thetag(indt)
      !integration is performed in the range theta:[the1, the2]
      !                                      phi  :[phi1, phi2] 
      call cal_limit_theta_phi(the1,the2,phi1,phi2,&
     phi,indt,nlat,nlon,thetag)
     
      do j = 1,nqla !in the theta direction
      tmpthe = 0.5*(xila(j)*(the2-the1)+(the1+the2))
      tmpwgt = 0
      tmp_xdotn = 0
      do i = 1,nqlo !in the phi   direction
      tmpphi = 0.5*(xilo(i)*(phi2-phi1)+(phi1+phi2))
!      write(*,'(5f12.8)') tmpthe,xi(j),the1,the2,cos(tmpthe)
      call surface_metric_xdotn(surfm,  xdotn,&
                               tmpthe,tmpphi,&
                               a,          b,&
                               nlat,    nlon)
      coef1 = wtlo(i)*surfm
      tmpwgt = tmpwgt + coef1
      tmp_xdotn = tmp_xdotn + coef1*xdotn
      end do
      coef2 = 0.5*(phi2-phi1)*wtla(j)
      weight = weight + coef2*tmpwgt
      int_xdotn = int_xdotn + coef2*tmp_xdotn
      end do
      weight = weight*0.5*(the2-the1)
      int_xdotn = int_xdotn*0.5*(the2-the1)
!      write(*,*) weight
      return
      end subroutine cal_weight_vol1


      subroutine cal_weight_vol_spectral(weight   ,&
                                        vol      ,&
                                        a,      b,&
                                        thetag   ,&
                                        nlat,nlon)
      real :: weight(nlat,nlon)
      real :: vol(nlat,nlon)
! volume can be interpreted by surface integration, vol(i,j) is the surface integration of the patch surrounding that node(i,j).
      real :: a(nlat,nlat,3)
      real :: b(nlat,nlat,3)
      real :: thetag(nlat)
      real :: xi(3)
      real :: wgt(3)
      real :: dphi
      real :: phi
      real :: vol1
!      real :: tmp
!      real :: vol
      tmp = 0
      vol1 = 0
      xi(1) = -sqrt(0.6)
      xi(2) = 0
      xi(3) = -xi(1)
      wgt(1) = 5./9
      wgt(2) = 8./9
      wgt(3) = wgt(1)

      dphi = 2*pi/nlon
      nq = 3
       do j= 1,nlon
         phi = (j-1)*dphi
       do i= 1,nlat
       call cal_weight_vol1(weight(i,j),&
                           vol   (i,j),&
                           phi      ,&
                           a,      b,&
                           i        ,&
                           thetag   ,&
                           nlat,nlon,&
                           evlxlat,wgtlat,nqlat,&
                           evlxlon,wgtlon,nqlon) 


      end do
      end do
      vol = 1./3*vol
!      write(*,*) tmp,vol1
      return
      end subroutine cal_weight_vol_spectral


      subroutine sum_area_vol_spectral(area,volumn,&
                                               weight,vol,&   
                                                nlat,nlon)
!      include 'SIZE'
      real :: area,volumn
      real :: weight(nlat,nlon)
      real :: vol(nlat,nlon)

      area = 0
      volumn = 0
      do j= 1,nlon
      do i= 1,nlat
      area = area + weight(i,j)
      volumn = volumn + vol(i,j)
!      write(*,*) nid,weight(i,j)
      end do
      end do
      return
      end subroutine sum_area_vol_spectral


      subroutine cal_area_vol_postupdate(area,volumn,&
                                              x,y,z,&
                                             thetag,&
                                          nlat,nlon)

      !after update the cell position, calculating its area and volumn for volumn correction.
      real :: area,volumn
      real :: x(nlat,nlon), y(nlat,nlon), z(nlat,nlon)
      real :: a(nlat,nlat,3)
      real :: b(nlat,nlat,3)
      real :: weight(nlat,nlon)
      real :: vol(nlat,nlon)
      real :: thetag(nlat)
      call cal_3dcrd_spectral_coeff(a,b,x,y,z,nlat,nlon)
      call cal_weight_vol_spectral(      weight   ,&
                                        vol      ,&
                                        a,      b,&
                                        thetag   ,&
                                        nlat,nlon)

      call sum_area_vol_spectral(area,volumn,&
                                 weight,vol,&   
                                  nlat,nlon)

      return
      end subroutine cal_area_vol_postupdate

      subroutine calvolcor_spectral(xcor,ycor,zcor,&
                                   unitnormal, &
                                   cvol,carea,&
                                   ivol,&
                                   nlat,nlon)

      real :: xcor(nlat,nlon)
      real :: ycor(nlat,nlon)
      real :: zcor(nlat,nlon)
      real :: unitnormal(3,nlat,nlon)
      real :: cvol,carea,ivol
      real :: coef
      coef = -1.*(cvol-ivol)/carea
      xcor = coef*unitnormal(1,:,:)
      ycor = coef*unitnormal(2,:,:)
      zcor = coef*unitnormal(3,:,:)
      return
      end subroutine calvolcor_spectral






      subroutine update_vol_cor_spectral
!      use addvar
      real :: xcor(nnlat2,nnlon2)
      real :: ycor(nnlat2,nnlon2)
      real :: zcor(nnlat2,nnlon2)
      real :: area,volumn !ivol: volumn in the reference status

      call cal_area_vol_postupdate(area,volumn,&
                            crdxD,crdyD,crdzD,&
                                       thetaO,&
                                nnlat2,nnlon2)      
!!!!!!!!!!!!!must use the unit normal of the original points !!!!!!!!!!!    
!!!!!!!!!!!!unit is for the new points
      call calvolcor_spectral(xcor,ycor,zcor,&
                                      unitO,& 
                                volumn,area,&
                                      cvoli,&
                              nnlat2,nnlon2)
!!!!!!!!!!!!!must use the unit normal of the original points !!!!!!!!!!!    
      crdxD = crdxD + xcor
      crdyD = crdyD + ycor
      crdzD = crdzD + zcor
!      if(nid.eq.0) then
      write(*,'(i5,a,6f16.12)') nid,' maximum correction',maxval(xcor),&
     maxval(ycor),maxval(zcor),volumn-cvoli,volumn,cvoli

      return
      end subroutine update_vol_cor_spectral

      subroutine update_vol_cor_spectral2
!      use addvar
      real :: xcor(nnlat2,nnlon2)
      real :: ycor(nnlat2,nnlon2)
      real :: zcor(nnlat2,nnlon2)
      real :: area,volumn !ivol: volumn in the reference status
      real :: a(nnlat2,nnlon2,3),b(nnlat2,nnlon2,3)
      call   cal_3dcrd_spectral_coeff(a,b,crdxD,crdyD,crdzD,&
     nnlat2,nnlon2)



      call cal_cell_area_vol_new(area,volumn,a,b,wgtdsin,&
     crdxD,crdyD,crdzD,thetaO,phiO,nnlat2,nnlon2,nnhm)

      call calvolcor_spectral(xcor,ycor,zcor,&
                                      unitO,& 
                                volumn,area,&
                                      cvoli,&
                              nnlat2,nnlon2)
!!!!!!!!!!!!!must use the unit normal of the original points !!!!!!!!!!!    
      crdxD = crdxD + xcor
      crdyD = crdyD + ycor
      crdzD = crdzD + zcor
!      if(nid.eq.0) then
      write(*,'(i5,a,6f16.12)') nid,' maximum correction',maxval(xcor),&
     maxval(ycor),maxval(zcor),volumn-cvoli

      return
      end subroutine update_vol_cor_spectral2



      subroutine cal_unitO_filter(unitout,&
     unitin,&
     theout,phiout,&
     nltout,nlnout,&
     nltin, nlnin)
      real :: unitout(3,nltout,nlnout)
      real :: unitin (3,nltin ,nlnin)
      real :: theout(nltout)
      real :: phiout(nlnout)
      real :: rnorm
      call filter_modes_unitout&
                      (unitout(1,:,:),unitout(2,:,:),unitout(3,:,:),&
                       unitin(1,:,:), unitin(2,:,:), unitin(3,:,:),&
                       theout,phiout,&
                       nltout,nlnout,&
                       nltin, nlnin)
      do j = 1,nlnout
      do i = 1,nltout
         rnorm = rnorm2(unitout(1,i,j),unitout(2,i,j),unitout(3,i,j))
         unitout(:,i,j)=unitout(:,i,j)/rnorm

      end do
      end do
      return
      end subroutine cal_unitO_filter




      subroutine cal_limit_theta_phi(the1,the2,phi1,phi2,&
     phi,indt,nlat,nlon,thetag)
      real :: the1,the2,phi1,phi2
      real :: phi
      real :: thetag(nlat)
      real :: pi


      dphi = 2*pi/nlon
      phi1 = phi - 0.5*dphi
      phi2 = phi + 0.5*dphi
      if(indt.eq.1) then
      the1 = 0
      the2 = 0.5*(thetag(indt) + thetag(indt+1))
      else if (indt.eq.nlat) then
      the1 = 0.5*(thetag(indt) + thetag(indt-1))
      the2 = pi
      else
      the1 = 0.5*(thetag(indt) + thetag(indt-1))
      the2 = 0.5*(thetag(indt) + thetag(indt+1))
      end if

      return
      end subroutine cal_limit_theta_phi


      subroutine cal_cell_shape(center,distdiff,&
                               x,y,z,&
                               nlat,nlon)

      real :: center(3),distdiff
      real :: x(nlat,nlon),y(nlat,nlon),z(nlat,nlon)
      real :: tx(nlat,nlon),ty(nlat,nlon),tz(nlat,nlon)
      real :: dist(nlat,nlon)
!      write(*,*) 'kkkkkkkk',avedist
      avedist = 0
      center(1) = vec2ave(x,nlat,nlon)
      center(2) = vec2ave(y,nlat,nlon)
      center(3) = vec2ave(z,nlat,nlon)
      tx = x - center(1)
      ty = y - center(2)
      tz = z - center(3)
      do j = 1,nlon
      do i = 1,nlat
      dist(i,j) = rnorm2(tx(i,j),ty(i,j),tz(i,j))   
      end do
      end do
      distdiff = maxval(dist)-minval(dist)
      return
      end subroutine cal_cell_shape


      subroutine cal_cell_center(center,&
                               x,y,z,&
                               nlat,nlon)

      real :: center(3)
      real :: x(nlat,nlon),y(nlat,nlon),z(nlat,nlon)
      center(1) = vec2ave(x,nlat,nlon)
      center(2) = vec2ave(y,nlat,nlon)
      center(3) = vec2ave(z,nlat,nlon)
      return
      end subroutine cal_cell_center
!***************************************************************calculates surface area and volume of the capsule********************************************************
      subroutine cal_cell_area_vol(area,vol,a,b,wgtdsin,&
      x,y,z,the,phi,nlat,nlon)
!!!!!!!!when calculating the initial volume, we need unit normal, and no dealiasing is used for the unit normal.
!      use sph_coeff
      real :: area,vol
      real :: a(nlat,nlat,3),b(nlat,nlat,3)
      real :: x(nlat,nlon),y(nlat,nlon),z(nlat,nlon)
      real :: wgtdsin(nlat,nlon),surfmet(nlat,nlon),wgt(nlat,nlon)
      real :: the(nlat),phi(nlon)
      real :: vat1(3),vat2(3),van3(3),xnouse(3),xdotn(nlat,nlon)
      call   cal_surface_metric_sph__OPT(surfmet,&
                                    a,b,the,phi,&
                                      nlat,nlon)
      wgt = surfmet*wgtdsin
      area = sum(wgt)
      do j = 1, nlon
         do i = 1, nlat
            call cal_co_vat1(vat1,vat2,&
                            xnouse   ,&
                        the(i),phi(j),&
                            a,      b,&
                            nlat,nlon)
            call vcross(van3(1),van3(2),van3(3),&
                       vat1(1),vat1(2),vat1(3),&
                       vat2(1),vat2(2),vat2(3),&
                                             1)
            van3 = van3/rnorm2(van3(1),van3(2),van3(3))
            call vdot3(xdotn(i,j),x(i,j),y(i,j),z(i,j),&
                              van3(1),van3(2),van3(3),&
                                                    1) 
         end do
      end do
      vol = sum(wgt*xdotn)/3.0


      

      end subroutine cal_cell_area_vol


      subroutine cal_crd_cell_peak(crdx,crdy,crdz,pbar,a,b,nlat)
      !!!!we assume that phi=0, theta=0 and pi correspond to 2 points
      real :: crdx(2),crdy(2),crdz(2)
      real :: pbar(2,nlat+1,0:nlat-1)
      real :: a(nlat,nlat,3),b(nlat,nlat,3)
      real :: tmp1(3)
      crdx=0; crdy=0; crdz=0
      do m=0,nlat-1
         do n=m,nlat-1
            ind=n-m+1
            crdx=crdx+pbar(:,ind,m)*a(m+1,n+1,1)
            crdy=crdy+pbar(:,ind,m)*a(m+1,n+1,2)
            crdz=crdz+pbar(:,ind,m)*a(m+1,n+1,3)
         end do
      end do
      end subroutine cal_crd_cell_peak



      subroutine cal_mean_curvature(H,conat,b,nlat,nlon)
      real :: H(nlat,nlon),conat(2,2,nlat,nlon),b(2,2,nlat,nlon)
      do j=1,nlon
         do i=1,nlat
            H(i,j)=0.5*(conat(1,1,i,j)*b(1,1,i,j)+&
                       conat(2,2,i,j)*b(2,2,i,j))+&
                       conat(1,2,i,j)*b(1,2,i,j)
         end do
      end do
      end subroutine cal_mean_curvature

      subroutine cal_deform_energy(eshr,edil,eben,I1,I2,H,&
     gs,csk,ebd,&
     nlat,nlon,ibend)
      real :: eshr(nlat,nlon),edil(nlat,nlon),eben(nlat,nlon)
      real :: I1(nlat,nlon),I2(nlat,nlon),H(nlat,nlon)
      integer ibend
      real :: gs,csk,ebd
      do j=1,nlon
         do i=1,nlat
            call calenergy(eshr(i,j),edil(i,j),gs,csk,I1(i,j),I2(i,j)) 
         end do
      end do

      if(ibend.eq.1) then
         eben = 2*ebd*H**2
      else
         eben = 0
      end if

      end subroutine cal_deform_energy









      end module mod_sph1



