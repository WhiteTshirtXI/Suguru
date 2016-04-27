     module mod_sph
     use mod_common
     use mod_common_mpi
     use mod_param
     use mod_sph1
     use mod_sph2
     use mod_sph3
     use mod_math
     public 
     contains

      



      subroutine capsule_init_cal_coeff  !(crdxR,crdyR,crdzR)      
      parameter (alpha_rbc=1.3858189)
      !real,intent(in):: crdxR(nnlat2,nnlon2),crdyR(nnlat2,nnlon2),crdzR(nnlat2,nnlon2)
      real   thetag(nnlat2),dwtsO(nnlat2),dwtsN(nnlat),dwtsOf(n_capf)
      real  dphi,dphif,rad
      real  dworkO(lldworkO),dworkOf(lldworkOf)
      real  dworkN(lldworkN)
      real  pi
      integer ldworkO,ldworkOf,ldworkN
      real  fac
      
      real,allocatable :: fcrd(:,:,:),fcrd2(:,:,:)
      real rm(3,3)
      real crdtmp(nnlat2,nnlon2,2)





        hf_gscell = 0.5*gscell




!me      fac = 1
      pi   = 4.0*atan(1.0)
      dphiO = 2*pi/nnlon2
      dphi = dphiO
      dphif=2*pi/m_capf
!me      facdfmx = 1
!me      facdfmy = 1
!me      facdfmz = 1
      ldworkOf = lldworkOf
      ldworkO = lldworkO
      ldworkN = lldworkN

!****************************************************this subroutine computes Gaussian distribution for theta ***************************
      call gaqd(nnlat2,thetaO,dwtsO,dworkO,ldworkO,ier)
      call gaqd(n_capf,thetaOf,dwtsOf,dworkOf,ldworkOf,ier)
      thetag = thetaO
      wtsO = dwtsO 
      wtsOf=dwtsOf
     if (ifshape=="s") then

         do i = 1,nnlat2
	    theta=thetag(i)
	    cost = cos(theta)
	    sint = sin(theta)
         do j = 1,nnlon2
            phiO(j) = (j-1)*dphi
	    phi  = phiO(j)
	    sinp = sin(phi)
	    cosp = cos(phi)
            crdzR(i,j) = cost
            crdyR(i,j) = sint*cosp 
            crdxR(i,j) = sint*sinp     
         end do
      end do

        do i = 1,n_capf
         theta=thetaOf(i)
	 cost = cos(theta)
	 sint = sin(theta)
         do j = 1,m_capf
          phiOf(j) = (j-1)*dphif
	  phi  = phiOf(j)
	  sinp = sin(phi)
	  cosp = cos(phi)
            crdzRf(i,j) = cost
            crdyRf(i,j) = sint*cosp 
            crdxRf(i,j) = sint*sinp     
         end do
      end do 
     elseif (ifshape=="b") then

         do i = 1,nnlat2
	    theta=thetag(i)
	    cost = cos(theta)
	    sint = sin(theta)
         do j = 1,nnlon2
            phiO(j) = (j-1)*dphi
	    phi  = phiO(j)
	    sinp = sin(phi)
	    cosp = cos(phi)            
            call cal_rad_rbc(theta,x_rad,y_rad,z_rad)
            if (i>int(floor((nnlat2+1.)/2.))) z_rad=-z_rad
            crdzR(i,j) = z_rad
            crdyR(i,j) = y_rad*cosp 
            crdxR(i,j) = x_rad*sinp 
         end do
      end do

        do i = 1,n_capf
         theta=thetaOf(i)
	 cost = cos(theta)
	 sint = sin(theta)
         do j = 1,m_capf
          phiOf(j) = (j-1)*dphif
	  phi  = phiOf(j)
	  sinp = sin(phi)
	  cosp = cos(phi)            
            call cal_rad_rbc(theta,x_rad,y_rad,z_rad)
            if (i>int(floor((n_capf+1.)/2.))) z_rad=-z_rad
            crdzRf(i,j) = z_rad
            crdyRf(i,j) = y_rad*cosp 
            crdxRf(i,j) = x_rad*sinp  
         end do
      end do
     endif


      crdxD=crdxR;crdyD=crdyR;crdzD=crdzR

      forall(j=1:nnlon2) wgtdsin(:,j) = wtsO/sin(thetaO)*dphiO
      forall(j=1:m_capf) wgtdsinf(:,j) = wtsOf/sin(thetaOf)*dphif
!      return

      if(idealias.eq.1) then
      call gaqd(nnlat,thetaN,dwtsN,dworkN,ldworkN,ier)
      do j = 1,nnlon
         phiN(j) = (j-1)*2*pi/nnlon
      end do
      else
      do j = 1,nnlat
      thetaN(j) = thetaO(j)
      end do
      do i = 1,nnlon
      phiN(i) = phiO(i)
      end do
      end if
!*********************************************************subroutine for finding a_mn and b_mn for legender polynomials***********************************
      call cal_3dcrd_spectral_coeff(acrdR,bcrdR,crdxR,crdyR,crdzR,&
     nnlat2,nnlon2)
      



!***************************************************subroutine calculating surface metrics(it is in the module spherharm3)********************************
      call cal_contra_metric_deri(conAR,&
                                 conARd1,conARd2,&
                                 BR, &
                                 BRd1,BRd2,&
                                 acrdR,bcrdR,&
                                 thetaN,phiN,&
                                 nnlat2,nnlon2,&
                                 nnlat,nnlon,&
                                 ibend)



!!!!!!!!!!!!!!!!!!!to get the volume and area of the preinflated capsule!!!!!!!!!!!!!!!!

!      call   cal_3dcrd_spectral_coeff(acrdD,bcrdD,crdxD,crdyD,crdzD,
!     $nnlat2,nnlon2)

!****************************************************************subroutine calculating surface metrics****************************************************

      call   cal_surface_metric_sph__OPT(surfmetD,&
                         acrdD,bcrdD,thetaO,phiO,&
                                   nnlat2,nnlon2)

          
!********************************************************************wgtd is area of each element**********************************************************
      wgtD = surfmetD*wgtdsin
!******************************************************computing surface and volume at the begining(located in spherharm1.f********************************
!!!!!!!!!!!!!!!!!!!to get the volume and area of the preinflated capsule!!!!!!!!!!!!!!!!
!      call cal_cell_area_vol(careai,cvoli,acrdD,bcrdD,wgtdsin,
!     $ crdxD,crdyD,crdzD,thetaO,phiO,nnlat2,nnlon2)


      fxL   = 0;fyL   = 0;fzL   = 0
!      call allocate_prinstress
!      call allocate_deform_energy

      return
      end subroutine capsule_init_cal_coeff
 
      subroutine cal_contra_metric_deri(at,&
                                    atd1,atd2,&
                                              bcur,&
                                   bcurde1,bcurde2,&
                                               a,b,&
                                       thetag,phig,&
                                         nlat2,nlon2,&
                                         nlat,nlon,&
                                             ibend)
!a11,a22,a12 are the contravarient tensor.
!ad11 is the deriavtive of a11, with respect to alpha and beta.
      real a   (nlat2,nlat2,3)
      real b   (nlat2,nlat2,3)

      real at(2,2,nlat,nlon)
      real atd1(2,2,2,nlat,nlon)
      real atd2(2,2,2,2,nlat,nlon)
      real cvat1(3,nlat,nlon)
      real cvat2(3,nlat,nlon)
      real cvan3(3,nlat,nlon)
      real surfmet(nlat,nlon)
      real cderi1(3,2,nlat,nlon)
      real cderi2(3,2,nlat,nlon)
      real cderi3(3,2,nlat,nlon)
      real convat1(3,nlat,nlon)
      real convat2(3,nlat,nlon)
      real convan3(3,nlat,nlon)
      real conderi1(3,2,nlat,nlon)
      real conderi2(3,2,nlat,nlon)
      real condeS1(3,3,nlat,nlon)
      real condeS2(3,3,nlat,nlon)
      real cdeSvat(3,4,nlat,nlon) ! S: second derivative covarient
      real cdeTvat(3,5,nlat,nlon) ! T: third derivative  covarient
      real condeSvat(3,4,nlat,nlon) ! S: second derivative contravarient
      real bcur(2,2,nlat,nlon)!curvature tensor.
      real bcurde1(2,2,2,nlat,nlon)
      real bcurde2(2,2,2,2,nlat,nlon)
      real cderiS3(3,3,nlat,nlon) ! \frac{\partial \vec{n}}{\partial \xi^{\alpha} \partial \xi^{\beta}}
      real thetag(nlat),phig(nlon)
!***********************************************************************************calculates covariant vectors***************************************************************
      call cal_co_all_basis(cvat1,cvat2,cvan3,&
      surfmet,&
      cderi1,cderi2,&
      cderi3,&
      cdeSvat,cdeTvat,&
      a,b,thetag,phig,&
      nlat2,nlon2,&
      nlat,nlon,&
      ibend)

      



      call cal_curvature(bcur,&
      cvan3,&
      cderi3,&
      cderi1,cderi2,&
      nlat,nlon,&
      ibend)

!      write(*,*) cderi1(:,1,:,:)
!********************************************************************************calculates contravariant vectors**************************************************************
      call cal_contra_all_basis(convat1,convat2,convan3,&
      conderi1,conderi2,&
      cvat1,cvat2,cvan3,&
      cderi1,cderi2,cderi3,&
      nlat,nlon) 

       call cal_curvature_derivative(bcurde1,bcurde2,&
      cderiS3,&
      bcur,&
      cvan3,&
      cderi3,&
      cderi1,cderi2,&
      cdeSvat,cdeTvat,&
      convat1,convat2,&
      conderi1,conderi2,&
      nlat,nlon,&
      ibend)

       call cal_contra_second_derivative(&
      condeS1,condeS2,&     
      cvat1,cvat2,cvan3,&
      cderi1,cderi2,cderi3,& 
      cdeSvat,&
      cderiS3,&
      nlat,nlon,&
      ibend)

      call cal_metric_deri(at,&
                          atd1,&
                          convat1,convat2,&
                          conderi1,conderi2,&                
                          nlat,nlon,&
                          ibend)


      call cal_contra_metric_second_deri(atd2,&
      convat1,convat2,&
      conderi1,conderi2,&
      condeS1,condeS2,&
      nlat,nlon,&
      ibend)      

      return
      end  subroutine cal_contra_metric_deri


      subroutine cal_co_contra_metric_deri(cvat1,cvat2,cvan3,&
                                                    surfmet,&
                                                        cat,&
                                                      catd1,&
                                                      conat,&
                                            conatd1,conatd2,&
                                            cdeSvat,cdeTvat,&
                                                   christof,&
                                                  chrisderi,&
                                                       bcur,&
                                            bcurde1,bcurde2,&
                                                        a,b,&
                                                thetag,phig,&
                                                nlat2,nlon2,&
                                                nlat,nlon,&
                                                      ibend)
!a11,a22,a12 are the contravarient tensor.
!ad11 is the deriavtive of a11, with respect to alpha and beta.
      real a(nlat2,nlat2,3),b(nlat2,nlat2,3)
      real cvat1(3,nlat,nlon)
      real cvat2(3,nlat,nlon)
      real cvan3(3,nlat,nlon)
      real surfmet(nlat,nlon)
      real cat(2,2,nlat,nlon),conat(2,2,nlat,nlon)
      real catd1(2,2,2,nlat,nlon)
!      real catd2(2,2,2,2,nlat,nlon)
      real conatd1(2,2,2,nlat,nlon)
      real conatd2(2,2,2,2,nlat,nlon)
      real cderi1(3,2,nlat,nlon)
      real cderi2(3,2,nlat,nlon)
      real cderi3(3,2,nlat,nlon)
      real convat1(3,nlat,nlon)
      real convat2(3,nlat,nlon)
      real convan3(3,nlat,nlon)
      real conderi1(3,2,nlat,nlon)
      real conderi2(3,2,nlat,nlon)
      real condeS1(3,3,nlat,nlon)
      real condeS2(3,3,nlat,nlon)
      real christof(2,2,2,nlat,nlon)
      real chrisderi(2,2,2,2,nlat,nlon)
      real bcur(2,2,nlat,nlon)!curvature tensor.
      real bcurde1(2,2,2,nlat,nlon)
      real bcurde2(2,2,2,2,nlat,nlon)
      real thetag(nlat),phig(nlon)
      real cdeSvat(3,4,nlat,nlon) ! S: second derivative
      real cdeTvat(3,5,nlat,nlon) ! T: third derivative      
      real cderiS3(3,3,nlat,nlon) ! \frac{\partial \vec{n}}{\partial \xi^{\alpha} \partial \xi^{\beta}}
      real condeSvat(3,4,nlat,nlon) ! contravarient



      call cal_co_all_basis(cvat1,cvat2,cvan3,&
      surfmet,&
      cderi1,cderi2,&
      cderi3,&
      cdeSvat,cdeTvat,&
      a,b,thetag,phig,&
      nlat2,nlon2,&
      nlat,nlon,&
      ibend)
!      write(*,*) cderi1(:,1,:,:)



!********************************************************calculates curvature******************************************************************
      
      call cal_curvature(bcur,&
      cvan3,&
      cderi3,&
      cderi1,cderi2,&
      nlat,nlon,&
      ibend)
!      write(*,*) b11
!      if(nid.eq.0) write(*,*) bcur




!*********************************************************calculation of contravariant vectors(in spherharm1.f)********************************
      call cal_contra_all_basis(convat1,convat2,convan3,&
      conderi1,conderi2,&
      cvat1,cvat2,cvan3,&
      cderi1,cderi2,cderi3,&
      nlat,nlon) 

      call cal_curvature_derivative(bcurde1,bcurde2,&
      cderiS3,&
      bcur,&
      cvan3,&
      cderi3,&
      cderi1,cderi2,&
      cdeSvat,cdeTvat,&
      convat1,convat2,&
      conderi1,conderi2,&
      nlat,nlon,&
      ibend)

      call cal_contra_second_derivative(&
      condeS1,condeS2,&     
      cvat1,cvat2,cvan3,&
      cderi1,cderi2,cderi3,& 
      cdeSvat,&
      cderiS3,&
      nlat,nlon,&
      ibend)
      
      call cal_christof(christof,&
      chrisderi,&
      convat1,convat2,&
      cderi1,cderi2,&
      cdeSvat,&
      conderi1,conderi2,&
      nlat,nlon,&
      ibend)

      call cal_metric_deri(cat,&
                          catd1,&
                          cvat1,cvat2,&
                          cderi1,cderi2,&                
                          nlat,nlon,&
                          ibend)

      call cal_metric_deri(conat,&
                          conatd1,&
                          convat1,convat2,&
                          conderi1,conderi2,&                
                          nlat,nlon,&
                          ibend)

      call cal_contra_metric_second_deri(conatd2,&
      convat1,convat2,&
      conderi1,conderi2,&
      condeS1,condeS2,&
      nlat,nlon,&
      ibend)

      return
      end subroutine cal_co_contra_metric_deri

!**********************************************************calculates curvature bcur is Hessian matrix***************************************************
      subroutine cal_curvature(bcur,&
      cvan3,&
      cderi3,&
      cderi1,cderi2,&
      nlat,nlon,&
      ibend)

      real bcur(2,2,nlat,nlon)
      real bcurde1(2,2,2,nlat,nlon)
      real bcurde2(2,2,2,2,nlat,nlon)
      !bcur(b,a,:,:) = b_{ab}
      !bcurde1(c,b,a,:,:) = \frac{\partial b_{ab}}{\partial \xi^{c}}
      !bcurde2(d,c,b,a,:,:) = \frac{\partial^{2g} b_{ab}}{\partial \xi^{c} \partial \xi^{d}}      
      real cvan3(3,nlat,nlon)
      real cderi1(3,2,nlat,nlon)
      real cderi2(3,2,nlat,nlon)
      real cderi3(3,2,nlat,nlon)
      do j = 1,nlon
      do i = 1,nlat
      call vdot3(bcur(1,1,i,j),&
     cvan3(1,i,j),cvan3(2,i,j),cvan3(3,i,j),&
     cderi1(1,1,i,j),cderi1(2,1,i,j),cderi1(3,1,i,j),1)
!      write(*,*)cderi1(1,1,i,j),cderi1(2,1,i,j),cderi1(3,1,i,j)
      call vdot3(bcur(2,2,i,j),&
     cvan3(1,i,j),cvan3(2,i,j),cvan3(3,i,j),&
     cderi2(1,2,i,j),cderi2(2,2,i,j),cderi2(3,2,i,j),1)
      call vdot3(bcur(1,2,i,j),&
     cvan3(1,i,j),cvan3(2,i,j),cvan3(3,i,j),&
     cderi1(1,2,i,j),cderi1(2,2,i,j),cderi1(3,2,i,j),1)
      bcur(2,1,i,j) = bcur(1,2,i,j)
      end do
      end do

      return
      end subroutine cal_curvature


      subroutine cal_curvature_derivative(&
      bcurde1,bcurde2,&
      cderiS3,&
      bcur,&
      cvan3,&
      cderi3,&
      cderi1,cderi2,&
      cdeSvat,cdeTvat,&
      convat1,convat2,&
      conderi1,conderi2,&
      nlat,nlon,&
      ibend)

      real convat1(3,nlat,nlon)
      real convat2(3,nlat,nlon)
      real conderi1(3,2,nlat,nlon)
      real conderi2(3,2,nlat,nlon)
      real bcur(2,2,nlat,nlon)
      real bcurde1(2,2,2,nlat,nlon)
      real bcurde2(2,2,2,2,nlat,nlon)
      !bcur(b,a,:,:) = b_{ab}
      !bcurde1(c,b,a,:,:) = \frac{\partial b_{ab}}{\partial \xi^{c}}
      !bcurde2(d,c,b,a,:,:) = \frac{\partial^{2g} b_{ab}}{\partial \xi^{c} \partial \xi^{d}}      
      real cvan3(3,nlat,nlon)
      real cderi1(3,2,nlat,nlon)
      real cderi2(3,2,nlat,nlon)
      real cderi3(3,2,nlat,nlon)
      real cdeSvat(3,4,nlat,nlon) ! S: second derivative
      real cdeTvat(3,5,nlat,nlon) ! T: third derivative      
      real cderiS3(3,3,nlat,nlon) ! \frac{\partial \vec{n}}{\partial \xi^{\alpha} \partial \xi^{\beta}}
      !cderiS3(:,1,i,j) for \frac{\partial^2 \vec{n}}{\partial \xi^{1}^2}}
      !cderiS3(:,2,i,j) for \frac{\partial^2 \vec{n}}{\partial \xi^{2}^2}}
      !cderiS3(:,3,i,j) for \frac{\partial^2 \vec{n}}{\partial \xi^{1} \partial \xi^{2}}}

!!!!!!!!!with bending moment!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(ibend.eq.1) then
      do j = 1,nlon
      do i = 1,nlat
!!!!!!!!!!!!!!for bcurde1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!alpha=1,beta=1,gamma=1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call vdot3(tmp1,&
     cdeSvat(1,1,i,j),cdeSvat(2,1,i,j),cdeSvat(3,1,i,j),&
     cvan3(1,i,j),cvan3(2,i,j),cvan3(3,i,j),1)
      call vdot3(tmp2,&
     cderi1(1,1,i,j),cderi1(2,1,i,j),cderi1(3,1,i,j),&
     cderi3(1,1,i,j),cderi3(2,1,i,j),cderi3(3,1,i,j),1)
      bcurde1(1,1,1,i,j) = tmp1 + tmp2
!!!!!!!!alpha=1,beta=1,gamma=2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call vdot3(tmp1,&
     cdeSvat(1,2,i,j),cdeSvat(2,2,i,j),cdeSvat(3,2,i,j),&
     cvan3(1,i,j),cvan3(2,i,j),cvan3(3,i,j),1)
      call vdot3(tmp2,&
     cderi1(1,1,i,j),cderi1(2,1,i,j),cderi1(3,1,i,j),&
     cderi3(1,2,i,j),cderi3(2,2,i,j),cderi3(3,2,i,j),1)
      bcurde1(2,1,1,i,j) = tmp1 + tmp2
!!!!!!!!alpha=2,beta=2,gamma=1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call vdot3(tmp1,&
     cdeSvat(1,3,i,j),cdeSvat(2,3,i,j),cdeSvat(3,3,i,j),&
     cvan3(1,i,j),cvan3(2,i,j),cvan3(3,i,j),1)
      call vdot3(tmp2,&
     cderi2(1,2,i,j),cderi2(2,2,i,j),cderi2(3,2,i,j),&
     cderi3(1,1,i,j),cderi3(2,1,i,j),cderi3(3,1,i,j),1)
      bcurde1(1,2,2,i,j) = tmp1 + tmp2
!!!!!!!!alpha=2,beta=2,gamma=2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call vdot3(tmp1,&
     cdeSvat(1,4,i,j),cdeSvat(2,4,i,j),cdeSvat(3,4,i,j),&
     cvan3(1,i,j),cvan3(2,i,j),cvan3(3,i,j),1)
      call vdot3(tmp2,&
     cderi2(1,2,i,j),cderi2(2,2,i,j),cderi2(3,2,i,j),&
     cderi3(1,2,i,j),cderi3(2,2,i,j),cderi3(3,2,i,j),1)
      bcurde1(2,2,2,i,j) = tmp1 + tmp2
!!!!!!!!alpha=2,beta=1,gamma=1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call vdot3(tmp1,&
     cdeSvat(1,2,i,j),cdeSvat(2,2,i,j),cdeSvat(3,2,i,j),&
     cvan3(1,i,j),cvan3(2,i,j),cvan3(3,i,j),1)
      call vdot3(tmp2,&
     cderi2(1,1,i,j),cderi2(2,1,i,j),cderi2(3,1,i,j),&
     cderi3(1,1,i,j),cderi3(2,1,i,j),cderi3(3,1,i,j),1)
      bcurde1(1,1,2,i,j) = tmp1 + tmp2
      bcurde1(1,2,1,i,j) = bcurde1(1,1,2,i,j)
!!!!!!!!alpha=2,beta=1,gamma=2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call vdot3(tmp1,&
     cdeSvat(1,3,i,j),cdeSvat(2,3,i,j),cdeSvat(3,3,i,j),&
     cvan3(1,i,j),cvan3(2,i,j),cvan3(3,i,j),1)
      call vdot3(tmp2,&
     cderi2(1,1,i,j),cderi2(2,1,i,j),cderi2(3,1,i,j),&
     cderi3(1,2,i,j),cderi3(2,2,i,j),cderi3(3,2,i,j),1)
      bcurde1(2,1,2,i,j) = tmp1 + tmp2
      bcurde1(2,2,1,i,j) = bcurde1(2,1,2,i,j)
!!!!!!!!!!!!!!for bcurde1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!second derivative of unit normal!!!!!!!!!!!!!!!!!!!!
      cderiS3(:,1,i,j)=-(bcurde1(1,1,1,i,j)*convat1(:,i,j)+&
                        bcurde1(1,2,1,i,j)*convat2(:,i,j)+&
                        bcur(1,1,i,j)*conderi1(:,1,i,j)+&
                        bcur(1,2,i,j)*conderi2(:,1,i,j)) 
      cderiS3(:,2,i,j)=-(bcurde1(2,1,2,i,j)*convat1(:,i,j)+&
                        bcurde1(2,2,2,i,j)*convat2(:,i,j)+&
                        bcur(2,1,i,j)*conderi1(:,2,i,j)+&
                        bcur(2,2,i,j)*conderi2(:,2,i,j))
      cderiS3(:,3,i,j)=-(bcurde1(2,1,1,i,j)*convat1(:,i,j)+&
                        bcurde1(2,2,1,i,j)*convat2(:,i,j)+&
                        bcur(1,1,i,j)*conderi1(:,2,i,j)+&
                        bcur(1,2,i,j)*conderi2(:,2,i,j)) 
!!!!!!!!!!!!!!second derivative of unit normal!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!for bcurde2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call vdot3(tmp1,&
     cdeTvat(1,1,i,j),cdeTvat(2,1,i,j),cdeTvat(3,1,i,j),&
     cvan3(1,i,j),cvan3(2,i,j),cvan3(3,i,j),1)
      call vdot3(tmp2,&
     cdeSvat(1,1,i,j),cdeSvat(2,1,i,j),cdeSvat(3,1,i,j),&
     cderi3(1,1,i,j),cderi3(2,1,i,j),cderi3(3,1,i,j),1)
      call vdot3(tmp3,&
     cdeSvat(1,1,i,j),cdeSvat(2,1,i,j),cdeSvat(3,1,i,j),&
     cderi3(1,1,i,j),cderi3(2,1,i,j),cderi3(3,1,i,j),1)
      call vdot3(tmp4,&
     cderi1(1,1,i,j),cderi1(2,1,i,j),cderi1(3,1,i,j),&
     cderiS3(1,1,i,j),cderiS3(2,1,i,j),cderiS3(3,1,i,j),1)
      bcurde2(1,1,1,1,i,j)=tmp1+tmp2+tmp3+tmp4

      call vdot3(tmp1,&
     cdeTvat(1,2,i,j),cdeTvat(2,2,i,j),cdeTvat(3,2,i,j),&
     cvan3(1,i,j),cvan3(2,i,j),cvan3(3,i,j),1)
      call vdot3(tmp2,&
     cdeSvat(1,1,i,j),cdeSvat(2,1,i,j),cdeSvat(3,1,i,j),&
     cderi3(1,2,i,j),cderi3(2,2,i,j),cderi3(3,2,i,j),1)
      call vdot3(tmp3,&
     cdeSvat(1,2,i,j),cdeSvat(2,2,i,j),cdeSvat(3,2,i,j),&
     cderi3(1,1,i,j),cderi3(2,1,i,j),cderi3(3,1,i,j),1)
      call vdot3(tmp4,&
     cderi1(1,1,i,j),cderi1(2,1,i,j),cderi1(3,1,i,j),&
     cderiS3(1,3,i,j),cderiS3(2,3,i,j),cderiS3(3,3,i,j),1)
      bcurde2(2,1,1,1,i,j)=tmp1+tmp2+tmp3+tmp4
      bcurde2(1,2,1,1,i,j) = bcurde2(2,1,1,1,i,j)

      call vdot3(tmp1,&
     cdeTvat(1,3,i,j),cdeTvat(2,3,i,j),cdeTvat(3,3,i,j),&
     cvan3(1,i,j),cvan3(2,i,j),cvan3(3,i,j),1)
      call vdot3(tmp2,&
     cdeSvat(1,2,i,j),cdeSvat(2,2,i,j),cdeSvat(3,2,i,j),&
     cderi3(1,2,i,j),cderi3(2,2,i,j),cderi3(3,2,i,j),1)
      call vdot3(tmp3,&
     cdeSvat(1,2,i,j),cdeSvat(2,2,i,j),cdeSvat(3,2,i,j),&
     cderi3(1,2,i,j),cderi3(2,2,i,j),cderi3(3,2,i,j),1)
      call vdot3(tmp4,&
     cderi1(1,1,i,j),cderi1(2,1,i,j),cderi1(3,1,i,j),&
     cderiS3(1,2,i,j),cderiS3(2,2,i,j),cderiS3(3,2,i,j),1)
      bcurde2(2,2,1,1,i,j)=tmp1+tmp2+tmp3+tmp4


      call vdot3(tmp1,&
     cdeTvat(1,3,i,j),cdeTvat(2,3,i,j),cdeTvat(3,3,i,j),&
     cvan3(1,i,j),cvan3(2,i,j),cvan3(3,i,j),1)
      call vdot3(tmp2,&
     cdeSvat(1,3,i,j),cdeSvat(2,3,i,j),cdeSvat(3,3,i,j),&
     cderi3(1,1,i,j),cderi3(2,1,i,j),cderi3(3,1,i,j),1)
      call vdot3(tmp3,&
     cdeSvat(1,3,i,j),cdeSvat(2,3,i,j),cdeSvat(3,3,i,j),&
     cderi3(1,1,i,j),cderi3(2,1,i,j),cderi3(3,1,i,j),1)
      call vdot3(tmp4,&
     cderi2(1,2,i,j),cderi2(2,2,i,j),cderi2(3,2,i,j),&
     cderiS3(1,1,i,j),cderiS3(2,1,i,j),cderiS3(3,1,i,j),1)
      bcurde2(1,1,2,2,i,j)=tmp1+tmp2+tmp3+tmp4

      call vdot3(tmp1,&
     cdeTvat(1,4,i,j),cdeTvat(2,4,i,j),cdeTvat(3,4,i,j),&
     cvan3(1,i,j),cvan3(2,i,j),cvan3(3,i,j),1)
      call vdot3(tmp2,&
     cdeSvat(1,3,i,j),cdeSvat(2,3,i,j),cdeSvat(3,3,i,j),&
     cderi3(1,2,i,j),cderi3(2,2,i,j),cderi3(3,2,i,j),1)
      call vdot3(tmp3,&
     cdeSvat(1,4,i,j),cdeSvat(2,4,i,j),cdeSvat(3,4,i,j),&
     cderi3(1,1,i,j),cderi3(2,1,i,j),cderi3(3,1,i,j),1)
      call vdot3(tmp4,&
     cderi2(1,2,i,j),cderi2(2,2,i,j),cderi2(3,2,i,j),&
     cderiS3(1,3,i,j),cderiS3(2,3,i,j),cderiS3(3,3,i,j),1)
      bcurde2(2,1,2,2,i,j)=tmp1+tmp2+tmp3+tmp4
      bcurde2(1,2,2,2,i,j)=bcurde2(2,1,2,2,i,j)

      call vdot3(tmp1,&
     cdeTvat(1,5,i,j),cdeTvat(2,5,i,j),cdeTvat(3,5,i,j),&
     cvan3(1,i,j),cvan3(2,i,j),cvan3(3,i,j),1)
      call vdot3(tmp2,&
     cdeSvat(1,4,i,j),cdeSvat(2,4,i,j),cdeSvat(3,4,i,j),&
     cderi3(1,2,i,j),cderi3(2,2,i,j),cderi3(3,2,i,j),1)
      call vdot3(tmp3,&
     cdeSvat(1,4,i,j),cdeSvat(2,4,i,j),cdeSvat(3,4,i,j),&
     cderi3(1,2,i,j),cderi3(2,2,i,j),cderi3(3,2,i,j),1)
      call vdot3(tmp4,&
     cderi2(1,2,i,j),cderi2(2,2,i,j),cderi2(3,2,i,j),&
     cderiS3(1,2,i,j),cderiS3(2,2,i,j),cderiS3(3,2,i,j),1)
      bcurde2(2,2,2,2,i,j)=tmp1+tmp2+tmp3+tmp4



      call vdot3(tmp1,&
     cdeTvat(1,2,i,j),cdeTvat(2,2,i,j),cdeTvat(3,2,i,j),&
     cvan3(1,i,j),cvan3(2,i,j),cvan3(3,i,j),1)
      call vdot3(tmp2,&
     cdeSvat(1,2,i,j),cdeSvat(2,2,i,j),cdeSvat(3,2,i,j),&
     cderi3(1,1,i,j),cderi3(2,1,i,j),cderi3(3,1,i,j),1)
      call vdot3(tmp3,&
     cdeSvat(1,2,i,j),cdeSvat(2,2,i,j),cdeSvat(3,2,i,j),&
     cderi3(1,1,i,j),cderi3(2,1,i,j),cderi3(3,1,i,j),1)
      call vdot3(tmp4,&
     cderi1(1,2,i,j),cderi1(2,2,i,j),cderi1(3,2,i,j),&
     cderiS3(1,1,i,j),cderiS3(2,1,i,j),cderiS3(3,1,i,j),1)
      bcurde2(1,1,2,1,i,j)=tmp1+tmp2+tmp3+tmp4
      bcurde2(1,1,1,2,i,j)=bcurde2(1,1,2,1,i,j)


      call vdot3(tmp1,&
     cdeTvat(1,3,i,j),cdeTvat(2,3,i,j),cdeTvat(3,3,i,j),&
     cvan3(1,i,j),cvan3(2,i,j),cvan3(3,i,j),1)
      call vdot3(tmp2,&
     cdeSvat(1,2,i,j),cdeSvat(2,2,i,j),cdeSvat(3,2,i,j),&
     cderi3(1,2,i,j),cderi3(2,2,i,j),cderi3(3,2,i,j),1)
      call vdot3(tmp3,&
     cdeSvat(1,3,i,j),cdeSvat(2,3,i,j),cdeSvat(3,3,i,j),&
     cderi3(1,1,i,j),cderi3(2,1,i,j),cderi3(3,1,i,j),1)
      call vdot3(tmp4,&
     cderi1(1,2,i,j),cderi1(2,2,i,j),cderi1(3,2,i,j),&
     cderiS3(1,3,i,j),cderiS3(2,3,i,j),cderiS3(3,3,i,j),1)
      bcurde2(2,1,2,1,i,j)=tmp1+tmp2+tmp3+tmp4
      bcurde2(1,2,2,1,i,j)=bcurde2(2,1,2,1,i,j)
      bcurde2(2,1,1,2,i,j)=bcurde2(2,1,2,1,i,j)
      bcurde2(1,2,1,2,i,j)=bcurde2(2,1,2,1,i,j)

      call vdot3(tmp1,&
     cdeTvat(1,4,i,j),cdeTvat(2,4,i,j),cdeTvat(3,4,i,j),&
     cvan3(1,i,j),cvan3(2,i,j),cvan3(3,i,j),1)
      call vdot3(tmp2,&
     cdeSvat(1,3,i,j),cdeSvat(2,3,i,j),cdeSvat(3,3,i,j),&
     cderi3(1,2,i,j),cderi3(2,2,i,j),cderi3(3,2,i,j),1)
      call vdot3(tmp3,&
     cdeSvat(1,3,i,j),cdeSvat(2,3,i,j),cdeSvat(3,3,i,j),&
     cderi3(1,2,i,j),cderi3(2,2,i,j),cderi3(3,2,i,j),1)
      call vdot3(tmp4,&
     cderi1(1,2,i,j),cderi1(2,2,i,j),cderi1(3,2,i,j),&
     cderiS3(1,2,i,j),cderiS3(2,2,i,j),cderiS3(3,2,i,j),1)
      bcurde2(2,2,2,1,i,j)=tmp1+tmp2+tmp3+tmp4
      bcurde2(2,2,1,2,i,j)=bcurde2(2,2,2,1,i,j)
!!!!!!!!!!!!!!for bcurde2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do
      end do

      end if
!!!!!!!!!with bending moment!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      return
      end subroutine cal_curvature_derivative



      subroutine cal_christof(christof,&
      chrisderi,&
      convat1,convat2,&
      cderi1,cderi2,&
      cdeSvat,&
      conderi1,conderi2,&
      nlat,nlon,&
      ibend)
      !christof(\gamma,\beta,\alpha) = \Tau_{\alpha \gamma}^{\beta} = \vec{a}_{\alpha,\gamma}\cdot\vec{a^{\beta}}
      real christof(2,2,2,nlat,nlon)
      real chrisderi(2,2,2,2,nlat,nlon)
      real convat1(3,nlat,nlon)
      real convat2(3,nlat,nlon)
      real cderi1(3,2,nlat,nlon)
      real cderi2(3,2,nlat,nlon)
      real cdeSvat(3,4,nlat,nlon)
      real conderi1(3,2,nlat,nlon)
      real conderi2(3,2,nlat,nlon)
      real tmp1,tmp2
      do j = 1,nlon
      do i = 1,nlat
      call vdot3(christof(1,1,1,i,j),&
     cderi1(1,1,i,j),cderi1(2,1,i,j),cderi1(3,1,i,j),&
     convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),1)
      call vdot3(christof(2,1,1,i,j),&
     cderi1(1,2,i,j),cderi1(2,2,i,j),cderi1(3,2,i,j),&
     convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),1)
      call vdot3(christof(1,2,1,i,j),&
     cderi1(1,1,i,j),cderi1(2,1,i,j),cderi1(3,1,i,j),&
     convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),1)
      call vdot3(christof(2,2,1,i,j),&
     cderi1(1,2,i,j),cderi1(2,2,i,j),cderi1(3,2,i,j),&
     convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),1)



      call vdot3(christof(2,1,2,i,j),&
     cderi2(1,2,i,j),cderi2(2,2,i,j),cderi2(3,2,i,j),&
     convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),1)


      call vdot3(christof(2,2,2,i,j),&
     cderi2(1,2,i,j),cderi2(2,2,i,j),cderi2(3,2,i,j),&
     convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),1)

      end do
      end do
      christof(1,1,2,:,:) = christof(2,1,1,:,:)
      christof(1,2,2,:,:) = christof(2,2,1,:,:)

      if(ibend.eq.1) then

      do j=1,nlon
      do i=1,nlat
      !alpha=1,beta=1,lambda=1,l=1 
      call vdot3(tmp1,cdeSvat(1,1,i,j),cdeSvat(2,1,i,j),cdeSvat(3,1,i,j)&
     ,          convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),1 )
      call vdot3(tmp2,cderi1(1,1,i,j),cderi1(2,1,i,j),cderi1(3,1,i,j)&
     ,conderi1(1,1,i,j),conderi1(2,1,i,j),conderi1(3,1,i,j),1 )
      chrisderi(1,1,1,1,i,j) = tmp1 + tmp2
      !alpha=1,beta=1,lambda=1,l=2 
      call vdot3(tmp1,cdeSvat(1,2,i,j),cdeSvat(2,2,i,j),cdeSvat(3,2,i,j)&
     ,          convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),1 )
      call vdot3(tmp2,cderi1(1,1,i,j),cderi1(2,1,i,j),cderi1(3,1,i,j)&
     ,conderi1(1,2,i,j),conderi1(2,2,i,j),conderi1(3,2,i,j),1 )
      chrisderi(2,1,1,1,i,j) = tmp1 + tmp2
      !alpha=1,beta=1,lambda=2,l=1 
      call vdot3(tmp1,cdeSvat(1,2,i,j),cdeSvat(2,2,i,j),cdeSvat(3,2,i,j)&
     ,          convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),1 )
      call vdot3(tmp2,cderi1(1,2,i,j),cderi1(2,2,i,j),cderi1(3,2,i,j)&
     ,conderi1(1,1,i,j),conderi1(2,1,i,j),conderi1(3,1,i,j),1 )
      chrisderi(1,2,1,1,i,j) = tmp1 + tmp2
      !alpha=1,beta=1,lambda=2,l=2 
      call vdot3(tmp1,cdeSvat(1,3,i,j),cdeSvat(2,3,i,j),cdeSvat(3,3,i,j)&
     ,          convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),1 )
      call vdot3(tmp2,cderi1(1,2,i,j),cderi1(2,2,i,j),cderi1(3,2,i,j)&
     ,conderi1(1,2,i,j),conderi1(2,2,i,j),conderi1(3,2,i,j),1 )
      chrisderi(2,2,1,1,i,j) = tmp1 + tmp2
      !alpha=1,beta=2,lambda=1,l=1 
      call vdot3(tmp1,cdeSvat(1,1,i,j),cdeSvat(2,1,i,j),cdeSvat(3,1,i,j)&
     ,          convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),1 )
      call vdot3(tmp2,cderi1(1,1,i,j),cderi1(2,1,i,j),cderi1(3,1,i,j)&
     ,conderi2(1,1,i,j),conderi2(2,1,i,j),conderi2(3,1,i,j),1 )
      chrisderi(1,1,2,1,i,j) = tmp1 + tmp2
      !alpha=1,beta=2,lambda=1,l=2
      call vdot3(tmp1,cdeSvat(1,2,i,j),cdeSvat(2,2,i,j),cdeSvat(3,2,i,j)&
     ,          convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),1 )
      call vdot3(tmp2,cderi1(1,1,i,j),cderi1(2,1,i,j),cderi1(3,1,i,j)&
     ,conderi2(1,2,i,j),conderi2(2,2,i,j),conderi2(3,2,i,j),1 )
      chrisderi(2,1,2,1,i,j) = tmp1 + tmp2
      !alpha=1,beta=2,lambda=2,l=1 
      call vdot3(tmp1,cdeSvat(1,2,i,j),cdeSvat(2,2,i,j),cdeSvat(3,2,i,j)&
     ,          convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),1 )
      call vdot3(tmp2,cderi1(1,2,i,j),cderi1(2,2,i,j),cderi1(3,2,i,j)&
     ,conderi2(1,1,i,j),conderi2(2,1,i,j),conderi2(3,1,i,j),1 )
      chrisderi(1,2,2,1,i,j) = tmp1 + tmp2
      !alpha=1,beta=2,lambda=2,l=2 
      call vdot3(tmp1,cdeSvat(1,3,i,j),cdeSvat(2,3,i,j),cdeSvat(3,3,i,j)&
     ,          convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),1 )
      call vdot3(tmp2,cderi1(1,2,i,j),cderi1(2,2,i,j),cderi1(3,2,i,j)&
     ,conderi2(1,2,i,j),conderi2(2,2,i,j),conderi2(3,2,i,j),1 )
      chrisderi(2,2,2,1,i,j) = tmp1 + tmp2


!c$$$      !alpha=2,beta=1,lambda=1,l=1 
!c$$$      call vdot3(tmp1,cdeSvat(1,2,i,j),cdeSvat(2,2,i,j),cdeSvat(3,2,i,j)
!c$$$     $,          convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),1 )
!c$$$      call vdot3(tmp2,cderi2(1,1,i,j),cderi2(2,1,i,j),cderi2(3,1,i,j)
!c$$$     $,conderi1(1,1,i,j),conderi1(2,1,i,j),conderi1(3,1,i,j),1 )
!c$$$      chrisderi(1,1,1,2,i,j) = tmp1 + tmp2
!c$$$      !alpha=2,beta=1,lambda=1,l=2 
!c$$$      call vdot3(tmp1,cdeSvat(1,3,i,j),cdeSvat(2,3,i,j),cdeSvat(3,3,i,j)
!c$$$     $,          convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),1 )
!c$$$      call vdot3(tmp2,cderi2(1,1,i,j),cderi2(2,1,i,j),cderi2(3,1,i,j)
!c$$$     $,conderi1(1,2,i,j),conderi1(2,2,i,j),conderi1(3,2,i,j),1 )
!c$$$      chrisderi(2,1,1,2,i,j) = tmp1 + tmp2
!      !alpha=2,beta=1,lambda=2,l=1 
      call vdot3(tmp1,cdeSvat(1,3,i,j),cdeSvat(2,3,i,j),cdeSvat(3,3,i,j)&
     ,          convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),1 )
      call vdot3(tmp2,cderi2(1,2,i,j),cderi2(2,2,i,j),cderi2(3,2,i,j)&
     ,conderi1(1,1,i,j),conderi1(2,1,i,j),conderi1(3,1,i,j),1 )
      chrisderi(1,2,1,2,i,j) = tmp1 + tmp2
      !alpha=2,beta=1,lambda=2,l=2 
      call vdot3(tmp1,cdeSvat(1,4,i,j),cdeSvat(2,4,i,j),cdeSvat(3,4,i,j)&
     ,          convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),1 )
      call vdot3(tmp2,cderi1(1,2,i,j),cderi1(2,2,i,j),cderi1(3,2,i,j)&
     ,conderi2(1,2,i,j),conderi2(2,2,i,j),conderi2(3,2,i,j),1 )
      chrisderi(2,2,1,2,i,j) = tmp1 + tmp2
!c$$$      !alpha=2,beta=2,lambda=1,l=1 
!c$$$      call vdot3(tmp1,cdeSvat(1,2,i,j),cdeSvat(2,2,i,j),cdeSvat(3,2,i,j)
!c$$$     $,          convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),1 )
!c$$$      call vdot3(tmp2,cderi2(1,1,i,j),cderi2(2,1,i,j),cderi2(3,1,i,j)
!c$$$     $,conderi2(1,1,i,j),conderi2(2,1,i,j),conderi2(3,1,i,j),1 )
!c$$$      chrisderi(1,1,2,2,i,j) = tmp1 + tmp2
!c$$$      !alpha=2,beta=2,lambda=1,l=2
!c$$$      call vdot3(tmp1,cdeSvat(1,3,i,j),cdeSvat(2,3,i,j),cdeSvat(3,3,i,j)
!c$$$     $,          convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),1 )
!c$$$      call vdot3(tmp2,cderi2(1,1,i,j),cderi2(2,1,i,j),cderi2(3,1,i,j)
!c$$$     $,conderi2(1,2,i,j),conderi2(2,2,i,j),conderi2(3,2,i,j),1 )
!c$$$      chrisderi(2,1,2,2,i,j) = tmp1 + tmp2
      !alpha=2,beta=2,lambda=2,l=1 
      call vdot3(tmp1,cdeSvat(1,3,i,j),cdeSvat(2,3,i,j),cdeSvat(3,3,i,j)&
     ,          convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),1 )
      call vdot3(tmp2,cderi2(1,2,i,j),cderi2(2,2,i,j),cderi2(3,2,i,j)&
     ,conderi2(1,1,i,j),conderi2(2,1,i,j),conderi2(3,1,i,j),1 )
      chrisderi(1,2,2,2,i,j) = tmp1 + tmp2
      !alpha=2,beta=2,lambda=2,l=2 
      call vdot3(tmp1,cdeSvat(1,4,i,j),cdeSvat(2,4,i,j),cdeSvat(3,4,i,j)&
     ,          convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),1 )
      call vdot3(tmp2,cderi2(1,2,i,j),cderi2(2,2,i,j),cderi2(3,2,i,j)&
     ,conderi2(1,2,i,j),conderi2(2,2,i,j),conderi2(3,2,i,j),1 )
      chrisderi(2,2,2,2,i,j) = tmp1 + tmp2
      end do
      end do
      chrisderi(:,1,:,2,:,:) = chrisderi(:,2,:,1,:,:)
      end if

      return
      end subroutine cal_christof

      
      subroutine cal_metric_deri(at,&
                                atd1,&
                                vat1,vat2,&
                                deri1,deri2,&                
                                nlat,nlon,&
                                ibend)

      real at(2,2,nlat,nlon)
      real atd1(2,2,2,nlat,nlon)
      real vat1(3,nlat,nlon)
      real vat2(3,nlat,nlon)
      real deri1(3,2,nlat,nlon)
      real deri2(3,2,nlat,nlon)
      integer al
      do j = 1,nlon
      do i = 1,nlat
      at(1,1,i,j)=vat1(1,i,j)**2+&
              vat1(2,i,j)**2+&
              vat1(3,i,j)**2

      at(2,2,i,j)=vat2(1,i,j)**2+&
              vat2(2,i,j)**2+&
              vat2(3,i,j)**2

      call vdot3(at(1,2,i,j),&
                vat1(1,i,j),vat1(2,i,j),vat1(3,i,j),&
                vat2(1,i,j),vat2(2,i,j),vat2(3,i,j),&
                1)
      
      at(2,1,i,j) = at(1,2,i,j)
      do al = 1,2
      call vdot3(atd1(al,1,1,i,j),&
          vat1(1,i,j),vat1(2,i,j),vat1(3,i,j),&
          deri1(1,al,i,j),deri1(2,al,i,j),deri1(3,al,i,j),&
          1)
      atd1(al,1,1,i,j) = 2*atd1(al,1,1,i,j)

      call vdot3(atd1(al,2,2,i,j),&
          vat2(1,i,j),vat2(2,i,j),vat2(3,i,j),&
          deri2(1,al,i,j),deri2(2,al,i,j),deri2(3,al,i,j),&
          1)
      atd1(al,2,2,i,j) = 2*atd1(al,2,2,i,j)


      call vdot3(tmp1,&
          vat1(1,i,j),vat1(2,i,j),vat1(3,i,j),&
          deri2(1,al,i,j),deri2(2,al,i,j),deri2(3,al,i,j),&
          1)

      call vdot3(tmp2,&
          vat2(1,i,j),vat2(2,i,j),vat2(3,i,j),&
          deri1(1,al,i,j),deri1(2,al,i,j),deri1(3,al,i,j),&
          1)
      atd1(al,1,2,i,j) = tmp1 + tmp2
      atd1(al,2,1,i,j) = atd1(al,1,2,i,j)
      end do
      
      end do
      end do


      return
      end subroutine cal_metric_deri



      subroutine cal_contra_metric_second_deri(conatd2,&
      convat1,convat2,&
      conderi1,conderi2,&
      condeS1,condeS2,&
      nlat,nlon,&
      ibend)
! conatd2(l,lambda,beta,alpha) =
! \frac{\partial^2 a^{\alpha \beta}}{\partial \xi^{lambda} \partial \xi^{l}}
      real conatd2(2,2,2,2,nlat,nlon)
      real convat1(3,nlat,nlon)
      real convat2(3,nlat,nlon)
      real conderi1(3,2,nlat,nlon)
      real conderi2(3,2,nlat,nlon)
      real condeS1(3,3,nlat,nlon)
      real condeS2(3,3,nlat,nlon)
      real tmp1,tmp2,tmp3,tmp4
      if(ibend.eq.1) then
      do j=1,nlon
      do i=1,nlat
! alpha=1 beta=1 lambda=1 l=1   
      call vdot3(tmp1,condeS1(1,1,i,j),condeS1(2,1,i,j),condeS1(3,1,i,j)&
     ,convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),1)
      call vdot3(tmp2,  conderi1(1,1,i,j),conderi1(2,1,i,j),&
     conderi1(3,1,i,j),conderi1(1,1,i,j),conderi1(2,1,i,j),&
     conderi1(3,1,i,j),1)
      call vdot3(tmp3,  conderi1(1,1,i,j),conderi1(2,1,i,j),&
     conderi1(3,1,i,j),conderi1(1,1,i,j),conderi1(2,1,i,j),&
     conderi1(3,1,i,j),1)
      call vdot3(tmp4,convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),&
     condeS1(1,1,i,j),condeS1(2,1,i,j),condeS1(3,1,i,j),1)
      conatd2(1,1,1,1,i,j) = tmp1 + tmp2 + tmp3 + tmp4


! alpha=1 beta=1 lambda=2 l=2   
      call vdot3(tmp1,condeS1(1,2,i,j),condeS1(2,2,i,j),condeS1(3,2,i,j)&
     ,convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),1)
      call vdot3(tmp2,  conderi1(1,2,i,j),conderi1(2,2,i,j),&
     conderi1(3,2,i,j),conderi1(1,2,i,j),conderi1(2,2,i,j),&
     conderi1(3,2,i,j),1)
      call vdot3(tmp3,  conderi1(1,2,i,j),conderi1(2,2,i,j),&
     conderi1(3,2,i,j),conderi1(1,2,i,j),conderi1(2,2,i,j),&
     conderi1(3,2,i,j),1)
      call vdot3(tmp4,convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),&
     condeS1(1,2,i,j),condeS1(2,2,i,j),condeS1(3,2,i,j),1)
      conatd2(2,2,1,1,i,j) = tmp1 + tmp2 + tmp3 + tmp4

! alpha=1 beta=1 lambda=1(2) l=2(1)   
      call vdot3(tmp1,condeS1(1,3,i,j),condeS1(2,3,i,j),condeS1(3,3,i,j)&
     ,convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),1)
      call vdot3(tmp2,  conderi1(1,1,i,j),conderi1(2,1,i,j),&
     conderi1(3,1,i,j),conderi1(1,2,i,j),conderi1(2,2,i,j),&
     conderi1(3,2,i,j),1)
      call vdot3(tmp3,  conderi1(1,2,i,j),conderi1(2,2,i,j),&
     conderi1(3,2,i,j),conderi1(1,1,i,j),conderi1(2,1,i,j),&
     conderi1(3,1,i,j),1)
      call vdot3(tmp4,convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),&
     condeS1(1,3,i,j),condeS1(2,3,i,j),condeS1(3,3,i,j),1)
      conatd2(2,1,1,1,i,j) = tmp1 + tmp2 + tmp3 + tmp4
      conatd2(1,2,1,1,i,j) = conatd2(2,1,1,1,i,j)

! alpha=2 beta=2 lambda=1 l=1   
      call vdot3(tmp1,condeS2(1,1,i,j),condeS2(2,1,i,j),condeS2(3,1,i,j)&
     ,convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),1)
      call vdot3(tmp2,  conderi2(1,1,i,j),conderi2(2,1,i,j),&
     conderi2(3,1,i,j),conderi2(1,1,i,j),conderi2(2,1,i,j),&
     conderi2(3,1,i,j),1)
      call vdot3(tmp3,  conderi2(1,1,i,j),conderi2(2,1,i,j),&
     conderi2(3,1,i,j),conderi2(1,1,i,j),conderi2(2,1,i,j),&
     conderi2(3,1,i,j),1)
      call vdot3(tmp4,convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),&
     condeS2(1,1,i,j),condeS2(2,1,i,j),condeS2(3,1,i,j),1)
      conatd2(1,1,2,2,i,j) = tmp1 + tmp2 + tmp3 + tmp4


! alpha=2 beta=2 lambda=2 l=2   
      call vdot3(tmp1,condeS2(1,2,i,j),condeS2(2,2,i,j),condeS2(3,2,i,j)&
     ,convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),1)
      call vdot3(tmp2,  conderi2(1,2,i,j),conderi2(2,2,i,j),&
     conderi2(3,2,i,j),conderi2(1,2,i,j),conderi2(2,2,i,j),&
     conderi2(3,2,i,j),1)
      call vdot3(tmp3,  conderi2(1,2,i,j),conderi2(2,2,i,j),&
     conderi2(3,2,i,j),conderi2(1,2,i,j),conderi2(2,2,i,j),&
     conderi2(3,2,i,j),1)
      call vdot3(tmp4,convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),&
     condeS2(1,2,i,j),condeS2(2,2,i,j),condeS2(3,2,i,j),1)
      conatd2(2,2,2,2,i,j) = tmp1 + tmp2 + tmp3 + tmp4

! alpha=2 beta=2 lambda=1(2) l=2(1)   
      call vdot3(tmp1,condeS2(1,3,i,j),condeS2(2,3,i,j),condeS2(3,3,i,j)&
     ,convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),1)
      call vdot3(tmp2,  conderi2(1,1,i,j),conderi2(2,1,i,j),&
     conderi2(3,1,i,j),conderi2(1,2,i,j),conderi2(2,2,i,j),&
     conderi2(3,2,i,j),1)
      call vdot3(tmp3,  conderi2(1,2,i,j),conderi2(2,2,i,j),&
     conderi2(3,2,i,j),conderi2(1,1,i,j),conderi2(2,1,i,j),&
     conderi2(3,1,i,j),1)
      call vdot3(tmp4,convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),&
     condeS2(1,3,i,j),condeS2(2,3,i,j),condeS2(3,3,i,j),1)
      conatd2(2,1,2,2,i,j) = tmp1 + tmp2 + tmp3 + tmp4
      conatd2(1,2,2,2,i,j) = conatd2(2,1,2,2,i,j)


! alpha=1(2) beta=2(1) lambda=1(2) l=2(1)   
      call vdot3(tmp1,condeS1(1,3,i,j),condeS1(2,3,i,j),condeS1(3,3,i,j)&
     ,convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),1)
      call vdot3(tmp2,  conderi1(1,1,i,j),conderi1(2,1,i,j),&
     conderi1(3,1,i,j),conderi2(1,2,i,j),conderi2(2,2,i,j),&
     conderi2(3,2,i,j),1)
      call vdot3(tmp3,  conderi1(1,2,i,j),conderi1(2,2,i,j),&
     conderi1(3,2,i,j),conderi2(1,1,i,j),conderi2(2,1,i,j),&
     conderi2(3,1,i,j),1)
      call vdot3(tmp4,convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),&
     condeS2(1,3,i,j),condeS2(2,3,i,j),condeS2(3,3,i,j),1)
      conatd2(2,1,2,1,i,j) = tmp1 + tmp2 + tmp3 + tmp4
      conatd2(1,2,2,1,i,j) = conatd2(2,1,2,1,i,j)
      conatd2(2,1,1,2,i,j) = conatd2(2,1,2,1,i,j)
      conatd2(1,2,1,2,i,j) = conatd2(2,1,2,1,i,j)

! alpha=1(2) beta=2(1) lambda=1 l=1   
      call vdot3(tmp1,condeS1(1,1,i,j),condeS1(2,1,i,j),condeS1(3,1,i,j)&
     ,convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),1)
      call vdot3(tmp2,  conderi1(1,1,i,j),conderi1(2,1,i,j),&
     conderi1(3,1,i,j),conderi2(1,1,i,j),conderi2(2,1,i,j),&
     conderi2(3,1,i,j),1)
      call vdot3(tmp3,  conderi1(1,1,i,j),conderi1(2,1,i,j),&
     conderi1(3,1,i,j),conderi2(1,1,i,j),conderi2(2,1,i,j),&
     conderi2(3,1,i,j),1)
      call vdot3(tmp4,convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),&
     condeS2(1,1,i,j),condeS2(2,1,i,j),condeS2(3,1,i,j),1)
      conatd2(1,1,2,1,i,j) = tmp1 + tmp2 + tmp3 + tmp4
      conatd2(1,1,1,2,i,j) =conatd2(1,1,2,1,i,j)

! alpha=1(2) beta=2(1) lambda=2 l=2   
      call vdot3(tmp1,condeS1(1,2,i,j),condeS1(2,2,i,j),condeS1(3,2,i,j)&
     ,convat2(1,i,j),convat2(2,i,j),convat2(3,i,j),1)
      call vdot3(tmp2,  conderi1(1,2,i,j),conderi1(2,2,i,j),&
     conderi1(3,2,i,j),conderi2(1,2,i,j),conderi2(2,2,i,j),&
     conderi2(3,2,i,j),1)
      call vdot3(tmp3,  conderi1(1,2,i,j),conderi1(2,2,i,j),&
     conderi1(3,2,i,j),conderi2(1,2,i,j),conderi2(2,2,i,j),&
     conderi2(3,2,i,j),1)
      call vdot3(tmp4,convat1(1,i,j),convat1(2,i,j),convat1(3,i,j),&
     condeS2(1,2,i,j),condeS2(2,2,i,j),condeS2(3,2,i,j),1)
      conatd2(2,2,2,1,i,j) = tmp1 + tmp2 + tmp3 + tmp4
      conatd2(2,2,1,2,i,j) = conatd2(2,2,2,1,i,j)
      end do
      end do
      end if
      return
      end subroutine cal_contra_metric_second_deri



      subroutine cal_Js_I1_I2_deri(Js,I1,I2,&
     Jsde,I1de,I2de,&
     conAR,&
     cat,&
     conARd1,&
     catd1)
      real Js, I1, I2
      real Jsde(2),I1de(2),I2de(2)
      real conAR(2,2)
      real cat(2,2)
      real conARd1(2,2,2)
      real catd1(2,2,2)
      integer al
      real abs_conA, abs_ca
      real abs_conAd(2),&
          abs_cad(2)
      abs_conA = conAR(1,1)*conAR(2,2)-conAR(1,2)**2
      abs_ca   = cat(1,1)*cat(2,2)-cat(1,2)**2
      I1 = conAR(1,1)*cat(1,1) + conAR(2,2)*cat(2,2) +&
      2*conAR(1,2)*cat(1,2) - 2
      I2 = abs_conA*abs_ca - 1
      Js = sqrt(1 + I2)
      
      
      do al = 1,2
      I1de(al)=conARd1(al,1,1)*cat(1,1)+conAR(1,1)*catd1(al,1,1)+&
              conARd1(al,2,2)*cat(2,2)+conAR(2,2)*catd1(al,2,2)+&
           2*(conARd1(al,1,2)*cat(1,2)+conAR(1,2)*catd1(al,1,2))

      abs_conAd(al)=conARd1(al,1,1)*conAR(2,2)&
                +conAR(1,1)*conARd1(al,2,2)&
                -2*conAR(1,2)*conARd1(al,1,2)
      abs_cad(al)  =catd1(al,1,1)*cat(2,2)+cat(1,1)*catd1(al,2,2)&
                -2*cat(1,2)*catd1(al,1,2)

      I2de(al)=abs_conAd(al)*abs_ca+&
              abs_cad(al)  *abs_conA
      end do

      Jsde = 0.5/Js*I2de
      return
      end subroutine cal_Js_I1_I2_deri


      subroutine cal_pwspI(pwspI1,pwspI2,&
     pwspI1de,pwspI2de,&
     I1,I2,&
     I1de,I2de,&
     GsCELL,CSK)
      real pwspI1, pwspI2
      real pwspI1de(2), pwspI2de(2)
      real I1, I2
      real I1de(2), I2de(2)
      real GsCELL, CSK
      call calpwspi_deri(pwspI1,pwspI2,&
     pwspI1de,pwspI2de,&
     GsCELL,CSK,&
     I1,I2,&
     I1de,I2de)

      return
      end subroutine cal_pwspI


      subroutine calculate_cauchy_tensor_deri(T11,T22,T12,&
     Td11,Td22,Td12,&
     Js,pwspI1,pwspI2,&
     Jsde,pwspI1de,pwspI2de,&
     conAR,&
     conARd1,&
     conaD,&
     conaDd1)
      real T11,T22,T12
      real Td11(2),Td22(2),Td12(2)
      real Js,pwspI1,pwspI2
      real Jsde(2),pwspI1de(2),pwspI2de(2)
      real conAR(2,2)
      real conARd1(2,2,2)
      real conaD(2,2)
      real conaDd1(2,2,2)
      real coef1, coef2, vcoef(6)
      integer al
      coef1 = 2./Js*pwspI1
      coef2 = 2.*Js*pwspI2
      T11 = coef1*conAR(1,1) + coef2*conaD(1,1)
      T22 = coef1*conAR(2,2) + coef2*conaD(2,2)
      T12 = coef1*conAR(1,2) + coef2*conaD(1,2)

      do al = 1,2
      vcoef(1) = -2.0/Js**2   *Jsde(al)*pwspI1
      vcoef(2) =  2.0/Js      *pwspI1de(al)  
      vcoef(3) =  2.0/Js      *pwspI1        
      vcoef(4) =  2.0*Jsde(al)*pwspI2
      vcoef(5) =  2.0*Js      *pwspI2de(al)
      vcoef(6) =  2.0*Js      *pwspI2

      Td11(al) = vcoef(1)*conAR(1,1)      +& 
                vcoef(2)*conAR(1,1)      +&
                vcoef(3)*conARd1(al,1,1) +&
                vcoef(4)*conaD(1,1)      +&
                vcoef(5)*conaD(1,1)      +&
                vcoef(6)*conaDd1(al,1,1)

      Td22(al) = vcoef(1)*conAR(2,2)      +& 
                vcoef(2)*conAR(2,2)      +&
                vcoef(3)*conARd1(al,2,2) +& 
                vcoef(4)*conaD(2,2)      +&
                vcoef(5)*conaD(2,2)      +&
                vcoef(6)*conaDd1(al,2,2)

      Td12(al) = vcoef(1)*conAR(1,2)      +& 
                vcoef(2)*conAR(1,2)      +&
                vcoef(3)*conARd1(al,1,2) +& 
                vcoef(4)*conaD(1,2)      +&
                vcoef(5)*conaD(1,2)      +&
                vcoef(6)*conaDd1(al,1,2) 
      end do

      return
      end subroutine calculate_cauchy_tensor_deri

      subroutine cal_moment_deri(mom,&
                              momdF,&
                          momdS_sum,&
              conat,conatd1,conatd2,&
              conAR,conARd1,conARd2,&
              bcur ,bcurde1,bcurde2,&
              BR   ,BRd1   ,BRd2   ,&
                                 EB,&
                          nlat,nlon,&
                              ibend)
      real mom(2,2,nlat,nlon) !moment mom(beta,alpha,i,j) to M^{\alpha \beta}
      real momdF(2,2,2,nlat,nlon)!first derivative of moment, momdF(\lambda,\beta,\alpha) 
!                                \frac{M^{\alpha \beta}}{\partial \xi^{\lambda}}
      real momdS_sum(nlat,nlon)     
      real conat(2,2,nlat,nlon)
      real conatd1(2,2,2,nlat,nlon)
      real conatd2(2,2,2,2,nlat,nlon)
      real conAR(2,2,nlat,nlon)
      real conARd1(2,2,2,nlat,nlon)
      real conARd2(2,2,2,2,nlat,nlon)
      real bcur(2,2,nlat,nlon)!curvature tensor.
      real bcurde1(2,2,2,nlat,nlon)
      real bcurde2(2,2,2,2,nlat,nlon)
      real BR(2,2,nlat,nlon)
      real BRd1(2,2,2,nlat,nlon)
      real BRd2(2,2,2,2,nlat,nlon)     
      real EB
      integer gam,gamp
      integer al,be

      real a(2,2)
      real ad1(2,2,2)
      real ad2(2,2,2,2)
      real AR(2,2)
      real ARd1(2,2,2)
      real ARd2(2,2,2,2)
      real b(2,2)
      real bd1(2,2,2)
      real bd2(2,2,2,2)
      real BRI(2,2)
      real BRId1(2,2,2)
      real BRId2(2,2,2,2)
      if(ibend.eq.1) then
      mom=0
      momdF=0
      momdS_sum=0
!      if(nid.eq.0) write(*,*) bcur(:,:,:,:)
      do j=1,nlon
      do i=1,nlat
         a   = conat(:,:,i,j)
         ad1 = conatd1(:,:,:,i,j)
         ad2 = conatd2(:,:,:,:,i,j)
         AR  = conAR(:,:,i,j)
         ARd1= conARd1(:,:,:,i,j)
         ARd2= conARd2(:,:,:,:,i,j)
         b   = bcur(:,:,i,j)
         bd1 = bcurde1(:,:,:,i,j)
         bd2 = bcurde2(:,:,:,:,i,j)
         BRI = BR(:,:,i,j)
         BRId1 = BRd1(:,:,:,i,j)
         BRId2 = BRd2(:,:,:,:,i,j)
         do al = 1,2
         do be = 1,2
!c$$$         if(nid.eq.0.and.i.eq.1.and.j.eq.1.and.al.eq.1.and.be.eq.1) then
!c$$$         write(*,*) maxval(b),maxval(bd1),maxval(bd2)
!c$$$         write(*,*) maxval(BRI),maxval(BRId1),maxval(BRId2)
!c$$$         end if
         !!!!!!!!!!!!!!!to get the sum of second derivative!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do gam=1,2
            do gamp=1,2
            momdS_sum(i,j)=momdS_sum(i,j)+&
           ad2(be,al,be,gam)*(AR(gamp,al)*BRI(gam,gamp)&
           -a(gamp,al)*b(gam,gamp))+&
           ad1(al,be,gam)*&
           (ARd1(be,gamp,al)*BRI(gam,gamp)+&
           AR(gamp,al)*BRId1(be,gam,gamp)-&
           ad1(be,gamp,al)*b(gam,gamp)-&
           a(gamp,al)*bd1(be,gam,gamp))&
           +ad1(be,be,gam)*&
           (ARd1(al,gamp,al)*BRI(gam,gamp)+&
           AR(gamp,al)*BRId1(al,gam,gamp)-&
           ad1(al,gamp,al)*b(gam,gamp)-&
           a(gamp,al)*bd1(al,gam,gamp))+a(be,gam)*(&
           ARd2(be,al,gamp,al)*BRI(gam,gamp)+&
           ARd1(al,gamp,al)*BRId1(be,gam,gamp)+&
           ARd1(be,gamp,al)*BRId1(al,gam,gamp)+&
           AR(gamp,al)*BRId2(be,al,gam,gamp)-&
           ad2(be,al,gamp,al)*b(gam,gamp)-&
           ad1(al,gamp,al)*bd1(be,gam,gamp)-&
           ad1(be,gamp,al)*bd1(al,gam,gamp)-&
           a(gamp,al)*bd2(be,al,gam,gamp))
            end do
            end do
         !!!!!!!!!!!!!!!to get the sum of second derivative!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
       
         !!!!!!!!!!!!!!!to get the moment!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if(al.eq.2.and.be.eq.1) CYCLE
               do gam=1,2
               do gamp=1,2
               mom(be,al,i,j)=mom(be,al,i,j)+a(be,gam)*&
              (AR(gamp,al)*BRI  (gam,gamp)-&
               a(gamp,al)*b(gam,gamp))   
               end do
               end do
         !!!!!!!!!!!!!!!to get the moment!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!!!!!!!!!!!!!to get the first derivative of moment!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do l=1,2
         do gam=1,2
         do gamp=1,2
            momdF(l,be,al,i,j)=momdF(l,be,al,i,j)+&
           ad1(l,be,gam)*(AR(gamp,al)*BRI(gam,gamp)-&
           a(gamp,al)*b(gam,gamp))+a(be,gam)*&
           (ARd1(l,gamp,al)*BRI(gam,gamp)&
           +AR(gamp,al)*BRId1(l,gam,gamp)&
           -ad1(l,gamp,al)*b(gam,gamp)&
           -a(gamp,al)*bd1(l,gam,gamp))
         end do
         end do
         end do
         !!!!!!!!!!!!!!!to get the first derivative of moment!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         end do
         end do  

      end do
      end do
      mom(1,2,:,:) = mom(2,1,:,:)
      momdF(:,1,2,:,:) = momdF(:,2,1,:,:)

      mom = mom*EB
      momdF = momdF*EB  
      momdS_sum = momdS_sum*EB
      end if

      return
      end subroutine cal_moment_deri


      subroutine cal_Q(Q,&
                  Qcode,&
                    mom,&
                  momdF,&
              momdS_sum,&
               christof,&
              chrisderi,&
              nlat,nlon,&
                  ibend)

      real Q(2,nlat,nlon)
      real Qcode(nlat,nlon)
      real mom(2,2,nlat,nlon) 
      real momdF(2,2,2,nlat,nlon)
      real momdS_sum(nlat,nlon)
      real christof(2,2,2,nlat,nlon)
      real chrisderi(2,2,2,2,nlat,nlon)
      integer al,be,lm
!!!!!!Qcode means $Q^{\beta}|_{\beta}$
      if(ibend.eq.1) then
      Q=0
      Qcode=momdS_sum
      do j=1,nlon
      do i=1,nlat

      do be = 1,2
         do al = 1,2
         Q(be,i,j)=Q(be,i,j)+momdF(al,be,al,i,j)
         do lm = 1,2
            Q(be,i,j)=Q(be,i,j)+christof(lm,al,al,i,j)*mom(be,lm,i,j)&
                              +christof(lm,be,al,i,j)*mom(lm,al,i,j)
      Qcode(i,j)=Qcode(i,j)+chrisderi(be,lm,al,al,i,j)*mom(be,lm,i,j)&
     +christof(lm,al,al,i,j)*momdF(be,be,lm,i,j)&
     +chrisderi(be,lm,be,al,i,j)*mom(lm,al,i,j)&
     +christof(lm,be,al,i,j)*momdF(be,lm,al,i,j)
         end do
         end do
      end do

      do be=1,2
      do lm=1,2
      Qcode(i,j)=Qcode(i,j)+Q(lm,i,j)*christof(lm,be,be,i,j)
      end do
      end do
      
      end do
      end do
      end if
      return
      end subroutine cal_Q


      subroutine cal_bcurmix(bcurmix,&
                               bcur,&
                              conat,&
                          nlat,nlon,&
                              ibend)

      real bcurmix(2,2,nlat,nlon)
      real bcur(2,2,nlat,nlon)
      real conat(2,2,nlat,nlon)
      integer al,be,lm
      if(ibend.eq.1) then
      bcurmix = 0
      do j=1,nlon
      do i=1,nlat
         do al=1,2
         do be=1,2
            do lm = 1,2
            bcurmix(be,al,i,j)=bcurmix(be,al,i,j)+&
           conat(lm,al,i,j)*bcur(be,lm,i,j)  
            end do
         end do
         end do
      end do
      end do
      end if

      return
      end subroutine cal_bcurmix

      subroutine cal_bcurcon(bcurcon,&
                               bcur,&
                              conat,&
                          nlat,nlon,&
                              ibend)

      real bcurcon(2,2,nlat,nlon)
      real bcur(2,2,nlat,nlon)
      real conat(2,2,nlat,nlon)
      integer al,be,lm,ga
      if(ibend.eq.1) then
      bcurcon = 0
      do j=1,nlon
      do i=1,nlat
         do al=1,2
         do be=1,2
            do lm = 1,2
            do ga = 1,2
            bcurcon(be,al,i,j)=bcurcon(be,al,i,j)+&
           conat(lm,al,i,j)*conat(ga,be,i,j)*&
           bcur(lm,ga,i,j)  
            end do
            end do
         end do
         end do
      end do
      end do
      end if

      return
      end subroutine cal_bcurcon




      subroutine cal_bcur(bcur,&
                            bcurmix,&
                              cat,&
                          nlat,nlon,&
                              ibend)

      real bcurmix(2,2,nlat,nlon)
      real bcur(2,2,nlat,nlon)
      real cat(2,2,nlat,nlon)
      integer al,be,lm
      if(ibend.eq.1) then
      bcur = 0
      do j=1,nlon
      do i=1,nlat
         do al=1,2
         do be=1,2
            do lm = 1,2
            bcur(be,al,i,j)=bcur(be,al,i,j)+&
           cat(lm,al,i,j)*bcurmix(be,lm,i,j)  
            end do
         end do
         end do
      end do
      end do
      end if

      return
      end subroutine cal_bcur




      subroutine cal_spectral_force(x,y,z,fx,fy,fz)
      real a(nnlat2,nnlat2,3)
      real b(nnlat2,nnlat2,3)
      real x(nnlat2,nnlon2), y(nnlat2,nnlon2), z(nnlat2,nnlon2)
      real fx(nnlat2,nnlon2),fy(nnlat2,nnlon2),fz(nnlat2,nnlon2)  
      real cat(2,2,nnlat,nnlon)
      real conat(2,2,nnlat,nnlon)
      real catd1(2,2,2,nnlat,nnlon)
      real conatd1(2,2,2,nnlat,nnlon)
      real conatd2(2,2,2,2,nnlat,nnlon)
      real cvat1(3,nnlat,nnlon)
      real cvat2(3,nnlat,nnlon)
      real cvan3(3,nnlat,nnlon)
      real surfmetN(nnlat,nnlon)
      real christof(2,2,2,nnlat,nnlon)
      real chrisderi(2,2,2,2,nnlat,nnlon)
      real bcur(2,2,nnlat,nnlon)!curvature tensor.
      real Hcur(2,2,nnlat,nnlon)!H=0.5*g^{\alpha\beta}b_{\alpha\beta}!!see pp581 of zhao's paper 'The dynamics of a vesicle in simple shear flow'.
      real bcurde1(2,2,2,nnlat,nnlon)
      real bcurde2(2,2,2,2,nnlat,nnlon)
      real bcurmix(2,2,nnlat,nnlon)!bcurmix(beta,alpha,:,:) to b_{beta}^{alpha}
!      real bcurcon(2,2,nnlat,nnlon)

      real chf(2,2,2)
      real Js, I1(nnlat,nnlon), I2(nnlat,nnlon)
      real Jsde(2),I1de(2),I2de(2)
      real pwspI1(nnlat,nnlon),pwspI2(nnlat,nnlon)
      real pwspI1de(2),pwspI2de(2)
      real T11,T22,T12
      real Td11(2),Td22(2),Td12(2)
      real q1,q2,q3
      real cdeSvat(3,4,nnlat,nnlon) ! S: second derivative
      real cdeTvat(3,5,nnlat,nnlon) ! T: third derivative 
      real mom(2,2,nnlat,nnlon) !moment mom(beta,alpha,i,j) to M^{\alpha \beta}
      real momdF(2,2,2,nnlat,nnlon)!first derivative of moment, momdF(\lambda,\beta,\alpha) 
!                                \frac{M^{\alpha \beta}}{\partial \xi^{\lambda}}
      real momdS_sum(nnlat,nnlon)!second derivative of moment, momdS(l,\lambda,\beta,\alpha)
!                                \frac{\partial^{2} M^{\alpha \beta}}{\partial \xi^{\alpha} \partial \xi^{\beta}}     ! summer over \alpha and \beta
      real Q(2,nnlat,nnlon) !surface transverse tensor, see poz book, Zhao, page 89
      real Qcode(nnlat,nnlon) ! Q^{\alpha}_{| \alpha} sum over \alpha, see page 89, 3.4.11 covarient derivative
      real fxN(nnlat,nnlon)
      real fyN(nnlat,nnlon)
      real fzN(nnlat,nnlon)  
      real Eshr(nnlat,nnlon),Edil(nnlat,nnlon),Eben(nnlat,nnlon)
      
      !!!!!!singular integration!!!!!!!!!!
!      real vxtsig(ntmppsi),vytsig(ntmppsi),vztsig(ntmppsi)
!      integer,save :: icald; data icald /0/
      !!!!!!singular integration!!!!!!!!!!

!c$$$      if(ifsing.eq.1.and.icald.eq.0) then
!c$$$         call singular_allocate2
!c$$$         icald = 1
!c$$$      end if

      call cal_3dcrd_spectral_coeff(a,b,x,y,z,nnlat2,nnlon2)
      acrdD = a
      bcrdD = b




      call cal_co_contra_metric_deri(cvat1,cvat2,cvan3,&
                                             surfmetN,&
                                                  cat,&
                                                catd1,&
                                                conat,&
                                      conatd1,conatd2,&
                                      cdeSvat,cdeTvat,&
                                             christof,&
                                            chrisderi,&
                                                 bcur,&
                                      bcurde1,bcurde2,&
                                                  a,b,&
                                          thetaN,phiN,&
                                        nnlat2,nnlon2,&
                                          nnlat,nnlon,&
                                                ibend)

         

!******************************************************************unitN is normal vector*****************************************************************************
!!!!!!!!!!!STORE THE UNIT OUT NORMAL VECTOR IN COMMON DATA!!!!!!!
      unitN = cvan3
!!!!!!!!!!!STORE THE UNIT OUT NORMAL VECTOR IN COMMON DATA!!!!!!!
!      if(nid.eq.0) write(*,*) cvat1
!******************************************************************computes mean curvature from curvature tensor******************************************************
             
      call cal_mean_curvature(Hcur,conat,bcur,nnlat,nnlon)

      if(idealias.eq.1) then
         call cal_unitO_filter(unitO,&
             unitN,&
             thetaO,phiO,&
             nnlat2,nnlon2,&
             nnlat,nnlon)



         call cal_surface_metric_sph__OPT(surfmetD,&
                             a,b,thetaO,phiO,&
                               nnlat2,nnlon2)
      else
      do j = 1,nnlon2
      do i = 1,nnlat2
         unitO(:,i,j)  = unitN(:,i,j) 
         surfmetD(i,j) = surfmetN(i,j)




      end do
      end do
      end if
      wgtD = surfmetD*wgtdsin
      call cal_moment_deri(mom,&
                        momdF,&
                    momdS_sum,&
        conat,conatd1,conatd2,&
        conAR,conARd1,conARd2,&
        bcur ,bcurde1,bcurde2,&
        BR   ,BRd1   ,BRd2   ,&
                        EBEND,&
                  nnlat,nnlon,&
                        ibend)


      call cal_Q(Q,&
            Qcode,&
              mom,&
            momdF,&
        momdS_sum,&
         christof,&
        chrisderi,&
      nnlat,nnlon,&
            ibend)

      call cal_bcurmix(bcurmix,&
      bcur,&
      conat,&
      nnlat,nnlon,&
      ibend)

      

      do j = 1,nnlon
      do i = 1,nnlat
      call cal_Js_I1_I2_deri(Js,I1(i,j),I2(i,j),&
     Jsde,I1de,I2de,&
     conAR(:,:,i,j),&
     cat(:,:,i,j),&
     conARd1(:,:,:,i,j),&
     catd1(:,:,:,i,j))



      call cal_pwspI(pwspI1(i,j),pwspI2(i,j),&
     pwspI1de,pwspI2de,&
     I1(i,j),I2(i,j),&
     I1de,I2de,&
     GsCELL,CSK)


      call calculate_cauchy_tensor_deri(T11,T22,T12,&
     Td11,Td22,Td12,&
     Js,pwspI1(i,j),pwspI2(i,j),&
     Jsde,pwspI1de,pwspI2de,&
     conAR(:,:,i,j),&
     conARd1(:,:,:,i,j),&
     conat(:,:,i,j),&
     conatd1(:,:,:,i,j))

      chf = christof(:,:,:,i,j)
      q1  = Td11(1)+Td12(2)+&
     chf(1,1,1)*T11+chf(2,1,1)*T12+chf(1,2,2)*T11+chf(2,2,2)*T12+&
     chf(1,1,1)*T11+chf(2,1,1)*T12+chf(1,1,2)*T12+chf(2,1,2)*T22
      q2  = Td12(1)+Td22(2)+&
     chf(1,1,1)*T12+chf(2,1,1)*T22+chf(1,2,2)*T12+chf(2,2,2)*T22+&
     chf(1,2,1)*T11+chf(2,2,1)*T12+chf(1,2,2)*T12+chf(2,2,2)*T22

      q3  = T11*bcur(1,1,i,j)+T22*bcur(2,2,i,j)+2*T12*bcur(1,2,i,j) 

      if(ibend.eq.1) then
      q1=q1-bcurmix(1,1,i,j)*Q(1,i,j)-bcurmix(2,1,i,j)*Q(2,i,j)
      q2=q2-bcurmix(1,2,i,j)*Q(1,i,j)-bcurmix(2,2,i,j)*Q(2,i,j) 
      q3=q3+Qcode(i,j)
      end if

      fxN(i,j) = q1*cvat1(1,i,j)+q2*cvat2(1,i,j)+q3*cvan3(1,i,j)
      fyN(i,j) = q1*cvat1(2,i,j)+q2*cvat2(2,i,j)+q3*cvan3(2,i,j)
      fzN(i,j) = q1*cvat1(3,i,j)+q2*cvat2(3,i,j)+q3*cvan3(3,i,j)
      end do
      end do

!!!!!!!!!!!!!!!calculate deforming energy!!!!!!!!!!!!!!!!!!!!!!!!!
      if(ifenergy.eq.1) then
         call cal_deform_energy(Eshr,Edil,Eben,I1,I2,Hcur,&
     GsCELL,CSK,EBEND,nnlat,nnlon,ibend)
      end if
!!!!!!!!!!!!!!!calculate deforming energy!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!calculate principal stress on higher modes nnlat*nnlon
      if(ifprinstress.eq.1) then
         call cal_prinstress(prstsmax,prstsmin,I1,I2,pwspI1,pwspI2,&
                            GsCELL,nnlat,nnlon)
      end if
!!!!!!calculate principal stress on higher modes nnlat*nnlon

      if(idealias.eq.1) then
            call filter_modes4__OPT(fx,fy,fz,&
                              prstsmax2,prstsmin2,&
                              E_shr,E_dil,E_ben,&
                              prstsmaxpk,prstsminpk,&
                              E_shrpk,E_dilpk,E_benpk,&
                              fxN,fyN,fzN,&
                              prstsmax,prstsmin,&
                              Eshr,Edil,Eben,&
                              thetaO,phiO,&
                              nnlat2,nnlon2,&
                              nnlat,nnlon,&
                              ifprinstress,ifenergy)
      else
         fx = fxN(1:nnlat2,1:nnlon2)  
         fy = fyN(1:nnlat2,1:nnlon2)  
         fz = fzN(1:nnlat2,1:nnlon2)  
      end if

!      write(*,*) 'For nid',nid, ' is: ', maxval(fx),minval(fx)
!      write(*,*) 'For nid',nid, ' is: ', maxval(fy),minval(fy)
!      write(*,*) 'For nid',nid, ' is: ', maxval(fz),minval(fz)

      return
      end subroutine cal_spectral_force

      subroutine calforcegram_spectral
      call cal_spectral_force(crdxD,crdyD,crdzD,fxL,fyL,fzL)
      fxL = fxL*wgtD
      fyL = fyL*wgtD
      fzL = fzL*wgtD
      return
      end subroutine calforcegram_spectral

      subroutine sum_loc_glb_c_spectral2(uxtt,uytt,uztt,&
                                        uxlc,uylc,uzlc,&
                                        uxgl,uygl,uzgl,&
                                                   ntt)
      real uxtt(ntt),uytt(ntt),uztt(ntt)
      real uxlc(ntt),uylc(ntt),uzlc(ntt)
      real uxgl(ntt),uygl(ntt),uzgl(ntt)
      uxtt = uxlc + uxgl
      uytt = uylc + uygl
      uztt = uzlc + uzgl      
      return
      end subroutine sum_loc_glb_c_spectral2


      subroutine updatecell_spectral(dt,it)
      integer ittag
      if(ifpop_cellvel.eq.123456) then
         ittag = it + 3 !!if restart is activated, alwasy go to default case
      else
         ittag = it
      end if

      select case(ittag)
      case (1)
         uxmul3 = uxtot;uymul3 = uytot;uzmul3 = uztot
         CRDXD=CRDXD+DT*UXMUL3
         CRDYD=CRDYD+DT*UYMUL3
         CRDZD=CRDZD+DT*UZMUL3
      case (2)
         uxmul2 = uxmul3;uymul2 = uymul3;uzmul2 = uzmul3
         uxmul3 = uxtot; uymul3 = uytot; uzmul3 = uztot
         CRDXD=CRDXD+DT*(1.5*UXMUL3-0.5*UXMUL2)
         CRDYD=CRDYD+DT*(1.5*UYMUL3-0.5*UYMUL2)
         CRDZD=CRDZD+DT*(1.5*UZMUL3-0.5*UZMUL2)
      case default
         uxmul1 = uxmul2;uymul1 = uymul2;uzmul1 = uzmul2
         uxmul2 = uxmul3;uymul2 = uymul3;uzmul2 = uzmul3
         uxmul3 = uxtot; uymul3 = uytot; uzmul3 = uztot
         CRDXD=CRDXD+DT*(23./12.*UXMUL3-4./3.*UXMUL2&
                  +5./12.*UXMUL1)
         CRDYD=CRDYD+DT*(23./12.*UYMUL3-4./3.*UYMUL2&
                  +5./12.*UYMUL1)
         CRDZD=CRDZD+DT*(23./12.*UZMUL3-4./3.*UZMUL2&
                  +5./12.*UZMUL1)
      end select

!      print * , nid,'sum x vel',sum(uxmul1),sum(uxmul2),sum(uxmul3)
!      print * , nid,'sum y vel',sum(uymul1),sum(uymul2),sum(uymul3)
!      print * , nid,'sum z vel',sum(uzmul1),sum(uzmul2),sum(uzmul3)
!      call calvolcor()
     
!      write(*,'(a,i5,f10.5)') 'max_X',nid,maxval(X3)
!      write (*,*) maxval(uy3g),minval(uy3g)
!      write (*,*) maxval(z3),minval(z3)
      end subroutine updatecell_spectral









!      function inid_tar(i,j,limit)
!      inid_tar = i + j
!      if(inid_tar.lt.1)inid_tar=inid_tar+limit
!      if(inid_tar.gt.limit)inid_tar=inid_tar-limit
!      inid_tar = nidcl(inid_tar)
!      return
!      end function inid_tar




      subroutine cal_sin_normal(rn1,rn2,rn3,x,y,z,a,b,c)
      common /myconstants/ pimy
      ax=2*x/a**2
      ay=2*y/b**2
      sz=1-0.7*cos(pimy/2.*(0.2+z/c))
      dsz=0.35*pimy/c*sin(pimy/2.*(0.2+z/c))
      az=2*z/c**2*sz**2+(z**2/c**2-1)*2*sz*dsz
      rnorm = sqrt(ax**2+ay**2+az**2)
      rn1 = ax/rnorm
      rn2 = ay/rnorm
      rn3 = az/rnorm
      end subroutine cal_sin_normal



!      subroutine cal_crd_triexpt
!      real xtmp(nn3,3)
!      call cal_sph_3dfun_vector(xtmp,thetri,phitri,acrdD,bcrdD,&
!     pbartri,nn3,nnlat2,nnlon2)
!      x3=xtmp(:,1);y3=xtmp(:,2);z3=xtmp(:,3)
!!      print *,x3-(x3_eq+ccen(1)),y3-(y3_eq+ccen(2)),
!!     $z3-(z3_eq+ccen(3))
!      end subroutine cal_crd_triexpt



!!!!!!!!!this one works for any models according skalak
!!!!!!!!!STRAIN ENERGY FUNCTION OF RED BLOOD CELL MEMBRANES
      subroutine cal_prinstress(prinmax,prinmin,I1,I2,pwpI1,pwpI2,&
      gs,nlat,nlon)
      real prinmax(nlat,nlon),prinmin(nlat,nlon)
      real I1(nlat,nlon),I2(nlat,nlon)
      real pwpI1(nlat,nlon),pwpI2(nlat,nlon)
      real lam1(nlat,nlon),lam2(nlat,nlon)
      real root(nlat,nlon)
      root = sqrt(I1**2+4*I1-4*I2)
      lam1 = sqrt(0.5*(2+I1+root))
      lam2 = sqrt(0.5*(2+I1-root))
      prinmax=2*lam1/lam2*(pwpI1+lam2*lam2*pwpI2)
      prinmin=2*lam2/lam1*(pwpI1+lam1*lam1*pwpI2)
      prinmax = prinmax / gs
      prinmin = prinmin / gs
      end subroutine cal_prinstress
!!!!!!!!!this one works for any models according skalak

      
      subroutine cal_rad_rbc(theta,x_rad,y_rad,z_rad)
      real, intent(in)::theta
      real,intent(out)::x_rad,y_rad,z_rad
      integer::k
      real::rd,rdO,f,df,rb
      rd=.02
      rdO=1.
      k=0
      do while ((abs(rd-rdO)/abs(rdO))>.00001)
      k=k+1
      f=(rd**4.)+(2.*((arbc*rd)**2.)*cos(2.*theta))+(arbc**4.)-(crbc**4.)
      df=(4.*(rd**3.))+(4.*rd*(arbc**2.)*cos(2.*theta))
      rdO=rd
      rd=rd-(f/df)
!      rb=rd*sin(theta)
!      f=(.5*sqrt(1.-(rb**2.))*(c0+(c2*(rb**2.))+(c4*(rb**4.))))-(rd*cos(theta))              
!      df=(-.5*rd*(sin(theta)**2.)*((1.-(rb**2.))**-.5)*(c0+(c2*(rb**2.))+(c4*(rb**4.))))+&
!      ((.5*sqrt(1.-(rb**2.)))*((2.*c2*(sin(theta)**2.)*rd)+(4.*c4*(sin(theta)**4.)*(rd**3.))))-cos(theta)
!      rdO=rd
!      rd=rd-(f/df)
      enddo
      x_rad=sqrt((((arbc**2.)+(rd**2.))**2.)-(crbc**4.))/(2.*arbc)
      y_rad=x_rad
      z_rad=sqrt(-(((arbc**2.)-(rd**2.))**2.)+(crbc**4.))/(2.*arbc)    
      return
      end subroutine cal_rad_rbc
      end module mod_sph


