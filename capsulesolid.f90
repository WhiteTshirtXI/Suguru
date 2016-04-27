      module mod_capsulesolid
      use mod_common
      use mod_param
      use mod_math
      public 
      contains

      
      subroutine calvol_area(vol,area,xnode,ynode,znode,triarea,nout)
      real crdcen(3)
      real xnode(nn3),ynode(nn3),znode(nn3)
      real triarea(ne3),nout(ne3,3)
      real vol,area
      vol = 0
      area = 0
      do ie = 1,ne3
         crdcen(1)=1.0/3*(xnode(in3(ie,1))+xnode(in3(ie,2))+&
                         xnode(in3(ie,3)))
         crdcen(2)=1.0/3*(ynode(in3(ie,1))+ynode(in3(ie,2))+&
                         ynode(in3(ie,3)))
         crdcen(3)=1.0/3*(znode(in3(ie,1))+znode(in3(ie,2))+&
                         znode(in3(ie,3)))
         vol=vol+1.0/3*(nout(ie,1)*crdcen(1)+&
                       nout(ie,2)*crdcen(2)+&
                       nout(ie,3)*crdcen(3))*triarea(ie)
         area=area+triarea(ie)
      end do
      return
      end subroutine calvol_area

      subroutine calavenodenormal!calculate average normal vector at the deformed status
      real vecnorm
      integer indelem
      nnvd = 0
      do i = 1,nn3
         naed(i) = 0
         do ie = 1,ie3(i,1)
            indelem = ie3(i,ie+1)
            nnvd(i,:)=nnvd(i,:)+aed(indelem)*nvd(indelem,:)
            naed(i)=naed(i)+aed(indelem)
         end do
         nnvd(i,:)=nnvd(i,:)/naed(i)
         vecnorm = rnorm2(nnvd(i,1),nnvd(i,2),nnvd(i,3))
         nnvd(i,:)=nnvd(i,:)/vecnorm
      end do
      
      naed = naed/3.0
      return
      end subroutine calavenodenormal

      function NtriF(eta1,eta2)!F means first order
      real eta1,eta2
      real NtriF(3)
      NtriF(1) = 1 - eta1 - eta2
      NtriF(2) = eta1
      NtriF(3) = eta2
      return
      end function NtriF

      function NDtriF(alpha)!\frac{d N^{1}}{d \eta^{\alpha}}
                            !\alpha = 1 or 2, here in both cases, value is -1
      integer alpha
      real NDtriF(3)
      NDtriF(1) = -1
      NDtriF(2) = 2 - alpha
      NDtriF(3) = alpha - 1      
      return
      end function NDtriF

      subroutine cobase(vec_a,x1,y1,z1,& !co: covarient,  base: base vector \vec{a}
                             x2,y2,z2,&
                             x3,y3,z3) !calculate local covarient base \vec{a}_{\alpha}
      !! vec_a is a array vec_a(2,3)
      real vec_a(2,3) ! 3 means three components x, y, z 2 means alpha = 1, 2
      integer alpha
      real val(3)
      do alpha = 1, 2
      val = NDtriF(alpha)
      vec_a(alpha,1) = val(1)*x1 &
                    + val(2)*x2 &
                    + val(3)*x3 

      vec_a(alpha,2) = val(1)*y1&
                    + val(2)*y2&
                    + val(3)*y3 

      vec_a(alpha,3) = val(1)*z1&
                    + val(2)*z2&
                    + val(3)*z3 
      end do
      return
      end subroutine cobase

      function hybprod(a1x,a1y,a1z,&
                      a2x,a2y,a2z,&
                      a3x,a3y,a3z)
      real hybprod
      real atmpx,atmpy,atmpz,temp
      atmpx = a2y*a3z-a3y*a2z
      atmpy = a2z*a3x-a3z*a2x
      atmpz = a2x*a3y-a3x*a2y
      hybprod = a1x*atmpx+a1y*atmpy+a1z*atmpz
      return
      end function hybprod

      subroutine co2contra(co1x,co1y,co1z,&
                          co2x,co2y,co2z,&
                          co3x,co3y,co3z,&
                          tr1x,tr1y,tr1z,&
                          tr2x,tr2y,tr2z,&
                          tr3x,tr3y,tr3z)
      real atmpx,atmpy,atmpz,hybtmp
      hybtmp = hybprod(co1x,co1y,co1z,&
                      co2x,co2y,co2z,&
                      co3x,co3y,co3z)

      atmpx = co2y*co3z-co3y*co2z
      atmpy = co2z*co3x-co3z*co2x
      atmpz = co2x*co3y-co3x*co2y
      tr1x = atmpx / hybtmp
      tr1y = atmpy / hybtmp
      tr1z = atmpz / hybtmp     


      hybtmp = hybprod(co2x,co2y,co2z,&
                      co3x,co3y,co3z,&
                      co1x,co1y,co1z)

      atmpx = co3y*co1z-co1y*co3z
      atmpy = co3z*co1x-co1z*co3x
      atmpz = co3x*co1y-co1x*co3y
      tr2x = atmpx / hybtmp
      tr2y = atmpy / hybtmp
      tr2z = atmpz / hybtmp     

      hybtmp = hybprod(co3x,co3y,co3z,&
                      co1x,co1y,co1z,&
                      co2x,co2y,co2z)
      atmpx = co1y*co2z-co2y*co1z
      atmpy = co1z*co2x-co2z*co1x
      atmpz = co1x*co2y-co2x*co1y
      tr3x = atmpx / hybtmp
      tr3y = atmpy / hybtmp
      tr3z = atmpz / hybtmp     
      return
      end subroutine co2contra


      subroutine calmetrtensor(mett,&
                              a1x,a1y,a1z,&
                              a2x,a2y,a2z)
      real mett(4) 
!     can be used for covarient and contravarient case
!     also   used for reference status and deformed status
!     mett(1:4) stores a_{\alpha \beta} or a^{\alpha \beta}
!     \alpha = 1 \beta = 1 mett(1)
!     \alpha = 2 \beta = 2 mett(2) 
!     \alpha = 1 \beta = 2 mett(3)
!     mett(4) = mett(1)*mett(2) - mett(3)**2
      mett(1) = a1x**2  + a1y**2  + a1z**2
      mett(2) = a2x**2  + a2y**2  + a2z**2
      mett(3) = a1x*a2x + a1y*a2y + a1z*a2z
      mett(4) = mett(1)*mett(2) - mett(3)**2
      return
      end subroutine calmetrtensor

      subroutine pWspI_NH(pWspI,GsCELL,CSK,I1,I2)
      real I1,I2,GsCELL,CSK
      real pWspI(2)
      pWspI(1) = 0.5*GsCELL
      pWspI(2) = -0.5*GsCELL*(I2+1)**(-2)    
!      write(*,'(5f12.8)') pWspI(2),GsCELL,CSK,I1,I2
      return
      end subroutine pWspI_NH 

      subroutine pWspI_SK(pWspI,GsCELL,CSK,I1,I2)
      real I1,I2,GsCELL,CSK
      real pWspI(2)
      pWspI(1) = 0.5*GsCELL*(1+I1)
      pWspI(2) = 0.5*GsCELL*(CSK*I2-1)    
      return
      end subroutine pWspI_SK  


      subroutine calI_Js(UAT,lac,I1,I2,Js) ! calculate I1, I2, element wise
      real UAT(4),lac(4) 
      real I1,I2,Js
      ! UAT upper case A contravariant
      ! lac lower case a covariant&
      I1 = UAT(1)*lac(1)+UAT(2)*lac(2)+&
          2*UAT(3)*lac(3) - 2
      I2 = UAT(4)*lac(4)-1
      Js = sqrt(1+I2) 
!      write (*,*) Js
      return
      end subroutine calI_Js

      subroutine calT_NH(T,UAT,lat,lac,GsCELL,CSK) ! calculate Cauthy tensor T element wise   
      real I1,I2,Js,GsCELL,CSK,T(3),pWspI(2)
      real UAT(4),lat(4),lac(4)
      call calI_Js(UAT,lac,I1,I2,Js)
!      call pWspI_NH(pWspI,GsCELL,CSK,I1,I2)
      call calpwspi(pWspI,GsCELL,CSK,I1,I2)
      do i = 1,3
         T(i) = 2.0/Js*pWspI(1)*UAT(i) +&
               2.0*Js*pWspI(2)*lat(i)
      end do
      return
      end subroutine calT_NH
      
      subroutine calchi(chi,vec_a) ! defined in Eq.2.5.14, only the first two terms are non-zero for first order elements. Element wise
      real chi(3,3,3)
      ! first  3, 3 vertices on the element
      ! second 3, 3 components x, y, z
      ! third  3, 1,\alpha=1 \beta=1 
      !           2,\alpha=2 \beta=2
      !           3,\alpha=1 \beta=2
      integer p,j
      real Npb(2,3),vec_a(2,3) ! 2 means alpha and beta
                               ! 3 means three vertices
      do i = 1,2
         Npb(i,:) = NDtriF(i)
      end do

      do j = 1,3
         do p = 1,3
            chi(p,j,1) = Npb(1,p)*vec_a(1,j)
            chi(p,j,2) = Npb(2,p)*vec_a(2,j)
            chi(p,j,3) = 0.5*Npb(2,p)*vec_a(1,j)+&
                        0.5*Npb(1,p)*vec_a(2,j)
          end do
      end do
      return
      end subroutine calchi

      subroutine MulchiT(chiMT,chi,T)!product of chi and T element wise
      real chiMT(3,3)
      real chi(3,3,3),T(3)
!     chiMT(3,3) first 3: 3 vertices, second 3: 3 components
      chiMT=chi(:,:,1)*T(1)+&
           chi(:,:,2)*T(2)+&
         2*chi(:,:,3)*T(3)   
      return
      end subroutine MulchiT

      subroutine calchiMT_NH(chiMT,vec_a,UAT,lat,lac,GsCELL,CSK)
      real chiMT(3,3)
      real vec_a(2,3),UAT(4),lat(4),lac(4),GsCELL,CSK
      real T(3),chi(3,3,3)
      call calT_NH(T,UAT,lat,lac,GsCELL,CSK)
      call calchi(chi,vec_a)
      call MulchiT(chiMT,chi,T)
      return
      end subroutine calchiMT_NH
      
      subroutine metric_ref()!reference status decide UAT
      integer nind(3) !global node index of the three vertices on the element.
      real baco(2,3),mett(4),batra(3,3),mett2(4)
      real prod(3,3)
      integer k
      real tmpa 
      tmpa= 0.
      do ie = 1, NE3
          prod = 0
          nind(:) = in3(ie,1:3)
          call cobase(baco,x3_eq(nind(1)),y3_eq(nind(1)),z3_eq(nind(1)),&
                          x3_eq(nind(2)),y3_eq(nind(2)),z3_eq(nind(2)),&
                          x3_eq(nind(3)),y3_eq(nind(3)),z3_eq(nind(3)))

          call calmetrtensor(mett2,&
                           baco(1,1),baco(1,2),baco(1,3),&
                           baco(2,1),baco(2,2),baco(2,3))
          tmpa=tmpa+0.5*sqrt(mett2(4))
!          baco(3,:) = nvr(ie,:)
          call co2contra(baco(1,1), baco(1,2), baco(1,3),&
                        baco(2,1), baco(2,2), baco(2,3),&
                        nvr(ie,1), nvr(ie,2), nvr(ie,3),&
                        batra(1,1),batra(1,2),batra(1,3),&
                        batra(2,1),batra(2,2),batra(2,3),&
                        batra(3,1),batra(3,2),batra(3,3))

          call calmetrtensor(mett,&
                           batra(1,1),batra(1,2),batra(1,3),&
                           batra(2,1),batra(2,2),batra(2,3))
          uat(ie,:) = mett(1:4)

      end do
!      write (*,'(a,f10.5)') 'area refffffffffffffffffff', tmpa
      return
      end subroutine metric_ref      



      subroutine ReadTriGaussConstants()
      if (NTGL.eq.1) then
         xiT(1) = 1.0/3
         etaT(1) = 1.0/3
         wT(1) = 1.0
      return
      else if (NTGL.eq.3) then
         xiT(1) = 1.0/6
         etaT(1) = 1.0/6
         wT(1) = 1.0/3
         xiT(2) = 2.0/3
         etaT(2) = 1.0/6
         wT(2) = wT(0)
         xiT(3) = 1.0/6
         etaT(3) = 2.0/3
         wT(3) = wT(0)
         return
      else if (NTGL.eq.4) then
         xiT(1) = 1.0/3
         etaT(1) = 1.0/3
         wT(1) = -27.0/48
         xiT(2) = 1.0/5
         etaT(2) = 1.0/5
         wT(2) = 25.0/48
         xiT(3) = 3.0/5
         etaT(3) = 1.0/5
         wT(3) = 25.0/48
         xiT(4) = 1.0/5
         etaT(4) = 3.0/5
         wT(4) = 25.0/48
         return
      else if (NTGL.eq.9) then
         xiT(1)  = 0.437525248383384
         xiT(2)  = 0.124949503233232
         xiT(3)  = 0.437525248383384
         xiT(4)  = 0.165409927389841
         xiT(5)  = 0.037477420750088     
         xiT(6)  = 0.797112651860071
         xiT(7)  = 0.165409927389841
         xiT(8)  = 0.037477420750088
         xiT(9)  = 0.797112651860071

         etaT(1)  = 0.437525248383384
         etaT(2)  = 0.437525248383384
         etaT(3)  = 0.124949503233232
         etaT(4)  = 0.037477420750088
         etaT(5)  = 0.165409927389841     
         etaT(6)  = 0.165409927389841
         etaT(7)  = 0.797112651860071
         etaT(8)  = 0.797112651860071
         etaT(9)  = 0.037477420750088

         wT(1)  = 0.205950504760887
         wT(2)  = 0.205950504760887
         wT(3)  = 0.205950504760887
         wT(4)  = 0.063691414286223
         wT(5)  = 0.063691414286223     
         wT(6)  = 0.063691414286223
         wT(7)  = 0.063691414286223
         wT(8)  = 0.063691414286223
         wT(9)  = 0.063691414286223
         return
      end if
      return
      end subroutine ReadTriGaussConstants

      end module mod_capsulesolid
