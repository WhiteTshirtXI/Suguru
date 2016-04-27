module mod_surf_vol
use mpi
use mod_common
use mod_common_mpi
use mod_interp_spread
use mod_param
use mod_sph
implicit none
private
public surf_vol
contains
subroutine surf_vol
implicit none
integer:: p,i,j,l,cap
real :: vol,vol2,areatot,x_cap_mean,y_cap_mean,z_cap_mean
real::  x_n_cap(n_capf,m_capf),y_n_cap(n_capf,m_capf),z_n_cap(n_capf,m_capf)
real::  sidee(n_capf,m_capf),sidew(n_capf,m_capf),siden(n_capf,m_capf),sides(n_capf,m_capf),rad_cap(n_capf,m_capf)
real::  cosinus(n_capf,m_capf),radyz(n_capf,m_capf),dr(n_capf,m_capf)
real::  inxx,inxy,inxz,inyy,inyz,inzz,i1,i2,i3,r,q,ip1,ip2,ip3,teta,lp1,lp2,lp3,mass,ic1,ic2,ic3,ma,na,radmax,radmin
do p=1,pmax

if (ap(p)%mslv.gt.0) then


 do i=1,n_capf
 do j=1,m_capf
 cap=((i-1)*m_capf)+j
 x_cap(i,j)=ap(p)%xfp(cap)
 y_cap(i,j)=ap(p)%yfp(cap)
 z_cap(i,j)=ap(p)%zfp(cap)
 x_n_cap(i,j)=ap(p)%xn(cap)
 y_n_cap(i,j)=ap(p)%yn(cap)
 z_n_cap(i,j)=ap(p)%zn(cap)
 area_cap(i,j)=ap(p)%area(cap)
 enddo
 enddo

 !**************************************computing some variables needed for calculation of volume *******************************************
x_cap_mean=0.
y_cap_mean=0.
z_cap_mean=0.
do i=1,n_capf
do j=1,m_capf
x_cap_mean=x_cap_mean+x_cap(i,j)
y_cap_mean=y_cap_mean+y_cap(i,j)
z_cap_mean=z_cap_mean+z_cap(i,j)
enddo
enddo
x_cap_mean=x_cap_mean/(m_capf*n_capf)
y_cap_mean=y_cap_mean/(m_capf*n_capf)
z_cap_mean=z_cap_mean/(m_capf*n_capf)


do i=1,n_capf
do j=1,m_capf
rad_cap(i,j)=sqrt(((x_cap(i,j)-x_cap_mean)**2.)+((y_cap(i,j)-y_cap_mean)**2.)+((z_cap(i,j)-z_cap_mean)**2.))
dr(i,j)=rad_cap(i,j)/1000.
cosinus(i,j)=(((x_cap(i,j)-x_cap_mean)*x_n_cap(i,j))+((y_cap(i,j)-y_cap_mean)*y_n_cap(i,j))+ &
((z_cap(i,j)-z_cap_mean)*z_n_cap(i,j)))/rad_cap(i,j)
enddo
enddo

!*********************************************************computing moment of inertia ********************************************
 mass=(4./3.)*pi



inxx=0.
inxy=0.
inxz=0.
inyy=0.
inyz=0.
inzz=0.
vol=0.
do l=1,1000
do i=1,n_capf
do j=1,m_capf
ma=(l-1)/1000.
na=(l-.5)/1000.
vol=vol+(area_cap(i,j)*dr(i,j)*(ma**2.)*abs(cosinus(i,j)))
inxx=inxx+(area_cap(i,j)*dr(i,j)*(ma**2.)*(na**2.)*abs(cosinus(i,j))*(((y_cap(i,j)-y_cap_mean)**2.)+ &
((z_cap(i,j)-z_cap_mean)**2.)))
inxy=inxy-(area_cap(i,j)*dr(i,j)*(ma**2.)*(na**2.)*abs(cosinus(i,j))*((x_cap(i,j)-x_cap_mean)*(y_cap(i,j)-y_cap_mean)))
inxz=inxz-(area_cap(i,j)*dr(i,j)*(ma**2.)*(na**2.)*abs(cosinus(i,j))*((x_cap(i,j)-x_cap_mean)*(z_cap(i,j)-z_cap_mean)))
inyy=inyy+(area_cap(i,j)*dr(i,j)*(ma**2.)*(na**2.)*abs(cosinus(i,j))*(((x_cap(i,j)-x_cap_mean)**2.)+ &
((z_cap(i,j)-z_cap_mean)**2.)))
inyz=inyz-(area_cap(i,j)*dr(i,j)*(ma**2.)*(na**2.)*abs(cosinus(i,j))*((y_cap(i,j)-y_cap_mean)*(z_cap(i,j)-z_cap_mean)))
inzz=inzz+(area_cap(i,j)*dr(i,j)*(ma**2.)*(na**2.)*abs(cosinus(i,j))*(((x_cap(i,j)-x_cap_mean)**2.)+ &
((y_cap(i,j)-y_cap_mean)**2.)))
enddo
enddo
enddo
i1=inxx+inyy+inzz
i2=(inxx*inyy)+(inyy*inzz)+(inxx*inzz)-(inxy**2.)-(inyz**2.)-(inxz**2.)
i3=(inxx*inyy*inzz)-(inxx*(inyz**2.))-(inyy*(inxz**2.))-(inzz*(inxy**2.))+(2.*inxy*inxz*inyz)
r=((2.*(i1**3.))-(9.*i1*i2)+(27.*i3))/54.
q=((3.*i2)-(i1**2.))/9.
teta=acos(r/sqrt(-(q**3.)))
ip1=(i1/3.)+(2.*sqrt(-q)*cos(teta/3.))
ip2=(i1/3.)+(2.*sqrt(-q)*cos((teta+(2.*pi))/3.))
ip3=(i1/3.)+(2.*sqrt(-q)*cos((teta+(4.*pi))/3.))
ic1=ip1/(mass/5.)
ic2=ip2/(mass/5.)
ic3=ip3/(mass/5.)
lp2=sqrt(.5*(ic1+ic3-ic2))
lp1=sqrt(ic3-(lp2**2.))
lp3=sqrt(ic1-(lp2**2.))
if (abs(ip1-inxx)<0.0001) then
exy1=abs(lp2-lp3)/(lp2+lp3)
elseif (abs(ip2-inxx)<0.0001) then
exy1=abs(lp1-lp3)/(lp1+lp3)
elseif (abs(ip3-inxx)<0.0001) then
exy1=abs(lp1-lp2)/(lp1+lp2)
endif


radmax=rad_cap(1,1)
radmin=rad_cap(1,1)
do i=1,n_capf
do j=1,m_capf
radmax=max(radmax,rad_cap(i,j))
radmin=min(radmin,rad_cap(i,j))
enddo
enddo
exy2=abs(radmax-radmin)/(radmax+radmin)
!****************************************************end of solving procedure************************************************************
endif
enddo

return
end subroutine surf_vol
end module mod_surf_vol
