module mod_interpolation
use mpi
use mod_common
use mod_common_mpi
use mod_param
implicit none
private
public interpolation
contains
subroutine interpolation
implicit none
integer::i,j,k,im,jm
real::xp,yp,x1,x2,y1,y2,q11,q12,q21,q22
real:: coeff11,coeff12,coeff21,coeff22
do i=0,i1c
 do j=0,j1c
 xp=(i*dxc)+(.5*(1.-(1./nr))*dx)
 yp=(j*dyc)+(.5*(1.-(1./nr))*dy)
 x1=floor(xp/dx)*dx
 x2=(floor(xp/dx)+1.)*dx
 y1=floor(yp/dy)*dy
 y2=(floor(yp/dy)+1.)*dy
 coeff11=((x2-xp)*(y2-yp))/((x2-x1)*(y2-y1))
 coeff12=((x2-xp)*(yp-y1))/((x2-x1)*(y2-y1))
 coeff21=((xp-x1)*(y2-yp))/((x2-x1)*(y2-y1))
 coeff22=((xp-x1)*(yp-y1))/((x2-x1)*(y2-y1))

 im=floor(xp/dx)
 jm=floor(yp/dy)
 q11=.5*(unew(im,jm,nf)+unew(im,jm,nf+1))
 q12=.5*(unew(im,jm+1,nf)+unew(im,jm+1,nf+1))
 q21=.5*(unew(im+1,jm,nf)+unew(im+1,jm,nf+1)) 
 q22=.5*(unew(im+1,jm+1,nf)+unew(im+1,jm+1,nf+1))

 ub(i,j)=(coeff11*q11)+(coeff12*q12)+(coeff21*q21)+(coeff22*q22)

 q11=.5*(vnew(im,jm,nf)+vnew(im,jm,nf+1))
 q12=.5*(vnew(im,jm+1,nf)+vnew(im,jm+1,nf+1))
 q21=.5*(vnew(im+1,jm,nf)+vnew(im+1,jm,nf+1)) 
 q22=.5*(vnew(im+1,jm+1,nf)+vnew(im+1,jm+1,nf+1))

 vb(i,j)=(coeff11*q11)+(coeff12*q12)+(coeff21*q21)+(coeff22*q22)

 q11= wnew(im,jm,nf+1)
 q12= wnew(im,jm+1,nf+1)
 q21= wnew(im+1,jm,nf+1)
 q22= wnew(im+1,jm+1,nf+1)

 wb(i,j)=(coeff11*q11)+(coeff12*q12)+(coeff21*q21)+(coeff22*q22)

 q11=.5*(pnew(im,jm,nf)+pnew(im,jm,nf+1))
 q12=.5*(pnew(im,jm+1,nf)+pnew(im,jm+1,nf+1))
 q21=.5*(pnew(im+1,jm,nf)+pnew(im+1,jm,nf+1)) 
 q22=.5*(pnew(im+1,jm+1,nf)+pnew(im+1,jm+1,nf+1))

 pb(i,j)=(coeff11*q11)+(coeff12*q12)+(coeff21*q21)+(coeff22*q22)


 q11=.5*(indic(im,jm,nf)+indic(im,jm,nf+1))
 q12=.5*(indic(im,jm+1,nf)+indic(im,jm+1,nf+1))
 q21=.5*(indic(im+1,jm,nf)+indic(im+1,jm,nf+1)) 
 q22=.5*(indic(im+1,jm+1,nf)+indic(im+1,jm+1,nf+1))
 
 indb(i,j)=(coeff11*q11)+(coeff12*q12)+(coeff21*q21)+(coeff22*q22)
 
 enddo
enddo
return
end subroutine interpolation
end module mod_interpolation
