      module mod_param
      implicit none
      integer imax,jmax,kmax,i1,j1,k1
      parameter(imax=32,jmax=32,kmax=32,i1=imax+1,j1=jmax+1,k1=kmax+1)

      real lx,ly,lz
      parameter(lx=2.,ly=2.,lz=2.)
      real dxi,dyi,dzi
      parameter(dxi=imax/lx,dyi=jmax/ly,dzi=kmax/lz)

      real radius,picon
      parameter(radius=0.5,picon=3.1415926535897932)

      real retraction
      parameter(retraction=0./dxi) !lfps retracted from interface into interior of obstacle
!     Number of Lagrangian force point over surface sphere (at r=R-0.5/dxi !)
      integer NL
      parameter(NL =nint((picon/3.)*(12.*(((radius-retraction       )*dxi)**2) + 1. )))
!     nr lfp's of second shell
      integer NL2
      parameter(NL2=nint((picon/3.)*(12.*(((radius-retraction-1./dxi)*dxi)**2) + 1. )))
!     nr lfp's of third shell
      integer NL3
      parameter(NL3=nint((picon/3.)*(12.*(((radius-retraction-2./dxi)*dxi)**2) + 1. )))
!     nr lfp's of fourth shell
      integer NL4
      parameter(NL4=nint((picon/3.)*(12.*(((radius-retraction-3./dxi)*dxi)**2) + 1. )))
!     nr lfp's of 1st-4th shell
      integer NLtot
      parameter(NLtot=NL+NL2+NL3+NL4)

      character*5 datadir
      parameter(datadir = 'data/')
      end module mod_param
