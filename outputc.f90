module mod_outputc
use mod_common
use mod_param
use mod_common_mpi
use decomp_2d
use decomp_2d_io
use mod_post
implicit none
private
public post2dc,write2dplanec
contains
!

subroutine post2dc(istep)
implicit none
integer, intent(in) :: istep
real, dimension(imaxc,jmaxc,kmaxc) :: aux
real    :: u4(1:jmaxc,1:kmaxc)
real    :: v4(1:jmaxc,1:kmaxc)
real    :: w4(1:jmaxc,1:kmaxc)
real    :: p4(1:jmaxc,1:kmaxc)
real    :: viscc4(1:jmaxc,1:kmaxc)
real    :: forceyc4(1:jmaxc,1:kmaxc)
real    :: forcezc4(1:jmaxc,1:kmaxc)
real    :: indich4(1:jmaxc,1:kmaxc)

real    :: u4r(0:(jtotc+1),0:(kmaxc+1))
real    :: v4r(0:(jtotc+1),0:(kmaxc+1))
real    :: w4r(0:(jtotc+1),0:(kmaxc+1))
real    :: p4r(0:(jtotc+1),0:(kmaxc+1))
real    :: viscc4r(0:(jtotc+1),0:(kmaxc+1))
real    :: forceyc4r(0:(jtotc+1),0:(kmaxc+1))
real    :: forcezc4r(0:(jtotc+1),0:(kmaxc+1))
real    :: indich4r(0:(jtotc+1),0:(kmaxc+1))

integer :: l,q,rank,sub_rank,sub_coordsc(ndims)
integer :: npoints
integer :: i,j,k,im,jm,km,it,jt,iloc,jloc
 character*3 number
 character*7 filenumber
! character*7 filenumber
aux(:,:,:) = 0.5*(vnewc(1:imaxc,1+1:jmaxc+1,1:kmaxc) + vnewc(1:imaxc,1:jmaxc,1:kmaxc))

call write2dplanec(aux,1,itotc/2,'yzi_fine',istep)


 
      i = itotc/2
      q=-1
  114 q=q+1
      if ( (1+q*imaxc) .lt. i ) go to 114
      q=q-1
      iloc = i - q*imaxc
      sub_coordsc(1) = q
      sub_coordsc(2) = 0
      call MPI_CART_RANK(comm_cart,sub_coordsc,rank,error) !rank of process to which data has to be sent
      do l=1,dims(2)-1
        sub_coordsc(1) = q
        sub_coordsc(2) = l
        call MPI_CART_RANK(comm_cart,sub_coordsc,sub_rank,error)
        if (myid .eq. sub_rank) then
          do k=1,kmaxc
            do j=1,jmaxc
              u4(j,k)=0.5*(unewc(iloc,j,k)+unewc(-1+iloc,j,k))
              v4(j,k)=vnewc(iloc,j,k)
              w4(j,k)=wnewc(iloc,j,k)
              p4(j,k)=pnewc(iloc,j,k)
              viscc4(j,k)=viscc(iloc,j,k)
              forceyc4(j,k)=forceyc(iloc,j,k)
              forcezc4(j,k)=forcezc(iloc,j,k)
              indich4(j,k)=indich(iloc,j,k)
            enddo
          enddo
          call MPI_SSEND(u4,jmaxc*kmaxc,MPI_REAL8,rank,7,comm_cart,error)
          call MPI_SSEND(v4,jmaxc*kmaxc,MPI_REAL8,rank,8,comm_cart,error)
          call MPI_SSEND(w4,jmaxc*kmaxc,MPI_REAL8,rank,9,comm_cart,error)
          call MPI_SSEND(p4,jmaxc*kmaxc,MPI_REAL8,rank,10,comm_cart,error)
          call MPI_SSEND(viscc4,jmaxc*kmaxc,MPI_REAL8,rank,11,comm_cart,error)
          call MPI_SSEND(forceyc4,jmaxc*kmaxc,MPI_REAL8,rank,12,comm_cart,error)
          call MPI_SSEND(forcezc4,jmaxc*kmaxc,MPI_REAL8,rank,13,comm_cart,error)
          call MPI_SSEND(indich4,jmaxc*kmaxc,MPI_REAL8,rank,14,comm_cart,error)
        endif
        if (myid .eq. rank) then
          call MPI_RECV(u4,jmaxc*kmaxc,MPI_REAL8,sub_rank,7,  &
                        comm_cart,status,error)
          call MPI_RECV(v4,jmaxc*kmaxc,MPI_REAL8,sub_rank,8,  &
                        comm_cart,status,error)
          call MPI_RECV(w4,jmaxc*kmaxc,MPI_REAL8,sub_rank,9,  &
                        comm_cart,status,error)
          call MPI_RECV(p4,jmaxc*kmaxc,MPI_REAL8,sub_rank,10, &
                        comm_cart,status,error)
          call MPI_RECV(viscc4,jmaxc*kmaxc,MPI_REAL8,sub_rank,11, &
                        comm_cart,status,error)
          call MPI_RECV(forceyc4,jmaxc*kmaxc,MPI_REAL8,sub_rank,12, &
                        comm_cart,status,error)
          call MPI_RECV(forcezc4,jmaxc*kmaxc,MPI_REAL8,sub_rank,13, &
                        comm_cart,status,error)
          call MPI_RECV(indich4,jmaxc*kmaxc,MPI_REAL8,sub_rank,14, &
                        comm_cart,status,error)
          do k=1,kmaxc
            do j=1,jmaxc
              u4r(j+l*jmaxc,k)=u4(j,k)
              v4r(j+l*jmaxc,k)=v4(j,k)
              w4r(j+l*jmaxc,k)=w4(j,k)
              p4r(j+l*jmaxc,k)=p4(j,k)
              viscc4r(j+l*jmaxc,k)=viscc4(j,k)
              forceyc4r(j+l*jmaxc,k)=forceyc4(j,k)
              forcezc4r(j+l*jmaxc,k)=forcezc4(j,k)
              indich4r(j+l*jmaxc,k)=indich4(j,k)
            enddo
          enddo
        endif
      enddo
      if (myid .eq. rank) then
        do k=1,kmaxc
          do j=1,jmaxc
            u4r(j,k)=0.5*(unewc(iloc,j,k)+unewc(-1+iloc,j,k))
            v4r(j,k)=vnewc(iloc,j,k)
            w4r(j,k)=wnewc(iloc,j,k)
            p4r(j,k)=pnewc(iloc,j,k)
            viscc4r(j,k)=viscc(iloc,j,k)
            forceyc4r(j,k)=forceyc(iloc,j,k)
            forcezc4r(j,k)=forcezc(iloc,j,k)
            indich4r(j,k)=indich(iloc,j,k)
          enddo
          !periodic b.c.
          u4r(0,k)     =u4r(jtotc,k)
        enddo
        !b.c.
        do j=1,jtotc
          w4r(j,0)=0.
        enddo

         

          write(filenumber,'(i7.7)') istep
          open(433,file=datadir//'yz_mid_fine'//filenumber//'.txt')
          write(433,*) 'VARIABLES = "y","z","u","v","w","p","visc","forcey","forcez","indicator"'
          write(433,*) 'ZONE T="Zone1"',' I=',jtotc+2*npoints,' J=',kmaxc,', F=POINT'
!'
          do k=1,kmaxc
            km=k-1
            do jt=-npoints+1,jtotc+npoints
              j=jt
              if (j .lt. 1)    j = j+jtotc
              if (j .gt. jtotc) j = j-jtotc
              jm=jt-1
              if (jm .lt. 1)    jm = jm+jtotc
              if (jm .gt. jtotc) jm = jm-jtotc
              write(433,'(10E15.7)') ((jt-0.5)/dyic),((k-0.5)/dzic),u4r(j,k), &
                                   0.5*(v4r(j,k)+v4r(jm,k)),                  &
                                   0.5*(w4r(j,k)+w4r(j,km)),p4r(j,k),&
                                   viscc4r(j,k),forceyc4r(j,k),forcezc4r(j,k),indich4r(j,k)                                      
                                  
            enddo
          enddo
          close(433)


        

      endif

!
return
end subroutine post2dc

subroutine write2dplanec(var,norm,islice,name,istep)
implicit none
real, intent(in), dimension(imaxc,jmaxc,kmaxc) :: var
integer, intent(in) :: norm,islice,istep
 character(len=3) :: name
 character(len=4) :: slicechar
 character(len=7) :: fldnum

write(fldnum,'(i7.7)') istep
write(slicechar,'(i4.4)') islice
select case(norm)
case(1) !normal to x --> yz plane
  call decomp_2d_write_plane(3,var,norm,islice,datadir//name//'_xeq_'//slicechar//'_fld_'//fldnum//'.out')
case(2) !normal to y --> zx plane
  call decomp_2d_write_plane(3,var,norm,islice,datadir//name//'_yeq_'//slicechar//'_fld_'//fldnum//'.out')
case(3) !normal to z --> xy plane
  call decomp_2d_write_plane(3,var,norm,islice,datadir//name//'_zeq_'//slicechar//'_fld_'//fldnum//'.out')
end select

return
end subroutine write2dplanec
!
end module mod_outputc

