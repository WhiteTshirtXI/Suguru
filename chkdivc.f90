module mod_chkdivc
implicit none
private
public chkdivc
contains
!
subroutine chkdivc
use mpi
use mod_param
use mod_common
use mod_common_mpi
implicit none
real :: divc,divctot,divctot_all,divcmax(2),divcmax_all(2)
integer :: i,j,k,im,jm,km
!
divcmax(1) = 0.
divcmax(2) = 1.*myid
divctot = 0.
im = 0
jm = 0
km = 0
!
do k=1,kmaxc
  do j=1,jmaxc
    do i=1,imaxc
      divc = (wnewc(i,j,k)-wnewc(i,j,k-1))*dzic + &
            (vnewc(i,j,k)-vnewc(i,j-1,k))*dyic + &
            (unewc(i,j,k)-unewc(i-1,j,k))*dxic
      divctot = divctot+divc
      divc = abs(divc)
      if(divc.gt.divcmax(1)) then
        divcmax(1) = divc
        im = i
        jm = j
        km = k
      endif
    enddo
  enddo
enddo
!
call mpi_allreduce(divctot,divctot_all,1,mpi_real8,mpi_sum,MPI_COMM_WORLD,error)
call mpi_allreduce(divcmax,divcmax_all,1,mpi_2double_precision,mpi_maxloc,MPI_COMM_WORLD,error)
!
if (myid.eq.int(divcmax_all(2))) then
  write(6,111) zstartc(1)-1+im, zstartc(2)-1+jm,km
  write(6,222) divctot_all,divcmax_all(1),int(divcmax_all(2))
  111 format('Maximal divcergence for fine mesh at i = ',I5,' j = ', I5,' k = ',I5)
  222 format('divcergence for fine mesh: Tot = ',e13.6,' Max = ',e13.6,' Rank = ',I3)
endif
!
if (divctot_all .gt. 1.e-3) then
  if (myid .eq. divcmax_all(2)) then
    write(6,*) 'Fatal error: total divcergence for fine mesh > 1.e-3!'
    write(6,*) 'Program aborted...'
  endif
  call mpi_finalize(error)
  stop
endif
!
return
end subroutine chkdivc
!
end module mod_chkdivc

